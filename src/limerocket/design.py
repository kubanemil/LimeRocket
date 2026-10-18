"""The design pipeline: spec -> chemistry -> thrust chamber -> feed -> tanks -> burn -> flight."""

from __future__ import annotations

import dataclasses
import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from limerocket import atmosphere, thermo
from limerocket.chamber import ThrustChamber, design_thrust_chamber
from limerocket.feed import (
    BurnResult,
    Engine,
    FeedMode,
    Pressurisation,
    TankLoad,
    blowdown_ullage_m3,
    bottle_volume_m3,
    simulate_burn,
)
from limerocket.hydraulics import Branch, Fluid, Orifices
from limerocket.nozzle import IdealGasCheck, divergence_efficiency, ideal_gas_check, rao_angles_deg
from limerocket.spec import DesignSpec
from limerocket.tanks import Vessel, size_vessel
from limerocket.thermo import Performance
from limerocket.trajectory import Flight, Vehicle, fly
from limerocket.units import G0_M_S2, m_to_mm, mm_to_m, pa_to_bar

INTERNAL_FIT_FRACTION = 0.9  # tanks and bottle sit inside the airframe with 10 % radial clearance
SUMMERFIELD_SEPARATION_RATIO = 0.4  # pe/pa below which nozzle flow separates (Sutton ch. 3)
MIN_LIFTOFF_THRUST_TO_WEIGHT = 5.0  # amateur rocketry rule of thumb
MIN_RAIL_EXIT_VELOCITY_M_S = 15.0  # amateur rocketry rule of thumb for fin authority
MIN_ORIFICE_DIAMETER_M = mm_to_m(0.3)
MIN_INJECTOR_STIFFNESS = 0.10  # injector dp / Pc below which chugging becomes likely
MAX_BURN_TIME_FACTOR = 10.0


class DesignError(ValueError):
    pass


@dataclass(frozen=True)
class InjectorSide:
    orifices: Orifices
    mass_flow_kg_s: float
    pressure_drop_pa: float
    velocity_m_s: float


@dataclass(frozen=True)
class Injector:
    oxidizer: InjectorSide
    fuel: InjectorSide

    @property
    def momentum_ratio(self) -> float:
        ox, fuel = self.oxidizer, self.fuel
        return (ox.mass_flow_kg_s * ox.velocity_m_s) / (fuel.mass_flow_kg_s * fuel.velocity_m_s)


@dataclass(frozen=True)
class MassItem:
    name: str
    mass_kg: float
    group: str


@dataclass(frozen=True)
class ApogeeSearch:
    target_m: float
    reached: bool
    burn_time_s: float
    evaluations: int


@dataclass(frozen=True)
class Design:
    spec: DesignSpec
    mixture_ratio: float
    ideal: Performance
    chamber: ThrustChamber
    engine: Engine
    mass_flow_kg_s: float
    c_star_m_s: float
    cf_vacuum: float
    divergence_efficiency: float
    ideal_gas: IdealGasCheck
    injector: Injector
    oxidizer_line_drop_pa: float
    fuel_line_drop_pa: float
    pressurisation: Pressurisation
    oxidizer_load: TankLoad
    fuel_load: TankLoad
    oxidizer_tank: Vessel
    fuel_tank: Vessel
    bottle: Vessel | None
    pressurant_mass_kg: float
    masses: tuple[MassItem, ...]
    vehicle_length_m: float
    burn: BurnResult
    flight: Flight
    warnings: tuple[str, ...]
    apogee_search: ApogeeSearch | None = None

    @property
    def thrust_n(self) -> float:
        return self.spec.engine.thrust_n

    @property
    def ambient_pa(self) -> float:
        return self.spec.engine.design_ambient_pressure_pa

    @property
    def isp_s(self) -> float:
        return self.thrust_n / (self.mass_flow_kg_s * G0_M_S2)

    @property
    def isp_vacuum_s(self) -> float:
        return self.c_star_m_s * self.cf_vacuum / G0_M_S2

    @property
    def cf(self) -> float:
        return self.thrust_n / (self.spec.engine.chamber_pressure_pa * self.chamber.throat_area_m2)

    @property
    def oxidizer_flow_kg_s(self) -> float:
        return self.mass_flow_kg_s * self.mixture_ratio / (1.0 + self.mixture_ratio)

    @property
    def fuel_flow_kg_s(self) -> float:
        return self.mass_flow_kg_s / (1.0 + self.mixture_ratio)

    @property
    def dry_mass_kg(self) -> float:
        return sum(m.mass_kg for m in self.masses if m.group in ("structure", "other"))

    @property
    def liftoff_mass_kg(self) -> float:
        return sum(m.mass_kg for m in self.masses)

    @property
    def propellant_mass_kg(self) -> float:
        return self.oxidizer_load.loaded_mass_kg + self.fuel_load.loaded_mass_kg


@dataclass(frozen=True)
class _Feed:
    injector: Injector
    oxidizer_branch: Branch
    fuel_branch: Branch
    tank_pressure_pa: float
    oxidizer_line_drop_pa: float
    fuel_line_drop_pa: float


def _feed(spec: DesignSpec, oxidizer_flow: float, fuel_flow: float) -> _Feed:
    eng, feed = spec.engine, spec.feed
    pc = eng.chamber_pressure_pa
    ox_fluid = Fluid(eng.oxidizer.density_kg_m3, eng.oxidizer.viscosity_pa_s)
    fuel_fluid = Fluid(eng.fuel.density_kg_m3, eng.fuel.viscosity_pa_s)
    ox_line_dp = feed.oxidizer_line.pressure_drop_pa(oxidizer_flow, ox_fluid)
    fuel_line_dp = feed.fuel_line.pressure_drop_pa(fuel_flow, fuel_fluid)
    # One regulator (or one blowdown start pressure) serves both tanks, so the branch with the
    # smaller line loss takes a stiffer injector to meet the same tank pressure.
    tank_pressure = pc * (1.0 + eng.injector.pressure_drop_fraction) + max(ox_line_dp, fuel_line_dp)

    def side(flow: float, fluid: Fluid, line_dp: float, cd: float, count: int) -> InjectorSide:
        dp = tank_pressure - pc - line_dp
        orifices = Orifices.sized_for(flow, dp, fluid, count, cd)
        velocity = flow / (fluid.density_kg_m3 * orifices.area_m2)
        return InjectorSide(orifices, flow, dp, velocity)

    ox_side = side(
        oxidizer_flow,
        ox_fluid,
        ox_line_dp,
        eng.injector.oxidizer.discharge_coefficient,
        eng.injector.oxidizer.orifices,
    )
    fuel_side = side(
        fuel_flow,
        fuel_fluid,
        fuel_line_dp,
        eng.injector.fuel.discharge_coefficient,
        eng.injector.fuel.orifices,
    )
    return _Feed(
        injector=Injector(ox_side, fuel_side),
        oxidizer_branch=Branch(ox_fluid, feed.oxidizer_line, ox_side.orifices),
        fuel_branch=Branch(fuel_fluid, feed.fuel_line, fuel_side.orifices),
        tank_pressure_pa=tank_pressure,
        oxidizer_line_drop_pa=ox_line_dp,
        fuel_line_drop_pa=fuel_line_dp,
    )


def _loads(
    spec: DesignSpec, f: _Feed, oxidizer_usable: float, fuel_usable: float
) -> tuple[TankLoad, TankLoad]:
    feed = spec.feed

    def load(branch: Branch, usable: float) -> TankLoad:
        residual = feed.residual_fraction * usable
        liquid_volume = (usable + residual) / branch.fluid.density_kg_m3
        if feed.mode is FeedMode.REGULATED:
            volume = liquid_volume / (1.0 - feed.ullage_fraction)
        else:
            usable_volume = usable / branch.fluid.density_kg_m3
            volume = liquid_volume + blowdown_ullage_m3(
                usable_volume, feed.blowdown_ratio, feed.polytropic_exponent
            )
        return TankLoad(branch, volume, usable, residual)

    return load(f.oxidizer_branch, oxidizer_usable), load(f.fuel_branch, fuel_usable)


def _pressurisation(
    spec: DesignSpec, tank_pressure: float, oxidizer: TankLoad, fuel: TankLoad
) -> Pressurisation:
    feed = spec.feed
    base = Pressurisation(
        mode=feed.mode,
        gas=feed.pressurant,
        gas_temperature_k=feed.pressurant_temperature_k,
        tank_pressure_pa=tank_pressure,
        polytropic_exponent=feed.polytropic_exponent,
    )
    if feed.mode is FeedMode.BLOWDOWN:
        return base
    dropout = tank_pressure + feed.regulator_dropout_margin_pa
    if feed.bottle_pressure_pa <= dropout:
        raise DesignError(
            f"bottle pressure {pa_to_bar(feed.bottle_pressure_pa):.0f} bar must exceed the "
            f"regulator dropout pressure {pa_to_bar(dropout):.1f} bar"
        )
    expelled = sum(t.usable_mass_kg / t.branch.fluid.density_kg_m3 for t in (oxidizer, fuel))
    volume = feed.bottle_volume_m3 or (1.0 + feed.pressurant_margin) * bottle_volume_m3(
        feed.pressurant, tank_pressure, expelled, feed.bottle_pressure_pa, dropout
    )
    return dataclasses.replace(
        base,
        bottle_volume_m3=volume,
        bottle_pressure_pa=feed.bottle_pressure_pa,
        regulator_dropout_pa=dropout,
    )


def _pressurant_mass_kg(p: Pressurisation, oxidizer: TankLoad, fuel: TankLoad) -> float:
    r_t = p.gas.gas_constant_j_kg_k * p.gas_temperature_k
    ullage = oxidizer.ullage_m3(oxidizer.loaded_mass_kg) + fuel.ullage_m3(fuel.loaded_mass_kg)
    return (p.tank_pressure_pa * ullage + p.bottle_pressure_pa * p.bottle_volume_m3) / r_t


def _vessels(
    spec: DesignSpec, p: Pressurisation, oxidizer: TankLoad, fuel: TankLoad
) -> tuple[Vessel, Vessel, Vessel | None]:
    st = spec.structure
    common = {
        "yield_factor": st.yield_factor,
        "ultimate_factor": st.ultimate_factor,
        "min_wall_m": st.min_wall_m,
    }
    radius = INTERNAL_FIT_FRACTION * 0.5 * st.max_diameter_m
    ox_tank = size_vessel(
        "Oxidizer tank",
        volume_m3=oxidizer.volume_m3,
        max_outer_radius_m=radius,
        pressure_pa=p.tank_pressure_pa,
        material=st.tank_material,
        **common,
    )
    fuel_tank = size_vessel(
        "Fuel tank",
        volume_m3=fuel.volume_m3,
        max_outer_radius_m=radius,
        pressure_pa=p.tank_pressure_pa,
        material=st.tank_material,
        **common,
    )
    bottle = None
    if p.mode is FeedMode.REGULATED:
        bottle = size_vessel(
            "Pressurant bottle",
            volume_m3=p.bottle_volume_m3,
            max_outer_radius_m=radius,
            pressure_pa=p.bottle_pressure_pa,
            material=st.bottle_material,
            **common,
        )
    return ox_tank, fuel_tank, bottle


def _warnings(design: Design) -> list[str]:
    spec, burn, flight = design.spec, design.burn, design.flight
    out: list[str] = []
    nozzle = design.chamber.nozzle
    if not nozzle.angles_in_chart_range:
        out.append(
            f"Area ratio {nozzle.area_ratio:.2f} is outside Rao's 80 % bell chart (4-100); "
            "wall angles are clamped to the nearest chart value"
        )
    site_pa = atmosphere.state(spec.vehicle.launch_altitude_m).pressure_pa
    pe_ratio = float(np.min(burn.exit_pressure_pa)) / site_pa
    if pe_ratio < SUMMERFIELD_SEPARATION_RATIO:
        out.append(
            f"Exit pressure falls to {pe_ratio:.2f} of ambient during the burn; flow separation is "
            "likely (Summerfield criterion pe/pa < 0.4) and sea-level thrust is overestimated"
        )
    stiffness = _min_injector_stiffness(design)
    if stiffness < MIN_INJECTOR_STIFFNESS:
        out.append(
            f"Injector pressure drop falls to {stiffness:.0%} of Pc during the burn; "
            "low-frequency instability (chugging) becomes likely"
        )
    for label, side in (("oxidizer", design.injector.oxidizer), ("fuel", design.injector.fuel)):
        if side.orifices.diameter_m < MIN_ORIFICE_DIAMETER_M:
            out.append(
                f"{label} orifices are {m_to_mm(side.orifices.diameter_m):.2f} mm; fewer, larger "
                "orifices are easier to drill"
            )
    if flight.liftoff_thrust_to_weight < MIN_LIFTOFF_THRUST_TO_WEIGHT:
        out.append(
            f"Lift-off thrust-to-weight is {flight.liftoff_thrust_to_weight:.1f}; "
            f"{MIN_LIFTOFF_THRUST_TO_WEIGHT:.0f} or more is the usual minimum"
        )
    if not flight.rail_exit_velocity_m_s >= MIN_RAIL_EXIT_VELOCITY_M_S:
        out.append(
            f"Rail exit velocity is {flight.rail_exit_velocity_m_s:.1f} m/s; aim for at least "
            f"{MIN_RAIL_EXIT_VELOCITY_M_S:.0f} m/s (longer rail or more thrust)"
        )
    if burn.dropout_time_s is not None:
        out.append(
            f"Regulator drops out at {burn.dropout_time_s:.2f} s; the rest of the burn is blowdown"
        )
    leftover = (
        burn.oxidizer_mass_kg[-1] - design.oxidizer_load.residual_mass_kg,
        burn.fuel_mass_kg[-1] - design.fuel_load.residual_mass_kg,
    )
    if max(leftover) > 0.01 * design.propellant_mass_kg:
        which = "oxidizer" if leftover[0] > leftover[1] else "fuel"
        out.append(
            f"O/F drift leaves {max(leftover):.3f} kg of usable {which} at burnout; "
            "trim the load or the injector"
        )
    if design.apogee_search and not design.apogee_search.reached:
        out.append(
            f"Apogee target {design.apogee_search.target_m:.0f} m is out of reach; the burn time "
            f"was set to the best found ({design.apogee_search.burn_time_s:.2f} s)"
        )
    return out


def _min_injector_stiffness(design: Design) -> float:
    burn = design.burn
    ox = design.oxidizer_load.branch
    fuel = design.fuel_load.branch
    pc = burn.chamber_pressure_pa
    dp_ox = np.array([ox.orifices.pressure_drop_pa(m, ox.fluid) for m in burn.oxidizer_flow_kg_s])
    dp_fuel = np.array([fuel.orifices.pressure_drop_pa(m, fuel.fluid) for m in burn.fuel_flow_kg_s])
    return float(np.min(np.minimum(dp_ox, dp_fuel) / pc))


def _design_once(spec: DesignSpec, mixture_ratio: float, ideal: Performance) -> Design:
    eng, st, veh = spec.engine, spec.structure, spec.vehicle
    pc, pa = eng.chamber_pressure_pa, eng.design_ambient_pressure_pa

    divergence = divergence_efficiency(rao_angles_deg(ideal.area_ratio)[1])
    cf_efficiency = divergence * eng.cf_efficiency
    cf_vacuum = cf_efficiency * ideal.cf_vacuum
    c_star = eng.c_star_efficiency * ideal.c_star_m_s
    throat_area = eng.thrust_n / (cf_vacuum * pc - pa * ideal.area_ratio)
    mass_flow = pc * throat_area / c_star

    chamber = design_thrust_chamber(
        throat_area_m2=throat_area,
        area_ratio=ideal.area_ratio,
        l_star_m=eng.l_star_m,
        contraction_ratio=eng.contraction_ratio,
        convergence_half_angle_deg=eng.convergence_half_angle_deg,
        chamber_pressure_pa=pc,
        material=eng.wall_material,
        yield_factor=st.yield_factor,
        ultimate_factor=st.ultimate_factor,
        min_wall_m=st.min_wall_m,
    )
    if chamber.outer_diameter_m > st.max_diameter_m:
        raise DesignError(
            f"thrust chamber is {m_to_mm(chamber.outer_diameter_m):.0f} mm across but the airframe "
            f"allows {m_to_mm(st.max_diameter_m):.0f} mm; lower the contraction ratio or raise Pc"
        )

    oxidizer_flow = mass_flow * mixture_ratio / (1.0 + mixture_ratio)
    fuel_flow = mass_flow - oxidizer_flow
    f = _feed(spec, oxidizer_flow, fuel_flow)
    engine = Engine(
        throat_area_m2=throat_area,
        area_ratio=ideal.area_ratio,
        c_star_efficiency=eng.c_star_efficiency,
        cf_efficiency=cf_efficiency,
        performance=thermo.performance_map(
            eng.oxidizer, eng.fuel, ideal.area_ratio, mixture_ratio, pc
        ),
    )

    def run(
        ox_usable: float, fuel_usable: float
    ) -> tuple[TankLoad, TankLoad, Pressurisation, BurnResult]:
        ox_load, fuel_load = _loads(spec, f, ox_usable, fuel_usable)
        p = _pressurisation(spec, f.tank_pressure_pa, ox_load, fuel_load)
        burn = simulate_burn(engine, p, ox_load, fuel_load, MAX_BURN_TIME_FACTOR * eng.burn_time_s)
        return ox_load, fuel_load, p, burn

    ox_usable, fuel_usable = oxidizer_flow * eng.burn_time_s, fuel_flow * eng.burn_time_s
    ox_load, fuel_load, p, burn = run(ox_usable, fuel_usable)
    if spec.feed.mode is FeedMode.BLOWDOWN:
        # With ullage proportional to the load, burn time scales exactly with the load.
        scale = eng.burn_time_s / burn.burn_time_s
        ox_load, fuel_load, p, burn = run(ox_usable * scale, fuel_usable * scale)

    ox_tank, fuel_tank, bottle = _vessels(spec, p, ox_load, fuel_load)
    pressurant = _pressurant_mass_kg(p, ox_load, fuel_load)
    masses = [
        MassItem("Thrust chamber and injector plate", chamber.mass_kg, "structure"),
        MassItem(ox_tank.name, ox_tank.mass_kg, "structure"),
        MassItem(fuel_tank.name, fuel_tank.mass_kg, "structure"),
    ]
    if bottle:
        masses.append(MassItem(bottle.name, bottle.mass_kg, "structure"))
    masses += [
        MassItem(name.replace("_", " ").capitalize(), kg, "other")
        for name, kg in veh.other_masses_kg
    ]
    masses += [
        MassItem(f"Pressurant ({spec.feed.pressurant.name.lower()})", pressurant, "pressurant"),
        MassItem(f"Oxidizer ({eng.oxidizer.name})", ox_load.loaded_mass_kg, "propellant"),
        MassItem(f"Fuel ({eng.fuel.name})", fuel_load.loaded_mass_kg, "propellant"),
    ]
    liftoff = sum(m.mass_kg for m in masses)
    length = (
        chamber.length_m
        + ox_tank.length_m
        + fuel_tank.length_m
        + (bottle.length_m if bottle else 0.0)
        + veh.other_length_m
    )
    vehicle = Vehicle(
        liftoff_mass_kg=liftoff,
        reference_area_m2=math.pi * (0.5 * st.max_diameter_m) ** 2,
        cd_table=veh.cd_table,
        rail_length_m=veh.rail_length_m,
        launch_altitude_m=veh.launch_altitude_m,
    )
    flight = fly(vehicle, burn)
    return Design(
        spec=spec,
        mixture_ratio=mixture_ratio,
        ideal=ideal,
        chamber=chamber,
        engine=engine,
        mass_flow_kg_s=mass_flow,
        c_star_m_s=c_star,
        cf_vacuum=cf_vacuum,
        divergence_efficiency=divergence,
        ideal_gas=ideal_gas_check(
            ideal.chamber.gamma,
            ideal.chamber.molar_mass_kg_mol,
            ideal.chamber.temperature_k,
            ideal.area_ratio,
        ),
        injector=f.injector,
        oxidizer_line_drop_pa=f.oxidizer_line_drop_pa,
        fuel_line_drop_pa=f.fuel_line_drop_pa,
        pressurisation=p,
        oxidizer_load=ox_load,
        fuel_load=fuel_load,
        oxidizer_tank=ox_tank,
        fuel_tank=fuel_tank,
        bottle=bottle,
        pressurant_mass_kg=pressurant,
        masses=tuple(masses),
        vehicle_length_m=length,
        burn=burn,
        flight=flight,
        warnings=(),
    )


def _with_burn_time(spec: DesignSpec, burn_time_s: float) -> DesignSpec:
    return dataclasses.replace(
        spec, engine=dataclasses.replace(spec.engine, burn_time_s=burn_time_s)
    )


BURN_TIME_BOUNDS = (0.2, 5.0)  # search range as multiples of the spec burn time
BRACKET_STEP = 1.5


def _search_apogee(
    spec: DesignSpec, mixture_ratio: float, ideal: Performance, target: float
) -> tuple[Design, ApogeeSearch]:
    """Bracket the burn time that reaches `target`, then refine with Brent's method."""
    designs: dict[float, Design] = {}

    def trial(burn_time: float) -> Design:
        if burn_time not in designs:
            designs[burn_time] = _design_once(
                _with_burn_time(spec, burn_time), mixture_ratio, ideal
            )
        return designs[burn_time]

    def apogee(burn_time: float) -> float:
        return trial(burn_time).flight.apogee_m

    t0 = spec.engine.burn_time_s
    lo_limit, hi_limit = BURN_TIME_BOUNDS[0] * t0, BURN_TIME_BOUNDS[1] * t0
    step = BRACKET_STEP if apogee(t0) < target else 1.0 / BRACKET_STEP
    a, b = t0, t0
    reached = False
    while lo_limit <= b * step <= hi_limit:
        a, b = b, b * step
        if (apogee(b) - target) * (apogee(a) - target) <= 0.0:
            reached = True
            break
        if step > 1.0 and apogee(b) < apogee(a):
            break
    if reached:
        lo, hi = sorted((a, b))
        burn_time = float(brentq(lambda t: apogee(t) - target, lo, hi, xtol=1e-3 * t0))
    else:
        burn_time = max(designs, key=apogee) if step > 1.0 else min(designs)
    final = trial(burn_time)
    return final, ApogeeSearch(target, reached, burn_time, len(designs))


def design(spec: DesignSpec) -> Design:
    eng = spec.engine
    pc, pa = eng.chamber_pressure_pa, eng.design_ambient_pressure_pa
    mixture_ratio = eng.mixture_ratio or thermo.optimal_mixture_ratio(
        eng.oxidizer, eng.fuel, pc, pa
    )
    ideal = thermo.performance(eng.oxidizer, eng.fuel, mixture_ratio, pc, exit_pressure_pa=pa)
    search = None
    if spec.vehicle.apogee_target_m is not None:
        result, search = _search_apogee(spec, mixture_ratio, ideal, spec.vehicle.apogee_target_m)
        result = dataclasses.replace(result, apogee_search=search)
    else:
        result = _design_once(spec, mixture_ratio, ideal)
    return dataclasses.replace(result, warnings=tuple(_warnings(result)))
