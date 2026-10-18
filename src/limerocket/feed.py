"""Pressure-fed operation over the burn: regulated or blowdown pressurisation, solved as an ODE.

Assumptions: ideal-gas pressurant, adiabatic tanks and bottle (no heat exchange with walls or
liquid), incompressible liquids, quasi-steady injector and line flow, no start/stop transients.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from enum import StrEnum

import numpy as np
import numpy.typing as npt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from limerocket.hydraulics import Branch
from limerocket.thermo import PerformanceMap
from limerocket.units import G0_M_S2, R_UNIVERSAL_J_MOL_K

FloatArray = npt.NDArray[np.float64]


class FeedError(RuntimeError):
    pass


class FeedMode(StrEnum):
    REGULATED = "regulated"
    BLOWDOWN = "blowdown"


@dataclass(frozen=True)
class Gas:
    key: str
    name: str
    molar_mass_kg_mol: float
    gamma: float

    @property
    def gas_constant_j_kg_k(self) -> float:
        return R_UNIVERSAL_J_MOL_K / self.molar_mass_kg_mol


GASES = {
    "helium": Gas("helium", "Helium", 4.002602e-3, 5.0 / 3.0),
    "nitrogen": Gas("nitrogen", "Nitrogen", 28.0134e-3, 1.40),
}


@dataclass(frozen=True)
class Engine:
    throat_area_m2: float
    area_ratio: float
    c_star_efficiency: float
    cf_efficiency: float
    performance: PerformanceMap

    @property
    def exit_area_m2(self) -> float:
        return self.throat_area_m2 * self.area_ratio


@dataclass(frozen=True)
class TankLoad:
    branch: Branch
    volume_m3: float
    usable_mass_kg: float
    residual_mass_kg: float

    @property
    def loaded_mass_kg(self) -> float:
        return self.usable_mass_kg + self.residual_mass_kg

    def ullage_m3(self, liquid_mass_kg: float) -> float:
        return self.volume_m3 - liquid_mass_kg / self.branch.fluid.density_kg_m3


@dataclass(frozen=True)
class Pressurisation:
    mode: FeedMode
    gas: Gas
    gas_temperature_k: float
    tank_pressure_pa: float
    polytropic_exponent: float
    bottle_volume_m3: float = 0.0
    bottle_pressure_pa: float = 0.0
    regulator_dropout_pa: float = 0.0


@dataclass(frozen=True)
class OperatingState:
    chamber_pressure_pa: float
    oxidizer_flow_kg_s: float
    fuel_flow_kg_s: float
    c_star_m_s: float
    thrust_vacuum_n: float
    exit_pressure_pa: float

    @property
    def mass_flow_kg_s(self) -> float:
        return self.oxidizer_flow_kg_s + self.fuel_flow_kg_s

    @property
    def mixture_ratio(self) -> float:
        return (
            self.oxidizer_flow_kg_s / self.fuel_flow_kg_s if self.fuel_flow_kg_s > 0 else math.inf
        )


def operate(
    engine: Engine, oxidizer: Branch, fuel: Branch, p_ox_pa: float, p_fuel_pa: float
) -> OperatingState:
    """Chamber pressure where injected flow matches nozzle flow: mdot * eta * c* = Pc * At."""

    def flows(pc: float) -> tuple[float, float]:
        return oxidizer.mass_flow_kg_s(p_ox_pa - pc), fuel.mass_flow_kg_s(p_fuel_pa - pc)

    def c_star(m_ox: float, m_fuel: float, pc: float) -> float:
        of = m_ox / m_fuel if m_fuel > 0 else 1e9
        return engine.c_star_efficiency * engine.performance(of, pc).c_star_m_s

    def residual(pc: float) -> float:
        m_ox, m_fuel = flows(pc)
        return (m_ox + m_fuel) * c_star(m_ox, m_fuel, pc) - pc * engine.throat_area_m2

    upper = min(p_ox_pa, p_fuel_pa)
    if upper <= 0.0 or residual(upper) > 0.0:
        raise FeedError("tank pressures cannot sustain flow from both tanks")
    pc = float(brentq(residual, 1e-9 * upper, upper, xtol=1e-9 * upper, rtol=1e-12))
    m_ox, m_fuel = flows(pc)
    of = m_ox / m_fuel if m_fuel > 0 else 1e9
    point = engine.performance(of, pc)
    return OperatingState(
        chamber_pressure_pa=pc,
        oxidizer_flow_kg_s=m_ox,
        fuel_flow_kg_s=m_fuel,
        c_star_m_s=engine.c_star_efficiency * point.c_star_m_s,
        thrust_vacuum_n=engine.cf_efficiency * point.cf_vacuum * pc * engine.throat_area_m2,
        exit_pressure_pa=point.exit_pressure_ratio * pc,
    )


def bottle_volume_m3(
    gas: Gas,
    tank_pressure_pa: float,
    expelled_volume_m3: float,
    bottle_pressure_pa: float,
    final_bottle_pressure_pa: float,
) -> float:
    """Adiabatic bottle sizing from an energy balance: Vb = gamma p dV / (pb0 - pbf)."""
    return (
        gas.gamma
        * tank_pressure_pa
        * expelled_volume_m3
        / (bottle_pressure_pa - final_bottle_pressure_pa)
    )


def blowdown_ullage_m3(usable_volume_m3: float, blowdown_ratio: float, exponent: float) -> float:
    """Initial ullage that makes p fall by `blowdown_ratio` over the burn, with p V^n = const."""
    return float(usable_volume_m3 / (blowdown_ratio ** (1.0 / exponent) - 1.0))


@dataclass(frozen=True)
class BurnResult:
    time_s: FloatArray
    chamber_pressure_pa: FloatArray
    oxidizer_flow_kg_s: FloatArray
    fuel_flow_kg_s: FloatArray
    thrust_vacuum_n: FloatArray
    exit_pressure_pa: FloatArray
    oxidizer_tank_pressure_pa: FloatArray
    fuel_tank_pressure_pa: FloatArray
    bottle_pressure_pa: FloatArray
    oxidizer_mass_kg: FloatArray
    fuel_mass_kg: FloatArray
    pressurant_mass_kg: FloatArray
    exit_area_m2: float
    dropout_time_s: float | None
    depleted: str

    @property
    def burn_time_s(self) -> float:
        return float(self.time_s[-1])

    @property
    def mass_flow_kg_s(self) -> FloatArray:
        return self.oxidizer_flow_kg_s + self.fuel_flow_kg_s

    @property
    def mixture_ratio(self) -> FloatArray:
        return self.oxidizer_flow_kg_s / self.fuel_flow_kg_s

    def thrust_n(self, ambient_pa: float) -> FloatArray:
        return self.thrust_vacuum_n - ambient_pa * self.exit_area_m2

    def isp_s(self, ambient_pa: float) -> FloatArray:
        return self.thrust_n(ambient_pa) / (self.mass_flow_kg_s * G0_M_S2)

    def total_impulse_ns(self, ambient_pa: float) -> float:
        return float(np.trapezoid(self.thrust_n(ambient_pa), self.time_s))

    @property
    def expelled_mass_kg(self) -> float:
        return float(np.trapezoid(self.mass_flow_kg_s, self.time_s))


SAMPLES_PER_PHASE = 200
ODE_RTOL = 1e-8
ODE_ATOL = 1e-10


@dataclass(frozen=True)
class _Phase:
    time: FloatArray
    m_ox: FloatArray
    m_fuel: FloatArray
    p_ox: FloatArray
    p_fuel: FloatArray
    p_bottle: FloatArray
    gas: FloatArray
    end: FloatArray
    depleted: str | None


def _depletion_event(tank: TankLoad, index: int) -> Callable[[float, FloatArray], float]:
    def event(_t: float, y: FloatArray) -> float:
        return float(y[index] - tank.residual_mass_kg)

    event.terminal = True  # type: ignore[attr-defined]
    event.direction = -1  # type: ignore[attr-defined]
    return event


def _which_depleted(t_events: list[FloatArray]) -> str | None:
    if t_events[0].size:
        return "oxidizer"
    if t_events[1].size:
        return "fuel"
    return None


def _bottle_state(p: Pressurisation) -> tuple[float, Callable[[float], tuple[float, float]]]:
    """Initial bottle gas mass and (pressure, temperature) as a function of remaining mass.

    The gas left in an adiabatic bottle expands isentropically, whatever leaves it.
    """
    gamma = p.gas.gamma
    t0 = p.gas_temperature_k
    mb0 = p.bottle_pressure_pa * p.bottle_volume_m3 / (p.gas.gas_constant_j_kg_k * t0)

    def state(mb: float) -> tuple[float, float]:
        ratio = max(mb, 0.0) / mb0
        return p.bottle_pressure_pa * ratio**gamma, t0 * ratio ** (gamma - 1.0)

    return mb0, state


def _regulated_phase(
    engine: Engine, p: Pressurisation, oxidizer: TankLoad, fuel: TankLoad, max_time_s: float
) -> _Phase:
    """Tanks held at the set pressure while the bottle stays above the dropout pressure."""
    r_gas = p.gas.gas_constant_j_kg_k
    t0 = p.gas_temperature_k
    p_set = p.tank_pressure_pa
    mb0, bottle = _bottle_state(p)
    op = operate(engine, oxidizer.branch, fuel.branch, p_set, p_set)
    ullage_rate_ox = op.oxidizer_flow_kg_s / oxidizer.branch.fluid.density_kg_m3
    ullage_rate_fuel = op.fuel_flow_kg_s / fuel.branch.fluid.density_kg_m3

    def rhs(_t: float, y: FloatArray) -> FloatArray:
        _, tb = bottle(y[2])
        # Throttling an ideal gas is isenthalpic, so gas enters the ullage at bottle temperature;
        # the energy balance at constant ullage pressure then gives mdot_gas = p dV/dt / (R Tb).
        gas_ox = p_set * ullage_rate_ox / (r_gas * tb)
        gas_fuel = p_set * ullage_rate_fuel / (r_gas * tb)
        return np.array(
            [-op.oxidizer_flow_kg_s, -op.fuel_flow_kg_s, -(gas_ox + gas_fuel), gas_ox, gas_fuel]
        )

    def dropout(_t: float, y: FloatArray) -> float:
        return bottle(y[2])[0] - p.regulator_dropout_pa

    dropout.terminal = True  # type: ignore[attr-defined]
    dropout.direction = -1  # type: ignore[attr-defined]

    y0 = np.array(
        [
            oxidizer.loaded_mass_kg,
            fuel.loaded_mass_kg,
            mb0,
            p_set * oxidizer.ullage_m3(oxidizer.loaded_mass_kg) / (r_gas * t0),
            p_set * fuel.ullage_m3(fuel.loaded_mass_kg) / (r_gas * t0),
        ]
    )
    sol = solve_ivp(
        rhs,
        (0.0, max_time_s),
        y0,
        events=[_depletion_event(oxidizer, 0), _depletion_event(fuel, 1), dropout],
        dense_output=True,
        rtol=ODE_RTOL,
        atol=ODE_ATOL,
    )
    times = np.linspace(0.0, float(sol.t[-1]), SAMPLES_PER_PHASE)
    assert sol.sol is not None
    y = sol.sol(times)
    return _Phase(
        time=times,
        m_ox=y[0],
        m_fuel=y[1],
        p_ox=np.full_like(times, p_set),
        p_fuel=np.full_like(times, p_set),
        p_bottle=np.array([bottle(mb)[0] for mb in y[2]]),
        gas=y[2] + y[3] + y[4],
        end=sol.y[:, -1],
        depleted=_which_depleted(sol.t_events),
    )


def _dropout_phase(
    engine: Engine,
    p: Pressurisation,
    oxidizer: TankLoad,
    fuel: TankLoad,
    *,
    start: FloatArray,
    t_start: float,
    max_time_s: float,
) -> _Phase:
    """Regulator wide open: tank pressure trails the bottle by the regulator's minimum drop.

    Gas keeps flowing from the bottle, so tank pressure falls smoothly from the set point and
    never rises above it. Holding p_tank = p_bottle - dp in the adiabatic energy balances of
    the bottle and the shared ullage gives mdot_gas = p_u dVu/dt / (R Tb (1 + Vu / Vb)).
    """
    r_gas = p.gas.gas_constant_j_kg_k
    vb = p.bottle_volume_m3
    margin = p.regulator_dropout_pa - p.tank_pressure_pa
    _, bottle = _bottle_state(p)
    rho_ox = oxidizer.branch.fluid.density_kg_m3
    rho_fuel = fuel.branch.fluid.density_kg_m3

    def tank_pressure(mb: float) -> float:
        return bottle(mb)[0] - margin

    def rhs(_t: float, y: FloatArray) -> FloatArray:
        m_ox, m_fuel, mb, _ = y
        pb, tb = bottle(mb)
        pu = pb - margin
        op = operate(engine, oxidizer.branch, fuel.branch, pu, pu)
        vu = oxidizer.ullage_m3(m_ox) + fuel.ullage_m3(m_fuel)
        dvu = op.oxidizer_flow_kg_s / rho_ox + op.fuel_flow_kg_s / rho_fuel
        gas = pu * dvu / (r_gas * tb * (1.0 + vu / vb))
        return np.array([-op.oxidizer_flow_kg_s, -op.fuel_flow_kg_s, -gas, gas])

    m_ox, m_fuel, mb, m_ull_ox, m_ull_fuel = start
    sol = solve_ivp(
        rhs,
        (t_start, max_time_s),
        np.array([m_ox, m_fuel, mb, m_ull_ox + m_ull_fuel]),
        events=[_depletion_event(oxidizer, 0), _depletion_event(fuel, 1)],
        dense_output=True,
        rtol=ODE_RTOL,
        atol=ODE_ATOL,
    )
    times = np.linspace(t_start, float(sol.t[-1]), SAMPLES_PER_PHASE)
    assert sol.sol is not None
    y = sol.sol(times)
    pu = np.array([tank_pressure(m) for m in y[2]])
    return _Phase(
        time=times,
        m_ox=y[0],
        m_fuel=y[1],
        p_ox=pu,
        p_fuel=pu,
        p_bottle=pu + margin,
        gas=y[2] + y[3],
        end=sol.y[:, -1],
        depleted=_which_depleted(sol.t_events),
    )


def _blowdown_phase(
    engine: Engine, p: Pressurisation, oxidizer: TankLoad, fuel: TankLoad, max_time_s: float
) -> _Phase:
    """Each tank's ullage expands polytropically, p V^n = const, as its liquid leaves."""
    n = p.polytropic_exponent
    p0 = p.tank_pressure_pa
    v_ox0 = oxidizer.ullage_m3(oxidizer.loaded_mass_kg)
    v_fuel0 = fuel.ullage_m3(fuel.loaded_mass_kg)
    gas_mass = p0 * (v_ox0 + v_fuel0) / (p.gas.gas_constant_j_kg_k * p.gas_temperature_k)

    def pressures(y: FloatArray) -> tuple[float, float]:
        return (
            p0 * (v_ox0 / oxidizer.ullage_m3(y[0])) ** n,
            p0 * (v_fuel0 / fuel.ullage_m3(y[1])) ** n,
        )

    def rhs(_t: float, y: FloatArray) -> FloatArray:
        op = operate(engine, oxidizer.branch, fuel.branch, *pressures(y))
        return np.array([-op.oxidizer_flow_kg_s, -op.fuel_flow_kg_s])

    sol = solve_ivp(
        rhs,
        (0.0, max_time_s),
        np.array([oxidizer.loaded_mass_kg, fuel.loaded_mass_kg]),
        events=[_depletion_event(oxidizer, 0), _depletion_event(fuel, 1)],
        dense_output=True,
        rtol=ODE_RTOL,
        atol=ODE_ATOL,
    )
    times = np.linspace(0.0, float(sol.t[-1]), SAMPLES_PER_PHASE)
    assert sol.sol is not None
    y = sol.sol(times)
    p_pairs = np.array([pressures(y[:, k]) for k in range(times.size)])
    return _Phase(
        time=times,
        m_ox=y[0],
        m_fuel=y[1],
        p_ox=p_pairs[:, 0],
        p_fuel=p_pairs[:, 1],
        p_bottle=np.full_like(times, np.nan),
        gas=np.full_like(times, gas_mass),
        end=sol.y[:, -1],
        depleted=_which_depleted(sol.t_events),
    )


def simulate_burn(
    engine: Engine,
    pressurisation: Pressurisation,
    oxidizer: TankLoad,
    fuel: TankLoad,
    max_time_s: float,
) -> BurnResult:
    """Integrate the burn from ignition until a tank empties or `max_time_s` passes."""
    p = pressurisation
    dropout_time: float | None = None
    if p.mode is FeedMode.REGULATED:
        regulated = _regulated_phase(engine, p, oxidizer, fuel, max_time_s)
        phases = [regulated]
        if regulated.depleted is None and regulated.time[-1] < max_time_s:
            dropout_time = float(regulated.time[-1])
            phases.append(
                _dropout_phase(
                    engine,
                    p,
                    oxidizer,
                    fuel,
                    start=regulated.end,
                    t_start=dropout_time,
                    max_time_s=max_time_s,
                )
            )
    else:
        phases = [_blowdown_phase(engine, p, oxidizer, fuel, max_time_s)]

    fields = ("time", "m_ox", "m_fuel", "p_ox", "p_fuel", "p_bottle", "gas")
    data = {
        f: np.concatenate([getattr(ph, f)[1 if i else 0 :] for i, ph in enumerate(phases)])
        for f in fields
    }
    states = [
        operate(engine, oxidizer.branch, fuel.branch, p_ox, p_fuel)
        for p_ox, p_fuel in zip(data["p_ox"], data["p_fuel"], strict=True)
    ]
    return BurnResult(
        time_s=data["time"],
        chamber_pressure_pa=np.array([s.chamber_pressure_pa for s in states]),
        oxidizer_flow_kg_s=np.array([s.oxidizer_flow_kg_s for s in states]),
        fuel_flow_kg_s=np.array([s.fuel_flow_kg_s for s in states]),
        thrust_vacuum_n=np.array([s.thrust_vacuum_n for s in states]),
        exit_pressure_pa=np.array([s.exit_pressure_pa for s in states]),
        oxidizer_tank_pressure_pa=data["p_ox"],
        fuel_tank_pressure_pa=data["p_fuel"],
        bottle_pressure_pa=data["p_bottle"],
        oxidizer_mass_kg=data["m_ox"],
        fuel_mass_kg=data["m_fuel"],
        pressurant_mass_kg=data["gas"],
        exit_area_m2=engine.exit_area_m2,
        dropout_time_s=dropout_time,
        depleted=phases[-1].depleted or "none (time limit)",
    )
