"""Equilibrium rocket performance from NASA CEA (the `cea` wheel, Gordon & McBride RP-1311)."""

from __future__ import annotations

import functools
from dataclasses import dataclass
from enum import IntEnum

import cea
import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import minimize_scalar

from limerocket.propellants import Propellant
from limerocket.units import G0_M_S2, KG_PER_G, bar_to_pa, pa_to_bar

BACKEND = f"NASA CEA {cea.__version__} (Python bindings)"

cea.set_log_level(cea.LOG_NONE)


class ThermoError(RuntimeError):
    pass


class Freeze(IntEnum):
    """Station at which composition freezes; NONE is shifting equilibrium."""

    NONE = 0
    CHAMBER = 1
    THROAT = 2


@dataclass(frozen=True)
class Reactants:
    names: tuple[str, ...]
    oxidizer_weights: tuple[float, ...]
    fuel_weights: tuple[float, ...]
    temperatures_k: tuple[float, ...]

    @classmethod
    def of(cls, oxidizer: Propellant, fuel: Propellant) -> Reactants:
        rows: dict[str, list[float]] = {}
        for propellant, column in ((oxidizer, 0), (fuel, 1)):
            for ingredient in propellant.ingredients:
                row = rows.setdefault(ingredient.cea_name, [0.0, 0.0, propellant.temperature_k])
                if row[2] != propellant.temperature_k:
                    raise ThermoError(
                        f"{ingredient.cea_name} appears in both propellants at different "
                        "temperatures; CEA needs one temperature per reactant"
                    )
                row[column] += ingredient.mass_fraction
        names = tuple(rows)
        return cls(
            names=names,
            oxidizer_weights=tuple(rows[n][0] for n in names),
            fuel_weights=tuple(rows[n][1] for n in names),
            temperatures_k=tuple(rows[n][2] for n in names),
        )


@dataclass(frozen=True)
class FlowState:
    pressure_pa: float
    temperature_k: float
    molar_mass_kg_mol: float
    gamma: float
    density_kg_m3: float
    sonic_velocity_m_s: float
    mach: float


@dataclass(frozen=True)
class Performance:
    """Ideal (100 % efficient) performance at one O/F, chamber pressure and area ratio."""

    mixture_ratio: float
    chamber_pressure_pa: float
    area_ratio: float
    chamber: FlowState
    throat: FlowState
    exit: FlowState
    c_star_m_s: float
    cf_vacuum: float
    freeze: Freeze

    @property
    def exit_pressure_pa(self) -> float:
        return self.exit.pressure_pa

    @property
    def isp_vacuum_s(self) -> float:
        return self.c_star_m_s * self.cf_vacuum / G0_M_S2

    def cf(self, ambient_pa: float) -> float:
        return self.cf_vacuum - ambient_pa * self.area_ratio / self.chamber_pressure_pa

    def isp_s(self, ambient_pa: float) -> float:
        return self.c_star_m_s * self.cf(ambient_pa) / G0_M_S2


def unknown_species(names: list[str]) -> list[str]:
    missing = []
    for name in names:
        try:
            cea.Mixture([name])
        except Exception:  # the binding raises a bare Exception when a name is unknown
            missing.append(name)
    return missing


@functools.lru_cache(maxsize=32)
def _solver(
    reactants: Reactants, products: tuple[str, ...] | None
) -> tuple[cea.Mixture, cea.RocketSolver]:
    try:
        mixture = cea.Mixture(list(reactants.names))
        product_mix = (
            cea.Mixture(list(products))
            if products
            else cea.Mixture(list(reactants.names), products_from_reactants=True)
        )
    except Exception as exc:
        raise ThermoError(f"CEA rejected the reactants {reactants.names}: {exc}") from exc
    return mixture, cea.RocketSolver(product_mix, reactants=mixture)


@functools.lru_cache(maxsize=16384)
def solve(
    reactants: Reactants,
    mixture_ratio: float,
    chamber_pressure_pa: float,
    *,
    exit_pressure_pa: float | None = None,
    area_ratio: float | None = None,
    freeze: Freeze = Freeze.NONE,
    products: tuple[str, ...] | None = None,
) -> Performance:
    """Infinite-area-combustor rocket problem; give exactly one of exit pressure or area ratio."""
    if (exit_pressure_pa is None) == (area_ratio is None):
        raise ValueError("give exactly one of exit_pressure_pa or area_ratio")
    mixture, solver = _solver(reactants, products)
    weights = mixture.of_ratio_to_weights(
        np.array(reactants.oxidizer_weights), np.array(reactants.fuel_weights), mixture_ratio
    )
    enthalpy = (
        mixture.calc_property(cea.ENTHALPY, weights, np.array(reactants.temperatures_k)) / cea.R
    )
    pc_bar = pa_to_bar(chamber_pressure_pa)
    station: dict[str, object] = {"hc": enthalpy, "iac": True}
    if exit_pressure_pa is not None:
        station["pi_p"] = [chamber_pressure_pa / exit_pressure_pa]
    else:
        station["supar"] = [area_ratio]
    if freeze is not Freeze.NONE:
        station["n_frz"] = int(freeze)
    solution = cea.RocketSolution(solver)
    try:
        solver.solve(solution, weights, pc_bar, **station)
    except Exception as exc:
        raise ThermoError(f"CEA failed at O/F={mixture_ratio:.3f}, Pc={pc_bar:.2f} bar") from exc
    if not solution.converged:
        raise ThermoError(f"CEA did not converge at O/F={mixture_ratio:.3f}, Pc={pc_bar:.2f} bar")

    def state(i: int) -> FlowState:
        return FlowState(
            pressure_pa=bar_to_pa(float(solution.P[i])),
            temperature_k=float(solution.T[i]),
            molar_mass_kg_mol=float(solution.MW[i]) * KG_PER_G,
            gamma=float(solution.gamma_s[i]),
            density_kg_m3=float(solution.density[i]),
            sonic_velocity_m_s=float(solution.sonic_velocity[i]),
            mach=float(solution.Mach[i]),
        )

    c_star = float(solution.c_star[-1])
    return Performance(
        mixture_ratio=mixture_ratio,
        chamber_pressure_pa=chamber_pressure_pa,
        area_ratio=float(solution.ae_at[-1]),
        chamber=state(0),
        throat=state(1),
        exit=state(-1),
        c_star_m_s=c_star,
        cf_vacuum=float(solution.Isp_vacuum[-1]) / c_star,
        freeze=freeze,
    )


def performance(
    oxidizer: Propellant,
    fuel: Propellant,
    mixture_ratio: float,
    chamber_pressure_pa: float,
    *,
    exit_pressure_pa: float | None = None,
    area_ratio: float | None = None,
) -> Performance:
    return solve(
        Reactants.of(oxidizer, fuel),
        mixture_ratio,
        chamber_pressure_pa,
        exit_pressure_pa=exit_pressure_pa,
        area_ratio=area_ratio,
    )


MIXTURE_RATIO_SEARCH = (0.3, 12.0)


def optimal_mixture_ratio(
    oxidizer: Propellant, fuel: Propellant, chamber_pressure_pa: float, exit_pressure_pa: float
) -> float:
    """O/F that maximises Isp with the nozzle optimally expanded to `exit_pressure_pa`."""

    def isp(of: float) -> float:
        try:
            p = performance(
                oxidizer, fuel, of, chamber_pressure_pa, exit_pressure_pa=exit_pressure_pa
            )
        except ThermoError:
            return 0.0
        return p.isp_s(exit_pressure_pa)

    grid = np.geomspace(*MIXTURE_RATIO_SEARCH, 48)
    values = [isp(float(of)) for of in grid]
    best = int(np.argmax(values))
    lo, hi = grid[max(best - 1, 0)], grid[min(best + 1, len(grid) - 1)]
    result = minimize_scalar(
        lambda of: -isp(float(of)), bounds=(lo, hi), method="bounded", options={"xatol": 1e-4}
    )
    return float(result.x)


@dataclass(frozen=True)
class OperatingPoint:
    c_star_m_s: float
    cf_vacuum: float
    exit_pressure_ratio: float


class PerformanceMap:
    """Off-design performance at a fixed nozzle area ratio, interpolated over O/F and log Pc."""

    def __init__(
        self,
        oxidizer: Propellant,
        fuel: Propellant,
        area_ratio: float,
        mixture_ratios: np.ndarray,
        chamber_pressures_pa: np.ndarray,
    ) -> None:
        self.area_ratio = area_ratio
        self.mixture_ratios = mixture_ratios
        self.chamber_pressures_pa = chamber_pressures_pa
        shape = (mixture_ratios.size, chamber_pressures_pa.size)
        c_star, cf_vac, pe_ratio = np.empty(shape), np.empty(shape), np.empty(shape)
        for i, of in enumerate(mixture_ratios):
            for j, pc in enumerate(chamber_pressures_pa):
                p = performance(oxidizer, fuel, float(of), float(pc), area_ratio=area_ratio)
                c_star[i, j] = p.c_star_m_s
                cf_vac[i, j] = p.cf_vacuum
                pe_ratio[i, j] = p.exit_pressure_pa / p.chamber_pressure_pa
        log_pc = np.log(chamber_pressures_pa)
        self._c_star = RectBivariateSpline(mixture_ratios, log_pc, c_star, kx=1, ky=1)
        self._cf_vac = RectBivariateSpline(mixture_ratios, log_pc, cf_vac, kx=1, ky=1)
        self._pe_ratio = RectBivariateSpline(mixture_ratios, log_pc, pe_ratio, kx=1, ky=1)

    def __call__(self, mixture_ratio: float, chamber_pressure_pa: float) -> OperatingPoint:
        of = float(np.clip(mixture_ratio, self.mixture_ratios[0], self.mixture_ratios[-1]))
        pc = float(
            np.clip(
                chamber_pressure_pa, self.chamber_pressures_pa[0], self.chamber_pressures_pa[-1]
            )
        )
        x, y = np.array([of]), np.array([np.log(pc)])
        return OperatingPoint(
            c_star_m_s=float(self._c_star.ev(x, y)[0]),
            cf_vacuum=float(self._cf_vac.ev(x, y)[0]),
            exit_pressure_ratio=float(self._pe_ratio.ev(x, y)[0]),
        )


@functools.lru_cache(maxsize=16)
def performance_map(
    oxidizer: Propellant,
    fuel: Propellant,
    area_ratio: float,
    design_mixture_ratio: float,
    design_chamber_pressure_pa: float,
) -> PerformanceMap:
    mixture_ratios = design_mixture_ratio * np.linspace(0.5, 1.6, 23)
    chamber_pressures = design_chamber_pressure_pa * np.geomspace(0.1, 1.25, 12)
    return PerformanceMap(oxidizer, fuel, area_ratio, mixture_ratios, chamber_pressures)
