"""Ideal-rocket relations (Sutton & Biblarz, 8th ed., ch. 3) and the Rao bell contour."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from scipy.optimize import brentq

from limerocket.units import R_UNIVERSAL_J_MOL_K

FloatArray = npt.NDArray[np.float64]


def area_ratio_from_mach(gamma: float, mach: float) -> float:
    """A/A* for isentropic flow (Sutton & Biblarz eq. 3-14)."""
    g = gamma
    term = (2.0 / (g + 1.0)) * (1.0 + 0.5 * (g - 1.0) * mach**2)
    return float(term ** ((g + 1.0) / (2.0 * (g - 1.0))) / mach)


def pressure_ratio_from_mach(gamma: float, mach: float) -> float:
    """p/p0 for isentropic flow (Sutton & Biblarz eq. 3-13)."""
    return float((1.0 + 0.5 * (gamma - 1.0) * mach**2) ** (-gamma / (gamma - 1.0)))


def mach_from_area_ratio(gamma: float, area_ratio: float, *, supersonic: bool = True) -> float:
    if area_ratio < 1.0:
        raise ValueError("area ratio must be >= 1")
    if area_ratio == 1.0:
        return 1.0
    lo, hi = (1.0, 100.0) if supersonic else (1e-9, 1.0)
    return float(brentq(lambda m: area_ratio_from_mach(gamma, m) - area_ratio, lo, hi, xtol=1e-14))


def thrust_coefficient(
    gamma: float, exit_to_chamber: float, ambient_to_chamber: float, area_ratio: float
) -> float:
    """Ideal thrust coefficient (Sutton & Biblarz eq. 3-30)."""
    g = gamma
    momentum = math.sqrt(
        2.0
        * g**2
        / (g - 1.0)
        * (2.0 / (g + 1.0)) ** ((g + 1.0) / (g - 1.0))
        * (1.0 - exit_to_chamber ** ((g - 1.0) / g))
    )
    return momentum + (exit_to_chamber - ambient_to_chamber) * area_ratio


def characteristic_velocity(gamma: float, molar_mass_kg_mol: float, chamber_temp_k: float) -> float:
    """Ideal c* (Sutton & Biblarz eq. 3-32)."""
    g = gamma
    gas_constant = R_UNIVERSAL_J_MOL_K / molar_mass_kg_mol
    return math.sqrt(g * gas_constant * chamber_temp_k) / (
        g * math.sqrt((2.0 / (g + 1.0)) ** ((g + 1.0) / (g - 1.0)))
    )


@dataclass(frozen=True)
class IdealGasCheck:
    """Ideal-rocket values from the chamber state at constant gamma, to set beside CEA's."""

    c_star_m_s: float
    exit_mach: float
    exit_pressure_ratio: float
    cf_vacuum: float


def ideal_gas_check(
    gamma: float, molar_mass_kg_mol: float, chamber_temp_k: float, area_ratio: float
) -> IdealGasCheck:
    mach = mach_from_area_ratio(gamma, area_ratio)
    pe_pc = pressure_ratio_from_mach(gamma, mach)
    return IdealGasCheck(
        c_star_m_s=characteristic_velocity(gamma, molar_mass_kg_mol, chamber_temp_k),
        exit_mach=mach,
        exit_pressure_ratio=pe_pc,
        cf_vacuum=thrust_coefficient(gamma, pe_pc, 0.0, area_ratio),
    )


def divergence_efficiency(exit_angle_deg: float) -> float:
    """Momentum lost to a non-axial exit flow, (1 + cos alpha) / 2 (Sutton & Biblarz eq. 3-34).

    Taken at the bell's exit angle theta_e, as is usual for parabolic bells.
    """
    return 0.5 * (1.0 + math.cos(math.radians(exit_angle_deg)))


# Rao's optimum-contour wall angles for an 80 % bell, digitized from NASA SP-8120 (1976)
# as reproduced by Newlands, "The Thrust Optimised Parabolic nozzle" (2017).
RAO_80_AREA_RATIOS = np.array([4.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 100.0])
RAO_80_THETA_N_DEG = np.array([21.5, 23.0, 26.3, 28.8, 30.0, 31.0, 31.5, 33.5])
RAO_80_THETA_E_DEG = np.array([14.0, 13.0, 11.0, 9.0, 8.5, 8.0, 7.5, 7.0])

BELL_FRACTION = 0.8
CONE_REFERENCE_HALF_ANGLE_DEG = 15.0
ENTRANT_ARC_RADIUS_RT = 1.5  # Rao throat upstream radius, in throat radii
EXIT_ARC_RADIUS_RT = 0.382  # Rao throat downstream radius, in throat radii


def rao_angles_deg(area_ratio: float) -> tuple[float, float, bool]:
    """(theta_n, theta_e, in_range) by log-linear interpolation of the 80 % bell chart."""
    log_eps = np.log(RAO_80_AREA_RATIOS)
    x = math.log(area_ratio)
    in_range = bool(log_eps[0] <= x <= log_eps[-1])
    theta_n = float(np.interp(x, log_eps, RAO_80_THETA_N_DEG))
    theta_e = float(np.interp(x, log_eps, RAO_80_THETA_E_DEG))
    return theta_n, theta_e, in_range


@dataclass(frozen=True)
class BellNozzle:
    throat_radius_m: float
    area_ratio: float
    theta_n_deg: float
    theta_e_deg: float
    angles_in_chart_range: bool
    x_m: FloatArray
    r_m: FloatArray

    @property
    def exit_radius_m(self) -> float:
        return self.throat_radius_m * math.sqrt(self.area_ratio)

    @property
    def length_m(self) -> float:
        return float(self.x_m[-1])


def bell_length_m(throat_radius_m: float, area_ratio: float) -> float:
    """80 % of the length of a 15 degree cone with the same area ratio (Newlands eq. 3)."""
    cone = (math.sqrt(area_ratio) - 1.0) * throat_radius_m
    return BELL_FRACTION * cone / math.tan(math.radians(CONE_REFERENCE_HALF_ANGLE_DEG))


def rao_bell(throat_radius_m: float, area_ratio: float, points: int = 72) -> BellNozzle:
    """Throat-downstream arc plus a quadratic Bezier bell; x = 0 at the throat."""
    if area_ratio <= 1.0:
        raise ValueError("area ratio must be > 1")
    rt = throat_radius_m
    theta_n_deg, theta_e_deg, in_range = rao_angles_deg(area_ratio)
    theta_n, theta_e = math.radians(theta_n_deg), math.radians(theta_e_deg)

    arc_points = max(points // 7, 8)
    angle = np.linspace(-math.pi / 2.0, theta_n - math.pi / 2.0, arc_points)
    radius = EXIT_ARC_RADIUS_RT * rt
    arc_x = radius * np.cos(angle)
    arc_r = radius * np.sin(angle) + radius + rt

    nx, ny = float(arc_x[-1]), float(arc_r[-1])
    ex, ey = bell_length_m(rt, area_ratio), rt * math.sqrt(area_ratio)
    m1, m2 = math.tan(theta_n), math.tan(theta_e)
    c1, c2 = ny - m1 * nx, ey - m2 * ex
    qx = (c2 - c1) / (m1 - m2)
    qy = (m1 * c2 - m2 * c1) / (m1 - m2)

    t = np.linspace(0.0, 1.0, points - arc_points + 1)[1:]
    bell_x = (1 - t) ** 2 * nx + 2 * (1 - t) * t * qx + t**2 * ex
    bell_r = (1 - t) ** 2 * ny + 2 * (1 - t) * t * qy + t**2 * ey

    return BellNozzle(
        throat_radius_m=rt,
        area_ratio=area_ratio,
        theta_n_deg=theta_n_deg,
        theta_e_deg=theta_e_deg,
        angles_in_chart_range=in_range,
        x_m=np.concatenate([arc_x, bell_x]),
        r_m=np.concatenate([arc_r, bell_r]),
    )
