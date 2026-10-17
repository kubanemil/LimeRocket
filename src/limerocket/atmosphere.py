"""U.S. Standard Atmosphere 1976 (NOAA-S/T 76-1562), 0-86 km geometric altitude.

Above 86 km the top layer is extended isothermally; the 80-86 km molecular-weight correction
(under 0.04 %) is neglected.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

G0 = 9.80665
R_STAR = 8.31432e3  # J/(kmol K), the value the 1976 standard uses
M0 = 28.9644  # kg/kmol, sea-level mean molecular weight
EARTH_RADIUS_M = 6_356_766.0
GAMMA_AIR = 1.4
SUTHERLAND_BETA = 1.458e-6
SUTHERLAND_S_K = 110.4

T0_K = 288.15
P0_PA = 101_325.0

LAYER_BASES_M = (0.0, 11_000.0, 20_000.0, 32_000.0, 47_000.0, 51_000.0, 71_000.0, 84_852.0)
LAPSE_RATES_K_M = (-0.0065, 0.0, 0.001, 0.0028, 0.0, -0.0028, -0.002, 0.0)


def _layer_bases() -> tuple[tuple[float, ...], tuple[float, ...]]:
    temps, pressures = [T0_K], [P0_PA]
    for i in range(len(LAYER_BASES_M) - 1):
        dh = LAYER_BASES_M[i + 1] - LAYER_BASES_M[i]
        t, p, lapse = temps[-1], pressures[-1], LAPSE_RATES_K_M[i]
        t_next = t + lapse * dh
        if lapse == 0.0:
            p_next = p * math.exp(-G0 * M0 * dh / (R_STAR * t))
        else:
            p_next = p * (t / t_next) ** (G0 * M0 / (R_STAR * lapse))
        temps.append(t_next)
        pressures.append(p_next)
    return tuple(temps), tuple(pressures)


LAYER_TEMPS_K, LAYER_PRESSURES_PA = _layer_bases()


def geopotential_m(geometric_m: float) -> float:
    return EARTH_RADIUS_M * geometric_m / (EARTH_RADIUS_M + geometric_m)


@dataclass(frozen=True)
class AirState:
    temperature_k: float
    pressure_pa: float
    density_kg_m3: float
    speed_of_sound_m_s: float
    viscosity_pa_s: float


def temperature_pressure(geometric_m: float) -> tuple[float, float]:
    h = geopotential_m(max(geometric_m, -5_000.0))
    i = len(LAYER_BASES_M) - 1
    while i > 0 and h < LAYER_BASES_M[i]:
        i -= 1
    hb, tb, pb, lapse = (
        LAYER_BASES_M[i],
        LAYER_TEMPS_K[i],
        LAYER_PRESSURES_PA[i],
        LAPSE_RATES_K_M[i],
    )
    t = tb + lapse * (h - hb)
    if lapse == 0.0:
        p = pb * math.exp(-G0 * M0 * (h - hb) / (R_STAR * tb))
    else:
        p = pb * (tb / t) ** (G0 * M0 / (R_STAR * lapse))
    return t, p


def state(geometric_m: float) -> AirState:
    t, p = temperature_pressure(geometric_m)
    return AirState(
        temperature_k=t,
        pressure_pa=p,
        density_kg_m3=p * M0 / (R_STAR * t),
        speed_of_sound_m_s=math.sqrt(GAMMA_AIR * R_STAR * t / M0),
        viscosity_pa_s=SUTHERLAND_BETA * t**1.5 / (t + SUTHERLAND_S_K),
    )


def gravity_m_s2(geometric_m: float) -> float:
    return G0 * (EARTH_RADIUS_M / (EARTH_RADIUS_M + geometric_m)) ** 2
