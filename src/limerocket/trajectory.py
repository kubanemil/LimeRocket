"""Vertical flight to apogee: variable-mass point, USSA 1976 air, Cd(Mach), 1/r^2 gravity."""

from __future__ import annotations

import functools
import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from scipy.integrate import solve_ivp

from limerocket import atmosphere
from limerocket.feed import BurnResult

FloatArray = npt.NDArray[np.float64]

# Generic power-off drag of a slender fin-stabilised rocket: flat subsonic, transonic rise
# peaking near Mach 1.1, supersonic decay. Replace with OpenRocket/RASAero data for a real airframe.
GENERIC_CD_TABLE: tuple[tuple[float, float], ...] = (
    (0.0, 0.45),
    (0.6, 0.45),
    (0.8, 0.48),
    (0.95, 0.60),
    (1.1, 0.72),
    (1.3, 0.66),
    (1.6, 0.57),
    (2.0, 0.49),
    (3.0, 0.40),
    (5.0, 0.32),
)

COAST_TIME_LIMIT_S = 1_000.0
OUTPUT_POINTS_PER_PHASE = 300


@dataclass(frozen=True)
class Vehicle:
    liftoff_mass_kg: float
    reference_area_m2: float
    cd_table: tuple[tuple[float, float], ...]
    rail_length_m: float
    launch_altitude_m: float

    @functools.cached_property
    def _cd_columns(self) -> tuple[FloatArray, FloatArray]:
        table = np.array(self.cd_table, dtype=float)
        return table[:, 0], table[:, 1]

    def drag_coefficient(self, mach: float) -> float:
        machs, cds = self._cd_columns
        return float(np.interp(mach, machs, cds))


@dataclass(frozen=True)
class Flight:
    time_s: FloatArray
    altitude_m: FloatArray
    velocity_m_s: FloatArray
    acceleration_m_s2: FloatArray
    mass_kg: FloatArray
    thrust_n: FloatArray
    drag_n: FloatArray
    mach: FloatArray
    dynamic_pressure_pa: FloatArray
    burnout_time_s: float
    rail_exit_velocity_m_s: float
    liftoff_thrust_to_weight: float

    @property
    def apogee_m(self) -> float:
        return float(self.altitude_m.max())

    @property
    def apogee_time_s(self) -> float:
        return float(self.time_s[int(np.argmax(self.altitude_m))])

    def _at_burnout(self, values: FloatArray) -> float:
        return float(np.interp(self.burnout_time_s, self.time_s, values))

    @property
    def burnout_velocity_m_s(self) -> float:
        return self._at_burnout(self.velocity_m_s)

    @property
    def burnout_altitude_m(self) -> float:
        return self._at_burnout(self.altitude_m)

    @property
    def max_q_pa(self) -> float:
        return float(self.dynamic_pressure_pa.max())

    @property
    def max_q_time_s(self) -> float:
        return float(self.time_s[int(np.argmax(self.dynamic_pressure_pa))])

    @property
    def max_acceleration_m_s2(self) -> float:
        return float(self.acceleration_m_s2.max())

    @property
    def max_mach(self) -> float:
        return float(self.mach.max())


def fly(vehicle: Vehicle, burn: BurnResult) -> Flight:
    t_burn = burn.time_s
    f_vac = burn.thrust_vacuum_n
    mdot = burn.mass_flow_kg_s
    exit_area = burn.exit_area_m2
    t_end = burn.burn_time_s
    area = vehicle.reference_area_m2
    site = vehicle.launch_altitude_m

    def forces(t: float, h: float, v: float, burning: bool) -> tuple[float, float, float]:
        air = atmosphere.state(site + h)
        thrust = (
            max(float(np.interp(t, t_burn, f_vac)) - air.pressure_pa * exit_area, 0.0)
            if burning
            else 0.0
        )
        mach = abs(v) / air.speed_of_sound_m_s
        drag = 0.5 * air.density_kg_m3 * v * abs(v) * vehicle.drag_coefficient(mach) * area
        return thrust, drag, atmosphere.gravity_m_s2(site + h)

    def rhs(burning: bool) -> Callable[[float, FloatArray], FloatArray]:
        def f(t: float, y: FloatArray) -> FloatArray:
            h, v, m = y
            thrust, drag, g = forces(t, h, v, burning)
            a = (thrust - drag) / m - g
            if h <= 0.0 and v <= 0.0 and a < 0.0:
                a = 0.0
            dm = -float(np.interp(t, t_burn, mdot)) if burning else 0.0
            return np.array([v, a, dm])

        return f

    def rail_exit(_t: float, y: FloatArray) -> float:
        return float(y[0] - vehicle.rail_length_m)

    rail_exit.direction = 1  # type: ignore[attr-defined]

    def apogee(_t: float, y: FloatArray) -> float:
        return float(y[1])

    apogee.terminal = True  # type: ignore[attr-defined]
    apogee.direction = -1  # type: ignore[attr-defined]

    y0 = np.array([0.0, 0.0, vehicle.liftoff_mass_kg])
    powered = solve_ivp(
        rhs(burning=True),
        (0.0, t_end),
        y0,
        events=[rail_exit],
        dense_output=True,
        max_step=t_end / 200.0,
        rtol=1e-8,
        atol=1e-8,
    )
    coast = solve_ivp(
        rhs(burning=False),
        (t_end, t_end + COAST_TIME_LIMIT_S),
        powered.y[:, -1],
        events=[apogee, rail_exit],
        dense_output=True,
        rtol=1e-8,
        atol=1e-8,
    )

    t1 = np.linspace(0.0, t_end, OUTPUT_POINTS_PER_PHASE)
    t2 = np.linspace(t_end, float(coast.t[-1]), OUTPUT_POINTS_PER_PHASE)[1:]
    times = np.concatenate([t1, t2])
    assert powered.sol is not None and coast.sol is not None
    states = np.concatenate([powered.sol(t1), coast.sol(t2)], axis=1)
    burning = times <= t_end

    rows = [
        forces(float(t), float(h), float(v), bool(b))
        for t, h, v, b in zip(times, states[0], states[1], burning, strict=True)
    ]
    thrust = np.array([r[0] for r in rows])
    drag = np.array([r[1] for r in rows])
    gravity = np.array([r[2] for r in rows])
    accel = (thrust - drag) / states[2] - gravity
    accel[(states[0] <= 0.0) & (accel < 0.0)] = 0.0
    air = [atmosphere.state(site + float(h)) for h in states[0]]
    mach = np.abs(states[1]) / np.array([a.speed_of_sound_m_s for a in air])
    q = 0.5 * np.array([a.density_kg_m3 for a in air]) * states[1] ** 2

    rail_hits = np.concatenate([powered.t_events[0], coast.t_events[1]])
    if rail_hits.size:
        t_rail = float(rail_hits[0])
        rail_velocity = float(np.interp(t_rail, times, states[1]))
    else:
        rail_velocity = math.nan

    g0 = atmosphere.gravity_m_s2(site)
    return Flight(
        time_s=times,
        altitude_m=states[0],
        velocity_m_s=states[1],
        acceleration_m_s2=accel,
        mass_kg=states[2],
        thrust_n=thrust,
        drag_n=drag,
        mach=mach,
        dynamic_pressure_pa=q,
        burnout_time_s=t_end,
        rail_exit_velocity_m_s=rail_velocity,
        liftoff_thrust_to_weight=float(thrust[0]) / (vehicle.liftoff_mass_kg * g0),
    )
