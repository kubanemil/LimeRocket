import math

import numpy as np
import pytest

from limerocket import atmosphere
from limerocket.feed import BurnResult
from limerocket.trajectory import Vehicle, fly

G = 9.80665


def constant_burn(thrust: float, mdot: float, burn_time: float) -> BurnResult:
    t = np.linspace(0.0, burn_time, 200)
    ones = np.ones_like(t)
    return BurnResult(
        time_s=t,
        chamber_pressure_pa=ones,
        oxidizer_flow_kg_s=0.5 * mdot * ones,
        fuel_flow_kg_s=0.5 * mdot * ones,
        thrust_vacuum_n=thrust * ones,
        exit_pressure_pa=ones,
        oxidizer_tank_pressure_pa=ones,
        fuel_tank_pressure_pa=ones,
        bottle_pressure_pa=ones,
        oxidizer_mass_kg=ones,
        fuel_mass_kg=ones,
        pressurant_mass_kg=ones,
        exit_area_m2=0.0,
        dropout_time_s=None,
        depleted="fuel",
    )


def vehicle(mass: float, area: float = 0.0) -> Vehicle:
    return Vehicle(mass, area, ((0.0, 0.5), (5.0, 0.5)), 5.0, 0.0)


def test_vacuum_constant_thrust_flight_matches_the_analytic_solution(monkeypatch):
    # Tsiolkovsky with gravity loss: v = ve ln(m0/m) - g t,
    # h = ve (m0/mdot) (1 - u + u ln u) - g t^2 / 2 with u = m/m0, then coast h + v^2 / 2g.
    monkeypatch.setattr(atmosphere, "gravity_m_s2", lambda _h: G)
    thrust, mdot, tb, m0 = 2000.0, 1.0, 8.0, 20.0
    flight = fly(vehicle(m0), constant_burn(thrust, mdot, tb))
    ve, u = thrust / mdot, (m0 - mdot * tb) / m0
    v_bo = ve * math.log(1 / u) - G * tb
    h_bo = ve * (m0 / mdot) * (1 - u + u * math.log(u)) - 0.5 * G * tb**2
    assert flight.burnout_velocity_m_s == pytest.approx(v_bo, rel=1e-6)
    assert flight.burnout_altitude_m == pytest.approx(h_bo, rel=1e-6)
    assert flight.apogee_m == pytest.approx(h_bo + v_bo**2 / (2 * G), rel=1e-6)
    assert flight.apogee_time_s == pytest.approx(tb + v_bo / G, rel=1e-5)
    assert flight.max_acceleration_m_s2 == pytest.approx(thrust / (m0 - mdot * tb) - G, rel=1e-5)


def test_rocket_heavier_than_its_thrust_stays_on_the_pad():
    flight = fly(vehicle(100.0), constant_burn(500.0, 0.2, 3.0))
    assert flight.apogee_m == pytest.approx(0.0, abs=1e-9)
    assert flight.liftoff_thrust_to_weight < 1
    assert math.isnan(flight.rail_exit_velocity_m_s)


def test_drag_lowers_apogee_and_dynamic_pressure_peaks_near_burnout():
    burn = constant_burn(2000.0, 1.0, 8.0)
    clean = fly(vehicle(20.0), burn)
    draggy = fly(vehicle(20.0, area=0.02), burn)
    assert draggy.apogee_m < clean.apogee_m
    assert 0.5 * burn.burn_time_s < draggy.max_q_time_s <= burn.burn_time_s + 0.5
