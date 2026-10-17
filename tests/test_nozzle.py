import math

import numpy as np
import pytest

from limerocket import nozzle
from limerocket.units import G0_M_S2

# NACA Report 1135 (1953), isentropic flow table for gamma = 1.4.
NACA_1135 = [(2.0, 1.6875, 0.12780), (3.0, 4.2346, 0.02722), (0.5, 1.3398, 0.84302)]


@pytest.mark.parametrize(("mach", "area_ratio", "pressure_ratio"), NACA_1135)
def test_isentropic_relations_match_naca_1135(mach, area_ratio, pressure_ratio):
    # The table prints five significant figures.
    assert nozzle.area_ratio_from_mach(1.4, mach) == pytest.approx(area_ratio, rel=2e-4)
    assert nozzle.pressure_ratio_from_mach(1.4, mach) == pytest.approx(pressure_ratio, rel=2e-4)
    supersonic = mach > 1
    assert nozzle.mach_from_area_ratio(1.4, area_ratio, supersonic=supersonic) == pytest.approx(
        mach, rel=1e-4
    )


def test_area_ratio_from_pressure_ratio_agrees_with_mach_relations():
    gamma, mach = 1.22, 3.1
    pe_pc = nozzle.pressure_ratio_from_mach(gamma, mach)
    assert nozzle.area_ratio_from_pressure_ratio(gamma, pe_pc) == pytest.approx(
        nozzle.area_ratio_from_mach(gamma, mach), rel=1e-12
    )


def test_ideal_engine_matches_proptools_worked_example():
    # proptools nozzle tutorial (after Sutton & Biblarz ch. 3): pc 10 MPa, pe = pa = 100 kPa,
    # gamma 1.2, M 20 g/mol, Tc 3000 K, throat diameter 0.1 m.
    gamma, molar_mass, tc, pc, pe = 1.2, 20e-3, 3000.0, 10e6, 100e3
    throat_area = math.pi * 0.05**2
    eps = nozzle.area_ratio_from_pressure_ratio(gamma, pe / pc)
    c_star = nozzle.characteristic_velocity(gamma, molar_mass, tc)
    cf = nozzle.thrust_coefficient(gamma, pe / pc, pe / pc, eps)
    ve = nozzle.exhaust_velocity(gamma, molar_mass, tc, pe / pc)
    assert eps == pytest.approx(11.9, abs=0.05)
    assert c_star == pytest.approx(1722, abs=0.5)
    assert ve == pytest.approx(2832, abs=0.5)
    assert c_star * cf / G0_M_S2 == pytest.approx(288.7, abs=0.1)
    assert cf * pc * throat_area == pytest.approx(129.2e3, abs=50)
    assert pc * throat_area / c_star == pytest.approx(45.6, abs=0.05)


def test_thrust_coefficient_approaches_the_infinite_expansion_limit():
    gamma = 1.2
    limit = math.sqrt(2 * gamma**2 / (gamma - 1) * (2 / (gamma + 1)) ** ((gamma + 1) / (gamma - 1)))
    assert nozzle.thrust_coefficient(gamma, 0.0, 0.0, 1.0) == pytest.approx(limit, rel=1e-12)
    assert nozzle.thrust_coefficient(gamma, 1e-4, 0.0, 1.0) < limit
    assert limit == pytest.approx(2.2466, abs=1e-4)


def test_rao_angles_match_newlands_worked_example_at_area_ratio_70():
    # Newlands, "The Thrust Optimised Parabolic nozzle": 80 % bell, eps 70 -> 33 deg / 7 deg.
    theta_n, theta_e, in_range = nozzle.rao_angles_deg(70.0)
    assert in_range
    assert theta_n == pytest.approx(33.0, abs=1.0)
    assert theta_e == pytest.approx(7.0, abs=1.0)


def test_bell_starts_at_the_throat_and_ends_at_the_exit_with_rao_angles():
    rt, eps = 0.01, 10.0
    bell = nozzle.rao_bell(rt, eps)
    assert (bell.x_m[0], bell.r_m[0]) == pytest.approx((0.0, rt))
    assert bell.r_m[-1] == pytest.approx(rt * math.sqrt(eps))
    assert bell.x_m[-1] == pytest.approx(nozzle.bell_length_m(rt, eps))
    assert np.all(np.diff(bell.x_m) > 0)
    exit_slope = (bell.r_m[-1] - bell.r_m[-2]) / (bell.x_m[-1] - bell.x_m[-2])
    assert math.degrees(math.atan(exit_slope)) == pytest.approx(bell.theta_e_deg, abs=0.5)
    slopes = np.degrees(np.arctan(np.diff(bell.r_m) / np.diff(bell.x_m)))
    assert slopes.max() == pytest.approx(bell.theta_n_deg, abs=0.5)
