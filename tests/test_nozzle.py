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


def test_ideal_engine_matches_proptools_worked_example():
    # proptools nozzle tutorial (after Sutton & Biblarz ch. 3): pc 10 MPa, gamma 1.2,
    # M 20 g/mol, Tc 3000 K, throat diameter 0.1 m, Ae/At 11.9 gives pe = pa = 100 kPa.
    gamma, molar_mass, tc, pc, eps = 1.2, 20e-3, 3000.0, 10e6, 11.9
    check = nozzle.ideal_gas_check(gamma, molar_mass, tc, eps)
    pe_pc = check.exit_pressure_ratio
    cf = nozzle.thrust_coefficient(gamma, pe_pc, pe_pc, eps)
    throat_area = math.pi * 0.05**2
    assert pe_pc * pc == pytest.approx(100e3, abs=500)
    assert check.c_star_m_s == pytest.approx(1722, abs=0.5)
    assert check.c_star_m_s * cf == pytest.approx(2832, abs=3)
    assert check.c_star_m_s * cf / G0_M_S2 == pytest.approx(288.7, abs=0.3)
    assert cf * pc * throat_area == pytest.approx(129.2e3, abs=150)
    assert pc * throat_area / check.c_star_m_s == pytest.approx(45.6, abs=0.05)
    assert check.cf_vacuum == pytest.approx(cf + pe_pc * eps)


def test_divergence_factor_matches_sutton_table_for_a_15_degree_cone():
    # Sutton & Biblarz table 3-3: lambda = 0.9830 at 15 deg and 0.9698 at 20 deg.
    assert nozzle.divergence_efficiency(15.0) == pytest.approx(0.9830, abs=1e-4)
    assert nozzle.divergence_efficiency(20.0) == pytest.approx(0.9698, abs=1e-4)


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
