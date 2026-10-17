import math

import numpy as np
import pytest

from limerocket import chamber
from limerocket.materials import CATALOG

SS316 = CATALOG["ss316"]


def make(**overrides):
    args = {
        "throat_area_m2": math.pi * 0.01**2,
        "area_ratio": 4.0,
        "l_star_m": 1.2,
        "contraction_ratio": None,
        "convergence_half_angle_deg": 30.0,
        "chamber_pressure_pa": 20e5,
        "material": SS316,
        "yield_factor": 1.5,
        "ultimate_factor": 2.0,
        "min_wall_m": 1e-3,
    } | overrides
    return chamber.design_thrust_chamber(**args)


def test_volume_up_to_the_throat_equals_l_star_times_throat_area():
    c = make()
    upstream = c.x_m <= 0.0
    volume = np.trapezoid(math.pi * c.r_m[upstream] ** 2, c.x_m[upstream])
    assert volume == pytest.approx(c.l_star_m * c.throat_area_m2, rel=1e-3)


def test_contraction_ratio_defaults_to_the_humble_correlation():
    # Ac/At = 8.0 Dt^-0.6 + 1.25 with Dt = 2 cm.
    assert chamber.humble_contraction_ratio(0.02) == pytest.approx(8.0 * 2.0**-0.6 + 1.25)
    c = make()
    assert c.contraction_ratio == pytest.approx(chamber.humble_contraction_ratio(0.02))
    assert c.chamber_radius_m == pytest.approx(0.01 * math.sqrt(c.contraction_ratio))


def test_contour_is_continuous_from_injector_face_to_exit():
    c = make(contraction_ratio=5.0)
    assert np.all(np.diff(c.x_m) > 0)
    assert c.r_m.min() == pytest.approx(c.throat_radius_m)
    assert c.r_m[0] == pytest.approx(c.chamber_radius_m)
    assert c.r_m[-1] == pytest.approx(c.throat_radius_m * 2.0)
    steps = np.hypot(np.diff(c.x_m), np.diff(c.r_m))
    assert steps.max() < c.cylinder_length_m + 1e-12


def test_wall_is_sized_by_hoop_stress_above_the_minimum_gauge():
    c = make(chamber_pressure_pa=150e5, min_wall_m=1e-4)
    allowable = SS316.allowable_stress_pa(1.5, 2.0)
    assert c.wall_thickness_m == pytest.approx(150e5 * c.chamber_radius_m / allowable)


def test_injector_plate_matches_roark_hand_calculation():
    # sigma = 3 (3 + nu) q a^2 / (8 t^2) with a = 30 mm, q = 20 bar, sigma = 100 MPa, nu = 0.3
    t = chamber.plate_thickness_m(0.03, 20e5, 100e6, 0.3)
    assert t == pytest.approx(0.03 * math.sqrt(3 * 3.3 * 20e5 / (8 * 100e6)))
    assert t == pytest.approx(4.720e-3, rel=1e-3)


def test_too_small_l_star_is_rejected():
    with pytest.raises(chamber.GeometryError, match="L\\*"):
        make(l_star_m=0.01, contraction_ratio=12.0)
