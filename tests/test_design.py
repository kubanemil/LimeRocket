import dataclasses

import pytest

from limerocket import design, spec
from limerocket.design import DesignError
from limerocket.units import G0_M_S2
from tests.conftest import EXAMPLES


def test_design_point_delivers_the_target_thrust(lime1):
    burn = lime1.burn
    assert burn.thrust_n(lime1.ambient_pa)[0] == pytest.approx(lime1.thrust_n, rel=1e-5)
    assert burn.chamber_pressure_pa[0] == pytest.approx(20e5, rel=1e-6)
    assert burn.mixture_ratio[0] == pytest.approx(lime1.mixture_ratio, rel=1e-6)


def test_delivered_vacuum_isp_applies_both_efficiencies(lime1):
    e = lime1.spec.engine
    expected = lime1.ideal.isp_vacuum_s * e.c_star_efficiency * e.cf_efficiency
    assert lime1.isp_vacuum_s == pytest.approx(expected, rel=1e-12)
    assert lime1.mass_flow_kg_s * lime1.isp_s * G0_M_S2 == pytest.approx(lime1.thrust_n)


def test_mass_budget_adds_up_and_matches_the_flight(lime1):
    total = sum(m.mass_kg for m in lime1.masses)
    assert lime1.liftoff_mass_kg == pytest.approx(total)
    assert lime1.flight.mass_kg[0] == pytest.approx(total)
    final = lime1.flight.mass_kg[-1]
    assert total - final == pytest.approx(lime1.burn.expelled_mass_kg, rel=1e-4)


def test_tanks_hold_the_load_with_the_requested_ullage(lime1):
    for load in (lime1.oxidizer_load, lime1.fuel_load):
        liquid = load.loaded_mass_kg / load.branch.fluid.density_kg_m3
        assert 1 - liquid / load.volume_m3 == pytest.approx(lime1.spec.feed.ullage_fraction)
    for vessel in (lime1.oxidizer_tank, lime1.fuel_tank, lime1.bottle):
        assert vessel is not None
        assert vessel.margin_of_safety >= -1e-12
        assert vessel.outer_diameter_m <= lime1.spec.structure.max_diameter_m


def test_apogee_search_hits_the_target(peroxide):
    search = peroxide.apogee_search
    assert search is not None and search.reached
    assert peroxide.flight.apogee_m == pytest.approx(search.target_m, abs=5.0)
    assert peroxide.burn.burn_time_s == pytest.approx(search.burn_time_s, rel=5e-3)


def test_engine_wider_than_the_airframe_is_rejected():
    s = spec.load(EXAMPLES / "lime-1.yaml")
    narrow = dataclasses.replace(s, structure=dataclasses.replace(s.structure, max_diameter_m=0.05))
    with pytest.raises(DesignError, match="thrust chamber is"):
        design.design(narrow)


def test_bottle_below_the_dropout_pressure_is_rejected():
    s = spec.load(EXAMPLES / "lime-1.yaml")
    weak = dataclasses.replace(s, feed=dataclasses.replace(s.feed, bottle_pressure_pa=30e5))
    with pytest.raises(DesignError, match="bottle pressure"):
        design.design(weak)
