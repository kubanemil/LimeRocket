import dataclasses

import numpy as np
import pytest

from limerocket.feed import bottle_volume_m3, operate, simulate_burn


def test_operating_point_reproduces_the_design_flows(lime1):
    p_tank = lime1.pressurisation.tank_pressure_pa
    op = operate(lime1.engine, lime1.oxidizer_load.branch, lime1.fuel_load.branch, p_tank, p_tank)
    assert op.chamber_pressure_pa == pytest.approx(lime1.spec.engine.chamber_pressure_pa, rel=1e-6)
    assert op.oxidizer_flow_kg_s == pytest.approx(lime1.oxidizer_flow_kg_s, rel=1e-6)
    assert op.fuel_flow_kg_s == pytest.approx(lime1.fuel_flow_kg_s, rel=1e-6)


def test_regulated_burn_holds_chamber_pressure_and_conserves_pressurant(lime1):
    burn = lime1.burn
    pc = lime1.spec.engine.chamber_pressure_pa
    assert np.allclose(burn.chamber_pressure_pa, pc, rtol=1e-6)
    assert np.allclose(burn.pressurant_mass_kg, burn.pressurant_mass_kg[0], rtol=1e-9)
    assert burn.burn_time_s == pytest.approx(lime1.spec.engine.burn_time_s, rel=1e-4)
    assert burn.dropout_time_s is None


@pytest.mark.parametrize("which", ["lime1", "peroxide"])
def test_expelled_propellant_equals_the_mass_that_left_the_tanks(which, request):
    d = request.getfixturevalue(which)
    burn = d.burn
    t = burn.time_s
    ox_out = np.trapezoid(burn.oxidizer_flow_kg_s, t)
    fuel_out = np.trapezoid(burn.fuel_flow_kg_s, t)
    assert ox_out == pytest.approx(burn.oxidizer_mass_kg[0] - burn.oxidizer_mass_kg[-1], rel=1e-4)
    assert fuel_out == pytest.approx(burn.fuel_mass_kg[0] - burn.fuel_mass_kg[-1], rel=1e-4)


def test_blowdown_ullage_follows_the_polytropic_law(peroxide):
    burn, p = peroxide.burn, peroxide.pressurisation
    n = p.polytropic_exponent
    for load, pressure, mass in (
        (peroxide.oxidizer_load, burn.oxidizer_tank_pressure_pa, burn.oxidizer_mass_kg),
        (peroxide.fuel_load, burn.fuel_tank_pressure_pa, burn.fuel_mass_kg),
    ):
        volume = np.array([load.ullage_m3(m) for m in mass])
        invariant = pressure * volume**n
        assert np.allclose(invariant, invariant[0], rtol=1e-10)
    assert np.allclose(burn.pressurant_mass_kg, burn.pressurant_mass_kg[0], rtol=1e-12)
    start = burn.oxidizer_tank_pressure_pa[0]
    end = p.tank_pressure_pa / peroxide.spec.feed.blowdown_ratio
    assert burn.oxidizer_tank_pressure_pa[-1] == pytest.approx(end, rel=0.02)
    assert start == pytest.approx(p.tank_pressure_pa)


def test_closed_form_bottle_size_lands_on_the_dropout_pressure_at_burnout(lime1):
    p = lime1.pressurisation
    ox, fuel = lime1.oxidizer_load, lime1.fuel_load
    expelled = sum(t.usable_mass_kg / t.branch.fluid.density_kg_m3 for t in (ox, fuel))
    exact = bottle_volume_m3(
        p.gas, p.tank_pressure_pa, expelled, p.bottle_pressure_pa, p.regulator_dropout_pa
    )
    sized = dataclasses.replace(p, bottle_volume_m3=1.001 * exact)
    burn = simulate_burn(lime1.engine, sized, ox, fuel, 60.0)
    assert burn.dropout_time_s is None
    assert burn.bottle_pressure_pa[-1] == pytest.approx(p.regulator_dropout_pa, rel=0.01)


def undersized(lime1, fraction=0.5):
    p = lime1.pressurisation
    return dataclasses.replace(p, bottle_volume_m3=fraction * p.bottle_volume_m3)


def test_undersized_bottle_drops_out_and_the_burn_finishes_in_blowdown(lime1):
    burn = simulate_burn(
        lime1.engine, undersized(lime1), lime1.oxidizer_load, lime1.fuel_load, 60.0
    )
    assert burn.dropout_time_s is not None
    assert burn.dropout_time_s < burn.burn_time_s
    after = burn.time_s > burn.dropout_time_s
    assert np.all(np.diff(burn.chamber_pressure_pa[after]) < 0)
    assert burn.burn_time_s > lime1.burn.burn_time_s
    assert np.allclose(burn.pressurant_mass_kg, burn.pressurant_mass_kg[0], rtol=1e-6)


def test_regulator_dropout_never_lifts_tank_pressure_above_the_set_point(lime1):
    # Regression: dropout used to mix bottle and ullage gas at once, a 1.2 bar step up.
    p_set = lime1.pressurisation.tank_pressure_pa
    burn = simulate_burn(
        lime1.engine, undersized(lime1, 0.3), lime1.oxidizer_load, lime1.fuel_load, 60.0
    )
    assert burn.dropout_time_s is not None
    tank = burn.oxidizer_tank_pressure_pa
    assert tank.max() <= p_set * (1 + 1e-12)
    assert burn.chamber_pressure_pa.max() <= lime1.spec.engine.chamber_pressure_pa * (1 + 1e-6)
    assert np.all(np.diff(tank) <= 1e-6 * p_set)
    margin = lime1.pressurisation.regulator_dropout_pa - p_set
    after = burn.time_s >= burn.dropout_time_s
    assert np.allclose(burn.bottle_pressure_pa[after] - tank[after], margin, rtol=1e-9)


def test_burn_stops_at_the_time_limit_even_after_dropout(lime1):
    # Regression: the post-dropout phase used to get the whole time limit again.
    burn = simulate_burn(lime1.engine, undersized(lime1), lime1.oxidizer_load, lime1.fuel_load, 4.0)
    assert burn.dropout_time_s is not None and burn.dropout_time_s < 4.0
    assert burn.burn_time_s == pytest.approx(4.0)
    assert burn.depleted.startswith("none")
