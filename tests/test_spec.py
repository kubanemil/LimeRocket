import copy

import pytest
import yaml

from limerocket import spec
from limerocket.feed import FeedMode
from tests.conftest import EXAMPLES


def load_example(name: str) -> dict:
    return yaml.safe_load((EXAMPLES / name).read_text())


@pytest.mark.parametrize("path", sorted(EXAMPLES.glob("*.yaml")), ids=lambda p: p.name)
def test_examples_parse(path):
    parsed = spec.load(path)
    assert parsed.engine.thrust_n > 0


def test_values_are_converted_to_si():
    s = spec.parse(load_example("lime-1.yaml"))
    assert s.engine.chamber_pressure_pa == pytest.approx(22e5)
    assert s.structure.max_diameter_m == pytest.approx(0.152)
    assert s.feed.oxidizer_line.inner_diameter_m == pytest.approx(7.7e-3)
    assert s.feed.oxidizer_line.roughness_m == pytest.approx(1.5e-6)
    assert s.engine.mixture_ratio is None
    assert s.feed.mode is FeedMode.REGULATED
    assert s.feed.polytropic_exponent == pytest.approx(5 / 3)


def test_every_problem_is_reported_at_once():
    data = load_example("lime-1.yaml")
    data["engine"]["thrust_n"] = "lots"
    data["engine"]["fuel"] = "unobtainium"
    data["engine"]["injector"]["oxidizer"]["cd"] = 0.7
    data["feed"]["blowdown_ratio"] = 2.0
    data["structure"]["tank_material"] = "cheese"
    data["vehicle"]["rail_length_m"] = -1
    data["colour"] = "lime"
    with pytest.raises(spec.SpecError) as err:
        spec.parse(data)
    text = str(err.value)
    for expected in (
        "engine.thrust_n: expected a number",
        "engine.fuel: unknown fuel: 'unobtainium'",
        "engine.injector.oxidizer.cd: unknown key",
        "feed.blowdown_ratio: applies only to mode: blowdown",
        "structure.tank_material: unknown value 'cheese'",
        "vehicle.rail_length_m: must be >= 0",
        "colour: unknown key",
    ):
        assert expected in text
    assert len(err.value.problems) == 7


def test_role_mismatch_is_rejected():
    data = load_example("lime-1.yaml")
    data["engine"]["oxidizer"] = "ethanol"
    with pytest.raises(spec.SpecError, match="'ethanol' is a fuel, not a oxidizer"):
        spec.parse(data)


def test_custom_propellant_from_cea_species():
    data = load_example("lime-1.yaml")
    data["engine"]["fuel"] = {
        "name": "E90",
        "ingredients": {"C2H5OH(L)": 0.9, "H2O(L)": 0.1},
        "temperature_k": 298.15,
        "density_kg_m3": 818,
        "viscosity_pa_s": 1.4e-3,
    }
    s = spec.parse(data)
    assert s.engine.fuel.name == "E90"
    assert [i.cea_name for i in s.engine.fuel.ingredients] == ["C2H5OH(L)", "H2O(L)"]


def test_custom_propellant_fractions_within_tolerance_are_normalised():
    # Regression: 3 x 0.3333333 passed the spec check, then Propellant raised a bare ValueError.
    data = load_example("lime-1.yaml")
    data["engine"]["fuel"] = {
        "name": "thirds",
        "ingredients": {"C2H5OH(L)": 0.3333333, "CH3OH(L)": 0.3333333, "H2O(L)": 0.3333333},
        "temperature_k": 298.15,
        "density_kg_m3": 850,
        "viscosity_pa_s": 1.4e-3,
    }
    fuel = spec.parse(data).engine.fuel
    assert sum(i.mass_fraction for i in fuel.ingredients) == pytest.approx(1.0, abs=1e-12)


def test_unknown_cea_species_is_a_spec_error():
    data = load_example("lime-1.yaml")
    data["engine"]["fuel"] = {
        "name": "mystery",
        "ingredients": {"Unobtainium(L)": 1.0},
        "temperature_k": 298.15,
        "density_kg_m3": 850,
        "viscosity_pa_s": 1.4e-3,
    }
    with pytest.raises(spec.SpecError, match="not in the CEA database: Unobtainium"):
        spec.parse(data)


@pytest.mark.parametrize(("pc_bar", "ambient_bar"), [(1.005, 1.01325), (1.5, 1.0), (3.0, 1.6)])
def test_chamber_pressure_must_choke_the_nozzle(pc_bar, ambient_bar):
    # Regression: Pc barely above ambient parsed fine and crashed the burn integrator.
    data = load_example("lime-1.yaml")
    data["engine"]["chamber_pressure_bar"] = pc_bar
    data["engine"]["design_ambient_pressure_bar"] = ambient_bar
    with pytest.raises(spec.SpecError, match="chamber_pressure_bar: must be at least 2"):
        spec.parse(data)


def test_custom_propellant_fractions_must_sum_to_one():
    data = load_example("lime-1.yaml")
    data["engine"]["fuel"] = {
        "name": "bad",
        "ingredients": {"C2H5OH(L)": 0.9},
        "temperature_k": 298.15,
        "density_kg_m3": 818,
        "viscosity_pa_s": 1.4e-3,
    }
    with pytest.raises(spec.SpecError, match="must sum to 1"):
        spec.parse(data)


def test_drag_table_must_increase_in_mach():
    data = load_example("lime-1.yaml")
    data["vehicle"]["drag_coefficient"] = [[0.0, 0.4], [0.0, 0.5]]
    with pytest.raises(spec.SpecError, match="increasing Mach"):
        spec.parse(data)
    data["vehicle"]["drag_coefficient"] = 0.5
    assert spec.parse(data).vehicle.cd_table == ((0.0, 0.5), (10.0, 0.5))


def test_fixed_mixture_ratio_and_auto_values():
    data = copy.deepcopy(load_example("lime-1.yaml"))
    data["engine"]["mixture_ratio"] = 1.2
    data["engine"]["contraction_ratio"] = 6
    s = spec.parse(data)
    assert s.engine.mixture_ratio == 1.2
    assert s.engine.contraction_ratio == 6
    del data["engine"]["mixture_ratio"]
    assert spec.parse(data).engine.mixture_ratio is None
    data["engine"]["mixture_ratio"] = "best"
    with pytest.raises(spec.SpecError, match="mixture_ratio: expected a number or 'auto'"):
        spec.parse(data)


def test_invalid_yaml_is_a_spec_error(tmp_path):
    path = tmp_path / "bad.yaml"
    path.write_text("engine: [unclosed")
    with pytest.raises(spec.SpecError, match="not valid YAML"):
        spec.load(path)
