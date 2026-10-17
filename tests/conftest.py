from pathlib import Path

import pytest

from limerocket import design, spec
from limerocket.design import Design
from limerocket.propellants import Ingredient, Propellant, Role

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"

LIQUID_HYDROGEN = Propellant(
    key="lh2",
    name="Liquid hydrogen",
    role=Role.FUEL,
    ingredients=(Ingredient("H2(L)", 1.0),),
    temperature_k=20.27,
    density_kg_m3=70.8,
    viscosity_pa_s=1.3e-5,
    source="test",
)
LIQUID_OXYGEN = Propellant(
    key="lox",
    name="Liquid oxygen",
    role=Role.OXIDIZER,
    ingredients=(Ingredient("O2(L)", 1.0),),
    temperature_k=90.17,
    density_kg_m3=1141.0,
    viscosity_pa_s=1.96e-4,
    source="test",
)


@pytest.fixture(scope="session")
def lime1() -> Design:
    return design.design(spec.load(EXAMPLES / "lime-1.yaml"))


@pytest.fixture(scope="session")
def peroxide() -> Design:
    return design.design(spec.load(EXAMPLES / "peroxide-blowdown.yaml"))
