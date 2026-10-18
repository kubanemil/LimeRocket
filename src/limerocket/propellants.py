"""Propellant catalog: CEA composition plus the liquid properties the feed system needs."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

FRACTION_TOLERANCE = 1e-6


class Role(StrEnum):
    OXIDIZER = "oxidizer"
    FUEL = "fuel"


@dataclass(frozen=True)
class Ingredient:
    cea_name: str
    mass_fraction: float


@dataclass(frozen=True)
class Propellant:
    key: str
    name: str
    role: Role
    ingredients: tuple[Ingredient, ...]
    temperature_k: float
    density_kg_m3: float
    viscosity_pa_s: float
    source: str

    def __post_init__(self) -> None:
        total = sum(i.mass_fraction for i in self.ingredients)
        if abs(total - 1.0) > FRACTION_TOLERANCE:
            raise ValueError(f"{self.key}: ingredient mass fractions sum to {total}, not 1")


def _single(cea_name: str) -> tuple[Ingredient, ...]:
    return (Ingredient(cea_name, 1.0),)


T_ROOM_K = 298.15
T_LOX_NBP_K = 90.17  # CEA thermo.lib reference temperature for O2(L)
T_LCH4_NBP_K = 111.643  # CEA thermo.lib reference temperature for CH4(L)

CATALOG: dict[str, Propellant] = {
    p.key: p
    for p in (
        Propellant(
            key="lox",
            name="Liquid oxygen",
            role=Role.OXIDIZER,
            ingredients=_single("O2(L)"),
            temperature_k=T_LOX_NBP_K,
            density_kg_m3=1141.0,
            viscosity_pa_s=1.96e-4,
            source="NIST Chemistry WebBook, saturated liquid at 90.2 K",
        ),
        Propellant(
            key="htp90",
            name="90% high-test peroxide",
            role=Role.OXIDIZER,
            ingredients=(Ingredient("H2O2(L)", 0.90), Ingredient("H2O(L)", 0.10)),
            temperature_k=T_ROOM_K,
            density_kg_m3=1387.0,
            viscosity_pa_s=1.2e-3,
            source="Schumb, Satterfield & Wentworth, Hydrogen Peroxide (1955), 25 °C; "
            "viscosity approximate",
        ),
        Propellant(
            key="ethanol",
            name="Ethanol",
            role=Role.FUEL,
            ingredients=_single("C2H5OH(L)"),
            temperature_k=T_ROOM_K,
            density_kg_m3=785.1,
            viscosity_pa_s=1.074e-3,
            source="CRC Handbook of Chemistry and Physics, 25 °C",
        ),
        Propellant(
            key="ethanol75",
            name="75% ethanol, 25% water (by mass)",
            role=Role.FUEL,
            ingredients=(Ingredient("C2H5OH(L)", 0.75), Ingredient("H2O(L)", 0.25)),
            temperature_k=T_ROOM_K,
            density_kg_m3=852.0,
            viscosity_pa_s=2.0e-3,
            source="CRC Handbook aqueous ethanol tables, 25 °C; viscosity approximate",
        ),
        Propellant(
            key="ipa",
            name="Isopropyl alcohol",
            role=Role.FUEL,
            ingredients=_single("C3H8O,2propanol"),
            temperature_k=T_ROOM_K,
            density_kg_m3=781.3,
            viscosity_pa_s=2.04e-3,
            source="CRC Handbook of Chemistry and Physics, 25 °C",
        ),
        Propellant(
            key="rp1",
            name="RP-1 kerosene",
            role=Role.FUEL,
            ingredients=_single("RP-1"),
            temperature_k=T_ROOM_K,
            density_kg_m3=810.0,
            viscosity_pa_s=1.6e-3,
            source="MIL-DTL-25576 density range 799-815 kg/m3; viscosity approximate",
        ),
        Propellant(
            key="methane",
            name="Liquid methane",
            role=Role.FUEL,
            ingredients=_single("CH4(L)"),
            temperature_k=T_LCH4_NBP_K,
            density_kg_m3=422.4,
            viscosity_pa_s=1.17e-4,
            source="NIST Chemistry WebBook, saturated liquid at 111.7 K",
        ),
    )
}


def get(key: str, role: Role) -> Propellant:
    try:
        propellant = CATALOG[key]
    except KeyError:
        known = ", ".join(k for k, p in CATALOG.items() if p.role is role)
        raise KeyError(f"unknown {role}: {key!r} (known: {known})") from None
    if propellant.role is not role:
        raise KeyError(f"{key!r} is a {propellant.role}, not a {role}")
    return propellant
