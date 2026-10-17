import math

import pytest

from limerocket import thermo
from limerocket.propellants import CATALOG, Ingredient, Propellant, Role
from limerocket.thermo import Freeze, Reactants
from limerocket.units import G0_M_S2, M_PER_FT, psi_to_pa
from tests.conftest import LIQUID_HYDROGEN, LIQUID_OXYGEN

LOX_LH2 = Reactants.of(LIQUID_OXYGEN, LIQUID_HYDROGEN)


# NASA RP-1311 Part II example 8, CEA2 reference output (nasa/cea test/main_interface).
CEA2_EX8_PC = 53.3172e5
CEA2_EX8_OF = 5.55157


def test_lox_lh2_chamber_and_c_star_match_cea2_example_8():
    p = thermo.solve(LOX_LH2, CEA2_EX8_OF, CEA2_EX8_PC, area_ratio=25.0)
    assert p.chamber.temperature_k == pytest.approx(3383.8446, rel=1e-5)
    assert p.throat.temperature_k == pytest.approx(3185.6732, rel=1e-5)
    assert p.c_star_m_s == pytest.approx(2332.3356, rel=1e-5)


@pytest.mark.parametrize(
    ("area_ratio", "ivac_m_s", "cf"),
    [(25.0, 4348.5138, 1.76836), (50.0, 4487.3021, 1.84756)],
)
def test_lox_lh2_vacuum_isp_at_area_ratio_matches_cea2_example_8(area_ratio, ivac_m_s, cf):
    p = thermo.solve(LOX_LH2, CEA2_EX8_OF, CEA2_EX8_PC, area_ratio=area_ratio)
    assert p.isp_vacuum_s * G0_M_S2 == pytest.approx(ivac_m_s, rel=1e-5)
    assert p.cf(p.exit_pressure_pa) == pytest.approx(cf, rel=1e-5)


def test_lox_lh2_isp_at_pressure_ratio_1000_matches_cea2_example_8():
    p = thermo.solve(LOX_LH2, CEA2_EX8_OF, CEA2_EX8_PC, exit_pressure_pa=CEA2_EX8_PC / 1000)
    assert p.area_ratio == pytest.approx(68.7533, rel=1e-5)
    assert p.isp_s(p.exit_pressure_pa) * G0_M_S2 == pytest.approx(4380.8113, rel=1e-5)


def test_nto_mmh_frozen_at_throat_matches_cea2_example_12():
    # RP-1311 example 12: 1000 psia, O/F 2.5, frozen at the throat, restricted product list.
    mmh = Propellant("mmh", "MMH", Role.FUEL, (Ingredient("CH6N2(L)", 1.0),), 298.15, 874, 1e-3, "")
    nto = Propellant(
        "nto", "NTO", Role.OXIDIZER, (Ingredient("N2O4(L)", 1.0),), 298.15, 1431, 1e-3, ""
    )
    only = [
        "CO",
        "CO2",
        "H",
        "HNO",
        "HO2",
        "H2",
        "H2O",
        "H2O2",
        "N",
        "NO",
        "NO2",
        "N2",
        "N2O",
        "O",
        "OH",
        "O2",
        "HCO",
        "NH",
        "CH4",
        "NH2",
        "NH3",
        "H2O(L)",
        "C(gr)",
    ]
    pc = psi_to_pa(1000)
    p = thermo.solve(
        Reactants.of(nto, mmh),
        2.5,
        pc,
        exit_pressure_pa=pc / 68.0457,
        freeze=Freeze.THROAT,
        products=tuple(only),
    )
    assert p.c_star_m_s == pytest.approx(1707.9214, rel=1e-5)
    assert p.area_ratio == pytest.approx(8.34153, rel=1e-5)
    assert p.isp_s(p.exit_pressure_pa) * G0_M_S2 == pytest.approx(2679.357, rel=1e-5)
    assert p.isp_vacuum_s * G0_M_S2 == pytest.approx(2888.726, rel=1e-5)


def test_lox_lh2_matches_rocketcea_documented_example():
    # RocketCEA docs, standard examples: LOX/LH2, Pc 1000 psia, MR 6, eps 40 (CEA2 Fortran).
    p = thermo.solve(LOX_LH2, 6.0, psi_to_pa(1000), area_ratio=40.0)
    assert p.chamber.temperature_k == pytest.approx(3483.35, rel=1e-4)
    assert p.c_star_m_s / M_PER_FT == pytest.approx(7560.0, rel=1e-4)
    assert p.cf(p.exit_pressure_pa) == pytest.approx(1.8351, rel=1e-4)
    assert p.isp_vacuum_s == pytest.approx(451.7, rel=2e-4)
    assert p.isp_s(p.exit_pressure_pa) == pytest.approx(431.2, rel=2e-4)


def test_optimal_mixture_ratio_sits_at_the_isp_peak():
    lox, fuel = CATALOG["lox"], CATALOG["ethanol75"]
    pc, pa = 20e5, 101_325.0
    best = thermo.optimal_mixture_ratio(lox, fuel, pc, pa)

    def isp(of: float) -> float:
        return thermo.performance(lox, fuel, of, pc, exit_pressure_pa=pa).isp_s(pa)

    assert isp(best) > isp(0.97 * best)
    assert isp(best) > isp(1.03 * best)


def test_performance_map_reproduces_cea_between_grid_points():
    lox, fuel = CATALOG["lox"], CATALOG["rp1"]
    design = thermo.performance(lox, fuel, 2.3, 30e5, exit_pressure_pa=101_325.0)
    table = thermo.performance_map(lox, fuel, design.area_ratio, 2.3, 30e5)
    of, pc = 2.13, 17.3e5
    direct = thermo.performance(lox, fuel, of, pc, area_ratio=design.area_ratio)
    point = table(of, pc)
    assert point.c_star_m_s == pytest.approx(direct.c_star_m_s, rel=2e-3)
    assert point.cf_vacuum == pytest.approx(direct.cf_vacuum, rel=2e-3)


def test_shared_species_in_both_propellants_are_merged():
    reactants = Reactants.of(CATALOG["htp90"], CATALOG["ethanol75"])
    water = reactants.names.index("H2O(L)")
    assert reactants.names.count("H2O(L)") == 1
    assert reactants.oxidizer_weights[water] == pytest.approx(0.10)
    assert reactants.fuel_weights[water] == pytest.approx(0.25)
    p = thermo.performance(CATALOG["htp90"], CATALOG["ethanol75"], 3.5, 20e5, exit_pressure_pa=1e5)
    assert math.isfinite(p.c_star_m_s)
