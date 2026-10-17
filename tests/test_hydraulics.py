import math

import numpy as np
import pytest

from limerocket.hydraulics import (
    Branch,
    Fluid,
    Line,
    Orifices,
    colebrook,
    darcy_friction_factor,
    haaland,
)


@pytest.mark.parametrize("reynolds", np.geomspace(1e4, 1e8, 9))
@pytest.mark.parametrize("roughness", [0.0, 1e-5, 1e-4, 1e-3, 1e-2, 5e-2])
def test_colebrook_agrees_with_haaland_within_its_stated_accuracy(reynolds, roughness):
    # Haaland (1983) quotes about 1.5 % from Colebrook-White.
    assert colebrook(reynolds, roughness) == pytest.approx(haaland(reynolds, roughness), rel=0.02)


def test_colebrook_matches_the_moody_chart():
    # Moody (1944): Re 1e5 with e/D 0.001 reads f = 0.022; fully rough e/D 0.01 reads 0.038.
    assert colebrook(1e5, 1e-3) == pytest.approx(0.0222, abs=3e-4)
    assert colebrook(1e8, 1e-2) == pytest.approx(0.0379, abs=3e-4)


@pytest.mark.parametrize("reynolds", [1e4, 3e4, 1e5])
def test_smooth_pipe_colebrook_agrees_with_blasius(reynolds):
    # Blasius (1913): f = 0.3164 Re^-0.25, valid for smooth pipes below Re ~ 1e5.
    assert colebrook(reynolds, 0.0) == pytest.approx(0.3164 * reynolds**-0.25, rel=0.03)


def test_laminar_flow_uses_64_over_re_and_transition_is_continuous():
    assert darcy_friction_factor(1000.0, 1e-3) == pytest.approx(0.064)
    below = darcy_friction_factor(2299.999, 0.0)
    above = darcy_friction_factor(2300.001, 0.0)
    assert above == pytest.approx(below, rel=1e-5)
    assert darcy_friction_factor(4000.0, 0.0) == pytest.approx(colebrook(4000.0, 0.0))


def test_orifice_flow_matches_hand_calculation():
    # 12 orifices of 1.2 mm, Cd 0.65, liquid oxygen 1141 kg/m3, 4 bar drop:
    # A = 12 * pi/4 * (1.2e-3)^2 = 1.35717e-5 m2
    # mdot = 0.65 * A * sqrt(2 * 1141 * 4e5) = 0.26652 kg/s
    orifices = Orifices(12, 1.2e-3, 0.65)
    lox = Fluid(1141.0, 1.96e-4)
    assert orifices.area_m2 == pytest.approx(1.35717e-5, rel=1e-5)
    assert orifices.mass_flow_kg_s(4e5, lox) == pytest.approx(0.26652, rel=1e-4)
    assert orifices.pressure_drop_pa(0.26652, lox) == pytest.approx(4e5, rel=2e-4)


def test_line_pressure_drop_matches_darcy_weisbach_by_hand():
    # Water-like fluid, 0.2 kg/s through 1 m of 8 mm tube: v = 3.979 m/s, Re = 31 831.
    water = Fluid(1000.0, 1e-3)
    line = Line(1.0, 8e-3, 1.5e-6, 2.0)
    v = 0.2 / (1000.0 * math.pi * 0.004**2)
    f = colebrook(4 * 0.2 / (math.pi * 8e-3 * 1e-3), 1.5e-6 / 8e-3)
    expected = (f * 1.0 / 8e-3 + 2.0) * 0.5 * 1000.0 * v**2
    assert line.pressure_drop_pa(0.2, water) == pytest.approx(expected, rel=1e-12)
    assert v == pytest.approx(3.979, rel=1e-3)


def test_branch_flow_inverts_branch_pressure_drop():
    branch = Branch(Fluid(800.0, 1.5e-3), Line(1.5, 6e-3, 1.5e-6, 5.0), Orifices(8, 1e-3, 0.7))
    dp = branch.pressure_drop_pa(0.15)
    assert branch.mass_flow_kg_s(dp) == pytest.approx(0.15, rel=1e-9)
    assert branch.mass_flow_kg_s(-1.0) == 0.0
