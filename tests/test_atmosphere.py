import pytest

from limerocket import atmosphere

# U.S. Standard Atmosphere 1976 (NOAA-S/T 76-1562), Table I, geometric altitude.
USSA_1976 = [
    (0.0, 288.150, 101_325.0, 1.2250),
    (1_000.0, 281.651, 89_876.0, 1.1117),
    (5_000.0, 255.676, 54_048.0, 0.73643),
    (10_000.0, 223.252, 26_500.0, 0.41351),
    (20_000.0, 216.650, 5_529.3, 0.088910),
    (30_000.0, 226.509, 1_197.0, 0.018410),
    (50_000.0, 270.650, 79.779, 1.0269e-3),
    (80_000.0, 198.639, 1.0524, 1.8458e-5),
]


@pytest.mark.parametrize(("altitude", "temperature", "pressure", "density"), USSA_1976)
def test_matches_ussa_1976_table(altitude, temperature, pressure, density):
    air = atmosphere.state(altitude)
    assert air.temperature_k == pytest.approx(temperature, rel=5e-5)
    assert air.pressure_pa == pytest.approx(pressure, rel=2e-4)
    assert air.density_kg_m3 == pytest.approx(density, rel=2e-4)


def test_sea_level_speed_of_sound_and_viscosity_match_the_standard():
    air = atmosphere.state(0.0)
    assert air.speed_of_sound_m_s == pytest.approx(340.294, rel=1e-5)
    assert air.viscosity_pa_s == pytest.approx(1.7894e-5, rel=1e-4)


def test_gravity_falls_with_inverse_square_of_radius():
    r0 = atmosphere.EARTH_RADIUS_M
    assert atmosphere.gravity_m_s2(r0) == pytest.approx(atmosphere.G0 / 4)
