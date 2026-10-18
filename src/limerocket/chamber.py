from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from limerocket.materials import Material
from limerocket.nozzle import ENTRANT_ARC_RADIUS_RT, BellNozzle, rao_bell

FloatArray = npt.NDArray[np.float64]

CHAMBER_FILLET_RADIUS_RC = 1.0  # cylinder-to-cone blend radius, in chamber radii (design choice)


class GeometryError(ValueError):
    pass


def humble_contraction_ratio(throat_diameter_m: float) -> float:
    """Ac/At = 8.0 Dt^-0.6 + 1.25, Dt in cm (Humble, Henry & Larson, Space Propulsion Analysis
    and Design, 1995)."""
    return float(8.0 * (throat_diameter_m * 100.0) ** -0.6 + 1.25)


@dataclass(frozen=True)
class ThrustChamber:
    throat_radius_m: float
    chamber_radius_m: float
    contraction_ratio: float
    l_star_m: float
    chamber_volume_m3: float
    cylinder_length_m: float
    convergent_length_m: float
    convergence_half_angle_deg: float
    nozzle: BellNozzle
    x_m: FloatArray
    r_m: FloatArray
    wall_thickness_m: float
    injector_plate_thickness_m: float
    shell_mass_kg: float
    injector_plate_mass_kg: float

    @property
    def throat_area_m2(self) -> float:
        return math.pi * self.throat_radius_m**2

    @property
    def exit_radius_m(self) -> float:
        return self.nozzle.exit_radius_m

    @property
    def length_m(self) -> float:
        return float(self.x_m[-1] - self.x_m[0])

    @property
    def outer_diameter_m(self) -> float:
        return 2.0 * (max(self.chamber_radius_m, self.exit_radius_m) + self.wall_thickness_m)

    @property
    def mass_kg(self) -> float:
        return self.shell_mass_kg + self.injector_plate_mass_kg


def _convergent(rt: float, rc: float, half_angle: float) -> tuple[FloatArray, FloatArray, float]:
    """Contour from the cylinder end to the throat (x = 0); returns x, r and the cylinder-end x."""
    r1 = ENTRANT_ARC_RADIUS_RT * rt
    p1x = -r1 * math.sin(half_angle)
    p1r = rt + r1 * (1.0 - math.cos(half_angle))
    if rc <= p1r:
        raise GeometryError(
            "contraction ratio too small for the throat entrant arc; raise it or the "
            "convergence angle"
        )
    fillet = min(CHAMBER_FILLET_RADIUS_RC * rc, (rc - p1r) / (1.0 - math.cos(half_angle)))
    e2r = rc - fillet * (1.0 - math.cos(half_angle))
    e2x = p1x - (e2r - p1r) / math.tan(half_angle)
    cyl_end_x = e2x - fillet * math.sin(half_angle)

    phi = np.linspace(0.0, half_angle, 12)
    fillet_x = cyl_end_x + fillet * np.sin(phi)
    fillet_r = rc - fillet * (1.0 - np.cos(phi))
    theta = np.linspace(-(math.pi / 2.0 + half_angle), -math.pi / 2.0, 16)
    arc_x = r1 * np.cos(theta)
    arc_r = r1 * np.sin(theta) + r1 + rt
    return (
        np.concatenate([fillet_x, arc_x]),
        np.concatenate([fillet_r, arc_r]),
        cyl_end_x,
    )


def _volume_m3(x: FloatArray, r: FloatArray) -> float:
    return float(np.trapezoid(math.pi * r**2, x))


def _surface_m2(x: FloatArray, r: FloatArray) -> float:
    ds = np.hypot(np.diff(x), np.diff(r))
    return float(np.sum(2.0 * math.pi * 0.5 * (r[1:] + r[:-1]) * ds))


def plate_thickness_m(radius_m: float, pressure_pa: float, allowable_pa: float, nu: float) -> float:
    """Simply supported flat plate, uniform load: s = 3(3+v)qa^2/8t^2 (Roark 7th ed. 11.2-10a)."""
    return radius_m * math.sqrt(3.0 * (3.0 + nu) * pressure_pa / (8.0 * allowable_pa))


def design_thrust_chamber(
    *,
    throat_area_m2: float,
    area_ratio: float,
    l_star_m: float,
    contraction_ratio: float | None,
    convergence_half_angle_deg: float,
    chamber_pressure_pa: float,
    material: Material,
    yield_factor: float,
    ultimate_factor: float,
    min_wall_m: float,
) -> ThrustChamber:
    rt = math.sqrt(throat_area_m2 / math.pi)
    eps_c = contraction_ratio or humble_contraction_ratio(2.0 * rt)
    rc = rt * math.sqrt(eps_c)
    half_angle = math.radians(convergence_half_angle_deg)

    conv_x, conv_r, cyl_end_x = _convergent(rt, rc, half_angle)
    chamber_volume = l_star_m * throat_area_m2
    convergent_volume = _volume_m3(conv_x, conv_r)
    cylinder_length = (chamber_volume - convergent_volume) / (math.pi * rc**2)
    if cylinder_length <= 0.0:
        raise GeometryError(
            f"L* = {l_star_m:.2f} m is smaller than the convergent section alone; raise L* or "
            "lower the contraction ratio"
        )
    injector_x = cyl_end_x - cylinder_length
    nozzle = rao_bell(rt, area_ratio)

    x = np.concatenate([[injector_x], conv_x, nozzle.x_m[1:]])
    r = np.concatenate([[rc], conv_r, nozzle.r_m[1:]])

    allowable = material.allowable_stress_pa(yield_factor, ultimate_factor)
    wall = max(min_wall_m, chamber_pressure_pa * rc / allowable)
    shell_mass = _surface_m2(x, r + 0.5 * wall) * wall * material.density_kg_m3
    plate = max(
        min_wall_m, plate_thickness_m(rc, chamber_pressure_pa, allowable, material.poisson_ratio)
    )
    plate_mass = math.pi * rc**2 * plate * material.density_kg_m3

    return ThrustChamber(
        throat_radius_m=rt,
        chamber_radius_m=rc,
        contraction_ratio=eps_c,
        l_star_m=l_star_m,
        chamber_volume_m3=chamber_volume,
        cylinder_length_m=cylinder_length,
        convergent_length_m=-cyl_end_x,
        convergence_half_angle_deg=convergence_half_angle_deg,
        nozzle=nozzle,
        x_m=x,
        r_m=r,
        wall_thickness_m=wall,
        injector_plate_thickness_m=plate,
        shell_mass_kg=shell_mass,
        injector_plate_mass_kg=plate_mass,
    )
