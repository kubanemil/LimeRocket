"""Thin-wall pressure vessels: cylinders with hemispherical heads, or spheres when they fit."""

from __future__ import annotations

import math
from dataclasses import dataclass

from limerocket.materials import Material


@dataclass(frozen=True)
class Vessel:
    name: str
    material: Material
    pressure_pa: float
    volume_m3: float
    inner_radius_m: float
    cylinder_length_m: float
    cylinder_wall_m: float
    head_wall_m: float
    allowable_stress_pa: float
    mass_kg: float

    @property
    def is_sphere(self) -> bool:
        return self.cylinder_length_m == 0.0

    @property
    def outer_diameter_m(self) -> float:
        return 2.0 * (self.inner_radius_m + max(self.cylinder_wall_m, self.head_wall_m))

    @property
    def length_m(self) -> float:
        return self.cylinder_length_m + 2.0 * (self.inner_radius_m + self.head_wall_m)

    @property
    def hoop_stress_pa(self) -> float:
        """sigma_h = p r / t (cylinder); zero for a sphere."""
        if self.is_sphere:
            return 0.0
        return self.pressure_pa * self.inner_radius_m / self.cylinder_wall_m

    @property
    def longitudinal_stress_pa(self) -> float:
        """sigma_l = p r / 2t (cylinder); zero for a sphere."""
        return 0.5 * self.hoop_stress_pa

    @property
    def head_stress_pa(self) -> float:
        """sigma = p r / 2t (hemisphere or sphere)."""
        return self.pressure_pa * self.inner_radius_m / (2.0 * self.head_wall_m)

    @property
    def margin_of_safety(self) -> float:
        """Allowable (already divided by the safety factors) over the peak stress, minus one."""
        peak = max(self.hoop_stress_pa, self.head_stress_pa)
        return self.allowable_stress_pa / peak - 1.0


def size_vessel(
    name: str,
    *,
    volume_m3: float,
    max_outer_radius_m: float,
    pressure_pa: float,
    material: Material,
    yield_factor: float,
    ultimate_factor: float,
    min_wall_m: float,
) -> Vessel:
    allowable = material.allowable_stress_pa(yield_factor, ultimate_factor)
    sphere_radius = (3.0 * volume_m3 / (4.0 * math.pi)) ** (1.0 / 3.0)
    head_wall_sphere = max(min_wall_m, pressure_pa * sphere_radius / (2.0 * allowable))
    if sphere_radius + head_wall_sphere <= max_outer_radius_m:
        radius, length, cyl_wall, head_wall = sphere_radius, 0.0, 0.0, head_wall_sphere
    else:
        cyl_wall = max(min_wall_m, pressure_pa * max_outer_radius_m / allowable)
        radius = max_outer_radius_m - cyl_wall
        head_wall = max(min_wall_m, pressure_pa * radius / (2.0 * allowable))
        length = (volume_m3 - 4.0 / 3.0 * math.pi * radius**3) / (math.pi * radius**2)
    shell = (
        2.0 * math.pi * (radius + 0.5 * cyl_wall) * length * cyl_wall
        + 4.0 * math.pi * (radius + 0.5 * head_wall) ** 2 * head_wall
    )
    return Vessel(
        name=name,
        material=material,
        pressure_pa=pressure_pa,
        volume_m3=volume_m3,
        inner_radius_m=radius,
        cylinder_length_m=length,
        cylinder_wall_m=cyl_wall,
        head_wall_m=head_wall,
        allowable_stress_pa=allowable,
        mass_kg=shell * material.density_kg_m3,
    )
