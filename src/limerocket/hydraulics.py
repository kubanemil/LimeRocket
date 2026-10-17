"""Incompressible feed hydraulics: Darcy-Weisbach lines and injector orifices."""

from __future__ import annotations

import math
from dataclasses import dataclass

from scipy.optimize import brentq

LAMINAR_RE = 2300.0
TURBULENT_RE = 4000.0


def haaland(reynolds: float, relative_roughness: float) -> float:
    """Explicit Darcy friction factor, Haaland (1983), within ~2 % of Colebrook."""
    inv_sqrt = -1.8 * math.log10((relative_roughness / 3.7) ** 1.11 + 6.9 / reynolds)
    return 1.0 / inv_sqrt**2


def colebrook(reynolds: float, relative_roughness: float) -> float:
    """Darcy friction factor from the Colebrook-White equation, solved by Newton iteration."""
    a = relative_roughness / 3.7
    b = 2.51 / reynolds
    x = 1.0 / math.sqrt(haaland(reynolds, relative_roughness))
    for _ in range(50):
        g = x + 2.0 * math.log10(a + b * x)
        dg = 1.0 + 2.0 * b / ((a + b * x) * math.log(10.0))
        step = g / dg
        x -= step
        if abs(step) < 1e-13 * x:
            break
    return 1.0 / x**2


def darcy_friction_factor(reynolds: float, relative_roughness: float) -> float:
    """64/Re when laminar, Colebrook when turbulent, linear blend across 2300 < Re < 4000."""
    if reynolds <= 0.0:
        return 0.0
    if reynolds < LAMINAR_RE:
        return 64.0 / reynolds
    if reynolds < TURBULENT_RE:
        w = (reynolds - LAMINAR_RE) / (TURBULENT_RE - LAMINAR_RE)
        return (1.0 - w) * 64.0 / LAMINAR_RE + w * colebrook(TURBULENT_RE, relative_roughness)
    return colebrook(reynolds, relative_roughness)


@dataclass(frozen=True)
class Fluid:
    density_kg_m3: float
    viscosity_pa_s: float


@dataclass(frozen=True)
class Line:
    length_m: float
    inner_diameter_m: float
    roughness_m: float
    minor_loss_k: float

    @property
    def area_m2(self) -> float:
        return math.pi * self.inner_diameter_m**2 / 4.0

    def velocity_m_s(self, mass_flow_kg_s: float, fluid: Fluid) -> float:
        return mass_flow_kg_s / (fluid.density_kg_m3 * self.area_m2)

    def reynolds(self, mass_flow_kg_s: float, fluid: Fluid) -> float:
        return 4.0 * mass_flow_kg_s / (math.pi * self.inner_diameter_m * fluid.viscosity_pa_s)

    def friction_factor(self, mass_flow_kg_s: float, fluid: Fluid) -> float:
        return darcy_friction_factor(
            self.reynolds(mass_flow_kg_s, fluid), self.roughness_m / self.inner_diameter_m
        )

    def pressure_drop_pa(self, mass_flow_kg_s: float, fluid: Fluid) -> float:
        """Darcy-Weisbach: dp = (f L/D + K) rho v^2 / 2."""
        v = self.velocity_m_s(mass_flow_kg_s, fluid)
        f = self.friction_factor(mass_flow_kg_s, fluid)
        k = f * self.length_m / self.inner_diameter_m + self.minor_loss_k
        return k * 0.5 * fluid.density_kg_m3 * v**2


@dataclass(frozen=True)
class Orifices:
    count: int
    diameter_m: float
    discharge_coefficient: float

    @property
    def area_m2(self) -> float:
        return self.count * math.pi * self.diameter_m**2 / 4.0

    def mass_flow_kg_s(self, pressure_drop_pa: float, fluid: Fluid) -> float:
        """Orifice equation: mdot = Cd A sqrt(2 rho dp)."""
        if pressure_drop_pa <= 0.0:
            return 0.0
        return (
            self.discharge_coefficient
            * self.area_m2
            * math.sqrt(2.0 * fluid.density_kg_m3 * pressure_drop_pa)
        )

    def pressure_drop_pa(self, mass_flow_kg_s: float, fluid: Fluid) -> float:
        return (mass_flow_kg_s / (self.discharge_coefficient * self.area_m2)) ** 2 / (
            2.0 * fluid.density_kg_m3
        )

    @classmethod
    def sized_for(
        cls,
        mass_flow_kg_s: float,
        pressure_drop_pa: float,
        fluid: Fluid,
        count: int,
        discharge_coefficient: float,
    ) -> Orifices:
        area = mass_flow_kg_s / (
            discharge_coefficient * math.sqrt(2.0 * fluid.density_kg_m3 * pressure_drop_pa)
        )
        return cls(count, math.sqrt(4.0 * area / (math.pi * count)), discharge_coefficient)


@dataclass(frozen=True)
class Branch:
    """Tank outlet to chamber: one feed line in series with the injector orifices."""

    fluid: Fluid
    line: Line
    orifices: Orifices

    def pressure_drop_pa(self, mass_flow_kg_s: float) -> float:
        return self.line.pressure_drop_pa(
            mass_flow_kg_s, self.fluid
        ) + self.orifices.pressure_drop_pa(mass_flow_kg_s, self.fluid)

    def mass_flow_kg_s(self, pressure_drop_pa: float) -> float:
        upper = self.orifices.mass_flow_kg_s(pressure_drop_pa, self.fluid)
        if upper <= 0.0:
            return 0.0
        return float(
            brentq(
                lambda m: self.pressure_drop_pa(m) - pressure_drop_pa,
                0.0,
                upper,
                xtol=1e-12 * upper,
                rtol=1e-12,
            )
        )
