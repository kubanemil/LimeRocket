"""Room-temperature structural properties (typical values, ASM Handbook / MatWeb)."""

from __future__ import annotations

from dataclasses import dataclass

from limerocket.units import PA_PER_MPA as MPA


@dataclass(frozen=True)
class Material:
    key: str
    name: str
    density_kg_m3: float
    yield_strength_pa: float
    ultimate_strength_pa: float
    poisson_ratio: float

    def allowable_stress_pa(self, yield_factor: float, ultimate_factor: float) -> float:
        return min(
            self.yield_strength_pa / yield_factor, self.ultimate_strength_pa / ultimate_factor
        )


CATALOG: dict[str, Material] = {
    m.key: m
    for m in (
        Material("al6061-t6", "Aluminium 6061-T6", 2700.0, 276 * MPA, 310 * MPA, 0.33),
        Material("al7075-t6", "Aluminium 7075-T6", 2810.0, 503 * MPA, 572 * MPA, 0.33),
        Material("ss304", "Stainless steel 304, annealed", 8000.0, 215 * MPA, 505 * MPA, 0.29),
        Material("ss316", "Stainless steel 316, annealed", 8000.0, 205 * MPA, 515 * MPA, 0.30),
        Material("steel4130", "Steel 4130, normalized", 7850.0, 435 * MPA, 670 * MPA, 0.29),
        Material("ti6al4v", "Titanium Ti-6Al-4V, annealed", 4430.0, 880 * MPA, 950 * MPA, 0.34),
        Material("copper", "Copper C10100, annealed", 8940.0, 69 * MPA, 220 * MPA, 0.34),
    )
}


def get(key: str) -> Material:
    try:
        return CATALOG[key]
    except KeyError:
        raise KeyError(f"unknown material {key!r} (known: {', '.join(CATALOG)})") from None
