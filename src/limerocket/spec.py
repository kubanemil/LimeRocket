"""Design spec: YAML in, validated frozen dataclasses (SI) out, with every problem reported."""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from limerocket import materials, propellants, thermo
from limerocket.feed import GASES, FeedMode, Gas
from limerocket.hydraulics import Line
from limerocket.materials import Material
from limerocket.propellants import FRACTION_TOLERANCE, Ingredient, Propellant, Role
from limerocket.trajectory import GENERIC_CD_TABLE
from limerocket.units import M3_PER_L, M_PER_UM, P_STANDARD_PA, bar_to_pa, mm_to_m, pa_to_bar


class SpecError(ValueError):
    def __init__(self, problems: list[str]) -> None:
        self.problems = problems
        super().__init__("invalid design spec:\n  - " + "\n  - ".join(problems))


@dataclass(frozen=True)
class InjectorSide:
    discharge_coefficient: float
    orifices: int


@dataclass(frozen=True)
class InjectorSpec:
    pressure_drop_fraction: float
    oxidizer: InjectorSide
    fuel: InjectorSide


@dataclass(frozen=True)
class EngineSpec:
    thrust_n: float
    burn_time_s: float
    chamber_pressure_pa: float
    design_ambient_pressure_pa: float
    oxidizer: Propellant
    fuel: Propellant
    mixture_ratio: float | None
    l_star_m: float
    contraction_ratio: float | None
    convergence_half_angle_deg: float
    c_star_efficiency: float
    cf_efficiency: float
    wall_material: Material
    injector: InjectorSpec


@dataclass(frozen=True)
class FeedSpec:
    mode: FeedMode
    pressurant: Gas
    pressurant_temperature_k: float
    polytropic_exponent: float
    residual_fraction: float
    oxidizer_line: Line
    fuel_line: Line
    ullage_fraction: float = 0.0
    bottle_pressure_pa: float = 0.0
    bottle_volume_m3: float | None = None
    regulator_dropout_margin_pa: float = 0.0
    pressurant_margin: float = 0.0
    blowdown_ratio: float = 0.0


@dataclass(frozen=True)
class StructureSpec:
    max_diameter_m: float
    tank_material: Material
    bottle_material: Material
    yield_factor: float
    ultimate_factor: float
    min_wall_m: float


@dataclass(frozen=True)
class VehicleSpec:
    other_masses_kg: tuple[tuple[str, float], ...]
    other_length_m: float
    cd_table: tuple[tuple[float, float], ...]
    rail_length_m: float
    launch_altitude_m: float
    apogee_target_m: float | None


@dataclass(frozen=True)
class DesignSpec:
    name: str
    description: str
    engine: EngineSpec
    feed: FeedSpec
    structure: StructureSpec
    vehicle: VehicleSpec


_MISSING: Any = object()
# Choked flow needs Pc/pa above about 1.8 for gamma 1.1 to 1.4; 2 leaves a little margin.
MIN_CHAMBER_TO_AMBIENT = 2.0
AUTO = "auto"


class _Section:
    def __init__(self, data: object, path: str, problems: list[str]) -> None:
        if data is None:
            data = {}
        if not isinstance(data, dict):
            problems.append(f"{path}: expected a mapping, got {type(data).__name__}")
            data = {}
        self.data: dict[str, Any] = data
        self.path = path
        self.problems = problems
        self.used: set[str] = set()

    def key(self, key: str) -> str:
        return f"{self.path}.{key}" if self.path else key

    def raw(self, key: str) -> object:
        self.used.add(key)
        return self.data.get(key, _MISSING)

    def error(self, key: str, message: str) -> None:
        self.problems.append(f"{self.key(key)}: {message}")

    def number(
        self,
        key: str,
        default: float | None = _MISSING,
        *,
        gt: float | None = None,
        ge: float | None = None,
        lt: float | None = None,
        le: float | None = None,
        auto: bool = False,
    ) -> float:
        value = self.raw(key)
        if value is _MISSING or (auto and value == AUTO):
            if default is _MISSING:
                self.error(key, "required")
                return math.nan
            return math.nan if default is None else default
        if isinstance(value, bool) or not isinstance(value, int | float):
            hint = f" or '{AUTO}'" if auto else ""
            self.error(key, f"expected a number{hint}, got {value!r}")
            return math.nan
        value = float(value)
        checks = (
            (gt, "> ", gt is None or value > gt),
            (ge, ">= ", ge is None or value >= ge),
            (lt, "< ", lt is None or value < lt),
            (le, "<= ", le is None or value <= le),
        )
        for bound, text, ok in checks:
            if not ok:
                self.error(key, f"must be {text}{bound:g}, got {value:g}")
        return value

    def optional(self, key: str, **bounds: float) -> float | None:
        value = self.number(key, None, auto=True, **bounds)
        return None if math.isnan(value) else value

    def integer(self, key: str, default: int, *, ge: int = 1) -> int:
        value = self.raw(key)
        if value is _MISSING:
            return default
        if isinstance(value, bool) or not isinstance(value, int):
            self.error(key, f"expected a whole number, got {value!r}")
            return default
        if value < ge:
            self.error(key, f"must be >= {ge}, got {value}")
        return value

    def text(self, key: str, default: str | None = _MISSING) -> str:
        value = self.raw(key)
        if value is _MISSING:
            if default is _MISSING:
                self.error(key, "required")
                return ""
            return default or ""
        if not isinstance(value, str):
            self.error(key, f"expected text, got {value!r}")
            return ""
        return value

    def choice[T](self, key: str, options: dict[str, T], default: str) -> T:
        name = self.text(key, default)
        if name not in options:
            self.error(key, f"unknown value {name!r} (choose from: {', '.join(options)})")
            return options[default]
        return options[name]

    def section(self, key: str) -> _Section:
        value = self.raw(key)
        return _Section(None if value is _MISSING else value, self.key(key), self.problems)

    def forbid(self, keys: tuple[str, ...], reason: str) -> None:
        for key in keys:
            if key in self.data:
                self.used.add(key)
                self.error(key, reason)

    def finish(self) -> None:
        for key in self.data:
            if key not in self.used:
                known = ", ".join(sorted(self.used))
                self.error(key, f"unknown key (known here: {known})")


def _propellant(section: _Section, key: str, role: Role) -> Propellant | None:
    value = section.raw(key)
    if value is _MISSING:
        section.error(key, "required")
        return None
    if isinstance(value, str):
        try:
            return propellants.get(value, role)
        except KeyError as exc:
            section.error(key, str(exc.args[0]))
            return None
    return _custom_propellant(_Section(value, section.key(key), section.problems), role)


def _custom_propellant(custom: _Section, role: Role) -> Propellant | None:
    name = custom.text("name")
    ingredients_raw = custom.raw("ingredients")
    ingredients: list[Ingredient] = []
    if not isinstance(ingredients_raw, dict) or not ingredients_raw:
        custom.error("ingredients", "expected a mapping of CEA species name to mass fraction")
    else:
        for species, fraction in ingredients_raw.items():
            if isinstance(fraction, bool) or not isinstance(fraction, int | float) or fraction <= 0:
                custom.error("ingredients", f"{species}: mass fraction must be a positive number")
            else:
                ingredients.append(Ingredient(str(species), float(fraction)))
    temperature = custom.number("temperature_k", gt=0)
    density = custom.number("density_kg_m3", gt=0)
    viscosity = custom.number("viscosity_pa_s", gt=0)
    custom.finish()
    total = sum(i.mass_fraction for i in ingredients)
    if ingredients and abs(total - 1.0) > FRACTION_TOLERANCE:
        custom.error("ingredients", f"mass fractions must sum to 1, got {total:g}")
        return None
    missing = thermo.unknown_species([i.cea_name for i in ingredients])
    if missing:
        custom.error("ingredients", f"not in the CEA database: {', '.join(missing)}")
        return None
    if not ingredients or any(math.isnan(v) for v in (temperature, density, viscosity)):
        return None
    return Propellant(
        key=name or "custom",
        name=name or "custom",
        role=role,
        ingredients=tuple(Ingredient(i.cea_name, i.mass_fraction / total) for i in ingredients),
        temperature_k=temperature,
        density_kg_m3=density,
        viscosity_pa_s=viscosity,
        source="user spec",
    )


def _injector_side(section: _Section) -> InjectorSide:
    side = InjectorSide(
        discharge_coefficient=section.number("discharge_coefficient", 0.7, gt=0, le=1),
        orifices=section.integer("orifices", 8),
    )
    section.finish()
    return side


def _line(section: _Section) -> Line:
    line = Line(
        length_m=section.number("length_m", ge=0),
        inner_diameter_m=mm_to_m(section.number("inner_diameter_mm", gt=0)),
        roughness_m=section.number("roughness_um", 1.5, ge=0) * M_PER_UM,
        minor_loss_k=section.number("minor_loss_k", 0.0, ge=0),
    )
    section.finish()
    return line


def _engine(s: _Section) -> EngineSpec | None:
    oxidizer = _propellant(s, "oxidizer", Role.OXIDIZER)
    fuel = _propellant(s, "fuel", Role.FUEL)
    mixture = s.optional("mixture_ratio", gt=0)
    material = s.choice("wall_material", materials.CATALOG, "ss316")
    inj = s.section("injector")
    injector = InjectorSpec(
        pressure_drop_fraction=inj.number("pressure_drop_fraction", 0.20, gt=0, lt=2),
        oxidizer=_injector_side(inj.section("oxidizer")),
        fuel=_injector_side(inj.section("fuel")),
    )
    inj.finish()
    thrust = s.number("thrust_n", gt=0)
    burn_time = s.number("burn_time_s", gt=0)
    chamber_pressure = s.number("chamber_pressure_bar", gt=0)
    ambient = s.number("design_ambient_pressure_bar", pa_to_bar(P_STANDARD_PA), gt=0)
    if chamber_pressure < MIN_CHAMBER_TO_AMBIENT * ambient:
        s.error(
            "chamber_pressure_bar",
            f"must be at least {MIN_CHAMBER_TO_AMBIENT:g} x design_ambient_pressure_bar "
            f"({MIN_CHAMBER_TO_AMBIENT * ambient:g} bar) for the nozzle to choke, "
            f"got {chamber_pressure:g}",
        )
    l_star = s.number("l_star_m", 1.0, gt=0)
    contraction = s.optional("contraction_ratio", gt=1)
    half_angle = s.number("convergence_half_angle_deg", 30.0, ge=15, le=60)
    c_star_eff = s.number("c_star_efficiency", 0.92, gt=0, le=1)
    cf_eff = s.number("thrust_coefficient_efficiency", 0.98, gt=0, le=1)
    s.finish()
    if oxidizer is None or fuel is None:
        return None
    return EngineSpec(
        thrust_n=thrust,
        burn_time_s=burn_time,
        chamber_pressure_pa=bar_to_pa(chamber_pressure),
        design_ambient_pressure_pa=bar_to_pa(ambient),
        oxidizer=oxidizer,
        fuel=fuel,
        mixture_ratio=mixture,
        l_star_m=l_star,
        contraction_ratio=contraction,
        convergence_half_angle_deg=half_angle,
        c_star_efficiency=c_star_eff,
        cf_efficiency=cf_eff,
        wall_material=material,
        injector=injector,
    )


REGULATED_KEYS = (
    "bottle_pressure_bar",
    "bottle_volume_l",
    "regulator_dropout_margin_bar",
    "pressurant_margin",
    "ullage_fraction",
)
BLOWDOWN_KEYS = ("blowdown_ratio", "polytropic_exponent")


def _feed(s: _Section) -> FeedSpec:
    mode = s.choice("mode", {m.value: m for m in FeedMode}, FeedMode.REGULATED.value)
    gas: Gas = s.choice("pressurant", GASES, "helium")
    lines = s.section("lines")
    common: dict[str, Any] = {
        "mode": mode,
        "pressurant": gas,
        "pressurant_temperature_k": s.number("pressurant_temperature_k", 293.15, gt=0),
        "polytropic_exponent": gas.gamma,
        "residual_fraction": s.number("residual_fraction", 0.03, ge=0, lt=1),
        "oxidizer_line": _line(lines.section("oxidizer")),
        "fuel_line": _line(lines.section("fuel")),
    }
    lines.finish()
    if mode is FeedMode.BLOWDOWN:
        s.forbid(REGULATED_KEYS, "applies only to mode: regulated")
        exponent = s.optional("polytropic_exponent", ge=1, le=1.7)
        if exponent is not None:
            common["polytropic_exponent"] = exponent
        spec = FeedSpec(**common, blowdown_ratio=s.number("blowdown_ratio", 2.5, gt=1))
    else:
        s.forbid(BLOWDOWN_KEYS, "applies only to mode: blowdown")
        volume_l = s.optional("bottle_volume_l", gt=0)
        spec = FeedSpec(
            **common,
            ullage_fraction=s.number("ullage_fraction", 0.05, gt=0, lt=1),
            bottle_pressure_pa=bar_to_pa(s.number("bottle_pressure_bar", 300.0, gt=0)),
            bottle_volume_m3=None if volume_l is None else volume_l * M3_PER_L,
            regulator_dropout_margin_pa=bar_to_pa(
                s.number("regulator_dropout_margin_bar", 10.0, ge=0)
            ),
            pressurant_margin=s.number("pressurant_margin", 0.20, ge=0),
        )
    s.finish()
    return spec


def _structure(s: _Section) -> StructureSpec:
    spec = StructureSpec(
        max_diameter_m=mm_to_m(s.number("max_diameter_mm", gt=0)),
        tank_material=s.choice("tank_material", materials.CATALOG, "al6061-t6"),
        bottle_material=s.choice("bottle_material", materials.CATALOG, "ti6al4v"),
        yield_factor=s.number("yield_safety_factor", 1.5, ge=1),
        ultimate_factor=s.number("ultimate_safety_factor", 2.0, ge=1),
        min_wall_m=mm_to_m(s.number("min_wall_mm", 1.0, gt=0)),
    )
    s.finish()
    return spec


def _cd_table(s: _Section) -> tuple[tuple[float, float], ...]:
    raw = s.raw("drag_coefficient")
    if raw is _MISSING:
        return GENERIC_CD_TABLE
    if isinstance(raw, int | float) and not isinstance(raw, bool) and raw > 0:
        return ((0.0, float(raw)), (10.0, float(raw)))

    def is_pair(row: object) -> bool:
        return (
            isinstance(row, list)
            and len(row) == 2
            and all(isinstance(v, int | float) and not isinstance(v, bool) for v in row)
        )

    raw_rows = raw if isinstance(raw, list) else []
    rows = [(float(row[0]), float(row[1])) for row in raw_rows if is_pair(row)]
    if (
        not rows
        or len(rows) != len(raw_rows)
        or any(b[0] <= a[0] for a, b in itertools.pairwise(rows))
    ):
        s.error(
            "drag_coefficient",
            "expected a number or a list of [mach, cd] pairs with increasing Mach",
        )
        return GENERIC_CD_TABLE
    return tuple(rows)


def _vehicle(s: _Section) -> VehicleSpec:
    masses_raw = s.raw("other_masses_kg")
    masses: list[tuple[str, float]] = []
    if masses_raw is not _MISSING:
        if not isinstance(masses_raw, dict):
            s.error("other_masses_kg", "expected a mapping of item name to mass in kg")
        else:
            for item, mass in masses_raw.items():
                if isinstance(mass, bool) or not isinstance(mass, int | float) or mass < 0:
                    s.error("other_masses_kg", f"{item}: mass must be a number >= 0")
                else:
                    masses.append((str(item), float(mass)))
    spec = VehicleSpec(
        other_masses_kg=tuple(masses),
        other_length_m=s.number("other_length_m", 0.0, ge=0),
        cd_table=_cd_table(s),
        rail_length_m=s.number("rail_length_m", 6.0, ge=0),
        launch_altitude_m=s.number("launch_altitude_m", 0.0, ge=-400, le=8000),
        apogee_target_m=s.optional("apogee_target_m", gt=0),
    )
    s.finish()
    return spec


def parse(data: object) -> DesignSpec:
    problems: list[str] = []
    root = _Section(data, "", problems)
    name = root.text("name")
    description = root.text("description", "")
    engine = _engine(root.section("engine"))
    feed = _feed(root.section("feed"))
    structure = _structure(root.section("structure"))
    vehicle = _vehicle(root.section("vehicle"))
    root.finish()
    if problems or engine is None:
        raise SpecError(problems)
    return DesignSpec(name, description, engine, feed, structure, vehicle)


def load(path: Path) -> DesignSpec:
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise SpecError([f"{path}: not valid YAML ({exc})"]) from exc
    return parse(data)
