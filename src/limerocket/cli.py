"""Command line: `limerocket design spec.yaml -o out/`."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Callable
from pathlib import Path

from limerocket import __version__, export, materials, propellants
from limerocket.chamber import GeometryError
from limerocket.design import Design, DesignError, design
from limerocket.feed import FeedError
from limerocket.plots import overview_figure, to_svg
from limerocket.report import html_report, terminal_summary
from limerocket.spec import SpecError, load
from limerocket.thermo import ThermoError
from limerocket.units import PA_PER_MPA


def write_outputs(result: Design, out: Path) -> list[str]:
    out.mkdir(parents=True, exist_ok=True)
    stem = result.spec.name.lower().replace(" ", "-")
    writers: dict[str, Callable[[Path], object]] = {
        f"{stem}-contour.csv": lambda p: export.write_contour_csv(p, result.chamber),
        f"{stem}-contour.dxf": lambda p: export.write_dxf(p, result.chamber),
        f"{stem}-chamber.stl": lambda p: export.write_stl(p, result.chamber),
        f"{stem}.eng": lambda p: export.write_eng(p, result),
        f"{stem}-thrust-rocketpy.csv": lambda p: export.write_rocketpy_csv(p, result),
        f"{stem}-burn.csv": lambda p: export.write_burn_csv(p, result),
        f"{stem}-flight.csv": lambda p: export.write_flight_csv(p, result),
        "overview.svg": lambda p: p.write_text(to_svg(overview_figure(result)), encoding="utf-8"),
    }
    for name, write in writers.items():
        write(out / name)
    files = list(writers)
    (out / "report.html").write_text(html_report(result, files), encoding="utf-8")
    return ["report.html", *files]


def _design(args: argparse.Namespace) -> int:
    result = design(load(args.spec))
    files = write_outputs(result, args.output)
    print(terminal_summary(result))
    print(f"\nWrote {len(files)} files to {args.output}/ (open report.html)")
    return 0


def _catalog(_args: argparse.Namespace) -> int:
    print("Propellants (key: name, role, density):")
    for p in propellants.CATALOG.values():
        print(f"  {p.key:10s} {p.name}, {p.role}, {p.density_kg_m3:g} kg/m3")
    print("\nMaterials (key: name, yield / ultimate):")
    for m in materials.CATALOG.values():
        strength = (
            f"{m.yield_strength_pa / PA_PER_MPA:g} / {m.ultimate_strength_pa / PA_PER_MPA:g} MPa"
        )
        print(f"  {m.key:10s} {m.name}, {strength}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="limerocket", description=__doc__)
    parser.add_argument("--version", action="version", version=f"limerocket {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("design", help="design a rocket from a YAML spec")
    run.add_argument("spec", type=Path)
    run.add_argument("-o", "--output", type=Path, default=Path("out"))
    run.set_defaults(func=_design)
    cat = sub.add_parser("catalog", help="list built-in propellants and materials")
    cat.set_defaults(func=_catalog)
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except (SpecError, DesignError, GeometryError, FeedError, ThermoError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except FileNotFoundError as exc:
        print(f"error: {exc.strerror}: {exc.filename}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
