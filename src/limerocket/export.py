"""Exports: chamber contour (CSV, DXF, STL), thrust curves (.eng, RocketPy CSV), time series."""

from __future__ import annotations

import csv
import math
import struct
from pathlib import Path

import numpy as np
import numpy.typing as npt

from limerocket import __version__
from limerocket.chamber import ThrustChamber
from limerocket.design import Design
from limerocket.units import P_STANDARD_PA

FloatArray = npt.NDArray[np.float64]
MM_PER_M = 1000.0
STL_SEGMENTS = 64
ENG_POINTS = 60
ENG_RAMP_S = 0.01  # RASP curves start and end at zero thrust; 10 ms edges stand in for transients


def outer_wall(chamber: ThrustChamber) -> tuple[FloatArray, FloatArray]:
    """Inner contour offset outward along its normal by the wall thickness."""
    x, r = chamber.x_m, chamber.r_m
    dx, dr = np.gradient(x), np.gradient(r)
    norm = np.hypot(dx, dr)
    xo = x - chamber.wall_thickness_m * dr / norm
    ro = r + chamber.wall_thickness_m * dx / norm
    keep = np.concatenate([[True], np.diff(xo) > 0.0])
    return xo[keep], ro[keep]


def write_contour_csv(path: Path, chamber: ThrustChamber) -> None:
    xo, ro = outer_wall(chamber)
    outer = np.interp(chamber.x_m, xo, ro)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["x_mm", "inner_radius_mm", "outer_radius_mm"])
        for x, r, o in zip(chamber.x_m, chamber.r_m, outer, strict=True):
            writer.writerow([f"{x * MM_PER_M:.4f}", f"{r * MM_PER_M:.4f}", f"{o * MM_PER_M:.4f}"])


def _dxf_polyline(points: list[tuple[float, float]], layer: str, closed: bool = False) -> list[str]:
    out = ["0", "POLYLINE", "8", layer, "66", "1", "70", "1" if closed else "0"]
    for x, y in points:
        out += ["0", "VERTEX", "8", layer, "10", f"{x:.4f}", "20", f"{y:.4f}", "30", "0.0"]
    out += ["0", "SEQEND", "8", layer]
    return out


def write_dxf(path: Path, chamber: ThrustChamber) -> None:
    """AutoCAD R12 ASCII DXF of the half-section in millimetres: wall profile plus centreline."""
    xo, ro = outer_wall(chamber)
    inner = list(zip(chamber.x_m * MM_PER_M, chamber.r_m * MM_PER_M, strict=True))
    outer = list(zip(xo * MM_PER_M, ro * MM_PER_M, strict=True))
    profile = inner + outer[::-1]
    x0, x1 = float(chamber.x_m[0] * MM_PER_M), float(chamber.x_m[-1] * MM_PER_M)
    lines = ["0", "SECTION", "2", "HEADER", "9", "$ACADVER", "1", "AC1009", "0", "ENDSEC"]
    lines += ["0", "SECTION", "2", "ENTITIES"]
    lines += _dxf_polyline(profile, "WALL", closed=True)
    lines += _dxf_polyline([(x0, 0.0), (x1, 0.0)], "CENTERLINE")
    lines += ["0", "ENDSEC", "0", "EOF"]
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def revolve(x: FloatArray, r: FloatArray, segments: int = STL_SEGMENTS) -> FloatArray:
    """Triangles (n, 3, 3) of a closed (x, r) profile revolved about the x axis."""
    theta = np.linspace(0.0, 2.0 * math.pi, segments, endpoint=False)
    ring = np.stack(
        [
            np.broadcast_to(x[:, None], (x.size, segments)),
            r[:, None] * np.cos(theta)[None, :],
            r[:, None] * np.sin(theta)[None, :],
        ],
        axis=-1,
    )
    nxt_i = np.roll(np.arange(x.size), -1)
    nxt_j = np.roll(np.arange(segments), -1)
    i, j = np.meshgrid(np.arange(x.size), np.arange(segments), indexing="ij")
    a, b = ring[i, j], ring[nxt_i[i], j]
    c, d = ring[nxt_i[i], nxt_j[j]], ring[i, nxt_j[j]]
    tris = np.concatenate([np.stack([a, b, c], -2), np.stack([a, c, d], -2)]).reshape(-1, 3, 3)
    area = np.linalg.norm(np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0]), axis=1)
    return tris[area > 0.0]


def wall_solid(chamber: ThrustChamber, segments: int = STL_SEGMENTS) -> FloatArray:
    """Watertight wall: inner contour downstream, outer contour back upstream, revolved."""
    xo, ro = outer_wall(chamber)
    x = np.concatenate([chamber.x_m, xo[::-1]]) * MM_PER_M
    r = np.concatenate([chamber.r_m, ro[::-1]]) * MM_PER_M
    return revolve(x, r, segments)


def write_stl(path: Path, chamber: ThrustChamber, segments: int = STL_SEGMENTS) -> None:
    """Binary STL of the chamber wall in millimetres."""
    tris = wall_solid(chamber, segments).astype(np.float32)
    normals = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    normals /= np.linalg.norm(normals, axis=1, keepdims=True)
    record = np.zeros(len(tris), dtype=[("n", "<f4", (3,)), ("v", "<f4", (3, 3)), ("attr", "<u2")])
    record["n"], record["v"] = normals, tris
    header = b"LimeRocket thrust chamber wall, mm".ljust(80, b" ")
    with path.open("wb") as fh:
        fh.write(header)
        fh.write(struct.pack("<I", len(tris)))
        fh.write(record.tobytes())


def _thrust_samples(
    design: Design, ambient_pa: float, points: int
) -> tuple[FloatArray, FloatArray]:
    burn = design.burn
    t = np.linspace(0.0, burn.burn_time_s, points)
    f = np.interp(t, burn.time_s, burn.thrust_n(ambient_pa))
    return t, np.maximum(f, 0.0)


def write_eng(path: Path, design: Design) -> None:
    """RASP .eng thrust curve (OpenRocket, RockSim); masses cover the whole propulsion system."""
    t, f = _thrust_samples(design, P_STANDARD_PA, ENG_POINTS)
    t = t + ENG_RAMP_S
    t[0] = ENG_RAMP_S
    expelled = design.burn.expelled_mass_kg
    propulsion = (
        design.chamber.mass_kg
        + design.oxidizer_tank.mass_kg
        + design.fuel_tank.mass_kg
        + (design.bottle.mass_kg if design.bottle else 0.0)
        + design.pressurant_mass_kg
        + design.propellant_mass_kg
    )
    length = (
        design.chamber.length_m
        + design.oxidizer_tank.length_m
        + design.fuel_tank.length_m
        + (design.bottle.length_m if design.bottle else 0.0)
    )
    spec = design.spec
    code = "".join(ch for ch in spec.name if ch.isalnum() or ch in "-_") or "LIME"
    lines = [
        f"; {spec.name}: {spec.engine.oxidizer.name} / {spec.engine.fuel.name}",
        f"; LimeRocket {__version__}, sea-level thrust, liquid pressure-fed",
        (
            f"{code} {spec.structure.max_diameter_m * MM_PER_M:.0f} {length * MM_PER_M:.0f} P "
            f"{expelled:.4f} {propulsion:.4f} LimeRocket"
        ),
    ]
    lines += [f"   {ti:.4f} {fi:.2f}" for ti, fi in zip(t, f, strict=True)]
    lines.append(f"   {t[-1] + ENG_RAMP_S:.4f} 0.00")
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def write_rocketpy_csv(path: Path, design: Design) -> None:
    """Two columns, time [s] and thrust [N], no header: RocketPy's thrust_source format."""
    burn = design.burn
    thrust = np.maximum(burn.thrust_n(P_STANDARD_PA), 0.0)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        for t, f in zip(burn.time_s, thrust, strict=True):
            writer.writerow([f"{t:.5f}", f"{f:.3f}"])


def write_burn_csv(path: Path, design: Design) -> None:
    burn = design.burn
    pa = design.ambient_pa
    columns = {
        "time_s": burn.time_s,
        "chamber_pressure_bar": burn.chamber_pressure_pa / 1e5,
        "thrust_n": burn.thrust_n(pa),
        "thrust_vacuum_n": burn.thrust_vacuum_n,
        "isp_s": burn.isp_s(pa),
        "oxidizer_flow_kg_s": burn.oxidizer_flow_kg_s,
        "fuel_flow_kg_s": burn.fuel_flow_kg_s,
        "mixture_ratio": burn.mixture_ratio,
        "oxidizer_tank_bar": burn.oxidizer_tank_pressure_pa / 1e5,
        "fuel_tank_bar": burn.fuel_tank_pressure_pa / 1e5,
        "bottle_bar": burn.bottle_pressure_pa / 1e5,
        "oxidizer_mass_kg": burn.oxidizer_mass_kg,
        "fuel_mass_kg": burn.fuel_mass_kg,
    }
    _write_columns(path, columns)


def write_flight_csv(path: Path, design: Design) -> None:
    fl = design.flight
    columns = {
        "time_s": fl.time_s,
        "altitude_m": fl.altitude_m,
        "velocity_m_s": fl.velocity_m_s,
        "acceleration_m_s2": fl.acceleration_m_s2,
        "mach": fl.mach,
        "dynamic_pressure_pa": fl.dynamic_pressure_pa,
        "mass_kg": fl.mass_kg,
        "thrust_n": fl.thrust_n,
        "drag_n": fl.drag_n,
    }
    _write_columns(path, columns)


def _write_columns(path: Path, columns: dict[str, FloatArray]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(columns)
        for row in zip(*columns.values(), strict=True):
            writer.writerow([f"{v:.6g}" for v in row])
