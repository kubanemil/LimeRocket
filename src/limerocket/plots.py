from __future__ import annotations

import io
from collections.abc import Iterator
from contextlib import contextmanager
from typing import Any

import matplotlib as mpl
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.typing import RcKeyType

from limerocket.design import Design
from limerocket.export import outer_wall
from limerocket.units import G0_M_S2, PA_PER_KPA, m_to_mm, pa_to_bar

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SERIES = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100")
WALL = "#9ec5f4"

STYLE: dict[RcKeyType, Any] = {
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "axes.edgecolor": AXIS,
    "axes.labelcolor": INK_2,
    "axes.titlecolor": INK,
    "axes.titlesize": 10,
    "axes.titleweight": "bold",
    "axes.titlelocation": "left",
    "axes.labelsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.8,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "lines.linewidth": 2.0,
    "lines.solid_capstyle": "round",
    "legend.frameon": False,
    "legend.fontsize": 8,
    "font.family": "sans-serif",
    "font.sans-serif": ["Inter", "Helvetica", "Arial", "DejaVu Sans"],
    "svg.fonttype": "none",
    "svg.hashsalt": "limerocket",
}


@contextmanager
def _style() -> Iterator[None]:
    with mpl.rc_context(STYLE):
        yield


def to_svg(fig: Figure) -> str:
    buf = io.StringIO()
    with _style():
        fig.savefig(buf, format="svg", bbox_inches="tight", metadata={"Date": None})
    return buf.getvalue()


def _engine_axes(ax: Axes, design: Design, *, annotate: bool) -> None:
    ch = design.chamber
    x, r = m_to_mm(ch.x_m), m_to_mm(ch.r_m)
    xo, ro = outer_wall(ch)
    xo, ro = m_to_mm(xo), m_to_mm(ro)
    for sign in (1, -1):
        ax.fill(
            np.concatenate([x, xo[::-1]]),
            sign * np.concatenate([r, ro[::-1]]),
            color=WALL,
            linewidth=0,
        )
        ax.plot(x, sign * r, color=SERIES[0], linewidth=1.5)
    ax.axhline(0.0, color=MUTED, linewidth=0.8, linestyle=(0, (6, 3, 1, 3)))
    ax.set_aspect("equal")
    ax.set_xlabel("Axial position from throat [mm]")
    ax.set_ylabel("Radius [mm]")
    if annotate:
        ax.annotate(
            f"Dt {m_to_mm(2 * ch.throat_radius_m):.1f} mm",
            (0.0, 0.0),
            xytext=(0, 4),
            textcoords="offset points",
            ha="center",
            fontsize=8,
            color=INK_2,
        )
        ax.annotate(
            f"Dc {m_to_mm(2 * ch.chamber_radius_m):.1f} mm",
            (x[0], m_to_mm(ch.chamber_radius_m)),
            xytext=(4, 6),
            textcoords="offset points",
            fontsize=8,
            color=INK_2,
        )
        ax.annotate(
            f"De {m_to_mm(2 * ch.exit_radius_m):.1f} mm",
            (x[-1], m_to_mm(ch.exit_radius_m)),
            xytext=(-4, 6),
            textcoords="offset points",
            ha="right",
            fontsize=8,
            color=INK_2,
        )
    ax.margins(x=0.04, y=0.25)


def engine_figure(design: Design) -> Figure:
    with _style():
        fig = Figure(figsize=(7.5, 3.0), layout="constrained")
        ax = fig.add_subplot()
        _engine_axes(ax, design, annotate=True)
        ax.set_title(
            f"Thrust chamber: L* {design.chamber.l_star_m:.2f} m, "
            f"Ac/At {design.chamber.contraction_ratio:.2f}, Ae/At {design.ideal.area_ratio:.2f}, "
            f"80 % Rao bell"
        )
        return fig


def burn_figure(design: Design) -> Figure:
    burn = design.burn
    t = burn.time_s
    pa = design.ambient_pa
    with _style():
        fig = Figure(figsize=(7.5, 5.6), layout="constrained")
        axes = fig.subplots(2, 2, sharex=True)
        ax = axes[0, 0]
        ax.plot(t, burn.thrust_n(pa), color=SERIES[0], label="Design ambient")
        ax.plot(t, burn.thrust_vacuum_n, color=SERIES[1], label="Vacuum")
        ax.set_title("Thrust [N]")
        ax.legend(loc="lower left")
        ax.set_ylim(bottom=0)

        ax = axes[0, 1]
        ax.plot(t, pa_to_bar(burn.chamber_pressure_pa), color=SERIES[0], label="Chamber")
        same = np.allclose(burn.oxidizer_tank_pressure_pa, burn.fuel_tank_pressure_pa, rtol=1e-3)
        ax.plot(
            t,
            pa_to_bar(burn.oxidizer_tank_pressure_pa),
            color=SERIES[1],
            label="Tanks" if same else "Oxidizer tank",
        )
        if not same:
            ax.plot(t, pa_to_bar(burn.fuel_tank_pressure_pa), color=SERIES[2], label="Fuel tank")
        ax.set_title("Pressure [bar]")
        ax.legend(loc="lower left")
        ax.set_ylim(bottom=0)

        ax = axes[1, 0]
        ax.plot(t, burn.mixture_ratio, color=SERIES[0])
        ax.axhline(design.mixture_ratio, color=MUTED, linewidth=0.8)
        ax.set_title("Mixture ratio O/F [-]")
        ax.set_xlabel("Time [s]")
        span = max(0.05 * design.mixture_ratio, float(np.ptp(burn.mixture_ratio)))
        ax.set_ylim(design.mixture_ratio - span, design.mixture_ratio + span)

        ax = axes[1, 1]
        if np.all(np.isnan(burn.bottle_pressure_pa)):
            ax.plot(t, burn.isp_s(pa), color=SERIES[0])
            ax.set_title("Specific impulse, design ambient [s]")
        else:
            ax.plot(t, pa_to_bar(burn.bottle_pressure_pa), color=SERIES[0])
            ax.axhline(
                pa_to_bar(design.pressurisation.regulator_dropout_pa), color=MUTED, linewidth=0.8
            )
            ax.set_title("Pressurant bottle [bar]")
            ax.set_ylim(bottom=0)
        ax.set_xlabel("Time [s]")
        return fig


def flight_figure(design: Design) -> Figure:
    fl = design.flight
    t = fl.time_s
    with _style():
        fig = Figure(figsize=(7.5, 5.6), layout="constrained")
        axes = fig.subplots(2, 2, sharex=True)
        panels = (
            (axes[0, 0], fl.altitude_m, "Altitude above pad [m]"),
            (axes[0, 1], fl.velocity_m_s, "Velocity [m/s]"),
            (axes[1, 0], fl.acceleration_m_s2 / G0_M_S2, "Acceleration [g]"),
            (axes[1, 1], fl.dynamic_pressure_pa / PA_PER_KPA, "Dynamic pressure [kPa]"),
        )
        for ax, values, title in panels:
            ax.plot(t, values, color=SERIES[0])
            ax.axvline(fl.burnout_time_s, color=MUTED, linewidth=0.8)
            ax.set_title(title)
        ax = axes[0, 0]
        ax.plot(
            [fl.apogee_time_s],
            [fl.apogee_m],
            "o",
            color=SERIES[0],
            markersize=7,
            markeredgecolor=SURFACE,
            markeredgewidth=2,
        )
        ax.annotate(
            f"apogee {fl.apogee_m:,.0f} m",
            (fl.apogee_time_s, fl.apogee_m),
            xytext=(-6, -14),
            textcoords="offset points",
            ha="right",
            fontsize=8,
            color=INK,
        )
        ax = axes[1, 1]
        ax.annotate(
            f"max q {fl.max_q_pa / PA_PER_KPA:.1f} kPa",
            (fl.max_q_time_s, fl.max_q_pa / PA_PER_KPA),
            xytext=(6, -4),
            textcoords="offset points",
            fontsize=8,
            color=INK,
        )
        for ax in axes[1]:
            ax.set_xlabel("Time [s]")
        axes[0, 1].annotate(
            "burnout",
            (fl.burnout_time_s, 0),
            xytext=(3, 3),
            textcoords="offset points",
            fontsize=8,
            color=MUTED,
        )
        return fig


def overview_figure(design: Design) -> Figure:
    """Three-panel summary for the README: chamber contour, thrust curve, altitude."""
    burn, fl = design.burn, design.flight
    with _style():
        fig = Figure(figsize=(11.0, 3.0), layout="constrained")
        grid = fig.add_gridspec(1, 3, width_ratios=[1.35, 1, 1])
        ax = fig.add_subplot(grid[0])
        _engine_axes(ax, design, annotate=True)
        ax.set_title(f"{design.spec.name} thrust chamber")

        ax = fig.add_subplot(grid[1])
        ax.plot(burn.time_s, burn.thrust_n(design.ambient_pa), color=SERIES[0])
        ax.set_ylim(bottom=0)
        ax.set_title("Thrust at design ambient [N]")
        ax.set_xlabel("Time [s]")

        ax = fig.add_subplot(grid[2])
        ax.plot(fl.time_s, fl.altitude_m, color=SERIES[0])
        ax.plot(
            [fl.apogee_time_s],
            [fl.apogee_m],
            "o",
            color=SERIES[0],
            markersize=7,
            markeredgecolor=SURFACE,
            markeredgewidth=2,
        )
        ax.annotate(
            f"{fl.apogee_m:,.0f} m",
            (fl.apogee_time_s, fl.apogee_m),
            xytext=(-6, -14),
            textcoords="offset points",
            ha="right",
            fontsize=8,
            color=INK,
        )
        ax.set_title("Altitude above pad [m]")
        ax.set_xlabel("Time [s]")
        return fig


def vehicle_figure(design: Design) -> Figure:
    """Side view of the stack, nose left, to scale."""
    radius = m_to_mm(design.spec.structure.max_diameter_m) / 2
    blocks: list[tuple[str, str, float]] = []
    if design.spec.vehicle.other_length_m > 0:
        blocks.append(("nose", "Nose, recovery, avionics", design.spec.vehicle.other_length_m))
    if design.bottle:
        blocks.append(("bottle", "Bottle", design.bottle.length_m))
    blocks += [
        ("oxidizer", "Oxidizer", design.oxidizer_tank.length_m),
        ("fuel", "Fuel", design.fuel_tank.length_m),
    ]
    fills = {"nose": "#f0efec", "bottle": GRID, "oxidizer": "#b7d3f6", "fuel": "#f7c6b0"}
    with _style():
        fig = Figure(figsize=(9.0, 1.9), layout="constrained")
        ax = fig.add_subplot()
        ax.grid(False)
        x = 0.0
        for kind, label, length in blocks:
            w = m_to_mm(length)
            ax.add_patch(
                Rectangle(
                    (x, -radius),
                    w,
                    2 * radius,
                    facecolor=fills[kind],
                    edgecolor=SURFACE,
                    linewidth=2,
                )
            )
            ax.text(
                x + w / 2,
                0,
                f"{label}\n{w:.0f} mm",
                ha="center",
                va="center",
                fontsize=7.5,
                color=INK,
            )
            x += w
        ch = design.chamber
        ex, er = m_to_mm(ch.x_m - ch.x_m[0]) + x, m_to_mm(ch.r_m)
        ax.fill_between(ex, -er, er, color=WALL, linewidth=0)
        ax.plot(ex, er, ex, -er, color=SERIES[0], linewidth=1.2)
        ax.set_xlim(-20, ex[-1] + 20)
        ax.set_ylim(-radius * 1.6, radius * 1.6)
        ax.set_aspect("equal")
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.set_xlabel("Length [mm]")
        ax.set_title(
            f"Vehicle stack: {design.vehicle_length_m:.2f} m long, "
            f"{m_to_mm(design.spec.structure.max_diameter_m):.0f} mm diameter, "
            f"{design.liftoff_mass_kg:.1f} kg at lift-off"
        )
        return fig
