from __future__ import annotations

import html
import math

import numpy as np

from limerocket import __version__, atmosphere, thermo
from limerocket.design import Design
from limerocket.feed import FeedMode
from limerocket.plots import burn_figure, engine_figure, flight_figure, to_svg, vehicle_figure
from limerocket.tanks import Vessel
from limerocket.units import (
    G0_M_S2,
    KG_PER_G,
    M3_PER_CM3,
    M3_PER_L,
    M_PER_UM,
    PA_PER_KPA,
    PA_PER_MPA,
    m_to_mm,
    pa_to_bar,
)

Row = tuple[str, str, str]


def fmt(value: float, digits: int = 4) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "-"
    if value == 0:
        return "0"
    magnitude = math.floor(math.log10(abs(value)))
    decimals = max(digits - 1 - magnitude, 0)
    return f"{value:,.{decimals}f}"


def _bar(pa: float) -> str:
    return fmt(pa_to_bar(pa))


def _mm(m: float) -> str:
    return fmt(m_to_mm(m))


def _table(rows: list[Row], head: tuple[str, ...] = ("Quantity", "Value", "Unit")) -> str:
    cells = "".join(
        f"<tr><td>{html.escape(a)}</td><td class=num>{html.escape(b)}</td>"
        f"<td>{html.escape(c)}</td></tr>"
        for a, b, c in rows
    )
    header = "".join(f"<th>{h}</th>" for h in head)
    return f"<table><thead><tr>{header}</tr></thead><tbody>{cells}</tbody></table>"


def _grid(headers: list[str], rows: list[list[str]]) -> str:
    head = "".join(f"<th>{html.escape(h)}</th>" for h in headers)
    body = "".join(
        "<tr>"
        + "".join(
            f"<td{' class=num' if i else ''}>{html.escape(c)}</td>" for i, c in enumerate(row)
        )
        + "</tr>"
        for row in rows
    )
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def input_rows(design: Design) -> list[Row]:
    s = design.spec
    e, f, st, v = s.engine, s.feed, s.structure, s.vehicle
    rows: list[Row] = [
        ("Target thrust at design ambient", fmt(e.thrust_n), "N"),
        (
            "Burn time",
            fmt(e.burn_time_s),
            "s, set by the apogee search" if design.apogee_search else "s",
        ),
        ("Chamber pressure", _bar(e.chamber_pressure_pa), "bar"),
        ("Design ambient pressure", _bar(e.design_ambient_pressure_pa), "bar"),
        ("Oxidizer", e.oxidizer.name, ""),
        ("Fuel", e.fuel.name, ""),
        (
            "Mixture ratio",
            "optimise for Isp" if e.mixture_ratio is None else fmt(e.mixture_ratio),
            "",
        ),
        ("Characteristic length L*", fmt(e.l_star_m), "m"),
        (
            "Contraction ratio",
            "auto" if e.contraction_ratio is None else fmt(e.contraction_ratio),
            "",
        ),
        ("Convergence half-angle", fmt(e.convergence_half_angle_deg), "deg"),
        ("c* efficiency", fmt(e.c_star_efficiency), ""),
        ("Thrust coefficient efficiency", fmt(e.cf_efficiency), ""),
        ("Chamber wall material", e.wall_material.name, ""),
        (
            "Injector pressure drop (minimum)",
            fmt(100 * e.injector.pressure_drop_fraction),
            "% of Pc",
        ),
        ("Feed mode", f.mode.value, ""),
        ("Pressurant", f.pressurant.name, ""),
        ("Pressurant temperature", fmt(f.pressurant_temperature_k), "K"),
        ("Polytropic exponent (expansion)", fmt(f.polytropic_exponent), ""),
        ("Unusable residual", fmt(100 * f.residual_fraction), "% of usable"),
    ]
    if f.mode is FeedMode.REGULATED:
        rows += [
            ("Bottle material", st.bottle_material.name, ""),
            ("Bottle pressure", _bar(f.bottle_pressure_pa), "bar"),
            ("Regulator dropout margin", _bar(f.regulator_dropout_margin_pa), "bar"),
            ("Pressurant margin", fmt(100 * f.pressurant_margin), "%"),
            ("Initial ullage", fmt(100 * f.ullage_fraction), "% of tank"),
        ]
    else:
        rows.append(("Blowdown ratio", fmt(f.blowdown_ratio), ""))
    for label, line in (("Oxidizer line", f.oxidizer_line), ("Fuel line", f.fuel_line)):
        rows.append(
            (
                label,
                (
                    f"{fmt(line.length_m)} m x {_mm(line.inner_diameter_m)} mm ID, "
                    f"k = {fmt(line.minor_loss_k)}, roughness {fmt(line.roughness_m / M_PER_UM)} um"
                ),
                "",
            )
        )
    rows += [
        ("Maximum diameter", _mm(st.max_diameter_m), "mm"),
        ("Tank material", st.tank_material.name, ""),
        (
            "Safety factor on yield / ultimate",
            f"{fmt(st.yield_factor)} / {fmt(st.ultimate_factor)}",
            "",
        ),
        ("Minimum wall", _mm(st.min_wall_m), "mm"),
        ("Other masses", fmt(sum(m for _, m in v.other_masses_kg)), "kg"),
        ("Other length", fmt(v.other_length_m), "m"),
        ("Launch rail", fmt(v.rail_length_m), "m"),
        ("Launch site altitude", fmt(v.launch_altitude_m), "m"),
        ("Apogee target", "-" if v.apogee_target_m is None else fmt(v.apogee_target_m), "m"),
    ]
    return rows


def performance_rows(design: Design) -> list[Row]:
    p = design.ideal
    pa = design.ambient_pa
    return [
        (
            "Mixture ratio O/F",
            fmt(design.mixture_ratio),
            "optimised" if design.spec.engine.mixture_ratio is None else "given",
        ),
        ("Chamber temperature", fmt(p.chamber.temperature_k), "K"),
        ("Chamber molar mass", fmt(p.chamber.molar_mass_kg_mol / KG_PER_G), "kg/kmol"),
        ("Chamber isentropic exponent", fmt(p.chamber.gamma), ""),
        (
            "Throat temperature / pressure",
            f"{fmt(p.throat.temperature_k)} / {_bar(p.throat.pressure_pa)}",
            "K / bar",
        ),
        (
            "Exit temperature / pressure",
            f"{fmt(p.exit.temperature_k)} / {_bar(p.exit.pressure_pa)}",
            "K / bar",
        ),
        ("Exit Mach number", fmt(p.exit.mach), ""),
        ("Area ratio Ae/At (optimum expansion)", fmt(p.area_ratio), ""),
        ("Bell divergence factor (1 + cos theta_e) / 2", fmt(design.divergence_efficiency), ""),
        ("c*, ideal / delivered", f"{fmt(p.c_star_m_s)} / {fmt(design.c_star_m_s)}", "m/s"),
        ("Cf vacuum, ideal / delivered", f"{fmt(p.cf_vacuum)} / {fmt(design.cf_vacuum)}", ""),
        ("Cf at design ambient, ideal / delivered", f"{fmt(p.cf(pa))} / {fmt(design.cf)}", ""),
        (
            "Isp at design ambient, ideal / delivered",
            f"{fmt(p.isp_s(pa))} / {fmt(design.isp_s)}",
            "s",
        ),
        (
            "Isp vacuum, ideal / delivered",
            f"{fmt(p.isp_vacuum_s)} / {fmt(design.isp_vacuum_s)}",
            "s",
        ),
        ("Total mass flow", fmt(design.mass_flow_kg_s), "kg/s"),
        (
            "Oxidizer / fuel flow",
            f"{fmt(design.oxidizer_flow_kg_s)} / {fmt(design.fuel_flow_kg_s)}",
            "kg/s",
        ),
    ]


def ideal_gas_rows(design: Design) -> list[list[str]]:
    p, check = design.ideal, design.ideal_gas
    pairs = (
        ("c* [m/s]", p.c_star_m_s, check.c_star_m_s),
        ("Exit Mach number", p.exit.mach, check.exit_mach),
        (
            "Exit pressure / chamber pressure",
            p.exit_pressure_pa / p.chamber_pressure_pa,
            check.exit_pressure_ratio,
        ),
        ("Cf vacuum", p.cf_vacuum, check.cf_vacuum),
    )
    return [
        [name, fmt(cea), fmt(ideal), f"{100 * (ideal / cea - 1):+.1f} %"]
        for name, cea, ideal in pairs
    ]


def chamber_rows(design: Design) -> list[Row]:
    c = design.chamber
    return [
        ("Throat diameter", _mm(2 * c.throat_radius_m), "mm"),
        ("Chamber diameter", _mm(2 * c.chamber_radius_m), "mm"),
        ("Exit diameter", _mm(2 * c.exit_radius_m), "mm"),
        (
            "Contraction ratio Ac/At",
            fmt(c.contraction_ratio),
            "auto" if design.spec.engine.contraction_ratio is None else "given",
        ),
        ("Chamber volume (L* x At)", fmt(c.chamber_volume_m3 / M3_PER_CM3), "cm3"),
        ("Cylinder length", _mm(c.cylinder_length_m), "mm"),
        ("Convergent length", _mm(c.convergent_length_m), "mm"),
        ("Bell length (80 % of 15 deg cone)", _mm(c.nozzle.length_m), "mm"),
        ("Injector face to exit", _mm(c.length_m), "mm"),
        (
            "Bell angles theta_n / theta_e",
            f"{fmt(c.nozzle.theta_n_deg, 3)} / {fmt(c.nozzle.theta_e_deg, 3)}",
            "deg",
        ),
        ("Wall thickness (hoop at Pc)", _mm(c.wall_thickness_m), "mm"),
        ("Injector plate thickness", _mm(c.injector_plate_thickness_m), "mm"),
        ("Shell / plate mass", f"{fmt(c.shell_mass_kg)} / {fmt(c.injector_plate_mass_kg)}", "kg"),
    ]


def feed_rows(design: Design) -> list[list[str]]:
    pc = design.spec.engine.chamber_pressure_pa
    rows = []
    for label, side, load, line_dp in (
        ("Oxidizer", design.injector.oxidizer, design.oxidizer_load, design.oxidizer_line_drop_pa),
        ("Fuel", design.injector.fuel, design.fuel_load, design.fuel_line_drop_pa),
    ):
        branch = load.branch
        flow = side.mass_flow_kg_s
        rows.append(
            [
                label,
                fmt(flow),
                fmt(branch.line.velocity_m_s(flow, branch.fluid), 3),
                fmt(branch.line.reynolds(flow, branch.fluid), 3),
                fmt(branch.line.friction_factor(flow, branch.fluid), 3),
                _bar(line_dp),
                f"{_bar(side.pressure_drop_pa)} ({fmt(100 * side.pressure_drop_pa / pc, 3)} %)",
                f"{side.orifices.count} x {fmt(m_to_mm(side.orifices.diameter_m), 3)} mm",
                fmt(side.orifices.discharge_coefficient, 2),
                fmt(side.velocity_m_s, 3),
            ]
        )
    return rows


def pressurisation_rows(design: Design) -> list[Row]:
    p = design.pressurisation
    rows: list[Row] = [
        (
            "Tank pressure (regulated set point)"
            if p.mode is FeedMode.REGULATED
            else "Initial tank pressure",
            _bar(p.tank_pressure_pa),
            "bar",
        ),
        ("Injector momentum ratio (ox/fuel)", fmt(design.injector.momentum_ratio, 3), ""),
        ("Pressurant mass", fmt(design.pressurant_mass_kg / KG_PER_G), "g"),
    ]
    if p.mode is FeedMode.REGULATED:
        rows += [
            ("Bottle volume", fmt(p.bottle_volume_m3 / M3_PER_L), "L"),
            ("Bottle pressure", _bar(p.bottle_pressure_pa), "bar"),
            ("Regulator dropout below", _bar(p.regulator_dropout_pa), "bar"),
        ]
    return rows


def vessel_rows(vessels: list[Vessel]) -> list[list[str]]:
    return [
        [
            v.name,
            fmt(v.volume_m3 / M3_PER_L),
            "sphere" if v.is_sphere else "cylinder",
            _mm(2 * v.inner_radius_m),
            _mm(v.length_m),
            "-" if v.is_sphere else _mm(v.cylinder_wall_m),
            _mm(v.head_wall_m),
            _bar(v.pressure_pa),
            fmt(max(v.hoop_stress_pa, v.head_stress_pa) / PA_PER_MPA, 3),
            fmt(v.allowable_stress_pa / PA_PER_MPA, 3),
            fmt(v.margin_of_safety, 2),
            fmt(v.mass_kg),
        ]
        for v in vessels
    ]


def mass_rows(design: Design) -> list[Row]:
    rows: list[Row] = [(m.name, fmt(m.mass_kg), m.group) for m in design.masses]
    rows += [
        ("Dry mass", fmt(design.dry_mass_kg), "total"),
        ("Lift-off mass", fmt(design.liftoff_mass_kg), "total"),
        (
            "Propellant mass fraction",
            fmt(design.propellant_mass_kg / design.liftoff_mass_kg, 3),
            "",
        ),
    ]
    return rows


def burn_rows(design: Design) -> list[Row]:
    b = design.burn
    pa = design.ambient_pa
    thrust = b.thrust_n(pa)
    impulse = b.total_impulse_ns(pa)
    rows: list[Row] = [
        ("Burn time", fmt(b.burn_time_s), "s"),
        ("Total impulse (design ambient)", fmt(impulse), "N s"),
        ("Average thrust", fmt(impulse / b.burn_time_s), "N"),
        ("Thrust start / end", f"{fmt(thrust[0])} / {fmt(thrust[-1])}", "N"),
        (
            "Chamber pressure start / end",
            f"{_bar(b.chamber_pressure_pa[0])} / {_bar(b.chamber_pressure_pa[-1])}",
            "bar",
        ),
        ("O/F start / end", f"{fmt(b.mixture_ratio[0])} / {fmt(b.mixture_ratio[-1])}", ""),
        ("Average delivered Isp", fmt(impulse / (b.expelled_mass_kg * G0_M_S2)), "s"),
        ("Propellant expelled", fmt(b.expelled_mass_kg), "kg"),
        ("First tank to empty", b.depleted, ""),
        ("Regulator dropout", "-" if b.dropout_time_s is None else fmt(b.dropout_time_s), "s"),
        (
            "Minimum exit / site ambient pressure",
            fmt(
                float(np.min(b.exit_pressure_pa))
                / atmosphere.state(design.spec.vehicle.launch_altitude_m).pressure_pa,
                3,
            ),
            "",
        ),
    ]
    return rows


def flight_rows(design: Design) -> list[Row]:
    f = design.flight
    rows: list[Row] = [
        ("Apogee above pad", fmt(f.apogee_m), "m"),
        ("Time to apogee", fmt(f.apogee_time_s), "s"),
        (
            "Burnout velocity / altitude",
            f"{fmt(f.burnout_velocity_m_s)} / {fmt(f.burnout_altitude_m)}",
            "m/s / m",
        ),
        ("Maximum Mach number", fmt(f.max_mach, 3), ""),
        (
            "Maximum dynamic pressure",
            fmt(f.max_q_pa / PA_PER_KPA),
            f"kPa at {fmt(f.max_q_time_s, 3)} s",
        ),
        ("Maximum acceleration", fmt(f.max_acceleration_m_s2 / G0_M_S2, 3), "g"),
        ("Lift-off thrust-to-weight", fmt(f.liftoff_thrust_to_weight, 3), ""),
        ("Rail exit velocity", fmt(f.rail_exit_velocity_m_s, 3), "m/s"),
        ("Vehicle length", fmt(design.vehicle_length_m), "m"),
    ]
    search = design.apogee_search
    if search:
        rows.append(
            (
                "Apogee search",
                (
                    f"target {fmt(search.target_m)} m, burn time {fmt(search.burn_time_s)} s, "
                    f"{'reached' if search.reached else 'not reachable'} "
                    f"after {search.evaluations} designs"
                ),
                "",
            )
        )
    return rows


CSS = """
:root { color-scheme: light; --ink:#0b0b0b; --ink2:#52514e; --muted:#898781; --line:#e1e0d9;
  --surface:#fcfcfb; --page:#f9f9f7; --accent:#2a78d6; --warn:#fab219; }
* { box-sizing: border-box; }
body { margin:0; background:var(--page); color:var(--ink);
  font:15px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif; }
main { max-width: 900px; margin: 0 auto; padding: 32px 16px 64px; }
h1 { font-size: 28px; margin: 0 0 4px; }
h2 { font-size: 19px; margin: 40px 0 12px; padding-top: 8px; border-top: 1px solid var(--line); }
p.sub { color: var(--ink2); margin: 0 0 20px; }
.tiles { display:grid; grid-template-columns: repeat(auto-fit, minmax(130px,1fr)); gap: 10px; }
.tile { background: var(--surface); border: 1px solid var(--line); border-radius: 8px;
  padding: 10px 12px; }
.tile .label { color: var(--ink2); font-size: 12px; }
.tile .value { font-size: 22px; font-weight: 600; }
.tile .unit { color: var(--muted); font-size: 13px; margin-left: 3px; font-weight: 400; }
.warn { background: #fff7e6; border-left: 4px solid var(--warn); padding: 10px 14px;
  margin: 16px 0; border-radius: 4px; }
.warn ul { margin: 4px 0 0; padding-left: 18px; }
table { width: 100%; border-collapse: collapse; background: var(--surface); margin: 8px 0 16px;
  font-size: 13.5px; }
th, td { text-align: left; padding: 5px 10px; border-bottom: 1px solid var(--line);
  vertical-align: top; }
th { color: var(--ink2); font-weight: 600; font-size: 12.5px; }
td.num { font-variant-numeric: tabular-nums; white-space: nowrap; }
.scroll { overflow-x: auto; }
figure { margin: 12px 0 20px; background: var(--surface); border: 1px solid var(--line);
  border-radius: 8px; padding: 8px; }
figure svg { width: 100%; height: auto; display:block; }
code { font-size: 13px; }
footer { color: var(--muted); font-size: 12.5px; margin-top: 40px; }
"""

ASSUMPTIONS = (
    (
        "Chemistry: NASA CEA shifting equilibrium, infinite-area combustor, ideal expansion; "
        "real losses enter only through the c* and Cf efficiencies."
    ),
    (
        "Thrust chamber walls are sized for hoop stress at room-temperature allowables; there is "
        "no heat-transfer or cooling analysis."
    ),
    (
        "Feed: incompressible liquids, quasi-steady orifice and Darcy-Weisbach line flow, "
        "ideal-gas pressurant, adiabatic tanks and bottle, no start/stop transients, "
        "no acceleration head."
    ),
    (
        "Flight: vertical 1-DOF point mass, U.S. Standard Atmosphere 1976, drag from the Cd(Mach) "
        "table, no wind, no recovery phase."
    ),
)


def html_report(design: Design, files: list[str]) -> str:
    s = design.spec
    f = design.flight
    tiles = [
        ("Thrust (design ambient)", fmt(design.thrust_n), "N"),
        ("Delivered Isp", fmt(design.isp_s, 3), "s"),
        ("Burn time", fmt(design.burn.burn_time_s, 3), "s"),
        ("Lift-off mass", fmt(design.liftoff_mass_kg, 3), "kg"),
        ("Apogee", fmt(f.apogee_m, 3), "m"),
        ("Max Mach", fmt(f.max_mach, 3), ""),
    ]
    tile_html = "".join(
        f"<div class=tile><div class=label>{label}</div><div class=value>{value}"
        f"<span class=unit>{unit}</span></div></div>"
        for label, value, unit in tiles
    )
    warnings = ""
    if design.warnings:
        items = "".join(f"<li>{html.escape(w)}</li>" for w in design.warnings)
        warnings = f"<div class=warn><strong>Check these</strong><ul>{items}</ul></div>"
    vessels = [v for v in (design.oxidizer_tank, design.fuel_tank, design.bottle) if v]
    file_items = "".join(
        f"<li><a href='{html.escape(n)}'><code>{html.escape(n)}</code></a></li>" for n in files
    )
    assumptions = "".join(f"<li>{html.escape(a)}</li>" for a in ASSUMPTIONS)

    def figure(svg: str) -> str:
        return f"<figure>{svg[svg.index('<svg') :]}</figure>"

    sections = [
        f"<h1>{html.escape(s.name)}</h1>",
        f"<p class=sub>{html.escape(s.description)}</p>" if s.description else "",
        f"<div class=tiles>{tile_html}</div>",
        warnings,
        figure(to_svg(vehicle_figure(design))),
        "<h2>1. Inputs</h2>",
        _table(input_rows(design)),
        "<h2>2. Propellant performance</h2>",
        (
            f"<p>Backend: {html.escape(thermo.BACKEND)}. Ideal values are CEA's; delivered "
            "values apply the c* and Cf efficiencies and the bell's divergence factor.</p>"
        ),
        _table(performance_rows(design)),
        (
            "<p>Cross-check: the ideal-rocket relations (Sutton &amp; Biblarz ch. 3) at the "
            "chamber's gamma and molar mass, held constant. The gap is what shifting "
            "equilibrium and a changing gamma add.</p>"
        ),
        _grid(["Quantity", "CEA", "Ideal gas", "Difference"], ideal_gas_rows(design)),
        "<h2>3. Thrust chamber and nozzle</h2>",
        figure(to_svg(engine_figure(design))),
        _table(chamber_rows(design)),
        "<h2>4. Injector and feed system</h2>",
        "<div class=scroll>"
        + _grid(
            [
                "Branch",
                "Flow [kg/s]",
                "Line v [m/s]",
                "Line Re",
                "Darcy f",
                "Line dp [bar]",
                "Injector dp [bar]",
                "Orifices",
                "Cd",
                "Jet v [m/s]",
            ],
            feed_rows(design),
        )
        + "</div>",
        _table(pressurisation_rows(design)),
        "<h2>5. Tanks and pressure vessels</h2>",
        "<div class=scroll>"
        + _grid(
            [
                "Vessel",
                "Volume [L]",
                "Form",
                "ID [mm]",
                "Length [mm]",
                "Cyl. wall [mm]",
                "Head wall [mm]",
                "MEOP [bar]",
                "Peak stress [MPa]",
                "Allowable [MPa]",
                "Margin",
                "Mass [kg]",
            ],
            vessel_rows(vessels),
        )
        + "</div>",
        (
            "<p>Allowable stress is min(yield / SF<sub>y</sub>, ultimate / SF<sub>u</sub>); "
            "a margin of 0 means the wall sits exactly at the allowable.</p>"
        ),
        "<h2>6. Mass budget</h2>",
        _table(mass_rows(design), ("Item", "Mass [kg]", "Group")),
        "<h2>7. Burn</h2>",
        figure(to_svg(burn_figure(design))),
        _table(burn_rows(design)),
        "<h2>8. Trajectory</h2>",
        figure(to_svg(flight_figure(design))),
        _table(flight_rows(design)),
        "<h2>9. Files</h2>",
        f"<ul>{file_items}</ul>",
        "<h2>10. Assumptions</h2>",
        f"<ul>{assumptions}</ul>",
        f"<footer>Generated by LimeRocket {__version__}</footer>",
    ]
    body = "\n".join(x for x in sections if x)
    return (
        "<!doctype html><html lang=en><head><meta charset=utf-8>"
        "<meta name=viewport content='width=device-width,initial-scale=1'>"
        f"<title>{html.escape(s.name)} design report</title><style>{CSS}</style></head>"
        f"<body><main>{body}</main></body></html>\n"
    )


def terminal_summary(design: Design) -> str:
    s = design.spec
    b, f, c = design.burn, design.flight, design.chamber
    pa = design.ambient_pa
    lines = [
        f"{s.name}: {s.engine.oxidizer.name} / {s.engine.fuel.name}, {s.feed.mode.value} feed",
        "",
        (
            f"  O/F {fmt(design.mixture_ratio)}   Pc {_bar(s.engine.chamber_pressure_pa)} bar   "
            f"Tc {fmt(design.ideal.chamber.temperature_k)} K   Ae/At {fmt(design.ideal.area_ratio)}"
        ),
        (
            f"  Isp {fmt(design.isp_s)} s (vac {fmt(design.isp_vacuum_s)} s)   "
            f"mdot {fmt(design.mass_flow_kg_s)} kg/s   c* {fmt(design.c_star_m_s)} m/s"
        ),
        (
            f"  Throat {_mm(2 * c.throat_radius_m)} mm   "
            f"chamber {_mm(2 * c.chamber_radius_m)} mm   "
            f"exit {_mm(2 * c.exit_radius_m)} mm   length {_mm(c.length_m)} mm"
        ),
        (
            f"  Tank pressure {_bar(design.pressurisation.tank_pressure_pa)} bar   "
            f"injector dp ox/fuel {_bar(design.injector.oxidizer.pressure_drop_pa)}/"
            f"{_bar(design.injector.fuel.pressure_drop_pa)} bar"
        ),
        (
            f"  Burn {fmt(b.burn_time_s)} s   impulse {fmt(b.total_impulse_ns(pa))} N s   "
            f"thrust {fmt(b.thrust_n(pa)[0])} -> {fmt(b.thrust_n(pa)[-1])} N"
        ),
        (
            f"  Mass dry {fmt(design.dry_mass_kg)} kg   "
            f"lift-off {fmt(design.liftoff_mass_kg)} kg   "
            f"length {fmt(design.vehicle_length_m)} m"
        ),
        (
            f"  Apogee {fmt(f.apogee_m)} m   max Mach {fmt(f.max_mach, 3)}   "
            f"max q {fmt(f.max_q_pa / PA_PER_KPA)} kPa   "
            f"max accel {fmt(f.max_acceleration_m_s2 / G0_M_S2, 3)} g"
        ),
    ]
    if design.warnings:
        lines.append("")
        lines += [f"  ! {w}" for w in design.warnings]
    return "\n".join(lines)
