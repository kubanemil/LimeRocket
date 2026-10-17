# Spec reference

A design spec is one YAML file. Every key carries its unit in the name. Anything not marked
required has the default shown. `auto` (or leaving the key out) lets LimeRocket choose.
Unknown keys are errors, and all problems are reported together.

Run `limerocket catalog` to list the built-in propellants and materials.

## Top level

| Key | Default | Meaning |
|---|---|---|
| `name` | required | Design name, also used for output file names |
| `description` | `""` | One line shown under the report title |

## `engine`

| Key | Default | Meaning |
|---|---|---|
| `thrust_n` | required | Thrust at the design ambient pressure, at the start of the burn |
| `burn_time_s` | required | Burn time; the propellant load is sized to it |
| `chamber_pressure_bar` | required | Chamber (injector-end stagnation) pressure |
| `design_ambient_pressure_bar` | `1.01325` | The nozzle is optimally expanded to this pressure |
| `oxidizer` | required | Catalog key (`lox`, `htp90`) or a custom propellant, see below |
| `fuel` | required | Catalog key (`ethanol`, `ethanol75`, `ipa`, `rp1`, `methane`) or custom |
| `mixture_ratio` | `optimize` | O/F by mass, or `optimize` for maximum Isp at the design point |
| `l_star_m` | `1.0` | Characteristic length, chamber volume / throat area |
| `contraction_ratio` | `auto` | Ac/At; `auto` uses Ac/At = 8.0 Dt^-0.6 + 1.25 (Dt in cm) |
| `convergence_half_angle_deg` | `30` | Convergent cone half-angle, 15 to 60 |
| `c_star_efficiency` | `0.92` | Delivered c* / ideal c* |
| `thrust_coefficient_efficiency` | `0.98` | Delivered vacuum Cf / ideal vacuum Cf |
| `wall_material` | `ss316` | Chamber and nozzle wall material |
| `injector.pressure_drop_fraction` | `0.20` | Minimum injector drop as a fraction of Pc |
| `injector.oxidizer.discharge_coefficient` | `0.7` | Orifice Cd |
| `injector.oxidizer.orifices` | `8` | Number of oxidizer orifices |
| `injector.fuel.*` | same | As for the oxidizer |

A custom propellant is a mapping of CEA species to mass fractions plus liquid properties:

```yaml
fuel:
  name: E90
  ingredients: {"C2H5OH(L)": 0.9, "H2O(L)": 0.1}
  temperature_k: 298.15
  density_kg_m3: 818
  viscosity_pa_s: 1.4e-3
```

## `feed`

| Key | Default | Meaning |
|---|---|---|
| `mode` | `regulated` | `regulated` (bottle and regulator) or `blowdown` (ullage only) |
| `pressurant` | `helium` | `helium` or `nitrogen` |
| `pressurant_temperature_k` | `293.15` | Initial gas temperature |
| `residual_fraction` | `0.03` | Unusable propellant, as a fraction of the usable load |
| `lines.oxidizer.length_m` | required | Tank-to-injector line length |
| `lines.oxidizer.inner_diameter_mm` | required | Line inner diameter |
| `lines.oxidizer.roughness_um` | `1.5` | Absolute roughness (drawn tubing ~1.5) |
| `lines.oxidizer.minor_loss_k` | `0` | Sum of fitting and valve loss coefficients |
| `lines.fuel.*` | same | As for the oxidizer |

Regulated mode only:

| Key | Default | Meaning |
|---|---|---|
| `bottle_pressure_bar` | `300` | Initial bottle pressure |
| `bottle_volume_l` | `auto` | Bottle volume; `auto` sizes it from an energy balance |
| `regulator_dropout_margin_bar` | `10` | Bottle must stay this far above the set pressure |
| `pressurant_margin` | `0.20` | Extra bottle volume over the exact ideal-gas requirement |
| `ullage_fraction` | `0.05` | Initial gas volume in each propellant tank |

Blowdown mode only:

| Key | Default | Meaning |
|---|---|---|
| `blowdown_ratio` | `2.5` | Initial over final tank pressure |
| `polytropic_exponent` | `auto` | Ullage expansion exponent, 1 (isothermal) to gamma; `auto` is the gas's gamma (adiabatic) |

In regulated mode the bottle and ullage are adiabatic. Below the dropout pressure the
regulator is wide open and the tanks trail the bottle by the dropout margin, so tank pressure
falls smoothly from the set point.

## `structure`

| Key | Default | Meaning |
|---|---|---|
| `max_diameter_mm` | required | Airframe diameter; tanks get 90 % of it |
| `tank_material` | `al6061-t6` | Propellant tank material |
| `bottle_material` | `ti6al4v` | Pressurant bottle material |
| `yield_safety_factor` | `1.5` | Allowable = min(yield / this, ultimate / next) |
| `ultimate_safety_factor` | `2.0` | |
| `min_wall_mm` | `1.0` | Minimum gauge for every wall |

## `vehicle`

| Key | Default | Meaning |
|---|---|---|
| `other_masses_kg` | none | Mapping of item to mass: airframe, recovery, avionics, plumbing |
| `other_length_m` | `0` | Length of everything that is not tanks or engine |
| `drag_coefficient` | generic table | A number, or `[[mach, cd], ...]` with increasing Mach |
| `rail_length_m` | `6` | Launch rail length |
| `launch_altitude_m` | `0` | Launch site altitude |
| `apogee_target_m` | none | If set, the burn time is searched to reach it |
