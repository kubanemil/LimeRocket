# LimeRocket

LimeRocket turns one YAML file into a first-cut design for a small pressure-fed liquid rocket.
You give it thrust, burn time, propellants, chamber pressure, feed lines and airframe diameter.
It gets the combustion from NASA CEA, draws the chamber and an 80 % Rao bell, sizes the
injector, tanks and pressurant, simulates the burn and flies the result to apogee.

![Peroxide-B: chamber contour, blowdown thrust curve and altitude](examples/out/peroxide-blowdown/overview.svg)

Peroxide-B is one of the two examples: a 2 kN HTP-90/RP-1 blowdown rocket whose burn time
LimeRocket searched to reach 5 km. Its full report is
[`examples/out/peroxide-blowdown/report.html`](examples/out/peroxide-blowdown/report.html).

```sh
git clone https://github.com/kubanemil/LimeRocket && cd LimeRocket
uv run limerocket design examples/lime-1.yaml -o out
```

That prints a summary and writes `out/report.html`, the chamber contour as CSV, DXF and STL,
and thrust curves for OpenRocket (`.eng`) and RocketPy. Every spec key, with its unit and
default, is in [docs/spec.md](docs/spec.md). To get `limerocket` on your path instead, run
`uv tool install git+https://github.com/kubanemil/LimeRocket`.

The burn is integrated in time. At each step LimeRocket finds the chamber pressure at which
the flow through the lines and injector orifices equals what the throat passes, using CEA
results mapped over O/F and pressure. Thrust decay in blowdown, O/F drift, regulator dropout
and which tank empties first all come out of that balance.

## Checked against published numbers

| Model | Reference | Worst error |
|---|---|---|
| Combustion | NASA RP-1311 examples 8 and 12, RocketCEA's LOX/LH2 example | 0.007 % |
| Isentropic flow | NACA Report 1135 tables | 0.014 % |
| Pipe friction | Moody chart | 0.12 % |
| Atmosphere | U.S. Standard Atmosphere 1976, 0 to 80 km | 0.007 % |

`./do test` runs these, plus closed-form checks such as mass conservation through the burn and
the rocket equation for a vacuum flight. `./do check` runs ruff and mypy, and `./do examples`
regenerates `examples/out`.

## What it leaves out

There is no heat transfer. Walls are sized for pressure at room-temperature strength, so a real
chamber still needs cooling or an ablative liner. Flight is vertical with one degree of freedom,
a generic drag curve and no wind; use OpenRocket or RASAero drag data for a real airframe.

MIT licensed.
