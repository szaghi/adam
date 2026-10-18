#!/usr/bin/env python3
"""Write the inputs of the FLUME radiative-cooling verification (issue #93) from the sod-x regression input.

    make_co.py <out.ini> [--case box|shock] [--model euler|euler-de|mhd-ideal] [--u U] [--p P] [--ni N]
                         [--time-max T] [--it-save N] [--cfl C] [--dt-factor F]
                         [--temperatures "T1 ... TN"] [--lambdas "L1 ... LN"] [--particle-mass M|none]

The case is a uniform box: the sod-x grid without its initial refinement (50 cells along x, the other directions
null, extrapolation boundaries), both regions set to rho = 1, u = U (default 0), p = P (default 1), gamma = 5/3 (so
R = 1 and T = p / rho = P), a conservation history row every step, and a [cooling] section with the given table
(default the four nodes T = 0.1 0.3 1 3, Lambda = 0.2 1 0.5 1.5: segment exponents 1.46, -0.58 and 1, the last one
continued above T = 3, the floor at T = 0.1). A uniform state has no residual, so every cell follows the isochoric
cooling ODE alone. `--particle-mass none` leaves the key out (the contract leg).

`--case shock` is the radiative shock: two cold streams rho = 1, u = +-U (default 1) at p = P (default 0.1, T = 0.1,
Mach 2.45) colliding at x = 0.5 on 400 cells (`--ni`), to t = 3 with a checkpoint every 200 steps, and the power law
Lambda = 0.18 (T / 0.1)^1.5 (the nodes T = 0.1 10, Lambda = 0.18 180): the upstream gas sits on the floor, the shocked
gas cools back to it over a layer of about 20 cells and piles into a dense layer at rest at the centre. The exponent
1.5 keeps the planar radiative shock clear of its overstability (Chevalier & Imamura 1982).
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

BASE = Path(__file__).resolve().parents[2] / "regression" / "sod-x" / "input.ini"
TEMPERATURES = "0.1 0.3 1.0 3.0"
LAMBDAS = "0.2 1.0 0.5 1.5"
# the radiative shock (co7): defaults of the case, each overridable
SHOCK = {"u": 1.0, "p": 0.1, "ni": 400, "time_max": 3.0, "it_save": 200,
         "temperatures": "0.1 10.0", "lambdas": "0.18 180.0"}
BOX = {"u": 0.0, "p": 1.0, "ni": 50, "time_max": 4.0, "it_save": 100000,
       "temperatures": TEMPERATURES, "lambdas": LAMBDAS}


def case_ini(a: argparse.Namespace) -> str:
    """Return the INI of a case."""
    s = BASE.read_text()
    s = re.sub(r"(?m)^physical_model\s*=.*$", f"physical_model = {a.model}", s)
    s = re.sub(r"(?m)^cp\s*=.*\n^cv\s*=.*\n", "gamma = 1.6666666666666667\n", s)
    regions = s.split("[initial_conditions_region_")
    for k in (1, 2):
        r = regions[k]
        u = -a.u if (a.case == "shock" and k == 2) else a.u  # the shock: the streams collide at x = 0.5
        r = re.sub(r"(?m)^r\s*=.*$", "r      = 1.0", r)
        r = re.sub(r"(?m)^u\s*=.*$", f"u      = {u!r}", r)
        r = re.sub(r"(?m)^p\s*=.*$", f"p      = {a.p!r}", r)
        regions[k] = r
    s = "[initial_conditions_region_".join(regions)
    s = re.sub(r"(?m)^(\s*time_max\s*=).*$", rf"\1 {a.time_max!r}", s)
    s = re.sub(r"(?m)^(\s*CFL\s*=).*$", rf"\1 {a.cfl!r}", s)
    s = re.sub(r"(?m)^(\s*it_save\s*=).*$", rf"\1 {a.it_save}", s)
    s = re.sub(r"(?m)^(\s*iu_ref_levels\s*=).*$", r"\1 0", s)  # the base cells only: no initial refinement
    s = re.sub(r"(?m)^(ni\s*=).*$", rf"\1 {a.ni}", s)
    if a.model == "mhd-ideal":
        s += "\n[mhd]\ndivergence_control = none\ndivb_tol = 0.0\ndivb_error = .false.\n"
        s += "rho_floor = 0.0\np_floor = 0.0\n"
    cooling = f"\n[cooling]\ntemperatures  = {a.temperatures}\nlambdas       = {a.lambdas}\n"
    if a.particle_mass != "none":
        cooling += f"particle_mass = {a.particle_mass}\n"
    if a.dt_factor is not None:
        cooling += f"dt_factor     = {a.dt_factor!r}\n"
    return s + cooling


def main() -> int:
    """Write the requested case."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("out", type=Path)
    ap.add_argument("--case", choices=("box", "shock"), default="box")
    ap.add_argument("--model", choices=("euler", "euler-de", "mhd-ideal"), default="euler")
    ap.add_argument("--u", type=float, default=None, help="the velocity along x (shock: the stream speed)")
    ap.add_argument("--p", type=float, default=None, help="the pressure (T = p, rho = 1)")
    ap.add_argument("--ni", type=int, default=None, help="cells along x")
    ap.add_argument("--time-max", type=float, default=None)
    ap.add_argument("--it-save", type=int, default=None, help="checkpoint cadence (steps)")
    ap.add_argument("--cfl", type=float, default=0.5)
    ap.add_argument("--dt-factor", type=float, default=None, help="[cooling].(dt_factor) (default: absent)")
    ap.add_argument("--temperatures", default=None)
    ap.add_argument("--lambdas", default=None)
    ap.add_argument("--particle-mass", default="1.0", help="a number, or none to leave the key out")
    a = ap.parse_args()
    for key, value in (SHOCK if a.case == "shock" else BOX).items():
        if getattr(a, key) is None:
            setattr(a, key, value)
    a.out.write_text(case_ini(a))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
