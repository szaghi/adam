#!/usr/bin/env python3
"""Write the inputs of the FLUME radiative-cooling verification (issue #93) from the sod-x regression input.

    make_co.py <out.ini> [--model euler|euler-de|mhd-ideal] [--u U] [--p P] [--time-max T] [--cfl C] [--dt-factor F]
                         [--temperatures "T1 ... TN"] [--lambdas "L1 ... LN"] [--particle-mass M|none]

The case is a uniform box: the sod-x grid without its initial refinement (50 cells along x, the other directions
null, extrapolation boundaries), both regions set to rho = 1, u = U (default 0), p = P (default 1), gamma = 5/3 (so
R = 1 and T = p / rho = P), a conservation history row every step, and a [cooling] section with the given table
(default the four nodes T = 0.1 0.3 1 3, Lambda = 0.2 1 0.5 1.5: segment exponents 1.46, -0.58 and 1, the last one
continued above T = 3, the floor at T = 0.1). A uniform state has no residual, so every cell follows the isochoric
cooling ODE alone. `--particle-mass none` leaves the key out (the contract leg).
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

BASE = Path(__file__).resolve().parents[2] / "regression" / "sod-x" / "input.ini"
TEMPERATURES = "0.1 0.3 1.0 3.0"
LAMBDAS = "0.2 1.0 0.5 1.5"


def case_ini(a: argparse.Namespace) -> str:
    """Return the INI of a case."""
    s = BASE.read_text()
    s = re.sub(r"(?m)^physical_model\s*=.*$", f"physical_model = {a.model}", s)
    s = re.sub(r"(?m)^cp\s*=.*\n^cv\s*=.*\n", "gamma = 1.6666666666666667\n", s)
    regions = s.split("[initial_conditions_region_")
    for k in (1, 2):
        r = regions[k]
        r = re.sub(r"(?m)^r\s*=.*$", "r      = 1.0", r)
        r = re.sub(r"(?m)^u\s*=.*$", f"u      = {a.u!r}", r)
        r = re.sub(r"(?m)^p\s*=.*$", f"p      = {a.p!r}", r)
        regions[k] = r
    s = "[initial_conditions_region_".join(regions)
    s = re.sub(r"(?m)^(\s*time_max\s*=).*$", rf"\1 {a.time_max!r}", s)
    s = re.sub(r"(?m)^(\s*CFL\s*=).*$", rf"\1 {a.cfl!r}", s)
    s = re.sub(r"(?m)^(\s*it_save\s*=).*$", r"\1 100000", s)
    s = re.sub(r"(?m)^(\s*iu_ref_levels\s*=).*$", r"\1 0", s)  # the 50 base cells: a uniform box needs no refinement
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
    ap.add_argument("--model", choices=("euler", "euler-de", "mhd-ideal"), default="euler")
    ap.add_argument("--u", type=float, default=0.0, help="the uniform velocity along x")
    ap.add_argument("--p", type=float, default=1.0, help="the uniform pressure (T = p, rho = 1)")
    ap.add_argument("--time-max", type=float, default=4.0)
    ap.add_argument("--cfl", type=float, default=0.5)
    ap.add_argument("--dt-factor", type=float, default=None, help="[cooling].(dt_factor) (default: absent)")
    ap.add_argument("--temperatures", default=TEMPERATURES)
    ap.add_argument("--lambdas", default=LAMBDAS)
    ap.add_argument("--particle-mass", default="1.0", help="a number, or none to leave the key out")
    a = ap.parse_args()
    a.out.write_text(case_ini(a))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
