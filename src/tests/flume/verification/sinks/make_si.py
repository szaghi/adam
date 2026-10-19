#!/usr/bin/env python3
"""Write the inputs of the FLUME sink-particle verification (issue #94) from the selfgravity-star regression input.

    make_si.py <out.ini> [--case rest|binary] [--cells N] [--rho R] [--p P] [--time-max T] [--it-max N]
                         [--x "x1 ..."] [--y ...] [--z ...] [--u ...] [--v ...] [--w ...] [--mass "m1 ..."]
                         [--restart-save N] [--restart] [--no-sinks-section] [--no-gravity]
                         [--model euler|euler-de|mhd-ideal]

The box: [-1, 1]^3, N^3 cells (default 32: blocks of N/2 cells, one uniform refinement), isolated self-gravity
(G = 1, the James step), inviscid walls, gamma = 2, the positivity limiter, a uniform background rho = R, p = P at
rest, a conservation and a sinks history row every step. Cases:
    rest    one sink of mass 1 at rest at (0.15, -0.1, 0.05), off every symmetry plane, in a background rho = 1e-3, the
            box periodic (the boundaries and the Poisson solve): no wall force, and a uniform background pulls nowhere
    binary  the equal-mass circular binary: masses 1 at (-+0.25, 0, 0), velocities (0, -+1, 0) (v^2 = G m / (2 d),
            d = 0.5: 8 cells at N = 32), period 2 pi 0.25 = 1.5708, in a near-vacuum background rho = 1e-6, p = 1e-6
    accrete  (P2) a sink of mass 0.1 at rest in a dense cold periodic background (rho = 1, p = 1e-3): it accretes
    collapse (P2) the selfgravity-star polytrope (radius 0.5, central density 1) at a tenth of its equilibrium
             pressure, kept isothermal at T = 0.01 by strong cooling (M11; adiabatic gamma = 2 never crosses the
             threshold: rho_thr grows with rho), closed walls, creation on, no initial sink: it collapses and must form
             one sink at the centre
    stable   (P2) a static periodic box, gamma = 1.4, rho = 1, p = 0.0185 (1.2 rho_thr at 2.5 cells), creation on: above
             the threshold but Jeans-stable (|E_grav| / 2 E_th = 0.66 on a control volume), no sink may form
    merge    (P2) two sinks of mass 1 at -+0.04 (1.28 cells apart) approaching at -+0.1: bound, converging, they merge
The sink keys override the case's sinks; --no-sinks-section leaves the section out (the reference runs).
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

BASE = Path(__file__).resolve().parents[2] / "regression" / "selfgravity-star" / "input.ini"
CASES = {
    "rest": {"rho": 1.0e-3, "p": 1.0e-3, "time_max": 0.5, "x": "0.15", "y": "-0.1", "z": "0.05",
             "u": "0.0", "v": "0.0", "w": "0.0", "mass": "1.0"},
    "binary": {"rho": 1.0e-6, "p": 1.0e-6, "time_max": 3.1416, "x": "-0.25 0.25", "y": "0.0 0.0", "z": "0.0 0.0",
               "u": "0.0 0.0", "v": "-1.0 1.0", "w": "0.0 0.0", "mass": "1.0 1.0"},
    # P2: accretion from a dense cold periodic background (rho_thr about 0.06 at r_acc = 2.5 cells)
    "accrete": {"rho": 1.0, "p": 1.0e-3, "time_max": 0.3, "x": "0.15", "y": "-0.1", "z": "0.05",
                "u": "0.0", "v": "0.0", "w": "0.0", "mass": "0.1"},
    # P2: a collapsing cold polytrope (the selfgravity-star case at a tenth of its pressure), creation on, no sink
    "collapse": {"rho": 1.0e-2, "p": 1.0e-4, "time_max": 0.6, "x": "", "y": "", "z": "", "u": "", "v": "", "w": "",
                 "mass": ""},
    # P2: a static periodic gamma = 1.4 box at 1.2 rho_thr: above the threshold but Jeans-stable, creation on
    "stable": {"rho": 1.0, "p": 0.0185, "time_max": 100.0, "x": "", "y": "", "z": "", "u": "", "v": "", "w": "",
               "mass": ""},
    # P2: two approaching sinks 1.28 cells apart: bound and converging, they merge
    "merge": {"rho": 1.0e-6, "p": 1.0e-6, "time_max": 100.0, "x": "-0.04 0.04", "y": "0.0 0.0", "z": "0.0 0.0",
              "u": "0.1 -0.1", "v": "0.0 0.0", "w": "0.0 0.0", "mass": "1.0 1.0"},
}


def case_ini(a: argparse.Namespace) -> str:
    """Return the INI of a case."""
    s = BASE.read_text()
    s = re.sub(r"(?m)^output_basename\s*=.*$", "output_basename = sinks", s)
    s = re.sub(r"(?m)^restart_basename\s*=.*$", "restart_basename = sinks-restart", s)
    s = re.sub(r"(?m)^restart_save\s*=.*$", f"restart_save = {a.restart_save}", s)
    s = re.sub(r"(?m)^restart\s*=.*$", f"restart = {'.true.' if a.restart else '.false.'}", s)
    for key in ("ni", "nj", "nk"):
        s = re.sub(rf"(?m)^{key}\s*=.*$", f"{key} = {a.cells // 2}", s)
    s = re.sub(r"(?m)^physical_model\s*=.*$", f"physical_model = {a.model}", s)
    s = re.sub(r"(?m)^it_max\s*=.*$", f"it_max = {a.it_max}", s)
    # the gas falls into the sinks' wells: the positivity limiter keeps a near-vacuum background admissible
    s = re.sub(r"(?m)^positivity_limiter\s*=.*$", "positivity_limiter = cell", s)
    s = re.sub(r"(?m)^time_max\s*=.*$", f"time_max = {a.time_max!r}", s)
    if a.case in ("rest", "accrete", "stable"):  # a periodic box: no wall forces, a uniform background pulls nowhere
        s = re.sub(r"(?m)^type\s*=\s*wall-inviscid\s*$", "type = periodic", s)
        s = re.sub(r"(?m)^poisson_bc\s*=.*$", "poisson_bc = periodic", s)
    if a.case == "stable":
        s = re.sub(r"(?m)^gamma\s*=.*$", "gamma = 1.4", s)
    if a.case == "collapse":  # the polytrope stays, at a tenth of its equilibrium pressure
        s = re.sub(r"(?m)^star_pressure_factor\s*=.*$", "star_pressure_factor = 0.1", s)
    if a.case != "collapse":  # the uniform background in place of the polytrope: one region over the whole domain
        i = s.index("[initial_conditions]")
        s = s[:i] + (
            "[initial_conditions]\ntype = riemann-problem\nregions_number = 1\namr_iterations = 0\n\n"
            f"[initial_conditions_region_1]\nr = {a.rho!r}\nu = 0.0\nv = 0.0\nw = 0.0\np = {a.p!r}\n"
            "emin_x = -1.0\nemin_y = -1.0\nemin_z = -1.0\nemax_x = 1.0\nemax_y = 1.0\nemax_z = 1.0\n"
        )
    if a.model == "mhd-ideal":
        s += "bx = 0.0\nby = 0.0\nbz = 0.0\n"  # the region needs the field (the refusal comes later, at [sinks])
        s += "\n[mhd]\ndivergence_control = none\ndivb_tol = 0.0\ndivb_error = .false.\n"
        s += "rho_floor = 0.0\np_floor = 0.0\n"
    if a.no_gravity:  # the contract leg: sinks without self-gravity
        i = s.index("[gravity]")
        j = s.index("[initial_conditions]")
        s = s[:i] + s[j:]
    if a.case == "collapse":  # isothermal at T = 0.01 (below the centre's 0.016): strong cooling with that floor
        s += "\n[cooling]\ntemperatures = 0.01 1.0\nlambdas = 1000.0 1000.0\nparticle_mass = 1.0\n"
    if not a.no_sinks_section:
        s += "\n[sinks]\n" + "".join(f"{k} = {getattr(a, k)}\n" for k in ("x", "y", "z", "u", "v", "w", "mass")
                                      if getattr(a, k))
        s += f"creation = {'.true.' if a.case in ('collapse', 'stable') or a.creation else '.false.'}\n"
    return s


def main() -> int:
    """Write the requested case."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("out", type=Path)
    ap.add_argument("--case", choices=tuple(CASES), default="binary")
    ap.add_argument("--model", choices=("euler", "euler-de", "mhd-ideal"), default="euler")
    ap.add_argument("--cells", type=int, default=32)
    ap.add_argument("--rho", type=float, default=None)
    ap.add_argument("--p", type=float, default=None)
    ap.add_argument("--time-max", type=float, default=None)
    ap.add_argument("--it-max", type=int, default=-1)
    ap.add_argument("--restart-save", type=int, default=0)
    ap.add_argument("--restart", action="store_true")
    ap.add_argument("--no-sinks-section", action="store_true")
    ap.add_argument("--no-gravity", action="store_true", help="leave the [gravity] section out (the contract leg)")
    ap.add_argument("--creation", action="store_true", help="[sinks] creation = .true. (default for collapse, stable)")
    for key in ("x", "y", "z", "u", "v", "w", "mass"):
        ap.add_argument(f"--{key}", default=None)
    a = ap.parse_args()
    for key, value in CASES[a.case].items():
        if getattr(a, key) is None:
            setattr(a, key, value)
    a.out.write_text(case_ini(a))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
