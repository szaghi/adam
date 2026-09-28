#!/usr/bin/env python3
"""Write a magnetised vortex input (issue #41, M2-P5, MV-7; Balsara 2004, ApJS 151).

Why: the magnetised vortex is an exact steady equilibrium (the radial balance of pressure, centrifugal force, magnetic
tension and magnetic pressure) convected by the free stream, with a divergence-free field: after one crossing of the
doubly periodic box the exact solution is the initial state. Unlike the linear waves of MV-5 and the CPAW of MV-6 it is
a 2-D structure of every variable (velocity, field and pressure all vary in both directions), and its Gaussian tails
make the periodic box an approximation whose error, exp((1 - (L/radius)^2) / 2) at the box edge, must stay below the
measured error. Background: rho = 1, p = 1, v = (1, 1, 0), B = 0, kappa = mu = 1, radius 1, gamma = 5/3.

Divergence control GLM (the production formulation, issue #41 D-1): with `none` the discrete div(B) of the vortex
grows exponentially, at a rate proportional to 1/h (measured, M2-P5c: the run stops on a negative pressure at t = 13.4,
8.2, 4.3 for N = 64, 128, 256 on [-7, 7]^2, whatever the CFL), which GLM removes. c_h = 3 above max(|u_d| + c_f) ~ 2.5
(glm_ch_check = error), damping alpha 0.18 over the vortex radius.

Geometry: the vortex grid of V2 (verification/vortex, 4x4 blocks) stretched to [-L, L]^2, N cells per side, doubly
periodic. WENO-5 characteristic, SSP-54, one crossing t = 2 L (the free stream moves by (2L, 2L)).

Usage:
    make_mhd_vortex.py <base.ini> <out.ini> --cells N [--half-width L] [--cfl C] [--time-max T]
"""

from __future__ import annotations

import argparse
import configparser
from pathlib import Path


def main() -> None:
    """Write the input."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("base", type=Path, help="verification/vortex/vortex-n064.ini")
    parser.add_argument("out", type=Path)
    parser.add_argument("--cells", type=int, required=True, help="cells per side")
    parser.add_argument("--half-width", type=float, default=10.0, help="the box is [-L, L]^2")
    parser.add_argument("--cfl", default="0.1")
    parser.add_argument("--time-max", type=float, default=None, help="default one crossing, 2 L")
    args = parser.parse_args()
    half = args.half_width
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.optionxform = str
    ini.read(args.base)
    ini["physics"].update({"physical_model": "mhd-ideal", "cp": "2.5", "cv": "1.5"})  # gamma = 5/3
    ini["mhd"] = {"divergence_control": "glm", "divb_tol": "0.0", "divb_error": ".false.", "rho_floor": "0.0",
                  "p_floor": "0.0", "glm_ch": "3.0", "glm_alpha": "0.18", "glm_damping_length": "1.0",
                  "glm_ch_check": "error"}
    ini["field"]["nv"] = "9"
    ini["grid"].update({"ni": str(args.cells // 4), "nj": str(args.cells // 4), "emin_x": repr(-half),
                        "emin_y": repr(-half), "emax_x": repr(half), "emax_y": repr(half)})
    for f in ("x", "y"):
        ini[f"bc_{f}_min"]["type"] = "periodic"
        ini[f"bc_{f}_max"]["type"] = "periodic"
    for sec in [k for k in ini.sections() if k.startswith("initial_conditions_region_")]:
        ini.remove_section(sec)
    ini["initial_conditions"] = {"type": "mhd-vortex", "amr_iterations": "0", "x0": "0.0", "y0": "0.0",
                                 "radius": "1.0", "kappa": "1.0", "mu": "1.0"}
    ini["initial_conditions_region_1"] = {"r": "1.0", "u": "1.0", "v": "1.0", "w": "0.0", "p": "1.0", "bx": "0.0",
                                          "by": "0.0", "bz": "0.0"}
    ini["runge_kutta"]["scheme"] = "runge-kutta-ssp-54"
    time_max = 2.0 * half if args.time_max is None else args.time_max
    ini["time"].update({"it_max": "-1", "time_max": repr(time_max), "CFL": args.cfl})
    ini["IO"].update({"output_basename": "mhd-vortex", "it_save": "1000000"})
    with args.out.open("w") as out:
        ini.write(out)


if __name__ == "__main__":
    main()
