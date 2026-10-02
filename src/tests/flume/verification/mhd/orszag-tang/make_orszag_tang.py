#!/usr/bin/env python3
"""Write an Orszag-Tang input (issue #41, M2-P6, MV-12; Stone et al. 2008, section 8.4).

Why: the Orszag-Tang vortex is the standard 2-D nonlinear MHD test: interacting shocks, a current sheet and low-beta
regions from a smooth, divergence-free, doubly periodic start. It has no exact solution, so what is checked is what the
equations guarantee: the 180 degrees rotation about (1/2, 1/2) maps the initial state to itself (scalars even, vectors
odd), so the discrete solution must keep that symmetry; the periodic box conserves every integral but psi's; with the
floors disabled a non-positive density or pressure stops the run, so a completed run has zero floored cells.

Grid: the vortex grid of V2 (verification/vortex, [0, 1]^2, 4x4 blocks), N cells per side, doubly periodic. GLM
(c_h above the fastest signal, glm_ch_check = error), WENO-5 characteristic, SSP-54, gamma = 5/3, t = 0.5.

Usage:
    make_orszag_tang.py <base.ini> <out.ini> --cells N [--cfl C] [--time-max T] [--glm-ch C] [--it-save N]
                        [--divergence-control glm|eglm]
                        [--refine-box XMIN YMIN XMAX YMAX]
"""

from __future__ import annotations

import argparse
import configparser
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from amr_box import refine_box  # noqa: E402


def main() -> None:
    """Write the input."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("base", type=Path, help="verification/vortex/vortex-n064.ini")
    parser.add_argument("out", type=Path)
    parser.add_argument("--cells", type=int, required=True, help="cells per side")
    parser.add_argument("--cfl", default="0.4")
    parser.add_argument("--time-max", default="0.5")
    parser.add_argument("--glm-ch", default="4.0")
    parser.add_argument("--divergence-control", choices=("glm", "eglm"), default="glm")
    parser.add_argument("--it-save", default="1000000", help="checkpoint period (the first and last are always saved)")
    parser.add_argument("--refine-box", type=float, nargs=4, default=None, metavar=("XMIN", "YMIN", "XMAX", "YMAX"),
                        help="refine the blocks whose centroid lies in the box by one 2:1 level (AMR variant)")
    args = parser.parse_args()
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.optionxform = str
    ini.read(args.base)
    ini["physics"].update({"physical_model": "mhd-ideal", "cp": "2.5", "cv": "1.5"})  # gamma = 5/3
    ini["mhd"] = {"divergence_control": args.divergence_control, "divb_tol": "0.0", "divb_error": ".false.",
                  "rho_floor": "0.0", "p_floor": "0.0", "glm_ch": args.glm_ch, "glm_alpha": "0.18",
                  "glm_damping_length": "1.0", "glm_ch_check": "error"}
    ini["field"]["nv"] = "9"
    ini["grid"].update({"ni": str(args.cells // 4), "nj": str(args.cells // 4), "emin_x": "0.0", "emin_y": "0.0",
                        "emax_x": "1.0", "emax_y": "1.0"})
    for f in ("x", "y"):
        ini[f"bc_{f}_min"]["type"] = "periodic"
        ini[f"bc_{f}_max"]["type"] = "periodic"
    for sec in [k for k in ini.sections() if k.startswith("initial_conditions_region_")]:
        ini.remove_section(sec)
    ini["initial_conditions"] = {"type": "orszag-tang", "amr_iterations": "0"}
    ini["runge_kutta"]["scheme"] = "runge-kutta-ssp-54"
    ini["time"].update({"it_max": "-1", "time_max": args.time_max, "CFL": args.cfl})
    ini["IO"].update({"output_basename": "orszag-tang", "it_save": args.it_save})
    ini["diagnostics"]["conservation_history_save"] = "1"
    if args.refine_box is not None:
        refine_box(ini, args.refine_box)
    with args.out.open("w") as out:
        ini.write(out)


if __name__ == "__main__":
    main()
