#!/usr/bin/env python3
"""Write an MHD rotor input (issue #41, M2-P6, MV-13; Balsara and Spicer 1999, Toth 2000 first rotor).

Why: a dense disc spun up in a magnetised ambient launches torsional Alfven waves that brake it, compresses the field
and drives strong shocks: a stringent 2-D test of positivity and of the multi-dimensional field evolution. No exact
solution: checked are the 180 degrees symmetry of the rotation composed with B -> -B (the ambient field is uniform)
and positivity with the floors disabled. Setup of Toth: gamma = 1.4, p = 1, ambient rho = 1 at rest,
B = (5 / sqrt(4 pi), 0, 0), disc rho = 10 for r < 0.1 in rigid rotation at speed 2 at r = 0.1, linear taper to
r = 0.115, t = 0.15, outflow boundaries.

Grid: the vortex grid of V2 (verification/vortex, [0, 1]^2, 4x4 blocks), N cells per side (256 in MV-13, not the 200
of Toth: with a power-of-two N the cell centres are exact binary fractions, so the initial state is bitwise
symmetric), extrapolation boundaries. GLM (c_h above the fastest signal, glm_ch_check = error), WENO-5
characteristic, SSP-54.

Usage:
    make_rotor.py <base.ini> <out.ini> --cells N [--cfl C] [--time-max T] [--glm-ch C] [--it-save N]
"""

from __future__ import annotations

import argparse
import configparser
import math
from pathlib import Path


def main() -> None:
    """Write the input."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("base", type=Path, help="verification/vortex/vortex-n064.ini")
    parser.add_argument("out", type=Path)
    parser.add_argument("--cells", type=int, required=True, help="cells per side")
    parser.add_argument("--cfl", default="0.4")
    parser.add_argument("--time-max", default="0.15")
    parser.add_argument("--glm-ch", default="6.0")
    parser.add_argument("--it-save", default="1000000", help="checkpoint period (the first and last are always saved)")
    args = parser.parse_args()
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.optionxform = str
    ini.read(args.base)
    ini["physics"].update({"physical_model": "mhd-ideal", "cp": "3.5", "cv": "2.5"})  # gamma = 1.4
    ini["mhd"] = {"divergence_control": "glm", "divb_tol": "0.0", "divb_error": ".false.", "rho_floor": "0.0",
                  "p_floor": "0.0", "glm_ch": args.glm_ch, "glm_alpha": "0.18", "glm_damping_length": "1.0",
                  "glm_ch_check": "error"}
    ini["field"]["nv"] = "9"
    ini["grid"].update({"ni": str(args.cells // 4), "nj": str(args.cells // 4), "emin_x": "0.0", "emin_y": "0.0",
                        "emax_x": "1.0", "emax_y": "1.0"})
    for f in ("x", "y"):
        ini[f"bc_{f}_min"]["type"] = "extrapolation"
        ini[f"bc_{f}_max"]["type"] = "extrapolation"
    for sec in [k for k in ini.sections() if k.startswith("initial_conditions_region_")]:
        ini.remove_section(sec)
    ini["initial_conditions"] = {"type": "mhd-rotor", "amr_iterations": "0", "x0": "0.5", "y0": "0.5", "r0": "0.1",
                                 "r1": "0.115", "rho_in": "10.0", "v0": "2.0"}
    ini["initial_conditions_region_1"] = {"r": "1.0", "u": "0.0", "v": "0.0", "w": "0.0", "p": "1.0",
                                          "bx": repr(5.0 / math.sqrt(4.0 * math.pi)), "by": "0.0", "bz": "0.0"}
    ini["runge_kutta"]["scheme"] = "runge-kutta-ssp-54"
    ini["time"].update({"it_max": "-1", "time_max": args.time_max, "CFL": args.cfl})
    ini["IO"].update({"output_basename": "rotor", "it_save": args.it_save})
    ini["diagnostics"]["conservation_history_save"] = "1"
    with args.out.open("w") as out:
        ini.write(out)


if __name__ == "__main__":
    main()
