#!/usr/bin/env python3
"""Write a field loop advection input (issue #41, M2-P6, MV-9; Gardiner and Stone 2005, Mignone and Tzeferacos 2010).

Why: a weak magnetic loop (A_z = A0 (R - r), A0 = 1e-3, R = 0.3, plasma beta ~ 1e6) advected obliquely across a doubly
periodic box: the field is divergence-free pointwise but discontinuous at r = R, and a collocated scheme cannot keep
the discrete div(B) zero there. With the out-of-plane velocity w = 1 of Mignone and Tzeferacos (2010, section 4.4.1)
the equations give dB_z/dt = w div(B), so the generated <|B_z|> measures the divergence error that GLM has to keep at
truncation level; the magnetic energy decay measures the numerical dissipation of the loop.

Setup: rho = 1, p = 1, v = (2, 1, 1), gamma = 5/3, [-1, 1] x [-0.5, 0.5], t = 2 (two crossings along each axis). Grid:
the vortex grid of V2 (verification/vortex, 4x4 blocks) stretched to the box, N x N/2 cells (square cells). GLM
(c_h above the fastest signal, glm_ch_check = error), WENO-5 characteristic, SSP-54.

Usage:
    make_field_loop.py <base.ini> <out.ini> --cells N [--cfl C] [--time-max T] [--glm-ch C] [--vz W]
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
    parser.add_argument("--cells", type=int, required=True, help="cells along x (half of them along y)")
    parser.add_argument("--cfl", default="0.4")
    parser.add_argument("--time-max", default="2.0")
    parser.add_argument("--glm-ch", default="4.0")
    parser.add_argument("--vz", default="1.0", help="out-of-plane velocity w")
    parser.add_argument("--refine-box", type=float, nargs=4, default=None, metavar=("XMIN", "YMIN", "XMAX", "YMAX"),
                        help="refine the blocks whose centroid lies in the box by one 2:1 level (AMR variant)")
    args = parser.parse_args()
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.optionxform = str
    ini.read(args.base)
    ini["physics"].update({"physical_model": "mhd-ideal", "cp": "2.5", "cv": "1.5"})  # gamma = 5/3
    ini["mhd"] = {"divergence_control": "glm", "divb_tol": "0.0", "divb_error": ".false.", "rho_floor": "0.0",
                  "p_floor": "0.0", "glm_ch": args.glm_ch, "glm_alpha": "0.18", "glm_damping_length": "1.0",
                  "glm_ch_check": "error"}
    ini["field"]["nv"] = "9"
    ini["grid"].update({"ni": str(args.cells // 4), "nj": str(args.cells // 8), "emin_x": "-1.0", "emin_y": "-0.5",
                        "emax_x": "1.0", "emax_y": "0.5"})
    for f in ("x", "y"):
        ini[f"bc_{f}_min"]["type"] = "periodic"
        ini[f"bc_{f}_max"]["type"] = "periodic"
    for sec in [k for k in ini.sections() if k.startswith("initial_conditions_region_")]:
        ini.remove_section(sec)
    ini["initial_conditions"] = {"type": "field-loop", "amr_iterations": "0", "x0": "0.0", "y0": "0.0",
                                 "loop_radius": "0.3", "loop_amplitude": "1.0e-3"}
    ini["initial_conditions_region_1"] = {"r": "1.0", "u": "2.0", "v": "1.0", "w": args.vz, "p": "1.0", "bx": "0.0",
                                          "by": "0.0", "bz": "0.0"}
    ini["runge_kutta"]["scheme"] = "runge-kutta-ssp-54"
    ini["time"].update({"it_max": "-1", "time_max": args.time_max, "CFL": args.cfl})
    ini["IO"].update({"output_basename": "field-loop", "it_save": "1000000"})
    ini["diagnostics"]["conservation_history_save"] = "1"
    if args.refine_box is not None:
        refine_box(ini, args.refine_box)
    with args.out.open("w") as out:
        ini.write(out)


if __name__ == "__main__":
    main()
