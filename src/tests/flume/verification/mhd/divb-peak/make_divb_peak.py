#!/usr/bin/env python3
"""Write the Dedner et al. (2002) div(B) peak MHD input (issue #41, M2-P4, MV-10).

Why: a peak in B_x on a uniform magnetised flow has a non-zero div(B) from t = 0. Without cleaning the error is only
advected with the flow; mixed GLM transports it away at c_h and damps it, so the integral of |div B| must decay
(Dedner et al. 2002, section 5.1). The input is derived from the periodic AMR conservation grid of V3
(verification/conservation/amr-periodic.ini) made 2-D (z null) on [-0.5, 1.5]^2, periodic, 64^2 cells (uniform level
2, 16^2 cells per block; nk = 4: the 2:1 seam ghost fill needs 4 cells per axis, null ones included); --amr refines
the blocks of [-0.5, 0.5]^2 once more (a 2:1 seam the peak crosses near t = 0.375). State of Dedner: rho = 1,
u = v = 1, w = 0, p = 6, B = (peak, 0, 1/sqrt(4 pi)), gamma = 5/3; the peak B_x += (1 - (r / R)^2)^2 / sqrt(4 pi),
R = 1/8, at the origin.

Usage:
    make_divb_peak.py <amr-periodic.ini> <out.ini> --divergence-control none|glm [--amr] [--time T] [--it-max N]
                      [--alpha A] [--divb-tol TOL] [--divb-error] [--save-aux] [--cells-per-block N]
"""

from __future__ import annotations

import argparse
import configparser
import math
from pathlib import Path

S4PI = math.sqrt(4.0 * math.pi)
CH = 5.0  # cleaning speed, above the fastest wave (|u_d| + c_f <= 1 + 3.2)


def main() -> None:
    """Write the input."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("base", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--divergence-control", choices=("none", "glm"), required=True)
    parser.add_argument("--amr", action="store_true", help="refine [-0.5, 0.5]^2 once more (2:1 seam)")
    parser.add_argument("--time", default="0.5", help="final time")
    parser.add_argument("--it-max", default="-1", help="iterations (-1: up to the final time)")
    parser.add_argument("--alpha", default="0.4", help="glm_alpha (damping length min-cell)")
    parser.add_argument("--divb-tol", default="0.0", help="[mhd].(divb_tol)")
    parser.add_argument("--divb-error", action="store_true", help="[mhd].(divb_error) = .true.")
    parser.add_argument("--save-aux", action="store_true", help="save the auxiliary and derived fields")
    parser.add_argument("--cells-per-block", default="16", help="ni = nj of a block (the AMR variant uses 8)")
    args = parser.parse_args()
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.optionxform = str
    ini.read(args.base)
    ini["physics"].update({"physical_model": "mhd-ideal", "cp": "2.5", "cv": "1.5"})  # gamma = 5/3
    ini["mhd"] = {"divergence_control": args.divergence_control, "divb_tol": args.divb_tol,
                  "divb_error": ".true." if args.divb_error else ".false.", "rho_floor": "0.0", "p_floor": "0.0"}
    if args.divergence_control == "glm":
        ini["mhd"].update({"glm_ch": repr(CH), "glm_alpha": args.alpha, "glm_damping_length": "min-cell",
                           "glm_ch_check": "error"})
    ini["field"]["nv"] = "9" if args.divergence_control == "glm" else "8"
    ini["grid"].update({"ni": args.cells_per_block, "nj": args.cells_per_block, "nk": "4", "emin_x": "-0.5", "emin_y": "-0.5", "emax_x": "1.5",
                        "emax_y": "1.5", "null_z": ".true."})
    ini["amr"].update({"iu_ref_levels": "2", "max_level": "3", "markers_number": "1" if args.amr else "0"})
    # the box spans the whole null z extent: a partial one would put 2:1 seams across z
    ini["amr_marker_1"].update({"box_xmin": "-0.5", "box_ymin": "-0.5", "box_zmin": "-1.0", "box_xmax": "0.5",
                                "box_ymax": "0.5", "box_zmax": "2.0", "target_level": "3"})
    ini["initial_conditions"] = {"type": "divb-peak", "amr_iterations": "1" if args.amr else "0", "peak_x0": "0.0",
                                 "peak_y0": "0.0", "peak_radius": "0.125", "peak_amplitude": repr(1.0 / S4PI)}
    ini["initial_conditions_region_1"] = {"r": "1.0", "u": "1.0", "v": "1.0", "w": "0.0", "p": "6.0", "bx": "0.0",
                                          "by": "0.0", "bz": repr(1.0 / S4PI)}
    ini["time"].update({"it_max": args.it_max, "time_max": args.time, "CFL": "0.5"})
    ini["IO"].update({"output_basename": "divb-peak", "it_save": "100000",
                      "save_auxiliary_fields": ".true." if args.save_aux else ".false."})
    ini["diagnostics"]["conservation_history_save"] = "1"
    with args.out.open("w") as out:
        ini.write(out)


if __name__ == "__main__":
    main()
