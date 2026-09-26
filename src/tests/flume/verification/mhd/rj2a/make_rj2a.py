#!/usr/bin/env python3
"""Write an RJ2a MHD Riemann-problem input along one axis (issue #41, M2-P3, MV-4).

Why: RJ2a (Ryu & Jones 1995 Fig. 2a; exact 7-state solution in Dai & Woodward 1994, Tables Ia and Ib) is run along x,
y and z on the Sod grid of the M1 verification (verification/sod/sod-{x,y,z}.ini: four blocks along the axis, the
discontinuity at 0.5), refined to N cells. The states are given in the frame (normal, tangent 1, tangent 2) with the
tangents cyclic as in FLUME (and Athena++): axis y -> (y, z, x), axis z -> (z, x, y); B in code units, B_SI/sqrt(mu0),
i.e. the Gaussian values divided by sqrt(4 pi).

Usage:
    make_rj2a.py <sod-axis.ini> <out.ini> --axis x|y|z --cells N --divergence-control none|glm
"""

from __future__ import annotations

import argparse
import configparser
import math
from pathlib import Path

S4PI = math.sqrt(4.0 * math.pi)
# (rho, u_n, u_t1, u_t2, p, b_n, b_t1, b_t2), left and right of the discontinuity at 0.5
LEFT = (1.08, 1.2, 0.01, 0.5, 0.95, 2.0 / S4PI, 3.6 / S4PI, 2.0 / S4PI)
RIGHT = (1.0, 0.0, 0.0, 0.0, 1.0, 2.0 / S4PI, 4.0 / S4PI, 2.0 / S4PI)
FRAME = {"x": (0, 1, 2), "y": (1, 2, 0), "z": (2, 0, 1)}  # (normal, tangent 1, tangent 2) as x/y/z indexes


def main() -> None:
    """Write the input."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("sod", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--axis", choices=("x", "y", "z"), required=True)
    parser.add_argument("--cells", type=int, required=True, help="cells along the axis (a multiple of 8)")
    parser.add_argument("--divergence-control", choices=("none", "glm"), required=True)
    args = parser.parse_args()
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.optionxform = str
    ini.read(args.sod)
    ini["physics"].update({"physical_model": "mhd-ideal", "cp": "2.5", "cv": "1.5"})  # gamma = 5/3
    ini["mhd"] = {"divergence_control": args.divergence_control, "divb_tol": "0.0", "divb_error": ".false.",
                  "rho_floor": "0.0", "p_floor": "0.0"}
    if args.divergence_control == "glm":
        # c_h below the fastest wave (the left state alone has |u| + c_f = 2.89 and persists at x = 0 up to t = 0.2): the
        # fluid sets dt, so the GLM run takes the steps of the run without cleaning (the pair check is bitwise); the
        # c_h check warns, as it should
        ini["mhd"].update({"glm_ch": "2.0", "glm_alpha": "0.18", "glm_damping_length": "1.0",
                           "glm_ch_check": "warning"})
    ini["field"]["nv"] = "9" if args.divergence_control == "glm" else "8"
    ini["grid"]["n" + "ijk"["xyz".index(args.axis)]] = str(args.cells // 4)  # four blocks along the axis
    frame = FRAME[args.axis]
    for region, state in (("initial_conditions_region_1", LEFT), ("initial_conditions_region_2", RIGHT)):
        vel, mag = [0.0] * 3, [0.0] * 3
        for k in range(3):
            vel[frame[k]] = state[1 + k]
            mag[frame[k]] = state[5 + k]
        sec = ini[region]
        sec["r"], sec["p"] = repr(state[0]), repr(state[4])
        for key, val in zip(("u", "v", "w"), vel, strict=True):
            sec[key] = repr(val)
        for key, val in zip(("bx", "by", "bz"), mag, strict=True):
            sec[key] = repr(val)
    ini["time"].update({"it_max": "-1", "time_max": "0.2", "CFL": "0.5"})
    ini["IO"]["it_save"] = "100000"
    with args.out.open("w") as out:
        ini.write(out)


if __name__ == "__main__":
    main()
