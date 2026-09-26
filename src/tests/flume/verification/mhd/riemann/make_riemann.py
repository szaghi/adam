#!/usr/bin/env python3
"""Write a 1-D MHD Riemann-problem input along x (issue #41, M2-P3c): Brio-Wu, RJ4d, or a negative-pressure state.

Why: the M2-P3 exit asks Brio-Wu and RJ4d to run clean (no floored cell) and the positivity floors to be exercised.
The inputs are derived from the Sod input of the M1 verification (verification/sod/sod-x.ini: four blocks along x, the
discontinuity at 0.5), refined to N cells. States (x frame) from the Athena test suites: Brio-Wu (Athena++
inputs/mhd/athinput.bw: gamma = 2, t = 0.1), RJ4d (Ryu & Jones 1995, Athena C tst/1D-mhd/athinput.rj4d: gamma = 5/3,
t = 0.16). The `negative` case is the Brio-Wu left state against a right state of pressure -0.1: a non-positive state
from the first stage, which the floors must fix or, disabled, report.

Usage:
    make_riemann.py <sod-x.ini> <out.ini> --case bw|rj4d|negative --cells N --divergence-control none|glm
                    [--rho-floor F] [--p-floor F]
"""

from __future__ import annotations

import argparse
import configparser
from pathlib import Path

# (gamma, t_end, left, right), states (rho, u, v, w, p, bx, by, bz)
CASES = {
    "bw": (2.0, 0.1, (1.0, 0.0, 0.0, 0.0, 1.0, 0.75, 1.0, 0.0), (0.125, 0.0, 0.0, 0.0, 0.1, 0.75, -1.0, 0.0)),
    "rj4d": (5.0 / 3.0, 0.16, (1.0, 0.0, 0.0, 0.0, 1.0, 0.7, 0.0, 0.0), (0.3, 0.0, 0.0, 1.0, 0.2, 0.7, 1.0, 0.0)),
    "negative": (2.0, 0.02, (1.0, 0.0, 0.0, 0.0, 1.0, 0.75, 1.0, 0.0), (0.125, 0.0, 0.0, 0.0, -0.1, 0.75, -1.0, 0.0)),
}


def main() -> None:
    """Write the input."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("sod", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--case", choices=tuple(CASES), required=True)
    parser.add_argument("--cells", type=int, required=True, help="cells along x (a multiple of 8)")
    parser.add_argument("--divergence-control", choices=("none", "glm"), required=True)
    parser.add_argument("--rho-floor", type=float, default=0.0)
    parser.add_argument("--p-floor", type=float, default=0.0)
    args = parser.parse_args()
    gamma, t_end, left, right = CASES[args.case]
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.optionxform = str
    ini.read(args.sod)
    ini["physics"].update({"physical_model": "mhd-ideal", "cp": repr(gamma), "cv": "1.0"})
    ini["mhd"] = {"divergence_control": args.divergence_control, "divb_tol": "0.0", "divb_error": ".false.",
                  "rho_floor": repr(args.rho_floor), "p_floor": repr(args.p_floor)}
    if args.divergence_control == "glm":
        ini["mhd"].update({"glm_ch": "4.0", "glm_alpha": "0.18", "glm_damping_length": "1.0",
                           "glm_ch_check": "warning"})
    ini["field"]["nv"] = "9" if args.divergence_control == "glm" else "8"
    ini["grid"]["ni"] = str(args.cells // 4)
    ini["IO"]["save_auxiliary_fields"] = ".true."
    ini["IO"]["it_save"] = "100000"
    for region, state in (("initial_conditions_region_1", left), ("initial_conditions_region_2", right)):
        for key, val in zip(("r", "u", "v", "w", "p", "bx", "by", "bz"), state, strict=True):
            ini[region][key] = repr(val)
    ini["time"].update({"it_max": "-1", "time_max": repr(t_end), "CFL": "0.5"})
    with args.out.open("w") as out:
        ini.write(out)


if __name__ == "__main__":
    main()
