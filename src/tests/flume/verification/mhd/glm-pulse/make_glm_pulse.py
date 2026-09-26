#!/usr/bin/env python3
"""Write a GLM pulse MHD input along x (issue #41, M2-P4, MV-3).

Why: in 1-D the (B_x, psi) pair of the mixed GLM system obeys the telegraph equations
dB_x/dt + dpsi/dx = 0, dpsi/dt + c_h^2 dB_x/dx = -(c_h^2/c_p^2) psi, decoupled from the fluid (the flux of B_x along x
is psi, the flux of psi is c_h^2 B_x): a Gaussian pulse in B_x on a fluid at rest has an exact solution (d'Alembert
without damping) computed by glm_pulse_oracle.py. The input is derived from the Sod grid of V1
(verification/sod/sod-x.ini: four blocks along x, the other directions null), refined to N cells, with x periodic or
bounded by reflecting walls (B_x odd, psi even at the wall, D-12).

Usage:
    make_glm_pulse.py <sod-x.ini> <out.ini> --cells N --divergence-control none|glm [--alpha A]
                      [--damping-length L|min-cell] [--bc periodic|wall] [--time T]
"""

from __future__ import annotations

import argparse
import configparser
from pathlib import Path

CH = 2.0  # cleaning speed, above the fastest wave (the fast speed of the background is sqrt(5/3) = 1.29)
PULSE = {"pulse_axis": "x", "pulse_center": "0.5", "pulse_width": "0.1", "pulse_amplitude": "0.1"}
BACKGROUND = {"r": "1.0", "u": "0.0", "v": "0.0", "w": "0.0", "p": "1.0", "bx": "0.0", "by": "0.0", "bz": "0.0"}


def main() -> None:
    """Write the input."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("sod", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--cells", type=int, required=True, help="cells along x (a multiple of 8)")
    parser.add_argument("--divergence-control", choices=("none", "glm"), required=True)
    parser.add_argument("--alpha", default="0.0", help="glm_alpha (0: no damping)")
    parser.add_argument("--damping-length", default="1.0", help="glm_damping_length: a positive length or min-cell")
    parser.add_argument("--bc", choices=("periodic", "wall"), default="periodic")
    parser.add_argument("--time", default="0.3", help="final time")
    args = parser.parse_args()
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.optionxform = str
    ini.read(args.sod)
    ini["physics"].update({"physical_model": "mhd-ideal", "cp": "2.5", "cv": "1.5"})  # gamma = 5/3
    ini["mhd"] = {"divergence_control": args.divergence_control, "divb_tol": "0.0", "divb_error": ".false.",
                  "rho_floor": "0.0", "p_floor": "0.0"}
    if args.divergence_control == "glm":
        # glm_ch_check = error: the fastest wave never outruns c_h here, so a trip is a defect
        ini["mhd"].update({"glm_ch": repr(CH), "glm_alpha": args.alpha, "glm_damping_length": args.damping_length,
                           "glm_ch_check": "error"})
    ini["field"]["nv"] = "9" if args.divergence_control == "glm" else "8"
    ini["grid"].update({"ni": str(args.cells // 4), "nj": "2", "nk": "2"})  # four blocks along x, thin transverse
    bc = "periodic" if args.bc == "periodic" else "wall-inviscid"
    ini["bc_x_min"]["type"] = bc
    ini["bc_x_max"]["type"] = bc
    ini["initial_conditions"] = {"type": "glm-pulse", "amr_iterations": "0", **PULSE}
    ini.remove_section("initial_conditions_region_2")
    ini["initial_conditions_region_1"] = dict(BACKGROUND)
    ini["time"].update({"it_max": "-1", "time_max": args.time, "CFL": "0.5"})
    ini["IO"].update({"output_basename": "glm-pulse", "it_save": "100000"})
    with args.out.open("w") as out:
        ini.write(out)


if __name__ == "__main__":
    main()
