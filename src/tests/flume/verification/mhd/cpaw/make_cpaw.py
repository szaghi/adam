#!/usr/bin/env python3
"""Write a circularly polarised Alfven wave input (issue #41, M2-P5, MV-6; Toth 2000, J. Comput. Phys. 161).

Why: the circularly polarised Alfven wave is an exact nonlinear solution of ideal MHD at any amplitude (|B_perp| is
uniform, so is the total pressure), so after one period the exact solution is the initial state. At the finite
amplitude of Toth (B_perp = 0.1, ten per cent of B_par) the WENO nonlinear weights, which the 1e-7 linear waves of MV-5
leave at their linear values, are exercised on a problem with a known answer. Background of Toth: rho = 1, p = 0.1,
B_par = 1, v_par = 0, gamma = 5/3, so the Alfven speed and the period (wavelength 1) are 1.

Geometry: 2-D inclined (tan a = 2, a = 63.43 deg, Mignone and Tzeferacos 2010) on the vortex grid of V2
(verification/vortex, 4x4 blocks) stretched to [0, sqrt(5)] x [0, sqrt(5)/2] with ni = 2 nj (square cells), one
wavelength along each axis, doubly periodic; or 1-D along x on the Sod grid of V1. WENO-5 characteristic, SSP-54,
CFL 0.1 (the spatial error dominates), t = one period.

Usage:
    make_cpaw.py <base.ini> <out.ini> --polarisation right|left --geometry 1d|2d --cells N [--amplitude A]
"""

from __future__ import annotations

import argparse
import configparser
import math
from pathlib import Path

TAN_A = 2.0


def main() -> None:
    """Write the input."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("base", type=Path, help="verification/sod/sod-x.ini (1d) or verification/vortex/*.ini (2d)")
    parser.add_argument("out", type=Path)
    parser.add_argument("--polarisation", choices=("right", "left"), required=True)
    parser.add_argument("--geometry", choices=("1d", "2d"), required=True)
    parser.add_argument("--cells", type=int, required=True, help="1d: cells along x; 2d: cells along x (2 x along y)")
    parser.add_argument("--amplitude", default="0.1", help="B_perp")
    args = parser.parse_args()
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.optionxform = str
    ini.read(args.base)
    ini["physics"].update({"physical_model": "mhd-ideal", "cp": "2.5", "cv": "1.5"})  # gamma = 5/3
    ini["mhd"] = {"divergence_control": "none", "divb_tol": "0.0", "divb_error": ".false.", "rho_floor": "0.0",
                  "p_floor": "0.0"}
    ini["field"]["nv"] = "8"
    angle = math.degrees(math.atan(TAN_A)) if args.geometry == "2d" else 0.0
    if args.geometry == "1d":
        ini["grid"]["ni"] = str(args.cells // 4)  # four blocks along x
        faces = ("x",)
    else:
        ini["grid"].update({"ni": str(args.cells // 4), "nj": str(args.cells // 8), "emin_x": "0.0", "emin_y": "0.0",
                            "emax_x": repr(math.sqrt(5.0)), "emax_y": repr(math.sqrt(5.0) / 2.0)})
        faces = ("x", "y")
    for f in faces:
        ini[f"bc_{f}_min"]["type"] = "periodic"
        ini[f"bc_{f}_max"]["type"] = "periodic"
    region = {"r": "1.0", "u": "0.0", "v": "0.0", "w": "0.0", "p": "0.1", "bx": "0.0", "by": "0.0", "bz": "0.0"}
    for sec in [k for k in ini.sections() if k.startswith("initial_conditions_region_")]:
        ini.remove_section(sec)
    ini["initial_conditions"] = {"type": "mhd-cpaw", "amr_iterations": "0", "polarisation": args.polarisation,
                                 "b_par": "1.0", "wave_angle": repr(angle), "wave_amplitude": args.amplitude,
                                 "wavelength": "1.0"}
    ini["initial_conditions_region_1"] = region
    ini["runge_kutta"]["scheme"] = "runge-kutta-ssp-54"
    ini["time"].update({"it_max": "-1", "time_max": "1.0", "CFL": "0.1"})
    ini["IO"].update({"output_basename": "cpaw", "it_save": "1000000"})
    with args.out.open("w") as out:
        ini.write(out)


if __name__ == "__main__":
    main()
