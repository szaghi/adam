#!/usr/bin/env python3
"""Write an MHD linear wave input (issue #41, M2-P5, MV-5; Stone et al. 2008, section 8.2).

Why: a small-amplitude eigenmode of the conservative MHD system travels unchanged at its wave speed, so after one
period the exact solution is the initial state: the error of every smooth-flow ingredient (eigen-decomposition, WENO
reconstruction, fluxes, periodic ghost fill, time integration) is measured on a problem with a known answer, family
by family. Background of Stone et al.: rho = 1, p = 3/5, B = (1, sqrt(2), 1/2) in the wave frame, gamma = 5/3, so
c_f = 2, c_a = 1, c_s = 1/2; the entropy wave rides a background flow u_n = 1. Amplitude 1e-6, wavelength 1.

Geometries: 1-D along x on the Sod grid of V1 (verification/sod/sod-x.ini, four blocks along x, x periodic), N cells
per wavelength; 2-D inclined (tan a = 2, a = 63.43 deg) on the vortex grid of V2 (verification/vortex, 4x4 blocks)
stretched to [0, sqrt(5)] x [0, sqrt(5)/2] with ni = 2 nj (square cells), one wavelength along each axis, doubly
periodic. WENO-5 characteristic, SSP-54, CFL 0.1 (the spatial error dominates), t = one period.

Usage:
    make_linear_wave.py <base.ini> <out.ini> --wave fast|alfven|slow|entropy --geometry 1d|2d --cells N [--amplitude A]
                        [--divergence-control none|glm|eglm]

With GLM or EGLM (issue #47, M3-P4c) the mode is the same (psi = 0, div B = 0), c_h = 1.1 (|u_n| + c_f) just above the
fastest wave (it sets the time step: 2.2, entropy 3.3), no damping.
"""

from __future__ import annotations

import argparse
import configparser
import math
from pathlib import Path

SPEED = {"fast": 2.0, "alfven": 1.0, "slow": 0.5, "entropy": 1.0}  # |u_n + c| of the right-going wave
TAN_A = 2.0


def main() -> None:
    """Write the input."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("base", type=Path, help="verification/sod/sod-x.ini (1d) or verification/vortex/*.ini (2d)")
    parser.add_argument("out", type=Path)
    parser.add_argument("--wave", choices=tuple(SPEED), required=True)
    parser.add_argument("--geometry", choices=("1d", "2d"), required=True)
    parser.add_argument("--cells", type=int, required=True, help="1d: cells along x; 2d: cells along x (2 x along y)")
    parser.add_argument("--amplitude", default="1.0e-6", help="wave amplitude")
    parser.add_argument("--divergence-control", choices=("none", "glm", "eglm"), default="none")
    args = parser.parse_args()
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.optionxform = str
    ini.read(args.base)
    ini["physics"].update({"physical_model": "mhd-ideal", "cp": "2.5", "cv": "1.5"})  # gamma = 5/3
    ini["mhd"] = {"divergence_control": args.divergence_control, "divb_tol": "0.0", "divb_error": ".false.",
                  "rho_floor": "0.0", "p_floor": "0.0"}
    ini["field"]["nv"] = "8" if args.divergence_control == "none" else "9"
    angle = math.degrees(math.atan(TAN_A)) if args.geometry == "2d" else 0.0
    c, s = math.cos(math.radians(angle)), math.sin(math.radians(angle))
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
    un = 1.0 if args.wave == "entropy" else 0.0
    if args.divergence_control != "none":
        ini["mhd"].update({"glm_ch": repr(1.1 * (un + 2.0)), "glm_alpha": "0.0", "glm_damping_length": "1.0",
                           "glm_ch_check": "error"})
    bn, bt1, bt2 = 1.0, math.sqrt(2.0), 0.5  # wave frame
    region = {"r": "1.0", "u": repr(un * c), "v": repr(un * s), "w": "0.0", "p": repr(0.6),
              "bx": repr(c * bn - s * bt1), "by": repr(s * bn + c * bt1), "bz": repr(bt2)}
    for sec in [k for k in ini.sections() if k.startswith("initial_conditions_region_")]:
        ini.remove_section(sec)
    ini["initial_conditions"] = {"type": "mhd-linear-wave", "amr_iterations": "0", "wave": args.wave,
                                 "wave_angle": repr(angle), "wave_amplitude": args.amplitude, "wavelength": "1.0"}
    ini["initial_conditions_region_1"] = region
    ini["runge_kutta"]["scheme"] = "runge-kutta-ssp-54"
    ini["time"].update({"it_max": "-1", "time_max": repr(1.0 / SPEED[args.wave]), "CFL": "0.1"})
    ini["IO"].update({"output_basename": "linear-wave", "it_save": "1000000"})
    with args.out.open("w") as out:
        ini.write(out)


if __name__ == "__main__":
    main()
