#!/usr/bin/env python3
"""Write a rotated shock tube input (issue #41, M2-P6, MV-8; Toth 2000, J. Comput. Phys. 161, section 6.3.2).

Why: a one-dimensional MHD Riemann problem rotated in the x-y plane is solved by the dimension-by-dimension scheme
with truly two-dimensional data, while its exact solution keeps the normal field uniform, B_par = 5 / sqrt(4 pi): the
deviation of B_par measures the divergence error of the base scheme plus GLM, and the whole solution must match the
same problem solved in 1-D.

Problem: Ryu and Jones (1995) 1a in the frame of the normal (rho, v_par, v_perp, p, B_par, B_perp) = (1, 10, 0, 20,
5/sqrt(4 pi), 5/sqrt(4 pi)) on the left, (1, -10, 0, 1, 5/sqrt(4 pi), 5/sqrt(4 pi)) on the right, gamma = 5/3. FLUME has
no shifted-periodic boundary, so the problem is periodic along its normal: IC rotated-riemann, s = n_x x + n_y y,
period 1, the shock tube jump at s = 1/4 + 1/(4N) (off every cell centre). On a periodic line the velocity jump must
be undone: a second jump would be the mirror problem, an expansion close to vacuum (density ~2e-4, Alfven speed ~100:
max(|u| + c_f) grows with the resolution of the fan and no c_h holds), so the return to the left state is a linear ramp
of the primitive variables over s in (1/2, 9/10] (a resolved expansion, density ~1/3 at the final time). B_par is
uniform in both states, so the initial field is divergence-free and B_par = 5/sqrt(4 pi) stays exact; the 1-D reference
solves the same problem, ramp and wave interactions included. The final time t = 0.04 / |n| keeps the fan of the jump
(fast speeds ~5, Toth: fast shocks at 0.1 and 0.9 at t = 0.08) within 0.2 / |n| of it, as in Toth (t = 0.08 cos a on
his strip). The fastest signal, ~20-24 in the ramp, sets c_h.

--angle 45 (normal (1, 1)) or 63 (normal (1, 2), Toth's 63.4 degrees). --reference writes the 1-D problem along x
with normal (|n|, 0) on [0, 1/|n|], so that s is the same coordinate, 4 N cells per period (4 times finer along the
normal than the 2-D grid, whose cell centres are 1/N apart in s).

Grid: 2-D, the vortex grid of V2 (verification/vortex, 4x4 blocks), [0, 1]^2, N cells per side, periodic; 1-D, the Sod
grid of V1 (verification/sod, 4 blocks along x), periodic. GLM (c_h above the fastest signal, glm_ch_check = error),
WENO-5 characteristic, SSP-54, floors off.

Usage:
    make_rotated_shock_tube.py <base.ini> <out.ini> --angle {45,63} --cells N [--reference] [--cfl C] [--glm-ch C]
"""

from __future__ import annotations

import argparse
import configparser
import math
from pathlib import Path

NORMAL = {"45": (1, 1), "63": (1, 2)}
B0 = 5.0 / math.sqrt(4.0 * math.pi)
LEFT = {"r": 1.0, "u": 10.0, "v": 0.0, "w": 0.0, "p": 20.0, "bx": B0, "by": B0, "bz": 0.0}
RIGHT = {"r": 1.0, "u": -10.0, "v": 0.0, "w": 0.0, "p": 1.0, "bx": B0, "by": B0, "bz": 0.0}


def main() -> None:
    """Write the input."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("base", type=Path, help="verification/vortex/vortex-n064.ini (2-D), sod/sod-x.ini (1-D)")
    parser.add_argument("out", type=Path)
    parser.add_argument("--angle", choices=sorted(NORMAL), required=True)
    parser.add_argument("--cells", type=int, required=True, help="2-D cells per side N (a multiple of 4)")
    parser.add_argument("--reference", action="store_true", help="the 1-D reference problem, 4 N cells per period")
    parser.add_argument("--cfl", default="0.4")
    parser.add_argument("--glm-ch", default="24.0")
    args = parser.parse_args()
    nx, ny = NORMAL[args.angle]
    norm = math.hypot(nx, ny)
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.optionxform = str
    ini.read(args.base)
    ini["physics"].update({"physical_model": "mhd-ideal", "cp": "2.5", "cv": "1.5"})  # gamma = 5/3
    ini["mhd"] = {"divergence_control": "glm", "divb_tol": "0.0", "divb_error": ".false.", "rho_floor": "0.0",
                  "p_floor": "0.0", "glm_ch": args.glm_ch, "glm_alpha": "0.18", "glm_damping_length": "1.0",
                  "glm_ch_check": "error"}
    ini["field"]["nv"] = "9"
    if args.reference:
        # 4 blocks along x: 4 N cells; one cell on the null axes (their extent is irrelevant, the cost is not)
        ini["grid"].update({"ni": str(args.cells), "nj": "1", "nk": "1", "emin_x": "0.0", "emax_x": repr(1.0 / norm)})
        normal = (norm, 0.0)
        axes = ("x",)
    else:
        ini["grid"].update({"ni": str(args.cells // 4), "nj": str(args.cells // 4), "emin_x": "0.0", "emin_y": "0.0",
                            "emax_x": "1.0", "emax_y": "1.0"})
        normal = (float(nx), float(ny))
        axes = ("x", "y")
    for f in axes:
        ini[f"bc_{f}_min"]["type"] = "periodic"
        ini[f"bc_{f}_max"]["type"] = "periodic"
    for sec in [k for k in ini.sections() if k.startswith("initial_conditions_region_")]:
        ini.remove_section(sec)
    shift = 0.25 / args.cells
    ini["initial_conditions"] = {"type": "rotated-riemann", "amr_iterations": "0", "normal_x": repr(normal[0]),
                                 "normal_y": repr(normal[1]), "interface_1": repr(0.25 + shift),
                                 "interface_2": repr(0.5 + shift), "interface_2_width": "0.4", "period": "1.0"}
    for n, state in ((1, LEFT), (2, RIGHT)):
        ini[f"initial_conditions_region_{n}"] = {k: repr(v) for k, v in state.items()}
    ini["runge_kutta"]["scheme"] = "runge-kutta-ssp-54"
    ini["time"].update({"it_max": "-1", "time_max": repr(0.04 / norm), "CFL": args.cfl})
    ini["IO"].update({"output_basename": "rst", "it_save": "1000000"})
    if "diagnostics" in ini:
        ini["diagnostics"]["conservation_history_save"] = "1"
    with args.out.open("w") as out:
        ini.write(out)


if __name__ == "__main__":
    main()
