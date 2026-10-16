#!/usr/bin/env python3
"""Write the inputs of the AD (AMR divergence) cases: D_flux across static 2:1 seams (issue #68 A2).

Why: the div(B) history column of every MHD check is the centred divergence of the cell values, which is not the
divergence GLM and EGLM act on (#68 F4). With `[diagnostics] divb_flux` FLUME also writes D_flux, the divergence of the
interface normal field read from the scheme's own psi face fluxes, the coarse 2:1 faces taking the fine fluxes (the
reflux), so that it sums exactly (Gauss) over a composite grid. These cases measure it where the seams matter:
  cpaw  the circularly polarised Alfven wave of MV-6 (an exact, source-free nonlinear solution, div B = 0), inclined
        at tan a = 2 on [0, sqrt 5] x [0, sqrt 5 / 2], the two middle block columns refined (x in [0.6, 1.6]: two x seams
        that the wave crosses all the time), GLM, many periods: is the seam D_flux bounded or does it grow?
  loop  the field loop of MV-9 (#78), GLM or EGLM, the box [0, 1] x [-0.5, 0.5] refined, t = 1.
Both on a quadtree (null z, nk = 1); --uniform writes the same case without the box (the control); the div(B) history
is written every --history steps.

The seam studies of #68 vary the transfer operators and the reflux on these cases (issue #85 R2b, so that they run from
the script rather than from inputs edited by hand): --seam-ghost-fill F and --seam-restriction R write `[amr]
seam_ghost_fill` and `seam_restriction` (D, G3a), --reflux-policy P `[numerics] reflux_policy` (E2), --no-reflux turns
the reflux off (G3a), --glm-damping-length L `[mhd] glm_damping_length` (a length or `min-cell`, G2), --it-save K
saves a checkpoint every K steps (the figures). Absent, each key keeps the default above.

Usage:
    make_amr_divb.py cpaw <out.ini> --cells N [--uniform] [--periods T] [--history S] [options]
    make_amr_divb.py loop <out.ini> --cells N --cleaning glm|eglm [--uniform] [--time-max T] [--history S] [options]
"""

from __future__ import annotations

import argparse
import configparser
import subprocess
import sys
import tempfile
from pathlib import Path

VERIF = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(VERIF / "mhd"))
from amr_box import refine_box  # noqa: E402

BOXES = {"cpaw": [0.6, 0.0, 1.6, 2.0], "loop": [0.0, -0.5, 1.0, 0.5]}
GLM = {"divb_tol": "0.0", "divb_error": ".false.", "rho_floor": "0.0", "p_floor": "0.0", "glm_alpha": "0.18",
       "glm_damping_length": "1.0", "glm_ch_check": "error"}


def generate(script: Path, args: list[str]) -> configparser.ConfigParser:
    """Run a sibling generator into a temporary file and return its input."""
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "case.ini"
        subprocess.run([sys.executable, str(script), str(VERIF / "vortex" / "vortex-n064.ini"), str(out), *args],
                       check=True)
        ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
        ini.optionxform = str
        ini.read(out)
        return ini


def main() -> None:
    """Write the input."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("case", choices=tuple(BOXES))
    parser.add_argument("out", type=Path)
    parser.add_argument("--cells", type=int, required=True, help="base cells along x (cpaw >= 64: blocks >= 2 ngc)")
    parser.add_argument("--cleaning", choices=("glm", "eglm"), default="glm", help="loop only (cpaw: glm)")
    parser.add_argument("--uniform", action="store_true", help="no refined box (the control)")
    parser.add_argument("--periods", type=float, default=10.0, help="cpaw: final time, in wave periods")
    parser.add_argument("--time-max", default="1.0", help="loop: final time")
    parser.add_argument("--history", type=int, default=10, help="div(B) history cadence, steps")
    parser.add_argument("--seam-ghost-fill", default=None,
                        choices=("injection", "restriction-compatible", "tricubic", "quintic", "weno6", "quintic-mp"),
                        help="[amr] seam_ghost_fill (default: the key absent, tricubic)")
    parser.add_argument("--seam-restriction", choices=("mean", "lagrange6", "weno6"), default=None,
                        help="[amr] seam_restriction (default: the key absent, mean)")
    parser.add_argument("--reflux-policy", choices=("always", "troubled"), default=None,
                        help="[numerics] reflux_policy (default: the key absent, always)")
    parser.add_argument("--no-reflux", action="store_true", help="turn the reflux off")
    parser.add_argument("--glm-damping-length", default=None, help="[mhd] glm_damping_length: a length or min-cell")
    parser.add_argument("--it-save", type=int, default=1000000, help="checkpoint cadence in steps")
    args = parser.parse_args()
    if args.case == "cpaw":
        ini = generate(VERIF / "mhd" / "cpaw" / "make_cpaw.py",
                       ["--polarisation", "right", "--geometry", "2d", "--cells", str(args.cells)])
        ini["mhd"] = {"divergence_control": "glm", "glm_ch": "2.0", **GLM}  # c_h above the fast speed, ~1.08
        ini["field"]["nv"] = "9"
        ini["time"].update({"it_max": "-1", "time_max": repr(args.periods), "CFL": "0.4"})
    else:
        ini = generate(VERIF / "mhd" / "field-loop" / "make_field_loop.py",
                       ["--cells", str(args.cells), "--time-max", args.time_max])
        ini["mhd"]["divergence_control"] = args.cleaning
    if not args.uniform:
        refine_box(ini, BOXES[args.case], ratio=4)
    ini["diagnostics"].update({"conservation_history_save": str(args.history), "divb_flux": ".true."})
    if args.seam_ghost_fill is not None:
        ini["amr"]["seam_ghost_fill"] = args.seam_ghost_fill
    if args.seam_restriction is not None:
        ini["amr"]["seam_restriction"] = args.seam_restriction
    if args.reflux_policy is not None:
        ini["numerics"]["reflux_policy"] = args.reflux_policy
    if args.no_reflux:
        ini["numerics"]["reflux"] = ".false."
    if args.glm_damping_length is not None:
        ini["mhd"]["glm_damping_length"] = args.glm_damping_length
    ini["IO"]["it_save"] = str(args.it_save)
    with args.out.open("w") as out:
        ini.write(out)


if __name__ == "__main__":
    main()
