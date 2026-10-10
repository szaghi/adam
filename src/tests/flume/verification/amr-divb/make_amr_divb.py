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

Usage:
    make_amr_divb.py cpaw <out.ini> --cells N [--uniform] [--periods T] [--history S]
    make_amr_divb.py loop <out.ini> --cells N --cleaning glm|eglm [--uniform] [--time-max T] [--history S]
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
    ini["IO"]["it_save"] = "1000000"
    with args.out.open("w") as out:
        ini.write(out)


if __name__ == "__main__":
    main()
