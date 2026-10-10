#!/usr/bin/env python3
"""Write the inputs of the AO (AMR order) ladders: a smooth vortex centred on the corner of a refined quadrant (#68 A1).

Why: every smooth order oracle of FLUME (V2, MV-5..7) runs on a uniform grid; no committed test measures the accuracy
across a static 2:1 seam, where the mean restriction of the coarse ghosts and the Berger-Colella reflux of point-value
fluxes are each expected to cap the composite error at second order (#68 F1, Appendix A). Here the vortex core sits on
the corner of the refined quadrant, so it straddles an x seam, a y seam and their corner at once, and it moves into the
fine quadrant during the run (free stream (1, 1)): the error of the seam band and of the interior are measured
separately by amr_order_oracle.py.

Cases:
  euler  the V2 isentropic vortex (verification/vortex, [0, 1]^2, radius 0.07, centre (0.5, 0.5)), t = 0.1 as V2,
         CFL 0.4 (V2 runs 0.1: the uniform N = 64 errors at 0.1 and 0.4 agree to 5 digits, L1 8.2184e-5 both, so
         the time error is negligible and the ladder costs 4x less); refined quadrant [0.5, 1]^2.
  mhd    the MV-7 magnetised vortex (verification/mhd/vortex, [-7, 7]^2, radius 1, centre (0, 0)), t = 0.5, GLM;
         refined quadrant [0, 7]^2.
Both: 4x4 base blocks, the quadrant refines four of them by one 2:1 level (mhd/amr_box.py, init-time primitive box,
no regridding); --ratio 4 is a quadtree (null z, nk = 1), 8 an octree (nk = 4); --uniform writes the same case
without the box (the uniform control at the coarse spacing).

Usage:
    make_amr_order.py euler <out.ini> --cells N [--ratio 4|8] [--uniform]
    make_amr_order.py mhd <out.ini> --cells N [--ratio 4|8] [--uniform] [--time-max T]
"""

from __future__ import annotations

import argparse
import configparser
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
VERIF = HERE.parent
sys.path.insert(0, str(VERIF / "mhd"))
from amr_box import refine_box  # noqa: E402

BOXES = {"euler": [0.5, 0.5, 1.0, 1.0], "mhd": [0.0, 0.0, 7.0, 7.0]}


def read(path: Path) -> configparser.ConfigParser:
    """Read an input keeping the case of its keys."""
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.optionxform = str
    ini.read(path)
    return ini


def base_input(case: str, cells: int, time_max: float | None) -> configparser.ConfigParser:
    """Return the uniform input of a case at N cells per side."""
    if case == "euler":
        base = VERIF / "vortex" / f"vortex-n{cells:03d}.ini"
        if not base.exists():
            sys.exit(f"make_amr_order: no V2 input {base.name} (N must be 64, 128 or 256)")
        ini = read(base)
        ini["IO"]["output_basename"] = "vortex"
        ini["time"]["CFL"] = "0.4"
        return ini
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "mhd-vortex.ini"
        cmd = [sys.executable, str(VERIF / "mhd" / "vortex" / "make_mhd_vortex.py"),
               str(VERIF / "vortex" / "vortex-n064.ini"), str(out), "--cells", str(cells), "--half-width", "7.0",
               "--cfl", "0.4", "--time-max", repr(0.5 if time_max is None else time_max)]
        subprocess.run(cmd, check=True)
        return read(out)


def main() -> None:
    """Write the input."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("case", choices=("euler", "mhd"))
    parser.add_argument("out", type=Path)
    parser.add_argument("--cells", type=int, required=True, help="base (coarse) cells per side")
    parser.add_argument("--ratio", type=int, default=4, choices=(4, 8), help="tree ratio: 4 quadtree, 8 octree")
    parser.add_argument("--uniform", action="store_true", help="no refined quadrant (the uniform control)")
    parser.add_argument("--time-max", type=float, default=None, help="mhd only, default 0.5")
    args = parser.parse_args()
    ini = base_input(args.case, args.cells, args.time_max)
    if args.uniform:
        ini["amr"]["ratio"] = str(args.ratio)
        if args.ratio == 8:
            ini["grid"]["nk"] = str(max(4, int(ini["grid"]["nk"])))
    else:
        refine_box(ini, BOXES[args.case], ratio=args.ratio)
    ini["IO"]["it_save"] = "1000000"
    with args.out.open("w") as out:
        ini.write(out)


if __name__ == "__main__":
    main()
