#!/usr/bin/env python3
"""Write the inputs of the AS (AMR shock) runs: a shock crossing a static 2:1 face (#68 A3).

Why: every committed check of a shock across a 2:1 face compares two FLUME runs (multirealm legs, regression goldens);
none compares one with the exact solution, so a conservation defect, a seam reflection or an entropy glitch common to
every run would pass. Here the 1-D Sod (V1) and RJ2a (MV-4) problems run on the V1 grid (four blocks along x, the
discontinuity at 0.5) with one 2:1 face at x = 0.75: the right shock crosses it at t = 0.143 (Sod) and t = 0.110
(RJ2a, the fast shock) and is past it by t = 0.2, so the face then sits in the constant state behind the shock, where
every deviation from the exact solution is error that amr_shock_oracle.py measures.

Sides:
  c2f  the block x > 0.75 is refined: the shock passes from the coarse to the fine side;
  f2c  the blocks x < 0.75 are refined: the shock passes from the fine to the coarse side.
--uniform writes the same case without refinement (the control; --cells 2N gives the control at the fine spacing).
--no-reflux turns the Berger-Colella reflux off (a measurement of what it buys, not a supported configuration);
--it-save K saves a checkpoint every K steps (default only the last; the A4 sensor logs of verification/amr-log).
--seam-restriction R writes `[amr] seam_restriction = R` (#68 C3 measurement); the octree refines the null y and z
too, so their 4-cell blocks are widened to 8 (the 6-point footprint needs 6), the same for both values of R so that
a comparison runs on one grid (the 1-D solution is the same on every transverse copy).

Usage:
    make_amr_shock.py sod|rj2a <out.ini> --cells N [--side c2f|f2c] [--uniform] [--no-reflux] [--it-save K]
                      [--seam-restriction mean|lagrange6|weno6]
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

SEAM = 0.75
BOXES = {"c2f": [SEAM, -1.0e30, 1.0e30, 1.0e30], "f2c": [-1.0e30, -1.0e30, SEAM, 1.0e30]}


def read(path: Path) -> configparser.ConfigParser:
    """Read an input keeping the case of its keys."""
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.optionxform = str
    ini.read(path)
    return ini


def base_input(case: str, cells: int) -> configparser.ConfigParser:
    """Return the uniform input of a case, N cells along x on four blocks."""
    if cells % 8:
        sys.exit("make_amr_shock: --cells must be a multiple of 8 (four blocks, an even block for the 2:1 level)")
    sod = VERIF / "sod" / "sod-x.ini"
    if case == "sod":
        ini = read(sod)
        ini["grid"]["ni"] = str(cells // 4)
        return ini
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "rj2a.ini"
        cmd = [sys.executable, str(VERIF / "mhd" / "rj2a" / "make_rj2a.py"), str(sod), str(out), "--axis", "x",
               "--cells", str(cells), "--divergence-control", "none"]
        subprocess.run(cmd, check=True)
        return read(out)


def main() -> None:
    """Write the input."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("case", choices=("sod", "rj2a"))
    parser.add_argument("out", type=Path)
    parser.add_argument("--cells", type=int, required=True, help="base (coarse) cells along x")
    parser.add_argument("--side", choices=tuple(BOXES), default="c2f", help="which side of x = 0.75 is refined")
    parser.add_argument("--uniform", action="store_true", help="no refinement (the control)")
    parser.add_argument("--no-reflux", action="store_true", help="turn the reflux off")
    parser.add_argument("--it-save", type=int, default=1000000, help="checkpoint cadence in steps")
    parser.add_argument("--seam-restriction", choices=("mean", "lagrange6", "weno6"), default=None,
                        help="[amr] seam_restriction, null axes widened to 8 cells (default: the key absent)")
    args = parser.parse_args()
    ini = base_input(args.case, args.cells)
    if not args.uniform:
        refine_box(ini, BOXES[args.side], ratio=8)
    if args.no_reflux:
        ini["numerics"]["reflux"] = ".false."
    if args.seam_restriction is not None:
        ini["amr"]["seam_restriction"] = args.seam_restriction
        for axis in "jk":
            if int(ini["grid"][f"n{axis}"]) == 4:
                ini["grid"][f"n{axis}"] = "8"
    ini["IO"]["output_basename"] = args.case
    ini["IO"]["it_save"] = str(args.it_save)
    with args.out.open("w") as out:
        ini.write(out)


if __name__ == "__main__":
    main()
