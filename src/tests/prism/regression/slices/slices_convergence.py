#!/usr/bin/env python3
"""PRISM [slices] convergence against an analytic solution (slices_restauration.md, F5).

Why: slices_oracle.py proves that every slice value is the trilinear (or inverse distance) combination of the right
eight cells; this check proves the end-to-end accuracy against physics. The case is an adimensional plane wave along +y
(AdimEM, kx = kz = 0, ky = 1): Dz = Bx = cos(2 pi / lambda (y - t)) in code units (light speed 1). After one step the
solver error is far below the interpolation error of a line slice off the cell centres, so the slice error against the
analytic wave is the trilinear error, O(h^2): halving h must divide it by ~4.

Usage:
    slices_convergence.py <work-coarse> <work-fine> [--ini NAME] [--p-min P]
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from slices_oracle import length_unit, read_ini  # noqa: E402

DZ, BX = 3 + 2, 3 + 3  # .mat columns: x, y, z, then Dx Dy Dz Bx ...


def slice_error(work: Path, ini_name: str) -> tuple[float, float, int]:
    """Return the max error of Dz and Bx of slice_1 against the analytic wave, the grid spacing and the iteration."""
    ini = read_ini(work / ini_name)
    unit = length_unit(ini)
    wavelength = float(ini["initial_conditions"]["lambda"]) / unit
    log = (work / "log.txt").read_text(errors="replace")
    times = re.findall(r"\[mpi-00000\]time:\s+(\S+)", log)
    time = float(times[-1])
    it = int(ini["time"]["it_max"])
    raw = next(work.glob(f"{ini['IO']['output_basename']}-slice_01-{it:09d}.mat")).read_bytes()
    nvar = int(np.frombuffer(raw[:16], dtype=np.int32)[3])
    data = np.frombuffer(raw[16:], dtype=np.float64).reshape(-1, nvar)
    exact = np.cos(2.0 * np.pi / wavelength * (data[:, 1] - time))
    err = max(float(np.max(np.abs(data[:, DZ] - exact))), float(np.max(np.abs(data[:, BX] - exact))))
    span = (float(ini["grid"]["emax_y"]) - float(ini["grid"]["emin_y"])) / unit
    h = span / (2 ** int(ini["amr"]["iu_ref_levels"]) * int(ini["grid"]["nj"]))
    return err, h, it


def main() -> int:
    """Run the check, print a report, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("coarse", type=Path)
    parser.add_argument("fine", type=Path)
    parser.add_argument("--ini", default="slices.ini")
    parser.add_argument("--p-min", type=float, default=1.8, help="minimum observed order")
    args = parser.parse_args()
    e_c, h_c, _ = slice_error(args.coarse, args.ini)
    e_f, h_f, _ = slice_error(args.fine, args.ini)
    p = np.log(e_c / e_f) / np.log(h_c / h_f)
    ok = p >= args.p_min
    print(f"plane wave slice vs analytic: h {h_c:.4g} error {e_c:.3e}, h {h_f:.4g} error {e_f:.3e}, observed order "
          f"{p:.2f}  {'PASS' if ok else 'FAIL'} (>= {args.p_min})")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
