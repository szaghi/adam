#!/usr/bin/env python3
"""Oracle of the Orszag-Tang vortex (issue #41, M2-P6, MV-12; Stone et al. 2008, section 8.4).

Why: there is no exact solution, but the equations guarantee two properties the discrete solution must keep:
* the 180 degrees rotation about (1/2, 1/2) maps the initial state to itself, so the solution stays symmetric: at the
  last checkpoint q_v(x, y) = P_v q_v(1 - x, 1 - y), P_v = +1 for rho, E and the z components (rw, bz), -1 for the
  in-plane vectors (ru, rv, bx, by); measured as max over the cells of the difference, relative to the scale of the
  quantity: max|q_v| for rho and E, the largest component of the vector for the momentum and field components
  (--sym-tol). rw and bz are zero in exact arithmetic: on the CPU they stay exactly zero, on FNL they carry round-off
  noise, whose relative defect against its own maximum would be meaningless;
* the periodic box conserves the integrals of rho, rho u, E and B (conservation history; psi is damped, so excluded):
  the drift from the first row, relative to max(|first|, 1), stays below --cons-tol.
Positivity is the run's own: with the floors disabled a non-positive density or pressure stops it.

Usage:
    orszag_tang_oracle.py <work> [--sym-tol T] [--cons-tol T] [--ngc N]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "linear-wave"))
from linear_wave_oracle import NAMES, profile  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "divb-peak"))
from divb_peak_oracle import history  # noqa: E402

PARITY = {"r": 1.0, "ru": -1.0, "rv": -1.0, "rw": 1.0, "rE": 1.0, "bx": -1.0, "by": -1.0, "bz": 1.0}
CONSERVED = ("r", "ru", "rv", "rw", "rE", "bx", "by", "bz")
GROUP = {"r": ("r",), "rE": ("rE",), "ru": ("ru", "rv", "rw"), "rv": ("ru", "rv", "rw"), "rw": ("ru", "rv", "rw"),
         "bx": ("bx", "by", "bz"), "by": ("bx", "by", "bz"), "bz": ("bx", "by", "bz")}


def grid(work: Path, ngc: int) -> np.ndarray:
    """Return the last checkpoint as an array [8, nx, ny] (z null), cells ordered by their centres."""
    cells = profile(work, "last", ngc)
    xs = sorted({k[0] for k in cells})
    ys = sorted({k[1] for k in cells})
    q = np.full((len(NAMES), len(xs), len(ys)), np.nan)
    ix = {x: n for n, x in enumerate(xs)}
    iy = {y: n for n, y in enumerate(ys)}
    for key, value in cells.items():
        q[:, ix[key[0]], iy[key[1]]] = value
    if np.isnan(q).any():
        sys.exit(f"orszag_tang_oracle: {work} is not a full uniform 2-D grid")
    return q


def main() -> int:
    """Run the MV-12 checks, print a report, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work", type=Path)
    parser.add_argument("--sym-tol", type=float, default=None, help="bound of the relative symmetry defect")
    parser.add_argument("--cons-tol", type=float, default=None, help="bound of the relative conservation drift")
    parser.add_argument("--ngc", type=int, default=3)
    args = parser.parse_args()
    ok = True
    q = grid(args.work, args.ngc)
    rotated = q[:, ::-1, ::-1]
    worst = 0.0
    for v, name in enumerate(NAMES):
        scale = max(float(np.abs(q[NAMES.index(g)]).max()) for g in GROUP[name])
        defect = float(np.abs(q[v] - PARITY[name] * rotated[v]).max() / scale) if scale > 0.0 else 0.0
        worst = max(worst, defect)
        print(f"   symmetry {name:3s}: {defect:.3e} (max|{name}| {np.abs(q[v]).max():.3e}, scale {scale:.3e})")
    line = f"180-degree symmetry, worst relative defect {worst:.3e}"
    if args.sym_tol is not None:
        good = worst <= args.sym_tol
        ok &= good
        line += f"  {'PASS' if good else 'FAIL'} (max {args.sym_tol:.1e})"
    print(line)
    hist = history(args.work, "conservation")
    drift = 0.0
    for name in CONSERVED:
        col = hist[f"int_{name}"]
        d = float(np.abs(col - col[0]).max() / max(abs(col[0]), 1.0))
        drift = max(drift, d)
        print(f"   conservation int_{name:3s}: first {col[0]:+.15e}, drift {d:.3e}")
    line = f"conservation, worst relative drift {drift:.3e} over {len(hist['it'])} rows (t = {hist['time'][-1]:.6g})"
    if args.cons_tol is not None:
        good = drift <= args.cons_tol
        ok &= good
        line += f"  {'PASS' if good else 'FAIL'} (max {args.cons_tol:.1e})"
    print(line)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
