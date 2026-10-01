#!/usr/bin/env python3
"""Oracle of the Orszag-Tang vortex and of the rotor (issue #41, M2-P6, MV-12, MV-13): 180 degrees symmetry.

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

The rotor (MV-13) has a uniform field, which the rotation maps to its opposite: its symmetry is the rotation composed
with the sign flip B -> -B (the equations are invariant under it), so B_x, B_y are even and B_z odd
(--field-parity even); its outflow boundaries do not conserve, so it passes no --cons-tol (the check is skipped).

Usage:
    orszag_tang_oracle.py <work> [--sym-tol T] [--cons-tol T] [--field-parity odd|even] [--step last|N] [--ngc N]

--step N measures the symmetry on the checkpoint of step N instead of the last one (the rotor under weno-riemann checks
it early, before the flow amplifies the round-off seed).
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


def pairs(work: Path, ngc: int, step: str = "last") -> tuple[np.ndarray, np.ndarray]:
    """Return a checkpoint (the last by default) as two arrays [cells, 8]: every cell and its image under the 180 degrees rotation.

    The image of the centre (x, y, z) is (x_min + x_max - x, y_min + y_max - y, z), the extremes over the cell centres:
    any grid symmetric under the rotation (uniform, or with a symmetric 2:1 refinement, M2-P7b), any number of cells
    along the null z axis. A cell without an image means an asymmetric grid, which is fatal."""
    cells = {tuple(round(c, 10) for c in key): value for key, value in profile(work, step, ngc).items()}
    xs = [k[0] for k in cells]
    ys = [k[1] for k in cells]
    sx, sy = min(xs) + max(xs), min(ys) + max(ys)
    keys = sorted(cells)
    images = [(round(sx - k[0], 10), round(sy - k[1], 10), k[2]) for k in keys]
    missing = [k for k, i in zip(keys, images, strict=True) if i not in cells]
    if missing:
        sys.exit(f"orszag_tang_oracle: {work}: {len(missing)} cells have no image under the rotation "
                 f"(e.g. {missing[0]})")
    return np.array([cells[k] for k in keys]), np.array([cells[i] for i in images])


def main() -> int:
    """Run the MV-12 checks, print a report, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work", type=Path)
    parser.add_argument("--sym-tol", type=float, default=None, help="bound of the relative symmetry defect")
    parser.add_argument("--cons-tol", type=float, default=None, help="bound of the relative conservation drift")
    parser.add_argument("--field-parity", choices=("odd", "even"), default="odd", help="in-plane field parity")
    parser.add_argument("--step", default="last", help="checkpoint of the symmetry check: last (default) or a step")
    parser.add_argument("--ngc", type=int, default=3)
    args = parser.parse_args()
    parity = dict(PARITY)
    if args.field_parity == "even":
        parity.update({"bx": 1.0, "by": 1.0, "bz": -1.0})
    ok = True
    q, rotated = pairs(args.work, args.ngc, args.step)
    worst = 0.0
    for v, name in enumerate(NAMES):
        scale = max(float(np.abs(q[:, NAMES.index(g)]).max()) for g in GROUP[name])
        defect = float(np.abs(q[:, v] - parity[name] * rotated[:, v]).max() / scale) if scale > 0.0 else 0.0
        worst = max(worst, defect)
        print(f"   symmetry {name:3s}: {defect:.3e} (max|{name}| {np.abs(q[:, v]).max():.3e}, scale {scale:.3e})")
    at = "" if args.step == "last" else f" at step {args.step}"
    line = f"180-degree symmetry{at}, worst relative defect {worst:.3e}"
    if args.sym_tol is not None:
        good = worst <= args.sym_tol
        ok &= good
        line += f"  {'PASS' if good else 'FAIL'} (max {args.sym_tol:.1e})"
    print(line)
    if args.cons_tol is None:
        return 0 if ok else 1
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
