#!/usr/bin/env python3
"""FLUME MHD verification MV-4 oracle (issue #41, M2-P3): RJ2a against the exact solution.

Why this oracle exists: RJ2a (Ryu & Jones 1995 Fig. 2a) develops all seven MHD waves (fast and slow shocks, rotational
discontinuities, a contact) from one discontinuity; its exact solution is tabulated in Dai & Woodward (1994), Tables Ia
and Ib, as used by the Athena++ problem generator shock_tube.cpp (five significant digits). Comparing every conservative
variable, in the frame (normal, tangent 1, tangent 2), against it checks the whole MHD path: eigensystem, splitting,
fluxes, the rotation of the transverse components. Checks, on the last checkpoint of each run:

* 1-D consistency: every transverse copy of a column is bitwise identical, and blocks agree on shared columns;
* accuracy: L1 of every conservative variable against the exact solution; the sum over the 8 variables below --l1-max;
* direction invariance (several runs): the runs equal the first one in the rotated frame within --dir-tol (0:
  bitwise, the check.sh setting): the MHD library sums |u|^2, |B|^2 and u.B independently of the order of the terms
  (mhd_sum3) and projects in the frame order of the direction, so a cyclic rotation reproduces the x run exactly;
* --pair NONE GLM: the GLM run equals the run without cleaning BITWISE on the 8 shared variables and psi is exactly
  zero (in 1-D B_n is uniform and psi starts at zero, the (B_n, psi) block is inert and block-diagonal).

Usage:
    rj2a_oracle.py <work> [<work> ...] [--l1-max L] [--dir-tol T] [--ngc N]
    rj2a_oracle.py --pair <work-none> <work-glm> [--ngc N]
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import h5py
import numpy as np

S4PI = math.sqrt(4.0 * math.pi)
GAMMA = 5.0 / 3.0
T = 0.2
X0 = 0.5
NAMES = ("r", "mn", "mt1", "mt2", "E", "bn", "bt1", "bt2")
FRAME = ((0, 1, 2), (1, 2, 0), (2, 0, 1))  # (normal, tangent 1, tangent 2) of the axis x, y, z


def exact(xs: np.ndarray) -> np.ndarray:
    """Return the exact RJ2a conservative state (frame components) at the cell centres xs, time T."""
    xfp = 2.2638 * T
    xrp = (0.53432 + 1.0 / math.sqrt(math.pi * 1.309)) * T
    xsp = (0.53432 + 0.48144 / 1.309) * T
    xc = 0.57538 * T
    xsm = (0.60588 - 0.51594 / 1.4903) * T
    xrm = (0.60588 - 1.0 / math.sqrt(math.pi * 1.4903)) * T
    xfm = (1.2 - 2.3305 / 1.08) * T
    # (lower bound of r, rho, u_n, u_t1, u_t2, p, b_t1, b_t2) from right to left; b_n = 2/sqrt(4 pi) everywhere
    states = ((xfp, 1.0, 0.0, 0.0, 0.0, 1.0, 4.0, 2.0),
              (xrp, 1.3090, 0.53432, -0.094572, -0.047286, 1.5844, 5.3452, 2.6726),
              (xsp, 1.3090, 0.53432, -0.18411, 0.17554, 1.5844, 5.7083, 1.7689),
              (xc, 1.4735, 0.57538, 0.047601, 0.24734, 1.9317, 5.0074, 1.5517),
              (xsm, 1.6343, 0.57538, 0.047601, 0.24734, 1.9317, 5.0074, 1.5517),
              (xrm, 1.4903, 0.60588, 0.22157, 0.30125, 1.6558, 5.5713, 1.7264),
              (xfm, 1.4903, 0.60588, 0.11235, 0.55686, 1.6558, 5.0987, 2.8326),
              (-math.inf, 1.08, 1.2, 0.01, 0.5, 0.95, 3.6, 2.0))
    q = np.zeros((8, xs.size))
    bn = 2.0 / S4PI
    for n, x in enumerate(xs - X0):
        _, d, un, ut1, ut2, p, bt1, bt2 = next(s for s in states if x > s[0])
        bt1, bt2 = bt1 / S4PI, bt2 / S4PI
        e = p / (GAMMA - 1.0) + 0.5 * d * (un * un + ut1 * ut1 + ut2 * ut2) + 0.5 * (bn * bn + bt1 * bt1 + bt2 * bt2)
        q[:, n] = (d, d * un, d * ut1, d * ut2, e, bn, bt1, bt2)
    return q


def load(work: Path, ngc: int) -> tuple[int, np.ndarray, np.ndarray, np.ndarray | None]:
    """Return the active axis, the cell centres, the frame-rotated conservative profile and psi (or None)."""
    files = sorted(p for p in work.glob("*-proc*.h5") if "restart" not in p.name)
    if not files:
        sys.exit(f"rj2a_oracle: no checkpoint in {work}")
    last = max(int(p.name.split("-")[-2]) for p in files)
    cells: dict[float, np.ndarray] = {}
    axis, has_psi = -1, False
    for path in (p for p in files if int(p.name.split("-")[-2]) == last):
        with h5py.File(path, "r") as h5:
            for blk in sorted({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")}):
                origin, dxyz = h5[f"{blk}-origin"][()][::-1], h5[f"{blk}-dxdydz"][()][::-1]
                has_psi = f"{blk}-psi" in h5
                names = ("r", "ru", "rv", "rw", "rE", "bx", "by", "bz") + (("psi",) if has_psi else ())
                q = np.stack([h5[f"{blk}-{v}"][()].transpose(2, 1, 0) for v in names])
                if ngc:
                    q = q[:, ngc:-ngc, ngc:-ngc, ngc:-ngc]
                if axis < 0:
                    axis = int(np.argmax(q.shape[1:]))
                n = q.shape[1 + axis]
                qa = np.moveaxis(q, 1 + axis, 1).reshape(q.shape[0], n, -1)
                if not np.all(qa == qa[:, :, :1]):
                    sys.exit(f"rj2a_oracle: transverse copies differ in {path.name}:{blk}")
                centres = origin[axis] + (np.arange(n) + ngc + 0.5) * dxyz[axis]
                for c, col in zip(centres, qa[:, :, 0].T, strict=True):
                    if c in cells and not np.array_equal(cells[c], col):
                        sys.exit(f"rj2a_oracle: two blocks disagree on the column at {c}")
                    cells[c] = col
    xs = np.array(sorted(cells))
    q = np.array([cells[c] for c in xs]).T
    fn, f1, f2 = FRAME[axis]
    rot = np.vstack([q[0], q[1 + fn], q[1 + f1], q[1 + f2], q[4], q[5 + fn], q[5 + f1], q[5 + f2]])
    return axis, xs, rot, (q[8] if has_psi else None)


def main() -> int:
    """Run the MV-4 checks, print a report, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work", type=Path, nargs="*")
    parser.add_argument("--l1-max", type=float, default=None, help="bound on the L1 sum over the 8 variables")
    parser.add_argument("--dir-tol", type=float, default=0.0, help="direction comparison tolerance (0: bitwise)")
    parser.add_argument("--pair", type=Path, nargs=2, metavar=("NONE", "GLM"))
    parser.add_argument("--ngc", type=int, default=3)
    args = parser.parse_args()
    status = 0
    if args.pair is not None:
        _, _, q_none, _ = load(args.pair[0], args.ngc)
        _, _, q_glm, psi = load(args.pair[1], args.ngc)
        diff = float(np.max(np.abs(q_glm - q_none)))
        psi_max = float(np.max(np.abs(psi))) if psi is not None else math.inf
        ok = diff == 0.0 and psi_max == 0.0
        print(f"{args.pair[1].name} vs {args.pair[0].name}: max |difference| {diff:.3e} (8 shared variables), "
              f"max |psi| {psi_max:.3e}  {'PASS (bitwise, psi exactly zero)' if ok else 'FAIL'}")
        return 0 if ok else 1
    profiles = []
    for work in args.work:
        axis, xs, q, _ = load(work, args.ngc)
        dx = xs[1] - xs[0]
        l1 = np.sum(np.abs(q - exact(xs)), axis=1) * dx
        total = float(np.sum(l1))
        verdict = ""
        if args.l1_max is not None:
            ok = total <= args.l1_max
            status |= 0 if ok else 1
            verdict = "  PASS" if ok else f"  FAIL (bound {args.l1_max:.4e})"
        print(f"{work.name}: axis {'xyz'[axis]}, cells {xs.size}, L1 " +
              " ".join(f"{n} {v:.3e}" for n, v in zip(NAMES, l1, strict=True)) + f", sum {total:.6e}{verdict}")
        profiles.append((work.name, q))
    for name, q in profiles[1:]:
        diff = float(np.max(np.abs(q - profiles[0][1])))
        ok = diff <= args.dir_tol
        status |= 0 if ok else 1
        kind = "bitwise" if args.dir_tol == 0.0 else f"tol {args.dir_tol:.1e}"
        print(f"{name} vs {profiles[0][0]} (rotated frame): max |difference| {diff:.3e}  "
              f"{'PASS' if ok else 'FAIL'} ({kind})")
    return status


if __name__ == "__main__":
    sys.exit(main())
