#!/usr/bin/env python3
"""Oracle of the MHD linear waves (issue #41, M2-P5, MV-5; Stone et al. 2008, section 8.2).

Why: after one period the exact solution of a linear eigenmode is the initial state, so the error of a run is the
difference between its last and its first checkpoint: delta_v = mean |q_v(T) - q_v(0)| over the cells, and the error
norm of Stone et al., eps = sqrt(sum_v delta_v^2), over the 8 conservative variables. A wrong eigenvector splits the
mode into the other families and leaves an O(amplitude) error that does not converge; a wrong ingredient of the smooth
path lowers the observed order p = log2(eps_N / eps_2N), asserted on the finest pair (--order-min).

Usage:
    linear_wave_oracle.py <work> [<work> ...] [--order-min P] [--eps-max E [E ...]] [--ngc N]
"""

from __future__ import annotations

import argparse
import configparser
import math
import sys
from pathlib import Path

import h5py
import numpy as np

NAMES = ("r", "ru", "rv", "rw", "rE", "bx", "by", "bz")


def amplitude(work: Path) -> float:
    """Return the wave amplitude of the input of a run (its only .ini file)."""
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.read(next(work.glob("*.ini")))
    return float(ini["initial_conditions"]["wave_amplitude"])


def profile(work: Path, step: str, ngc: int) -> dict[tuple[float, ...], np.ndarray]:
    """Return the conservative variables of every interior cell of the first or last checkpoint (or of the checkpoint
    of a given step number), keyed by the cell centre."""
    files = sorted(p for p in work.glob("*-proc*.h5") if "restart" not in p.name)
    if not files:
        sys.exit(f"linear_wave_oracle: no checkpoint in {work}")
    steps = sorted({int(p.name.split("-")[-2]) for p in files})
    if len(steps) < 2:
        sys.exit(f"linear_wave_oracle: {work} needs the first and the last checkpoint")
    pick = steps[0] if step == "first" else steps[-1] if step == "last" else int(step)
    if pick not in steps:
        sys.exit(f"linear_wave_oracle: {work} has no checkpoint of step {pick}")
    cells: dict[tuple[float, ...], np.ndarray] = {}
    for path in (p for p in files if int(p.name.split("-")[-2]) == pick):
        with h5py.File(path, "r") as h5:
            for blk in sorted({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")}):
                origin, dxyz = h5[f"{blk}-origin"][()][::-1], h5[f"{blk}-dxdydz"][()][::-1]
                q = np.stack([h5[f"{blk}-{v}"][()].transpose(2, 1, 0) for v in NAMES])
                if ngc:
                    q = q[:, ngc:-ngc, ngc:-ngc, ngc:-ngc]
                idx = np.indices(q.shape[1:]).reshape(3, -1).T
                for i, j, k in idx:
                    key = tuple(round(float(origin[a] + (n + ngc + 0.5) * dxyz[a]), 12) for a, n in enumerate((i, j, k)))
                    cells[key] = q[:, i, j, k]
    return cells


def volumes(work: Path, step: str, ngc: int) -> dict[tuple[float, ...], float]:
    """Return the volume of every interior cell of the first or last checkpoint, keyed as profile() keys it (for the
    volume-weighted means of AMR runs, whose cells differ in size)."""
    files = sorted(p for p in work.glob("*-proc*.h5") if "restart" not in p.name)
    steps = sorted({int(p.name.split("-")[-2]) for p in files})
    pick = steps[0] if step == "first" else steps[-1]
    out: dict[tuple[float, ...], float] = {}
    for path in (p for p in files if int(p.name.split("-")[-2]) == pick):
        with h5py.File(path, "r") as h5:
            for blk in sorted({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")}):
                origin, dxyz = h5[f"{blk}-origin"][()][::-1], h5[f"{blk}-dxdydz"][()][::-1]
                shape = h5[f"{blk}-r"][()].shape[::-1]
                vol = float(np.prod(dxyz))
                for i, j, k in np.ndindex(*(n - 2 * ngc for n in shape)):
                    key = tuple(round(float(origin[a] + (n + ngc + 0.5) * dxyz[a]), 12)
                                for a, n in enumerate((i, j, k)))
                    out[key] = vol
    return out


def error(work: Path, ngc: int) -> float:
    """Return eps = sqrt(sum_v mean|q_v(T) - q_v(0)|^2) of a run."""
    first, last = profile(work, "first", ngc), profile(work, "last", ngc)
    if first.keys() != last.keys():
        sys.exit(f"linear_wave_oracle: the first and last checkpoints of {work} have different cells")
    keys = sorted(first)
    diff = np.abs(np.array([last[k] for k in keys]) - np.array([first[k] for k in keys]))
    return float(math.sqrt(np.sum(np.mean(diff, axis=0) ** 2)))


def main() -> int:
    """Run the MV-5 checks, print a report, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work", type=Path, nargs="+")
    parser.add_argument("--order-min", type=float, default=None, help="minimum observed order of the finest pair")
    parser.add_argument("--eps-max", type=float, nargs="+", default=None, help="bound of eps (one, or one per run)")
    parser.add_argument("--ngc", type=int, default=3)
    args = parser.parse_args()
    if args.eps_max is not None and len(args.eps_max) not in (1, len(args.work)):
        sys.exit("linear_wave_oracle: --eps-max takes one bound or one per run")
    ok = True
    eps = []
    for n, w in enumerate(args.work):
        e = error(w, args.ngc)
        eps.append(e)
        line = f"{w.name}: eps {e:.4e} (eps / amplitude {e / amplitude(w):.3e})"
        if len(eps) > 1:
            line += f"  order {math.log2(eps[-2] / e):+.2f}"
        if args.eps_max is not None:
            bound = args.eps_max[0] if len(args.eps_max) == 1 else args.eps_max[n]
            good = e <= bound
            ok &= good
            line += f"  {'PASS' if good else 'FAIL'} (max {bound:.3e})"
        print(line)
    if args.order_min is not None and len(eps) > 1:
        p = math.log2(eps[-2] / eps[-1])
        good = p >= args.order_min
        ok &= good
        print(f"   finest-pair order {p:+.2f} {'>=' if good else '<'} {args.order_min}: {'PASS' if good else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
