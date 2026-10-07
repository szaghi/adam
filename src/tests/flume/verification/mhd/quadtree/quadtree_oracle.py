#!/usr/bin/env python3
"""FLUME quadtree AMR oracle (issue #46, MV-15): 2-D runs on different trees must give the same (x, y) solution.

Why this oracle exists: a quadtree (ratio 4) refines x and y only, so its 2:1 seams are 1:1 in z; until #46 the seam
machinery assumed an octree and the coarse cells beside a seam picked up a spurious z dependence. The same 2-D problem
run on an octree with null z (the reference, the layout the AMR variants used before #46), on a quadtree with nk = 1
and on a quadtree with an active z axis (nk = 4) holds the same cells in (x, y). Two checks, per variable:

* z invariance: within each run, the cells of one (x, y) column differ by at most --z-tol (relative to the maximum
  magnitude of the variable in the reference): the extrapolated z boundaries keep a z-invariant start z-invariant;
* agreement: the column value of every run against the reference at the same (x, y), within
  max(--tol, --spread-factor * the reference's own z spread of the variable). The trees order the cells differently (an
  octree restriction averages 8 fine cells, 2 identical layers, a quadtree 4), so the runs agree to round-off, not bit
  for bit. Even the reference is z-invariant to round-off only (its tricubic coarse-fine weights differ by z
  sub-position); on a problem whose limiter switches amplify round-off, the reference's own spread across identical
  layers is the floor of any comparison, and --spread-factor ties the bound to it.

Every compared run must hold the same (x, y) cell set and the same last saved step (the same time step sequence). With
no run under test only the z invariance of the reference is checked (a run whose time steps differ, e.g. an active z
axis entering the CFL condition).

Geometry (XH5F): origin and dxdydz stored (z, y, x), the origin at the corner of the first ghost cell; field datasets
(z, y, x) with ghosts.

Usage:
    quadtree_oracle.py <reference-work-dir> [<work-dir> ...] --ngc N [--tol T] [--z-tol T] [--spread-factor F]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np

GEOMETRY = ("origin", "dxdydz", "time_iteration")


def checkpoints(work: Path) -> list[Path]:
    """Return the field checkpoints of a work directory (restart files excluded)."""
    return [p for p in work.glob("*-[0-9]*-proc*.h5") if "-restart-" not in p.name]


def last_step(work: Path) -> int:
    """Return the last saved step of a work directory."""
    steps = {int(p.name.split("-")[-2]) for p in checkpoints(work)}
    if not steps:
        sys.exit(f"quadtree_oracle: no checkpoint in {work}")
    return max(steps)


def load_columns(work: Path, step: int, ngc: int) -> tuple[tuple[str, ...], dict[tuple[float, float], np.ndarray]]:
    """Return the field names and, per (x, y) column, the interior cells of the column (shape: cells, variables)."""
    columns: dict[tuple[float, float], list[np.ndarray]] = {}
    variables: tuple[str, ...] = ()
    for path in sorted(p for p in checkpoints(work) if p.name.split("-")[-2] == f"{step:09d}"):
        with h5py.File(path, "r") as h5:
            if not len(h5):  # a rank that owns no block
                continue
            found = tuple(sorted({k.rsplit("-", 1)[1] for k in h5} - set(GEOMETRY)))
            if not variables:
                variables = found
            elif found != variables:
                sys.exit(f"quadtree_oracle: {path} holds fields {found}, not {variables}")
            for blk in {k.rsplit("-", 1)[0] for k in h5}:
                dx = h5[f"{blk}-dxdydz"][()][::-1]
                lo = h5[f"{blk}-origin"][()][::-1] + ngc * dx
                q = np.stack([h5[f"{blk}-{v}"][()].transpose(2, 1, 0) for v in variables])
                q = q[:, ngc:-ngc, ngc:-ngc, ngc:-ngc]
                for i in range(q.shape[1]):
                    for j in range(q.shape[2]):
                        key = (round(float(lo[0] + (i + 0.5) * dx[0]), 12), round(float(lo[1] + (j + 0.5) * dx[1]), 12))
                        columns.setdefault(key, []).extend(q[:, i, j, k] for k in range(q.shape[3]))
    return variables, {key: np.array(cells) for key, cells in columns.items()}


def per_variable(variables: tuple[str, ...], values: np.ndarray) -> str:
    """Return `name value` pairs, one per variable."""
    return " ".join(f"{v} {x:.1e}" for v, x in zip(variables, values, strict=True))


def main() -> int:
    """Run the checks, print a report, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("reference", type=Path, help="work directory of the reference run (octree, null z)")
    parser.add_argument("works", type=Path, nargs="*", help="work directories of the runs under test")
    parser.add_argument("--ngc", type=int, required=True, help="ghost cells number")
    parser.add_argument("--tol", type=float, default=0.0, help="relative tolerance against the reference")
    parser.add_argument("--z-tol", type=float, default=0.0, help="relative tolerance of the z invariance")
    parser.add_argument("--spread-factor", type=float, default=0.0,
                        help="agreement bound as a multiple of the reference z spread (per variable)")
    args = parser.parse_args()

    step = last_step(args.reference)
    variables, ref = load_columns(args.reference, step, args.ngc)
    scale = np.max(np.abs(np.concatenate(list(ref.values()))), axis=0)
    scale = np.where(scale > 0.0, scale, 1.0)
    keys = sorted(ref)
    status = 0
    ref_spread = np.zeros(len(variables))
    for work in [args.reference, *args.works]:
        name = work.name
        if last_step(work) != step:
            print(f"{name}: last saved step {last_step(work)}, reference {step} (different time steps)  FAIL")
            status = 1
            continue
        found, cols = (variables, ref) if work == args.reference else load_columns(work, step, args.ngc)
        if found != variables or set(cols) != set(ref):
            print(f"{name}: fields {found} or (x, y) cell set ({len(cols)} columns, reference {len(ref)}) differ  FAIL")
            status = 1
            continue
        spread = np.max([np.max(np.abs(cols[key] - cols[key][0]), axis=0) for key in keys], axis=0) / scale
        ok = bool(np.all(spread <= args.z_tol))
        status |= 0 if ok else 1
        nz = sorted({len(cols[key]) for key in keys})
        print(f"{name}: step {step}, {len(keys)} columns of {nz} cells, z spread {per_variable(variables, spread)}"
              f"  {'PASS' if ok else 'FAIL'}")
        if work == args.reference:
            ref_spread = spread
            continue
        diff = np.max([np.abs(cols[key][0] - ref[key][0]) for key in keys], axis=0) / scale
        ok = bool(np.all(diff <= np.maximum(args.tol, args.spread_factor * ref_spread)))
        status |= 0 if ok else 1
        print(f"{name}: against the reference, relative difference {per_variable(variables, diff)}"
              f"  {'PASS' if ok else 'FAIL'}")
    return status


if __name__ == "__main__":
    sys.exit(main())
