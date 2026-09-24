#!/usr/bin/env python3
"""FLUME multi-realm oracle (issue #37): a multi-realm run against the single-realm run on the same cell set.

Why this oracle exists: the realms of a forest manifest partition a domain whose union can be the domain of a
single-realm run (same cell size, same cell centres). With a mirror seam filled at every Runge-Kutta stage (beta
cadence) the seam is a block interface like any other, so the union must reproduce the single-realm run. Two checks:

* fields: the interior cells of every realm, keyed by their centre coordinates, against the single-realm cells at the
  same step: the maximum relative difference per variable (scale: the maximum magnitude of the reference) must not
  exceed the tolerance (0 = bitwise); the worst cell is reported by coordinates;
* conservation: the per-step sum over the realms of the volume integrals (`<basename>-conservation_history.dat`)
  against the single-realm integrals, within CONSERVATION_TOL: the realms sum their cells separately, so even
  bitwise-identical fields give integrals that differ by the summation order (measured 2e-13 on the 2-realm Sod).

Geometry (XH5F): origin and dxdydz stored (z, y, x), the origin at the corner of the first ghost cell; field datasets
(z, y, x) with ghosts.

Usage:
    multirealm_oracle.py <multi-work-dir> <single-work-dir> --ngc N [--tol T] [--step S]
"""

from __future__ import annotations

import argparse
import configparser
import sys
from pathlib import Path

import h5py
import numpy as np

VARIABLES = ("r", "ru", "rv", "rw", "rE")
CONSERVATION_TOL = 1.0e-12


def basenames(work: Path) -> list[str]:
    """Return the output basenames of the INI files of a work directory (manifests have none)."""
    names = []
    for ini_path in sorted(work.glob("*.ini")):
        ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None, strict=False)
        ini.read(ini_path)
        if ini.has_option("IO", "output_basename"):
            names.append(ini["IO"]["output_basename"].strip())
    return names


def last_step(work: Path, names: list[str]) -> int:
    """Return the last saved step common to every basename."""
    steps = None
    for name in names:
        found = {int(p.name.split("-")[-2]) for p in work.glob(f"{name}-*-proc*.h5")}
        steps = found if steps is None else steps & found
    if not steps:
        sys.exit(f"multirealm_oracle: no common checkpoint step in {work}")
    return max(steps)


def load_cells(work: Path, names: list[str], step: int, ngc: int) -> dict[tuple[float, float, float], np.ndarray]:
    """Return the interior cells of every basename at `step`, keyed by their rounded centre coordinates."""
    cells: dict[tuple[float, float, float], np.ndarray] = {}
    for name in names:
        for path in sorted(work.glob(f"{name}-{step:09d}-proc*.h5")):
            with h5py.File(path, "r") as h5:
                for blk in {k.rsplit("-", 1)[0] for k in h5}:
                    dx = h5[f"{blk}-dxdydz"][()][::-1]
                    lo = h5[f"{blk}-origin"][()][::-1] + ngc * dx
                    q = np.stack([h5[f"{blk}-{v}"][()].transpose(2, 1, 0) for v in VARIABLES])
                    q = q[:, ngc:-ngc, ngc:-ngc, ngc:-ngc]
                    for i in range(q.shape[1]):
                        for j in range(q.shape[2]):
                            for k in range(q.shape[3]):
                                key = (round(float(lo[0] + (i + 0.5) * dx[0]), 12),
                                       round(float(lo[1] + (j + 0.5) * dx[1]), 12),
                                       round(float(lo[2] + (k + 0.5) * dx[2]), 12))
                                cells[key] = q[:, i, j, k]
    return cells


def conservation_sum(work: Path, names: list[str]) -> np.ndarray:
    """Return the per-step sum over the basenames of the conservation histories (rows: it, time, integrals)."""
    total = None
    for name in names:
        data = np.loadtxt(work / f"{name}-conservation_history.dat", skiprows=1, ndmin=2)
        if total is None:
            total = data.copy()
        else:
            n = min(len(total), len(data))
            total = total[:n]
            total[:, 2:] += data[:n, 2:]
    return total


def main() -> int:
    """Run the checks, print a report, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("multi", type=Path, help="work directory of the multi-realm run")
    parser.add_argument("single", type=Path, help="work directory of the single-realm run")
    parser.add_argument("--ngc", type=int, required=True, help="ghost cells number")
    parser.add_argument("--tol", type=float, default=0.0, help="relative tolerance (0 = bitwise)")
    parser.add_argument("--step", type=int, default=None, help="step to compare (default: the last common one)")
    args = parser.parse_args()

    multi_names, single_names = basenames(args.multi), basenames(args.single)
    last_multi, last_single = last_step(args.multi, multi_names), last_step(args.single, single_names)
    if args.step is None and last_multi != last_single:
        print(f"last saved steps differ: multi-realm {last_multi}, single-realm {last_single} (different time steps)  FAIL")
        return 1
    step = args.step if args.step is not None else last_single
    multi = load_cells(args.multi, multi_names, step, args.ngc)
    single = load_cells(args.single, single_names, step, args.ngc)
    status = 0

    missing = set(single) ^ set(multi)
    if missing:
        print(f"cell sets differ: {len(missing)} cells not in both runs  FAIL")
        return 1
    keys = sorted(single)
    a = np.array([multi[key] for key in keys])
    b = np.array([single[key] for key in keys])
    scale = np.max(np.abs(b), axis=0)
    scale = np.where(scale > 0.0, scale, 1.0)
    rel = np.abs(a - b) / scale
    for v, name in enumerate(VARIABLES):
        worst = int(np.argmax(rel[:, v]))
        ok = rel[worst, v] <= args.tol
        status |= 0 if ok else 1
        where = "" if rel[worst, v] == 0.0 else f", worst at x = {keys[worst][0]:.6f}"
        print(f"step {step} {name:>2}: {len(keys)} cells, max relative difference {rel[worst, v]:.3e}{where}"
              f"  {'PASS' if ok else 'FAIL'}")

    cm, cs = conservation_sum(args.multi, multi_names), conservation_sum(args.single, single_names)
    n = min(len(cm), len(cs))
    cscale = np.max(np.abs(cs[:n, 2:]), axis=0)
    cscale = np.where(cscale > 0.0, cscale, 1.0)
    cdiff = np.max(np.abs(cm[:n, 2:] - cs[:n, 2:]) / cscale, axis=0)
    ok = bool(np.all(cdiff <= max(args.tol, CONSERVATION_TOL)))
    status |= 0 if ok else 1
    print(f"conservation: sum over realms vs single realm, {n} steps, max relative difference per variable "
          f"{' '.join(f'{d:.2e}' for d in cdiff)}  {'PASS' if ok else 'FAIL'}")
    return status


if __name__ == "__main__":
    sys.exit(main())
