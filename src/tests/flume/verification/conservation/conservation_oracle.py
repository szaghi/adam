#!/usr/bin/env python3
"""FLUME verification V3 oracle: conservation across AMR coarse-fine faces (issue #35, section 11).

Why this oracle exists: a finite-volume update conserves only if every interior face flux leaves one cell and enters
its neighbour. At a 2:1 coarse-fine face the two sides compute different fluxes (different stencils, different ghost
fills), so the grid leaks unless the Berger-Colella reflux replaces the coarse flux with the restricted fine one. In a
periodic box nothing crosses the boundary: every volume integral must stay constant to round-off.

The run without reflux is the negative control: it must drift well above round-off, otherwise the case does not
exercise the seams and a passing reflux leg proves nothing (af54afee pattern).

Usage:
    conservation_oracle.py --conserved <work-dir> [--max-drift D] [--leaky <work-dir> --min-drift M]
    conservation_oracle.py --compare <work-dir-a> <work-dir-b> --tol T

--compare checks two runs of the same case (e.g. CPU and FNL, V4) cell by cell on the last checkpoint, blocks matched
by their origin (rank-independent); the difference of each variable is relative to its largest magnitude, so that one
tolerance fits density and energy alike.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np

VARIABLES = ("r", "ru", "rv", "rw", "rE")


def drift(work: Path) -> np.ndarray:
    """Return the maximum relative drift of each volume integral of the conservation history of one run."""
    hist = next(work.glob("*-conservation_history.dat"))
    rows = [line.split() for line in hist.open() if line.split() and line.split()[0][0] in "+-0123456789"]
    h = np.array([[float(v) for v in r] for r in rows])
    if len(h) < 2:
        sys.exit(f"conservation_oracle: {hist} holds fewer than two rows")
    return np.max(np.abs(h[:, 2:] - h[0, 2:]), axis=0) / np.abs(h[0, 2:])


def last_fields(work: Path) -> dict[tuple[float, ...], np.ndarray]:
    """Return the conservative fields of the last checkpoint of one run, keyed by block origin."""
    files = sorted(p for p in work.glob("*-proc*.h5") if "restart" not in p.name)
    last = max(int(p.name.split("-")[-2]) for p in files)
    out = {}
    for path in (p for p in files if int(p.name.split("-")[-2]) == last):
        with h5py.File(path, "r") as h5:
            for blk in {k.rsplit("-", 1)[0] for k in h5}:
                key = tuple(float(x) for x in np.round(h5[f"{blk}-origin"][()], 12))
                out[key] = np.stack([h5[f"{blk}-{v}"][()] for v in VARIABLES])
    return out


def report(work: Path, d: np.ndarray) -> str:
    """Return the drift line of one run."""
    return f"{work.name}: relative drift " + " ".join(f"{n} {x:.2e}" for n, x in zip(VARIABLES, d, strict=True))


def main() -> int:
    """Run the requested checks, print a report, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--conserved", type=Path, help="run with reflux: every integral must be conserved")
    parser.add_argument("--max-drift", type=float, default=1.0e-13, help="bound of the conserved run")
    parser.add_argument("--leaky", type=Path, help="run without reflux (negative control): it must drift")
    parser.add_argument("--min-drift", type=float, default=1.0e-10, help="lower bound of the negative control")
    parser.add_argument("--compare", type=Path, nargs=2, help="two runs to compare cell by cell")
    parser.add_argument("--tol", type=float, default=0.0, help="relative comparison tolerance (default 0: bitwise)")
    args = parser.parse_args()

    status = 0
    if args.conserved is not None:
        d = drift(args.conserved)
        ok = bool(np.all(d <= args.max_drift))
        status |= 0 if ok else 1
        print(f"{report(args.conserved, d)}  {'PASS' if ok else 'FAIL'} (<= {args.max_drift:.1e})")
    if args.leaky is not None:
        d = drift(args.leaky)
        ok = bool(np.max(d) >= args.min_drift)
        status |= 0 if ok else 1
        print(f"{report(args.leaky, d)}  {'PASS' if ok else 'FAIL'} (negative control, max >= {args.min_drift:.1e})")
    if args.compare is not None:
        a, b = (last_fields(w) for w in args.compare)
        if a.keys() != b.keys():
            sys.exit("conservation_oracle: the two runs hold different block sets")
        scale = np.max([np.max(np.abs(a[k]), axis=(1, 2, 3)) for k in a], axis=0)
        delta = np.max([np.max(np.abs(a[k] - b[k]), axis=(1, 2, 3)) for k in a], axis=0)
        diff = float(np.max(delta / scale))
        ok = diff <= args.tol
        status |= 0 if ok else 1
        kind = "bitwise" if args.tol == 0.0 else f"tol {args.tol:.1e}"
        print(f"{args.compare[0].name} vs {args.compare[1].name}: {len(a)} blocks, max relative difference {diff:.3e}"
              f"  {'PASS' if ok else 'FAIL'} ({kind})")
    return status


if __name__ == "__main__":
    sys.exit(main())
