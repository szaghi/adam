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
    conservation_oracle.py --compare <work-dir-a> <work-dir-b> --tol T [--ngc N] [--nonzero NAME ...]

--compare checks two runs of the same case (e.g. CPU and FNL, V4) cell by cell on the last checkpoint, blocks matched
by their origin (rank-independent), every field the checkpoints hold (Euler r, ru, rv, rw, rE; MHD adds bx, by, bz,
psi); the difference of each variable is relative to its largest magnitude, so that one tolerance fits density and
energy alike (an identically zero variable is compared absolutely). `--ngc N` strips N ghost layers first (the digest
semantics): edge and corner ghost cells that no map fills keep stale values, which the directional stencils never
read. `--nonzero NAME` fails the comparison when a field is identically zero in the first run (a restart that must
carry psi proves nothing on a zero psi).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np

VARIABLES = ("r", "ru", "rv", "rw", "rE")
GEOMETRY = ("origin", "dxdydz", "time_iteration")


def drift(work: Path) -> tuple[tuple[str, ...], np.ndarray]:
    """Return (names, maximum relative drift) of each volume integral of the conservation history of one run.

    The names come from the history header (`int_r` -> `r`; MHD adds `bx, by, bz [, psi]`). An integral that starts at
    exactly zero (psi) has no scale of its own: its drift is relative to the largest initial integral."""
    hist = next(work.glob("*-conservation_history.dat"))
    lines = hist.read_text().splitlines()
    rows = [line.split() for line in lines if line.split() and line.split()[0][0] in "+-0123456789"]
    h = np.array([[float(v) for v in r] for r in rows])
    if len(h) < 2:
        sys.exit(f"conservation_oracle: {hist} holds fewer than two rows")
    header = [line for line in lines if line.startswith("VARIABLES=")]
    names = VARIABLES
    if header:
        names = tuple(n.strip('"').removeprefix("int_") for n in header[0].split("=", 1)[1].split()[2:])
    h0 = np.abs(h[0, 2:])
    scale = np.where(h0 > 0.0, h0, np.max(h0))
    return names, np.max(np.abs(h[:, 2:] - h[0, 2:]), axis=0) / scale


def last_fields(work: Path) -> tuple[tuple[str, ...], dict[tuple[float, ...], np.ndarray]]:
    """Return the field names and the fields of the last checkpoint of one run, keyed by block origin."""
    files = sorted(p for p in work.glob("*-proc*.h5") if "restart" not in p.name)
    last = max(int(p.name.split("-")[-2]) for p in files)
    out = {}
    names: tuple[str, ...] = ()
    for path in (p for p in files if int(p.name.split("-")[-2]) == last):
        with h5py.File(path, "r") as h5:
            if not len(h5):  # a rank that owns no block (issue #42)
                continue
            found = tuple(sorted({k.rsplit("-", 1)[1] for k in h5} - set(GEOMETRY)))
            if names and found != names:
                sys.exit(f"conservation_oracle: {path} holds fields {found}, not {names}")
            names = found
            for blk in {k.rsplit("-", 1)[0] for k in h5}:
                key = tuple(float(x) for x in np.round(h5[f"{blk}-origin"][()], 12))
                out[key] = np.stack([h5[f"{blk}-{v}"][()] for v in names])
    return names, out


def report(work: Path, names: tuple[str, ...], d: np.ndarray) -> str:
    """Return the drift line of one run."""
    return f"{work.name}: relative drift " + " ".join(f"{n} {x:.2e}" for n, x in zip(names, d, strict=True))


def main() -> int:
    """Run the requested checks, print a report, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--conserved", type=Path, help="run with reflux: every integral must be conserved")
    parser.add_argument("--max-drift", type=float, default=1.0e-13, help="bound of the conserved run")
    parser.add_argument("--leaky", type=Path, help="run without reflux (negative control): it must drift")
    parser.add_argument("--min-drift", type=float, default=1.0e-10, help="lower bound of the negative control")
    parser.add_argument("--exclude", nargs="+", default=[], help="integrals reported but not bounded (e.g. psi)")
    parser.add_argument("--compare", type=Path, nargs=2, help="two runs to compare cell by cell")
    parser.add_argument("--tol", type=float, default=0.0, help="relative comparison tolerance (default 0: bitwise)")
    parser.add_argument("--ngc", type=int, default=0, help="ghost layers stripped before the comparison")
    parser.add_argument("--nonzero", nargs="+", default=[], help="fields that must not be identically zero (--compare)")
    args = parser.parse_args()

    status = 0
    if args.conserved is not None:
        names, d = drift(args.conserved)
        bounded = np.array([n not in args.exclude for n in names])
        ok = bool(np.all(d[bounded] <= args.max_drift))
        status |= 0 if ok else 1
        excluded = f", {' '.join(args.exclude)} not bounded" if args.exclude else ""
        print(f"{report(args.conserved, names, d)}  {'PASS' if ok else 'FAIL'} (<= {args.max_drift:.1e}{excluded})")
    if args.leaky is not None:
        names, d = drift(args.leaky)
        ok = bool(np.max(d) >= args.min_drift)
        status |= 0 if ok else 1
        print(f"{report(args.leaky, names, d)}  {'PASS' if ok else 'FAIL'} "
              f"(negative control, max >= {args.min_drift:.1e})")
    if args.compare is not None:
        (names_a, a), (names_b, b) = (last_fields(w) for w in args.compare)
        if names_a != names_b:
            sys.exit(f"conservation_oracle: the two runs hold different fields, {names_a} and {names_b}")
        if args.ngc > 0:
            g = args.ngc
            a = {k: v[:, g:-g, g:-g, g:-g] for k, v in a.items()}
            b = {k: v[:, g:-g, g:-g, g:-g] for k, v in b.items()}
        if a.keys() != b.keys():
            sys.exit("conservation_oracle: the two runs hold different block sets")
        scale = np.max([np.max(np.abs(a[k]), axis=(1, 2, 3)) for k in a], axis=0)
        scale = np.where(scale > 0.0, scale, 1.0)  # an identically zero variable is compared absolutely
        delta = np.max([np.max(np.abs(a[k] - b[k]), axis=(1, 2, 3)) for k in a], axis=0)
        diff = float(np.max(delta / scale))
        ok = diff <= args.tol
        status |= 0 if ok else 1
        kind = "bitwise" if args.tol == 0.0 else f"tol {args.tol:.1e}"
        print(f"{args.compare[0].name} vs {args.compare[1].name}: {len(a)} blocks, {len(names_a)} fields "
              f"({' '.join(names_a)}), max relative difference {diff:.3e}  {'PASS' if ok else 'FAIL'} ({kind})")
        for name in args.nonzero:
            if name not in names_a:
                sys.exit(f"conservation_oracle: --nonzero {name}: no such field in {args.compare[0]}")
            top = float(np.max([np.max(np.abs(a[k][names_a.index(name)])) for k in a]))
            ok = top > 0.0
            status |= 0 if ok else 1
            print(f"   max |{name}| {top:.3e}  {'PASS' if ok else 'FAIL'} (must not be identically zero)")
    return status


if __name__ == "__main__":
    sys.exit(main())
