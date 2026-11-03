#!/usr/bin/env python3
"""Check the channel runs of issue #119 T3. `laminar <work>...`, coarse to fine: the driving force against 3 mu U / h^2 and its order, the mass conserved; the bulk momentum is held by construction, so not checked.
`same a b --tol`: the force of two runs. `agree a b`: the histories of two runs. `energy a b --tol`: the change of the total energy of two runs."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np


def force(work: Path, last: float = 0.25) -> tuple[float, float]:
    """Return the mean driving force over the last part of the history, and the drift of the mass over the run."""
    f = np.loadtxt(work / "channel-forcing_history.dat", ndmin=2)
    tail = f[int(len(f) * (1.0 - last)) :, 2]
    path = work / "channel-conservation_history.dat"
    with open(path) as fh:
        names = [w.strip('"') for w in fh.readline().split("=")[1].split()]
    c = np.loadtxt(path, skiprows=1, ndmin=2)[:, names.index("int_r")]
    return float(tail.mean()), float(np.max(np.abs(c - c[0])) / abs(c[0]))


def laminar(works: list[str], reynolds: float, tol: float, min_order: float, mass_tol: float, height: float = 0.0) -> int:
    """Laminar Poiseuille: f = 3 mu U / h^2, here 3 / Re. With a wall model matched at a fixed height h (a laminar layer, linear below h: tau = mu u(h) / h),
    the steady flow is a parabola A (1 - y^2) plus a slip U_s; tau = 2 mu A, U_s = A h^2 and the bulk U_s + 2A/3 = 1 give f = 3 / Re / (1 + 3 h^2 / 2)."""
    exact = 3.0 / reynolds / (1.0 + 1.5 * height * height)
    errs = []
    ok = True
    for w in works:
        f, drift = force(Path(w))
        errs.append(abs(f - exact) / exact)
        ok = ok and drift < mass_tol
        print(f"{w.rsplit('-', 1)[-1]}: f {f:.6e} against {exact:.6e}, relative {errs[-1]:.2e}; mass drift {drift:.2e} (bound {mass_tol:.0e})")
    order = float(np.log2(errs[-2] / errs[-1])) if len(errs) > 1 else 0.0
    ok = ok and errs[-1] < tol and (len(errs) == 1 or order > min_order)
    print(f"order {order:.2f} (min {min_order}), finest {errs[-1]:.2e} (bound {tol:.0e})  {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def history(path: Path, skip: int) -> np.ndarray:
    """Load a history."""
    return np.loadtxt(path, skiprows=skip, ndmin=2)


def agree(a: str, b: str) -> int:
    """Two runs on the same step sequence: forcing history within 1e-9 (a mean differenced to round-off over dt), conservation within 1e-11; a floor of 1, the transverse momenta being round-off about 0."""
    ok = True
    for name, skip, tol in (("channel-forcing_history.dat", 0, 1e-9), ("channel-conservation_history.dat", 1, 1e-11)):
        x, y = history(Path(a) / name, skip), history(Path(b) / name, skip)
        d = float(np.max(np.abs(x - y)[:, 2:] / np.maximum(np.abs(x[:, 2:]).max(axis=0), 1.0)))
        good = d < tol
        ok = ok and good
        print(f"{name}: {d:.2e} (bound {tol:.0e})  {'PASS' if good else 'FAIL'}")
    return 0 if ok else 1


def same_force(a: str, b: str, tol: float) -> int:
    """The force of a run against a reference run."""
    fa, da = force(Path(a))
    fb, _ = force(Path(b))
    rel = abs(fa - fb) / abs(fb)
    ok = rel < tol and da < 1e-12
    print(f"force {fa:.6e} against {fb:.6e}, relative {rel:.2e} (bound {tol:.0e}); mass drift {da:.2e} (bound 1e-12)  {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def differ_force(a: str, b: str, tol: float) -> int:
    """Two runs whose forces must differ by more than tol, relative: an option that acts."""
    fa, _ = force(Path(a))
    fb, _ = force(Path(b))
    rel = abs(fa - fb) / abs(fb)
    ok = rel > tol
    print(f"force {fa:.6e} against {fb:.6e}, relative {rel:.2e} (min {tol:.0e})  {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def between(a: str, coarse: str, fine: str) -> int:
    """The force of a run with half the wall refined, strictly between those of the coarse and the fine uniform runs; its mass conserved."""
    fa, da = force(Path(a))
    fc, _ = force(Path(coarse))
    ff, _ = force(Path(fine))
    ok = min(fc, ff) < fa < max(fc, ff) and da < 1e-12
    print(f"force {fa:.6e}, coarse {fc:.6e}, fine {ff:.6e}, mean of the two {(fc + ff) / 2:.6e}; mass drift {da:.2e} (bound 1e-12)  {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def energy(a: str, b: str, tol: float) -> int:
    """The change of the total energy of a run against a reference run."""
    out = []
    for w in (a, b):
        path = Path(w) / "channel-conservation_history.dat"
        with open(path) as fh:
            names = [v.strip('"') for v in fh.readline().split("=")[1].split()]
        c = np.loadtxt(path, skiprows=1, ndmin=2)[:, names.index("int_rE")]
        out.append(float(c[-1] - c[0]))
    rel = abs(out[0] - out[1]) / abs(out[1])
    ok = rel < tol and out[0] * out[1] > 0.0
    print(f"change of the total energy {out[0]:.4e} against {out[1]:.4e}, relative {rel:.2e} (bound {tol:.0e})  {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def main() -> int:
    """Check."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("what", choices=["laminar", "same", "agree", "energy", "between", "differ"])
    ap.add_argument("works", nargs="+")
    ap.add_argument("--reynolds", type=float, default=100.0)
    ap.add_argument("--tol", type=float, default=0.02)
    ap.add_argument("--min-order", type=float, default=1.5)
    ap.add_argument("--mass-tol", type=float, default=1e-12)
    ap.add_argument("--height", type=float, default=0.0)
    args = ap.parse_args()
    if args.what == "same":
        return same_force(args.works[0], args.works[1], args.tol)
    if args.what == "differ":
        return differ_force(args.works[0], args.works[1], args.tol)
    if args.what == "between":
        return between(args.works[0], args.works[1], args.works[2])
    if args.what == "agree":
        return agree(args.works[0], args.works[1])
    if args.what == "energy":
        return energy(args.works[0], args.works[1], args.tol)
    return laminar(args.works, args.reynolds, args.tol, args.min_order, args.mass_tol, args.height)


if __name__ == "__main__":
    sys.exit(main())
