#!/usr/bin/env python3
"""FLUME NV-6 oracle (issue #49, N2d): Sod in interstellar units, limiting active at the shock.

Why: the classic WENO weights `d / (zeps + IS)^p` carry an absolute `zeps = 1e-6` (issue #48, A3). `IS` scales as the
square of the reconstructed variable, so with a density of 1e-21 kg/m^3 it is ~1e-42, `zeps` dominates, the weights
collapse to the linear ones and the scheme stops limiting; the `weno-riemann` sensor reads the same weights and goes
blind too. Two cures exist since #49: the scale-invariant weights (`[weno] weights = si`, N1) and the reference layer
(`[reference]`, N2), which hands the solver O(1) numbers.

Every run is the Sod problem; `physical` runs are written in SI interstellar units (`scaling.py physical`). The
profiles are compared non-dimensionally (x / L0, rho / rho0) after the last step:

* `--same A B`: run B equals run A within `--tol` (maximum difference of rho / rho0): limiting is the one of A;
* `--collapse A B`: the total variation of B exceeds the one of A by more than `--tv-ratio` (negative control: the case
  detects a scheme that does not limit).

Usage:
    interstellar_oracle.py --ini <sod.ini> --references L0 U0 RHO0 (--same | --collapse) <A> <B> [--tol T]
                           [--tv-ratio R]

A work directory whose name contains `-physical` is in SI units (code units otherwise: output_units = code).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "sod"))
from sod_oracle import exact_riemann_density, load_profile, read_ini


def profile(work: Path, ngc: int, references: tuple[float, float, float]) -> tuple[int, np.ndarray, np.ndarray]:
    """Last iteration, x / L0 and rho / rho0 of a run (code-unit runs are already non-dimensional)."""
    basename = next(work.glob("*-residuals.dat")).name.removesuffix("-residuals.dat")
    it, _, xs, q = load_profile(work, basename, ngc)
    length, _, density = references if "-physical" in work.name else (1.0, 1.0, 1.0)
    return it, xs / length, q[0] / density


def main() -> int:
    """Run one NV-6 assertion, print a report, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--ini", type=Path, required=True, help="the code-unit Sod input (physics, time, IC)")
    parser.add_argument("--references", type=float, nargs=3, required=True, metavar=("L0", "U0", "RHO0"))
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--same", nargs=2, type=Path, metavar=("A", "B"))
    mode.add_argument("--collapse", nargs=2, type=Path, metavar=("A", "B"))
    parser.add_argument("--tol", type=float, default=1.0e-12, help="--same: max |rho_B - rho_A| / rho0")
    parser.add_argument("--tv-ratio", type=float, default=1.1, help="--collapse: TV_B / TV_A lower bound")
    args = parser.parse_args()

    ini = read_ini(args.ini)
    gamma = float(ini["physics"]["cp"]) / float(ini["physics"]["cv"])
    t = float(ini["time"]["time_max"])
    ngc = int(ini["grid"]["ngc"])
    r1, r2 = ini["initial_conditions_region_1"], ini["initial_conditions_region_2"]
    left = (float(r1["r"]), float(r1["u"]), float(r1["p"]))
    right = (float(r2["r"]), float(r2["u"]), float(r2["p"]))
    works = args.same or args.collapse
    report = []
    for work in works:
        it, x, rho = profile(work, ngc, tuple(args.references))
        exact = exact_riemann_density(x, t, 0.5, gamma, left, right)
        l1 = float(np.sum(np.abs(rho - exact)) * (x[1] - x[0]))
        tv = float(np.sum(np.abs(np.diff(rho))))
        new_extrema = max(float(rho.max()) - left[0], 0.0) + max(right[0] - float(rho.min()), 0.0)
        report.append((work.name, it, x, rho, l1, tv, new_extrema))
    for name, it, _, _, l1, tv, ext in report:
        print(f"   {name}: it {it}, L1(rho/rho0) {l1:.6e}, TV {tv:.6f}, new extrema {ext:.2e}")
    (_, _, xa, ra, _, tva, _), (_, _, xb, rb, _, tvb, _) = report
    if args.same:
        diff = float(np.max(np.abs(rb - ra))) if xa.shape == xb.shape else float("inf")
        ok = diff <= args.tol
        verdict = "PASS" if ok else "FAIL"
        print(f"   same limiting: max |rho_B - rho_A| / rho0 = {diff:.2e}  {verdict} (tol {args.tol:.0e})")
    else:
        ratio = tvb / tva
        ok = ratio > args.tv_ratio
        print(f"   collapse detected: TV_B / TV_A = {ratio:.3f}  {'PASS' if ok else 'FAIL'} (> {args.tv_ratio})")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
