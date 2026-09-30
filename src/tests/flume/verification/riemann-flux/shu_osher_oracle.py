#!/usr/bin/env python3
"""FLUME verification RV-3 oracle: the Shu-Osher shock-density wave interaction along x, y, z (issue #47).

Why this oracle exists: the Shu-Osher problem (Shu & Osher 1989, example 8) has no closed-form solution, and it is the
standard test of whether a shock-capturing scheme keeps the high-frequency entropy waves behind the shock. RV-3 asks
the `weno-riemann` scheme to be no worse than the flux-splitting scheme there. The reference is the splitting scheme of
`hybrid_proto.py` (FLUME's `weno` scheme in NumPy) on REF_FACTOR times finer cells, averaged to the grid of the runs;
it is computed once and cached in the temporary directory (a few minutes). A split-scheme reference biases the
comparison toward the split scheme by the reference's own error, about 1/REF_FACTOR of the coarse error at the shocks.

Checks, per run: L1(rho) against the reference (bound `--l1-max`), positive density and pressure. The y and z runs
must equal the x run bitwise after the axes permutation, every transverse copy bitwise identical (as V1).

Usage:
    shu_osher_oracle.py <case.ini> <work-dir> [<work-dir> ...] [--l1-max L]
"""

from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "sod"))
import hybrid_proto  # noqa: E402
from sod_oracle import load_profile, read_ini  # noqa: E402

REF_FACTOR = 16


def reference_density(cells: int, t_end: float) -> np.ndarray:
    """Return the reference density averaged to `cells` cells of [-5, 5] at `t_end` (cached)."""
    n = REF_FACTOR * cells
    cache = Path(tempfile.gettempdir()) / f"flume-shu-osher-ref-{n}-t{t_end:g}.npy"
    try:
        rho = np.load(cache)
    except OSError:
        u0, h, _, bc = hybrid_proto.problem("shu-osher", n)
        u = hybrid_proto.run(u0.copy(), h, t_end, {"scheme": "split", "corr_off": 0}, bc)
        if u is None:
            sys.exit("shu_osher_oracle: the reference run failed")
        rho = u[0]
        np.save(cache, rho)
    return rho.reshape(cells, REF_FACTOR).mean(axis=1)


def main() -> int:
    """Run the checks, print a report, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("ini", type=Path, help="the x input of the runs (time, grid)")
    parser.add_argument("work", type=Path, nargs="+", help="work directories of the x (and y, z) runs")
    parser.add_argument("--l1-max", type=float, default=None, help="L1(rho) bound (per run)")
    args = parser.parse_args()

    ini = read_ini(args.ini)
    gamma = float(ini["physics"]["cp"]) / float(ini["physics"]["cv"])
    if abs(gamma - hybrid_proto.GAMMA) > 1e-6:
        sys.exit(f"shu_osher_oracle: gamma {gamma} differs from the reference's {hybrid_proto.GAMMA}")
    t = float(ini["time"]["time_max"])
    ngc = int(ini["grid"]["ngc"])

    status, profiles, rho_ref = 0, [], None
    for work in args.work:
        basename = next(work.glob("*-residuals.dat")).name.removesuffix("-residuals.dat")
        it, axis, xs, q = load_profile(work, basename, ngc)
        if rho_ref is None:
            rho_ref = reference_density(xs.size, t)
        dx = xs[1] - xs[0]
        l1 = float(np.sum(np.abs(q[0] - rho_ref)) * dx)
        rho = q[0]
        p = (gamma - 1.0) * (q[4] - 0.5 * (q[1] ** 2 + q[2] ** 2 + q[3] ** 2) / rho)
        ok = bool(np.all(rho > 0.0) and np.all(p > 0.0))
        verdict = ""
        if args.l1_max is not None:
            ok = ok and l1 <= args.l1_max
            verdict = f"  (bound {args.l1_max:.3e})"
        status |= 0 if ok else 1
        print(f"{work.name}/{basename}: axis {'xyz'[axis]}, it {it}, cells {xs.size}, L1(rho) = {l1:.6e}, "
              f"max rho {rho.max():.4f} (ref {rho_ref.max():.4f}), min p {p.min():.4e}  "
              f"{'PASS' if ok else 'FAIL'}{verdict}")
        mom = np.roll(q[1:4], -axis, axis=0)
        profiles.append((f"{work.name}/{basename}", np.vstack([q[0:1], mom, q[4:5]])))
    for name, prof in profiles[1:]:
        diff = float(np.max(np.abs(prof - profiles[0][1])))
        ok = diff == 0.0
        status |= 0 if ok else 1
        print(f"{name} vs {profiles[0][0]}: max |difference| = {diff:.3e}  {'PASS' if ok else 'FAIL'} (bitwise)")
    return status


if __name__ == "__main__":
    sys.exit(main())
