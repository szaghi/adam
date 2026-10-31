#!/usr/bin/env python3
"""Check the subgrid-scale model (issue #119, T1).

formula   the Fortran formulas (test_sgs, the executable given) against an independent computation: Vreman and WALE from the
          tensor expressions by numpy, Sigma from a singular value decomposition, on random tensors; and the zeros each model is
          built for: Vreman where a_ij has rank 1, WALE in pure shear, Sigma in pure shear, 2-D flow and solid rotation.
budget    dK/dt = pd - eps - eps_sgs of a ladder (argument n:work, coarse to fine): what is left must fall with n, and at the finest
          the budget must close better with eps_sgs than without it or with its sign flipped, so the column measures what the
          fluxes do.
inviscid  a run with no viscosity: eps exactly zero, eps_sgs positive, mass and total energy conserved.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np


def vreman(c: float, dl: np.ndarray, g: np.ndarray) -> float:
    """Vreman (2004), a_ij = d u_j / d x_i."""
    a = g.T
    b = (a.T * dl**2) @ a
    bb = b[0, 0] * b[1, 1] - b[0, 1] ** 2 + b[0, 0] * b[2, 2] - b[0, 2] ** 2 + b[1, 1] * b[2, 2] - b[1, 2] ** 2
    aa = float(np.sum(a * a))
    return c * np.sqrt(bb / aa) if aa > 0.0 and bb > 0.0 else 0.0


def wale(c: float, delta: float, g: np.ndarray) -> float:
    """WALE (Nicoud & Ducros 1999), the traceless symmetric part of g^2."""
    s = 0.5 * (g + g.T)
    g2 = g @ g
    sd = 0.5 * (g2 + g2.T) - np.trace(g2) / 3.0 * np.eye(3)
    ss, dd = float(np.sum(s * s)), float(np.sum(sd * sd))
    den = ss**2.5 + dd**1.25
    return (c * delta) ** 2 * dd**1.5 / den if den > 0.0 else 0.0


def sigma(c: float, delta: float, g: np.ndarray) -> float:
    """Sigma (Nicoud et al. 2011), from the singular values."""
    s = np.linalg.svd(g, compute_uv=False)
    return (c * delta) ** 2 * s[2] * (s[0] - s[1]) * (s[1] - s[2]) / s[0] ** 2 if s[0] > 0.0 else 0.0


def formula(args: argparse.Namespace) -> int:
    """Compare the Fortran formulas with the independent ones, and check the zeros."""
    rng = np.random.default_rng(119)
    cases = []
    for _ in range(200):
        dl = rng.uniform(0.01, 1.0, 3)
        cases.append((rng.uniform(0.01, 2.0), dl, float(np.prod(dl) ** (1 / 3)), rng.normal(size=(3, 3)) * 10.0 ** rng.uniform(-3, 3)))
    shear = np.zeros((3, 3))
    shear[0, 1] = 3.0
    planar = np.array([[1.0, 2.0, 0.0], [-3.0, -1.0, 0.0], [0.0, 0.0, 0.0]])
    rotation = np.array([[0.0, -2.0, 0.0], [2.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    rank_one = np.outer([1.0, 2.0, -1.0], [3.0, 1.0, 2.0])
    zero = np.zeros((3, 3))
    unit = np.ones(3)
    # name, gradient, which models must vanish (Vreman, WALE, Sigma)
    zeros = [("pure shear", shear, (1, 1, 1)), ("2-D", planar, (0, 0, 1)), ("solid rotation", rotation, (0, 0, 1)),
             ("rank one", rank_one, (1, 0, 1)), ("no gradient", zero, (1, 1, 1))]
    text = "".join(f"{c:.17e} " + " ".join(f"{x:.17e}" for x in [*dl, delta, *g.flatten(order="F")]) + "\n" for c, dl, delta, g in cases)
    text += "".join("1.0 1.0 1.0 1.0 1.0 " + " ".join(f"{x:.17e}" for x in g.flatten(order="F")) + "\n" for _, g, _ in zeros)
    res = subprocess.run([args.exe], input=text, capture_output=True, text=True, check=True)
    got = np.loadtxt(res.stdout.splitlines())
    want = np.array([[vreman(c, dl, g), wale(c, delta, g), sigma(c, delta, g)] for c, dl, delta, g in cases])
    head = got[: len(cases)]
    scale = np.maximum(np.abs(want), 1e-300)
    errs = np.max(np.abs(head - want) / scale, axis=0)
    ok = True
    for name, err, tol in zip(("Vreman", "WALE", "Sigma"), errs, (1e-12, 1e-12, args.sigma_tol)):
        good = err < tol
        ok = ok and good
        print(f"{name:7s} largest relative difference on {len(cases)} random tensors {err:.2e} (bound {tol:.0e})  {'PASS' if good else 'FAIL'}")
    for (name, g, must), row in zip(zeros, got[len(cases) :]):
        scale_g = max(float(np.max(np.abs(g))), 1.0)
        bad = [m for m, flag, v in zip(("Vreman", "WALE", "Sigma"), must, row) if flag and abs(v) > 1e-10 * scale_g]
        ok = ok and not bad
        print(f"{name:15s} nu = {row[0]:.2e} {row[1]:.2e} {row[2]:.2e}: {'FAIL ' + ' '.join(bad) if bad else 'the expected zeros PASS'}")
    return 0 if ok else 1


def history(work: str) -> np.ndarray:
    """Return the history rows."""
    return np.loadtxt(Path(work) / "tgv-turbulence_history.dat", skiprows=1, ndmin=2)


def budget(args: argparse.Namespace) -> int:
    """Check the budget with the subgrid-scale dissipation."""
    left = []
    for arg in args.runs:
        n, work = arg.split(":")
        h = history(work)
        h = h[h[:, 1] <= args.time]
        t, kin, eps, pd, sgs = h[:, 1], h[:, 2], h[:, 6], h[:, 7], h[:, 8]
        drop, phys, model = kin[-1] - kin[0], np.trapezoid(pd - eps, t), np.trapezoid(sgs, t)
        good, without, flipped = (abs(drop - phys + c * model) / abs(drop) for c in (1.0, 0.0, -1.0))
        left.append(good)
        print(f"n {int(n):4d}: integral(eps_sgs) over the K drop {model / abs(drop):.3e}; budget left {good:.3e}, without it {without:.3e}, flipped {flipped:.3e}")
    ok = good < without and good < flipped and model > 0.0
    ratios = [left[i] / left[i + 1] for i in range(len(left) - 1)]
    ok = ok and all(r >= args.ratio_min for r in ratios)
    print(f"ratios per doubling {' '.join(f'{r:.2f}' for r in ratios)} (min {args.ratio_min})  {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def inviscid(args: argparse.Namespace) -> int:
    """Check a run with no viscosity."""
    h = history(args.work)
    eps_max, sgs_min = float(np.max(np.abs(h[:, 6]))), float(np.max(h[:, 8]))
    with open(Path(args.work) / "tgv-conservation_history.dat") as fh:
        names = [w.strip('"') for w in fh.readline().split("=")[1].split()]
    c = np.loadtxt(Path(args.work) / "tgv-conservation_history.dat", skiprows=1, ndmin=2)
    drift = {v: float(np.max(np.abs(c[:, names.index(v)] - c[0, names.index(v)])) / abs(c[0, names.index(v)])) for v in ("int_r", "int_rE")}
    ok = eps_max == 0.0 and sgs_min > 0.0 and max(drift.values()) < args.tol
    print(f"eps {eps_max:.2e}, largest eps_sgs {sgs_min:.3e}, drift of mass {drift['int_r']:.2e}, of energy {drift['int_rE']:.2e} (bound {args.tol:.0e})  {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def main() -> int:
    """Run a check."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("formula")
    p.add_argument("exe")
    p.add_argument("--sigma-tol", type=float, default=1e-6, help="the closed form loses digits on nearly equal singular values")
    p.set_defaults(func=formula)
    p = sub.add_parser("budget")
    p.add_argument("runs", nargs="+")
    p.add_argument("--time", type=float, default=2.0)
    p.add_argument("--ratio-min", type=float, default=4.0)
    p.set_defaults(func=budget)
    p = sub.add_parser("inviscid")
    p.add_argument("work")
    p.add_argument("--tol", type=float, default=1e-12)
    p.set_defaults(func=inviscid)
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
