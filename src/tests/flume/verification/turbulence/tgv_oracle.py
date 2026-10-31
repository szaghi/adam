#!/usr/bin/env python3
"""Oracles of the turbulence history (issue #119, T0.1, T0.2), multirealm-free, numpy only.

t0      the first row of each run (argument n:work) against the exact Taylor-Green integrals: K, whose trigonometric
        integrand the midpoint rule integrates exactly, within --kin-tol; the gradient integrals Z and eps (4th-order
        central gradients) with observed order at least --order-min; eps_d and pd, zero by symmetry, within --zero-tol.
agree   two turbulence histories (two rank counts, two backends) within --tol relative to each column's largest value.
budget  dK/dt = pd - eps: the integral of `-(dK/dt) + pd - eps`, the scheme's numerical dissipation, over the
        K drop, from t = 0 to --time, for a resolution ladder n:work; it must fall by at least --ratio-min per doubling. Its size
        is a measurement, not a contract: at Re 1600, M 0.1, t = 2, WENO-5, 70 % of the drop at 32^3, 8.8 % at 64^3.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

COLUMNS = ("it", "time", "K", "Z", "eps_s", "eps_d", "eps", "pd")


def history(work: str) -> np.ndarray:
    """Return the rows of the turbulence history of a work directory."""
    return np.loadtxt(Path(work) / "tgv-turbulence_history.dat", skiprows=1, ndmin=2)


def exact(n: int, mach: float, reynolds: float) -> dict[str, float]:
    """Return the exact integrals at t = 0, by the midpoint rule on the run's own cell centres."""
    gamma = 1.4
    h = 2.0 * np.pi / n
    x = (np.arange(n) + 0.5) * h
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    p0 = 1.0 / (gamma * mach**2)
    p = p0 + (np.cos(2 * X) + np.cos(2 * Y)) * (np.cos(2 * Z) + 2.0) / 16.0
    rho = p / p0
    u = np.sin(X) * np.cos(Y) * np.cos(Z)
    v = -np.cos(X) * np.sin(Y) * np.cos(Z)
    wx = -np.cos(X) * np.sin(Y) * np.sin(Z)  # -dv/dz
    wy = -np.sin(X) * np.cos(Y) * np.sin(Z)  # du/dz
    wz = 2.0 * np.sin(X) * np.sin(Y) * np.cos(Z)  # dv/dx - du/dy
    w2 = wx**2 + wy**2 + wz**2
    vol = h**3
    mu = 1.0 / reynolds if reynolds > 0.0 else 0.0
    return {"K": float(np.sum(0.5 * rho * (u**2 + v**2)) * vol), "Z": float(np.sum(0.5 * rho * w2) * vol),
            "eps": float(mu * np.sum(w2) * vol)}


def order(ns: list[int], errs: list[float]) -> list[float]:
    """Return the observed orders of a ladder."""
    return [float(np.log2(errs[i] / errs[i + 1]) / np.log2(ns[i + 1] / ns[i])) for i in range(len(errs) - 1)]


def t0(args: argparse.Namespace) -> int:
    """Check the first rows against the exact integrals."""
    ok = True
    ladder = [(int(a.split(":")[0]), a.split(":")[1]) for a in args.runs]
    errs: dict[str, list[float]] = {"Z": [], "eps": []}
    for n, work in ladder:
        row = dict(zip(COLUMNS, history(work)[0]))
        ref = exact(n, args.mach, args.reynolds)
        kin = abs(row["K"] - ref["K"]) / ref["K"]
        zero = max(abs(row["eps_d"]), abs(row["pd"])) / ref["K"]
        errs["Z"].append(abs(row["Z"] - ref["Z"]) / ref["Z"])
        errs["eps"].append(abs(row["eps"] - ref["eps"]) / ref["eps"])
        split = abs(row["eps"] - row["eps_s"] - row["eps_d"]) / ref["eps"]
        good = kin <= args.kin_tol and zero <= args.zero_tol and split <= args.zero_tol
        ok = ok and good
        print(f"n {n:4d}: K rel {kin:.2e}, eps_d and pd over K {zero:.2e}, Z rel {errs['Z'][-1]:.3e}, "
              f"eps rel {errs['eps'][-1]:.3e}, eps - eps_s - eps_d {split:.2e}  {'PASS' if good else 'FAIL'}")
    ns = [n for n, _ in ladder]
    for name, err in errs.items():
        orders = order(ns, err)
        good = min(orders) >= args.order_min
        ok = ok and good
        print(f"{name}: observed orders {' '.join(f'{o:.2f}' for o in orders)} (min {args.order_min})  "
              f"{'PASS' if good else 'FAIL'}")
    return 0 if ok else 1


def agree(args: argparse.Namespace) -> int:
    """Compare two histories."""
    a, b = history(args.work_a), history(args.work_b)
    if a.shape != b.shape:
        print(f"agree: shapes differ {a.shape} {b.shape}  FAIL")
        return 1
    ok = True
    for c, name in enumerate(COLUMNS[1:], start=1):
        scale = max(np.max(np.abs(a[:, c])), np.max(np.abs(a[:, 2])) if name in ("eps_d", "pd") else 0.0, 1e-300)
        diff = float(np.max(np.abs(a[:, c] - b[:, c])) / scale)
        good = diff <= args.tol
        ok = ok and good
        print(f"{name:6s} largest difference over its scale {diff:.2e} (bound {args.tol:.0e})  {'PASS' if good else 'FAIL'}")
    return 0 if ok else 1


def budget(args: argparse.Namespace) -> int:
    """Measure the numerical dissipation of a ladder, as a share of the K drop."""
    res = []
    for arg in args.runs:
        n, work = arg.split(":")
        h = history(work)
        h = h[h[:, 1] <= args.time]
        t, kin, eps, pd = h[:, 1], h[:, 2], h[:, 6], h[:, 7]
        work_phys = np.trapezoid(pd - eps, t)
        res.append(abs(kin[-1] - kin[0] - work_phys) / abs(kin[-1] - kin[0]))
        print(f"n {int(n):4d}: t {t[-1]:.3f}, {len(t)} rows, K drop {kin[0] - kin[-1]:.5e}, |K drop - integral(eps - pd)| "
              f"over the drop {res[-1]:.3e}")
    ratios = [res[i] / res[i + 1] for i in range(len(res) - 1)]
    ok = bool(ratios) and min(ratios) >= args.ratio_min
    print(f"ratios per doubling {' '.join(f'{r:.2f}' for r in ratios)} (min {args.ratio_min}), finest share {res[-1]:.3e}  "
          f"{'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def pressure_work(args: argparse.Namespace) -> int:
    """Check pd where it weighs: the budget closes better with it than without it or with its sign flipped."""
    h = history(args.work)
    h = h[h[:, 1] <= args.time]
    t, kin, eps, pd = h[:, 1], h[:, 2], h[:, 6], h[:, 7]
    drop, viscous, work = kin[-1] - kin[0], np.trapezoid(-eps, t), np.trapezoid(pd, t)
    good, dropped, flipped = (abs(drop - viscous - c * work) / abs(drop) for c in (1.0, 0.0, -1.0))
    ok = good < dropped and good < flipped
    print(f"integral(pd) over the K drop {work / abs(drop):.3e}; budget left {good:.3e}, without pd {dropped:.3e}, "
          f"pd flipped {flipped:.3e}  {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def main() -> int:
    """Dispatch the oracles."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("t0")
    p.add_argument("runs", nargs="+", help="n:work")
    p.add_argument("--mach", type=float, default=0.1)
    p.add_argument("--reynolds", type=float, default=1600.0)
    p.add_argument("--kin-tol", type=float, default=1e-12)
    p.add_argument("--zero-tol", type=float, default=1e-9)
    p.add_argument("--order-min", type=float, default=3.5)
    p.set_defaults(func=t0)
    p = sub.add_parser("agree")
    p.add_argument("work_a")
    p.add_argument("work_b")
    p.add_argument("--tol", type=float, default=1e-11)
    p.set_defaults(func=agree)
    p = sub.add_parser("pd")
    p.add_argument("work")
    p.add_argument("--time", type=float, default=2.0)
    p.set_defaults(func=pressure_work)

    p = sub.add_parser("budget")
    p.add_argument("runs", nargs="+", help="n:work, coarse to fine")
    p.add_argument("--time", type=float, default=2.0)
    p.add_argument("--ratio-min", type=float, default=4.0)
    p.set_defaults(func=budget)
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
