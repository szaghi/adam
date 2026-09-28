#!/usr/bin/env python3
"""Oracle of the magnetised vortex (issue #41, M2-P5, MV-7; Balsara 2004, ApJS 151).

Why: the vortex is an exact steady equilibrium convected by the free stream, so at the final time T (the last row of
the div(B) history) the exact solution is the initial vortex translated by v0 T, periodically wrapped: every interior
cell of the last checkpoint is compared with it pointwise. Errors: L1 of Stone et al. (eps = sqrt(sum_v mean|dq_v|^2))
and Linf (max over the cells and the 8 conservative variables); observed orders p = log2(e_N / e_2N), asserted on the
finest pair (--order-min for L1, --linf-order-min for Linf: Linf reaches the asymptotic range later). The exact
solution is evaluated here independently of the Fortran IC (same formulas, re-derived from the radial balance
dp/dr = rho v^2/r - B^2/r - d(B^2/2)/dr), so a wrong equilibrium in either shows up as a non-converging error. The
initial checkpoint must match it to round-off (--ic-tol). The final max|div B| of the history, the discrete divergence
of a div-free analytic field (second-order operator, fdv_order = 2), must converge too (--divb-order-min).

Usage:
    mhd_vortex_oracle.py <work> [<work> ...] [--order-min P] [--linf-order-min P] [--divb-order-min P]
                         [--l1-max E [E ...]] [--linf-max E [E ...]] [--ic-tol T] [--ngc N]
"""

from __future__ import annotations

import argparse
import configparser
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "linear-wave"))
from linear_wave_oracle import profile  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "divb-peak"))
from divb_peak_oracle import history  # noqa: E402


def read_ini(work: Path) -> configparser.ConfigParser:
    """Return the input of a run (its only .ini file)."""
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.read(next(work.glob("*.ini")))
    return ini


def exact(ini: configparser.ConfigParser, xy: np.ndarray, t: float) -> np.ndarray:
    """Return the exact conservative state [8, cells] at the cell centres xy [cells, 2] and time t."""
    ic, fs, grid = ini["initial_conditions"], ini["initial_conditions_region_1"], ini["grid"]
    gamma = float(ini["physics"]["cp"]) / float(ini["physics"]["cv"])
    x0, y0, a, kappa, mu = (float(ic[k]) for k in ("x0", "y0", "radius", "kappa", "mu"))
    r0, u0, v0, w0, p0, bz0 = (float(fs[k]) for k in ("r", "u", "v", "w", "p", "bz"))
    lo = np.array([float(grid["emin_x"]), float(grid["emin_y"])])
    size = np.array([float(grid["emax_x"]), float(grid["emax_y"])]) - lo
    d = xy - np.array([x0 + u0 * t, y0 + v0 * t])
    d = (d + 0.5 * size) % size - 0.5 * size  # nearest periodic image
    dx, dy = d[:, 0] / a, d[:, 1] / a
    rr = dx * dx + dy * dy
    e = np.exp(0.5 * (1.0 - rr))
    u = u0 - kappa / (2.0 * math.pi) * dy * e
    v = v0 + kappa / (2.0 * math.pi) * dx * e
    bx = -mu / (2.0 * math.pi) * dy * e
    by = mu / (2.0 * math.pi) * dx * e
    p = p0 + (mu**2 * (1.0 - rr) - r0 * kappa**2) / (8.0 * math.pi**2) * e * e
    r = np.full_like(u, r0)
    w = np.full_like(u, w0)
    bz = np.full_like(u, bz0)
    energy = p / (gamma - 1.0) + 0.5 * r * (u * u + v * v + w * w) + 0.5 * (bx * bx + by * by + bz * bz)
    return np.stack([r, r * u, r * v, r * w, energy, bx, by, bz])


def errors(work: Path, step: str, t: float, ngc: int) -> tuple[float, float]:
    """Return (L1 eps, Linf) of the first or last checkpoint of a run against the exact solution at t."""
    cells = profile(work, step, ngc)
    keys = sorted(cells)
    q = np.array([cells[k] for k in keys]).T
    xy = np.array([[k[0], k[1]] for k in keys])
    diff = np.abs(q - exact(read_ini(work), xy, t))
    return float(math.sqrt(np.sum(np.mean(diff, axis=1) ** 2))), float(diff.max())


def bound(values: list[float] | None, n: int) -> float | None:
    """Return the n-th bound of a per-run list (one value applies to every run)."""
    if values is None:
        return None
    return values[0] if len(values) == 1 else values[n]


def main() -> int:
    """Run the MV-7 checks, print a report, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work", type=Path, nargs="+", help="runs, coarse first, N doubling")
    parser.add_argument("--order-min", type=float, default=None, help="minimum L1 order of the finest pair")
    parser.add_argument("--linf-order-min", type=float, default=None, help="minimum Linf order of the finest pair")
    parser.add_argument("--divb-order-min", type=float, default=None, help="minimum order of the final max|div B|")
    parser.add_argument("--l1-max", type=float, nargs="+", default=None, help="bound of L1 (one, or one per run)")
    parser.add_argument("--linf-max", type=float, nargs="+", default=None, help="bound of Linf (one, or one per run)")
    parser.add_argument("--ic-tol", type=float, default=None, help="bound of Linf of the initial checkpoint")
    parser.add_argument("--ngc", type=int, default=3)
    args = parser.parse_args()
    for name in ("l1_max", "linf_max"):
        values = getattr(args, name)
        if values is not None and len(values) not in (1, len(args.work)):
            sys.exit(f"mhd_vortex_oracle: --{name.replace('_', '-')} takes one bound or one per run")
    ok = True
    l1, linf, divb = [], [], []
    for n, w in enumerate(args.work):
        hist = history(w, "divb")
        t_end = float(hist["time"][-1])
        divb.append(float(hist["max_divb"][-1]))
        _, ic_inf = errors(w, "first", 0.0, args.ngc)
        e1, einf = errors(w, "last", t_end, args.ngc)
        l1.append(e1)
        linf.append(einf)
        line = f"{w.name}: t {t_end:.6g}  L1 {e1:.4e}  Linf {einf:.4e}  max|div B| {divb[-1]:.3e}  IC Linf {ic_inf:.1e}"
        if n > 0:
            line += (f"  orders L1 {math.log2(l1[-2] / e1):+.2f} Linf {math.log2(linf[-2] / einf):+.2f}"
                     f" divB {math.log2(divb[-2] / divb[-1]):+.2f}")
        print(line)
        for label, value, b in (("L1", e1, bound(args.l1_max, n)), ("Linf", einf, bound(args.linf_max, n)),
                                ("IC Linf", ic_inf, args.ic_tol)):
            if b is not None:
                good = value <= b
                ok &= good
                print(f"   {label} {value:.4e} {'<=' if good else '>'} {b:.3e}: {'PASS' if good else 'FAIL'}")
    if len(args.work) > 1:
        for label, values, minimum in (("L1", l1, args.order_min), ("Linf", linf, args.linf_order_min),
                                       ("max|div B|", divb, args.divb_order_min)):
            if minimum is None:
                continue
            p = math.log2(values[-2] / values[-1])
            good = p >= minimum
            ok &= good
            verdict = "PASS" if good else "FAIL"
            print(f"   {label} finest-pair order {p:+.2f} {'>=' if good else '<'} {minimum}: {verdict}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
