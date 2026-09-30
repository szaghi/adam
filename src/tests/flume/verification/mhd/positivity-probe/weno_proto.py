#!/usr/bin/env python3
"""WENO5 + SSP-RK3 on the Balsara-Spicer blast with three positivity limiters (issue #41, M2 stretch S1).

Scheme: component-wise WENO5-JS (Jiang & Shu 1996) of the Lax-Friedrichs split fluxes f+- = (f +- alpha q) / 2, alpha
the maximum of |u_n| + c_f (and c_h) over the 6-cell stencil of each face (FLUME's conservative-variable mode); FLUME's
GLM model (lf_proto.py); SSP-RK3 in Shu-Osher form, every stage a forward-Euler step of dt; periodic 2-D grid.
Limiters, applied to the face fluxes of every stage:
  none      no limiting;
  onesided  Hu, Adams & Shu (2013): per face, theta keeps the one-sided partial states q_L - 2 D lambda (F - f(q_L)) and
            q_R + 2 D lambda (F - f(q_R)) admissible (the current Fortran design);
  cell      Xu (2014) parametrised: per cell, Lambda = the largest value with every corner of [0, Lambda]^(2 D) of the
            face thetas keeping the cell update admissible (needs only the first-order update admissible);
            theta_face = min(Lambda_L, Lambda_R).
Reports the first stage with a non-positive density or pressure, or the final time, and the limiter counters.
--powell adds the Powell source -(div B)(0, B, u.B, u, 0) to both updates; --alpha-factor scales the first-order
Lax-Friedrichs speed.
Usage: weno_proto.py <limiter> [--cells N] [--b0 B] [--p-in P] [--cfl C] [--t-end T] [--powell] [--alpha-factor F]
"""

from __future__ import annotations

import argparse
import itertools
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from lf_proto import GAMMA, flux, pressure  # noqa: E402

MODEL, CH, EPS = "glm", 60.0, 1.0e-13
ALPHA_FACTOR = [1.0]  # set from --alpha-factor
SHIFTS = (2, 1, 0, -1, -2, -3)  # np.roll shifts giving cells i-2 .. i+3 at index i


def admissible(u: np.ndarray) -> np.ndarray:
    """Cells with density and pressure at least EPS."""
    return (u[0] >= EPS) & (pressure(u, MODEL) >= EPS)


def weno5(vm2: np.ndarray, vm1: np.ndarray, v0: np.ndarray, vp1: np.ndarray, vp2: np.ndarray) -> np.ndarray:
    """WENO5-JS value at the face between v0 and vp1 (upwind side v0)."""
    b0 = 13 / 12 * (vm2 - 2 * vm1 + v0) ** 2 + 0.25 * (vm2 - 4 * vm1 + 3 * v0) ** 2
    b1 = 13 / 12 * (vm1 - 2 * v0 + vp1) ** 2 + 0.25 * (vm1 - vp1) ** 2
    b2 = 13 / 12 * (v0 - 2 * vp1 + vp2) ** 2 + 0.25 * (3 * v0 - 4 * vp1 + vp2) ** 2
    a0, a1, a2 = 0.1 / (1e-6 + b0) ** 2, 0.6 / (1e-6 + b1) ** 2, 0.3 / (1e-6 + b2) ** 2
    s = a0 + a1 + a2
    q0, q1, q2 = 2 * vm2 - 7 * vm1 + 11 * v0, -vm1 + 5 * v0 + 2 * vp1, 2 * v0 + 5 * vp1 - vp2
    return (a0 * q0 + a1 * q1 + a2 * q2) / (6 * s)


def fluxes(u: np.ndarray, d: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """High-order and Lax-Friedrichs face fluxes at i+1/2 along d, and the physical flux of the cells."""
    f, speed = flux(u, d, MODEL, CH)
    speed = np.maximum(speed, CH)
    alpha = np.max([np.roll(speed, s, axis=d) for s in SHIFTS], axis=0)
    fh = np.empty_like(u)
    for m in range(u.shape[0]):
        fs = {s: np.roll(f[m], s, axis=d) for s in SHIFTS}
        us = {s: np.roll(u[m], s, axis=d) for s in SHIFTS}
        fp = {s: 0.5 * (fs[s] + alpha * us[s]) for s in SHIFTS}  # right-going, left-biased: cells i-2 .. i+2
        fm = {s: 0.5 * (fs[s] - alpha * us[s]) for s in SHIFTS}  # left-going, right-biased: cells i+3 .. i-1
        fh[m] = weno5(fp[2], fp[1], fp[0], fp[-1], fp[-2]) + weno5(fm[-3], fm[-2], fm[-1], fm[0], fm[1])
    ur = np.roll(u, -1, axis=1 + d)
    fr = np.roll(f, -1, axis=1 + d)
    a_lf = ALPHA_FACTOR[0] * np.maximum(speed, np.roll(speed, -1, axis=d))
    flf = 0.5 * (f + fr) - 0.5 * a_lf * (ur - u)
    return fh, flf, f


def theta_max(v0: np.ndarray, dv: np.ndarray, ok0: np.ndarray, iters: int = 30) -> np.ndarray:
    """Per cell, the largest theta in [0, 1] with v0 + theta dv admissible (0 where v0 is not)."""
    full = admissible(v0 + dv)
    need = ok0 & ~full
    lo, hi = np.zeros(full.shape), np.ones(full.shape)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        ok = admissible(v0 + mid * dv)
        lo = np.where(need & ok, mid, lo)
        hi = np.where(need & ~ok, mid, hi)
    return np.where(ok0, np.where(full, 1.0, lo), 0.0)


def powell(u: np.ndarray, h: float) -> np.ndarray:
    """Powell (Godunov) source -(div B)(0, B, u.B, u, 0), div B by central differences (Wu & Shu 2019 form)."""
    vel = u[1:4] / u[0]
    divb = sum((np.roll(u[5 + d], -1, axis=d) - np.roll(u[5 + d], 1, axis=d)) / (2 * h) for d in (0, 1))
    s = np.zeros_like(u)
    s[1:4] = -divb * u[5:8]
    s[4] = -divb * (vel * u[5:8]).sum(axis=0)
    s[5:8] = -divb * vel
    return s


def limit(u: np.ndarray, fh: list, flf: list, f: list, lam: float, kind: str,
          src: np.ndarray) -> tuple[list, int, int]:
    """Limited face fluxes of the two directions (src: dt times the source, part of the first-order update)."""
    if kind == "none":
        return fh, 0, 0
    if kind == "onesided":
        out, lim, unres = [], 0, 0
        step = 2 * 2 * lam
        for d in (0, 1):
            df = fh[d] - flf[d]
            ur, fr = np.roll(u, -1, axis=1 + d), np.roll(f[d], -1, axis=1 + d)
            vl0, vr0 = u - step * (flf[d] - f[d]), ur + step * (flf[d] - fr)
            okl, okr = admissible(vl0), admissible(vr0)
            th = np.minimum(theta_max(vl0, -step * df, okl), theta_max(vr0, step * df, okr))
            unres += int(np.sum(~(okl & okr)))
            lim += int(np.sum(th < 1.0))
            out.append(flf[d] + th * df)
        return out, lim, unres
    ul = u - lam * sum(flf[d] - np.roll(flf[d], 1, axis=1 + d) for d in (0, 1)) + src
    contrib = []
    for d in (0, 1):
        df = fh[d] - flf[d]
        contrib += [-lam * df, lam * np.roll(df, 1, axis=1 + d)]  # face i+1/2, face i-1/2 of each cell
    ok0 = admissible(ul)
    lam_cell = np.where(ok0, 1.0, 0.0)
    for r in range(1, 5):
        for sub in itertools.combinations(range(4), r):
            lam_cell = np.minimum(lam_cell, theta_max(ul, sum(contrib[k] for k in sub), ok0))
    out, lim = [], 0
    for d in (0, 1):
        th = np.minimum(lam_cell, np.roll(lam_cell, -1, axis=d))
        lim += int(np.sum(th < 1.0))
        out.append(flf[d] + th * (fh[d] - flf[d]))
    return out, lim, int(np.sum(~ok0))


def main() -> int:
    """Run the blast with one limiter."""
    ap = argparse.ArgumentParser()
    ap.add_argument("limiter", choices=("none", "onesided", "cell"))
    ap.add_argument("--cells", type=int, default=128)
    ap.add_argument("--b0", type=float, default=100.0 / np.sqrt(4.0 * np.pi))
    ap.add_argument("--p-in", type=float, default=1000.0)
    ap.add_argument("--cfl", type=float, default=0.25)
    ap.add_argument("--t-end", type=float, default=0.01)
    ap.add_argument("--powell", action="store_true", help="add the Powell source term")
    ap.add_argument("--alpha-factor", type=float, default=1.0, help="scale of the Lax-Friedrichs speed")
    a = ap.parse_args()
    ALPHA_FACTOR[0] = a.alpha_factor
    n = a.cells
    h = 1.0 / n
    x = (np.arange(n) + 0.5) * h
    xx, yy = np.meshgrid(x, x, indexing="ij")
    u = np.zeros((9, n, n))
    u[0] = 1.0
    u[5] = u[6] = a.b0 / np.sqrt(2.0)
    p0 = np.where((xx - 0.5) ** 2 + (yy - 0.5) ** 2 < 0.01, a.p_in, 0.1)
    u[4] = p0 / (GAMMA - 1.0) + 0.5 * (u[5] ** 2 + u[6] ** 2)
    t, it, counters = 0.0, 0, {"lim": 0, "unres": 0}

    def fe(v: np.ndarray, dt: float) -> np.ndarray:
        res = [fluxes(v, d) for d in (0, 1)]
        src = dt * powell(v, h) if a.powell else np.zeros_like(v)
        fl, lim, unres = limit(v, [r[0] for r in res], [r[1] for r in res], [r[2] for r in res], dt / h, a.limiter,
                               src)
        counters["lim"] += lim
        counters["unres"] = max(counters["unres"], unres)
        return v - dt / h * sum(fl[d] - np.roll(fl[d], 1, axis=1 + d) for d in (0, 1)) + src

    def report(what: str) -> str:
        name = a.limiter + ("+powell" if a.powell else "")
        return f"{name:15s}: {what}; limited faces {counters['lim']}, max unresolved per stage {counters['unres']}"

    while t < a.t_end - 1e-15:
        speeds = [np.maximum(flux(u, d, MODEL, CH)[1], CH) for d in (0, 1)]
        dt = min(a.cfl / np.max(speeds[0] / h + speeds[1] / h), a.t_end - t)
        it += 1
        stage = fe(u, dt)
        for s, (cu, cs) in enumerate(((0.75, 0.25), (1.0 / 3.0, 2.0 / 3.0), (0.0, 0.0)), start=1):
            if not admissible(stage).all():
                print(report(f"INADMISSIBLE at step {it} stage {s}, min p {pressure(stage, MODEL).min():.3e}"))
                return 1
            if s == 3:
                break
            stage = cu * u + cs * fe(stage, dt)
        u = stage
        t += dt
    print(report(f"reached t {t:.3e} in {it} steps, min rho {u[0].min():.3e}, min p {pressure(u, MODEL).min():.3e}"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
