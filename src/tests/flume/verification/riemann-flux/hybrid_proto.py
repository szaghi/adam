#!/usr/bin/env python3
"""Hybrid face-correction scheme prototype (issue #47, M3-P0): 1-D Euler, WENO5 interpolation, Riemann flux, correction.

Face flux at i+1/2 (Chen, Toth & Gombosi 2016 form, coefficients derived with sympy):
  6th: F^ = 64/45 F - 13/60 (f_i + f_i+1) + 1/180 (f_i-1 + f_i+2)
  4th: F^ = 4/3 F - 1/6 (f_i + f_i+1)
  none: F^ = F
with F = RS(q_L, q_R) and q_L, q_R the WENO5 interpolation (linear combination (3, -20, 90, 60, -5)/128) of point
values, in characteristic (Roe eigenvectors at the face) or primitive variables. The correction C = F^ - F is
multiplied by a per-face sensor s in [0, 1]:
  on      s = 1;
  weno    s = 1 where every characteristic field's nonlinear weights stay within a factor of the linear weights
          (min_k w_k / d_k >= TAU), else 0;
  jump    s = 1 where the pressure and density jumps across the 6-cell stencil are below TAU_J (relative), else 0.
Baseline 'split': FLUME's current scheme (characteristic WENO5 reconstruction of the local Lax-Friedrichs split fluxes).
SSP-RK3, CFL 0.4 on max(|u| + a).

Tests: 'order' (density wave, periodic, N = 20..320, t = 1, dt ~ h^(5/3) so that the SSP-RK3 error stays below
the space error), 'sod', 'lax', 'shu-osher' (L1 of rho against an 8x finer 'split' reference averaged to the grid and
cached in the temporary directory; overshoot, min p).
Usage: hybrid_proto.py <test> [--scheme split|hybrid] [--rs llf|hll|hllc] [--corr 6|4|0] [--sensor on|weno|jump]
[--vars char|prim] [--cells N]
"""

from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path

import numpy as np

GAMMA, G = 1.4, 4  # ratio of specific heats, ghost cells
TAU, TAU_J = 0.2, 0.1


def prim(u: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Density, velocity, pressure."""
    r = u[0]
    v = u[1] / r
    return r, v, (GAMMA - 1.0) * (u[2] - 0.5 * r * v * v)


def cons(r: np.ndarray, v: np.ndarray, p: np.ndarray) -> np.ndarray:
    """Conservative state."""
    return np.array([r, r * v, p / (GAMMA - 1.0) + 0.5 * r * v * v])


def flux(u: np.ndarray) -> np.ndarray:
    """Physical flux."""
    r, v, p = prim(u)
    return np.array([r * v, r * v * v + p, v * (u[2] + p)])


def roe(ul: np.ndarray, ur: np.ndarray) -> tuple:
    """Roe-averaged velocity, enthalpy, sound speed and the right/left eigenvector matrices (3, 3, n)."""
    rl, vl, pl = prim(ul)
    rr, vr, pr = prim(ur)
    sl, sr = np.sqrt(rl), np.sqrt(rr)
    v = (sl * vl + sr * vr) / (sl + sr)
    h = (sl * (ul[2] + pl) / rl + sr * (ur[2] + pr) / rr) / (sl + sr)
    a = np.sqrt((GAMMA - 1.0) * (h - 0.5 * v * v))
    one = np.ones_like(v)
    rm = np.array([[one, one, one], [v - a, v, v + a], [h - v * a, 0.5 * v * v, h + v * a]])
    b1 = (GAMMA - 1.0) / (a * a)
    b2 = 0.5 * v * v * b1
    lm = np.array([[0.5 * (b2 + v / a), -0.5 * (b1 * v + 1.0 / a), 0.5 * b1],
                   [1.0 - b2, b1 * v, -b1],
                   [0.5 * (b2 - v / a), -0.5 * (b1 * v - 1.0 / a), 0.5 * b1]])
    return v, h, a, (rm, lm)


def smoothness(vm2, vm1, v0, vp1, vp2) -> tuple:
    """Jiang-Shu smoothness indicators of the three candidate stencils."""
    b0 = 13 / 12 * (vm2 - 2 * vm1 + v0) ** 2 + 0.25 * (vm2 - 4 * vm1 + 3 * v0) ** 2
    b1 = 13 / 12 * (vm1 - 2 * v0 + vp1) ** 2 + 0.25 * (vm1 - vp1) ** 2
    b2 = 13 / 12 * (v0 - 2 * vp1 + vp2) ** 2 + 0.25 * (3 * v0 - 4 * vp1 + vp2) ** 2
    return b0, b1, b2


def weno5_interp(vm2, vm1, v0, vp1, vp2) -> tuple[np.ndarray, np.ndarray]:
    """WENO5 interpolation at the face between v0 and vp1 (upwind side v0), and min_k w_k / d_k."""
    d = (1 / 16, 10 / 16, 5 / 16)
    a = [d[k] / (1e-6 + b) ** 2 for k, b in enumerate(smoothness(vm2, vm1, v0, vp1, vp2))]
    s = a[0] + a[1] + a[2]
    w = [ak / s for ak in a]
    q0 = (3 * vm2 - 10 * vm1 + 15 * v0) / 8
    q1 = (-vm1 + 6 * v0 + 3 * vp1) / 8
    q2 = (3 * v0 + 6 * vp1 - vp2) / 8
    return w[0] * q0 + w[1] * q1 + w[2] * q2, np.minimum.reduce([w[k] / d[k] for k in range(3)])


def weno5_recon(vm2, vm1, v0, vp1, vp2) -> np.ndarray:
    """WENO5-JS reconstruction (cell averages -> face value), upwind side v0."""
    b0, b1, b2 = smoothness(vm2, vm1, v0, vp1, vp2)
    a0, a1, a2 = 0.1 / (1e-6 + b0) ** 2, 0.6 / (1e-6 + b1) ** 2, 0.3 / (1e-6 + b2) ** 2
    q0, q1, q2 = 2 * vm2 - 7 * vm1 + 11 * v0, -vm1 + 5 * v0 + 2 * vp1, 2 * v0 + 5 * vp1 - vp2
    return (a0 * q0 + a1 * q1 + a2 * q2) / (6 * (a0 + a1 + a2))


def riemann(ul: np.ndarray, ur: np.ndarray, kind: str) -> np.ndarray:
    """LLF, HLL (Einfeldt speeds) or HLLC (Batten speeds) flux."""
    rl, vl, pl = prim(ul)
    rr, vr, pr = prim(ur)
    al, ar = np.sqrt(GAMMA * pl / rl), np.sqrt(GAMMA * pr / rr)
    fl, fr = flux(ul), flux(ur)
    if kind == "llf":
        s = np.maximum(np.abs(vl) + al, np.abs(vr) + ar)
        return 0.5 * (fl + fr) - 0.5 * s * (ur - ul)
    v, _, a, _ = roe(ul, ur)
    sl = np.minimum(vl - al, v - a)
    sr = np.maximum(vr + ar, v + a)
    if kind == "hll":
        fh = (sr * fl - sl * fr + sl * sr * (ur - ul)) / (sr - sl)
        return np.where(sl >= 0, fl, np.where(sr <= 0, fr, fh))
    sm = (pr - pl + rl * vl * (sl - vl) - rr * vr * (sr - vr)) / (rl * (sl - vl) - rr * (sr - vr))

    def star(u, r, vv, p, s):
        c = r * (s - vv) / (s - sm)
        return np.array([c, c * sm, c * (u[2] / r + (sm - vv) * (sm + p / (r * (s - vv))))])

    usl, usr = star(ul, rl, vl, pl, sl), star(ur, rr, vr, pr, sr)
    return np.where(sl >= 0, fl, np.where(sm >= 0, fl + sl * (usl - ul), np.where(sr > 0, fr + sr * (usr - ur), fr)))


def rhs(u: np.ndarray, h: float, cfg: dict, bc: str) -> np.ndarray:
    """-(F^_{i+1/2} - F^_{i-1/2}) / h on the interior cells."""
    up = np.pad(u, ((0, 0), (G, G)), mode="wrap" if bc == "periodic" else "edge")
    n = u.shape[1]
    f = flux(up)
    idx = np.arange(G - 1, G + n)  # faces i+1/2 for i = G-1 .. G+n-1

    def at(k):
        return up[:, idx + k]

    if cfg["scheme"] == "split":
        _, _, _, (rm, lm) = roe(at(0), at(1))
        speeds = []
        for k in range(-2, 4):
            r, v, p = prim(at(k))
            speeds.append(np.abs(v) + np.sqrt(GAMMA * p / r))
        alpha = np.max(speeds, axis=0)
        fp = {k: np.einsum("ijn,jn->in", lm, 0.5 * (f[:, idx + k] + alpha * at(k))) for k in range(-2, 4)}
        fm = {k: np.einsum("ijn,jn->in", lm, 0.5 * (f[:, idx + k] - alpha * at(k))) for k in range(-2, 4)}
        wc = weno5_recon(fp[-2], fp[-1], fp[0], fp[1], fp[2]) + weno5_recon(fm[3], fm[2], fm[1], fm[0], fm[-1])
        fh = np.einsum("ijn,jn->in", rm, wc)
    else:
        if cfg["vars"] == "char":
            _, _, _, (rm, lm) = roe(at(0), at(1))
            w = {k: np.einsum("ijn,jn->in", lm, at(k)) for k in range(-2, 4)}
        else:
            w = {k: np.array(prim(at(k))) for k in range(-2, 4)}
        wl, sml = weno5_interp(w[-2], w[-1], w[0], w[1], w[2])
        wr, smr = weno5_interp(w[3], w[2], w[1], w[0], w[-1])
        if cfg["vars"] == "char":
            ql, qr = np.einsum("ijn,jn->in", rm, wl), np.einsum("ijn,jn->in", rm, wr)
        else:
            ql, qr = cons(*wl), cons(*wr)
        fr = riemann(ql, qr, cfg["rs"])
        c = np.zeros_like(fr)
        if cfg["corr"] == 6:
            c = 19 / 45 * fr - 13 / 60 * (f[:, idx] + f[:, idx + 1]) + 1 / 180 * (f[:, idx - 1] + f[:, idx + 2])
        elif cfg["corr"] == 4:
            c = 1 / 3 * fr - 1 / 6 * (f[:, idx] + f[:, idx + 1])
        if cfg["sensor"] == "weno":
            s = (np.minimum(sml.min(axis=0), smr.min(axis=0)) >= TAU).astype(float)
        elif cfg["sensor"] == "jump":
            r6 = np.array([prim(at(k))[0] for k in range(-2, 4)])
            p6 = np.array([prim(at(k))[2] for k in range(-2, 4)])
            jr = (r6.max(axis=0) - r6.min(axis=0)) / r6.min(axis=0)
            jp = (p6.max(axis=0) - p6.min(axis=0)) / p6.min(axis=0)
            s = ((jr < TAU_J) & (jp < TAU_J)).astype(float)
        else:
            s = np.ones(fr.shape[1])
        cfg["corr_off"] += int(np.sum(s < 1.0))
        fh = fr + s * c
    return -(fh[:, 1:] - fh[:, :-1]) / h


def run(u: np.ndarray, h: float, t_end: float, cfg: dict, bc: str, cfl: float = 0.4) -> np.ndarray | None:
    """SSP-RK3 to t_end; None on a non-positive density or pressure."""
    t = 0.0
    while t < t_end - 1e-14:
        r, v, p = prim(u)
        dt = min(cfl * h / np.max(np.abs(v) + np.sqrt(GAMMA * p / r)), t_end - t)
        u1 = u + dt * rhs(u, h, cfg, bc)
        u2 = 0.75 * u + 0.25 * (u1 + dt * rhs(u1, h, cfg, bc))
        u = u / 3.0 + 2.0 / 3.0 * (u2 + dt * rhs(u2, h, cfg, bc))
        t += dt
        r, _, p = prim(u)
        if np.any(r <= 0) or np.any(p <= 0) or not np.all(np.isfinite(u)):
            return None
    return u


def problem(name: str, n: int) -> tuple[np.ndarray, float, float, str]:
    """Initial state, spacing, final time, boundary kind."""
    if name == "order":
        h = 1.0 / n
        x = (np.arange(n) + 0.5) * h
        return cons(1.0 + 0.2 * np.sin(2 * np.pi * x), np.ones(n), np.ones(n)), h, 1.0, "periodic"
    h = 1.0 / n if name in ("sod", "lax") else 10.0 / n
    x = (np.arange(n) + 0.5) * h
    if name == "sod":
        m = x < 0.5
        return cons(np.where(m, 1.0, 0.125), np.zeros(n), np.where(m, 1.0, 0.1)), h, 0.2, "edge"
    if name == "lax":
        m = x < 0.5
        return cons(np.where(m, 0.445, 0.5), np.where(m, 0.698, 0.0), np.where(m, 3.528, 0.571)), h, 0.14, "edge"
    x = x - 5.0
    m = x < -4.0
    r = np.where(m, 3.857143, 1.0 + 0.2 * np.sin(5 * x))
    return cons(r, np.where(m, 2.629369, 0.0), np.where(m, 10.33333, 1.0)), h, 1.8, "edge"


def main() -> int:
    """Run one test."""
    ap = argparse.ArgumentParser()
    ap.add_argument("test", choices=("order", "sod", "lax", "shu-osher"))
    ap.add_argument("--scheme", choices=("split", "hybrid"), default="hybrid")
    ap.add_argument("--rs", choices=("llf", "hll", "hllc"), default="hllc")
    ap.add_argument("--corr", type=int, choices=(0, 4, 6), default=6)
    ap.add_argument("--sensor", choices=("on", "weno", "jump"), default="on")
    ap.add_argument("--vars", choices=("char", "prim"), default="char")
    ap.add_argument("--cells", type=int, default=256)
    a = ap.parse_args()
    cfg = {"scheme": a.scheme, "rs": a.rs, "corr": a.corr, "sensor": a.sensor, "vars": a.vars, "corr_off": 0}
    tag = "split" if a.scheme == "split" else f"hybrid/{a.rs}/c{a.corr}/{a.sensor}/{a.vars}"
    if a.test == "order":
        prev, out = None, []
        for n in (20, 40, 80, 160, 320):
            u0, h, t, bc = problem("order", n)
            u = run(u0.copy(), h, t, cfg, bc, cfl=0.4 * (20 / n) ** (2 / 3))  # dt ~ h^(5/3): time error below space
            err = np.mean(np.abs(u[0] - u0[0]))  # t = 1: one period
            out.append(f"N {n} L1 {err:.3e}" + (f" p {np.log2(prev / err):.2f}" if prev else ""))
            prev = err
        print(tag, "order:", "; ".join(out))
        return 0
    u0, h, t, bc = problem(a.test, a.cells)
    u = run(u0.copy(), h, t, cfg, bc)
    if u is None:
        print(tag, a.test, "FAILED: non-positive state")
        return 1
    ref_n = 8 * a.cells
    cache = Path(tempfile.gettempdir()) / f"flume-hybrid-proto-ref-{a.test}-{ref_n}.npy"
    try:
        uref = np.load(cache)
    except OSError:
        uref0, href, _, _ = problem(a.test, ref_n)
        uref = run(uref0.copy(), href, t, {"scheme": "split", "corr_off": 0}, bc)
        np.save(cache, uref)
    rref = uref[0].reshape(a.cells, -1).mean(axis=1)
    r, _, p = prim(u)
    l1 = np.mean(np.abs(r - rref)) * (h * a.cells)
    print(f"{tag} {a.test} N {a.cells}: L1(rho) {l1:.4e}, max rho {r.max():.4f} (ref {uref[0].max():.4f}), "
          f"min p {p.min():.3e}, corrections off {cfg['corr_off']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
