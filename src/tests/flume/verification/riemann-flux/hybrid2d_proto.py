#!/usr/bin/env python3
"""2-D MHD gate for the M3 scheme (issue #47, M3-P0): hybrid face flux + EGLM + cell-based positivity limiter.

Face flux at i+1/2 along d: WENO5 interpolation of the primitive variables (rho, u, v, w, p, B, psi) to the two face
states; the linear (B_n, psi) subsystem solved exactly at the face (GLM: B_n flux psi, psi flux c_h^2 B_n; EGLM, Derigs
et al. 2018: B_n flux c_h psi, psi flux c_h B_n) and B~_n, psi~ imposed on both states (Mignone & Tzeferacos 2010);
HLL (Davis speeds, min/max of u_n -+ c_f) or LLF of the full system; correction 6th/4th/none (Chen, Toth & Gombosi
2016) times the WENO-weight sensor (off where min_k w_k / d_k < 0.2 in any variable on either side).
Positivity: the cell-based limiter (Xu 2014; Christlieb et al. 2015) toward the first-order LF backbone with the
2nd-order EGLM sources; the high-order update uses the sources of order --src-order (2 or 4). SSP-RK3.
Tests:
  blast  Balsara-Spicer (beta 2.5e-4), t = 0.01; reports the first inadmissible stage or the end.
  cpaw   circularly polarised Alfven wave along (1, 2) on [0, 1]^2 (Toth 2000: rho 1, p 0.1, B_par 1, amplitude 0.1),
         one period; L1 of B and the order over N = 32, 64, 128 (dt ~ h^(5/3)).
  brio-wu periodic double Brio-Wu (use --gamma 2), t = 0.1; L1 and overshoot against LLF on 8x the cells.
Usage: hybrid2d_proto.py <test> [--model glm|eglm|none] [--rs hll|llf] [--corr 6|4|0] [--sensor on|weno]
       [--limiter none|cell] [--src-order 2|4] [--cells N] [--cfl C] [--ch C] [--b0 B]
       [--gamma G]
"""

from __future__ import annotations

import argparse
import itertools
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "mhd" / "positivity-probe"))
import lf_proto  # noqa: E402
from lf_proto import flux, pressure  # noqa: E402

EPS, TAU = 1.0e-13, 0.2
CFG = {"model": "eglm", "ch": 60.0, "rs": "hll", "corr": 6, "sensor": "weno", "limiter": "cell", "src": 2}
STATS = {"backbone_bad": 0, "limited": 0, "corr_off": 0}


def admissible(u: np.ndarray) -> np.ndarray:
    """Cells with density and pressure at least EPS."""
    return (u[0] >= EPS) & (pressure(u, CFG["model"]) >= EPS)


def to_prim(u: np.ndarray) -> np.ndarray:
    """(rho, u, v, w, p, Bx, By, Bz, psi)."""
    w = u.copy()
    w[1:4] = u[1:4] / u[0]
    w[4] = pressure(u, CFG["model"])
    return w


def to_cons(w: np.ndarray) -> np.ndarray:
    """Inverse of to_prim."""
    u = w.copy()
    u[1:4] = w[0] * w[1:4]
    pe = 0.5 * w[8] ** 2 if CFG["model"] == "eglm" else 0.0
    u[4] = w[4] / (lf_proto.GAMMA - 1.0) + 0.5 * w[0] * (w[1:4] ** 2).sum(axis=0) + 0.5 * (w[5:8] ** 2).sum(axis=0) + pe
    return u


def weno5_interp(vm2, vm1, v0, vp1, vp2) -> tuple[np.ndarray, np.ndarray]:
    """WENO5 interpolation at the face between v0 and vp1 (upwind side v0), and min_k w_k / d_k."""
    b0 = 13 / 12 * (vm2 - 2 * vm1 + v0) ** 2 + 0.25 * (vm2 - 4 * vm1 + 3 * v0) ** 2
    b1 = 13 / 12 * (vm1 - 2 * v0 + vp1) ** 2 + 0.25 * (vm1 - vp1) ** 2
    b2 = 13 / 12 * (v0 - 2 * vp1 + vp2) ** 2 + 0.25 * (3 * v0 - 4 * vp1 + vp2) ** 2
    d = (1 / 16, 10 / 16, 5 / 16)
    a = [d[k] / (1e-6 + b) ** 2 for k, b in enumerate((b0, b1, b2))]
    s = a[0] + a[1] + a[2]
    q0 = (3 * vm2 - 10 * vm1 + 15 * v0) / 8
    q1 = (-vm1 + 6 * v0 + 3 * vp1) / 8
    q2 = (3 * v0 + 6 * vp1 - vp2) / 8
    return (a[0] * q0 + a[1] * q1 + a[2] * q2) / s, np.minimum.reduce([a[k] / s / d[k] for k in range(3)])


def fast_speed(w: np.ndarray, d: int) -> np.ndarray:
    """Fast magnetosonic speed along d from primitives."""
    a2 = lf_proto.GAMMA * np.maximum(w[4], 0.0) / w[0]
    bb = (w[5:8] ** 2).sum(axis=0) / w[0]
    return np.sqrt(0.5 * (a2 + bb + np.sqrt(np.maximum((a2 + bb) ** 2 - 4.0 * a2 * w[5 + d] ** 2 / w[0], 0.0))))


def riemann(wl: np.ndarray, wr: np.ndarray, d: int) -> np.ndarray:
    """HLL or LLF flux with the (B_n, psi) subsystem solved exactly."""
    ch, model = CFG["ch"], CFG["model"]
    wl, wr = wl.copy(), wr.copy()
    if model in ("glm", "eglm"):
        bl, br, pl, pr = wl[5 + d], wr[5 + d], wl[8], wr[8]
        if model == "glm":
            bt, pt = 0.5 * (bl + br) - (pr - pl) / (2 * ch), 0.5 * (pl + pr) - 0.5 * ch * (br - bl)
        else:
            bt, pt = 0.5 * (bl + br) - 0.5 * (pr - pl), 0.5 * (pl + pr) - 0.5 * (br - bl)
        wl[5 + d] = wr[5 + d] = bt
        wl[8] = wr[8] = pt
    ul, ur = to_cons(wl), to_cons(wr)
    fl, _ = flux(ul, d, model, ch)
    fr, _ = flux(ur, d, model, ch)
    cl, cr = fast_speed(wl, d), fast_speed(wr, d)
    if CFG["rs"] == "llf":
        s = np.maximum(np.abs(wl[1 + d]) + cl, np.abs(wr[1 + d]) + cr)
        return 0.5 * (fl + fr) - 0.5 * s * (ur - ul)
    sl = np.minimum(wl[1 + d] - cl, wr[1 + d] - cr)
    sr = np.maximum(wl[1 + d] + cl, wr[1 + d] + cr)
    fh = (sr * fl - sl * fr + sl * sr * (ur - ul)) / (sr - sl)
    return np.where(sl >= 0, fl, np.where(sr <= 0, fr, fh))


def hybrid_flux(u: np.ndarray, d: int) -> np.ndarray:
    """Hybrid face flux at i+1/2 along d."""
    ax = 1 + d
    w = to_prim(u)
    ws = {k: np.roll(w, -k, axis=ax) for k in range(-2, 4)}  # cell i+k at index i
    wl, sl = weno5_interp(ws[-2], ws[-1], ws[0], ws[1], ws[2])
    wr, sr = weno5_interp(ws[3], ws[2], ws[1], ws[0], ws[-1])
    fr = riemann(wl, wr, d)
    f, _ = flux(u, d, CFG["model"], CFG["ch"])
    fs = {k: np.roll(f, -k, axis=ax) for k in (-1, 0, 1, 2)}
    if CFG["corr"] == 6:
        c = 19 / 45 * fr - 13 / 60 * (fs[0] + fs[1]) + 1 / 180 * (fs[-1] + fs[2])
    elif CFG["corr"] == 4:
        c = 1 / 3 * fr - 1 / 6 * (fs[0] + fs[1])
    else:
        c = np.zeros_like(fr)
    if CFG["sensor"] == "weno":
        s = (np.minimum(sl.min(axis=0), sr.min(axis=0)) >= TAU).astype(float)
        STATS["corr_off"] += int(np.sum(s < 1.0))
        c = s * c
    return fr + c


def lf_flux(u: np.ndarray, d: int) -> tuple[np.ndarray, np.ndarray]:
    """First-order LF face flux at i+1/2 along d, and |u_n| + c_f (at least c_h) of the cells."""
    f, s = flux(u, d, CFG["model"], CFG["ch"])
    if CFG["model"] in ("glm", "eglm"):
        s = np.maximum(s, CFG["ch"])
    ax = 1 + d
    a = np.maximum(s, np.roll(s, -1, axis=d))
    return 0.5 * (f + np.roll(f, -1, axis=ax)) - 0.5 * a * (np.roll(u, -1, axis=ax) - u), s


def deriv(v: np.ndarray, d: int, h: float, order: int) -> np.ndarray:
    """Central first derivative of order 2 or 4 along d (periodic)."""
    if order == 2:
        return (np.roll(v, -1, axis=d) - np.roll(v, 1, axis=d)) / (2 * h)
    num = -np.roll(v, -2, axis=d) + 8 * np.roll(v, -1, axis=d) - 8 * np.roll(v, 1, axis=d) + np.roll(v, 2, axis=d)
    return num / (12 * h)


def source(u: np.ndarray, h: float, order: int) -> np.ndarray:
    """EGLM nonconservative source (zero for the other models)."""
    s = np.zeros_like(u)
    if CFG["model"] != "eglm":
        return s
    vel = u[1:4] / u[0]
    divb = sum(deriv(u[5 + d], d, h, order) for d in (0, 1))
    ugp = sum(vel[d] * deriv(u[8], d, h, order) for d in (0, 1))
    s[1:4] = -divb * u[5:8]
    s[4] = -divb * (vel * u[5:8]).sum(axis=0) - ugp * u[8]
    s[5:8] = -divb * vel
    s[8] = -ugp
    return s


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


def fe(u: np.ndarray, dt: float, h: float) -> np.ndarray:
    """One forward-Euler step with the limited hybrid fluxes."""
    lam = dt / h
    fh = [hybrid_flux(u, d) for d in (0, 1)]
    sh = dt * source(u, h, CFG["src"])
    if CFG["limiter"] == "none":
        return u - lam * sum(fh[d] - np.roll(fh[d], 1, axis=1 + d) for d in (0, 1)) + sh
    fl = [lf_flux(u, d)[0] for d in (0, 1)]
    sl = dt * source(u, h, 2)
    ul = u - lam * sum(fl[d] - np.roll(fl[d], 1, axis=1 + d) for d in (0, 1)) + sl
    ok0 = admissible(ul)
    STATS["backbone_bad"] = max(STATS["backbone_bad"], int(np.sum(~ok0)))
    contrib = [sh - sl]
    for d in (0, 1):
        df = fh[d] - fl[d]
        contrib += [-lam * df, lam * np.roll(df, 1, axis=1 + d)]
    lam_cell = np.where(ok0, 1.0, 0.0)
    for r in range(1, 6):
        for sub in itertools.combinations(range(5), r):
            lam_cell = np.minimum(lam_cell, theta_max(ul, sum(contrib[k] for k in sub), ok0))
    out = []
    for d in (0, 1):
        th = np.minimum(lam_cell, np.roll(lam_cell, -1, axis=d))
        STATS["limited"] += int(np.sum(th < 1.0))
        out.append(fl[d] + th * (fh[d] - fl[d]))
    return u - lam * sum(out[d] - np.roll(out[d], 1, axis=1 + d) for d in (0, 1)) + sl + lam_cell * (sh - sl)


def stable_dt(u: np.ndarray, h: float, cfl: float) -> float:
    """dt = cfl / max(sum_d (|u_d| + c_f, at least c_h) / h)."""
    rate = sum(lf_flux(u, d)[1] for d in (0, 1)) / h
    return cfl / rate.max()


def rk3(u: np.ndarray, h: float, t_end: float, cfl: float) -> tuple[np.ndarray | None, str]:
    """SSP-RK3 to t_end; stops at the first inadmissible stage."""
    t, it = 0.0, 0
    while t < t_end - 1e-15:
        dt = min(stable_dt(u, h, cfl), t_end - t)
        it += 1
        s1 = fe(u, dt, h)
        if not admissible(s1).all():
            return None, f"INADMISSIBLE at step {it} stage 1, min p {pressure(s1, CFG['model']).min():.3e}"
        s2 = 0.75 * u + 0.25 * fe(s1, dt, h)
        if not admissible(s2).all():
            return None, f"INADMISSIBLE at step {it} stage 2, min p {pressure(s2, CFG['model']).min():.3e}"
        u = u / 3.0 + 2.0 / 3.0 * fe(s2, dt, h)
        if not admissible(u).all():
            return None, f"INADMISSIBLE at step {it} stage 3, min p {pressure(u, CFG['model']).min():.3e}"
        t += dt
    return u, f"reached t {t:.3e} in {it} steps"


def blast(n: int, cfl: float, b0: float) -> str:
    """Balsara-Spicer blast."""
    h = 1.0 / n
    x = (np.arange(n) + 0.5) * h
    xx, yy = np.meshgrid(x, x, indexing="ij")
    w = np.zeros((9, n, n))
    w[0] = 1.0
    w[4] = np.where((xx - 0.5) ** 2 + (yy - 0.5) ** 2 < 0.01, 1000.0, 0.1)
    w[5] = w[6] = b0 / np.sqrt(2.0)
    u, msg = rk3(to_cons(w), h, 0.01, cfl)
    tail = "" if u is None else f", min rho {u[0].min():.3e}, min p {pressure(u, CFG['model']).min():.3e}"
    return (f"blast {n}^2: {msg}{tail}; max backbone-inadmissible cells {STATS['backbone_bad']}, limited faces "
            f"{STATS['limited']}, corrections off {STATS['corr_off']}")


def brio_wu(n: int, cfl: float) -> str:
    """Periodic double Brio-Wu problem on [0, 2) (left state for 0.5 <= x < 1.5), 8 uniform cells in y, t = 0.1.

    L1(rho) and L1(B_y) against LLF without correction on 8x the cells averaged to the grid, and the overshoot of rho and B_y
    above the reference maxima (the oscillation measure of risk R-1).
    """

    def state(m: int) -> tuple[np.ndarray, float]:
        h = 2.0 / m
        x = (np.arange(m) + 0.5) * h
        left = (x >= 0.5) & (x < 1.5)
        w = np.zeros((9, m, 8))
        w[0] = np.where(left, 1.0, 0.125)[:, None]
        w[4] = np.where(left, 1.0, 0.1)[:, None]
        w[5] = 0.75
        w[6] = np.where(left, 1.0, -1.0)[:, None]
        return to_cons(w), h

    u0, h = state(n)
    u, msg = rk3(u0, h, 0.1, cfl)
    if u is None:
        return f"brio-wu N {n}: {msg}"
    r0, rh = state(8 * n)
    saved = dict(CFG)
    CFG.update(rs="llf", corr=0, sensor="on", limiter="none")  # one non-oscillatory reference for every variant
    ref, _ = rk3(r0, rh, 0.1, cfl)
    CFG.update(saved)
    rho, by = u[0, :, 0], u[6, :, 0]
    rho_ref, by_ref = ref[0, :, 0].reshape(n, 8).mean(axis=1), ref[6, :, 0].reshape(n, 8).mean(axis=1)
    return (f"brio-wu N {n}: {msg}; L1(rho) {np.mean(np.abs(rho - rho_ref)):.4e}, L1(By) "
            f"{np.mean(np.abs(by - by_ref)):.4e}, overshoot rho {rho.max() - ref[0].max():+.3e}, "
            f"By {by.max() - ref[6].max():+.3e}, corrections off {STATS['corr_off']}")


def cpaw_state(n: int, t: float) -> tuple[np.ndarray, float]:
    """Circularly polarised Alfven wave along (1, 2)/sqrt 5 at time t (v_A = 1, travelling toward -x_par).

    Not along a diagonal: there, B depends on x + y only and the discrete div B cancels exactly, which hides the source.
    """
    h = 1.0 / n
    x = (np.arange(n) + 0.5) * h
    xx, yy = np.meshgrid(x, x, indexing="ij")
    phase = 2 * np.pi * ((xx + 2 * yy) + np.sqrt(5.0) * t)  # 2 pi (x_par + t) / lambda, lambda = 1/sqrt 5
    e_par = np.array([1.0, 2.0]) / np.sqrt(5.0)
    e_perp = np.array([-2.0, 1.0]) / np.sqrt(5.0)
    b_perp, b_3 = 0.1 * np.sin(phase), 0.1 * np.cos(phase)
    w = np.zeros((9, n, n))
    w[0], w[4] = 1.0, 0.1
    for c in (0, 1):
        w[5 + c] = e_par[c] + e_perp[c] * b_perp
        w[1 + c] = e_perp[c] * b_perp
    w[7], w[3] = b_3, b_3
    return w, h


def cpaw(cfl: float) -> str:
    """Order of the oblique CPAW over one period."""
    period = 1.0 / np.sqrt(5.0)
    out, prev = [], None
    for n in (32, 64, 128):
        w0, h = cpaw_state(n, 0.0)
        u, msg = rk3(to_cons(w0), h, period, cfl * (32 / n) ** (2 / 3))
        if u is None:
            return f"cpaw N {n}: {msg}"
        err = np.mean(np.abs(to_prim(u)[5:8] - cpaw_state(n, period)[0][5:8]))
        out.append(f"N {n} L1(B) {err:.3e}" + (f" p {np.log2(prev / err):.2f}" if prev else ""))
        prev = err
    return "cpaw: " + "; ".join(out)


def main() -> int:
    """Run one test."""
    ap = argparse.ArgumentParser()
    ap.add_argument("test", choices=("blast", "cpaw", "brio-wu"))
    ap.add_argument("--model", choices=("glm", "eglm", "none"), default="eglm")
    ap.add_argument("--rs", choices=("hll", "llf"), default="hll")
    ap.add_argument("--corr", type=int, choices=(0, 4, 6), default=6)
    ap.add_argument("--sensor", choices=("on", "weno"), default="weno")
    ap.add_argument("--limiter", choices=("none", "cell"), default="cell")
    ap.add_argument("--src-order", type=int, choices=(2, 4), default=2)
    ap.add_argument("--cells", type=int, default=64)
    ap.add_argument("--cfl", type=float, default=0.4)
    ap.add_argument("--ch", type=float, default=None)
    ap.add_argument("--gamma", type=float, default=1.4)
    ap.add_argument("--b0", type=float, default=100.0 / np.sqrt(4.0 * np.pi), help="blast |B|")
    a = ap.parse_args()
    CFG.update(model=a.model, rs=a.rs, corr=a.corr, sensor=a.sensor, limiter=a.limiter, src=a.src_order)
    CFG["ch"] = a.ch if a.ch is not None else {"cpaw": 2.0, "brio-wu": 3.0}.get(a.test, 60.0)
    lf_proto.GAMMA = a.gamma
    tag = f"{a.model}/{a.rs}/c{a.corr}/{a.sensor}/{a.limiter}/src{a.src_order} cfl {a.cfl}"
    if a.test == "blast":
        print(tag, blast(a.cells, a.cfl, a.b0))
    elif a.test == "cpaw":
        print(tag, cpaw(a.cfl))
    else:
        print(tag, brio_wu(a.cells, a.cfl))
    return 0


if __name__ == "__main__":
    sys.exit(main())
