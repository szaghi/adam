#!/usr/bin/env python3
"""Wu-Shu positivity backbone test (M3 gate): first-order LF with the Wu (2018) speed and the Godunov-Powell source.

Wu (2018, SINUM) / Wu & Shu (2018, SISC): the multi-D first-order Lax-Friedrichs scheme for ideal MHD is positivity
preserving for ANY admissible input (no discrete divergence-free condition) when (a) the Powell source
-(div B)(0, B, u.B, u) is added with div B by central differences of the cell values, and (b) the LF speed exceeds
  sigma(U, V) = max(|u_n| + c_f (U), |u_n| + c_f (V), |rho-weighted u_n| + max c_f)
                + |B_U - B_V| / (sqrt rho_U + sqrt rho_V),
under dt * sum_d sigma_d / dx_d <= CFL. This script checks that claim before any Fortran:
  random   one step on random admissible states with huge cell-to-cell B jumps (the theorem's setting);
  blast    WENO5 + SSP-RK3 + cell-based (Xu 2014) limiter on the Balsara-Spicer blast, whose first-order backbone is
           the scheme above.
Options: --speed std|wu (LF speed), --global (one sigma per direction, as in the proof), --no-powell,
         --model none|glm|eglm (eglm: Derigs et al. 2018 with its sources, psi in B units),
         --first-order (the backbone alone).
"""

from __future__ import annotations

import argparse
import itertools
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from lf_proto import GAMMA, flux, pressure  # noqa: E402
from weno_proto import weno5  # noqa: E402

EPS = 1.0e-13
SHIFTS = (2, 1, 0, -1, -2, -3)
CFG = {"model": "none", "ch": 60.0, "speed": "wu", "global": False, "powell": True}


def admissible(u: np.ndarray) -> np.ndarray:
    """Cells with density and pressure at least EPS."""
    return (u[0] >= EPS) & (pressure(u, CFG["model"]) >= EPS)


def cell_flux(u: np.ndarray, d: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Physical flux along d, |u_n| + c_f, and c_f."""
    f, s = flux(u, d, CFG["model"], CFG["ch"])
    return f, s, s - np.abs(u[1 + d] / u[0])


def face_sigma(u: np.ndarray, d: int, s: np.ndarray, cf: np.ndarray) -> np.ndarray:
    """LF speed at the face i+1/2 along d."""
    sr, cfr, ur = np.roll(s, -1, axis=d), np.roll(cf, -1, axis=d), np.roll(u, -1, axis=1 + d)
    sig = np.maximum(s, sr)
    if CFG["speed"] == "wu":
        ql, qr = np.sqrt(u[0]), np.sqrt(ur[0])
        unw = (ql * u[1 + d] / u[0] + qr * ur[1 + d] / ur[0]) / (ql + qr)
        sig = np.maximum(sig, np.abs(unw) + np.maximum(cf, cfr))
        sig = sig + np.sqrt(((u[5:8] - ur[5:8]) ** 2).sum(axis=0)) / (ql + qr)
    if CFG["model"] in ("glm", "eglm"):
        sig = np.maximum(sig, CFG["ch"])
    if CFG["global"]:
        sig = np.full_like(sig, sig.max())
    return sig


def powell(u: np.ndarray, h: float) -> np.ndarray:
    """Godunov-Powell source -(div B)(0, B, u.B, u, 0), div B by central differences."""
    vel = u[1:4] / u[0]
    divb = sum((np.roll(u[5 + d], -1, axis=d) - np.roll(u[5 + d], 1, axis=d)) / (2 * h) for d in (0, 1))
    s = np.zeros_like(u)
    s[1:4] = -divb * u[5:8]
    s[4] = -divb * (vel * u[5:8]).sum(axis=0)
    s[5:8] = -divb * vel
    if CFG["model"] == "eglm":  # Derigs et al. 2018: -(grad psi).(0, 0, u psi, 0, u)
        ugp = sum(vel[d] * (np.roll(u[8], -1, axis=d) - np.roll(u[8], 1, axis=d)) / (2 * h) for d in (0, 1))
        s[4] -= ugp * u[8]
        s[8] -= ugp
    return s if CFG["powell"] else np.zeros_like(u)


def lf_fluxes(u: np.ndarray) -> tuple[list, list, list]:
    """LF face fluxes, face speeds and cell fluxes of both directions."""
    flf, sig, fc = [], [], []
    for d in (0, 1):
        f, s, cf = cell_flux(u, d)
        sg = face_sigma(u, d, s, cf)
        flf.append(0.5 * (f + np.roll(f, -1, axis=1 + d)) - 0.5 * sg * (np.roll(u, -1, axis=1 + d) - u))
        sig.append(sg)
        fc.append(f)
    return flf, sig, fc


def stable_dt(u: np.ndarray, h: float, cfl: float) -> float:
    """dt with dt * sum_d max(sigma at the two faces of the cell) / h <= cfl."""
    _, sig, _ = lf_fluxes(u)
    rate = sum(np.maximum(sig[d], np.roll(sig[d], 1, axis=d)) for d in (0, 1)) / h
    return cfl / rate.max()


def weno_flux(u: np.ndarray, d: int) -> np.ndarray:
    """WENO5 LF-split face flux at i+1/2 along d (alpha over the 6-cell stencil)."""
    f, s, _ = cell_flux(u, d)
    if CFG["model"] in ("glm", "eglm"):
        s = np.maximum(s, CFG["ch"])
    alpha = np.max([np.roll(s, k, axis=d) for k in SHIFTS], axis=0)
    fh = np.empty_like(u)
    for m in range(u.shape[0]):
        fs = {k: np.roll(f[m], k, axis=d) for k in SHIFTS}
        us = {k: np.roll(u[m], k, axis=d) for k in SHIFTS}
        fp = {k: 0.5 * (fs[k] + alpha * us[k]) for k in SHIFTS}
        fm = {k: 0.5 * (fs[k] - alpha * us[k]) for k in SHIFTS}
        fh[m] = weno5(fp[2], fp[1], fp[0], fp[-1], fp[-2]) + weno5(fm[-3], fm[-2], fm[-1], fm[0], fm[1])
    return fh


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


def fe(u: np.ndarray, dt: float, h: float, high: bool, stats: dict) -> np.ndarray:
    """One forward-Euler step: first-order backbone, or WENO blended by the cell-based limiter."""
    lam = dt / h
    flf, _, _ = lf_fluxes(u)
    src = dt * powell(u, h)
    ul = u - lam * sum(flf[d] - np.roll(flf[d], 1, axis=1 + d) for d in (0, 1)) + src
    ok0 = admissible(ul)
    stats["backbone_bad"] = max(stats["backbone_bad"], int(np.sum(~ok0)))
    if not high:
        return ul
    fh = [weno_flux(u, d) for d in (0, 1)]
    contrib = []
    for d in (0, 1):
        df = fh[d] - flf[d]
        contrib += [-lam * df, lam * np.roll(df, 1, axis=1 + d)]
    lam_cell = np.where(ok0, 1.0, 0.0)
    for r in range(1, 5):
        for sub in itertools.combinations(range(4), r):
            lam_cell = np.minimum(lam_cell, theta_max(ul, sum(contrib[k] for k in sub), ok0))
    out = []
    for d in (0, 1):
        th = np.minimum(lam_cell, np.roll(lam_cell, -1, axis=d))
        stats["limited"] += int(np.sum(th < 1.0))
        out.append(flf[d] + th * (fh[d] - flf[d]))
    return u - lam * sum(out[d] - np.roll(out[d], 1, axis=1 + d) for d in (0, 1)) + src


def random_test(n: int, cfl: float, trials: int, seed: int, plog: tuple, vs: float) -> str:
    """Single first-order steps on random admissible states; count inadmissible outputs."""
    rng = np.random.default_rng(seed)
    bad_total, worst = 0, 0.0
    for _ in range(trials):
        rho = np.exp(rng.uniform(-3.0, 3.0, (n, n)))
        vel = rng.normal(0.0, vs, (3, n, n))
        b = rng.normal(0.0, 30.0, (3, n, n))
        p = np.exp(rng.uniform(plog[0], plog[1], (n, n)))
        u = np.zeros((9, n, n))
        u[0], u[1:4], u[5:8] = rho, rho * vel, b
        u[4] = p / (GAMMA - 1.0) + 0.5 * rho * (vel**2).sum(axis=0) + 0.5 * (b**2).sum(axis=0)
        h = 1.0 / n
        v = fe(u, stable_dt(u, h, cfl), h, False, {"backbone_bad": 0, "limited": 0})
        bad = ~admissible(v)
        bad_total += int(bad.sum())
        if bad.any():
            worst = min(worst, float(pressure(v, CFG["model"])[bad].min()))
    return f"random: {trials} trials of {n}^2, inadmissible cells {bad_total}, worst p {worst:.3e}"


def blast(n: int, cfl: float, t_end: float, b0: float, high: bool) -> str:
    """Balsara-Spicer blast, WENO5 + cell limiter over the backbone (or the backbone alone)."""
    h = 1.0 / n
    x = (np.arange(n) + 0.5) * h
    xx, yy = np.meshgrid(x, x, indexing="ij")
    u = np.zeros((9, n, n))
    u[0] = 1.0
    u[5] = u[6] = b0 / np.sqrt(2.0)
    p0 = np.where((xx - 0.5) ** 2 + (yy - 0.5) ** 2 < 0.01, 1000.0, 0.1)
    u[4] = p0 / (GAMMA - 1.0) + 0.5 * (u[5] ** 2 + u[6] ** 2)
    t, it, stats = 0.0, 0, {"backbone_bad": 0, "limited": 0}
    while t < t_end - 1e-15:
        dt = min(stable_dt(u, h, cfl), t_end - t)
        it += 1
        stage = fe(u, dt, h, high, stats)
        for s, (cu, cs) in enumerate(((0.75, 0.25), (1.0 / 3.0, 2.0 / 3.0), (0.0, 0.0)), start=1):
            if not admissible(stage).all():
                return (f"blast: INADMISSIBLE at step {it} stage {s}, min p {pressure(stage, CFG['model']).min():.3e}; "
                        f"max backbone-inadmissible cells {stats['backbone_bad']}, limited faces {stats['limited']}")
            if s == 3:
                break
            stage = cu * u + cs * fe(stage, dt, h, high, stats)
        u = stage
        t += dt
    return (f"blast: reached t {t:.3e} in {it} steps, min rho {u[0].min():.3e}, "
            f"min p {pressure(u, CFG['model']).min():.3e}; max backbone-inadmissible cells {stats['backbone_bad']}, "
            f"limited faces {stats['limited']}")


def main() -> int:
    """Run one test."""
    ap = argparse.ArgumentParser()
    ap.add_argument("test", choices=("random", "blast"))
    ap.add_argument("--model", choices=("none", "glm", "eglm"), default="none")
    ap.add_argument("--speed", choices=("std", "wu"), default="wu")
    ap.add_argument("--global", dest="glob", action="store_true")
    ap.add_argument("--no-powell", action="store_true")
    ap.add_argument("--first-order", action="store_true", help="blast with the backbone alone")
    ap.add_argument("--cells", type=int, default=64)
    ap.add_argument("--cfl", type=float, default=0.4)
    ap.add_argument("--trials", type=int, default=20)
    ap.add_argument("--plog-min", type=float, default=-8.0)
    ap.add_argument("--plog-max", type=float, default=1.0)
    ap.add_argument("--vel", type=float, default=10.0)
    ap.add_argument("--t-end", type=float, default=0.01)
    ap.add_argument("--b0", type=float, default=100.0 / np.sqrt(4.0 * np.pi))
    a = ap.parse_args()
    CFG.update(model=a.model, speed=a.speed, powell=not a.no_powell)
    CFG["global"] = a.glob
    tag = f"{a.model}/{a.speed}{'/global' if a.glob else ''}{'' if CFG['powell'] else '/no-powell'} cfl {a.cfl}"
    if a.test == "random":
        print(tag, random_test(a.cells, a.cfl, a.trials, 1, (a.plog_min, a.plog_max), a.vel))
    else:
        print(tag, blast(a.cells, a.cfl, a.t_end, a.b0, not a.first_order))
    return 0


if __name__ == "__main__":
    sys.exit(main())
