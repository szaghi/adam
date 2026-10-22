#!/usr/bin/env python3
"""FLUME isothermal and barotropic MHD oracle (issue #100, M13 P5): the eigensystem and the HLLD solver, checked.

    mhd_oracle.py selftest        the oracle checks its two literature sources against properties independent of them

The system (Mignone 2007, eqs. 1-6; Kim et al. 1999, eqs. 2-1..2-8) is ideal MHD without an energy equation, the
total pressure p_T = P(rho) + |B|^2 / 2 with P from the [eos] law of eos_oracle.Eos (isothermal P = c0^2 rho, or the
barotropic piecewise polytrope). FLUME keeps its MHD state layout (D-M13-1): U = (rho, rho u, rho v, rho w, slot5, Bx,
By, Bz[, psi]), slot 5 inert (zero flux, its eigenvector the unit vector at speed 0), the normal field a parameter of
the 1-D core (its row decoupled at speed 0, or the GLM pair at -+c_h), as the adiabatic library does.

The 1-D core (normal n, tangential t1, t2) is the six-wave system on (rho, rho u_n, rho u_t1, rho u_t2, B_t1, B_t2).
Its Jacobian holds the gas law only through dP/drho (Kim's row 2: a^2 - u^2), so the isothermal eigensystem of Kim et
al. (1999, eqs. 2-9..2-53) is the barotropic one with a^2 -> c^2 = dP/drho at the state's density: the oracle checks
both gases. The wave speeds (u_n -+ c_f, u_n -+ c_a, u_n -+ c_s) are those of adiabatic MHD with a^2 in place of
gamma p / rho; there is no entropy mode.

Kim et al. normalise the fast pair with alpha_f = sqrt(c_f^2 - b_x^2) / sqrt(c_f^2 - c_s^2) (eq. 2-24), which is
Ryu & Jones's sqrt(a^2 - c_s^2) / sqrt(c_f^2 - c_s^2) times c_f / a; the left vectors carry the matching theta_1,
theta_2 (eqs. 2-51, 2-52). The self-test decides whether the printed set is consistent (L R = I and R Lambda L = A).

The HLLD solver (Mignone 2007, section 3): three states between the outer fast speeds S_L, S_R and the rotational
speeds S*_{L,R} = u* -+ |B_x| / sqrt(rho*) (eq. 29); rho*, m_x* and the first two flux components are the HLL averages
(eqs. 20-23), u* = F^hll_rho / rho^hll (the advective velocity), the transverse components from eqs. 30-37, the flux
sampled by eq. 38. Eq. 39 is printed with S_R = min(lambda_f+(U_L), lambda_f+(U_R)): it must be the max (the outer
speed bounds both states' fast speeds); the oracle, as the FLUME libraries, takes the wider Miyoshi & Kusano (2005,
eq. 67) speeds, which bound Davis's. The degenerate case (S_a - S*_L)(S_a - S*_R) -> 0
(no tangential field, B_x^2 > a^2 rho) takes zero transverse jumps (Mignone, section 3, first bullet). The speeds
can be out of order (S*_R >= S_R or S*_L <= S_L): with u* = F_rho^hll / rho^hll a strong density jump drives the
advective velocity far off, and on Mignone's own tube 1 (rho_L / rho_R = 10) the first face has S*_R = 3.90 > S_R =
3.37 (Davis's estimate is the same there), a star region of negative width that the paper does not discuss; there,
as the adiabatic FLUME HLLD does, the flux is HLL's.

The self-test (all to round-off unless stated):
  kim       L R = I and R Lambda L = A (A by central differences of the core flux, 1e-7) at random states, isothermal
            and barotropic, and at the degenerate points: b_t = 0 (b_x^2 above and below c^2), b_x = 0, b_t = 0 with
            b_x^2 = c^2 (the triple umbilic, eq. 2-29); both as printed and in the FLUME form (kim_eigensystem);
  embed     the 8- and 9-component eigensystems (slot 5, the decoupled B_n, the GLM pair) satisfy L R = I in all
            three directions;
  speeds    the eigenvalues are the fast, Alfven and slow speeds of eqs. 2-15..2-17, ordered;
  consist   HLLD(U, U) = HLL(U, U) = LLF(U, U) = F(U), subsonic and supersonic, both signs of B_x;
  rotation  an isolated rotational discontinuity (rho, u_n, |B_t| continuous, the Walen relation) is resolved exactly:
            the HLLD flux is the upwind flux of the exact solution; HLL is not;
  hydro     with B = 0 the HLLD flux is the isothermal HLLC flux of eos_oracle (the attribution P1 deferred);
  positive  the HLLD star density is positive on Mignone's tubes 1-3 (Table 1, H = sqrt(4 pi) B);
  cpaw      the circularly polarised Alfven wave is an exact solution: its residual dU/dt + dF/dx vanishes pointwise.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from eos_oracle import POLY, Eos, check, flux_hllc, isothermal  # noqa: E402

NV = 8  # rho, rho u, rho v, rho w, slot5, Bx, By, Bz
IB = 5  # first magnetic slot
FRAME = {0: (1, 2, 3), 1: (2, 3, 1), 2: (3, 1, 2)}  # normal, tangential 1, tangential 2 (momentum slots)


# --------------------------------------------------------------------------------------------------------------------
# the full state (FLUME layout)
# --------------------------------------------------------------------------------------------------------------------
def state(rho: float, v: tuple[float, float, float], b: tuple[float, float, float]) -> np.ndarray:
    """U = (rho, rho v, slot5 = 0, B)."""
    return np.array([rho, rho * v[0], rho * v[1], rho * v[2], 0.0, b[0], b[1], b[2]])


def flux(eos: Eos, q: np.ndarray, d: int = 0) -> np.ndarray:
    """The ideal MHD flux in direction d, total pressure P(rho) + |B|^2 / 2; slot 5 and B_n carry no flux."""
    rho = q[0]
    vel = q[1:4] / rho
    bf = q[IB : IB + 3]
    pt = eos.p(rho) + 0.5 * float(bf @ bf)
    f = np.zeros(NV)
    f[0] = q[1 + d]
    f[1:4] = q[1:4] * vel[d] - bf[d] * bf
    f[1 + d] += pt
    f[IB : IB + 3] = bf * vel[d] - bf[d] * vel
    f[IB + d] = 0.0
    return f


def fast_speed(eos: Eos, q: np.ndarray, d: int) -> float:
    rho = q[0]
    c2 = eos.c2(rho)
    b2 = float(q[IB : IB + 3] @ q[IB : IB + 3]) / rho
    bn2 = q[IB + d] ** 2 / rho
    return math.sqrt(0.5 * (c2 + b2 + math.sqrt(max((c2 + b2) ** 2 - 4.0 * c2 * bn2, 0.0))))


# --------------------------------------------------------------------------------------------------------------------
# the 1-D core of Kim et al. (1999): (rho, rho u, rho v, rho w, By, Bz), x normal, Bx a parameter
# --------------------------------------------------------------------------------------------------------------------
def core_flux(eos: Eos, w: np.ndarray, bx: float) -> np.ndarray:
    rho, mx, my, mz, by, bz = w
    u, v, ww = mx / rho, my / rho, mz / rho
    pt = eos.p(rho) + 0.5 * (bx * bx + by * by + bz * bz)
    return np.array([mx, mx * u + pt - bx * bx, my * u - bx * by, mz * u - bx * bz, by * u - bx * v, bz * u - bx * ww])


def core_jacobian(eos: Eos, w: np.ndarray, bx: float, h: float = 1e-7) -> np.ndarray:
    a = np.zeros((6, 6))
    for k in range(6):
        e = np.zeros(6)
        e[k] = h * max(abs(w[k]), 1.0)
        a[:, k] = (core_flux(eos, w + e, bx) - core_flux(eos, w - e, bx)) / (2.0 * e[k])
    return a


def kim_eigensystem(eos: Eos, w: np.ndarray, bx: float, flume: bool = False
                    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(lambda, R, L) of Kim et al. (1999), eqs. 2-9..2-53, ordered u-cf, u-ca, u-cs, u+cs, u+ca, u+cf, a^2 = dP/drho.

    `flume` selects the form the Fortran library implements: the same vectors with alpha_f = alf c_f / a, alpha_s =
    als, (alf, als) the cancellation-free Ryu-Jones pair of the adiabatic library (issue #79), and its degenerate values:
    the transverse direction (beta_y, beta_z) = (1, 0) at b_t = 0, (alf, als) = (1, 0) at the triple umbilic."""
    rho, mx, my, mz, by, bz = w
    u, v, ww = mx / rho, my / rho, mz / rho
    sr = math.sqrt(rho)
    a2 = eos.c2(rho)
    a = math.sqrt(a2)
    bxs, bys, bzs = bx / sr, by / sr, bz / sr
    b2 = bxs * bxs + bys * bys + bzs * bzs
    disc = math.sqrt(max((a2 + b2) ** 2 - 4.0 * a2 * bxs * bxs, 0.0))
    cf = math.sqrt(0.5 * (a2 + b2 + disc))
    cs = math.sqrt(max(0.5 * (a2 + b2 - disc), 0.0))
    ca = abs(bxs)
    sbx = 1.0 if bx >= 0.0 else -1.0
    bt = math.hypot(bys, bzs)
    if bt > 1e-12 * max(cf, 1.0):
        betay, betaz = bys / bt, bzs / bt
    elif flume:
        betay, betaz = 1.0, 0.0
    else:
        betay = betaz = 1.0 / math.sqrt(2.0)  # eq. 2-28
    if flume:
        ct2, tdif = bys * bys + bzs * bzs, bxs * bxs + bys * bys + bzs * bzs - a2
        cf2_cs2 = math.sqrt(tdif * tdif + 4.0 * a2 * ct2)
        if cf2_cs2 <= 1e-12 * cf * cf:
            alf, als = 1.0, 0.0
        else:
            if tdif >= 0.0:
                amc = 2.0 * a2 * ct2 / (cf2_cs2 + tdif)
                cma = cf2_cs2 - amc
            else:
                cma = 2.0 * a2 * ct2 / (cf2_cs2 - tdif)
                amc = cf2_cs2 - cma
            alf, als = min(1.0, math.sqrt(amc / cf2_cs2)), min(1.0, math.sqrt(cma / cf2_cs2))
        alphaf, alphas = alf * cf / a, als
    elif cf * cf - cs * cs > 1e-12 * cf * cf:
        alphaf = math.sqrt(max(cf * cf - bxs * bxs, 0.0)) / math.sqrt(cf * cf - cs * cs)  # eq. 2-24
        alphas = math.sqrt(max(cf * cf - a2, 0.0)) / math.sqrt(cf * cf - cs * cs)  # eq. 2-25
    else:
        alphaf = alphas = 1.0  # eq. 2-29
    # eq. 2-53: sign(b_T) on the slow pair if a^2 > c_a^2, on the fast pair if a^2 < c_a^2
    if bys > 0.0 or (bys == 0.0 and bzs > 0.0):
        sbt = 1.0
    elif bys < 0.0 or (bys == 0.0 and bzs < 0.0):
        sbt = -1.0
    else:
        sbt = 1.0
    lam = np.array([u - cf, u - ca, u - cs, u + cs, u + ca, u + cf])
    r = np.zeros((6, 6))
    lv = np.zeros((6, 6))
    th1 = 2.0 * (alphaf * alphaf * a2 + alphas * alphas * cf * cf)  # eq. 2-51
    th2 = 2.0 * (alphaf * alphaf * cf * a + alphas * alphas * ca * cs)  # eq. 2-52
    vt = betay * v + betaz * ww
    for col, sg in ((0, -1.0), (5, 1.0)):  # fast, eqs. 2-21, 2-31..2-36
        r[:, col] = [alphaf, alphaf * (u + sg * cf), alphaf * v - sg * alphas * betay * bxs,
                     alphaf * ww - sg * alphas * betaz * bxs, alphas * betay * cf / sr, alphas * betaz * cf / sr]
        lv[col] = [alphaf * a2 / th1 + sg * (-alphaf * a * u + alphas * cs * vt * sbx) / th2, sg * alphaf * a / th2,
                   -sg * alphas * betay * cs * sbx / th2, -sg * alphas * betaz * cs * sbx / th2,
                   alphas * betay * cf * sr / th1, alphas * betaz * cf * sr / th1]
    for col, sg in ((1, -1.0), (4, 1.0)):  # Alfven, eqs. 2-22, 2-38..2-43
        r[:, col] = [0.0, 0.0, -sg * betaz * sbx, sg * betay * sbx, betaz / sr, -betay / sr]
        lv[col] = [sg * 0.5 * (betaz * v - betay * ww) * sbx, 0.0, -sg * 0.5 * betaz * sbx, sg * 0.5 * betay * sbx,
                   0.5 * betaz * sr, -0.5 * betay * sr]
    for col, sg in ((2, -1.0), (3, 1.0)):  # slow, eqs. 2-23, 2-45..2-50
        r[:, col] = [alphas, alphas * (u + sg * cs), alphas * v + sg * alphaf * betay * a * sbx,
                     alphas * ww + sg * alphaf * betaz * a * sbx, -alphaf * betay * a2 / (cf * sr),
                     -alphaf * betaz * a2 / (cf * sr)]
        lv[col] = [alphas * cf * cf / th1 - sg * (alphas * ca * u + alphaf * cf * vt * sbx) / th2,
                   sg * alphas * ca / th2, sg * alphaf * betay * cf * sbx / th2, sg * alphaf * betaz * cf * sbx / th2,
                   -alphaf * betay * cf * sr / th1, -alphaf * betaz * cf * sr / th1]
    if a2 > ca * ca:
        r[:, 2:4] *= sbt
        lv[2:4] *= sbt
    elif a2 < ca * ca:
        r[:, [0, 5]] *= sbt
        lv[[0, 5]] *= sbt
    return lam, r, lv


def to_core(q: np.ndarray, d: int) -> tuple[np.ndarray, float]:
    n, t1, t2 = FRAME[d]
    return np.array([q[0], q[n], q[t1], q[t2], q[IB + t1 - 1], q[IB + t2 - 1]]), q[IB + d]


# --------------------------------------------------------------------------------------------------------------------
# the embedded eigensystem (FLUME layout)
# --------------------------------------------------------------------------------------------------------------------
def eigensystem(eos: Eos, q: np.ndarray, d: int, ch: float | None = None, flume: bool = False
                ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(lambda, R, L) on the 8-component state (no cleaning: B_n decoupled at speed 0) or 9 (GLM: (B_n, psi) at
    -+c_h); slot 5 the unit vector at speed 0."""
    w, bx = to_core(q, d)
    lam6, r6, l6 = kim_eigensystem(eos, w, bx, flume)
    n, t1, t2 = FRAME[d]
    mp = [0, n, t1, t2, IB + t1 - 1, IB + t2 - 1]
    nv = NV if ch is None else NV + 1
    r = np.zeros((nv, nv))
    lv = np.zeros((nv, nv))
    lam = np.zeros(nv)
    for k in range(6):
        lam[k] = lam6[k]
        for i in range(6):
            r[mp[i], k] = r6[i, k]
            lv[k, mp[i]] = l6[k, i]
    r[4, 6] = lv[6, 4] = 1.0
    lam[6] = 0.0
    if ch is None:
        r[IB + d, 7] = lv[7, IB + d] = 1.0
    else:
        r[IB + d, 7], r[NV, 7] = 1.0, -1.0
        r[IB + d, 8], r[NV, 8] = 1.0, 1.0
        lv[7, IB + d], lv[7, NV] = 0.5, -0.5
        lv[8, IB + d], lv[8, NV] = 0.5, 0.5
        lam[7], lam[8] = -ch, ch
    return lam, r, lv


# --------------------------------------------------------------------------------------------------------------------
# Riemann solvers (frame of direction d)
# --------------------------------------------------------------------------------------------------------------------
def outer_speeds(eos: Eos, ql: np.ndarray, qr: np.ndarray, d: int) -> tuple[float, float]:
    """The outer speeds of the FLUME MHD libraries, S_L = min(u_L, u_R) - max(c_fL, c_fR), S_R = max(u_L, u_R) +
    max(c_fL, c_fR) (Miyoshi & Kusano 2005, eq. 67): they bound Mignone's Davis estimates (eq. 39, printed with S_R the
    min: it must be the max)."""
    ul, ur = ql[1 + d] / ql[0], qr[1 + d] / qr[0]
    cf = max(fast_speed(eos, ql, d), fast_speed(eos, qr, d))
    return min(ul, ur) - cf, max(ul, ur) + cf


def flux_llf(eos: Eos, ql: np.ndarray, qr: np.ndarray, d: int = 0) -> np.ndarray:
    """LLF with each state's own B_n (the B_n jump is dissipated, as in the FLUME MHD libraries)."""
    a = max(abs(ql[1 + d] / ql[0]) + fast_speed(eos, ql, d), abs(qr[1 + d] / qr[0]) + fast_speed(eos, qr, d))
    return 0.5 * (flux(eos, ql, d) + flux(eos, qr, d)) - 0.5 * a * (qr - ql)


def flux_hll(eos: Eos, ql: np.ndarray, qr: np.ndarray, d: int = 0) -> np.ndarray:
    """HLL with each state's own B_n (as LLF)."""
    sl, sr = outer_speeds(eos, ql, qr, d)
    if sl >= 0.0:
        return flux(eos, ql, d)
    if sr <= 0.0:
        return flux(eos, qr, d)
    return (sr * flux(eos, ql, d) - sl * flux(eos, qr, d) + sl * sr * (qr - ql)) / (sr - sl)


def flux_hlld(eos: Eos, ql: np.ndarray, qr: np.ndarray, d: int = 0) -> tuple[np.ndarray, float]:
    """Mignone (2007) isothermal HLLD in direction d, both states with the average B_n (as the FLUME libraries; the B_n
    flux is then zero); returns the flux and the star density."""
    n, t1, t2 = FRAME[d]
    ql, qr = ql.copy(), qr.copy()
    ql[IB + d] = qr[IB + d] = 0.5 * (ql[IB + d] + qr[IB + d])
    sl, sr = outer_speeds(eos, ql, qr, d)
    fl, fr = flux(eos, ql, d), flux(eos, qr, d)
    uh = (sr * qr - sl * ql - fr + fl) / (sr - sl)  # eq. 15
    fh = (sr * fl - sl * fr + sl * sr * (qr - ql)) / (sr - sl)  # eq. 17
    rs = uh[0]
    if sl >= 0.0:
        return fl, rs
    if sr <= 0.0:
        return fr, rs
    if rs <= 0.0:
        return flux_hll(eos, ql, qr, d), rs
    us = fh[0] / rs  # the advective velocity
    bx = ql[IB + d]
    sqr = math.sqrt(rs)
    ssl, ssr = us - abs(bx) / sqr, us + abs(bx) / sqr  # eq. 29
    if not (sl < ssl and ssr < sr):  # out of order: a star region of negative width (the FLUME fallback)
        return flux_hll(eos, ql, qr, d), rs

    def star(q: np.ndarray, s: float) -> np.ndarray:  # eqs. 20-21, 30-33
        rho = q[0]
        u = q[n] / rho
        den = (s - ssl) * (s - ssr)
        qs = np.zeros(NV)
        qs[0], qs[n] = rs, uh[n]
        qs[IB + d] = bx
        scale = abs(den) > 1e-12 * max(rho * (s - u) ** 2, bx * bx, 1e-300)
        for t, bt in ((t1, IB + t1 - 1), (t2, IB + t2 - 1)):
            if scale:
                qs[t] = rs * q[t] / rho - bx * q[bt] * (us - u) / den
                qs[bt] = q[bt] / rs * (rho * (s - u) ** 2 - bx * bx) / den
            else:  # degenerate: zero transverse jumps
                qs[t] = rs * q[t] / rho
                qs[bt] = q[bt]
        return qs

    usl, usr = star(ql, sl), star(qr, sr)
    if ssl >= 0.0:
        f = fl + sl * (usl - ql)
    elif ssr <= 0.0:
        f = fr + sr * (usr - qr)
    else:  # the centre state, eqs. 34-37, and its flux, eq. 24
        x = sqr * (1.0 if bx >= 0.0 else -1.0)
        uc = usl.copy()
        for t, bt in ((t1, IB + t1 - 1), (t2, IB + t2 - 1)):
            uc[t] = 0.5 * (usl[t] + usr[t]) + 0.5 * x * (usr[bt] - usl[bt])
            uc[bt] = 0.5 * (usl[bt] + usr[bt]) + 0.5 * (usr[t] - usl[t]) / x
        f = np.zeros(NV)
        f[0], f[n] = fh[0], fh[n]
        for t, bt in ((t1, IB + t1 - 1), (t2, IB + t2 - 1)):
            f[t] = uc[t] * us - bx * uc[bt]
            f[bt] = uc[bt] * us - bx * uc[t] / rs
    f[4] = f[IB + d] = 0.0
    return f, rs


# --------------------------------------------------------------------------------------------------------------------
# exact solutions
# --------------------------------------------------------------------------------------------------------------------
def cpaw(rho0: float, b0: float, eps: float, k: float, x: np.ndarray, t: float, direction: float = 1.0) -> np.ndarray:
    """The circularly polarised Alfven wave along x (an exact nonlinear solution of ideal MHD with any P(rho)):
    rho = rho0, u = 0, B = (b0, eps b0 sin(phi), eps b0 cos(phi)), v_t = -direction B_t / sqrt(rho0),
    phi = k (x - direction c_a t), c_a = b0 / sqrt(rho0). Returns U[NV, x]."""
    ca = b0 / math.sqrt(rho0)
    phi = k * (x - direction * ca * t)
    by, bz = eps * b0 * np.sin(phi), eps * b0 * np.cos(phi)
    q = np.zeros((NV, x.size))
    q[0] = rho0
    q[2] = -direction * rho0 * by / math.sqrt(rho0)
    q[3] = -direction * rho0 * bz / math.sqrt(rho0)
    q[IB] = b0
    q[IB + 1], q[IB + 2] = by, bz
    return q


MIGNONE_TUBES = {  # Mignone (2007), Table 1: rho, u, v, w, Hx, Hy, Hz (H = sqrt(4 pi) B), t_s; a = 1
    1: ((1.0, 0.0, 0.0, 0.0, 3.0, 5.0, 0.0), (0.1, 0.0, 0.0, 0.0, 3.0, 2.0, 0.0), 0.1),
    2: ((1.08, 1.2, 0.01, 0.5, 2.0, 3.6, 2.0), (1.0, 0.0, 0.0, 0.0, 2.0, 4.0, 2.0), 0.2),
    3: ((0.1, 5.0, 0.0, 0.0, 0.0, -1.0, -2.0), (0.1, -5.0, 0.0, 0.0, 0.0, 1.0, 2.0), 0.25),
}


def mignone_state(row: tuple[float, ...]) -> np.ndarray:
    s = math.sqrt(4.0 * math.pi)
    return state(row[0], (row[1], row[2], row[3]), (row[4] / s, row[5] / s, row[6] / s))


# --------------------------------------------------------------------------------------------------------------------
# self-test
# --------------------------------------------------------------------------------------------------------------------
GASES = (("isothermal", isothermal(1.0)), ("barotropic", POLY))


def random_state(rng: np.random.Generator, eos: Eos) -> np.ndarray:
    rho = float(rng.uniform(0.3, 12.0)) if eos is POLY else float(rng.uniform(0.2, 3.0))
    return state(rho, tuple(rng.uniform(-2.0, 2.0, 3)), tuple(rng.uniform(-2.0, 2.0, 3)))


def check_kim() -> bool:
    rng = np.random.default_rng(100)
    ok = True
    for name, eos in GASES:
        cases = [random_state(rng, eos) for _ in range(40)]
        rho = 1.0 if name == "isothermal" else 5.0
        c = eos.c(rho)
        cases += [state(rho, (0.3, -0.2, 0.1), (2.0 * c * math.sqrt(rho), 0.0, 0.0)),  # b_t = 0, b_x^2 > c^2
                  state(rho, (0.3, -0.2, 0.1), (0.5 * c * math.sqrt(rho), 0.0, 0.0)),  # b_t = 0, b_x^2 < c^2
                  state(rho, (0.3, -0.2, 0.1), (0.0, 0.7, -0.4)),  # b_x = 0
                  state(rho, (0.3, -0.2, 0.1), (c * math.sqrt(rho), 0.0, 0.0)),  # triple umbilic
                  state(rho, (0.3, -0.2, 0.1), (-1.3, 0.4, 0.9))]  # negative b_x
        for flume in (False, True):
            worst_i = worst_a = 0.0
            for q in cases:
                w, bx = to_core(q, 0)
                lam, r, lv = kim_eigensystem(eos, w, bx, flume)
                worst_i = max(worst_i, float(np.abs(lv @ r - np.eye(6)).max()))
                a = core_jacobian(eos, w, bx)
                worst_a = max(worst_a, float(np.abs(r @ np.diag(lam) @ lv - a).max() / max(np.abs(a).max(), 1.0)))
            ok &= check("kim", worst_i < 1e-10 and worst_a < 1e-6,
                        f"{name}, {'FLUME form' if flume else 'as printed'}: {len(cases)} states (5 degenerate), "
                        f"max |L R - I| {worst_i:.1e}, max |R Lambda L - A| / max|A| {worst_a:.1e} (A by central "
                        "differences)")
    return ok


def check_embed() -> bool:
    rng = np.random.default_rng(7)
    worst = 0.0
    for _, eos in GASES:
        for d in range(3):
            for ch in (None, 3.0):
                q = random_state(rng, eos)
                if ch is not None:
                    q = np.append(q, 0.2)
                _, r, lv = eigensystem(eos, q[:NV], d, ch)
                worst = max(worst, float(np.abs(lv @ r - np.eye(r.shape[0])).max()))
    return check("embed", worst < 1e-10, f"8 and 9 components, x, y, z: max |L R - I| {worst:.1e}")


def check_speeds() -> bool:
    rng = np.random.default_rng(3)
    ok = True
    for _ in range(30):
        q = random_state(rng, isothermal(1.0))
        w, bx = to_core(q, 0)
        lam, _, _ = kim_eigensystem(isothermal(1.0), w, bx)
        ok &= bool(np.all(np.diff(lam) >= -1e-14)) and abs(lam[5] - w[1] / w[0] - fast_speed(isothermal(1.0), q, 0)) < 1e-12
    return check("speeds", ok, "30 states: ordered, the fast speed of eq. 2-15")


def check_consist() -> bool:
    rng = np.random.default_rng(11)
    worst = 0.0
    for _, eos in GASES:
        for _ in range(30):
            q = random_state(rng, eos)
            for d in range(3):
                f = flux(eos, q, d)
                for g in (flux_hlld(eos, q, q, d)[0], flux_hll(eos, q, q, d), flux_llf(eos, q, q, d)):
                    worst = max(worst, float(np.abs(g - f).max() / max(np.abs(f).max(), 1.0)))
    return check("consist", worst < 1e-13, f"HLLD, HLL, LLF (U, U) = F(U), 60 states x 3 directions: {worst:.1e}")


def rotational_pair(eos: Eos, rho: float, u: float, bx: float, bt: float, th0: float, th1: float, family: float
                    ) -> tuple[np.ndarray, np.ndarray, float]:
    """Two states joined by one rotational discontinuity of family -1 (u - c_a) or +1 (u + c_a): rho, u, |B_t|
    continuous, B_t rotated, v_t from the Walen relation v_t = v0 - family sign(B_x) B_t / sqrt(rho) (the
    B_t row of Rankine-Hugoniot: (s - u) dB_t = -B_x dv_t, s - u = family |B_x| / sqrt(rho))."""
    sr = math.sqrt(rho)
    sgn = 1.0 if bx >= 0.0 else -1.0

    def side(th: float) -> np.ndarray:
        by, bz = bt * math.cos(th), bt * math.sin(th)
        v = 0.1 - family * sgn * by / sr
        w = -0.2 - family * sgn * bz / sr
        return state(rho, (u, v, w), (bx, by, bz))

    return side(th0), side(th1), u + family * abs(bx) / sr


def check_rotation() -> bool:
    ok = True
    for name, eos in GASES:
        rho = 1.0 if name == "isothermal" else 5.0
        worst_d = 0.0
        worst_l = 0.0
        rh = 0.0
        for family, u in ((-1.0, 1.5), (-1.0, 0.3), (1.0, -1.5), (1.0, -0.3)):
            for bx in (0.8, -0.8):
                ql, qr, s = rotational_pair(eos, rho, u, bx, 0.6, 0.3, 2.1, family)
                rh = max(rh, float(np.abs(s * (qr - ql) - (flux(eos, qr, 0) - flux(eos, ql, 0))).max()))
                exact = flux(eos, ql, 0) if s > 0.0 else flux(eos, qr, 0)
                worst_d = max(worst_d, float(np.abs(flux_hlld(eos, ql, qr, 0)[0] - exact).max()))
                worst_l = max(worst_l, float(np.abs(flux_hll(eos, ql, qr, 0) - exact).max()))
        ok &= check("rotation", rh < 1e-13 and worst_d < 1e-12 and worst_l > 1e-3,
                    f"{name}: Rankine-Hugoniot of the pair {rh:.1e}; HLLD against the exact flux {worst_d:.1e}, "
                    f"HLL {worst_l:.1e} (diffusive)")
    return ok


def check_hydro() -> bool:
    rng = np.random.default_rng(5)
    worst = 0.0
    eos = isothermal(1.0)
    for _ in range(40):
        ql = state(float(rng.uniform(0.2, 3.0)), tuple(rng.uniform(-2.0, 2.0, 3)), (0.0, 0.0, 0.0))
        qr = state(float(rng.uniform(0.2, 3.0)), tuple(rng.uniform(-2.0, 2.0, 3)), (0.0, 0.0, 0.0))
        f = flux_hlld(eos, ql, qr, 0)[0]
        g = flux_hllc(eos, ql[:5], qr[:5], 0)
        worst = max(worst, float(np.abs(f[:4] - g[:4]).max()))
    return check("hydro", worst < 1e-13, f"B = 0: HLLD against the isothermal HLLC of eos_oracle, 40 pairs: {worst:.1e}")


def check_positive() -> bool:
    worst = math.inf
    for _, (left, right, _) in MIGNONE_TUBES.items():
        _, rs = flux_hlld(isothermal(1.0), mignone_state(left), mignone_state(right), 0)
        worst = min(worst, rs)
    return check("positive", worst > 0.0, f"Mignone's tubes 1-3: smallest HLLD star density {worst:.4f}")


def check_cpaw() -> bool:
    eos = isothermal(1.0)
    k = 2.0 * math.pi
    x = np.linspace(0.0, 1.0, 9)
    h = 1e-5
    worst = 0.0
    for direction in (1.0, -1.0):
        q = cpaw(1.3, 0.9, 0.1, k, x, 0.37, direction)
        dqdt = (cpaw(1.3, 0.9, 0.1, k, x, 0.37 + h, direction) - cpaw(1.3, 0.9, 0.1, k, x, 0.37 - h, direction)) / (2 * h)
        fp = np.stack([flux(eos, c, 0) for c in cpaw(1.3, 0.9, 0.1, k, x + h, 0.37, direction).T], axis=1)
        fm = np.stack([flux(eos, c, 0) for c in cpaw(1.3, 0.9, 0.1, k, x - h, 0.37, direction).T], axis=1)
        res = dqdt + (fp - fm) / (2 * h)
        worst = max(worst, float(np.abs(res).max() / np.abs(q).max()))
    return check("cpaw", worst < 1e-8, f"residual of the exact wave, both directions: {worst:.1e} (central differences)")


def selftest() -> bool:
    results = [check_kim(), check_embed(), check_speeds(), check_consist(), check_rotation(), check_hydro(),
               check_positive(), check_cpaw()]
    ok = all(results)
    print(f"MHD oracle selftest: {'PASS' if ok else 'FAIL'}")
    return ok


# --------------------------------------------------------------------------------------------------------------------
# runs (the P5 legs, eos/mhd.sh)
# --------------------------------------------------------------------------------------------------------------------
NAMES = ("r", "ru", "rv", "rw", "rE", "bx", "by", "bz", "psi")


def tube_line(work: Path, it: int | None = None) -> tuple[np.ndarray, dict[str, np.ndarray], str]:
    """The cells of a 1-D run along its active axis (first or last checkpoint, or step `it`): FLUME checkpoints store
    the arrays, the origin and the spacing in (z, y, x) order; the other two coordinates are null (one row)."""
    import configparser  # noqa: PLC0415

    import h5py  # noqa: PLC0415

    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None, strict=False)
    ini.read(work / "input.ini")
    ng = int(ini["grid"]["ngc"])
    axis = next(ax for ax in "xyz" if ini["grid"][f"null_{ax}"].strip().lower() == ".false.")
    k = "xyz".index(axis)
    steps = sorted({int(f.name.rsplit("-", 2)[-2]) for f in work.glob("*-proc*.h5")})
    it = steps[-1] if it is None else it
    cells: dict[float, dict[str, float]] = {}
    for path in work.glob(f"*-{it:09d}-proc*.h5"):
        with h5py.File(path, "r") as h5:
            for b in {name.rsplit("-", 1)[0] for name in h5}:
                o, d = h5[f"{b}-origin"][:][::-1], h5[f"{b}-dxdydz"][:][::-1]
                arr = {n: h5[f"{b}-{n}"][:] for n in NAMES if f"{b}-{n}" in h5}
                n = arr["r"].shape[2 - k] - 2 * ng
                for i in range(n):
                    idx = [ng, ng, ng]
                    idx[2 - k] = ng + i
                    cells[round(float(o[k] + (ng + i + 0.5) * d[k]), 10)] = {v: float(a[tuple(idx)]) for v, a in arr.items()}
    xs = np.array(sorted(cells))
    return xs, {v: np.array([cells[x][v] for x in xs]) for v in cells[xs[0]]}, axis


def to_x_frame(f: dict[str, np.ndarray], axis: str) -> dict[str, np.ndarray]:
    """The fields of a run along `axis` in the frame of the run along x (the inverse of make_eos_mhd.py's AXES)."""
    perm = {"x": (0, 1, 2), "y": (2, 0, 1), "z": (1, 2, 0)}[axis]
    out = {"r": f["r"], "rE": f["rE"]}
    for comps in (("ru", "rv", "rw"), ("bx", "by", "bz")):
        for j in range(3):
            out[comps[j]] = f[comps[perm.index(j)]]
    return out


def compare_rotation(works: list[Path]) -> bool:
    """The same tube along x, y and z: every field of the last checkpoint, mapped to the x frame, bitwise equal."""
    ref = to_x_frame(*tube_line(works[0])[1:])
    worst = 0.0
    for w in works[1:]:
        f = to_x_frame(*tube_line(w)[1:])
        worst = max(worst, max(float(np.abs(f[v] - ref[v]).max()) for v in ref))
    return check("rotation", worst == 0.0, f"{len(works)} axes, the last checkpoints in the x frame: max |diff| {worst:.1e}")


def compare_agree(a: Path, b: Path, tol: float) -> bool:
    """Two runs of the same problem (the one-segment table against the isothermal gas): every field of the last
    checkpoint within `tol` of its largest magnitude."""
    fa, fb = tube_line(a)[1], tube_line(b)[1]
    worst = max(float(np.abs(fa[v] - fb[v]).max() / max(np.abs(fa[v]).max(), 1e-300)) for v in fa if v != "rE")
    same_e = float(np.abs(fa["rE"]).max()) == 0.0 == float(np.abs(fb["rE"]).max())
    return check("agree", worst <= tol and same_e, f"{a.name} vs {b.name}: max relative difference {worst:.1e} "
                 f"(bound {tol:g}), rE 0 in both: {same_e}")


def coarsen(x: np.ndarray, y: np.ndarray, n: int) -> np.ndarray:
    """Cell averages of a fine profile on n equal cells (the fine cells nest)."""
    return y.reshape(n, -1).mean(axis=1) if y.size % n == 0 else np.interp((np.arange(n) + 0.5) / n, x, y)


def compare_selfconv(works: list[Path], ref: Path, order_min: float) -> bool:
    """Runs coarse to fine against the finest: the L1 distance of rho (cell averages of the reference) shrinking at an
    order >= order_min over the whole ladder, first run to last (a discontinuous solution: first order at best), rho > 0,
    the energy slot exactly 0. The pair orders are printed, not asserted: where the solution is shocks only (tube 3) the
    L1 distance is set by each shock's position within its cell, which scatters the pairs (0.17 and 1.15 on 200, 400,
    800 cells; 0.66 over the ladder)."""
    xr, fr, _ = tube_line(ref)
    ok, dist, cells = True, [], []
    for w in works:
        x, f, _ = tube_line(w)
        cells.append(x.size)
        dist.append(float(np.abs(f["r"] - coarsen(xr, fr["r"], x.size)).mean()))
        ok &= check("positive", float(f["r"].min()) > 0.0 and float(np.abs(f["rE"]).max()) == 0.0,
                    f"{w.name}: {x.size} cells, min rho {f['r'].min():.4f}, max |rE| {np.abs(f['rE']).max():.1e}")
    orders = [math.log2(dist[i] / dist[i + 1]) for i in range(len(dist) - 1)]
    ladder = math.log(dist[0] / dist[-1]) / math.log(cells[-1] / cells[0]) if len(dist) > 1 else math.inf
    detail = ", ".join(f"{d:.3e}" for d in dist) + "; pair orders " + ", ".join(f"{p:+.2f}" for p in orders)
    if len(dist) > 1:
        detail += f"; ladder order {ladder:+.2f}"
    return check("selfconv", ladder >= order_min, f"L1 rho to {ref.name}: {detail} (min {order_min})") and ok


def compare_gain(hlld: Path, hll: Path, ref: Path) -> bool:
    """HLLD resolves the rotational discontinuities that HLL smears: its L1 distance of B_y to the reference is smaller."""
    xr, fr, _ = tube_line(ref)
    dist = {}
    for name, w in (("hlld", hlld), ("hll", hll)):
        x, f, _ = tube_line(w)
        dist[name] = float(np.abs(f["by"] - coarsen(xr, fr["by"], x.size)).mean())
    return check("hlld-gain", dist["hlld"] < dist["hll"],
                 f"L1 B_y to {ref.name}: HLLD {dist['hlld']:.3e}, HLL {dist['hll']:.3e} "
                 f"(ratio {dist['hlld'] / dist['hll']:.2f})")


def last_fields(work: Path) -> dict[tuple[float, ...], dict[str, np.ndarray]]:
    """Every block of the last checkpoint of a run, by its origin (the block names follow the decomposition, so two runs
    on different rank counts name the same block differently): the interior of each field (arrays stored (z, y, x) with
    ngc ghosts on every side)."""
    import configparser  # noqa: PLC0415

    import h5py  # noqa: PLC0415

    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None, strict=False)
    ini.read(work / "input.ini")
    ng = int(ini["grid"]["ngc"])
    steps = sorted({int(f.name.rsplit("-", 2)[-2]) for f in work.glob("*-proc*.h5")})
    out: dict[tuple[float, ...], dict[str, np.ndarray]] = {}
    for path in work.glob(f"*-{steps[-1]:09d}-proc*.h5"):
        with h5py.File(path, "r") as h5:
            for b in {name.rsplit("-", 1)[0] for name in h5}:
                key = tuple(round(float(o), 12) for o in h5[f"{b}-origin"][:])
                out[key] = {n: h5[f"{b}-{n}"][:][ng:-ng, ng:-ng, ng:-ng] for n in NAMES if f"{b}-{n}" in h5}
    return out


def compare_backends(cpu: Path, fnl: Path, tol: float) -> bool:
    """P5b: the last checkpoint of a CPU run against the FNL run of the same input. The same blocks; every field within
    `tol` of its scale; the energy slot exactly 0 on both. Scales: rho its largest magnitude; the momenta at least
    max(rho) c_ref and B and psi at least max|B|, c_ref = sqrt(max(rho c0^2, |B|^2) / min(rho)) the largest wave speed
    of the background (a wave of amplitude 1e-7 carries momenta of that order, whose round-off is that of the O(1)
    fluxes)."""
    a, b = last_fields(cpu), last_fields(fnl)
    if a.keys() != b.keys():
        return check("backends", False, f"{cpu.name} vs {fnl.name}: different blocks ({len(a)} against {len(b)})")
    big = lambda w, v: max(float(np.abs(w[k][v]).max()) for k in w) if v in next(iter(w.values())) else 0.0  # noqa: E731
    rho_max, rho_min = big(a, "r"), min(float(a[k]["r"].min()) for k in a)
    b_max = math.sqrt(max(float((a[k]["bx"] ** 2 + a[k]["by"] ** 2 + a[k]["bz"] ** 2).max()) for k in a))
    c_ref = math.sqrt(max(rho_max * run_c0(cpu) ** 2, b_max**2) / rho_min)
    report, ok = [], True
    for v in [n for n in NAMES if n != "rE" and n in next(iter(a.values()))]:
        d = max(float(np.abs(a[k][v] - b[k][v]).max()) for k in a)
        scale = big(a, v)
        if v in ("ru", "rv", "rw"):
            scale = max(scale, rho_max * c_ref)
        elif v in ("bx", "by", "bz", "psi"):
            scale = max(scale, b_max)
        gap = d / scale if scale > 0.0 else d
        report.append(f"{v} {gap:.1e}")
        ok = ok and gap <= tol
    inert = max(float(np.abs(w[k]["rE"]).max()) for w in (a, b) for k in w)
    return check("backends", ok and inert == 0.0,
                 f"{cpu.name} vs {fnl.name}: {', '.join(report)} (bound {tol:.0e}); max |rE| {inert:.1e}")


def run_c0(work: Path) -> float:
    """The [eos] sound_speed of a run's input."""
    import configparser  # noqa: PLC0415

    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None, strict=False)
    ini.read(work / "input.ini")
    return float(ini["eos"]["sound_speed"])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest", help="check the eigensystem and the HLLD solver")
    p = sub.add_parser("rotation", help="the same tube along x, y, z: bitwise equal in the x frame")
    p.add_argument("works", type=Path, nargs="+")
    p = sub.add_parser("selfconv", help="tube runs coarse to fine against the finest")
    p.add_argument("works", type=Path, nargs="+")
    p.add_argument("--ref", type=Path, required=True)
    p.add_argument("--order-min", type=float, default=0.5)
    p = sub.add_parser("agree", help="two runs of the same problem agree within --tol")
    p.add_argument("a", type=Path)
    p.add_argument("b", type=Path)
    p.add_argument("--tol", type=float, default=1e-10)
    p = sub.add_parser("backends", help="a CPU run against the FNL run of the same input (P5b)")
    p.add_argument("cpu", type=Path)
    p.add_argument("fnl", type=Path)
    p.add_argument("--tol", type=float, default=1e-10)
    p = sub.add_parser("hlld-gain", help="HLLD closer to the reference than HLL (B_y)")
    p.add_argument("hlld", type=Path)
    p.add_argument("hll", type=Path)
    p.add_argument("--ref", type=Path, required=True)
    args = ap.parse_args()
    if args.cmd == "selftest":
        return 0 if selftest() else 1
    if args.cmd == "rotation":
        return 0 if compare_rotation(args.works) else 1
    if args.cmd == "selfconv":
        return 0 if compare_selfconv(args.works, args.ref, args.order_min) else 1
    if args.cmd == "agree":
        return 0 if compare_agree(args.a, args.b, args.tol) else 1
    if args.cmd == "backends":
        return 0 if compare_backends(args.cpu, args.fnl, args.tol) else 1
    if args.cmd == "hlld-gain":
        return 0 if compare_gain(args.hlld, args.hll, args.ref) else 1
    return 2


if __name__ == "__main__":
    sys.exit(main())
