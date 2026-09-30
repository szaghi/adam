#!/usr/bin/env python3
"""First-order Rusanov (the positivity limiter's theta = 0 backbone) on the Balsara-Spicer blast (#41, M2 S2).

Models (state rho, m, E, B, psi on a periodic 2-D grid, gamma = 1.4):
  glm    FLUME's current mixed GLM: B_n flux psi, psi flux c_h^2 B_n, psi not in E, no sources (Dedner 2002);
  eglm   ideal GLM-MHD, Derigs et al. 2018 eqs. (3.16)-(3.18): psi in B units, B_n flux c_h psi, psi flux c_h B_n,
         energy flux + c_h psi B_n, E includes psi^2/2,
         source -[(div B)(0, B, u.B, u, 0) + (grad psi).(0, 0, u psi, 0, u)];
  none   no cleaning, no sources (B_n flux 0);
  powell 8-wave: none plus -(div B)(0, B, u.B, u) (Powell 1994).
Forward Euler with dt = CFL / max(sum_d alpha_d / dx_d); div B and grad psi by central differences. Reports the first
step with a non-positive density or pressure, or the final time.
Usage: lf_proto.py <model> --cells N --b0 B --p-in P [--cfl C] [--ch C] [--t-end T]
"""

from __future__ import annotations

import argparse

import numpy as np

GAMMA = 1.4


def pressure(u: np.ndarray, model: str) -> np.ndarray:
    """Thermal pressure of the state (psi energy only in EGLM)."""
    ke = 0.5 * (u[1] ** 2 + u[2] ** 2 + u[3] ** 2) / u[0]
    me = 0.5 * (u[5] ** 2 + u[6] ** 2 + u[7] ** 2)
    pe = 0.5 * u[8] ** 2 if model == "eglm" else 0.0
    return (GAMMA - 1.0) * (u[4] - ke - me - pe)


def flux(u: np.ndarray, d: int, model: str, ch: float) -> tuple[np.ndarray, np.ndarray]:
    """Physical flux along axis d (0 = x, 1 = y) and the normal fast speed plus |u_n|."""
    rho = u[0]
    vel = u[1:4] / rho
    b = u[5:8]
    p = pressure(u, model)
    b2 = (b**2).sum(axis=0)
    ptot = p + 0.5 * b2
    un, bn = vel[d], b[d]
    f = np.empty_like(u)
    f[0] = u[1 + d]
    for k in range(3):
        f[1 + k] = u[1 + k] * un - bn * b[k] + (ptot if k == d else 0.0)
    f[4] = un * (u[4] + ptot) - bn * (vel * b).sum(axis=0)
    for k in range(3):
        f[5 + k] = un * b[k] - bn * vel[k]
    f[5 + d] = 0.0
    f[8] = 0.0
    if model == "glm":
        f[5 + d] = u[8]
        f[8] = ch**2 * bn
    elif model == "eglm":
        f[5 + d] = ch * u[8]
        f[8] = ch * bn
        f[4] = f[4] + ch * u[8] * bn
    a2 = GAMMA * np.maximum(p, 0.0) / rho
    bb = b2 / rho
    cf = np.sqrt(0.5 * (a2 + bb + np.sqrt(np.maximum((a2 + bb) ** 2 - 4.0 * a2 * bn**2 / rho, 0.0))))
    return f, np.abs(un) + cf


def run(model: str, n: int, b0: float, p_in: float, cfl: float, ch: float, t_end: float) -> str:
    """Integrate until t_end or the first inadmissible state."""
    h = 1.0 / n
    x = (np.arange(n) + 0.5) * h
    xx, yy = np.meshgrid(x, x, indexing="ij")
    u = np.zeros((9, n, n))
    u[0] = 1.0
    p0 = np.where((xx - 0.5) ** 2 + (yy - 0.5) ** 2 < 0.01, p_in, 0.1)
    u[5] = u[6] = b0 / np.sqrt(2.0)
    u[4] = p0 / (GAMMA - 1.0) + 0.5 * (u[5] ** 2 + u[6] ** 2)
    use_ch = model in ("glm", "eglm")
    t, it = 0.0, 0
    while t < t_end - 1e-15:
        fs, speeds = zip(*(flux(u, d, model, ch) for d in (0, 1)), strict=True)
        alpha = [np.maximum(s, ch) if use_ch else s for s in speeds]
        dt = min(cfl / np.max(alpha[0] / h + alpha[1] / h), t_end - t)
        du = np.zeros_like(u)
        for d in (0, 1):
            ur, fr, ar = np.roll(u, -1, axis=1 + d), np.roll(fs[d], -1, axis=1 + d), np.roll(alpha[d], -1, axis=d)
            face = 0.5 * (fs[d] + fr) - 0.5 * np.maximum(alpha[d], ar) * (ur - u)
            du -= (face - np.roll(face, 1, axis=1 + d)) / h
        if model in ("eglm", "powell"):
            vel = u[1:4] / u[0]
            divb = sum((np.roll(u[5 + d], -1, axis=d) - np.roll(u[5 + d], 1, axis=d)) / (2 * h) for d in (0, 1))
            du[1:4] -= divb * u[5:8]
            du[4] -= divb * (vel * u[5:8]).sum(axis=0)
            du[5:8] -= divb * vel
            if model == "eglm":
                ugp = sum(vel[d] * (np.roll(u[8], -1, axis=d) - np.roll(u[8], 1, axis=d)) / (2 * h) for d in (0, 1))
                du[4] -= ugp * u[8]
                du[8] -= ugp
        u = u + dt * du
        t += dt
        it += 1
        p = pressure(u, model)
        if np.any(u[0] <= 0.0) or np.any(p <= 0.0):
            return (f"{model:6s}: INADMISSIBLE at step {it} t {t:.3e}, {int(np.sum(p <= 0))} cells p <= 0, "
                    f"min p {p.min():.3e}")
    p = pressure(u, model)
    return f"{model:6s}: reached t {t:.3e} in {it} steps, min rho {u[0].min():.3e}, min p {p.min():.3e}"


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("model", choices=("glm", "eglm", "none", "powell"))
    ap.add_argument("--cells", type=int, default=128)
    ap.add_argument("--b0", type=float, default=100.0 / np.sqrt(4.0 * np.pi))
    ap.add_argument("--p-in", type=float, default=1000.0)
    ap.add_argument("--cfl", type=float, default=0.25)
    ap.add_argument("--ch", type=float, default=60.0)
    ap.add_argument("--t-end", type=float, default=0.01)
    a = ap.parse_args()
    print(run(a.model, a.cells, a.b0, a.p_in, a.cfl, a.ch, a.t_end))
