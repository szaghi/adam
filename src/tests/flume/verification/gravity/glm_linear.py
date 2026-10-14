#!/usr/bin/env python3
"""Linear stability of a magnetised shear layer under the GLM and EGLM cleaned MHD equations (issue #83).

Why this script exists: in the magnetised Rayleigh-Taylor runs (rt2d-mhd) the mirror asymmetry grows at a rate that
rises with the cleaning speed glm_ch, by about 1.5 glm_ch, and the growth starts on the shear flanks between bubble
and spike. This asks whether the cleaned PDEs themselves have such a mode: the linearised GLM and EGLM systems about
an exact, divergence-free, sheared and magnetised steady state. If their growth rates do not depend on c_h, the
c_h-proportional growth of FLUME is not a property of the cleaned equations but of their discretisation.

Base state (exact steady ideal MHD, no gravity): rho_0 = 1, u_0 = (0, V(x)), B_0 = (0, b(x)), p_0 = P - b^2 / 2,
psi_0 = 0, on the periodic x interval [-1, 1) holding a double shear layer V = V_0 (tanh((x + a) / w) - tanh((x - a) / w)
- 1), a = 0.5, and b = b_0 + b_s (V / V_0). Perturbations f(x) exp(i k y + sigma t), Fourier collocation in x: the
eigenvalues of a 7 N x 7 N matrix in (rho, u, v, p, B_x, B_y, psi). With s = sigma + i k V (the advective derivative):

EGLM (Derigs et al. 2018: the Powell-type sources -(div B)(0, B, u.B, u) and -(u.grad psi)(0, 0, psi, 0, 1), psi in
field units, advected):
    s rho' + rho_0 (u'_x + i k v') = 0
    rho_0 s u' = -d_x (p' + b B_y') + i k b B_x'
    rho_0 (s v' + V_x u') = -i k p' + b_x B_x'
    s p' + p_0x u' + gamma p_0 (u'_x + i k v') = 0
    s B_x' = i k b u' - c_h psi'_x
    s B_y' = V_x B_x' - b_x u' - b u'_x - i k c_h psi'
    s psi' + c_h (B_x'_x + i k B_y') = -damp psi'
GLM (Dedner et al. 2002, mixed; conservative MHD, psi not advected): the same plus, with D' = B_x'_x + i k B_y', the
terms b D' / rho_0 in y momentum, V D' in B_y and -(gamma - 1) V b D' in the pressure, and sigma in place of s in the
psi equation. At c_h = 0 GLM is the conservative MHD system without the Powell terms, which is not hyperbolic when
div B is free: its growth rate rises in proportion to k (the reported negative control).

Usage:
    glm_linear.py [--n N ...] [--k K ...] [--ch C ...] [--models eglm glm] [--v0 V] [--b0 B] [--bshear B] [--w W]
                  [--damp D]
prints the largest growth rate per model, cleaning speed and wavenumber, at each resolution N.
"""

from __future__ import annotations

import argparse
import math

import numpy as np

LENGTH = 2.0  # the periodic x interval [-1, 1)
HALF_GAP = 0.5  # the shear layers at x = -0.5 and x = 0.5


def diff_matrix(n: int, length: float) -> np.ndarray:
    """Fourier first-derivative matrix on n equispaced points of a periodic interval of the given length."""
    h = 2.0 * math.pi / n
    col = np.zeros(n)
    m = np.arange(1, n)
    col[1:] = 0.5 * (-1.0) ** m / np.tan(m * h / 2.0)
    idx = (np.arange(n)[:, None] - np.arange(n)[None, :]) % n
    return col[idx] * (2.0 * math.pi / length)


def spectrum(model: str, n: int, k: float, ch: float, damp: float = 0.0, v0: float = 0.5, b0: float = 0.3,
             bshear: float = 0.0, w: float = 0.05, gamma: float = 1.4, ptot: float = 1.0) -> np.ndarray:
    """Eigenvalues sigma of the linearised `model` ("eglm" or "glm") about the double shear layer, wavenumber k along
    the flow, cleaning speed ch and psi damping rate damp (FLUME's glm_alpha c_h / glm_damping_length)."""
    x = -1.0 + LENGTH * np.arange(n) / n
    prof = np.tanh((x + HALF_GAP) / w) - np.tanh((x - HALF_GAP) / w) - 1.0  # +1 inside |x| < a, -1 outside
    dprof = (1.0 / np.cosh((x + HALF_GAP) / w) ** 2 - 1.0 / np.cosh((x - HALF_GAP) / w) ** 2) / w
    v, vx = v0 * prof, v0 * dprof
    b, bx_ = b0 + bshear * prof, bshear * dprof
    p0 = ptot - 0.5 * b**2
    p0x = -b * bx_
    dm = diff_matrix(n, LENGTH)
    eye = np.eye(n)
    ik = 1j * k
    dg = np.diag
    mat = np.zeros((7 * n, 7 * n), dtype=complex)

    def put(r: int, c: int, m: np.ndarray) -> None:
        mat[r * n:(r + 1) * n, c * n:(c + 1) * n] += m

    put(0, 1, -dm)                                         # mass (rho_0 = 1)
    put(0, 2, -ik * eye)
    put(1, 3, -dm)                                         # x momentum
    put(1, 5, -(dm @ dg(b)))
    put(1, 4, ik * dg(b))
    put(2, 1, -dg(vx))                                     # y momentum
    put(2, 3, -ik * eye)
    put(2, 4, dg(bx_))
    put(3, 1, -dg(p0x) - gamma * dg(p0) @ dm)              # pressure
    put(3, 2, -gamma * ik * dg(p0))
    put(4, 1, ik * dg(b))                                  # B_x
    put(4, 6, -ch * dm)
    put(5, 4, dg(vx))                                      # B_y
    put(5, 1, -dg(bx_) - dg(b) @ dm)
    put(5, 6, -ik * ch * eye)
    put(6, 4, -ch * dm)                                    # psi
    put(6, 5, -ik * ch * eye)
    put(6, 6, -damp * eye)
    if model == "glm":  # the divergence terms that the Powell-type sources of EGLM remove
        put(2, 4, dg(b) @ dm)
        put(2, 5, dg(b) * ik)
        put(5, 4, dg(v) @ dm)
        put(5, 5, dg(v) * ik)
        put(3, 4, -(gamma - 1.0) * dg(v * b) @ dm)
        put(3, 5, -(gamma - 1.0) * dg(v * b) * ik)
    for r in range(7):  # sigma = s - i k V on every advected row (psi is not advected under glm)
        if not (model == "glm" and r == 6):
            put(r, r, -dg(ik * v))
    return np.linalg.eigvals(mat)


def max_growth(model: str, n: int, k: float, ch: float, **base: float) -> float:
    """The largest growth rate Re(sigma) of `spectrum`."""
    return float(np.max(spectrum(model, n, k, ch, **base).real))


def main() -> int:
    """Print the largest growth rate per model, cleaning speed and wavenumber, at each resolution."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--n", type=int, nargs="+", default=[64, 128])
    ap.add_argument("--k", type=float, nargs="+", default=[2 * math.pi, 4 * math.pi])
    ap.add_argument("--ch", type=float, nargs="+", default=[0.0, 1.5, 3.0, 6.0])
    ap.add_argument("--v0", type=float, default=0.5)
    ap.add_argument("--b0", type=float, default=0.3)
    ap.add_argument("--bshear", type=float, default=0.0)
    ap.add_argument("--w", type=float, default=0.05)
    ap.add_argument("--damp", type=float, default=0.0, help="psi damping rate glm_alpha c_h / L")
    ap.add_argument("--models", nargs="+", choices=("eglm", "glm"), default=["eglm", "glm"])
    a = ap.parse_args()
    base = {"damp": a.damp, "v0": a.v0, "b0": a.b0, "bshear": a.bshear, "w": a.w}
    for model in a.models:
        for ch in a.ch:
            cells = [f"k = {k:6.2f}: " + " / ".join(f"{max_growth(model, n, k, ch, **base):7.3f}" for n in a.n)
                     for k in a.k]
            print(f"{model:4s} c_h {ch:4.1f}  " + "  ".join(cells), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
