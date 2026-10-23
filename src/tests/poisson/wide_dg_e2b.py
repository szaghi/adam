#!/usr/bin/env python3
"""#84 E2b prototype: the zero-ghost wide D(G) next to a reflected (PEC) axis, and the exact E2b preconditioner.

2-D, x with zero ghosts at both faces (PRISM's Dirichlet), y reflected at both faces (odd: PEC on the scalar, or even:
the PEC normal component). The operator is the wide D(G) of order p; the preconditioner is the wide p = 2 operator
with the same faces, solved exactly by separable transforms: along x the two parity sub-lattices, each a compact
Laplacian at spacing 2h with zero at the ghost centres; along y the p = 2 operator with the reflection. Asserted by
measurement: PCG iterations to 1e-10 against N, for p = 2, 4, 6, 8 and both y parities.
Usage: python3 wide_dg_e2b.py
"""

import numpy as np
from wide_dg_zero_ghost import pcg, x_matrix


def y_matrix(p, n, h, odd):
    """1-D wide D(G) on n cells, ghosts reflected about the faces: odd (Dirichlet-like) or even (Neumann-like)."""
    from wide_dg_proto import dg_stencil

    c = dg_stencil(p)
    a = np.zeros((n, n))
    sgn = -1.0 if odd else 1.0
    for i in range(n):
        for o, co in enumerate(c):
            j = i + o - p
            if co == 0.0:
                continue
            if 0 <= j < n:
                a[i, j] += co
            else:
                jm = -1 - j if j < 0 else 2 * n - 1 - j
                a[i, jm] += sgn * co
    return a / h**2


def main():
    rng = np.random.default_rng(84284)
    print(
        f"   {'p':>2} {'N':>4} {'y':>5} {'min eig':>10} {'PCG, E2b preconditioner':>24}"
    )
    for odd in (True, False):
        for p in (2, 4, 6, 8):
            for n in (32, 64, 128):
                h = 1.0 / n
                ax = -x_matrix(p, n, h, reflect=False)
                ay = -y_matrix(p, n, h, odd)
                # the preconditioner: x the p = 2 zero-ghost operator (two parity sub-lattices, compact at 2h, zero at
                # the ghost centres), y the p = 2 operator with the reflection; separable, inverted by eigenvectors
                mx = -x_matrix(2, n, h, reflect=False)
                my = -y_matrix(2, n, h, odd)
                wx, vx = np.linalg.eigh(mx)
                wy, vy = np.linalg.eigh(my)
                lam = wx[:, None] + wy[None, :]
                tiny = 1e-9 * lam.max()
                # the operator's own null space (y even and x ... ): A = kron(ax, I) + kron(I, ay)
                ex, ux = np.linalg.eigh(ax)
                ey, uy = np.linalg.eigh(ay)
                lam_a = ex[:, None] + ey[None, :]
                min_a = np.abs(lam_a).min()

                def a_op(u, ax=ax, ay=ay):
                    return ax @ u + u @ ay.T

                def m_op(r, vx=vx, vy=vy, lam=lam, tiny=tiny):
                    rh = vx.T @ r @ vy
                    ok = lam > tiny
                    yh = np.where(ok, rh / np.where(ok, lam, 1.0), 0.0)
                    return vx @ yh @ vy.T

                b = rng.standard_normal((n, n))
                null = np.abs(lam_a) < 1e-9 * np.abs(lam_a).max()
                if null.any():
                    bh = ux.T @ b @ uy
                    bh[null] = 0.0
                    b = ux @ bh @ uy.T
                it = pcg(a_op, m_op, b)
                print(
                    f"   {p:>2} {n:>4} {'odd' if odd else 'even':>5} {min_a:10.3e} {it:24d}"
                )


if __name__ == "__main__":
    main()
