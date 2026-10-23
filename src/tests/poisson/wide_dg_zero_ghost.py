#!/usr/bin/env python3
"""#84 E2 prototype, part 2: PRISM's `Dirichlet` elliptic face writes 0 in every ghost layer (not a reflection), so
the wide D(G) operator with it is the infinite stencil truncated: symmetric, but no FFT diagonalises it. Measured on a
2-D grid bounded along x (zero ghosts at both faces) and periodic along y:
  - the smallest eigenvalues of the zero-ghost operator (is it singular, like the reflected one?);
  - PCG iterations to 1e-10 with the reflected operator (odd about the faces: what the FFT solver inverts exactly) as
    the preconditioner, its null modes either dropped (pseudo-inverse) or regularised (the zero eigenvalues replaced by
    the smallest nonzero one); and, as the reference, with the wide p = 2 zero-ghost operator inverted exactly.
Usage: python3 wide_dg_zero_ghost.py
"""

import numpy as np
from wide_dg_proto import dg_stencil, s_p


def x_matrix(p, n, h, reflect):
    """1-D wide D(G) on n cells: ghosts zero (reflect=False) or odd about the faces (reflect=True)."""
    c = dg_stencil(p)
    a = np.zeros((n, n))
    for i in range(n):
        for o, co in enumerate(c):
            j = i + o - p
            if co == 0.0:
                continue
            if 0 <= j < n:
                a[i, j] += co
            elif reflect:  # odd about the face: ghost j < 0 mirrors -1 - j, ghost j >= n mirrors 2n - 1 - j
                jm = -1 - j if j < 0 else 2 * n - 1 - j
                a[i, jm] -= co
    return a / h**2


def pcg(apply_a, apply_m, b, tol=1e-10, maxit=5000):
    x = np.zeros_like(b)
    r = b.copy()
    z = apply_m(r)
    d = z.copy()
    rz = np.vdot(r, z)
    nb = np.linalg.norm(b)
    for it in range(1, maxit + 1):
        ad = apply_a(d)
        al = rz / np.vdot(d, ad)
        x += al * d
        r -= al * ad
        if np.linalg.norm(r) <= tol * nb:
            return it
        z = apply_m(r)
        rz_new = np.vdot(r, z)
        d = z + (rz_new / rz) * d
        rz = rz_new
    return maxit


def main():
    rng = np.random.default_rng(842)
    print(
        f"   {'p':>2} {'N':>4} {'min|eig| zero-ghost':>20} {'PCG reflected (pinv)':>21} "
        f"{'PCG reflected (reg)':>20} {'PCG zero-ghost p=2':>19}"
    )
    for p in (2, 4, 6, 8):
        for n in (32, 64, 128):
            h = 1.0 / n
            ty = 2 * np.pi * np.fft.fftfreq(n)
            ky = np.fft.fftfreq(n) * n
            ly = s_p(p, ty) ** 2
            ly[np.abs(np.abs(ky) - n / 2) < 0.5] = 0.0
            ly[ky == 0] = 0.0
            ly = ly / h**2  # -lambda_y, >= 0
            ax = -x_matrix(p, n, h, reflect=False)  # SPD sign: -D(G)
            mx = -x_matrix(p, n, h, reflect=True)
            wx, vx = np.linalg.eigh(ax)
            mw, mv = np.linalg.eigh(mx)
            # the zero-ghost operator's spectrum: eigenvalues wx_i + ly_j
            lam_a = wx[:, None] + ly[None, :]
            min_a = np.abs(lam_a).min()

            def a_op(u):  # u (n, n): x first index
                ux = ax @ u
                uy = np.fft.ifft(ly[None, :] * np.fft.fft(u, axis=1), axis=1).real
                return ux + uy

            def m_inv(reg, mwx=mw, mvx=mv):
                lam = mwx[:, None] + ly[None, :]
                tiny = 1e-9 * lam.max()
                if reg:
                    lam = np.where(lam > tiny, lam, lam[lam > tiny].min())

                def m(r):
                    rh = np.fft.fft(mvx.T @ r, axis=1)
                    ok = lam > tiny
                    yh = np.where(ok, rh / np.where(ok, lam, 1.0), 0.0)
                    return mvx @ np.fft.ifft(yh, axis=1).real

                return m

            # the wide p = 2 zero-ghost operator, inverted exactly (the reference preconditioner)
            a2 = -x_matrix(2, n, h, reflect=False)
            w2, v2 = np.linalg.eigh(a2)
            ly2 = np.sin(ty) ** 2
            ly2[np.abs(np.abs(ky) - n / 2) < 0.5] = 0.0
            ly2 = ly2 / h**2
            lam2 = w2[:, None] + ly2[None, :]
            tiny2 = 1e-9 * lam2.max()

            def m2(r):
                rh = np.fft.fft(v2.T @ r, axis=1)
                ok = lam2 > tiny2
                yh = np.where(ok, rh / np.where(ok, lam2, 1.0), 0.0)
                return v2 @ np.fft.ifft(yh, axis=1).real

            b = rng.standard_normal((n, n))
            # project b off the zero-ghost operator's null space, if any
            null = np.abs(lam_a) < 1e-9 * np.abs(lam_a).max()
            if null.any():
                bh = np.fft.fft(vx.T @ b, axis=1)
                bh[null] = 0.0
                b = vx @ np.fft.ifft(bh, axis=1).real
            it_p = pcg(a_op, m_inv(False), b)
            it_r = pcg(a_op, m_inv(True), b)
            it_2 = pcg(a_op, m2, b)
            print(f"   {p:>2} {n:>4} {min_a:20.3e} {it_p:21d} {it_r:20d} {it_2:19d}")


if __name__ == "__main__":
    main()
