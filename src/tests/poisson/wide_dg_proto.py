#!/usr/bin/env python3
"""#84 E2 prototype: the PRISM-consistent wide operator D(G) of order p, and how to solve it.

D(G) composes the order-p centred first derivative with itself per axis (adam_flail_legacy_object, the
compute_smoothing_gauss_seidel_{2nd,4th,6th,8th} stencils). On a uniform periodic grid every operator here is
diagonal in Fourier space, so each preconditioner is applied as its EXACT inverse: the best case of any V-cycle that
approximates it. Measured:
  - the stencil built in real space against the symbol -s_p(theta)^2 / h^2 (the operator is the one PRISM uses);
  - the direct FFT solve: one transform, residual at round-off;
  - PCG iterations to 1e-10 for: none (plain CG), the compact 5-point Laplacian (the "compact-MG" hypothesis of the
    #84 plan), and the wide p = 2 operator (the 2^d sub-lattice compact Laplacians at spacing 2h);
  - the spectral bounds of the two preconditioned operators.
Usage: python3 wide_dg_proto.py
"""

import numpy as np

# centred first-derivative coefficients a_k (k = 1..p/2): D1 q_i = sum_k a_k (q_{i+k} - q_{i-k}) / h
D1 = {
    2: [1 / 2],
    4: [2 / 3, -1 / 12],
    6: [3 / 4, -3 / 20, 1 / 60],
    8: [4 / 5, -1 / 5, 4 / 105, -1 / 280],
}


def s_p(p, th):
    """Symbol of h D1 (times -i): sum_k 2 a_k sin(k theta)."""
    return sum(2 * a * np.sin((k + 1) * th) for k, a in enumerate(D1[p]))


def dg_stencil(p):
    """1-D coefficients of D1 o D1 (h^2 units), offsets -p..p."""
    a = np.zeros(2 * (p // 2) + 1)
    for k, ak in enumerate(D1[p]):
        a[p // 2 + k + 1] = ak
        a[p // 2 - k - 1] = -ak
    return np.convolve(a, a)  # offsets -p..p


def apply_dg_real(p, q, h):
    """D(G) q on a periodic 2-D grid by the real-space stencil."""
    c = dg_stencil(p)
    out = np.zeros_like(q)
    for o, co in enumerate(c):
        off = o - p
        if co == 0.0:
            continue
        out += co * (np.roll(q, -off, axis=0) + np.roll(q, -off, axis=1)) / h**2
    return out


def thetas(n):
    return 2 * np.pi * np.fft.fftfreq(n)


def pcg(apply_a, apply_m, b, tol=1e-10, maxit=20000):
    x = np.zeros_like(b)
    r = b.copy()
    z = apply_m(r)
    p_ = z.copy()
    rz = np.vdot(r, z).real
    nb = np.linalg.norm(b)
    for it in range(1, maxit + 1):
        ap = apply_a(p_)
        al = rz / np.vdot(p_, ap).real
        x += al * p_
        r -= al * ap
        if np.linalg.norm(r) <= tol * nb:
            return it
        z = apply_m(r)
        rz_new = np.vdot(r, z).real
        p_ = z + (rz_new / rz) * p_
        rz = rz_new
    return maxit


def main():
    rng = np.random.default_rng(84)
    print(
        "1) the real-space stencil against the symbol, N = 32, max |D(G) q - IFFT(lambda FFT q)| / max |D(G) q|"
    )
    n = 32
    h = 1.0 / n
    q = rng.standard_normal((n, n))
    t = thetas(n)
    for p in (2, 4, 6, 8):
        lam = -(s_p(p, t)[:, None] ** 2 + s_p(p, t)[None, :] ** 2) / h**2
        ref = np.fft.ifft2(lam * np.fft.fft2(q)).real
        got = apply_dg_real(p, q, h)
        print(
            f"   p = {p}: {np.abs(got - ref).max() / np.abs(got).max():.2e}   stencil {np.round(dg_stencil(p) * 1, 6)}"
        )
    print(
        "\n2) spectral bounds over theta in (0, pi]: s_p^2 / (4 sin^2(theta/2)) (compact preconditioner) and "
        "s_p^2 / sin^2(theta) (wide p=2)"
    )
    th = np.linspace(1e-6, np.pi - 1e-6, 200001)
    for p in (2, 4, 6, 8):
        rc = s_p(p, th) ** 2 / (4 * np.sin(th / 2) ** 2)
        rw = s_p(p, th) ** 2 / np.sin(th) ** 2
        print(
            f"   p = {p}: compact [{rc.min():.3e}, {rc.max():.3f}]  (-> 0 at pi: unbounded condition)   "
            f"wide-2 [{rw.min():.3f}, {rw.max():.3f}]  (kappa {rw.max() / rw.min():.2f})"
        )
    print(
        "\n3) direct FFT solve (the null modes, theta in {0, pi}^2, dropped), residual / |f|, and PCG iterations to "
        "1e-10, periodic 2-D, random source projected onto the range"
    )
    print(
        f"   {'p':>2} {'N':>5} {'FFT residual':>13} {'CG':>6} {'PCG compact':>12} {'PCG wide-2':>11}"
    )
    for p in (2, 4, 6, 8):
        for n in (32, 64, 128, 256):
            h = 1.0 / n
            t = thetas(n)
            k = np.fft.fftfreq(n) * n
            lam1 = s_p(p, t) ** 2
            lam1[np.abs(np.abs(k) - n / 2) < 0.5] = (
                0.0  # theta = pi: exactly 0 (sin(pi) is 1e-16, not 0)
            )
            lam1[k == 0] = 0.0
            lam = -(lam1[:, None] + lam1[None, :]) / h**2
            null = lam == 0.0
            f = rng.standard_normal((n, n))
            fh = np.fft.fft2(f)
            fh[null] = 0.0
            f = np.fft.ifft2(fh).real
            # direct
            uh = np.zeros_like(fh)
            uh[~null] = fh[~null] / lam[~null]
            u = np.fft.ifft2(uh).real
            res = np.abs(apply_dg_real(p, u, h) - f).max() / np.abs(f).max()
            # preconditioners (exact inverses, on the range)
            lc1 = 4 * np.sin(t / 2) ** 2
            lc = -(lc1[:, None] + lc1[None, :]) / h**2
            lw1 = np.sin(t) ** 2
            lw1[np.abs(np.abs(k) - n / 2) < 0.5] = 0.0
            lw = -(lw1[:, None] + lw1[None, :]) / h**2

            def a_op(x, lam=lam):
                return np.fft.ifft2(lam * np.fft.fft2(x)).real

            def inv(lmat):
                def m(x):
                    xh = np.fft.fft2(x)
                    yh = np.zeros_like(xh)
                    ok = lmat != 0.0
                    yh[ok] = xh[ok] / lmat[ok]
                    yh[null] = 0.0
                    return -np.fft.ifft2(yh).real  # -L^-1: SPD for the negated operator

                return m

            def neg_a(x, lam=lam):
                return -a_op(x)

            it_cg = pcg(neg_a, lambda x: x, -f)
            it_c = pcg(neg_a, inv(lc), -f)
            it_w = pcg(neg_a, inv(lw), -f)
            print(f"   {p:>2} {n:>5} {res:13.2e} {it_cg:6d} {it_c:12d} {it_w:11d}")


if __name__ == "__main__":
    main()
