#!/usr/bin/env python3
"""#84 E5c prototype: multigrid levels inside the blocks below the coarsest tree level, and their transfers.

The coarsest level of poisson_mg_object is the whole uniform grid at the minimum leaf level, solved exactly by the FFT
solver in every V-cycle (on rank 0). E5c coarsens it further inside every block, halving the cells per axis down to one
cell per block, and solves only that root grid exactly. This prototype measures, on a uniform 3-D grid of N^3 cells
(blocks of B^3 cells), the FCG iterations to 1e-10 (relative max-norm residual, as the solver's stop) with the V-cycle
of the solver (red-black Gauss-Seidel, 2 sweeps down and up, the V-weighted 7-point operator rediscretised per level,
restriction the transpose of the prolongation), for:
  - "exact": the coarsest level solved exactly (the present solver: 1 iteration on a uniform grid);
  - "pc":    sub-levels down to one cell per block, piecewise-constant prolongation (the solver's tree-level transfer);
  - "tl":    the same with trilinear (cell-centred) prolongation;
on periodic and ghost-Dirichlet faces, for N growing at fixed B (the root grows) and at fixed root (B grows).
Usage: python3 mg_subblock_e5c.py
"""

import numpy as np


def apply_a(x, h, periodic):
    """y = A x: the V-weighted 7-point Laplacian, w = h (face area over spacing, 3-D), zero ghosts or periodic."""
    if periodic:
        s = sum(np.roll(x, sh, ax) for ax in range(3) for sh in (-1, 1))
    else:
        p = np.pad(x, 1)
        s = (p[2:, 1:-1, 1:-1] + p[:-2, 1:-1, 1:-1] + p[1:-1, 2:, 1:-1] + p[1:-1, :-2, 1:-1]
             + p[1:-1, 1:-1, 2:] + p[1:-1, 1:-1, :-2])
    return h * (s - 6.0 * x)


def smooth(e, rhs, h, periodic, order, sweeps=2):
    """Red-black Gauss-Seidel on A e = rhs (diagonal -6 h)."""
    n = e.shape[0]
    i, j, k = np.indices((n, n, n))
    colour = (i + j + k) % 2
    for _ in range(sweeps):
        for c in (order, 1 - order):
            off = apply_a(e, h, periodic) + 6.0 * h * e
            m = colour == c
            e[m] = (rhs[m] - off[m]) / (-6.0 * h)
    return e


def exact(rhs, h, periodic):
    """A e = rhs exactly: FFT (periodic, mean dropped) or the odd extension about the ghost centres (zero ghosts)."""
    n = rhs.shape[0]
    f = rhs / h**3  # A = h^3 L_h: solve L_h e = rhs / V
    if periodic:
        k = np.arange(n)
        lam1 = (2.0 * np.cos(2.0 * np.pi * k / n) - 2.0) / h**2
        lam = lam1[:, None, None] + lam1[None, :, None] + lam1[None, None, :]
        fh = np.fft.fftn(f)
        lam[0, 0, 0] = 1.0
        fh /= lam
        fh[0, 0, 0] = 0.0
        return np.real(np.fft.ifftn(fh))
    m = 2 * (n + 1)
    g = np.zeros((m, m, m))
    g[1:n + 1, 1:n + 1, 1:n + 1] = f
    g[n + 2:, :, :] = -g[n:0:-1, :, :]
    g[:, n + 2:, :] = -g[:, n:0:-1, :]
    g[:, :, n + 2:] = -g[:, :, n:0:-1]
    k = np.arange(m)
    lam1 = (2.0 * np.cos(2.0 * np.pi * k / m) - 2.0) / h**2
    lam = lam1[:, None, None] + lam1[None, :, None] + lam1[None, None, :]
    lam[0, 0, 0] = 1.0
    gh = np.fft.fftn(g) / lam
    gh[0, 0, 0] = 0.0
    return np.real(np.fft.ifftn(gh))[1:n + 1, 1:n + 1, 1:n + 1]


def prolong(c, kind, periodic):
    """Coarse -> fine: piecewise constant, or trilinear cell-centred (3/4, 1/4 per axis; zero beyond a bounded face)."""
    f = c
    for ax in range(3):
        f = np.repeat(f, 2, axis=ax)
        if kind == "tl":
            n = f.shape[ax]
            lo = np.roll(f, 2, ax) if periodic else np.concatenate(
                [np.zeros_like(np.take(f, [0, 1], ax)), np.take(f, range(n - 2), ax)], ax)
            hi = np.roll(f, -2, ax) if periodic else np.concatenate(
                [np.take(f, range(2, n), ax), np.zeros_like(np.take(f, [0, 1], ax))], ax)
            even = (np.arange(n) % 2 == 0).reshape([-1 if a == ax else 1 for a in range(3)])
            f = 0.75 * f + 0.25 * np.where(even, lo, hi)
    return f


def restrict(r, kind, periodic):
    """Fine -> coarse: the transpose of prolong (the sum of the children for piecewise constant)."""
    n = r.shape[0] // 2
    # transpose by linearity: R = P^T, built column by column is too slow; use the adjoint stencil per axis
    out = r
    for ax in range(3):
        m = out.shape[ax]
        if kind == "tl":
            lo = np.roll(out, -2, ax) if periodic else np.concatenate(
                [np.take(out, range(2, m), ax), np.zeros_like(np.take(out, [0, 1], ax))], ax)
            hi = np.roll(out, 2, ax) if periodic else np.concatenate(
                [np.zeros_like(np.take(out, [0, 1], ax)), np.take(out, range(m - 2), ax)], ax)
            even = (np.arange(m) % 2 == 0).reshape([-1 if a == ax else 1 for a in range(3)])
            out = 0.75 * out + 0.25 * np.where(even, lo, hi)
        shape = list(out.shape)
        shape[ax] = m // 2
        shape.insert(ax + 1, 2)
        out = out.reshape(shape).sum(axis=ax + 1)
    assert out.shape[0] == n
    return out


def vcycle(r, h, levels, kind, periodic):
    """One V-cycle from the finest level (h) down `levels` coarsenings, the last solved exactly."""
    if levels == 0:
        return exact(r, h, periodic)
    e = smooth(np.zeros_like(r), r, h, periodic, 0)
    res = r - apply_a(e, h, periodic)
    ec = vcycle(restrict(res, kind, periodic), 2.0 * h, levels - 1, kind, periodic)
    e = e + prolong(ec, kind, periodic)
    return smooth(e, r, h, periodic, 1)


def fcg(n, levels, kind, periodic, rng, tol=1e-10, itmax=200):
    """FCG(1) on A u = f (f of zero mean when periodic), preconditioned by the V-cycle; the iterations to tol."""
    h = 1.0 / n
    f = rng.standard_normal((n, n, n)) * h**3
    if periodic:
        f -= f.mean()
    u = np.zeros_like(f)
    r = f.copy()
    fmax = np.abs(f).max()
    pp = qp = None
    for it in range(1, itmax + 1):
        z = vcycle(r, h, levels, kind, periodic)
        if periodic:
            z -= z.mean()
        p = z.copy()
        if pp is not None:
            p -= np.vdot(z, qp) / np.vdot(pp, qp) * pp
        q = apply_a(p, h, periodic)
        a = np.vdot(p, r) / np.vdot(p, q)
        u += a * p
        r -= a * q
        pp, qp = p, q
        if np.abs(r).max() / fmax <= tol:
            return it
    return itmax


def main():
    rng = np.random.default_rng(845)
    print(f"{'faces':>9} {'N':>4} {'B':>3} {'root':>4} {'sub':>3} {'exact':>5} {'pc':>4} {'tl':>4}")
    for periodic in (True, False):
        for n, b in ((32, 8), (64, 8), (128, 8), (64, 16), (128, 16)):
            sub = int(np.log2(b))
            its = [fcg(n, 0, "pc", periodic, rng)]
            its += [fcg(n, sub, kind, periodic, rng) for kind in ("pc", "tl")]
            print(f"{'periodic' if periodic else 'dirichlet':>9} {n:4d} {b:3d} {n // b:4d} {sub:3d} "
                  f"{its[0]:5d} {its[1]:4d} {its[2]:4d}")


if __name__ == "__main__":
    main()
