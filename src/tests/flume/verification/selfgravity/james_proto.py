#!/usr/bin/env python3
"""The isolated Poisson solve on an AMR grid by a symmetric James method: the P0 prototype of issue #87.

Why: self-gravity conserves E + W exactly only if the Poisson solve is a linear map phi = M (V rho) with V M V
symmetric (#81 P4). The multigrid solver of #81 P5a is periodic only; isolated boundaries on AMR need the free-space
solution of the composite operator. The textbook James method (screening charge on the boundary surface, potential
evaluated at the ghost centres) is not symmetric. The discrete form is: let A be the composite operator with zero
Dirichlet ghosts (V-weighted, symmetric, P5a's mimetic mean flux at the 2:1 faces) and B the map of the ghost values
onto the right-hand side of the boundary rows. The zero-ghost solve phi0 = A^-1 (V rho), extended by zero outside,
satisfies the infinite composite operator exactly with an extra source -B^T phi0 on the ghost cells (the screening
charge), so the ghost values of the free-space solution are g = -G B^T phi0, with G the inverse of the infinite
operator from ghost cells to ghost cells. Replacing G by the coarse lattice Green's function (exact when every cell
touching the boundary is at the coarsest level and the grid is uniform) gives

    phi = A^-1 (V rho - B g) = (A^-1 + A^-1 B G B^T A^-1) V rho,

symmetric because A and G are. On a uniform grid it is the exact discrete free-space solution; on AMR G misses the
perturbation of the refined interior, a truncation-order error.

The model: the unit cube, N coarse cells per side (H = 1/N), the box [1/4, 3/4]^3 refined one 2:1 level (h = H/2),
the second-order Laplacian in finite-volume form (T = face area / centre distance; at a 2:1 face the mimetic mean
flux of P5a: every fine face carries (2/3) h (mean of the coarse face's 4 fine cells - phi_C)), one coarse ghost
layer outside each face. G is the lattice Green's function of the 7-point Laplacian (L_1 G = delta, scaled 1/H),
tabulated here by the near-field DST of scripts/lattice_green.py with the G0 + G1 expansion on the faces (numpy only;
its face error, ~1e-9, is far below every measure here and identical in the two solves the uniform check compares).

Checks (exit 1 on a failure):
  1. symmetry: |rho1 . V phi2 - rho2 . V phi1| / |rho1 . V phi2| at round-off (bound 1e-10, the CG tolerance), on the
     AMR grid, random densities;
  2. uniform grid: phi against the direct lattice convolution sum_q G(n_p - n_q) / H V_q rho_q (bound 1e-9);
  3. accuracy on AMR: a Gaussian blob (width 0.1) in the refined box against its exact potential, N = 8, 16, 24, 32:
     the observed order of max |phi - phi_exact| on the finest pair at least 1.8 (a narrower blob, 0.06, is not yet
     resolved at N = 24: its order there is 1.79, and the uniform control's is 0.74);
  4. negative control: the textbook form (the screening charge at the boundary faces, the continuum kernel from the
     faces to the ghost centres) breaks the symmetry by many orders (at least 1e-6).

Usage: james_proto.py [--cells 8 16 24 32]
"""

from __future__ import annotations

import argparse
import math
import sys
from collections.abc import Callable

import numpy as np

BOX = (0.25, 0.75)


def dst1(x: np.ndarray, axis: int) -> np.ndarray:
    """Return the type-I discrete sine transform along `axis` (scripts/lattice_green.py)."""
    x = np.moveaxis(x, axis, -1)
    k = x.shape[-1]
    zero = np.zeros(x.shape[:-1] + (1,))
    y = -np.fft.fft(np.concatenate([zero, x, zero, -x[..., ::-1]], axis=-1), axis=-1).imag[..., 1:k + 1] / 2.0
    return np.moveaxis(y, -1, axis)


def lattice_green(radius: int) -> np.ndarray:
    """Return G on [-R, R]^3 (index n + R), L_1 G = delta, G0 + G1 on the faces, solved by DST-I."""
    m = 2 * radius + 1
    n = np.arange(m, dtype=float) - radius
    n1, n2, n3 = np.meshgrid(n, n, n, indexing="ij")
    face = (np.abs(n1) == radius) | (np.abs(n2) == radius) | (np.abs(n3) == radius)
    g = np.zeros((m, m, m))
    a, b, c = n1[face], n2[face], n3[face]
    r2 = a * a + b * b + c * c
    r = np.sqrt(r2)
    g[face] = -1.0 / (4.0 * math.pi * r) - 5.0 * (a**4 + b**4 + c**4 - 0.6 * r2 * r2) / (32.0 * math.pi * r**7)
    rhs = np.zeros((m - 2,) * 3)
    rhs[radius - 1, radius - 1, radius - 1] = 1.0
    rhs[0] -= g[0, 1:-1, 1:-1]
    rhs[-1] -= g[-1, 1:-1, 1:-1]
    rhs[:, 0] -= g[1:-1, 0, 1:-1]
    rhs[:, -1] -= g[1:-1, -1, 1:-1]
    rhs[:, :, 0] -= g[1:-1, 1:-1, 0]
    rhs[:, :, -1] -= g[1:-1, 1:-1, -1]
    k = np.arange(1, m - 1, dtype=float)
    lam1 = 2.0 * np.cos(math.pi * k / (m - 1)) - 2.0
    lam = lam1[:, None, None] + lam1[None, :, None] + lam1[None, None, :]
    t = rhs
    for ax in range(3):
        t = dst1(t, ax)
    t = t / lam
    for ax in range(3):
        t = dst1(t, ax)
    g[1:-1, 1:-1, 1:-1] = t * (2.0 / (m - 1)) ** 3
    return g


class Grid:
    """The composite grid: coarse cells outside the box, fine cells inside, the coarse ghost layer, the operator."""

    def __init__(self, n: int, refined: bool) -> None:
        self.n, self.big, self.small = n, 1.0 / n, 0.5 / n
        lo, hi = round(BOX[0] * n), round(BOX[1] * n)
        self.idx: dict[tuple, int] = {}
        xs, vol = [], []
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    if refined and lo <= i < hi and lo <= j < hi and lo <= k < hi:
                        for a in range(2):
                            for b in range(2):
                                for c in range(2):
                                    self.idx[("f", 2 * i + a, 2 * j + b, 2 * k + c)] = len(xs)
                                    xs.append(((2 * i + a + 0.5) * self.small, (2 * j + b + 0.5) * self.small,
                                               (2 * k + c + 0.5) * self.small))
                                    vol.append(self.small**3)
                    else:
                        self.idx[("c", i, j, k)] = len(xs)
                        xs.append(((i + 0.5) * self.big, (j + 0.5) * self.big, (k + 0.5) * self.big))
                        vol.append(self.big**3)
        self.x, self.vol, self.m = np.array(xs), np.array(vol), len(xs)
        rows: list[int] = []
        cols: list[int] = []
        vals: list[float] = []

        def add(p: int, q: int, v: float) -> None:
            rows.append(p)
            cols.append(q)
            vals.append(v)

        def couple(p: int, q: int, t: float) -> None:
            add(p, q, t)
            add(p, p, -t)
            add(q, p, t)
            add(q, q, -t)

        # same-level faces along +d; from every coarse cell, the 2:1 faces on both sides (a group: the coarse cell and
        # the 4 fine cells across its face)
        for key, p in self.idx.items():
            lev, *ijk = key
            for d in range(3):
                nb = list(ijk)
                nb[d] += 1
                same = (lev, *nb)
                if same in self.idx:
                    couple(p, self.idx[same], self.big if lev == "c" else self.small)
            if lev == "c" and refined:
                for d in range(3):
                    for side in (-1, 1):
                        self._seam(p, ijk, d, side, add)
        # the ghost layer: the coarse cells just outside each face of the cube, B its coupling to the boundary rows
        self.ghosts: list[tuple[int, int, int]] = []
        brow, bcol = [], []
        for key, p in self.idx.items():
            lev, *ijk = key
            if lev != "c":
                continue
            for d in range(3):
                for side in (-1, 1):
                    if (side < 0 and ijk[d] == 0) or (side > 0 and ijk[d] == n - 1):
                        g = list(ijk)
                        g[d] += side
                        self.ghosts.append((g[0], g[1], g[2]))
                        add(p, p, -self.big)
                        brow.append(p)
                        bcol.append(len(self.ghosts) - 1)
        self.rows, self.cols, self.vals = np.array(rows), np.array(cols), np.array(vals)
        self.brow, self.bcol = np.array(brow), np.array(bcol)
        self.ng = len(self.ghosts)
        on_diag = self.rows == self.cols
        self.diag = np.bincount(self.rows[on_diag], self.vals[on_diag], minlength=self.m)

    def _seam(self, p: int, ijk: list[int], d: int, side: int, add: Callable[[int, int, float], None]) -> None:
        """Add the mimetic mean flux between coarse cell p and the 4 fine cells across its face (d, side), if any."""
        nb = list(ijk)
        nb[d] += side
        f0 = [2 * c for c in nb]
        f0[d] = 2 * nb[d] + (0 if side > 0 else 1)
        others = [e for e in range(3) if e != d]
        fine = []
        for a in range(2):
            for b in range(2):
                f = list(f0)
                f[others[0]] += a
                f[others[1]] += b
                fine.append(self.idx.get(("f", *f)))
        if any(q is None for q in fine):
            return
        w = 2.0 / 3.0 * self.small
        for q in fine:
            add(q, p, w)
            for r in fine:
                add(q, r, -w / 4.0)
            add(p, q, w)
            add(p, p, -w)

    def apply(self, phi: np.ndarray) -> np.ndarray:
        """Return A phi (zero Dirichlet ghosts)."""
        return np.bincount(self.rows, self.vals * phi[self.cols], minlength=self.m)

    def solve(self, rhs: np.ndarray, tol: float = 1e-13) -> np.ndarray:
        """Return A^-1 rhs by Jacobi-preconditioned CG on -A (symmetric positive definite)."""
        b = -rhs
        x = np.zeros(self.m)
        r = b.copy()
        z = r / -self.diag
        p = z.copy()
        rz = r @ z
        nb = math.sqrt(b @ b)
        for _ in range(20000):
            ap = -self.apply(p)
            alpha = rz / (p @ ap)
            x += alpha * p
            r -= alpha * ap
            if math.sqrt(r @ r) <= tol * nb:
                return x
            z = r / -self.diag
            rz_new = r @ z
            p = z + rz_new / rz * p
            rz = rz_new
        sys.exit("james_proto: CG did not converge")

    def b_apply(self, g: np.ndarray) -> np.ndarray:
        """Return B g: the ghost values onto the boundary rows (coefficient H)."""
        return np.bincount(self.brow, self.big * g[self.bcol], minlength=self.m)

    def bt_apply(self, phi: np.ndarray) -> np.ndarray:
        """Return B^T phi: the boundary cells onto their ghosts."""
        return np.bincount(self.bcol, self.big * phi[self.brow], minlength=self.ng)


def ghost_green(grid: Grid, table: np.ndarray, radius: int) -> np.ndarray:
    """Return the ghost-to-ghost block of the inverse of the infinite coarse operator, G(n_q - n_q') / H."""
    gh = np.array(grid.ghosts)
    dn = gh[:, None, :] - gh[None, :, :] + radius
    return table[dn[..., 0], dn[..., 1], dn[..., 2]] / grid.big


def james(grid: Grid, gg: np.ndarray, rho: np.ndarray) -> np.ndarray:
    """Return the symmetric (ghost-charge) James solution phi = A^-1 (V rho - B g), g = -G B^T A^-1 V rho."""
    f = grid.vol * rho
    phi0 = grid.solve(f)
    g = -gg @ grid.bt_apply(phi0)
    return grid.solve(f - grid.b_apply(g))


def james_textbook(grid: Grid, rho: np.ndarray) -> np.ndarray:
    """Return the textbook James solution: the screening charge -B^T phi0 at the boundary faces, the continuum kernel
    -1/(4 pi r) from the faces to the ghost centres (not symmetric: sources and evaluation points differ)."""
    f = grid.vol * rho
    phi0 = grid.solve(f)
    sigma = -grid.bt_apply(phi0)
    gh = (np.array(grid.ghosts) + 0.5) * grid.big
    owner = np.empty((grid.ng, 3))
    owner[grid.bcol] = grid.x[grid.brow]
    faces = 0.5 * (gh + owner)
    dist = np.linalg.norm(gh[:, None, :] - faces[None, :, :], axis=2)
    g = (1.0 / (4.0 * math.pi * dist)) @ sigma
    return grid.solve(f - grid.b_apply(g))


def gaussian(x: np.ndarray, width: float = 0.1) -> tuple[np.ndarray, np.ndarray]:
    """Return a Gaussian blob centred in the cube and its exact free-space potential (lap phi = rho)."""
    r = np.linalg.norm(x - 0.5, axis=1)
    rho = np.exp(-((r / width) ** 2))
    mass = math.pi**1.5 * width**3
    erf = np.array([math.erf(v / width) for v in r])
    safe = np.where(r > 1e-12, r, 1.0)
    centre = -mass / (4.0 * math.pi) * 2.0 / (math.sqrt(math.pi) * width)
    phi = np.where(r > 1e-12, -mass * erf / (4.0 * math.pi * safe), centre)
    return rho, phi


def main() -> int:
    """Run the four checks; return the exit status."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cells", type=int, nargs="+", default=[8, 16, 24, 32])
    args = ap.parse_args()
    radius = max(args.cells) + 4
    table = lattice_green(radius)
    ok = True
    rng = np.random.default_rng(87)

    # 1 and 4: the symmetry of the AMR solve, the textbook form as the negative control
    n = args.cells[1] if len(args.cells) > 1 else args.cells[0]
    grid = Grid(n, refined=True)
    gg = ghost_green(grid, table, radius)
    r1, r2 = rng.standard_normal(grid.m), rng.standard_normal(grid.m)
    p1, p2 = james(grid, gg, r1), james(grid, gg, r2)
    a12, a21 = r1 @ (grid.vol * p2), r2 @ (grid.vol * p1)
    asym = abs(a12 - a21) / abs(a12)
    print(f"1. symmetry on AMR, N = {n} ({grid.m} cells, {grid.ng} ghosts): |r1.V phi2 - r2.V phi1| / |r1.V phi2| = "
          f"{asym:.2e} (bound 1e-10): {'PASS' if asym <= 1e-10 else 'FAIL'}")
    ok &= asym <= 1e-10
    t1, t2 = james_textbook(grid, r1), james_textbook(grid, r2)
    tasym = abs(r1 @ (grid.vol * t2) - r2 @ (grid.vol * t1)) / abs(r1 @ (grid.vol * t2))
    print(f"4. textbook James (surface charge, continuum kernel): asymmetry {tasym:.2e} (at least 1e-6): "
          f"{'PASS' if tasym >= 1e-6 else 'FAIL'}")
    ok &= tasym >= 1e-6

    # 2: a uniform grid against the direct lattice free-space solution
    nu = args.cells[0]
    uni = Grid(nu, refined=False)
    rho = rng.standard_normal(uni.m)
    phi = james(uni, ghost_green(uni, table, radius), rho)
    ijk = np.rint(uni.x / uni.big - 0.5).astype(int)
    dn = ijk[:, None, :] - ijk[None, :, :] + radius
    direct = (table[dn[..., 0], dn[..., 1], dn[..., 2]] / uni.big) @ (uni.vol * rho)
    err = np.max(np.abs(phi - direct)) / np.max(np.abs(direct))
    print(f"2. uniform grid N = {nu}: max |James - direct lattice convolution| / max = {err:.2e} (bound 1e-9): "
          f"{'PASS' if err <= 1e-9 else 'FAIL'}")
    ok &= err <= 1e-9

    # 3: the accuracy on AMR against the exact potential of a Gaussian blob
    print("3. accuracy on AMR, a Gaussian blob (width 0.1) in the refined box:")
    errs = []
    for n in args.cells:
        g = Grid(n, refined=True)
        rho, exact = gaussian(g.x)
        e = float(np.max(np.abs(james(g, ghost_green(g, table, radius), rho) - exact)))
        u = Grid(n, refined=False)
        rho_u, exact_u = gaussian(u.x)
        eu = float(np.max(np.abs(james(u, ghost_green(u, table, radius), rho_u) - exact_u)))
        errs.append(e)
        print(f"   N = {n:3d}: max |phi - phi_exact| = {e:.3e} (uniform coarse grid {eu:.3e})")
    order = math.log(errs[-2] / errs[-1]) / math.log(args.cells[-1] / args.cells[-2])
    print(f"   observed order (finest pair) {order:.2f} (at least 1.8): {'PASS' if order >= 1.8 else 'FAIL'}")
    ok &= order >= 1.8
    print("JAMES PROTOTYPE PASSED" if ok else "JAMES PROTOTYPE FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
