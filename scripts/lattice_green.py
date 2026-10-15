#!/usr/bin/env python3
"""Generator and validator of the lattice Green's function of the 7-point Laplacian (issue #81, M8-P4).

The isolated FFT Poisson solver (``src/lib/common/adam_poisson_fft_object.F90``) convolves the source with the
**lattice** Green's function G of the second-order 7-point Laplacian on the unit cubic lattice,

    L_1 G(n) = sum_d [G(n + e_d) - 2 G(n) + G(n - e_d)] = delta_{n,0},   G(n) -> -1 / (4 pi |n|) as |n| -> oo,

so that the discrete Poisson equation holds inside the domain (Hockney and Eastwood's zero padding with the lattice
kernel instead of the sampled continuum one). On a lattice of spacing h the kernel is h^2 G(n).

G is computed here once and emitted as the Fortran module ``src/lib/common/adam_lattice_green_library.F90``:

- **Far field, max |n_d| > R:** the asymptotic expansion G0 + G1 + G2 + G3 of the Taylor series of L_1 = Delta +
  (1/12) sum_d d^4/dx_d^4 + (1/360) sum_d d^6/dx_d^6 + (2/8!) sum_d d^8/dx_d^8 + ...: G0 = -1 / (4 pi r),
  G1 = -5 (n1^4 + n2^4 + n3^4 - 3 r^4 / 5) / (32 pi r^7) (Delta G1 = -(1/12) sum d^4 G0; unique, since no
  cubic-invariant harmonic of degree 2 exists), G2 = P8(n) / r^13 and G3 = P12(n) / r^19, each solving its order of
  the series (`derive_g2`, `derive_g3`). G2 is unique up to c3 Y4 r^4 / r^13 and G3 up to d6 Y6 r^6 / r^19 (Y4, Y6
  the cubic harmonics of degree 4 and 6, one each): the Taylor series cannot fix c3 and d6, the lattice does, so they
  are fitted to the table (least squares over 12 <= r <= 26, iterated with the table, see `fit_lattice`). The next
  term is O(r^-9).
- **Near field, max |n_d| <= R:** the Dirichlet problem L_1 G = delta inside the cube [-R, R]^3, with G = G0 + G1 +
  G2 on its faces, solved exactly by a type-I discrete sine transform (it diagonalises L_1 with Dirichlet faces).
  The table therefore satisfies L_1 G = delta to round-off inside the cube; the face data's truncation error enters
  the table only as a discrete-harmonic perturbation (the discrete maximum principle bounds it by the face error),
  which the residual does not see. It does show where the table meets the expansion: hence the composite check.

Stored: the wedge n1 >= n2 >= n3 >= 0, n1 <= R (the cubic symmetry gives the rest), index
n1 (n1 + 1) (n1 + 2) / 6 + n2 (n2 + 1) / 2 + n3 + 1.

Checks (exit 1 on a failure):
- the expansion: G1, G2 and G3 against their orders of the Taylor series of L_1, symbolically;
- the residual of L_1 on the table, every interior point, against delta;
- the residual of L_1 on the composite kernel (table and expansion, as the Fortran function evaluates it) on
  [-2R, 2R]^3, across the faces where they meet;
- G(0) against Watson's closed form for the simple cubic lattice (mpmath), which the Dirichlet problem knows nothing
  about, so it measures the face error that reaches the centre;
- G(1, 0, 0) - G(0) = 1/6 (the residual at the origin, written out);
- with ``--check``, the committed module against a fresh emission, byte for byte.

Usage: ``python3 scripts/lattice_green.py [--radius R] [--out PATH | --check]`` (needs numpy, sympy, mpmath).
"""

from __future__ import annotations

import argparse
import itertools
import math
import sys
from pathlib import Path

import mpmath
import numpy as np
import sympy as sp

REPO = Path(__file__).resolve().parent.parent
MODULE = REPO / "src" / "lib" / "common" / "adam_lattice_green_library.F90"
RADIUS = 32  # half side of the tabulated cube
PER_LINE = 3  # values per Fortran line
PER_DATA = 300  # values per DATA statement (100 lines, under the 255 continuation lines of the standard)


def check_expansion() -> bool:
    """Return True when Delta G1 = -(1/12) sum_d d^4 G0 / dx_d^4 for r > 0 (the h^2 term of L_h's Taylor series)."""
    x, y, z = sp.symbols("x y z", positive=True)
    r = sp.sqrt(x**2 + y**2 + z**2)
    g0 = -1 / (4 * sp.pi * r)
    g1 = -sp.Rational(5, 32) / sp.pi * (x**4 + y**4 + z**4 - sp.Rational(3, 5) * r**4) / r**7
    lap = sum(sp.diff(g1, v, 2) for v in (x, y, z))
    rhs = -sp.Rational(1, 12) * sum(sp.diff(g0, v, 4) for v in (x, y, z))
    point = {x: sp.Rational(3, 7), y: sp.Rational(5, 11), z: sp.Rational(2, 13)}
    return sp.simplify((lap - rhs).subs(point)) == 0


def g2_coefficients(c3: float) -> tuple[float, float, float, float]:
    """Return the coefficients of P8 = c0 sum n^8 + c1 sum n^6 n'^2 + c2 sum n^4 n'^4 + c3 n1^2 n2^2 n3^2 r^2
    (`derive_g2` solves for c0, c1, c2 given c3)."""
    return (-c3 / 13 - 71 / (1664 * math.pi), c3 / 13 + 23 / (13 * math.pi), 4 * c3 / 13 - 8985 / (1664 * math.pi), c3)


def derive_g2() -> bool:
    """Return True when the coefficients of `g2_coefficients` solve Delta G2 = -(1/12) sum d^4 G1 - (1/360) sum d^6 G0
    for any c3 (symbolically, at a rational point)."""
    x, y, z, c3 = sp.symbols("x y z c3", positive=True)
    r2 = x**2 + y**2 + z**2
    r = sp.sqrt(r2)
    g0 = -1 / (4 * sp.pi * r)
    g1 = -sp.Rational(5, 32) / sp.pi * (x**4 + y**4 + z**4 - sp.Rational(3, 5) * r2**2) / r**7
    c0, c1, c2 = (-c3 / 13 - sp.Rational(71, 1664) / sp.pi, c3 / 13 + sp.Rational(23, 13) / sp.pi,
                  4 * c3 / 13 - sp.Rational(8985, 1664) / sp.pi)
    p8 = (c0 * (x**8 + y**8 + z**8) + c1 * (x**6 * (y**2 + z**2) + y**6 * (x**2 + z**2) + z**6 * (x**2 + y**2))
          + c2 * (x**4 * y**4 + y**4 * z**4 + z**4 * x**4) + c3 * x**2 * y**2 * z**2 * r2)
    g2 = p8 / r**13

    def d(f: sp.Expr, k: int) -> sp.Expr:
        return sum(sp.diff(f, v, k) for v in (x, y, z))

    defect = d(g2, 2) + sp.Rational(1, 12) * d(g1, 4) + sp.Rational(1, 360) * d(g0, 6)
    point = {x: sp.Rational(3, 7), y: sp.Rational(5, 11), z: sp.Rational(2, 13)}
    return sp.simplify(defect.subs(point)) == 0


G3_PARTS = ((12, 0, 0), (10, 2, 0), (8, 4, 0), (8, 2, 2), (6, 6, 0), (6, 4, 2), (4, 4, 4))  # P12's monomial orbits


def g3_coefficients(c3: float, d6: float) -> tuple[float, ...]:
    """Return the coefficients of P12 on the orbits `G3_PARTS` (`derive_g3` solves for d0..d5 given c3 and d6)."""
    pi = math.pi
    return (d6 / 405 + (-5280 * pi * c3 + 29165) / (11232 * pi),
            -d6 / 90 + (21120 * pi * c3 + 19105) / (3072 * pi),
            -d6 / 15 + (-63360 * pi * c3 - 1815365) / (6656 * pi),
            17 * d6 / 135 + (638880 * pi * c3 + 1630735) / (14976 * pi),
            -43 * d6 / 405 + (-1515360 * pi * c3 + 12498115) / (44928 * pi),
            13 * d6 / 30 + (147840 * pi * c3 + 8735195) / (13312 * pi),
            d6)


def orbit_sum(n: tuple, part: tuple) -> object:
    """Return the sum of n1^a n2^b n3^c over the distinct permutations (a, b, c) of `part` (numpy or sympy)."""
    out = 0
    for e in sorted(set(itertools.permutations(part))):
        out = out + n[0] ** e[0] * n[1] ** e[1] * n[2] ** e[2]
    return out


def derive_g3() -> bool:
    """Return True when `g3_coefficients` solve Delta G3 = -(1/12) sum d^4 G2 - (1/360) sum d^6 G1 - (2/8!) sum d^8 G0
    for any c3, d6 (symbolically, at a rational point)."""
    x, y, z, c3, d6 = sp.symbols("x y z c3 d6")
    r2 = x**2 + y**2 + z**2
    r = sp.sqrt(r2)
    pi = sp.pi
    g0 = -1 / (4 * pi * r)
    g1 = -sp.Rational(5, 32) / pi * (orbit_sum((x, y, z), (4, 0, 0)) - sp.Rational(3, 5) * r2**2) / r**7
    c = (-c3 / 13 - sp.Rational(71, 1664) / pi, c3 / 13 + sp.Rational(23, 13) / pi,
         4 * c3 / 13 - sp.Rational(8985, 1664) / pi, c3)
    p8 = sum(ci * orbit_sum((x, y, z), part) for ci, part in zip(c, ((8, 0, 0), (6, 2, 0), (4, 4, 0), (4, 2, 2)),
                                                               strict=True))
    g2 = p8 / r**13
    d = (d6 / 405 + (-5280 * pi * c3 + 29165) / (11232 * pi), -d6 / 90 + (21120 * pi * c3 + 19105) / (3072 * pi),
         -d6 / 15 + (-63360 * pi * c3 - 1815365) / (6656 * pi),
         sp.Rational(17, 135) * d6 + (638880 * pi * c3 + 1630735) / (14976 * pi),
         -sp.Rational(43, 405) * d6 + (-1515360 * pi * c3 + 12498115) / (44928 * pi),
         sp.Rational(13, 30) * d6 + (147840 * pi * c3 + 8735195) / (13312 * pi), d6)
    g3 = sum(di * orbit_sum((x, y, z), part) for di, part in zip(d, G3_PARTS, strict=True)) / r**19

    def dd(f: sp.Expr, k: int) -> sp.Expr:
        return sum(sp.diff(f, v, k) for v in (x, y, z))

    defect = (dd(g3, 2) + sp.Rational(1, 12) * dd(g2, 4) + sp.Rational(1, 360) * dd(g1, 6)
              + sp.Rational(2, 40320) * dd(g0, 8))
    point = {x: sp.Rational(3, 7), y: sp.Rational(5, 11), z: sp.Rational(2, 13), c3: sp.Rational(4, 7),
             d6: sp.Rational(-3, 5)}
    return sp.simplify(defect.subs(point)) == 0


def g3_term(n1: np.ndarray, n2: np.ndarray, n3: np.ndarray, c3: float, d6: float) -> np.ndarray:
    """Return G3 at the lattice points (n1, n2, n3), none at the origin."""
    coef = g3_coefficients(c3, d6)
    r2 = n1 * n1 + n2 * n2 + n3 * n3
    p12 = sum(ci * orbit_sum((n1, n2, n3), part) for ci, part in zip(coef, G3_PARTS, strict=True))
    return p12 / np.sqrt(r2) ** 19


def g2_term(n1: np.ndarray, n2: np.ndarray, n3: np.ndarray, c3: float) -> np.ndarray:
    """Return G2 at the lattice points (n1, n2, n3), none at the origin."""
    c0, c1, c2, c3 = g2_coefficients(c3)
    r2 = n1 * n1 + n2 * n2 + n3 * n3
    m8 = n1**8 + n2**8 + n3**8
    m62 = n1**6 * (n2**2 + n3**2) + n2**6 * (n1**2 + n3**2) + n3**6 * (n1**2 + n2**2)
    m44 = n1**4 * n2**4 + n2**4 * n3**4 + n3**4 * n1**4
    m422 = n1**2 * n2**2 * n3**2 * r2
    return (c0 * m8 + c1 * m62 + c2 * m44 + c3 * m422) / np.sqrt(r2) ** 13


def asymptotic(n1: np.ndarray, n2: np.ndarray, n3: np.ndarray, c3: float | None = None,
               d6: float | None = None) -> np.ndarray:
    """Return G0 + G1 (+ G2 when `c3` is given, + G3 when `d6` is too) at the lattice points, none at the origin."""
    r2 = n1 * n1 + n2 * n2 + n3 * n3
    r = np.sqrt(r2)
    quartic = n1**4 + n2**4 + n3**4 - 0.6 * r2 * r2
    g = -1.0 / (4.0 * math.pi * r) - 5.0 * quartic / (32.0 * math.pi * r**7)
    if c3 is not None:
        g = g + g2_term(n1, n2, n3, c3)
        if d6 is not None:
            g = g + g3_term(n1, n2, n3, c3, d6)
    return g


def fit_lattice(g: np.ndarray, radius: int, rmin: float = 12.0, rmax: float = 26.0) -> tuple[float, float]:
    """Return the least-squares (c3, d6) of G - (G0 + G1 + G2 + G3)(0, 0) = c3 d(G2 + G3)/dc3 + d6 dG3/dd6 over
    rmin <= r <= rmax in the table (the expansion is linear in both)."""
    n = np.arange(2 * radius + 1, dtype=float) - radius
    n1, n2, n3 = np.meshgrid(n, n, n, indexing="ij")
    rr = np.sqrt(n1 * n1 + n2 * n2 + n3 * n3)
    sel = (rr >= rmin) & (rr <= rmax)
    a, b, c = n1[sel], n2[sel], n3[sel]
    zero = asymptotic(a, b, c, 0.0, 0.0)
    rem = g[sel] - zero
    basis = np.stack([asymptotic(a, b, c, 1.0, 0.0) - zero, asymptotic(a, b, c, 0.0, 1.0) - zero], axis=1)
    coef = np.linalg.lstsq(basis, rem, rcond=None)[0]
    return float(coef[0]), float(coef[1])


def dst1(x: np.ndarray, axis: int) -> np.ndarray:
    """Return the type-I discrete sine transform along `axis`, y_k = sum_j x_j sin(pi j k / (K + 1)), j, k = 1..K,
    from the FFT of the odd extension [0, x, 0, -reversed x] (exact, numpy only)."""
    x = np.moveaxis(x, axis, -1)
    k = x.shape[-1]
    zero = np.zeros(x.shape[:-1] + (1,))
    ext = np.concatenate([zero, x, zero, -x[..., ::-1]], axis=-1)
    y = -np.fft.fft(ext, axis=-1).imag[..., 1:k + 1] / 2.0
    return np.moveaxis(y, -1, axis)


def dst1n(x: np.ndarray) -> np.ndarray:
    """Return the type-I sine transform along every axis (its own inverse up to prod (2 / (K_d + 1)))."""
    for axis in range(x.ndim):
        x = dst1(x, axis)
    return x


def near_field(radius: int, c3: float | None, d6: float | None = None) -> np.ndarray:
    """Return G on [-R, R]^3 (index n + R): the Dirichlet problem with the asymptotic faces, solved by DST-I."""
    m = 2 * radius + 1
    n = np.arange(m, dtype=float) - radius
    n1, n2, n3 = np.meshgrid(n, n, n, indexing="ij")
    face = (np.abs(n1) == radius) | (np.abs(n2) == radius) | (np.abs(n3) == radius)
    g = np.zeros((m, m, m))
    g[face] = asymptotic(n1[face], n2[face], n3[face], c3, d6)
    rhs = np.zeros((m - 2, m - 2, m - 2))
    rhs[radius - 1, radius - 1, radius - 1] = 1.0
    # move the face values to the right-hand side of the interior equations
    rhs[0, :, :] -= g[0, 1:-1, 1:-1]
    rhs[-1, :, :] -= g[-1, 1:-1, 1:-1]
    rhs[:, 0, :] -= g[1:-1, 0, 1:-1]
    rhs[:, -1, :] -= g[1:-1, -1, 1:-1]
    rhs[:, :, 0] -= g[1:-1, 1:-1, 0]
    rhs[:, :, -1] -= g[1:-1, 1:-1, -1]
    k = np.arange(1, m - 1, dtype=float)
    lam1 = 2.0 * np.cos(math.pi * k / (m - 1)) - 2.0
    lam = lam1[:, None, None] + lam1[None, :, None] + lam1[None, None, :]
    g[1:-1, 1:-1, 1:-1] = dst1n(dst1n(rhs) / lam) * (2.0 / (m - 1)) ** 3
    return g


def residual(g: np.ndarray) -> float:
    """Return max |L_1 G - delta| over the interior of the tabulated cube."""
    c = g[1:-1, 1:-1, 1:-1]
    lap = (g[2:, 1:-1, 1:-1] + g[:-2, 1:-1, 1:-1] + g[1:-1, 2:, 1:-1] + g[1:-1, :-2, 1:-1] + g[1:-1, 1:-1, 2:]
           + g[1:-1, 1:-1, :-2] - 6.0 * c)
    centre = (g.shape[0] - 2) // 2
    lap[centre, centre, centre] -= 1.0
    return float(np.max(np.abs(lap)))


def composite(g: np.ndarray, radius: int, extent: int, c3: float, d6: float) -> np.ndarray:
    """Return the kernel as the Fortran function evaluates it on [-E, E]^3: the table for max |n_d| <= R, the
    asymptotic expansion beyond (the origin is always in the table)."""
    m = 2 * extent + 1
    n = np.arange(m, dtype=float) - extent
    n1, n2, n3 = np.meshgrid(n, n, n, indexing="ij")
    out = np.empty((m, m, m))
    inner = (np.abs(n1) <= radius) & (np.abs(n2) <= radius) & (np.abs(n3) <= radius)
    out[~inner] = asymptotic(n1[~inner], n2[~inner], n3[~inner], c3, d6)
    lo, hi = extent - radius, extent + radius + 1
    out[lo:hi, lo:hi, lo:hi] = g
    return out


def watson_g0() -> float:
    """Return G(0) of the simple cubic lattice from Watson's closed form, -W / 6, to double precision."""
    mpmath.mp.dps = 30
    w = mpmath.sqrt(6) / (32 * mpmath.pi**3) * (mpmath.gamma(mpmath.mpf(1) / 24) * mpmath.gamma(mpmath.mpf(5) / 24)
                                               * mpmath.gamma(mpmath.mpf(7) / 24)
                                               * mpmath.gamma(mpmath.mpf(11) / 24))
    return float(-w / 6)


def wedge(g: np.ndarray, radius: int) -> list[float]:
    """Return the wedge n1 >= n2 >= n3 >= 0 of the table, in the storage order."""
    vals = []
    for a in range(radius + 1):
        for b in range(a + 1):
            for c in range(b + 1):
                vals.append(float(g[a + radius, b + radius, c + radius]))
    return vals


def emit(vals: list[float], radius: int, c3: float, d6: float) -> str:
    """Return the Fortran module."""
    count = len(vals)
    k0, k1, k2, k3 = g2_coefficients(c3)
    q = g3_coefficients(c3, d6)
    head = f"""!< ADAM, lattice Green's function of the 7-point Laplacian (generated: scripts/lattice_green.py, do not edit).
module adam_lattice_green_library
!< ADAM, lattice Green's function of the 7-point Laplacian (generated: scripts/lattice_green.py, do not edit).
!<
!< `lattice_green(n)` returns G(n) on the unit cubic lattice, `L_1 G = delta_(n,0)`, `G -> -1 / (4 pi |n|)`: the
!< kernel of the isolated FFT Poisson solver (issue #81, M8-P4; on a lattice of spacing h the kernel is h^2 G).
!< For max |n_d| <= {radius} it is the table (the Dirichlet problem on the cube, round-off residual inside); beyond, the
!< asymptotic expansion G0 + G1 + G2 + G3, G0 = -1 / (4 pi r), G1 = -5 (n1^4 + n2^4 + n3^4 - 3 r^4 / 5) / (32 pi r^7),
!< G2 = P8(n) / r^13, G3 = P12(n) / r^19 (their lattice coefficients fitted to the table), error O(r^-9).
!< Regenerate and validate with `python3 scripts/lattice_green.py` (`--check` compares this file).
use :: penf, only : I4P, R8P

implicit none
private
public :: lattice_green

integer(I4P), parameter :: LGF_RADIUS={radius}_I4P !< Half side of the tabulated cube.
integer(I4P), parameter :: LGF_COUNT={count}_I4P !< Tabulated values (the wedge n1 >= n2 >= n3 >= 0).
real(R8P),    parameter :: PI=acos(-1._R8P) !< Pi.
real(R8P),    parameter :: G2_C(4)=[{k0:+.17E}_R8P, {k1:+.17E}_R8P, &
                                    {k2:+.17E}_R8P, {k3:+.17E}_R8P] !< G2 = P8 / r^13 coefficients.
real(R8P),    parameter :: G3_C(7)=[{q[0]:+.17E}_R8P, {q[1]:+.17E}_R8P, &
                                    {q[2]:+.17E}_R8P, {q[3]:+.17E}_R8P, &
                                    {q[4]:+.17E}_R8P, {q[5]:+.17E}_R8P, &
                                    {q[6]:+.17E}_R8P] !< G3 = P12 / r^19 coefficients.
real(R8P)               :: LGF_TABLE(LGF_COUNT) !< The wedge, index n1 (n1+1) (n1+2) / 6 + n2 (n2+1) / 2 + n3 + 1.
"""
    data = []
    for s in range(0, count, PER_DATA):
        chunk = vals[s:s + PER_DATA]
        lines = []
        for t in range(0, len(chunk), PER_LINE):
            row = ", ".join(f"{v:+.17E}_R8P" for v in chunk[t:t + PER_LINE])
            lines.append("   " + row)
        body = ", &\n".join(lines)
        data.append(f"data LGF_TABLE({s + 1}:{s + len(chunk)}) / &\n{body} /")
    tail = """
contains
   pure function lattice_green(n) result(g)
   !< Return the lattice Green's function of the 7-point Laplacian at the lattice offset `n`.
   integer(I4P), intent(in) :: n(3)    !< Lattice offset.
   real(R8P)                :: g       !< G(n).
   integer(I4P)             :: a, b, c !< |n| sorted, a >= b >= c.
   integer(I4P)             :: t       !< Swap.
   real(R8P)                :: x(3)    !< |n| sorted, real.
   real(R8P)                :: r2, r   !< |n|^2, |n|.
   real(R8P)                :: p8      !< G2's polynomial.
   real(R8P)                :: p12     !< G3's polynomial.
   real(R8P)                :: s2(3)   !< x^2.

   a = abs(n(1)) ; b = abs(n(2)) ; c = abs(n(3))
   if (a < b) then ; t = a ; a = b ; b = t ; endif
   if (b < c) then ; t = b ; b = c ; c = t ; endif
   if (a < b) then ; t = a ; a = b ; b = t ; endif
   if (a <= LGF_RADIUS) then
      g = LGF_TABLE(a * (a + 1) * (a + 2) / 6 + b * (b + 1) / 2 + c + 1)
   else
      x = [real(a, R8P), real(b, R8P), real(c, R8P)]
      r2 = sum(x**2)
      r = sqrt(r2)
      p8 = G2_C(1) * sum(x**8) + &
           G2_C(2) * (x(1)**6 * (x(2)**2 + x(3)**2) + x(2)**6 * (x(1)**2 + x(3)**2) + x(3)**6 * (x(1)**2 + x(2)**2)) + &
           G2_C(3) * (x(1)**4 * x(2)**4 + x(2)**4 * x(3)**4 + x(3)**4 * x(1)**4) + &
           G2_C(4) * x(1)**2 * x(2)**2 * x(3)**2 * r2
      s2 = x**2
      p12 = G3_C(1) * sum(s2**6) + &
            G3_C(2) * (s2(1)**5 * (s2(2) + s2(3)) + s2(2)**5 * (s2(1) + s2(3)) + s2(3)**5 * (s2(1) + s2(2))) + &
            G3_C(3) * (s2(1)**4 * (s2(2)**2 + s2(3)**2) + s2(2)**4 * (s2(1)**2 + s2(3)**2) + &
                       s2(3)**4 * (s2(1)**2 + s2(2)**2)) + &
            G3_C(4) * (s2(1)**4 * s2(2) * s2(3) + s2(2)**4 * s2(1) * s2(3) + s2(3)**4 * s2(1) * s2(2)) + &
            G3_C(5) * (s2(1)**3 * s2(2)**3 + s2(2)**3 * s2(3)**3 + s2(3)**3 * s2(1)**3) + &
            G3_C(6) * s2(1) * s2(2) * s2(3) * (s2(1)**2 * (s2(2) + s2(3)) + s2(2)**2 * (s2(1) + s2(3)) + &
                                               s2(3)**2 * (s2(1) + s2(2))) + &
            G3_C(7) * (s2(1) * s2(2) * s2(3))**2
      g = -1._R8P / (4._R8P * PI * r) - 5._R8P * (sum(x**4) - 0.6_R8P * r2 * r2) / (32._R8P * PI * r**7) + p8 / r**13 + &
          p12 / r**19
   endif
   endfunction lattice_green
endmodule adam_lattice_green_library
"""
    return head + "\n".join(data) + "\n" + tail


def main() -> int:
    """Generate, validate and write (or check) the Fortran module."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--radius", type=int, default=RADIUS, help="half side of the tabulated cube")
    ap.add_argument("--out", type=Path, default=MODULE, help="output module")
    ap.add_argument("--check", action="store_true", help="compare the committed module instead of writing it")
    args = ap.parse_args()
    ok = True
    expansion = check_expansion() and derive_g2() and derive_g3()
    print(f"expansion: G1, G2, G3 against the Taylor series of L_1: {'PASS' if expansion else 'FAIL'}")
    ok = ok and expansion
    # c3, d6: the table with G0 + G1 faces gives a first fit; passes with the full expansion on the faces converge it
    c3, d6 = fit_lattice(near_field(args.radius, None), args.radius)
    for _ in range(3):
        g = near_field(args.radius, c3, d6)
        c3, d6 = fit_lattice(g, args.radius)
    g = near_field(args.radius, c3, d6)
    print(f"c3 = {c3:.16e}, d6 = {d6:.16e} (the lattice coefficients of G2 and G3, fitted over 12 <= r <= 26)")
    res = residual(g)
    print(f"table residual: max |L_1 G - delta| = {res:.2e} (bound 1e-13): {'PASS' if res <= 1e-13 else 'FAIL'}")
    ok = ok and res <= 1e-13
    c = args.radius
    g0, g0_exact = g[c, c, c], watson_g0()
    err0 = abs(g0 - g0_exact) / abs(g0_exact)
    print(f"G(0) = {g0:.16e}, Watson {g0_exact:.16e}, relative {err0:.2e} (bound 1e-12): "
          f"{'PASS' if err0 <= 1e-12 else 'FAIL'}")
    ok = ok and err0 <= 1e-12
    step = g[c + 1, c, c] - g0 - 1.0 / 6.0
    print(f"G(1,0,0) - G(0) - 1/6 = {step:.2e} (bound 1e-14): {'PASS' if abs(step) <= 1e-14 else 'FAIL'}")
    ok = ok and abs(step) <= 1e-14
    ext = 2 * args.radius
    cres = residual(composite(g, args.radius, ext, c3, d6))
    print(f"kernel residual on [-{ext}, {ext}]^3 (table and expansion): max |L_1 G - delta| = {cres:.2e} "
          f"(bound 1e-12): {'PASS' if cres <= 1e-12 else 'FAIL'}")
    ok = ok and cres <= 1e-12
    for rr in (8, 16, 24):
        one = (np.array(float(rr)), np.array(0.0), np.array(0.0))
        d1 = abs(g[c + rr, c, c] - asymptotic(*one))
        d2 = abs(g[c + rr, c, c] - asymptotic(*one, c3))
        d3 = abs(g[c + rr, c, c] - asymptotic(*one, c3, d6))
        print(f"   G({rr},0,0): - (G0 + G1) = {float(d1):.2e}, - (... + G2) = {float(d2):.2e}, - (... + G3) = "
              f"{float(d3):.2e}")
    text = emit(wedge(g, args.radius), args.radius, c3, d6)
    if args.check:
        same = args.out.read_text() == text
        print(f"{args.out}: {'identical' if same else 'DIFFERS from a fresh emission'}")
        ok = ok and same
    elif ok:
        args.out.write_text(text)
        print(f"wrote {args.out} ({len(text.splitlines())} lines)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
