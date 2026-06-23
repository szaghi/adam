#!/usr/bin/env python3
"""Exact-rational generator and validator for ADAM FDV coefficient tables.

This is the offline source of truth for the finite-difference (FD) and
finite-volume (FV) coefficient tables hard-coded in
``src/lib/common/adam_fdv_operators_library.F90``.

All coefficients are derived in *exact rational arithmetic* (``sympy.Rational``)
from first principles, then validated bit-for-bit against the committed Fortran
tables. Phase 1 of the FDV SOTA rebuild is "validation-first": prove this script
reproduces the existing hand-entered tables before extending to higher orders.

Conventions (matching the Fortran module):
- ``s`` is the *half stencil length*; classical order of accuracy is ``2*s`` for
  the centered families.
- Tables store only half the stencil, exploiting symmetry/antisymmetry:
    * FD derivative ``n`` centered: weight ``c_m`` with
      ``d^n q/ds^n ~= (1/ds^n) * sum_{m=1..s} FDn_CC(m,s) * (q(1+m) +/- q(1-m))``
      (``+`` for even ``n`` plus a central ``c_0`` term, ``-`` for odd ``n``).
    * FV centered reconstruction at the right face ``i+1/2``:
      ``q_{i+1/2} ~= sum_{m=1..s} FV1_CC(m,s) * (q(m) + q(1-m))``  (symmetric).
    * FV right-upwind reconstruction:
      ``q_{i+1/2} ~= sum_{m=1..s} FV1_UR(m,s) * q(m)``  (one-sided, nodes 0..s-1
      i.e. cells i..i+s-1; stored shifted so index 1 <-> node 0).

Run ``python3 scripts/fdv_coefficients.py`` to validate. Non-zero exit on any
mismatch.
"""

from __future__ import annotations

import sys

import sympy as sp

S_MAX = 5  # current maximum half-stencil length in the Fortran module


def fd_weights(deriv: int, nodes: list[int]) -> list[sp.Rational]:
    """Exact finite-difference weights for ``deriv``-th derivative at x=0.

    Solves the Taylor-matching linear system exactly: find weights ``w_j`` on the
    integer ``nodes`` (in units of the grid spacing) such that
    ``sum_j w_j f(node_j) = f^{(deriv)}(0) + O(h^p)``. Exact rational solve; no
    conditioning concern at these sizes.
    """
    n = len(nodes)
    # Vandermonde-style moment matrix M[k, j] = node_j^k / k! is unnecessary;
    # match raw moments: sum_j w_j node_j^k = deriv! if k==deriv else 0.
    matrix = sp.Matrix([[sp.Integer(node) ** k for node in nodes] for k in range(n)])
    rhs = sp.Matrix([sp.factorial(deriv) if k == deriv else 0 for k in range(n)])
    weights = matrix.solve(rhs)
    return [sp.nsimplify(w) for w in weights]


def _stencil_too_small(deriv: int, s: int) -> bool:
    """A centered stencil of half-width ``s`` (2s+1 points) cannot represent a
    derivative of order ``deriv`` if ``deriv > 2s`` (not enough points). The
    Fortran tables store all-zero rows in those cases (e.g. FD3 for s=1, FD5 for
    s<=2). Mirror that convention exactly."""
    return deriv > 2 * s


def fd_centered_halfstencil(deriv: int, s: int) -> list[sp.Rational]:
    """Half-stencil FD centered weights for the ODD-derivative layout.

    Returns ``[c_1, ..., c_s]`` applied to the antisymmetric pairs
    ``q(1+m) - q(1-m)`` (for odd ``deriv``) or symmetric pairs ``q(1+m) + q(1-m)``
    (for ``deriv``=0). Zero row when the stencil is too small.
    """
    if _stencil_too_small(deriv, s):
        return [sp.Integer(0)] * s
    nodes = list(range(-s, s + 1))  # -s .. s about the center
    w = fd_weights(deriv, nodes)
    wmap = {node: wj for node, wj in zip(nodes, w, strict=True)}
    return [wmap[m] for m in range(1, s + 1)]


def fd_centered_even_layout(deriv: int, s: int) -> list[sp.Rational]:
    """Even-derivative FD centered layout: ``[c_0, c_1, ..., c_s]`` (length s+1),
    where ``c_0`` is the central weight and ``c_m`` multiplies the symmetric pair
    ``q(1+m) + q(1-m)``. Zero row when the stencil is too small."""
    if _stencil_too_small(deriv, s):
        return [sp.Integer(0)] * (s + 1)
    nodes = list(range(-s, s + 1))
    w = fd_weights(deriv, nodes)
    wmap = {node: wj for node, wj in zip(nodes, w, strict=True)}
    return [wmap[0]] + [wmap[m] for m in range(1, s + 1)]


def fd0_reconstruction_halfstencil(s: int) -> list[sp.Rational]:
    """Half-stencil FD reconstruction (deriv 0) of the face value q_{i+1/2} from
    POINT values (Lagrange), symmetric. Returns ``[a_1..a_s]`` on pairs
    ``q(m) + q(1-m)``. Distinct from the FV reconstruction (cell averages):
    here the data are point samples, matched against the point value at +1/2."""
    nodes = list(range(1 - s, s + 1))  # point nodes straddling the face
    n = len(nodes)
    half = sp.Rational(1, 2)
    # point value of monomial x^k at integer node j is just j^k
    matrix = sp.Matrix([[sp.Integer(node) ** k for node in nodes] for k in range(n)])
    rhs = sp.Matrix([half**k for k in range(n)])
    a = matrix.solve(rhs)
    amap = {node: aj for node, aj in zip(nodes, a, strict=True)}
    return [amap[m] for m in range(1, s + 1)]


def fv_lupwind_deriv_halfstencil(m: int, s: int) -> list[sp.Rational]:
    """Left-upwind reconstruction of the m-th derivative at the face i+1/2 from the
    one-sided cell averages on offsets {0, +1, ..., +(s-1)} (cells at and ahead of
    the face) -- the mirror of fv_rupwind_deriv. m=0 is FV1_UL. Returns weights in
    Fortran storage order matching FV1_UL (applied as a_k * q(1-k), k=1..s, i.e.
    node order [0, -? ...]); we return [a_node0, a_node+1, ...] then the caller maps.
    Stored to match the existing FV1_UL row (reverse of FV1_UR for m=0)."""
    nodes = list(range(s, 0, -1))  # cells i+s, ..., i+2, i+1 (descending), matches FV1_UL storage
    n = len(nodes)
    half = sp.Rational(1, 2)

    def avg(node: int, k: int) -> sp.Rational:
        j = sp.Integer(node)
        return ((j + half) ** (k + 1) - (j - half) ** (k + 1)) / (k + 1)

    def dmonom_at_half(k: int) -> sp.Rational:
        if k < m:
            return sp.Integer(0)
        return sp.factorial(k) / sp.factorial(k - m) * half ** (k - m)

    matrix = sp.Matrix([[avg(node, k) for node in nodes] for k in range(n)])
    rhs = sp.Matrix([dmonom_at_half(k) for k in range(n)])
    a = matrix.solve(rhs)
    return list(a)


def fv_lupwind_halfstencil(s: int) -> list[sp.Rational]:
    """Left-upwind FV reconstruction: the mirror of the right-upwind family. The
    Fortran storage is the reversed right-upwind row. We derive it directly by the
    mirrored node set: cells i+1, i+2, ..., i+s reconstructing the face i+1/2 from
    the downwind side, stored as ``[a_{1-s}, ..., a_0]`` (reversed UR order)."""
    ur = fv_rupwind_halfstencil(s)
    return list(reversed(ur))


def fv_centered_face_halfstencil(s: int) -> list[sp.Rational]:
    """Half-stencil FV centered reconstruction weights at the face i+1/2.

    Primitive-function construction: the face value is the derivative of the
    interpolated primitive F (whose face point-values are exact running sums of
    cell averages). Equivalently, the reconstruction weights are the FD weights
    that map cell averages on nodes (-(s-1)..s) to the point value at +1/2,
    obtained by differentiating the degree-2s primitive interpolant once.

    We compute it directly: find weights a_j on cell-center offsets j = 1-s..s
    such that sum_j a_j <q>_{i+j} = q(i+1/2) + O(h^{2s}), where <q> is the cell
    average. On a uniform grid the cell-average of a monomial x^k over the cell
    centered at j is (1/h) * integral_{j-1/2}^{j+1/2} x^k dx (h=1). Match these
    moments against the point value (1/2)^k at the face.
    """
    nodes = list(range(1 - s, s + 1))  # 2s cells straddling the face i+1/2
    n = len(nodes)
    # cell-average of monomial degree k over unit cell centered at node j:
    #   avg_k(j) = ((j+1/2)^{k+1} - (j-1/2)^{k+1}) / (k+1)
    half = sp.Rational(1, 2)

    def avg(node: int, k: int) -> sp.Rational:
        j = sp.Integer(node)
        return ((j + half) ** (k + 1) - (j - half) ** (k + 1)) / (k + 1)

    matrix = sp.Matrix([[avg(node, k) for node in nodes] for k in range(n)])
    # target: point value of x^k at the face x = +1/2
    rhs = sp.Matrix([half**k for k in range(n)])
    a = matrix.solve(rhs)
    amap = {node: aj for node, aj in zip(nodes, a, strict=True)}
    # fold symmetric pair: weight on (q(m) + q(1-m)) is a[m] (== a[1-m]); m=1..s
    # node m corresponds to offset m, node 1-m to offset 1-m.
    return [amap[m] for m in range(1, s + 1)]


def fv_face_deriv_weights(m: int, s: int) -> dict[int, sp.Rational]:
    """Full node->weight map for reconstructing the m-th derivative of q at the
    face x=+1/2 from CELL AVERAGES on the 2s cells straddling the face.

    This is the primitive-function construction generalized to derivative m:
    interpolate the primitive F=int q (degree-2s through 2s face running-sums),
    differentiate (m+1) times, evaluate at the face. Equivalently, solve the
    moment system on cell averages targeting d^m/dx^m[x^k] at x=1/2.

    For m=0 this is exactly the face value reconstruction (FV1_CC). The FV
    derivative-n ladder is built as a flux difference of the (n-1)th face
    reconstruction:  dn|_i = (g_{i+1/2} - g_{i-1/2})/h  with g = m=(n-1) recon.
    """
    nodes = list(range(1 - s, s + 1))  # 2s cells straddling the face i+1/2
    n = len(nodes)
    half = sp.Rational(1, 2)

    def avg(node: int, k: int) -> sp.Rational:
        j = sp.Integer(node)
        return ((j + half) ** (k + 1) - (j - half) ** (k + 1)) / (k + 1)

    def dmonom_at_half(k: int) -> sp.Rational:
        # d^m/dx^m (x^k) at x = 1/2  =  k!/(k-m)! * (1/2)^(k-m)  for k>=m else 0
        if k < m:
            return sp.Integer(0)
        return sp.factorial(k) / sp.factorial(k - m) * half ** (k - m)

    matrix = sp.Matrix([[avg(node, k) for node in nodes] for k in range(n)])
    rhs = sp.Matrix([dmonom_at_half(k) for k in range(n)])
    w = matrix.solve(rhs)
    return {node: wj for node, wj in zip(nodes, w, strict=True)}


def fv_face_deriv_halfstencil(m: int, s: int) -> list[sp.Rational]:
    """Half-stencil form of the m-th-derivative face reconstruction.

    Parity about the face i+1/2: even m -> symmetric (w[m_node] == w[1-m_node]),
    odd m -> antisymmetric (w[m_node] == -w[1-m_node]). Returns the positive-side
    weights [a_1..a_s] applied to the (anti)symmetric pairs q(m_)+/-q(1-m_).
    """
    wmap = fv_face_deriv_weights(m, s)
    return [wmap[j] for j in range(1, s + 1)]


def deconvolution_avg_to_point_weights(s: int) -> dict[int, sp.Rational]:
    """Weights c_j on cell averages <q>_{i+j}, j=-s..s, recovering the cell-center
    POINT value q_i to order 2s+1:  q_i = sum_j c_j <q>_{i+j} + O(h^{2s+2}).

    Symmetric (even) operator. For s=1 this is the classic
    q_i = <q>_i - (1/24)(<q>_{i+1} - 2<q>_i + <q>_{i-1}) + O(h^4)
    i.e. c = {-1:-1/24, 0:13/12, +1:-1/24} (McCorquodale-Colella / Mignone Eq.52).

    Moment match: sum_j c_j avg_k(j) = (point value of x^k at x=0) = delta_{k,0},
    for k=0..2s.
    """
    nodes = list(range(-s, s + 1))
    n = len(nodes)
    half = sp.Rational(1, 2)

    def avg(node: int, k: int) -> sp.Rational:
        j = sp.Integer(node)
        return ((j + half) ** (k + 1) - (j - half) ** (k + 1)) / (k + 1)

    matrix = sp.Matrix([[avg(node, k) for node in nodes] for k in range(n)])
    rhs = sp.Matrix([sp.Integer(1) if k == 0 else sp.Integer(0) for k in range(n)])
    c = matrix.solve(rhs)
    return {node: cj for node, cj in zip(nodes, c, strict=True)}


def deconvolution_point_to_avg_weights(s: int) -> dict[int, sp.Rational]:
    """Inverse: weights recovering the cell AVERAGE <q>_i from POINT values q_{i+j}
    to order 2s+1.  <q>_i = sum_j c_j q_{i+j} + O(h^{2s+2}).

    For s=1: <q>_i = q_i + (1/24)(q_{i+1} - 2 q_i + q_{i-1}) + O(h^4).

    Moment match: sum_j c_j (j^k) = avg_k(0) (cell-average of x^k over cell i),
    for k=0..2s.
    """
    nodes = list(range(-s, s + 1))
    n = len(nodes)
    half = sp.Rational(1, 2)

    def avg0(k: int) -> sp.Rational:
        return (half ** (k + 1) - (-half) ** (k + 1)) / (k + 1)

    matrix = sp.Matrix([[sp.Integer(node) ** k for node in nodes] for k in range(n)])
    rhs = sp.Matrix([avg0(k) for k in range(n)])
    c = matrix.solve(rhs)
    return {node: cj for node, cj in zip(nodes, c, strict=True)}


def fv_rupwind_halfstencil(s: int) -> list[sp.Rational]:
    """Half-stencil FV right-upwind reconstruction weights at the face i+1/2.

    One-sided *upwind* reconstruction of the downwind face i+1/2 from the s cells
    at and behind the face: cells i, i-1, ..., i-(s-1) (offsets 0, -1, ..., -(s-1)).
    This is the classical upwind-biased FV reconstruction (e.g. for a right-moving
    wave the i+1/2 flux uses upwind-side cells <= i). Stored as [a_0, a_-1, ...].
    """
    return fv_rupwind_deriv_halfstencil(0, s)


def fv_rupwind_deriv_halfstencil(m: int, s: int) -> list[sp.Rational]:
    """Right-upwind reconstruction of the m-th derivative at the face i+1/2 from the
    one-sided cell averages on offsets {0, -1, ..., -(s-1)} (cells at and behind the
    face). m=0 is FV1_UR. Used for the upwind derivative ladder exactly as the
    centered FVRm tables: dn|_i = (g_{i+1/2} - g_{i-1/2})/ds^n, g = m=(n-1) recon.
    Returns [a_0, a_-1, ..., a_-(s-1)] in node order (no symmetry to exploit)."""
    nodes = list(range(0, -s, -1))  # cells i, i-1, ..., i-(s-1)
    n = len(nodes)
    half = sp.Rational(1, 2)

    def avg(node: int, k: int) -> sp.Rational:
        j = sp.Integer(node)
        return ((j + half) ** (k + 1) - (j - half) ** (k + 1)) / (k + 1)

    def dmonom_at_half(k: int) -> sp.Rational:
        if k < m:
            return sp.Integer(0)
        return sp.factorial(k) / sp.factorial(k - m) * half ** (k - m)

    matrix = sp.Matrix([[avg(node, k) for node in nodes] for k in range(n)])
    rhs = sp.Matrix([dmonom_at_half(k) for k in range(n)])
    a = matrix.solve(rhs)
    return list(a)


# committed Fortran tables (numerator list, common divisor) for validation.
# Half-stencil layout (length s after slicing) unless marked even-layout (length s+1).
COMMITTED = {
    "FD0_CC": {  # deriv-0 reconstruction from point values, symmetric, length s
        1: ([1, 0, 0, 0, 0], 2),
        2: ([9, -1, 0, 0, 0], 16),
        3: ([150, -25, 3, 0, 0], 256),
        4: ([1225, -245, 49, -5, 0], 2048),
        5: ([39690, -8820, 2268, -405, 35], 65536),
    },
    "FD1_CC": {  # deriv 1, antisymmetric, length s
        1: ([1, 0, 0, 0, 0], 2),
        2: ([8, -1, 0, 0, 0], 12),
        3: ([45, -9, 1, 0, 0], 60),
        4: ([672, -168, 32, -3, 0], 840),
        5: ([2100, -600, 150, -25, 2], 2520),
    },
    "FD2_CC": {  # deriv 2, even-layout [c0, c1..cs], length s+1
        1: ([-2, 1, 0, 0, 0, 0], 1),
        2: ([-30, 16, -1, 0, 0, 0], 12),
        3: ([-490, 270, -27, 2, 0, 0], 180),
        4: ([-14350, 8064, -1008, 128, -9, 0], 5040),
        5: ([-73766, 42000, -6000, 1000, -125, 8], 25200),
    },
    "FD3_CC": {  # deriv 3, antisymmetric, length s; zero rows where too small
        1: ([0, 0, 0, 0, 0], 1),
        2: ([-2, 1, 0, 0, 0], 2),
        3: ([-13, 8, -1, 0, 0], 8),
        4: ([-488, 338, -72, 7, 0], 240),
        5: ([-70098, 52428, -14607, 2522, -205], 30240),
    },
    "FD4_CC": {  # deriv 4, even-layout, length s+1
        1: ([0, 0, 0, 0, 0, 0], 1),
        2: ([6, -4, 1, 0, 0, 0], 1),
        3: ([56, -39, 12, -1, 0, 0], 6),
        4: ([2730, -1952, 676, -96, 7, 0], 240),
        5: ([192654, -140196, 52428, -9738, 1261, -82], 15120),
    },
    "FD5_CC": {  # deriv 5, antisymmetric, length s; zero rows where too small
        1: ([0, 0, 0, 0, 0], 1),
        2: ([0, 0, 0, 0, 0], 1),
        3: ([5, -4, 1, 0, 0], 2),
        4: ([29, -26, 9, -1, 0], 6),
        5: ([1938, -1872, 783, -152, 13], 288),
    },
    "FD6_CC": {  # deriv 6, even-layout, length s+1; zero rows where too small
        1: ([0, 0, 0, 0, 0, 0], 1),
        2: ([0, 0, 0, 0, 0, 0], 1),
        3: ([-20, 15, -6, 1, 0, 0], 1),
        4: ([-150, 116, -52, 12, -1, 0], 4),
        5: ([-12276, 9690, -4680, 1305, -190, 13], 240),
    },
    "FV1_CC": {  # FV centered reconstruction from cell averages, length s
        1: ([1, 0, 0, 0, 0], 2),
        2: ([7, -1, 0, 0, 0], 12),
        3: ([37, -8, 1, 0, 0], 60),
        4: ([533, -139, 29, -3, 0], 840),
        5: ([1627, -473, 127, -23, 2], 2520),
    },
    "FV1_UR": {  # FV right-upwind reconstruction, length s
        1: ([1, 0, 0, 0, 0], 1),
        2: ([3, -1, 0, 0, 0], 2),
        3: ([11, -7, 2, 0, 0], 6),
        4: ([25, -23, 13, -3, 0], 12),
        5: ([137, -163, 137, -63, 12], 60),
    },
    "FV1_UL": {  # FV left-upwind reconstruction (reversed UR), length s
        1: ([1, 0, 0, 0, 0], 1),
        2: ([-1, 3, 0, 0, 0], 2),
        3: ([2, -7, 11, 0, 0], 6),
        4: ([-3, 13, -23, 25, 0], 12),
        5: ([12, -63, 137, -163, 137], 60),
    },
}

# tables stored with the even-derivative [c0, c1..cs] layout (length s+1)
EVEN_LAYOUT = {"FD2_CC", "FD4_CC", "FD6_CC"}


def committed_row(table: str, s: int) -> list[sp.Rational]:
    nums, div = COMMITTED[table][s]
    length = s + 1 if table in EVEN_LAYOUT else s
    return [sp.Rational(num, div) for num in nums[:length]]


def validate() -> int:
    failures = 0
    checks = [
        ("FD0_CC", fd0_reconstruction_halfstencil),
        ("FD1_CC", lambda s: fd_centered_halfstencil(1, s)),
        ("FD2_CC", lambda s: fd_centered_even_layout(2, s)),
        ("FD3_CC", lambda s: fd_centered_halfstencil(3, s)),
        ("FD4_CC", lambda s: fd_centered_even_layout(4, s)),
        ("FD5_CC", lambda s: fd_centered_halfstencil(5, s)),
        ("FD6_CC", lambda s: fd_centered_even_layout(6, s)),
        ("FV1_CC", fv_centered_face_halfstencil),
        ("FV1_UR", fv_rupwind_halfstencil),
        ("FV1_UL", fv_lupwind_halfstencil),
    ]
    # Phase 2/3 families (no committed counterpart): assert the two pinned-down
    # invariants instead of bit-for-bit. (a) m=0 face-derivative recon == FV1_CC;
    # (b) s=1 avg->point deconvolution == canonical {-1/24, 13/12, -1/24}.
    for s in range(1, S_MAX + 1):
        got = fv_face_deriv_halfstencil(0, s)
        want = committed_row("FV1_CC", s)
        ok = got == want
        if not ok:
            failures += 1
        print(f"[{'OK ' if ok else 'FAIL'}] FVface(m=0) s={s} == FV1_CC: {got}")
    canon = deconvolution_avg_to_point_weights(1)
    want_canon = {-1: sp.Rational(-1, 24), 0: sp.Rational(13, 12), 1: sp.Rational(-1, 24)}
    ok = canon == want_canon
    if not ok:
        failures += 1
    print(f"[{'OK ' if ok else 'FAIL'}] deconv(s=1) avg->point == Mignone Eq.52: {dict(canon)}")
    for table, gen in checks:
        for s in range(1, S_MAX + 1):
            got = gen(s)
            want = committed_row(table, s)
            ok = got == want
            status = "OK " if ok else "FAIL"
            if not ok:
                failures += 1
            print(f"[{status}] {table} s={s}: generated={got} committed={want}")
    print()
    if failures:
        print(f"VALIDATION FAILED: {failures} mismatch(es)")
    else:
        print("VALIDATION PASSED: all generated tables match committed Fortran bit-for-bit")
    return 1 if failures else 0


MAX_COL = 132  # free-form Fortran hard line-length limit (house style)


def _common_divisor_row(rats: list[sp.Rational]) -> tuple[list[int], int]:
    """Express a list of Rationals as integer numerators over a single common
    divisor (the house style of the Fortran tables: ``[ints]/divisor``)."""
    dens = [r.q for r in rats]
    lcm = sp.ilcm(*dens) if dens else sp.Integer(1)
    nums = [int(r * lcm) for r in rats]
    return nums, int(lcm)


def _emit_param_row(name: str, dim: str, nums: list[int], div: int) -> list[str]:
    """Emit one ``real(R8P), parameter :: NAME(dim)=[...]/div`` declaration,
    wrapping the bracketed body across continuation lines to stay <= MAX_COL.
    Each element keeps its own line when wrapping, so no single line overflows."""
    elems = [f"{n}._R8P" for n in nums]
    suffix = f"/{div}._R8P" if div != 1 else ""
    head = f"real(R8P), parameter :: {name}({dim})=["
    oneline = f"{head}{','.join(elems)}]{suffix}"
    if len(oneline) <= MAX_COL:
        return [oneline]
    # wrap: one element per continuation line, aligned under the bracket
    indent = " " * len(head)
    out = [f"{head}{elems[0]},   &"]
    for e in elems[1:-1]:
        out.append(f"{indent}{e},   &")
    out.append(f"{indent}{elems[-1]}]{suffix}")
    # guard: even a single wrapped element must fit
    for ln in out:
        if len(ln) > MAX_COL:
            raise ValueError(f"line still exceeds {MAX_COL} cols: {ln!r}")
    return out


def _emit_halfstencil(name: str, rows: list[list[sp.Rational]], doc: str) -> str:
    """Emit a half-stencil parameter table block (length S_MAX rows) + reshape."""
    lines = [f"!< {doc}"]
    for s, rats in enumerate(rows, start=1):
        padded = list(rats) + [sp.Integer(0)] * (S_MAX - len(rats))
        nums, div = _common_divisor_row(padded)
        lines += _emit_param_row(f"{name}_S{s}", "S_MAX", nums, div)
    lines.append(f"real(R8P), parameter :: {name}(S_MAX,S_MAX)=reshape([{name}_S1, &")
    for s in range(2, S_MAX + 1):
        cont = ", &" if s < S_MAX else "],&"
        lines.append(f"                                                     {name}_S{s}{cont}")
    lines.append(f"                                                    [S_MAX,S_MAX]) !< {name} table.")
    accnames = ",".join(f"{name}_S{s}" for s in range(1, S_MAX + 1)) + f",{name}"
    lines.append(f"!$acc declare copyin({accnames})")
    return "\n".join(lines)


def _emit_deconv(name: str, weights_fn, doc: str) -> str:
    """Emit a symmetric deconvolution table: [c0, c1..cs] layout (S_MAX+1)."""
    lines = [f"!< {doc}"]
    for s in range(1, S_MAX + 1):
        w = weights_fn(s)
        rats = [w[0]] + [w[m] for m in range(1, s + 1)] + [sp.Integer(0)] * (S_MAX - s)
        nums, div = _common_divisor_row(rats)
        lines += _emit_param_row(f"{name}_S{s}", "S_MAX+1", nums, div)
    lines.append(f"real(R8P), parameter :: {name}(S_MAX+1,S_MAX)=reshape([{name}_S1, &")
    for s in range(2, S_MAX + 1):
        cont = ", &" if s < S_MAX else "],&"
        lines.append(f"                                                         {name}_S{s}{cont}")
    lines.append(f"                                                        [S_MAX+1,S_MAX]) !< {name} table.")
    accnames = ",".join(f"{name}_S{s}" for s in range(1, S_MAX + 1)) + f",{name}"
    lines.append(f"!$acc declare copyin({accnames})")
    return "\n".join(lines)


def _emit_onesided(name: str, deriv_fn, doc: str) -> str:
    """Emit a one-sided (upwind) full-row table: S_MAX-length rows, no symmetry
    folding (each row holds all s weights directly in node order)."""
    lines = [f"!< {doc}"]
    for s in range(1, S_MAX + 1):
        rats = list(deriv_fn(s)) + [sp.Integer(0)] * (S_MAX - s)
        nums, div = _common_divisor_row(rats)
        lines += _emit_param_row(f"{name}_S{s}", "S_MAX", nums, div)
    lines.append(f"real(R8P), parameter :: {name}(S_MAX,S_MAX)=reshape([{name}_S1, &")
    for s in range(2, S_MAX + 1):
        cont = ", &" if s < S_MAX else "],&"
        lines.append(f"                                                     {name}_S{s}{cont}")
    lines.append(f"                                                    [S_MAX,S_MAX]) !< {name} table.")
    accnames = ",".join(f"{name}_S{s}" for s in range(1, S_MAX + 1)) + f",{name}"
    lines.append(f"!$acc declare copyin({accnames})")
    return "\n".join(lines)


def emit_fortran_upwind() -> None:
    """Emit the upwind-family derivative reconstruction tables (m=1..4).

    FVUR{m}_CC: right-upwind reconstruction of the m-th derivative at face i+1/2 from
    one-sided cells {0,-1,..,-(s-1)} (node order). FVUL{m}_CC: left-upwind, cells
    {s,..,2,1} (descending, matching FV1_UL storage). m=0 are the existing FV1_UR/UL.
    No d5 table (m=4 covers up to d5; d6 would need m=5 = 6-cell stencil, unrepresentable
    at S_MAX=5, so upwind d6 stays on recursion). Used as the centered ladder:
    dn|_i = (g_{i+1/2} - g_{i-1/2})/ds^n, g = FVU?{n-1} reconstruction.
    """
    for m in range(1, S_MAX):  # m=1..4 -> supports upwind d2..d5
        doc = f"FV right-upwind face-reconstruction of derivative {m} (one-sided, cells <= i)."
        print(_emit_onesided(f"FVUR{m}_CC", lambda s, m=m: fv_rupwind_deriv_halfstencil(m, s), doc))
        print()
    for m in range(1, S_MAX):
        doc = f"FV left-upwind face-reconstruction of derivative {m} (one-sided, cells >= i+1)."
        print(_emit_onesided(f"FVUL{m}_CC", lambda s, m=m: fv_lupwind_deriv_halfstencil(m, s), doc))
        print()


def emit_fortran() -> None:
    """Emit the Phase-2 centered-family Fortran coefficient blocks to stdout.

    - FVR{m}_CC: half-stencil face reconstruction of the m-th derivative from cell
      averages, m=1..5 (m=0 is the existing FV1_CC). Parity: even m symmetric,
      odd m antisymmetric about the face i+1/2 (same half-stencil pair layout as
      FV1_CC). Used as: dn|_i = (g_{i+1/2} - g_{i-1/2})/ds, g = FVR{n-1} recon.
    - DECONV_A2P / DECONV_P2A: symmetric cell-average<->point deconvolution.
    """
    for m in range(1, S_MAX + 1):
        parity = "symmetric" if m % 2 == 0 else "antisymmetric"
        doc = f"FV face-reconstruction of derivative {m} (cell averages -> face i+1/2), {parity}."
        rows = [fv_face_deriv_halfstencil(m, s) for s in range(1, S_MAX + 1)]
        print(_emit_halfstencil(f"FVR{m}_CC", rows, doc))
        print()
    print(_emit_deconv("DECONV_A2P", deconvolution_avg_to_point_weights,
                       "Deconvolution cell-average -> cell-center point value (symmetric)."))
    print()
    print(_emit_deconv("DECONV_P2A", deconvolution_point_to_avg_weights,
                       "Deconvolution cell-center point value -> cell-average (symmetric)."))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--emit-fortran":
        emit_fortran()
        sys.exit(0)
    if len(sys.argv) > 1 and sys.argv[1] == "--emit-fortran-upwind":
        emit_fortran_upwind()
        sys.exit(0)
    sys.exit(validate())
