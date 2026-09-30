#!/usr/bin/env python3
"""Generate and check the WENO interpolation tables of `adam_weno_object` (issue #47, M3-P1).

Interpolation (Jiang, Shu & Zhang 2013): face values from POINT values, unlike the library's reconstruction tables
(face values from cell averages). Same layout as the reconstruction tables `a` and `p`:
  p(f, s2, s1, S)  coefficient of cell (s1 - s2) in the candidate polynomial of stencil s1 (s1, s2 = 0 .. S-1);
  a(f, s1, S)      linear (optimal) weight of stencil s1;
f = 1 is the left interface of cell 0 (x = -1/2), f = 2 the right one (x = +1/2), cell centres at the integers.
Candidates: the Lagrange interpolant of degree S-1 on the cells of the stencil; weights: those whose combination is the
Lagrange interpolant of degree 2S-2 on the union of the stencils, cells -(S-1) .. S-1.

Prints the Fortran assignments and checks: each candidate exact to degree S-1, the weighted combination exact to degree
2S-2, weights positive and summing to 1. Usage: weno_interpolation_tables.py [S_MAX] (default 3).
"""

from __future__ import annotations

import sys

import sympy as sp


def lagrange(nodes: list[int], node: int, at: sp.Rational) -> sp.Rational:
    """Lagrange basis polynomial of `node` over `nodes`, evaluated at `at`."""
    out = sp.Rational(1)
    for m in nodes:
        if m != node:
            out *= (at - m) / sp.Rational(node - m)
    return sp.nsimplify(out)


def tables(s: int) -> tuple[dict, dict]:
    """Candidate coefficients p[(f, s2, s1)] and linear weights a[(f, s1)] for S = s."""
    p, a = {}, {}
    for f, at in ((1, sp.Rational(-1, 2)), (2, sp.Rational(1, 2))):
        for s1 in range(s):
            nodes = [s1 - s2 for s2 in range(s)]
            for s2 in range(s):
                p[(f, s2, s1)] = lagrange(nodes, s1 - s2, at)
        union = list(range(-(s - 1), s))
        target = {m: lagrange(union, m, at) for m in union}
        w = sp.symbols(f"w0:{s}")
        eqs = [sum(w[s1] * p[(f, s2, s1)] for s1 in range(s) for s2 in range(s) if s1 - s2 == m) - target[m]
               for m in union]
        sol = sp.solve(eqs, w, dict=True)
        if len(sol) != 1:
            raise SystemExit(f"S = {s}, f = {f}: no unique linear weights")
        for s1 in range(s):
            a[(f, s1)] = sp.nsimplify(sol[0][w[s1]])
    return p, a


def check(s: int, p: dict, a: dict) -> None:
    """Exactness on monomials and weight positivity."""
    for f, at in ((1, sp.Rational(-1, 2)), (2, sp.Rational(1, 2))):
        if sum(a[(f, s1)] for s1 in range(s)) != 1 or any(a[(f, s1)] <= 0 for s1 in range(s)):
            raise SystemExit(f"S = {s}, f = {f}: weights not positive or not summing to 1")
        for deg in range(2 * s - 1):
            comb = sum(a[(f, s1)] * p[(f, s2, s1)] * sp.Rational(s1 - s2) ** deg for s1 in range(s) for s2 in range(s))
            if comb != at**deg:
                raise SystemExit(f"S = {s}, f = {f}: combination not exact at degree {deg}")
            if deg < s:
                for s1 in range(s):
                    cand = sum(p[(f, s2, s1)] * sp.Rational(s1 - s2) ** deg for s2 in range(s))
                    if cand != at**deg:
                        raise SystemExit(f"S = {s}, f = {f}, stencil {s1}: candidate not exact at degree {deg}")


def fortran(v: sp.Rational) -> str:
    """Fortran R8P literal of a rational."""
    num, den = v.p, v.q
    return f"{num}._R8P" if den == 1 else f"{num}._R8P/{den}._R8P"


def main() -> int:
    """Print the tables of S = 1 .. S_MAX."""
    s_max = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    for s in range(1, s_max + 1):
        p, a = tables(s)
        check(s, p, a)
        print(f"! S = {s}: checked (candidates exact to degree {s - 1}, combination to degree {2 * s - 2})")
        for f in (1, 2):
            print("     " + " ; ".join(f"a({f},{s1},S) = {fortran(a[(f, s1)])}" for s1 in range(s)))
        for f in (1, 2):
            for s1 in range(s):
                print("     " + " ; ".join(f"p({f},{s2},{s1},S) = {fortran(p[(f, s2, s1)])}" for s2 in range(s)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
