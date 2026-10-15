#!/usr/bin/env python3
"""The composite Laplacian across 2:1 faces: three couplings, measured (issue #81, P5a; the choice behind SV-9).

Self-gravity on AMR needs a composite operator A (finite-volume form, V-weighted: (A phi)_i = sum over faces of
T (phi_nbr - phi_i)) that is
  - **symmetric**, so that phi = K f has V K symmetric: the exact E + W exchange of #81 P4;
  - **consistent**, so that the truncation error tau = A phi_exact / V - Lap phi_exact stays bounded at the 2:1 faces
    and the face gradient (the gravity the momentum flux uses there) converges.
The 2-D model: the doubly periodic unit square, coarse cells (H = 1/N) on x < 1/2, fine cells (h = H/2) on x > 1/2, so
two 2:1 seams (x = 1/2 and x = 0 = 1). A coarse cell C faces two fine cells f0, f1 across a seam; the couplings:
  B   Tomida & Stone (2023) flux matching: the fine ghost phi_g = phi_f / 3 + 2/3 (phi_C +- (h/2) D), D the central
      tangential difference of the coarse cells, the coarse face flux the sum of the fine ones. Consistent, not
      symmetric.
  S1  (B + B^T) / 2. Symmetric; the tangential D terms of B land, transposed, in the coarse rows as O(1/h)
      truncation: inconsistent.
  M   the mimetic mean flux (adam_poisson_mg_object): every fine face flux (2/3)(mean(phi_f0, phi_f1) - phi_C), the
      coarse cell their sum. Symmetric and consistent (exact for linear phi face by face).
Measured per N: max |phi - phi_exact| (both mean-free), max |tau|, the seam face gradient error on the fine side (each
operator's own: its seam flux over the face area), and the asymmetry ||A - A^T|| / ||A||.

Usage: seam_operators.py [N ...]   (default 8 16 32 48; dense solves, numpy only)
"""

from __future__ import annotations

import math
import sys

import numpy as np

KINDS = ("B", "S1", "M")
TWO_PI = 2.0 * math.pi


def exact(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Periodic exact potential."""
    return np.sin(TWO_PI * x) * np.cos(TWO_PI * y) + 0.5 * np.cos(TWO_PI * x + 2.0 * TWO_PI * y)


def lap_exact(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Laplacian of the exact potential."""
    return (-2.0 * TWO_PI**2 * np.sin(TWO_PI * x) * np.cos(TWO_PI * y)
            - 0.5 * 5.0 * TWO_PI**2 * np.cos(TWO_PI * x + 2.0 * TWO_PI * y))


def grad_x_exact(x: float, y: float) -> float:
    """x-derivative of the exact potential."""
    return (TWO_PI * math.cos(TWO_PI * x) * math.cos(TWO_PI * y)
            - 0.5 * TWO_PI * math.sin(TWO_PI * x + 2.0 * TWO_PI * y))


def build(n: int, kind: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, list]:
    """Assemble the V-weighted composite operator.

    Parameters
    ----------
    n : int
        Coarse cells per direction (even).
    kind : str
        Seam coupling, one of KINDS.

    Returns
    -------
    tuple
        The operator A, cell volumes, cell centres x and y, and the seam fine cells as (fine index, side, seam
        couplings), side +1 when the coarse cell lies at -x of the fine one, the couplings a {column: weight} map of
        the seam part of the fine cell's row.
    """
    big, small = 1.0 / n, 0.5 / n
    idx: dict = {}
    xs, ys, vol = [], [], []
    for i in range(n // 2):
        for j in range(n):
            idx[("c", i, j)] = len(xs)
            xs.append((i + 0.5) * big)
            ys.append((j + 0.5) * big)
            vol.append(big * big)
    for i in range(n):
        for j in range(2 * n):
            idx[("f", i, j)] = len(xs)
            xs.append(0.5 + (i + 0.5) * small)
            ys.append((j + 0.5) * small)
            vol.append(small * small)
    m = len(xs)
    a = np.zeros((m, m))

    def couple(p: int, q: int) -> None:
        a[p, q] += 1.0
        a[p, p] -= 1.0
        a[q, p] += 1.0
        a[q, q] -= 1.0

    for i in range(n // 2):  # coarse-coarse
        for j in range(n):
            couple(idx[("c", i, j)], idx[("c", i, (j + 1) % n)])
            if i < n // 2 - 1:
                couple(idx[("c", i, j)], idx[("c", i + 1, j)])
    for i in range(n):  # fine-fine
        for j in range(2 * n):
            couple(idx[("f", i, j)], idx[("f", i, (j + 1) % (2 * n))])
            if i < n - 1:
                couple(idx[("f", i, j)], idx[("f", i + 1, j)])
    b = np.zeros((m, m))  # the seam part, B or M
    seam = []
    for ic, ifn, side in ((n // 2 - 1, 0, 1), (0, n - 1, -1)):  # x = 1/2 (coarse at -x), x = 0 = 1 (coarse at +x)
        for jc in range(n):
            pc = idx[("c", ic, jc)]
            fine = [idx[("f", ifn, 2 * jc + t)] for t in (0, 1)]
            for t, pf in enumerate(fine):
                if kind == "M":
                    w = {pc: 2.0 / 3.0}
                    for q in fine:
                        w[q] = w.get(q, 0.0) - 1.0 / 3.0
                    for q, cq in w.items():
                        b[pf, q] += cq
                    b[pc, pf] += 2.0 / 3.0
                    b[pc, pc] -= 2.0 / 3.0
                else:
                    sgn = -1.0 if t == 0 else 1.0
                    up, dn = idx[("c", ic, (jc + 1) % n)], idx[("c", ic, (jc - 1) % n)]
                    c = 2.0 / 3.0 * sgn * (small / 2) / (2 * big)
                    w = {pc: 2.0 / 3.0, up: c, dn: -c, pf: -2.0 / 3.0}
                    for q, cq in w.items():
                        b[pf, q] += cq
                        if q != pf:
                            b[pc, q] -= cq
                    b[pc, pf] += 2.0 / 3.0
                seam.append((pf, side))
    if kind == "S1":
        b = 0.5 * (b + b.T)
    seam = [(pf, side, {int(q): b[pf, q] for q in np.nonzero(b[pf])[0]}) for pf, side in seam]
    return a + b, np.array(vol), np.array(xs), np.array(ys), seam


def measure(n: int, kind: str) -> tuple[float, float, float, float]:
    """Errors of one coupling at one N: phi, truncation, seam face gradient, asymmetry."""
    a, vol, xs, ys, seam = build(n, kind)
    f = vol * lap_exact(xs, ys)
    f -= vol * (f.sum() / vol.sum())
    # A is singular (the constants): pin the mean with a rank-one term, then compare V-weighted mean-free fields
    phi = np.linalg.solve(a - np.outer(vol, vol) / vol.sum(), f)
    pe = exact(xs, ys)
    phi -= (vol @ phi) / vol.sum()
    pe -= (vol @ pe) / vol.sum()
    e_phi = float(np.max(np.abs(phi - pe)))
    tau = float(np.max(np.abs(a @ exact(xs, ys) / vol - lap_exact(xs, ys))))
    small = 0.5 / n
    e_face = 0.0
    for pf, side, w in seam:  # seam flux = area * outward gradient = small * (-side) * dphi/dx
        flux = sum(cq * phi[q] for q, cq in w.items())
        g = -side * flux / small
        e_face = max(e_face, abs(g - grad_x_exact(xs[pf] - side * small / 2, ys[pf])))
    asym = float(np.linalg.norm(a - a.T) / np.linalg.norm(a))
    return e_phi, tau, e_face, asym


def ladder(ns: list[int]) -> dict[str, list[tuple[int, float, float, float, float]]]:
    """Every coupling over the N ladder: rows (N, phi error, truncation, seam face gradient error, asymmetry)."""
    return {k: [(n, *measure(n, k)) for n in ns] for k in KINDS}


def main() -> None:
    """Print the ladder and the orders."""
    ns = [int(v) for v in sys.argv[1:]] or [8, 16, 32, 48]
    for kind, rows in ladder(ns).items():
        print(f"coupling {kind}")
        prev = None
        for n, e_phi, tau, e_face, asym in rows:
            line = (f"   N = {n:3d}: |phi| {e_phi:.3e}, |tau| {tau:.3e}, seam face gradient {e_face:.3e}, "
                    f"asymmetry {asym:.1e}")
            if prev:
                r = math.log(n / prev[0])
                line += "  orders " + ", ".join(f"{math.log(p / q) / r:+.2f}"
                                                for p, q in zip(prev[1:4], (e_phi, tau, e_face), strict=True))
            print(line)
            prev = (n, e_phi, tau, e_face)


if __name__ == "__main__":
    main()
