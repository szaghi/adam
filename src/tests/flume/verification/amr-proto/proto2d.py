#!/usr/bin/env python3
"""The 2-D NumPy prototype of #68 Phase B: a composite grid with one 2:1 box, its transfer operators and its reflux.

Why: #68 Appendix A measured, in 1-D, that the mean restriction and the reflux of point-value fluxes each cap the
order of the composite solution at 2, and that a 6-point restriction with a quintic prolongation and no reflux
recovers order 5. A 1-D model has no tangential direction: no tangential averaging in the register (the coarse face
takes the mean of two fine point fluxes, one more O(H^2) term, #68 F1), no tangential footprint of the prolongation,
no seam corners. B1 repeats the measurement in 2-D, where all three are present, and locates the maximum error
(corner, edge or interior) to test whether the seam corners explain the L-infinity order 1.46 of A1's Euler vortex.

Model (B1): u_t + a u_x + b u_y = 0 on the periodic unit square, conservative finite differences of point values with
the 5th-order linear upwind flux (the WENO5 optimal weights) dimension by dimension, 3 ghost layers, one global time
step, classical RK4 with dt = 0.4 h (64 h)^(1/4) / max(a, b) so that the spatial error dominates (Appendix A). The
coarse grid (N cells per side) covers the domain; the box [0.25, 0.75]^2 is refined 2:1 (spacing h = H/2) and the
coarse cells it covers are overwritten by the restriction of the fine ones before every evaluation, so the coarse
stencils near the box read restricted values (FLUME's coarse ghosts). Arms:
  restriction  mean (the 2x2 mean of the fine point values, FLUME's) | lag6 (6x6 tensor Lagrange at the coarse
               centre, its footprint shifted inside the fine box at the box boundary);
  prolongation cubic (4x4 tensor Lagrange, FLUME's tricubic) | quintic (6x6), from the coarse cells outside the box
               only: along an axis on which the ghost lies outside the box the stencil takes the nearest coarse
               centres on that side (shift-inward, FLUME's donor-block footprint), along the other axis the nearest
               ones (centred);
  reflux       on: a coarse face on the box boundary takes the mean of the two fine face fluxes it covers (the 2-D
               register); off: each grid keeps its own flux.

Usage:
    proto2d.py b1 [--cells 32 64 128] [--velocity 1 1] [--assert ARM:NORM:MIN ...]
    --assert lag6-quintic-off:linf:4.5 fails unless that arm's finest-pair order of that norm is at least MIN
"""

from __future__ import annotations

import argparse
import math
import sys

import numpy as np

G = 3  # ghost layers
BOX = (0.25, 0.75)  # the refined box, the same interval on both axes
C5 = np.array([2.0, -13.0, 47.0, 27.0, -3.0]) / 60.0  # 5th-order upwind face value from cells i-2 .. i+2


def lagrange(xs: np.ndarray, x: float) -> np.ndarray:
    """Return the Lagrange weights of the nodes xs at x."""
    w = np.ones(len(xs))
    for j in range(len(xs)):
        for m in range(len(xs)):
            if m != j:
                w[j] *= (x - xs[m]) / (xs[j] - xs[m])
    return w


def face_flux(p: np.ndarray, axis: int) -> np.ndarray:
    """Return the upwind face values (velocity > 0) of an array padded by G along axis: n + 1 faces for n cells."""
    p = np.moveaxis(p, axis, 0)
    n = p.shape[0] - 2 * G
    f = sum(C5[k] * p[k:k + n + 1] for k in range(5))  # face f between padded cells f + 2 and f + 3
    return np.moveaxis(f, 0, axis)


class Composite:
    """A periodic coarse grid with the box refined 2:1, and the transfer operators of one arm."""

    def __init__(self, n: int, restriction: str, prolongation: str, reflux: bool, refined: bool = True) -> None:
        self.n, self.H = n, 1.0 / n
        self.h = self.H / 2.0
        self.refined, self.reflux = refined, reflux
        self.xc = (np.arange(n) + 0.5) * self.H
        self.i0, self.i1 = round(BOX[0] * n), round(BOX[1] * n)  # covered coarse cells i0 .. i1 - 1 (both axes)
        self.nf = 2 * (self.i1 - self.i0)
        self.xf = BOX[0] + (np.arange(self.nf) + 0.5) * self.h
        if refined:
            self._restriction_stencils(restriction)
            self._prolongation_stencils(4 if prolongation == "cubic" else 6)

    def _restriction_stencils(self, kind: str) -> None:
        """Per covered coarse index (one axis): the fine indices and weights (the 2-D operator is their tensor)."""
        idx, wts = [], []
        for i in range(self.i0, self.i1):
            k = 2 * (i - self.i0)  # the fine cells k, k + 1 hold the coarse cell i
            if kind == "mean":
                idx.append([k, k + 1])
                wts.append([0.5, 0.5])
            else:
                j0 = min(max(k - 2, 0), self.nf - 6)
                idx.append(list(range(j0, j0 + 6)))
                wts.append(list(lagrange(self.xf[j0:j0 + 6], self.xc[i])))
        self.r_idx, self.r_wts = np.array(idx), np.array(wts)

    def _axis_stencil(self, x: float, npts: int) -> tuple[np.ndarray, np.ndarray]:
        """Return the coarse indices (periodic) and weights of a 1-D stencil at x: one-sided away from the box when x
        lies outside it on this axis, the nearest centres otherwise."""
        cand = np.arange(-self.n, 2 * self.n)  # periodic images of the coarse centres
        xs = (cand + 0.5) * self.H
        if x < BOX[0]:
            cand, xs = cand[xs < BOX[0]], xs[xs < BOX[0]]
        elif x > BOX[1]:
            cand, xs = cand[xs > BOX[1]], xs[xs > BOX[1]]
        near = np.sort(np.argsort(np.abs(xs - x), kind="stable")[:npts])
        return cand[near] % self.n, lagrange(xs[near], x)

    def _prolongation_stencils(self, npts: int) -> None:
        """Per fine ghost-frame index (one axis, padded index 0 .. nf + 2G - 1): the coarse stencil and weights."""
        xs = BOX[0] + (np.arange(-G, self.nf + G) + 0.5) * self.h
        st = [self._axis_stencil(x, npts) for x in xs]
        self.p_idx = np.array([s[0] for s in st])
        self.p_wts = np.array([s[1] for s in st])
        pad = np.zeros((self.nf + 2 * G, self.nf + 2 * G), dtype=bool)
        pad[:G, :], pad[-G:, :], pad[:, :G], pad[:, -G:] = True, True, True, True
        self.ghost = np.nonzero(pad)

    def restrict(self, uc: np.ndarray, uf: np.ndarray) -> None:
        """Overwrite the covered coarse cells with the restriction of the fine ones."""
        block = uf[self.r_idx[:, None, :, None], self.r_idx[None, :, None, :]]  # (I, J, a, b)
        uc[self.i0:self.i1, self.i0:self.i1] = np.einsum("ijab,ia,jb->ij", block, self.r_wts, self.r_wts)

    def pad_fine(self, uc: np.ndarray, uf: np.ndarray) -> np.ndarray:
        """Return the fine array with its ghost frame prolonged from the coarse cells outside the box."""
        p = np.zeros((self.nf + 2 * G, self.nf + 2 * G))
        p[G:-G, G:-G] = uf
        gi, gj = self.ghost
        block = uc[self.p_idx[gi][:, :, None], self.p_idx[gj][:, None, :]]
        p[gi, gj] = np.einsum("gab,ga,gb->g", block, self.p_wts[gi], self.p_wts[gj])
        return p

    def rhs(self, uc: np.ndarray, uf: np.ndarray, vel: tuple[float, float]) -> tuple[np.ndarray, np.ndarray]:
        """Return the time derivatives of the coarse (covered cells: zero) and fine grids."""
        if self.refined:
            uc = uc.copy()
            self.restrict(uc, uf)
        pc = np.pad(uc, G, mode="wrap")
        fcx = vel[0] * face_flux(pc[:, G:-G], 0)  # (n + 1, n)
        fcy = vel[1] * face_flux(pc[G:-G, :], 1)  # (n, n + 1)
        if not self.refined:
            return -(np.diff(fcx, axis=0) + np.diff(fcy, axis=1)) / self.H, np.zeros_like(uf)
        pf = self.pad_fine(uc, uf)
        ffx = vel[0] * face_flux(pf[:, G:-G], 0)  # (nf + 1, nf)
        ffy = vel[1] * face_flux(pf[G:-G, :], 1)
        if self.reflux:  # each coarse face on the box boundary takes the mean of its two fine faces
            s = slice(self.i0, self.i1)
            fcx[self.i0, s] = 0.5 * (ffx[0, 0::2] + ffx[0, 1::2])
            fcx[self.i1, s] = 0.5 * (ffx[-1, 0::2] + ffx[-1, 1::2])
            fcy[s, self.i0] = 0.5 * (ffy[0::2, 0] + ffy[1::2, 0])
            fcy[s, self.i1] = 0.5 * (ffy[0::2, -1] + ffy[1::2, -1])
        dc = -(np.diff(fcx, axis=0) + np.diff(fcy, axis=1)) / self.H
        dc[self.i0:self.i1, self.i0:self.i1] = 0.0
        return dc, -(np.diff(ffx, axis=0) + np.diff(ffy, axis=1)) / self.h

    def covered(self) -> np.ndarray:
        """Return the mask of the covered coarse cells."""
        m = np.zeros((self.n, self.n), dtype=bool)
        if self.refined:
            m[self.i0:self.i1, self.i0:self.i1] = True
        return m


def exact(x: np.ndarray, y: np.ndarray, t: float, vel: tuple[float, float]) -> np.ndarray:
    """Return the advected initial condition sin(2 pi x) sin(2 pi y) at time t."""
    return np.sin(2.0 * np.pi * (x - vel[0] * t)) * np.sin(2.0 * np.pi * (y - vel[1] * t))


def run(n: int, arm: tuple[str, str, bool] | None, vel: tuple[float, float], t_end: float = 1.0) -> dict[str, float]:
    """Advance one arm (None: the uniform coarse control) to t_end; return its error norms and its mass drift."""
    grid = Composite(n, *(arm or ("mean", "cubic", False)), refined=arm is not None)
    X, Y = np.meshgrid(grid.xc, grid.xc, indexing="ij")
    XF, YF = np.meshgrid(grid.xf, grid.xf, indexing="ij")
    uc, uf = exact(X, Y, 0.0, vel), exact(XF, YF, 0.0, vel)
    live = ~grid.covered()
    h = grid.h if arm is not None else grid.H
    dt = 0.4 * h * (64.0 * h) ** 0.25 / max(vel)
    steps = math.ceil(t_end / dt)
    dt = t_end / steps

    def mass(uc: np.ndarray, uf: np.ndarray) -> float:
        return float(np.sum(uc[live]) * grid.H**2 + (np.sum(uf) * grid.h**2 if arm is not None else 0.0))

    m0 = mass(uc, uf)
    for _ in range(steps):
        k1 = grid.rhs(uc, uf, vel)
        k2 = grid.rhs(uc + 0.5 * dt * k1[0], uf + 0.5 * dt * k1[1], vel)
        k3 = grid.rhs(uc + 0.5 * dt * k2[0], uf + 0.5 * dt * k2[1], vel)
        k4 = grid.rhs(uc + dt * k3[0], uf + dt * k3[1], vel)
        uc = uc + dt / 6.0 * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0])
        uf = uf + dt / 6.0 * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1])
    ec = np.abs(uc - exact(X, Y, t_end, vel))[live]
    xs, ys, hs = X[live], Y[live], np.full(ec.size, grid.H)
    err = ec
    if arm is not None:
        ef = np.abs(uf - exact(XF, YF, t_end, vel)).ravel()
        err = np.concatenate([ec, ef])
        xs, ys = np.concatenate([xs, XF.ravel()]), np.concatenate([ys, YF.ravel()])
        hs = np.concatenate([hs, np.full(ef.size, grid.h)])
    out = {"l1": float(np.sum(err * hs**2)), "linf": float(err.max()), "mass": mass(uc, uf) - m0}
    # where the error lives: per region, its maximum; the corner region is within 3 coarse cells of a box corner (both
    # axes), the edge region within 3 coarse cells of a box edge otherwise, the interior the rest
    dx = np.minimum(np.abs(xs - BOX[0]), np.abs(xs - BOX[1]))
    dy = np.minimum(np.abs(ys - BOX[0]), np.abs(ys - BOX[1]))
    near = 3.0 * grid.H
    in_x = (ys > BOX[0] - near) & (ys < BOX[1] + near)  # within the box's y span: an x edge is near
    in_y = (xs > BOX[0] - near) & (xs < BOX[1] + near)
    corner = (dx < near) & (dy < near)
    edge = ~corner & (((dx < near) & in_x) | ((dy < near) & in_y))
    if arm is None:
        corner[:], edge[:] = False, False
    for region, mask in (("corner", corner), ("edge", edge), ("interior", ~corner & ~edge)):
        out[f"linf_{region}"] = float(err[mask].max()) if mask.any() else float("nan")
    return out


ARMS = {"mean-cubic-on": ("mean", "cubic", True), "mean-cubic-off": ("mean", "cubic", False),
        "lag6-cubic-on": ("lag6", "cubic", True), "lag6-cubic-off": ("lag6", "cubic", False),
        "lag6-quintic-on": ("lag6", "quintic", True), "lag6-quintic-off": ("lag6", "quintic", False),
        "mean-quintic-off": ("mean", "quintic", False)}
NORMS = ("l1", "linf", "linf_corner", "linf_edge", "linf_interior")


def b1(cells: list[int], vel: tuple[float, float], asserts: list[str]) -> int:
    """Run the B1 table; return the exit status."""
    bounds = {}
    for spec in asserts:
        arm, norm, value = spec.split(":")
        bounds[(arm, norm)] = float(value)
    print(f"B1: u_t + {vel[0]} u_x + {vel[1]} u_y = 0, box {BOX} refined 2:1, t = 1, N = {cells}")
    names = " ".join(f"{n:>13s}" for n in NORMS)
    print(f"{'arm':18s} {'L1 (finest)':>11s} {'mass':>9s} | finest-pair orders: {names}")
    status = 0
    for name, arm in [("uniform", None), *ARMS.items()]:
        res = [run(n, arm, vel) for n in cells]
        orders = {k: math.log2(res[-2][k] / res[-1][k]) for k in NORMS if not math.isnan(res[-1][k])}
        all_orders = " ".join(f"{orders.get(k, float('nan')):13.2f}" for k in NORMS)
        print(f"{name:18s} {res[-1]['l1']:11.3e} {res[-1]['mass']:9.1e} | {'':20s}{all_orders}", flush=True)
        for (a, k), b in bounds.items():
            if a == name and orders.get(k, -math.inf) < b:
                status = 1
                print(f"  FAIL {name}: {k} order {orders.get(k, float('nan')):.2f} < {b}")
    return status


def main() -> int:
    """Dispatch the subcommand."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("b1")
    p.add_argument("--cells", type=int, nargs="+", default=[32, 64, 128])
    p.add_argument("--velocity", type=float, nargs=2, default=[1.0, 1.0])
    p.add_argument("--assert", dest="asserts", action="append", default=[], metavar="ARM:NORM:MIN")
    args = parser.parse_args()
    return b1(args.cells, tuple(args.velocity), args.asserts)


if __name__ == "__main__":
    sys.exit(main())
