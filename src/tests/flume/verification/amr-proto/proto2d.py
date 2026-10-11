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

B2 (the go/no-go of #68 E, D-M6-4): the same composite grid with WENO5-JS (Lax-Friedrichs splitting, the weights of
adam_weno_object) and a third reflux policy, troubled: a coarse box face takes the mean of its fine fluxes only where
the WENO sensor r = min_k w_k / d_k of the coarse face or of either fine face is below tau (0.2, the weno-riemann
threshold). Arms: today's mean/cubic/always, and lag6/quintic under always, troubled and never. Legs:
  smooth  the B1 advection: orders, and the share of seam faces the sensor flags (false positives);
  c2f     a Burgers pulse (u = 1 on [0, 0.15]) whose shock crosses x = 0.25 coarse to fine;
  f2c     a pulse on [0.45, 0.6] whose shock crosses x = 0.75 fine to coarse;
  weak    a 1e-3 pulse riding on u = 1 (crossing x = 0.25): its smoothness indicators are near zeps, so the JS sensor
          misses it, the missed-detection control.
The pulses start from exact cell averages (the discrete mass is the exact one) and are measured against the exact
solution: the mass defect over the pulse's mass, and the mass-equivalent shock position of every row, in local cells.

Usage:
    proto2d.py b1 [--cells 32 64 128] [--velocity 1 1] [--assert ARM:NORM:MIN ...]
    proto2d.py b2 [--cells 32 64 128] [--shock-cells 64 128] [--leg smooth|c2f|f2c|weak ...] [--tau 0.2] [--check]
    --assert lag6-quintic-off:linf:4.5 fails unless that arm's finest-pair order of that norm is at least MIN
    --check asserts the B2 criteria (see b2)
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


def weno5(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray,
          e: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return the WENO5-JS value at i + 1/2 from the cells i-2 .. i+2 and its sensor r = min_k w_k / d_k.

    The weights are those of adam_weno_object: zeps = 1e-6, exponent 3 (weno_weights_exponent(S = 3)), the
    Jiang-Shu smoothness indicators (its table expands to them)."""
    is0 = 13.0 / 12.0 * (a - 2.0 * b + c) ** 2 + 0.25 * (a - 4.0 * b + 3.0 * c) ** 2
    is1 = 13.0 / 12.0 * (b - 2.0 * c + d) ** 2 + 0.25 * (b - d) ** 2
    is2 = 13.0 / 12.0 * (c - 2.0 * d + e) ** 2 + 0.25 * (3.0 * c - 4.0 * d + e) ** 2
    lin = (0.1, 0.6, 0.3)
    alpha = [lin[k] / (1.0e-6 + s) ** 3 for k, s in enumerate((is0, is1, is2))]
    tot = alpha[0] + alpha[1] + alpha[2]
    w = [al / tot for al in alpha]
    value = (w[0] * (2.0 * a - 7.0 * b + 11.0 * c) + w[1] * (-b + 5.0 * c + 2.0 * d)
             + w[2] * (2.0 * c + 5.0 * d - e)) / 6.0
    return value, np.minimum.reduce([w[k] / lin[k] for k in range(3)])


def weno_flux(p: np.ndarray, axis: int, f, alpha: float) -> tuple[np.ndarray, np.ndarray]:
    """Return the WENO5-JS Lax-Friedrichs face fluxes of the flux function f (padded array, n + 1 faces) and their
    sensor, the minimum of the two split reconstructions."""
    p = np.moveaxis(p, axis, 0)
    n = p.shape[0] - 2 * G
    fp, fm = 0.5 * (f(p) + alpha * p), 0.5 * (f(p) - alpha * p)
    s = [slice(k, k + n + 1) for k in range(6)]  # face j between padded cells j + 2 and j + 3
    plus, rp = weno5(fp[s[0]], fp[s[1]], fp[s[2]], fp[s[3]], fp[s[4]])
    minus, rm = weno5(fm[s[5]], fm[s[4]], fm[s[3]], fm[s[2]], fm[s[1]])
    return np.moveaxis(plus + minus, 0, axis), np.moveaxis(np.minimum(rp, rm), 0, axis)


class Law:
    """A scalar conservation law u_t + f(u)_x + g(u)_y = 0 and the Lax-Friedrichs speeds of its splitting."""

    def __init__(self, kind: str, vel: tuple[float, float] = (1.0, 1.0), umax: float = 1.0) -> None:
        self.kind = kind
        if kind == "advect":
            self.fx, self.fy = (lambda u: vel[0] * u), (lambda u: vel[1] * u)
            self.ax, self.ay = abs(vel[0]), abs(vel[1])
        else:  # burgers along x: f = u^2 / 2, g = 0 (the data vary along x only)
            self.fx, self.fy = (lambda u: 0.5 * u * u), (lambda u: 0.0 * u)
            self.ax, self.ay = umax, 0.0
        self.speed = max(self.ax, self.ay)


class Composite:
    """A periodic coarse grid with the box refined 2:1, and the transfer operators of one arm.

    reflux: True or "always" (every coarse box face takes the mean of its two fine fluxes), False or "never", or
    "troubled" (only where the sensor of the coarse face or of either fine face is below tau; scheme weno only)."""

    def __init__(self, n: int, restriction: str, prolongation: str, reflux: bool | str, refined: bool = True,
                 scheme: str = "linear", law: Law | None = None, tau: float = 0.2) -> None:
        self.n, self.H = n, 1.0 / n
        self.h = self.H / 2.0
        self.refined = refined
        self.policy = {True: "always", False: "never"}.get(reflux, reflux)
        self.scheme, self.law, self.tau = scheme, law, tau
        self.trouble = [0, 0]  # troubled and total coarse box faces met by the policy
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

    def fluxes(self, p: np.ndarray, vel: tuple[float, float]) -> tuple[np.ndarray, ...]:
        """Return the x and y face fluxes of a padded array and their sensors (None for the linear scheme)."""
        px, py = p[:, G:-G], p[G:-G, :]
        if self.scheme == "linear":
            return vel[0] * face_flux(px, 0), vel[1] * face_flux(py, 1), None, None
        fx, rx = weno_flux(px, 0, self.law.fx, self.law.ax)
        fy, ry = weno_flux(py, 1, self.law.fy, self.law.ay)
        return fx, fy, rx, ry

    def rhs(self, uc: np.ndarray, uf: np.ndarray, vel: tuple[float, float]) -> tuple[np.ndarray, np.ndarray]:
        """Return the time derivatives of the coarse (covered cells: zero) and fine grids."""
        if self.refined:
            uc = uc.copy()
            self.restrict(uc, uf)
        fcx, fcy, rcx, rcy = self.fluxes(np.pad(uc, G, mode="wrap"), vel)  # (n + 1, n), (n, n + 1)
        if not self.refined:
            return -(np.diff(fcx, axis=0) + np.diff(fcy, axis=1)) / self.H, np.zeros_like(uf)
        ffx, ffy, rfx, rfy = self.fluxes(self.pad_fine(uc, uf), vel)  # (nf + 1, nf), (nf, nf + 1)
        if self.policy != "never":  # a coarse face on the box boundary takes the mean of its two fine faces
            s = slice(self.i0, self.i1)
            faces = ((fcx, rcx, (self.i0, s), ffx, rfx, (0, slice(None))),
                     (fcx, rcx, (self.i1, s), ffx, rfx, (-1, slice(None))),
                     (fcy, rcy, (s, self.i0), ffy, rfy, (slice(None), 0)),
                     (fcy, rcy, (s, self.i1), ffy, rfy, (slice(None), -1)))
            for fc, rc, cidx, ff, rf, fidx in faces:
                fine = ff[fidx]
                mean = 0.5 * (fine[0::2] + fine[1::2])
                if self.policy == "always":
                    fc[cidx] = mean
                    continue
                rfine = rf[fidx]
                bad = (rc[cidx] < self.tau) | (rfine[0::2] < self.tau) | (rfine[1::2] < self.tau)
                fc[cidx] = np.where(bad, mean, fc[cidx])
                self.trouble[0] += int(bad.sum())
                self.trouble[1] += bad.size
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


def advance(grid: Composite, uc: np.ndarray, uf: np.ndarray, vel: tuple[float, float], speed: float,
            t_end: float, refined: bool) -> tuple[np.ndarray, np.ndarray]:
    """Advance both grids to t_end with classical RK4, dt = 0.4 h (64 h)^(1/4) / speed (Appendix A)."""
    h = grid.h if refined else grid.H
    dt = 0.4 * h * (64.0 * h) ** 0.25 / speed
    steps = math.ceil(t_end / dt)
    dt = t_end / steps
    for _ in range(steps):
        k1 = grid.rhs(uc, uf, vel)
        k2 = grid.rhs(uc + 0.5 * dt * k1[0], uf + 0.5 * dt * k1[1], vel)
        k3 = grid.rhs(uc + 0.5 * dt * k2[0], uf + 0.5 * dt * k2[1], vel)
        k4 = grid.rhs(uc + dt * k3[0], uf + dt * k3[1], vel)
        uc = uc + dt / 6.0 * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0])
        uf = uf + dt / 6.0 * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1])
    return uc, uf


def run(n: int, arm: tuple[str, str, bool | str] | None, vel: tuple[float, float], t_end: float = 1.0,
        scheme: str = "linear", tau: float = 0.2) -> dict[str, float]:
    """Advance one arm (None: the uniform coarse control) to t_end; return its error norms and its mass drift."""
    grid = Composite(n, *(arm or ("mean", "cubic", False)), refined=arm is not None, scheme=scheme,
                     law=Law("advect", vel), tau=tau)
    X, Y = np.meshgrid(grid.xc, grid.xc, indexing="ij")
    XF, YF = np.meshgrid(grid.xf, grid.xf, indexing="ij")
    uc, uf = exact(X, Y, 0.0, vel), exact(XF, YF, 0.0, vel)
    live = ~grid.covered()

    def mass(uc: np.ndarray, uf: np.ndarray) -> float:
        return float(np.sum(uc[live]) * grid.H**2 + (np.sum(uf) * grid.h**2 if arm is not None else 0.0))

    m0 = mass(uc, uf)
    uc, uf = advance(grid, uc, uf, vel, max(vel), t_end, arm is not None)
    ec = np.abs(uc - exact(X, Y, t_end, vel))[live]
    xs, ys, hs = X[live], Y[live], np.full(ec.size, grid.H)
    err = ec
    if arm is not None:
        ef = np.abs(uf - exact(XF, YF, t_end, vel)).ravel()
        err = np.concatenate([ec, ef])
        xs, ys = np.concatenate([xs, XF.ravel()]), np.concatenate([ys, YF.ravel()])
        hs = np.concatenate([hs, np.full(ef.size, grid.h)])
    out = {"l1": float(np.sum(err * hs**2)), "linf": float(err.max()), "mass": mass(uc, uf) - m0,
           "trouble": grid.trouble[0] / grid.trouble[1] if grid.trouble[1] else float("nan")}
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


B2_SMOOTH = {"mean-cubic-always": ("mean", "cubic", "always"), "lag6-quintic-always": ("lag6", "quintic", "always"),
             "lag6-quintic-troubled": ("lag6", "quintic", "troubled"),
             "lag6-quintic-never": ("lag6", "quintic", "never")}
# Burgers pulses u = ub + amp on [x0, x0 + w] (ub elsewhere), y-independent; the shock (the pulse's right edge) crosses
# a box face: c2f x = 0.25 at t = 0.2 (the rarefaction catches it at t = 0.3); f2c x = 0.75 at t = 0.3; weak a jump
# of 1e-3 riding on u = 1, crossing x = 0.25 at t = 0.1, whose smoothness indicators are near zeps: the sensor's blind
# spot (the missed-detection case)
PULSES = {"c2f": {"ub": 0.0, "amp": 1.0, "x0": 0.0, "w": 0.15, "t_end": 0.5},
          "f2c": {"ub": 0.0, "amp": 1.0, "x0": 0.45, "w": 0.15, "t_end": 0.6},
          "weak": {"ub": 1.0, "amp": 1.0e-3, "x0": 0.0, "w": 0.15, "t_end": 0.5}}


def pulse(x: np.ndarray, t: float, ub: float, amp: float, x0: float, w: float, **_: float) -> np.ndarray:
    """Return the exact Burgers solution of the pulse at time t: in xi = x - ub t and tau = amp t, the unit pulse on
    [x0, x0 + w] (a rarefaction from x0, a shock at x0 + w + tau/2, which the rarefaction catches at tau = 2w, after
    which the shock sits at x0 + sqrt(2 w tau))."""
    xi, tau = x - ub * t, amp * t
    if tau == 0.0:
        v = ((xi >= x0) & (xi < x0 + w)).astype(float)
    elif tau <= 2.0 * w:
        v = np.where((xi >= x0) & (xi < x0 + tau), (xi - x0) / tau, 0.0)
        v = np.where((xi >= x0 + tau) & (xi < x0 + w + 0.5 * tau), 1.0, v)
    else:
        v = np.where((xi >= x0) & (xi < x0 + math.sqrt(2.0 * w * tau)), (xi - x0) / tau, 0.0)
    return ub + amp * v


def pulse_average(x: np.ndarray, h: float, ub: float, amp: float, x0: float, w: float, **_: float) -> np.ndarray:
    """Return the cell averages of the initial pulse: the discrete mass is then the exact one, so a shock position read
    from the mass in a window measures the scheme, not where the initial step falls among the cell centres."""
    overlap = np.clip(np.minimum(x + 0.5 * h, x0 + w) - np.maximum(x - 0.5 * h, x0), 0.0, None)
    return ub + amp * overlap / h


def shock_at(t: float, ub: float, amp: float, x0: float, w: float, **_: float) -> tuple[float, float]:
    """Return the exact shock position and jump of the pulse at time t."""
    tau = amp * t
    if tau <= 2.0 * w:
        return x0 + w + 0.5 * tau + ub * t, amp
    xs = x0 + math.sqrt(2.0 * w * tau)
    return xs + ub * t, amp * (xs - x0) / tau


def run_shock(n: int, arm: tuple[str, str, str] | None, case: str, tau: float = 0.2) -> dict[str, float]:
    """Advance a Burgers pulse on one arm; return the mass defect, the shock position errors and the troubled share.

    The shock position of a row is mass-equivalent: the cells within W = 6 H of the exact shock hold the exact
    integral plus jump * dx, so dx = (sum u h - integral of the exact u) / jump, in cells of the spacing that holds the
    exact shock in that row."""
    pp = PULSES[case]
    law = Law("burgers", umax=pp["ub"] + pp["amp"])
    grid = Composite(n, *(arm or ("mean", "cubic", "never")), refined=arm is not None, scheme="weno", law=law,
                     tau=tau)
    X, _ = np.meshgrid(grid.xc, grid.xc, indexing="ij")
    XF, _ = np.meshgrid(grid.xf, grid.xf, indexing="ij")
    uc, uf = pulse_average(X, grid.H, **pp), pulse_average(XF, grid.h, **pp)
    live = ~grid.covered()

    def mass(uc: np.ndarray, uf: np.ndarray) -> float:
        return float(np.sum(uc[live]) * grid.H**2 + (np.sum(uf) * grid.h**2 if arm is not None else 0.0))

    m0 = mass(uc, uf)
    uc, uf = advance(grid, uc, uf, (1.0, 1.0), law.speed, pp["t_end"], arm is not None)
    xs, jump = shock_at(pp["t_end"], **pp)
    width = 6.0 * grid.H
    errors = []
    rows = [(None, j) for j in range(grid.n)]
    if arm is not None:  # the rows inside the box's y span are fine rows: coarse cells outside, fine cells inside
        rows = [(None, j) for j in range(grid.n) if not grid.i0 <= j < grid.i1] + [(k, None) for k in range(grid.nf)]
    for k, j in rows:
        if k is None:
            x, u, h = grid.xc, uc[:, j], np.full(grid.n, grid.H)
        else:
            out = (grid.xc < BOX[0]) | (grid.xc > BOX[1])
            jc = grid.i0 + k // 2
            x = np.concatenate([grid.xc[out], grid.xf])
            u = np.concatenate([uc[out, jc], uf[:, k]])
            h = np.concatenate([np.full(out.sum(), grid.H), np.full(grid.nf, grid.h)])
        win = np.abs(x - xs) <= width
        a, b = float(np.min(x[win] - 0.5 * h[win])), float(np.max(x[win] + 0.5 * h[win]))
        sample = np.linspace(a, b, 400001)
        exact_int = float(np.trapezoid(pulse(sample, pp["t_end"], **pp), sample))
        local = h[np.argmin(np.abs(x - xs))]
        errors.append((float(np.sum(u[win] * h[win])) - exact_int) / jump / local)
    errors = np.array(errors)
    # the mass defect relative to the pulse's own mass (amp w), not to the total, which the weak case's u = 1 dominates
    return {"mass": (mass(uc, uf) - m0) / (pp["amp"] * pp["w"]), "shock_max": float(np.max(np.abs(errors))),
            "shock_mean": float(np.mean(errors)),
            "trouble": grid.trouble[0] / grid.trouble[1] if grid.trouble[1] else float("nan")}


def b2(cells: list[int], shock_cells: list[int], vel: tuple[float, float], legs: list[str], tau: float,
       check: bool) -> int:
    """Run the B2 tables (smooth order, shock crossings); return the exit status.

    With check: the smooth order of troubled (and of the uniform control) at least 4.7 in L1 and Linf; on the strong
    shocks (c2f, f2c) troubled within 0.05 cells (mean shock position) of always and its mass defect at most 1e-3 of
    the pulse mass, while never's is at least 1e-3 (the negative control: the metrics see a missing reflux)."""
    status = 0

    def verdict(ok: bool, what: str) -> None:
        nonlocal status
        status |= 0 if ok else 1
        print(f"  {'PASS' if ok else 'FAIL'} {what}")

    if "smooth" in legs:
        print(f"B2 smooth: WENO5-JS, u_t + {vel[0]} u_x + {vel[1]} u_y = 0, box {BOX} refined 2:1, t = 1, "
              f"N = {cells}, tau = {tau}")
        names = " ".join(f"{n:>13s}" for n in NORMS)
        print(f"{'arm':24s} {'L1 (finest)':>11s} {'mass':>9s} {'troubled':>9s} | finest-pair orders: {names}")
        for name, arm in [("uniform", None), *B2_SMOOTH.items()]:
            res = [run(n, arm, vel, scheme="weno", tau=tau) for n in cells]
            orders = {k: math.log2(res[-2][k] / res[-1][k]) for k in NORMS if not math.isnan(res[-1][k])}
            all_orders = " ".join(f"{orders.get(k, float('nan')):13.2f}" for k in NORMS)
            print(f"{name:24s} {res[-1]['l1']:11.3e} {res[-1]['mass']:9.1e} {res[-1]['trouble']:9.3f} | {'':20s}"
                  f"{all_orders}", flush=True)
            if check and name in ("uniform", "lag6-quintic-troubled"):
                verdict(min(orders["l1"], orders["linf"]) >= 4.7, f"{name}: L1 and Linf orders at least 4.7")
    for case in (c for c in PULSES if c in legs):
        print(f"B2 shock {case}: Burgers pulse {PULSES[case]}, WENO5-JS, tau = {tau}; shock position error in local "
              "cells (max |.| and mean over the rows), mass defect relative")
        for n in shock_cells:
            res = {}
            for name, arm in [("uniform", None), *B2_SMOOTH.items()]:
                r = res[name] = run_shock(n, arm, case, tau=tau)
                print(f"  N = {n:4d} {name:24s} mass {r['mass']:10.2e}  shock max {r['shock_max']:7.3f}  mean "
                      f"{r['shock_mean']:7.3f}  troubled {r['trouble']:6.3f}", flush=True)
            if check and case != "weak":
                t, a, v = res["lag6-quintic-troubled"], res["lag6-quintic-always"], res["lag6-quintic-never"]
                verdict(abs(t["shock_mean"] - a["shock_mean"]) <= 0.05, "troubled places the shock as always (0.05)")
                verdict(abs(t["mass"]) <= 1.0e-3, "troubled mass defect at most 1e-3 of the pulse")
                verdict(abs(v["mass"]) >= 1.0e-3, "never's mass defect at least 1e-3 (negative control)")
    return status


def main() -> int:
    """Dispatch the subcommand."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("b1")
    p.add_argument("--cells", type=int, nargs="+", default=[32, 64, 128])
    p.add_argument("--velocity", type=float, nargs=2, default=[1.0, 1.0])
    p.add_argument("--assert", dest="asserts", action="append", default=[], metavar="ARM:NORM:MIN")
    q = sub.add_parser("b2")
    q.add_argument("--cells", type=int, nargs="+", default=[32, 64, 128])
    q.add_argument("--shock-cells", type=int, nargs="+", default=[64, 128])
    q.add_argument("--velocity", type=float, nargs=2, default=[1.0, 1.0])
    q.add_argument("--leg", dest="legs", action="append", default=[], choices=("smooth", *PULSES))
    q.add_argument("--tau", type=float, default=0.2)
    q.add_argument("--check", action="store_true", help="assert the B2 go/no-go criteria")
    args = parser.parse_args()
    if args.cmd == "b2":
        return b2(args.cells, args.shock_cells, tuple(args.velocity), args.legs or ["smooth", *PULSES], args.tau,
                  args.check)
    return b1(args.cells, tuple(args.velocity), args.asserts)


if __name__ == "__main__":
    sys.exit(main())
