#!/usr/bin/env python3
"""FLUME isothermal and barotropic equation-of-state oracle (issue #100, M13): the exact solutions P1-P4 check against.

    eos_oracle.py selftest                    the oracle checks itself against properties independent of its construction
    eos_oracle.py riemann <work> [--case tube|collide] [--c0 C0] [--mach M]
                                              a run against the exact solution: the star plateau, the shock and contact
                                              positions, the inert energy slot, the L1 density error (reported)
    eos_oracle.py wave <work> ... [--c0 C0]   acoustic-wave runs, coarse to fine: the L1 error after one period, the
                                              observed orders, the finest above ORDER_MIN
    eos_oracle.py inert <work>                every checkpoint's energy slot exactly 0, the total mass constant

The equation of state is barotropic, P = P(rho), so the hydro system is mass and momentum only; FLUME keeps its state
layout and the energy slot is inert (D-M13-1): zero flux, zero source, its eigenvector the unit vector. Two forms:

  isothermal   P = c0^2 rho, c^2 = c0^2;
  barotropic   a continuous piecewise polytrope P = K_k rho^G_k on rho_{k-1} <= rho < rho_k (rho_0 = 0, the last
               segment unbounded), K_0 = c0^2 rho_ref^(1 - G_0) and K_{k+1} = K_k rho_k^(G_k - G_{k+1}) (continuity);
               c^2 = dP/drho = G_k P / rho jumps at a break where G changes, P does not.

The speed of sound is c^2 = dP/drho everywhere (D-M13-4): the Jeans dispersion relation and the characteristic speeds
use it. The Riemann invariants of the hydro system are u + phi(rho) (carried by the left fan) and u - phi(rho) (the
right fan), phi(rho) = int c / rho drho, in closed form on each segment: phi = 2 c / (G - 1) for G != 1 and
c0 ln rho for G = 1, offsets added so that phi is continuous.

The exact Riemann solver (Toro's convention): the star density rho* satisfies u_L - f_L(rho*) = u_R + f_R(rho*), with,
for each side K,
  a shock (rho* > rho_K):        f = sqrt((P* - P_K) (1 / rho_K - 1 / rho*)) > 0, the barotropic Rankine-Hugoniot jump;
  a rarefaction (rho* <= rho_K): f = phi(rho*) - phi(rho_K) <= 0, the Riemann invariant.
f_L + f_R increases with rho*, so the root is bracketed and found by bisection in ln rho. The transverse velocities
jump only at the contact u*. There is no energy equation: an isothermal shock compresses by rho_2 / rho_1 = M_1^2, M_1
the upstream Mach number in the shock frame, without the adiabatic bound (gamma + 1) / (gamma - 1).

The eigensystem of the x flux F = (rho u, rho u^2 + P, rho u v, rho u w, 0) in U = (rho, rho u, rho v, rho w, slot5):
eigenvalues u - c, u, u, u + c, 0; right vectors (1, u - c, v, w, 0), (0, 0, 1, 0, 0), (0, 0, 0, 1, 0),
(1, u + c, v, w, 0), (0, 0, 0, 0, 1); the y and z systems by permuting the momentum components.

The self-test (all to round-off unless stated):
  rh        every shock of the sampled problems satisfies the Rankine-Hugoniot mass and momentum conditions in its frame;
  ratio     an isothermal shock compresses by M_1^2 (two colliding streams at Mach 1.5 to 50);
  fan       inside every rarefaction the sampled state lies on its own characteristic, u -+ c(rho) = x / t, and the fan
            meets the star state at its tail (1e-10, the fan bisection);
  phi       the closed-form phi against Gauss-Legendre quadrature of c / rho, across the breaks;
  closed    the isothermal wave functions against their closed forms, c0 (rho* - rho) / sqrt(rho* rho) for a shock and
            c0 ln(rho* / rho) for a rarefaction;
  mirror    the mirrored problem has the mirrored solution;
  toro      a polytropic two-rarefaction problem matches Toro's ideal-gas isentropic star pressure (Toro 2009,
            eq. 4.46), the barotropic invariant being the isentropic one (1e-11);
  eigen     L R = I and R diag(lambda) L = dF/dU at random states, in x, y and z;
  flux      the reference LLF, HLL and isothermal HLLC fluxes are consistent (F(U, U) = F(U)) and HLLC resolves an
            isolated contact (u = 0, a jump in the transverse velocities only) exactly;
  jeans     the Jeans dispersion relation omega^2 = c^2 k^2 - 4 pi G rho_0 vanishes at k_J, growth below, oscillation above;
  default   the Mathew & Federrath 2021 table: P continuous at every break, c^2 / (P / rho) = G on each segment;
  convex    the breaks where G decreases upward are flagged, and the exact solver refuses such a law (below).

Convexity. Inside a segment the system is genuinely nonlinear, (rho c)' = (G + 1) c / 2 > 0. At a break where G
increases upward c jumps up, and a rarefaction through it splits into two fans with a plateau at the break density
(sample() returns it: the fan bisection converges to the break). Where G decreases upward, rho c jumps down: the
system is not genuinely nonlinear there, and the exact Riemann solutions crossing that break contain composite waves
(shock-rarefaction), which this solver does not construct, so it refuses such a law. The Mathew & Federrath 2021
table is non-convex at its 1.4 -> 1.1 break, rho = 3.84e-8 g cm^-3: the exact-solution tests use convex tables, and
the code's schemes (entropy-satisfying LLF/HLL fluxes) do not depend on convexity.
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

TOL = 1.0e-12
TOL_FAN = 1.0e-10
TOL_TAIL = 1.0e-7
TOL_TORO = 1.0e-11
BISECTION_ITERATIONS = 200
GAUSS_POINTS = 40

# Mathew & Federrath 2021 (MNRAS 507, 2448, eqs. 5-7, after Larson 1969 and Masunaga & Inutsuka 2000): the breaks in
# g cm^-3 and the polytropic exponents; c_s0 = 0.2 km/s.
MF2021_BREAKS = (2.50e-16, 3.84e-13, 3.84e-8, 3.84e-3)
MF2021_EXPONENTS = (1.0, 1.1, 1.4, 1.1, 5.0 / 3.0)


@dataclass
class Eos:
    """A continuous piecewise polytrope; one segment with G = 1 is the isothermal gas."""

    c0: float
    breaks: tuple[float, ...] = ()
    exponents: tuple[float, ...] = (1.0,)
    rho_ref: float = 1.0
    k: list[float] = field(init=False)
    phi0: list[float] = field(init=False)

    def __post_init__(self) -> None:
        if len(self.exponents) != len(self.breaks) + 1:
            raise ValueError("one exponent per segment: len(exponents) == len(breaks) + 1")
        if any(b <= a for a, b in zip(self.breaks[:-1], self.breaks[1:], strict=True)):
            raise ValueError("the breaks must increase")
        self.k = [self.c0**2 * self.rho_ref ** (1.0 - self.exponents[0])]
        for i, rho_b in enumerate(self.breaks):
            self.k.append(self.k[i] * rho_b ** (self.exponents[i] - self.exponents[i + 1]))
        # phi continuous: each segment's closed form plus an offset fixed at the break below it.
        self.phi0 = [0.0]
        for i, rho_b in enumerate(self.breaks):
            self.phi0.append(self.phi0[i] + self._phi_raw(rho_b, i) - self._phi_raw(rho_b, i + 1))

    def segment(self, rho: float) -> int:
        """Segment index of rho (a break belongs to the segment above it)."""
        return sum(1 for b in self.breaks if rho >= b)

    def p(self, rho: float) -> float:
        i = self.segment(rho)
        return self.k[i] * rho ** self.exponents[i]

    def c2(self, rho: float) -> float:
        """c^2 = dP/drho."""
        i = self.segment(rho)
        return self.exponents[i] * self.k[i] * rho ** (self.exponents[i] - 1.0)

    def c(self, rho: float) -> float:
        return math.sqrt(self.c2(rho))

    def c_in(self, rho: float, i: int) -> float:
        """c of segment i's law at rho (a one-sided limit at a break)."""
        return math.sqrt(self.exponents[i] * self.k[i] * rho ** (self.exponents[i] - 1.0))

    def nonconvex_breaks(self) -> list[float]:
        """The breaks where G decreases upward: rho c (so the characteristic speed) jumps down there and the system
        loses genuine nonlinearity (inside a segment (rho c)' = (G + 1) rho c / (2 rho) > 0)."""
        return [b for i, b in enumerate(self.breaks) if self.exponents[i + 1] < self.exponents[i]]

    def _phi_raw(self, rho: float, i: int) -> float:
        g = self.exponents[i]
        if g == 1.0:
            return math.sqrt(self.k[i]) * math.log(rho)
        return 2.0 * math.sqrt(g * self.k[i] * rho ** (g - 1.0)) / (g - 1.0)

    def phi(self, rho: float) -> float:
        """The Riemann invariant integral int c / rho drho, continuous across the breaks."""
        i = self.segment(rho)
        return self.phi0[i] + self._phi_raw(rho, i)


def isothermal(c0: float) -> Eos:
    return Eos(c0=c0)


@dataclass
class State:
    rho: float
    u: float
    v: float = 0.0
    w: float = 0.0


def wave_function(eos: Eos, side: State, rho_star: float) -> float:
    """f_K(rho*), u* = u_L - f_L = u_R + f_R: positive across a shock, negative across a rarefaction."""
    if rho_star > side.rho:
        return math.sqrt((eos.p(rho_star) - eos.p(side.rho)) * (1.0 / side.rho - 1.0 / rho_star))
    return eos.phi(rho_star) - eos.phi(side.rho)


def star_density(eos: Eos, left: State, right: State) -> float:
    """rho* from f_L(rho*) + f_R(rho*) + u_R - u_L = 0, by bisection in ln rho (the residual increases with rho*)."""

    def residual(lnrho: float) -> float:
        r = math.exp(lnrho)
        return wave_function(eos, left, r) + wave_function(eos, right, r) + right.u - left.u

    lo = math.log(min(left.rho, right.rho)) - 1.0
    while residual(lo) > 0.0:
        lo -= 1.0
        if lo < -700.0:
            raise ValueError("vacuum: the star density underflows")
    hi = math.log(max(left.rho, right.rho)) + 1.0
    while residual(hi) < 0.0:
        hi += 1.0
    for _ in range(BISECTION_ITERATIONS):
        mid = 0.5 * (lo + hi)
        if residual(mid) > 0.0:
            hi = mid
        else:
            lo = mid
    return math.exp(0.5 * (lo + hi))


@dataclass
class Solution:
    eos: Eos
    left: State
    right: State
    rho_star: float = field(init=False)
    u_star: float = field(init=False)

    def __post_init__(self) -> None:
        if self.eos.nonconvex_breaks():
            raise ValueError(
                f"non-convex barotropic law (G decreases at {self.eos.nonconvex_breaks()}): its exact "
                "Riemann solutions contain composite waves, which this solver does not construct"
            )
        self.rho_star = star_density(self.eos, self.left, self.right)
        self.u_star = self.left.u - wave_function(self.eos, self.left, self.rho_star)

    def shock_speed(self, side: str) -> float:
        """Speed of the shock on 'left' or 'right' (only meaningful when that wave is a shock)."""
        s = self.left if side == "left" else self.right
        m = math.sqrt((self.eos.p(self.rho_star) - self.eos.p(s.rho)) / (1.0 / s.rho - 1.0 / self.rho_star))
        return s.u - m / s.rho if side == "left" else s.u + m / s.rho

    def _fan_density(self, side: State, xi: float, sign: float) -> float:
        """Inside a fan (sign +1 left, -1 right): u - sign c = xi, u = side.u - sign (phi(rho) - phi(rho_K)).

        g(rho) = u - sign c - xi decreases with rho in the left fan and increases in the right one: g > 0 means the
        root lies above rho in the left fan and below it in the right one.
        """
        lo, hi = math.log(self.rho_star), math.log(side.rho)
        for _ in range(BISECTION_ITERATIONS):
            mid = 0.5 * (lo + hi)
            rho = math.exp(mid)
            u = side.u - sign * (self.eos.phi(rho) - self.eos.phi(side.rho))
            g = u - sign * self.eos.c(rho) - xi
            if (g > 0.0) == (sign > 0.0):
                lo = mid
            else:
                hi = mid
        return math.exp(0.5 * (lo + hi))

    def sample(self, xi: float) -> State:
        """The self-similar solution at x / t = xi."""
        eos, left, right = self.eos, self.left, self.right
        if xi <= self.u_star:  # left of the contact
            if self.rho_star > left.rho:
                if xi < self.shock_speed("left"):
                    return left
                return State(self.rho_star, self.u_star, left.v, left.w)
            head, tail = left.u - eos.c(left.rho), self.u_star - eos.c(self.rho_star)
            if xi <= head:
                return left
            if xi >= tail:
                return State(self.rho_star, self.u_star, left.v, left.w)
            rho = self._fan_density(left, xi, +1.0)
            return State(rho, left.u - (eos.phi(rho) - eos.phi(left.rho)), left.v, left.w)
        if self.rho_star > right.rho:
            if xi > self.shock_speed("right"):
                return right
            return State(self.rho_star, self.u_star, right.v, right.w)
        head, tail = right.u + eos.c(right.rho), self.u_star + eos.c(self.rho_star)
        if xi >= head:
            return right
        if xi <= tail:
            return State(self.rho_star, self.u_star, right.v, right.w)
        rho = self._fan_density(right, xi, -1.0)
        return State(rho, right.u + (eos.phi(rho) - eos.phi(right.rho)), right.v, right.w)


def conserved(s: State) -> np.ndarray:
    """U = (rho, rho u, rho v, rho w, slot5); slot 5 holds its frozen initial value, 0."""
    return np.array([s.rho, s.rho * s.u, s.rho * s.v, s.rho * s.w, 0.0])


DIRECTION = {0: (1, 2, 3), 1: (2, 3, 1), 2: (3, 1, 2)}  # normal, then the two tangential momentum slots


def flux(eos: Eos, q: np.ndarray, d: int = 0) -> np.ndarray:
    n, t1, t2 = DIRECTION[d]
    rho, un = q[0], q[n] / q[0]
    f = np.zeros(5)
    f[0] = q[n]
    f[n] = q[n] * un + eos.p(rho)
    f[t1] = q[t1] * un
    f[t2] = q[t2] * un
    return f


def jacobian(eos: Eos, q: np.ndarray, d: int = 0) -> np.ndarray:
    n, t1, t2 = DIRECTION[d]
    rho = q[0]
    un, ut1, ut2 = q[n] / rho, q[t1] / rho, q[t2] / rho
    a = np.zeros((5, 5))
    a[0, n] = 1.0
    a[n, 0], a[n, n] = eos.c2(rho) - un * un, 2.0 * un
    a[t1, 0], a[t1, n], a[t1, t1] = -un * ut1, ut1, un
    a[t2, 0], a[t2, n], a[t2, t2] = -un * ut2, ut2, un
    return a


def eigensystem(eos: Eos, q: np.ndarray, d: int = 0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(lambda, R, L): the columns of R are the right vectors, the rows of L the left ones, L R = I."""
    n, t1, t2 = DIRECTION[d]
    rho = q[0]
    un, ut1, ut2 = q[n] / rho, q[t1] / rho, q[t2] / rho
    c = eos.c(rho)
    lam = np.array([un - c, un, un, un + c, 0.0])
    r = np.zeros((5, 5))
    r[0, 0], r[n, 0], r[t1, 0], r[t2, 0] = 1.0, un - c, ut1, ut2
    r[t1, 1] = 1.0
    r[t2, 2] = 1.0
    r[0, 3], r[n, 3], r[t1, 3], r[t2, 3] = 1.0, un + c, ut1, ut2
    r[4, 4] = 1.0
    lv = np.zeros((5, 5))
    lv[0, 0], lv[0, n] = (un + c) / (2.0 * c), -1.0 / (2.0 * c)
    lv[1, 0], lv[1, t1] = -ut1, 1.0
    lv[2, 0], lv[2, t2] = -ut2, 1.0
    lv[3, 0], lv[3, n] = -(un - c) / (2.0 * c), 1.0 / (2.0 * c)
    lv[4, 4] = 1.0
    return lam, r, lv


def flux_llf(eos: Eos, ql: np.ndarray, qr: np.ndarray, d: int = 0) -> np.ndarray:
    n = DIRECTION[d][0]
    a = max(abs(ql[n] / ql[0]) + eos.c(ql[0]), abs(qr[n] / qr[0]) + eos.c(qr[0]))
    return 0.5 * (flux(eos, ql, d) + flux(eos, qr, d)) - 0.5 * a * (qr - ql)


def signal_speeds(eos: Eos, ql: np.ndarray, qr: np.ndarray, d: int) -> tuple[float, float]:
    n = DIRECTION[d][0]
    ul, ur = ql[n] / ql[0], qr[n] / qr[0]
    return min(ul - eos.c(ql[0]), ur - eos.c(qr[0])), max(ul + eos.c(ql[0]), ur + eos.c(qr[0]))


def flux_hll(eos: Eos, ql: np.ndarray, qr: np.ndarray, d: int = 0) -> np.ndarray:
    sl, sr = signal_speeds(eos, ql, qr, d)
    if sl >= 0.0:
        return flux(eos, ql, d)
    if sr <= 0.0:
        return flux(eos, qr, d)
    return (sr * flux(eos, ql, d) - sl * flux(eos, qr, d) + sl * sr * (qr - ql)) / (sr - sl)


def flux_hllc(eos: Eos, ql: np.ndarray, qr: np.ndarray, d: int = 0) -> np.ndarray:
    """Isothermal HLLC: HLL for mass and normal momentum; the transverse momenta are the HLL mass flux times the
    upwind transverse velocity. Without an energy equation the contact carries only the tangential velocities, so
    resolving it needs no star pressure or density (the attribution to the hydro limit of an isothermal HLLD is to be
    checked against Mignone 2007 in P5)."""
    _, t1, t2 = DIRECTION[d]
    f = flux_hll(eos, ql, qr, d)
    up = ql if f[0] >= 0.0 else qr
    f[t1] = f[0] * up[t1] / up[0]
    f[t2] = f[0] * up[t2] / up[0]
    return f


def jeans_omega2(eos: Eos, rho0: float, k: float, grav: float) -> float:
    """omega^2 = c^2 k^2 - 4 pi G rho_0, c^2 = P'(rho_0): oscillation above k_J, growth (omega^2 < 0) below."""
    return eos.c2(rho0) * k * k - 4.0 * math.pi * grav * rho0


def jeans_wavenumber(eos: Eos, rho0: float, grav: float) -> float:
    return math.sqrt(4.0 * math.pi * grav * rho0 / eos.c2(rho0))


def linear_wave(eos: Eos, rho0: float, amplitude: float, k: float, x: np.ndarray, t: float) -> tuple[np.ndarray, ...]:
    """The right-going acoustic wave of small amplitude: rho, u (linear theory, valid for t << 1 / (amplitude k c))."""
    c = eos.c(rho0)
    phase = np.sin(k * (x - c * t))
    return rho0 * (1.0 + amplitude * phase), amplitude * c * phase


# --------------------------------------------------------------------------------------------------------------------
# self-test
# --------------------------------------------------------------------------------------------------------------------
POLY = Eos(c0=1.0, breaks=(2.0, 8.0), exponents=(1.0, 1.4, 5.0 / 3.0))


def problems() -> list[tuple[str, Eos, State, State]]:
    iso = isothermal(1.0)
    return [
        ("iso sod-like", iso, State(1.0, 0.0, 0.3), State(0.125, 0.0, -0.2)),
        ("iso collision", iso, State(1.0, 3.0), State(1.0, -3.0)),
        ("iso expansion", iso, State(1.0, -1.0, 0.0, 0.5), State(1.0, 1.0)),
        ("iso strong", iso, State(100.0, 0.0), State(1.0, 0.0)),
        ("poly across breaks", POLY, State(10.0, 0.0), State(0.5, 0.0, 0.7)),
        ("poly collision", POLY, State(1.0, 4.0), State(1.5, -2.0)),
    ]


def check(name: str, ok: bool, detail: str) -> bool:
    print(f"{'PASS' if ok else 'FAIL'} [{name}] {detail}")
    return ok


def check_waves() -> bool:
    worst_rh, worst_fan, worst_tail = 0.0, 0.0, 0.0
    for _, eos, left, right in problems():
        sol = Solution(eos, left, right)
        for side, s in (("left", left), ("right", right)):
            if sol.rho_star > s.rho:  # a shock: mass and momentum flux continuous in its frame
                sp = sol.shock_speed(side)
                m1, m2 = s.rho * (s.u - sp), sol.rho_star * (sol.u_star - sp)
                p1 = eos.p(s.rho) + s.rho * (s.u - sp) ** 2
                p2 = eos.p(sol.rho_star) + sol.rho_star * (sol.u_star - sp) ** 2
                worst_rh = max(worst_rh, abs(m1 - m2) / abs(m1), abs(p1 - p2) / abs(p1))
            else:  # a rarefaction: the sampled state lies on its own characteristic, u - sign c = xi (self-similarity),
                # and the fan joins the star state at its tail (the invariant itself is how sample() builds u)
                sign = 1.0 if side == "left" else -1.0
                head = s.u - sign * eos.c(s.rho)
                tail = sol.u_star - sign * eos.c(sol.rho_star)
                for xi in np.linspace(min(head, tail), max(head, tail), 41)[1:-1]:
                    st = sol.sample(float(xi))
                    err = abs(st.u - sign * eos.c(st.rho) - xi) / max(1.0, abs(xi))
                    i = eos.segment(st.rho)  # a plateau at a convex break: xi between its two one-sided speeds
                    if i > 0 and abs(st.rho - eos.breaks[i - 1]) <= TOL_FAN * st.rho:
                        lo_speed = st.u - sign * eos.c_in(st.rho, i - 1)
                        hi_speed = st.u - sign * eos.c_in(st.rho, i)
                        if min(lo_speed, hi_speed) - TOL_FAN <= xi <= max(lo_speed, hi_speed) + TOL_FAN:
                            err = 0.0
                    worst_fan = max(worst_fan, err)
                eps = 1.0e-9 * max(1.0, abs(tail - head))
                inner = sol.sample(tail - math.copysign(eps, tail - head))
                worst_tail = max(worst_tail, abs(inner.rho - sol.rho_star) / sol.rho_star)
    ok = check("rh", worst_rh < TOL, f"worst Rankine-Hugoniot residual {worst_rh:.1e}")
    return (
        check(
            "fan",
            worst_fan < TOL_FAN and worst_tail < TOL_TAIL,
            f"off-characteristic {worst_fan:.1e} (plateaus at convex breaks accepted), tail {worst_tail:.1e}",
        )
        and ok
    )


def check_isothermal() -> bool:
    iso = isothermal(1.0)
    worst = 0.0
    for mach in (1.5, 3.0, 10.0, 50.0):
        sol = Solution(iso, State(1.0, mach), State(1.0, -mach))
        sp = sol.shock_speed("left")
        worst = max(worst, abs(sol.rho_star - (mach - sp) ** 2) / sol.rho_star)  # c0 = 1: M_1 = u_1 - S
    ok = check("ratio", worst < TOL, f"isothermal compression against M_1^2, worst {worst:.1e} (Mach 1.5 to 50)")
    c0, worst = 1.7, 0.0
    iso = isothermal(c0)
    for rho_k in (0.1, 1.0, 9.0):
        for rho_s in (0.05, 0.5, 2.0, 30.0):
            exact = c0 * (rho_s - rho_k) / math.sqrt(rho_s * rho_k) if rho_s > rho_k else c0 * math.log(rho_s / rho_k)
            worst = max(worst, abs(wave_function(iso, State(rho_k, 0.0), rho_s) - exact) / max(1.0, abs(exact)))
    return (
        check("closed", worst < TOL, f"isothermal wave functions against their closed forms, worst {worst:.1e}") and ok
    )


def check_invariant() -> bool:
    xg, wg = np.polynomial.legendre.leggauss(GAUSS_POINTS)
    worst = 0.0
    for a, b in ((0.3, 1.9), (1.0, 5.0), (3.0, 20.0), (0.5, 30.0)):
        cuts = [a] + [br for br in POLY.breaks if a < br < b] + [b]
        quad = 0.0
        for lo, hi in zip(cuts[:-1], cuts[1:], strict=True):
            x = 0.5 * (hi - lo) * xg + 0.5 * (hi + lo)
            quad += 0.5 * (hi - lo) * sum(wi * POLY.c(xi) / xi for xi, wi in zip(x, wg, strict=True))
        worst = max(worst, abs(POLY.phi(b) - POLY.phi(a) - quad) / abs(quad))
    ok = check("phi", worst < TOL, f"closed-form invariant against quadrature, worst {worst:.1e}")
    worst = 0.0
    for _, eos, left, right in problems():
        a = Solution(eos, left, right)
        b = Solution(eos, State(right.rho, -right.u, right.v, right.w), State(left.rho, -left.u, left.v, left.w))
        worst = max(worst, abs(a.rho_star - b.rho_star) / a.rho_star, abs(a.u_star + b.u_star))
    ok = check("mirror", worst < TOL, f"mirrored problems, worst {worst:.1e}") and ok
    gamma = 1.4
    kpoly = Eos(c0=1.0, exponents=(gamma,))
    left, right = State(1.0, -0.5), State(0.8, 0.6)  # a two-rarefaction problem: both waves isentropic
    sol = Solution(kpoly, left, right)
    al, ar = kpoly.c(left.rho), kpoly.c(right.rho)
    pl, pr = kpoly.p(left.rho), kpoly.p(right.rho)
    z = (gamma - 1.0) / (2.0 * gamma)
    p_star = ((al + ar - 0.5 * (gamma - 1.0) * (right.u - left.u)) / (al / pl**z + ar / pr**z)) ** (1.0 / z)
    rel = abs(kpoly.p(sol.rho_star) - p_star) / p_star
    two_rarefactions = sol.rho_star < min(left.rho, right.rho)
    return (
        check(
            "toro",
            two_rarefactions and rel < TOL_TORO,
            f"two-rarefaction star pressure against Toro eq. 4.46, {rel:.1e}",
        )
        and ok
    )


def check_linear_algebra() -> bool:
    rng = np.random.default_rng(100)
    worst_i, worst_a = 0.0, 0.0
    for _ in range(50):
        q = conserved(State(rng.uniform(0.2, 20.0), *rng.uniform(-3.0, 3.0, 3)))
        for d in range(3):
            lam, r, lv = eigensystem(POLY, q, d)
            worst_i = max(worst_i, float(np.abs(lv @ r - np.eye(5)).max()))
            a = jacobian(POLY, q, d)
            worst_a = max(worst_a, float(np.abs(r @ np.diag(lam) @ lv - a).max() / max(1.0, np.abs(a).max())))
    ok = check("eigen", worst_i < TOL and worst_a < TOL, f"L R - I {worst_i:.1e}, R lambda L - A {worst_a:.1e}")
    worst = 0.0
    for _ in range(20):
        q = conserved(State(rng.uniform(0.2, 20.0), *rng.uniform(-3.0, 3.0, 3)))
        for d in range(3):
            for f in (flux_llf, flux_hll, flux_hllc):
                worst = max(worst, float(np.abs(f(POLY, q, q, d) - flux(POLY, q, d)).max()))
    iso = isothermal(1.0)
    ql, qr = conserved(State(1.0, 0.0, 0.7, -0.2)), conserved(State(1.0, 0.0, -0.4, 0.9))
    contact = float(np.abs(flux_hllc(iso, ql, qr) - flux(iso, ql)).max())  # zero mass flux: the exact flux (P only)
    return (
        check("flux", worst < TOL and contact < TOL, f"consistency {worst:.1e}, isolated contact {contact:.1e}") and ok
    )


def check_tables() -> bool:
    grav, rho0 = 0.7, 2.5
    kj = jeans_wavenumber(POLY, rho0, grav)
    marginal = abs(jeans_omega2(POLY, rho0, kj, grav)) / (4.0 * math.pi * grav * rho0)
    sides = jeans_omega2(POLY, rho0, 0.5 * kj, grav) < 0.0 < jeans_omega2(POLY, rho0, 2.0 * kj, grav)
    ok = check(
        "jeans", marginal < TOL and sides, f"omega^2(k_J) {marginal:.1e} relative; growth below k_J, oscillation above"
    )
    mf = Eos(c0=1.0, breaks=MF2021_BREAKS, exponents=MF2021_EXPONENTS, rho_ref=1.0e-20)
    worst_p, worst_g = 0.0, 0.0
    for i, b in enumerate(MF2021_BREAKS):
        below = mf.k[i] * b ** MF2021_EXPONENTS[i]
        worst_p = max(worst_p, abs(mf.p(b) - below) / below)
    for i, g in enumerate(MF2021_EXPONENTS):
        rho = 1.5 * (MF2021_BREAKS[i - 1] if i > 0 else 1.0e-20)
        worst_g = max(worst_g, abs(mf.c2(rho) / (mf.p(rho) / rho) - g))
    ok = (
        check(
            "default",
            worst_p < TOL and worst_g < TOL,
            f"Mathew & Federrath 2021: P continuity {worst_p:.1e}, c^2 / (P / rho) - G {worst_g:.1e}",
        )
        and ok
    )
    try:
        Solution(mf, State(1.0e-7, 0.0), State(1.0e-8, 0.0))
        refused = False
    except ValueError:
        refused = True
    flagged = mf.nonconvex_breaks() == [3.84e-8] and not POLY.nonconvex_breaks()
    return (
        check(
            "convex",
            flagged and refused,
            f"non-convex breaks of the default table {mf.nonconvex_breaks()} (G 1.4 -> 1.1), refused by the "
            "exact solver; the test table is convex",
        )
        and ok
    )


def selftest() -> bool:
    results = [check_waves(), check_isothermal(), check_invariant(), check_linear_algebra(), check_tables()]
    return all(results)


# --------------------------------------------------------------------------------------------------------------------
# runs (P1): the FLUME checkpoints against the exact solutions
# --------------------------------------------------------------------------------------------------------------------
TOL_PLATEAU = 1.0e-3  # star plateau mean density against rho*, relative
TOL_WAVE_CELLS = 1.5  # shock position (mid-density crossing) against the exact one, in cells
TOL_CONTACT_CELLS = 2.0  # contact position (mid transverse-velocity crossing), in cells
TOL_MASS = 1.0e-13  # total mass of a periodic box, relative drift
ORDER_MIN = 4.5  # observed order of the finest pair of the acoustic wave (WENO-5, SSP-54 at CFL 0.1)
EDGE_CELLS = 5  # cells kept away from every wave edge when the plateau is measured


def run_ngc(work: Path) -> int:
    for line in (work / "input.ini").read_text().splitlines():
        if line.split("=")[0].strip() == "ngc":
            return int(line.split("=")[1].split(";")[0])
    raise ValueError(f"{work}: no [grid] ngc")


def checkpoints(work: Path) -> list[int]:
    return sorted({int(p.name.rsplit("-proc", 1)[0].rsplit("-", 1)[1]) for p in work.glob("*-proc*.h5")})


def checkpoint_line(work: Path, it: int) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Cell centres along x and r, ru, rv, rE along x of a checkpoint (the first cell in y and z)."""
    import h5py  # noqa: PLC0415 (only the run comparisons need it)

    ngc = run_ngc(work)
    cells: dict[float, list[float]] = {}
    names = ("r", "ru", "rv", "rE")
    for path in sorted(work.glob(f"*-{it:09d}-proc*.h5")):
        with h5py.File(path, "r") as h5:
            for blk in sorted({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")}):
                origin, dxyz = h5[f"{blk}-origin"][()][::-1], h5[f"{blk}-dxdydz"][()][::-1]
                q = [h5[f"{blk}-{v}"][()][ngc:-ngc, ngc:-ngc, ngc:-ngc][0, 0, :] for v in names]
                for i in range(q[0].size):
                    cells[round(float(origin[0] + (i + ngc + 0.5) * dxyz[0]), 10)] = [float(a[i]) for a in q]
    x = np.array(sorted(cells))
    return x, {v: np.array([cells[c][k] for c in x]) for k, v in enumerate(names)}


def history(work: Path) -> tuple[list[str], np.ndarray]:
    path = next(work.glob("*conservation_history*"))
    lines = path.read_text().splitlines()
    names = [n.strip('"') for n in lines[0].split("=", 1)[1].split()]
    return names, np.array([[float(v) for v in ln.split()] for ln in lines[1:] if ln.strip()])


def checkpoint_time(work: Path, it: int) -> float:
    names, rows = history(work)
    row = rows[rows[:, names.index("it")] == it]
    assert row.size, f"{work}: no history row for it = {it}"
    return float(row[0, names.index("time")])


def crossing(x: np.ndarray, y: np.ndarray, level: float, lo: float, hi: float) -> float:
    """The first x in [lo, hi] where y crosses level (linear interpolation between the bracketing cells)."""
    for i in range(x.size - 1):
        if lo <= x[i] <= hi and (y[i] - level) * (y[i + 1] - level) <= 0.0 and y[i] != y[i + 1]:
            return float(x[i] + (level - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i]))
    return math.nan


def compare_riemann(work: Path, case: str, c0: float, mach: float) -> bool:
    """A shock tube or a collision against the exact solution: plateau, waves, inert slot, L1."""
    eos = isothermal(c0)
    if case == "tube":
        left, right = State(1.0, 0.0, 0.3), State(0.125, 0.0, -0.2)
    else:
        left, right = State(1.0, mach * c0), State(1.0, -mach * c0)
    sol = Solution(eos, left, right)
    it = checkpoints(work)[-1]
    t = checkpoint_time(work, it)
    x, q = checkpoint_line(work, it)
    dx = float(x[1] - x[0])
    rho_ex = np.array([sol.sample((xi - 0.5) / t).rho for xi in x])
    l1 = float(np.mean(np.abs(q["r"] - rho_ex)) / np.mean(rho_ex))
    edges = [0.5 + t * s for s in sol_edges(sol)]
    inner = [
        xi
        for xi in x
        if all(abs(xi - e) > EDGE_CELLS * dx for e in edges)
        and abs(sol.sample((xi - 0.5) / t).rho - sol.rho_star) == 0.0
    ]
    mask = np.isin(x, inner)
    # the plateau mean: a strong shock leaves a decaying, zero-mean oscillation train behind it (WENO, Arora & Roe
    # 1997), reported by its largest deviation, which the mean does not see
    plateau = float(abs(np.mean(q["r"][mask]) - sol.rho_star) / sol.rho_star) if mask.any() else math.inf
    ripple = float(np.max(np.abs(q["r"][mask] - sol.rho_star)) / sol.rho_star) if mask.any() else math.inf
    ok = check(
        "plateau",
        plateau < TOL_PLATEAU,
        f"{work.name}: star density {sol.rho_star:.6f}, plateau mean error {plateau:.1e} ({mask.sum()} cells, "
        f"largest deviation {ripple:.1e})",
    )
    for side in ("left", "right"):
        s = left if side == "left" else right
        if sol.rho_star <= s.rho:
            continue
        xs = 0.5 + t * sol.shock_speed(side)
        level = 0.5 * (sol.rho_star + s.rho)
        lo, hi = (0.0, 0.5 + t * sol.u_star) if side == "left" else (0.5 + t * sol.u_star, 1.0)
        xm = crossing(x, q["r"], level, lo, hi)
        err = abs(xm - xs) / dx
        ok = (
            check(
                "shock",
                err < TOL_WAVE_CELLS,
                f"{work.name}: {side} shock at {xm:.5f}, exact {xs:.5f} ({err:.2f} cells)",
            )
            and ok
        )
    if case == "tube":
        xc = 0.5 + t * sol.u_star
        v = q["rv"] / q["r"]
        xm = crossing(x, v, 0.5 * (left.v + right.v), xc - 0.1, xc + 0.1)
        err = abs(xm - xc) / dx
        ok = (
            check(
                "contact",
                err < TOL_CONTACT_CELLS,
                f"{work.name}: contact at {xm:.5f}, exact {xc:.5f} ({err:.2f} cells)",
            )
            and ok
        )
    inert = float(np.max(np.abs(q["rE"])))
    ok = check("inert", inert == 0.0, f"{work.name}: max |rE| {inert:.1e} (the energy slot is inert)") and ok
    print(f"INFO [l1] {work.name}: L1 density error {l1:.3e} relative ({x.size} cells)")
    if case == "tube":
        v_ex = np.array([sol.sample((xi - 0.5) / t).v for xi in x])
        l1v = float(np.mean(np.abs(q["rv"] / q["r"] - v_ex)))
        print(f"INFO [l1] {work.name}: L1 transverse-velocity error {l1v:.3e} (the contact)")
    return ok


def sol_edges(sol: Solution) -> list[float]:
    """The speeds of every wave edge: shocks, fan heads and tails, the contact (when it carries a jump)."""
    eos = sol.eos
    out = [sol.u_star] if (sol.left.v, sol.left.w) != (sol.right.v, sol.right.w) else []  # a contact only if it jumps
    for side, s, sign in (("left", sol.left, 1.0), ("right", sol.right, -1.0)):
        if sol.rho_star > s.rho:
            out.append(sol.shock_speed(side))
        else:
            out += [s.u - sign * eos.c(s.rho), sol.u_star - sign * eos.c(sol.rho_star)]
    return out


def compare_wave(works: list[Path], c0: float, amplitude: float) -> bool:
    """The acoustic wave after one period against its initial profile: the L1 error at each resolution, the orders."""
    rows = []
    for work in works:
        it = checkpoints(work)[-1]
        x, q = checkpoint_line(work, it)
        t = checkpoint_time(work, it)
        rho_ex, _ = linear_wave(isothermal(c0), 1.0, amplitude, 2.0 * math.pi, x, t)
        rows.append((x.size, float(np.mean(np.abs(q["r"] - rho_ex)))))
    ok = True
    for (n0, e0), (n1, e1) in zip(rows[:-1], rows[1:], strict=True):
        order = math.log(e0 / e1) / math.log(n1 / n0)
        print(f"INFO [order] N {n0} -> {n1}: L1 {e0:.3e} -> {e1:.3e}, order {order:.2f}")
    (n0, e0), (n1, e1) = rows[-2], rows[-1]
    finest = math.log(e0 / e1) / math.log(n1 / n0)
    return check("order", finest > ORDER_MIN, f"finest pair N {n0} -> {n1}: order {finest:.2f} (> {ORDER_MIN})") and ok


def compare_inert(work: Path) -> bool:
    """Every checkpoint's energy slot exactly 0; the total mass of a periodic box constant (TOL_MASS)."""
    worst = max(float(np.max(np.abs(checkpoint_line(work, it)[1]["rE"]))) for it in checkpoints(work))
    names, rows = history(work)
    mass = rows[:, names.index("int_r")]
    drift = float(np.max(np.abs(mass - mass[0])) / abs(mass[0]))
    ok = check("inert", worst == 0.0, f"{work.name}: max |rE| over every checkpoint {worst:.1e}")
    return check("mass", drift < TOL_MASS, f"{work.name}: total mass drift {drift:.1e} ({len(mass)} rows)") and ok


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest", help="check the oracle against properties independent of its construction")
    p = sub.add_parser("riemann", help="a shock tube or collision run against the exact solution")
    p.add_argument("work", type=Path)
    p.add_argument("--case", choices=("tube", "collide"), default="tube")
    p.add_argument("--c0", type=float, default=1.0)
    p.add_argument("--mach", type=float, default=5.0)
    p = sub.add_parser("wave", help="acoustic-wave runs (coarse to fine): L1 errors and orders")
    p.add_argument("works", type=Path, nargs="+")
    p.add_argument("--c0", type=float, default=1.0)
    p.add_argument("--amplitude", type=float, default=1.0e-6)
    p = sub.add_parser("inert", help="the energy slot of every checkpoint and the total mass of a run")
    p.add_argument("work", type=Path)
    args = ap.parse_args()
    if args.cmd == "selftest":
        return 0 if selftest() else 1
    if args.cmd == "riemann":
        return 0 if compare_riemann(args.work, args.case, args.c0, args.mach) else 1
    if args.cmd == "wave":
        return 0 if compare_wave(args.works, args.c0, args.amplitude) else 1
    if args.cmd == "inert":
        return 0 if compare_inert(args.work) else 1
    return 2


if __name__ == "__main__":
    sys.exit(main())
