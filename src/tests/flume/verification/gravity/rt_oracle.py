#!/usr/bin/env python3
"""The linear Rayleigh-Taylor growth rate of the FLUME `rayleigh-taylor` initial condition (issue #73, M7-G3, GV-4).

Why this oracle exists: the single-mode Rayleigh-Taylor test is compared visually in the literature (its nonlinear
stage has no converged solution), so GV-4 asserts the one quantitative property it has, the growth rate of the linear
stage, against the linearised compressible Euler equations of the exact equilibrium profile. Nothing here comes from
FLUME: the equilibrium is the closed form below, the eigenproblem is solved by shooting.

The equilibrium (gravity -g along the height h, two isothermal layers, cold and dense above):
    T(h) = a + b tanh(z), z = (h - h_0) / delta, a = (T_l + T_h) / 2, b = (T_h - T_l) / 2, T_l = p_1 / rho_1,
    T_h = T_l / R (R = density_ratio: heavy over light at equal pressure), p(h_0) = p_1,
    p(h) = p_1 exp(-g int_{h_0}^h dh' / T) with int_0^z dz' / (a + b tanh z') = (a z - b ln(a cosh z + b sinh z) +
    b ln a) / (a^2 - b^2), rho = p / T.

The perturbations ~ exp(i k x + sigma t), adiabatic (ideal gas, gamma), eliminating u (x momentum) and rho_1 (mass):
    v'  = rho g / (gamma p) v - (sigma / (gamma p) + k^2 / (sigma rho)) p_1
    p_1' = -sigma rho v + (g / sigma) (k^2 / sigma p_1 + rho' v + rho v')
with v = 0 at both ends of the column (walls, or `hydrostatic` faces: the normal velocity is odd there). For a given
sigma the system is integrated from the bottom (v = 0, p_1 = 1) by RK4; the growth rate is the largest sigma at which
v vanishes at the top (the fundamental mode, v without a node), found by a downward scan and bisection.

Usage:
    rt_oracle.py sigma <input.ini> [--steps N]       growth rate of an input of make_atmosphere.py (rt cases)
    rt_oracle.py selftest                            the closed-form profile against quadrature, and the
                                                     incompressible sharp-interface limit against
                                                     sigma^2 = g k (rho_h - rho_l) / (rho_h coth(k H_h) + rho_l coth(k H_l))
"""

from __future__ import annotations

import argparse
import configparser
import math
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Column:
    """The equilibrium column of the `rayleigh-taylor` initial condition."""

    g: float  # gravity magnitude
    gamma: float  # specific heats ratio
    rho1: float  # light layer: T_l = p_1 / rho_1
    p1: float  # pressure at the interface height
    ratio: float  # density ratio heavy / light at equal pressure, T_l / T_h
    h0: float  # interface height
    delta: float  # interface width
    lo: float  # bottom of the column
    hi: float  # top of the column

    def temperature(self, h: float) -> tuple[float, float]:
        """T(h) and dT/dh."""
        tl = self.p1 / self.rho1
        th = tl / self.ratio
        a, b = 0.5 * (tl + th), 0.5 * (th - tl)
        t = math.tanh((h - self.h0) / self.delta)
        return a + b * t, b * (1.0 - t * t) / self.delta

    def state(self, h: float) -> tuple[float, float, float]:
        """rho, p, drho/dh at the height h."""
        tl = self.p1 / self.rho1
        th = tl / self.ratio
        a, b = 0.5 * (tl + th), 0.5 * (th - tl)
        z = (h - self.h0) / self.delta
        # ln(a cosh z + b sinh z), stable for large |z|
        s = 1.0 if z >= 0.0 else -1.0
        lncosh = abs(z) + math.log(0.5 * (a + b * s) + 0.5 * (a - b * s) * math.exp(-2.0 * abs(z)))
        integral = self.delta * (a * z - b * lncosh + b * math.log(a)) / (a * a - b * b)
        p = self.p1 * math.exp(-self.g * integral)
        t, dt = self.temperature(h)
        rho = p / t
        return rho, p, -self.g * rho / t - rho * dt / t


def rhs(col: Column, k: float, sigma: float, h: float, v: float, q: float) -> tuple[float, float]:
    """The linearised system at h: (v', p_1')."""
    rho, p, drho = col.state(h)
    dv = rho * col.g / (col.gamma * p) * v - (sigma / (col.gamma * p) + k * k / (sigma * rho)) * q
    dq = -sigma * rho * v + col.g / sigma * (k * k / sigma * q + drho * v + rho * dv)
    return dv, dq


def shoot(col: Column, k: float, sigma: float, steps: int) -> float:
    """v at the top of the column, integrated from v = 0, p_1 = 1 at the bottom (normalised by |p_1| there)."""
    dh = (col.hi - col.lo) / steps
    v, q = 0.0, 1.0
    for n in range(steps):
        h = col.lo + n * dh
        k1 = rhs(col, k, sigma, h, v, q)
        k2 = rhs(col, k, sigma, h + 0.5 * dh, v + 0.5 * dh * k1[0], q + 0.5 * dh * k1[1])
        k3 = rhs(col, k, sigma, h + 0.5 * dh, v + 0.5 * dh * k2[0], q + 0.5 * dh * k2[1])
        k4 = rhs(col, k, sigma, h + dh, v + dh * k3[0], q + dh * k3[1])
        v += dh / 6.0 * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0])
        q += dh / 6.0 * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1])
        m = max(abs(v), abs(q))
        if m > 1e100:  # the solution grows like exp(k h): rescale (the sign of v is what matters)
            v, q = v / m, q / m
    return v / max(abs(q), 1e-300)


def eigenfunction(col: Column, k: float, sigma: float, steps: int = 4000) -> tuple[list[float], list[float]]:
    """The vertical velocity of the mode of rate sigma along the column, normalised to max |v| = 1."""
    dh = (col.hi - col.lo) / steps
    v, q = 0.0, 1.0
    hs, vs = [col.lo], [0.0]
    for n in range(steps):
        h = col.lo + n * dh
        k1 = rhs(col, k, sigma, h, v, q)
        k2 = rhs(col, k, sigma, h + 0.5 * dh, v + 0.5 * dh * k1[0], q + 0.5 * dh * k1[1])
        k3 = rhs(col, k, sigma, h + 0.5 * dh, v + 0.5 * dh * k2[0], q + 0.5 * dh * k2[1])
        k4 = rhs(col, k, sigma, h + dh, v + dh * k3[0], q + dh * k3[1])
        v += dh / 6.0 * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0])
        q += dh / 6.0 * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1])
        hs.append(h + dh)
        vs.append(v)
    m = max(abs(x) for x in vs)
    return hs, [x / m for x in vs]


def growth_rate(col: Column, k: float, steps: int = 4000, scan: int = 200) -> float:
    """The largest sigma with v(top) = 0: scan down from sqrt(g k) (the rate of A = 1), then bisect."""
    top = math.sqrt(col.g * k)
    prev_s, prev_v = top, shoot(col, k, top, steps)
    for n in range(1, scan + 1):
        s = top * (1.0 - n / (scan + 1))
        v = shoot(col, k, s, steps)
        if (v > 0.0) != (prev_v > 0.0):
            lo, hi, flo = s, prev_s, v
            for _ in range(50):
                mid = 0.5 * (lo + hi)
                fm = shoot(col, k, mid, steps)
                if (fm > 0.0) == (flo > 0.0):
                    lo, flo = mid, fm
                else:
                    hi = mid
            return 0.5 * (lo + hi)
        prev_s, prev_v = s, v
    raise ValueError("no unstable mode found below sqrt(g k)")


def column_of(ini_path: Path) -> tuple[Column, float]:
    """The column and the wavenumber of an input of make_atmosphere.py (gravity along -y, x periodic)."""
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None, strict=False)
    ini.read(ini_path)

    def val(sec: str, key: str) -> float:
        return float(ini[sec][key].split(";")[0])

    ic = "initial_conditions"
    col = Column(g=-val("gravity", "g_y"), gamma=val("physics", "gamma"), rho1=val("initial_conditions_region_1", "r"),
                 p1=val("initial_conditions_region_1", "p"), ratio=val(ic, "rt_density_ratio"),
                 h0=val(ic, "rt_height"), delta=val(ic, "rt_width"), lo=val("grid", "emin_y"),
                 hi=val("grid", "emax_y"))
    return col, 2.0 * math.pi / val(ic, "rt_wavelength")


def selftest() -> int:
    """Check the profile against quadrature and the solver against the incompressible sharp-interface rate."""
    status = 0
    col = Column(g=0.1, gamma=1.4, rho1=1.0, p1=1.0, ratio=2.0, h0=0.5, delta=0.05, lo=0.0, hi=1.0)
    # the closed-form p against a midpoint quadrature of dp/dh = -rho g, from h_0 to the ends
    for h_end in (0.0, 1.0):
        n = 100000
        dh = (h_end - col.h0) / n
        lnp = math.log(col.p1)
        for i in range(n):
            t, _ = col.temperature(col.h0 + (i + 0.5) * dh)
            lnp -= col.g * dh / t
        err = abs(math.exp(lnp) / col.state(h_end)[1] - 1.0)
        ok = err < 1e-10
        status |= 0 if ok else 1
        print(f"   profile p({h_end:g}), closed form against quadrature: relative error {err:.1e}  "
              f"{'PASS' if ok else 'FAIL'}")
    # incompressible (huge p_1: c -> infinity, the layers uniform) and sharp (delta -> 0): the two-layer rate
    k = 4.0 * math.pi
    rl, rh = 1.0, 2.0
    exact = math.sqrt(0.1 * k * (rh - rl) / (rh / math.tanh(k * 0.5) + rl / math.tanh(k * 0.5)))
    for delta in (4e-3, 2e-3):
        sharp = Column(g=0.1, gamma=1.4, rho1=1.0, p1=1e6, ratio=2.0, h0=0.5, delta=delta, lo=0.0, hi=1.0)
        s = growth_rate(sharp, k, steps=8000)
        err = abs(s / exact - 1.0)
        ok = err < k * delta  # the finite width lowers the rate by O(k delta)
        status |= 0 if ok else 1
        print(f"   incompressible, delta {delta:g}: sigma {s:.6f}, sharp two-layer {exact:.6f}, relative {err:.1e}  "
              f"{'PASS' if ok else 'FAIL'}")
    # control: light over heavy is stable, no growing mode
    stable = Column(g=0.1, gamma=1.4, rho1=1.0, p1=1.0, ratio=0.5, h0=0.5, delta=0.05, lo=0.0, hi=1.0)
    try:
        growth_rate(stable, k, steps=1000)
        ok = False
    except ValueError:
        ok = True
    status |= 0 if ok else 1
    print(f"   control, light over heavy: no growing mode  {'PASS' if ok else 'FAIL'}")
    # the solver converged in its steps
    s1, s2 = growth_rate(col, k, steps=2000), growth_rate(col, k, steps=4000)
    err = abs(s1 / s2 - 1.0)
    ok = err < 1e-8
    status |= 0 if ok else 1
    print(f"   delta 0.05, g 0.1: sigma {s2:.10f}, 2000 against 4000 steps {err:.1e}  {'PASS' if ok else 'FAIL'}")
    print(f"rt_oracle selftest: {'PASS' if status == 0 else 'FAIL'}")
    return status


def main() -> int:
    """Parse the command line."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="check", required=True)
    p = sub.add_parser("sigma")
    p.add_argument("ini", type=Path)
    p.add_argument("--steps", type=int, default=4000)
    sub.add_parser("selftest")
    args = parser.parse_args()
    if args.check == "selftest":
        return selftest()
    col, k = column_of(args.ini)
    print(f"{growth_rate(col, k, steps=args.steps):.10f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
