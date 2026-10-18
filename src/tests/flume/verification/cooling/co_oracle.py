#!/usr/bin/env python3
"""FLUME radiative-cooling oracle (issue #93): a uniform box against the isochoric cooling ODE.

    co_oracle.py exact <work> [--floor]      the temperature of every history row against the ODE solution, to TOL_T
                                             relative, until the floor; with --floor, the floor reached when the ODE
                                             reaches it and held (T = T_1 to TOL_FLOOR, E bitwise constant); E plus the
                                             radiated energy constant to TOL_E; with [cooling].(dt_factor), every step
                                             within dt_factor t_cool and the limit reached
    co_oracle.py entropy <work>              (euler-de) the entropy of every row on the adiabat of its energy, to TOL_S
    co_oracle.py same <work a> <work b>      the density, momentum and energy columns of two runs bitwise equal

The box is uniform (make_co.py), so the volume integrals of the history are the state times the volume and every
row is one state. The ODE is dT/dt = -(gamma - 1) rho Lambda(T) / (m^2 R): its solution is written as the time to cool
from T_0 to T, t(T) = m^2 R / ((gamma - 1) rho) int_T^T_0 dT' / Lambda(T'), by Gauss-Legendre quadrature in ln T on
each piece of the table (an integrand smooth on the piece: the result does not use Townsend's closed form), and
inverted by bisection.
"""

from __future__ import annotations

import argparse
import configparser
import sys
from pathlib import Path

import numpy as np

TOL_T = 1.0e-12
TOL_FLOOR = 1.0e-14
TOL_E = 1.0e-13
TOL_S = 1.0e-13
GAUSS_POINTS = 40


class Case:
    """The input and the history of a run."""

    def __init__(self, work: Path) -> None:
        ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
        ini.read(work / "input.ini")
        phys, cool, grid = ini["physics"], ini["cooling"], ini["grid"]
        if "gamma" in phys:
            self.gamma, self.R = float(phys["gamma"]), 1.0
        else:
            cp, cv = float(phys["cp"]), float(phys["cv"])
            self.gamma, self.R = cp / cv, cp - cv
        self.t_nodes = np.array([float(v) for v in cool["temperatures"].split()])
        self.l_nodes = np.array([float(v) for v in cool["lambdas"].split()])
        self.alpha = np.log(self.l_nodes[1:] / self.l_nodes[:-1]) / np.log(self.t_nodes[1:] / self.t_nodes[:-1])
        self.m = float(cool["particle_mass"])
        self.dt_factor = float(cool.get("dt_factor", "0"))
        volume = np.prod([float(grid[f"emax_{d}"]) - float(grid[f"emin_{d}"]) for d in "xyz"])
        hist = sorted(work.glob("*-conservation_history.dat"))[0]
        names = hist.read_text().splitlines()[0].split("=", 1)[1].replace('"', "").split()
        data = np.loadtxt(hist, skiprows=1, ndmin=2)
        self.col = {n: data[:, k] for k, n in enumerate(names)}
        self.time = self.col["time"]
        self.r = self.col["int_r"] / volume
        ke = (self.col["int_ru"] ** 2 + self.col["int_rv"] ** 2 + self.col["int_rw"] ** 2) / (2.0 * self.col["int_r"])
        self.e = (self.col["int_rE"] - ke) / volume
        self.temperature = (self.gamma - 1.0) * self.e / (self.r * self.R)
        self.rs = self.col["int_rs"] / volume if "int_rs" in self.col else None

    def lam(self, t: np.ndarray) -> np.ndarray:
        """Return Lambda(T) (the table extended above T_N by its last segment)."""
        k = np.clip(np.searchsorted(self.t_nodes, t, side="right") - 1, 0, self.t_nodes.size - 2)
        return self.l_nodes[k] * (t / self.t_nodes[k]) ** self.alpha[k]

    def inv_lambda_integral(self, lo: float, hi: float) -> float:
        """Return int_lo^hi dT / Lambda(T), lo <= hi, by Gauss-Legendre in ln T on each piece of the table."""
        x, w = np.polynomial.legendre.leggauss(GAUSS_POINTS)
        cuts = [lo, *[t for t in self.t_nodes if lo < t < hi], hi]
        total = 0.0
        for a, b in zip(cuts[:-1], cuts[1:], strict=True):
            ua, ub = np.log(a), np.log(b)
            u = 0.5 * (ub - ua) * x + 0.5 * (ub + ua)
            t = np.exp(u)
            total += 0.5 * (ub - ua) * float(np.sum(w * t / self.lam(t)))
        return total

    def cooling_time_to(self, t0: float, rho: float, t: float) -> float:
        """Return the time to cool from t0 to t at the density rho."""
        return self.m**2 * self.R / ((self.gamma - 1.0) * rho) * self.inv_lambda_integral(t, t0)

    def temperature_at(self, t0: float, rho: float, time: float) -> float:
        """Return the ODE temperature at a time (the floor T_1 after it is reached), by bisection in ln T."""
        floor = self.t_nodes[0]
        if time >= self.cooling_time_to(t0, rho, floor):
            return floor
        lo, hi = np.log(floor), np.log(t0)
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if self.cooling_time_to(t0, rho, float(np.exp(mid))) > time:
                lo = mid
            else:
                hi = mid
            if hi - lo < 1.0e-17:
                break
        return float(np.exp(0.5 * (lo + hi)))


def exact(work: Path, with_floor: bool) -> bool:
    """Check the temperature history against the ODE."""
    c = Case(work)
    ok = True
    t0, rho = float(c.temperature[0]), float(c.r[0])
    t_floor = c.cooling_time_to(t0, rho, float(c.t_nodes[0]))
    before = c.time < t_floor * (1.0 - 1.0e-9)
    oracle = np.array([c.temperature_at(t0, rho, float(t)) for t in c.time[before]])
    err = np.abs(c.temperature[before] - oracle) / oracle
    status = "PASS" if err.max() <= TOL_T else "FAIL"
    ok &= status == "PASS"
    print(f"   T(t): {before.sum()} rows before the floor (t < {t_floor:.6f}), T from {t0:.3f} to "
          f"{c.temperature[before].min():.4f}, max relative error {err.max():.2e} (bound {TOL_T:.0e})  {status}")
    e_tot = c.col["E_plus_radiated"]
    drift = float(np.abs(e_tot - e_tot[0]).max() / abs(e_tot[0]))
    status = "PASS" if drift <= TOL_E else "FAIL"
    ok &= status == "PASS"
    print(f"   E + radiated: drift {drift:.2e} (bound {TOL_E:.0e}), radiated {c.col['radiated'][-1]:.6f}  {status}")
    if with_floor:
        after = ~before
        reached = after & (np.abs(c.temperature / c.t_nodes[0] - 1.0) <= TOL_FLOOR)
        if not after.any():
            print("   floor: the run ends before the ODE reaches it  FAIL")
            return False
        first = int(np.argmax(after))
        # the step that crosses the ODE floor time ends on the floor: every row from the first one after it is at T_1
        held = bool(reached[first:].all()) and bool((c.col["int_rE"][first:] == c.col["int_rE"][first]).all())
        late = float(np.abs(c.temperature[after] / c.t_nodes[0] - 1.0).max())
        status = "PASS" if held and after.sum() > 1 else "FAIL"
        ok &= status == "PASS"
        print(f"   floor: ODE at t = {t_floor:.6f}, the run at T_1 on all {after.sum()} rows after it (from t = "
              f"{c.time[first]:.6f}), max |T/T_1 - 1| {late:.1e} (bound {TOL_FLOOR:.0e}), E bitwise constant: {held}  "
              f"{status}")
    if c.dt_factor > 0.0:
        dt = np.diff(c.time)[:-1]  # the last step is clamped to time_max
        tn, en, rn = c.temperature[:-2], c.e[:-2], c.r[:-2]
        cooling = tn > c.t_nodes[0] * (1.0 + 1.0e-12)
        t_cool = en[cooling] * c.m**2 / (rn[cooling] ** 2 * c.lam(tn[cooling]))
        ratio = dt[cooling] / (c.dt_factor * t_cool)
        status = "PASS" if ratio.max() <= 1.0 + 1.0e-12 and ratio.max() >= 1.0 - 1.0e-9 else "FAIL"
        ok &= status == "PASS"
        print(f"   dt_factor {c.dt_factor}: max dt / (dt_factor t_cool) {ratio.max():.12f} over {cooling.sum()} "
              f"cooling steps (bound 1, and reached), {np.mean(ratio > 1.0 - 1.0e-9) * 100:.0f} % of them limited  "
              f"{status}")
    return ok


def entropy(work: Path) -> bool:
    """Check that the entropy of every row is on the adiabat of its energy."""
    c = Case(work)
    if c.rs is None:
        print("   no int_rs column: not an euler-de run  FAIL")
        return False
    adiabat = (c.gamma - 1.0) * c.e * c.r ** (1.0 - c.gamma)
    err = float(np.abs(c.rs - adiabat).max() / np.abs(adiabat).max())
    status = "PASS" if err <= TOL_S else "FAIL"
    print(f"   entropy on the cooled adiabat: max relative gap {err:.2e} over {c.rs.size} rows (bound {TOL_S:.0e})  "
          f"{status}")
    return status == "PASS"


def same(work_a: Path, work_b: Path) -> bool:
    """Check that the density, momentum and energy histories of two runs are bitwise equal."""
    a, b = Case(work_a), Case(work_b)
    ok = True
    for v in ("time", "int_r", "int_ru", "int_rv", "int_rw", "int_rE", "radiated"):
        equal = a.col[v].shape == b.col[v].shape and bool((a.col[v] == b.col[v]).all())
        ok &= equal
        print(f"   {v}: {'bitwise equal' if equal else 'DIFFERENT'}")
    print(f"   {'PASS' if ok else 'FAIL'}")
    return ok


def main() -> int:
    """Run the requested check."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mode", choices=("exact", "entropy", "same"))
    ap.add_argument("works", type=Path, nargs="+")
    ap.add_argument("--floor", action="store_true", help="exact: also check the floor")
    a = ap.parse_args()
    if a.mode == "exact":
        ok = exact(a.works[0], a.floor)
    elif a.mode == "entropy":
        ok = entropy(a.works[0])
    else:
        ok = same(a.works[0], a.works[1])
    print(f"CO {a.mode}: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
