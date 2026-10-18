#!/usr/bin/env python3
"""FLUME radiative-cooling oracle (issue #93): a uniform box against the isochoric cooling ODE.

    co_oracle.py exact <work> [--floor]      the temperature of every history row against the ODE solution, to TOL_T
                                             relative, until the floor; with --floor, the floor reached when the ODE
                                             reaches it and held (T = T_1 to TOL_FLOOR, E bitwise constant); E plus the
                                             radiated energy constant to TOL_E; with [cooling].(dt_factor), every step
                                             within dt_factor t_cool and the limit reached
    co_oracle.py entropy <work>              (euler-de) the entropy of every row on the adiabat of its energy, to TOL_S
    co_oracle.py same <work a> <work b>      the density, momentum and energy columns of two runs bitwise equal
    co_oracle.py agree <cpu work> <fnl work> CPU and FNL agree: the same steps, every row's time, temperature and
                                             radiated energy within TOL_AGREE relative
    co_oracle.py shock <work>                the radiative shock (make_co.py --case shock) against its steady structure:
                                             the shock speed, the cooling-layer profile and thickness, the dense layer

The box is uniform (make_co.py), so the volume integrals of the history are the state times the volume and every
row is one state. The ODE is dT/dt = -(gamma - 1) rho Lambda(T) / (m^2 R): its solution is written as the time to cool
from T_0 to T, t(T) = m^2 R / ((gamma - 1) rho) int_T^T_0 dT' / Lambda(T'), by Gauss-Legendre quadrature in ln T on
each piece of the table (an integrand smooth on the piece: the result does not use Townsend's closed form), and
inverted by bisection.

The radiative shock: two streams rho_0, -+u_0 at the floor temperature T_1 collide; the shocked gas cools back to T_1
and piles up at rest. Taken whole, shock plus cooling layer is an isothermal jump at T_1, so in the frame of the
structure (speed V_s outward) mass and momentum give V_s (u_0 + V_s) = R T_1 and the dense layer rho_d = rho_0 (u_0 +
V_s) / V_s, with no energy equation. The cooling layer is steady in that frame: rho w = m, p + rho w^2 = P, and the
total enthalpy flux falls at the cooling rate, m dH/dxi = -n^2 Lambda(T), H = gamma R T / (gamma - 1) + w^2 / 2, from
the Rankine-Hugoniot state behind the shock (xi = 0) to T_1, where w reaches V_s (the dense layer). The profile rho(xi)
is integrated in ln T on a fine grid. Checked on the run: the outward speed of the shock (rho crossing the mean of
rho_0 and the post-shock density), fitted over the checkpoints from t = SHOCK_T_FIT, against V_s; the density of the
layer against the ODE profile aligned at the shock (TOL_PROFILE, L1 relative), and the point where rho = SHOCK_RHO_MARK
against the ODE's (one cell); the median density and temperature of the outer half of the dense layer (the gas cooled
last: the centre keeps the start-up transient of the collision) against rho_d and T_1. The planar radiative shock
is overstable for small cooling exponents (Chevalier & Imamura 1982): the case uses alpha = 1.5, where the structure
settles after t ~ 1 (at alpha = 0.5 it still oscillates at t = 1.5, the speed 17 % off).
"""

from __future__ import annotations

import argparse
import configparser
import sys
from pathlib import Path

import h5py
import numpy as np

TOL_T = 1.0e-12
TOL_FLOOR = 1.0e-14
TOL_E = 1.0e-13
TOL_S = 1.0e-13
TOL_AGREE = 1.0e-12
TOL_SPEED = 1.0e-2
TOL_DENSE = 2.0e-2
TOL_PROFILE = 5.0e-2
SHOCK_T_FIT = 1.5
SHOCK_RHO_MARK = 6.0
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


def agree(work_cpu: Path, work_fnl: Path) -> bool:
    """Check that the CPU and FNL histories agree row by row."""
    a, b = Case(work_cpu), Case(work_fnl)
    if a.time.size != b.time.size:
        print(f"   rows: CPU {a.time.size}, FNL {b.time.size}  FAIL")
        return False
    ok = True
    for name, x, y in (("time", a.time[1:], b.time[1:]), ("T", a.temperature, b.temperature),
                       ("radiated", a.col["radiated"][1:], b.col["radiated"][1:])):
        gap = float(np.abs(x - y).max() / np.abs(x).max())
        status = "PASS" if gap <= TOL_AGREE else "FAIL"
        ok &= status == "PASS"
        print(f"   {name}: max relative gap {gap:.2e} over {a.time.size} rows (bound {TOL_AGREE:.0e})  {status}")
    return ok


def checkpoint_line(work: Path, it: int, ngc: int) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Return the cell centres along x and r, ru, rE along x of a checkpoint (the first cell in y and z)."""
    cells: dict[float, list[float]] = {}
    for path in sorted(work.glob(f"*-{it:09d}-proc*.h5")):
        with h5py.File(path, "r") as h5:
            for blk in sorted({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")}):
                origin, dxyz = h5[f"{blk}-origin"][()][::-1], h5[f"{blk}-dxdydz"][()][::-1]
                q = [h5[f"{blk}-{v}"][()][ngc:-ngc, ngc:-ngc, ngc:-ngc][0, 0, :] for v in ("r", "ru", "rE")]
                for i in range(q[0].size):
                    cells[round(float(origin[0] + (i + ngc + 0.5) * dxyz[0]), 10)] = [float(a[i]) for a in q]
    x = np.array(sorted(cells))
    return x, {v: np.array([cells[c][k] for c in x]) for k, v in enumerate(("r", "ru", "rE"))}


def shock(work: Path) -> bool:
    """Check the radiative shock against its steady structure."""
    c = Case(work)
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.read(work / "input.ini")
    ngc = int(ini["grid"]["ngc"])
    reg = ini["initial_conditions_region_1"]
    rho0, u0, p0 = float(reg["r"]), float(reg["u"]), float(reg["p"])
    g, gas, t1 = c.gamma, c.R, float(c.t_nodes[0])
    vs = (-u0 + np.sqrt(u0**2 + 4.0 * gas * t1)) / 2.0
    w0 = u0 + vs
    m, pm = rho0 * w0, rho0 * w0**2 + p0
    rho_d = m / vs
    mach2 = w0**2 / (g * p0 / rho0)
    rho1 = rho0 * (g + 1.0) * mach2 / ((g - 1.0) * mach2 + 2.0)
    p1 = p0 * (2.0 * g * mach2 - (g - 1.0)) / (g + 1.0)
    # the steady cooling layer, integrated in ln T from the post-shock temperature to the floor
    temps = np.geomspace(p1 / (rho1 * gas), t1 * (1.0 + 1.0e-12), 200001)
    # the subsonic branch of p + rho w^2 = P: w = m / rho1 behind the shock, V_s at the floor
    w = (pm / m - np.sqrt((pm / m) ** 2 - 4.0 * gas * temps)) / 2.0
    h = g / (g - 1.0) * gas * temps + 0.5 * w**2
    mid = np.sqrt(temps[1:] * temps[:-1])
    rho_mid = m / ((pm / m - np.sqrt((pm / m) ** 2 - 4.0 * gas * mid)) / 2.0)
    xi = np.concatenate([[0.0], np.cumsum(-m * np.diff(h) / ((rho_mid / c.m) ** 2 * c.lam(mid)))])
    rho_ode = m / w
    print(f"   theory: V_s {vs:.6f}, dense layer rho {rho_d:.4f}, post-shock rho {rho1:.4f}, "
          f"cooling layer {xi[-1]:.4f}")
    its = sorted(int(f.name.split("-")[-2]) for f in work.glob("*-proc000000.h5"))
    time_of = dict(zip(c.col["it"].astype(int), c.time, strict=True))  # the checkpoints' times, from the history

    def front(x: np.ndarray, r: np.ndarray) -> float:
        """Return the shock position on the left half, where rho crosses the mean of rho_0 and rho1."""
        xl, rl = x[x < 0.5], r[x < 0.5]
        thr = 0.5 * (rho0 + rho1)
        i = int(np.argmax(rl > thr))
        return float(xl[i - 1] + (thr - rl[i - 1]) * (xl[i] - xl[i - 1]) / (rl[i] - rl[i - 1]))

    ok = True
    times, dist = [], []
    for it in its:
        if it == 0 or time_of[it] < SHOCK_T_FIT:
            continue
        x, q = checkpoint_line(work, it, ngc)
        times.append(time_of[it])
        dist.append(0.5 - front(x, q["r"]))
    speed = float(np.polyfit(times, dist, 1)[0])
    err = abs(speed - vs) / vs
    status = "PASS" if err <= TOL_SPEED and len(times) >= 3 else "FAIL"
    ok &= status == "PASS"
    print(f"   shock speed {speed:.6f} over {len(times)} checkpoints from t = {times[0]:.3f} to {times[-1]:.3f}: "
          f"{err:.1e} relative (bound {TOL_SPEED:.0e})  {status}")
    x, q = checkpoint_line(work, its[-1], ngc)
    r = q["r"]
    temp = (g - 1.0) * (q["rE"] - 0.5 * q["ru"] ** 2 / r) / (r * gas)
    xs = front(x, r)
    dx = float(x[1] - x[0])
    layer = (x > xs) & (x - xs < xi[-1])
    expected = np.interp(x[layer] - xs, xi, rho_ode)
    l1 = float(np.sum(np.abs(r[layer] - expected)) / np.sum(expected))
    status = "PASS" if l1 <= TOL_PROFILE else "FAIL"
    ok &= status == "PASS"
    print(f"   cooling layer: {layer.sum()} cells, density against the ODE {l1:.2e} L1 relative (bound "
          f"{TOL_PROFILE:.0e})  {status}")
    xl, rl = x[(x > xs) & (x < 0.5)], r[(x > xs) & (x < 0.5)]
    i = int(np.argmax(rl > SHOCK_RHO_MARK))
    mark = float(xl[i - 1] + (SHOCK_RHO_MARK - rl[i - 1]) * (xl[i] - xl[i - 1]) / (rl[i] - rl[i - 1])) - xs
    mark_ode = float(np.interp(SHOCK_RHO_MARK, rho_ode, xi))
    status = "PASS" if abs(mark - mark_ode) <= dx else "FAIL"
    ok &= status == "PASS"
    print(f"   rho = {SHOCK_RHO_MARK} at {mark:.5f} behind the shock, ODE {mark_ode:.5f}: "
          f"{abs(mark - mark_ode) / dx:.2f} cells (bound 1)  {status}")
    start = xs + xi[-1]
    dense = (x > start) & (x < start + 0.5 * (0.5 - start))
    rho_med, t_med = float(np.median(r[dense])), float(np.median(temp[dense]))
    status = "PASS" if abs(rho_med / rho_d - 1.0) <= TOL_DENSE and abs(t_med / t1 - 1.0) <= TOL_DENSE else "FAIL"
    ok &= status == "PASS"
    print(f"   dense layer (outer half, {dense.sum()} cells): rho {rho_med:.4f} ({rho_med / rho_d - 1.0:+.1e}), "
          f"T {t_med:.5f} ({t_med / t1 - 1.0:+.1e}) (bound {TOL_DENSE:.0e})  {status}")
    return ok


def main() -> int:
    """Run the requested check."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mode", choices=("exact", "entropy", "same", "agree", "shock"))
    ap.add_argument("works", type=Path, nargs="+")
    ap.add_argument("--floor", action="store_true", help="exact: also check the floor")
    a = ap.parse_args()
    if a.mode == "exact":
        ok = exact(a.works[0], a.floor)
    elif a.mode == "entropy":
        ok = entropy(a.works[0])
    elif a.mode == "same":
        ok = same(a.works[0], a.works[1])
    elif a.mode == "agree":
        ok = agree(a.works[0], a.works[1])
    else:
        ok = shock(a.works[0])
    print(f"CO {a.mode}: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
