#!/usr/bin/env python3
"""Write the inputs of the FLUME barotropic-gravity verification (issue #100, M13 P4): the isothermal and barotropic
gases under gravity, self-gravity, the Jeans marker and sinks.

    make_eos_gravity.py <out.ini> <case> [--eos isothermal|barotropic] [--cells N] [--cfl C] [--time-max T]
                        [--it-save N] [--it-max N] [--amplitude A] [--box] [--solver auto|fft|multigrid] [--amr L]
                        [--xi XI] [--overdensity F] [--set SECTION:KEY=VALUE ...] [--unset SECTION:KEY ...]

The gas: `--eos isothermal` (default), c0 of the case, or `--eos barotropic`, the convex test table POLY of
eos_oracle.py (c0 = 1, densities 2 8, exponents 1 1.4 5/3) on a background in its Gamma = 1.4 or 5/3 segment, where
c^2 = dP/drho differs from P / rho. The cases:

  freefall  make_eos.py's periodic box at rest (sod-x, x periodic, 2 blocks of 24 cells), uniform, under the constant
            gravity g_x = 0.5 to t = 0.4: rho u = rho g t, rho uniform, rE = 0 (barotropic: rho0 = 5);
  jeans     make_jeans.py jeans2d (the stable wave, 4x4 blocks of --cells, k = 2 pi (1, 1)) with the gas in place of
            the gamma-law one, G set so that 4 pi G rho0 = c^2 k^2 / 4 (lambda = lambda_J / 2, w = sqrt(3) c k / 2)
            and one period; isothermal c0 = 1, rho0 = 1 (the gamma-law case's numbers), barotropic rho0 = 5
            (c^2 = P'(5) = 2.0198, P / rho = 1.4427); --box refines the central box (2:1 seams, multigrid);
  jeansu    make_jeans.py jeansu2d, 4 pi G rho0 = 9/4 c^2 k^2 (lambda = 1.5 lambda_J, sigma = sqrt(5/4) c k), the growth
            to t = 0.5 / c; with --amplitude above 1e-3, t = 0.4 / c, the nonlinear collapse (the momentum, the mass);
  marker    a static periodic uniform box ([-1, 1]^3, 32^3, the sinks' box without sinks), a Jeans marker of 5.4 cells
            per Jeans length, a regrid at step 4: barotropic rho0 = 10 (Gamma = 5/3, c^2 = 3.37) with G such that
            lambda_J / h = 6 by c^2 = dP/drho (no refinement; by P / rho it would be 4.65: refined); isothermal, c0^2
            such that lambda_J / h = 5 (refined);
  shu       make_si.py shu (Federrath et al. 2010, A = 29.3, a seed sink) with the isothermal gas c0 = 0.5 (the run's
            a^2 = 0.25);
  stable    make_si.py stable (a static periodic box at 1.2 rho_thr, creation on) with c0^2 = 0.0259, the gamma = 1.4
            case's c^2: Jeans-stable by the virial thermal term (|E_grav| / 2 E_th = 0.79 with E_th = 3/2 int P dV),
            no sink may form (E_th = 0, the energy slot, would form one);
  accrete   make_si.py accrete (a sink of mass 0.1 in a dense cold periodic background) with c0^2 = 2e-3;
  collapse  the sinks' box (walls, isolated self-gravity), a supercritical Bonnor-Ebert sphere: xi = --xi (default 10)
            at the walls, c0 = 0.2, `be_overdensity` --overdensity (default 2), creation on, to t = 3: one sink at
            the centre;
  be        the same box, the sub-critical Bonnor-Ebert sphere, xi = --xi (default 3) at the walls, c0 = 1, no sinks,
            to t = 2 (3 free-fall times of the centre): an equilibrium;
  bb        the Boss-Bodenheimer cloud (Burkert & Bodenheimer 1993 in G = 1 units: R = 1, M = 1, alpha = 0.26,
            beta = 0.16, A = 0.1) in [-2, 2]^3 with an ambient of a hundredth of its density, isothermal, a Jeans
            marker (8 cells) over --amr levels, sinks with creation: the morphology, to t = 1.25 t_ff.

`--cells` sets the cells per block side (jeans) or per box side (the 3-D boxes), `--amr L` runtime AMR levels (shu,
bb). The 3-D boxes keep make_si.py's positivity limiter, for this gas the density-only one (P4): the accretion
drains the cells within the accretion radius to the threshold density. `--set`/`--unset` edit the result (the
contract cases).
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "selfgravity"))
sys.path.insert(0, str(HERE.parent / "sinks"))

import make_eos  # noqa: E402
import make_jeans  # noqa: E402
import make_si  # noqa: E402
from eos_oracle import Eos  # noqa: E402
from make_eos import set_key, unset_key  # noqa: E402

POLY = {"densities": "2.0 8.0", "exponents": "1.0 1.4 1.6666666666666667"}
POLY_EOS = Eos(c0=1.0, breaks=(2.0, 8.0), exponents=(1.0, 1.4, 5.0 / 3.0))
K2 = 8.0 * math.pi**2  # the Jeans cases' k^2: k = 2 pi (1, 1)


def drop_section(s: str, name: str) -> str:
    """Remove the section `[name]` (to the next section or the end)."""
    i = s.index(f"\n[{name}]")
    j = s.find("\n[", i + 1)
    return s[:i] + (s[j:] if j >= 0 else "\n")


def replace_ic(s: str, ic: str) -> str:
    """Replace the [initial_conditions] and [initial_conditions_region_1] sections (the others stay where they are)."""
    i = s.index("[initial_conditions]")
    j = s.find("\n[", s.index("[initial_conditions_region_1]") + 1)
    return s[:i] + ic + (s[j:] if j >= 0 else "")


def gas(s: str, eos: str, c0: float) -> str:
    """Put the gas in place of the gamma-law one: no gamma/cp/cv, no region pressure, the [eos] section."""
    for key in ("gamma", "cp", "cv"):
        s = unset_key(s, "physics", key)
    for k in (1, 2):
        if f"[initial_conditions_region_{k}]" in s:
            s = unset_key(s, f"initial_conditions_region_{k}", "p")
    s = set_key(s, "eos", "type", eos)
    s = set_key(s, "eos", "sound_speed", repr(c0))
    if eos == "barotropic":
        for key, value in POLY.items():
            s = set_key(s, "eos", key, value)
    return s


def jeans(a: argparse.Namespace) -> str:
    """The Jeans cases on the isothermal or the barotropic gas."""
    unstable = a.case == "jeansu"
    rho0 = 5.0 if a.eos == "barotropic" else 1.0
    c2 = POLY_EOS.c2(rho0) if a.eos == "barotropic" else 1.0
    t_end = a.time_max
    if t_end is None:
        if unstable:
            t_end = (0.4 if a.amplitude and a.amplitude > 1.0e-3 else 0.5) / math.sqrt(c2)
        else:
            t_end = 2.0 * math.pi / math.sqrt(0.75 * c2 * K2)
    s = make_jeans.case_ini("jeansu2d" if unstable else "jeans2d", a.cells, a.cfl, t_end, a.it_save, a.amplitude,
                            None, a.box, a.solver, None, None, 0)
    s = gas(s, a.eos, 1.0)
    s = set_key(s, "initial_conditions_region_1", "r", repr(rho0))
    grav = (2.25 if unstable else 0.25) * c2 * K2 / (4.0 * math.pi * rho0)
    return set_key(s, "gravity", "G", repr(grav))


def si_namespace(case: str, cells: int, amr: int, time_max: float | None, it_max: int) -> argparse.Namespace:
    """make_si.py's arguments of a case."""
    a = argparse.Namespace(case=case, model="euler", cells=cells, rho=None, p=None, time_max=time_max, it_max=it_max,
                           restart_save=0, restart=False, no_sinks_section=False, no_gravity=False, creation=False,
                           base_levels=1, amr=amr, overdensity=29.3, x=None, y=None, z=None, u=None, v=None, w=None,
                           mass=None)
    for key, value in make_si.CASES[case].items():
        if getattr(a, key) is None:
            setattr(a, key, value)
    return a


def sinks_box(a: argparse.Namespace, case: str, c0: float) -> str:
    """A make_si.py case with the isothermal gas (and make_si.py's positivity limiter: density-only, P4)."""
    s = make_si.case_ini(si_namespace(case, a.cells, a.amr, a.time_max, a.it_max))
    return gas(s, "isothermal", c0)


def sphere(a: argparse.Namespace, xi: float, c0: float, overdensity: float, sinks: bool) -> str:
    """The sinks' box with a Bonnor-Ebert sphere of xi at the walls (G = 1, half-width 1: r_0 = 1 / xi)."""
    s = make_si.case_ini(si_namespace("collapse", a.cells, a.amr, a.time_max, a.it_max))
    s = s[: s.index("\n[cooling]")] + s[s.index("\n[sinks]"):]  # the collapse case's cooling: not needed here
    s = gas(s, "isothermal", c0)
    rho_c = c0**2 * xi**2 / (4.0 * math.pi)
    s = replace_ic(s, "[initial_conditions]\ntype = bonnor-ebert\nregions_number = 1\namr_iterations = 0\n"
                      f"be_density = {rho_c!r}\nbe_overdensity = {overdensity!r}\n\n"
                      "[initial_conditions_region_1]\nr = 1.0\nu = 0.0\nv = 0.0\nw = 0.0\n")
    if not sinks:
        s = drop_section(s, "sinks")
    return s


def marker(a: argparse.Namespace) -> str:
    """The Jeans marker on a static periodic box: refined by c^2 = dP/drho, or not."""
    s = make_si.case_ini(si_namespace("stable", a.cells, 1, 100.0, a.it_max if a.it_max > 0 else 4))
    s = drop_section(s, "sinks")
    s = s.replace("jeans_cells = 8.0", "jeans_cells = 5.4\njeans_cells_derefine = 1000.0")
    h = 2.0 / a.cells
    if a.eos == "barotropic":
        rho0, cells = 10.0, 6.0
        grav = math.pi * POLY_EOS.c2(rho0) / (rho0 * (cells * h) ** 2)
        s = gas(s, "barotropic", 1.0)
    else:
        rho0, cells, grav = 1.0, 5.0, 1.0
        s = gas(s, "isothermal", math.sqrt((cells * h) ** 2 * grav * rho0 / math.pi))
    s = set_key(s, "initial_conditions_region_1", "r", repr(rho0))
    return set_key(s, "gravity", "G", repr(grav))


def boss_bodenheimer(a: argparse.Namespace) -> str:
    """The Boss-Bodenheimer cloud in [-2, 2]^3, isothermal, Jeans-refined, with sinks."""
    mass, radius = 1.0, 1.0
    rho0 = 3.0 * mass / (4.0 * math.pi * radius**3)
    w = 0.6 * mass**2 / radius  # |W| of the uniform sphere, G = 1
    c0 = math.sqrt(0.26 * w / (1.5 * mass))  # alpha = E_th / |W| = 0.26
    omega = math.sqrt(2.0 * 0.16 * w / (0.4 * mass * radius**2))  # beta = E_rot / |W| = 0.16
    t_ff = math.sqrt(3.0 * math.pi / (32.0 * rho0))
    s = make_si.case_ini(si_namespace("collapse", a.cells, a.amr, a.time_max or 1.25 * t_ff, a.it_max))
    s = s[: s.index("\n[cooling]")] + s[s.index("\n[sinks]"):]
    s = gas(s, "isothermal", c0)
    for x in "xyz":
        s = set_key(s, "grid", f"emin_{x}", "-2.0")
        s = set_key(s, "grid", f"emax_{x}", "2.0")
    s = replace_ic(s, "[initial_conditions]\ntype = boss-bodenheimer\nregions_number = 1\namr_iterations = 0\n"
                      f"bb_radius = {radius!r}\nbb_density = {rho0!r}\nbb_omega = {omega!r}\nbb_amplitude = 0.1\n\n"
                      f"[initial_conditions_region_1]\nr = {rho0 / 100.0!r}\nu = 0.0\nv = 0.0\nw = 0.0\n")
    return s


def freefall(a: argparse.Namespace) -> str:
    """make_eos.py's periodic box, uniform, under a constant gravity along x."""
    ns = argparse.Namespace(case="wave", ni=48, time_max=0.4, cfl=0.5, c0=1.0, mach=5.0, scheme="weno",
                            recon="characteristic", solver="hllc", rk=None, keep_gas=False, null_y=False, level=None,
                            eos=a.eos, densities=None, exponents=None, rho_ref=None,
                            rho0=5.0 if a.eos == "barotropic" else 1.0, left=None, right=None, set=[], unset=[])
    if a.eos == "barotropic":
        ns.densities, ns.exponents = POLY["densities"], POLY["exponents"]
    s = make_eos.case_ini(ns)
    s = set_key(s, "initial_conditions", "wave_amplitude", "0.0")
    return set_key(s, "gravity", "g_x", "0.5")


def case_ini(a: argparse.Namespace) -> str:
    """Return the INI of a case."""
    xi, over = a.xi, a.overdensity
    builders = {
        "freefall": lambda: freefall(a),
        "jeans": lambda: jeans(a),
        "jeansu": lambda: jeans(a),
        "marker": lambda: marker(a),
        "shu": lambda: sinks_box(a, "shu", 0.5),
        "stable": lambda: sinks_box(a, "stable", math.sqrt(0.0259)),
        "accrete": lambda: sinks_box(a, "accrete", math.sqrt(2.0e-3)),
        "collapse": lambda: sphere(a, xi or 10.0, 0.2, over or 2.0, sinks=True),
        "be": lambda: sphere(a, xi or 3.0, 1.0, over or 1.0, sinks=False),
        "bb": lambda: boss_bodenheimer(a),
    }
    s = builders[a.case]()
    if a.case not in ("jeans", "jeansu"):  # make_jeans.py writes it_save itself
        s = set_key(s, "IO", "it_save", str(a.it_save))
    if a.case in ("be", "collapse") and a.time_max is None:
        s = set_key(s, "time", "time_max", "2.0" if a.case == "be" else "3.0")
    for item in a.set:
        sk, value = item.split("=", 1)
        section, key = sk.split(":", 1)
        s = set_key(s, section, key, value)
    for item in a.unset:
        section, key = item.split(":", 1)
        s = unset_key(s, section, key)
    return s


def main() -> int:
    """Parse the command line and write the input."""
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", type=Path)
    ap.add_argument("case", choices=("freefall", "jeans", "jeansu", "marker", "shu", "stable", "accrete", "collapse",
                                     "be", "bb"))
    ap.add_argument("--eos", choices=("isothermal", "barotropic"), default="isothermal")
    ap.add_argument("--cells", type=int, default=None)
    ap.add_argument("--cfl", type=float, default=0.4)
    ap.add_argument("--time-max", type=float, default=None)
    ap.add_argument("--it-save", type=int, default=1000000)
    ap.add_argument("--it-max", type=int, default=-1)
    ap.add_argument("--amplitude", type=float, default=None)
    ap.add_argument("--box", action="store_true")
    ap.add_argument("--solver", choices=("auto", "fft", "multigrid"), default=None)
    ap.add_argument("--amr", type=int, default=0)
    ap.add_argument("--xi", type=float, default=None)
    ap.add_argument("--overdensity", type=float, default=None)
    ap.add_argument("--set", action="append", default=[])
    ap.add_argument("--unset", action="append", default=[])
    a = ap.parse_args()
    if a.cells is None:
        a.cells = 8 if a.case in ("jeans", "jeansu") else 32
    a.out.write_text(case_ini(a))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
