#!/usr/bin/env python3
"""Oracle of the FLUME self-gravity polytrope cases (issue #81, M8-P4, SV-7), written by make_polytrope.py.

    star_oracle.py equilibrium <work> [--tol-centre C] [--tol-mach M] [--tol-energy E]
    star_oracle.py convergence <coarse work> <fine work> [--min-reduction R]
    star_oracle.py unstable <work> [--min-departure D]

The n = 1 polytrope star of [initial_conditions] (`polytrope`), with `star_pressure_factor` 1, is an equilibrium of
the Euler equations under its own isolated gravity; the discrete scheme holds it up to its truncation error. Per
checkpoint:

- the central density: the mean density of the cells within one cell of the star's centre, against its initial
  value (the drift of the core);
- the largest Mach number inside the star (cells denser than a tenth of the central density): the spurious velocity;

and from the conservation history the drift of E + W, max |(E + W)(t) - (E + W)(0)| / |(E + W)(0)| (the walls close
the box: no boundary term). Exit 1 when the largest central-density drift exceeds `--tol-centre`, the largest Mach
number `--tol-mach`, or the E + W drift `--tol-energy`.

`convergence` (SV-7) takes the same star on two grids, the second twice as fine: the largest central-density drift and
the largest Mach number inside 0.8 a of each, and how much they shrink (coarse over fine: 2 at first order, 4 at
second). Exit 1 when either shrinks by less than `--min-reduction`. Printed, not asserted: the largest Mach number
over the whole star, set at its surface (the density's kink there, and the positivity limiter acting), which
converges slowest, and the E + W drift (with the limiter it is not exact).

The SV-8 collapse (`star_pressure_factor` < 1) is checked by `jeans_oracle.py conservation` (E + W and momentum
against dt).

The `bonnor-ebert` sphere of the isothermal gas (issue #100, P4) is read the same way: its central density
`be_density`, its "star" the sphere inscribed in the box around `be_x0, be_y0, be_z0` (radius a), the Mach number
over c0. An isothermal gas conserves no energy, so the E + W drift is printed, not asserted (`--tol-energy inf`).
In a closed cube the Lane-Emden profile is no exact equilibrium (the gas in the corners pulls through its cubic
multipoles, a few per cent of the gravity), so the sphere is checked for stability: a sub-critical one stays within
`--tol-centre` of its central density (`equilibrium`), a super-critical one departs from it by at least
`--min-departure`, in either direction (`unstable`).
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "gravity"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from atmosphere_oracle import checkpoints, interior_steps, read_ini  # noqa: E402
from jeans_oracle import history  # noqa: E402


def star_rows(work: Path) -> list[tuple[float, float, float, float]]:
    """Return (time, central density, largest Mach number in the star, largest Mach number inside 0.8 a) of every
    checkpoint."""
    ini = read_ini(work)
    sec = ini["initial_conditions"]
    bonnor_ebert = sec["type"].split(";")[0].strip() == "bonnor-ebert"  # issue #100, P4
    pre = "be" if bonnor_ebert else "star"
    centre = np.array([float(sec.get(f"{pre}_{x}0", "0.0").split(";")[0]) for x in "xyz"])
    rho_c = float(sec[f"{pre}_density"].split(";")[0])
    if bonnor_ebert:  # the inscribed sphere, every cell inside it counted; the sound speed c0
        lo = np.array([float(ini["grid"][f"emin_{x}"].split(";")[0]) for x in "xyz"])
        hi = np.array([float(ini["grid"][f"emax_{x}"].split(";")[0]) for x in "xyz"])
        radius = float(np.min(np.minimum(centre - lo, hi - centre)))
        c0 = float(ini["eos"]["sound_speed"].split(";")[0])
    else:
        radius = float(sec["star_radius"].split(";")[0])
        gamma = float(ini["physics"]["gamma"].split(";")[0])
    rows = []
    for _, (t, paths) in sorted(checkpoints(work).items()):
        core, mach, mach_in = [], 0.0, 0.0
        for c, dxyz, _, q in interior_steps(ini, paths):
            rho = q["r"]
            if np.linalg.norm(c - centre) < float(np.max(dxyz)):
                core.append(rho)
            inside = np.linalg.norm(c - centre) < radius if bonnor_ebert else rho > 0.1 * rho_c
            if inside:
                u = np.array([q["ru"], q["rv"], q["rw"]]) / rho
                if bonnor_ebert:
                    m = float(np.linalg.norm(u)) / c0
                else:
                    p = (gamma - 1.0) * (q["rE"] - 0.5 * rho * float(np.dot(u, u)))
                    m = float(np.linalg.norm(u)) / math.sqrt(gamma * p / rho)
                mach = max(mach, m)
                if np.linalg.norm(c - centre) < 0.8 * radius:
                    mach_in = max(mach_in, m)
        rows.append((t, float(np.mean(core)), mach, mach_in))
    return rows


def equilibrium(args: argparse.Namespace) -> int:
    """SV-7: the central-density drift, the spurious Mach number and the E + W drift of the equilibrium star."""
    rows = star_rows(args.work)
    rho0 = rows[0][1]
    drift_c = max(abs(r[1] / rho0 - 1.0) for r in rows)
    mach = max(r[2] for r in rows)
    h = history(args.work)
    ew = h["E_plus_W"]
    drift_e = float(np.max(np.abs(ew - ew[0])) / abs(ew[0]))
    for t, rc, m, _ in rows:
        print(f"   t = {t:7.4f}: central density / initial - 1 {rc / rho0 - 1.0:+.3e}, max Mach in the star {m:.3e}")
    ok = drift_c <= args.tol_centre and mach <= args.tol_mach and drift_e <= args.tol_energy
    print(f"   {args.work.name}: {len(rows)} checkpoints to t = {rows[-1][0]:.4f}; largest central-density drift "
          f"{drift_c:.2e}, largest Mach {mach:.2e}, E + W drift {drift_e:.2e}")
    print(f"Equilibrium: {'PASS' if ok else 'FAIL'} (central {args.tol_centre:g}, Mach {args.tol_mach:g}, "
          f"E + W {args.tol_energy:g})")
    return 0 if ok else 1


def unstable(args: argparse.Namespace) -> int:
    """The super-critical Bonnor-Ebert sphere (issue #100, P4): its central density departs from the initial one."""
    rows = star_rows(args.work)
    rho0 = rows[0][1]
    dep = max(abs(r[1] / rho0 - 1.0) for r in rows)
    last = rows[-1][1] / rho0 - 1.0
    ok = dep >= args.min_departure
    print(f"   {args.work.name}: {len(rows)} checkpoints to t = {rows[-1][0]:.4f}; largest central-density departure "
          f"{dep:.2e} (last {last:+.2e})")
    print(f"Unstable: {'PASS' if ok else 'FAIL'} (min departure {args.min_departure:g})")
    return 0 if ok else 1


def summary(work: Path) -> tuple[float, float, float, float, float]:
    """Return (last time, largest central-density drift, largest Mach number inside 0.8 a, largest Mach number in the
    star, E + W drift) of a run."""
    rows = star_rows(work)
    rho0 = rows[0][1]
    ew = history(work)["E_plus_W"]
    return (rows[-1][0], max(abs(r[1] / rho0 - 1.0) for r in rows), max(r[3] for r in rows),
            max(r[2] for r in rows), float(np.max(np.abs(ew - ew[0])) / abs(ew[0])))


def convergence(args: argparse.Namespace) -> int:
    """SV-7: the spurious motion of the equilibrium star shrinks with the grid."""
    coarse, fine = summary(args.works[0]), summary(args.works[1])
    for work, (t, dc, mi, ma, de) in zip(args.works, (coarse, fine), strict=True):
        print(f"   {work.name}: to t = {t:.4f}, largest central-density drift {dc:.3e}, largest Mach inside 0.8 a "
              f"{mi:.3e} (whole star {ma:.3e}), E + W drift {de:.2e}")
    red_c, red_m, red_s = coarse[1] / fine[1], coarse[2] / fine[2], coarse[3] / fine[3]
    ok = red_c >= args.min_reduction and red_m >= args.min_reduction
    print(f"   reduction, coarse over fine: central density {red_c:.2f}, Mach inside 0.8 a {red_m:.2f} (whole star "
          f"{red_s:.2f}, not asserted)")
    print(f"Equilibrium convergence: {'PASS' if ok else 'FAIL'} (min reduction {args.min_reduction:g})")
    return 0 if ok else 1


def main() -> int:
    """Parse the command line and run the check."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="check", required=True)
    p = sub.add_parser("equilibrium")
    p.add_argument("work", type=Path)
    p.add_argument("--tol-centre", type=float, default=1.0e-2)
    p.add_argument("--tol-mach", type=float, default=1.0e-1)
    p.add_argument("--tol-energy", type=float, default=1.0e-6)
    p = sub.add_parser("convergence")
    p.add_argument("works", type=Path, nargs=2)
    p.add_argument("--min-reduction", type=float, default=2.0)
    p = sub.add_parser("unstable")
    p.add_argument("work", type=Path)
    p.add_argument("--min-departure", type=float, default=0.5)
    args = parser.parse_args()
    return {"equilibrium": equilibrium, "convergence": convergence, "unstable": unstable}[args.check](args)


if __name__ == "__main__":
    raise SystemExit(main())
