#!/usr/bin/env python3
"""Oracle of the AS runs (#68 A3): a shock across a static 2:1 face, against the exact solution.

Why: a seam can fail a shock in three ways a run-to-run comparison does not see: conservation (the shock lands in the
wrong place), a reflection (acoustic waves left in the state behind the shock) and an entropy glitch (a density
spike, at constant pressure, advected away from the face by the flow). At t = 0.2 the tracked right shock (Sod: the
shock; RJ2a: the fast shock) is past the face at x = 0.75, which then sits in the constant state behind the shock, the
plateau, so the exact solution there is one state and each failure has its own number:

* mass: the relative defect of the mass integral against its exact value, the initial mass plus t times the net
  inflow at the two ends (their states never change by t = 0.2);
* shock: the mass-equivalent shock position (the position of the step that holds the mass of the cells within W of
  the exact shock, so a smeared profile counts as its integral) minus the exact one, in cells of the local spacing;
  W = --margin, a length (default 0.03, 6 cells of the Sod N = 200 grid) the same for every run: a margin in cells of
  the run would shrink on the fine controls, whose wider-in-cells contact would then enter the plateau;
* plateau: the cells between the left wave of the plateau (Sod: the contact; RJ2a: the rotational discontinuity)
  and the shock, W away from both: max |q - q2| / |q1 - q2| for the density, the normal velocity and the pressure (q1
  ahead of the shock, q2 behind it), max |s / s2 - 1| of the entropy function s = p / rho^gamma, and its mean;
* band: the same maxima over the cells within ngc cells (of their own spacing) of the face, the seam-local part;
* spread: the largest relative difference among the transverse copies of a column (round-off, see load).

Each run is reported; a run given as --control (the uniform run at the spacing of the region that holds the shock at
t = 0.2) is reported beside it, with the ratios refined / control.

Usage:
    amr_shock_oracle.py <work> [<work> ...] [--control <work>] [--seam X] [--max METRIC:VALUE ...]
    --max plateau_p:2e-3 fails unless that metric of every positional run is at most VALUE (shock is |shock|)
"""

from __future__ import annotations

import argparse
import configparser
import math
import sys
from pathlib import Path

import h5py
import numpy as np

VERIF = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(VERIF / "sod"))
sys.path.insert(0, str(VERIF / "mhd" / "rj2a"))
import rj2a_oracle  # noqa: E402
from sod_oracle import exact_riemann_density, star_state  # noqa: E402

EULER = ("r", "ru", "rv", "rw", "rE")
MHD = EULER + ("bx", "by", "bz")
METRICS = ("spread", "mass", "shock", "l1_r", "plateau_r", "plateau_u", "plateau_p", "plateau_s", "plateau_s_mean",
           "band_r", "band_u", "band_p", "band_s")


def read_ini(work: Path) -> configparser.ConfigParser:
    """Return the input of a run (its only .ini file)."""
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.read(next(work.glob("*.ini")))
    return ini


def load(work: Path, names: tuple[str, ...], ngc: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Return (centres, spacings, conservative variables [v, n], spread) of the last checkpoint along x.

    The transverse copies of a column (and the blocks holding the same column) agree only to round-off on a refined
    grid: the 2:1 restriction sums the fine cells of each coarse cell in an order that depends on the transverse
    position. The first copy is kept; spread is the largest relative difference among the copies, reported (it must
    stay at round-off: a transverse dependence of the solution is a defect, M2-P7b).
    """
    files = sorted(p for p in work.glob("*-proc*.h5") if "restart" not in p.name)
    if not files:
        sys.exit(f"amr_shock_oracle: no checkpoint in {work}")
    last = max(int(p.name.split("-")[-2]) for p in files)
    cells: dict[float, tuple[float, np.ndarray]] = {}
    spread = 0.0
    for path in (p for p in files if int(p.name.split("-")[-2]) == last):
        with h5py.File(path, "r") as h5:
            for blk in sorted({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")}):
                # geometry and datasets are stored (z, y, x) with ghosts: reverse to (x, y, z), strip the ghosts
                origin, dxyz = h5[f"{blk}-origin"][()][::-1], h5[f"{blk}-dxdydz"][()][::-1]
                q = np.stack([h5[f"{blk}-{v}"][()].transpose(2, 1, 0)[ngc:-ngc, ngc:-ngc, ngc:-ngc] for v in names])
                qa = q.reshape(len(names), q.shape[1], -1)
                scale = float(np.max(np.abs(qa)))  # one scale: the transverse momenta are zero
                spread = max(spread, float(np.max(np.abs(qa - qa[:, :, :1]))) / scale)
                centres = origin[0] + (np.arange(q.shape[1]) + ngc + 0.5) * dxyz[0]
                for c, col in zip(centres, qa[:, :, 0].T, strict=True):
                    if c in cells:
                        spread = max(spread, float(np.max(np.abs(cells[c][1] - col))) / scale)
                    else:
                        cells[c] = (float(dxyz[0]), col)
    xs = np.array(sorted(cells))
    return xs, np.array([cells[c][0] for c in xs]), np.array([cells[c][1] for c in xs]).T, spread


def primitives(q: np.ndarray, gamma: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return density, x velocity and pressure of conservative columns (Euler or MHD, by their count)."""
    rho, u = q[0], q[1] / q[0]
    kinetic = 0.5 * (q[1] ** 2 + q[2] ** 2 + q[3] ** 2) / q[0]
    magnetic = 0.5 * np.sum(q[5:8] ** 2, axis=0) if q.shape[0] >= 8 else 0.0
    return rho, u, (gamma - 1.0) * (q[4] - kinetic - magnetic)


class Problem:
    """The exact solution of a run: density, the tracked shock, the plateau behind it and the end fluxes."""

    def __init__(self, ini: configparser.ConfigParser) -> None:
        self.gamma = float(ini["physics"]["cp"]) / float(ini["physics"]["cv"])
        self.t = float(ini["time"]["time_max"])
        self.mhd = ini["physics"]["physical_model"].strip().startswith("mhd")
        r1, r2 = ini["initial_conditions_region_1"], ini["initial_conditions_region_2"]
        left = (float(r1["r"]), float(r1["u"]), float(r1["p"]))
        right = (float(r2["r"]), float(r2["u"]), float(r2["p"]))
        # mass: initial integral on [0, 1] (discontinuity at 0.5) plus the inflow at x = 0 minus the outflow at x = 1
        self.mass = 0.5 * (left[0] + right[0]) + self.t * (left[0] * left[1] - right[0] * right[1])
        if self.mhd:
            if abs(self.t - rj2a_oracle.T) > 1e-12:
                sys.exit("amr_shock_oracle: the RJ2a exact solution is tabulated at t = 0.2 only")
            # the wave positions of rj2a_oracle.exact (Dai & Woodward 1994): fast shock, rotational discontinuity
            self.x_shock = rj2a_oracle.X0 + 2.2638 * self.t
            self.x_left = rj2a_oracle.X0 + (0.53432 + 1.0 / math.sqrt(math.pi * 1.309)) * self.t
            ahead = rj2a_oracle.exact(np.array([self.x_shock + 1e-9]))
            behind = rj2a_oracle.exact(np.array([self.x_shock - 1e-9]))
            self.q1 = np.array(primitives(ahead, self.gamma))[:, 0]
            self.q2 = np.array(primitives(behind, self.gamma))[:, 0]
        else:
            p, u = star_state(self.gamma, left, right)
            if p <= right[2]:
                sys.exit("amr_shock_oracle: the right wave of the Sod problem is not a shock")
            g1 = (self.gamma - 1.0) / (self.gamma + 1.0)
            ar = math.sqrt(self.gamma * right[2] / right[0])
            speed = right[1] + ar * math.sqrt((self.gamma + 1.0) / (2.0 * self.gamma) * p / right[2]
                                              + (self.gamma - 1.0) / (2.0 * self.gamma))
            self.x_shock, self.x_left = 0.5 + speed * self.t, 0.5 + u * self.t
            self.q1 = np.array(right)
            self.q2 = np.array([right[0] * (p / right[2] + g1) / (g1 * p / right[2] + 1.0), u, p])
            self.left, self.right = left, right

    def density(self, xs: np.ndarray) -> np.ndarray:
        """Return the exact density at the cell centres xs."""
        if self.mhd:
            return rj2a_oracle.exact(xs)[0]
        return exact_riemann_density(xs, self.t, 0.5, self.gamma, self.left, self.right)


def measure(work: Path, seam: float, width: float) -> dict[str, float]:
    """Return the metrics of one run."""
    ini = read_ini(work)
    prob = Problem(ini)
    ngc = int(ini["grid"]["ngc"])
    xs, h, q, spread = load(work, MHD if prob.mhd else EULER, ngc)
    rho, u, p = primitives(q, prob.gamma)
    out = {"spread": spread,
           "mass": abs(float(np.sum(rho * h)) - prob.mass) / prob.mass,
           "l1_r": float(np.sum(np.abs(rho - prob.density(xs)) * h))}
    win = np.abs(xs - prob.x_shock) <= width
    a, b = float(np.min(xs[win] - 0.5 * h[win])), float(np.max(xs[win] + 0.5 * h[win]))
    r1, r2 = prob.q1[0], prob.q2[0]
    x_num = (float(np.sum(rho[win] * h[win])) - r1 * b + r2 * a) / (r2 - r1)
    out["shock"] = (x_num - prob.x_shock) / float(h[np.argmin(np.abs(xs - prob.x_shock))])
    jumps = np.abs(prob.q1 - prob.q2)
    s2 = prob.q2[2] / prob.q2[0] ** prob.gamma
    ds = np.abs(p / rho ** prob.gamma / s2 - 1.0)
    for region, sel in (("plateau", (xs >= prob.x_left + width) & (xs <= prob.x_shock - width)),
                        ("band", np.abs(xs - seam) <= ngc * h)):
        if not np.any(sel):
            sys.exit(f"amr_shock_oracle: no cell in the {region} of {work.name}")
        for k, (name, v) in enumerate((("r", rho), ("u", u), ("p", p))):
            out[f"{region}_{name}"] = float(np.max(np.abs(v[sel] - prob.q2[k]))) / jumps[k]
        out[f"{region}_s"] = float(np.max(ds[sel]))
        if region == "plateau":
            out["plateau_s_mean"] = float(np.sum(ds[sel] * h[sel]) / np.sum(h[sel]))
    return out


def main() -> int:
    """Measure the runs, print a report, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work", type=Path, nargs="+")
    parser.add_argument("--control", type=Path, default=None, help="the uniform run to compare with")
    parser.add_argument("--seam", type=float, default=0.75, help="the 2:1 face position")
    parser.add_argument("--margin", type=float, default=0.03, help="W, the distance kept from the waves")
    parser.add_argument("--max", action="append", default=[], metavar="METRIC:VALUE")
    args = parser.parse_args()
    bounds = {}
    for spec in args.max:
        name, value = spec.split(":")
        if name not in METRICS:
            sys.exit(f"amr_shock_oracle: unknown metric {name} (accepted: {', '.join(METRICS)})")
        bounds[name] = float(value)
    control = measure(args.control, args.seam, args.margin) if args.control is not None else None
    print(f"{'run':44s} " + " ".join(f"{m:>14s}" for m in METRICS))
    if control is not None:
        print(f"{args.control.name + ' (control)':44s} " + " ".join(f"{control[m]:14.3e}" for m in METRICS))
    status = 0
    for work in args.work:
        res = measure(work, args.seam, args.margin)
        print(f"{work.name:44s} " + " ".join(f"{res[m]:14.3e}" for m in METRICS))
        if control is not None:
            ratios = [f"{res[m] / control[m]:14.2f}" if control[m] != 0.0 else f"{'-':>14s}" for m in METRICS]
            print(f"{'  / control':44s} " + " ".join(ratios))
        for name, bound in bounds.items():
            value = abs(res[name])
            if value > bound:
                status = 1
                print(f"  FAIL {work.name}: {name} {value:.3e} > {bound:.3e}")
    return status


if __name__ == "__main__":
    sys.exit(main())
