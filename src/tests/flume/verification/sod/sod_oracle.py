#!/usr/bin/env python3
"""FLUME verification V1 oracle: Sod shock tube along x, y, z (issue #35, section 11).

Why this oracle exists: CHASE ran only x-aligned 1-D Riemann problems, which is how its transposed characteristic
projection and its singular y/z eigenvectors went unnoticed. V1 closes that gap with two independent checks:

* accuracy: the density of each run is compared with the exact Riemann solution (Toro, "Riemann Solvers and
  Numerical Methods for Fluid Dynamics", 3rd ed., section 4.2-4.5), the L1 error must stay below a recorded bound;
* direction invariance: the x, y and z runs must be identical after the permutation of the axes, to the last bit
  (cyclic tangents make the three directional operators the same arithmetic on permuted data).

Every copy of the 1-D solution along the null (transverse) directions must also be bitwise identical.

Usage:
    sod_oracle.py <case.ini> <work-dir> [<work-dir> ...] [--l1-max L] [--tol T] [--mirror M]

Every work directory after the first is compared with the first one: bitwise by default (the three directions of one
backend), within `--tol` (maximum absolute difference) across backends, whose compilers contract differently into
FMAs.

`--mirror M` switches to the reflecting-wall double Sod (sod-wall-*.ini): instead of the exact solution, each run must be
mirror-symmetric about the centre, q(x) = S q(1 - x) with the normal momentum negated, within M. Symmetry alone does
not prove the walls reflect (extrapolation at both ends is symmetric too): `--closed C` also requires the box to be
closed, the relative drift of the mass and energy integrals of the conservation history below C.
"""

from __future__ import annotations

import argparse
import configparser
import math
import sys
from pathlib import Path

import h5py
import numpy as np


def read_ini(path: Path) -> configparser.ConfigParser:
    """Read a FLUME INI file (`;` comments, no interpolation)."""
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.read(path)
    return ini


def exact_riemann_density(x: np.ndarray, t: float, x0: float, gamma: float,
                          left: tuple[float, float, float], right: tuple[float, float, float]) -> np.ndarray:
    """Return the exact density of a 1-D Riemann problem at time `t` (Toro, chapter 4).

    Parameters
    ----------
    x : np.ndarray
        Sampling abscissae.
    t : float
        Time.
    x0 : float
        Initial discontinuity position.
    gamma : float
        Specific heats ratio.
    left, right : tuple[float, float, float]
        Primitive states (density, velocity, pressure).

    Returns
    -------
    np.ndarray
        Exact density at `x`.
    """
    rl, ul, pl = left
    rr, ur, pr = right
    al, ar = math.sqrt(gamma * pl / rl), math.sqrt(gamma * pr / rr)

    def f_k(p: float, rk: float, pk: float, ak: float) -> tuple[float, float]:
        """Pressure function of one side and its derivative (Toro, eq. 4.6-4.7 and 4.37)."""
        if p > pk:
            ak_, bk_ = 2.0 / ((gamma + 1.0) * rk), (gamma - 1.0) / (gamma + 1.0) * pk
            sq = math.sqrt(ak_ / (p + bk_))
            return (p - pk) * sq, sq * (1.0 - 0.5 * (p - pk) / (p + bk_))
        pr_ = p / pk
        return (2.0 * ak / (gamma - 1.0) * (pr_ ** ((gamma - 1.0) / (2.0 * gamma)) - 1.0),
                1.0 / (rk * ak) * pr_ ** (-(gamma + 1.0) / (2.0 * gamma)))

    p = 0.5 * (pl + pr)
    for _ in range(100):
        fl, dfl = f_k(p, rl, pl, al)
        fr, dfr = f_k(p, rr, pr, ar)
        p_new = max(1e-12, p - (fl + fr + ur - ul) / (dfl + dfr))
        if abs(p_new - p) < 1e-15 * p:
            p = p_new
            break
        p = p_new
    fl, _ = f_k(p, rl, pl, al)
    fr, _ = f_k(p, rr, pr, ar)
    u = 0.5 * (ul + ur) + 0.5 * (fr - fl)
    g1, g2 = (gamma - 1.0) / (gamma + 1.0), (gamma - 1.0) / (2.0 * gamma)

    rho = np.empty_like(x)
    for n, xn in enumerate(x):
        s = (xn - x0) / t
        if s <= u:  # left of the contact
            if p > pl:  # left shock
                sl = ul - al * math.sqrt((gamma + 1.0) / (2.0 * gamma) * p / pl + g2)
                rho[n] = rl if s <= sl else rl * (p / pl + g1) / (g1 * p / pl + 1.0)
            else:  # left rarefaction
                shl, stl = ul - al, u - al * (p / pl) ** g2
                if s <= shl:
                    rho[n] = rl
                elif s >= stl:
                    rho[n] = rl * (p / pl) ** (1.0 / gamma)
                else:
                    rho[n] = rl * (2.0 / (gamma + 1.0) + g1 / al * (ul - s)) ** (2.0 / (gamma - 1.0))
        else:  # right of the contact
            if p > pr:  # right shock
                sr = ur + ar * math.sqrt((gamma + 1.0) / (2.0 * gamma) * p / pr + g2)
                rho[n] = rr if s >= sr else rr * (p / pr + g1) / (g1 * p / pr + 1.0)
            else:  # right rarefaction
                shr, str_ = ur + ar, u + ar * (p / pr) ** g2
                if s >= shr:
                    rho[n] = rr
                elif s <= str_:
                    rho[n] = rr * (p / pr) ** (1.0 / gamma)
                else:
                    rho[n] = rr * (2.0 / (gamma + 1.0) - g1 / ar * (s - ur)) ** (2.0 / (gamma - 1.0))
    return rho


def load_profile(work: Path, basename: str, ngc: int) -> tuple[int, int, np.ndarray, np.ndarray]:
    """Return the last saved iteration, the active axis and its 1-D profiles (centres, conservative variables).

    The active axis is the longest block extent. Every transverse copy of a 1-D column must be bitwise identical to
    the first one, and two blocks holding the same column must agree bitwise, otherwise the oracle stops.
    """
    files = sorted(work.glob(f"{basename}-*-proc*.h5"))
    if not files:
        sys.exit(f"sod_oracle: no {basename}-*-proc*.h5 in {work}")
    last = max(int(f.name.split("-")[-2]) for f in files)
    cells: dict[float, np.ndarray] = {}
    axis = -1
    for path in (f for f in files if int(f.name.split("-")[-2]) == last):
        with h5py.File(path, "r") as h5:
            blocks = sorted({k.rsplit("-", 1)[0] for k in h5})
            for blk in blocks:
                # geometry is stored (z, y, x), as XDMF ORIGIN_DXDYDZ expects: reverse it to (x, y, z)
                origin, dxyz = h5[f"{blk}-origin"][()][::-1], h5[f"{blk}-dxdydz"][()][::-1]
                # datasets are stored (z, y, x) with ghosts: transpose to (x, y, z), strip the ghosts
                q = np.stack([h5[f"{blk}-{v}"][()].transpose(2, 1, 0) for v in ("r", "ru", "rv", "rw", "rE")])
                if ngc:
                    q = q[:, ngc:-ngc, ngc:-ngc, ngc:-ngc]
                if axis < 0:
                    axis = int(np.argmax(q.shape[1:]))
                n = q.shape[1 + axis]
                qa = np.moveaxis(q, 1 + axis, 1).reshape(5, n, -1)
                if not np.all(qa == qa[:, :, :1]):
                    sys.exit(f"sod_oracle: transverse copies differ in {path.name}:{blk}")
                # checkpoints are saved with ghosts, so the origin is the corner of the first ghost cell
                centres = origin[axis] + (np.arange(n) + ngc + 0.5) * dxyz[axis]
                for c, col in zip(centres, qa[:, :, 0].T, strict=True):
                    if c in cells and not np.array_equal(cells[c], col):
                        sys.exit(f"sod_oracle: two blocks disagree on the column at {c}")
                    cells[c] = col
    xs = np.array(sorted(cells))
    return last, axis, xs, np.array([cells[c] for c in xs]).T


def main() -> int:
    """Run the V1 checks, print a report, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("ini", type=Path, help="one of the sod-{x,y,z}.ini inputs (physics, time, IC)")
    parser.add_argument("work", type=Path, nargs="+", help="work directories of the x (and y, z) runs")
    parser.add_argument("--l1-max", type=float, default=None, help="L1(rho) bound (per run)")
    parser.add_argument("--tol", type=float, default=0.0, help="comparison tolerance (default 0: bitwise)")
    parser.add_argument("--mirror", type=float, default=None, help="mirror-symmetry tolerance (wall case)")
    parser.add_argument("--closed", type=float, default=None, help="mass/energy relative drift bound (wall case)")
    args = parser.parse_args()

    ini = read_ini(args.ini)
    gamma = float(ini["physics"]["cp"]) / float(ini["physics"]["cv"])
    t = float(ini["time"]["time_max"])
    ngc = int(ini["grid"]["ngc"])
    r1, r2 = ini["initial_conditions_region_1"], ini["initial_conditions_region_2"]
    left = (float(r1["r"]), 0.0, float(r1["p"]))
    right = (float(r2["r"]), 0.0, float(r2["p"]))

    status, profiles = 0, []
    for work in args.work:
        basename = next(work.glob("*-residuals.dat")).name.removesuffix("-residuals.dat")
        it, axis, xs, q = load_profile(work, basename, ngc)
        if args.mirror is not None:
            if np.max(np.abs(xs + xs[::-1] - 1.0)) > 1e-12:
                sys.exit(f"sod_oracle: the cell centres of {work.name} are not symmetric about 0.5")
            mirrored = q[:, ::-1].copy()
            mirrored[1 + axis] = -mirrored[1 + axis]
            asym = float(np.max(np.abs(q - mirrored)))
            ok = asym <= args.mirror
            status |= 0 if ok else 1
            print(f"{work.name}/{basename}: axis {'xyz'[axis]}, it {it}, cells {xs.size}, max |q(x) - S q(1-x)| = "
                  f"{asym:.3e}  {'PASS' if ok else 'FAIL'} (tol {args.mirror:.1e})")
            if args.closed is not None:
                rows = [line.split() for line in open(work / f"{basename}-conservation_history.dat")]
                hist = np.array([[float(v) for v in r] for r in rows if r and r[0][0] in "+-0123456789"])
                drift = [float(np.max(np.abs(hist[:, c] - hist[0, c])) / abs(hist[0, c])) for c in (2, 6)]
                ok = max(drift) <= args.closed
                status |= 0 if ok else 1
                print(f"{work.name}/{basename}: closed box, relative drift mass {drift[0]:.2e} energy {drift[1]:.2e}  "
                      f"{'PASS' if ok else 'FAIL'} (tol {args.closed:.1e})")
        else:
            dx = xs[1] - xs[0]
            rho_exact = exact_riemann_density(xs, t, 0.5, gamma, left, right)
            l1 = float(np.sum(np.abs(q[0] - rho_exact)) * dx)
            verdict = ""
            if args.l1_max is not None:
                ok = l1 <= args.l1_max
                status |= 0 if ok else 1
                verdict = "  PASS" if ok else f"  FAIL (bound {args.l1_max:.3e})"
            print(f"{work.name}/{basename}: axis {'xyz'[axis]}, it {it}, cells {xs.size}, L1(rho) = {l1:.6e}{verdict}")
        # permute the momentum so that the active-axis component always comes first
        mom = np.roll(q[1:4], -axis, axis=0)
        profiles.append((f"{work.name}/{basename}", np.vstack([q[0:1], mom, q[4:5]])))
    for name, prof in profiles[1:]:
        diff = float(np.max(np.abs(prof - profiles[0][1])))
        ok = diff <= args.tol
        status |= 0 if ok else 1
        kind = "bitwise" if args.tol == 0.0 else f"tol {args.tol:.1e}"
        print(f"{name} vs {profiles[0][0]}: max |difference| = {diff:.3e}  {'PASS' if ok else 'FAIL'} ({kind})")
    return status


if __name__ == "__main__":
    sys.exit(main())
