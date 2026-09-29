#!/usr/bin/env python3
"""Oracle of the rotated shock tube (issue #41, M2-P6, MV-8; Toth 2000, J. Comput. Phys. 161, section 6.3.2).

Why: the exact solution of the rotated Ryu-Jones 1a problem keeps the normal field uniform, B_par = 5 / sqrt(4 pi), and
is the 1-D solution along the normal. From the last checkpoint of a 2-D run, in the frame of the normal n:
  1. dB_par = sum |B_par - 5/sqrt(4 pi)| / sum 5/sqrt(4 pi) (Toth eq. 45 against the analytic value; his base scheme,
     second order, minmod: 0.0037 at 63.4 degrees, N = 256), at most --dbpar-max;
  2. the relative L1 errors d(rho, v_par, v_perp, p, B_par, B_perp) against the 1-D reference run (the same problem
     along x, 4 times finer along the normal, linearly interpolated in s = n_x x + n_y y, periodic), and their mean,
     at most --dmean-max (Toth, 63.4 degrees, N = 256: 0.0238 for his base scheme);
  3. over a ladder of runs (coarse first), dB_par and the mean error decrease (discontinuities: first order).

Usage:
    rotated_shock_tube_oracle.py <work> [<work> ...] --reference <work-1d> [--dbpar-max E [E ...]]
                                 [--dmean-max E [E ...]] [--ngc N]
"""

from __future__ import annotations

import argparse
import configparser
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "linear-wave"))
from linear_wave_oracle import NAMES, profile  # noqa: E402

B0 = 5.0 / math.sqrt(4.0 * math.pi)
PRIMS = ("rho", "v_par", "v_perp", "p", "B_par", "B_perp")


def setup(work: Path) -> tuple[float, float, float, float]:
    """Return (normal_x, normal_y, period, gamma) of the input of a run (its only .ini file)."""
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.read(next(work.glob("*.ini")))
    ic = ini["initial_conditions"]
    gamma = float(ini["physics"]["cp"]) / float(ini["physics"]["cv"])
    return float(ic["normal_x"]), float(ic["normal_y"]), float(ic["period"]), gamma


def primitives(work: Path, ngc: int) -> tuple[np.ndarray, np.ndarray]:
    """Return (s modulo the period, primitives in the frame of the normal [cells, 6]) of the last checkpoint."""
    nx, ny, period, gamma = setup(work)
    cells = profile(work, "last", ngc)
    centres = np.array(list(cells))
    q = np.array(list(cells.values()))
    norm = math.hypot(nx, ny)
    n = (nx / norm, ny / norm)
    r = q[:, NAMES.index("r")]
    mx, my, mz = (q[:, NAMES.index(v)] for v in ("ru", "rv", "rw"))
    bx, by, bz = (q[:, NAMES.index(v)] for v in ("bx", "by", "bz"))
    p = (gamma - 1.0) * (q[:, NAMES.index("rE")] - 0.5 * (mx**2 + my**2 + mz**2) / r - 0.5 * (bx**2 + by**2 + bz**2))
    prim = np.stack([r, (mx * n[0] + my * n[1]) / r, (-mx * n[1] + my * n[0]) / r, p, bx * n[0] + by * n[1],
                     -bx * n[1] + by * n[0]], axis=1)
    s = np.mod(nx * centres[:, 0] + ny * centres[:, 1], period)
    return s, prim


def bound(values: list[float] | None, n: int) -> float | None:
    """Return the n-th bound of a per-run list (one value applies to every run)."""
    if values is None:
        return None
    return values[0] if len(values) == 1 else values[n]


def measures(work: Path, s_ref: np.ndarray, prim_ref: np.ndarray, ngc: int) -> tuple[float, float, float, np.ndarray]:
    """Return (dB_par, max |B_par - B0| / B0, mean error, errors per primitive) of a 2-D run against the reference."""
    s, prim = primitives(work, ngc)
    period = setup(work)[2]
    ref = np.stack([np.interp(s, s_ref, prim_ref[:, v], period=period) for v in range(len(PRIMS))], axis=1)
    dbpar = float(np.sum(np.abs(prim[:, 4] - B0)) / (B0 * len(s)))
    dmax = float(np.max(np.abs(prim[:, 4] - B0)) / B0)
    errors = np.sum(np.abs(prim - ref), axis=0) / np.sum(np.abs(ref), axis=0)
    return dbpar, dmax, float(np.mean(errors)), errors


def main() -> int:
    """Run the MV-8 checks, print a report, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work", type=Path, nargs="+", help="the 2-D runs, coarse first")
    parser.add_argument("--reference", type=Path, required=True, help="the 1-D reference run")
    parser.add_argument("--dbpar-max", type=float, nargs="+", default=None, help="bound of dB_par (one or per run)")
    parser.add_argument("--dmean-max", type=float, nargs="+", default=None, help="bound of the mean error (idem)")
    parser.add_argument("--ngc", type=int, default=3)
    args = parser.parse_args()
    for name in ("dbpar_max", "dmean_max"):
        values = getattr(args, name)
        if values is not None and len(values) not in (1, len(args.work)):
            sys.exit(f"rotated_shock_tube_oracle: --{name.replace('_', '-')} takes one bound or one per run")
    s_ref, prim_ref = primitives(args.reference, args.ngc)
    order = np.argsort(s_ref)
    s_ref, prim_ref = s_ref[order], prim_ref[order]
    ok = True
    previous = None
    for n, work in enumerate(args.work):
        dbpar, dmax, dmean, errors = measures(work, s_ref, prim_ref, args.ngc)
        print(f"{work.name} (reference {args.reference.name}):")
        print(f"   dB_par (analytic) {dbpar:.4e}, max |B_par - B0| / B0 {dmax:.4e}")
        print("   vs reference: " + ", ".join(f"d{name} {e:.4e}" for name, e in zip(PRIMS, errors, strict=True)) +
              f"; mean {dmean:.4e}")
        for name, value, limit in (("dB_par", dbpar, bound(args.dbpar_max, n)),
                                   ("mean error", dmean, bound(args.dmean_max, n))):
            if limit is not None:
                good = value <= limit
                ok &= good
                print(f"   {name} {value:.4e} {'<=' if good else '>'} {limit:.3e}: {'PASS' if good else 'FAIL'}")
        if previous is not None:
            for name, coarse, fine in (("dB_par", previous[0], dbpar), ("mean error", previous[1], dmean)):
                good = fine < coarse
                ok &= good
                print(f"   {name} rate {math.log2(coarse / fine):+.2f} (decreasing): {'PASS' if good else 'FAIL'}")
        previous = (dbpar, dmean)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
