#!/usr/bin/env python3
"""FLUME verification V2 oracle: isentropic vortex convergence (issue #35, section 11).

Why this oracle exists: V1 (Sod) proves the scheme is correct across discontinuities, where WENO drops to first
order; it says nothing about the design order on smooth flows. The isentropic vortex (Shu 1998, ICASE 97-65, section
5.1) is an exact steady solution convected by the free stream, so the error of each run is measured pointwise against
it, and the observed order between successive resolutions, p = log2(E_N / E_2N), must reach the design order.

The exact solution uses the same formulas as the FLUME initial condition (adam_flume_ic_object), with the vortex
centre moved by the free stream and wrapped periodically to the nearest image.

Usage:
    vortex_oracle.py <work-dir> [<work-dir> ...] [--order-min P]

Each work directory holds one run (its INI and final checkpoint); directories are sorted by resolution.
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


def exact_density(x: np.ndarray, y: np.ndarray, ini: configparser.ConfigParser) -> np.ndarray:
    """Return the exact vortex density at the final time, at points (x, y) of the periodic box."""
    gamma = float(ini["physics"]["cp"]) / float(ini["physics"]["cv"])
    t = float(ini["time"]["time_max"])
    ic, fs = ini["initial_conditions"], ini["initial_conditions_region_1"]
    x0, y0, radius, strength = (float(ic[k]) for k in ("x0", "y0", "radius", "strength"))
    r0, u0, v0, p0 = (float(fs[k]) for k in ("r", "u", "v", "p"))
    lx = float(ini["grid"]["emax_x"]) - float(ini["grid"]["emin_x"])
    ly = float(ini["grid"]["emax_y"]) - float(ini["grid"]["emin_y"])
    dx = x - (x0 + u0 * t)
    dy = y - (y0 + v0 * t)
    dx -= lx * np.round(dx / lx)
    dy -= ly * np.round(dy / ly)
    rr = (dx * dx + dy * dy) / (radius * radius)
    t0 = p0 / r0
    temp = t0 - (gamma - 1.0) / gamma * strength**2 / (8.0 * math.pi**2) * np.exp(1.0 - rr)
    return r0 * (temp / t0) ** (1.0 / (gamma - 1.0))


def run_error(work: Path) -> tuple[int, float, float]:
    """Return cells per side, L1 and Linf density errors of the last checkpoint of one run."""
    ini = read_ini(next(work.glob("vortex-*.ini")))
    ngc = int(ini["grid"]["ngc"])
    basename = ini["IO"]["output_basename"]
    files = sorted(work.glob(f"{basename}-*-proc*.h5"))
    if not files:
        sys.exit(f"vortex_oracle: no {basename}-*-proc*.h5 in {work}")
    last = max(int(f.name.split("-")[-2]) for f in files)
    err1, errinf, area = 0.0, 0.0, 0.0
    for path in (f for f in files if int(f.name.split("-")[-2]) == last):
        with h5py.File(path, "r") as h5:
            for blk in sorted({k.rsplit("-", 1)[0] for k in h5}):
                # geometry is stored (z, y, x) and the origin is the corner of the first ghost cell
                origin, dxyz = h5[f"{blk}-origin"][()][::-1], h5[f"{blk}-dxdydz"][()][::-1]
                rho = h5[f"{blk}-r"][()].transpose(2, 1, 0)[ngc:-ngc, ngc:-ngc, ngc:-ngc]
                if not np.all(rho == rho[:, :, :1]):
                    sys.exit(f"vortex_oracle: copies along the null z direction differ in {path.name}:{blk}")
                rho = rho[:, :, 0]
                xc = origin[0] + (np.arange(rho.shape[0]) + ngc + 0.5) * dxyz[0]
                yc = origin[1] + (np.arange(rho.shape[1]) + ngc + 0.5) * dxyz[1]
                xx, yy = np.meshgrid(xc, yc, indexing="ij")
                diff = np.abs(rho - exact_density(xx, yy, ini))
                err1 += float(diff.sum()) * dxyz[0] * dxyz[1]
                errinf = max(errinf, float(diff.max()))
                area += diff.size * dxyz[0] * dxyz[1]
    # the blocks stacked along the null z direction repeat the same 2-D field: they weigh err1 and area alike
    lx = float(ini["grid"]["emax_x"]) - float(ini["grid"]["emin_x"])
    return round(lx / dxyz[0]), err1 / area, errinf


def main() -> int:
    """Measure the errors of every run and the observed orders, print a report, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work", type=Path, nargs="+", help="work directories, one per resolution")
    parser.add_argument("--order-min", type=float, default=None, help="minimum L1 order of the finest pair")
    args = parser.parse_args()

    results = sorted(run_error(w) for w in args.work)
    status = 0
    prev = None
    for n, e1, einf in results:
        line = f"N {n:4d}: L1(rho) = {e1:.6e}  Linf(rho) = {einf:.6e}"
        if prev is not None:
            p1 = math.log(prev[1] / e1) / math.log(n / prev[0])
            pinf = math.log(prev[2] / einf) / math.log(n / prev[0])
            line += f"  order L1 {p1:+.2f}  Linf {pinf:+.2f}"
        print(line)
        prev = (n, e1, einf)
    if args.order_min is not None and len(results) > 1:
        (n0, e0, _), (n1, e1, _) = results[-2], results[-1]
        p = math.log(e0 / e1) / math.log(n1 / n0)
        ok = p >= args.order_min
        status = 0 if ok else 1
        print(f"finest-pair L1 order {p:+.2f} {'>=' if ok else '<'} {args.order_min}: {'PASS' if ok else 'FAIL'}")
    return status


if __name__ == "__main__":
    sys.exit(main())
