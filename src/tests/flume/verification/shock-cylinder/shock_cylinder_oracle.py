#!/usr/bin/env python3
"""FLUME verification oracle of the shock-cylinder case: immersed boundary with init-time AMR (issue #35, P6).

Why this oracle exists: the shock-cylinder flow has no closed-form solution, but its setup (cylinder centred at
y = 0.5, walls at y = 0 and 1) is mirror-symmetric about y = 0.5, and so is the exact solution. Mirror symmetry is
blind to defects that are themselves symmetric (a wrong wall law gives a wrong but symmetric flow); it catches the
ones that are not: indexing and ordering errors in the distance function, the eikonal extrapolation, the cut spacing
or the solid AMR marker, stencils reading the wrong ghost, data races. The checks, on the last checkpoint:

* the refinement: the blocks crossed by the cylinder surface are one level finer than the others;
* mirror symmetry: every block equals its mirror block (y -> 1 - y, v -> -v) within the tolerance;
* sanity: density and pressure are positive and finite in the fluid (outside the cylinder);
* --compare: two runs (e.g. CPU and FNL) agree cell by cell within a relative tolerance.

Usage:
    shock_cylinder_oracle.py <work-dir> [--mirror-tol T] [--compare <work-dir-b> --tol T]
"""

from __future__ import annotations

import argparse
import configparser
import math
import sys
from pathlib import Path

import h5py
import numpy as np

VARIABLES = ("r", "ru", "rv", "rw", "rE")


def read_ini(path: Path) -> configparser.ConfigParser:
    """Read a FLUME INI file (`;` comments, no interpolation)."""
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.read(path)
    return ini


def last_blocks(work: Path) -> dict[tuple[float, ...], tuple[np.ndarray, np.ndarray]]:
    """Return {interior origin (x, y, z): (spacing, q[v, x, y, z] interior)} of the last checkpoint of one run."""
    ini = read_ini(next(work.glob("*.ini")))
    ngc = int(ini["grid"]["ngc"])
    files = sorted(p for p in work.glob("*-proc*.h5") if "restart" not in p.name)
    last = max(int(p.name.split("-")[-2]) for p in files)
    out = {}
    for path in (p for p in files if int(p.name.split("-")[-2]) == last):
        with h5py.File(path, "r") as h5:
            for blk in {k.rsplit("-", 1)[0] for k in h5}:
                # geometry stored (z, y, x); the origin is the corner of the first ghost cell
                dx = h5[f"{blk}-dxdydz"][()][::-1]
                lo = h5[f"{blk}-origin"][()][::-1] + ngc * dx
                q = np.stack([h5[f"{blk}-{v}"][()].transpose(2, 1, 0) for v in VARIABLES])
                out[tuple(float(x) for x in np.round(lo, 12))] = (dx, q[:, ngc:-ngc, ngc:-ngc, ngc:-ngc])
    return out


def main() -> int:
    """Run the checks, print a report, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work", type=Path, help="work directory of the run")
    parser.add_argument("--mirror-tol", type=float, default=1.0e-10, help="mirror symmetry tolerance (relative)")
    parser.add_argument("--compare", type=Path, default=None, help="second run to compare with")
    parser.add_argument("--tol", type=float, default=1.0e-10, help="comparison tolerance (relative)")
    args = parser.parse_args()

    ini = read_ini(next(args.work.glob("*.ini")))
    sol = ini["solid_1"]
    cx, cy, radius = float(sol["circle_center_x"]), float(sol["circle_center_y"]), float(sol["circle_radius"])
    gamma = float(ini["physics"]["cp"]) / float(ini["physics"]["cv"])
    blocks = last_blocks(args.work)
    status = 0

    # refinement: the blocks crossed by the surface are finer
    coarse = max(float(dx[0]) for dx, _ in blocks.values())
    crossed, fine, wrong = 0, 0, 0
    for lo, (dx, q) in blocks.items():
        xs = lo[0] + (np.arange(q.shape[1]) + 0.5) * dx[0]
        ys = lo[1] + (np.arange(q.shape[2]) + 0.5) * dx[1]
        xx, yy = np.meshgrid(xs, ys, indexing="ij")
        dist = np.hypot(xx - cx, yy - cy) - radius
        is_fine = bool(dx[0] < coarse)
        is_crossed = bool(dist.min() < 0.0 < dist.max())
        crossed += is_crossed
        fine += is_fine
        wrong += is_crossed and not is_fine
    ok = wrong == 0 and fine > 0
    status |= 0 if ok else 1
    print(f"refinement: {len(blocks)} blocks, {fine} fine, {crossed} crossed by the surface, {wrong} crossed but coarse"
          f"  {'PASS' if ok else 'FAIL'}")

    # mirror symmetry about y = 0.5
    scale = np.max([np.max(np.abs(q), axis=(1, 2, 3)) for _, q in blocks.values()], axis=0)
    scale = np.where(scale > 0.0, scale, 1.0)  # an identically zero variable (w in 2-D) is compared absolutely
    asym = 0.0
    for lo, (dx, q) in blocks.items():
        ly = q.shape[2] * dx[1]
        key = (lo[0], round(1.0 - lo[1] - ly, 12), lo[2])
        if key not in blocks:
            sys.exit(f"shock_cylinder_oracle: block {lo} has no mirror block {key}")
        m = blocks[key][1][:, :, ::-1, :].copy()
        m[2] = -m[2]
        asym = max(asym, float(np.max(np.max(np.abs(q - m), axis=(1, 2, 3)) / scale)))
    ok = asym <= args.mirror_tol
    status |= 0 if ok else 1
    print(f"mirror symmetry about y = 0.5: max relative asymmetry {asym:.3e}  {'PASS' if ok else 'FAIL'}"
          f" (tol {args.mirror_tol:.1e})")

    # sanity in the fluid
    bad = 0
    rmin, pmin = math.inf, math.inf
    for lo, (dx, q) in blocks.items():
        xs = lo[0] + (np.arange(q.shape[1]) + 0.5) * dx[0]
        ys = lo[1] + (np.arange(q.shape[2]) + 0.5) * dx[1]
        xx, yy = np.meshgrid(xs, ys, indexing="ij")
        fluid = np.broadcast_to((np.hypot(xx - cx, yy - cy) > radius)[:, :, None], q.shape[1:])
        r = q[0][fluid]
        p = (gamma - 1.0) * (q[4][fluid] - 0.5 * (q[1][fluid]**2 + q[2][fluid]**2 + q[3][fluid]**2) / r)
        bad += int(np.sum(~np.isfinite(r)) + np.sum(~np.isfinite(p)))
        rmin, pmin = min(rmin, float(r.min())), min(pmin, float(p.min()))
    ok = bad == 0 and rmin > 0.0 and pmin > 0.0
    status |= 0 if ok else 1
    print(f"fluid sanity: min rho {rmin:.4f}, min p {pmin:.4f}, non-finite {bad}  {'PASS' if ok else 'FAIL'}")

    if args.compare is not None:
        other = last_blocks(args.compare)
        if other.keys() != blocks.keys():
            sys.exit("shock_cylinder_oracle: the two runs hold different block sets")
        delta = np.max([np.max(np.abs(blocks[k][1] - other[k][1]), axis=(1, 2, 3)) for k in blocks], axis=0)
        diff = float(np.max(delta / scale))
        ok = diff <= args.tol
        status |= 0 if ok else 1
        print(f"{args.work.name} vs {args.compare.name}: max relative difference {diff:.3e}  {'PASS' if ok else 'FAIL'}"
              f" (tol {args.tol:.1e})")
    return status


if __name__ == "__main__":
    sys.exit(main())
