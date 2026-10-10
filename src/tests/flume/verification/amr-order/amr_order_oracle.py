#!/usr/bin/env python3
"""Oracle of the AO ladders (#68 A1): the error of a smooth vortex across static 2:1 seams, per region, and its orders.

Why: a composite-grid error norm mixes the seam band, where the transfer operators and the reflux act, with the
interior, where the uniform scheme acts; the #68 audit predicts second order in the band and design order elsewhere,
and only separate norms can show it. Each run's last checkpoint is compared cell by cell with the exact convected
vortex (euler: the V2 density; mhd: the MV-7 conservative state, its 8 variables combined as MV-7 does,
sqrt(sum_v mean_v^2)). The cells are split by the periodic distance of their centre to the boundary of the refined
box: the seam band holds the cells within `ngc` of their own spacing from it (where a stencil reads seam ghosts or
the reflux acts), the interior holds the rest; a uniform run (no box) is all interior. For each region: L1 (the
volume-weighted mean of |error|) and Linf, and the observed orders of successive runs (N doubling).

Usage:
    amr_order_oracle.py <case> <work> [<work> ...] [--box xmin ymin xmax ymax] [--assert REGION:NORM:MIN ...]
    --assert interior:l1:4.5 fails unless the finest-pair order of that region and norm is at least MIN
    --max band:l1:1.6e-5 fails unless that error of the finest run is at most MAX
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
sys.path.insert(0, str(VERIF / "vortex"))
sys.path.insert(0, str(VERIF / "mhd" / "vortex"))
sys.path.insert(0, str(VERIF / "mhd" / "linear-wave"))
from mhd_vortex_oracle import exact as mhd_exact  # noqa: E402
from vortex_oracle import exact_density  # noqa: E402

VARIABLES = {"euler": ("r",), "mhd": ("r", "ru", "rv", "rw", "rE", "bx", "by", "bz")}
REGIONS = ("all", "band", "interior")


def read_ini(work: Path) -> configparser.ConfigParser:
    """Return the input of a run (its only .ini file)."""
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.read(next(work.glob("*.ini")))
    return ini


def last_cells(work: Path, names: tuple[str, ...]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (centres [n, 2], spacings [n, 2], values [n, v]) of the interior cells of the last checkpoint; the copies
    along the null z direction are one cell (they must be equal)."""
    ini = read_ini(work)
    ngc = int(ini["grid"]["ngc"])
    files = sorted(p for p in work.glob("*-proc*.h5") if "restart" not in p.name)
    if not files:
        sys.exit(f"amr_order_oracle: no checkpoint in {work}")
    last = max(int(p.name.split("-")[-2]) for p in files)
    xy, hh, qq = [], [], []
    for path in (p for p in files if int(p.name.split("-")[-2]) == last):
        with h5py.File(path, "r") as h5:
            for blk in sorted({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")}):
                origin, dxyz = h5[f"{blk}-origin"][()][::-1], h5[f"{blk}-dxdydz"][()][::-1]
                q = np.stack([h5[f"{blk}-{v}"][()].transpose(2, 1, 0)[ngc:-ngc, ngc:-ngc, ngc:-ngc] for v in names])
                if not np.all(q == q[:, :, :, :1]):
                    sys.exit(f"amr_order_oracle: copies along the null z direction differ in {path.name}:{blk}")
                if origin[2] + ngc * dxyz[2] > float(ini["grid"]["emin_z"]) + 0.5 * dxyz[2]:
                    continue  # a block stacked along the null z direction repeats the bottom one
                q = q[:, :, :, 0]
                xc = origin[0] + (np.arange(q.shape[1]) + ngc + 0.5) * dxyz[0]
                yc = origin[1] + (np.arange(q.shape[2]) + ngc + 0.5) * dxyz[1]
                xx, yy = np.meshgrid(xc, yc, indexing="ij")
                xy.append(np.stack([xx.ravel(), yy.ravel()], axis=1))
                hh.append(np.tile(dxyz[:2], (xx.size, 1)))
                qq.append(q.reshape(len(names), -1).T)
    return np.concatenate(xy), np.concatenate(hh), np.concatenate(qq)


def box_distance(xy: np.ndarray, box: list[float], lo: np.ndarray, size: np.ndarray) -> np.ndarray:
    """Return the distance of each point to the boundary of the box, over the periodic images of the point."""
    best = np.full(len(xy), np.inf)
    for sx in (-1, 0, 1):
        for sy in (-1, 0, 1):
            p = xy + np.array([sx, sy]) * size
            inside = (p[:, 0] >= box[0]) & (p[:, 0] <= box[2]) & (p[:, 1] >= box[1]) & (p[:, 1] <= box[3])
            d_in = np.minimum.reduce([p[:, 0] - box[0], box[2] - p[:, 0], p[:, 1] - box[1], box[3] - p[:, 1]])
            ox = np.maximum.reduce([box[0] - p[:, 0], np.zeros(len(p)), p[:, 0] - box[2]])
            oy = np.maximum.reduce([box[1] - p[:, 1], np.zeros(len(p)), p[:, 1] - box[3]])
            best = np.minimum(best, np.where(inside, d_in, np.hypot(ox, oy)))
    return best


def region_errors(case: str, work: Path, box: list[float] | None) -> tuple[int, dict[str, tuple[float, float, int]]]:
    """Return the base cells per side and, per region, (L1, Linf, cells) of one run."""
    ini = read_ini(work)
    ngc = int(ini["grid"]["ngc"])
    names = VARIABLES[case]
    xy, h, q = last_cells(work, names)
    t = float(ini["time"]["time_max"])
    if case == "euler":
        err = np.abs(q[:, 0] - exact_density(xy[:, 0], xy[:, 1], ini))[:, None]
    else:
        err = np.abs(q - mhd_exact(ini, xy, t).T)
    lo = np.array([float(ini["grid"]["emin_x"]), float(ini["grid"]["emin_y"])])
    size = np.array([float(ini["grid"]["emax_x"]), float(ini["grid"]["emax_y"])]) - lo
    vol = h[:, 0] * h[:, 1]
    band = np.zeros(len(xy), dtype=bool)
    if box is not None:
        band = box_distance(xy, box, lo, size) < ngc * h.max(axis=1)
    out = {}
    for region, mask in (("all", np.ones(len(xy), dtype=bool)), ("band", band), ("interior", ~band)):
        if not mask.any():
            out[region] = (math.nan, math.nan, 0)
            continue
        means = (err[mask] * vol[mask, None]).sum(axis=0) / vol[mask].sum()
        out[region] = (float(math.sqrt(np.sum(means**2))), float(err[mask].max()), int(mask.sum()))
    return round(size[0] / h.max()), out


def main() -> int:
    """Measure every run, print the per-region errors and orders, check the assertions, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("case", choices=tuple(VARIABLES))
    parser.add_argument("work", type=Path, nargs="+", help="runs, coarse first, N doubling")
    parser.add_argument("--box", type=float, nargs=4, default=None, help="refined box xmin ymin xmax ymax")
    parser.add_argument("--assert", dest="asserts", nargs="+", default=[], help="REGION:NORM:MIN, norm l1 or linf")
    parser.add_argument("--max", dest="maxima", nargs="+", default=[], help="REGION:NORM:MAX on the finest run")
    args = parser.parse_args()
    rows = [region_errors(args.case, w, args.box) for w in args.work]
    for w, (n, reg) in zip(args.work, rows, strict=True):
        print(f"{w.name}: N {n}  " + "  ".join(
            f"{r} L1 {reg[r][0]:.4e} Linf {reg[r][1]:.4e} ({reg[r][2]} cells)" for r in REGIONS if reg[r][2]))
    orders: dict[tuple[str, str], list[float]] = {}
    for r in REGIONS:
        for k, norm in enumerate(("l1", "linf")):
            vals = [reg[r][k] for _, reg in rows]
            if any(not v > 0.0 for v in vals):
                continue
            orders[(r, norm)] = [math.log2(a / b) for a, b in zip(vals[:-1], vals[1:], strict=True)]
            print(f"   order {r:>8} {norm:>4}: " + " ".join(f"{p:+.2f}" for p in orders[(r, norm)]))
    status = 0
    for spec in args.asserts:
        region, norm, low = spec.split(":")
        p = orders.get((region, norm), [math.nan])[-1]
        ok = p >= float(low)
        status |= 0 if ok else 1
        print(f"   {region} {norm} finest-pair order {p:+.2f} >= {float(low):.2f}: {'PASS' if ok else 'FAIL'}")
    for spec in args.maxima:
        region, norm, high = spec.split(":")
        e = rows[-1][1][region][("l1", "linf").index(norm)]
        ok = e <= float(high)
        status |= 0 if ok else 1
        print(f"   {region} {norm} of the finest run {e:.4e} <= {float(high):.3e}: {'PASS' if ok else 'FAIL'}")
    return status


if __name__ == "__main__":
    sys.exit(main())
