"""RIPPLE oracle (issue #59): read the runs' XH5F output and check them against the exact solution.

The exact solution of u_t + a . grad(u) = 0 on the periodic box is the initial profile translated by a t,
u0 = prod_d sin(2 pi (x_d - emin_d) / L_d) over the non-null axes (adam_ripple_common_object). The checkpoints hold
point values at the cell centres, one HDF5 file per rank, every array in (z, y, x) order.

Subcommands:
  convergence WORK...  the L1 and Linf errors of the last checkpoint of each run (coarse to fine) and the observed
                       orders of each pair; --min-order bounds the L1 order of the finest pair, --max-l1 the finest L1
  conservation WORK    the drift of the total of u dV in the history file, relative to max|u0| times the domain
                       volume, bounded by --tol
  compare WORK WORK    the largest difference of u between two runs' last checkpoints, block by block (keyed by the
                       block origin), bounded by --tol
  refined WORK         the base blocks a static AMR run refined (those holding blocks finer than the coarsest) must be
                       those where u0 at the base resolution rises above --threshold (issue #110: the crest criterion)
"""

from __future__ import annotations

import argparse
import configparser
import math
import sys
from pathlib import Path

import h5py
import numpy as np


def read_input(work: Path) -> configparser.ConfigParser:
    """The run's input.ini, inline `;` comments stripped."""
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",))
    ini.read(work / "input.ini")
    return ini


def last_iteration(work: Path) -> int:
    """The iteration of the last checkpoint (the largest `-<it>-proc` file name)."""
    return max(int(p.name.rsplit("-proc", 1)[0].rsplit("-", 1)[1]) for p in work.glob("*-proc*.h5"))


def blocks(work: Path) -> dict[tuple[float, ...], tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Every block of the last checkpoint: origin (x, y, z) -> (spacing (x, y, z), u in (x, y, z) order, origin)."""
    it = last_iteration(work)
    out = {}
    for path in sorted(work.glob(f"*-{it:09d}-proc*.h5")):
        with h5py.File(path, "r") as h5:
            for blk in sorted({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")}):
                origin = h5[f"{blk}-origin"][()][::-1]
                dxyz = h5[f"{blk}-dxdydz"][()][::-1]
                u = h5[f"{blk}-u"][()].transpose(2, 1, 0)
                out[tuple(np.round(origin, 12))] = (dxyz, u, origin)
    return out


def final_time(work: Path) -> float:
    """The time of the last history row."""
    rows = [ln.split() for ln in (work / f"{basename(work)}-history.dat").read_text().splitlines() if ln[:1] != "#"]
    return float(rows[-1][1])


def basename(work: Path) -> str:
    return read_input(work)["IO"]["output_basename"].strip()


def domain(work: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Velocity, domain corners and null axes of a run."""
    ini = read_input(work)
    a = np.array([float(ini["ripple"][f"velocity_{c}"]) for c in "xyz"])
    emin = np.array([float(ini["grid"][f"emin_{c}"]) for c in "xyz"])
    emax = np.array([float(ini["grid"][f"emax_{c}"]) for c in "xyz"])
    null = np.array([ini["grid"][f"null_{c}"].strip().lower() in (".true.", "true") for c in "xyz"])
    return np.where(null, 0.0, a), emin, emax, null


def errors(work: Path) -> tuple[float, float, int]:
    """L1 (volume-weighted mean) and Linf errors of the last checkpoint against the exact solution, cells number."""
    a, emin, emax, null = domain(work)
    t = final_time(work)
    e1, vol, einf, n = 0.0, 0.0, 0.0, 0
    for dxyz, u, origin in blocks(work).values():
        x = [origin[d] + (np.arange(u.shape[d]) + 0.5) * dxyz[d] for d in range(3)]
        xx = np.meshgrid(*x, indexing="ij")
        exact = np.ones_like(u)
        for d in range(3):
            if not null[d]:
                exact *= np.sin(2.0 * math.pi * (xx[d] - a[d] * t - emin[d]) / (emax[d] - emin[d]))
        dv = float(np.prod(np.where(null, 1.0, dxyz)))
        e = np.abs(u - exact)
        e1 += float(e.sum()) * dv
        vol += e.size * dv
        einf = max(einf, float(e.max()))
        n += e.size
    return e1 / vol, einf, n


def convergence(args: argparse.Namespace) -> int:
    rows = [(w, *errors(Path(w))) for w in args.work]
    print(f"{'cells':>8} {'L1':>12} {'Linf':>12} {'order L1':>9} {'order Linf':>10}")
    orders = []
    for k, (_, e1, einf, n) in enumerate(rows):
        if k == 0:
            print(f"{n:8d} {e1:12.4e} {einf:12.4e}")
            continue
        _, p1, pinf, pn = rows[k - 1]
        r = (n / pn) ** (1.0 / args.dims)
        o1, oinf = math.log(p1 / e1) / math.log(r), math.log(pinf / einf) / math.log(r)
        orders.append(o1)
        print(f"{n:8d} {e1:12.4e} {einf:12.4e} {o1:9.3f} {oinf:10.3f}")
    ok = True
    for w in args.work:
        t, t_max = final_time(Path(w)), float(read_input(Path(w))["time"]["time_max"])
        if abs(t - t_max) > 1e-12 * t_max:
            print(f"FAIL: {w} ended at t = {t!r}, not at time_max = {t_max!r}")
            ok = False
    if args.min_order is not None and orders[-1] < args.min_order:
        print(f"FAIL: finest-pair L1 order {orders[-1]:.3f} < {args.min_order}")
        ok = False
    if args.max_l1 is not None and rows[-1][1] > args.max_l1:
        print(f"FAIL: finest L1 error {rows[-1][1]:.3e} > {args.max_l1:.3e}")
        ok = False
    return 0 if ok else 1


def conservation(args: argparse.Namespace) -> int:
    work = Path(args.work)
    _, emin, emax, null = domain(work)
    rows = [ln.split() for ln in (work / f"{basename(work)}-history.dat").read_text().splitlines() if ln[:1] != "#"]
    total = np.array([float(r[2]) for r in rows])
    scale = float(np.prod(np.where(null, 1.0, emax - emin)))  # max|u0| = 1 times the volume
    drift = float(np.max(np.abs(total - total[0]))) / scale
    print(f"total of u dV: {len(rows)} rows, largest drift {drift:.3e} of max|u0| times the volume")
    if drift > args.tol:
        print(f"FAIL: drift {drift:.3e} > {args.tol:.1e}")
        return 1
    return 0


def compare(args: argparse.Namespace) -> int:
    b1, b2 = blocks(Path(args.work[0])), blocks(Path(args.work[1]))
    if set(b1) != set(b2):
        print(f"FAIL: the runs hold different blocks ({len(b1)} and {len(b2)})")
        return 1
    diff = max(float(np.max(np.abs(b1[k][1] - b2[k][1]))) for k in b1)
    print(f"{len(b1)} blocks, largest |u_1 - u_2| = {diff:.3e}")
    if diff > args.tol:
        print(f"FAIL: {diff:.3e} > {args.tol:.1e}")
        return 1
    return 0


def refined(args: argparse.Namespace) -> int:
    """The base blocks the run refined against those where u0 rises above the threshold at the base resolution."""
    work = Path(args.work)
    ini = read_input(work)
    _, emin, emax, null = domain(work)
    nb = args.base_blocks
    ncell = [int(ini["grid"][f"n{c}"]) for c in "ij"]
    size = (emax - emin) / nb
    blks = blocks(work)
    coarsest = max(float(v[0][0]) for v in blks.values())
    actual = set()
    for dxyz, _, origin in blks.values():
        if float(dxyz[0]) < coarsest:
            actual.add(tuple(int(np.floor((origin[d] - emin[d]) / size[d] + 1e-9)) for d in range(2)))
    expected = set()
    for ix in range(nb):
        for iy in range(nb):
            xs = emin[0] + size[0] * (ix + (np.arange(ncell[0]) + 0.5) / ncell[0])
            ys = emin[1] + size[1] * (iy + (np.arange(ncell[1]) + 0.5) / ncell[1])
            u0 = np.outer(np.sin(2.0 * math.pi * (xs - emin[0]) / (emax[0] - emin[0])),
                          np.sin(2.0 * math.pi * (ys - emin[1]) / (emax[1] - emin[1])))
            if u0.max() > args.threshold:
                expected.add((ix, iy))
    print(f"refined base blocks: {sorted(actual)}; expected (u0 > {args.threshold}): {sorted(expected)}")
    if actual != expected:
        print("FAIL: the refined blocks are not those where u0 rises above the threshold")
        return 1
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("convergence")
    p.add_argument("work", nargs="+")
    p.add_argument("--dims", type=int, default=2, help="non-null axes (cells ratio -> spacing ratio)")
    p.add_argument("--min-order", type=float)
    p.add_argument("--max-l1", type=float)
    p.set_defaults(func=convergence)
    p = sub.add_parser("conservation")
    p.add_argument("work")
    p.add_argument("--tol", type=float, default=1e-13)
    p.set_defaults(func=conservation)
    p = sub.add_parser("compare")
    p.add_argument("work", nargs=2)
    p.add_argument("--tol", type=float, default=1e-12)
    p.set_defaults(func=compare)
    p = sub.add_parser("refined")
    p.add_argument("work")
    p.add_argument("--threshold", type=float, required=True)
    p.add_argument("--base-blocks", type=int, default=4, help="base blocks per non-null axis")
    p.set_defaults(func=refined)
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
