#!/usr/bin/env python3
"""FLUME verification oracle of the output channels: slices and auxiliary fields (issue #35, P7).

Why this oracle exists: the slices (library slices_object, MPI-IO .mat files) and the auxiliary fields of the XH5F
checkpoints are derived outputs: nothing in a run reads them back, so a wrong value would go unnoticed. Two checks on
the last checkpoint of a Sod-x run:

* slices: the slice points lie on the x cell centres, where the trilinear interpolation must return the cell values
  of the checkpoint (the transverse copies are identical, so the transverse interpolation is exact too);
* auxiliary fields: the saved u, v, w, p, H, a must equal the values recomputed from the saved conservative fields.

.mat layout (adam_object%save_slice): 4 int32 (ni, nj, nk, nvar = 3 + nq), then per point 3 float64 coordinates and
nq float64 variables, points in (i, j, k) order, i fastest.

Usage:
    io_oracle.py <work-dir> [--tol T]
"""

from __future__ import annotations

import argparse
import configparser
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


def main() -> int:
    """Run the checks, print a report, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work", type=Path, help="work directory of the Sod-x run with slices and auxiliary fields")
    parser.add_argument("--tol", type=float, default=1.0e-12, help="relative tolerance")
    args = parser.parse_args()

    ini = read_ini(next(args.work.glob("*.ini")))
    ngc = int(ini["grid"]["ngc"])
    gamma = float(ini["physics"]["cp"]) / float(ini["physics"]["cv"])
    basename = ini["IO"]["output_basename"]
    files = sorted(p for p in args.work.glob(f"{basename}-*-proc*.h5"))
    last = max(int(p.name.split("-")[-2]) for p in files)
    status = 0

    # cell values along x (one transverse copy) and auxiliary fields check
    cells: dict[float, np.ndarray] = {}
    aux_err = 0.0
    for path in (p for p in files if int(p.name.split("-")[-2]) == last):
        with h5py.File(path, "r") as h5:
            for blk in {k.rsplit("-", 1)[0] for k in h5}:
                dx = h5[f"{blk}-dxdydz"][()][::-1]
                lo = h5[f"{blk}-origin"][()][::-1] + ngc * dx
                q = np.stack([h5[f"{blk}-{v}"][()].transpose(2, 1, 0) for v in VARIABLES])[:, ngc:-ngc, ngc:-ngc, ngc:-ngc]
                aux = {v: h5[f"{blk}-{v}"][()].transpose(2, 1, 0)[ngc:-ngc, ngc:-ngc, ngc:-ngc]
                       for v in ("u", "v", "w", "p", "H", "a")}
                r = q[0]
                u, v, w = q[1] / r, q[2] / r, q[3] / r
                p = (gamma - 1.0) * (q[4] - 0.5 * r * (u * u + v * v + w * w))
                ref = {"u": u, "v": v, "w": w, "p": p, "H": (q[4] + p) / r, "a": np.sqrt(gamma * p / r)}
                for name, val in ref.items():
                    scale = max(float(np.max(np.abs(val))), 1.0e-300)
                    aux_err = max(aux_err, float(np.max(np.abs(aux[name] - val))) / scale)
                for i in range(q.shape[1]):
                    cells[round(float(lo[0] + (i + 0.5) * dx[0]), 12)] = q[:, i, 0, 0]
    ok = aux_err <= args.tol
    status |= 0 if ok else 1
    print(f"auxiliary fields (u, v, w, p, H, a) vs recomputed: max relative difference {aux_err:.3e}"
          f"  {'PASS' if ok else 'FAIL'} (tol {args.tol:.1e})")

    # slices
    mats = sorted(args.work.glob(f"{basename}-slice_01-*.mat"))
    if not mats:
        sys.exit(f"io_oracle: no {basename}-slice_01-*.mat in {args.work}")
    mat = mats[-1]
    raw = mat.read_bytes()
    ni, nj, nk, nvar = (int(x) for x in np.frombuffer(raw[:16], dtype=np.int32))
    data = np.frombuffer(raw[16:], dtype=np.float64).reshape(ni * nj * nk, nvar)
    scale = np.max([np.abs(c) for c in cells.values()], axis=0)
    scale = np.where(scale > 0.0, scale, 1.0)
    err, missing = 0.0, 0
    for row in data:
        key = round(float(row[0]), 12)
        if key not in cells:
            missing += 1
            continue
        err = max(err, float(np.max(np.abs(row[3:3 + len(VARIABLES)] - cells[key]) / scale)))
    ok = missing == 0 and err <= args.tol
    status |= 0 if ok else 1
    print(f"slice {mat.name}: {ni}x{nj}x{nk} points, {missing} off the cell centres, max relative difference {err:.3e}"
          f"  {'PASS' if ok else 'FAIL'} (tol {args.tol:.1e})")
    return status


if __name__ == "__main__":
    sys.exit(main())
