#!/usr/bin/env python3
"""Check that the restart snapshot places its blocks where the regular field files do (issue #109).

Every restart also writes `<restart_basename>.xdmf` + `-procNNNNNN.h5`, the interior cells only. Its block origin must
be the corner of the first stored cell: the origin of the same block in a regular field file of the same step (ghost
cells included, origin at the ghost corner) plus `ngc` cells. The arrays are in (z, y, x) order, and so are `origin`
and `dxdydz`. Before #109 the snapshot also wrote the ghost-corner origin, so every reader placed it `ngc` cells off.

Usage: snapshot_oracle.py WORK RESTART_BASENAME FIELD_BASENAME NGC
       (FIELD_BASENAME is the regular save of the same step, e.g. sod-x-000000020)
"""

from __future__ import annotations

import sys
from pathlib import Path

import h5py
import numpy as np


def blocks(work: Path, basename: str) -> dict[str, tuple[np.ndarray, np.ndarray, tuple[int, ...]]]:
    """Every block of a save: name -> (origin, dxdydz, shape of the first variable)."""
    out = {}
    for path in sorted(work.glob(f"{basename}-proc*.h5")):
        with h5py.File(path, "r") as h5:
            for key in h5:
                if not key.endswith("-origin"):
                    continue
                name = key[: -len("-origin")]
                var = next(
                    k
                    for k in h5
                    if k.startswith(name + "-") and k.rsplit("-", 1)[1] not in ("origin", "dxdydz", "time_iteration")
                )
                out[name] = (h5[key][()], h5[name + "-dxdydz"][()], h5[var].shape)
    return out


def main() -> int:
    work, restart, field, ngc = Path(sys.argv[1]), sys.argv[2], sys.argv[3], int(sys.argv[4])
    snap, full = blocks(work, restart), blocks(work, field)
    if not snap or snap.keys() != full.keys():
        print(f"FAIL: the snapshot ({len(snap)} blocks) and the field save ({len(full)} blocks) differ")
        return 1
    worst = 0.0
    for name, (origin, dxyz, shape) in snap.items():
        f_origin, f_dxyz, f_shape = full[name]
        if any(s + 2 * ngc != f for s, f in zip(shape, f_shape)):
            print(f"FAIL: {name}: snapshot shape {shape}, field shape {f_shape}: not the interior of it")
            return 1
        worst = max(worst, float(np.max(np.abs(origin - (f_origin + ngc * f_dxyz)) / f_dxyz)))
    verdict = "PASS" if worst < 1e-9 else "FAIL"
    print(f"   restart snapshot: {len(snap)} blocks, origin off the first stored cell by {worst:.3e} cells  {verdict}")
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
