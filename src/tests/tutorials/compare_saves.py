"""Compare two saves: max |q_b - q_a| / max |q_a| per variable, over every block (tutorials, issue #60).

Usage: python compare_saves.py WORK BASENAME IT_A IT_B [WORK_B]

WORK holds the run's <BASENAME>-<it:9>-procNNNNNN.h5 field files (FLUME, PRISM and RIPPLE write them alike). Save IT_B
is read from WORK_B when given (two runs of the same input, e.g. CPU and FNL), else from WORK. Blocks are matched by
their origin, so the two saves must hold the same grid.
"""

import glob
import sys

import h5py
import numpy as np

META = ("-origin", "-dxdydz", "-time_iteration")


def load(work: str, base: str, it: int) -> dict[tuple, np.ndarray]:
    """Every variable of every block of one save, keyed by (block origin, variable name)."""
    out = {}
    for path in sorted(glob.glob(f"{work}/{base}-{it:09d}-proc*.h5")):
        with h5py.File(path, "r") as h5:
            for key in h5:
                if key.endswith(META):
                    continue
                block, var = key.rsplit("-", 1)
                origin = tuple(np.round(h5[f"{block}-origin"][()], 12))
                out[(origin, var)] = h5[key][()]
    return out


def main() -> int:
    work, base, it_a, it_b = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    work_b = sys.argv[5] if len(sys.argv) > 5 else work
    a, b = load(work, base, it_a), load(work_b, base, it_b)
    if not a or a.keys() != b.keys():
        sys.exit("compare_saves.py: the two saves are missing or hold different grids")
    for var in sorted({v for _, v in a}):
        num = max(np.abs(b[k] - a[k]).max() for k in a if k[1] == var)
        den = max(np.abs(a[k]).max() for k in a if k[1] == var)
        if den > 0:
            print(f"{var:8s} max|q({it_b}) - q({it_a})| / max|q({it_a})| = {num / den:.3e}")
        else:  # zero in the first save: report the absolute difference, so a spurious field is not skipped
            print(f"{var:8s} max|q({it_b}) - q({it_a})| = {num:.3e} (zero in q({it_a}))")
    return 0


if __name__ == "__main__":
    sys.exit(main())
