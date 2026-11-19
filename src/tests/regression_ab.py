#!/usr/bin/env python3
"""Compare two runs of the same regression case bit for bit: the A/B check of a refactoring (SELAH M7, ADAM #132).

Why this script exists: the regression harnesses (run.sh) compare a run with its goldens within rtol = 1e-6, so they
cannot tell an identical result from one that moved at round-off. A refactoring that must not change the arithmetic
(ADAM's integrators replaced by SELAH's, #132) needs the stronger statement: the same case run by the old executable
(A) and by the new one (B) gives the same bits. This script compares the outputs of the two run directories:
- every history log `*.dat` both runs wrote (residuals, conservation, div(B), turbulence): byte for byte; on a
  difference, the first row that differs and the largest absolute and relative deviation of its numbers;
- every checkpoint `<basename>-<step>-proc<rank>.h5` both runs wrote: every dataset, bit for bit (NaN equal to NaN);
  on a difference, the dataset and the largest deviations.
A file only one run wrote is a difference. Exit 0 when everything is identical, 1 otherwise, 2 on a usage error.

With `--ngc N` a block dataset (rank 3, ghosts included) that differs only in its N ghost layers is not a difference:
it is reported as a note. Two runs of the same FNL executable store different values in some ghost cells of a null
axis (measured on develop 61fa3344, #132 A1: 3 of 56 FLUME cases, histories and interiors bitwise; q and q_gpu are
zeroed at allocation, so these cells are written, with run-dependent values: cause not identified); the interior and
the histories are what a run computes.

Usage:
    regression_ab.py RUN_A RUN_B [--ngc N]
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

CHECKPOINT = re.compile(r"-[0-9]+-proc[0-9]+\.h5$")


def numbers(line: str) -> list[float]:
    """Return the numbers of a history row (an unparsable token is skipped)."""
    out = []
    for tok in line.split():
        try:
            out.append(float(tok))
        except ValueError:
            pass
    return out


def compare_log(a: Path, b: Path) -> bool:
    """Compare two history logs byte for byte; report the first differing row and the largest deviations."""
    ta, tb = a.read_bytes(), b.read_bytes()
    if ta == tb:
        return True
    la, lb = ta.decode(errors="replace").splitlines(), tb.decode(errors="replace").splitlines()
    first = None
    max_abs = max_rel = 0.0
    for n, (ra, rb) in enumerate(zip(la, lb)):
        if ra == rb:
            continue
        if first is None:
            first = n
        for x, y in zip(numbers(ra), numbers(rb)):
            d = abs(x - y)
            max_abs = max(max_abs, d)
            if d > 0.0:
                max_rel = max(max_rel, d / max(abs(x), abs(y)))
    rows = f"{len(la)} vs {len(lb)} rows" if len(la) != len(lb) else f"{len(la)} rows"
    where = f"first differing row {first} ({rows})" if first is not None else f"row counts differ ({rows})"
    print(f"DIFF {a.name}: {where}; largest deviation abs {max_abs:.3e}, rel {max_rel:.3e}")
    if first is not None:
        print(f"     A: {la[first][:160]}")
        print(f"     B: {lb[first][:160]}")
    return False


def compare_h5(a: Path, b: Path, ngc: int | None = None) -> tuple[bool, int]:
    """Compare every dataset of two checkpoints bit for bit; with ngc, a rank-3 dataset that differs in its ghost layers
    only is counted apart. Return (identical, datasets differing in the ghosts only)."""
    import h5py
    import numpy as np

    same = True
    ghost_only = 0
    with h5py.File(a, "r") as fa, h5py.File(b, "r") as fb:
        names_a, names_b = [], []
        fa.visititems(lambda n, o: names_a.append(n) if isinstance(o, h5py.Dataset) else None)
        fb.visititems(lambda n, o: names_b.append(n) if isinstance(o, h5py.Dataset) else None)
        if sorted(names_a) != sorted(names_b):
            print(f"DIFF {a.name}: the datasets differ ({len(names_a)} vs {len(names_b)})")
            return False, 0
        worst = (0.0, 0.0, "")
        ndiff = 0
        for name in sorted(names_a):
            xa, xb = fa[name][()], fb[name][()]
            if np.shape(xa) != np.shape(xb):
                print(f"DIFF {a.name}:{name}: shapes {np.shape(xa)} vs {np.shape(xb)}")
                same = False
                continue
            if np.asarray(xa).dtype.kind in "fc":
                eq = np.array_equal(xa, xb, equal_nan=True)
            else:
                eq = np.array_equal(xa, xb)
            if eq:
                continue
            if ngc and np.ndim(xa) == 3 and min(np.shape(xa)) > 2 * ngc:
                core = (slice(ngc, -ngc),) * 3
                if np.array_equal(np.asarray(xa)[core], np.asarray(xb)[core], equal_nan=True):
                    ghost_only += 1
                    continue
            ndiff += 1
            same = False
            if np.asarray(xa).dtype.kind in "fc":
                d = np.abs(np.asarray(xa, dtype=float) - np.asarray(xb, dtype=float))
                m = float(np.nanmax(d))
                scale = float(np.nanmax(np.abs(np.asarray(xa, dtype=float)))) or 1.0
                if m > worst[0]:
                    worst = (m, m / scale, name)
        if ndiff:
            print(f"DIFF {a.name}: {ndiff} of {len(names_a)} datasets differ; largest abs {worst[0]:.3e} "
                  f"(rel to the dataset max {worst[1]:.3e}) in {worst[2]}")
        if ghost_only:
            print(f"NOTE {a.name}: {ghost_only} datasets differ in their ghost cells only (--ngc {ngc})")
    return same, ghost_only


def main(argv: list[str]) -> int:
    """Compare the two run directories; print a one-line verdict."""
    ngc = None
    if len(argv) == 5 and argv[3] == "--ngc":
        ngc = int(argv[4])
        argv = argv[:3]
    if len(argv) != 3:
        print(__doc__.strip().splitlines()[-1], file=sys.stderr)
        return 2
    ra, rb = Path(argv[1]), Path(argv[2])
    if not (ra.is_dir() and rb.is_dir()):
        print(f"regression_ab: {ra} and {rb} must be run directories", file=sys.stderr)
        return 2
    files_a = {p.name for p in ra.iterdir() if p.suffix == ".dat" or CHECKPOINT.search(p.name)}
    files_b = {p.name for p in rb.iterdir() if p.suffix == ".dat" or CHECKPOINT.search(p.name)}
    ok = True
    for name in sorted(files_a ^ files_b):
        print(f"DIFF {name}: written by {'A' if name in files_a else 'B'} only")
        ok = False
    common = sorted(files_a & files_b)
    if not common:
        print("regression_ab: no history log and no checkpoint to compare")
        return 1
    ghosts = 0
    for name in common:
        if name.endswith(".dat"):
            same = compare_log(ra / name, rb / name)
        else:
            same, g = compare_h5(ra / name, rb / name, ngc)
            ghosts += g
        ok = ok and same
    nlog = sum(1 for n in common if n.endswith(".dat"))
    note = f" (ghost cells differ in {ghosts} datasets)" if ghosts else ""
    print(f"{'BITWISE' if ok else 'DIFFERENT'}: {nlog} history logs, {len(common) - nlog} checkpoints compared{note}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
