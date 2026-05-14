#!/usr/bin/env python3
"""PRISM regression field digest.

Reduces every HDF5 checkpoint to a compact, committable text fingerprint and
compares a freshly produced fingerprint against a committed golden one.

Why a digest instead of the raw HDF5: a full PRISM checkpoint is ~100 MB per
rank — far too large to commit as a regression golden. The per-iteration
residuals log catches global-norm regressions, but a residual is an aggregate
and can mask a sign-cancelling regional error. This digest closes that gap: it
records point-wise reductions (count, min, max, sum, sum of squares) per field
variable, which shift if any single cell changes, yet stay in the kilobyte
range.

Granularity is per-variable-per-file: the 32 blocks of each variable (e.g. all
``block_*-proc*-Bx`` datasets) are reduced together. This is intentionally
invariant to block-to-rank ownership — redistributing blocks across ranks is
not a physics regression — while remaining sensitive to any change in the
field values themselves.

Usage:
    digest.py write  <out.txt> <checkpoint.h5> [<checkpoint.h5> ...]
    digest.py compare <produced.txt> <golden.txt> [--rtol R] [--atol A]

``write`` emits a deterministic, sorted digest. ``compare`` exits 0 when every
row matches within tolerance, non-zero otherwise, printing the offending rows.
"""

from __future__ import annotations

import sys
from pathlib import Path

import h5py
import numpy as np

# Dataset name layout: block_<NNN>-proc<NNN>-<variable>. The variable suffix is
# everything after the second '-'. Reductions are aggregated over all blocks
# and ranks sharing a variable suffix.
_NAME_SEP = "-"

# Default comparison tolerances. sum / sum_sq accumulate floating-point
# association noise across compilers, MPI rank counts and FMA contraction;
# rtol absorbs that while staying tight enough to catch a real regression.
# count / min / max are exact-ish (min/max get a tiny atol for denormal noise).
_DEFAULT_RTOL = 1.0e-11
_DEFAULT_ATOL = 1.0e-13


def _variable_of(dataset_name: str) -> str | None:
    """Return the variable suffix of a ``block_*-proc*-<var>`` dataset name.

    Returns None for any dataset that does not match the expected layout, so
    unexpected top-level datasets are skipped rather than crashing the digest.
    """
    parts = dataset_name.split(_NAME_SEP, 2)
    if len(parts) != 3:
        return None
    return parts[2]


def _reduce_file(path: Path) -> dict[str, np.ndarray]:
    """Accumulate per-variable reductions over every block dataset in one file.

    Each variable maps to a float64 array [count, min, max, sum, sum_sq].
    """
    acc: dict[str, np.ndarray] = {}
    with h5py.File(path, "r") as h5:
        for name, obj in h5.items():
            if not isinstance(obj, h5py.Dataset):
                continue
            var = _variable_of(name)
            if var is None:
                continue
            data = np.asarray(obj[()], dtype=np.float64).ravel()
            if data.size == 0:
                continue
            stats = np.array(
                [
                    float(data.size),
                    float(data.min()),
                    float(data.max()),
                    float(data.sum()),
                    float(np.square(data).sum()),
                ],
                dtype=np.float64,
            )
            if var in acc:
                prev = acc[var]
                acc[var] = np.array(
                    [
                        prev[0] + stats[0],
                        min(prev[1], stats[1]),
                        max(prev[2], stats[2]),
                        prev[3] + stats[3],
                        prev[4] + stats[4],
                    ],
                    dtype=np.float64,
                )
            else:
                acc[var] = stats
    return acc


def _format_row(checkpoint: str, var: str, stats: np.ndarray) -> str:
    count = int(stats[0])
    return (
        f"{checkpoint}\t{var}\t{count}\t"
        f"{stats[1]:.16e}\t{stats[2]:.16e}\t{stats[3]:.16e}\t{stats[4]:.16e}"
    )


def cmd_write(out_path: str, h5_paths: list[str]) -> int:
    if not h5_paths:
        print("ERROR: write needs at least one HDF5 file", file=sys.stderr)
        return 2
    rows: list[str] = []
    for h5_path in sorted(h5_paths):
        p = Path(h5_path)
        if not p.is_file():
            print(f"ERROR: not a file: {h5_path}", file=sys.stderr)
            return 2
        # Key the digest on the basename only: the work directory path is
        # volatile, the checkpoint file name is the stable identity.
        checkpoint = p.name
        acc = _reduce_file(p)
        for var in sorted(acc):
            rows.append(_format_row(checkpoint, var, acc[var]))
    header = "# checkpoint\tvariable\tcount\tmin\tmax\tsum\tsum_sq"
    Path(out_path).write_text(header + "\n" + "\n".join(rows) + "\n")
    print(f">> wrote digest: {out_path} ({len(rows)} rows)")
    return 0


def _parse_digest(path: Path) -> dict[tuple[str, str], list[str]]:
    """Parse a digest file into {(checkpoint, variable): [count,min,max,sum,sumsq]}."""
    table: dict[tuple[str, str], list[str]] = {}
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 7:
            raise ValueError(f"malformed digest row in {path}: {line!r}")
        checkpoint, var = fields[0], fields[1]
        table[(checkpoint, var)] = fields[2:]
    return table


def _values_match(produced: list[str], golden: list[str], rtol: float, atol: float) -> str | None:
    """Return None if all five reductions match within tolerance, else a reason."""
    # count is integral and must match exactly.
    if produced[0] != golden[0]:
        return f"count {produced[0]} != {golden[0]}"
    labels = ("min", "max", "sum", "sum_sq")
    for idx, label in enumerate(labels, start=1):
        p = float(produced[idx])
        g = float(golden[idx])
        if not np.isclose(p, g, rtol=rtol, atol=atol):
            return f"{label} {p:.16e} != {g:.16e}"
    return None


def cmd_compare(produced_path: str, golden_path: str, rtol: float, atol: float) -> int:
    produced_p = Path(produced_path)
    golden_p = Path(golden_path)
    if not produced_p.is_file():
        print(f"ERROR: produced digest missing: {produced_path}", file=sys.stderr)
        return 2
    if not golden_p.is_file():
        print(f"ERROR: golden digest missing: {golden_path}", file=sys.stderr)
        return 2

    produced = _parse_digest(produced_p)
    golden = _parse_digest(golden_p)

    failures: list[str] = []

    missing = sorted(set(golden) - set(produced))
    for key in missing:
        failures.append(f"MISSING  {key[0]} / {key[1]} — in golden, not produced")

    extra = sorted(set(produced) - set(golden))
    for key in extra:
        failures.append(f"EXTRA    {key[0]} / {key[1]} — produced, not in golden")

    for key in sorted(set(produced) & set(golden)):
        reason = _values_match(produced[key], golden[key], rtol, atol)
        if reason is not None:
            failures.append(f"MISMATCH {key[0]} / {key[1]} — {reason}")

    if failures:
        print(f"FAIL: digest comparison ({len(failures)} discrepancies)", file=sys.stderr)
        for line in failures:
            print("  " + line, file=sys.stderr)
        return 1

    print(f">> digest match: {len(produced)} rows within rtol={rtol:g} atol={atol:g}")
    return 0


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(__doc__, file=sys.stderr)
        return 2
    mode = argv[1]
    if mode == "write":
        return cmd_write(argv[2], argv[3:]) if len(argv) >= 4 else cmd_write("", [])
    if mode == "compare":
        if len(argv) < 4:
            print("ERROR: compare needs <produced.txt> <golden.txt>", file=sys.stderr)
            return 2
        rtol = _DEFAULT_RTOL
        atol = _DEFAULT_ATOL
        rest = argv[4:]
        i = 0
        while i < len(rest):
            if rest[i] == "--rtol" and i + 1 < len(rest):
                rtol = float(rest[i + 1])
                i += 2
            elif rest[i] == "--atol" and i + 1 < len(rest):
                atol = float(rest[i + 1])
                i += 2
            else:
                print(f"ERROR: unknown compare argument {rest[i]!r}", file=sys.stderr)
                return 2
        return cmd_compare(argv[2], argv[3], rtol, atol)
    print(f"ERROR: unknown mode {mode!r} (expected 'write' or 'compare')", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv))
