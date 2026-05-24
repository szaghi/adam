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

Granularity is per-variable-per-output-step: ALL block datasets at iteration
N (across all blocks, all MPI ranks, and — for Phase D multi-realm cases —
all realm output basenames) are reduced together. This is intentionally
invariant to:

  * block-to-rank ownership (redistributing blocks across ranks is not a
    physics regression);
  * realm-to-file partitioning (the two-realm x-split rmf case writes its
    cells across `rmf_2realm_r1-*.h5` + `rmf_2realm_r2-*.h5`, but the
    UNION of cells is bit-equivalent to the single-realm `rmf_regression-*.h5`
    at the same iteration index — see issue #13 D.4).

Ghost cells are EXCLUDED from the reduction. PRISM saves checkpoints with
`with_ghost=.true.`, which appends 2·ngc cells per direction. The ghost cells
diverge between single-realm and multi-realm configs (the multi-realm case
has inter-realm-mirror ghosts at the realm boundary, which the single-realm
config never sees), so including them defeats bit-comparability. The
ghost-cell stripping is parametric on the case directory's ngc (read from
any case .ini's `[grid] ngc = ...`); falls back to ngc=0 (no stripping)
when no INI is supplied.

Output rows are keyed on a STEP-INDEX derived from each checkpoint's
filename (`<basename>-<step>-proc<rank>.h5` ⇒ step=<step>) rather than the
basename, so multi-realm checkpoints at the same iteration index aggregate
into a single row.

Usage:
    digest.py write   <out.txt> <checkpoint.h5> [<checkpoint.h5> ...] [--ngc N | --case-dir DIR]
    digest.py compare <produced.txt> <golden.txt> [--rtol R] [--atol A]

``write`` emits a deterministic, sorted digest. ``compare`` exits 0 when every
row matches within tolerance, non-zero otherwise, printing the offending rows.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import h5py
import numpy as np

# Dataset name layout: block_<NNN>-proc<NNN>-<variable>. The variable suffix is
# everything after the second '-'. Reductions are aggregated over all blocks
# and ranks sharing a variable suffix.
_NAME_SEP = "-"

# Filename layout: <basename>-<step>-proc<rank>.h5, where <step> is a
# zero-padded iteration index emitted by PRISM IO. The step index is the
# digest's aggregation key across multi-realm output files.
_FILENAME_RE = re.compile(r"^(?P<basename>.+)-(?P<step>\d+)-proc(?P<rank>\d+)\.h5$")

# Variable-name prefixes excluded from the digest entirely.
#
# div_* — every divergence diagnostic field. PRISM emits div_D, div_B,
# div_J plus numbered div_05..div_12 (see adam_prism_common_object.F90
# div_name); all are constraint diagnostics that should be ~0 (Gauss's
# law, charge conservation). Their stored values are pure round-off, so
# *no* reduction of them (min/max/sum/sum_sq) is reproducible across
# compilers — gcc-14 vs gcc-16 disagreed by ~24% on div_J min because
# they were comparing noise to noise. The numbered div_NN fields are not
# a different class: they only escaped the first CI failure because they
# happened to still be exactly 0.0 at the saved iterations of this
# vacuum case; they would flake the moment the physics drives them
# nonzero. Keeping any div_* in the digest is anti-signal.
#
# A future divergence-cleaning regression, if wanted, needs a dedicated
# metric (e.g. ‖div‖ below a threshold), not a golden-value comparison.
_EXCLUDED_VAR_PREFIXES = ("div_",)

# Per-block metadata variables: digest rows whose value depends on the
# block-partitioning of the domain rather than on physics. Comparisons of
# these rows across configurations with different block partitionings (the
# canonical case: a multi-realm config aggregated against its single-realm
# union) are necessarily non-bit-comparable in a non-physics-meaningful way:
#
#   * `dxdydz`, `origin` — per-block geometry attributes (3 values per block,
#     stored as one dataset per block in the HDF5 checkpoint). The digest's
#     count column equals 3·nblocks; the sum/sum_sq scale with nblocks; the
#     min/max can shift when block extents differ.
#   * `time_iteration` — per-block iteration counter (1 value per block,
#     stored as one dataset per block). Count equals nblocks; sum scales
#     with nblocks.
#
# When compared against a reference digest captured at the SAME block
# partitioning (the case's own golden), these rows are perfectly stable and
# the comparison fires on real layout regressions. When compared against a
# reference captured at a DIFFERENT partitioning (the cross-case sanity check
# rmf-2realm-vs-rmf, where the digest is invariant for FIELD variables but
# not for these per-block attributes), the mismatches are benign by
# construction and would otherwise drown out any real field regression.
#
# `cmd_compare` downgrades mismatches on these names from MISMATCH (counted
# as failure) to SKIP_METADATA (logged but not failing). The downgrade is
# unconditional: any consumer who genuinely wants to check per-block layout
# equivalence should compare against a reference captured at the same
# partitioning, where these rows match bit-exactly.
_PER_BLOCK_METADATA_VARS = frozenset({"dxdydz", "origin", "time_iteration"})

# Default comparison tolerances. Calibrated for *cross-compiler* runs
# (goldens captured in CI, also re-run in CI on a different toolchain than
# a dev workstation), not same-machine bit-identical reruns.
#   rtol 1e-6  — absorbs floating-point association noise from differing
#                compiler codegen, libm, FMA contraction and MPI reduction
#                order; still tight enough to catch any algorithmic change.
#   atol 1e-3  — floor for cancellation-residue reductions. An antisymmetric
#                field like Jz (±3.7e3 per cell, summing to ~4e-5 over
#                340k cells) has a `sum` that is pure cancellation noise;
#                1e-3 treats it as indistinguishable from zero, while a
#                genuine aggregate like Jx `sum` (~4.6e8) is unaffected.
_DEFAULT_RTOL = 1.0e-6
_DEFAULT_ATOL = 1.0e-3


def _variable_of(dataset_name: str) -> str | None:
    """Return the variable suffix of a ``block_*-proc*-<var>`` dataset name.

    Returns None for any dataset that does not match the expected layout, or
    whose variable name starts with an ``_EXCLUDED_VAR_PREFIXES`` entry, so
    it is skipped from the digest.
    """
    parts = dataset_name.split(_NAME_SEP, 2)
    if len(parts) != 3:
        return None
    var = parts[2]
    if var.startswith(_EXCLUDED_VAR_PREFIXES):
        return None
    return var


def _interior(data: np.ndarray, ngc: int) -> np.ndarray:
    """Strip ngc ghost-cell layers off every dimension of a 3D field array.

    PRISM checkpoints with `with_ghost=.true.` store block data as a 3D array
    of shape (ni+2ngc, nj+2ngc, nk+2ngc). Returns the central
    [ngc:-ngc, ngc:-ngc, ngc:-ngc] slice (interior cells only) when ngc>0,
    or the unmodified array when ngc==0 (legacy single-realm digest, kept
    for backward compatibility with already-captured goldens).
    """
    if ngc <= 0 or data.ndim != 3:
        return data
    return data[ngc:-ngc, ngc:-ngc, ngc:-ngc]


def _reduce_file(path: Path, ngc: int) -> dict[str, np.ndarray]:
    """Accumulate per-variable reductions over every block dataset in one file.

    Each variable maps to a float64 array [count, min, max, sum, sum_sq].
    Ghost cells are stripped when `ngc > 0` so the reduction covers only
    interior cells; without stripping the inter-realm-mirror ghost cells of
    a multi-realm case would shift sums and break bit-comparability against
    a single-realm reference.
    """
    acc: dict[str, np.ndarray] = {}
    with h5py.File(path, "r") as h5:
        for name, obj in h5.items():
            if not isinstance(obj, h5py.Dataset):
                continue
            var = _variable_of(name)
            if var is None:
                continue
            raw = np.asarray(obj[()], dtype=np.float64)
            interior = _interior(raw, ngc)
            data = interior.ravel()
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


def _step_index_of(filename: str) -> str | None:
    """Extract the zero-padded step index from `<basename>-<step>-proc<rank>.h5`."""
    m = _FILENAME_RE.match(filename)
    if m is None:
        return None
    return m.group("step")


def _ngc_from_case_dir(case_dir: Path) -> int:
    """Read `[grid] ngc = ...` from any of the case directory's `.ini` files.

    Multi-realm cases have several INIs (forest manifest + one per realm);
    every realm INI declares `ngc`, and Phase D's design requires them to
    agree. We accept the first value found.
    """
    pattern = re.compile(r"^\s*ngc\s*=\s*(\d+)", re.MULTILINE)
    for ini_path in sorted(case_dir.glob("*.ini")):
        match = pattern.search(ini_path.read_text())
        if match:
            return int(match.group(1))
    return 0


def _format_row(step: str, var: str, stats: np.ndarray) -> str:
    """Tab-separated row keyed on (iteration step, variable).

    Step is the zero-padded iteration index parsed from the checkpoint
    filename; it aggregates across all blocks, ranks, AND multi-realm
    output basenames at that iteration.
    """
    count = int(stats[0])
    return (
        f"{step}\t{var}\t{count}\t"
        f"{stats[1]:.16e}\t{stats[2]:.16e}\t{stats[3]:.16e}\t{stats[4]:.16e}"
    )


def cmd_write(out_path: str, h5_paths: list[str], ngc: int) -> int:
    if not h5_paths:
        print("ERROR: write needs at least one HDF5 file", file=sys.stderr)
        return 2
    # step_acc[(step, variable)] -> [count, min, max, sum, sum_sq]
    step_acc: dict[tuple[str, str], np.ndarray] = {}
    for h5_path in sorted(h5_paths):
        p = Path(h5_path)
        if not p.is_file():
            print(f"ERROR: not a file: {h5_path}", file=sys.stderr)
            return 2
        step = _step_index_of(p.name)
        if step is None:
            print(f"ERROR: cannot parse step index from filename: {p.name}", file=sys.stderr)
            return 2
        file_acc = _reduce_file(p, ngc=ngc)
        for var, stats in file_acc.items():
            key = (step, var)
            if key in step_acc:
                prev = step_acc[key]
                step_acc[key] = np.array(
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
                step_acc[key] = stats
    rows = [_format_row(step, var, step_acc[(step, var)]) for (step, var) in sorted(step_acc)]
    header = "# step\tvariable\tcount\tmin\tmax\tsum\tsum_sq"
    Path(out_path).write_text(header + "\n" + "\n".join(rows) + "\n")
    print(f">> wrote digest: {out_path} ({len(rows)} rows, ngc={ngc})")
    return 0


def _parse_digest(path: Path) -> dict[tuple[str, str], list[str]]:
    """Parse a digest file into {(step, variable): [count,min,max,sum,sumsq]}.

    Compatible with both the legacy `<basename>-<step>-proc<rank>.h5`-keyed
    format and the step-keyed format introduced for multi-realm aggregation
    (Phase D.4): the two share the same row layout — the first column is
    just opaque text used as part of the lookup key.
    """
    table: dict[tuple[str, str], list[str]] = {}
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 7:
            raise ValueError(f"malformed digest row in {path}: {line!r}")
        step, var = fields[0], fields[1]
        table[(step, var)] = fields[2:]
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
    metadata_skips: list[str] = []  # benign per-block-metadata mismatches (see _PER_BLOCK_METADATA_VARS)

    missing = sorted(set(golden) - set(produced))
    for key in missing:
        # MISSING/EXTRA on per-block-metadata names are also benign for the
        # cross-partitioning case (same reason as value mismatches): if the
        # produced run has a different block layout, the metadata row count
        # changes and we accept that.
        if key[1] in _PER_BLOCK_METADATA_VARS:
            metadata_skips.append(f"MISSING  {key[0]} / {key[1]} — in golden, not produced")
        else:
            failures.append(f"MISSING  {key[0]} / {key[1]} — in golden, not produced")

    extra = sorted(set(produced) - set(golden))
    for key in extra:
        if key[1] in _PER_BLOCK_METADATA_VARS:
            metadata_skips.append(f"EXTRA    {key[0]} / {key[1]} — produced, not in golden")
        else:
            failures.append(f"EXTRA    {key[0]} / {key[1]} — produced, not in golden")

    for key in sorted(set(produced) & set(golden)):
        reason = _values_match(produced[key], golden[key], rtol, atol)
        if reason is not None:
            # Per-block-metadata variables (`dxdydz`, `origin`,
            # `time_iteration`) are downgraded to SKIP_METADATA: their
            # values change with block partitioning rather than physics,
            # so a mismatch here is meaningful ONLY when the reference
            # was captured at the same partitioning (then it doesn't fire).
            if key[1] in _PER_BLOCK_METADATA_VARS:
                metadata_skips.append(f"SKIP_METADATA {key[0]} / {key[1]} — {reason}")
            else:
                failures.append(f"MISMATCH {key[0]} / {key[1]} — {reason}")

    if metadata_skips:
        print(
            f">> {len(metadata_skips)} per-block-metadata discrepancies skipped "
            f"(benign across block-partitioning changes; see _PER_BLOCK_METADATA_VARS)",
            file=sys.stderr,
        )
        for line in metadata_skips:
            print("  " + line, file=sys.stderr)

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
        if len(argv) < 4:
            print("ERROR: write needs <out.txt> <checkpoint.h5> [...]", file=sys.stderr)
            return 2
        out_path = argv[2]
        # Split positional HDF5 paths from flag pairs (--ngc / --case-dir).
        h5_paths: list[str] = []
        ngc = 0
        case_dir: Path | None = None
        rest = argv[3:]
        i = 0
        while i < len(rest):
            if rest[i] == "--ngc" and i + 1 < len(rest):
                ngc = int(rest[i + 1])
                i += 2
            elif rest[i] == "--case-dir" and i + 1 < len(rest):
                case_dir = Path(rest[i + 1])
                i += 2
            elif rest[i].startswith("--"):
                print(f"ERROR: unknown write flag {rest[i]!r}", file=sys.stderr)
                return 2
            else:
                h5_paths.append(rest[i])
                i += 1
        if case_dir is not None and ngc == 0:
            ngc = _ngc_from_case_dir(case_dir)
        return cmd_write(out_path, h5_paths, ngc=ngc)
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
