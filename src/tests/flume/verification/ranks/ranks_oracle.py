#!/usr/bin/env python3
"""Oracle of the FLUME empty-rank leg DV-1 (issue #82): block ownership per rank, and np N against np 1 bitwise.

    ranks_oracle.py ownership <work> --empty-at-start [--filled-at-end]
    ranks_oracle.py same <work_np1> <work_npN>

`ownership` counts the blocks of every rank in the first and the last checkpoint of a run (the blocks of the
`<basename>-<step>-proc<rank>.h5` file of each rank). `--empty-at-start` asserts that every rank but rank 0 owns no
block in the first checkpoint, so the leg exercises an empty rank and cannot become vacuous if the distribution changes;
`--filled-at-end` asserts that every rank owns a block in the last one (a runtime regrid repopulated them).

`same` compares two runs of the same input on 1 and on N ranks at their last common checkpoint, cell by cell (matched by
centre) and field by field, and passes only when every value is equal (bitwise). With every block on rank 0 the
arithmetic of the two runs is identical; the reductions add exact zeros from the empty ranks. After a regrid the blocks
are spread, and the scheme's operations are still local to a cell's stencil, which the ghost exchange fills with the
same values whatever rank owns the neighbour.

Exit 0 on pass, 1 on fail.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "gravity"))  # the checkpoint readers
from atmosphere_oracle import FIELDS, checkpoints, interior, read_ini  # noqa: E402


def blocks_per_rank(paths: list[Path]) -> dict[int, int]:
    """Blocks of each rank in one checkpoint (rank -> count), from the per-rank files."""
    out = {}
    for path in paths:
        rank = int(path.stem.rsplit("-proc", 1)[1])
        with h5py.File(path, "r") as h5:
            out[rank] = len({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")})
    return out


def ownership(args: argparse.Namespace) -> int:
    """Report the blocks per rank at the first and last checkpoints; assert the requested emptiness and filling."""
    cps = checkpoints(args.work)
    first, last = min(cps), max(cps)
    start, end = blocks_per_rank(cps[first][1]), blocks_per_rank(cps[last][1])
    print(f"   {args.work.name}: blocks per rank at step {first}: {dict(sorted(start.items()))}, at step {last}: "
          f"{dict(sorted(end.items()))}")
    ok = True
    if len(start) < 2:
        print(f"   {args.work.name}: a single rank, nothing to check  FAIL")
        ok = False
    if args.empty_at_start and any(n for r, n in start.items() if r > 0):
        print(f"   {args.work.name}: a rank other than 0 owns blocks at step {first}: no empty rank exercised  FAIL")
        ok = False
    if args.filled_at_end and not all(end.values()):
        print(f"   {args.work.name}: a rank still owns no block at step {last}: the regrid did not repopulate it  FAIL")
        ok = False
    return 0 if ok else 1


def same(args: argparse.Namespace) -> int:
    """Assert np N = np 1 bitwise at the last common checkpoint."""
    ca, cb = checkpoints(args.a), checkpoints(args.b)
    step = max(set(ca) & set(cb))
    ini = read_ini(args.a)
    fa = {tuple(c.round(12)): q for c, _, q in interior(ini, ca[step][1])}
    fb = {tuple(c.round(12)): q for c, _, q in interior(ini, cb[step][1])}
    if set(fa) != set(fb):
        print(f"{args.a.parent.name}: the two runs hold different cells at step {step} ({len(fa)} and {len(fb)})  FAIL")
        return 1
    worst = 0.0
    for name in FIELDS:
        if name not in next(iter(fa.values())):
            continue
        diff = max(abs(fa[c][name] - fb[c][name]) for c in fa)
        worst = max(worst, diff)
        if diff > 0.0:
            print(f"   {name:>3}: max |np1 - npN| = {diff:.3e}")
    ok = worst == 0.0
    print(f"{args.a.parent.name}: {args.a.name} vs {args.b.name} at step {step}, {len(fa)} cells: "
          f"{'bitwise' if ok else f'max difference {worst:.3e}'}  {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def main() -> int:
    """Parse the command line and run the check."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="check", required=True)
    p = sub.add_parser("ownership")
    p.add_argument("work", type=Path)
    p.add_argument("--empty-at-start", action="store_true", help="ranks > 0 own no block at the first checkpoint")
    p.add_argument("--filled-at-end", action="store_true", help="every rank owns a block at the last checkpoint")
    p = sub.add_parser("same")
    p.add_argument("a", type=Path, help="the run on 1 rank")
    p.add_argument("b", type=Path, help="the run on N ranks")
    args = parser.parse_args()
    return {"ownership": ownership, "same": same}[args.check](args)


if __name__ == "__main__":
    raise SystemExit(main())
