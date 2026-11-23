#!/usr/bin/env python3
"""Repartition a FLUME restart written on N MPI ranks for a restart on M ranks (stdlib only, issue #133 phase 1).

Why: a restart needs the rank count it was saved with (`tree%load_nodes` stops on a node of a rank that does not
exist, `field%load_blocks` on a missing per-rank file). The speedup matrix warms each leg up once and starts every
backend and rank count from the same bits (#133 D3), so the per-rank files are redistributed here, offline, exactly as
ADAM would have laid the blocks out on M ranks.

The restart files of one realm (`<basename>` = its `[IO] restart_basename`; writers in parentheses):

    <basename>.time               int32 step, float64 time (adam_object%save_restart_files): copied
    <basename>.tnd                int64 n, then n x (int64 code, int32 owner rank, int64 block index into the owner's
                                  .fbd, 1-based), sorted by Morton position (tree_object%save_nodes): REWRITTEN
    <basename>-proc<r:06d>.fbd    int32 nv, ni, nj, nk, ngc, nv_pic, np; int32 blocks; per block int64 code,
                                  int32 i, j, k, level, float64 q(nv, 1-ngc:ni+ngc, ...) with its ghosts
                                  (field_object%save_blocks): SPLIT into M files
    <basename>-proc<r:06d>.wmf    int32 blocks, ni, nj, nk, fresh; per block int64 code and the wall-model filter
                                  averages (flume_common_object%save_wm_filter): SPLIT like the .fbd, when present
    <basename>.reference          the references of the files (flume_reference_object%save_restart): copied
    <basename>.sinks              the replicated sink list (flume_sinks_object%save_restart): copied when present
    <basename>.dff                the digital-filter inflow plane (flume_inflow_df_object%save_state): copied when
                                  present

The split rule (`--rule`) reproduces `tree_object%mpi_redistribute` (adam_tree_object.F90) on the codes sorted by Morton
position (each code compared through its finest descendant, `is_lower_by_pos`):

    units     (default) the units are the complete sibling sets (the `ratio` children of one parent, consecutive in the
              Morton order) and the single blocks; `contiguous_min_max_split` cuts them into M contiguous runs
              minimising the busiest rank (#123). This is the layout of every runtime regrid, and the only one a run
              that derefines may restart from (`tree%derefine` reads a sibling set on one rank): use it for runtime
              AMR (`[amr] frequency > 0`), and wherever in doubt.
    siblings  the sorted codes cut into M contiguous chunks of floor or ceiling of n/M blocks, sibling sets split
              (`split_siblings = .true.`): the layout of a static grid's cold start (`[amr] frequency = 0` with
              `balance_static`, the default; adam_object%rebalance). Only for static grids.

On each rank the blocks are in Morton order with block indexes 1, 2, ... (as `mpi_redistribute` numbers them), so a
repartition to the rank count and rule a dump was written with reproduces its files byte for byte (`--info` says
which rule reproduces a dump's layout). A restart never rebalances (adam_realm_object.F90), so the restarted run keeps
the layout written here until its first regrid. A dump with particles (PIC, nv_pic or np > 0) is refused.

Multi-realm forests: every realm has its own basename; pass them all, each is repartitioned on its own.

Usage:
    repartition.py <basename> [<basename> ...] --to M --out DIR [--rule units|siblings] [--ratio R] [--force]
    repartition.py <basename> [<basename> ...] --info [--ratio R]

`--out` must not be the dump's own directory; with `--force` the restart files of each basename already in DIR are
removed first (without it they stop the tool). `--ratio` (2, 4 or 8: binary tree, quadtree, octree) is inferred from
the codes and levels of the blocks when omitted, and must be given when they cannot tell (e.g. a single level-1 set).
Runs under `python3 -I` (no third-party package: medooza's python3 has none).
"""

from __future__ import annotations

import argparse
import re
import shutil
import struct
import sys
from dataclasses import dataclass
from pathlib import Path

COPIED = (".time", ".reference", ".sinks", ".dff")  # global files copied verbatim (.time is required)
FBD_HEAD = struct.Struct("=7i")  # nv, ni, nj, nk, ngc, nv_pic, np
WMF_HEAD = struct.Struct("=5i")  # blocks, ni, nj, nk, fresh
TND_REC = struct.Struct("=qiq")  # code, rank, block index
I4 = struct.Struct("=i")
I8 = struct.Struct("=q")
BLOCK_HEAD = struct.Struct("=q4i")  # code, i, j, k, level


class RepartitionError(Exception):
    """A restart that cannot be repartitioned (missing, inconsistent or unsupported files)."""


@dataclass
class Block:
    """One block of the dump: its code and level, and where its data sit in the source files."""

    code: int
    level: int
    fbd: Path  # source .fbd
    offset: int  # offset of the block record (code) in the source .fbd
    wmf: Path | None = None  # source .wmf
    wmf_offset: int = 0  # offset of the block record (code) in the source .wmf


@dataclass
class Dump:
    """The restart of one realm as read from disk."""

    basename: Path
    ranks: int
    header: tuple[int, ...]  # the .fbd header
    record: int  # bytes of a block record in the .fbd (code, coordinates, q)
    blocks: list[Block]  # in the order of the source files (rank, then file position)
    layout: list[tuple[int, int, int]]  # the .tnd records as read (code, rank, block index)
    wmf_record: int  # bytes of a block record in the .wmf (code, averages); 0 without .wmf
    wmf_head: tuple[int, ...] | None  # (ni, nj, nk, fresh) of the .wmf


# Morton arithmetic of tree_object (adam_tree_object.F90): codes >= 0, the root (the ancestor of all) is -1


def level_of(code: int, ratio: int) -> int:
    """Return the refinement level of a code (`tree_object%level`)."""
    if code <= -1:
        return 0
    level, c = 1, code
    while c >= ratio:
        c = (c - ratio) // ratio
        if c >= 0:
            level += 1
    return level


def parent_of(code: int, ratio: int) -> int:
    """Return the parent of a code, -2 for the level-1 codes and the root (`tree_object%parent`)."""
    return (code - ratio) // ratio if code > ratio - 1 else -2


def finest_of(code: int, level: int, ratio: int, depth: int) -> int:
    """Return the last descendant of a code at a depth (`tree_object%finest_at_level`): its Morton position key."""
    for _ in range(level, depth):
        code = ratio * code + ratio + (ratio - 1)
    return code


def morton_sorted(codes: list[int], ratio: int) -> list[int]:
    """Return the codes sorted by Morton position (`is_lower_by_pos`, through the finest descendants).

    Leaves never overlap, so the order does not depend on the depth of comparison once it reaches the deepest code;
    ADAM compares at `max_level`, which is that deep or deeper.
    """
    levels = {c: level_of(c, ratio) for c in codes}
    depth = max(levels.values(), default=0)
    return sorted(codes, key=lambda c: finest_of(c, levels[c], ratio, depth))


def contiguous_min_max_split(weights: list[int], parts: int) -> list[int]:
    """Return the part of each unit, a contiguous split minimising the largest part (port of the Fortran function).

    Bisection on the capacity over [max(weights), sum(weights)] for the smallest one the greedy fill fits in `parts`
    parts; the fill also closes a part when the units left equal the parts left, so no part is empty when there are at
    least as many units as parts.
    """
    n = len(weights)
    part = [0] * n
    if n == 0 or parts <= 1:
        return part

    def parts_needed(capacity: int) -> int:
        needed, load = 1, 0
        for w in weights:
            if load + w > capacity:
                needed += 1
                load = 0
            load += w
        return needed

    lo, hi = max(weights), sum(weights)
    while lo < hi:
        cap = lo + (hi - lo) // 2
        if parts_needed(cap) <= parts:
            hi = cap
        else:
            lo = cap + 1
    cap, p, load = lo, 0, 0
    for u, w in enumerate(weights):
        if load > 0 and p < parts - 1 and (load + w > cap or n - u <= parts - 1 - p):
            p += 1
            load = 0
        part[u] = p
        load += w
    return part


def split(codes_sorted: list[int], ranks: int, ratio: int, rule: str) -> list[int]:
    """Return the rank of each code (in Morton order) under a rule of `tree_object%mpi_redistribute`."""
    n = len(codes_sorted)
    if rule == "siblings":
        return [(c * ranks) // n for c in range(n)]
    sizes = []  # the units: a complete sibling set (ratio codes), else one code
    c = 0
    while c < n:
        size = 1
        first = codes_sorted[c]
        if (ratio >= 2 and c + ratio <= n and first >= 0 and first % ratio == 0
                and all(parent_of(codes_sorted[c + k], ratio) == parent_of(first, ratio) for k in range(1, ratio))):
            size = ratio
        sizes.append(size)
        c += size
    owner = []
    for size, part in zip(sizes, contiguous_min_max_split(sizes, ranks), strict=True):
        owner += [part] * size
    return owner


def infer_ratio(blocks: list[Block]) -> int | None:
    """Return the only ratio (2, 4, 8) whose levels of the codes equal the stored block levels, or None."""
    fits = [r for r in (2, 4, 8) if all(level_of(b.code, r) == b.level for b in blocks)]
    return fits[0] if len(fits) == 1 else None


# reading


def rank_files(basename: Path, suffix: str) -> dict[int, Path]:
    """Return the per-rank files `<basename>-proc<rank>.<suffix>` by rank."""
    pattern = re.compile(re.escape(basename.name) + r"-proc(\d{6})" + re.escape(suffix) + "$")
    out = {}
    for f in basename.parent.glob(f"{basename.name}-proc*{suffix}"):
        m = pattern.match(f.name)
        if m:
            out[int(m.group(1))] = f
    return out


def read_dump(basename: Path) -> Dump:
    """Read the layout of the restart of one realm (headers, block codes and offsets; the data stay on disk)."""
    fbds = rank_files(basename, ".fbd")
    if not fbds:
        raise RepartitionError(f"{basename}: no {basename.name}-proc*.fbd file")
    ranks = len(fbds)
    if sorted(fbds) != list(range(ranks)):
        raise RepartitionError(f"{basename}: the .fbd files are not those of ranks 0..{ranks - 1}: {sorted(fbds)}")
    for suffix in (".time", ".tnd"):
        if not Path(f"{basename}{suffix}").is_file():
            raise RepartitionError(f"{basename}: {basename.name}{suffix} is missing")
    header = None
    blocks: list[Block] = []
    record = 0
    for rank in range(ranks):
        with fbds[rank].open("rb") as f:
            head = FBD_HEAD.unpack(f.read(FBD_HEAD.size))
            if header is None:
                header = head
                nv, ni, nj, nk, ngc, nv_pic, np_ = head
                if nv_pic != 0 or np_ != 0:
                    raise RepartitionError(f"{fbds[rank]}: a dump with particles (nv_pic {nv_pic}, np {np_}) is not "
                                           "supported (PIC restarts are PRISM's, not FLUME's)")
                record = BLOCK_HEAD.size + 8 * nv * (ni + 2 * ngc) * (nj + 2 * ngc) * (nk + 2 * ngc)
            elif head != header:
                raise RepartitionError(f"{fbds[rank]}: header {head} differs from rank 0's {header}")
            (nb,) = I4.unpack(f.read(I4.size))
            size = fbds[rank].stat().st_size
            if size != FBD_HEAD.size + I4.size + nb * record:
                raise RepartitionError(f"{fbds[rank]}: {size} bytes, {nb} blocks of {record} bytes expected "
                                       f"{FBD_HEAD.size + I4.size + nb * record}")
            offset = FBD_HEAD.size + I4.size
            for _ in range(nb):
                f.seek(offset)
                code, _i, _j, _k, level = BLOCK_HEAD.unpack(f.read(BLOCK_HEAD.size))
                blocks.append(Block(code=code, level=level, fbd=fbds[rank], offset=offset))
                offset += record
    assert header is not None
    layout = read_tnd(Path(f"{basename}.tnd"))
    check_layout(basename, ranks, blocks, layout)
    wmf_record, wmf_head = read_wmf(basename, ranks, blocks, header)
    return Dump(basename=basename, ranks=ranks, header=header, record=record, blocks=blocks, layout=layout,
                wmf_record=wmf_record, wmf_head=wmf_head)


def read_tnd(path: Path) -> list[tuple[int, int, int]]:
    """Return the (code, rank, block index) records of a .tnd (trailing bytes of an older, longer file ignored)."""
    data = path.read_bytes()
    if len(data) < I8.size:
        raise RepartitionError(f"{path}: empty")
    (n,) = I8.unpack_from(data, 0)
    if len(data) < I8.size + n * TND_REC.size:
        raise RepartitionError(f"{path}: {n} records announced, {(len(data) - I8.size) // TND_REC.size} present")
    return [TND_REC.unpack_from(data, I8.size + k * TND_REC.size) for k in range(n)]


def check_layout(basename: Path, ranks: int, blocks: list[Block], layout: list[tuple[int, int, int]]) -> None:
    """Refuse a dump whose .tnd does not describe its .fbd files (owner rank and position of every block)."""
    where = {}
    position: dict[Path, int] = {}
    for b in blocks:
        position[b.fbd] = position.get(b.fbd, 0) + 1
        rank = int(b.fbd.name[-10:-4])
        if b.code in where:
            raise RepartitionError(f"{basename}: block code {b.code} stored twice")
        where[b.code] = (rank, position[b.fbd])
    tnd = {code: (rank, index) for code, rank, index in layout}
    if len(tnd) != len(layout):
        raise RepartitionError(f"{basename}.tnd: a code is listed twice")
    if set(tnd) != set(where):
        raise RepartitionError(f"{basename}: the .tnd lists {len(tnd)} codes, the .fbd files hold {len(where)}, and "
                               "the two sets differ")
    for code, (rank, index) in tnd.items():
        if rank >= ranks:
            raise RepartitionError(f"{basename}.tnd: block {code} owned by rank {rank}, only {ranks} .fbd files")
        if where[code] != (rank, index):
            raise RepartitionError(f"{basename}.tnd: block {code} at rank {rank} index {index}, the .fbd files hold it "
                                   f"at rank {where[code][0]} index {where[code][1]}")


def read_wmf(basename: Path, ranks: int, blocks: list[Block],
             header: tuple[int, ...]) -> tuple[int, tuple[int, ...] | None]:
    """Attach the wall-model filter records to the blocks; return the record size and the (ni, nj, nk, fresh)."""
    wmfs = rank_files(basename, ".wmf")
    if not wmfs:
        return 0, None
    if sorted(wmfs) != list(range(ranks)):
        raise RepartitionError(f"{basename}: .wmf files for ranks {sorted(wmfs)}, the .fbd for 0..{ranks - 1}")
    by_code = {b.code: b for b in blocks}
    record = 0
    head = None
    for rank in range(ranks):
        with wmfs[rank].open("rb") as f:
            nb, ni, nj, nk, fresh = WMF_HEAD.unpack(f.read(WMF_HEAD.size))
            if (ni, nj, nk) != tuple(header[1:4]):
                raise RepartitionError(f"{wmfs[rank]}: blocks of {ni}x{nj}x{nk} cells, the .fbd {header[1:4]}")
            if head is None:
                head = (ni, nj, nk, fresh)
            elif head != (ni, nj, nk, fresh):
                raise RepartitionError(f"{wmfs[rank]}: header (ni, nj, nk, fresh) {(ni, nj, nk, fresh)} differs from "
                                       f"rank 0's {head}")
            size = wmfs[rank].stat().st_size
            if nb:
                if (size - WMF_HEAD.size) % nb:
                    raise RepartitionError(f"{wmfs[rank]}: {size} bytes do not hold {nb} equal block records")
                this = (size - WMF_HEAD.size) // nb
                if record and this != record:
                    raise RepartitionError(f"{wmfs[rank]}: block records of {this} bytes, other ranks {record}")
                record = this
            elif size != WMF_HEAD.size:
                raise RepartitionError(f"{wmfs[rank]}: no block and {size} bytes")
            for k in range(nb):
                offset = WMF_HEAD.size + k * record
                f.seek(offset)
                (code,) = I8.unpack(f.read(I8.size))
                blk = by_code.get(code)
                if blk is None or int(blk.fbd.name[-10:-4]) != rank:
                    raise RepartitionError(f"{wmfs[rank]}: block {code} is not one of rank {rank}'s .fbd blocks")
                blk.wmf, blk.wmf_offset = wmfs[rank], offset
    missing = [b.code for b in blocks if b.wmf is None]
    if missing:
        raise RepartitionError(f"{basename}: {len(missing)} blocks have no .wmf record (first {missing[0]})")
    return record, head


# writing


def copy_record(dst, src: Path, offset: int, size: int, cache: dict[Path, object]) -> None:
    """Append `size` bytes of `src` from `offset` to the open file `dst` (source files kept open in `cache`)."""
    f = cache.get(src)
    if f is None:
        f = cache[src] = src.open("rb")
    f.seek(offset)
    left = size
    while left:
        chunk = f.read(min(left, 1 << 26))
        if not chunk:
            raise RepartitionError(f"{src}: truncated at offset {offset}")
        dst.write(chunk)
        left -= len(chunk)


def restart_files(basename: Path) -> list[Path]:
    """Return the restart files of a basename present in its directory."""
    out = [Path(f"{basename}{s}") for s in (*COPIED, ".tnd")]
    out = [p for p in out if p.is_file()]
    for suffix in (".fbd", ".wmf"):
        out += list(rank_files(basename, suffix).values())
    return out


def plan(dump: Dump, ratio: int, ranks: int, rule: str) -> tuple[list[int], list[int]]:
    """Return the codes in Morton order and the new owner rank of each."""
    codes = morton_sorted([b.code for b in dump.blocks], ratio)
    return codes, split(codes, ranks, ratio, rule)


def write_dump(dump: Dump, out: Path, ratio: int, ranks: int, rule: str) -> list[int]:
    """Write the restart of `dump` for `ranks` ranks under `out`; return the blocks of each new rank."""
    codes, owner = plan(dump, ratio, ranks, rule)
    by_code = {b.code: b for b in dump.blocks}
    target = out / dump.basename.name
    per_rank: list[list[int]] = [[] for _ in range(ranks)]
    for code, rank in zip(codes, owner, strict=True):
        per_rank[rank].append(code)
    cache: dict[Path, object] = {}
    try:
        for rank, mine in enumerate(per_rank):
            with Path(f"{target}-proc{rank:06d}.fbd").open("wb") as f:
                f.write(FBD_HEAD.pack(*dump.header))
                f.write(I4.pack(len(mine)))
                for code in mine:
                    copy_record(f, by_code[code].fbd, by_code[code].offset, dump.record, cache)
            if dump.wmf_head is not None:
                with Path(f"{target}-proc{rank:06d}.wmf").open("wb") as f:
                    f.write(WMF_HEAD.pack(len(mine), *dump.wmf_head))
                    for code in mine:
                        blk = by_code[code]
                        assert blk.wmf is not None
                        copy_record(f, blk.wmf, blk.wmf_offset, dump.wmf_record, cache)
    finally:
        for f in cache.values():
            f.close()
    index = {code: k + 1 for mine in per_rank for k, code in enumerate(mine)}
    with Path(f"{target}.tnd").open("wb") as f:
        f.write(I8.pack(len(codes)))
        for code, rank in zip(codes, owner, strict=True):
            f.write(TND_REC.pack(code, rank, index[code]))
    for suffix in COPIED:
        src = Path(f"{dump.basename}{suffix}")
        if src.is_file():
            shutil.copyfile(src, f"{target}{suffix}")
    return [len(mine) for mine in per_rank]


def matching_rules(dump: Dump, ratio: int) -> list[str]:
    """Return the rules whose layout on the dump's own rank count is the dump's (owners and block indexes)."""
    written = {code: (rank, index) for code, rank, index in dump.layout}
    out = []
    for rule in ("units", "siblings"):
        codes, owner = plan(dump, ratio, dump.ranks, rule)
        index, seen = {}, [0] * dump.ranks
        for code, rank in zip(codes, owner, strict=True):
            seen[rank] += 1
            index[code] = (rank, seen[rank])
        if index == written:
            out.append(rule)
    return out


def main() -> int:
    """Parse the command line, repartition every basename (or report on it with --info)."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("basenames", nargs="+", type=Path, help="restart basenames (one per realm), paths allowed")
    parser.add_argument("--to", type=int, help="the rank count of the restart to write")
    parser.add_argument("--out", type=Path, help="the directory to write into (not the dump's own)")
    parser.add_argument("--rule", choices=("units", "siblings"), default="units",
                        help="units: sibling sets whole, the runtime-regrid layout (default); siblings: even chunks, "
                             "the static cold-start layout")
    parser.add_argument("--ratio", type=int, choices=(2, 4, 8), help="tree ratio (default: inferred)")
    parser.add_argument("--force", action="store_true", help="remove the basename's restart files in --out first")
    parser.add_argument("--info", action="store_true", help="report the dump's layout, write nothing")
    args = parser.parse_args()
    if not args.info and (args.to is None or args.out is None):
        parser.error("--to and --out are required (or --info)")
    if args.to is not None and args.to < 1:
        parser.error("--to must be at least 1")
    try:
        for basename in args.basenames:
            dump = read_dump(basename)
            ratio = args.ratio or infer_ratio(dump.blocks)
            if args.ratio and any(level_of(b.code, args.ratio) != b.level for b in dump.blocks):
                raise RepartitionError(f"{basename}: --ratio {args.ratio} does not give the stored block levels "
                                       "(the input's [amr] ratio is another)")
            if ratio is None:
                raise RepartitionError(f"{basename}: the tree ratio cannot be inferred from the block codes and "
                                       "levels; pass --ratio (2, 4 or 8, the input's [amr] ratio)")
            nb = [0] * dump.ranks
            for b in dump.blocks:
                nb[int(b.fbd.name[-10:-4])] += 1
            extra = [s for s in (".reference", ".sinks", ".dff") if Path(f"{basename}{s}").is_file()]
            extra += [".wmf"] if dump.wmf_head is not None else []
            step, time = struct.unpack("=id", Path(f"{basename}.time").read_bytes()[:12])
            print(f"{basename}: step {step}, time {time!r}, {len(dump.blocks)} blocks on {dump.ranks} ranks {nb}, "
                  f"ratio {ratio}, also {' '.join(extra) or 'nothing'}")
            if args.info:
                rules = matching_rules(dump, ratio)
                print(f"   layout reproduced by rule: {' and '.join(rules) if rules else 'none (not a fresh split)'}")
                continue
            out = args.out.resolve()
            if out == basename.parent.resolve():
                raise RepartitionError(f"{basename}: --out is the dump's own directory; write elsewhere")
            out.mkdir(parents=True, exist_ok=True)
            stale = restart_files(out / basename.name)
            if stale and not args.force:
                raise RepartitionError(f"{out}: holds restart files of {basename.name} ({stale[0].name}, ...); "
                                       "--force removes them")
            for p in stale:
                p.unlink()
            counts = write_dump(dump, out, ratio, args.to, args.rule)
            print(f"   -> {out / basename.name}: {args.to} ranks, rule {args.rule}, blocks per rank {counts}")
    except RepartitionError as err:
        print(f"repartition: {err}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
