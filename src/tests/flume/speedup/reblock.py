#!/usr/bin/env python3
"""Re-block a FLUME restart one tree level deeper: every block becomes its children (stdlib only, issue #133 leg 3).

Why: `field%load_blocks` reads only a dump whose blocks have the input's `[grid] ni, nj, nk`, so a layout of smaller
blocks needs the dump rewritten. Leg 3 of the speedup matrix restarts from the developed T3 channel dump (L2, blocks
of 64 x 48 x 32 cells), and its Kolmogorov marker can coarsen the channel core only with blocks of 32 x 24 x 16 at L3
(an L1 parent of the L2 layout spans wall to centreline). This tool writes that layout from the dump, offline.

Every block is replaced by its `ratio` children at the next level (octree 2 x 2 x 2; quadtree 2 x 2 in x, y, z never
split; binary tree 2 in x): the same cells, the same spacing h, no interpolation. A child takes its octant of the
parent's cells, ghosts included: the child's ghost cells are parent interior cells at the faces between siblings and
the parent's own ghosts at the parent's faces, so every value of the child's record is one of the parent's record.
The restart refreshes every ghost before the first step anyway (`post_initial_state_forest`).

The restart files of one realm (formats in repartition.py's docstring; the writers there):

    <basename>.time               copied
    <basename>.tnd                REWRITTEN: the children, sorted by Morton position, with their new owners and indexes
    <basename>-proc<r:06d>.fbd    REWRITTEN: header (nv, ni/2, nj/2, nk/2 on the split axes, ngc, 0, 0); per child
                                  int64 code, int32 i, j, k, level (the parent's doubled plus the octant bits, level
                                  + 1), float64 q(nv, 1-ngc:ni/2+ngc, ...) cut from the parent's q
    <basename>.reference, .sinks, .dff
                                  copied (references, the replicated sink list, the inflow plane: no block in them)
    <basename>-proc<r:06d>.wmf    REFUSED: the wall-model filter averages live on the block faces, and a child's face
                                  between siblings has no parent face to take them from (delete the files to restart
                                  the filter from the next sample, not bitwise)
    <basename>-proc<r:06d>.h5, <basename>.xdmf
                                  NOT WRITTEN: the restart's visualisation mirror (save_xh5f); no restart reads it

The Morton arithmetic (adam_tree_object.F90, MORTIF): the children of a code c are `ratio * c + ratio + o`,
o = 0 .. ratio - 1, and the bits of o are the child's offsets, bit 0 along x (i), bit 1 along y (j), bit 2 along z (k)
(`demorton3D`/`demorton2D` of the local child index); a block's coordinates (i, j, k) are its integer position at its
level, so a child is at (2 i + bit 0, 2 j + bit 1, 2 k + bit 2) on the split axes. The parents' stored coordinates are
checked against their codes before anything is written.

The new layout is a fresh split of the children over `--to` ranks (default: the dump's rank count) by the rule of
repartition.py (`--rule units`, the default, keeps every sibling set on one rank: the layout of a regrid, the only one
a run that derefines may restart from). The input that reads the result needs `[grid] ni, nj, nk` halved on the split
axes (ngc unchanged), `[amr] max_level` at least the deepest level + 1 and `[amr] iu_ref_levels` + 1 where it names the
uniform level (make_channel.py `--levels L + 1` for a dump written with `--levels L`). A dump with particles is refused.

Usage:
    reblock.py <src> [<src> ...] --out DIR [--to M] [--rule units|siblings] [--ratio R] [--force]

A <src> is a restart basename (`dir/channel-restart`) or a directory, which stands for every basename with a `.tnd` in
it. `--out` must not be a source directory; with `--force` the restart files of each basename already in DIR are removed
first. Runs under `python3 -I` (no third-party package: medooza's python3 has none).
"""

from __future__ import annotations

import argparse
import importlib.util
import shutil
import struct
import sys
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Any


def _load_repartition() -> ModuleType:
    """Import repartition.py from this script's directory (`python3 -I` keeps that directory off sys.path)."""
    path = Path(__file__).resolve().parent / "repartition.py"
    spec = importlib.util.spec_from_file_location("repartition", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules["repartition"] = module
    spec.loader.exec_module(module)
    return module


rp = _load_repartition()

R8 = 8  # bytes of a float64


class ReblockError(Exception):
    """A restart that cannot be re-blocked (missing, inconsistent or unsupported files)."""


@dataclass(frozen=True)
class Geometry:
    """The cell counts of the parents and of the children, and the axes the tree splits."""

    nv: int
    ngc: int
    parent: tuple[int, int, int]  # ni, nj, nk of a parent block
    child: tuple[int, int, int]  # ni, nj, nk of a child block
    split: tuple[bool, bool, bool]  # the axes the tree splits (x; x, y; x, y, z)


def split_axes(ratio: int) -> tuple[bool, bool, bool]:
    """Return the axes a tree of `ratio` splits (grid_object%initialize: binary x, quadtree x and y, octree all)."""
    return {2: (True, False, False), 4: (True, True, False), 8: (True, True, True)}[ratio]


def octant_offsets(o: int) -> tuple[int, int, int]:
    """Return the offsets (along x, y, z) of the local child index `o` (bits 0, 1, 2: demorton of MORTIF)."""
    return o & 1, (o >> 1) & 1, (o >> 2) & 1


def coordinates_of(code: int, ratio: int) -> tuple[int, int, int, int]:
    """Return the block coordinates (i, j, k, level) of a code (`tree_object%morton_to_coordinates`)."""
    level = rp.level_of(code, ratio)
    ijk = [0, 0, 0]
    c, p = code, 0
    while True:
        local = (c + ratio) % ratio if c > 0 else 0  # child_local
        for axis, bit in enumerate(octant_offsets(local)):
            ijk[axis] += bit << p
        if rp.level_of(c, ratio) <= 1:
            break
        c = rp.parent_of(c, ratio)
        p += 1
    return ijk[0], ijk[1], ijk[2], level


def geometry(header: tuple[int, ...], ratio: int) -> Geometry:
    """Return the geometry of the re-blocked dump; refuse block sizes that do not halve into usable children."""
    nv, ni, nj, nk, ngc, _nv_pic, _np = header
    split = split_axes(ratio)
    parent = (ni, nj, nk)
    child = []
    for axis, (n, s) in enumerate(zip(parent, split, strict=True)):
        if not s:
            child.append(n)
            continue
        if n % 2:
            raise ReblockError(f"the blocks have {n} cells along axis {'xyz'[axis]}: odd, no halving")
        if n // 2 < ngc:
            raise ReblockError(f"children of {n // 2} cells along axis {'xyz'[axis]} are thinner than the {ngc} ghost "
                               "layers (the ghost exchange reads ngc interior cells of a neighbour)")
        if n // 2 < 2 * ngc:
            print(f"   warning: children of {n // 2} cells along {'xyz'[axis]}, fewer than 2 ngc = {2 * ngc}: a 2:1 "
                  "face (runtime AMR) restricts them to fewer coarse cells than the ghost layers", file=sys.stderr)
        child.append(n // 2)
    return Geometry(nv=nv, ngc=ngc, parent=parent, child=(child[0], child[1], child[2]), split=split)


def children_records(record: bytes, geo: Geometry, ratio: int) -> list[tuple[int, bytes]]:
    """Return (code, record) of the children of one parent block record, in Morton order (local index 0 .. ratio-1)."""
    code, i, j, k, level = rp.BLOCK_HEAD.unpack_from(record, 0)
    if (i, j, k, level) != coordinates_of(code, ratio):
        raise ReblockError(f"block {code}: stored coordinates {(i, j, k, level)}, its code gives "
                           f"{coordinates_of(code, ratio)}")
    q = memoryview(record)[rp.BLOCK_HEAD.size:]
    g = geo.ngc
    pni, pnj, _pnk = (n + 2 * g for n in geo.parent)  # parent extents with ghosts
    cni, cnj, cnk = (n + 2 * g for n in geo.child)  # child extents with ghosts
    row = R8 * geo.nv * cni  # one child i-run (all variables, contiguous in the parent too)
    out = []
    for o in range(ratio):
        off = octant_offsets(o)
        # child cell c (1-ngc based) is parent cell off * n_child + c: the 0-based parent index of the child's first
        # ghost cell is off * n_child (the child range never leaves the parent's ghosted range)
        start = [off[a] * geo.child[a] if geo.split[a] else 0 for a in range(3)]
        ci, cj, ck = (2 * v + off[a] if geo.split[a] else v for a, v in enumerate((i, j, k)))
        child_code = ratio * code + ratio + o
        if rp.level_of(child_code, ratio) != level + 1:
            raise ReblockError(f"block {code}: child {child_code} is not at level {level + 1}")
        parts = [rp.BLOCK_HEAD.pack(child_code, ci, cj, ck, level + 1)]
        for kk in range(cnk):
            pk = start[2] + kk
            for jj in range(cnj):
                pj = start[1] + jj
                first = R8 * geo.nv * ((pk * pnj + pj) * pni + start[0])
                parts.append(q[first:first + row])
        out.append((child_code, b"".join(parts)))
    return out


def sources(args: list[Path]) -> list[Path]:
    """Return the basenames named by the command line (a directory stands for every `.tnd` basename in it)."""
    out = []
    for src in args:
        if src.is_dir():
            found = sorted(p.with_suffix("") for p in src.glob("*.tnd"))
            if not found:
                raise ReblockError(f"{src}: no restart (no *.tnd file)")
            out += found
        else:
            out.append(src)
    return out


def reblock(dump: Any, geo: Geometry, out: Path, ratio: int, ranks: int, rule: str) -> tuple[list[int], int]:
    """Write the re-blocked restart of `dump` for `ranks` ranks under `out`; return the blocks per rank and the depth."""
    parents = rp.morton_sorted([b.code for b in dump.blocks], ratio)
    children = [ratio * c + ratio + o for c in parents for o in range(ratio)]
    if children != rp.morton_sorted(children, ratio):
        raise ReblockError(f"{dump.basename}: the children of the Morton-sorted parents are not Morton-sorted")
    owner = rp.split(children, ranks, ratio, rule)
    index, seen = {}, [0] * ranks
    for code, rank in zip(children, owner, strict=True):
        seen[rank] += 1
        index[code] = seen[rank]
    owner_of = dict(zip(children, owner, strict=True))
    target = out / dump.basename.name
    by_code = {b.code: b for b in dump.blocks}
    files = [Path(f"{target}-proc{r:06d}.fbd").open("wb") for r in range(ranks)]
    cache: dict[Path, Any] = {}
    try:
        for r, f in enumerate(files):
            f.write(rp.FBD_HEAD.pack(geo.nv, *geo.child, geo.ngc, 0, 0))
            f.write(rp.I4.pack(seen[r]))
        for code in parents:
            blk = by_code[code]
            src = cache.get(blk.fbd)
            if src is None:
                src = cache[blk.fbd] = blk.fbd.open("rb")
            src.seek(blk.offset)
            record = src.read(dump.record)
            if len(record) != dump.record:
                raise ReblockError(f"{blk.fbd}: truncated at offset {blk.offset}")
            for child_code, child_record in children_records(record, geo, ratio):
                files[owner_of[child_code]].write(child_record)
    finally:
        for f in files:
            f.close()
        for f in cache.values():
            f.close()
    with Path(f"{target}.tnd").open("wb") as f:
        f.write(rp.I8.pack(len(children)))
        for code, rank in zip(children, owner, strict=True):
            f.write(rp.TND_REC.pack(code, rank, index[code]))
    for suffix in rp.COPIED:
        src_file = Path(f"{dump.basename}{suffix}")
        if src_file.is_file():
            shutil.copyfile(src_file, f"{target}{suffix}")
    depth = max(rp.level_of(c, ratio) for c in children)
    return seen, depth


def main() -> int:
    """Parse the command line and re-block every basename."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("sources", nargs="+", type=Path, help="restart basenames or directories holding them")
    parser.add_argument("--out", type=Path, required=True, help="the directory to write into (not a source's own)")
    parser.add_argument("--to", type=int, help="the rank count of the restart to write (default: the dump's)")
    parser.add_argument("--rule", choices=("units", "siblings"), default="units",
                        help="units: sibling sets whole, the runtime-regrid layout (default); siblings: even chunks, "
                             "the static cold-start layout")
    parser.add_argument("--ratio", type=int, choices=(2, 4, 8), help="tree ratio (default: inferred)")
    parser.add_argument("--force", action="store_true", help="remove the basename's restart files in --out first")
    args = parser.parse_args()
    if args.to is not None and args.to < 1:
        parser.error("--to must be at least 1")
    try:
        for basename in sources(args.sources):
            try:
                dump = rp.read_dump(basename)
            except rp.RepartitionError as err:
                raise ReblockError(str(err)) from err
            if dump.wmf_head is not None:
                raise ReblockError(f"{basename}: wall-model filter averages (.wmf) cannot be re-blocked (a child face "
                                   "between siblings has no parent face); delete the .wmf files to restart the filter "
                                   "from the next sample")
            ratio = args.ratio or rp.infer_ratio(dump.blocks)
            if args.ratio and any(rp.level_of(b.code, args.ratio) != b.level for b in dump.blocks):
                raise ReblockError(f"{basename}: --ratio {args.ratio} does not give the stored block levels")
            if ratio is None:
                raise ReblockError(f"{basename}: the tree ratio cannot be inferred from the block codes and levels; "
                                   "pass --ratio (2, 4 or 8, the input's [amr] ratio)")
            geo = geometry(dump.header, ratio)  # refuses unusable children before anything is written
            out = args.out.resolve()
            if out == basename.parent.resolve():
                raise ReblockError(f"{basename}: --out is the dump's own directory; write elsewhere")
            out.mkdir(parents=True, exist_ok=True)
            stale = rp.restart_files(out / basename.name)
            if stale and not args.force:
                raise ReblockError(f"{out}: holds restart files of {basename.name} ({stale[0].name}, ...); --force "
                                   "removes them")
            for p in stale:
                p.unlink()
            ranks = args.to or dump.ranks
            levels = sorted({b.level for b in dump.blocks})
            counts, depth = reblock(dump, geo, out, ratio, ranks, args.rule)
            step, time = struct.unpack("=id", Path(f"{basename}.time").read_bytes()[:12])
            print(f"{basename}: step {step}, time {time!r}, {len(dump.blocks)} blocks of {'x'.join(map(str, geo.parent))} "
                  f"cells at L{levels[0]}..L{levels[-1]} on {dump.ranks} ranks, ratio {ratio}")
            print(f"   -> {out / basename.name}: {len(dump.blocks) * ratio} blocks of "
                  f"{'x'.join(map(str, geo.child))} cells at L{levels[0] + 1}..L{depth} on {ranks} ranks, rule "
                  f"{args.rule}, blocks per rank {counts}")
            print(f"   the input: [grid] ni = {geo.child[0]}, nj = {geo.child[1]}, nk = {geo.child[2]}; [amr] "
                  f"max_level >= {depth}, iu_ref_levels one higher")
    except ReblockError as err:
        print(f"reblock: {err}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
