# AMR ghost-fill benchmark (2:1 neighbourhoods)

**Question** (2026-11-25): the uniform benchmark (`../README.md`) showed ADAM's local ghost fill moves its data at
10-16 % of a GPU's copy rate, the per-cell index map being the cost. On the frozen AMR sphere of #133 the same kernel
costs 5.6 times more per launch than on a uniform grid of the same row count (`speedup/README.md`, "Where the AMR
overhead goes": 14.6 against 2.6 ms on an A30), and how that split among the interpolation arithmetic, the branch
divergence and the donor reads was not measured. What does a 2:1 row cost, why, and what would a redesign gain?

`test_amr_ghost_bench.F90` (FoBiS mode `test-amr-ghost-bench-fnl-nvf`, linked to ADAM) builds a real 2:1 grid with
ADAM's own tree, grid and maps objects: the unit cube, an octree of blocks of b^3 cells (5 variables, `ngc` 3), refined
uniformly to a base level, then one more level on the blocks inside the box [0,s]^3: 2:1 faces, edges and corners
around the box. The FLUME density guard of the seam fill (#101) is on, the field smooth and positive, the seam regimes
those of `[amr] seam_ghost_fill` and `seam_restriction`, every high-order one included.

ADAM describes the ghost fill with a table, `maps%local_map_ghost(f, 1:13)`, one descriptor per (receiver block, donor
block, region): ranges, index shifts and the refinement relation (same level, finer donor: restriction, coarser
donor: prolongation). It expands the table into `maps%local_map_ghost_cell`, one row per ghost cell with its seam
metadata, deduplicated so that the last descriptor covering a ghost wins (#22 F1); only the rows reach the GPU, where
`update_ghost_local_gpu_dev` runs every kind of row in one launch, a branch per row.

| Variant | Kernel | Isolates |
|---|---|---|
| V0 | `update_ghost_local_gpu_dev` on ADAM's rows | today |
| V0s | the same kernel on the rows split by kind (copy / restriction / prolongation), one launch each | branch divergence; the per-kind costs |
| D1 | one kernel per kind over the descriptors: work items decoded from the ranges, the seam metadata computed on the device by ADAM's pure routines, the values by ADAM's operators (`interp_seam_ghost_gpu_dev`, `restrict_seam_high_order_gpu_dev`, the mean) | the per-cell map, on AMR |
| D2 | D1 batched in shape classes (same kind, ranges, shifts), the member blocks innermost | D1 with the uniform benchmark's coalescing |
| V0ng | V0 with the density guard off | the guard's cost; the reference of D3 |
| D3 | prolongation per coarse anchor (tricubic, guard off): the union of its 8 children's 4^3 footprints gathered once, each child interpolated by ADAM's `seam_interpolate_tricubic` on its slice | the reads the children of an anchor share |

Overlapping descriptors: an owner mask (one int32 per stored cell of every block, from the descriptor order) gives each
ghost its winning descriptor; D1, D2 and D3 write a ghost only from its owner, ADAM's last-writer rule. **Every variant
is checked bit for bit after its timing: V0s, D1 and D2 against V0, D3 against V0ng's prolongation; every run, every
regime pair, both GPUs: identical.**

## Runs (2026-11-25)

`run.sh <varset> <out>` builds the mode and runs the matrix on one GPU (median of 30 fills, 20 at base level 3):

1. the **production-like grid**: blocks of 24^3, base level 2, box 0.75: 253 blocks (levels 2-3), 3.5 M cells, 2.87 M
   rows (88.1 % copy, 2.3 % restriction, 9.7 % prolongation), the row count per GPU of the frozen sphere 2c on 4 A30
   (254 blocks of 24^3, 2.86 M rows); every prolongation (tricubic, quintic, quintic-mp, weno6) x restriction (mean,
   lagrange6, weno6) pair;
2. the refined fraction: box 0.25 and 0.5 (4.2 % and 9.7 % prolongation rows);
3. the block edge: 16 and 32 at box 0.5;
4. a larger load: base level 3, box 0.5 (960 blocks, 13 M cells; `fields` 40 so the capacity holds the block indices).

Jobs: ESA SpaceHPC, one NVIDIA H100 (nvfortran 26.3, `spacehpc_nvf`): 184161 (V0, V0s, D1, D2), 184174 (with V0ng,
D3), 184175 (the larger load); medooza, one NVIDIA A30 (nvfortran 26.9, `iac_nvf`): 1369803, 1369806, 1369807 (the larger load). Results in
`results/`, figures by `plot.py` (in the docs page, below).

**Repeatability.** Every case ran twice per GPU (the two passes). On the H100 every time reproduces within 0.1 %. On
the A30, V0s reproduces within 4.2 % and V0 within 0.3 % in 15 of 20 cases; in 5 cases one pass of V0 ran 3-34 %
slower (tricubic / mean 7.83 and 10.45 ms, the first case of its job; 32^3 tricubic 5.86 and 7.86 ms): the
figures take the lower V0 of the two passes, the divergence figures below quote only cases stable in both.

## Results (H100, the reference)

**The fill on the production-like grid** (ms; V0s's launches by kind):

| Prolongation / restriction | V0 today | copy (88 % of rows) | restriction (2.3 %) | prolongation (9.7 %) | V0ng, guard off (its prolongation) |
|---|---|---|---|---|---|
| tricubic / mean (the defaults) | 2.54 | 0.92 | 0.10 | 1.38 | 1.15 (0.23) |
| tricubic / lagrange6 | 2.67 | 0.92 | 0.32 | 1.37 | 1.29 (0.24) |
| quintic / lagrange6 | 3.94 | 0.92 | 0.32 | 2.63 | 1.89 (0.88) |
| quintic-mp / mean | 5.18 | 0.92 | 0.10 | 4.12 | 2.42 (1.52) |
| weno6 / mean | 7.53 | 0.92 | 0.10 | 6.38 | 3.77 (2.84) |
| weno6 / weno6 | 8.57 | 0.92 | 1.02 | 6.38 | 4.68 (2.84) |

**The cost of one row** (all 5 variables, guard on): copy **0.36 ns**; restriction mean 1.5, lagrange6 4.9, weno6
15.5 ns; prolongation tricubic **5.0**, quintic 9.5, quintic-mp 14.9, weno6 **23.0 ns**. A prolongation row costs 14 to
64 copy rows.

- **The 2:1 cost is the prolongation arithmetic.** Under 10 % of the rows, the prolongation takes 54 % (tricubic) to
  85 % (weno6) of the fill; a weno6 restriction costs 1.0 ms, as much as all the copies. The copies cost what they cost on a
  uniform grid.
- **The density guard is most of it.** Without the guard the tricubic prolongation takes 0.23 ms instead of 1.38 (the
  guard is 83 % of it), the quintic 0.88 instead of 2.63, the weno6 2.84 instead of 6.38; the whole fill 45-59 % of its
  guarded time. The guard re-interpolates the density for each of the other four variables and re-gathers its
  footprint, per fine ghost: the same density work 5 times.
- **Branch divergence is minor on the H100**: splitting the rows by kind (V0s) changes the fill by -10 % to +7 %. On the
  A30 it gains up to 38 % (tricubic / lagrange6 11.5 -> 7.9 ms; 16^3 quintic / lagrange6 4.44 -> 2.74 ms).
- **The descriptor kernels do not transfer the uniform gain to AMR.** D2 makes the copies 18 % faster (0.75 against
  0.92 ms), far from the uniform benchmark's 3.8x: AMR descriptors are many, small and irregular (5,310 here, 154 shape
  classes), and the owner mask and the search are paid per cell. On the prolongation D1 and D2 are 1.2-2.1x **slower**
  than V0: they recompute the seam metadata per cell that the rows carry precomputed.
- **Sharing a footprint per anchor in one thread loses** (D3, tricubic, guard off): 0.58 ms against 0.23 for V0ng's
  prolongation, 2.5x slower. Eight children computed serially by one thread, with a 2.7 KB private footprint, lose more
  parallelism than the shared reads save; sharing has to be cooperative (threads of a block over shared memory), not
  serial.
- **The same at 4 times the load** (base level 3, 960 blocks, 13.3 M cells, 11.3 M rows of which 4.4 % prolongation;
  job 184175): tricubic / mean V0 6.18 ms (copy 3.33, restriction 0.15, prolongation 2.60), copy 0.31 ns per row,
  tricubic prolongation 5.2 ns, weno6 22.8 ns; guard off, the tricubic prolongation 0.43 ms (the guard 84 %); D3 0.94
  against 0.43, 2.2x slower. The A30 at the same load (job 1369807): V0 20.2 ms (copy 9.55, restriction 0.55,
  prolongation 9.36), guard off 1.20 (the guard 87 %), D3 2.60, 2.2x slower; every variant bitwise.
- The A30 shows the same split, 2.6-4.3x slower per row (row costs and second-pass numbers in `results/`).

**Production check.** The copy rows reproduce FLUME: on the A30 the copy launch costs 2.44 ms for 2.53 M rows, FLUME's
uniform 2a launch 2.6 ms for a similar row count. FLUME's AMR launch on 2c costs 14.6 ms against 7.8 ms here at the same
total row count: the difference is the mix. With the per-row costs above, 2c's launch implies ~26 % prolongation rows,
against 9.7 % in the box: plausible for a refinement that hugs a sphere over three levels (more seam surface per volume
than a box), but 2c's row mix is not measured.

## What it says for the redesign

1. **The density guard is the first target of the 2:1 fill**, before any data-structure change: compute the density,
   its bound and the guard decision once per fine ghost and reuse them for the other variables (bitwise: the same
   operations, done once). On this grid it is worth up to half of the fill.
2. **The prolongation operators are compute, not data movement**: 5-23 ns per row, 14-64 copy rows. They belong to the
   first track (arithmetic) more than to the second; a block-cooperative kernel (the children of a set of anchors
   sharing their footprints through shared memory) is the path the per-thread D3 rules out.
3. **Group the rows by kind**: free with today's maps, bitwise, up to 38 % on the A30.
4. **Descriptors need structure, not just compression**: on AMR the gain of the uniform slabs comes back only if the
   descriptors are few and regular (non-overlapping by construction, no owner mask, member blocks contiguous) and the
   seam metadata stays precomputed per class rather than per cell.
