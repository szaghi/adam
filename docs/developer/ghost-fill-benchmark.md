# Ghost fill on the GPU: a microbenchmark

A standalone measurement of how fast ADAM fills the ghost cells of its blocks on a GPU, against what the same GPU can
copy, and against the alternatives an AMR redesign could take. It is the first piece of evidence for that redesign
(2026-11-25): **today's fill moves its data at 10–16 % of the GPU's copy rate, and structured slab copies, with no
per-cell index map, reach 43–57 % and are 3.0–5.5 times faster, bit for bit identical.** On a 2:1 AMR grid the picture
changes ([AMR neighbourhoods](#amr-neighbourhoods-2-1-seams)): **the cost moves to the prolongation rows, under 10 % of
the rows and 54–85 % of the fill, and most of that is the density guard of the seam fill.** The code, the job scripts
and the raw results are in
[`src/tests/perf/ghost-fill/`](https://github.com/szaghi/adam/tree/develop/src/tests/perf/ghost-fill) (uniform) and its
`amr/` (2:1 neighbourhoods).

## Why it was done

The scaling campaigns of [#133](https://github.com/szaghi/adam/issues/133) timed a FLUME step phase by phase on A30 and
H100 nodes ([performance](/applications/flume/performance)). Three findings set the question:

1. **The arithmetic dominates.** The face-flux kernel takes 60–85 % of a step; a cell's cost depends weakly on the block
   size (on one H100, blocks of 32³ cost 8 % more per cell than 64³; AthenaK, with a lighter cell update, loses 1.7× over
   the same range). FLUME's high-order cell update is heavy, by design.
2. **The face kernel is not yet compute-bound.** On the H100 it runs at roughly 3–7 % of the FP64 peak, held back by
   local memory and latency (its private stencil arrays spill; profile in `speedup/profile-couette3d/README.md`). Making
   it faster is the first track of the next optimisation round.
3. **Everything around the ghost cells is data movement.** Without an immersed body, the local fill, the boundary
   conditions and the MPI exchange take up to ~20 % of a step on the H100 (~30 % with the body). With a faster face kernel
   that share grows (Amdahl): at a face kernel twice as fast the ghost work of the body case passes a third of the step.
   So the ghost machinery has to reach its own floor: the second track.

Before redesigning the AMR data structures for the second track, the question is how far today's fill is from that
floor, and why. A whole FLUME step mixes too many effects to answer it; a microbenchmark isolates the fill.

## What it measures

`ghost_fill_bench.F90` (plain OpenACC, no ADAM dependency) builds a uniform periodic brick of blocks of b³ cells with
`ngc` ghost layers and 5 variables, and fills every ghost cell (faces, edges and corners) from its same-level
neighbour: the rows that ADAM's local ghost map holds. Each variant changes one thing:

| Variant | Kernel | Isolates |
|---|---|---|
| **V0** | ADAM's `update_ghost_local_gpu_dev` (`adam_fnl_field_kernels.F90`) replicated: one row per ghost cell, 10 `integer(I8P)` columns of which 9 are read, `collapse(2)` with the variable loop outermost, the branch on the row's kind, layout `q(b,i,j,k,v)` (block index fastest), rows in the order `adam_maps_object` builds them | today |
| V1 | V0 with `integer(int32)` columns | the bytes of the indices |
| V2 | one linear offset for the donor and one for the receiver per row, the variables inside the row | the index traffic and the map read once per variable |
| V3 | **structured slabs**: no per-cell map; one kernel per neighbour direction (26) over (variable, k, j, i, block), the donor block from a per-block neighbour table, the donor cell by a constant shift; layout `q(b,i,j,k,v)` | the per-cell map |
| V4 | V3 on the layout `q(i,j,k,v,b)`, cell index fastest | the layout |
| V5 | V3 **fused into one launch**: the ghost cells of a block numbered direction by direction, the direction and cell decoded from a 26-entry table | the launch count |
| V6 | V5 on the layout `q(i,j,k,v,b)` | the layout, one launch |
| ref | a contiguous copy of as many values | the floor the GPU reaches |

**Method.** Median of 50 fills after a warm-up, each followed by a device wait. Sizes: ~7 M cells per GPU in blocks of
16³, 24³, 32³, 38³ and 48³ (the load of the 7 M-cell weak pairs of #133), ~56 M in 24³, 38³ and 48³, and 38³ with
`ngc` 4 (WENO-7). The bandwidth counts 16 bytes per ghost value (one read, one write); the index reads of the map
variants are not counted, so their figure is what the data alone would need. **Every variant, at every size, on both
GPUs, is compared bit for bit with a fill on the host and is identical**: a copy is exact, so any difference would be a
bug.

**Machines.** ESA SpaceHPC, one NVIDIA H100 (nvfortran 26.3, `-gpu=cc90,cuda12.9`, PBS job 184030), the reference; and
medooza, one NVIDIA A30 (nvfortran 26.9, `-gpu=cc80,cuda12.9`, Slurm job 1369794), for comparison. Each run printed the
GPU's clock and throttle bits before and after: no throttling.

## The benchmark reproduces ADAM's fill

V0 against the local ghost fill timed inside FLUME steps by the phase timer, at the same load and block size (per
call: the step does 3 fills without the body, 12 with it):

| | V0 | FLUME, phase timer |
|---|---|---|
| H100, 56 M cells per GPU, 38³ | 9.09 ms | 9.5 ms |
| H100, 7 M cells per GPU, 38³ | 1.10 ms | 1.12 ms |
| A30, 7 M cells per GPU, 38³ | 3.27 ms | 3.1 ms |

Within 5 %: what the benchmark says about V0 is what ADAM does.

## Results

![The fraction of the copy floor every variant reaches against the block edge, on one H100 and one A30](/developer/ghost-fill-floor.png)

**Today's fill (V0) moves its data at 10–16 % of what the GPU copies**, on both GPUs, at every block size and both loads.
**Structured slabs (V3, V5) reach 43–57 %.** Narrowing the map to 32-bit indices (V1) helps a little; precomputed offsets
with the variables inside the row (V2) are worse than today, because each thread then walks five addresses a whole
variable apart.

![Today's fill time over the best structured variant's, and over the 32-bit map's, against the block edge](/developer/ghost-fill-speedup.png)

**The structured fill is 3.0–5.5 times faster than today's**, more at larger blocks and loads (H100: 3.4× at 16³, 7 M
cells, to 4.5× at 48³, 56 M; A30: up to 5.5×). The 32-bit map gives 1.3–1.4× at every point.

Time of one fill, the reference sizes (all in `src/tests/perf/ghost-fill/README.md`):

| | copy floor | V0 today | V1 I4P | best slabs | speedup |
|---|---|---|---|---|---|
| H100, 7 M cells, 38³ | 0.164 ms | 1.096 ms (15 %) | 0.809 ms | 0.291 ms (56 %, V5) | 3.8× |
| H100, 56 M cells, 38³ | 1.156 ms | 9.094 ms (13 %) | 6.773 ms | 2.187 ms (53 %, V5) | 4.2× |
| A30, 7 M cells, 38³ | 0.394 ms | 3.269 ms (12 %) | 2.312 ms | 0.905 ms (44 %, V5) | 3.6× |
| A30, 56 M cells, 38³ | 3.031 ms | 24.35 ms (12 %) | 19.00 ms | 5.559 ms (55 %, V3) | 4.4× |

![A FLUME step on 4 H100 split into its phases: today, with the local fill at the slab speed, and with that plus a face kernel twice as fast](/developer/ghost-fill-step.png)

**What it means for a FLUME step** (4 H100, the uniform sphere in 38³ blocks, phase timer; the local fill scaled by the
measured slab/V0 ratio; the right-hand bar of each group also halves the face fluxes, a **scenario** for the first track,
not a measurement). Today the slab fill alone would take 2.8 % (7 M cells, no body) to 10 % (56 M cells, with the body)
off a step. With a face kernel twice as fast, the ghost work, at today's boundary-condition and exchange cost, is
13 % (no body) and 34 % (body) of the step at 7 M cells: the ghost machinery would again be a third of the step where
the immersed boundary is present.

## What it proves

- **The cost of today's fill is the per-cell index map, not the data.** Every copied value carries a row of nine 8-byte
  indices; the copies are gathers and scatters through them. The data alone would move 6–10 times faster.
- **Removing the map is worth three times more than shrinking it.** 32-bit indices gain 1.3–1.4×; computing the indices
  from a per-block neighbour table and the slab geometry gains 3.0–5.5×. A map rewrite does not reach the floor; a
  structural change does.
- **ADAM's block-fastest layout costs the ghost machinery nothing.** The cell-fastest layout is never faster (V6 against
  V5, V4 against V3). The layout chosen for the face kernel's coalescing can stay; the first track is free to keep it.
- **The launch strategy depends on the GPU and the load.** One fused launch wins on the H100 nearly everywhere; on the
  A30, 26 per-direction launches win at 56 M cells and at 16³. The design should allow both.
- **The structured fill is about a factor 2 from the copy floor**: thin slabs (3–4 cells) gathered from 26 neighbours.
  That factor is the fill's own floor, not a defect.

## AMR neighbourhoods (2:1 seams)

The second benchmark, `src/tests/perf/ghost-fill/amr/`, measures the fill where the rows are not copies. It is linked
to ADAM: the grid is built by ADAM's own tree, grid and maps objects (the unit cube, an octree of blocks of b³ cells
refined to a base level, then one more level inside the box [0,s]³, so 2:1 faces, edges and corners around the box);
the values come from ADAM's own seam operators, every high-order regime of `[amr] seam_ghost_fill` and
`seam_restriction`; FLUME's density guard of the seam fill ([#101](https://github.com/szaghi/adam/issues/101)) is on.

ADAM describes a fill with a descriptor table (`maps%local_map_ghost`, one entry per receiver block, donor block and
region, with its ranges, shifts and refinement relation), expands it into one row per ghost cell with its seam
metadata (`maps%local_map_ghost_cell`, deduplicated so the last descriptor covering a ghost wins), and runs every row in
one launch, a branch per row (V0). The variants:

| Variant | Kernel | Isolates |
|---|---|---|
| V0 | `update_ghost_local_gpu_dev` on ADAM's rows | today |
| V0s | the same kernel on the rows split by kind (copy, restriction, prolongation), one launch each | branch divergence; the cost of each kind |
| D1 | one kernel per kind over the descriptors, the seam metadata computed on the device by ADAM's routines | the per-cell map, on AMR |
| D2 | D1 batched in shape classes, the member blocks innermost | D1 with the uniform benchmark's coalescing |
| V0ng | V0 with the density guard off | the guard's cost |
| D3 | prolongation per coarse anchor (tricubic, guard off): the 8 children's footprints gathered once, each child by ADAM's `seam_interpolate_tricubic` | the reads the children share |
| V0p, V0png | since [#157](https://github.com/szaghi/adam/issues/157): the kernel of b8641a59 (its guard function copied verbatim into the harness), guard on and off | the before of #157 in the same run |
| V0k, V0kng | since #157: ADAM's kernel on the rows sorted by kind, guard on and off | grouping the rows by kind |

An owner mask (one integer per stored cell) reproduces the last-writer rule for the descriptor kernels. **Every variant,
every regime pair, both GPUs: bit for bit identical to its reference** (V0, and V0ng for D3).

The production-like grid has blocks of 24³ at levels 2–3, 253 blocks, 3.5 M cells and 2.87 M rows (88.1 % copy, 2.3 %
restriction, 9.7 % prolongation): the row count per GPU of the frozen AMR sphere of [#133](https://github.com/szaghi/adam/issues/133)
on 4 A30.

![The fill on the 2:1 grid for every prolongation and restriction regime, split by kind of row, H100 and A30](/developer/amr-ghost-kinds.png)

**The 2:1 cost is the prolongation.** Under 10 % of the rows, the prolongation takes 54 % (tricubic) to 85 % (weno6) of
the fill on the H100; the copies cost what they cost on a uniform grid; a weno6 restriction costs as much as all the
copies.

![The cost of one row of each kind, with the density guard on (bars) and off (dots)](/developer/amr-ghost-row-cost.png)

**One row** (5 variables, H100): a copy 0.36 ns; a restriction 1.5 (mean) to 15.5 ns (weno6); a prolongation 5.0
(tricubic) to 23.0 ns (weno6): 14 to 64 copy rows. The A30 is 2.6–4.3 times slower per row, with the same proportions.
**The density guard is most of the prolongation**: without it a tricubic row costs 0.84 ns instead of 5.0 (83 % of the
row is the guard; 85 % on the A30), a weno6 row 10 instead of 23. The guard re-interpolates the density, and re-gathers
its footprint, for each of the other four variables of every fine ghost: the same density work five times.

![The variants on the 2:1 grid, the default pair and the WENO pair, H100 and A30](/developer/amr-ghost-variants.png)

- **Branch divergence is minor on the H100** (V0s: −10 % to +7 %), larger on the A30 (up to −38 % where both passes
  agree).
- **The descriptor kernels do not carry the uniform gain to AMR.** D2 makes the copies only 18 % faster (5,310 small,
  irregular descriptors in 154 shape classes; the owner mask and the search paid per cell) and the prolongation
  1.2–2.1 times slower than V0 (the seam metadata the rows carry precomputed is recomputed per cell).
- **Sharing a footprint per anchor in one thread loses**: D3 is 2.5 times slower than V0ng's prolongation on the H100,
  2.2 times on the A30; eight children computed serially with a 2.7 KB private footprint lose more parallelism than the
  shared reads save.
- At four times the load (960 blocks, 13.3 M cells, 11.3 M rows) the per-row costs and the guard's share are the same.

**Production check.** On the A30 the copy rows cost 2.44 ms for 2.53 M rows, FLUME's uniform launch 2.6 ms for a
similar count: the copies reproduce FLUME. FLUME's AMR launch on the frozen sphere costs 14.6 ms against 7.8 ms here at
the same total row count; with the per-row costs above that implies about 26 % prolongation rows in the sphere's grid
(three levels hugging the body), against 9.7 % in the box. The sphere's row mix itself was not measured.

**What it says for the 2:1 fill.**

1. **The density guard is the first target**, before any change of data structure: compute the density, its bound and
   the guard decision once per fine ghost and reuse them for the other variables (the same operations, done once:
   bitwise). It is worth up to half of the 2:1 fill. **Applied** ([below](#applied-the-density-guard-once-per-ghost)).
2. **The prolongation operators are arithmetic** (5–23 ns per row): they belong to the first track. Sharing the
   children's footprints has to be cooperative (threads of a block through shared memory); the serial per-anchor thread
   is ruled out.
3. **Group the rows by kind**: free with today's maps, bitwise, up to 38 % on the A30. **Measured with the guard
   once: slower on the H100, not applied** ([below](#applied-the-density-guard-once-per-ghost)).
4. **Descriptors need structure, not only compression**: few, regular, non-overlapping by construction (no owner mask),
   member blocks contiguous, the seam metadata precomputed per class. **Superseded for the copies by C4**
   ([below](#applied-the-copy-rows-in-their-own-launch)): the rows themselves, packed and sorted with the receiving block
   fastest, give the copies their 4x without descriptors.

### Applied: the density guard once per ghost

[#157](https://github.com/szaghi/adam/issues/157) B2 (2026-11-25) splits the device guard as the CPU's: a first
kernel over the coarse-to-fine rows computes each fine ghost's clamped density once (`seam_guard_density_gpu_dev`) into
a scratch value per row, and the fill applies it to each variable (`seam_guard_value_gpu_dev`); with the guard off the
fill runs a kernel without the guard's code. The benchmark checks it in every run: bit for bit against the b8641a59
kernel (V0p) on every regime pair and both GPUs, and against the CPU fill.

![#157 on the production-like grid: the b8641a59 fill, the guard once per ghost, and the rows sorted by kind; guard on and off, H100 and A30](/developer/amr-ghost-i157.png)

| Production-like grid (ms) | H100 before | H100 #157 | A30 before | A30 #157 |
|---|---|---|---|---|
| tricubic / mean (the defaults) | 2.54 | 1.67 (−34 %) | 7.83 | 6.77 (−14 %) |
| quintic / lagrange6 | 3.94 | 2.37 (−40 %) | 13.58 | 7.16 (−47 %) |
| quintic-mp / mean | 5.18 | 2.94 (−43 %) | 17.26 | 8.96 (−48 %) |
| weno6 / weno6 | 8.57 | 5.22 (−39 %) | 31.61 | 19.76 (−37 %) |

- **Every grid and regime pair is faster: 20–43 % on the H100, 3–50 % on the A30.** The tricubic fill on the H100 is
  now 0.46 ms above its guard-off floor, against 1.39 ms before.
- **With the guard off** (PRISM, and FLUME's non-state fields), against the b8641a59 kernel in the same run, the 6-point
  regimes gain 8–20 % on the H100 and 5–19 % on the A30; tricubic moves by −0.3 to +6.3 % on the H100 (−12 to +2 % on
  the A30): the b8641a59 kernel carried the guard's code unused, and that compiled shape happened to suit the
  copy-dominated tricubic fill on the H100.
- **Grouping the rows by kind does not pay with the guard once.** The same kernel on the rows sorted by kind is slower
  on the H100 in 18 of 22 runs (guard off up to +16 %): the heavy prolongation threads become the tail of the launch,
  while in ADAM's order they run beside the copies. The A30 still prefers the sorted rows (up to −33 %). The H100 is the
  reference, and the guard-once kernel does not depend on the order, so the rows keep ADAM's order.

Per-case numbers, the nvfortran hazards met on the way and the repeatability notes are in the benchmark's README
(`src/tests/perf/ghost-fill/amr/README.md`, "#157").

### Applied: the copy rows in their own launch

[#157](https://github.com/szaghi/adam/issues/157) C (2026-11-26). The plan proposed structured descriptors for the copy
rows; the benchmark found a shorter way. What gave the uniform benchmark its 3–5x is **the receiving block running
fastest** (consecutive threads write consecutive memory of the block-fastest field), not the absence of a map. New
variants, every one bitwise against V0:

| Variant | What it runs | A30, production-like grid, tricubic / mean |
|---|---|---|
| K0 | the copy rows alone, ADAM's kernel and order | 2.64 ms |
| Ka | the copy rows as int32 donor and receiver offsets, ADAM's order | 1.54 |
| Kb / Kd | the same sorted by receiver / donor offset (the block fastest) | **0.53 / 0.51** |
| Kp | the copy rows packed as [block, i j k] int32 pairs, sorted as Kb, the field indexed directly | **0.65** |
| Ca–Cg | every kind of row in one compact launch, in ADAM's order or sorted (by receiver, by kind, heavy rows first or last) | 5.2–12.6 |
| Ci / Cj | ADAM's kernel on the restriction and prolongation rows in ADAM's order, then Kb / Kp | **3.71 / 3.84** |

- **Sorting the copies gains 4x, compressing them alone 1.7x.** Sorting the *prolongation* rows loses: the children of
  a coarse anchor stop being neighbours and their shared footprints leave the caches (Cb, every row sorted by receiver,
  is twice V0).
- **Mixing the kinds in one launch wastes the gain**, whatever the order: the copies need their own launch.
- **The dedupe never removes a row** on any FLUME or PRISM regression case nor on the #133 grids (its message never
  appears in the logs): the rows are non-overlapping already, so any split and order is bitwise.

**Applied as Cj** (C4): `maps_fnl_object%copy_cpu_gpu` derives, from the same CPU map, the heavy rows (restriction and
prolongation, the map's order, int64) and the copy rows (four int32 `[b_send, ijk_send, b_recv, ijk_recv]`, ijk packed
10 bits each, sorted by receiver cell then block); the local fill runs `update_ghost_local_gpu_dev` on the first and
`update_ghost_local_copy_gpu_dev` on the second. Bitwise against the b8641a59 kernel on the full map in every run, and
on every FLUME and PRISM regression case (A/B), on both GPUs.

| Production-like grid (ms) | A30 before C4 | A30 C4 | H100 before C4 | H100 C4 |
|---|---|---|---|---|
| tricubic / mean | 5.86 | **3.88** (−34 %) | 1.67 | **1.04** (−37 %) |
| weno6 / weno6 | 22.3 | 20.5 (−8 %) | 5.18 | 4.56 (−12 %) |
| guard off, tricubic | 3.52 (b8641a59) | 1.65 | 1.21 | 0.58 (−53 %) |

On the H100 (SpaceHPC job 184729, HEAD 80582712 against C4 alternated on one GPU) every case gains, the split launch
costs nothing there: quintic / lagrange6 2.37 → 1.88 ms (−21 %), tricubic on 16³ blocks 0.396 → 0.277 (−30 %), at 4
times the load (960 blocks) 4.58 → 2.32 (−49 %); the copies alone 0.90 → 0.26 ms (3.5x).
The A30 figures of the C4 fill are its two launches timed apart and the fill timed again at the end of the run
(3.24 + 0.65 = 3.88 ms): **the first variant a run times is not warm** (V0, 5.37 ms in the same run). The harness now
fills for a second before V0; the figures of the earlier sections, each run's first, are biased high by up to that much,
so the gains reported against V0 there were conservative.

**After C4 the 2:1 fill is the prolongation and its guard**: on the A30 the copies are 0.65 of 3.88 ms; the guard alone
(3.88 against 1.65 guard-off) is 2.2 ms. That is where the cooperative prolongation (#157 D) and the rest of the guard
have to act.

## What it does not cover

- **Real AMR layouts.** The 2:1 benchmark refines a box; the row mix of a grid that hugs a body (the sphere) is
  inferred, not measured, and three-level grids are not built.
- **The boundary-condition rows, the blend and the MPI pack/unpack.** They use the same per-row maps; they are expected to
  carry a similar factor and are the next measurements.
- **Several realms, empty ranks, an immersed body.** The fill here is the same-level local fill only.
- **Other vendors.** Plain OpenACC on NVIDIA here; the AMD path (FUNDAL's OpenMP offload) was measured by Giacomo Rossi
  on an MI300X (#157): there an array constructor with variable elements in device code was a device `malloc`, which
  made the 2:1 fill 7.8 times the H100's; fixed (element assignments, bitwise; nvfortran's quintic-mp fill gained 14 %
  on the A30 too), and `src/tests/lint/check-device-array-constructors.sh` keeps it out.

## Future paths it indicates

1. **Ghost exchange by slabs, not by cells.** The maps should describe (receiver block, donor block, direction, ranges,
   shift) per neighbour pair, a few dozen descriptors per block instead of thousands of rows, and the kernels compute the
   cell indices. For a uniform same-level neighbourhood this is V3/V5. The AMR redesign has to give every block this kind
   of description of its neighbourhood, including the 2:1 faces.
2. **2:1 seams: the guard and the operators first.** The 2:1 benchmark shows the seam rows are arithmetic (the
   prolongation operators and, above all, the density guard), not data movement: once the guard runs once per ghost
   (#157 B2) and the copies have their own sorted launch (#157 C4), what is left is the guard's velocity bound and the
   footprint reads that the cooperative prolongation (#157 D) would share.
3. **One pass per refresh.** Local fill, boundary conditions and the blend run as three passes over three maps today;
   with slab descriptors they can be one fused pass per ghost refresh, which matters most where the refreshes are many
   (12 per step with the immersed body).
4. **MPI pack and unpack by slabs.** The send and receive buffers are filled through the same kind of per-row maps; slab
   descriptors per neighbour rank give contiguous pack kernels.
5. **A selectable launch strategy** (fused or per direction), chosen once per GPU and load.
6. **The first track keeps the layout.** With the ghost machinery near its floor, the face kernel (local memory, register
   pressure, shared-memory stencils) is free to be rewritten for its own sake, without a layout change imposed by the
   ghosts.

## Reproduce

```bash
cd src/tests/perf/ghost-fill
./run.sh cc90,cuda12.9 results/h100-<id>.txt     # one GPU; nvfortran on PATH (inside jobs/spacehpc.pbs on SpaceHPC)
./run.sh cc80,cuda12.9 results/a30-<id>.txt      # jobs/medooza.slurm on medooza
cd ../../../..   # the repository root
exe/.regression-venv/bin/python src/tests/perf/ghost-fill/plot.py --out docs/public/developer --theme both
```

`plot.py` reads the newest `results/h100-*.txt` and `results/a30-*.txt` and the H100 phase reports of #133
(`src/tests/flume/speedup/results/2026-11-25-spacehpc/phase/reports/`), and writes the first three figures of this page
with their dark twins.

The 2:1 benchmark builds through FoBiS (mode `test-amr-ghost-bench-fnl-nvf`) inside its run script:

```bash
cd src/tests/perf/ghost-fill/amr
./run.sh spacehpc_nvf results/h100-<id>.txt      # jobs/spacehpc.pbs (RUN_ONLY=large: the larger load only)
./run.sh iac_nvf results/a30-<id>.txt            # jobs/medooza.slurm
cd ../../../../..   # the repository root
exe/.regression-venv/bin/python src/tests/perf/ghost-fill/amr/plot.py --out docs/public/developer --theme both
```
