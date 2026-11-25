# Ghost-fill microbenchmark

The documentation page, with the figures and the analysis: [Ghost fill on the GPU](../../../../docs/developer/ghost-fill-benchmark.md)
(published at `/developer/ghost-fill-benchmark`). `plot.py` draws the figures from `results/` and the H100 phase
reports of #133 (`exe/.regression-venv/bin/python src/tests/perf/ghost-fill/plot.py --out docs/public/developer --theme
both`, from the repository root).

**Question** (AMR redesign, 2026-11-25): how far ADAM's local ghost fill is from the memory floor of a GPU, and which
design choice explains the gap. Without an immersed body the ghost work of a FLUME step is at most ~20 % of it on an
H100 (`src/tests/flume/speedup/README.md`, "ESA SpaceHPC"), but it is pure data movement: as the face-flux kernel gets
faster (track 1), its share grows, so the ghost machinery has to reach its own floor (track 2).

`ghost_fill_bench.F90` is standalone (plain OpenACC, no ADAM dependency): a uniform periodic brick of blocks of b^3
cells, `ngc` ghost layers, 5 variables; every ghost cell (faces, edges, corners) filled from its same-level neighbour,
the rows ADAM's local map holds. Every variant is checked bit for bit against a host fill after its timing.

| Variant | Kernel | Isolates |
|---|---|---|
| V0 | ADAM's `update_ghost_local_gpu_dev` (`adam_fnl_field_kernels.F90`): one row per ghost cell, 10 `integer(I8P)` columns of which 9 are read, `collapse(2)` with the variable loop outermost, the branch on column 9, layout `q(b,i,j,k,v)`, rows in `adam_maps_object` order | today |
| V1 | V0 with `integer(int32)` columns | the bytes of the indices |
| V2 | one donor and one receiver linear offset per row, the variables inside the row | index traffic and the map re-read per variable |
| V3 | structured slabs: no per-cell map, one kernel per neighbour direction (26) over (v, k, j, i, block), the donor block from a per-block neighbour table, the donor cell by a constant shift; layout `q(b,i,j,k,v)` | the per-cell map |
| V4 | V3 on the layout `q(i,j,k,v,b)` (cell index fastest) | the layout |
| V5 | V3 fused into one launch: every ghost cell of a block numbered direction by direction, direction and cell decoded from a 26-entry table; layout `q(b,i,j,k,v)` | the launch count |
| V6 | V5 on the layout `q(i,j,k,v,b)` | the layout, one launch |
| ref | a contiguous copy of as many values | the floor the GPU reaches |

Bandwidth counts 16 B per ghost value (one read, one write); index reads are not counted, so V0's figure is what the
data alone would need.

## Runs (2026-11-25)

`run.sh <-gpu target> <out>` builds and runs the matrix: ~7 M cells per GPU at block edges 16, 24, 32, 38, 48; ~56 M
at 24, 38, 48; 38^3 with `ngc` 4 (WENO-7); median of 50 fills. Jobs: ESA SpaceHPC, one NVIDIA H100 (`gpu-007`,
nvfortran 26.3, `-gpu=cc90,cuda12.9`, job 184030, `jobs/spacehpc.pbs`); medooza, one NVIDIA A30 (`gnode23`, nvfortran
26.9, `-gpu=cc80,cuda12.9`, job 1369794, `jobs/medooza.slurm`). Raw output in `results/`. Every variant, every size,
both GPUs: bitwise identical to the host fill.

**V0 reproduces production.** The phase timer of the FLUME step (`speedup/results/2026-11-25-spacehpc/phase/`,
`speedup/results/2026-11-24/phase/`) gives the local fill per call: H100 9.5 ms at 56 M cells per GPU in 38^3 blocks
(V0 9.09 ms), 1.12 ms at 7 M (V0 1.10 ms); A30 3.1 ms at 7 M (V0 3.27 ms). The benchmark's V0 is ADAM's fill within 5 %.

**Time of one fill (ms), and fraction of the copy floor** (H100 first, the reference):

| Cells per GPU, b | copy floor | V0 (ADAM) | V1 I4P | V3 slabs | V5 slabs, 1 launch | V6 cell-fastest, 1 launch | best / V0 |
|---|---|---|---|---|---|---|---|
| **H100** 7.1 M, 16 | 0.436 | 2.775 (16 %) | 1.956 | 1.197 | **0.809 (54 %)** | 1.051 | 3.4x |
| 7.1 M, 24 | 0.268 | 1.691 (16 %) | 1.201 | 0.914 | **0.490 (55 %)** | 0.707 | 3.5x |
| 7.1 M, 32 | 0.196 | 1.313 (15 %) | 0.966 | 0.828 | **0.354 (55 %)** | 0.449 | 3.7x |
| 7.0 M, 38 | 0.164 | 1.096 (15 %) | 0.809 | 0.774 | **0.291 (56 %)** | 0.316 | 3.8x |
| 7.1 M, 48 | 0.131 | 0.898 (15 %) | 0.672 | 0.710 | **0.230 (57 %)** | 0.255 | 3.9x |
| 56.6 M, 24 | 2.001 | 13.85 (14 %) | 10.04 | **3.594 (56 %)** | 3.785 | 5.611 | 3.9x |
| 56.2 M, 38 | 1.156 | 9.094 (13 %) | 6.773 | 2.280 | **2.187 (53 %)** | 2.379 | 4.2x |
| 56.6 M, 48 | 0.900 | 7.702 (12 %) | 5.816 | 1.890 | **1.707 (53 %)** | 1.920 | 4.5x |
| 7.0 M, 38, ngc 4 | 0.221 | 1.640 (13 %) | 1.219 | 0.844 | **0.402 (55 %)** | 0.424 | 4.1x |
| **A30** 7.1 M, 16 | 1.130 | 7.686 (15 %) | 6.095 | **2.450 (46 %)** | 2.724 | 2.998 | 3.1x |
| 7.1 M, 24 | 0.671 | 4.755 (14 %) | 3.782 | 1.602 | **1.561 (43 %)** | 1.655 | 3.0x |
| 7.1 M, 32 | 0.480 | 4.122 (12 %) | 2.822 | 1.239 | **1.112 (43 %)** | 1.123 | 3.7x |
| 7.0 M, 38 | 0.394 | 3.269 (12 %) | 2.312 | 1.126 | **0.905 (44 %)** | 0.916 | 3.6x |
| 7.1 M, 48 | 0.309 | 2.605 (12 %) | 1.847 | 0.998 | **0.705 (44 %)** | 0.716 | 3.7x |
| 56.6 M, 24 | 5.272 | 38.30 (14 %) | 30.64 | **10.05 (52 %)** | 13.23 | 13.51 | 3.8x |
| 56.2 M, 38 | 3.031 | 24.35 (12 %) | 19.00 | **5.559 (55 %)** | 7.405 | 7.772 | 4.4x |
| 56.6 M, 48 | 2.372 | 23.90 (10 %) | 16.96 | **4.372 (54 %)** | 5.715 | 6.074 | 5.5x |
| 7.0 M, 38, ngc 4 | 0.544 | 4.667 (12 %) | 3.431 | 1.368 | **1.260 (43 %)** | 1.290 | 3.7x |

(The copy floor: H100 1.83-2.16 TB/s, A30 0.78-0.82 TB/s. V2 and V4 are slower than the best on every point and are
in `results/`: V2, the offsets with the variables inside the row, falls to 6-22 % of the floor, because each thread then
walks five addresses a whole variable apart.)

## What it says

- **Today's fill moves its data at 10-16 % of what the GPU can copy**, on both GPUs, at every block size and load: the
  per-cell index map is the cost, not the data.
- **Structured slabs reach 43-57 % of the floor and are 3.0-5.5 times faster than today**, bitwise identical. On the
  H100 one launch (V5) is best everywhere but at 24^3 / 56 M, where 26 launches (V3) are as good; on the A30 the 26
  launches win at 56 M (and at 16^3), one launch at 7 M. A design that can do either is right on both.
- **Narrowing the map (I4P) gains 1.3-1.4x**: index bytes matter, but removing the map gains three times more.
- **The layout is not the issue for the fill**: the cell-fastest layout (V6) is as fast as ADAM's block-fastest (V5) or
  slower. The block-fastest layout chosen for the face kernel's coalescing costs the ghost machinery nothing.
- **The remaining factor ~2 to the copy floor** is the shape of the work (ghost slabs 3-4 cells thin, gathers from 26
  neighbours), the floor of a fill, not of a copy.

What it means for a FLUME step (H100, phase timer of the uniform sphere, the local fill only, at today's face kernel):
38^3 blocks at 7 M cells per GPU, 3.4 -> ~0.9 ms per step without the body (-2.8 % of 90 ms), 13.4 -> ~3.5 ms with it
(-8 % of 119 ms); at 56 M, 28.4 -> ~6.8 ms (-3 % of 710 ms) and 113.5 -> ~27 ms (-10 % of 861 ms). The boundary rows
and the MPI pack/unpack use the same per-row maps and are not measured here. As track 1 shortens the face kernel the
same milliseconds weigh more.

**Not covered:** the 2:1 seam rows (interpolation, restriction), the boundary-condition rows, the MPI pack/unpack, an
unstructured (AMR) neighbourhood: here every block has 26 same-level neighbours.
