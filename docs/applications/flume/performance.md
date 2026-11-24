# Performance: GPU and CPU nodes, uniform-grid scaling

What a GPU node buys over a CPU node, what sets the cost of a cell on an A30, and how FLUME scales inside a node and
across nodes. **Scaling is measured and claimed on uniform 3-D grids only.** Runtime AMR (the host regrid, 2:1 seams,
reflux, the runtime rank layout) and the 2-D quadtree path (ghost layers along the null axis) are implemented for
accuracy and not yet optimised for performance: their timings are kept as baselines at the end of the page
([The accuracy-first paths](#the-accuracy-first-paths-runtime-amr-and-2-d)), with no scaling claim.

The measurements are the speedup matrix of [#133](https://github.com/szaghi/adam/issues/133) and its grid study; the
harness, the protocol and every table are in
[`src/tests/flume/speedup/README.md`](https://github.com/szaghi/adam/blob/develop/src/tests/flume/speedup/README.md),
the raw tables in its `results/2026-11-23/` and `results/2026-11-24/`. The kernel-level profile of one GPU (where a
step's time goes inside the A30) is the companion study
[`profile-couette3d`](https://github.com/szaghi/adam/blob/develop/src/tests/flume/speedup/profile-couette3d/README.md).

## The matrix

Three flows (*legs*), each on a uniform grid and on two AMR grids that refine to the uniform grid's spacing, so that
the solution where the AMR grid is finest and the global time step are the same: the per-step ratio is the
time-to-solution ratio.

| Leg | Flow | Grid | Scheme | Uniform cells |
|---|---|---|---|---|
| 1 | Woodward–Colella Mach 3 forward-facing step, three realms glued by `auto` seams | 2-D quadtrees, h = 1/1280 (L4), runtime AMR from L1 | WENO-7 characteristic, SSP-54 | 4.13 M |
| 2 | Mach 2 shock diffraction over a sphere, immersed boundary (Euler wall) | 3-D octree, 24³ blocks, h = 1/384 (L4), runtime AMR from L2 | WENO-5 characteristic, SSP-33 | 56.6 M |
| 3 | Coleman–Kim–Moser channel, M_b 1.5, Re_b 3000, isothermal walls, Re_τ ≈ 202 | 3-D octree, 32×24×16 blocks, L3, runtime AMR from L2 (`r = h/η` marker) | `hybrid-kep-weno`, viscous, SSP-54 | 6.29 M |

Each leg has three cases: **uniform** (U), **runtime AMR** (R: regrid every 10 steps) and **frozen AMR** (F: R's grid
and state with `[amr] frequency = 0`). Every backend and rank count starts from the same restart dump (written once,
repartitioned for the reading rank count). A step cost is the median over three repetitions of
$(T_{\rm loop}(110) - T_{\rm loop}(10))/100$, every output off; the grind time is
$g = \text{step cost} \times \text{devices} / N$ (ns per cell-step, devices = GPUs or cores).

The scaling cases are the uniform sphere (2a and the grid study's `wg-*` cases, the same flow on other grids) and
the uniform channel (3a).

**Hardware.** medooza `prod-gn`: a node holds 4 NVIDIA A30 (24 GB) and two Intel Xeon Gold 6338 (64 cores, 4 NUMA
nodes). A *GPU node* is 4 A30 with one MPI rank per GPU (FNL, OpenACC); a *CPU node* is 64 MPI ranks (CPU backend).

## GPU node against CPU node

T(CPU node) / T(GPU node) on the uniform grids:

| Case | Grid | Blocks | GPU node / CPU node |
|---|---|---|---|
| 2a sphere (3-D, WENO-5, IB, SSP-33) | 56.6 M cells, L4 | 4096 of 24³ | 16.3 |
| 3a channel (3-D, hybrid, viscous, SSP-54) | 6.29 M cells, L3 | 512 of 32×24×16 | 13.8 |
| sphere, grid study | 7.08 M cells, L3 | 512 of 24³ | 11.9 |
| sphere, grid study | 28.1 M cells, L3 | 512 of 38³ | 12.6 |

A GPU node runs the uniform 3-D flows 12 to 16 times faster than a CPU node; on the same flow the factor moves with
the grid. On the CPU node the sphere costs 14.4 µs per cell-step per core on 2a, 12.4 on the 24³ grid and 10.2 on the
38³ grid: 2a's 64 blocks per rank carry the immersed boundary's load imbalance, which the GPU node, with 1024 blocks
per A30, does not feel. The two grid-study grids load the GPU node with 1.8 and 7.0 M cells per A30, below 2a's
14.2 M ([what a cell costs](#what-a-cell-costs-on-one-a30)). Larger blocks make both nodes cheaper per cell.

## What a cell costs on one A30

![Grind time of the uniform sphere against the block edge: on 1 A30, and the same grids on 1, 2 and 4 A30](/flume/perf-grid-study.png)

The uniform sphere on grids of 2.1 to 28.1 M cells, 1 A30, from 20-step cold starts (the per-cell cost does not
depend on the flow state: 0.5 %, measured on 2a):

| Block edge b (L3, 512 blocks) | 16 | 20 | 24 | 28 | 32 | 36 | 38 | 40 |
|---|---|---|---|---|---|---|---|---|
| cells | 2.10 M | 4.10 M | 7.08 M | 11.2 M | 16.8 M | 23.9 M | 28.1 M | 32.8 M |
| g, ns per cell-step | 63.1 | 54.4 | 51.0 | 47.7 | 46.2 | 44.4 | 43.9 | does not fit |

- **The block size sets the cost of a cell.** The same 192³ grid costs 66.0 ns per cell-step as 4096 blocks of 12³,
  51.0 as 512 of 24³ and 44.0 as 64 of 48³: a block of b³ cells stores (b + 6)³ (WENO-5, 3 ghost layers), and the
  ghost cells are filled and partly computed every stage. Above a few million cells per A30 the load itself matters
  little: 7.1 M cells in 48³ blocks cost what 28.1 M cells in 38³ blocks cost.
- **The matrix's 24³ blocks leave ~16 % on the table**: 51.0 against 43.9–44.0 ns with 38³–48³ blocks.
- **An A30 holds 28.1 M cells of the sphere** (38³ blocks, 88 % of the blocks budget); 32.8 M is refused at
  initialisation (the budget, 0.85 of the 25.2 GB the A30 reports, holds 510 blocks of 46³ stored cells, 512 needed).

## Inside a GPU node

The same grids on 1, 2 and 4 A30 of one node (the right panel of the figure above):

| b | 16 | 20 | 24 | 28 | 32 | 36 | 38 |
|---|---|---|---|---|---|---|---|
| g, 2 A30 / g, 1 A30 | 1.09 | 1.05 | 1.04 | 1.03 | 1.02 | 1.02 | 1.02 |
| g, 4 A30 / g, 2 A30 | 1.30 | 1.25 | 1.22 | 1.18 | 1.19 | 1.16 | 1.13 |

Two A30 cost 2–9 % more per cell than one; four cost another 13–30 %. The node explains the step: its A30 are
NVLink-bridged in pairs (GPU 0–1 and GPU 2–3, 4 × 25 GB/s each way), while between the pairs, and from every GPU to
the InfiniBand card, the path is PCIe plus the inter-socket link. The four GPUs hang off NUMA nodes 0, 1, 3 and 3.

**Rank placement is not the cause.** The default launch leaves three of the four ranks on host cores away from their
GPU's NUMA node (ranks 1 and 3 on NUMA 2, which has no GPU). Pinning each rank to its GPU, to 8 cores of that GPU's
NUMA node and to its memory changes the grind time by −0.8 to 0 % on one and on two nodes (2a, and the 24³ and 38³
grids), inside the spreads. The halo exchange posts MPI on device buffers (no host staging), so the host core is off its
path, and the cost of the fourth A30 is the link between the NVLink pairs. Pinning does one useful thing: unpinned, every rank also opens a context on GPU 0, which then has
0.9 GB less free memory than the others; pinned, all four have the same.

## Scaling across nodes

**Strong scaling of the uniform sphere** (2a, 56.6 M cells, 24³ blocks; step cost in s, efficiency against one node):

| Node | 1 | 2 | 4 | 8 |
|---|---|---|---|---|
| GPU nodes (4 A30 each) | 0.780 | 0.491 | 0.260 | 0.167 |
| efficiency | 1 | 0.79 | 0.75 | 0.58 |
| cells per A30 | 14.2 M | 7.08 M | 3.54 M | 1.77 M |
| CPU nodes (64 cores each) | 12.7 | 6.48 | 3.24 | – |
| efficiency | 1 | 0.98 | 0.98 | – |

The GPU efficiency at 8 nodes is set by the load as much as by the network: 1.77 M cells per A30 in 24³ blocks is
where the one-node curve already rises (a 4-A30 node at that load costs 64.8 ns per cell-step against 55.1 at
14.2 M). The largest run is the L5 sphere, 453 M cells on 32 A30 (8 nodes), at 0.98 s per step.

**Weak scaling of the uniform sphere** (g in ns per cell-step; inter-node efficiency = g on one 4-A30 node / g on
8 nodes at the same cells per A30 and the same blocks: an octree level is 8 times the cells, so 4 A30 pair with 32):

| Cells per A30 | Blocks | One node, 4 A30 | 8 nodes, 32 A30 | Inter-node efficiency |
|---|---|---|---|---|
| 14.2 M | 24³ | 55.1 (2a, 56.6 M cells) | 69.3 (L5, 453 M cells) | 0.80 |
| 7.02 M | 38³ | 50.5 (L3, 28.1 M cells) | 67.0 (L4, 225 M cells) | 0.75 |

- **Leaving the node costs 20–25 % per cell on 8 nodes** at 7–14 M cells per A30. The pair at 38³ blocks loses more
  than the one at 24³ although its blocks carry fewer ghost cells: its cells are 15 % cheaper, and at half the load
  per A30 the same exchange weighs more. In throughput the two 8-node runs are equal (478 and 462 M cell-steps per
  second).
- **The weak efficiencies of the first version (0.72–0.73 at 3.5–7 M cells per A30) mixed two costs**: their 1- and
  2-A30 references left out the cross-pair cost of the fourth A30 (13–22 % at those loads), which belongs to the node,
  not to the network.
- **The CPU node scales almost perfectly** (0.98 on 2 and 4 nodes), but each CPU node is 12–16 times slower.

## Device memory

FNL allocates every block-sized device array at the block capacity that its budget computes from the device memory
(`memory_fraction` × memory / realms / (fields × 8 B × stored cells per block)), not at the blocks in use. Measured
with three values of `ADAM_MEMORY_FRACTION` per case on 4 A30 (`nvidia-smi` peaks; the peak is a straight line in the
stored-cell slots of the capacity, residuals below 5 MiB):

| Case | Bytes per stored-cell slot (doubles) | Budget's fields per slot | Memory outside the capacity, per GPU | Needed per A30 at its layout |
|---|---|---|---|---|
| step (3 realms, WENO-7, SSP-54) | 404.9 B (50.6) | 64 | 3.40 GiB | 7.7 GiB (1a, 4 A30) |
| sphere (IB, WENO-5, SSP-33) | 392.4 B (49.1) | 54 | 1.95 GiB | 12.1 GiB (2a, 4 A30) |
| channel (hybrid, viscous, SSP-54) | 493.4 B (61.7) | 69 | 0.75 GiB | 6.7 GiB (3a, 1 A30) |

The budget overestimates a slot by 10–26 % (it counts the face fluxes as full stored blocks), and at its default
fraction 0.85 an A30 keeps 2.5–4.9 GiB free. GPU 0 holds 0.87 GiB more than the others in every run (not identified).
The budget, not the device, sets the sphere's ceiling of 28.1 M cells per A30.

## The accuracy-first paths: runtime AMR and 2-D

Runtime AMR and the 2-D quadtrees are implemented for accuracy; their performance is not optimised yet, and the
numbers below are baselines for that work, not scaling results.

![Step cost of the AMR grids over the uniform grid's, split inside and outside the time steps](/flume/perf-amr-decomposition.png)

| Leg | Node | AMR speedup T(U)/T(R) | Ideal N(U)/N(R) | Frozen grid T(U)/T(F) | Grid efficiency | T(R)/T(F) − 1 |
|---|---|---|---|---|---|---|
| 1 | GPU (4 A30) | 1.04 | 1.59 | 1.40 | 0.89 | 34 % |
| 1 | CPU (64) | 1.06 | 1.59 | 1.42 | 0.90 | 34 % |
| 2 | GPU (4 A30) | 1.48 | 4.02 | 2.04 | 0.51 | 38 % |
| 2 | CPU (64) | 2.20 | 4.02 | 2.73 | 0.68 | 24 % |
| 3 | GPU (4 A30) | 1.25 | 1.78 | 1.60 | 0.90 | 29 % |
| 3 | CPU (64) | 1.03 | 1.78 | 1.43 | 0.81 | 39 % |

The grid efficiency is $[T(U)\,N(F)/N(U)]/T(F)$: 1 when an AMR cell costs what a uniform cell costs. Where the AMR
cost goes:

- **A cell of the sphere's AMR grid costs twice a uniform cell on the A30** (grid efficiency 0.51): at the same
  254–256 blocks per GPU the face fluxes cost the same on both grids, while the intra-rank ghost fill (one map-driven
  kernel that copies same-level ghosts and, on the AMR grid, also interpolates coarse-to-fine and restricts
  fine-to-coarse ones) costs 5.6 times as much: 96 % of the difference.
- **On the GPU node the runtime overhead is the regrid, which runs on the host** at every regrid step (0.26–1.4 s
  when the grid does not change, 2.0–2.5 s when it does); it scales with the blocks per GPU and stays 22–27 % of a
  runtime step.
- **On the CPU node the regrid is cheap** (1.3 % of the sphere's runtime step) **but the runtime layout is not**: a
  runtime regrid distributes whole sibling sets and splits every realm from rank 0, so the busiest rank keeps the
  uniform grid's load (the step's runtime case: 64,512 cells against a mean of 40,680). The runtime step then costs
  what the uniform step costs, and does not shrink with more ranks (1.19 s on 64, 128 and 256 ranks).
- **The 2-D step's cell costs 1.60 cells of the same step extruded in 3-D** at the same order and load (202 against
  126 ns per cell-step): the boundary-condition kernel refills the `ngc` ghost layers along the null axis at every
  stage (59 % of the 2-D kernel time). The same ghosts cap the step at about 1.7 M cells per A30. A GPU node runs the
  step 5.8 times faster than a CPU node (uniform 5.83, frozen 5.75, runtime 5.72), against 14–16 on the 3-D uniform
  grids.

Timings of these paths on more nodes (baselines; step cost in s on 4 / 8 / 16 / 32 A30 and 64 / 128 / 256 cores):

| Case | GPU nodes 1 / 2 / 4 / 8 | CPU nodes 1 / 2 / 4 |
|---|---|---|
| 1a step uniform (2-D) | 0.216 / 0.119 / 0.070 / – | 1.26 / 0.642 / 0.328 |
| 1b step frozen | 0.154 / 0.089 / 0.057 / – | 0.887 / 0.544 / 0.334 |
| 1c step runtime | 0.207 / 0.119 / 0.076 / – | 1.19 / 1.18 / 1.18 |
| 2b sphere runtime | 0.526 / 0.350 / 0.312 / 0.188 | 5.76 / 3.95 / 3.06 |
| 2c sphere frozen | 0.381 / 0.284 / 0.257 / 0.161 | 4.64 / 3.13 / 1.97 |

![Grind time against cells per device, every case of the matrix](/flume/perf-grind.png)

## What the measurements say about the hypotheses of the plan

The plan of [#133](https://github.com/szaghi/adam/issues/133) stated four hypotheses before the runs, three of them on
the accuracy-first paths:

| Hypothesis | Verdict | The measurement |
|---|---|---|
| H1: on the step, runtime AMR falls well short of its ideal on GPUs (launch floor), and comes close on CPUs | **reformulated** | short of the ideal on both nodes (1.04 and 1.06 of 1.59); no launch floor (the step scales to 0.26 M cells per A30); the losses are the AMR grid and the host regrid on the GPU node, the runtime layout on the CPU node |
| H2: a 2-D cell costs more than two thirds of a 3-D cell at the same order, by the null-axis ghosts | **confirmed** | at the same order and 1.03 M cells per A30 a 2-D cell costs 1.60 3-D cells (202 against 126 ns); the boundary-condition kernel, which refills the null axis's ghost layers every stage, takes 59 % of the 2-D step's kernel time against 10 % in 3-D |
| H3: on the sphere, the GPU runtime overhead is the regrid's host round trip; on the CPU it stays a few percent | **confirmed** | 38 % on 4 A30, all outside the steps, paid at every regrid step (1.4 s) whether or not the grid changes; 1.3 % on the CPU node, whose runtime overhead is inside the steps (half of it the runtime layout) |
| H4: on the channel, T(runtime) ≈ T(frozen) + regrid, the gain from the core alone | **confirmed on the GPU node, refuted on the CPU node** | 4 A30: 0.1432 + 0.0414 = 0.1846 s against 0.1842; the grid coarsens the core only. 64 cores: the runtime layout leaves the busiest rank with the uniform grid's load (gain 1.03) |

## Limits

- One machine (A30, Xeon Gold 6338); the AMD MI300X node of the plan is not measured yet.
- The matrix and its multi-node runs use 24³ blocks; the grid study shows 38³–48³ blocks are ~16 % cheaper per cell
  on one A30, so the multi-node numbers are not the best this hardware gives.
- Up to 8 nodes (32 A30), the most medooza schedules.
- The windows are short (110 steps); two residual biases were measured and are not in these numbers: printing the
  progress lines on every rank costs nothing measurable, and computing the conservation integrals every step costs
  0.9 % on 4 A30.

## Open issues found

| Cost | Where it shows | Proposal |
|---|---|---|
| The regrid runs on the host at every regrid step | 22–27 % of a runtime step on 4 A30, also when the grid does not change | per-phase timers, then marking on the device |
| The runtime layout keeps sibling sets whole and splits every realm from rank 0 | the CPU runtime step stops at 64 ranks; legs 1 and 3 gain nothing from runtime AMR on the CPU node | offset or joint split of the realms, weighted by cells; split sibling sets between regrids |
| A cell of the sphere's AMR grid costs twice a uniform cell on the A30 | grid efficiency 0.51 on 4 A30 | one kernel per kind of ghost row (copy, interpolation, restriction), then a cheaper interpolation |
| A 2-D block stores ghost layers along its null axis; realms share the device equally | the step fits ~1.7 M cells per A30 | no ghosts on a null axis; device share in proportion to stored cells |
| The FNL job's free-memory check sees the batch node only | multi-node jobs | check every node |

*Measured: ADAM `0bcdf870` + the night-2026-11-22 tree (since committed, `6f5e57a1`..`b06f02d3`), CPU gfortran
16.2.0 + Open MPI 5.0.11 (`iac_gnu`) on medooza (Intel Xeon Gold 6338, 64 ranks per node), FNL nvfortran 26.9, CUDA
12.9, LTO (`iac_nvf`) on NVIDIA A30 (4 per node): the matrix 2026-11-23 (jobs 1369400-1369597, 1369616), the grid study and the
rank placement 2026-11-24 (jobs 1369668-1369679, 1369687-1369691); `speedup/run.sh` windows, `analyze.py`, `matrix_plot.py`
and `grid_study_plot.py`.*
