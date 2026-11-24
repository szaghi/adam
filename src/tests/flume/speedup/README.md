# FLUME speedup matrix

**Status: every leg measured on medooza (2026-11-23), single node and multi-node; results, figures, the memory
calibration, the residual biases and the verdicts on H1-H4 below; docs page
[Performance](../../../../docs/applications/flume/performance.md)** ([issue
#133](https://github.com/szaghi/adam/issues/133), the plan's only source of truth). Phases 0-3 are the harness;
phases 5-7 the legs; phase 9 (owner, 2026-11-23) lifts the plan's "no scaling beyond one node": strong scaling on up to
8 nodes, weak scaling at fixed cells per GPU, grind time. Not done: the AMD node (fnl-amd, Giacomo Rossi).

The matrix answers two questions: how much faster a GPU node is than a CPU node on each flow, and how much runtime AMR
buys at equal finest resolution, with where its overhead goes (the AMR grid itself against the regrids).

## The matrix

Three **legs** (flows), each on a uniform grid and on AMR grids refined to the same finest spacing (D1), on every
**backend** (`cpu`, `fnl-nvf`, `fnl-amd`). Every number below is in `matrix.conf`; `make_cases.py` writes the inputs.

| Case | Role | Flow | Dims, tree | Scheme | Grid | Finest h | Cells |
|---|---|---|---|---|---|---|---|
| `1a-step-uniform` | uniform | Woodward-Colella Mach 3 step, 3-realm forest, aligned blocks, `auto` seams | 2-D, quadtree | WENO-7 characteristic, SSP-54 | uniform L4, blocks A 48x16, B 48x64, C 192x64 | 1/1280 | 4.13 M |
| `1b-step-fixed` | frozen | same | same | same | 1c's grid frozen | 1/1280 | measured |
| `1c-step-runtime` | runtime | same | same | same | Loehner density marker (refine 0.3, derefine 0.1), regrid every 10 steps, base L1 = 1/160 | 1/1280 | measured |
| `2a-sphere-uniform` | uniform | Mach 2 shock over a sphere (r 0.13), immersed boundary, Euler wall | 3-D, octree | WENO-5 characteristic, SSP-33 | uniform L4, 24^3 blocks | 1/384 | 56.6 M |
| `2b-sphere-runtime` | runtime | same | same | same | density gradient marker, tol 0.1 drho/h(base), regrid every 10 steps, base L2 | 1/384 | measured |
| `2c-sphere-fixed` | frozen | same | same | same | 2b's grid frozen | 1/384 | measured |
| `3a-channel-uniform` | uniform | CKM channel, M_b 1.5, Re_b 3000, isothermal walls | 3-D, octree, periodic x, z | `hybrid-kep-weno`, viscous, SSP-54 | uniform L3, 32x24x16 blocks (the T3 `resolved-hybrid` dump re-blocked) | 4pi/256 x 2/192 x (4pi/3)/128 | 6.29 M |
| `3b-channel-runtime` | runtime | same | same | same | `r = h/eta` marker (mode 7, refine 9.5, derefine 4.4: phase 7.2, a priori), regrid every 10 steps, base L2 (predicted: the core coarsens to L2, notes) | as 3a | measured |
| `3c-channel-fixed` | frozen | same | same | same | 3b's grid frozen | as 3a | measured |

Notes on the table, against the plan's first version (the re-analysis of 2026-11-22 on #133):

- **Leg 1 is WENO-7.** WENO-9 failed the feasibility gate 4.1 (non-finite density at t = 3.62 just behind the step
  corner, uniform h = 1/160, CPU, medooza job st-J1 1369311; no positivity limiter on multi-realm runs). WENO-7 is the
  plan's fallback and its gate is pending; the order is the `weno` key of `[leg.1]`. The warm-up ends at t_w = 2.0,
  well before that failure time.
- **Leg 1's marker is make_step.py's Loehner density marker** (`marker = lohner`), not a plain density gradient: the
  aligned layout of `verification/step/make_step.py` ships it, its thresholds are scale-free, and the plan gives no
  gradient tolerance for the step. `marker = gradient` with `tol_factor` switches it.
- **Leg 1's blocks**: one root block per realm and square cells fix the block shapes, A (3m, m), B (3m, 4m),
  C (12m, 4m) cells with m = 16 (`make_step.py --layout aligned`); `auto` seams need the same block cell counts along a
  seam. Every realm gets an equal share of the memory, so 1a at 1/1280 fits 4 A30 by only 6-21 % (re-analysis).
- **Leg 2 is generated afresh** (`make_sphere.py`), with the immersed-boundary sphere kind (`analytical_sphere`); the
  removed matrix's input survives as `templates/sphere-legacy.ini`, a record only. The rectangle kind is not used (its
  distance function is wrong, #143).
- **Leg 3 runs on the T3 dump re-blocked** (L3, 32x24x16 blocks, base L2), the plan's first layout: the dump's own
  layout (L2, 64x48x32 blocks) is the one `load_blocks` would read as it is, but **on it the runtime case cannot
  coarsen the core** (phase 7.2): an L1 block is a half
  channel, wall to centreline, and a sibling group derefines only when all 8 children vote derefine, so the 4 wall
  blocks of every L1 parent pin it (or, with derefine_tol above the wall's r_max, everything coarsens to a uniform
  L1). The marker itself separates the two regions: `kolmogorov_apriori.py` on the settle dump (step 55626) and the
  T3 dump (step 51626), Re_tau 202, h+ 5.2, gives block r_max 5.3-9.4 on the wall blocks (|y| > 0.5) and 2.0-4.1
  in the core; one level coarser (restricted field, 2h) 10.2-14.4 and 3.8-7.2. The thresholds in `[leg.3]`,
  `refine_tol` 9.5 and `derefine_tol` 4.4, derefine exactly the core and keep it coarse with blocks of 32x24x16
  at L3 (`max_level` 3, `base_level` 2): 256 L3 + 32 L2 blocks, 3.54 M cells, 56 % of 3a. That layout is the T3
  dump re-blocked by `reblock.py` (each 64x48x32 block into its 8 children, "Re-blocking the T3 dump" below) and
  doubles the ghost share of a block (2.04 vs 1.46 cells read per cell).

**Roles.** The plan's letters differ between legs (leg 1: b frozen, c runtime; legs 2 and 3: b runtime, c frozen), so
the harness names the cases by role: **U** uniform, **R** runtime, **F** frozen (`[amr] frequency = 0` on the runtime
grid, D2).

## Protocol

**The chain of a leg** (D3), `run.sh` stages:

1. **warmup**, once per leg, on the fastest backend: the uniform case from the cold start to t_w (leg 1 2.0, leg 2
   0.15, leg 3 150), restart dump in `<work>/dumps/leg<N>/warmup/`. Leg 3: the T3 `resolved-hybrid` dump, re-blocked
   (`reblock.py`, below), replaces the run (`run.sh warmup --cases 3 --use-dump DIR` links its files). Then run
   `settle` and `windows`, not `chain`: `chain` runs its own warm-up, and every stage moves an existing directory to `<work>/attic/` first, so a
   chain after `--use-dump` throws the linked dump away and warms up from the cold start (to t_w = 150).
2. **settle**, once per leg: the runtime case restarted from the warm-up dump (repartitioned for the settle's rank
   count, rule `units`), `settle_steps` steps (legs 1, 2: 100; leg 3: 4000, one eddy turnover); the first regrids
   derefine what the marker does not need; dump in `<work>/dumps/leg<N>/settle/`.
3. **windows**, every backend and rank count: U from the warm-up dump, R and F from the settle dump, each dump
   repartitioned for the reading rank count by `repartition.py` (rule `siblings` for the static grids U and F, the
   layout of a static cold start; `units` for R, the layout of a regrid).

**A window**: three times, interleaved over the cases, a run of 10 and a run of 110 steps from the window dump, every
output off (every `*_save` key of `[IO]` and `[diagnostics]` set to 0), `ADAM_SEAM_TIMING=1`. 10 and 110 are multiples of
the regrid period, so the two runs differ by exactly ten regrids.

**What `analyze.py` computes**, per backend, rank count and case:

- step cost = `median(T_loop(110) - T_loop(10)) / 100` from the forest timing line (no start-up, no restart read; the
  difference removes the first-step transients); the same difference of the steps, regrid and seam-fill timers splits
  it;
- spread = `(max - min) / median` of the three differences, flagged above 2 %;
- N = the mean cell count over the 100 steps the difference keeps: the blocks per realm of the window dump times the
  cells per block, updated by the `flume: regrid at step S: blocks A -> B` lines (a regrid at S acts from step S + 1;
  the lines do not name their realm, they come in realm order, and the realm of each is the one whose count matches;
  when realms of equal counts leave it ambiguous the summary says so);
- throughput = `step cost x ranks / N`, in ns per cell-step (ranks = devices): the grind time g of phase 9, with the
  cells per device `N / ranks` beside it, the strong and weak efficiencies ("Multi-node and weak scaling (phase 9)").

**Ratios** (per leg, backend and rank count; backend speedup per case):

| Ratio | Measures |
|---|---|
| T_cpu-node / T_backend | backend speedup, per case (the reference: the CPU run with the fewest ranks, one node; the multi-node CPU runs of phase 9 are rows like any other) |
| T(U) / T(R) | AMR speedup at equal finest h; ideal bound N(U) / N(R) |
| [T(U) N(F) / N(U)] / T(F) | AMR grid efficiency at constant cell work (ghosts, 2:1 seam fills, reflux) |
| T(R) / T(F) - 1 | regrid cost (marking, regrid, prolongation; the FNL host round trip) |
| step cost x devices / N | throughput, ns per cell-step |

The hypotheses H1-H4 these ratios test are in the issue.

## Results (medooza, 2026-11-23 and 2026-11-24)

The tables of `analyze.py` are in `results/2026-11-23/summary.{md,csv}` (the matrix and its multi-node runs) and
`results/2026-11-24/` (the grid study, the rank placement and the CPU node at two block sizes). The executables:
`build.sh` on the night tree of 2026-11-22 (develop `0bcdf870` and the night's steps, since committed in
`6f5e57a1`..`b06f02d3`: `src_modified=yes` in every `run.info`), `flume-cpu-gnu` (`iac_gnu`, gfortran 16.2.0 + Open
MPI 5.0.11) and `flume-fnl-nvf` (`iac_nvf`, nvfortran 26.9, CUDA 12.9, LTO, not fastlink). A CPU node is 64 ranks, a
GPU node 4 A30; step cost = median over 3 repetitions of (T_loop(110) - T_loop(10)) / 100, spreads below 1 % on FNL,
up to 7 % on the CPU node.

**Scope** (owner, 2026-11-24): scaling is measured and claimed on the **uniform 3-D grids** only (2a, 3a and the grid
study's `wg-*`, the uniform sphere on other grids). Runtime AMR (the host regrid, 2:1 seams, reflux, the runtime
layout) and the 2-D quadtree path (ghost layers along the null axis) are implemented for accuracy and not yet optimised
for performance; their timings are kept as baselines ("The accuracy-first paths" below), with no scaling claim.

**A GPU node against a CPU node, uniform grids** (T(cpu, 64 ranks) / T(fnl-nvf, 4 A30)):

| Case | Cells | Blocks | CPU node (s) | GPU node (s) | GPU node / CPU node |
|---|---|---|---|---|---|
| 2a sphere (3-D, WENO-5, IB, SSP-33) | 56.6 M | 4096 of 24^3 | 12.7 | 0.780 | 16.3 |
| 3a channel (3-D, hybrid-kep-weno, viscous, SSP-54) | 6.29 M | 512 of 32x24x16 | 3.17 | 0.230 | 13.8 |
| `wg-b24` sphere, L3 (2026-11-24) | 7.08 M | 512 of 24^3 | 1.366 | 0.1147 | 11.9 |
| `wg-b38` sphere, L3 (2026-11-24) | 28.1 M | 512 of 38^3 | 4.458 | 0.355 | 12.6 |

The factor is not one number: on the same flow it is 11.9-16.3 with the grid. On the CPU node 2a costs 14.4 us per
cell-step per core, against 12.4 for `wg-b24` (8 blocks per rank) and 10.2 for `wg-b38` (job 1369687; spreads 5.3 and
0.7 %). The 2a figure is inflated by subnormal arithmetic in its developed state, which slows the CPU ranks that hold it
and not the A30 ("Where a step's time goes" below): with flush to zero 2a costs 8.74 us per cell-step per core, the cold
`wg-b24` 8.81-9.54. (This paragraph first blamed the immersed boundary's imbalance: wrong, 2026-11-24.) On the GPU node
the two grid-study grids are at 1.8 and 7.0 M cells per A30, below 2a's 14.2 M.

**Strong scaling, uniform sphere** (2a, 56.6 M cells, 24^3 blocks; step cost in s, efficiency g(1 node) / g):

| | 1 node | 2 nodes | 4 nodes | 8 nodes |
|---|---|---|---|---|
| fnl-nvf (4 A30 per node) | 0.780 | 0.491 (0.79) | 0.260 (0.75) | 0.167 (0.58) |
| cells per A30 | 14.2 M | 7.08 M | 3.54 M | 1.77 M |
| cpu (64 ranks per node) | 12.7 | 6.48 (0.98) | 3.24 (0.98) | 2.31 (0.69; job 1369694) |

**Weak scaling, uniform sphere** (g = step cost x devices / N, ns per cell-step; 24^3 blocks). The weak efficiency of
the 2026-11-23 version mixed two costs: the cross-pair exchange inside a node (1 -> 4 A30) and the network (one node
-> several). Separated with the grid study ("Grid size" below):

| Cells per A30 | Blocks | One node of 4 A30 | Several nodes | Inter-node efficiency |
|---|---|---|---|---|
| 14.2 M | 24^3 | 55.1 (2a) | 69.3 (`w2-sphere-l5`, 453 M cells, 8 nodes) | **0.80**, measured |
| 7.02 M | 38^3 | 50.5 (`wg-b38`) | 67.0 (`wg-l4b38`, 225 M cells, 8 nodes; job 1369691) | **0.75**, measured |
| 7.08 M | 24^3 | ~58.7, not run (`wg-b38` x 51.0/43.9, the 1-A30 block factor) | 69.4 (2a, 2 nodes) | ~0.85, estimate |
| 7.02 M, no body | 38^3 | 37.1 (`wg-b38`, the sphere dropped) | 41.3 (`wg-l4b38`, the sphere dropped; jobs 1369710, 1369711) | **0.90**, measured |

**Without the immersed body the 8-node pair keeps 0.90** (the same grids and dumps, `[solids]` dropped from each
rank's input by `speedup-133/noib-fnl.sh`, a `run.sh --exe` wrapper; `results/2026-11-24/gpu-weak-noib/`). The cost of
leaving the node is the same per ghost exchange in both runs, 9.7-9.8 ms ((67.0 - 50.5) ns x 7.02 M / 12 exchanges per
step with the body, its eikonal extension exchanging the ghosts 4 times per stage; (41.3 - 37.1) x 7.02 M / 3 without):
on the GPU the multi-node loss follows the number of ghost exchanges. The body costs 36 % on one node (50.5 against
37.1 ns: the 9 extra exchanges and the eikonal kernels).

The 2026-11-23 weak efficiencies at 7.08 and 3.54 M cells per A30 (0.73, 0.72) took 1- and 2-A30 references, so
they held the cross-pair cost of the fourth A30 (13-22 % at those loads) as well as the network's. The 8-node pair at
38^3 blocks loses more than the one at 24^3 (0.75 against 0.80) although its blocks carry fewer ghost cells: a 38^3
step is 15 % cheaper per cell, and at half the load per A30 the same exchange weighs more (16.5 ns per cell-step on
top of one node's, against 14.2). In throughput the two 8-node runs are equal: 225 M cells at 0.470 s per step and
453 M at 0.981 s, 478 and 462 M cell-steps per second.

**Findings (uniform grids)**

- **On one A30 the block size sets the cost of a cell; the load matters little above a few M cells.** The same 192^3
  grid costs 66.0, 51.0 and 44.0 ns per cell-step in 12^3, 24^3 and 48^3 blocks; the matrix's 24^3 leaves ~16 % on
  the table. An A30 holds 28.1 M cells of the sphere (38^3 blocks).
- **Inside a node the cost is the fourth GPU, not the second**: 2 A30 cost 2-9 % more per cell than one (an NVLink
  pair), 4 A30 another 13-30 % (the exchange between the pairs). Pinning ranks to their GPU's NUMA node changes
  nothing (within 0.8 %).
- **Across nodes the uniform sphere keeps 0.75-0.80 of one node's throughput per A30 on 8 nodes** at 7-14 M cells per
  A30 with the immersed body, 0.90 without it: the loss is the MPI ghost exchange, GPU halos crossing nodes at 1.5 GB/s
  without GPUDirect RDMA ("Where a step's time goes"). Strong scaling of 2a from 1 to 8 nodes (0.79, 0.75, 0.58) runs
  down the load curve: at 8 nodes 2a has 1.77 M cells per A30 in 24^3 blocks, where one node alone already costs 18 %
  more per cell than at 14.2 M.
- **The CPU node scales to 4 nodes (2a 0.98 on 2 and 4) and, as measured, loses half on 8** (weak 0.53): subnormal
  arithmetic on 76 ranks, not the network (with flush to zero 0.93-1.0, "Where a step's time goes"); a GPU node is
  12-16 times faster than a CPU node
  on the uniform 3-D flows (11.9 and 12.6 on the grid-study grids, 16.3 on 2a, 13.8 on the channel).
- **The per-cell cost does not depend on the flow state** (0.5 %, "State independence" below): cold-started cases and
  developed ones are comparable.


**Figures** (`matrix_plot.py results/2026-11-23/summary.csv --out results/2026-11-23 --theme both`, phase 8.1;
`grid_study_plot.py`, 2026-11-24; copies in `docs/public/flume/perf-*.png` for the docs page
[Performance](../../../../docs/applications/flume/performance.md)): `perf-grid-study.png` (g against the block edge
on 1 A30, and the same grids on 1, 2 and 4 A30 with the rank-placement A/B), `perf-backend-speedup.png` (T(CPU node) /
T(GPU node) per case), `perf-amr-decomposition.png` (per leg and node, the step cost of U, F and R over U's, split
inside and outside the steps, with the ideal N/N(U)), `perf-grind.png` (g against the cells per device, every case);
`grind.png` is `scaling_plot.py`'s earlier figure.

### Grid size on the GPU node (2026-11-24)

Owner, 2026-11-24: some findings looked like the grids were too small on 8 and 16 A30; before anything else, do the
grid dimensions affect the results? The grid study (`[weak.wg-*]` in the study's conf, `matrix.conf` plus the cases
below; work tree `/ifs/hpc/home/szaghi/adam-tests/speedup-133/grid-study`, jobs 1369668-1369678) times the uniform
sphere of leg 2 (WENO-5, IB, SSP-33) from 20-step cold starts, the protocol's windows, on 1, 2 and 4 A30 of one node
with the matrix's launch. `grid_study_plot.py` draws `results/2026-11-24/perf-grid-study.png`.

**Block size** (1 A30; L3 = 512 blocks of b^3 cells; g in ns per cell-step):

| b | 16 | 20 | 24 | 28 | 32 | 36 | 38 | 40 |
|---|---|---|---|---|---|---|---|---|
| cells | 2.10 M | 4.10 M | 7.08 M | 11.2 M | 16.8 M | 23.9 M | 28.1 M | 32.8 M |
| g, 1 A30 | 63.1 | 54.4 | 51.0 | 47.7 | 46.2 | 44.4 | 43.9 | refused |

The same 192^3 grid cut three ways on 1 A30: 4096 blocks of 12^3 **66.0**, 512 of 24^3 **51.0**, 64 of 48^3 **44.0**.
At a fixed grid the block size alone moves g by 50 %; at a fixed block size (b 24) 7.1 M cells cost what the 48^3 cut
of the same grid costs at b 48, so above a few M cells per A30 the load matters little and the block size a lot. Over
the L3 series g = 15.7 + 18.0 ((b+6)/b)^3 ns fits every point within 0.9 ns (a cost per interior cell plus one per
stored cell, ngc 3); it is a description of the series, not a law of b: the 4096-block cut lies 23 ns below its
extrapolation to b 12. **Every matrix GPU number uses 24^3 blocks: ~16 % of its step cost is avoidable with 38^3-48^3
blocks** (51.0 against 43.9-44.0 on one A30).

**The memory ceiling of an A30 is 28.1 M cells** of this case: b 40 (32.8 M) is refused by the blocks budget (0.85 x
the 25.23 GB the A30 reports, 23.5 GiB, holds 510 blocks of 46^3 stored cells at 54 fields, 512 needed); b 38 fits
at 88 % of it (583 blocks of 44^3).

**GPUs inside a node** (the same grids on 2 and 4 A30: at one b the work is the same, the ratio is the GPUs'):

| b | 16 | 20 | 24 | 28 | 32 | 36 | 38 |
|---|---|---|---|---|---|---|---|
| g, 2 A30 | 68.8 | 57.1 | 53.0 | 49.2 | 47.0 | 45.4 | 44.7 |
| g, 4 A30 | 89.5 | 71.2 | 64.8 | 58.2 | 55.7 | 52.4 | 50.5 |
| g(2) / g(1) | 1.09 | 1.05 | 1.04 | 1.03 | 1.02 | 1.02 | 1.02 |
| g(4) / g(2) | 1.30 | 1.25 | 1.22 | 1.18 | 1.19 | 1.16 | 1.13 |

Spreads are below 1 %. Going from 1 to 2 A30 costs 2-9 %, from 2 to 4 A30 13-30 %, at every load. The node's
topology (`nvidia-smi topo -m`, job 1369679) pairs the A30 by NVLink (GPU0-GPU1, GPU2-GPU3, NV4, 4 x 25 GB/s); between
the pairs, and from every GPU to the InfiniBand card, the path is SYS (PCIe and the inter-socket link). The GPUs sit
on NUMA nodes 0, 1, 3 and 3 (cores 0-15, 16-31, 48-63). The 4-A30 penalty is the first cross-pair exchange; "Rank
placement" below measures how much of it the launch adds.

### Rank placement in a GPU node (2026-11-24)

The launch of the matrix (`srun --mpi=pmix -n 4`, `--cpus-per-task=8`, no binding asked) lets Slurm place the ranks:
ranks 0 and 2 run on cores 0-15 (NUMA 0), ranks 1 and 3 on cores 32-47 (NUMA 2, which has no GPU), and every rank sees
the four GPUs and picks `mod(local rank, 4)`. So only rank 0 runs next to its A30. `jobs/medooza-bind-fnl.sh` (passed to
`run.sh --exe`, launched with `srun --mpi=pmix -c 16 --cpu-bind=none -n`: the default step's cpuset holds cores 0-15 and
32-47 only) pins rank l to GPU l (`CUDA_VISIBLE_DEVICES=l`), 8 cores of that GPU's NUMA node and its memory (`numactl
--physcpubind --membind`; cores 0-7, 16-23, 48-55, 56-63). Jobs 1369688 (one node) and 1369689 (two nodes),
`--exclusive`, ran each case unbound then bound on the same nodes (work trees `bind/{off,on}-np{4,8}`):

| Case | Cells per A30 | 4 A30 unbound | 4 A30 bound | 8 A30 unbound | 8 A30 bound |
|---|---|---|---|---|---|
| 2a (L4, 24^3) | 14.2 M / 7.08 M | 55.4 | 55.2 | 69.5 | 69.3 |
| `wg-b38` (L3, 38^3) | 7.02 M / 3.51 M | 50.6 | 50.6 | 67.2 | 67.1 |
| `wg-b24` (L3, 24^3) | 1.77 M / 0.88 M | 64.7 | 64.3 | 99.6 | 98.8 |

(g in ns per cell-step; spreads at most 0.9 %.) The unbound runs reproduce the grid study and the matrix within 0.5 %
(2a 55.1 and 69.4, `wg-b38` 50.5, `wg-b24` 64.8).

**Pinning changes nothing measurable**: -0.4 to 0 % on one node, -0.8 to -0.1 % on two, inside the spreads. The
FNL halo exchange posts MPI on device buffers (GPU-direct by design, no host staging), which is consistent with the
host core being off its path (the transports UCX picks per link were not inspected); the 2 -> 4 A30 cost of "Grid
size" is then the cross-pair link, not the launch. Binding is therefore not made the default;
`jobs/medooza-bind-fnl.sh` is kept for the memory headroom and for other nodes.

**The placement explains GPU 0's extra memory** ("Device memory" below): unbound, rank 0 reports 23.99 GB of device
memory available at start and ranks 1-3 24.91 GB, because every rank touches GPU 0 before its own device (the uploads
the memory study saw); bound, every rank sees one GPU and all four report 24.91 GB. The blocks budget is computed from
the device's total memory, so the capacity does not change; the 0.9 GB is headroom on GPU 0.

### Where a step's time goes, and what does not scale (2026-11-24)

Owner, 2026-11-24: "I already suspect that the ghost exchange is the current bottleneck of ADAM: measure all relevant
parts (flux computation, ghost exchange, time integrator, etc.), rank them, find which is not scaling and why (maybe
not only the MPI barrier)."

**Method.** A measurement-only patch, not committed (`results/2026-11-24/phase/phase-timer.patch`, applied to HEAD
`a88d8367` by `apply_phase_timer.py`; the module is `adam_phase_timer.F90.txt`), adds cumulative `MPI_Wtime` timers to
22 phases of the FLUME step on both backends: dt (local and its `MPI_ALLREDUCE`), the immersed-boundary kernels, the
local ghost fill, the MPI ghost exchange split into pack / post / `MPI_WAITALL` / barrier 1 / unpack / barrier 2,
boundary conditions and the ghost blend, `q_aux`, face fluxes, flux difference, the integrator (SELAH), the forest's
post step, regrid and `is_done` (two `MPI_ALLREDUCE`), progress printing. FNL kernels are synchronous (no `async`
clause anywhere), so a host timer measures their device time. `ADAM_PHASE_TIMING=1` arms it, the first 10 steps are
not timed, the report gives min / mean / max over the ranks in ms per step; `ADAM_GHOST_BARRIER=0` skips the two
barriers of every MPI ghost exchange, `ADAM_FTZ=1` sets abrupt underflow (flush to zero) in the forest's time loop.
Runs: 110 steps from the window dumps of the weak pairs, CPU 110,592 cells per core (24^3 blocks; `wg-b24` on 1 node,
2a on 8), FNL 7.02 M cells per A30 (38^3; `wg-b38` on 1 node, `wg-l4b38` on 8), each with the immersed body and
without it, barriers on and off (jobs 1369712-1369714, 1369718, 1369721-1369723; reports in
`results/2026-11-24/phase/reports/`, ranked by `phase_compare.py`). Every barriers-on/off pair has a bitwise-identical
residual history.

**GPU: the MPI ghost exchange is the whole multi-node loss** (ms per step, means over the ranks):

| Phase | 1 node, no body | 8 nodes, no body | 1 node, body | 8 nodes, body |
|---|---|---|---|---|
| step | 260.5 | 289.0 | 355.3 | 467.0 |
| face fluxes | 208.5 | 208.0 | 207.4 | 208.4 |
| **ghost MPI** | 12.7 | **43.7** | 50.0 | **172.2** |
| ... `MPI_WAITALL` | 10.3 | 36.0 | 40.7 | 144.2 |
| ... both barriers | 0.3 | 3.8 | 0.8 | 12.5 |
| ghost local | 9.4 | 9.4 | 37.4 | 37.4 |
| BC + blend | 7.3 | 5.0 | 29.3 | 20.1 |
| IB kernels | - | - | 5.8 | 3.4 |
| flux diff, q_aux, integrator, dt, post step, reductions, progress | 22.6 | 22.9 | 25.4 | 25.6 |

From 1 to 8 nodes the ghost MPI exchange grows by 31.0 ms (no body) and 122.2 ms (body): 109 % of the step's growth
both times, 90-93 % of it in `MPI_WAITALL`; every compute phase is flat (face fluxes within 0.5 %). Skipping the
barriers saves 0.4 % (no body) and 1.2 % (body) on 8 nodes. **Why: the transport.** One HDR InfiniBand card per node
(`mlx5_0`, 200 Gb/s), reached from every A30 over SYS, and no GPUDirect RDMA (neither `nvidia_peermem` nor `gdrdrv` is
loaded; driver 545.23; the HPC-X UCX 1.21 of the FNL toolchain sees only `rc_verbs`/`ud_verbs` on the card): UCX moves
GPU buffers between nodes through host staging. `ucx_perftest`, 33 MB messages (jobs 1369716, 1369717, 1369719;
`results/2026-11-24/phase/network.txt`):

| Path | GB/s |
|---|---|
| between nodes, host to host | 23.5 |
| **between nodes, GPU to GPU** | **1.53** (1.75 at 8 MB; best UCX setting, 4 MB rendezvous fragments, 1.89) |
| in a node, GPU 0 to GPU 1 (NVLink pair) | 88.5 |
| in a node, GPU 0 to GPU 2 (across the pairs) | 17.7 |

With 38^3 blocks on 8 nodes each node holds an octant of the domain and shares 3 faces with other nodes: ~33 MB of
halo each way per exchange (an estimate from the layout), 12 ms at the measured waitall, ~2.8 GB/s effective. The
across-the-pairs 17.7 GB/s is the cost of the third and fourth A30 of a node ("Grid size" above). **The immersed body
multiplies it**: its eikonal extension refreshes every ghost of every variable 4 times per stage (12 exchanges per
step instead of 3), so on one node the body costs 94.8 ms per step, of which the ghost work is 87.4 (MPI +37.3, local
+28.0, BC +22.0) and its own kernels 5.8; on 8 nodes the ghost update (local + MPI + BC) is 230 ms, 49 % of the step,
more than the face fluxes.

**CPU: subnormal numbers, not the exchange.** On 8 CPU nodes the face fluxes, the same work on every rank, cost
838 / 980 / 1904 ms per step (min / mean / max): the slowest rank takes twice as long, and the others wait for it in
barrier 1 (425 ms), `MPI_WAITALL` (199) and the step's allreduces (310). Skipping the barriers moves the wait into the
allreduces (2067 against 2009 ms); pinning every rank to a core changes nothing (1987 against 1990; the default launch
binds each rank to a 16-core NUMA domain); another set of nodes gives the same straggler. The developed 2a state holds
745,771 subnormal doubles (below 2.2e-308) in its interior cells, all in rho v and rho w (transverse momenta decaying
ahead of the shock), on 76 of the 512 ranks (up to 11,720 each, ~2 % of a rank's values); the cold-started grid of the
1-node reference holds none (`subnormal_count.py`). x86 cores handle subnormal operands in microcode, many times slower;
the A30 at full speed. **With flush to zero** (`ADAM_FTZ=1`, the same dump and binary):

| CPU, no body | gradual underflow | flush to zero |
|---|---|---|
| 2a, 8 nodes, 512 ranks: step (face fluxes min / mean / max) | 1988 ms (836 / 976 / 1890) | **1044 ms** (836 / 848 / 962) |
| ... `MPI_WAITALL`, barrier 1, regrid + is_done | 192, 421, 304 | 18, 46, 33 |
| 2a, 1 node, 64 ranks: step (face fluxes min / mean / max) | 11.83 s (6715 / 7860 / 11143) | **7.73 s** (6755 / 6847 / 6986) |

The CPU weak pair becomes 0.93-1.0 instead of 0.51: 1044 ms on 8 nodes against, on one node, `wg-b24` (no subnormal
value) at 974 ms with gradual underflow and 1055 ms with flush to zero (job 1369722, another node: the CPU node's
spread; flush to zero changes nothing where there is no subnormal); two grids at the same load, so an estimate. On one
node 2a costs 8.74 us per cell-step per core with flush to zero, no more than the cold grid (8.81 us gradual, 9.54 with
flush to zero): the "14.4 against 11.2 us" of "9.2.4" was this. The CPU numbers of the matrix on 2a (and on every
developed state) carry the subnormal slowdown, so the GPU node / CPU node factor 16.3 on 2a is inflated; the grid-study
factors (11.9, 12.6, cold grids) are the clean ones. Skipping the barriers on the CPU brings nothing (one node: 1150
against 974 ms without the body, 1367 against 1266 with it; single runs, inside the CPU node's spread or worse).

**Ranked, with the measured or estimated gain** (proposals, none applied):

| # | Cost | Where | Gain | Proposal |
|---|---|---|---|---|
| 1 | subnormal arithmetic | CPU, developed states | -35 % (1 node) to -47 % (8 nodes) on 2a | flush to zero for the time loop (or `-ffast-math`'s crtfastmath at link), or floor values below a tiny threshold; changes results in the last bits: the regression goldens would need a refresh |
| 2 | GPU halo between nodes at 1.5 GB/s | FNL, multi-node | up to ~9 ms of the ~12 ms per exchange (estimate at 23.5 GB/s host path plus the two device copies) | stage the inter-node halo through pinned host buffers; or GPUDirect RDMA on medooza (`nvidia_peermem`, admins) |
| 3 | 12 ghost exchanges per step with the immersed body | both backends | up to 9 of the 12 exchanges per step | the eikonal iterations exchange only the variables and blocks they need, or iterate inside a wider ghost band between exchanges |
| 4 | the two `MPI_Barrier` per ghost exchange | both | 0.4-1.2 % (GPU), none (CPU) | remove (bitwise identical in every pair) |
| 5 | across-the-pairs link in a node (17.7 GB/s) | FNL, 4 A30 | the 13-30 % of "Grid size" | rank order that keeps neighbours in an NVLink pair |

### The accuracy-first paths: runtime AMR and 2-D (baselines)

Measured with the matrix; recorded as the baselines of the optimisation work, not as scaling results.

**A GPU node against a CPU node, every case** (T(cpu, 64 ranks) / T(fnl-nvf, 4 A30)):

| Leg | Uniform | Frozen AMR grid | Runtime AMR |
|---|---|---|---|
| 1 step (2-D, WENO-7, SSP-54) | 5.83 | 5.75 | 5.72 |
| 2 sphere (3-D, WENO-5, IB, SSP-33) | 16.3 | 12.2 | 10.9 |
| 3 channel (3-D, hybrid-kep-weno, viscous, SSP-54) | 13.8 | 15.5 | 16.7 |

**Runtime AMR** (AMR speedup T(U)/T(R) against its ideal N(U)/N(R); grid efficiency [T(U) N(F)/N(U)] / T(F); regrid
cost T(R)/T(F) - 1):

| Leg | Backend | Devices | AMR speedup | Ideal | Grid efficiency | Regrid cost |
|---|---|---|---|---|---|---|
| 1 | fnl-nvf | 4 | 1.04 | 1.59 | 0.89 | 34 % |
| 1 | cpu | 64 | 1.06 | 1.59 | 0.90 | 34 % |
| 2 | fnl-nvf | 4 | 1.48 | 4.02 | 0.51 | 38 % |
| 2 | cpu | 64 | 2.20 | 4.02 | 0.68 | 24 % |
| 3 | fnl-nvf | 4 | 1.25 | 1.78 | 0.90 | 29 % |
| 3 | cpu | 64 | 1.03 | 1.78 | 0.81 | 39 % |

**Step costs on more nodes** (s per step; efficiency against one node):

| Case | fnl-nvf 4 / 8 / 16 / 32 A30 | Efficiency 8 / 16 / 32 | cpu 64 / 128 / 256 ranks | Efficiency 128 / 256 |
|---|---|---|---|---|
| 1a step uniform (2-D) | 0.216 / 0.119 / 0.070 / - | 0.91 / 0.77 / - | 1.26 / 0.642 / 0.328 | 0.98 / 0.96 |
| 1b step frozen | 0.154 / 0.089 / 0.057 / - | 0.87 / 0.67 / - | 0.887 / 0.544 / 0.334 | 0.81 / 0.66 |
| 1c step runtime | 0.207 / 0.119 / 0.076 / - | 0.87 / 0.69 / - | 1.19 / 1.18 / 1.18 | 0.50 / 0.25 |
| 2b sphere runtime | 0.526 / 0.350 / 0.312 / 0.188 | 0.75 / 0.42 / 0.35 | 5.76 / 3.95 / 3.06 | 0.73 / 0.47 |
| 2c sphere frozen | 0.381 / 0.284 / 0.257 / 0.161 | 0.67 / 0.37 / 0.30 | 4.64 / 3.13 / 1.97 | 0.74 / 0.59 |

The 2-D step's weak series (1.03 M cells per A30, memory-capped by the null-axis ghosts): g 202 (1 A30), 209 (4),
212 (8, 2 nodes), 214 (16, 4 nodes).

- **The 2-D cell is the dearer one**: the 2-D step at 1.03 M cells per A30 spends 202 ns per cell-step, the same forest
  extruded in 3-D at the same order and load 126 (H2): the boundary-condition kernel refills the null-axis ghost
  layers at every stage. A GPU node runs the step 5.8 times faster than a CPU node, against 14-16 on the 3-D uniform
  grids.
- **Runtime AMR buys little on either node: 1.04 (step), 1.48 (sphere), 1.25 (channel) on 4 A30; 1.06, 2.20, 1.03 on
  64 cores**, against ideals of 1.59, 4.02 and 1.78. On the GPU node the losses are the AMR grid itself (grid
  efficiency 0.89, 0.51, 0.90) and the regrid, which runs on the host (22-27 % of the runtime step); on the CPU node
  the AMR grid (0.90, 0.68, 0.81) and the **rank layout of the runtime grid**, which leaves the busiest rank with the
  uniform case's load on legs 1 and 3; the regrid itself is 1-2 % there ("Where the AMR overhead goes" below).
- **The CPU's runtime-AMR step does not scale with ranks** (1c 1.18 s per step on 64, 128 and 256 ranks) **because of
  its layout, not its regrid**: a runtime regrid distributes whole sibling sets and splits every realm from rank 0, so
  1c's busiest rank holds 64,512 cells at every rank count.
- **The sphere's AMR grids stop scaling below about 1 M cells per A30**: 2c costs 0.26 s per step on 16 A30 at 0.88 M
  cells per A30, the same as 2a with four times the cells.
- **The step is not launch-bound at 1 M cells per A30** (the plan's reading of H1, refuted): it scales to 0.26 M cells
  per A30 (0.77).
- **Leg 3** needed the T3 dump re-blocked (`reblock.py`) for the Kolmogorov marker to coarsen the core; the a priori
  prediction (288 blocks) held exactly at the first regrid. Its margins (~1.08) rest on two snapshots.

### Results per leg (phases 5.4, 6.4, 7.6)

The case directories (`make_cases.py`, `<work>/cases/<case>/README.md`) are regenerated by every `run.sh` call; their
measured results are here. Step costs in s per step (median of 3; spreads in `summary.md`), N the window's mean cell
count; the dumps and logs are on medooza under `/ifs/hpc/home/szaghi/adam-tests/speedup-133/work`.

**Leg 1, Mach 3 step** (2-D, three realms, WENO-7, SSP-54). Warm-up on 4 A30 to t_w = 2.0 (34,672 steps, job 1369400;
its last checkpoint shows the bow shock, the Mach stem, the slip line with Kelvin-Helmholtz roll-ups and the reflected
shocks, nothing at the `auto` seams), settle 100 steps; the window's runtime grid changed once in ten regrids
(realm C 187 -> 184 blocks).

| Case | N | Blocks (A, B, C) | 4 / 8 / 16 A30 | 64 / 128 / 256 cores | GPU node / CPU node |
|---|---|---|---|---|---|
| 1a uniform | 4.13 M | 256, 256, 256 | 0.216 / 0.119 / 0.0704 | 1.26 / 0.642 / 0.328 | 5.83 |
| 1b frozen | 2.60 M | 106, 73, 187 | 0.154 / 0.0890 / 0.0573 | 0.887 / 0.544 / 0.334 | 5.75 |
| 1c runtime | 2.60 M | 106, 73, 187 | 0.207 / 0.119 / 0.0757 | 1.19 / 1.18 / 1.18 | 5.72 |

AMR speedup 1.04 (GPU node), 1.06 (CPU node), ideal 1.59; grid efficiency 0.89 / 0.90; T(R)/T(F) - 1 34 % / 34 %
(on the GPU node the regrid; on the CPU node the runtime layout, 1c's steps cost 1.17 s against 1b's 0.88).

**Leg 2, shock over a sphere** (3-D, IB, WENO-5, SSP-33). Warm-up on 4 A30 to t_w = 0.15 (663 steps), settle 100
steps (job 1369401); the window's runtime grid changed once in ten regrids (1016 -> 1072 blocks).

| Case | N | Blocks | 4 / 8 / 16 / 32 A30 | 64 / 128 / 256 cores | GPU node / CPU node |
|---|---|---|---|---|---|
| 2a uniform | 56.6 M | 4096 | 0.780 / 0.491 / 0.260 / 0.167 | 12.7 / 6.48 / 3.24 | 16.3 |
| 2b runtime | 14.1 M | 1016 -> 1072 | 0.526 / 0.350 / 0.312 / 0.188 | 5.76 / 3.95 / 3.06 | 10.9 |
| 2c frozen | 14.0 M | 1016 | 0.381 / 0.284 / 0.257 / 0.161 | 4.64 / 3.13 / 1.97 | 12.2 |

AMR speedup 1.48 (GPU node), 2.20 (CPU node), ideal 4.02; grid efficiency 0.51 / 0.68; T(R)/T(F) - 1 38 % (GPU: the
regrid, 0.144 s per step) / 24 % (CPU: inside the steps, 5.70 against 4.63 s; the runtime layout gives the busiest
rank 1.13x the mean load against 2c's 1.01x, about half of it; the rest is not located; the regrid timer is 1.3 %).

**Leg 3, turbulent channel** (3-D, `hybrid-kep-weno`, viscous, SSP-54; Re_tau 202). The T3 `resolved-hybrid` dump
(t = 150, step 51626) re-blocked to L3 blocks of 32x24x16 (`reblock.py`, job 1369523), settle 4000 steps on 4 A30
(job 1369525: one regrid, 512 -> 288 blocks, the core coarsened to L2 exactly as predicted a priori, then no grid
change); no grid change in the windows.

| Case | N | Blocks | 1 / 4 A30 | 64 cores | GPU node / CPU node |
|---|---|---|---|---|---|
| 3a uniform | 6.29 M | 512 (L3) | 0.800 / 0.230 | 3.17 | 13.8 |
| 3b runtime | 3.54 M | 288 (256 L3 + 32 L2) | - / 0.184 | 3.08 | 16.7 |
| 3c frozen | 3.54 M | 288 | - / 0.143 | 2.22 | 15.5 |

AMR speedup 1.25 (GPU node), 1.03 (CPU node), ideal 1.78; grid efficiency 0.90 / 0.81; T(R)/T(F) - 1 29 % (GPU: the
regrid steps, which change nothing here) / 39 % (CPU: the layout, busiest rank 98,304 cells = 3a's per-rank load).

### Where the AMR overhead goes

**Inside and outside the steps.** `step cost - steps` (both medians of the window differences) is the time outside the
time steps: the regrid. On the GPU node the steps of R cost what the steps of F cost (2b 42.03 s, 2c 41.94 s per 110
steps; 3b and 3c 0.1428 s per step), so R - F is the regrid. The forest's timers are maxima over the ranks, so on the
CPU the regrid timer also holds the wait of the ranks that reach its reduction first (2a's uniform grid shows 1.44 s
per step of "regrid" on 64 ranks: the wait for the slowest ranks, which 2026-11-24 found slowed by subnormal arithmetic,
not by the immersed boundary; "Where a step's time goes").

**The GPU regrid** (forest regrid timer and `flume: regrid at step` lines of the 110-step windows): 2b 16.1 s for 11
regrid steps, one of which changed the grid (logged host round trip 2.03 s), so ~1.4 s per regrid step that changes
nothing; 1c 5.1 s (one change, 2.46 s; ~0.26 s otherwise); 3b 4.5 s with no change (~0.41 s per regrid step). The cost
per regrid step scales with the blocks per GPU (2b: 14.7, 7.0, 4.8, 2.9 s per 100 steps on 4, 8, 16, 32 A30; 5.5-9 ms
per block), not with the changes. Its split among the device-to-host copy, the host marking and the tree work is not
measured (no per-phase timers): draft issue `issue-fnl-regrid-host-cost.md`.

**The CPU runtime layout** (the window dumps' `.tnd` owner ranks; `repartition.py --rule units` reproduces
`tree%mpi_redistribute`, the layout of every runtime regrid):

| Case, ranks | Busiest rank (cells) | Mean | Busiest / mean | Ranks with no block |
|---|---|---|---|---|
| 1b frozen (`siblings`), 64 / 128 / 256 | 44,544 / 28,416 / 16,128 | 40,680 / 20,340 / 10,170 | 1.09 / 1.40 / 1.59 | 0 / 0 / 29 |
| 1c runtime (`units`), 64 / 128 / 256 | 64,512 / 64,512 / 64,512 | 40,680 / 20,340 / 10,170 | 1.59 / 3.17 / 6.34 | 3 / 67 / 195 |
| 3b runtime / 3c frozen, 64 | 98,304 / 61,440 | 55,296 | 1.78 / 1.11 | 0 / 0 |
| 2b runtime, 64 / 128 / 256 | 248,832 / 124,416 / 110,592 | 219,456 / 109,728 / 54,864 | 1.13 / 1.13 / 2.02 | 0 / 0 / 52 |

On FNL the same layout is near balance up to 16 A30 (1c 1.05 on 4, 1.12 on 16; 2b 1.00 on 4, 1.13 on 32). **Measured directly** (job n3m-A): 1b's inputs (`[amr] frequency = 0`) on 1c's window dump (the `units` layout) on 64 ranks cost **1.184 s per step** (spread 1.2 %), against 0.887 s for 1b on its own (`siblings`) layout and 1.187 s for 1c: the layout alone accounts for the runtime step's cost on the CPU node; the regrid adds nothing measurable. Draft issue:
`issue-runtime-layout-sibling-units.md`.

**The sphere's AMR grid on the GPU** (2c: 109 ns per cell-step on 4 A30 at 3.5 M cells per A30, where the
uniform curve is flat at 51-55 ns; 254 blocks per GPU, layout balanced, reflux 2.4 % of the step). Located with
`NVCOMPILER_ACC_TIME=1` (job n3m-A,
10-step runs from the window dumps, device time per GPU, the reports of ranks whose output did not interleave): 2c
(254 blocks per GPU) against `w2-sphere-l3` on 2 A30 (256 uniform blocks per GPU):

| Kernel (per GPU, 10-step run) | 2c, AMR | `w2-sphere-l3`, uniform | 2a, uniform, 1024 blocks |
|---|---|---|---|
| all kernels | 3.57 s | 2.06 s | 7.12 s |
| face fluxes (`compute_face_fluxes_s3_dev`) | 1.09 | 1.11 | 4.36 |
| intra-rank ghost fill (`update_ghost_local_gpu_dev`) | **1.77** | 0.32 | 1.31 |
| ghost send buffers (`populate_send_buffer_ghost_gpu_dev`) | 0.11 | 0.01 | 0.05 |
| boundary conditions | 0.12 | 0.20 | 0.36 |

The ghost fill takes 1.45 s of the 1.51 s difference (96 %). Its 121 launches have the same size on both grids
(111,580 against 105,596 gangs of 128 rows x variables of the map) and take 14.6 ms against 2.6 ms each: on the AMR
grid the map holds the 2:1 rows, coarse-to-fine interpolation (tricubic, `seam_ghost_fill` regime 2) and fine-to-coarse
restriction, which the kernel (`adam_fnl_field_kernels.F90:388`) runs in the same loop as the copies, one branch per
row. How the 5.6x splits among the interpolation arithmetic, the branch divergence and the donor reads is not
measured. Draft issue:
`issue-sphere-amr-gpu-cell-cost.md`.

### Device memory (phase 4.3)

The FNL backend allocates every block-sized device array at the block **capacity** of the budget,
`nb = memory_fraction x device memory / realms / (fields_number x 8 B x block_weight)` (`adam_adam_object`
`compute_blocks_number`), not at the blocks in use; so the device memory a run holds is set by the budget, and the
question of 4.3 is what one stored-cell slot of that capacity really costs against the budget's `fields_number`
(the plan's "~55 doubles per stored cell"). Measured (job n3m-A 1369616, gnode11): the window inputs of 1a, 2a and 3a
on 4 A30 at three values of `ADAM_MEMORY_FRACTION` each, 10 steps, `nvidia-smi` sampled every 0.2 s; the peak per GPU
against the stored-cell slots per rank (sum over the realms of `nb x block_weight`, from the "blocks budget" lines) is
a straight line (residuals below 5 MiB):

| Case | Fractions | Slope: bytes per slot (doubles) | Budget `fields_number` | Intercept per GPU (GPU 0) | Peak at 0.85, GPU 1-3 (GPU 0) |
|---|---|---|---|---|---|
| 1a step, 3 realms, WENO-7, SSP-54 | 0.60, 0.72, 0.85 | 404.9 B (50.6) | 64 | 3.40 GiB (4.27) | 19.21 GiB (20.08) |
| 2a sphere, IB, WENO-5, SSP-33 | 0.55, 0.70, 0.85 | 392.4 B (49.1) | 54 | 1.95 GiB (2.81) | 20.09 GiB (20.96) |
| 3a channel, hybrid, viscous, SSP-54 | 0.30, 0.55, 0.85 | 493.4 B (61.7) | 69 | 0.75 GiB (1.61) | 18.60 GiB (19.46) |

- **The budget overestimates the per-slot cost** by 26 % (step), 10 % (sphere) and 12 % (channel): it counts the three
  face fluxes as full stored blocks, while they are allocated over the interior faces only (in 2-D a tenth of a stored
  block). The plan's ~55 doubles was within 12 % for the 3-D cases.
- **The memory not proportional to the capacity** is 0.75-3.4 GiB per GPU (largest for the three-realm step). GPU 0
  holds 0.87 GiB more than the other three in every run (the "device memory_avail" each rank prints at start already
  shows 0.93 GB less on rank 0). Every rank touches GPU 0 before its own device: the `NVCOMPILER_ACC_TIME` report of
  the ranks on GPUs 1-3 lists 33 uploads to device 0 before their own device's, consistent with three extra contexts
  on GPU 0 (3 x ~0.29 GiB; the context size is not measured). Confirmed 2026-11-24: with each rank seeing only its
  own GPU ("Rank placement" above) all four ranks report the same free memory.
- **At the default fraction 0.85 the A30 (23.5 GiB) keeps 2.5-4.9 GiB free**: no case is near an out-of-memory.

**The matrix table, corrected** (device memory the case needs at its matrix layout, slope x the stored slots of its
blocks in use + intercept, per GPU; GPU 0 + 0.87 GiB):

| Case | Plan estimate | Needed per A30 | Total |
|---|---|---|---|
| 1a step (64 blocks of each realm per GPU, 4 A30) | ~34 GB (WENO-9) | 4.3 GiB of blocks + 3.4 = 7.7 GiB | 31.6 GiB |
| 2a sphere (1024 blocks per GPU, 4 A30) | ~49 GB | 10.1 + 1.95 = 12.1 GiB | 49.1 GiB |
| 3a channel (512 blocks, 1 A30; 128 per GPU on 4) | ~6 GB | 5.9 + 0.75 = 6.7 GiB on one A30 | 6.7 GiB |
| `w2-sphere-l5` (1024 blocks per GPU, 32 A30) | - | 12.1 GiB | 387 GiB |

The step's ceiling ("The step's memory ceiling" in phase 9) is unchanged in kind: the capacity of its binding realm
is set by the budget, and the budget's 64 fields, the null-axis ghosts and the equal share per realm all shrink it;
draft issue `issue-step-memory-ceiling.md`.

### Residual biases (phase 4.4)

Measured in job n3m-A (one node, three interleaved windows per variant, from the matrix's window dumps). The progress
lines (`time_object%print_progress`: 7 lines per realm, per rank, per step) have no input key: a measurement-only patch
of `adam_time_object.F90` (`ADAM_PROGRESS=all|rank0|none`, built in a scratch tree on medooza, not committed) switched
them; `compute_conservation` every step is `[diagnostics] conservation_history_save = 1`.

| Variant | 2a, 4 A30 | 1a, 64 cores |
|---|---|---|
| matrix executable, every output off (the protocol) | 0.7815 s (spread 0.3 %) | 1.253 s (4.5 %) |
| patched executable, progress on every rank | 0.7851 s (0.0 %) | - |
| patched executable, no progress printing | 0.7870 s (0.1 %) | 1.273 s (1.3 %) |
| matrix executable, `compute_conservation` every step | 0.7888 s (0.4 %) | 1.252 s (4.7 %) |

- **Progress printing: no measurable cost.** Removing it does not make a step faster on either backend (FNL +0.2 %
  against the patched build that prints; CPU +1.6 %, inside the 4.5 % spread). Not changed; recorded.
- **`compute_conservation` every step: +0.9 % on 4 A30** (0.0073 s per step, spreads below 0.4 %); not resolved on the
  CPU node (-0.1 %, spreads 4.5-4.7 %).
- Neither is in the matrix's numbers: the windows set `conservation_history_save = 0`, and `time%is_to_save` gates the
  computation itself, not only the write. The warm-ups and settles ran the generated cadence (every step on legs 1 and
  2), outside any window.

### State independence (9.2.5) and a CPU weak point (9.2.4)

- **9.2.5, the state-independence premise: holds to 0.5 %.** The L4 sphere from a 20-step cold start (2a's input, the
  shock still far from the body) against 2a from its warm-up at t_w = 0.15, 4 A30, the same node, three interleaved
  windows (job n3m-A): **0.7776 s against 0.7815 s per step** (spreads 0.4 and 0.3 %). The weak cases, which start
  cold, and 1a/2a, which start developed, are comparable. (2a's window here, 0.7815 s, reproduces the matrix's 0.7799 s
  of job 1369401 to 0.2 %.)
- **9.2.4, the CPU weak point (2026-11-24): 0.53 on 8 nodes, and not because of the immersed boundary.** The uniform
  sphere at 110,592 cells per core, 24^3 blocks, 8 per rank (one sibling set): `wg-b24` (L3, 7.08 M cells) on one node,
  64 ranks, against 2a (L4, 56.6 M) on 8 nodes, 512 ranks, the same day (work tree `speedup-133/cpu-weak`, jobs 1369694,
  1369695; tables in `results/2026-11-24/cpu-weak/`). To test whether the immersed boundary's load imbalance is the
  cause (at 512 ranks the 8 ranks at the sphere's centre hold 59 % solid cells, at 64 ranks 7.4 %), the same pair ran
  without the body (`speedup-133/noib-cpu.sh` on medooza, a `run.sh --exe` wrapper that drops `[solids]` and `[solid_1]`
  from each rank's copy of the input; `cpu-weak-noib`, jobs 1369697, 1369698; `results/2026-11-24/cpu-weak-noib/`):

  | Pair | 1 node, 64 ranks: s per step (spread), g | 8 nodes, 512 ranks: s per step (spread), g | Weak efficiency | Wait at the reduction, 1 / 8 nodes |
  |---|---|---|---|---|
  | with the sphere | 1.234 (3.3 %), 11.2 us | 2.311 (1.2 %), 20.9 us | 0.53 | 0.060 / 0.366 s per step |
  | without the body | 1.035 (9.2 %), 9.36 us | 2.033 (2.8 %), 18.4 us | 0.51 | 0.040 / 0.363 s per step |

  The body costs ~17-19 % on both sides and leaves the efficiency where it is, so it is not the immersed boundary.
  2a's step on 64 / 128 / 256 / 512 ranks costs 12.7 / 6.48 / 3.24 / 2.31 s (g 14.4 / 14.7 / 14.7 / 20.9 us). **The
  cause, found by the phase breakdown ("Where a step's time goes" above), is subnormal arithmetic, not the exchange
  between nodes and not the barriers** (both suspected here first): the developed 2a state holds subnormal values on
  76 of the 512 ranks, which compute their face fluxes twice as slowly while every other rank waits at the barriers and
  reductions (the 0.36 s per step of "wait at the reduction"); with flush to zero the 8-node run drops from 1988 to
  1044 ms per step and the weak efficiency is 0.93-1.0. The cold-started `wg-b24` holds no subnormal value. The two
  1-node measurements of the same grid on different days differ by 7.5 % (1.270 and 1.366 s; 1.234 today): the CPU
  node's spread.

### Hypotheses (phase 8.2)

- **H1** (leg 1 on GPUs: the AMR speedup falls well short of N(a)/N(c), the 2-D runtime grid leaving each GPU near
  the launch floor; on the CPU it comes close). **Reformulated.** The first clause holds: 1.04 against 1.59 on 4 A30.
  The mechanism is refuted: the step is not at a launch floor (1a keeps scaling to 0.26 M cells per A30, strong
  efficiency 0.77 on 16 A30; 1b, the AMR grid frozen, 0.67 at 0.16 M); the loss is the AMR grid (efficiency 0.89) and
  the regrid on the host (0.048 s of the 0.207 s step, 23 %). The CPU clause is refuted: 1.06 against 1.59, the same as
  on the GPU node, for another reason: the runtime layout puts one sibling set of each realm on the busiest rank (1c's
  steps cost 1.17 s against 1b's 0.88, 1.59x the mean load), while the regrid is 1 % of the step.
  Corrected statement: *on the 2-D step, runtime AMR gains nothing at equal finest spacing on either node, the GPU
  node losing it to the host regrid, the CPU node to the layout of the regrid's output.*
- **H2** (leg 1: a 2-D cell costs more than two thirds of a 3-D cell at the same order, by the null-axis ghosts).
  **Confirmed, and located.** Measured (job n3m-A, the same node, three interleaved windows each): the step forest
  extruded in 3-D at the same order (`make_step.py --layout aligned --dims 3 --m 8 --nk 8 --weno 7`, octrees, L2,
  2.06 M cells, cold start of 20 steps) on 2 A30 costs 0.1304 s per step (g 126 ns per cell-step, spread 0.3 %); the
  2-D `w1-step-m8` on 1 A30 0.2083 s (g 202 ns, 0.2 %): both at 1.03 M cells per A30. **A 2-D cell costs 1.60 3-D
  cells**, not two thirds. Kernel times (`NVCOMPILER_ACC_TIME`, 10-step runs, per GPU): the face fluxes follow the
  dimension count (0.54 s in 2-D, 0.88 s in 3-D: 0.62); the boundary-condition kernel `set_boundary_conditions_dev`
  takes 1.34 s in 2-D against 0.16 s in 3-D, 59 % of the 2-D kernel time against 10 %, and its launches cover up to
  19,567 gangs of 128 rows against 901 (21x the rows at the same cells). The rest of the 2-D kernel time is 0.69 of the
  3-D step's. On a null axis each block's two z faces are domain faces, so the `ngc` z ghost layers of every block (8/9
  of a stored 2-D block at `ngc` 4) are boundary rows, refilled at every stage: the null-axis ghosts cost the 2-D step
  time as well as memory (issue draft `issue-step-memory-ceiling.md`, which carries both).
- **H3** (leg 2, FNL: T(R)/T(F) - 1 is dominated by the regrid host round trip; on the CPU it stays a few percent).
  **Confirmed, with a correction of what the round trip is.** On 4 A30, T(2b)/T(2c) - 1 = 38 %, and all of it is
  outside the steps (2b's steps 42.03 s, 2c's 41.94 s per 110 steps): the regrid, which runs on the host (the state
  copied to the host, the markers and the tree there, the state copied back). It is paid at every regrid step, not
  only when the grid changes: ~1.4 s per regrid step that changes nothing, 2.0 s for the one that changed the grid;
  it scales with the blocks per GPU (5.5-9 ms per block), and its split between the copies and the host work is not
  measured. On the CPU the regrid timer is 1.3 % of 2b's step (the steps are 11x longer); the CPU's T(R)/T(F) - 1
  (24 %) is inside the steps (5.70 against 4.63 s), half of it accounted for by the runtime layout (busiest rank
  1.13x the mean against 2c's 1.01x), the rest not located; it is not the regrid.
- **H4** (leg 3: T(3b) ≈ T(3c) + the regrid cost; the AMR speedup comes from the core alone). **Confirmed on the GPU
  node, refuted on the CPU node.** 4 A30: 3b and 3c both spend 0.1428 s per step in the steps, and 3b's 0.0414 s
  outside them is its regrid: 0.1432 + 0.0414 = 0.1846 against 3b's 0.1842. The runtime grid is the a priori
  prediction (256 L3 wall blocks, 32 L2 core blocks: only the core coarsened) and did not change in 4000 settle steps
  after its first regrid, so the AMR speedup (1.25; frozen 1.60; ideal 1.78) comes from the core alone; there was no
  churn to measure, only the regrid steps' fixed cost (0.41 s each, every 10 steps, 29 %). 64 cores: 3b's steps cost
  3.07 s against 3c's 2.21: the runtime layout (busiest rank 98,304 cells, the uniform 3a's per-rank load; 3a costs
  3.17 s) erases the gain (AMR speedup 1.03).

**Caveats**: the commit is not yet recorded (the night tree); the GPU memory preflight checks only the batch node
(draft issue `issue-gpu-preflight-batch-node.md`); one 16-A30 run of 1c was killed by Slurm (SIGKILL, no memory
record) and rerun alone (job 1369494); the leg-3 runs before the re-blocking (layout 64 x 48 x 32) are superseded and
not in the tables; the windows are 110 steps, in which the runtime grids changed once in ten regrids (legs 1, 2) or
never (leg 3): a flow with more churn pays more regrid.

## Files

| Path | What |
|---|---|
| `matrix.conf` | every parameter: per leg the grids, levels, block sizes, scheme, marker and tolerances, t_w, settle length, regrid period, scaling series; the window (10/110 steps, 3 repetitions, 2 % spread); `[leg.N.reduced]` the WSL plumbing sizes; `[weak.<case>]` the weak-scaling cases (phase 9) |
| `make_cases.py` | writes one directory per case (inputs, `case.conf` for run.sh, `README.md`) from `matrix.conf`; calls `verification/step/make_step.py`, `make_sphere.py`, `verification/turbulence/make_channel.py` |
| `make_sphere.py` | leg 2's input: shock over a sphere, immersed boundary, octree, uniform or runtime AMR |
| `run.sh` | the chain: `cases`, `warmup`, `settle`, `windows`, `chain` stages; a `run.info` per run |
| `repartition.py` | a restart written on N ranks rewritten for M ranks (phase 1) |
| `reblock.py` | a restart rewritten one tree level deeper, every block into its children (leg 3's T3 dump) |
| `ini_set.py` | rewrites `section.key=value` options of an input (the per-run changes of run.sh) |
| `analyze.py` | step costs, spreads, N, throughputs, the ratios, the grind times and the strong and weak efficiencies -> `results/summary.{md,csv}` |
| `scaling_plot.py` | the grind-time curves g(cells per device) from `summary.csv` (matplotlib, not needed by `analyze.py`) |
| `grid-study-cases.conf` | the grid-size study's cases (`[weak.wg-*]`), an overlay: `cat matrix.conf grid-study-cases.conf` is the study's `--conf` |
| `grid_study_plot.py` | the grid-size figure from the study's `summary.csv` (and the rank-placement A/B summaries): g against the block edge on 1 A30 and on 1, 2, 4 A30 |
| `phase_compare.py` | ranks the phases of the measurement-only phase-timer reports, one node against several (`results/2026-11-24/phase/phase-compare.md`) |
| `results/2026-11-24/phase/` | the phase breakdown: the timer patch (`phase-timer.patch`, `apply_phase_timer.py`, `adam_phase_timer.F90.txt`; not for commit to the source), every run's report, the network measurements, `subnormal_count.py` |
| `jobs/medooza-bind-fnl.sh` | the FNL executable with each rank pinned to its A30's NUMA node (`run.sh --exe` with an absolute path; "Rank placement") |
| `matrix_plot.py` | the summary figures from `summary.csv`: backend speedup per case, AMR decomposition per leg, grind time (matplotlib; `--theme both` writes the docs' dark twins) |
| `jobs/medooza-{cpu,fnl}.slurm` | the medooza jobs: one run.sh stage per job, outputs on `/ifs` |
| `build.sh` | builds one backend from a clean tree into `exe/<backend>/` (it still calls backends "legs") |
| `templates/sphere-legacy.ini` | the removed matrix's sphere input, a record (stale, does not run) |
| `profile-couette3d/` | the profiling case (outside the matrix; shares the executables of `build.sh`) |

The run tree (`--work`, default `runs/full` or `runs/reduced` here):

```
<work>/cases/<case>/                         inputs, case.conf, README.md (rewritten by every run.sh call)
<work>/dumps/leg<N>/{warmup,settle}/         the two dumps of a leg, with the logs of the runs that wrote them
<work>/dumps/<weak case>/warmup/             the dump of a weak case (phase 9: a cold start, no settle)
<work>/windows/<backend>-np<N>/<case>/dump/  the window dump repartitioned for N ranks
<work>/windows/<backend>-np<N>/<case>/r<rep>-n<steps>/   one window run: inputs, log.txt, wall.txt, run.info
<work>/attic/                                directories a stage found already there (moved, never deleted)
```

## How to run

**WSL: dry runs and the plumbing test only, never timings.** The reduced sizes (`--reduced`) run every case through
the whole chain in seconds per run (WSL: one run at a time, under the shared lock when agents run):

```bash
cd src/tests/flume/speedup
./build.sh cpu --varset local_gnu                                 # or --exe an existing executable
./run.sh chain --reduced --cases 1 --dry-run                      # the commands, nothing run
./run.sh warmup  --reduced --cases 1,2,3a --np 2                  # dumps written on 2 ranks...
./run.sh settle  --reduced --cases 1,2 --np 2
./run.sh windows --reduced --cases 1,2,3a --np 1                  # ...read on 1 (repartitioned)
./analyze.py --work runs/reduced --out runs/reduced/results       # a summary of meaningless timings
```

Tested this way on 2026-11-22 (a CPU executable of develop `0bcdf870`, without the `mode = 7` marker): legs 1 and 2 through warm-up and
settle on 2 ranks, their windows and 3a's on 1 rank (repartitioned), 2a also on 2 ranks; `analyze.py` built the summary.
3b and 3c need the `mode = 7` marker (phase 2) in the executable: their inputs are written, not run.

**medooza** (batch only: nothing but `sbatch` on the login node; at most 4 jobs at once; runs on `/ifs`, never in
the repository; prod-gn FNL needs UCC off, which the FNL job sets). Build in a batch job, then submit from
`src/tests/flume/speedup` of the clone; the jobs pass their arguments to `run.sh` and set `--backend`, `--np` (the
job's tasks) and `--work` (`$SPEEDUP_WORK`, default `/ifs/hpc/home/szaghi/adam-tests/speedup-133`):

```bash
# builds (inside batch jobs): ./build.sh cpu --varset iac_gnu ; ./build.sh fnl-nvf --varset iac_nvf
sbatch jobs/medooza-fnl.slurm warmup --cases 2                    # leg 2 warm-up, 4 A30 (the fastest backend)
sbatch --dependency=afterok:<job> jobs/medooza-fnl.slurm settle --cases 2
sbatch --dependency=afterok:<job> jobs/medooza-fnl.slurm windows --cases 2         # fnl-nvf, 4 A30
sbatch --dependency=afterok:<job> jobs/medooza-cpu.slurm windows --cases 2         # cpu, 64 ranks
sbatch --gres=gpu:1 --ntasks-per-node=1 jobs/medooza-fnl.slurm windows --cases 3a  # one A30, where the case fits
```

The windows and chain stages end with `analyze.py --work $SPEEDUP_WORK --out $SPEEDUP_WORK/results`, which reads
every backend's windows found there. The FNL job refuses a node whose GPUs have less than 20000 MiB free (exit 4:
resubmit with `-x <node>`). Leg 3 from the T3 dump: re-block it first ("Re-blocking the T3 dump" below), then `sbatch jobs/medooza-fnl.slurm
warmup --cases 3 --use-dump <out dir>` (the directory `reblock.py` wrote, holding `channel-restart.*`).

**The grid-size study and the rank placement** (2026-11-24; the study's own work tree, so the matrix's is untouched):

```bash
W=/ifs/hpc/home/szaghi/adam-tests/speedup-133; cat matrix.conf grid-study-cases.conf > $W/grid-study.conf
R() { ./run.sh "$@" --backend fnl-nvf --conf $W/grid-study.conf --work $W/grid-study --launcher "srun --mpi=pmix -n"; }
R warmup --np 1 --cases wg-b24          # in a 1-node FNL job, one case at a time (wg-b40 is refused: over the budget)
R windows --np 4 --cases wg-b24,wg-b38  # likewise --np 1, 2; then analyze.py --work $W/grid-study
R chain --np 32 --cases wg-l4b38        # in an 8-node job: cold start and windows on 32 A30
# rank placement, in an --exclusive job: the same windows bound, in another work tree
./run.sh windows --backend fnl-nvf --np 4 --conf $W/grid-study.conf --work $W/bind/on-np4 --cases wg-b38 \
   --exe $PWD/jobs/medooza-bind-fnl.sh --launcher "srun --mpi=pmix -c 16 --cpu-bind=none -n"
```

A work tree other than the study's needs the warm-up dumps linked: `ln -s $W/grid-study/dumps/<case>
<work>/dumps/<case>`. The CPU node: `./run.sh windows --backend cpu --np 64` with the same conf and cases.

**AMD node, `fnl-amd`** (Giacomo Rossi). The windows start from the dumps written on medooza, so every backend times
the same bits (D3):

1. Use the commit of the medooza runs (`run.info` `commit=`), so that `matrix.conf` and the generators are the same.
2. Build: `./build.sh fnl-amd --varset "thera_amd gfx942" --jobs 16` (the toolchain notes in
   `profile-couette3d/README.md`, "AMD node").
3. Copy the dumps of the legs to measure from medooza into the same layout,
   `<work>/dumps/leg<N>/{warmup,settle}/` (only the restart files `*-restart*` are needed: `.time`, `.tnd`,
   `.reference`, `-proc*.fbd`; a leg-2 L4 dump is ~5 GB).
4. Run the windows, the whole node and one GPU (one MPI rank per GPU):

   ```bash
   ./run.sh windows --backend fnl-amd --np 8 --work <work> --cases 2 --launcher "mpirun -x ADAM_SEAM_TIMING -np"
   ./run.sh windows --backend fnl-amd --np 1 --work <work> --cases 3a
   ```

   `run.sh` repartitions each dump for the rank count (`repartition.py`, stdlib Python 3), regenerates the inputs and
   writes `<work>/windows/fnl-amd-np<N>/`.
5. Send back `<work>/windows/fnl-amd-np*/` without the dumps (`rsync -a --exclude dump --exclude '*.fbd'
   --exclude '*-restart*' ...`): `run.info`, `log.txt` and `wall.txt` per run are what `analyze.py` reads.

### Re-blocking the T3 dump (leg 3)

`reblock.py` rewrites a restart one tree level deeper: every block becomes its `ratio` children (octree 2x2x2,
quadtree 2x2 in x and y, binary tree 2 in x), each child the octant of its parent's cells, ghosts included. No
interpolation: the same cells, the same h, the tree one level deeper; the children are laid out over the ranks by
`repartition.py`'s rules (`--to`, default the dump's rank count; `--rule units`, the default). It rewrites the `.fbd`
files and the `.tnd`, copies `.time`, `.reference`, `.sinks` and `.dff`, refuses `.wmf` (wall-model filter averages
live on block faces) and does not write the `.h5`/`.xdmf` mirror (no restart reads it). The input that reads the result
has `[grid] ni, nj, nk` halved on the split axes and `[amr] max_level`, `iu_ref_levels` one higher.

Leg 3 needs it: the T3 `resolved-hybrid` dump (L2, 64x48x32 blocks, 4 ranks) becomes L3 blocks of 32x24x16, the
layout of `[leg.3]` (`max_level` 3, `base_level` 2). On medooza, in a batch job (a few seconds, ~0.5 GB written):

```bash
sbatch -A iacroma -p dev-gn --mem 4G -t 00:10:00 \
       --wrap "python3 -I reblock.py <T3 dir> --out <out dir>"     # <T3 dir> holds channel-restart.*
sbatch jobs/medooza-fnl.slurm warmup --cases 3 --use-dump <out dir>
```

Tested on WSL (CPU, 2026-11-23): a 64x48x32 channel (`make_channel.py`, the leg-3 physics), 10 steps on 2 ranks, its
L2 dump (64 blocks of 16x12x8) continued 10 steps as it is and re-blocked (512 blocks of 8x6x4 at L3), on 2 ranks.
The 98,304 cells agree to 1.9e-14 (absolute, every variable). The difference comes from the mass-flow forcing, a
global sum whose order follows the blocks: without `[forcing]` the two are bitwise. The re-blocked dump repartitioned
to 3 ranks (`repartition.py`) and continued agrees to 1.5e-12, as the L2 dump repartitioned to 3 ranks does (1.2e-12,
the forcing sum again); without `[forcing]` it is bitwise too.

## Multi-node and weak scaling (phase 9)

Owner, 2026-11-23: the plan had "strong and weak scaling beyond one node" out of scope; it is now in. medooza `prod-gn`
nodes hold 4 A30 (24 GB) and 64 cores; up to 8 nodes (32 A30, 512 cores). The FNL job runs 4 ranks per node, one per
GPU (`--ntasks-per-node=4`), so `--nodes=N` gives `--np 4N` (`SLURM_NTASKS`); the CPU job 64 ranks per node.

### 9.1 Strong scaling: the legs on more nodes (no new case)

The windows of the existing cases, their dumps repartitioned for the rank count (`repartition.py`), exactly as on one
node. FNL np 8/16/32 = 2/4/8 nodes; CPU 2/4 nodes = 128/256 ranks. Leg 3 stays on one node (3a is ~6 GB: no multi-node
load).

```bash
sbatch --nodes=2 jobs/medooza-fnl.slurm windows --cases 2        # np 8
sbatch --nodes=4 jobs/medooza-fnl.slurm windows --cases 2        # np 16
sbatch --nodes=8 jobs/medooza-fnl.slurm windows --cases 2        # np 32
sbatch --nodes=2 jobs/medooza-fnl.slurm windows --cases 1        # likewise 4, 8 for leg 1
sbatch --nodes=2 jobs/medooza-cpu.slurm windows --cases 2        # 128 ranks
sbatch --nodes=4 jobs/medooza-cpu.slurm windows --cases 2        # 256 ranks; likewise --cases 1
```

`--cases 1` and `--cases 2` select the leg's matrix cases only (the weak cases are named `w...`). Two things to know:
the FNL job's free-memory preflight runs `nvidia-smi` on the first node only; and at 32 ranks a runtime grid (1c, 2b)
of a few hundred blocks leaves ranks with few or no blocks of a realm (an empty rank was a crash once, #82, MHD).

### 9.2 Weak scaling: uniform cases at fixed cells per GPU

New uniform cases (`[weak.<case>]` in `matrix.conf`, role U only: no runtime or frozen twin, no settle). Each warms up
from a **cold start of 20 steps** (`warmup_steps`, `[weak]`; the per-cell cost of a uniform WENO grid does not depend
on the flow state) into its own dump `<work>/dumps/<case>/warmup/`; the window is the protocol's, 3 x (10, 110) steps.
1a and 2a are reused as they are (their windows from the legs' warm-ups at t_w), joined to the series by the leg's
`series` key. The names start with `w`; select them with `--cases w`, `w1`, `w2` or a name.

**The sphere is the weak-scaling and grind-time series** (WENO-5, 24^3 blocks: the block shape is the same at every
level, only the block count changes):

| Load (cells / A30) | Case | Cells | A30 | Nodes | Role in the analysis |
|---|---|---|---|---|---|
| 14.2 M | `2a-sphere-uniform` (L4) | 56.6 M | 4 | 1 | weak reference |
| 14.2 M | `w2-sphere-l5` (L5, 768^3) | 453 M | 32 | 8 | weak point |
| 7.08 M | `w2-sphere-l3` (L3, 192^3) | 7.08 M | 1 | 1 | weak reference |
| 7.08 M | `2a-sphere-uniform` | 56.6 M | 8 | 2 | weak point (9.1's np 8) |
| 3.54 M | `w2-sphere-l3` | 7.08 M | 2 | 1 | grind curve (optional: pairs 2a at 16) |
| 3.54 M | `2a-sphere-uniform` | 56.6 M | 16 | 4 | grind curve (9.1's np 16) |
| 1.77 M | `w2-sphere-l3` | 7.08 M | 4 | 1 | grind curve, weak reference |
| 1.77 M | `2a-sphere-uniform` | 56.6 M | 32 | 8 | grind curve, weak point (9.1's np 32) |

**The step series is memory-capped 2-D (~1 M cells per A30; the ceiling below).** It was planned as launch-bound
evidence (H1), on the reading that an A30 is launch-bound below a few M cells (`profile-couette3d`: `small`, 2.1 M, is
the launch-overhead profile; `medium`, 16.8 M, loads an A30). The measurements refute that for the step: 1a keeps
scaling on 8 and 16 A30 (strong efficiency 0.91 and 0.77 at 0.52 and 0.26 M cells per A30, "Results"); the grind time
of the step (209 ns per cell-step) is higher than the sphere's (55) because its cell-step does more work (WENO-7
characteristic, SSP-54: 5 stages, against WENO-5 and SSP-33), not because the GPU waits. Its block shape also changes
with m (realm C stores 12.2 cells per cell at m 8, 9.8 at m 32: the ghost share changes along the series),
so its "weak efficiency" mixes block size with communication:

| Case | m | Cells (L4) | A30 | Nodes | Cells / A30 |
|---|---|---|---|---|---|
| `w1-step-m8` | 8 | 1.03 M | 1 | 1 | 1.03 M |
| `1a-step-uniform` | 16 | 4.13 M | 4 | 1 | 1.03 M |
| `w1-step-m22` | 22 | 7.81 M | 8 | 2 | 0.976 M (m^2 = 512 has no even root; 5.4 % below) |
| `w1-step-m32` | 32 | 16.5 M | 16 | 4 | 1.03 M |

Every m is even and >= 2 ngc = 8 (WENO-7), as `make_step.py --layout aligned` requires.

```bash
# sphere (the weak and grind series); 2a at np 4 is phase 5's, at 8/16/32 is 9.1's
sbatch --gres=gpu:1 --ntasks-per-node=1 jobs/medooza-fnl.slurm chain --cases w2-sphere-l3          # L3 on 1 A30
sbatch --dependency=afterok:<L3 job> --gres=gpu:2 --ntasks-per-node=2 jobs/medooza-fnl.slurm windows --cases w2-sphere-l3
sbatch --dependency=afterok:<L3 job> jobs/medooza-fnl.slurm windows --cases w2-sphere-l3          # L3 on 4 A30
sbatch --nodes=8 jobs/medooza-fnl.slurm chain --cases w2-sphere-l5                                 # L5 on 32 A30
# step (memory-capped); 1a at np 4 is phase 6's
sbatch --gres=gpu:1 --ntasks-per-node=1 jobs/medooza-fnl.slurm chain --cases w1-step-m8
sbatch --nodes=2 jobs/medooza-fnl.slurm chain --cases w1-step-m22
sbatch --nodes=4 jobs/medooza-fnl.slurm chain --cases w1-step-m32
```

The `chain` of a weak case is its warm-up then its windows (no settle). CPU, optional, the same per-node load:
`sbatch jobs/medooza-cpu.slurm chain --cases w2-sphere-l3` (1 node) against 2a on 8 CPU nodes
(`sbatch --nodes=8 jobs/medooza-cpu.slurm windows --cases 2a`).

**Memory** (FNL, the block budget of `adam_adam_object` `compute_blocks_number`: `0.85 x device memory / realms /
(fields x 8 B x stored cells per block)`, stored cells `(ni+2 ngc)(nj+2 ngc)(nk+2 ngc)` with the null axis too; fields
64 for the step at WENO-7/SSP-54 and 54 for the sphere at WENO-5/SSP-33, the counts the WSL "blocks budget" lines
print; A30 taken as 24 GB, the conservative end: the job log's "blocks budget" line gives the exact capacity). The share
is the binding realm's blocks per rank over its capacity:

| Case | A30 | Blocks per rank (binding realm) | Capacity | Used | Fits |
|---|---|---|---|---|---|
| `w1-step-m8` | 1 | 256 (C, 104 x 40 x 9 stored) | 355 | 72 % | yes |
| `1a-step-uniform` | 4 | 64 (C, 200 x 72 x 9) | 102 | 63 % | yes |
| `w1-step-m22` | 8 | 32 (C, 272 x 96 x 9) | 57 | 56 % | yes |
| `w1-step-m32` | 16 | 16 (C, 392 x 136 x 9) | 28 | 57 % | yes |
| `w2-sphere-l3` | 1 / 2 / 4 | 512 / 256 / 128 (30^3 stored) | 1749 | 29 / 15 / 7 % | yes |
| `2a-sphere-uniform` | 4 / 8 / 16 / 32 | 1024 / 512 / 256 / 128 | 1749 | 59 / 29 / 15 / 7 % | yes |
| `w2-sphere-l5` | 32 | 1024 | 1749 | 59 % (11.9 GB of fields per A30) | yes |

Disk (`/ifs`): a restart stores q with its ghosts, `blocks x 30^3 x 5 x 8 B` for the sphere: L3 0.55 GB, L4 4.4 GB,
**L5 35 GB**, written by the warm-up and again by the window's repartition (~71 GB for `w2-sphere-l5`, plus one 4.4 GB
copy per rank count of 2a); the step's dumps are 0.5-6.6 GB.

**The step's memory ceiling** (a follow-up for the owner, not changed here). Two causes, with their factors at m 16:

1. **Null-axis ghosts.** A 2-D quadtree block stores `ngc` ghost layers on each side of its single z cell too: realm C's
   192 x 64 x 1 block is stored as 200 x 72 x 9 (10.5 stored cells per cell; 1.17 without the z ghosts). Without them
   the step's memory per cell drops ~9x.
2. **FLUME's equal memory split per realm** (`adam_flume_fnl_object.F90` divides the device by the realm count): realm
   C holds 76 % of the cells with 1/3 of the memory and binds first.

The ceiling is then ~1.6 M cells per A30 today (realm C full; at 1.03 M it is 63 % full), ~3.6 M with a split
proportional to the cells, ~30 M with both fixed. (The ~34 GB of the plan's memory note for 1a was the WENO-9, ngc 5
estimate; at WENO-7 1a's blocks take ~23 GB over the 4 A30.)

### 9.3 Grind time and the efficiencies (`analyze.py`)

Per window (backend, case, ranks; devices = ranks: GPUs of an FNL run, cores of a CPU run):

- **grind time** `g = step cost x devices / N`, ns per cell-step (the step cost of the protocol: the median of
  `(T_loop(110) - T_loop(10)) / 100`); **cells per device** `N / devices`. g against cells per device is the grind
  curve; a multi-node point above the single-node curve at the same load is the multi-node overhead.
- **strong efficiency**, per backend and case: `g(ref) / g(np)`, ref the smallest rank count present (4 for the
  legs' FNL windows, 64 for CPU); 1 is ideal.
- **weak efficiency**, per backend and series (`series` of the leg in `matrix.conf`: the uniform case and the weak
  cases of the leg): the windows are grouped by cells per device (a group: loads within `load_tolerance` = 10 % of its
  largest, `[protocol]`), and `g(ref) / g(np)` with ref the group's smallest device count; a group of one window has
  none.

`summary.md` gains the table "Scaling and grind time" (backend, series, case, devices, N, cells per device, step
cost, g, strong efficiency and its reference, load group, weak efficiency and its reference), sorted by backend,
series and load, and the series labels; `summary.csv` the same columns (`kind`, `series`, `cells_per_device`,
`strong_ref_np`, `strong_efficiency`, `weak_group`, `weak_ref_np`, `weak_efficiency`). The AMR ratios take the
matrix cases only. `scaling_plot.py results/summary.csv --out scaling.png` draws g(cells per device) per backend
(filled: one node, open: several; matplotlib, outside the stdlib-only scripts).

The premise that the per-cell cost does not depend on the flow state (weak cases from a 20-step cold start, 1a and 2a
from t_w) holds to 0.5 % ("State independence" above).

Tested on WSL at reduced sizes (2026-11-23, CPU, the `[*.reduced]` sizes: step m 8/12/16 at L0, sphere blocks 8 and 10
at L1): warm-ups on 2 ranks, windows on 1, 2 and 4 ranks (repartitioned), `analyze.py` grouped w1-step-m8 on 1 rank
with 1a and w1-step-m32 on 4 (4,032 cells per rank), and w2-sphere-l3 on 1 with w2-sphere-l5 on 2 (4,096 / 4,000).

## Decisions

Owner, 2026-11-20:

- **D1, equal finest h.** Every AMR case refines to the uniform case's level with the same block size: same solution
  where the AMR grid is finest, same global dt. The AMR speedup is measured, not estimated, and the per-step ratio is
  the time-to-solution ratio.
- **D2, a fixed grid is the runtime grid frozen.** The frozen case restarts from the dump that starts the runtime
  case's window with `[amr] frequency = 0`: same grid, same state.
- **D3, one warm-up per leg.** Restart dumps are per rank (`<basename>-proc<rank>.fbd`); `repartition.py` makes a
  dump written on N ranks readable on M. Each leg warms up once, on the fastest backend, and every backend and rank
  count starts its windows from the same bits.
- **D4, the channel marker is the local Kolmogorov resolution** `r = h/eta` (phase 2, `[amr_marker] mode = 7`).
- **D5, words.** *Leg* = a flow of the matrix (as `--legs` in the verification scripts); *case* = one grid of a leg
  (1a, 1b, ...); *backend* = one build made by `build.sh`. `profile-couette3d` and `build.sh` keep calling builds
  "legs" until they are next touched.

Owner, 2026-11-22 (the re-analysis against `0bcdf870`, comment on #133):

1. **Leg 1 stays 2-D** (D1 holds on every backend); a **3-D demonstration step at h = 1/320 on 4 A30** is added for
   the figure of #142.
2. The paper commit is tagged after the step and #142 Part A.
3. The immersed-boundary rectangle distance bug goes to its own issue (#143).
