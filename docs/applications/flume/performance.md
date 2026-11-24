# Performance: GPU and CPU nodes, runtime AMR, scaling

What a GPU node buys over a CPU node on three flows, what runtime AMR buys at equal finest spacing and where its
overhead goes, and how the uniform and AMR grids scale across nodes. The measurements are the speedup matrix of
[#133](https://github.com/szaghi/adam/issues/133); the harness, the protocol and every table are in
[`src/tests/flume/speedup/README.md`](https://github.com/szaghi/adam/blob/develop/src/tests/flume/speedup/README.md),
the raw tables in its `results/2026-11-23/summary.{md,csv}`. The kernel-level profile of one GPU (where a step's time
goes inside the A30) is the companion study
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

**Hardware.** medooza `prod-gn`: a node holds 4 NVIDIA A30 (24 GB) and two Intel Xeon Gold 6338 (64 cores). A *GPU
node* is 4 A30 with one MPI rank per GPU (FNL, OpenACC); a *CPU node* is 64 MPI ranks (CPU backend).

## GPU node against CPU node

![T(CPU node) / T(GPU node) per case](/flume/perf-backend-speedup.png)

| Leg | Uniform | Frozen AMR grid | Runtime AMR |
|---|---|---|---|
| 1 step (2-D, WENO-7, SSP-54) | 5.83 | 5.75 | 5.72 |
| 2 sphere (3-D, WENO-5, IB, SSP-33) | 16.3 | 12.2 | 10.9 |
| 3 channel (3-D, hybrid, viscous, SSP-54) | 13.8 | 15.5 | 16.7 |

A GPU node runs a step 11 to 17 times faster than a CPU node on the 3-D flows and 5.8 times faster on the 2-D step. The per-cell work differs between
the legs, so their grind times are not comparable across legs; within the step the GPU's grind time is
209 ns per cell-step at 1.03 M cells per A30, the CPU's 19.5 µs per cell-step per core.

## Runtime AMR at equal finest spacing

![Step cost of the AMR grids over the uniform grid's, split inside and outside the time steps](/flume/perf-amr-decomposition.png)

| Leg | Node | AMR speedup T(U)/T(R) | Ideal N(U)/N(R) | Frozen grid T(U)/T(F) | Grid efficiency | T(R)/T(F) − 1 |
|---|---|---|---|---|---|---|
| 1 | GPU (4 A30) | 1.04 | 1.59 | 1.40 | 0.89 | 34 % |
| 1 | CPU (64) | 1.06 | 1.59 | 1.42 | 0.90 | 34 % |
| 2 | GPU (4 A30) | 1.48 | 4.02 | 2.04 | 0.51 | 38 % |
| 2 | CPU (64) | 2.20 | 4.02 | 2.73 | 0.68 | 24 % |
| 3 | GPU (4 A30) | 1.25 | 1.78 | 1.60 | 0.90 | 29 % |
| 3 | CPU (64) | 1.03 | 1.78 | 1.43 | 0.81 | 39 % |

The grid efficiency is $[T(U)\,N(F)/N(U)]/T(F)$: 1 when an AMR cell costs what a uniform cell costs. Reading the
table:

- **The frozen AMR grid keeps most of the cell saving on the step and the channel** (efficiency 0.81–0.90), and
  half of it on the sphere on the GPU (0.51): a cell of the sphere's AMR grid costs 109 ns per step on 4 A30, a cell of
  the uniform grid 55 ns, at a load (3.5 M cells per A30) where the uniform grid's grind time is flat
  (51–55 ns from 3.5 to 14 M cells per A30). Kernel times locate it: at the same 254–256 blocks per GPU, the face fluxes cost the same on
  the AMR and the uniform grid, while the intra-rank ghost fill (one map-driven kernel that copies same-level ghosts
  and, on the AMR grid, also interpolates coarse-to-fine and restricts fine-to-coarse ones) costs 5.6 times as much:
  96 % of the difference.
- **On the GPU node the runtime overhead is the regrid**: the steps of R cost what the steps of F cost (2b 42.03 s,
  2c 41.94 s per 110 steps), and the regrid, done on the host, costs 0.26–1.4 s per regrid step on 4 A30 when the grid
  does not change and 2.0–2.5 s when it does (in these windows it changed once in ten regrid steps on legs 1 and 2,
  never on leg 3). It scales with the blocks per GPU (sphere: 5.5–9 ms per block), so it shrinks with more GPUs but
  stays 22–27 % of a runtime step that costs 0.2–0.5 s.
- **On the CPU node the regrid is cheap** (1.3 % of the sphere's runtime step) **but the runtime layout is not**: a
  runtime regrid distributes whole sibling sets (4 or 8 blocks) and splits every realm from rank 0, so on 64 ranks the
  busiest rank of the step's runtime case holds 64,512 cells (one sibling set of each realm: the per-rank load of the
  *uniform* 1a) against a mean of 40,680, and the channel's 98,304 (the uniform 3a's) against 55,296. The runtime step
  then costs what the uniform step costs (1c 1.19 s against 1a 1.26 s; 3b 3.08 s against 3a 3.17 s). Run on the runtime layout, the frozen
  grid costs what the runtime grid costs (1b on 1c's layout, 64 ranks: 1.184 s against 1b's own 0.887 s and 1c's
  1.187 s): the layout, not the regrid, is the CPU's runtime overhead.

## Scaling across nodes

![Grind time against cells per device](/flume/perf-grind.png)

**Strong scaling** (step cost in s; efficiency against one node):

| Case | 4 / 8 / 16 / 32 A30 | Efficiency 8 / 16 / 32 | CPU 64 / 128 / 256 ranks | Efficiency 128 / 256 |
|---|---|---|---|---|
| 1a step uniform | 0.216 / 0.119 / 0.070 / – | 0.91 / 0.77 / – | 1.26 / 0.642 / 0.328 | 0.98 / 0.96 |
| 1b step frozen | 0.154 / 0.089 / 0.057 / – | 0.87 / 0.67 / – | 0.887 / 0.544 / 0.334 | 0.81 / 0.66 |
| 1c step runtime | 0.207 / 0.119 / 0.076 / – | 0.87 / 0.69 / – | 1.19 / 1.18 / 1.18 | 0.50 / 0.25 |
| 2a sphere uniform | 0.780 / 0.491 / 0.260 / 0.167 | 0.79 / 0.75 / 0.58 | 12.7 / 6.48 / 3.24 | 0.98 / 0.98 |
| 2b sphere runtime | 0.526 / 0.350 / 0.312 / 0.188 | 0.75 / 0.42 / 0.35 | 5.76 / 3.95 / 3.06 | 0.73 / 0.47 |
| 2c sphere frozen | 0.381 / 0.284 / 0.257 / 0.161 | 0.67 / 0.37 / 0.30 | 4.64 / 3.13 / 1.97 | 0.74 / 0.59 |
| 3a channel uniform | 0.800 (1 A30) / 0.230 (4) | 0.87 (1 → 4) | 3.17 (64) | – |

**Weak scaling and grind time** (g in ns per cell-step; weak efficiency $g(\text{ref})/g$ at the same cells per A30):

| Series | Cells per A30 | A30 (nodes): g | Weak efficiency |
|---|---|---|---|
| sphere | 14.2 M | 4 (1): 55.1; 32 (8): 69.3 | 0.80 |
| sphere | 7.08 M | 1 (1): 51.0; 8 (2): 69.4 | 0.73 |
| sphere | 3.54 M | 2 (1): 52.7; 16 (4): 73.4 | 0.72 |
| sphere | 1.77 M | 4 (1): 64.8; 32 (8): 94.6 | 0.68 |
| step | 1.03 M | 1: 202; 4: 209; 8 (2): 212; 16 (4): 214 | 0.97 / 0.96 / 0.95 |

- **Uniform grids scale**: 0.96–0.98 on 4 CPU nodes; 0.79–0.91 on 2 GPU nodes and 0.75–0.77 on 4, 0.58 for the sphere
  on 8 GPU nodes at 1.77 M cells per A30. The largest run is the L5 sphere, 453 M cells on 32 A30 (8 nodes), at 0.98 s per step.
- **On one GPU node the grind time is flat** at 51–55 ns from 3.5 to 14 M cells per A30 and rises below (65 ns at
  1.77 M). Leaving the node costs the sphere 20–32 % at equal load (weak efficiency 0.68–0.80); the 2-D step 3–5 %
  (why it loses less is not measured).
- **The AMR grids of the sphere stop scaling below about 1 M cells per A30** (2c costs 0.26 s per step on 16 A30, the
  same as 2a with four times the cells), and the CPU's runtime step stops at 64 ranks (the layout above).
- **The per-cell cost does not depend on the flow state**: a cold start of the L4 sphere (20 steps, the shock far
  from the body) costs 0.7776 s per step on 4 A30, the developed state of 2a (t = 0.15) 0.7815 s, on the same node:
  0.5 %. The weak cases, which start cold, and 2a, which starts developed, are comparable.

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
The 2-D step pays for its null axis: a quadtree block stores `ngc` ghost layers along its single $z$ cell (realm C's
192 × 64 × 1 block is stored as 200 × 72 × 9), and each realm of a forest gets an equal share of the device; together
they cap the step at about 1.7 M cells per A30.

## What the measurements say about the hypotheses of the plan

The plan of [#133](https://github.com/szaghi/adam/issues/133) stated four hypotheses before the runs:

| Hypothesis | Verdict | The measurement |
|---|---|---|
| H1: on the step, runtime AMR falls well short of its ideal on GPUs (launch floor), and comes close on CPUs | **reformulated** | short of the ideal on both nodes (1.04 and 1.06 of 1.59); no launch floor (the step scales to 0.26 M cells per A30); the losses are the AMR grid and the host regrid on the GPU node, the runtime layout on the CPU node |
| H2: a 2-D cell costs more than two thirds of a 3-D cell at the same order, by the null-axis ghosts | **confirmed** | at the same order and 1.03 M cells per A30 a 2-D cell costs 1.60 3-D cells (202 against 126 ns); the boundary-condition kernel, which refills the null axis's ghost layers every stage, takes 59 % of the 2-D step's kernel time against 10 % in 3-D |
| H3: on the sphere, the GPU runtime overhead is the regrid's host round trip; on the CPU it stays a few percent | **confirmed** | 38 % on 4 A30, all outside the steps, paid at every regrid step (1.4 s) whether or not the grid changes; 1.3 % on the CPU node, whose runtime overhead is inside the steps (half of it the runtime layout) |
| H4: on the channel, T(runtime) ≈ T(frozen) + regrid, the gain from the core alone | **confirmed on the GPU node, refuted on the CPU node** | 4 A30: 0.1432 + 0.0414 = 0.1846 s against 0.1842; the grid coarsens the core only. 64 cores: the runtime layout leaves the busiest rank with the uniform grid's load (gain 1.03) |

## Limits

- One machine (A30, Xeon Gold 6338); the AMD MI300X node of the plan is not measured yet.
- The windows are short (110 steps): a runtime grid changed once in ten regrids on legs 1 and 2 and never on leg 3, so
  a flow with more churn pays more regrid than shown here.
- Two residual biases were measured and are not in these numbers: printing the progress lines on every rank costs
  nothing measurable, and computing the conservation integrals every step costs 0.9 % on 4 A30; the windows run with
  every diagnostic output off.

## Open issues found

The measurements located five costs, recorded for follow-up issues (drafts of the #133 night, not yet filed):

| Cost | Where it shows | Proposal |
|---|---|---|
| The regrid runs on the host at every regrid step | 22–27 % of a runtime step on 4 A30, also when the grid does not change | per-phase timers, then marking on the device |
| The runtime layout keeps sibling sets whole and splits every realm from rank 0 | the CPU runtime step stops at 64 ranks; legs 1 and 3 gain nothing from runtime AMR on the CPU node | offset or joint split of the realms, weighted by cells; split sibling sets between regrids |
| A cell of the sphere's AMR grid costs twice a uniform cell on the A30 | grid efficiency 0.51 on 4 A30 | one kernel per kind of ghost row (copy, interpolation, restriction), then a cheaper interpolation |
| A 2-D block stores ghost layers along its null axis; realms share the device equally | the step fits ~1.7 M cells per A30 | no ghosts on a null axis; device share in proportion to stored cells |
| The FNL job's free-memory check sees the batch node only | multi-node jobs | check every node |

*Measured: ADAM `0bcdf870` + the night-2026-11-22 tree (commit pending), CPU gfortran 16.2.0 + Open MPI 5.0.11
(`iac_gnu`) on medooza (Intel Xeon Gold 6338, 64 ranks per node), FNL nvfortran 26.9, CUDA 12.9, LTO (`iac_nvf`) on
NVIDIA A30 (4 per node), 2026-11-23 (`speedup/run.sh` windows, `analyze.py`; jobs 1369400-1369597 of the matrix and
1369616 (n3m-A) of the memory, bias, state, 3-D step and kernel-time checks; `matrix_plot.py` draws the figures from `summary.csv`).*
