# FLUME speedup matrix

**Status: harness ready, no leg measured yet** ([issue #133](https://github.com/szaghi/adam/issues/133), the plan's only
source of truth: every step is ticked there with its commit, every phase ends with a comment carrying its results).
Phase 3 (this harness) is in place and tested on WSL at reduced sizes; the legs run in phases 5-7, the summary and the
figures come in phase 8. Phase 9 (owner, 2026-11-23) lifts the plan's "no scaling beyond one node": strong scaling of the
legs on up to 8 medooza nodes, weak scaling at fixed cells per GPU and the grind-time analysis, "Multi-node and weak
scaling (phase 9)" below.

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
| `3a-channel-uniform` | uniform | CKM channel, M_b 1.5, Re_b 3000, isothermal walls | 3-D, octree, periodic x, z | `hybrid-kep-weno`, viscous, SSP-54 | uniform L2, 64x48x32 blocks (the T3 `resolved-hybrid` dump) | 4pi/256 x 2/192 x (4pi/3)/128 | 6.29 M |
| `3b-channel-runtime` | runtime | same | same | same | `r = h/eta` marker (mode 7, refine 4, derefine 1.9: placeholders until 7.2), regrid every 10 steps, base L1 | as 3a | measured |
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
- **Leg 3 has the layout of the T3 dump** (L2, 64x48x32 blocks: `load_blocks` refuses another layout), not the plan's
  first L3, 32x24x16; the runtime case coarsens to L1. Its marker thresholds are fixed a priori on a 3a checkpoint in
  phase 7.2 (`refine_tol`, `derefine_tol` in `[leg.3]`; the parser needs `derefine_tol < refine_tol / 2`).

**Roles.** The plan's letters differ between legs (leg 1: b frozen, c runtime; legs 2 and 3: b runtime, c frozen), so
the harness names the cases by role: **U** uniform, **R** runtime, **F** frozen (`[amr] frequency = 0` on the runtime
grid, D2).

## Protocol

**The chain of a leg** (D3), `run.sh` stages:

1. **warmup**, once per leg, on the fastest backend: the uniform case from the cold start to t_w (leg 1 2.0, leg 2
   0.15, leg 3 150), restart dump in `<work>/dumps/leg<N>/warmup/`. Leg 3: a surviving T3 `resolved-hybrid` dump
   replaces the run (`run.sh warmup --cases 3 --use-dump DIR` links its files). Then run `settle` and `windows`, not
   `chain`: `chain` runs its own warm-up, and every stage moves an existing directory to `<work>/attic/` first, so a
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
| T_cpu-node / T_backend | backend speedup, per case (the reference: the CPU run with the most ranks) |
| T(U) / T(R) | AMR speedup at equal finest h; ideal bound N(U) / N(R) |
| [T(U) N(F) / N(U)] / T(F) | AMR grid efficiency at constant cell work (ghosts, 2:1 seam fills, reflux) |
| T(R) / T(F) - 1 | regrid cost (marking, regrid, prolongation; the FNL host round trip) |
| step cost x devices / N | throughput, ns per cell-step |

The hypotheses H1-H4 these ratios test are in the issue.

## Files

| Path | What |
|---|---|
| `matrix.conf` | every parameter: per leg the grids, levels, block sizes, scheme, marker and tolerances, t_w, settle length, regrid period, scaling series; the window (10/110 steps, 3 repetitions, 2 % spread); `[leg.N.reduced]` the WSL plumbing sizes; `[weak.<case>]` the weak-scaling cases (phase 9) |
| `make_cases.py` | writes one directory per case (inputs, `case.conf` for run.sh, `README.md`) from `matrix.conf`; calls `verification/step/make_step.py`, `make_sphere.py`, `verification/turbulence/make_channel.py` |
| `make_sphere.py` | leg 2's input: shock over a sphere, immersed boundary, octree, uniform or runtime AMR |
| `run.sh` | the chain: `cases`, `warmup`, `settle`, `windows`, `chain` stages; a `run.info` per run |
| `repartition.py` | a restart written on N ranks rewritten for M ranks (phase 1) |
| `ini_set.py` | rewrites `section.key=value` options of an input (the per-run changes of run.sh) |
| `analyze.py` | step costs, spreads, N, throughputs, the ratios, the grind times and the strong and weak efficiencies -> `results/summary.{md,csv}` |
| `scaling_plot.py` | the grind-time curves g(cells per device) from `summary.csv` (matplotlib, not needed by `analyze.py`) |
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
resubmit with `-x <node>`). Leg 3 with a surviving T3 dump: `sbatch jobs/medooza-fnl.slurm warmup --cases 3
--use-dump /ifs/.../resolved-hybrid` (the directory holding `channel-restart.*`).

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

**The step series is memory-capped 2-D, in the launch-bound regime (H1): not a grind measurement.** An A30 is
launch-bound below a few M cells (`profile-couette3d`: `small`, 2.1 M, is the launch-overhead profile; `medium`,
16.8 M, loads an A30), and the step cannot go above ~1.6 M cells per A30 (the memory ceiling below). Its block shape
also changes with m (realm C stores 12.2 cells per cell at m 8, 9.8 at m 32: the ghost share changes along the series),
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
# step (memory-capped, launch-bound); 1a at np 4 is phase 6's
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
from t_w) is untested here; a 4-A30 window of a cold-start L4 sphere next to 2a's would measure it in one job.

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
