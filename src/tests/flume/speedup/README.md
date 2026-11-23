# FLUME speedup matrix

**Status: harness ready, no leg measured yet** ([issue #133](https://github.com/szaghi/adam/issues/133), the plan's only
source of truth: every step is ticked there with its commit, every phase ends with a comment carrying its results).
Phase 3 (this harness) is in place and tested on WSL at reduced sizes; the legs run in phases 5-7, the summary and the
figures come in phase 8.

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
   replaces the run (`run.sh warmup --cases 3 --use-dump DIR` links its files).
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
- throughput = `step cost x ranks / N`, in ns per cell-step (ranks = devices).

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
| `matrix.conf` | every parameter: per leg the grids, levels, block sizes, scheme, marker and tolerances, t_w, settle length, regrid period; the window (10/110 steps, 3 repetitions, 2 % spread); `[leg.N.reduced]` the WSL plumbing sizes |
| `make_cases.py` | writes one directory per case (inputs, `case.conf` for run.sh, `README.md`) from `matrix.conf`; calls `verification/step/make_step.py`, `make_sphere.py`, `verification/turbulence/make_channel.py` |
| `make_sphere.py` | leg 2's input: shock over a sphere, immersed boundary, octree, uniform or runtime AMR |
| `run.sh` | the chain: `cases`, `warmup`, `settle`, `windows`, `chain` stages; a `run.info` per run |
| `repartition.py` | a restart written on N ranks rewritten for M ranks (phase 1) |
| `ini_set.py` | rewrites `section.key=value` options of an input (the per-run changes of run.sh) |
| `analyze.py` | step costs, spreads, N, throughputs and the ratios -> `results/summary.{md,csv}` |
| `jobs/medooza-{cpu,fnl}.slurm` | the medooza jobs: one run.sh stage per job, outputs on `/ifs` |
| `build.sh` | builds one backend from a clean tree into `exe/<backend>/` (it still calls backends "legs") |
| `templates/sphere-legacy.ini` | the removed matrix's sphere input, a record (stale, does not run) |
| `profile-couette3d/` | the profiling case (outside the matrix; shares the executables of `build.sh`) |

The run tree (`--work`, default `runs/full` or `runs/reduced` here):

```
<work>/cases/<case>/                         inputs, case.conf, README.md (rewritten by every run.sh call)
<work>/dumps/leg<N>/{warmup,settle}/         the two dumps of a leg, with the logs of the runs that wrote them
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
