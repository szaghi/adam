# FLUME speedup matrix

How fast a FLUME step is on each backend: CPU (gfortran + MPI), FNL OpenACC on NVIDIA GPUs, and FNL OpenMP offload on
AMD GPUs. The matrix also measures dynamic AMR against a fixed uniform grid, how the cost changes with the 2-D plane orientation and the tree, and how it changes from 2-D to 3-D. Every case is
timed over a window of at most 100 steps on a developed flow, so a measurement takes minutes. The method is designed
so that start-up and output do not bias the window.

**Leg**: one build of FLUME for one backend and compiler, named by a short word: `cpu` (gfortran + MPI), `fnl-nvf`
(OpenACC on NVIDIA GPUs) and `fnl-amd` (OpenMP offload on AMD GPUs). Every case is measured on every leg; the table of
the legs, with the command that builds each, is in [Legs and how to run them](#legs-and-how-to-run-them).

## Matrix

Two flows, both the Mach 2 shock over a body of the V6 verification (`verification/shock-cylinder`): the incident
shock starts at x = 0.2 and hits a body of radius 0.13 at the domain centre (an immersed boundary with an Euler
wall); compressible Euler, WENO-5 characteristic flux splitting, SSP-33, CFL 0.5, reflux on.

| Case directory | Dims | Plane, tree | AMR | Full simulation |
|---|---|---|---|---|
| `xy-quad-shock-cylinder-noamr` | 2-D | xy, quadtree | uniform at NOAMR_LEVEL | yes |
| `xy-quad-shock-cylinder-amr` | 2-D | xy, quadtree | dynamic, ∇ρ marker, BASE..MAX | yes |
| `xy-oct-shock-cylinder-noamr` | 2-D | xy, octree, z null | uniform | — |
| `xy-oct-shock-cylinder-amr` | 2-D | xy, octree, z null | dynamic | — |
| `xz-oct-shock-cylinder-noamr` | 2-D | xz, octree, y null | uniform | — |
| `xz-oct-shock-cylinder-amr` | 2-D | xz, octree, y null | dynamic | — |
| `yz-oct-shock-cylinder-noamr` | 2-D | yz, octree, x null | uniform | — |
| `yz-oct-shock-cylinder-amr` | 2-D | yz, octree, x null | dynamic | — |
| `xyz-shock-sphere-noamr` | 3-D | octree | uniform | yes |
| `xyz-shock-sphere-amr` | 3-D | octree | dynamic | yes |
| `xy-quad-heavy-shock-cylinder-noamr` | 2-D | xy, quadtree, NI 64 | uniform, 16.8 M cells | — |
| `xy-quad-heavy-shock-cylinder-amr` | 2-D | xy, quadtree, NI 64 | dynamic | — |
| `xy-quad-amd-shock-cylinder-noamr` | 2-D | xy, quadtree, NI 128 | uniform, 67 M cells; **fnl-amd only** | — |
| `xy-quad-amd-shock-cylinder-amr` | 2-D | xy, quadtree, NI 128 | dynamic; **fnl-amd only** | — |
| `xyz-amd-shock-sphere-noamr` | 3-D | octree, NI 32 | uniform, 134 M cells; **fnl-amd only** | — |
| `xyz-amd-shock-sphere-amr` | 3-D | octree, NI 32 | dynamic, ~330 M cells; **fnl-amd only** | — |

The first ten cases are the base matrix (every leg). The two `heavy` cases put the 2-D production path, the quadtree,
at a size that loads a GPU (the base quadtree is too small for that), on every leg. The four `amd` cases are sized
for the 192 GB AMD GPUs and run on `fnl-amd` only (`AMD_ONLY_FLOWS` in `matrix.conf`; `run.sh` skips them on the other
legs unless `--cases` names them).

Each axis plays the same role in every orientation: the shock runs along the first axis of the plane, the walls sit
on the second, and the null axis is extrapolated. Each case directory has its `README.md` (the case, its grid, its
toolchain), its `input.ini`, `data/<leg>-np<N>/` (the measurement after the warmup) and, for the four reference cases,
`full/<leg>-np<N>/` (the complete simulation to t = 0.25).

**What is compared, and why so:**

- **AMR:** the noamr level is independent of the amr levels (owner decision 2026-10-21). The noamr cases are a
  fixed-grid reference sized to run on every leg (L4: up to 16.8 M cells, which fits 4 A30); the amr cases go one
  level deeper (L2 → L5, up to about 15 M cells) so that they load the GPUs. The two do not compute the same
  solution, so the equal-resolution AMR speedup is an **estimate**: the noamr cost per cell-step times the cells of a
  uniform grid at the amr finest level, over the measured amr step. It assumes the throughput does not depend on the
  grid size, and on a GPU a larger grid usually runs faster per cell, so the estimate overstates the AMR gain there.
  The amr cases refine on one gradient-of-density marker (mode 2, conservative ρ). Its tolerance is absolute, so it is
  derived from the grid: `0.1·Δρ/h(base)`, with Δρ the incident jump. A regrid runs every 10 steps, with
  `iters = max - base` sweeps and conservative prolongation.
- **Orientation:** the quadtree refines x and y only (`adam_tree_object.F90`: z of a quadtree has one child), so a
  2-D run in the xz or yz plane must be an octree with a null axis. xy, xz and yz are therefore all run on octrees,
  with identical cell counts, so that they differ only in which array index is null. `xy-quad` is the extra point
  that prices the octree: an octree splits its null axis at every level too (4·2^L cells per column against 1).
- **2-D against 3-D:** the cost of a cell-step (ns per cell per step) of the sphere against the cylinder, on the same
  flow.
- **Backends:** the step cost of a reference label (by default the CPU run with the most ranks) divided by that of
  every other label, case by case.
- **Throughput against size:** the noamr ns per cell-step of the same flow at growing sizes (quadtree 0.26 M, 16.8 M,
  67 M cells; sphere 16.8 M, 134 M cells): where each backend saturates, and a direct measure of how far the
  equal-resolution AMR estimate (which assumes a size-independent throughput) is from the truth.

**Dropped:** the forward-facing step (EV-step). It is a three-realm forest, and FLUME refuses runtime regridding on a
multi-realm forest (`forest_object%check_runtime_amr`, issue #74 D-M5-4; refined blocks along an inter-realm seam are
refused too), so it has no dynamic-AMR variant. Owner decision 2026-10-21: drop every step case for now.

## Grid

`matrix.conf` holds every parameter; after editing it, run `./make_inputs.py inputs` (`run.sh` refuses an input that
still holds `TBD`). The grid (2026-10-21), sized on medooza (64 cores, 1 TB, 4 A30 of 24 GB per node) and Giacomo's
node (8 AMD GPUs of 192 GB):

| Flow | NI | noamr | amr | noamr size | amr size (estimate) |
|---|---|---|---|---|---|
| cylinder | 32 | L4, h = 1/512 | L2 → L5, h = 1/128 → 1/1024 | quad 256 blocks, 0.26 M cells; oct 4096 blocks, 16.8 M cells, ~26 GB | quad ~160 blocks, ~0.16 M cells; oct ~3600 blocks, ~15 M cells, ~23 GB |
| sphere | 16 | L4, h = 1/256 | L2 → L5, h = 1/64 → 1/512 | 4096 blocks, 16.8 M cells, ~19 GB | ~2900 blocks, ~12 M cells, ~13 GB |
| quad heavy | 64 | L6, h = 1/4096 | L3 → L7, h = 1/512 → 1/8192 | 4096 blocks, 16.8 M cells, ~62 GB | ~700 blocks, ~2.9 M cells, ~11 GB |
| quad AMD | 128 | L6, h = 1/8192 | L2 → L7, h = 1/512 → 1/16384 | 4096 blocks, 67 M cells, ~226 GB | ~650 blocks, ~10 M cells, ~36 GB |
| sphere AMD | 32 | L4, h = 1/512 | L2 → L6, h = 1/128 → 1/2048 | 4096 blocks, 134 M cells, ~99 GB | ~10000 blocks, ~330 M cells, ~240 GB |

The base and heavy cases fit 4 A30 (the heavy noamr at ~15 GB per A30) and one AMD GPU. The AMD-only cases need the
AMD node: the sphere noamr fits one GPU (61 % of its budget) or eight (17 M cells each); the quad noamr and the sphere
amr need the node (~28-30 GB per GPU).

**The warmup sets the ceiling, not the memory.** dt scales with h, so the warmup to t = 0.15 takes ~0.9/h steps: 460
at h = 1/512, ~3700 at 1/4096, ~7400 at 1/8192, ~14700 at 1/16384. The quad heavy noamr on 4 A30 is the case that
may not fit medooza's 1 h `dev-gn` limit: run it alone (`--cases xy-quad-heavy-shock-cylinder-noamr`). **A 2-D amr case
cannot be made heavy:** its refined region is lines (the shock, the ring around the body), so its cells grow as 1/h
while its steps grow as 1/h too; 16 M amr cells would need h ~ 1/50000 and ~45000 warmup steps. The 2-D amr cases
stay light by nature; the noamr cases are the GPU-loading 2-D points.

The memory model is about 55 block-sized arrays per block, ghosts included on
every axis (null axes too, which makes the quadtree's 1-cell-thick blocks 7 cells thick in memory); the amr sizes are
geometric estimates (the planar incident shock is a full slab of finest blocks, plus a shell around the body and the
2:1 transition), and the regrid lines of the logs print the real ones. Two consequences to read the results with:

- the base quadtree cases (0.26 M and ~0.16 M cells) are too small to load a GPU: they price the octree at the same
  h; the heavy pair measures the quadtree's GPU throughput;
- the base cases on 8 AMD GPUs are ~2 M cells per GPU: those np 8 runs are strong scaling, np 1 is the GPU-to-GPU
  comparison; the AMD-only cases are the ones that load the AMD node.

Before the campaign, calibrate: run `xyz-shock-sphere-noamr` once with `SPH_NOAMR_LEVEL=3` (2 M cells) on the CPU node
and on 4 A30, to measure the real ns per cell-step and check the memory model. On medooza the FNL jobs run on `dev-gn`
(1 h limit): split them with `--cases`.

## Measurement protocol

FLUME has no step timer, and a run's wall time also holds its start-up (MPI, device set-up, the initial AMR, the
immersed-body distance function) and its output. `run.sh` therefore measures by difference, with no code change:

1. **Warmup:** the case runs to t = 0.15, the incident shock at mid-body, every output off but one restart dump at
   the end. Every leg starts its window from the same physical state, each from its own dump, because restart files
   are per rank.
2. **Window:** three times, interleaved (short, long, short, long, …), a run restarted from the dump for 10 steps
   (short) and one for 110 (long). Every output is off: `it_save = restart_save = residuals_save =
   conservation_history_save = 0`, so not even the last step writes (`is_to_save` treats cadence 0 as never).
3. **Cost:** `D_r = T_long - T_short` is the cost of steps 11 to 110 after the restart. Start-up, the restart read
   and the first-step transients (first device launches, first touch) cancel. The step cost is `median(D_r)/100`, and
   the spread `(max - min)/median` is reported, with a warning above 5 %.
4. **AMR:** the regrid cadence keys on the step count, which the restart restores, and N1 = 10 and NWIN = 100 are
   multiples of it. The two runs therefore differ by exactly 10 regrids, which belong to the cost of a dynamic-AMR
   step. FNL prints the host round trip of each regrid; its share of the window is reported.
5. **Comparable work:** 100 steps after t = 0.15 are a few percent of the run, so the grid and the solution barely
   move across the window, and the same window is measured on every leg. The analysis records the window's t-range,
   its mean dt and the mean number of blocks it ran on (from the restart and regrid lines). The throughput
   (ns per cell-step) therefore normalises differences in grid size.

**Known residual biases** (not removable without a code change; each is the same for both runs, so it enters the
window only per step):

- every rank prints 7 progress lines per step (`print_progress`): a stdout cost that grows with np. It is small
  against a GPU step, but on many CPU ranks with small blocks it is visible. Remedy if needed: print on rank 0 only;
- `compute_conservation` (a global reduction) runs every step whatever `conservation_history_save` says: a genuine
  per-step cost of FLUME today, so it is kept;
- the FNL regrid goes through the host: a genuine cost of FNL dynamic AMR, reported separately.

## Legs and how to run them

| Leg | fobis mode | Who, where | Build |
|---|---|---|---|
| `cpu` | `flume-cpu-gnu` (MPI only) | medooza | `source /ifs/hpc/home/szaghi/gnu/gnu.sh && ./build.sh cpu --varset iac_gnu` |
| `fnl-nvf` | `flume-fnl-nvf` (OpenACC) | medooza, A30 | `source /ifs/hpc/home/szaghi/nvidia/nvhpc.sh && ./build.sh fnl-nvf --varset iac_nvf` |
| `fnl-amd` | `flume-fnl-omp-amd` (OpenMP offload) | Giacomo Rossi, AMD node (8 × 192 GB) | `./build.sh fnl-amd --varset "thera_amd gfx942"` (see [his section](#instructions-for-giacomo-fnl-amd-leg)) |

`build.sh` cleans the mode, builds it, and copies the executable to `exe/<leg>/`, because both FNL legs build
`exe/adam_flume_fnl`. It writes `builds/<leg>/build.txt` (the full fobis log) and `builds/<leg>/toolchain.md`
(compiler, MPI and fobis versions, varset, compile and link flags, commit, executable digest), which every case README
embeds. `--jobs N` compiles N files at a time. `--tag NAME` (with `--gpu OPTS` for extra nvfortran `-gpu=` options)
builds a variant into `exe/<leg>-<NAME>/` for profiling experiments; a variant is never a matrix leg (see
[profile-couette3d](profile-couette3d/README.md#build-variants)).

```bash
cd src/tests/flume/speedup
# medooza, CPU: one full node, all ranks
sbatch --ntasks=<cores per node> jobs/medooza-cpu.slurm                # measurements (data/)
sbatch --ntasks=<cores per node> jobs/medooza-cpu.slurm --what full    # the four full simulations
# medooza, FNL: 4 A30 (or --gres=gpu:1 --ntasks-per-node=1 for one GPU)
sbatch jobs/medooza-fnl.slurm --cases cylinder
sbatch jobs/medooza-fnl.slurm --cases sphere
# AMD (Giacomo): any launcher, one rank per GCD/GPU
./run.sh fnl-amd --np 8 --launcher "srun -n"          # or "mpirun -np"
# anywhere, after any run
./analyze.sh                                         # results/summary.md, data/*/result.md, README results
```

`run.sh` options: `--cases` (substrings: `cylinder`, `sphere`, `-amr`, `xz-oct`, …), `--skip` (substrings to leave out), `--what data|full|all`,
`--reuse-warmup` (re-measure on an existing dump), `--exe`. The label of a measurement is `<leg>-np<N>`, so several
rank counts sit side by side; for the backend comparison prefer whole nodes (CPU: all cores; FNL: all GPUs) and one
GPU.

The `fnl-amd` leg has its own step-by-step section below: [Instructions for Giacomo](#instructions-for-giacomo-fnl-amd-leg).

## Instructions for Giacomo (fnl-amd leg)

Giacomo, this section is everything you need to run the AMD leg on your node (8 AMD GPUs of 192 GB; MI300X assumed,
`gfx942`). You will build one executable, run every case of the matrix twice (on the whole node and on one GPU),
analyse, and send back a few small text files. Every command is run from `src/tests/flume/speedup/` unless stated.

### 1. Prerequisites

- **Compiler:** AMD flang with OpenMP offload (ROCm `amdflang`), and an **MPI wrapper (`mpif90`) built on it**.
- **GPU-aware MPI:** the FNL halo exchange posts MPI directly on device buffers (there is no host-staging fallback), so
  the MPI must be ROCm-aware (Open MPI + UCX with ROCm, or Cray MPICH with `MPICH_GPU_SUPPORT_ENABLED=1`). Without it
  every run with more than one rank fails at the first ghost exchange.
- **HDF5:** parallel HDF5 built with the same `amdflang` and MPI. Point a fobos varset at it: `thera_amd` if your
  node is thera, or `local_amd` after building it into `lib/hdf5/develop/amd/latest` with
  `scripts/hdf5_build.sh -build -sdk ... -hdf5 $PWD/lib/hdf5/develop/amd/latest` (from the repository root; see
  `docs/start/install.md`).
- **FoBiS** ≥ 3.8 (`pip install FoBiS.py`, the command is `fobis`) and **python3 ≥ 3.7** (standard library only).

### 2. Build

```bash
cd <adam>                       # repository root
git pull && fobis fetch         # current sources and dependencies
cd src/tests/flume/speedup
./build.sh fnl-amd --varset "thera_amd gfx942"    # or "local_amd gfx942"; gfx90a for MI250X
```

`build.sh` cleans the `flume-fnl-omp-amd` mode, builds it, copies the executable to `exe/fnl-amd/adam_flume_fnl`, and
writes `builds/fnl-amd/build.txt` (the full fobis log) and `builds/fnl-amd/toolchain.md` (compiler, MPI, FoBiS, commit,
executable digest). Check that `toolchain.md` lists your `amdflang` and `mpif90`. A failed build prints the path of the
log; send it to Stefano as it is.

### 3. GPUs and ranks

Run **one MPI rank per GPU**: FLUME picks the device `mod(node-local rank, visible GPUs)` by itself, so do not bind
devices per rank (if your launcher does, e.g. `--gpus-per-task=1`, that works too). Never put two ranks on one GPU:
each rank pre-allocates 85 % of its GPU's memory for its blocks pool, and two of them overflow it.

`run.sh` launches with `<launcher> <np> <executable> input.ini`; the default launcher is `mpirun -np`. Pass yours with
`--launcher`, e.g. `--launcher "srun -n"` or `--launcher "mpirun --bind-to numa -np"`.

### 4. Smoke test and calibration first (about 10 minutes)

The times below are estimates (no FLUME run has been measured on MI300X yet). Two short runs check the build and the
MPI, and measure the real throughput of one loaded GPU:

```bash
./run.sh fnl-amd --np 8 --cases xy-quad-shock-cylinder-noamr    # ~2 min: smoke test (8 ranks, GPU-aware MPI)
./run.sh fnl-amd --np 1 --cases xyz-shock-sphere-noamr          # ~7 min: calibration (16.8 M cells on one GPU)
./analyze.sh
cat xyz-shock-sphere-noamr/data/fnl-amd-np1/result.md
```

Each must end with `run.sh: fnl-amd-np<N> done`. The first case is too small to load a GPU, so take the throughput
from the second: `result.md` gives its step cost and its ns per cell-step. If that is far from the 20 ns assumed
below, rescale the table: a case costs about
`(warmup steps + 360) × step cost + 7 × start-up`, where the step cost is `cells × ns per cell-step / GPUs` (but not
below ~10 ms, the launch floor) and a start-up is 15 s to 2 min (it grows with the memory a case allocates).

### 5. The cases, what each is for, and what they cost

Each case runs 7 times: a warmup to t = 0.15 (it ends with a restart dump), then 3 × (10 + 110) restarted steps with
all output off. "Warmup steps" is about 0.9/h at the finest spacing h. The **estimated** wall times assume ~20 ns per
cell-step on a loaded MI300X, a ~10 ms per-step floor, and the host round trip of every regrid (every 10 steps); the
amr cell counts are geometric estimates at the peak of the run, so their times are upper bounds. Trust them within a
factor of ~3 until the calibration says otherwise.

| Case | What it is for | Cells | Warmup steps | np 8 | np 1 |
|---|---|---|---|---|---|
| `xy-quad-shock-cylinder-noamr` | base 2-D reference on the quadtree; small (launch-bound on a GPU) | 0.26 M | 460 | ~2 min | ~2 min |
| `xy-quad-shock-cylinder-amr` | base 2-D dynamic AMR, quadtree | ~0.16 M | 920 | ~2 min | ~2 min |
| `xy-oct-shock-cylinder-noamr` | 2-D on an octree (z null): the octree's cost against the quadtree | 16.8 M | 460 | ~4 min | ~8 min |
| `xy-oct-shock-cylinder-amr` | the same with dynamic AMR | ~15 M | 920 | ~5 min | ~11 min |
| `xz-oct-shock-cylinder-noamr` | orientation: the same flow in the xz plane (y null) | 16.8 M | 460 | ~4 min | ~8 min |
| `xz-oct-shock-cylinder-amr` | orientation, dynamic AMR | ~15 M | 920 | ~5 min | ~11 min |
| `yz-oct-shock-cylinder-noamr` | orientation: the yz plane (x null) | 16.8 M | 460 | ~4 min | ~8 min |
| `yz-oct-shock-cylinder-amr` | orientation, dynamic AMR | ~15 M | 920 | ~5 min | ~11 min |
| `xyz-shock-sphere-noamr` | base 3-D reference; the CPU/A30/AMD comparison | 16.8 M | 230 | ~4 min | ~7 min |
| `xyz-shock-sphere-amr` | base 3-D dynamic AMR | ~12 M | 460 | ~4 min | ~8 min |
| `xy-quad-heavy-shock-cylinder-noamr` | the quadtree at a size that loads a GPU (2-D throughput) | 16.8 M | 3700 | ~6 min | ~27 min |
| `xy-quad-heavy-shock-cylinder-amr` | its dynamic-AMR twin (light: 2-D AMR refines lines) | ~2.9 M | 7400 | ~4 min | ~11 min |
| `xy-quad-amd-shock-cylinder-noamr` | **AMD-only**: 2-D quadtree filling the node | 67 M | 7400 | ~29 min | does not fit |
| `xy-quad-amd-shock-cylinder-amr` | **AMD-only**: deep 2-D AMR (levels 2 → 7) | ~10 M | 14700 | ~12 min | ~66 min |
| `xyz-amd-shock-sphere-noamr` | **AMD-only**: 3-D, one GPU nearly full (99 GB) | 134 M | 460 | ~12 min | ~44 min |
| `xyz-amd-shock-sphere-amr` | **AMD-only**: deep 3-D AMR filling the node (~240 GB) | ~330 M | 1840 | ~51 min | does not fit |
| **Total** | | | | **~2.5 h** | **~3.7 h** (without the two that do not fit) |

What the numbers are used for: the base cases compare the three legs on identical work (medooza's 64 CPU cores, 4 A30,
your node); np 8 against np 1 is your strong scaling, np 1 the GPU-to-GPU comparison with one A30; the heavy and AMD
cases show the throughput at sizes that load a GPU (the "throughput against size" table of the results).

### 6. Run one case

Give its exact directory name to `--cases` (each full name selects exactly one case):

```bash
./run.sh fnl-amd --np 8 --cases xyz-amd-shock-sphere-noamr
./run.sh fnl-amd --np 1 --cases xyz-amd-shock-sphere-noamr
```

Several at once: `--cases name1,name2`. Substrings select groups: `--cases sphere` (the four sphere cases),
`--cases -amr` (every amr case), `--cases amd-` (the four AMD-only cases). A rerun of a case replaces its
`data/fnl-amd-np<N>/`; add `--reuse-warmup` to re-measure on the existing warmup dump without running the warmup again.

### 7. Run the whole matrix

```bash
nohup ./run.sh fnl-amd --np 8 > run-amd-np8.txt 2>&1 &                       # all 16 cases, ~2.5 h
# when it has finished (one run at a time: the two would share the GPUs):
nohup ./run.sh fnl-amd --np 1 --skip xy-quad-amd-shock-cylinder-noamr,xyz-amd-shock-sphere-amr \
   > run-amd-np1.txt 2>&1 &                                                    # 14 cases, ~3.7 h
./analyze.sh
```

Inside a batch system, put the same lines in a job script (one node, 8 GPUs, exclusive). The run continues after a
failed case and lists the failures at the end; a failed run keeps its `log.txt`. The most likely failure is a case that
does not fit the memory: the log then shows a FUNDAL device allocation error (`ierr = 101`). Full simulations
(`--what full`) are not needed from your side: the reference simulations are run on medooza.

### 8. What to send back

```bash
cd <adam>/src/tests/flume/speedup
tar czf flume-speedup-amd.tgz --exclude='*.fbd' --exclude='*.tnd' --exclude='*.time' --exclude='*.h5' \
   --exclude='*.xdmf' --exclude='*.reference' builds/fnl-amd */data/fnl-amd-np* results run-amd-np*.txt
```

The archive holds the build log and toolchain, every `timings.csv`, `run.info`, `result.md` and run log, and the
summary: a few MB. The excludes leave out the restart dumps under each `warmup/`, which reach hundreds of GB for the
AMD-only cases; delete them once the archive is made (`find . -path '*/warmup/*' \( -name '*.fbd' -o -name '*.h5' \)
-delete`), unless you want `--reuse-warmup` later. Or commit the same files on a branch: the restart dumps and
the executable are git-ignored.

## Files

```
matrix.conf        every parameter (grid, times, window, AMR cadence and tolerance factor)
make_inputs.py     writes <case>/input.ini and README.md; helpers for run.sh (stdlib only, Python >= 3.7)
build.sh           build one leg from clean, keep its executable, record its toolchain
run.sh             warmup + windows (data/) and full simulations (full/) of one leg
analyze.sh/.py     step costs, throughputs and every comparison -> results/summary.{md,csv}
jobs/              medooza Slurm templates (CPU, FNL)
builds/<leg>/      build.txt, toolchain.md (builds/<leg>-<tag>/ for a variant)
results/           summary.md, summary.csv
profile-couette3d/ the profiling case (not a speedup case): a 3-D Couette channel, see its README.md
```

The root `.gitignore` ignores `*.sh` and `*.log`: add the scripts with `git add -f`. The build logs are named
`build.txt` so that they are tracked. Restart dumps, checkpoints and executables are ignored (`.gitignore` here).

**First on record:** an immersed boundary on a quadtree (`xy-quad`). V6 and its runtime-AMR variants all ran on
octrees. Before trusting `xy-quad` timings, compare its full simulation with the `xy-oct` solution, for example with
`verification/mhd/quadtree/quadtree_oracle.py`.
