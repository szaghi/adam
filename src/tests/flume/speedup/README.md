# FLUME speedup matrix

How fast a FLUME step is on each backend: CPU (gfortran + MPI), FNL OpenACC on NVIDIA GPUs, and FNL OpenMP offload on
AMD GPUs. The matrix also measures what dynamic AMR saves against a uniform grid of the same finest resolution, how
the cost changes with the 2-D plane orientation and the tree, and how it changes from 2-D to 3-D. Every case is
timed over a window of at most 100 steps on a developed flow, so a measurement takes minutes. The method is designed
so that start-up and output do not bias the window.

## Matrix

Two flows, both the Mach 2 shock over a body of the V6 verification (`verification/shock-cylinder`): the incident
shock starts at x = 0.2 and hits a body of radius 0.13 at the domain centre (an immersed boundary with an Euler
wall); compressible Euler, WENO-5 characteristic flux splitting, SSP-33, CFL 0.5, reflux on.

| Case directory | Dims | Plane, tree | AMR | Full simulation |
|---|---|---|---|---|
| `xy-quad-shock-cylinder-noamr` | 2-D | xy, quadtree | uniform at MAX_LEVEL | yes |
| `xy-quad-shock-cylinder-amr` | 2-D | xy, quadtree | dynamic, ∇ρ marker | yes |
| `xy-oct-shock-cylinder-noamr` | 2-D | xy, octree, z null | uniform | — |
| `xy-oct-shock-cylinder-amr` | 2-D | xy, octree, z null | dynamic | — |
| `xz-oct-shock-cylinder-noamr` | 2-D | xz, octree, y null | uniform | — |
| `xz-oct-shock-cylinder-amr` | 2-D | xz, octree, y null | dynamic | — |
| `yz-oct-shock-cylinder-noamr` | 2-D | yz, octree, x null | uniform | — |
| `yz-oct-shock-cylinder-amr` | 2-D | yz, octree, x null | dynamic | — |
| `xyz-shock-sphere-noamr` | 3-D | octree | uniform | yes |
| `xyz-shock-sphere-amr` | 3-D | octree | dynamic | yes |

Each axis plays the same role in every orientation: the shock runs along the first axis of the plane, the walls sit
on the second, and the null axis is extrapolated. Each case directory has its `README.md` (the case, its grid, its
toolchain), its `input.ini`, `data/<leg>-np<N>/` (the measurement after the warmup) and, for the four reference cases,
`full/<leg>-np<N>/` (the complete simulation to t = 0.25).

**What is compared, and why so:**

- **AMR speedup:** each noamr case is uniform at the finest level of its amr twin, so the ratio is the cost of the
  same finest resolution with and without dynamic AMR. The amr cases refine on one gradient-of-density marker (mode 2,
  conservative ρ). Its tolerance is absolute, so it is derived from the grid: `0.1·Δρ/h(base)`, with Δρ the incident
  jump. A regrid runs every 10 steps, with `iters = max - base` sweeps and conservative prolongation.
- **Orientation:** the quadtree refines x and y only (`adam_tree_object.F90`: z of a quadtree has one child), so a
  2-D run in the xz or yz plane must be an octree with a null axis. xy, xz and yz are therefore all run on octrees,
  with identical cell counts, so that they differ only in which array index is null. `xy-quad` is the extra point
  that prices the octree: an octree splits its null axis at every level too (4·2^L cells per column against 1).
- **2-D against 3-D:** the cost of a cell-step (ns per cell per step) of the sphere against the cylinder, on the same
  flow.
- **Backends:** the step cost of a reference label (by default the CPU run with the most ranks) divided by that of
  every other label, case by case.

**Dropped:** the forward-facing step (EV-step). It is a three-realm forest, and FLUME refuses runtime regridding on a
multi-realm forest (`forest_object%check_runtime_amr`, issue #74 D-M5-4; refined blocks along an inter-realm seam are
refused too), so it has no dynamic-AMR variant. Owner decision 2026-10-21: drop every step case for now.

## The grid is not decided yet

`matrix.conf` holds every parameter. The block size and the levels are `TBD`, and so are the values derived from them
(`delta_fine`, `delta_coarse`, the marker tolerance, `iters`): fill `CYL_NI`, `CYL_BASE_LEVEL`, `CYL_MAX_LEVEL`,
`SPH_NI`, `SPH_BASE_LEVEL`, `SPH_MAX_LEVEL`, then run

```bash
./make_inputs.py inputs
```

`run.sh` refuses an input that still holds `TBD`. Sizing notes:

- cells of a noamr case: quadtree (NI·2^L)², octree 2-D (NI·2^L)²·4·2^L, sphere (NI·2^L)³ (the case README
  prints them);
- NI even and ≥ 6 (2 ngc along a refined axis, `check_amr_block_cells`);
- the window must not be dominated by launch overhead on the GPU: aim at ≥ 10⁶ cells per GPU for the noamr cases;
- on medooza the FNL jobs run on `dev-gn` (1 h limit): split them with `--cases` if the warmups are long.

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
| `fnl-amd` | `flume-fnl-omp-amd` (OpenMP offload) | Giacomo Rossi, AMD cluster | `./build.sh fnl-amd --varset "local_amd gfx90a"` (his HDF5 varset and GPU arch) |

`build.sh` cleans the mode, builds it, and copies the executable to `exe/<leg>/`, because both FNL legs build
`exe/adam_flume_fnl`. It writes `builds/<leg>/build.txt` (the full fobis log) and `builds/<leg>/toolchain.md`
(compiler, MPI and fobis versions, varset, first compile line, commit, executable digest), which every case README
embeds.

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

`run.sh` options: `--cases` (substrings: `cylinder`, `sphere`, `-amr`, `xz-oct`, …), `--what data|full|all`,
`--reuse-warmup` (re-measure on an existing dump), `--exe`. The label of a measurement is `<leg>-np<N>`, so several
rank counts sit side by side; for the backend comparison prefer whole nodes (CPU: all cores; FNL: all GPUs) and one
GPU.

**For Giacomo (fnl-amd):** build with your AMD varset, run `./run.sh fnl-amd --np <GPUs> --launcher "<your launcher>
-n"` (all cases, or by `--cases`), then `./analyze.sh`. Send back `builds/fnl-amd/` and `*/data/fnl-amd-np*/`
(timings, logs, run.info, result.md; the restart dumps are git-ignored).

## Files

```
matrix.conf        every parameter (grid TBD, times, window, AMR cadence and tolerance factor)
make_inputs.py     writes <case>/input.ini and README.md; helpers for run.sh (stdlib only, Python >= 3.7)
build.sh           build one leg from clean, keep its executable, record its toolchain
run.sh             warmup + windows (data/) and full simulations (full/) of one leg
analyze.sh/.py     step costs, throughputs and every comparison -> results/summary.{md,csv}
jobs/              medooza Slurm templates (CPU, FNL)
builds/<leg>/      build.txt, toolchain.md
results/           summary.md, summary.csv
```

The root `.gitignore` ignores `*.sh` and `*.log`: add the scripts with `git add -f`. The build logs are named
`build.txt` so that they are tracked. Restart dumps, checkpoints and executables are ignored (`.gitignore` here).

**First on record:** an immersed boundary on a quadtree (`xy-quad`). V6 and its runtime-AMR variants all ran on
octrees. Before trusting `xy-quad` timings, compare its full simulation with the `xy-oct` solution, for example with
`verification/mhd/quadtree/quadtree_oracle.py`.
