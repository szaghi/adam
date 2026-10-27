# profile-couette3d

A case for **profiling** FLUME and the ADAM library, not for speedup measurements or verification: one short,
steady, reproducible run that goes through the main regions of the code, on WSL, on medooza and on the AMD node. It is
outside the speedup matrix (`../run.sh` and `../analyze.sh` ignore it) and shares only its executables
(`../build.sh <leg>`).

## The flow

The exact compressible Couette flow of `regression/couette` (issue #65, VV-3), extended to 3-D: an isothermal no-slip
wall at y = 0 (T = 1), a no-slip adiabatic wall moving at u = 1 at y = 1, periodic in x and z; γ = 1.4, μ = 0.01,
Pr = 0.75. The initial condition is the exact steady solution, so the flow does not evolve: every profile of a given
size sees the same work, whatever the step.

## What a profile goes through

| Region | Where | How the case exercises it |
|---|---|---|
| Convective residual | WENO-5 flux splitting on characteristic variables (`compute_residuals_weno`, FLUME cpu/fnl objects and Euler kernels) | x, y and z (3-D octree, nothing null) |
| Dissipative residual | viscous stress and heat flux, 4th order (`adam_flume_dissipation_library`, cpu/fnl dissipation kernels) | μ = 0.01: the viscous time step bound is of the order of the convective one |
| Time integration | SSP Runge–Kutta 5 stages (`adam_rk_object`, FNL `rk_fnl`) | stage buffers, a ghost exchange per stage |
| Boundary conditions | `adam_flume_bc_object`, `set_boundary_conditions` | isothermal and moving no-slip walls, periodic x and z, edges and corners |
| Ghost exchange | `update_ghost` (local copies and MPI, GPU-direct on FNL) | every stage, across ranks when np > 1 |
| AMR seams and reflux | 2:1 seam ghost fill (restriction, prolongation), the forest flux register (`adam_forest_object`) | the block layer along each wall is one level finer (two box markers): 2:1 seams parallel to the walls, reflux every step |
| Runtime regrid | `regrid_forest` (marking, ghost refresh, the library regrid; FNL: the host round trip) | every 10 steps; the box markers vote "keep", so the grid does not change |
| Reductions | the global time step, the conservation integrals | every step |
| Output (optional) | HDF5 checkpoint (`save_xh5f`) | `--io`: the first and last step |

**Not exercised:** the immersed boundary (distance function, eikonal extrapolation, wall inversion), the prolongation
of a regrid that changes the grid, inter-realm (multi-realm) seams, MHD, gravity, the positivity limiter. Profile
those on the speedup cases (`../xyz-shock-sphere-amr`, for example, has the immersed boundary and real regrids).

## Two variants

- **Refined walls** (default): everything in the table above.
- **Uniform** (`--uniform`): the simplest run, one realm on a uniform grid with no AMR marker, no regrid and no 2:1
  seam. It keeps the convective and dissipative residuals, the stages, the boundaries, the ghost exchange and the
  reductions, so its profile is the solver alone; set against the refined-walls profile of a similar size, the
  difference prices the AMR machinery (seam fills, reflux, regrids).

| Uniform size | Grid | Cells | Memory (estimate) | For |
|---|---|---|---|---|
| `small` | NI 16, level 3: 512 blocks | 2.1 M | ~2.8 GB | WSL (CPU np 2, one GPU) |
| `medium` | NI 32, level 3: 512 blocks | 16.8 M | ~15 GB | one A30, the medooza CPU node |
| `large` | NI 32, level 4: 4096 blocks | 134 M | ~117 GB | one AMD GPU (192 GB) |

## Sizes of the refined-walls variant

| Size | Blocks | Cells | Memory (estimate) | For |
|---|---|---|---|---|
| `small` | NI 16, base level 2, walls level 3: 288 | 1.2 M | ~1.6 GB | WSL (CPU np 2, one GPU), one A30 |
| `medium` | NI 16, base level 3, walls level 4: 1408 | 5.8 M | ~8 GB | the medooza CPU node (np 64), one or four A30 |
| `large` | NI 32, base level 3, walls level 4: 1408 | 46 M | ~40 GB | one AMD GPU (192 GB); four A30 (~10 GB each) |

The size matters on a GPU: a case that is small for the device leaves it waiting on kernel launches, and its profile
then shows launch latency rather than kernel cost. Profile the kernels on a size that loads the GPU (`medium` on an
A30, `large` on an AMD GPU); `small` on a big GPU is a profile of the launch overhead. The committed
`input-<size>.ini` and `input-<size>-uniform.ini` are what `profile.sh` writes (20 steps, a regrid every 10 on the
refined-walls variant, no output).

Wall times of the 20 steps, start-up included, measured on medooza (commit `f46b56fb`): `medium` on one A30, 50 s
refined walls and 89 s uniform under Nsight Systems; `medium` uniform on 16 ranks of a CPU node, 427 s without a
profiler. On WSL `small` takes a few minutes on 2 CPU ranks and tens of seconds on one GPU (a planning figure: WSL is
not a benchmark machine). `large` has not been run yet. Profilers add their own overhead: Nsight Systems and rocprofv3
tens of percent; Nsight Compute replays each profiled kernel many times, so `profile.sh` limits it to 300 kernels of
rank 0.

## Run

Build the leg first (`cd .. && ./build.sh <leg> --varset ...`), then from this directory:

```bash
./profile.sh <leg> --size <small|medium|large> [--uniform] --np <N> --tool <none|nsys|ncu|rocprof|perf> [--steps N] [--io]
             [--tag NAME] [--exe PATH]
```

For example, the uniform variant on WSL: `./profile.sh cpu --size small --uniform --np 2 --tool perf`.

Every run goes to `runs/<leg>-<size>[-uniform]-np<N>-<tool>[-<tag>]/` (a previous one is renamed, not overwritten):
`input.ini`, `log.txt`, `wall.txt`, `run.info` (host, job, command, executable and its digest, commit, `src_modified`)
and the profiler files, one per rank. `ADAM_SEAM_TIMING=1` is always set, so the end of `log.txt` carries the forest
timing report (loop, steps, reflux, seam fills; max over the ranks), which `profile.sh` prints. `--dry-run` prints the
command without running it.

**WSL** (one run at a time; FNL at most np 2, the two GPUs; the UCX knobs of issue #12, as in
`regression/run-fnl-local.sh`):

```bash
./profile.sh cpu     --size small --np 2 --tool perf      # then: perf report -i runs/cpu-small-np2-perf/perf-r0.data
./profile.sh fnl-nvf --size small --np 1 --tool nsys      # then: nsys stats runs/fnl-nvf-small-np1-nsys/nsys-r0.nsys-rep
```

Nsight Compute on WSL needs the GPU performance counters enabled in the Windows NVIDIA Control Panel
(`ERR_NVGPUCTRPERM` otherwise); compute-sanitizer is blind there.

**medooza** (inside an allocation, the toolchain sourced; `srun` forwards the environment):

```bash
source /ifs/hpc/home/szaghi/nvidia/nvhpc.sh
salloc -A iacroma -p dev-gn --gres=gpu:4 --ntasks-per-node=4 --mem=96G -t 00:30:00
./profile.sh fnl-nvf --size medium --np 1 --tool nsys --launcher "srun --mpi=pmix -n"
./profile.sh fnl-nvf --size large  --np 4 --tool nsys --launcher "srun --mpi=pmix -n"
# CPU node: source /ifs/hpc/home/szaghi/gnu/gnu.sh; export OMPI_MCA_coll_hcoll_enable=0
./profile.sh cpu     --size medium --np 64 --tool none --launcher "srun --mpi=pmix --cpu-bind=cores -n"
```

Two tools do not run on medooza (October 2026): **Nsight Compute** stops with `ERR_NVGPUCTRPERM` (no permission to
access the GPU performance counters: the driver reserves them to the administrators, presumably through its
`NVreg_RestrictProfilingToAdminUsers` option), and **perf** is not installed (`which perf` finds nothing on the login
node nor on a compute node). Without counters, the registers and the stack (local memory) of every kernel are in the
executable (`cuobjdump --dump-resource-usage <exe>`) and in the Nsight Systems trace (`nsys stats --report
cuda_gpu_trace`, column `Reg/Trd`); for the counters themselves use Nsight Compute on WSL (second GPU), keeping in mind
that its consumer GPU runs FP64 at 1/64 of the FP32 rate, so its kernels are bound by FP64 throughput where an A30 is
not.

### Build variants

A profiling experiment that changes how the code is compiled builds a **variant** of a leg, so that the leg's
executable and toolchain record stay as they are and the two can be profiled side by side:

```bash
cd .. && ./build.sh fnl-nvf --varset iac_nvf --jobs 16 --tag lto --gpu lto   # -> exe/fnl-nvf-lto/, builds/fnl-nvf-lto/
cd profile-couette3d
./profile.sh fnl-nvf --size medium --np 1 --tool nsys --tag lto --launcher "srun --mpi=pmix -n"
```

`--gpu OPTS` (the `fnl-nvf` leg only, with `--tag` and a clean build) appends `-gpu=OPTS` to the compile and link flags
of the fobis mode: fobis appends `--cflags`/`--lflags` to the mode's own, and nvfortran merges the two `-gpu=` options.
The fobos files and the default build do not change; removing the variant is `rm -r ../exe/fnl-nvf-<tag>
../builds/fnl-nvf-<tag>`. A variant's `toolchain.md` records the full compile and link flags and is not embedded in the
speedup case READMEs. A source experiment (a patch) is built the same way, from a tree with the patch applied:
`toolchain.md` then marks the commit `src/ modified`.

**AMD node** (one rank per GPU; rocprofv3 from ROCm 6.2 or later):

```bash
./profile.sh fnl-amd --size large --np 1 --tool rocprof   # per-kernel statistics in runs/.../rocprof-r0/
./profile.sh fnl-amd --size large --np 8 --tool rocprof   # the same split over the node: communication shows up
```

## Reading the results

- `nsys stats <file>.nsys-rep` (or the Nsight Systems GUI): the kernel and memory-copy summaries, the OpenACC and
  MPI ranges; the timeline shows where the steps wait.
- `ncu-r0.ncu-rep` (Nsight Compute GUI or `ncu --import ... --page details`): occupancy, memory throughput, registers
  and spills of each kernel.
- `rocprof-r<N>/`: the `*_stats.csv` files rank the kernels by total time; the trace files hold every launch.
- `perf report -i perf-r<N>.data`: the CPU hot spots with their call graphs.
- The forest timing line: the wall time of the loop and of the steps; their difference is the regrid and post-step
  share, the reflux and seam entries price the AMR machinery.

## Results

Measured in October 2026 at commit `f46b56fb` (nvfortran 26.9 with CUDA 12.9, gfortran 16.2 with Open MPI 5.0.11). The
GPU numbers come from one A30 of medooza (`gnode01`, driver 575, 1440 MHz), one rank, the `medium` size, 20 steps, under
Nsight Systems; kernel times are the mean per launch of `nsys stats --report cuda_gpu_kern_sum`, and the loop is the
forest timing line (the steps plus the regrids, without start-up). Repeated runs of the same executable, in three jobs,
agree to 0.6% in loop time and kernel means; the wall time, start-up included, varies by up to 8%.

### Where the time goes

**FNL, A30, the build of the `fnl-nvf` leg:**

| | Refined walls | Uniform |
|---|---|---|
| `compute_face_fluxes_dev` (WENO-5 convective fluxes, one direction per launch) | 61.6 ms, 67% of the GPU time | 172.4 ms, 71% |
| `add_dissipative_fluxes_o4_dev` (viscous fluxes) | 21.0 ms, 23% | 55.3 ms, 23% |
| every other kernel | 10% | 6% |
| GPU busy (kernel time / loop) | 27.7 s / 35.6 s = 78% | 72.7 s / 72.6 s = 100% |

The two flux kernels are the whole cost of the uniform run. On the refined walls 8 s of the loop (22%) have no kernel
running (a lower bound: the kernel total includes the start-up kernels). At most 1.9 s of it is the regrid and
post-step share (loop minus steps); the other 6 s fall inside the steps and are not attributed yet. The 3708 launches
of `blend_inadmissible_ghosts_dev` (22 µs each on average, 82 ms in all) cannot fill them; an Nsight Systems timeline of
one refined step would locate them.

**CPU** (Nsight Systems sampling on WSL, 2 ranks; shares only, WSL is not a benchmark machine): the WENO weights 31%,
the dissipative fluxes about 20%, `malloc`/`free` 10.5% and the integer power `__powidf2` 2.9%, both inside the WENO
weights, the AMR machinery about 2%. The last two were avoidable: the locals of `adam_weno_object` had bounds depending
on `S`, so gfortran put them on the heap, and `(zeps + IS)**wexp` with a variable integer `wexp` called the integer
power of libgcc. Commit `f46b56fb` fixed both, bitwise (the products are formed in `__powidf2`'s order): `medium`
uniform on 16 ranks of a medooza CPU node went from 497 s to 427 s (three runs each, alternated; the "after" build
was the fix before its commit, so the run records carry the HEAD of the time, `ff80d8da`):

![CPU wall times before and after the WENO fix](figures/cpu-weno.png)

### The face-flux kernel on the A30

`cuobjdump --dump-resource-usage` gives `compute_face_fluxes_dev` 254 registers and a 3608-byte stack per thread: its
private arrays (`qs`, `qas`, `fsplit`, `er`, …, sized by `S_MAX = 5`) go to local memory, because their addresses are
passed to device routines of other modules (`face_split_fluxes`, the Roe average and eigenvectors, the WENO
reconstruction, the back-projection) that nvfortran compiles as calls (`-Minfo`: "Local memory used for qas, mu, flux,
fsplit, qs, vr_, vr, v, er"). Four variants were built (`../build.sh --tag`, see [Build variants](#build-variants);
Slurm job 1357634) and timed with the leg (built in job 1357614) in one job, 1357635:

| Variant | `-gpu=` added | Face kernel: registers, stack | Face kernel, uniform | Dissipative, uniform | Loop, uniform | Loop, refined |
|---|---|---|---|---|---|---|
| leg (base) | none | 254, 3608 B | 172.4 ms | 55.3 ms | 72.6 s | 35.6 s |
| `r128` | `maxregcount:128` | 128, 3704 B | 200.4 ms (+16%) | 53.9 ms | 81.2 s (+12%) | 39.1 s |
| `lto` | `lto` | 148, 2720 B | 103.1 ms (−40%) | 43.6 ms (−21%) | 47.3 s (−35%) | 26.2 s (−27%) |
| `lto-r128` | `lto,maxregcount:128` | 148, 2720 B | 103.6 ms | 44.2 ms | 47.6 s | 26.2 s |
| `B0-lto` | `lto`, and the patch below | 210, 2368 B | 56.3 ms (−67%) | 44.0 ms (−20%) | 33.4 s (−54%) | 21.2 s (−40%) |

and, from an earlier job (1357619), `B0` without LTO: 164.1 ms (−5%), loop 70.1 s uniform and 34.8 s refined. The figure
splits each loop into the two flux kernels, the other kernels and the time with no kernel running; its bars are the
means over every run of a variant (the leg ran four times in three jobs, `B0` twice), so they differ from the table,
which is job 1357635 alone, by at most 0.2%:

![The A30 loop of each build variant, split by kernel](figures/a30-variants.png)

- **Link-time optimisation is the largest single gain.** With `-gpu=lto` the device link inlines the device routines
  across modules (a two-module test: one `CALL` in the SASS without it, none with it), which removes the calls and much
  of the local memory they forced. Both flux kernels gain (the face one 40%, the dissipative one 21%), and the other
  kernels together 24% (4.37 s to 3.31 s, uniform). The `fnl-nvf` build takes longer: about 830 s against 540 s without
  LTO, on 16 cores of a GPU node.
- **Occupancy is not the limit.** Capping the registers at 128 (`r128`) raises the theoretical occupancy from 2 to 4
  blocks of 128 threads per SM and makes the face kernel 16% slower: the spills it adds cost more than the warps it
  gains. Under LTO `maxregcount:128` is not honoured (`lto-r128`: 148 and 134 registers, both above the cap), so that
  variant tests no cap and the conclusion rests on `r128` alone.
- **A compile-time stencil width pays only with LTO.** `patches/face-kernel-s3.patch` (experiment B0) makes the
  half-width `S` a constant, `SC = 3`, in `compute_face_fluxes_dev`: alone it saves 5%, because the arrays are still
  passed to calls; with LTO the calls are inlined, the constant bounds presumably let the stencil loops unroll (not
  checked in the SASS), and the face kernel drops to a third of the leg's. The patch is **valid only for WENO-5 inputs**
  and is not in the code: making it a feature means one kernel per supported `S`, chosen on the host (one `select case`,
  no configuration branch inside a kernel).
- **The compiler's own inlining does not work here.** `-Minline` inlines only routines of the same file (cross-file
  inlining needs extracted libraries and still refuses `weno_reconstruct_upwind_dev`: "missing prototype during
  crossing files"). With the library bodies included into the kernel file and `-Minline=reshape` everything inlines, but
  the inlined routines' local arrays are left out of the kernel's `private` clause and become implicit copies and
  reductions (`Generating implicit reduction(+:..inline)`): the launch faults in `__pgi_uacc_dataup`. `-gpu=lto`
  inlines after the OpenACC lowering, and does not have this problem.
- **On the refined walls the time outside the kernels stays.** The loop minus the kernel time is 7.95 s for the leg
  and 8.04 s for `B0-lto`; as the kernels shrink from 27.7 s to 13.2 s, that share grows from 22% to 38% of the loop.
  Most of it is not attributed yet (see above).

**Correctness of the variants.** Every variant but `B0` (timed only) ran `small`, refined walls, 10 steps with
checkpoints at steps 0 and 10 (`--steps 10 --io`); `h5diff` found every checkpoint bitwise identical to the leg's, and
so does the figure, cell by cell. The regrid at step 10 leaves the grid unchanged (the markers are geometric), so it
rebuilds the maps and moves no data: these runs do not check a regrid data transfer.

![u and T against the exact Couette solution, their error, and the variants against the leg](figures/fields.png)

The initial condition is the exact solution at the cell centres (step 0: errors of round-off, 1e-16). After 10 steps u
and T follow the exact profiles; the error stays at round-off (1e-16 to 1e-14) inside each level and grows only where
the scheme is not exact for this solution: at the isothermal wall (4e-7 on T: its mirror ghost is second order, as in
VV-3) and on both sides of the two 2:1 seams, y = 0.25 and 0.75, where the wall layers meet the base level (up to 1e-8;
the seam ghosts are prolonged and restricted). This case is still a weak check of a build change: the solution is
steady, so 10 steps move it only by these amounts. It shows that the variants run and agree to the last bit here; it is
not a regression suite. A build change adopted from these results passes the FLUME and PRISM FNL regression suites
first.

**Nsight Compute on WSL** (the second GPU, `small` uniform, one launch of the face kernel) agrees on the cause and not
on the cost: the leg's kernel writes 5.75 GB to DRAM for 89 MB of fluxes (the rest is local memory) and executes 139 M
local loads and stores against 33 M global loads; `lto` cuts the local traffic by half and `B0-lto` by 83%, with 65%
fewer instructions. Its duration changes by less than 10%: the RTX 4070 runs FP64 at 1/64 of the FP32 rate, and the
full Nsight Compute set of the leg's kernel names FP64 as its busiest pipeline (83%); the A30 runs FP64 at 1/2 of the
FP32 rate.

### Reproducing the results

The measurements are `results/runs.csv`, one row per run directory (variant, size, Slurm job, date, commit, executable,
wall time, forest loop and steps, the mean per launch of the two flux kernels, the total kernel time), and the figures
are drawn from it and from the checkpoints of the correctness runs:

```bash
./make_figures.py collect runs/            # runs/ here or a copy of the profiling machine's; -> results/runs.csv
./make_figures.py figures runs/            # -> figures/*.png (needs numpy, matplotlib, h5py)
```

On medooza, one A30 of a `dev-gn` node held exclusively, every run with `--launcher "srun --mpi=pmix -n"`:

- job 1357614: the leg built (`../build.sh fnl-nvf --varset iac_nvf`) and timed, `./profile.sh fnl-nvf --size medium
  [--uniform] --np 1 --tool nsys`;
- job 1357619: `B0` built from the patch (applied, built with `fobis build --mode flume-fnl-nvf --varset iac_nvf`, kept
  as `../exe/fnl-nvf-B0/`, reverted), then the leg and `B0` timed twice each, alternated (`--exe`; `--tag` did not
  exist yet);
- job 1357634: the variants built, `../build.sh fnl-nvf --varset iac_nvf --jobs 16 --tag T --gpu G` for `lto`, `r128`
  (`maxregcount:128`), `lto-r128` and, with `patches/face-kernel-s3.patch` applied and then reverted, `B0-lto`;
- job 1357635: for the leg and each of the four, `./profile.sh fnl-nvf --size small --np 1 --steps 10 --io [--tag T]`
  (correctness, compared with `h5diff`) and the two timings above `[--tag T]`.

The CPU before/after runs (job 1357589) are `./profile.sh cpu --size medium --uniform --np 16 --tool none --exe
<exe>`, alternated, three each. The figures average the nsys runs of one rank and 20 steps (`FIG_NP`, `FIG_STEPS` in
`make_figures.py`) and warn when a bar mixes commits; the CPU figure is the runs of job 1357589.
