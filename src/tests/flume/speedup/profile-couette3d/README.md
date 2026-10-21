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

Rough wall times of the 20 steps, without a profiler (estimates: nothing has been measured on this case yet): `small`
~2 min on 2 WSL CPU ranks, seconds on one GPU; `medium` ~20 s on the 64-core medooza node; `large` ~1 min on one AMD
GPU. Start-up adds seconds to a minute. Profilers add their own overhead: Nsight Systems and rocprofv3 tens of
percent; Nsight Compute replays each profiled kernel many times, so `profile.sh` limits it to 300 kernels of rank 0.

## Run

Build the leg first (`cd .. && ./build.sh <leg> --varset ...`), then from this directory:

```bash
./profile.sh <leg> --size <small|medium|large> [--uniform] --np <N> --tool <none|nsys|ncu|rocprof|perf> [--steps N] [--io]
```

For example, the uniform variant on WSL: `./profile.sh cpu --size small --uniform --np 2 --tool perf`.

Every run goes to `runs/<leg>-<size>[-uniform]-np<N>-<tool>/` (a previous one is renamed, not overwritten): `input.ini`,
`log.txt`, `wall.txt`, `run.info` and the profiler files, one per rank. `ADAM_SEAM_TIMING=1` is always set, so the
end of `log.txt` carries the forest timing report (loop, steps, reflux, seam fills; max over the ranks), which
`profile.sh` prints. `--dry-run` prints the command without running it.

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
./profile.sh fnl-nvf --size medium --np 1 --tool ncu  --launcher "srun --mpi=pmix -n"
./profile.sh fnl-nvf --size large  --np 4 --tool nsys --launcher "srun --mpi=pmix -n"
# CPU node: source /ifs/hpc/home/szaghi/gnu/gnu.sh; export OMPI_MCA_coll_hcoll_enable=0
./profile.sh cpu     --size medium --np 64 --tool perf --launcher "srun --mpi=pmix --cpu-bind=cores -n"
```

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
