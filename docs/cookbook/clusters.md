# MPI and clusters

Recipes for running on more than one rank, on a workstation and under a batch scheduler. The GPU-specific settings are
in [GPU](/cookbook/gpu).

## C1. Run on several ranks

**Problem.** Run a case with MPI.

**Command** (the FLUME Sod tube, `src/tests/flume/regression/sod-x`):

```bash
cd src/tests/flume/regression/sod-x
mpirun -np 4 ../../../../../exe/adam_flume_cpu input.ini
```

The executable takes the input file as its only argument (default `input.ini`); a file with a `[forest]` section is a
multi-realm manifest ([forest R0](/cookbook/forest#r0-the-driver-pattern)). Each rank writes its own
`-procNNNNNN.h5` field files and `-procNNNNNN.fbd` restart files.

**What you get.** The blocks are spread over the ranks by Morton order, in whole sibling families. A grid with fewer
families than ranks leaves ranks with **no block**: that is supported (they take part in every collective and do no
work), not an error, on both backends ([ranks without blocks](/applications/flume/verification/infrastructure/ranks-without-blocks)).
A restart needs the rank count it was saved with; more or fewer ranks stop the run ([IO4](/cookbook/io#io4-restart-a-run)).

**Reference.** [`[amr] iu_ref_levels`](/reference/library-input#amr) (the family distribution).

## C2. A batch job

**Problem.** Run under Slurm on a cluster node.

**Job** (the shape of `src/tests/flume/speedup/jobs/medooza-cpu.slurm` and `medooza-fnl.slurm`, with the site names
removed):

```bash
#!/bin/bash
#SBATCH --job-name=flume-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=96G
#SBATCH --time=02:00:00
set -euo pipefail
# load the toolchain the executable was built with (the varset's compiler and MPI)
srun --mpi=pmix --cpu-bind=cores -n "$SLURM_NTASKS" ./adam_flume_cpu input.ini
```

For FNL, one rank per GPU:

```bash
#SBATCH --gres=gpu:4
#SBATCH --ntasks-per-node=4
srun --mpi=pmix -n "$SLURM_NTASKS" ./adam_flume_fnl input.ini
```

**Rules.**
- **Same toolchain at run time as at build time.** Load the compiler and MPI of the varset the executable was built
  with; HDF5 and the MPI Fortran modules are compiler-specific.
- **Give `--mem`.** On the CPU backend a rank's block budget is the node's memory, capped by the job's cgroup limit
  (what `--mem` sets), divided by the ranks on the node (`mpih_object`, `src/lib/common/adam_mpih_object.F90`); on FNL
  it is the device memory.
- **GPU binding.** FUNDAL picks device `mod(local_rank, visible devices)` on each node. With one rank per GPU and
  `srun --gpu-bind` (or one visible device per task) every rank gets its own GPU; ranks that see the same device share
  it, and the block budget is divided among them (`count_device_sharers`), so they never oversubscribe it.
- **Open MPI 5 with hcoll on some clusters** segfaults in `MPI_Finalize` (exit 139) after a correct run. Disable it
  with `export OMPI_MCA_coll_hcoll_enable=0`; the regression harnesses fail a case whose run exits non-zero, so a crash
  at finalize does not pass unseen.

**Reference.** [GPU5](/cookbook/gpu#gpu5-run-on-a-development-box-with-wsl2-not-on-a-cluster) for the settings that
belong to a WSL2 workstation only.

## C3. Workstation (WSL2) and cluster settings are not interchangeable

**Problem.** Know which environment settings belong where.

| Setting | WSL2 workstation | Cluster with InfiniBand and GPUDirect |
|---|---|---|
| `UCX_RNDV_THRESH=inf` | required for FNL with 2+ ranks ([#12](https://github.com/szaghi/adam/issues/12)) | **harmful**: rendezvous is the fast GPUDirect path |
| `UCX_TLS=^cma` | required (CMA restricted) | not needed |
| `OMPI_MCA_coll_hcoll_enable=0` | required (no Mellanox hardware) | only where hcoll crashes at finalize (C2) |
| GPU-aware UCX transports, e.g. `UCX_TLS=rc_x,sm,cuda_copy,cuda_ipc` | — | the site's recommended tuning |

Timings measured on WSL2 are not benchmarks: its GPU and MPI stack is a correctness environment only.

**Reference.** [GPU5](/cookbook/gpu#gpu5-run-on-a-development-box-with-wsl2-not-on-a-cluster).
