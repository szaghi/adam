#!/bin/bash
# FLUME speedup matrix, CPU backend (cpu) on medooza: a PLACEHOLDER, it stops before running anything.
#
# The stale matrix this job ran (16 shock-over-a-body cases, run.sh) was removed in issue #133 phase 0. The job of
# the new matrix (warm-up, settle, repartition, windows) is written in phase 3 (step 3.5). Until then the file keeps
# its name and the shape of a medooza job (Slurm header, toolchain, environment), which docs/cookbook/clusters.md
# cites. See ../README.md and #133.
#SBATCH --job-name=flume-speedup-cpu
#SBATCH --account=iacroma
#SBATCH --partition=prod-gn
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=jobs/%x-%j.txt
set -euo pipefail
source /ifs/hpc/home/szaghi/gnu/gnu.sh
export OMPI_MCA_coll_hcoll_enable=0 # hcoll segfaults in MPI_Finalize on the compute nodes (exit 139)
echo "jobs/medooza-cpu.slurm: the speedup matrix harness does not exist yet (issue #133 phase 3)" >&2
exit 2
