#!/bin/bash
# FLUME speedup matrix, CPU leg on medooza: whole nodes, one MPI rank per core.
#
# Set the rank count at submission, e.g. one full node:
#   cd ~/adam/src/tests/flume/speedup && sbatch --ntasks=<cores of a node> jobs/medooza-cpu.slurm [run.sh options]
# Build first on the login node: source /ifs/hpc/home/szaghi/gnu/gnu.sh && ./build.sh cpu --varset iac_gnu
# Submit from the speedup directory: run.sh writes data/ and full/ next to the cases.
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
./run.sh cpu --np "$SLURM_NTASKS" --launcher "srun --mpi=pmix --cpu-bind=cores -n" "$@"
