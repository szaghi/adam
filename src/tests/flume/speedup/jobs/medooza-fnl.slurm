#!/bin/bash
# FLUME speedup matrix, FNL OpenACC leg on medooza: one node, one MPI rank per A30 GPU.
#
#   cd ~/adam/src/tests/flume/speedup && sbatch jobs/medooza-fnl.slurm [run.sh options]
#   one GPU instead of four: sbatch --gres=gpu:1 --ntasks-per-node=1 jobs/medooza-fnl.slurm
# Build first on the login node: source /ifs/hpc/home/szaghi/nvidia/nvhpc.sh && ./build.sh fnl-nvf --varset iac_nvf
# dev-gn until the prod-gn drivers are updated (gnode20 runs driver 545, which the NVHPC 26.9 HPC-X cannot use); its
# limit is 1 h, so split the matrix with --cases when the grid is large.
#SBATCH --job-name=flume-speedup-fnl
#SBATCH --account=iacroma
#SBATCH --partition=dev-gn
#SBATCH --nodes=1
#SBATCH --gres=gpu:4
#SBATCH --ntasks-per-node=4
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --output=jobs/%x-%j.txt
set -euo pipefail
source /ifs/hpc/home/szaghi/nvidia/nvhpc.sh
nvidia-smi --query-gpu=name,driver_version --format=csv,noheader
./run.sh fnl-nvf --np "$SLURM_NTASKS" --launcher "srun --mpi=pmix -n" "$@"
