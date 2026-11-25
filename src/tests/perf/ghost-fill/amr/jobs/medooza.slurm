#!/bin/bash
# The AMR ghost-fill benchmark on one A30 of medooza: build (iac_nvf) and run.
#   sbatch --chdir=<an ADAM tree>/src/tests/perf/ghost-fill/amr jobs/medooza.slurm
#SBATCH -J amr-ghost
#SBATCH -A iacroma
#SBATCH -p prod-gn
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH -t 03:00:00
#SBATCH -o amr-ghost-%j.out
source /ifs/hpc/home/szaghi/nvidia/nvhpc.sh ompi5 > /dev/null 2>&1 || exit 1
export OMPI_MCA_coll_ucc_enable=0 OMPI_MCA_coll=^ucc
./run.sh iac_nvf "results/a30-${SLURM_JOB_ID}.txt"
