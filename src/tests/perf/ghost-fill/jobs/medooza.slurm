#!/bin/bash
# ghost_fill_bench on one A30 of medooza (prod-gn).   sbatch --chdir=<bench dir> jobs/medooza.slurm
#SBATCH -J ghost-fill
#SBATCH -A iacroma
#SBATCH -p prod-gn
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH -t 01:00:00
#SBATCH -o ghost-fill-%j.out
source /ifs/hpc/home/szaghi/nvidia/nvhpc.sh ompi5 > /dev/null 2>&1 || exit 1
./run.sh cc80,cuda12.9 "results/a30-${SLURM_JOB_ID}.txt"
