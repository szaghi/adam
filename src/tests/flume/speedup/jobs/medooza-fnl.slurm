#!/bin/bash
# FLUME speedup matrix, FNL OpenACC backend (fnl-nvf) on medooza: a PLACEHOLDER, it stops before running anything.
#
# The stale matrix this job ran (16 shock-over-a-body cases, run.sh) was removed in issue #133 phase 0. The job of
# the new matrix (warm-up, settle, repartition, windows) is written in phase 3 (step 3.5). Until then the file keeps
# its name and the shape of a medooza job (Slurm header, toolchain, environment), which docs/cookbook/clusters.md
# cites. See ../README.md and #133.
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
echo "jobs/medooza-fnl.slurm: the speedup matrix harness does not exist yet (issue #133 phase 3)" >&2
exit 2
