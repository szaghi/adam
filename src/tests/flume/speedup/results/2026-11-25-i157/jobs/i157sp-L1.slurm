#!/bin/bash
# #157 night 2: leg 1 windows (1a-step-uniform 1b-step-fixed 1c-step-runtime) on 4 A30, executables A (b8641a59) and B (#157) alternated.
#SBATCH --job-name=i157sp-L1
#SBATCH --account=iacroma
#SBATCH --partition=prod-gn
#SBATCH --nodes=1
#SBATCH --gres=gpu:4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=10:00:00
#SBATCH --output=/ifs/hpc/home/szaghi/adam-tests/speedup-133/i157/%x.%j.out
set -uo pipefail
source /ifs/hpc/home/szaghi/adam-tests/speedup-133/i157/i157sp-common.sh
prepare
windows 4 1a-step-uniform 1b-step-fixed 1c-step-runtime
for t in A B; do python3 -I "$R/analyze.py" --work "$R/$t" --out "$R/$t/results-L1" > /dev/null 2>&1; log "analyze $t rc $?"; done
log "I157SP-L1-DONE"
