#!/bin/bash
# #157 C4: leg 2 windows (2a-sphere-uniform 2b-sphere-runtime 2c-sphere-fixed) on 4 A30, executables A (HEAD 80582712) and B (C4) alternated.
#SBATCH --job-name=i157c-L2
#SBATCH --account=iacroma
#SBATCH --partition=prod-gn
#SBATCH --nodes=1
#SBATCH --gres=gpu:4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/ifs/hpc/home/szaghi/adam-tests/speedup-133/i157c/%x.%j.out
set -uo pipefail
source /ifs/hpc/home/szaghi/adam-tests/speedup-133/i157c/i157c-common.sh
prepare
windows 4 2a-sphere-uniform 2b-sphere-runtime 2c-sphere-fixed
for t in A B; do python3 -I "$R/analyze.py" --work "$R/$t" --out "$R/$t/results-L2" > /dev/null 2>&1; log "analyze $t rc $?"; done
log "I157C-L2-DONE"
