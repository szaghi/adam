#!/bin/bash
# #157 night 2: device time per kernel (NVCOMPILER_ACC_TIME=1, 10-step runs from the window dumps, not timings), the
# runs of the #133 "Where the AMR overhead goes" table (2c on 4 A30, w2-sphere-l3 on 2, 2a on 4) and the AMR cases of
# legs 1 and 3 (1b, 3c on 4), executables A (b8641a59) and B (#157).
#SBATCH --job-name=i157sp-P
#SBATCH --account=iacroma
#SBATCH --partition=prod-gn
#SBATCH --nodes=1
#SBATCH --gres=gpu:4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=/ifs/hpc/home/szaghi/adam-tests/speedup-133/i157/%x.%j.out
set -uo pipefail
source /ifs/hpc/home/szaghi/adam-tests/speedup-133/i157/i157sp-common.sh
R=$R/profile
prepare
for tag in A B; do
   win $tag 4 2c-sphere-fixed 1 10 NVCOMPILER_ACC_TIME=1
   win $tag 2 w2-sphere-l3 1 10 NVCOMPILER_ACC_TIME=1
   win $tag 4 2a-sphere-uniform 1 10 NVCOMPILER_ACC_TIME=1
   win $tag 4 1b-step-fixed 1 10 NVCOMPILER_ACC_TIME=1
   win $tag 4 3c-channel-fixed 1 10 NVCOMPILER_ACC_TIME=1
done
log "I157SP-P-DONE"
