#!/bin/bash
# Night 2026-11-27-i161-cpu-channel, #161 G1 and 2.1: the 1a windows (i160 dumps, 64 ranks) A/B of 0f75ftz / b864 (G1:
# the guilty commit against its parent) / fixA / 82c7, 2 pairs; then the 2a windows (the #133 matrix's dumps) A/B of
# 82c7 / fixA, 2 pairs. fixA is built by the i161-fix job (waits for trees/fixA/READY).
#SBATCH -J i161-1a2a
#SBATCH -A iacroma
#SBATCH -p prod-gn
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=256G
#SBATCH -t 04:00:00
#SBATCH -o /ifs/hpc/home/szaghi/adam-tests/speedup-133/i161/%x.%j.out
set -uo pipefail
D=/ifs/hpc/home/szaghi/adam-tests/speedup-133/i161
n=0; until [[ -f $D/trees/fixA/READY ]] || (( n > 120 )); do sleep 60; n=$((n+1)); done
[[ -f $D/trees/fixA/READY ]] || { echo "fixA not built"; exit 4; }
LABELS="0f75ftz b864 fixA 82c7" PAIRS=2 CASE=1a-step-uniform TAGNAME=1a bash $D/i161-ab2.slurm
SRCROOT=$D/../work/windows/cpu-np64 LABELS="82c7 fixA" PAIRS=2 CASE=2a-sphere-uniform TAGNAME=2a bash $D/i161-ab2.slurm
echo "I161-1A2A-DONE $(date -Is)"
