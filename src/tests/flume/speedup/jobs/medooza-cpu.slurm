#!/bin/bash
# FLUME speedup matrix, CPU backend (cpu) on medooza: one stage of run.sh on one node (issue #133).
#
# Submit from src/tests/flume/speedup of the medooza clone (the job finds run.sh through $SLURM_SUBMIT_DIR, or
# $SPEEDUP_DIR), after ./build.sh cpu --varset iac_gnu (in a batch job: nothing but sbatch on the login node). The
# arguments are run.sh's: the stage, then its options; --backend cpu and --work are set here (a --work given in the
# arguments wins). Every output goes under $SPEEDUP_WORK (default /ifs/hpc/home/szaghi/adam-tests/speedup-133), never
# into the repository. The windows and chain stages end with analyze.py into $SPEEDUP_WORK/results/.
#
#   sbatch jobs/medooza-cpu.slurm chain --cases 2           # warm-up, settle and windows of leg 2 on 64 ranks
#   sbatch jobs/medooza-cpu.slurm windows --cases 1a,1b,1c  # windows only, from the dumps already on /ifs
#   sbatch --ntasks-per-node=32 jobs/medooza-cpu.slurm windows --cases 2 --np 32
#
# The node: 64 ranks (one per core), --mem sized for the largest case (the block pool is 0.85 x --mem / ranks; 2a at
# L4 holds ~48 GB of fields). Medooza rules: batch only, at most 4 jobs at once, the user-space toolchain (never the
# cluster modules), hcoll off (it segfaults in MPI_Finalize on the compute nodes, exit 139).
#SBATCH --job-name=flume-speedup-cpu
#SBATCH --account=iacroma
#SBATCH --partition=prod-gn
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --output=/ifs/hpc/home/szaghi/adam-tests/%x-%j.txt
set -euo pipefail
# shellcheck source=/dev/null
source /ifs/hpc/home/szaghi/gnu/gnu.sh
export OMPI_MCA_coll_hcoll_enable=0
SPEEDUP_DIR="${SPEEDUP_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
SPEEDUP_WORK="${SPEEDUP_WORK:-/ifs/hpc/home/szaghi/adam-tests/speedup-133}"
[[ -x $SPEEDUP_DIR/run.sh ]] || { echo "medooza-cpu.slurm: no run.sh in $SPEEDUP_DIR (submit from speedup/)" >&2 ; exit 2 ; }
case "$SPEEDUP_WORK" in /home/*) echo "medooza-cpu.slurm: $SPEEDUP_WORK is under /home: runs go on /ifs" >&2 ; exit 2 ;; esac
[[ $# -ge 1 ]] || { echo "medooza-cpu.slurm: give run.sh's stage and options" >&2 ; exit 2 ; }
echo "job ${SLURM_JOB_ID:-} on $(hostname), $(date -Is): run.sh $*"
lscpu | grep -E '^(Model name|Socket|Core|Thread)' || true
"$SPEEDUP_DIR/run.sh" "$1" --backend cpu --np "${SLURM_NTASKS:-64}" --work "$SPEEDUP_WORK" --launcher "srun --mpi=pmix -n" \
   "${@:2}"
if [[ $1 == windows || $1 == chain ]]; then
   python3 -I "$SPEEDUP_DIR/analyze.py" --work "$SPEEDUP_WORK" --out "$SPEEDUP_WORK/results"
fi
