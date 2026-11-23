#!/bin/bash
# FLUME speedup matrix, FNL OpenACC backend (fnl-nvf) on medooza: one stage of run.sh on one node of 4 A30 (#133).
#
# Submit from src/tests/flume/speedup of the medooza clone (the job finds run.sh through $SLURM_SUBMIT_DIR, or
# $SPEEDUP_DIR), after ./build.sh fnl-nvf --varset iac_nvf (in a batch job). The arguments are run.sh's: the stage, then
# its options; --backend fnl-nvf, --np (one rank per GPU) and --work are set here (a --work given in the arguments
# wins). Every output goes under $SPEEDUP_WORK (default /ifs/hpc/home/szaghi/adam-tests/speedup-133), never into the
# repository. The windows and chain stages end with analyze.py into $SPEEDUP_WORK/results/.
#
#   sbatch jobs/medooza-fnl.slurm warmup --cases 2          # the warm-up of leg 2 on 4 A30 (the fastest backend, D3)
#   sbatch jobs/medooza-fnl.slurm chain --cases 2           # warm-up, settle and windows of leg 2
#   sbatch --gres=gpu:1 --ntasks-per-node=1 jobs/medooza-fnl.slurm windows --cases 3a   # one A30
#
# Medooza rules: batch only, at most 4 jobs at once, the user-space toolchain (never the cluster modules); the prod-gn
# drivers (545) are older than CUDA 12.9 wants, which only HPC-X UCC trips on: UCC off. A node whose GPUs another
# process holds is refused before the run (the FNL block budget is sized on the device total): resubmit with -x <node>.
#SBATCH --job-name=flume-speedup-fnl
#SBATCH --account=iacroma
#SBATCH --partition=prod-gn
#SBATCH --nodes=1
#SBATCH --gres=gpu:4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=/ifs/hpc/home/szaghi/adam-tests/%x-%j.txt
set -euo pipefail
# shellcheck source=/dev/null
source /ifs/hpc/home/szaghi/nvidia/nvhpc.sh
export OMPI_MCA_coll_ucc_enable=0 OMPI_MCA_coll=^ucc
SPEEDUP_DIR="${SPEEDUP_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
SPEEDUP_WORK="${SPEEDUP_WORK:-/ifs/hpc/home/szaghi/adam-tests/speedup-133}"
[[ -x $SPEEDUP_DIR/run.sh ]] || { echo "medooza-fnl.slurm: no run.sh in $SPEEDUP_DIR (submit from speedup/)" >&2 ; exit 2 ; }
case "$SPEEDUP_WORK" in /home/*) echo "medooza-fnl.slurm: $SPEEDUP_WORK is under /home: runs go on /ifs" >&2 ; exit 2 ;; esac
[[ $# -ge 1 ]] || { echo "medooza-fnl.slurm: give run.sh's stage and options" >&2 ; exit 2 ; }
echo "job ${SLURM_JOB_ID:-} on $(hostname), $(date -Is): run.sh $*"
nvidia-smi --query-gpu=index,name,driver_version,memory.total,memory.free --format=csv,noheader
# preflight: every GPU of the job must have at least 20000 MiB free (another process holding a GPU fails the
# allocation at start, or worse shares the device under the timing)
while IFS= read -r free; do
   if [[ ${free// /} -lt 20000 ]]; then
      echo "medooza-fnl.slurm: a GPU of $(hostname) has only $free MiB free: resubmit with -x $(hostname)" >&2
      exit 4
   fi
done < <(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits)
"$SPEEDUP_DIR/run.sh" "$1" --backend fnl-nvf --np "${SLURM_NTASKS:-4}" --work "$SPEEDUP_WORK" \
   --launcher "srun --mpi=pmix -n" "${@:2}"
if [[ $1 == windows || $1 == chain ]]; then
   python3 -I "$SPEEDUP_DIR/analyze.py" --work "$SPEEDUP_WORK" --out "$SPEEDUP_WORK/results"
fi
