#!/usr/bin/env bash
# FLUME profiling case: run the 3-D Couette channel once, under a profiler or none, on any machine.
#
# Why: one short, steady, reproducible run that goes through the main regions of FLUME and of the ADAM library (WENO
# convective and 4th-order viscous fluxes, SSP-54, walls and periodic faces, MPI ghosts, 2:1 seams and reflux, runtime
# regrids, global reductions; make_profile.py lists them). The same script runs on WSL, on medooza and on the AMD node:
# only the leg, the size, the rank count and the profiler change.
#
# Each rank runs under its own profiler instance and writes its own file (the rank comes from OMPI_COMM_WORLD_RANK,
# PMI_RANK or SLURM_PROCID). ADAM_SEAM_TIMING=1 is always exported: the forest then prints, at the end of the log, the
# wall time of its loop, of the steps and of the reflux, and per realm of the seam fills (max over the ranks); it puts
# a barrier before the reflux reduction to separate the load imbalance, a small perturbation.
#
# Usage: ./profile.sh <leg> [--size small|medium|large] [--np N] [--tool T] [--steps N] [--io]
#                     [--launcher "CMD"] [--exe PATH] [--dry-run]
#   leg         cpu, fnl-nvf or fnl-amd; the executable is ../exe/<leg>/ (../build.sh <leg>)
#   --size      small (default; WSL), medium (medooza), large (one AMD GPU, or 4 A30); see make_profile.py
#   --uniform   the simplest variant: one realm on a uniform grid, no AMR marker, no regrid, no 2:1 seam
#   --np        MPI ranks (default 2, or $SLURM_NTASKS); FNL: one rank per GPU
#   --tool      none (default: wall time and the forest timing report only)
#               nsys     NVIDIA Nsight Systems timeline, every rank (fnl-nvf; CPU ranges too)
#               ncu      NVIDIA Nsight Compute kernel metrics, rank 0 only, the first 300 kernels after 200 skipped
#               rocprof  AMD rocprofv3 kernel trace and per-kernel statistics, every rank (fnl-amd)
#               perf     Linux perf record with call graphs, every rank (cpu; host side of the FNL legs)
#   --steps     time steps (default 20: two runtime regrids, at steps 10 and 20)
#   --io        also write the first and last checkpoint (default: no output at all, a pure compute profile)
#   --launcher  the launcher followed by the rank count (default "mpirun -x ADAM_SEAM_TIMING -np"; Open MPI
#               forwards only the variables -x names; srun forwards the environment: "srun -n")
#   --dry-run   print the command, run nothing
#
# Output: runs/<leg>-<size>[-uniform]-np<N>-<tool>/ with input.ini, log.txt, wall.txt (seconds), run.info and the profiler files.
# WSL: one run at a time, FNL np <= 2 (two GPUs), and the UCX knobs of issue #12 (CLAUDE.md, WSL2 GPU+MPI caveats);
# Nsight Compute there needs the GPU counters enabled in the Windows NVIDIA Control Panel.
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
SPEEDUP_DIR="$(cd "$CASE_DIR/.." && pwd)"
REPO_ROOT="$(cd "$SPEEDUP_DIR/../../../.." && pwd)"
LEG="${1:-}"
[[ $# -gt 0 ]] && shift
SIZE=small
NP="${SLURM_NTASKS:-2}"
TOOL=none
STEPS=20
IO=0
UNIFORM=""
LAUNCHER="mpirun -x ADAM_SEAM_TIMING -np"
EXE=""
DRY=0
while [[ $# -gt 0 ]]; do
   case "$1" in
      --size)     SIZE="$2" ; shift 2 ;;
      --np)       NP="$2" ; shift 2 ;;
      --tool)     TOOL="$2" ; shift 2 ;;
      --steps)    STEPS="$2" ; shift 2 ;;
      --io)       IO=1 ; shift ;;
      --uniform)  UNIFORM="-uniform" ; shift ;;
      --launcher) LAUNCHER="$2" ; shift 2 ;;
      --exe)      EXE="$2" ; shift 2 ;;
      --dry-run)  DRY=1 ; shift ;;
      *) echo "profile.sh: unknown argument '$1' (see the usage at the top of the script)" >&2 ; exit 2 ;;
   esac
done
case "$LEG" in
   cpu)             EXE="${EXE:-$SPEEDUP_DIR/exe/cpu/adam_flume_cpu}" ;;
   fnl-nvf|fnl-amd) EXE="${EXE:-$SPEEDUP_DIR/exe/$LEG/adam_flume_fnl}" ;;
   *) echo "profile.sh: the leg must be cpu, fnl-nvf or fnl-amd, got '$LEG'" >&2 ; exit 2 ;;
esac
case "$SIZE" in small|medium|large) ;; *) echo "profile.sh: --size must be small, medium or large" >&2 ; exit 2 ;; esac

# the profiler command each rank runs in front of the executable; RANK is set per rank by the wrapper below, so the
# single quotes are deliberate: $RANK and "$@" expand in each rank's shell, not here
# shellcheck disable=SC2016
case "$TOOL" in
   none)    PROF='' ;;
   nsys)    PROF='nsys profile --trace=cuda,openacc,mpi,osrt --sample=cpu --force-overwrite=true -o nsys-r$RANK' ;;
   ncu)     PROF='if [ "$RANK" = 0 ]; then set -- ncu --set full --launch-skip 200 --launch-count 300 -f -o ncu-r0 "$@"; fi;' ;;
   rocprof) PROF='rocprofv3 --kernel-trace --memory-copy-trace --stats --output-format csv -d rocprof-r$RANK --' ;;
   perf)    PROF='perf record -g -o perf-r$RANK.data --' ;;
   *) echo "profile.sh: --tool must be none, nsys, ncu, rocprof or perf" >&2 ; exit 2 ;;
esac
if [[ $TOOL == ncu ]]; then
   WRAP="RANK=\${OMPI_COMM_WORLD_RANK:-\${PMI_RANK:-\${SLURM_PROCID:-0}}}; $PROF exec \"\$@\""
elif [[ -n $PROF ]]; then
   WRAP="RANK=\${OMPI_COMM_WORLD_RANK:-\${PMI_RANK:-\${SLURM_PROCID:-0}}}; exec $PROF \"\$@\""
else
   WRAP='exec "$@"'
fi

WORK="$CASE_DIR/runs/$LEG-$SIZE$UNIFORM-np$NP-$TOOL"
CMD="$LAUNCHER $NP bash -c '$WRAP' rank $EXE input.ini"
if [[ $DRY -eq 1 ]]; then
   echo "cd $WORK && ADAM_SEAM_TIMING=1 $CMD"
   exit 0
fi
[[ -x $EXE ]] || { echo "profile.sh: executable '$EXE' not found (../build.sh $LEG, or --exe)" >&2 ; exit 2 ; }
if [[ -d $WORK ]]; then mv "$WORK" "$WORK.old-$(date +%Y%m%d-%H%M%S)"; fi
mkdir -p "$WORK"
python3 "$CASE_DIR/make_profile.py" "$SIZE" "$WORK/input.ini" --steps "$STEPS" ${UNIFORM:+--uniform}
if [[ $IO -eq 1 ]]; then
   python3 "$SPEEDUP_DIR/make_inputs.py" set "$WORK/input.ini" "IO.it_save=$STEPS"
fi
{
   echo "date=$(date -Is)"
   echo "host=$(hostname)"
   echo "slurm_job=${SLURM_JOB_ID:-}"
   echo "leg=$LEG"
   echo "size=$SIZE${UNIFORM}"
   echo "np=$NP"
   echo "tool=$TOOL"
   echo "steps=$STEPS"
   echo "io=$IO"
   echo "command=$CMD"
   echo "exe_sha256=$(sha256sum "$EXE" | cut -c1-16)"
   echo "commit=$(git -C "$REPO_ROOT" rev-parse HEAD 2> /dev/null || echo unknown)"
} > "$WORK/run.info"

echo ">> profile: $LEG, size $SIZE${UNIFORM}, np $NP, tool $TOOL, $STEPS steps -> $WORK"
export ADAM_SEAM_TIMING=1
t0=$(date +%s.%N)
status=0
(cd "$WORK" && eval "$CMD" > log.txt 2>&1) || status=$?
t1=$(date +%s.%N)
awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.3f\n", b - a}' > "$WORK/wall.txt"
echo "   wall $(cat "$WORK/wall.txt") s, exit $status"
grep -a "forest timing\|forest seam timing\|regrid at step" "$WORK/log.txt" | sed 's/^/   /' || true
if [[ $status -ne 0 ]]; then
   echo "profile.sh: the run exited $status, see $WORK/log.txt" >&2
   exit 1
fi
