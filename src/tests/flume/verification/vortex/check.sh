#!/usr/bin/env bash
# FLUME verification V2: isentropic vortex convergence (issue #35, section 11).
#
# Why: V1 (Sod) proves correctness across discontinuities, where WENO is first order; it cannot see a loss of design
# order on smooth flows (a wrong stencil weight, a wrong ghost fill across periodic block faces, a first-order
# boundary). The isentropic vortex is an exact convected solution in a doubly periodic box: every run is compared
# pointwise with it and the observed L1 order of the finest pair must reach ORDER_MIN.
#
# Also exercised: library periodicity across blocks and ranks (16 quadtree blocks, 8/8 over two ranks) and the
# `isentropic-vortex` initial condition.
#
# The ladder is 64/128/256 and the order is asserted on the finest pair only: 64 -> 128 is pre-asymptotic (the vortex
# radius spans 4.5 cells at 64, and WENO-JS weights are still far from optimal there; measured L1 order 4.15).
#
# Usage: ./check.sh [--build] [--np N]
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12);
# --build always builds the CPU default, never the override.
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
ORDER_MIN="4.5"
NP=2
BUILD=0
LADDER=(064 128 256)

while [[ $# -gt 0 ]]; do
   case "$1" in
      --build)       BUILD=1 ; shift ;;
      --np)          NP="$2" ; shift 2 ;;
      *)             echo "check.sh: unknown argument '$1' (accepted: --build, --np N)" >&2 ; exit 2 ;;
   esac
done

if [[ $BUILD -eq 1 ]]; then
   (cd "$REPO_ROOT" && fobis build --mode flume-cpu-gnu)
fi
if [[ ! -x "$EXE" ]]; then
   echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2
   exit 2
fi

VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
if ! "$VENV_PY" -c 'import h5py, numpy' 2>/dev/null; then
   echo ">> creating the oracle venv at $REPO_ROOT/exe/.regression-venv"
   python3 -m venv "$REPO_ROOT/exe/.regression-venv"
   "$VENV_PY" -m pip install --quiet --upgrade pip
   "$VENV_PY" -m pip install --quiet h5py numpy
fi

TAG="$(basename "$EXE")-np$NP"
WORK=()
for n in "${LADDER[@]}"; do
   work="$CASE_DIR/work-$TAG-n$n"
   rm -rf "$work"
   mkdir -p "$work"
   cp "$CASE_DIR/vortex-n$n.ini" "$work/"
   echo ">> vortex-n$n: mpirun -np $NP $(basename "$EXE")"
   start=$(date +%s)
   if ! (cd "$work" && mpirun -np "$NP" "$EXE" "vortex-n$n.ini" > log.txt 2>&1); then
      echo "check.sh: vortex-n$n run failed, see $work/log.txt" >&2
      exit 1
   fi
   echo "   done in $(( $(date +%s) - start )) s"
   WORK+=("$work")
done

"$VENV_PY" "$CASE_DIR/vortex_oracle.py" "${WORK[@]}" --order-min "$ORDER_MIN"
echo "V2 PASSED ($TAG, N = ${LADDER[*]})"
