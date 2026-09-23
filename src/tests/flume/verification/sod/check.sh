#!/usr/bin/env bash
# FLUME verification V1: Sod shock tube along x, y and z (issue #35, section 11).
#
# Why: CHASE only ever ran x-aligned 1-D Riemann problems, so its transposed characteristic projection and its
# singular y/z eigenvectors went unnoticed. V1 runs the same Sod problem along each axis (the other two directions
# null) and asserts, through sod_oracle.py:
#   1. accuracy: L1(rho) against the exact Riemann solution <= L1_MAX, for every direction;
#   2. direction invariance: the y and z runs equal the x run after the axes permutation, BITWISE;
#   3. 1-D consistency: every transverse copy of the solution, in every block, is bitwise identical.
#
# Baseline provenance (P3): WENO-5 characteristic, SSP-33, CFL 0.5, 200 cells, t = 0.2, np 2 (64 blocks, 32/32):
# CPU L1(rho) = 3.244140e-03, FNL L1(rho) = 3.075502e-03. The two backends differ because the library WENO weights
# use a different exponent on the host (1/(eps+IS)**S) and on the device (1/(eps+IS)**2); with the same exponent
# they agree to 2.4e-14. L1_MAX is the CPU value plus 2%, which bounds both backends.
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
L1_MAX="3.31e-03"
NP=2
BUILD=0

while [[ $# -gt 0 ]]; do
   case "$1" in
      --build) BUILD=1 ; shift ;;
      --np)    NP="$2" ; shift 2 ;;
      *)       echo "check.sh: unknown argument '$1' (accepted: --build, --np N)" >&2 ; exit 2 ;;
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
for axis in x y z; do
   work="$CASE_DIR/work-$TAG-$axis"
   rm -rf "$work"
   mkdir -p "$work"
   cp "$CASE_DIR/sod-$axis.ini" "$work/"
   echo ">> sod-$axis: mpirun -np $NP $(basename "$EXE")"
   if ! (cd "$work" && mpirun -np "$NP" "$EXE" "sod-$axis.ini" > log.txt 2>&1); then
      echo "check.sh: sod-$axis run failed, see $work/log.txt" >&2
      exit 1
   fi
   WORK+=("$work")
done

"$VENV_PY" "$CASE_DIR/sod_oracle.py" "$CASE_DIR/sod-x.ini" "${WORK[@]}" --l1-max "$L1_MAX"
echo "V1 PASSED ($TAG)"
