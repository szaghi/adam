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
# L1(rho) = 3.244140e-03 on CPU and FNL (the backends agree to 7.2e-14 since the host and device WENO weights share
# one exponent, weno_weights_exponent). L1_MAX is that value plus 2%.
#
# --wall runs the reflecting-wall double Sod instead (sod-wall-{x,y,z}.ini, issue #35 P4 exit): walls at both ends,
# high state in the middle, t = 0.4; each run must be mirror-symmetric about the centre within MIRROR_TOL (normal
# momentum negated), the box must be closed (mass and energy relative drift <= CLOSED_TOL; extrapolation ends would
# leak ~10%), and the three directions must agree bitwise. It exercises the wall-inviscid ghost mirror on all six faces.
#
# Usage: ./check.sh [--build] [--np N] [--wall]
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12);
# --build always builds the CPU default, never the override.
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
L1_MAX="3.31e-03"
MIRROR_TOL="1.0e-12"
CLOSED_TOL="1.0e-11"
NP=2
BUILD=0
CASE="sod"

while [[ $# -gt 0 ]]; do
   case "$1" in
      --build) BUILD=1 ; shift ;;
      --np)    NP="$2" ; shift 2 ;;
      --wall)  CASE="sod-wall" ; shift ;;
      *)       echo "check.sh: unknown argument '$1' (accepted: --build, --np N, --wall)" >&2 ; exit 2 ;;
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

TAG="$CASE-$(basename "$EXE")-np$NP"
WORK=()
for axis in x y z; do
   work="$CASE_DIR/work-$TAG-$axis"
   rm -rf "$work"
   mkdir -p "$work"
   cp "$CASE_DIR/$CASE-$axis.ini" "$work/"
   echo ">> $CASE-$axis: mpirun -np $NP $(basename "$EXE")"
   if ! (cd "$work" && mpirun -np "$NP" "$EXE" "$CASE-$axis.ini" > log.txt 2>&1); then
      echo "check.sh: $CASE-$axis run failed, see $work/log.txt" >&2
      exit 1
   fi
   WORK+=("$work")
done

if [[ "$CASE" == "sod" ]]; then
   "$VENV_PY" "$CASE_DIR/sod_oracle.py" "$CASE_DIR/sod-x.ini" "${WORK[@]}" --l1-max "$L1_MAX"
   echo "V1 PASSED ($TAG)"
else
   "$VENV_PY" "$CASE_DIR/sod_oracle.py" "$CASE_DIR/sod-wall-x.ini" "${WORK[@]}" --mirror "$MIRROR_TOL" --closed "$CLOSED_TOL"
   echo "reflecting-wall Sod PASSED ($TAG)"
fi
