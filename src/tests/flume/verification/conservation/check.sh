#!/usr/bin/env bash
# FLUME verification V3: conservation across AMR coarse-fine faces with reflux (issue #35, section 11).
#
# Why: at a 2:1 coarse-fine face the coarse and fine sides compute different fluxes, so the grid leaks mass, momentum
# and energy unless the Berger-Colella reflux corrects it. The case is a triply periodic box (nothing crosses the
# boundary) with one refined octant: 6 coarse-fine faces, periodic ones included, carrying a uniform flow with a 1%
# seeded perturbation. Two legs:
#   1. reflux on : the five volume integrals must be constant within MAX_DRIFT (round-off);
#   2. reflux off: the negative control must drift by at least MIN_DRIFT, which proves the seams are exercised.
#
# FLUME accumulates the seam fluxes of every Runge-Kutta stage weighted by its SSP coefficient, so the register holds
# the flux the committed step actually used; measured (P5): drift <= 6e-15 CPU, 2.2e-16 FNL, with reflux; 2.2e-5
# without, on np 1 and np 2.
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
MAX_DRIFT="1.0e-13"
MIN_DRIFT="1.0e-10"
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
for reflux in true false; do
   work="$CASE_DIR/work-$TAG-reflux-$reflux"
   rm -rf "$work"
   mkdir -p "$work"
   sed "s/^reflux                   = .true./reflux                   = .$reflux./" \
       "$CASE_DIR/amr-periodic.ini" > "$work/amr-periodic.ini"
   echo ">> amr-periodic, reflux .$reflux.: mpirun -np $NP $(basename "$EXE")"
   if ! (cd "$work" && mpirun -np "$NP" "$EXE" amr-periodic.ini > log.txt 2>&1); then
      echo "check.sh: amr-periodic (reflux .$reflux.) run failed, see $work/log.txt" >&2
      exit 1
   fi
   if ! grep -q "registered intra-realm AMR seam faces: +6" "$work/log.txt"; then
      echo "check.sh: the case must register 6 coarse-fine faces, see $work/log.txt" >&2
      exit 1
   fi
done

"$VENV_PY" "$CASE_DIR/conservation_oracle.py" --conserved "$CASE_DIR/work-$TAG-reflux-true" --max-drift "$MAX_DRIFT" \
                                              --leaky "$CASE_DIR/work-$TAG-reflux-false" --min-drift "$MIN_DRIFT"
echo "V3 PASSED ($TAG)"
