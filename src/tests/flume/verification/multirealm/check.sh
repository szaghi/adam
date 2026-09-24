#!/usr/bin/env bash
# FLUME verification of the multi-realm coupling (issue #37): 2-realm Sod against the single-realm Sod.
#
# Why: a forest manifest glues realms through inter-realm seams. sod-2realm.ini splits the sod-x domain at x = 0.5
# (the diaphragm, so the seam carries the Riemann problem from the first step) into two realms with the sod-x cell
# size; with a mirror seam filled at every Runge-Kutta stage (beta cadence) the seam is a block interface like any
# other, so the union of the two realms must reproduce the single-realm sod-x bit for bit. multirealm_oracle.py
# compares the interior cells by their centres (tolerance 0) and the sum of the realms' conservation histories with
# the single-realm one (round-off: the realms sum their cells separately).
#
# Measured (issue #37, np 2, t = 0.2, 174 steps): CPU and FNL bitwise on all 51200 cells, the summed conservation
# histories within 2e-13. The first runs found three defects: the CPU BC routine stopped on the forest's BC_SEAM crown
# rows; the fine side of the inter-realm reflux register was 2:1-restricted like an AMR seam, which wrote F_coarse - 0
# into three quarters of the realm-1 seam skin at every step; the FNL backend never copied the forest-built seam and
# BC maps to the device (illegal address in the seam fill kernel).
#
# Usage: ./check.sh [--build] [--np N]
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12);
# --build always builds the CPU default, never the override.
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
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

run() { # run <work-dir> <input> <files...>
   local work="$1" input="$2"
   shift 2
   rm -rf "$work" ; mkdir -p "$work" ; cp "$@" "$work/"
   echo ">> $(basename "$work"): mpirun -np $NP $(basename "$EXE") $input"
   if ! (cd "$work" && mpirun -np "$NP" "$EXE" "$input" > log.txt 2>&1); then
      echo "check.sh: run failed, see $work/log.txt" >&2
      exit 1
   fi
}

echo "== leg 1: 2-realm Sod vs sod-x"
single="$CASE_DIR/work-$TAG-single"
multi="$CASE_DIR/work-$TAG-2realm"
run "$single" sod-x.ini "$VERIF_DIR/sod/sod-x.ini"
run "$multi" sod-2realm.ini "$CASE_DIR/sod-2realm.ini" "$CASE_DIR/sod-2realm-r1.ini" "$CASE_DIR/sod-2realm-r2.ini"
"$VENV_PY" "$CASE_DIR/multirealm_oracle.py" "$multi" "$single" --ngc 3 --tol 0
echo "multi-realm verification PASSED ($TAG)"
