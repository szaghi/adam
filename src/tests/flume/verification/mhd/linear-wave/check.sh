#!/usr/bin/env bash
# FLUME MHD verification MV-5 (issue #41, M2-P5): linear waves, the design order of the smooth MHD path per family.
#
# Why: after one period a small-amplitude eigenmode returns to its initial state, so every smooth-flow ingredient
# (eigenvectors, WENO-5 characteristic reconstruction, fluxes, periodic ghost fill, SSP-54) is measured against an
# exact answer, per wave family (fast, Alfven, slow, entropy), 1-D along x and 2-D inclined (tan a = 2). Through
# linear_wave_oracle.py (eps of Stone et al. 2008, section 8.2, last vs first checkpoint), np 2:
#   - the error of every run below its bound (CPU baseline plus 2 %, M2-P5a);
#   - the observed order of each ladder at least ORDER_MIN (design 5 minus 0.3).
# Amplitude 1e-7, not the 1e-6 of Stone et al.: at 1e-6 the O(A^2) nonlinear error of the mode (measured: eps scales
# as A^2 at N = 128, 3.26e-10 at A = 1e-5, 3.26e-12 at 1e-6) floors a fifth-order scheme already at N = 64-128, and
# below ~1e-8 the accumulated round-off floors it (1.5e-14 at A = 1e-8). Ladders 16/32 (1-D) and 32/64 cells along x
# (2-D, twice the cells per block along x than y), where every family converges at 4.96-5.02; at 2-D N = 128 the slow
# and Alfven waves reach the floor (orders 2.6, 4.3 on 64/128).
# One run at a time; checkpoints deleted after use.
#
# Usage: ./check.sh [--np N]
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/../.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
NP=2
ORDER_MIN="4.7"
AMPLITUDE="1.0e-7"
declare -A EPS_MAX=( # CPU baseline (M2-P5a) plus 2 %, per geometry, wave and N
   [1d-fast-16]=1.512e-10    [1d-fast-32]=4.843e-12
   [1d-alfven-16]=8.693e-11  [1d-alfven-32]=2.784e-12
   [1d-slow-16]=9.621e-11    [1d-slow-32]=3.079e-12
   [1d-entropy-16]=9.220e-11 [1d-entropy-32]=2.951e-12
   [2d-fast-32]=1.309e-10    [2d-fast-64]=4.178e-12
   [2d-alfven-32]=1.227e-10  [2d-alfven-64]=3.940e-12
   [2d-slow-32]=6.826e-11    [2d-slow-64]=2.106e-12
   [2d-entropy-32]=7.502e-11 [2d-entropy-64]=2.386e-12
)

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np) NP="$2" ; shift 2 ;;
      *)    echo "check.sh: unknown argument '$1' (accepted: --np N)" >&2 ; exit 2 ;;
   esac
done
if [[ ! -x "$EXE" ]]; then
   echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2
   exit 2
fi
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
ORACLE="$CASE_DIR/linear_wave_oracle.py"
TAG="$(basename "$EXE")-np$NP"
FAILED=0

case_run() { # case_run <geometry> <wave> <cells>: run one case, print its work directory
   local base w="$CASE_DIR/work-$TAG-$1-$2-$3"
   if [[ $1 == 1d ]]; then base="$VERIF_DIR/sod/sod-x.ini" ; else base="$VERIF_DIR/vortex/vortex-n064.ini" ; fi
   rm -rf "$w" ; mkdir -p "$w"
   "$VENV_PY" "$CASE_DIR/make_linear_wave.py" "$base" "$w/linear-wave.ini" --wave "$2" --geometry "$1" --cells "$3" \
      --amplitude "$AMPLITUDE"
   if ! (cd "$w" && mpirun -np "$NP" "$EXE" linear-wave.ini > log.txt 2>&1); then
      echo "check.sh: linear wave $1 $2 N=$3 run failed, see $w/log.txt" >&2
      return 1
   fi
   echo "$w"
}

for geo in 1d 2d; do
   if [[ $geo == 1d ]]; then ladder=(16 32) ; else ladder=(32 64) ; fi
   for wave in fast alfven slow entropy; do
      echo ">> MV-5 $wave wave, $geo, N = ${ladder[*]} ($(basename "$EXE"), np $NP)"
      works=() ; bounds=()
      for n in "${ladder[@]}"; do
         w="$(case_run $geo $wave $n)" || exit 1
         works+=("$w") ; bounds+=("${EPS_MAX[$geo-$wave-$n]}")
      done
      "$VENV_PY" "$ORACLE" "${works[@]}" --eps-max "${bounds[@]}" --order-min "$ORDER_MIN" || FAILED=1
      for w in "${works[@]}"; do find "$w" -name '*.h5' -delete; done
   done
done

if [[ $FAILED -eq 0 ]]; then
   echo "MV-5 PASSED ($TAG)"
else
   echo "MV-5 FAILED ($TAG)"
   exit 1
fi
