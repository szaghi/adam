#!/usr/bin/env bash
# FLUME MHD verification MV-6 (issue #41, M2-P5): the circularly polarised Alfven wave, design order at finite amplitude.
#
# Why: the wave of Toth (2000) is an exact nonlinear solution at any amplitude, so after one period the exact solution
# is the initial state; at B_perp = 0.1 it exercises the WENO nonlinear weights that the 1e-7 linear waves of MV-5 leave
# linear. 2-D inclined (tan a = 2, Mignone and Tzeferacos 2010), both polarisations, np 2. Through cpaw_oracle.py
# (eps of Stone et al., last vs first checkpoint):
#   - the error of every run below its bound (CPU baseline plus 2 %, M2-P5b);
#   - the observed order of each ladder at least ORDER_MIN (design 5 minus 0.3);
#   - the left and right errors (mirror images under z -> -z) equal to a relative LR_TOL at every N: 0, bitwise, since
#     rounding is symmetric under a sign change (measured on CPU, M2-P5b).
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
LR_TOL="0.0"
LADDER=(32 64 128)
EPS_MAX=(1.049e-03 3.537e-05 1.112e-06) # CPU baseline (M2-P5b) plus 2 %, per N of the ladder

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
TAG="$(basename "$EXE")-np$NP"

case_run() { # case_run <polarisation> <cells>: run one case, print its work directory
   local w="$CASE_DIR/work-$TAG-$1-$2"
   rm -rf "$w" ; mkdir -p "$w"
   "$VENV_PY" "$CASE_DIR/make_cpaw.py" "$VERIF_DIR/vortex/vortex-n064.ini" "$w/cpaw.ini" --polarisation "$1" \
      --geometry 2d --cells "$2"
   if ! (cd "$w" && mpirun -np "$NP" "$EXE" cpaw.ini > log.txt 2>&1); then
      echo "check.sh: CPAW $1 N=$2 run failed, see $w/log.txt" >&2
      return 1
   fi
   echo "$w"
}

declare -A works=()
for pol in right left; do
   echo ">> MV-6 $pol polarisation, 2d, N = ${LADDER[*]} ($(basename "$EXE"), np $NP)"
   for n in "${LADDER[@]}"; do
      w="$(case_run $pol $n)" || exit 1
      works[$pol]+="$w "
   done
done
# shellcheck disable=SC2086
if "$VENV_PY" "$CASE_DIR/cpaw_oracle.py" --right ${works[right]} --left ${works[left]} --eps-max "${EPS_MAX[@]}" \
      --order-min "$ORDER_MIN" --lr-tol "$LR_TOL"; then
   STATUS=0
else
   STATUS=1
fi
for w in ${works[right]} ${works[left]}; do find "$w" -name '*.h5' -delete; done

if [[ $STATUS -eq 0 ]]; then
   echo "MV-6 PASSED ($TAG)"
else
   echo "MV-6 FAILED ($TAG)"
   exit 1
fi
