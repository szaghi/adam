#!/usr/bin/env bash
# FLUME MHD verification MV-7 (issue #41, M2-P5): the magnetised vortex of Balsara (2004), design order in L1 and Linf.
#
# Why: a steady 2-D equilibrium of every variable (pressure, centrifugal force, magnetic tension and pressure) convected
# by the free stream, with a divergence-free field: after one crossing of the doubly periodic box the exact solution is
# the initial state, and mhd_vortex_oracle.py compares every cell with the exact translated vortex. [-7, 7]^2 (tails
# e^-24 at the box edge), t = 14, GLM (with `none` the discrete div(B) grows exponentially, see make_mhd_vortex.py),
# CFL 0.4 (errors equal to 4-5 digits at CFL 0.1), np 2. Checks:
#   - the L1 and Linf errors of every run below their bounds (CPU baseline plus 2 %, M2-P5c);
#   - the observed L1 order of the finest pair at least ORDER_MIN (design 5 minus 0.3); the Linf order at least
#     LINF_ORDER_MIN (measured 3.41, 3.91, 4.41 on 32/64, 64/128, 128/256: still pre-asymptotic at 256, and 512 costs
#     2.5 h on the CPU);
#   - max|div B| at truncation level: converging at least at DIVB_ORDER_MIN (second-order operator, fdv_order = 2);
#   - the initial checkpoint equal to the exact vortex to IC_TOL.
# 64/128/256 (32 is pre-asymptotic, order 3.65); about 40 min on the CPU, 31 of them at 256. Checkpoints deleted after
# use.
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
LINF_ORDER_MIN="4.0"
DIVB_ORDER_MIN="1.7"
IC_TOL="1.0e-14"
LADDER=(64 128 256)
L1_MAX=(4.327e-04 2.224e-05 6.900e-07)   # CPU baseline (M2-P5c) plus 2 %, per N of the ladder
LINF_MAX=(7.812e-03 5.194e-04 2.439e-05) # idem

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

echo ">> MV-7 magnetised vortex, N = ${LADDER[*]} ($(basename "$EXE"), np $NP)"
works=()
for n in "${LADDER[@]}"; do
   w="$CASE_DIR/work-$TAG-n$n"
   rm -rf "$w" ; mkdir -p "$w"
   "$VENV_PY" "$CASE_DIR/make_mhd_vortex.py" "$VERIF_DIR/vortex/vortex-n064.ini" "$w/mhd-vortex.ini" --cells "$n" \
      --half-width 7.0 --cfl 0.4
   if ! (cd "$w" && mpirun -np "$NP" "$EXE" mhd-vortex.ini > log.txt 2>&1); then
      echo "check.sh: magnetised vortex N=$n run failed, see $w/log.txt" >&2
      exit 1
   fi
   works+=("$w")
done
if "$VENV_PY" "$CASE_DIR/mhd_vortex_oracle.py" "${works[@]}" --l1-max "${L1_MAX[@]}" --linf-max "${LINF_MAX[@]}" \
      --order-min "$ORDER_MIN" --linf-order-min "$LINF_ORDER_MIN" --divb-order-min "$DIVB_ORDER_MIN" \
      --ic-tol "$IC_TOL"; then
   STATUS=0
else
   STATUS=1
fi
for w in "${works[@]}"; do find "$w" -name '*.h5' -delete; done

if [[ $STATUS -eq 0 ]]; then
   echo "MV-7 PASSED ($TAG)"
else
   echo "MV-7 FAILED ($TAG)"
   exit 1
fi
