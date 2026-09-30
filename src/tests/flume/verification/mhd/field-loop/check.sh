#!/usr/bin/env bash
# FLUME MHD verification MV-9 (issue #41, M2-P6): field loop advection, divergence error and loop dissipation.
#
# Why: a weak loop (A0 = 1e-3, R = 0.3) advected obliquely with the out-of-plane velocity w = 1 (Gardiner and Stone
# 2005; Mignone and Tzeferacos 2010, section 4.4.1): dB_z/dt = w div(B), so the generated <|B_z|> measures the
# divergence error GLM leaves. [-1, 1] x [-0.5, 0.5], t = 2, GLM, np 2. Through field_loop_oracle.py:
#   - <|B_z|> / A0 below its bound (CPU baseline plus 2 %, M2-P6c) and decreasing under refinement (rate > 0). The
#     loop field jumps at r = R, so the discrete div(B) there is ~A0 / h over a ring ~h wide: the B_z source is O(1)
#     in h and the decrease is slow, no convergence order is expected (measured 1.02e-3, 5.89e-4, 5.13e-4 at
#     N = 64, 128, 256: rates 0.79, 0.20; the Mignone and Tzeferacos GLM level is ~1e-3);
#   - the magnetic energy E_B(T) / E_B(0) above its bound (CPU baseline minus 0.2 %; measured 0.872, 0.937, 0.967 at
#     N = 64, 128, 256).
# Ladder 64/128 (N x N/2 cells), about 10 min on the CPU (256 adds 26 min for no new information).
#
# --amr runs the AMR leg instead (M2-P7b): ladder 32/64 with [0, 1] x [-0.5, 0.5] refined once (2:1 seams at x = 0
# and x = 1 that the loop crosses twice, octree with null z, nk = 4, see ../amr_box.py), t = 1 (one crossing), the
# same checks with their own baselines (measured <|B_z|>/A0 2.23e-3, 1.10e-3, rate +1.03; E_B 0.802, 0.901, against
# 1.99e-3, 1.11e-3 and 0.759, 0.884 of the uniform 32/64 runs: the refined half dissipates less, the seams add B_z at
# N = 32 and none at N = 64), plus the seam div(B) (seam_divb_oracle.py): its final value at most SEAM_DECAY times its
# peak (measured 0.033, 0.011; the peak is the initial-data error of the loop edge at step 1, GLM removes it; the
# collocated PRISM seam runs away instead, issue #29). About 28 min on the CPU (the N = 64 AMR run takes 23).
#
# Usage: ./check.sh [--np N] [--amr]
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/../.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
NP=2
LADDER=(64 128)
BZ_MAX=(1.040e-03 6.003e-04)   # CPU baseline (M2-P6c) plus 2 %, per N of the ladder
BZ_RATE_MIN="0.0"
ENERGY_MIN=(0.87011 0.93544)   # CPU baseline (M2-P6c) minus 0.2 %, per N of the ladder
AMR=0
LADDER_AMR=(32 64)
BZ_MAX_AMR=(2.277e-03 1.118e-03)   # CPU baseline (M2-P7b) plus 2 %
ENERGY_MIN_AMR=(0.80066 0.89888)   # CPU baseline (M2-P7b) minus 0.2 %
SEAM_DECAY="0.1"

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np)  NP="$2" ; shift 2 ;;
      --amr) AMR=1 ; shift ;;
      *)     echo "check.sh: unknown argument '$1' (accepted: --np N, --amr)" >&2 ; exit 2 ;;
   esac
done
if [[ ! -x "$EXE" ]]; then
   echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2
   exit 2
fi
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
TAG="$(basename "$EXE")-np$NP"
MAKE_OPTS=()
LEG=""
if [[ $AMR -eq 1 ]]; then
   TAG="$TAG-amr" ; LEG=", right half refined 2:1"
   LADDER=("${LADDER_AMR[@]}") ; BZ_MAX=("${BZ_MAX_AMR[@]}") ; ENERGY_MIN=("${ENERGY_MIN_AMR[@]}")
   MAKE_OPTS=(--time-max 1.0 --refine-box 0.0 -0.5 1.0 0.5)
fi

echo ">> MV-9 field loop, N = ${LADDER[*]}$LEG ($(basename "$EXE"), np $NP)"
works=()
for n in "${LADDER[@]}"; do
   w="$CASE_DIR/work-$TAG-n$n"
   rm -rf "$w" ; mkdir -p "$w"
   "$VENV_PY" "$CASE_DIR/make_field_loop.py" "$VERIF_DIR/vortex/vortex-n064.ini" "$w/field-loop.ini" --cells "$n" \
      "${MAKE_OPTS[@]}"
   if ! (cd "$w" && mpirun -np "$NP" "$EXE" field-loop.ini > log.txt 2>&1); then
      echo "check.sh: field loop N=$n run failed, see $w/log.txt" >&2
      exit 1
   fi
   works+=("$w")
done
if "$VENV_PY" "$CASE_DIR/field_loop_oracle.py" "${works[@]}" --bz-max "${BZ_MAX[@]}" --bz-order-min "$BZ_RATE_MIN" \
      --energy-min "${ENERGY_MIN[@]}"; then
   STATUS=0
else
   STATUS=1
fi
if [[ $AMR -eq 1 ]]; then
   for w in "${works[@]}"; do
      "$VENV_PY" "$CASE_DIR/../seam_divb_oracle.py" "$w" --decay-max "$SEAM_DECAY" || STATUS=1
   done
fi
for w in "${works[@]}"; do find "$w" -name '*.h5' -delete; done

if [[ $STATUS -eq 0 ]]; then
   echo "MV-9 PASSED ($TAG)"
else
   echo "MV-9 FAILED ($TAG)"
   exit 1
fi
