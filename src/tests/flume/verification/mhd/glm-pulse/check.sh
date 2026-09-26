#!/usr/bin/env bash
# FLUME MHD verification MV-3 (issue #41, M2-P4): the GLM pulse against the exact solution of the (B_x, psi) pair.
#
# Why: in 1-D the (B_x, psi) pair of mixed GLM obeys the telegraph equations, decoupled from the fluid; a Gaussian pulse
# in B_x on a fluid at rest has an exact solution (d'Alembert without damping), computed mode by mode by
# glm_pulse_oracle.py. It checks the GLM block of the face flux, the damping source (glm_alpha, glm_damping_length,
# min-cell), c_h in dt (c_h = 2 is above the fastest wave, 1.29, so it sets dt; glm_ch_check = error must not trip) and
# the wall parity of psi. Pulse width 0.1 at x = 0.5, t = 0.3 (the pulses cross the periodic boundary), gamma = 5/3,
# WENO-5 characteristic, SSP-33, CFL 0.5, np 2, on the Sod grid of V1 (four blocks along x):
#   1. d'Alembert (glm_alpha = 0), periodic, N = 64, 128, 256: L1 sum of B_x and psi / c_h below L1_MAX_P, observed
#      order at least ORDER_MIN, int B_x constant to CONSERVED_TOL (conservation history);
#   2. damping on (glm_alpha = 0.5, glm_damping_length = 0.1, c_h^2/c_p^2 = 10), N = 64, 128, 256: below L1_MAX_D, order at
#      least ORDER_MIN (a wrong rate is an O(1) relative error: 1.23 with no damping at N = 128, against 9.2e-5);
#   3. damping length min-cell (glm_alpha = 0.1, rate 25.6 at N = 128): below L1_MAX_M;
#   4. reflecting walls (t = 0.4, both pulses reflected once): the odd/even extension, below L1_MAX_W;
#   5. divergence_control = none (N = 128): B_x of the last checkpoint equals the first one BITWISE (zero flux).
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
L1_MAX_P=(1.584e-04 9.342e-06 7.487e-07) # CPU baseline 1.553008e-04, 9.159216e-06, 7.340314e-07 (M2-P4a), plus 2%
L1_MAX_D=(3.393e-05 1.962e-06 1.572e-07) # CPU baseline 3.326672e-05, 1.923056e-06, 1.541121e-07 (M2-P4a), plus 2%
L1_MAX_M="1.487e-07"                     # CPU baseline 1.458195e-07 (M2-P4a), plus 2%
L1_MAX_W="1.364e-05"                     # CPU baseline 1.336939e-05 (M2-P4a), plus 2%
ORDER_MIN="3.0"                          # SSP-33 in time (measured 4.1 and 3.6: WENO-5 in space still visible)
CONSERVED_TOL="1.0e-14"                  # int B_x ~ 1.8e-2 (measured drift <= 5.8e-16)

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
ORACLE="$CASE_DIR/glm_pulse_oracle.py"
TAG="$(basename "$EXE")-np$NP"
FAILED=0

case_run() { # case_run <name> <make_glm_pulse.py options...>: run one case, print its work directory
   local w="$CASE_DIR/work-$TAG-$1"
   shift
   rm -rf "$w" ; mkdir -p "$w"
   "$VENV_PY" "$CASE_DIR/make_glm_pulse.py" "$VERIF_DIR/sod/sod-x.ini" "$w/glm-pulse.ini" "$@"
   if ! (cd "$w" && mpirun -np "$NP" "$EXE" glm-pulse.ini > log.txt 2>&1); then
      echo "check.sh: GLM pulse $(basename "$w") run failed, see $w/log.txt" >&2
      return 1
   fi
   echo "$w"
}

echo ">> MV-3 d'Alembert, periodic, N = 64, 128, 256 ($(basename "$EXE"), np $NP)"
WP=()
for n in 64 128 256; do w="$(case_run p$n --cells $n --divergence-control glm)" || exit 1 ; WP+=("$w"); done
"$VENV_PY" "$ORACLE" "${WP[@]}" --l1-max "${L1_MAX_P[@]}" --order-min "$ORDER_MIN" --conserved "$CONSERVED_TOL" \
   || FAILED=1
echo ">> MV-3 damping on (glm_alpha = 0.5, glm_damping_length = 0.1), N = 64, 128, 256"
WD=()
for n in 64 128 256; do
   w="$(case_run d$n --cells $n --divergence-control glm --alpha 0.5 --damping-length 0.1)" || exit 1 ; WD+=("$w")
done
"$VENV_PY" "$ORACLE" "${WD[@]}" --l1-max "${L1_MAX_D[@]}" --order-min "$ORDER_MIN" || FAILED=1
echo ">> MV-3 damping length min-cell (glm_alpha = 0.1), N = 128"
WM="$(case_run m128 --cells 128 --divergence-control glm --alpha 0.1 --damping-length min-cell)" || exit 1
echo "   $(grep -m1 'MHD GLM damping' "$WM/log.txt" | sed 's/^.*MHD GLM damping: //')"
"$VENV_PY" "$ORACLE" "$WM" --l1-max "$L1_MAX_M" || FAILED=1
echo ">> MV-3 reflecting walls (B_x odd, psi even), t = 0.4, N = 128"
WW="$(case_run w128 --cells 128 --divergence-control glm --bc wall --time 0.4)" || exit 1
"$VENV_PY" "$ORACLE" "$WW" --l1-max "$L1_MAX_W" || FAILED=1
echo ">> MV-3 divergence_control = none, N = 128: B_x frozen"
WN="$(case_run n128 --cells 128 --divergence-control none)" || exit 1
"$VENV_PY" "$ORACLE" "$WN" --unchanged || FAILED=1
for w in "${WP[@]}" "${WD[@]}" "$WM" "$WW" "$WN"; do find "$w" -name '*.h5' -delete; done

if [[ $FAILED -eq 0 ]]; then
   echo "MV-3 PASSED ($TAG)"
else
   echo "MV-3 FAILED ($TAG)"
   exit 1
fi
