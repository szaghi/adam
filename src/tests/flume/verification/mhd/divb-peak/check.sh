#!/usr/bin/env bash
# FLUME MHD verification MV-10 (issue #41, M2-P4): the Dedner et al. (2002) div(B) peak, the div(B) history and monitor.
#
# Why: the peak in B_x has a non-zero div(B) from t = 0; mixed GLM must remove it, while without cleaning it is only
# advected. It checks the div(B) history (max, L1, seam-local max; host-side FD of the library, null directions
# weighted zero), the divb_tol monitor (warning and stop), the derived output fields (pt, beta, bmag, divb) and, with
# a 2:1 AMR box, the seam-local column. Through divb_peak_oracle.py, gamma = 5/3, WENO-5 characteristic, SSP-33,
# CFL 0.5, np 2, 2-D periodic [-0.5, 1.5]^2, t = 0.5 (c_h = 5, glm_alpha = 0.4, damping length min-cell):
#   1. uniform 64^2, GLM (auxiliary fields saved) and none: the GLM L1 of div B never exceeds its initial value (not
#      monotone: a few-percent transient while the peak spreads into outgoing waves), ends below DECAY_MAX of it and
#      below RATIO_MAX of the run without cleaning; int B conserved to CONSERVED_TOL; max |psi| below PSI_MAX; the
#      saved pt, beta, bmag equal their definitions and the saved max|divb| equals the history (DERIVED_TOL);
#   2. AMR ([-0.5, 0.5]^2 refined once more, 8^2 cells per base block), GLM: L1 below DECAY_MAX_AMR of its initial
#      value, int B conserved, psi bounded, the seam-local column zero at step 0, never above the global maximum,
#      positive at the end;
#   3. monitor, 5 steps without cleaning, divb_tol = 1 (initial max|div B| = 3.0): divb_error = .false. warns,
#      divb_error = .true. stops the run with its message.
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
DECAY_MAX="0.1"       # GLM final / initial L1 of div B (CPU baseline 0.057, M2-P4b)
RATIO_MAX="0.1"       # GLM / none final L1 of div B (CPU baseline 0.027)
DECAY_MAX_AMR="0.2"   # AMR GLM final / initial (CPU baseline 0.091: the 2:1 seam injects div B, issue #29)
CONSERVED_TOL="1e-11" # int B_x, B_y, B_z (int B_z = 1.13; measured drift <= 8.3e-13)
PSI_MAX="0.141"       # c_h A / 10, A the peak amplitude (CPU baseline 3.6e-3 uniform, 1.3e-2 AMR)
DERIVED_TOL="1e-12"   # saved derived fields vs definitions and history (CPU: exact)

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
ORACLE="$CASE_DIR/divb_peak_oracle.py"
TAG="$(basename "$EXE")-np$NP"
FAILED=0

prepare() { # prepare <name> <make_divb_peak.py options...>: write the input of one case, print its work directory
   local w="$CASE_DIR/work-$TAG-$1"
   shift
   rm -rf "$w" ; mkdir -p "$w"
   "$VENV_PY" "$CASE_DIR/make_divb_peak.py" "$VERIF_DIR/conservation/amr-periodic.ini" "$w/divb-peak.ini" "$@"
   echo "$w"
}
case_run() { # case_run <name> <options...>: run one case, print its work directory
   local w
   w="$(prepare "$@")"
   if ! (cd "$w" && mpirun -np "$NP" "$EXE" divb-peak.ini > log.txt 2>&1); then
      echo "check.sh: div(B) peak $(basename "$w") run failed, see $w/log.txt" >&2
      return 1
   fi
   echo "$w"
}

echo ">> MV-10 Dedner peak, uniform 64^2, GLM and none ($(basename "$EXE"), np $NP)"
WG="$(case_run glm --divergence-control glm --save-aux)" || exit 1
WN="$(case_run none --divergence-control none)" || exit 1
"$VENV_PY" "$ORACLE" --decay "$WG" "$WN" --decay-max "$DECAY_MAX" --ratio-max "$RATIO_MAX" || FAILED=1
"$VENV_PY" "$ORACLE" "$WG" --conserved "$CONSERVED_TOL" --psi-max "$PSI_MAX" --derived --derived-tol "$DERIVED_TOL" \
   || FAILED=1
"$VENV_PY" "$ORACLE" "$WN" --conserved "$CONSERVED_TOL" || FAILED=1
echo ">> MV-10 Dedner peak, 2:1 AMR box, GLM"
WA="$(case_run amr --divergence-control glm --amr --cells-per-block 8)" || exit 1
"$VENV_PY" "$ORACLE" --decay "$WA" --decay-max "$DECAY_MAX_AMR" || FAILED=1
"$VENV_PY" "$ORACLE" "$WA" --conserved "$CONSERVED_TOL" --psi-max "$PSI_MAX" --seam || FAILED=1
for w in "$WG" "$WN" "$WA"; do find "$w" -name '*.h5' -delete; done

echo ">> MV-10 divb_tol monitor, 5 steps without cleaning, divb_tol = 1"
WW="$(case_run warn --divergence-control none --it-max 5 --divb-tol 1.0)" || exit 1
n=$(grep -c 'warning: max|div B|' "$WW/log.txt" || true)
if [[ "$n" -gt 0 ]]; then echo "   divb_error = .false.: $n warnings, the run completes  PASS"
else echo "   divb_error = .false.: no warning in $WW/log.txt  FAIL" ; FAILED=1 ; fi
WS="$(prepare stop --divergence-control none --it-max 5 --divb-tol 1.0 --divb-error)"
if (cd "$WS" && mpirun -np "$NP" "$EXE" divb-peak.ini > log.txt 2>&1); then
   echo "   divb_error = .true.: the run completed instead of stopping  FAIL" ; FAILED=1
elif grep -q 'max|div B| = .*divb_error' "$WS/log.txt"; then
   echo "   divb_error = .true.: stops with: $(grep -m1 'max|div B|' "$WS/log.txt" | sed 's/^.*error stop : //' | cut -c1-90)  PASS"
else
   echo "   divb_error = .true.: the run failed without the expected message, see $WS/log.txt  FAIL" ; FAILED=1
fi
for w in "$WW" "$WS"; do find "$w" -name '*.h5' -delete; done

if [[ $FAILED -eq 0 ]]; then
   echo "MV-10 PASSED ($TAG)"
else
   echo "MV-10 FAILED ($TAG)"
   exit 1
fi
