#!/usr/bin/env bash
# pic-pec-t0: the t0 electric field of a PIC plasma between PEC walls (issue #84, E3b).
#
# THE CLAIM UNDER TEST: FLAIL's t0 solve D(G) phi = -rho/eps, then D <- D - eps G phi, gives div D = rho with PRISM's own
# divergence, to round-off, up to the one component no D can carry: PRISM's centred div and grad of order p annihilate
# the checkerboard chi = (-1)^(i+j+k), so with PEC faces (the potential odd about them) D(G) is singular on chi and
# div D - rho = -c chi, c = <rho, chi> / N (t0_elliptic.py gauss). The legacy smoother left 3.4e-5 max|rho| beyond
# that (measured in #84 E3b); here it is ~1e-16.
#
# ACCEPTANCE, at -np 1 and -np 2:
#   1. the run completes, and both t0 solves ran on FLAIL (no legacy fallback);
#   2. max|div D - rho + c chi| <= 1e-12 max|rho| on the step-0 checkpoint;
#   3. the step-0 D and B of the two rank counts agree to 1e-12 (the solves gather to rank 0).
#
# Usage: ./check.sh [--build]
# PRISM_EXE selects the executable (default exe/adam_prism_cpu), PRISM_PY a python with h5py (default the regression
# venv exe/.regression-venv, which run.sh creates). mpirun must be on PATH.
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${PRISM_EXE:-$REPO_ROOT/exe/adam_prism_cpu}"
BACKEND="${EXE##*_}"
PY="${PRISM_PY:-$REPO_ROOT/exe/.regression-venv/bin/python}"
CHECK="$CASE_DIR/../t0_elliptic.py"
BASE=pic_pec_t0

if [[ "${1:-}" == "--build" ]]; then
   echo ">> building prism-cpu-gnu"
   (cd "$REPO_ROOT" && fobis build --mode prism-cpu-gnu)
fi
[[ -x "$EXE" ]] || { echo "ERROR: $EXE not found — run with --build" >&2; exit 2; }
[[ -x "$PY" ]] || { echo "ERROR: $PY not found (run run.sh once, or set PRISM_PY)" >&2; exit 2; }
command -v mpirun >/dev/null 2>&1 || { echo "ERROR: mpirun not on PATH" >&2; exit 2; }

fail=0
for np in 1 2; do
   wd="$CASE_DIR/work-check-$BACKEND-np$np"
   rm -rf "$wd" && mkdir -p "$wd"
   cp "$CASE_DIR/input.ini" "$wd/"
   if ! ( cd "$wd" && timeout 600 mpirun -np "$np" "$EXE" > run.log 2>&1 ); then
      echo "FAIL [pic-pec-t0/$BACKEND/np$np] the run exited non-zero, see $wd/run.log"; fail=1; continue
   fi
   n="$(grep -ac 'FLAIL t0 solve' "$wd/run.log" || true)"
   if [[ "$n" -ne 4 ]]; then
      echo "FAIL [pic-pec-t0/$BACKEND/np$np] expected 4 FLAIL t0 solves (D, then A x3), found $n"; fail=1
   fi
   "$PY" "$CHECK" gauss "$wd" "$BASE" || { echo "FAIL [pic-pec-t0/$BACKEND/np$np] Gauss law at t0"; fail=1; }
done
if [[ $fail -eq 0 ]]; then
   "$PY" "$CHECK" same "$CASE_DIR/work-check-$BACKEND-np1" "$CASE_DIR/work-check-$BACKEND-np2" "$BASE" Dx Dy Dz Bx By Bz \
      || { echo "FAIL [pic-pec-t0/$BACKEND] the t0 fields depend on the rank count"; fail=1; }
fi
[[ $fail -eq 0 ]] && echo "PASS [pic-pec-t0/$BACKEND] t0 Gauss law (np 1, 2) and rank independence" || echo "FAIL [pic-pec-t0/$BACKEND]"
exit $fail
