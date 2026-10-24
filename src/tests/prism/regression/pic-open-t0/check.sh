#!/usr/bin/env bash
# pic-open-t0: the t0 electric field of a PIC plasma inside exact/open faces (PRISM's Neumann; issue #84, E3c).
#
# THE CLAIM UNDER TEST: an exact/open face holds, in every ghost layer, the continuum free-space potential of the source,
# -sum f dV / (4 pi r) (FLAIL_BC_FREE_GHOSTS). FLAIL sums it on rank 0 over the gathered source, so the t0 fields do not
# depend on the rank count (the legacy smoother summed over each rank's own blocks; this nearly neutral plasma, far from
# the faces, cannot show the difference: rmf-open-t0 does). The operator with prescribed ghosts is nonsingular, so div D = rho to the CG tolerance on the
# cells that read no field-BC ghost (t0_elliptic.py open, margin = hs = 3).
#
# ACCEPTANCE, at -np 1 and -np 2:
#   1. the run completes, and both t0 solves ran on FLAIL (no legacy fallback);
#   2. max|div D - rho| <= 1e-11 max|rho| beyond 3 cells of the domain faces, on the step-0 checkpoint (measured 6.6e-13;
#      the legacy smoother, 10000 sweeps, reached 4.6e-16 on this source);
#   3. the step-0 D and B of the two rank counts agree to 1e-12.
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
BASE=pic_open_t0

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
      echo "FAIL [pic-open-t0/$BACKEND/np$np] the run exited non-zero, see $wd/run.log"; fail=1; continue
   fi
   n="$(grep -ac 'FLAIL t0 solve' "$wd/run.log" || true)"
   if [[ "$n" -ne 4 ]]; then
      echo "FAIL [pic-open-t0/$BACKEND/np$np] expected 4 FLAIL t0 solves (D, then A x3), found $n"; fail=1
   fi
   "$PY" "$CHECK" open "$wd" "$BASE" margin=3 || { echo "FAIL [pic-open-t0/$BACKEND/np$np] Gauss law at t0"; fail=1; }
done
if [[ $fail -eq 0 ]]; then
   "$PY" "$CHECK" same "$CASE_DIR/work-check-$BACKEND-np1" "$CASE_DIR/work-check-$BACKEND-np2" "$BASE" Dx Dy Dz Bx By Bz \
      || { echo "FAIL [pic-open-t0/$BACKEND] the t0 fields depend on the rank count"; fail=1; }
fi
[[ $fail -eq 0 ]] && echo "PASS [pic-open-t0/$BACKEND] t0 Gauss law (np 1, 2) and rank independence (free-space ghosts)" || echo "FAIL [pic-open-t0/$BACKEND]"
exit $fail
