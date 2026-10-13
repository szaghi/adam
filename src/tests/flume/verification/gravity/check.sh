#!/usr/bin/env bash
# FLUME verification GV (issue #73, M7): stratified atmospheres under constant gravity.
#
# M7-G0 (the inputs, before the gravity source of G1): every case of make_atmosphere.py runs one negligible step
# (CFL = 1e-30) with the ghosts poisoned, and atmosphere_oracle.py profile checks the checkpoint:
#   - the interior against the analytic atmosphere (isothermal or polytropic, region 1 at the origin);
#   - the ghosts beyond the `hydrostatic` faces against the boundary rule (the isothermal continuation from the first
#     interior cell along the normal; for the isothermal atmosphere the analytic profile itself, edges and corners too),
#     the velocity zero, the field and psi even;
#   - the periodic ghosts against the analytic profile at the wrapped position.
# Cases: iso2d, poly2d, tilted2d (every face hydrostatic, a refined block), tilted3d, iso2d-mhd (EGLM, a field oblique
# to gravity), froude2d (the number form of [gravity]), on each of --np.
# Negative control: iso2d run with extrapolation on the y faces and judged as if they were hydrostatic must FAIL (the
# oracle tells the boundary rule from a copy).
# NV-5 (--reference): ../scaling/reference.sh --inputs on iso2d, iso2d-mhd and froude2d: the [gravity] components
# convert as an acceleration (u0^2 / L0), froude does not, and the dimensional run equals the base one bit for bit.
#
# Usage: ./check.sh [--np "1 2"] [--cases "iso2d ..."] [--reference] [--keep]
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
NPS="1 2"
CASES="iso2d poly2d tilted2d tilted3d iso2d-mhd froude2d"
REFERENCE=0
KEEP=0
while [[ $# -gt 0 ]]; do
   case "$1" in
      --np)        NPS="$2" ; shift 2 ;;
      --cases)     CASES="$2" ; shift 2 ;;
      --reference) REFERENCE=1 ; shift ;;
      --keep)      KEEP=1 ; shift ;;
      *) echo "check.sh: unknown argument '$1' (accepted: --np, --cases, --reference, --keep)" >&2 ; exit 2 ;;
   esac
done
TAG="$(basename "$EXE")"
FAILED=0

run_case() { # run_case <work> <np>: run the input of <work>, report a failure
   if ! (cd "$1" && mpirun -np "$2" "$EXE" input.ini < /dev/null > log.txt 2>&1); then
      echo "   run failed: see $1/log.txt" ; FAILED=1 ; return 1
   fi
}

echo ">> GV: stratified atmospheres ($TAG)"
for np in $NPS; do
   for case in $CASES; do
      work="$CASE_DIR/work-$TAG-np$np-$case"
      rm -rf "${work:?}" ; mkdir -p "$work"
      python3 -I "$CASE_DIR/make_atmosphere.py" "$work/input.ini" "$case"
      echo "-- $case, np $np"
      run_case "$work" "$np" || continue
      "$VENV_PY" "$CASE_DIR/atmosphere_oracle.py" profile "$work" || FAILED=1
      [[ $KEEP -eq 1 ]] || find "$work" -name '*.h5' -delete
   done
done

echo "-- negative control: iso2d with extrapolation on the y faces, judged as hydrostatic (must FAIL)"
work="$CASE_DIR/work-$TAG-negative"
rm -rf "${work:?}" ; mkdir -p "$work"
python3 -I "$CASE_DIR/make_atmosphere.py" "$work/input.ini" iso2d
sed -i '/^\[bc_y_m[ai][nx]\]/{n;s/^type = hydrostatic/type = extrapolation/}' "$work/input.ini"
if run_case "$work" 1; then
   sed -i '/^\[bc_y_m[ai][nx]\]/{n;s/^type = extrapolation/type = hydrostatic/}' "$work/input.ini"
   if "$VENV_PY" "$CASE_DIR/atmosphere_oracle.py" profile "$work" > "$work/oracle.txt"; then
      echo "   the oracle passed copied ghosts: FAIL" ; FAILED=1
   else
      echo "   the oracle rejects copied ghosts ($(grep -a 'hydrostatic ghosts rho' "$work/oracle.txt" | sed 's/^ *//')): PASS"
   fi
fi
[[ $KEEP -eq 1 ]] || find "$work" -name '*.h5' -delete

if [[ $REFERENCE -eq 1 ]]; then
   echo "-- NV-5: [reference] on iso2d, iso2d-mhd, froude2d"
   inputs=()
   mkdir -p "$CASE_DIR/work-$TAG-ref"
   for case in iso2d iso2d-mhd froude2d; do
      python3 -I "$CASE_DIR/make_atmosphere.py" "$CASE_DIR/work-$TAG-ref/$case.ini" "$case"
      inputs+=("$CASE_DIR/work-$TAG-ref/$case.ini")
   done
   (cd "$CASE_DIR/../scaling" && FLUME_EXE="$EXE" ./reference.sh --inputs "${inputs[*]}") || FAILED=1
fi

if [[ $FAILED -eq 0 ]]; then
   echo "GV PASSED ($TAG)"
else
   echo "GV FAILED ($TAG)"
fi
exit $FAILED
