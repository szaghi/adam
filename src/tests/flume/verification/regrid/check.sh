#!/usr/bin/env bash
# FLUME verification RG: runtime regridding (issue #74, M5).
#
# Why: runtime AMR changes the grid during the run; every piece that was built once (the 2:1 flux register, the ghost
# and boundary maps, the device copies, the grid-dependent state) must follow, and a regrid must not change the
# conserved integrals. The legs grow with the phases of #74:
#   rg0  (P0, P1) the library regrid round trip, tests/amr/test_amr_regrid_roundtrip on 1, 2 and 3 ranks: refine a
#        block, derefine it back, with the linear and the conservative prolongation; the leaves and the blocks summed
#        over the ranks consistent;
#   rg1  (P0) the [amr] frequency contract on the sod-x regression input: 0 runs (no runtime regridding), a negative
#        value is fatal, n > 0 is refused until P2 lands the regrid (each with its message);
#   rg2  (P1) the [amr] regrid_prolongation contract on sod-x: absent resolves to conservative (the FLUME default),
#        linear and conservative are taken as given, any other value is refused.
#
# The P1 legs of rg0 also cover the conservative prolongation: the round trip exact and the integral unchanged on
# every field, linear data reproduced, a steep front kept in the range of its data (no negative child).
#
# Usage: ./check.sh [--leg rg0|rg1|rg2 ...]
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
UNIT="$REPO_ROOT/exe/test_amr_regrid_roundtrip"
LEGS=()
while [[ $# -gt 0 ]]; do
   case "$1" in
      --leg) LEGS+=("$2") ; shift 2 ;;
      *)     echo "check.sh: unknown argument '$1' (accepted: --leg L)" >&2 ; exit 2 ;;
   esac
done
[[ ${#LEGS[@]} -eq 0 ]] && LEGS=(rg0 rg1 rg2)
TAG="$(basename "$EXE")"
FAILED=0

expect_refused() { # expect_refused <work> <message regex>: the run must fail and print the message
   local work="$1" msg="$2"
   if (cd "$work" && mpirun -np 2 "$EXE" input.ini < /dev/null > log.txt 2>&1); then
      echo "   run succeeded, expected a refusal: see $work/log.txt" ; FAILED=1
   elif grep -aqE "$msg" "$work/log.txt"; then
      echo "   refused as expected"
   else
      echo "   failed without the expected message ($msg): see $work/log.txt" ; FAILED=1
   fi
}

echo ">> RG: runtime regridding ($TAG)"
for leg in "${LEGS[@]}"; do
   case "$leg" in
      rg0)
         echo "-- rg0: library regrid round trip"
         if [[ ! -x "$UNIT" ]]; then
            echo "   $UNIT not found (fobis build --mode test-amr-regrid-roundtrip-gnu)" ; FAILED=1 ; continue
         fi
         work="$CASE_DIR/work-unit" ; rm -rf "$work" ; mkdir -p "$work"
         for np in 1 2 3; do
            if (cd "$work" && mpirun -np "$np" --oversubscribe "$UNIT" < /dev/null > "log-np$np.txt" 2>&1) && \
               grep -q "test_amr_regrid_roundtrip: PASSED" "$work/log-np$np.txt"; then
               echo "   np $np: PASS"
            else
               echo "   np $np: FAIL, see $work/log-np$np.txt" ; FAILED=1
            fi
         done ;;
      rg1)
         base="$REPO_ROOT/src/tests/flume/regression/sod-x/input.ini"
         for f in 0 -1 5; do
            echo "-- rg1: [amr] frequency = $f"
            work="$CASE_DIR/work-$TAG-frequency$f" ; rm -rf "$work" ; mkdir -p "$work"
            sed -E "s/^(\s*frequency\s*=\s*)0\b/\1$f/; s/^(\s*time_max\s*=).*/\1 0.002/" "$base" > "$work/input.ini"
            case "$f" in
               0)  if (cd "$work" && mpirun -np 2 "$EXE" input.ini < /dev/null > log.txt 2>&1); then echo "   runs"
                   else echo "   run failed: see $work/log.txt" ; FAILED=1 ; fi ;;
               -1) expect_refused "$work" "must be 0 \(no runtime regridding\) or positive" ;;
               5)  expect_refused "$work" "asks for runtime regridding, which is not implemented yet" ;;
            esac
         done ;;
      rg2)
         base="$REPO_ROOT/src/tests/flume/regression/sod-x/input.ini"
         for p in absent linear conservative bogus; do
            echo "-- rg2: [amr] regrid_prolongation = $p"
            work="$CASE_DIR/work-$TAG-prolongation-$p" ; rm -rf "$work" ; mkdir -p "$work"
            if [[ "$p" == absent ]]; then
               sed -E "s/^(\s*time_max\s*=).*/\1 0.002/" "$base" > "$work/input.ini"
               want=conservative
            else
               sed -E "s/^(\s*frequency\s*=.*)$/\1\nregrid_prolongation = $p/; s/^(\s*time_max\s*=).*/\1 0.002/" \
                  "$base" > "$work/input.ini"
               want=$p
            fi
            if [[ "$p" == bogus ]]; then
               expect_refused "$work" "unknown \[amr\]\.\(regrid_prolongation\) \"bogus\""
            elif (cd "$work" && mpirun -np 2 "$EXE" input.ini < /dev/null > log.txt 2>&1); then
               if grep -aq "flume: \[amr\] regrid_prolongation = $want\$" "$work/log.txt"; then
                  echo "   runs, resolved to $want"
               else
                  echo "   runs, but not resolved to $want: see $work/log.txt" ; FAILED=1
               fi
            else
               echo "   run failed: see $work/log.txt" ; FAILED=1
            fi
         done ;;
      *) echo "check.sh: unknown leg '$leg'" >&2 ; exit 2 ;;
   esac
done

if [[ $FAILED -eq 0 ]]; then
   echo "RG PASSED ($TAG)"
else
   echo "RG FAILED ($TAG)"
   exit 1
fi
