#!/usr/bin/env bash
# FLUME verification RG: runtime regridding (issue #74, M5).
#
# Why: runtime AMR changes the grid during the run; every piece that was built once (the 2:1 flux register, the ghost
# and boundary maps, the device copies, the grid-dependent state) must follow, and a regrid must not change the
# conserved integrals. The legs grow with the phases of #74:
#   rg0  (P0, P1) the library regrid round trip, tests/amr/test_amr_regrid_roundtrip on 1, 2 and 3 ranks: refine a
#        block, derefine it back, with the linear and the conservative prolongation; the leaves and the blocks summed
#        over the ranks consistent;
#   rg1  (P0, P2) the [amr] frequency contract on the sod-x regression input: 0 runs (no runtime regridding), a
#        negative value is fatal, n > 0 without markers is fatal (nothing would mark a block), each with its message;
#   rg2  (P1) the [amr] regrid_prolongation contract on sod-x: absent resolves to conservative (the FLUME default),
#        linear and conservative are taken as given, any other value is refused;
#   rg3  (P2) conservation through regrids: the isentropic vortex travelling across a periodic quadtree, a Loehner
#        marker on the density, a regrid every 5 steps, reflux on (make_regrid.py vortex). The grid must refine and
#        derefine during the run, and with the conservative prolongation every volume integral must stay constant to
#        round-off; with the linear one they must drift (the negative control: else the regrids move no data).
#        The flux register is rebuilt at every regrid, so this leg also covers the reflux on seams that come and go;
#   rg4  (P2) restart across regrids: 30 steps against 20 + restart + 10, the restart saved at step 20, a regrid
#        step that changes the grid; the last fields bitwise and the histories identical (the forest regrids before
#        post_step, so the restart holds the grid the next step runs on);
#   rg5  (P2) immersed solid: the shock over the cylinder with the solid marker and a Loehner marker, a regrid every 5
#        steps (make_regrid.py cylinder): the run completes with at least one regrid that changes the grid, the
#        distance function recomputed on it, the state admissible (each regrid checks it) and finite.
#
# On the FNL backend rg3 checks the refusal of runtime regridding (it lands in #74 P3); rg4 and rg5 are CPU only.
#
# The P1 legs of rg0 also cover the conservative prolongation: the round trip exact and the integral unchanged on
# every field, linear data reproduced, a steep front kept in the range of its data (no negative child).
#
# Usage: ./check.sh [--leg rg0|rg1|rg2|rg3|rg4|rg5 ...]
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
[[ ${#LEGS[@]} -eq 0 ]] && LEGS=(rg0 rg1 rg2 rg3 rg4 rg5)
TAG="$(basename "$EXE")"
FAILED=0
MAKE=(python3 -I "$CASE_DIR/make_regrid.py")
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
ORACLE="$REPO_ROOT/src/tests/flume/verification/conservation/conservation_oracle.py"
IS_FNL=0 ; [[ "$TAG" == *fnl* ]] && IS_FNL=1

run_case() { # run_case <work> [log]: run input.ini on 2 ranks, report a failure
   local work="$1" log="${2:-log.txt}"
   if (cd "$work" && mpirun -np 2 "$EXE" input.ini < /dev/null > "$log" 2>&1); then
      return 0
   fi
   echo "   run failed: see $work/$log" ; FAILED=1 ; return 1
}
regrids() { # regrids <log>: print the regrid lines and the totals of blocks refined and families coarsened
   grep -a "flume: regrid at" "$1" | sed 's/^/   /'
   grep -a "flume: regrid at" "$1" | awk '{for (i=1;i<=NF;i++) {if ($i ~ /^refined/) r+=$(i-1); if ($i ~ /^coarsened/) d+=$(i-1)}} END {print r+0, d+0}'
}

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
               5)  expect_refused "$work" "nothing would mark a block" ;;
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
      rg3)
         if [[ $IS_FNL -eq 1 ]]; then
            echo "-- rg3: runtime regridding refused on FNL (lands in #74 P3)"
            work="$CASE_DIR/work-$TAG-vortex-refused" ; rm -rf "$work" ; mkdir -p "$work"
            "${MAKE[@]}" vortex "$work/input.ini"
            expect_refused "$work" "which the FNL backend does not do yet"
            continue
         fi
         for p in conservative linear; do
            echo "-- rg3: travelling vortex, regrid every 5 steps, $p prolongation"
            work="$CASE_DIR/work-$TAG-vortex-$p" ; rm -rf "$work" ; mkdir -p "$work"
            "${MAKE[@]}" vortex "$work/input.ini" --prolongation "$p"
            run_case "$work" || continue
            read -r nref nder < <(regrids "$work/log.txt" | tail -1)
            regrids "$work/log.txt" | sed '$d'
            if [[ $nref -gt 0 && $nder -gt 0 ]]; then
               echo "   $nref blocks refined, $nder families coarsened over the run"
            else
               echo "   the grid did not both refine and coarsen ($nref refined, $nder coarsened)" ; FAILED=1
            fi
         done
         "$VENV_PY" "$ORACLE" --conserved "$CASE_DIR/work-$TAG-vortex-conservative" \
                              --leaky "$CASE_DIR/work-$TAG-vortex-linear" || FAILED=1 ;;
      rg4)
         [[ $IS_FNL -eq 1 ]] && { echo "-- rg4: CPU only until #74 P3" ; continue ; }
         echo "-- rg4: restart across regrids, 30 steps vs 20 + restart + 10"
         a="$CASE_DIR/work-$TAG-restart-A" ; b="$CASE_DIR/work-$TAG-restart-B"
         rm -rf "$a" "$b" ; mkdir -p "$a" "$b"
         "${MAKE[@]}" vortex "$a/input.ini" --it-max 30
         run_case "$a" || continue
         "${MAKE[@]}" vortex "$b/input.ini" --it-max 20 --restart-save 20
         run_case "$b" log-1.txt || continue
         if grep -aq "flume: regrid at step 20: .* [1-9][0-9]* refined\|flume: regrid at step 20: .* [1-9][0-9]* coarsened" \
               "$b/log-1.txt"; then
            echo "   step 20 regrids and changes the grid"
         else
            echo "   step 20 does not change the grid: the restart does not test the ordering" ; FAILED=1
         fi
         "${MAKE[@]}" vortex "$b/input.ini" --it-max 30 --restart
         run_case "$b" log-2.txt || continue
         "$VENV_PY" "$ORACLE" --compare "$a" "$b" --tol 0 --ngc 3 || FAILED=1
         for f in $(cd "$a" && ls ./*.dat); do
            if cmp -s "$a/$f" "$b/$f"; then echo "   ${f#./}: identical"
            else echo "   ${f#./}: differs" ; FAILED=1 ; fi
         done ;;
      rg5)
         [[ $IS_FNL -eq 1 ]] && { echo "-- rg5: CPU only until #74 P3" ; continue ; }
         echo "-- rg5: shock over the cylinder, solid + Loehner markers, regrid every 5 steps"
         work="$CASE_DIR/work-$TAG-cylinder" ; rm -rf "$work" ; mkdir -p "$work"
         "${MAKE[@]}" cylinder "$work/input.ini" --it-max 60
         run_case "$work" || continue
         read -r nref nder < <(regrids "$work/log.txt" | tail -1)
         regrids "$work/log.txt" | sed '$d'
         if [[ $nref -gt 0 ]]; then echo "   completes, $nref blocks refined, $nder families coarsened"
         else echo "   no regrid changed the grid" ; FAILED=1 ; fi ;;
      *) echo "check.sh: unknown leg '$leg'" >&2 ; exit 2 ;;
   esac
done

if [[ $FAILED -eq 0 ]]; then
   echo "RG PASSED ($TAG)"
else
   echo "RG FAILED ($TAG)"
   exit 1
fi
