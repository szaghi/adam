#!/usr/bin/env bash
# FLUME verification AV: accuracy and robustness of runtime regridding (issue #74, M5-P4).
#
# Why: RG (check.sh) proves the regrid machinery (round trip, conservation, restart, CPU/FNL agreement on one case);
# AV asks whether a run that regrids as the solution moves is about as accurate as the uniform run at its finest level
# while holding fewer cells, on the problems FLUME solves, and whether a regrid keeps what the scheme guarantees.
# Bounds stated in the #74 P4 plan before the runs; a measured value outside one is reported, never moved. Inputs from
# make_regrid.py, checks by av_oracle.py (reusing the V1, V2, V6 and MV-9 oracles):
#   av2  isentropic vortex (V2), base N = 32 tracked to N = 128 (gradient marker on rho, tol 0.05), a regrid every 5
#        steps, t = 0.2: L1(rho) <= 1.5x the uniform N = 128 run and <= 1/8 of the uniform N = 32 run;
#   av3a Sod (V1) on a quadtree with y null, base 48 cells tracked to 192 (gradient OR Loehner on rho), on [0, 1.2] so
#        that the jump at 0.5 is not a block face (issue #76): L1(rho) <= 1.3x the uniform 192-cell run, the contact
#        and the shock on finest cells at t = 0.2, the null copies within 1e-12;
#   av3b Balsara-Spicer MHD blast (the blast-limiter golden: EGLM, positivity limiter), base N = 64 tracked to N = 128
#        (gradient OR Loehner on p): the outer shock radius on the four half-axes within 2 finest cells of the uniform
#        N = 128 run, rho conserved to 1e-13, rhoE and B drifting at most 1.1x the uniform run (the scheme itself does
#        not conserve them here); every regrid checks the new state admissible. Point symmetry is not asserted: the
#        uniform run breaks it too (issue #77);
#   av4  field loop (MV-9) with GLM and with EGLM, base N = 64 tracked to N = 128 (gradient on B_x OR Loehner on B_x
#        and B_y, floor 1e-4), t = 1: E_B(T)/E_B(0) >= the uniform N = 128 run's - 0.5 %, <|B_z|>/A0 <= 1.2x, max
#        |div B| at the end <= 2x (EGLM: <|B_z|> reported only, its seam level is 70x its uniform one even without
#        any regrid, issue #78); the jump of max |div B| at each regrid reported;
#   av6  rank invariance: the av2 tracked vortex on 1, 2 and 3 ranks (1 and 2 on FNL: the WSL development box has two
#        GPUs, and three ranks sharing them run ~170x slower): the same regrids, the fields within 1e-13;
#   av7  shock over the cylinder (V6) to t = 0.25, tracked (solid OR Loehner on rho): the bow shock stand-off within 2
#        finest cells of the uniform run at V6's finest level; the V6 oracle (blocks crossed by the surface refined,
#        mirror symmetry 1e-10, positivity) on the tracked run; the V6 init-AMR run's stand-off reported;
#   av8  CPU against FNL: the tracked runs of av2, av3a, av3b, av4 and av7 of both backends (run the legs with each
#        executable first) regrid alike, and per field their CPU-FNL difference is at most 2x the CPU-FNL difference
#        of the uniform-fine runs (or below 1e-10): the regrid adds no divergence between the backends. A fixed
#        tolerance would not test the regrid: on MHD the backends differ without any regrid (B by ~1e-6 on the field
#        loop, the blast at O(1), issues #77 and #79), on Euler they agree to round-off.
#
# Usage: ./accuracy.sh [--leg av2|av3a|av3b|av4|av6|av7|av8 ...]
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./accuracy.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
MAKE=(python3 -I "$CASE_DIR/make_regrid.py")
ORACLE=("$VENV_PY" "$CASE_DIR/av_oracle.py")
CONSERVATION="$VERIF_DIR/conservation/conservation_oracle.py"
LEGS=()
while [[ $# -gt 0 ]]; do
   case "$1" in
      --leg) LEGS+=("$2") ; shift 2 ;;
      *)     echo "accuracy.sh: unknown argument '$1' (accepted: --leg L)" >&2 ; exit 2 ;;
   esac
done
[[ ${#LEGS[@]} -eq 0 ]] && LEGS=(av2 av3a av3b av4 av6 av7)
TAG="$(basename "$EXE")"
FAILED=0
UNTIMED='s/, [^,]* s (host round trip)$//' # the FNL regrid log lines end with their wall time

case_run() { # case_run <work> <ini> [np]: run one case, report a failure, return its status
   local work="$1" ini="$2" np="${3:-2}" start
   start=$(date +%s)
   if (cd "$work" && mpirun -np "$np" --oversubscribe "$EXE" "$ini" < /dev/null > log.txt 2>&1); then
      echo "   $(basename "$work"): $(( $(date +%s) - start )) s, $(grep -ac 'flume: regrid at' "$work/log.txt") regrids"
      return 0
   fi
   echo "   run failed: see $work/log.txt" ; FAILED=1 ; return 1
}
prepare() { # prepare <name> <case> <ini> [make_regrid options]: a fresh work dir with its input, path on stdout
   local name="$1" kind="$2" ini="$3" work ; shift 3
   work="$CASE_DIR/work-$TAG-$name" ; rm -rf "${work:?}" ; mkdir -p "$work"
   "${MAKE[@]}" "$kind" "$work/$ini" "$@"
   echo "$work"
}

echo ">> AV: runtime AMR accuracy ($TAG)"
for leg in "${LEGS[@]}"; do
   case "$leg" in
      av2)
         echo "-- av2: isentropic vortex, N = 32 tracked to N = 128 against uniform 32 and 128"
         c=$(prepare av2-u32 vortex-av vortex-av.ini --base-level 1 --max-level 1 --frequency 0)
         f=$(prepare av2-u128 vortex-av vortex-av.ini --base-level 3 --max-level 3 --frequency 0)
         t=$(prepare av2-tracked vortex-av vortex-av.ini --base-level 1 --max-level 3 --frequency 5 --tol 0.05)
         case_run "$c" vortex-av.ini && case_run "$f" vortex-av.ini && case_run "$t" vortex-av.ini || continue
         "${ORACLE[@]}" vortex "$t" --fine "$f" --coarse "$c" || FAILED=1 ;;
      av3a)
         echo "-- av3a: Sod, 48 cells tracked to 192 against uniform 192"
         f=$(prepare av3a-u192 sod-av sod-av.ini --base-level 4 --max-level 4 --frequency 0)
         t=$(prepare av3a-tracked sod-av sod-av.ini --base-level 2 --max-level 4 --frequency 5 --tol 1.0 --refine 0.5 \
                                                    --derefine 0.2)
         case_run "$f" sod-av.ini && case_run "$t" sod-av.ini || continue
         "${ORACLE[@]}" sod "$t" --fine "$f" || FAILED=1 ;;
      av3b)
         echo "-- av3b: MHD blast, N = 64 tracked to N = 128 against uniform 128"
         f=$(prepare av3b-u128 blast-av blast-av.ini --base-level 3 --max-level 3 --frequency 0 --time-max 0.01)
         t=$(prepare av3b-tracked blast-av blast-av.ini --base-level 2 --max-level 3 --frequency 5 --time-max 0.01 \
                                                        --tol 100.0 --refine 0.5 --derefine 0.2)
         case_run "$f" blast-av.ini && case_run "$t" blast-av.ini || continue
         "${ORACLE[@]}" blast "$t" --fine "$f" --no-symmetry || FAILED=1
         "${ORACLE[@]}" drift "$t" --fine "$f" --conserved r --ratio 1.1 || FAILED=1 ;;
      av4)
         for cleaning in glm eglm; do
            extra=() ; bz=() ; [[ "$cleaning" == eglm ]] && { extra=(--eglm) ; bz=(--bz-report) ; }
            echo "-- av4: field loop ($cleaning), N = 64 tracked to N = 128 against uniform 128"
            f=$(prepare "av4-$cleaning-u128" loop-av field-loop.ini --cells 128 --frequency 0 --time-max 1.0 "${extra[@]}")
            t=$(prepare "av4-$cleaning-tracked" loop-av field-loop.ini --cells 64 --max-level 3 --frequency 5 \
                                                --tol 0.01 --time-max 1.0 "${extra[@]}")
            case_run "$f" field-loop.ini && case_run "$t" field-loop.ini || continue
            "${ORACLE[@]}" loop "$t" --fine "$f" "${bz[@]}" || FAILED=1
         done ;;
      av6)
         ranks=(1 2 3) ; [[ "$TAG" == *fnl* ]] && ranks=(1 2)
         echo "-- av6: the av2 tracked vortex on ${ranks[*]} ranks"
         for np in "${ranks[@]}"; do
            w=$(prepare "av6-np$np" vortex-av vortex-av.ini --base-level 1 --max-level 3 --frequency 5 --tol 0.05)
            case_run "$w" vortex-av.ini "$np" || continue 2
         done
         for np in "${ranks[@]:1}"; do
            if diff <(grep -a "flume: regrid at" "$CASE_DIR/work-$TAG-av6-np1/log.txt" | sed "$UNTIMED") \
                    <(grep -a "flume: regrid at" "$CASE_DIR/work-$TAG-av6-np$np/log.txt" | sed "$UNTIMED") > /dev/null; then
               echo "   np $np: the same regrids as np 1"
            else
               echo "   np $np: the regrids differ from np 1" ; FAILED=1
            fi
            "$VENV_PY" "$CONSERVATION" --compare "$CASE_DIR/work-$TAG-av6-np1" "$CASE_DIR/work-$TAG-av6-np$np" \
                                       --tol 1e-13 --ngc 3 || FAILED=1
         done ;;
      av7)
         echo "-- av7: shock over the cylinder to t = 0.25, tracked against uniform at V6's finest level"
         f=$(prepare av7-uniform cylinder-av shock-cylinder.ini --uniform)
         t=$(prepare av7-tracked cylinder-av shock-cylinder.ini --frequency 5)
         i="$CASE_DIR/work-$TAG-av7-init-amr" ; rm -rf "${i:?}" ; mkdir -p "$i"
         cp "$VERIF_DIR/shock-cylinder/shock-cylinder.ini" "$i/"
         case_run "$f" shock-cylinder.ini && case_run "$t" shock-cylinder.ini && case_run "$i" shock-cylinder.ini \
            || continue
         "${ORACLE[@]}" cylinder "$t" --fine "$f" || FAILED=1
         "${ORACLE[@]}" cylinder "$i" --fine "$f" --report-only
         "$VENV_PY" "$VERIF_DIR/shock-cylinder/shock_cylinder_oracle.py" "$t" --no-refinement-check || FAILED=1 ;;
      av8)
         echo "-- av8: CPU against FNL on the tracked runs"
         for pair in av2-tracked:av2-u128 av3a-tracked:av3a-u192 av3b-tracked:av3b-u128 \
                     av4-glm-tracked:av4-glm-u128 av4-eglm-tracked:av4-eglm-u128 av7-tracked:av7-uniform; do
            name=${pair%%:*} ; uniform=${pair#*:}
            c="$CASE_DIR/work-adam_flume_cpu-$name" ; f="$CASE_DIR/work-adam_flume_fnl-$name"
            if [[ ! -f "$c/log.txt" || ! -f "$f/log.txt" ]]; then
               echo "   $name: missing a backend's run (run the legs with both executables first)" ; FAILED=1 ; continue
            fi
            if diff <(grep -a "flume: regrid at" "$c/log.txt" | sed "$UNTIMED") \
                    <(grep -a "flume: regrid at" "$f/log.txt" | sed "$UNTIMED") > /dev/null; then
               echo "   $name: the same regrids"
            else
               echo "   $name: the regrids differ" ; FAILED=1
            fi
            "${ORACLE[@]}" agree "$c" "$f" "$CASE_DIR/work-adam_flume_cpu-$uniform" \
                                 "$CASE_DIR/work-adam_flume_fnl-$uniform" || FAILED=1
         done ;;
      *) echo "accuracy.sh: unknown leg '$leg'" >&2 ; exit 2 ;;
   esac
done

if [[ $FAILED -eq 0 ]]; then
   echo "AV PASSED ($TAG)"
else
   echo "AV FAILED ($TAG)"
   exit 1
fi
