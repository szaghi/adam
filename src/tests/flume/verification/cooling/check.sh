#!/usr/bin/env bash
# FLUME verification CO: optically thin radiative cooling, `[cooling]` (issue #93).
#
# Why: cooling is applied once per step by Townsend's exact integration (D-M11-1), so on a uniform box, where the flux
# divergence vanishes, every cell must follow the isochoric cooling ODE to round-off whatever the step, reach the floor
# T_1 exactly and stay there; the energy it removes is the radiated energy. make_co.py writes the box (rho = 1, p = 1,
# T = 1, the four-node table T = 0.1 0.3 1 3, Lambda = 0.2 1 0.5 1.5), co_oracle.py integrates the ODE by quadrature,
# independently of the closed form. The legs:
#   co1  the temperature of every step against the ODE to 1e-12 relative, and E plus the radiated energy constant to
#        1e-13, before the floor: with the CFL step (CFL 0.5) from T = 5 to t = 6.5 (the table continued above its
#        last node, then all three segments), and with the cooling limiter (dt_factor 0.01: every step within 0.01
#        t_cool, and the limit reached) from T = 1 to t = 2;
#   co2  the floor: to t = 3, T = T_1 exactly from the step that crosses the ODE's floor time on, E bitwise constant;
#   co3  the contract: [cooling] with mhd-ideal, non-increasing temperatures, lambdas of another length, more than 64
#        nodes, no particle_mass and a negative dt_factor are refused;
#   co4  the dual energy: the box moving at u = 10 (e = 0.03 E, below the dual-energy switch 0.05) with euler-de, to
#        t = 1.5: the temperature against the ODE, the entropy on the adiabat of the cooled energy every step (1e-13),
#        no energy added by the entropy floor (it may report lifts of round-off size, rho s -> p -> e not being an
#        exact round trip, each raising E by 0), and the flow history bitwise that of the same box with euler;
#   co5  the [reference] layer: the co1 limited box to t = 1 written in the units L0 = 2, u0 = 2, rho0 = 4
#        (scaling.py dimensionalize: the temperatures times u0^2, Lambda times L0^5 u0^3 rho0, the particle mass times
#        rho0 L0^3, the node lists value by value) runs bitwise as the code-units box (powers of two convert exactly);
#   co6  (P2) CPU and FNL agree on the co1 CFL and co4 euler-de runs (run co1 and co4 with each first): the same steps,
#        every row's time, temperature and radiated energy within 1e-12 relative (a uniform box: the runs without
#        cooling agree bitwise, so the AV av8 floor, 1e-10, is tightened to the oracle's 1e-12).
#
# Usage: ./check.sh [--legs LIST]   (legs co1 co2 co3 co4 co5 co6, all; default co1 to co5)
#
# FLUME_EXE overrides the executable under test (default exe/adam_flume_cpu); the caller owns the matching environment.
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
source "$VERIF_DIR/legs.sh"
legs_declare "co1 co2 co3 co4 co5 co6" "co1 co2 co3 co4 co5"
while [[ $# -gt 0 ]]; do
   case "$1" in
      --legs) legs_select "$2" ; shift 2 ;;
      *)      echo "check.sh: unknown argument '$1' (accepted: --legs LIST)" >&2 ; exit 2 ;;
   esac
done
legs_resolve check.sh
if [[ ! -x "$EXE" ]]; then echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2 ; exit 2 ; fi
EXE="$(realpath "$EXE")"
TAG="$(basename "$EXE")"
FAILED=0

run() { # run <work dir>: run the input of the work dir; 0 on success
   (cd "$1" && mpirun -np 2 "$EXE" input.ini < /dev/null > log.txt 2>&1)
}
make_case() { # make_case <name> [make_co.py options]: write one case, its dir in W
   local name=$1 ; shift
   W="$CASE_DIR/work-$TAG-$name"
   rm -rf "$W" ; mkdir -p "$W"
   python3 -I "$CASE_DIR/make_co.py" "$W/input.ini" "$@"
}
run_oracle() { # run_oracle <work dir> <co_oracle.py arguments>: run the case, then the oracle
   local w=$1 ; shift
   if run "$w"; then
      "$VENV_PY" "$CASE_DIR/co_oracle.py" "$@" || FAILED=1
   else
      echo "   the run failed: see $w/log.txt  FAIL" ; FAILED=1
   fi
}
expect_refused() { # expect_refused <work dir> <message regex>
   if run "$1"; then
      echo "   accepted, expected a refusal: see $1/log.txt  FAIL" ; FAILED=1
   elif grep -aqE "$2" "$1/log.txt"; then
      echo "   refused as expected"
   else
      echo "   failed without the expected message ($2): see $1/log.txt  FAIL" ; FAILED=1
   fi
}

echo ">> CO: radiative cooling ($TAG)"
for leg in $LEGS; do
   case $leg in
      co1)
         echo "-- co1: the box against the ODE, the CFL step, from T = 5"
         make_case co1-cfl --p 5.0 --time-max 6.5
         run_oracle "$W" exact "$W"
         echo "-- co1: the box against the ODE, the cooling limiter (dt_factor 0.01)"
         make_case co1-limited --time-max 2.0 --dt-factor 0.01
         run_oracle "$W" exact "$W" ;;
      co2)
         echo "-- co2: the floor"
         make_case co2-floor --time-max 3.0
         run_oracle "$W" exact "$W" --floor ;;
      co3)
         echo "-- co3: [cooling] with mhd-ideal"
         make_case co3-mhd --model mhd-ideal
         expect_refused "$W" "\[cooling\] needs a hydrodynamic physical_model"
         echo "-- co3: non-increasing temperatures"
         make_case co3-order --temperatures "0.1 0.1 1.0" --lambdas "0.2 1.0 0.5"
         expect_refused "$W" "\(temperatures\) must be strictly increasing"
         echo "-- co3: lambdas of another length"
         make_case co3-length --lambdas "0.2 1.0 0.5"
         expect_refused "$W" "\(lambdas\) has 3 values, \(temperatures\) 4"
         echo "-- co3: more than 64 nodes"
         nodes=$(python3 -c "print(' '.join(repr(0.1 * 1.05**k) for k in range(65)))")
         make_case co3-nodes --temperatures "$nodes" --lambdas "$nodes"
         expect_refused "$W" "\(temperatures\) has 65 nodes, at most 64 are accepted"
         echo "-- co3: no particle_mass"
         make_case co3-mass --particle-mass none
         expect_refused "$W" "failed to load \[cooling\]\.\(particle_mass\)"
         echo "-- co3: a negative dt_factor"
         make_case co3-dt --dt-factor -1.0
         expect_refused "$W" "\(dt_factor\) must be >= 0" ;;
      co4)
         echo "-- co4: the moving box, euler-de against the ODE and against euler"
         make_case co4-euler-de --model euler-de --u 10.0 --time-max 1.5 ; de="$W"
         make_case co4-euler --model euler --u 10.0 --time-max 1.5 ; eu="$W"
         if run "$de" && run "$eu"; then
            "$VENV_PY" "$CASE_DIR/co_oracle.py" exact "$de" || FAILED=1
            "$VENV_PY" "$CASE_DIR/co_oracle.py" entropy "$de" || FAILED=1
            # the floor may report round-off lifts (rho s -> p -> e is not an exact round trip), never an energy gain
            lifts=$(grep -ac "lifted to their adiabat" "$de/log.txt" || true)
            if grep -a "lifted to their adiabat" "$de/log.txt" | grep -avq "raised by +0.00000000000000000E+000"; then
               echo "   the entropy floor added energy to cooled cells: see $de/log.txt  FAIL" ; FAILED=1
            else
               echo "   the entropy floor added no energy ($lifts round-off lift reports, each raising E by 0)"
            fi
            "$VENV_PY" "$CASE_DIR/co_oracle.py" same "$de" "$eu" || FAILED=1
         else
            echo "   a run failed: see $de/log.txt, $eu/log.txt  FAIL" ; FAILED=1
         fi ;;
      co5)
         echo "-- co5: the box in the units L0 = 2, u0 = 2, rho0 = 4 through [reference]"
         make_case co5-base --time-max 1.0 --dt-factor 0.01 ; base="$W"
         make_case co5-reference --time-max 1.0 --dt-factor 0.01 ; ref="$W"
         "$VENV_PY" "$VERIF_DIR/scaling/scaling.py" dimensionalize "$base/input.ini" "$ref/input.ini" --j 1 --k 1 --m 1
         if run "$base" && run "$ref"; then
            "$VENV_PY" "$CASE_DIR/co_oracle.py" same "$base" "$ref" || FAILED=1
         else
            echo "   a run failed: see $base/log.txt, $ref/log.txt  FAIL" ; FAILED=1
         fi ;;
      co6)
         for case_ in co1-cfl co4-euler-de; do
            echo "-- co6: CPU and FNL agree, $case_"
            cpu="$CASE_DIR/work-adam_flume_cpu-$case_" ; fnl="$CASE_DIR/work-adam_flume_fnl-$case_"
            if [[ -f "$cpu/log.txt" && -f "$fnl/log.txt" ]]; then
               "$VENV_PY" "$CASE_DIR/co_oracle.py" agree "$cpu" "$fnl" || FAILED=1
            else
               echo "   missing a run: run co1 and co4 with both executables first" ; FAILED=1
            fi
         done ;;
   esac
done
if [[ $FAILED -ne 0 ]]; then echo "CO FAILED ($TAG, legs: $LEGS)" ; exit 1 ; fi
echo "CO PASSED ($TAG, legs: $LEGS)"
