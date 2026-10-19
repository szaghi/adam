#!/usr/bin/env bash
# FLUME verification SI: sink particles, `[sinks]` (issue #94).
#
# Why: sinks are a replicated list of point masses coupled to the gas through the self-gravity solve only (their mass
# CIC-deposited into the Poisson source, the mesh gravity CIC-interpolated back, D-M12-2), advanced by kick-drift-kick
# once per step. make_si.py writes the box ([-1, 1]^3, 32^3, gamma = 2, G = 1), si_oracle.py checks it. The legs:
#   si1  a sink at rest in a uniform periodic background (rho = 1e-3, the sink off every symmetry plane, t = 0.5): it
#        stays put within 1e-3 cells, and gas plus sink momentum is conserved within 1e-3 of the momentum they exchange
#        (the gas feels face gradients of the potential, the sink centred ones: equal and opposite to truncation);
#   si2  the equal-mass circular binary (masses 1, separation 0.5 = 8 cells, isolated, near-vacuum background, one
#        orbit): the separation within 3 %, the period against Kepler's at the mean separation within 2 %, the
#        orbital energy and the angular momentum within 1 %, the centre of mass within 1e-10 (the mesh force is
#        softened at about two cells: not round-off);
#   si3  the contract: [sinks] without self-gravity, with mhd-ideal, lists of different lengths, a non-positive mass
#        and a sink outside the domain are refused;
#   si4  restart: the si1 box 20 steps in one run, and 10 + 10 across a restart (<restart_basename>.sinks): the
#        conservation and sinks histories bitwise equal;
#   si5  (P2) accretion: a sink of mass 0.1 at rest in a dense cold periodic background (rho = 1, p = 1e-3, rho_thr
#        about 0.06), t = 0.3: gas plus sink mass conserved within 1e-13, the sink grows, gas plus sink momentum within
#        5 % of the momentum exchanged (the gas feels face gradients, the sink the CIC gather: near a point mass they
#        differ by about 1 %);
#   si6  (P2) creation: the selfgravity-star polytrope at a tenth of its pressure, kept isothermal at T = 0.01 by strong
#        cooling (M11), closed walls, to t = 0.6: exactly one sink, within a cell of the centre, mass conserved within
#        1e-13 (adiabatic gamma = 2 never crosses the threshold: rho_thr grows with rho);
#   si7  (P2) no creation: a static periodic gamma = 1.4 box at 1.2 rho_thr, 20 steps: above the threshold, every cell a
#        potential minimum, but Jeans-stable (|E_grav| / 2 E_th = 0.66): no sink;
#   si8  (P2) merger: two sinks of mass 1, 1.28 cells apart, approaching: one sink of mass 2, its momentum the sum
#        (zero) within 1e-10;
#   si9  (P3) CPU and FNL agree on the si2, si5, si6, si8, si10 and si11 runs (run them with both executables first):
#        the same steps and sinks, masses within 1e-8 of the largest, positions within 1e-8 cells, velocities within
#        1e-8 of the largest speed; si6 and the si10 collapse (a sink formed at a symmetric centre: which of the 8 tied
#        cells wins is round-off) and si11 (the sink at the centre: its velocity is round-off, 1e-11) the masses and the
#        final distance from the centre (1e-2 cells);
#   si10 (P4) runtime AMR (a Jeans marker, a regrid every 4 steps): the si1 sink on a 4^3-block base of 8 cells with one
#        level above it, 40 steps: the blocks within r_acc + 2 cells of the sink on the finest level, no other block
#        refined (the Jeans marker never fires in that background), the sink within 1e-3 cells, the mass conserved
#        within 1e-13 across the regrids; the si6 collapse at 16^3 with two levels above the base: one sink, formed on
#        the finest level (the only one it may form on with runtime regridding), at the centre, its neighbourhood on
#        that level, the mass within 1e-13;
#   si11 (P4) the overdense singular isothermal sphere of Federrath et al. (2010, section 3.4; A = 29.3, gamma = 1.001,
#        a seed sink at the centre, t = 0.25): the accretion rate over the last 60 % of the run within 6 % of Shu's
#        (1977) m_0(A) a^3 / G (4.1 % at 32^3; F10 about 3 %), the mass within 1e-13.
#
# Usage: ./check.sh [--legs LIST]   (legs si1 to si11, all; default every leg but si9)
#
# FLUME_EXE overrides the executable under test (default exe/adam_flume_cpu_fft: self-gravity needs the FFT build); the
# caller owns the matching environment.
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu_fft}"
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
source "$VERIF_DIR/legs.sh"
legs_declare "si1 si2 si3 si4 si5 si6 si7 si8 si9 si10 si11" "si1 si2 si3 si4 si5 si6 si7 si8 si10 si11"
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
make_case() { # make_case <name> [make_si.py options]: write one case, its dir in W
   local name=$1 ; shift
   W="$CASE_DIR/work-$TAG-$name"
   rm -rf "$W" ; mkdir -p "$W"
   python3 -I "$CASE_DIR/make_si.py" "$W/input.ini" "$@"
}
run_oracle() { # run_oracle <work dir> <si_oracle.py arguments>: run the case, then the oracle
   local w=$1 ; shift
   if run "$w"; then
      "$VENV_PY" "$CASE_DIR/si_oracle.py" "$@" || FAILED=1
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

echo ">> SI: sink particles ($TAG)"
for leg in $LEGS; do
   case $leg in
      si1)
         echo "-- si1: a sink at rest in a uniform periodic background"
         make_case si1-rest --case rest
         run_oracle "$W" rest "$W" ;;
      si2)
         echo "-- si2: the circular binary, one orbit"
         make_case si2-binary --case binary --time-max 1.5708
         run_oracle "$W" binary "$W" ;;
      si3)
         echo "-- si3: [sinks] without self-gravity"
         make_case si3-gravity --case rest --no-gravity
         expect_refused "$W" "\[sinks\] needs \[gravity\]\.\(self_gravity\)"
         echo "-- si3: [sinks] with mhd-ideal"
         make_case si3-mhd --case rest --model mhd-ideal
         expect_refused "$W" "\[sinks\] needs a hydrodynamic physical_model"
         echo "-- si3: lists of different lengths"
         make_case si3-length --case binary --x "0.1"
         expect_refused "$W" "\[sinks\]\.\(x\) has 1 values, \(mass\) 2"
         echo "-- si3: a non-positive mass"
         make_case si3-mass --case rest --mass "-1.0"
         expect_refused "$W" "\[sinks\]\.\(mass\) must all be > 0"
         echo "-- si3: a sink outside the domain"
         make_case si3-outside --case rest --x "1.5"
         expect_refused "$W" "lies outside the domain" ;;
      si4)
         echo "-- si4: restart, 20 steps against 10 + 10"
         make_case si4-A --case rest --it-max 20 --time-max 100.0 ; a="$W"
         make_case si4-B --case rest --it-max 10 --time-max 100.0 --restart-save 10 ; b="$W"
         if run "$a" && run "$b"; then
            python3 -I "$CASE_DIR/make_si.py" "$b/input.ini" --case rest --it-max 20 --time-max 100.0 --restart
            mv "$b/log.txt" "$b/log-first.txt"
            if run "$b"; then
               "$VENV_PY" "$CASE_DIR/si_oracle.py" same "$a" "$b" || FAILED=1
            else
               echo "   the restarted run failed: see $b/log.txt  FAIL" ; FAILED=1
            fi
         else
            echo "   a run failed: see $a/log.txt, $b/log.txt  FAIL" ; FAILED=1
         fi ;;
      si5)
         echo "-- si5: accretion from a dense cold background"
         make_case si5-accrete --case accrete
         run_oracle "$W" accrete "$W" ;;
      si6)
         echo "-- si6: the isothermal collapse forms one central sink"
         make_case si6-collapse --case collapse
         run_oracle "$W" collapse "$W" ;;
      si7)
         echo "-- si7: a Jeans-stable box above the threshold forms none"
         make_case si7-stable --case stable --it-max 20
         run_oracle "$W" none "$W" ;;
      si8)
         echo "-- si8: two bound, converging sinks merge"
         make_case si8-merge --case merge --it-max 5
         run_oracle "$W" merge "$W" ;;
      si9)
         for case_ in si2-binary si5-accrete si6-collapse si8-merge si10-rest-amr si10-collapse-amr si11-shu; do
            echo "-- si9: CPU and FNL agree, $case_"
            cpu="$CASE_DIR/work-adam_flume_cpu_fft-$case_" ; fnl="$CASE_DIR/work-adam_flume_fnl_fft-$case_"
            if [[ -f "$cpu/log.txt" && -f "$fnl/log.txt" ]]; then
               # a sink at a symmetric centre: a tie whose winner is round-off (si6, si10), or a velocity that is
               # round-off (si11, 1e-11): the masses and the final distance from the centre
               only="" ; [[ $case_ =~ ^(si6-collapse|si10-collapse-amr|si11-shu)$ ]] && only="--mass-only"
               "$VENV_PY" "$CASE_DIR/si_oracle.py" agree "$cpu" "$fnl" $only || FAILED=1
            else
               echo "   missing a run: run $case_ with both executables first" ; FAILED=1
            fi
         done ;;
      si10)
         echo "-- si10: a sink at rest keeps its neighbourhood on the finest level across the regrids"
         make_case si10-rest-amr --case rest --cells 32 --base-levels 2 --amr 1 --it-max 40 --time-max 100.0
         run_oracle "$W" amr "$W" --exclusive --rest
         echo "-- si10: the collapse on runtime AMR forms one sink on the finest level"
         make_case si10-collapse-amr --case collapse --cells 16 --amr 2
         run_oracle "$W" amr "$W"
         "$VENV_PY" "$CASE_DIR/si_oracle.py" collapse "$W" || FAILED=1 ;;
      si11)
         echo "-- si11: the overdense singular isothermal sphere accretes at Shu's rate"
         make_case si11-shu --case shu
         run_oracle "$W" shu "$W" ;;
   esac
done
if [[ $FAILED -ne 0 ]]; then echo "SI FAILED ($TAG, legs: $LEGS)" ; exit 1 ; fi
echo "SI PASSED ($TAG, legs: $LEGS)"
