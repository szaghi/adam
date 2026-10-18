#!/usr/bin/env bash
# FLUME verification DE: dual energy, the advected entropy of `physical_model = euler-de` (issue #91, P2).
#
# Why: in a cold supersonic flow the internal energy e = E - |rho u|^2 / (2 rho) is a small difference of large numbers,
# and the truncation error of the total-energy scheme drains it (the #88 collapse). euler-de carries rho s, s = p/rho^g,
# and at the end of every step lifts a cold cell (e < eta E) that fell below its adiabat back to it; everywhere else the
# entropy follows the energy. The legs check that the model changes nothing where it must not and keeps shocks heated:
#   de1  Sod (the sod-x input): no cell is cold, so euler-de gives the Euler fields of euler bitwise;
#   de2  the cold planar Noh shock (make_de.py noh, Mach 845, t = 0.6): every cell starts cold, the floor must not
#        remove the shock's heat: the post-shock pressure of both models within 2 % of the exact 1.2, the cold streams ahead at
#        p = 1e-6 within 1e-3, the two densities within 2e-2;
#   de3  the contract: [physics].(dual_energy_eta) outside (0, 1) is refused; on FNL euler-de itself is refused (the
#        FNL kernels are P3).
# The runtime-AMR collapse of #88 with euler-de is selfgravity/star.sh SV-17.
#
# Usage: ./check.sh [--legs LIST]   (legs de1 de2 de3, all; default all)
#
# FLUME_EXE overrides the executable under test (default exe/adam_flume_cpu); the caller owns the matching environment.
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
source "$VERIF_DIR/legs.sh"
legs_declare "de1 de2 de3"
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
FNL=0 ; [[ $TAG == *fnl* ]] && FNL=1
FAILED=0

run() { # run <work dir>: run the input of the work dir; 0 on success
   (cd "$1" && mpirun -np 2 "$EXE" input.ini < /dev/null > log.txt 2>&1)
}
make_run() { # make_run <leg> <case> <model> [make_de.py options]: write and run one case, its dir in W
   local leg=$1 case_=$2 model=$3 ; shift 3
   W="$CASE_DIR/work-$TAG-$leg-$model"
   rm -rf "$W" ; mkdir -p "$W"
   python3 -I "$CASE_DIR/make_de.py" "$W/input.ini" "$case_" --model "$model" "$@"
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

echo ">> DE: dual energy ($TAG)"
for leg in $LEGS; do
   if [[ $FNL -eq 1 && $leg != de3 ]]; then echo "-- $leg: CPU only until #91 P3 (skipped)" ; continue ; fi
   case $leg in
      de1|de2)
         case_=sod ; [[ $leg == de2 ]] && case_=noh
         echo "-- $leg: $case_, euler against euler-de"
         dirs=()
         for model in euler euler-de; do
            make_run "$leg" "$case_" "$model"
            if run "$W"; then dirs+=("$W") ; else echo "   $model: the run failed, see $W/log.txt  FAIL" ; FAILED=1 ; fi
         done
         if [[ ${#dirs[@]} -eq 2 ]]; then "$VENV_PY" "$CASE_DIR/de_oracle.py" "$case_" "${dirs[@]}" || FAILED=1 ; fi ;;
      de3)
         if [[ $FNL -eq 1 ]]; then
            echo "-- de3: euler-de on FNL"
            make_run de3 sod euler-de
            expect_refused "$W" "euler-de runs on the CPU backend only"
         else
            for eta in 0.0 1.0; do
               echo "-- de3: dual_energy_eta = $eta"
               make_run "de3-$eta" sod euler-de --eta "$eta"
               expect_refused "$W" "dual_energy_eta\) must lie in \(0, 1\)"
            done
         fi ;;
   esac
done
if [[ $FAILED -ne 0 ]]; then echo "DE FAILED ($TAG, legs: $LEGS)" ; exit 1 ; fi
echo "DE PASSED ($TAG, legs: $LEGS)"
