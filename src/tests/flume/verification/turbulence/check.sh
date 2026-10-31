#!/usr/bin/env bash
# FLUME verification of the turbulence diagnostics and of the Taylor-Green initial condition (issue #119, T0.1, T0.2).
#
# Why: the ILES baseline of #119 (T0.3-T0.6) rests on the integrals of `[diagnostics] turbulence_history_save` and on the
# `taylor-green` initial condition, so both are checked first against what is known exactly.
#
#   tt1  the first history row at M 0.1 and M 1.25, n = 16, 32, 64, against the exact integrals (tgv_oracle.py t0): K within
#        1e-12 (the midpoint rule is exact on its trigonometric integrand), Z and eps at the 4th order of the gradients, eps_d
#        and pd zero (the field is solenoidal, and the same difference of x and y cancels), eps = eps_s + eps_d.
#   tt2  the same 32^3 run, 10 steps, on one rank and on N: the histories agree to 1e-11 of each column's scale (the
#        sums are in another order, so not bitwise).
#   tt3  the energy budget, Re 1600, M 0.1, to t = 2 (laminar), n = 32, 64 (--budget-ladder; 128 on the GPU: about 14 h on two
#        CPU ranks): the K drop against the integral of eps - pd. What is left is the scheme's numerical dissipation: it
#        must fall with the resolution.
#   tt4  spectra.py on the initial field: Parseval to 1e-12, the energy of the single-mode field in shell 2.
#   tt6  the pressure work: at M 0.8, n = 64, pd is large against what the scheme dissipates, and the budget to t = 2 must close
#        better with it than without it or with its sign flipped (at M 0.1 neither tt1, zero by symmetry, nor tt3 sees it).
#   tt5  (--agree TAG) this executable's tt2 history against the executable TAG (tt2 run): backends within 1e-11.
#
# Usage: ./check.sh [--np N] [--legs LIST] [--budget-ladder "32 64"] [--agree EXE-NAME]   FLUME_EXE overrides the executable.
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
NP=2
AGREE=""
BUDGET_LADDER="32 64"
source "$VERIF_DIR/legs.sh"
legs_declare "tt1 tt2 tt3 tt4 tt6"
while [[ $# -gt 0 ]]; do
   case "$1" in
      --np)    NP="$2" ; shift 2 ;;
      --legs)  legs_select "$2" ; shift 2 ;;
      --agree) AGREE="$2" ; shift 2 ;;
      --budget-ladder) BUDGET_LADDER="$2" ; shift 2 ;;
      *)       echo "check.sh: unknown argument '$1' (accepted: --np N, --legs LIST, --agree EXE-NAME)" >&2 ; exit 2 ;;
   esac
done
legs_resolve check.sh
[[ -x "$EXE" ]] || { echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2 ; exit 2 ; }
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
TAG="$(basename "$EXE")"

run() { # run <name> <np> <make_tgv options>: echo the work directory
   local name="$1" np="$2" work="$CASE_DIR/work-$TAG-$1" ; shift 2
   rm -rf "${work:?}" ; mkdir -p "$work"
   python3 -I "$CASE_DIR/make_tgv.py" "$work/input.ini" "$@"
   (cd "$work" && mpirun -np "$np" "$EXE" input.ini > log.txt 2>&1) || { echo "check.sh: $name failed, see $work/log.txt" >&2 ; exit 1 ; }
   echo "$work"
}

if has_leg tt1; then
   for mach in 0.1 1.25; do
      echo "== tt1: first history row against the exact integrals, M $mach"
      ladder=()
      for n in 16 32 64; do
         w=$(run "tt1-m$mach-n$n" "$NP" --n $n --levels 1 --mach $mach --it-max 1 --cfl 0.1)
         ladder+=("$n:$w")
      done
      "$VENV_PY" "$CASE_DIR/tgv_oracle.py" t0 "${ladder[@]}" --mach $mach
   done
fi
if has_leg tt2; then
   echo "== tt2: one rank against $NP"
   a=$(run tt2-np1 1 --n 32 --levels 1 --it-max 10)
   b=$(run tt2-npN "$NP" --n 32 --levels 1 --it-max 10)
   "$VENV_PY" "$CASE_DIR/tgv_oracle.py" agree "$a" "$b"
fi
if has_leg tt3; then
   echo "== tt3: energy budget to t = 2, Re 1600, M 0.1"
   ladder=()
   for n in $BUDGET_LADDER; do
      w=$(run "tt3-n$n" "$NP" --n $n --levels 2 --time-max 2.0)
      ladder+=("$n:$w")
   done
   "$VENV_PY" "$CASE_DIR/tgv_oracle.py" budget "${ladder[@]}" --time 2.0
fi
if has_leg tt4; then
   echo "== tt4: spectrum of the initial field"
   w=$(run tt4 "$NP" --n 32 --levels 1 --it-max 1 --it-save 1)
   "$VENV_PY" "$CASE_DIR/spectra.py" "$w" --step 0 --peak-share 0.99 --peak-shell 2
fi
if has_leg tt6; then
   echo "== tt6: pressure work, M 0.8"
   w=$(run tt6 "$NP" --n 64 --levels 2 --mach 0.8 --time-max 2.0)
   "$VENV_PY" "$CASE_DIR/tgv_oracle.py" pd "$w" --time 2.0
fi
if [[ -n "$AGREE" ]]; then
   echo "== tt5: $TAG against $AGREE"
   "$VENV_PY" "$CASE_DIR/tgv_oracle.py" agree "$CASE_DIR/work-$TAG-tt2-npN" "$CASE_DIR/work-$AGREE-tt2-npN"
fi
echo "turbulence verification PASSED (legs: $LEGS) ($TAG)"
