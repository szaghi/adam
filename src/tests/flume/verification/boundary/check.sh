#!/usr/bin/env bash
# FLUME verification of the non-reflecting boundary (issue #121, S2): `[bc_*] type = non-reflecting`.
#
# What the kind is for, measured: in FLUME's characteristic-upwind WENO, extrapolation at an outflow and a prescribed state at an inflow
# already let linear waves out (a packet leaving at rest: 1.2e-4 of its amplitude left behind by extrapolation; through a subsonic inflow
# face, 1.4e-5 by `inflow`). What they cannot do is hold a mean pressure at a subsonic outflow: extrapolation keeps whatever pressure the
# box has. The non-reflecting kind relaxes the incoming acoustic wave toward the target pressure (Rudy & Strikwerda 1980, Poinsot &
# Lele 1992), at the price of a reflection of order sigma c / (omega L) of what passes through. The legs:
#
#   nr1  exact limits: sigma = 0 at a subsonic outflow, and any sigma at a supersonic one, equal extrapolation bit for bit.
#   nr2  the reflection of a packet leaving at rest (wavelength 0.2, width 0.1, c = 1, L = 2): below 4e-3 of its amplitude at
#        sigma = 0.28 (measured 2.86e-3) and below 1.1e-3 at sigma = 0.07 (7.56e-4), the ratio of the two within [3, 5] (linear in sigma).
#   nr3  the anchor: a uniform 2 % overpressure at Mach 0.5, both faces non-reflecting, decays below 0.15 of itself by t = 20 (measured
#        0.099); with extrapolation it stays (above 0.95).
#   nr4  a subsonic inflow face: a left-running packet leaving through it (Mach 0.5, entropy and tangential velocity from the target)
#        leaves below 6e-3 of its amplitude (measured 4.28e-3); `inflow`'s, printed, is the reference.
#
# Box [0, 2] x [0, 0.0625]^2, 256 x 8 x 8 cells (64 x 8 x 8 for nr3), 8 blocks over the ranks, periodic y and z (make_nrbc.py).
#
# Usage: ./check.sh [--np N] [--legs LIST]   LIST: nr1, nr2, nr3, nr4 or all (default all)
# FLUME_EXE overrides the executable under test (the caller owns the environment, as in the other verifications).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
NP=2
source "$VERIF_DIR/legs.sh"
legs_declare "nr1 nr2 nr3 nr4"

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np)   NP="$2" ; shift 2 ;;
      --legs) legs_select "$2" ; shift 2 ;;
      *)      echo "check.sh: unknown argument '$1' (accepted: --np N, --legs LIST)" >&2 ; exit 2 ;;
   esac
done
legs_resolve check.sh
[[ -x "$EXE" ]] || { echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2 ; exit 2 ; }
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
TAG="$(basename "$EXE")-np$NP"
ORACLE="$CASE_DIR/nrbc_oracle.py"

run() { # run <name> <make_nrbc options...>: a fresh work dir, the input, the run
   local name=$1 ; shift
   local work="$CASE_DIR/work-$TAG-$name"
   rm -rf "${work:?}"
   mkdir -p "$work"
   python3 -I "$CASE_DIR/make_nrbc.py" "$work/input.ini" "$@"
   echo ">> $name: mpirun -np $NP $(basename "$EXE")" >&2
   if ! (cd "$work" && mpirun -np "$NP" "$EXE" > log.txt 2>&1); then
      echo "check.sh: $name failed, see $work/log.txt" >&2
      exit 1
   fi
   echo "$work"
}

if has_leg nr1; then
   echo "== nr1: exact limits"
   a=$(run nr1-sub-nr --u 0.5 --sigma 0.0 --time-max 1.0 --xmin inflow)
   b=$(run nr1-sub-ex --u 0.5 --time-max 1.0 --xmin inflow --xmax extrapolation)
   "$VENV_PY" "$ORACLE" bitwise "$a" "$b"
   a=$(run nr1-sup-nr --u 1.5 --time-max 1.0 --xmin inflow)
   b=$(run nr1-sup-ex --u 1.5 --time-max 1.0 --xmin inflow --xmax extrapolation)
   "$VENV_PY" "$ORACLE" bitwise "$a" "$b"
fi
if has_leg nr2; then
   echo "== nr2: the reflection of a packet leaving at rest, linear in sigma"
   a=$(run nr2-s28 --time-max 1.5 --sigma 0.28)
   b=$(run nr2-s07 --time-max 1.5 --sigma 0.07)
   c=$(run nr2-ex --time-max 1.5 --xmin extrapolation --xmax extrapolation)
   "$VENV_PY" "$ORACLE" residual "$a" --amplitude 1e-3 --bound 4e-3
   "$VENV_PY" "$ORACLE" residual "$b" --amplitude 1e-3 --bound 1.1e-3
   "$VENV_PY" "$ORACLE" residual "$c" --amplitude 1e-3
   "$VENV_PY" "$ORACLE" ratio "$a" "$b" --amplitude 1e-3 --min 3 --max 5
fi
if has_leg nr3; then
   echo "== nr3: the pressure anchor"
   a=$(run nr3-nr --u 0.5 --ic overpressure --amplitude 0.02 --time-max 20 --nx 64)
   b=$(run nr3-ex --u 0.5 --ic overpressure --amplitude 0.02 --time-max 20 --nx 64 --xmin extrapolation --xmax extrapolation)
   "$VENV_PY" "$ORACLE" relax "$a" --deviation 0.02 --bound 0.15
   "$VENV_PY" "$ORACLE" relax "$b" --deviation 0.02 --above 0.95
fi
if has_leg nr4; then
   echo "== nr4: a packet leaving through a subsonic inflow face"
   a=$(run nr4-nr --u 0.5 --left --time-max 3.0)
   b=$(run nr4-in --u 0.5 --left --time-max 3.0 --xmin inflow)
   "$VENV_PY" "$ORACLE" residual "$a" --amplitude 1e-3 --bound 6e-3
   "$VENV_PY" "$ORACLE" residual "$b" --amplitude 1e-3
fi
echo "NR PASSED ($TAG, legs: $LEGS)"
