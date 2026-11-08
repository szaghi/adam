#!/usr/bin/env bash
# FLUME verification of the open boundaries of issue #121: the non-reflecting boundary (S2, `[bc_*] type = non-reflecting`, legs nr*)
# and the oblique-shock generator (S3, `[bc_y_*] type = shock-generator`, legs sg*).
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
# The shock generator: a Mach 2 free stream in [0, 2] x [0, 1] x [0, 0.25], 64 x 32 x 8 cells, inflow at x_min, extrapolation at x_max and
# y_min, periodic z (make_shockgen.py); the shock of deflection 10 degrees (beta 39.314, p2/p1 1.70658) enters through y = 1 at x = 0.5.
#   sg1  the uniform region behind the shock (16 cells off it, 4 below the streamline from its foot, above y = 0.2) against
#        Rankine-Hugoniot at t = 2.5 (steady): p and rho within 3e-3 (measured 1.47e-3, 1.22e-3), the flow angle within 0.05 degrees
#        (0.026); the shock angle fitted from the mid-pressure crossing of the rows within 0.3 degrees (39.322).
#   sg2  the same shock given by its angle (shock_angle 39.313932): the same bounds (the same digits).
#   sg3  a shock that enters beyond the box (shock_x 3) equals an extrapolated top bit for bit.
#   sg4  the regular reflection off a slip wall (S4): y_min `wall-inviscid`, 96 x 48 x 8 cells, the shock entering at x = 0.1 (it meets the
#        wall at x = 1.321), t = 5 (steady; at t = 2.5 the region behind the reflected shock is still 6 % low): behind the reflected shock
#        p and rho within 4e-3 of the exact p3 / p1 = 2.80319 (measured 2.1e-3, 1.5e-3), the flow parallel to the wall within 0.1 degrees
#        (0.062); the reflected shock's angle to the wall within 0.3 degrees of 39.384 (39.437).
#
# Usage: ./check.sh [--np N] [--legs LIST]   LIST: nr1, nr2, nr3, nr4, sg1, sg2, sg3, sg4 or all (default all)
# FLUME_EXE overrides the executable under test (the caller owns the environment, as in the other verifications).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
NP=2
source "$VERIF_DIR/legs.sh"
legs_declare "nr1 nr2 nr3 nr4 sg1 sg2 sg3 sg4"

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
SG_ORACLE="$CASE_DIR/shockgen_oracle.py"
SG_BOUND=3e-3
SG_ANGLE=0.05

run() { # run <name> <make_nrbc options...>: a fresh work dir, the input, the run; MAKE overrides the generator
   local name=$1 ; shift
   local work="$CASE_DIR/work-$TAG-$name"
   rm -rf "${work:?}"
   mkdir -p "$work"
   python3 -I "$CASE_DIR/${MAKE:-make_nrbc.py}" "$work/input.ini" "$@"
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
SG=(--mach 2 --shock-x 0.5)
if has_leg sg1; then
   echo "== sg1: the post-shock state and the shock angle"
   a=$(MAKE=make_shockgen.py run sg1-def --nx 64 --deflection 10 --time-max 2.5)
   "$VENV_PY" "$SG_ORACLE" post "$a" "${SG[@]}" --deflection 10 --bound "$SG_BOUND" --angle-bound "$SG_ANGLE"
   "$VENV_PY" "$SG_ORACLE" angle "$a" "${SG[@]}" --deflection 10 --bound 0.3
fi
if has_leg sg2; then
   echo "== sg2: the shock given by its angle"
   a=$(MAKE=make_shockgen.py run sg2-ang --nx 64 --angle 39.313932 --time-max 2.5)
   "$VENV_PY" "$SG_ORACLE" post "$a" "${SG[@]}" --angle 39.313932 --bound "$SG_BOUND" --angle-bound "$SG_ANGLE"
fi
if has_leg sg3; then
   echo "== sg3: a shock beyond the box is an extrapolated top"
   a=$(MAKE=make_shockgen.py run sg3-far --nx 64 --deflection 10 --shock-x 3.0 --time-max 0.5)
   b=$(MAKE=make_shockgen.py run sg3-ex --nx 64 --top extrapolation --time-max 0.5)
   "$VENV_PY" "$SG_ORACLE" bitwise "$a" "$b"
fi
if has_leg sg4; then
   echo "== sg4: the regular reflection off a slip wall"
   a=$(MAKE=make_shockgen.py run sg4-refl --nx 96 --deflection 10 --shock-x 0.1 --bottom wall-inviscid --time-max 5.0)
   "$VENV_PY" "$SG_ORACLE" reflect "$a" --mach 2 --deflection 10 --shock-x 0.1 --bound 4e-3 --angle-bound 0.1
fi
echo "NR/SG PASSED ($TAG, legs: $LEGS)"
