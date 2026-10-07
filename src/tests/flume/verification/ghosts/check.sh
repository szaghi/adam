#!/usr/bin/env bash
# FLUME verification GP (issue #65, P0): every face and edge ghost cell holds the value a stencil must read.
#
# Why: the dissipative terms of M4 take cross derivatives (the tangential derivatives of the viscous stress and of the
# current at a face), which read the edge ghosts: the cells outside a block along two axes. The directional WENO
# stencils never read them, so until P0 nothing checked what they held. Each case holds a linear field (exact under
# the same-level copy, the 2:1 restriction and the tricubic coarse->fine fill) and takes one negligible step
# (CFL = 1e-30), after which the fields are written with their ghosts filled in the order of a Runge-Kutta stage;
# ghost_probe.py checks every face and edge ghost against the value it must hold (the boundary condition of a face
# composed with the exchange). Cases (make_probe.py): box3d (octree, 2:1 seams meeting inflow and two walls; Euler
# and MHD), channel2d (quadtree, periodic x with a wall), mirror3d and refined3d (two realms, 1:1 and 2:1 inter-realm
# seams), walls3d (the no-slip and isothermal walls of issue #65 P1, moving and resting; Euler and MHD with EGLM),
# step (the Woodward-Colella three-realm forest of ../step at N = 80, refined boxes on), each on 1, 2 and 3 ranks.
#
# Every case sets [diagnostics] ghost_poison: each ghost is NaN before the fill, so a ghost read before its donor is
# written shows as NaN (the negligible step makes a stale ghost equal to a fresh one otherwise).
#
# Measured (CPU np 1-3, FNL np 1-2): every face and edge group within 1.3e-15 of the field scale. Before P0 the realm
# edges (rows beyond two realm faces: wall + wall, wall + seam, inflow + wall, ...) held an inward diagonal copy, off by
# up to 0.78; before P1 the extrapolation chains raced on FNL (block+extrapolation edges read NaN with the poison).
#
# Usage: ./check.sh [--np "1 2 3"] [--cases "box3d box3d-mhd channel2d mirror3d refined3d walls3d walls3d-mhd step"]
#        [--keep]
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
NPS="1 2 3"
CASES="box3d box3d-mhd channel2d mirror3d refined3d walls3d walls3d-mhd step"
KEEP=0

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np)    NPS="$2" ; shift 2 ;;
      --cases) CASES="$2" ; shift 2 ;;
      --keep)  KEEP=1 ; shift ;;
      *)       echo "check.sh: unknown argument '$1' (accepted: --np \"N ...\", --cases \"...\", --keep)" >&2 ; exit 2 ;;
   esac
done
if [[ ! -x "$EXE" ]]; then
   echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2
   exit 2
fi
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
TAG="$(basename "$EXE")"
FAILED=0

for np in $NPS; do
   for case in $CASES; do
      w="$CASE_DIR/work-$TAG-np$np-$case"
      rm -rf "$w" ; mkdir -p "$w"
      input=probe.ini
      case "$case" in
         step)      "$VENV_PY" "$CASE_DIR/../step/make_step.py" "$w" --cells 80 --refine > /dev/null
                    "$VENV_PY" "$CASE_DIR/make_probe.py" "$w" --linearize
                    input=step.ini ;;
         box3d-mhd)   "$VENV_PY" "$CASE_DIR/make_probe.py" "$w" box3d --model mhd ;;
         walls3d-mhd) "$VENV_PY" "$CASE_DIR/make_probe.py" "$w" walls3d --model mhd ;;
         *)         "$VENV_PY" "$CASE_DIR/make_probe.py" "$w" "$case" ;;
      esac
      echo ">> GP $case, np $np ($TAG)"
      if ! (cd "$w" && mpirun -np "$np" "$EXE" "$input" > log.txt 2>&1); then
         echo "check.sh: $case np $np failed, see $w/log.txt" >&2
         FAILED=1
         continue
      fi
      "$VENV_PY" "$CASE_DIR/ghost_probe.py" "$w" | grep -v ' info$' || FAILED=1
      [[ $KEEP -eq 1 ]] || find "$w" -name '*.h5' -delete
   done
done

if [[ $FAILED -eq 0 ]]; then
   echo "GP PASSED ($TAG)"
else
   echo "GP FAILED ($TAG)"
   exit 1
fi
