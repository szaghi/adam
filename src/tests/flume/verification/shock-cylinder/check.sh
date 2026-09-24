#!/usr/bin/env bash
# FLUME verification: Mach 2 shock over a cylinder, immersed boundary with init-time AMR on the solid (issue #35, P6).
#
# Why: exercises every immersed-boundary piece together: the distance function, the eikonal extrapolation into the
# solid (Jacobi sweeps) and the Euler-wall inversion, the cut spacing of the fluid cells next to the surface, the
# solid masks of the Runge-Kutta stages, and the solid AMR marker (with the coarse-fine reflux of the refined ring).
# The flow has no closed form; shock_cylinder_oracle.py asserts what the setup guarantees:
#   1. the blocks crossed by the cylinder surface are refined;
#   2. the solution is mirror-symmetric about y = 0.5 within MIRROR_TOL (catches every non-symmetric defect:
#      indexing, ordering, stencils reading the wrong ghost, races; a symmetric wrong wall law is not caught);
#   3. density and pressure are positive and finite in the fluid.
#
# Measured (P6, np 2): CPU asymmetry 1.2e-11, FNL 4.8e-12; CPU vs FNL fields 1.1e-11 (relative); 176 blocks, 96 of
# them crossed by the surface, 32 coarse-fine faces, 270 steps. The radius 0.13 makes the surface cross block faces at
# a shallow angle, so the case depends on the all-solids phi summary in the ghost cells (with the former interior-only
# summary the solution moves by 3.9e-2). The FNL run is ~4x slower than the CPU on the WSL box (per-face, per-stage
# seam copies).
#
# Usage: ./check.sh [--build] [--np N]
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12);
# --build always builds the CPU default, never the override.
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
MIRROR_TOL="1.0e-10"
NP=2
BUILD=0

while [[ $# -gt 0 ]]; do
   case "$1" in
      --build) BUILD=1 ; shift ;;
      --np)    NP="$2" ; shift 2 ;;
      *)       echo "check.sh: unknown argument '$1' (accepted: --build, --np N)" >&2 ; exit 2 ;;
   esac
done

if [[ $BUILD -eq 1 ]]; then
   (cd "$REPO_ROOT" && fobis build --mode flume-cpu-gnu)
fi
if [[ ! -x "$EXE" ]]; then
   echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2
   exit 2
fi

VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
if ! "$VENV_PY" -c 'import h5py, numpy' 2>/dev/null; then
   echo ">> creating the oracle venv at $REPO_ROOT/exe/.regression-venv"
   python3 -m venv "$REPO_ROOT/exe/.regression-venv"
   "$VENV_PY" -m pip install --quiet --upgrade pip
   "$VENV_PY" -m pip install --quiet h5py numpy
fi

work="$CASE_DIR/work-$(basename "$EXE")-np$NP"
rm -rf "$work"
mkdir -p "$work"
cp "$CASE_DIR/shock-cylinder.ini" "$work/"
echo ">> shock-cylinder: mpirun -np $NP $(basename "$EXE")"
if ! (cd "$work" && mpirun -np "$NP" "$EXE" shock-cylinder.ini > log.txt 2>&1); then
   echo "check.sh: shock-cylinder run failed, see $work/log.txt" >&2
   exit 1
fi
"$VENV_PY" "$CASE_DIR/shock_cylinder_oracle.py" "$work" --mirror-tol "$MIRROR_TOL"
echo "shock-cylinder PASSED ($(basename "$work"))"
