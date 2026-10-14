#!/usr/bin/env bash
# FLUME verification of self-gravity (issue #81, M8). SV-1 (P1): gravity from a cell potential.
#
# Why: self-gravity applies gravity through a cell potential Phi (face gravity -dPhi/dx, cell gravity the mean of the
# cell's two faces; momentum rho times the cell gravity, energy from the face mass fluxes times the face gravity).
# Before any Poisson solver exists, that path is checked on the external field it must reproduce: `[gravity] potential
# = linear` writes the constant g as Phi = -g . x and applies it through the potential path only (the constant-g
# kernels get zero). The M7 oracles then hold with their tolerances:
#   - iso2d, iso2d-mhd (one negligible step, ghosts poisoned): the hydrostatic ghosts, now filled from the potential
#     drop -(Phi_ghost - Phi_1), against the isothermal profile (atmosphere_oracle.py profile, 1e-13);
#   - freefall2d: rho u = rho g t, E = E_0 + rho g^2 t^2 / 2 (atmosphere_oracle.py freefall, 1e-12);
#   - box2d, box2d-mhd: a closed box, mass and E + rho Phi conserved (atmosphere_oracle.py energy, 1e-12; the
#     potential path is the second-order face form, whose budget closes on E + rho Phi exactly, without the wall term
#     of the fourth-order constant-g source);
#   - box2d-lim (the positivity limiter acting): the run completes, the energy drift is reported;
#   - wave2d (the M7-G2a falling entropy wave, periodic, exact solution) at 8, 16, 32 cells per block side to t = 0.2,
#     CFL scaled as GV-9: the L1 errors of rho and p, whose order is the potential path's, second order (the face form
#     without the fourth-order correction); asserted >= MIN_ORDER on the finest pair (atmosphere_oracle.py wave).
#     Measured (CPU and FNL): orders 2.12 (rho) and 2.00 (p); the constant-g path gives 5.0 and 4.6 (GV-9).
# The constant-g path itself is unchanged bitwise: the regression suite (src/tests/flume/regression) holds it.
#
# Usage: ./check.sh [--np N] [--cases "iso2d ..."] [--keep]
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
GRAVITY_DIR="$CASE_DIR/../gravity"
NP=2
MIN_ORDER=1.8
CASES="iso2d iso2d-mhd freefall2d box2d box2d-mhd box2d-lim"
KEEP=0

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np)    NP="$2" ; shift 2 ;;
      --cases) CASES="$2" ; shift 2 ;;
      --keep)  KEEP=1 ; shift ;;
      *)       echo "check.sh: unknown argument '$1' (accepted: --np N, --cases LIST, --keep)" >&2 ; exit 2 ;;
   esac
done
if [[ ! -x "$EXE" ]]; then
   echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2 ; exit 2
fi

VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
if ! "$VENV_PY" -c 'import h5py, numpy' 2>/dev/null; then
   echo ">> creating the oracle venv at $REPO_ROOT/exe/.regression-venv"
   python3 -m venv "$REPO_ROOT/exe/.regression-venv"
   "$VENV_PY" -m pip install --quiet --upgrade pip
   "$VENV_PY" -m pip install --quiet h5py numpy
fi

TAG="$(basename "$EXE")-np$NP"
FAILED=0
echo ">> SV-1: gravity from the potential Phi = -g . x ([gravity] potential = linear)"
for c in $CASES; do
   work="$CASE_DIR/work-$TAG-sv1-$c"
   rm -rf "$work" ; mkdir -p "$work"
   python3 -I "$GRAVITY_DIR/make_atmosphere.py" "$work/input.ini" "$c" --potential linear
   if ! (cd "$work" && mpirun -np "$NP" "$EXE" input.ini > log.txt 2>&1); then
      echo "   $c: the run failed, see $work/log.txt  FAIL"
      FAILED=1
      continue
   fi
   case $c in
      iso2d|iso2d-mhd)  "$VENV_PY" "$GRAVITY_DIR/atmosphere_oracle.py" profile "$work" || FAILED=1 ;;
      freefall2d)       "$VENV_PY" "$GRAVITY_DIR/atmosphere_oracle.py" freefall "$work" || FAILED=1 ;;
      box2d|box2d-mhd)  "$VENV_PY" "$GRAVITY_DIR/atmosphere_oracle.py" energy "$work" || FAILED=1 ;;
      *)                "$VENV_PY" "$GRAVITY_DIR/atmosphere_oracle.py" energy "$work" --report-total || FAILED=1 ;;
   esac
   if [[ $KEEP -eq 0 ]]; then rm -rf "$work" ; fi
done
works=()
for cells in 8 16 32; do
   work="$CASE_DIR/work-$TAG-sv1-wave-$cells"
   rm -rf "$work" ; mkdir -p "$work"
   python3 -I "$GRAVITY_DIR/make_atmosphere.py" "$work/input.ini" wave2d --cells "$cells" --time-max 0.2 \
      --cfl "$(python3 -c "print(0.4 * (8 / $cells) ** (2 / 3))")" --potential linear
   sed -i 's/^it_save = .*/it_save = 100000/' "$work/input.ini"
   if ! (cd "$work" && mpirun -np "$NP" "$EXE" input.ini > log.txt 2>&1); then
      echo "   wave2d $cells: the run failed, see $work/log.txt  FAIL"
      FAILED=1
   fi
   works+=("$work")
done
if [[ $FAILED -eq 0 ]]; then
   "$VENV_PY" "$GRAVITY_DIR/atmosphere_oracle.py" wave "${works[@]}" --min-order "$MIN_ORDER" || FAILED=1
fi
if [[ $KEEP -eq 0 ]]; then rm -rf "${works[@]}" ; fi
if [[ $FAILED -ne 0 ]]; then
   echo "SV-1 FAILED ($TAG)" ; exit 1
fi
echo "SV-1 PASSED ($TAG)"
