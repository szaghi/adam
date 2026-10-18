#!/usr/bin/env bash
# FLUME verification of self-gravity (issue #81, M8-P4): SV-7 and SV-8, the n = 1 polytrope star under isolated
# self-gravity.
#
# Why: isolated boundaries (zero padding with the lattice Green's function, SV-6) let a self-gravitating object sit in
# a box without its periodic images. The cases (make_polytrope.py) are an n = 1 polytrope (radius 1/2, rho_c = 1,
# G = 1, gamma = 2) centred in the closed box [-1, 1]^3 (inviscid walls: no mass or energy crosses it):
#   - SV-7, the equilibrium (`star`), in a light atmosphere (1e-4 rho_c, hydrostatic) with the positivity limiter,
#     which a light atmosphere needs at the star's surface: the drift of the central density and the largest Mach
#     number inside 0.8 a, the spurious motion of the discrete scheme around the continuum equilibrium, to t = 2; they
#     shrink with the grid (32^3, 64^3) by at least MIN_REDUCTION. The star's surface (the density's kink, the
#     limiter) converges at about first order and drives the pulsation: measured reductions 1.74 and 2.06, the
#     whole-star Mach number 1.50 (not asserted);
#   - SV-8, the collapse (`starc`, the pressure halved), in a 1e-2 atmosphere without the limiter (whose gravity energy
#     correction breaks the exchange where it acts), 32^3 to t = 1.5 at CFL 0.4, 0.2, 0.1: E and W exchange about 74 %
#     of |E + W|, while E + W drifts only by the time integration (its order in dt >= MIN_DRIFT_ORDER, SSP-RK3) and the
#     momentum stays at round-off (TOL_MOMENTUM).
#
# SV-14 (issue #87, P3): a smaller star on AMR (radius SV14_RADIUS, so that the 2:1 seam of the box [-1/2,1/2]^3 lies
# in the atmosphere, off the surface: a seam that cuts the surface of a star in a light atmosphere blows up, #88),
# isolated self-gravity by the multigrid's symmetric James step (every leaf on the domain boundary at the coarsest
# level):
#   - the equilibrium with the box refined once over the 32^3 base (make_polytrope.py --box: the star at the 64^3
#     spacing): the spurious motion shrinks against the same star on the uniform 32^3 grid by at least MIN_REDUCTION;
#   - the collapse on that grid with the per-stage reflux, CFL 0.4, 0.2, 0.1: E + W drifts only by the time integration
#     (order in dt >= MIN_DRIFT_ORDER), the momentum and the mass at round-off;
#   - runtime AMR: refused (D-87-4) at a regrid that refines a block on the domain boundary (--boundary-box). A
#     positive runtime case cannot change the grid in this box at a gate's cost: one level above the 32^3 base can only
#     ever refine the 8 inner blocks (0.5 wide), two levels carry the 2:1 balance to the boundary blocks (refused),
#     and 4-cell blocks are below 2 ngc.
#
# SV-15 (issue #88): the runtime-AMR collapse. Without --box the star of radius SV14_RADIUS sits on the 32^3 base and a
# regrid every 5 steps refines its surface (the Loehner marker; the 8 inner blocks, off the domain boundary): the grid
# changes at the first regrid. With `regrid_prolongation = conservative-state` the collapse reaches t = 1.5, a regrid
# changes the grid, and the mass and the momentum stay at round-off (asserted). The same run with FLUME's default
# prolongation, `conservative`, is reported: after the regrid its coldest cells (children limited variable by
# variable) keep cooling and the collapse drains one to a negative internal energy (non-finite at step 64).
#
# SV-16 (issue #91): the Jeans criterion. The star of radius SV16_RADIUS (3.5 cells per radius on the 32^3 base; the
# Jeans length of this collapsing n = 1, gamma = 2 polytrope is sqrt(2) times its radius, 5 base cells) with a Jeans
# marker (8 cells per Jeans length) next to the SV-15 Loehner marker, a regrid every 5 steps: the initial AMR refines
# the star before the collapse (more blocks than the 64 of the base at step 0, asserted), and the collapse reaches
# t = 1.5 with the default prolongation, mass and momentum at round-off (asserted). Without the Jeans marker the same
# star is refined only by the step-5 regrid, mid-collapse, and goes non-finite (reported).
#
# Usage: ./star.sh [--np N] [--legs LIST] [--keep]   LIST: sv7, sv8, sv14, sv15, sv16 or all (default all)
#
# FLUME_EXE overrides the executable under test (default exe/adam_flume_cpu_fft, mode flume-cpu-gnu-fft; FNL:
# exe/adam_flume_fnl_fft, mode flume-fnl-nvf-fft). The caller owns the matching environment (FNL: nvhpc mpirun on PATH
# and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu_fft}"
NP=2
source "$VERIF_DIR/legs.sh"
legs_declare "sv7 sv8 sv14 sv15 sv16"
KEEP=0
SV7_TIME=2.0
MIN_REDUCTION=1.5
TOL_MOMENTUM=1.0e-14
MIN_DRIFT_ORDER=2.5
TOL_MASS=1.0e-14
SV14_RADIUS=0.3
SV16_RADIUS=0.22

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np)   NP="$2" ; shift 2 ;;
      --legs) legs_select "$2" ; shift 2 ;;
      --keep) KEEP=1 ; shift ;;
      *)      echo "star.sh: unknown argument '$1' (accepted: --np N, --legs LIST, --keep)" >&2 ; exit 2 ;;
   esac
done
legs_resolve star.sh
if [[ ! -x "$EXE" ]]; then
   echo "star.sh: executable '$EXE' not found (build mode flume-cpu-gnu-fft or set FLUME_EXE)" >&2 ; exit 2
fi
EXE="$(realpath "$EXE")"

VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
if ! "$VENV_PY" -c 'import h5py, numpy' 2>/dev/null; then
   echo ">> creating the oracle venv at $REPO_ROOT/exe/.regression-venv"
   python3 -m venv "$REPO_ROOT/exe/.regression-venv"
   "$VENV_PY" -m pip install --quiet --upgrade pip
   "$VENV_PY" -m pip install --quiet h5py numpy
fi

TAG="$(basename "$EXE")-np$NP"
FAILED=0
WORKS=()
run() { # run <name> <make_polytrope arguments>: write and run one case, path in W
   local name="$1" ; shift
   W="$CASE_DIR/work-$TAG-$name"
   rm -rf "$W" ; mkdir -p "$W"
   python3 -I "$CASE_DIR/make_polytrope.py" "$W/input.ini" "$@"
   WORKS+=("$W")
   if ! (cd "$W" && mpirun -np "$NP" "$EXE" input.ini > log.txt 2>&1); then
      echo "   $name: the run failed, see $W/log.txt  FAIL"
      FAILED=1
      return 1
   fi
}
refuse() { # refuse <name> <expected message> <make_polytrope arguments>: the run must stop, naming the cause
   local name="$1" msg="$2" ; shift 2
   W="$CASE_DIR/work-$TAG-$name"
   rm -rf "$W" ; mkdir -p "$W"
   python3 -I "$CASE_DIR/make_polytrope.py" "$W/input.ini" "$@"
   WORKS+=("$W")
   if (cd "$W" && mpirun -np "$NP" "$EXE" input.ini > log.txt 2>&1); then
      echo "   refusal, $name: accepted  FAIL" ; FAILED=1
   elif grep -aq "$msg" "$W/log.txt"; then
      echo "   refusal, $name: stops the run ($msg)  PASS"
   else
      echo "   refusal, $name: the run failed without the expected message, see $W/log.txt  FAIL" ; FAILED=1
   fi
}

for leg in $LEGS; do
   case $leg in
      sv7)
         echo ">> SV-7: the equilibrium star (light atmosphere, positivity limiter)"
         ws=()
         for lev in 2 3; do
            run "sv7-l$lev" star --levels "$lev" --atmosphere 1e-4 --limiter --time-max "$SV7_TIME" --it-save 10 && \
               ws+=("$W")
         done
         if [[ ${#ws[@]} -eq 2 ]]; then
            "$VENV_PY" "$CASE_DIR/star_oracle.py" convergence "${ws[@]}" --min-reduction "$MIN_REDUCTION" || FAILED=1
         fi ;;
      sv8)
         echo ">> SV-8: the collapsing star (E + W against dt)"
         ws=()
         for cfl in 0.4 0.2 0.1; do run "sv8-collapse-$cfl" starc --cfl "$cfl" && ws+=("$W") ; done
         if [[ ${#ws[@]} -eq 3 ]]; then
            "$VENV_PY" "$CASE_DIR/jeans_oracle.py" conservation "${ws[@]}" --tol-momentum "$TOL_MOMENTUM" \
               --min-order "$MIN_DRIFT_ORDER" --label Polytrope || FAILED=1
         fi ;;
      sv14)
         echo ">> SV-14: the star on AMR, isolated self-gravity by the multigrid's James step (issue #87)"
         ws=()
         for grid in uniform box; do
            run "sv14-$grid" star --levels 2 --radius "$SV14_RADIUS" $([[ $grid == box ]] && echo --box) \
               --atmosphere 1e-4 --limiter --time-max "$SV7_TIME" --it-save 10 && ws+=("$W")
         done
         if [[ ${#ws[@]} -eq 2 ]]; then
            "$VENV_PY" "$CASE_DIR/star_oracle.py" convergence "${ws[@]}" --min-reduction "$MIN_REDUCTION" || FAILED=1
         fi
         ws=()
         for cfl in 0.4 0.2 0.1; do
            run "sv14-collapse-$cfl" starc --radius "$SV14_RADIUS" --box --reflux-cadence stage --cfl "$cfl" && \
               ws+=("$W")
         done
         if [[ ${#ws[@]} -eq 3 ]]; then
            "$VENV_PY" "$CASE_DIR/jeans_oracle.py" conservation "${ws[@]}" --tol-momentum "$TOL_MOMENTUM" \
               --min-order "$MIN_DRIFT_ORDER" --tol-mass "$TOL_MASS" --label "Polytrope AMR" || FAILED=1
         fi
         refuse sv14-boundary 'touches the domain boundary' starc --regrid 2 --boundary-box --time-max 0.2 ;;
      sv15)
         echo ">> SV-15: the runtime-AMR collapse, a regrid every 5 steps (issue #88)"
         if run sv15-state starc --radius "$SV14_RADIUS" --reflux-cadence stage --regrid 5 \
               --prolongation conservative-state; then
            n=$(grep -a "flume: regrid at" "$W/log.txt" | grep -ac " [1-9][0-9]* refined\| [1-9][0-9]* coarsened" || true)
            if [[ $n -gt 0 ]]; then
               echo "   conservative-state: reaches t = 1.5; regrids that change the grid: $n  PASS"
            else
               echo "   conservative-state: no regrid changed the grid  FAIL" ; FAILED=1
            fi
            "$VENV_PY" "$VERIF_DIR/conservation/conservation_oracle.py" --conserved "$W" --exclude rE W E_plus_W || FAILED=1
         fi
         W="$CASE_DIR/work-$TAG-sv15-conservative" ; rm -rf "$W" ; mkdir -p "$W" ; WORKS+=("$W")
         python3 -I "$CASE_DIR/make_polytrope.py" "$W/input.ini" starc --radius "$SV14_RADIUS" --reflux-cadence stage \
            --regrid 5
         if (cd "$W" && mpirun -np "$NP" "$EXE" input.ini > log.txt 2>&1); then
            echo "   conservative (the default): reaches t = 1.5 (reported)"
         else
            echo "   conservative (the default): stops, $(grep -a -m1 -o 'non-finite.*at step +[0-9]*' "$W/log.txt" || true)" \
                 "(reported, #88)"
         fi ;;
      sv16)
         echo ">> SV-16: the runtime-AMR collapse refined by the Jeans criterion (issue #91)"
         if run sv16-jeans starc --radius "$SV16_RADIUS" --reflux-cadence stage --regrid 5 --jeans 8; then
            n0=$("$VENV_PY" -c 'import sys, h5py, pathlib
print(sum(sum(1 for k in h5py.File(f, "r") if k.endswith("-origin"))
          for f in pathlib.Path(sys.argv[1]).glob("starc-000000000-proc*.h5")))' "$W")
            if [[ $n0 -gt 64 ]]; then
               echo "   jeans: the initial grid holds $n0 blocks (64 on the base), and the collapse reaches t = 1.5  PASS"
            else
               echo "   jeans: the initial grid holds $n0 blocks, the star is not refined before the collapse  FAIL"
               FAILED=1
            fi
            "$VENV_PY" "$VERIF_DIR/conservation/conservation_oracle.py" --conserved "$W" --exclude rE W E_plus_W || FAILED=1
         fi
         W="$CASE_DIR/work-$TAG-sv16-control" ; rm -rf "$W" ; mkdir -p "$W" ; WORKS+=("$W")
         python3 -I "$CASE_DIR/make_polytrope.py" "$W/input.ini" starc --radius "$SV16_RADIUS" --reflux-cadence stage \
            --regrid 5
         if (cd "$W" && mpirun -np "$NP" "$EXE" input.ini > log.txt 2>&1); then
            echo "   without the Jeans marker: reaches t = 1.5 (reported)"
         else
            echo "   without the Jeans marker: stops, $(grep -a -m1 -o 'non-finite.*at step +[0-9]*' "$W/log.txt" || true)" \
                 "(reported)"
         fi ;;
   esac
done
if [[ $KEEP -eq 0 && $FAILED -eq 0 ]]; then rm -rf "${WORKS[@]}" ; fi
if [[ $FAILED -ne 0 ]]; then
   echo "SV-7, SV-8, SV-14, SV-15, SV-16 ($LEGS) FAILED ($TAG)" ; exit 1
fi
echo "SV-7, SV-8, SV-14, SV-15, SV-16 ($LEGS) PASSED ($TAG)"
