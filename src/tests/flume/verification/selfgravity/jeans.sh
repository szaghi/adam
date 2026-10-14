#!/usr/bin/env bash
# FLUME verification of self-gravity (issue #81, M8-P3): SV-3, SV-4 and SV-5, the Jeans modes.
#
# Why: self-gravity couples the FFT Poisson solve (one per Runge-Kutta stage) to the potential path of SV-1. A linear
# Jeans mode obeys w^2 = c^2 k^2 - 4 pi G rho_0, which only a correct coupling reproduces: the factor 4 pi G, the sign
# of the force, the mean removed, the per-stage solve. The cases (make_jeans.py) are a uniform gas (rho_0 = 1, c = 1,
# gamma = 5/3, periodic) with the mode of k = 2 pi (1, 1):
#   - SV-3, the stable wave (lambda = lambda_J / 2, w = pi sqrt(6)) to one period: the L1 error of rho against the
#     exact wave over the amplitude, second order (>= MIN_ORDER on the finest pair), in 2-D (32^2, 64^2, 128^2) and in
#     3-D (16^3, 32^3, the wave in the x-y plane);
#   - SV-4, the unstable mode (lambda = 1.5 lambda_J, sigma = pi sqrt(10)): the growth rate of the linear mode (64^2,
#     amplitude 1e-6, t <= 0.5) against the exact one (TOL_GROWTH, relative), its mirror symmetry (the point reflection
#     about the density's centre of symmetry, an exact symmetry of the equations: TOL_SYMMETRY), and, into the
#     nonlinear collapse (amplitude 1e-2 to t = 0.4, CFL 0.4, 0.2, 0.1), the momentum (round-off: TOL_MOMENTUM) and the
#     drift of E + W, which only the time integration leaves: its order in dt >= MIN_DRIFT_ORDER (SSP-RK3);
#   - SV-5, scaling covariance: the stable 2-D case written in the dimensional units L0 = 2^2, u0 = 2^-1, rho0 = 4^-2
#     with a [reference] section (G converted with u0^2 / (rho0 L0^2)) equals the code-unit run bit for bit
#     (scaling/reference.sh, NV-5).
#
# Usage: ./jeans.sh [--np N] [--legs "sv3-2d sv3-3d sv4 sv5"] [--keep]
#
# FLUME_EXE overrides the executable under test (default exe/adam_flume_cpu_fft, mode flume-cpu-gnu-fft; FNL:
# exe/adam_flume_fnl_fft, mode flume-fnl-nvf-fft). The caller owns the matching environment (FNL: nvhpc mpirun on PATH
# and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu_fft}"
NP=2
LEGS="sv3-2d sv3-3d sv4 sv5"
KEEP=0
MIN_ORDER=1.8
TOL_GROWTH=1.0e-2
TOL_SYMMETRY=1.0e-9
TOL_MOMENTUM=1.0e-14
MIN_DRIFT_ORDER=2.5

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np)   NP="$2" ; shift 2 ;;
      --legs) LEGS="$2" ; shift 2 ;;
      --keep) KEEP=1 ; shift ;;
      *)      echo "jeans.sh: unknown argument '$1' (accepted: --np N, --legs LIST, --keep)" >&2 ; exit 2 ;;
   esac
done
if [[ ! -x "$EXE" ]]; then
   echo "jeans.sh: executable '$EXE' not found (build mode flume-cpu-gnu-fft or set FLUME_EXE)" >&2 ; exit 2
fi
EXE="$(realpath "$EXE")"

VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
if ! "$VENV_PY" -c 'import h5py, numpy' 2>/dev/null; then
   echo ">> creating the oracle venv at $REPO_ROOT/exe/.regression-venv"
   python3 -m venv "$REPO_ROOT/exe/.regression-venv"
   "$VENV_PY" -m pip install --quiet --upgrade pip
   "$VENV_PY" -m pip install --quiet h5py numpy
fi
ORACLE="$CASE_DIR/jeans_oracle.py"

TAG="$(basename "$EXE")-np$NP"
FAILED=0
WORKS=()
case_dir() { # case_dir <name> <make_jeans arguments>: a fresh work-$TAG-<name> with its input, path in W
   local name="$1" ; shift
   W="$CASE_DIR/work-$TAG-$name"
   rm -rf "$W" ; mkdir -p "$W"
   python3 -I "$CASE_DIR/make_jeans.py" "$W/input.ini" "$@"
   WORKS+=("$W")
}
run() { # run <name> <make_jeans arguments>: write and run one case
   local name="$1"
   case_dir "$@"
   if ! (cd "$W" && mpirun -np "$NP" "$EXE" input.ini > log.txt 2>&1); then
      echo "   $name: the run failed, see $W/log.txt  FAIL"
      FAILED=1
      return 1
   fi
}

for leg in $LEGS; do
   case $leg in
      sv3-2d)
         echo ">> SV-3: the stable Jeans wave, 2-D"
         ws=()
         for cells in 8 16 32; do run "sv3-2d-$cells" jeans2d --cells "$cells" && ws+=("$W") ; done
         if [[ ${#ws[@]} -eq 3 ]]; then "$VENV_PY" "$ORACLE" wave "${ws[@]}" --min-order "$MIN_ORDER" || FAILED=1 ; fi ;;
      sv3-3d)
         echo ">> SV-3: the stable Jeans wave, 3-D"
         ws=()
         for cells in 8 16; do run "sv3-3d-$cells" jeans3d --cells "$cells" && ws+=("$W") ; done
         if [[ ${#ws[@]} -eq 2 ]]; then "$VENV_PY" "$ORACLE" wave "${ws[@]}" --min-order "$MIN_ORDER" || FAILED=1 ; fi ;;
      sv4)
         echo ">> SV-4: the unstable Jeans mode"
         if run sv4-linear jeansu2d --cells 16 --it-save 5; then
            "$VENV_PY" "$ORACLE" growth "$W" --tol "$TOL_GROWTH" || FAILED=1
            "$VENV_PY" "$ORACLE" symmetry "$W" --tol "$TOL_SYMMETRY" || FAILED=1
         fi
         ws=()
         for cfl in 0.4 0.2 0.1; do run "sv4-collapse-$cfl" jeansc2d --cells 16 --cfl "$cfl" && ws+=("$W") ; done
         if [[ ${#ws[@]} -eq 3 ]]; then
            "$VENV_PY" "$ORACLE" conservation "${ws[@]}" --tol-momentum "$TOL_MOMENTUM" \
               --min-order "$MIN_DRIFT_ORDER" || FAILED=1
         fi ;;
      sv5)
         echo ">> SV-5: scaling covariance with G"
         case_dir sv5-input jeans2d --cells 8 --time-max 0.2
         FLUME_EXE="$EXE" "$CASE_DIR/../scaling/reference.sh" --np "$NP" --inputs "$W/input.ini" || FAILED=1 ;;
      *)
         echo "jeans.sh: unknown leg '$leg' (accepted: sv3-2d, sv3-3d, sv4, sv5)" >&2 ; exit 2 ;;
   esac
done
if [[ $KEEP -eq 0 && $FAILED -eq 0 ]]; then rm -rf "${WORKS[@]}" ; fi
if [[ $FAILED -ne 0 ]]; then
   echo "SV-3..5 FAILED ($TAG)" ; exit 1
fi
echo "SV-3..5 PASSED ($TAG)"
