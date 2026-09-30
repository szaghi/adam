#!/usr/bin/env bash
# FLUME verification RV-2, RV-3: the weno-riemann scheme on the V2 vortex and the V1 Sod problems (issue #47).
#
# Why: `scheme_space = weno-riemann` replaces the WENO flux splitting by the WENO interpolation of the face states, a
# Riemann flux and a high-order correction (Chen, Toth & Gombosi 2016). Its oracles are those of the splitting scheme,
# on the same inputs with the [numerics] block rewritten:
#   RV-2 vortex  the isentropic vortex ladder 64/128/256 (vortex/): L1 order of the finest pair >= ORDER_MIN;
#   RV-3 sod     Sod along x, y, z (sod/): L1(rho) against the exact solution <= L1_MAX, y and z BITWISE equal to x
#                after the permutation, every transverse copy bitwise identical.
#
# Usage: ./check.sh [--np N] [--leg vortex|sod|all] [--solver S] [--correction C] [--sensor S] [--recon R]
#                   [--order-min X] [--l1-max X]
# Defaults: all legs, llf, 6th, weno, characteristic; ORDER_MIN 4.5 (as V2), L1_MAX 3.31e-03 (as V1).
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VER_DIR="$(cd "$CASE_DIR/.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
NP=2
LEG=all
SOLVER=llf
CORRECTION=6th
SENSOR=weno
RECON=characteristic
ORDER_MIN="4.5"
L1_MAX="3.31e-03"

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np)         NP="$2" ; shift 2 ;;
      --leg)        LEG="$2" ; shift 2 ;;
      --solver)     SOLVER="$2" ; shift 2 ;;
      --correction) CORRECTION="$2" ; shift 2 ;;
      --sensor)     SENSOR="$2" ; shift 2 ;;
      --recon)      RECON="$2" ; shift 2 ;;
      --order-min)  ORDER_MIN="$2" ; shift 2 ;;
      --l1-max)     L1_MAX="$2" ; shift 2 ;;
      *) echo "check.sh: unknown argument '$1'" >&2 ; exit 2 ;;
   esac
done
if [[ ! -x "$EXE" ]]; then
   echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2
   exit 2
fi
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
TAG="$(basename "$EXE")-np$NP-$SOLVER-$CORRECTION-$SENSOR-$RECON"

set_numerics() { # rewrite the [numerics] block of an input for weno-riemann
   "$VENV_PY" - "$1" "$SOLVER" "$CORRECTION" "$SENSOR" "$RECON" <<'EOF'
import re, sys
path, solver, corr, sensor, recon = sys.argv[1:]
text = open(path).read()
new = (f"scheme_space             = weno-riemann\nreconstruction_variables = {recon}\nriemann_solver           = {solver}\n"
       f"flux_correction          = {corr}\nflux_correction_sensor   = {sensor}\n")
text, n = re.subn(r"^scheme_space\s*=.*\n^reconstruction_variables\s*=.*\n", new, text, flags=re.M)
if n != 1:
    sys.exit(f"{path}: [numerics] block not found")
open(path, "w").write(text)
EOF
}

run() { # case-dir input work
   rm -rf "$3"
   mkdir -p "$3"
   cp "$1/$2" "$3/"
   set_numerics "$3/$2"
   echo ">> $2: mpirun -np $NP $(basename "$EXE") ($SOLVER, $CORRECTION, $SENSOR, $RECON)"
   local start
   start=$(date +%s)
   if ! (cd "$3" && mpirun -np "$NP" "$EXE" "$2" > log.txt 2>&1); then
      echo "check.sh: $2 run failed, see $3/log.txt" >&2
      exit 1
   fi
   echo "   done in $(( $(date +%s) - start )) s"
}

if [[ $LEG == vortex || $LEG == all ]]; then
   WORK=()
   for n in 064 128 256; do
      run "$VER_DIR/vortex" "vortex-n$n.ini" "$CASE_DIR/work-$TAG-vortex-n$n"
      WORK+=("$CASE_DIR/work-$TAG-vortex-n$n")
   done
   "$VENV_PY" "$VER_DIR/vortex/vortex_oracle.py" "${WORK[@]}" --order-min "$ORDER_MIN"
   echo "RV-2 PASSED ($TAG)"
fi
if [[ $LEG == sod || $LEG == all ]]; then
   WORK=()
   for d in x y z; do
      run "$VER_DIR/sod" "sod-$d.ini" "$CASE_DIR/work-$TAG-sod-$d"
      WORK+=("$CASE_DIR/work-$TAG-sod-$d")
   done
   "$VENV_PY" "$VER_DIR/sod/sod_oracle.py" "$VER_DIR/sod/sod-x.ini" "${WORK[@]}" --l1-max "$L1_MAX"
   echo "RV-3 PASSED ($TAG)"
fi
