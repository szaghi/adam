#!/usr/bin/env bash
# FLUME verification NV-6 (issue #49, N2d): Sod in interstellar units, limiting active at the shock (issue #48, A3).
#
# Why: with the classic WENO weights, an absolute zeps = 1e-6 dominates the smoothness indicators of a density of
# 1e-21 kg/m^3 (IS ~ 1e-42): the weights collapse to the linear ones and the scheme stops limiting, and the
# weno-riemann sensor, which reads the same weights, goes blind. #49 gives two cures: the scale-invariant weights
# ([weno] weights = si, N1) and the reference layer ([reference], N2). The V1 Sod problem (sod/sod-x.ini) is written in
# SI units, L0 = 1 pc, u0 = 10 km/s, rho0 = 1e-21 kg/m^3 (scaling.py physical; the references are not powers of two),
# and run on both flux paths (weno, weno-riemann HLLC with the 6th-order correction and the weno sensor):
#   1. raw SI, weights si          == code units, weights si   (limiting unchanged by the units, round-off)
#   2. SI with [reference], js     == code units, js           (the layer hands the solver the code numbers)
#   3. raw SI, weights js: total variation > 1.1 x code units, js, or a stop on a non-finite state (negative control:
#      the collapse #48 A3 predicted, which shows the case detects a scheme that does not limit; the unlimited
#      weno-riemann run undershoots towards vacuum and, round-off deciding, survives on the CPU and stops on the FNL)
# Comparisons of rho / rho0 against x / L0, tolerance 1e-12 (interstellar_oracle.py).
#
# Measured at N2d (CPU, np 2), TV of rho / rho0 (code units js: 0.877 weno, 0.884 weno-riemann): raw SI js 1.162 and
# 2.455 (new extrema 2e-2 and 9e-2; weno-riemann L1 x20); legs 1 and 2 within 1.1e-13.
#
# Usage: ./interstellar.sh [--np N]
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./interstellar.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
TOOL="$CASE_DIR/scaling.py"
ORACLE="$CASE_DIR/interstellar_oracle.py"
SOD="$CASE_DIR/../sod/sod-x.ini"
L0=3.0857e16 ; U0=1.0e4 ; RHO0=1.0e-21
NP=2

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np) NP="$2" ; shift 2 ;;
      *) echo "interstellar.sh: unknown argument '$1'" >&2 ; exit 2 ;;
   esac
done
[[ -x "$EXE" ]] || { echo "interstellar.sh: executable '$EXE' not found" >&2 ; exit 2 ; }
TAG="$(basename "$EXE")-np$NP"

run() { # run [--may-stop] <work-dir> <scaling.py physical options...>: write the input and run it
   local may_stop=0
   [[ $1 == --may-stop ]] && { may_stop=1 ; shift ; }
   local work="$1" ; shift
   rm -rf "$work" ; mkdir -p "$work"
   "$VENV_PY" "$TOOL" physical "$SOD" "$work/input.ini" "$@"
   if [[ -n $SCHEME ]]; then
      "$VENV_PY" - "$work/input.ini" <<'PYEOF'
import configparser, sys
ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None, strict=False)
ini.optionxform = str
ini.read(sys.argv[1])
ini["numerics"].update({"scheme_space": "weno-riemann", "riemann_solver": "hllc", "flux_correction": "6th",
                        "flux_correction_sensor": "weno"})
with open(sys.argv[1], "w") as f:
    ini.write(f)
PYEOF
   fi
   if ! (cd "$work" && mpirun -np "$NP" "$EXE" input.ini > log.txt 2>&1); then
      [[ $may_stop -eq 1 ]] && grep -aq "non-finite (NaN or infinite) values in the state" "$work/log.txt" && return 0
      echo "interstellar.sh: run failed, see $work/log.txt" >&2
      exit 1
   fi
}

fails=0
echo ">> NV-6 interstellar Sod ($(basename "$EXE"), np $NP): L0 = $L0 m, u0 = $U0 m/s, rho0 = $RHO0 kg/m^3"
for SCHEME in "" weno-riemann; do
   path="${SCHEME:-weno}"
   w="$CASE_DIR/work-nv6-$TAG-$path"
   physical=(--length "$L0" --velocity "$U0" --density "$RHO0")
   run "$w-code-js"              --length 1 --velocity 1 --density 1 --weights js
   run "$w-code-si"              --length 1 --velocity 1 --density 1 --weights si
   run "$w-physical-si"          "${physical[@]}" --weights si
   run "$w-reference-js"         "${physical[@]}" --weights js --reference
   run --may-stop "$w-physical-js" "${physical[@]}" --weights js
   echo "   [$path] 1. raw SI, si weights == code units, si weights"
   "$VENV_PY" "$ORACLE" --ini "$SOD" --references "$L0" "$U0" "$RHO0" --same "$w-code-si" "$w-physical-si" || \
      fails=$(( fails + 1 ))
   echo "   [$path] 2. SI with [reference], js weights == code units, js weights"
   "$VENV_PY" "$ORACLE" --ini "$SOD" --references "$L0" "$U0" "$RHO0" --same "$w-code-js" "$w-reference-js" || \
      fails=$(( fails + 1 ))
   echo "   [$path] 3. raw SI, js weights: the collapse (negative control)"
   stop="$(grep -aom1 "non-finite (NaN or infinite) values in the state at step +[0-9]*" "$w-physical-js/log.txt" || true)"
   if [[ -n $stop ]]; then
      echo "   collapse detected: the run stopped on ${stop#*values in the }  PASS"
   else
      "$VENV_PY" "$ORACLE" --ini "$SOD" --references "$L0" "$U0" "$RHO0" --collapse "$w-code-js" "$w-physical-js" || \
         fails=$(( fails + 1 ))
   fi
   find "$w"-* -name '*.h5' -delete
done
if [[ $fails -gt 0 ]]; then
   echo "NV-6 FAILED: $fails assertions ($TAG)"
   exit 1
fi
echo "NV-6: limiting active in interstellar units with si weights and with [reference], collapse detected without ($TAG)"
