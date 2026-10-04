#!/usr/bin/env bash
# FLUME MHD verification MV-14 (issue #41, M2-P7c): multi-realm MHD and restart round trips.
#
# Why: the forest seam (#37) and the restart files carry the MHD state (B, psi) only if every channel was widened to
# nv; a missed one is silent. Three legs, np 2:
#   1. RJ2a (MV-4 input, GLM, N = 256 along x, t = 0.2) split at the diaphragm x = 0.5 into two realms (make_split.py:
#      mirror seam, beta cadence), against the single-realm run: the union must be BITWISE on all 9 fields
#      (multirealm_oracle.py, tolerance 0) and the summed conservation histories within round-off;
#   2. restart round trip of the 2-realm run: the straight run of leg 1 against a run stopped at half its steps with
#      restart files, then restarted to t = 0.2: interior cells BITWISE on all 9 fields of both realms, residuals,
#      conservation and div(B) histories byte-identical (the V7 pattern). In 1-D, GLM keeps psi exactly 0 (MV-4), so
#      this leg cannot see psi;
#   3. restart round trip of the MV-11 box (3-D, 2:1 AMR, reflux, GLM with damping, a seeded perturbation): psi is
#      non-zero there, and the round trip must be BITWISE on all 9 fields with max |psi| > 0 (--nonzero psi), histories
#      byte-identical;
#   4. (issue #52) RJ2a (GLM) on the cells of the Euler seam + AMR leg (verification/multirealm/sod-amr.ini: x > 0.75
#      refined 2:1), against the same cells split at the 2:1 face into a coarse and a fine realm glued by
#      `coupling = refined` (sod-amr-refined*.ini, the realm grids kept by make_rj2a.py): the 2:1 seam ghosts and the
#      2:1 inter-realm reflux carry all 9 fields, so the union must be BITWISE.
# Measured (CPU, M2-P7c): leg 1 bitwise over 297 steps, summed conservation within 1.7e-16; leg 2 (148 + restart)
# bitwise, 6 histories identical; leg 3 bitwise with max |psi| 1.7e-3. About 4.5 min on the CPU.
#
# Usage: ./check.sh [--np N]
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/../.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
NP=2

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np) NP="$2" ; shift 2 ;;
      *)    echo "check.sh: unknown argument '$1' (accepted: --np N)" >&2 ; exit 2 ;;
   esac
done
if [[ ! -x "$EXE" ]]; then
   echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2
   exit 2
fi
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
TAG="$(basename "$EXE")-np$NP"
ORACLE="$VERIF_DIR/conservation/conservation_oracle.py"
STATUS=0

run() { # run <work-dir> <ini> <log>
   if ! (cd "$1" && mpirun -np "$NP" "$EXE" "$2" > "$3" 2>&1); then
      echo "check.sh: run failed, see $1/$3" >&2
      exit 1
   fi
}
same_histories() { # same_histories <dir-a> <dir-b>: every history file of a byte-identical in b
   local f
   for f in $(cd "$1" && ls ./*.dat); do
      if cmp -s "$1/$f" "$2/$f"; then
         echo "   ${f#./}: identical  PASS"
      else
         echo "   ${f#./}: differs  FAIL"
         STATUS=1
      fi
   done
}
steps() { tail -1 "$1" | awk '{print $1+0}' ; } # last step of a residuals history

echo ">> MV-14 leg 1: RJ2a N = 256, 2 realms vs 1 ($(basename "$EXE"), np $NP)"
single="$CASE_DIR/work-$TAG-single" ; multi="$CASE_DIR/work-$TAG-2realm"
rm -rf "$single" "$multi" ; mkdir -p "$single" "$multi"
"$VENV_PY" "$VERIF_DIR/mhd/rj2a/make_rj2a.py" "$VERIF_DIR/sod/sod-x.ini" "$single/rj2a.ini" --axis x --cells 256 \
   --divergence-control glm
"$VENV_PY" "$CASE_DIR/make_split.py" "$single/rj2a.ini" "$multi" rj2a-2realm
run "$single" rj2a.ini log.txt
run "$multi" rj2a-2realm.ini log.txt
"$VENV_PY" "$VERIF_DIR/multirealm/multirealm_oracle.py" "$multi" "$single" --ngc 3 || STATUS=1

n=$(steps "$multi/rj2a-2realm-r1-residuals.dat") ; h=$(( n / 2 ))
echo ">> MV-14 leg 2: 2-realm restart round trip, $n steps vs $h + restart"
rst="$CASE_DIR/work-$TAG-2realm-restart"
rm -rf "$rst" ; mkdir -p "$rst"
"$VENV_PY" "$CASE_DIR/make_split.py" "$single/rj2a.ini" "$rst" rj2a-2realm --set "time.it_max=$h" \
   --set "IO.restart_save=$h"
run "$rst" rj2a-2realm.ini log-1.txt
"$VENV_PY" "$CASE_DIR/make_split.py" "$single/rj2a.ini" "$rst" rj2a-2realm --set "IO.restart=.true."
run "$rst" rj2a-2realm.ini log-2.txt
"$VENV_PY" "$ORACLE" --compare "$multi" "$rst" --tol 0 --ngc 3 || STATUS=1
same_histories "$multi" "$rst"

echo ">> MV-14 leg 3: MV-11 box (AMR, reflux, damped GLM, psi != 0) restart round trip, 20 steps vs 10 + restart"
box=(--divergence-control glm --b 0.8 0.5 -0.3 --set "mhd.glm_alpha=0.18")
make_box=("$VENV_PY" "$VERIF_DIR/mhd/zero-field/make_input.py" "$VERIF_DIR/conservation/amr-periodic.ini")
a="$CASE_DIR/work-$TAG-box-A" ; b="$CASE_DIR/work-$TAG-box-B"
rm -rf "$a" "$b" ; mkdir -p "$a" "$b"
"${make_box[@]}" "$a/input.ini" "${box[@]}" \
   --set "time.it_max=20"
run "$a" input.ini log.txt
"${make_box[@]}" "$b/input.ini" "${box[@]}" \
   --set "time.it_max=10" --set "IO.restart_save=10"
run "$b" input.ini log-1.txt
"${make_box[@]}" "$b/input.ini" "${box[@]}" \
   --set "time.it_max=20" --set "IO.restart=.true."
run "$b" input.ini log-2.txt
"$VENV_PY" "$ORACLE" --compare "$a" "$b" --tol 0 --ngc 3 --nonzero psi || STATUS=1
same_histories "$a" "$b"

echo ">> MV-14 leg 4: RJ2a on the sod-amr cells, refined (2:1) 2-realm split vs 1 realm (issue #52)"
mr="$VERIF_DIR/multirealm"
rs="$CASE_DIR/work-$TAG-amr-single" ; rm_="$CASE_DIR/work-$TAG-refined"
rm -rf "$rs" "$rm_" ; mkdir -p "$rs" "$rm_"
rj2a=("$VENV_PY" "$VERIF_DIR/mhd/rj2a/make_rj2a.py")
"${rj2a[@]}" "$mr/sod-amr.ini" "$rs/sod-amr.ini" --axis x --divergence-control glm
cp "$mr/sod-amr-refined.ini" "$rm_/"
for r in 1 2; do
   "${rj2a[@]}" "$mr/sod-amr-refined-r$r.ini" "$rm_/sod-amr-refined-r$r.ini" --axis x --divergence-control glm
done
run "$rs" sod-amr.ini log.txt
run "$rm_" sod-amr-refined.ini log.txt
"$VENV_PY" "$VERIF_DIR/multirealm/multirealm_oracle.py" "$rm_" "$rs" --ngc 3 || STATUS=1

for w in "$single" "$multi" "$rst" "$a" "$b" "$rs" "$rm_"; do find "$w" -name '*.h5' -delete; done
if [[ $STATUS -eq 0 ]]; then
   echo "MV-14 PASSED ($TAG)"
else
   echo "MV-14 FAILED ($TAG)"
   exit 1
fi
