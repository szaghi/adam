#!/usr/bin/env bash
# rmf-2realm-fv-pulse-blocks oracle (issue #51, inter-realm mirror seam between blocks that do not line up).
#
# THE CLAIM UNDER TEST: on the FV path, a 1:1 mirror seam whose two sides have different block sizes (the same cells)
# keeps the flux register exact. Each fine-side (realm_b) block face is scattered into every register face it overlaps
# (maps%seam_overlap, flux_register%accumulate_fine_overlaps); on a mirror seam the fluxes of the two sides are the
# same, so every register face must close: max|F_coarse-F_fine_sum| = 0, the line PRISM prints per face and step.
#
# ACCEPTANCE, for each way of misaligning (small blocks in realm 1, the register side, then in realm 2):
#   1. the run reaches 100% with no error/abort/NaN;
#   2. the reflux diagnostic is printed and its maximum is exactly 0;
#   3. the union of the realms matches the single-realm FV pulse within FIELD_TOL (relative): PRISM FV is not
#      invariant under a change of block layout at round-off (the single realm alone on 4 blocks per axis instead of 2
#      differs by 6.6e-16), so the fields are not compared bitwise; the register check is the exact one.
#
# Measured (issue #51, CPU and FNL, np 2, 20 steps): reflux mismatch 0 on all faces and steps both ways (320 and 80
# lines); fields within 3.7e-15 of the single realm.
#
# The same acceptance holds with one realm on a SINGLE block (`one1`, `one2`: iu_ref_levels = 0, 16 x 32 x 32 cells,
# max_level kept 1): the tree lookup of the seam peers then walks up from a level-1 code to the root, and aborted
# (error -111) before issue #54.
#
# Usage: ./check.sh [--build] [--np N]
#
# PRISM_EXE: override the executable under test, e.g.
#   PRISM_EXE=$REPO/exe/adam_prism_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH plus the
# WSL UCX knobs of issue #12); --build always builds the CPU default.
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${PRISM_EXE:-$REPO_ROOT/exe/adam_prism_cpu}"
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
ORACLE="$REPO_ROOT/src/tests/flume/verification/multirealm/multirealm_oracle.py"
FIELD_TOL="1.0e-13"
NP=2

do_build=0
while [[ $# -gt 0 ]]; do
   case "$1" in
      --build) do_build=1 ; shift ;;
      --np)    NP="$2" ; shift 2 ;;
      *) echo "ERROR: unknown flag $1 (use --build, --np N)" >&2; exit 2 ;;
   esac
done
if [[ $do_build -eq 1 ]]; then
   echo ">> building prism-cpu-gnu"
   (cd "$REPO_ROOT" && fobis build --mode prism-cpu-gnu)
fi
[[ -x "$EXE" ]] || { echo "ERROR: $EXE not found — run with --build" >&2; exit 2; }
command -v mpirun >/dev/null 2>&1 || { echo "ERROR: mpirun not on PATH" >&2; exit 2; }

TAG="$(basename "$EXE")-np$NP"
FV='s/^scheme_space             = fd_centered/scheme_space             = fv_centered/; s/^fdv_scheme = fd/fdv_scheme = fv/'
SMALL='s/^ni     = 8/ni     = 4/; s/^nj     = 16/nj     = 8/; s/^nk     = 16/nk     = 8/; s/^iu_ref_levels  = 1/iu_ref_levels  = 2/; s/^max_level      = 1/max_level      = 2/'
ONE='s/^ni     = 8/ni     = 16/; s/^nj     = 16/nj     = 32/; s/^nk     = 16/nk     = 32/; s/^iu_ref_levels  = 1/iu_ref_levels  = 0/'
fail=0

run_in() { # workdir
   ( cd "$1" && timeout 900 mpirun -np "$NP" "$EXE" > run.log 2>&1 ) || true
   if grep -qiE 'error|abort| nan |segfault' "$1/run.log"; then
      echo "FAIL [rmf-2realm-fv-pulse-blocks] $(basename "$1") reported an error/abort/NaN"; fail=1
   fi
   if ! grep -qE 'progress:[[:space:]]*100%' "$1/run.log"; then
      echo "FAIL [rmf-2realm-fv-pulse-blocks] $(basename "$1") did not reach 100%"; fail=1
   fi
}

single="$CASE_DIR/work-$TAG-single"
rm -rf "$single" ; mkdir -p "$single"
sed "$FV; s/^it_max   = 5$/it_max   = 20/; s/^it_save                = 5/it_save                = 20/; s/^markers_number = 1/markers_number = 0/; s/^max_level      = 2/max_level      = 1/; s/^amr_iterations  = 1 .*/amr_iterations  = 0/; s/^save_residual_fields   = .true./save_residual_fields   = .false./; s/^save_divergence_fields = .true./save_divergence_fields = .false./" \
   "$CASE_DIR/../rmf-amr-fd-pulse/input.ini" > "$single/input.ini"
echo ">> [rmf-2realm-fv-pulse-blocks] single-realm FV pulse ($TAG)"
run_in "$single"

for leg in small1 small2 one1 one2; do
   small="${leg: -1}"
   w="$CASE_DIR/work-$TAG-$leg"
   rm -rf "$w" ; mkdir -p "$w"
   cp "$CASE_DIR/input.ini" "$CASE_DIR/realm_1.ini" "$CASE_DIR/realm_2.ini" "$w/"
   if [[ $leg == small* ]]; then
      sed "$SMALL" "$CASE_DIR/realm_$small.ini" > "$w/realm_$small.ini"
      echo ">> [rmf-2realm-fv-pulse-blocks] split, realm $small on the small blocks ($TAG)"
   else
      sed "$ONE" "$CASE_DIR/realm_$small.ini" > "$w/realm_$small.ini"
      echo ">> [rmf-2realm-fv-pulse-blocks] split, realm $small on a single block (issue #54, $TAG)"
   fi
   run_in "$w"
   mism="$(grep -ah 'max|F_coarse-F_fine_sum|' "$w/run.log" | awk '{v=$NF+0; if (v<0) v=-v; if (v>m) m=v} END {printf "%.6E %d", m+0, NR}')"
   echo "   reflux max|F_coarse-F_fine_sum| = ${mism% *} over ${mism#* } face-step lines (must be 0)"
   if [[ "${mism#* }" == "0" ]] || ! awk "BEGIN{exit !(${mism% *} == 0)}"; then
      echo "FAIL [rmf-2realm-fv-pulse-blocks] the register does not close on the misaligned mirror seam"; fail=1
   fi
   "$VENV_PY" "$ORACLE" "$w" "$single" --ngc 3 --tol "$FIELD_TOL" --fields-only | sed 's/^/   /' || fail=1
done

for w in "$CASE_DIR"/work-"$TAG"-*; do
   find "$w" -type f \( -name '*.h5' -o -name '*.fbd' -o -name '*.xdmf' -o -name '*.tnd' \) -delete
done
if [[ $fail -eq 0 ]]; then
   echo "PASS [rmf-2realm-fv-pulse-blocks] misaligned mirror seam: register closed, fields within $FIELD_TOL ($TAG)"
   exit 0
fi
echo "FAIL [rmf-2realm-fv-pulse-blocks] ($TAG)"
exit 1
