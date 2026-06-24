#!/usr/bin/env bash
# rmf-amr M1 structural check (issue #13 §7.5, milestone M1).
#
# Unlike run.sh (digest/golden regression), this is a STRUCTURAL assertion: it
# drives the live IC-time AMR refine path from the AMR_GEO primitive-box marker
# under scheme_space=fv_centered and confirms a deterministic intra-realm 2:1
# coarse-fine jump forms. M1 has no golden — it de-risks the refine path before
# the M2 registration pass / M3 reflux depend on it.
#
# Asserts:
#   1. the run completes (no abort/NaN) and the fv_centered solve advances;
#   2. blocks_number > 1 (refinement happened);
#   3. BOTH a coarse level and a finer level are present (a 2:1 jump exists),
#      and a control run with the marker disabled produces NO such jump — so the
#      box marker, not the uniform iu_ref_levels, is what creates the seam.
#
# Usage: ./check.sh            (expects exe/adam_prism_cpu already built)
#        ./check.sh --build    (build prism-cpu-gnu first)
#
# mpirun and the GNU MPI toolchain must be on PATH (see run.sh header).

set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="$REPO_ROOT/exe/adam_prism_cpu"

if [[ "${1:-}" == "--build" ]]; then
   echo ">> building prism-cpu-gnu"
   (cd "$REPO_ROOT" && fobis build --mode prism-cpu-gnu)
fi
[[ -x "$EXE" ]] || { echo "ERROR: $EXE not found — run with --build" >&2; exit 2; }
command -v mpirun >/dev/null 2>&1 || { echo "ERROR: mpirun not on PATH" >&2; exit 2; }

# level(code) for refinement ratio 8: code<=0 -> 0 (ancestor), 1..7 -> 1, >=8 -> 2.
# Returns the sorted-unique set of levels present among the "code=+N" log lines.
levels_present() { grep -oE 'code=\+[0-9]+' "$1" \
   | awk -F+ '{c=$2; if(c<=0)l=0; else if(c<8)l=1; else l=2; print l}' | sort -un | tr '\n' ' '; }

run_in() { # workdir, ini-transform-sed
   local wd="$1" sed_expr="$2"
   rm -rf "$wd" && mkdir -p "$wd"
   sed "$sed_expr" "$CASE_DIR/input.ini" > "$wd/input.ini"
   ( cd "$wd" && timeout 300 mpirun -np 1 "$EXE" > run.log 2>&1 )
}

fail=0

# --- Main run: marker active. ---
WORK="$CASE_DIR/work-cpu"
echo ">> [rmf-amr] running marker-active case"
run_in "$WORK" 's/^$/&/'   # identity transform (copy as-is)

if grep -qiE 'error|abort| nan |segfault' "$WORK/run.log"; then
   echo "FAIL [rmf-amr] run reported an error/abort/NaN"; fail=1
fi
if ! grep -qE 'progress:[[:space:]]*100%' "$WORK/run.log"; then
   echo "FAIL [rmf-amr] fv_centered time loop did not reach 100%"; fail=1
fi
nblocks="$(grep -oE 'nodes number:[[:space:]]*[0-9]+' "$WORK/run.log" | grep -oE '[0-9]+' | tail -1)"
if [[ -z "$nblocks" || "$nblocks" -le 1 ]]; then
   echo "FAIL [rmf-amr] blocks_number not > 1 (got '${nblocks:-none}')"; fail=1
else
   echo ">> [rmf-amr] blocks_number = $nblocks (> 1)"
fi
main_levels="$(levels_present "$WORK/run.log")"
echo ">> [rmf-amr] refinement levels present (marker active): $main_levels"
if [[ "$main_levels" != *"1"* || "$main_levels" != *"2"* ]]; then
   echo "FAIL [rmf-amr] no coarse-fine jump: expected both level 1 and level 2"; fail=1
fi

# --- Control run: marker disabled. Proves the box (not iu_ref_levels) makes the jump. ---
CTRL="$CASE_DIR/work-cpu-ctrl"
echo ">> [rmf-amr] running marker-disabled control case"
run_in "$CTRL" 's/markers_number = 1/markers_number = 0/'
ctrl_levels="$(levels_present "$CTRL/run.log")"
echo ">> [rmf-amr] refinement levels present (marker off):   $ctrl_levels"
if [[ "$ctrl_levels" == *"2"* ]]; then
   echo "FAIL [rmf-amr] control produced a level-2 jump without the marker — fixture not isolating the box"; fail=1
fi

if [[ $fail -eq 0 ]]; then
   echo "PASS [rmf-amr] M1: deterministic intra-realm 2:1 jump forms under fv_centered, driven by the box marker"
   exit 0
else
   echo "FAIL [rmf-amr] M1 structural check failed"
   exit 1
fi
