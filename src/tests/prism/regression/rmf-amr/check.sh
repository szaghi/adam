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
#   4. (M2) the forest registers > 0 intra-realm AMR seam faces with the marker
#      active, and exactly 0 in the marker-disabled control — proving the
#      registration pass keys off the real coarse-fine tree adjacency.
#   5. (#28) register parity np1 vs np2: at -np 2 every registered face must be
#      printed (the issue-#23 matched `reflux face` diagnostics) exactly once
#      per step across the rank union, with values identical to the np1 run.
#      Registration aliasing, missing cross-rank fine-sum reduction, or an
#      un-gated reflux apply all break this assertion.
#
# Usage: ./check.sh            (expects exe/adam_prism_cpu already built)
#        ./check.sh --build    (build prism-cpu-gnu first)
#        PRISM_EXE=<path> ./check.sh   (check a different backend executable,
#                                       e.g. exe/adam_prism_fnl with the WSL
#                                       UCX env from run-fnl-local.sh)
#
# mpirun and the GNU MPI toolchain must be on PATH (see run.sh header).

set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${PRISM_EXE:-$REPO_ROOT/exe/adam_prism_cpu}"

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

run_in() { # workdir, ini-transform-sed [, nprocs (default 1)]
   local wd="$1" sed_expr="$2" np="${3:-1}"
   rm -rf "$wd" && mkdir -p "$wd"
   sed "$sed_expr" "$CASE_DIR/input.ini" > "$wd/input.ini"
   ( cd "$wd" && timeout 300 mpirun -np "$np" "$EXE" > run.log 2>&1 )
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
# M2: intra-realm AMR seam faces must be registered (> 0).
amr_faces="$(grep -oE 'registered intra-realm AMR seam faces: \+?[0-9]+' "$WORK/run.log" | grep -oE '[0-9]+' | tail -1)"
echo ">> [rmf-amr] registered intra-realm AMR seam faces (marker active): ${amr_faces:-none}"
if [[ -z "$amr_faces" || "$amr_faces" -le 0 ]]; then
   echo "FAIL [rmf-amr] M2 registration pass registered no intra-realm AMR seam faces"; fail=1
fi
# M3: the Berger-Colella reflux must fire with a NON-ZERO coarse/fine flux
# mismatch (the 2:1 jump makes F_coarse != F_fine_sum). A same-resolution mirror
# seam would stay round-off zero and emit nothing.
reflux_delta="$(grep -oE 'reflux max\|F_coarse-F_fine_sum\| = [0-9.E+-]+' "$WORK/run.log" | grep -oE '[0-9.E+-]+$' | tail -1)"
echo ">> [rmf-amr] reflux max|F_coarse-F_fine_sum| (marker active): ${reflux_delta:-none}"
if [[ -z "$reflux_delta" ]] || ! awk "BEGIN{exit !(${reflux_delta:-0} > 0)}"; then
   echo "FAIL [rmf-amr] M3 reflux did not fire a non-zero correction across the 2:1 seam"; fail=1
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
# M2: with no coarse-fine jump, the registration pass must register exactly 0 faces.
ctrl_faces="$(grep -oE 'registered intra-realm AMR seam faces: \+?[0-9]+' "$CTRL/run.log" | grep -oE '[0-9]+' | tail -1)"
echo ">> [rmf-amr] registered intra-realm AMR seam faces (marker off):   ${ctrl_faces:-none}"
if [[ "${ctrl_faces:-0}" -ne 0 ]]; then
   echo "FAIL [rmf-amr] control registered AMR seam faces without a coarse-fine jump"; fail=1
fi
# M3: with no jump, the reflux diagnostic must be silent (no non-zero correction).
if grep -q 'reflux max|F_coarse-F_fine_sum|' "$CTRL/run.log"; then
   echo "FAIL [rmf-amr] control fired a reflux correction without a coarse-fine jump"; fail=1
fi

# --- #28 register-parity leg (default-on): np2 register content must equal np1. ---
# The per-face `reflux face N coarse_block M max|F_coarse-F_fine_sum|` prints
# (issue #23 R3 matched diagnostics) are the acceptance instrument: field digests
# cannot see the reflux correction at this scale (#23 R2 finding), so register
# content is asserted directly. At np>1 each face must be printed EXACTLY ONCE
# per step across the rank union — by its coarse-owner rank — with per-step
# values identical to np1's (block groups are rank-complete at this
# decomposition, so the cross-rank reduce adds exact zeros; string equality is
# the intended strength). Duplicated faces, 0.0-stuck entries, or shifted
# values are the issue-#28 defect signatures.
REG1="$CASE_DIR/work-reg-np1"
REG2="$CASE_DIR/work-reg-np2"
echo ">> [rmf-amr] running register-parity np1 run"
run_in "$REG1" 's/^$/&/' 1
echo ">> [rmf-amr] running register-parity np2 run"
run_in "$REG2" 's/^$/&/' 2
for wd in "$REG1" "$REG2"; do
   if ! grep -qE 'progress:[[:space:]]*100%' "$wd/run.log"; then
      echo "FAIL [rmf-amr] register-parity run in $(basename "$wd") did not reach 100%"; fail=1
   fi
done
# Registration must be decomposition-invariant: same face count at np1 and np2.
reg_count() { grep -oE 'registered intra-realm AMR seam faces: \+?[0-9]+' "$1" | grep -oE '[0-9]+' | tail -1; }
np1_faces="$(reg_count "$REG1/run.log")"
np2_faces="$(reg_count "$REG2/run.log")"
if [[ "${np1_faces:-0}" -ne "${np2_faces:-1}" || "${np1_faces:-0}" -le 0 ]]; then
   echo "FAIL [rmf-amr] #28 registered face count differs: np1=${np1_faces:-none} np2=${np2_faces:-none}"; fail=1
fi
# Per-face per-step value sequences (log order = step order; each face's np2
# prints come from a single owner rank, so mpirun's per-rank ordering holds).
face_values() { # log, face-index -> newline list of printed max|...| values
   grep -E "reflux face $2 " "$1" | awk '{print $NF}'
}
for f in $(seq 1 "${np1_faces:-0}"); do
   v1="$(face_values "$REG1/run.log" "$f")"
   v2="$(face_values "$REG2/run.log" "$f")"
   n1="$(printf '%s\n' "$v1" | grep -c . || true)"
   n2="$(printf '%s\n' "$v2" | grep -c . || true)"
   if [[ "$n1" -ne "$n2" ]]; then
      echo "FAIL [rmf-amr] #28 face $f printed $n2 times at np2 vs $n1 at np1 (must be exactly once per step)"; fail=1
   elif [[ "$v1" != "$v2" ]]; then
      echo "FAIL [rmf-amr] #28 face $f np2 values differ from np1:"
      paste <(printf '%s\n' "$v1") <(printf '%s\n' "$v2") | sed 's/^/       np1 vs np2: /'
      fail=1
   else
      echo ">> [rmf-amr] #28 face $f: np2 == np1 ($n1 steps, exact match)"
   fi
done
# WSL disk discipline: the parity runs' checkpoints are not needed downstream.
rm -f "$REG1"/*.h5 "$REG2"/*.h5

if [[ $fail -eq 0 ]]; then
   echo "PASS [rmf-amr] M1-M3 structural + #28 np1/np2 register parity"
   exit 0
else
   echo "FAIL [rmf-amr] structural / register-parity check failed"
   exit 1
fi
