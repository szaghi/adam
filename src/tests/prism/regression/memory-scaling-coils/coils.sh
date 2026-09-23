#!/usr/bin/env bash
# memory-scaling-coils — check that the capacity budget counts coil memory.
#
# THE QUESTION: does coil memory enter compute_blocks_number?
#
# History: it did not. fields_number was hardcoded to 80 (adam_adam_object.F90,
# "remember to change") and coil%initialize ran after nb was fixed, so the coil
# arrays below grew outside the budget. On FNL that tipped the coil cases into a
# device out-of-memory once 3db858fb budgeted from total rather than free memory.
# Now prism_common_object%compute_fields_number counts them before nb is chosen:
#
#   J_vec        (3, ghosted-ijk, nb, nc)   adam_prism_coil_object.F90   [host]
#   j_vec_gpu    (nb, ghosted-ijk, 3, nc)   adam_prism_fnl_coil_object.F90 [device]
#
# 3 block-sized fields per coil, on the budgeted memory space of each backend.
#
# WHAT THIS ASSERTS (the two inputs differ in nothing but the coils)
#   1. fields_number grows by exactly 3 per coil: F(nc coils) = F(0) + 3*nc.
#   2. nb scales by the inverse ratio: nb(nc) = nb(0) * F(0)/F(nc), within 1 block
#      (each nb is rounded to the nearest integer independently).
#   3. both runs complete.
#
# and REPORTS the coil footprint against the budget.
#
# Usage:
#   ./coils.sh [--cpu|--fnl] [--np N] [--build]
#   PRISM_EXE=/path/to/exe ./coils.sh --fnl
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
SUITE_DIR="$(cd "$CASE_DIR/.." && pwd)"
BASE_DIR="$SUITE_DIR/memory-scaling"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"

backend="cpu" ; np=1 ; do_build=0
while [[ $# -gt 0 ]]; do
   case "$1" in
      --fnl)   backend="fnl" ;;
      --cpu)   backend="cpu" ;;
      --build) do_build=1 ;;
      --np)    shift; np="${1:?--np needs a value}" ;;
      *) echo "ERROR: unknown flag $1 (use --cpu / --fnl / --np N / --build)" >&2; exit 2 ;;
   esac
   shift
done

if [[ "$backend" == "fnl" ]]; then
   EXE="${PRISM_EXE:-$REPO_ROOT/exe/adam_prism_fnl}"
   BUDGET_RE='device memory_total \[GB\]'   # the budget is TOTAL; memory_avail is free, used only by the OOM diagnostic
else
   EXE="${PRISM_EXE:-$REPO_ROOT/exe/adam_prism_cpu}"
   BUDGET_RE='\]  memory_avail \[GB\]'
fi

if [[ $do_build -eq 1 ]]; then
   if [[ "$backend" == "fnl" ]]; then
      (cd "$REPO_ROOT" && fobis build --mode prism-fnl-nvf --varset local_nvf)
   else
      (cd "$REPO_ROOT" && fobis build --mode prism-cpu-gnu)
   fi
fi
[[ -x "$EXE" ]] || { echo "ERROR: $EXE not found — run with --build" >&2; exit 2; }
[[ -f "$BASE_DIR/input.ini" ]] || { echo "ERROR: base case missing: $BASE_DIR/input.ini" >&2; exit 2; }
command -v mpirun >/dev/null 2>&1 || { echo "ERROR: mpirun not on PATH" >&2; exit 2; }

nb_of()     { grep -a -m1 'all blocks number (nb)' "$1" | grep -oE '[+-]?[0-9]+$' | tr -d '+'; }
budget_of() { grep -a -m1 -E "$BUDGET_RE" "$1" | grep -oE '[+-]?0\.[0-9]+E[+-][0-9]+' | head -1; }
fields_of() { grep -a -m1 'initialize fields_number:' "$1" | grep -oE '[+-]?[0-9]+$' | tr -d '+'; }
ini_int()   { grep -m1 -E "^\s*$2\s*=" "$1" | grep -oE '[0-9]+' | head -1; }

run_case() { # label, case-dir -> echoes workdir
   local label="$1" src="$2"
   local wd="$CASE_DIR/work-${backend}-${label}-np${np}"
   rm -rf "$wd" && mkdir -p "$wd"
   cp "$src/input.ini" "$wd/"
   ( cd "$wd" && timeout 300 mpirun -np "$np" "$EXE" > run.log 2>&1 ) || true
   find "$wd" -type f \( -name '*.h5' -o -name '*.fbd' -o -name '*.xdmf' -o -name '*.tnd' \) -delete
   echo "$wd"
}

fail=0
echo ">> [memory-scaling-coils/$backend] np$np — comparing 0 coils vs 4 coils"

WD0="$(run_case nocoil "$BASE_DIR")"
WD4="$(run_case coils  "$CASE_DIR")"

for spec in "0:$WD0" "4:$WD4"; do
   n="${spec%%:*}" ; wd="${spec#*:}"
   if grep -aqiE 'error stop|segfault' "$wd/run.log"; then
      echo "FAIL [memory-scaling-coils/$backend] the ${n}-coil run failed:"
      grep -aiE 'error stop|segfault' "$wd/run.log" | head -2 | sed 's/^/       /'
      if [[ "$n" == "4" ]] && grep -aqi 'failed to allocate j_vec_gpu' "$wd/run.log"; then
         echo "       ^ device out-of-memory on the coil array: the budget undercounts what the"
         echo "         realm allocates (see prism_common_object%compute_fields_number)."
      fi
      fail=1
   fi
done
[[ $fail -ne 0 ]] && { echo "FAIL [memory-scaling-coils/$backend] runs did not complete"; exit 1; }

nb0="$(nb_of "$WD0/run.log")" ; nb4="$(nb_of "$WD4/run.log")"
fn0="$(fields_of "$WD0/run.log")" ; fn4="$(fields_of "$WD4/run.log")"
bud="$(budget_of "$WD4/run.log")"
ni="$(ini_int "$CASE_DIR/input.ini" ni)"  ; nj="$(ini_int "$CASE_DIR/input.ini" nj)"
nk="$(ini_int "$CASE_DIR/input.ini" nk)"  ; ngc="$(ini_int "$CASE_DIR/input.ini" ngc)"
nc="$(ini_int "$CASE_DIR/input.ini" rectangular_coils_number)"

echo ">> [memory-scaling-coils/$backend] nb without coils = $nb0"
echo ">> [memory-scaling-coils/$backend] nb with $nc coils  = $nb4   budget = ${bud:-?} GB"

echo ">> [memory-scaling-coils/$backend] fields_number without coils = $fn0, with $nc coils = $fn4"

if [[ -z "$fn0" || -z "$fn4" ]]; then
   echo "FAIL [memory-scaling-coils/$backend] fields_number not found in the run logs"
   fail=1
elif (( fn4 != fn0 + 3*nc )); then
   echo "FAIL [memory-scaling-coils/$backend] fields_number grew by $((fn4 - fn0)), expected 3*$nc = $((3*nc))"
   fail=1
else
   expected="$(awk -v n="$nb0" -v f0="$fn0" -v f4="$fn4" 'BEGIN{printf "%d", n*f0/f4 + 0.5}')"
   delta=$(( nb4 - expected )) ; delta=${delta#-}
   if (( delta > 1 )); then
      echo "FAIL [memory-scaling-coils/$backend] nb with coils = $nb4, expected $expected = nb0*$fn0/$fn4 (+-1)"
      fail=1
   else
      echo ">> [memory-scaling-coils/$backend] nb scales with the counted coils: $nb0 -> $nb4 (expected $expected)"
   fi
fi

# footprint: 3 components * ghosted block * 8 bytes, per block per coil
awk -v nb="$nb4" -v nc="$nc" -v ni="$ni" -v nj="$nj" -v nk="$nk" -v ngc="$ngc" -v bud="${bud:-0}" -v be="$backend" '
BEGIN {
   gb  = (ni+2*ngc) * (nj+2*ngc) * (nk+2*ngc)
   per = 3 * gb * 8                       # bytes per block, per coil
   tot = per * nb * nc
   MB = 1024*1024 ; GB = MB*1024
   printf ">> [footprint] ghosted block = %d cells; per block-coil = %d B (3 comps x %d x 8)\n", gb, per, gb
   printf ">> [footprint] per rank, %d coils at nb=%d:\n", nc, nb
   printf ">>   J_vec      (host)   %8.1f MB\n", tot/MB
   if (be == "fnl") {
      printf ">>   j_vec_gpu  (device) %8.1f MB\n", tot/MB
      printf ">>   buf_6D_R8P (host)   %8.1f MB\n", tot/MB
   }
   printf ">>   cost of ONE coil      %8.1f MB per array\n", (per*nb)/MB
   if (bud+0 > 0) {
      if (be == "fnl")
         printf ">> [footprint] device j_vec_gpu = %.2f GB of a %.2f GB budget (%.1f%%), counted in nb\n", tot/GB, bud, 100*(tot/GB)/bud
      else
         printf ">> [footprint] host J_vec = %.2f GB of a %.2f GB budget (%.1f%%), counted in nb\n", tot/GB, bud, 100*(tot/GB)/bud
   }
}'

if [[ $fail -eq 0 ]]; then
   echo "PASS [memory-scaling-coils/$backend] coil memory counted in the budget"
   exit 0
else
   echo "FAIL [memory-scaling-coils/$backend] see above"
   exit 1
fi
