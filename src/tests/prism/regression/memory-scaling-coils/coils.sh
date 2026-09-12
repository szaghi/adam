#!/usr/bin/env bash
# memory-scaling-coils — quantify the coil memory that the capacity budget ignores.
#
# THE QUESTION: does coil memory enter compute_blocks_number?
#
# It does not. compute_blocks_number runs inside realm_object%initialize
# (adam_prism_common_object.F90:335) with a hardcoded fields_number=80
# (adam_adam_object.F90:221, commented "remember to change"); coil%initialize runs
# later at :348, AFTER nb is fixed. Three arrays then grow linearly in coil count,
# each dimensioned on the per-rank nb:
#
#   J_vec        (3, ghosted-ijk, nb, nc)   adam_prism_coil_object.F90:141
#   j_vec_gpu    (nb, ghosted-ijk, 3, nc)   adam_prism_fnl_coil_object.F90:177
#   buf_6D_R8P   (nb, ghosted-ijk, 3, nc)   adam_prism_fnl_object.F90:565   [FNL only]
#
# The leading 3 is literal -- nv = size(coil%j_vec,dim=1) at fnl_coil_object:165,
# NOT the state-vector nv. The device array mirrors the host exactly.
#
# WHAT THIS ASSERTS
#   1. nb is IDENTICAL with and without coils, at the same rank count. That is the
#      coil-blindness itself: the budget cannot see them. The two inputs differ in
#      nothing else.
#   2. both runs complete -- i.e. the uncounted memory still fits inside the
#      save_factor=0.4 margin at this size. A failure of the 4-coil run is the
#      interesting one: it means the margin was consumed and the real ceiling
#      arrived before nb predicted.
#
# and REPORTS the measured footprint against the budget, so the omission is a
# number rather than an argument.
#
# This is a CALIBRATION defect, not a SCALABILITY one: it shifts the baseline, it
# does not break the doubling. memory-scaling/ covers the doubling.
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
         echo "       ^ THE INTERESTING FAILURE: the uncounted coil memory exhausted the"
         echo "         save_factor=0.4 margin. nb was chosen without it (fields_number=80,"
         echo "         adam_adam_object.F90:221), so the real ceiling arrived before nb predicted."
      fi
      fail=1
   fi
done
[[ $fail -ne 0 ]] && { echo "FAIL [memory-scaling-coils/$backend] runs did not complete"; exit 1; }

nb0="$(nb_of "$WD0/run.log")" ; nb4="$(nb_of "$WD4/run.log")"
bud="$(budget_of "$WD4/run.log")"
ni="$(ini_int "$CASE_DIR/input.ini" ni)"  ; nj="$(ini_int "$CASE_DIR/input.ini" nj)"
nk="$(ini_int "$CASE_DIR/input.ini" nk)"  ; ngc="$(ini_int "$CASE_DIR/input.ini" ngc)"
nc="$(ini_int "$CASE_DIR/input.ini" rectangular_coils_number)"

echo ">> [memory-scaling-coils/$backend] nb without coils = $nb0"
echo ">> [memory-scaling-coils/$backend] nb with $nc coils  = $nb4   budget = ${bud:-?} GB"

if [[ "$nb0" != "$nb4" ]]; then
   echo "FAIL [memory-scaling-coils/$backend] nb changed with coils: $nb0 -> $nb4."
   echo "                                     Either the budget now accounts for coil memory (good news,"
   echo "                                     but this oracle must then be re-derived), or the two inputs"
   echo "                                     differ in something other than coils — they must not."
   fail=1
else
   echo ">> [memory-scaling-coils/$backend] nb IDENTICAL — the budget is blind to coil memory (as expected)"
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
         printf ">> [footprint] device j_vec_gpu = %.2f GB of a %.2f GB budget (%.1f%%), UNCOUNTED when nb was chosen\n", tot/GB, bud, 100*(tot/GB)/bud
      else
         printf ">> [footprint] host J_vec = %.2f GB of a %.2f GB budget (%.1f%%), UNCOUNTED when nb was chosen\n", tot/GB, bud, 100*(tot/GB)/bud
   }
}'

echo ">> [memory-scaling-coils/$backend] the save_factor=0.4 margin (adam_adam_object.F90:164) is the only"
echo "                                   headroom absorbing this; fields_number=80 (:221) models block"
echo "                                   fields alone and is marked 'remember to change'."

if [[ $fail -eq 0 ]]; then
   echo "PASS [memory-scaling-coils/$backend] coil memory quantified; nb confirmed coil-blind"
   exit 0
else
   echo "FAIL [memory-scaling-coils/$backend] see above"
   exit 1
fi
