#!/usr/bin/env bash
# memory-scaling — capacity-scaling probe for the nb / nodes_number algorithm.
#
# Runs the minimal case at -np 1 and -np 2 and reports, per rank count:
#   budget      the memory handed to compute_blocks_number
#   nb          per-rank block capacity      (adam_field_object.F90:272)
#   aggregate   nb * np, i.e. nodes_number   (adam_adam_object.F90:166)
#
# and asserts the BACKEND-SPECIFIC expectation. The two backends must behave
# DIFFERENTLY, which is the whole point of the probe:
#
#   CPU  host RAM is SHARED by the ranks on a node, and mpih_object divides
#        /proc/meminfo's MemTotal by the MPI_COMM_TYPE_SHARED communicator size
#        (the ranks on THIS node). Correct -> nb HALVES as ranks-per-node double,
#        aggregate STAYS FLAT within a node. Extra ranks on one node buy no extra
#        grid; extra NODES must.
#
#   GPU  each rank owns its OWN device; adam_mpih_nvf_object.F90:63 does NOT
#        divide. Correct -> nb CONSTANT, aggregate DOUBLES. Two GPUs buy twice
#        the grid.
#
# A regression that conflates the two memory models — dividing the device budget
# by procs_number, or forgetting to divide the host one — flips one of these and
# is caught here.
#
# Usage:
#   ./scaling.sh                 # CPU backend (default)
#   ./scaling.sh --fnl           # GPU backend; caller owns nvhpc PATH + UCX knobs
#   ./scaling.sh --build         # build the selected backend first
#   PRISM_EXE=/path/to/exe ./scaling.sh --fnl   # explicit override
#
# EXIT 0 = observed scaling matches the backend's memory model.
#
# NOT a golden test: no field values are compared, and run.sh skips this case
# (no golden/<backend>/). Run it by hand, like the check.sh oracles.
#
# CAVEAT recorded, not asserted: the device budget is FREE memory, not total —
# FUNDAL's dev_init passes dev_memory_avail positionally into
# dev_get_device_memory_info(mem_free, mem_total) (fundal_dev_handling.F90:97),
# binding it to ACC_PROPERTY_FREE_MEMORY. GPU `nb` therefore depends on whatever
# else is resident at startup. The probe prints the budget so drift is visible,
# but asserts only the SCALING RATIO, never an absolute value.
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"

backend="cpu" ; do_build=0
for arg in "$@"; do
   case "$arg" in
      --fnl)   backend="fnl" ;;
      --cpu)   backend="cpu" ;;
      --build) do_build=1 ;;
      *) echo "ERROR: unknown flag $arg (use --cpu / --fnl / --build)" >&2; exit 2 ;;
   esac
done

if [[ "$backend" == "fnl" ]]; then
   EXE="${PRISM_EXE:-$REPO_ROOT/exe/adam_prism_fnl}"
   BUDGET_RE='device memory_avail \[GB\]'
else
   EXE="${PRISM_EXE:-$REPO_ROOT/exe/adam_prism_cpu}"
   BUDGET_RE='\]  memory_avail \[GB\]'
fi

if [[ $do_build -eq 1 ]]; then
   if [[ "$backend" == "fnl" ]]; then
      echo ">> building prism-fnl-nvf"
      (cd "$REPO_ROOT" && fobis build --mode prism-fnl-nvf --varset local_nvf)
   else
      echo ">> building prism-cpu-gnu"
      (cd "$REPO_ROOT" && fobis build --mode prism-cpu-gnu)
   fi
fi
[[ -x "$EXE" ]] || { echo "ERROR: $EXE not found — run with --build" >&2; exit 2; }
command -v mpirun >/dev/null 2>&1 || { echo "ERROR: mpirun not on PATH" >&2; exit 2; }

# nb is printed by field_object%description() as:
#     [mpi-00000]  all blocks number (nb):          +682
nb_of()     { grep -m1 'all blocks number (nb)' "$1" | grep -oE '[+-]?[0-9]+$'; }
nodeproc_of() { grep -m1 'node procs_number' "$1" | grep -oE '[+-]?[0-9]+$'; }
worldproc_of(){ grep -m1 -E '\]  procs_number:' "$1" | grep -oE '[+-]?[0-9]+$'; }
budget_of() { grep -m1 -E "$BUDGET_RE" "$1" | grep -oE '[+-]?0\.[0-9]+E[+-][0-9]+' | head -1; }

run_np() { # rank count -> echoes workdir
   local np="$1" wd="$CASE_DIR/work-${backend}-np${np}"
   rm -rf "$wd" && mkdir -p "$wd"
   cp "$CASE_DIR/input.ini" "$wd/"
   ( cd "$wd" && timeout 300 mpirun -np "$np" "$EXE" > run.log 2>&1 ) || true
   find "$wd" -type f \( -name '*.h5' -o -name '*.fbd' -o -name '*.xdmf' -o -name '*.tnd' \) -delete
   echo "$wd"
}

fail=0
declare -A NB BUDGET NODEP WORLDP

for np in 1 2; do
   wd="$(run_np "$np")"
   if grep -qiE 'error stop|abort|segfault' "$wd/run.log"; then
      echo "FAIL [memory-scaling/$backend] np$np run reported an error/abort"
      grep -iE 'error stop|abort|segfault' "$wd/run.log" | head -3 | sed 's/^/       /'
      fail=1 ; continue
   fi
   nb="$(nb_of "$wd/run.log" || true)"
   bud="$(budget_of "$wd/run.log" || true)"
   if [[ -z "$nb" ]]; then
      echo "FAIL [memory-scaling/$backend] np$np: could not parse 'all blocks number (nb)' from run.log"
      fail=1 ; continue
   fi
   NB[$np]="${nb#+}" ; BUDGET[$np]="${bud:-unparsed}"
   NODEP[$np]="$(nodeproc_of "$wd/run.log" | tr -d '+' || true)"
   WORLDP[$np]="$(worldproc_of "$wd/run.log" | tr -d '+' || true)"
   echo ">> [memory-scaling/$backend] np$np  budget = ${BUDGET[$np]} GB  nb = ${NB[$np]}  aggregate = $(( ${NB[$np]} * np ))"
done

if [[ $fail -ne 0 ]]; then
   echo "FAIL [memory-scaling/$backend] could not complete both runs"
   exit 1
fi

nb1="${NB[1]}" ; nb2="${NB[2]}"
agg1=$(( nb1 )) ; agg2=$(( nb2 * 2 ))
echo ">> [memory-scaling/$backend] nb: $nb1 -> $nb2   aggregate: $agg1 -> $agg2"

# 10% band: nint() rounding and a drifting free-memory budget both perturb nb.
within() { awk "BEGIN{d=($1-$2); if(d<0)d=-d; exit !(d <= $3*($2))}"; }

if [[ "$backend" == "cpu" ]]; then
   echo ">> [memory-scaling/cpu] expect nb to HALVE and aggregate to stay FLAT (shared host RAM)"
   if ! within "$nb2" "$(( nb1 / 2 ))" 0.10; then
      echo "FAIL [memory-scaling/cpu] nb did not halve: $nb1 -> $nb2 (expected ~$(( nb1 / 2 )))."
      echo "                          mpih_object.F90:136 divides host RAM by procs_number; if that"
      echo "                          division was removed, ranks on one node now overcommit shared RAM."
      fail=1
   fi
   if ! within "$agg2" "$agg1" 0.10; then
      echo "FAIL [memory-scaling/cpu] aggregate capacity changed: $agg1 -> $agg2 (expected flat)."
      echo "                          On shared host RAM, adding ranks cannot create capacity."
      fail=1
   fi
else
   echo ">> [memory-scaling/fnl] expect nb CONSTANT and aggregate to DOUBLE (one device per rank)"
   if ! within "$nb2" "$nb1" 0.10; then
      echo "FAIL [memory-scaling/fnl] nb changed with rank count: $nb1 -> $nb2 (expected ~constant)."
      echo "                          Each rank owns its own device, so the device budget must NOT be"
      echo "                          divided by procs_number (adam_mpih_nvf_object.F90:63). A halving"
      echo "                          here means the host memory model leaked into the GPU path, and"
      echo "                          doubling the GPUs no longer doubles the affordable grid."
      fail=1
   fi
   if ! within "$agg2" "$(( agg1 * 2 ))" 0.10; then
      echo "FAIL [memory-scaling/fnl] aggregate capacity did not double: $agg1 -> $agg2."
      echo "                          This is the weak-scaling contract: 2 GPUs must afford 2x the grid."
      fail=1
   fi
fi

# --- multi-node leg (self-disabling) -----------------------------------------
# On ONE node every rank shares the same /proc/meminfo, so node_procs_number ==
# procs_number and the host-memory division is trivially right. The interesting
# case is >1 node: get_memory_info reads THIS node's RAM (PENF
# penf_allocatable_memory.F90:3799) while procs_number counts the WHOLE world, so
# dividing by the world size under-estimates each rank's share BY THE NODE COUNT
# -- 2 nodes x 4 ranks gave mem/8 where the true share is mem/4, and adding nodes
# bought no capacity at all. adam_mpih_object.F90 now divides by the
# MPI_COMM_TYPE_SHARED communicator size instead.
#
# This leg cannot run on a single-node workstation, so it skips itself rather
# than passing vacuously. On a cluster, launch across >1 node to arm it.
np2_node="${NODEP[2]:-}" ; np2_world="${WORLDP[2]:-}"
if [[ -z "$np2_node" || -z "$np2_world" ]]; then
   echo ">> [memory-scaling/$backend] multi-node leg: SKIPPED (could not parse node/world rank counts)"
elif [[ "$np2_node" == "$np2_world" ]]; then
   echo ">> [memory-scaling/$backend] multi-node leg: SKIPPED — all $np2_world rank(s) on one node"
   echo "                             (node_procs_number == procs_number; run across >1 node to arm it)"
else
   echo ">> [memory-scaling/$backend] multi-node leg: ARMED — $np2_world ranks over $(( np2_world / np2_node )) nodes"
   if [[ "$backend" == "cpu" ]]; then
      # Each rank's budget must reflect ITS OWN node's RAM divided by the ranks on
      # THAT node -- not the world size. Equivalently nb must match what the same
      # per-node rank count yields on one node.
      if [[ "${NODEP[1]:-}" == "${NODEP[2]:-}" ]] && ! within "${NB[2]}" "${NB[1]}" 0.10; then
         echo "FAIL [memory-scaling/cpu] same ranks-per-node but nb changed: ${NB[1]} -> ${NB[2]}."
         echo "                          The host budget is still being divided by the WORLD size"
         echo "                          (adam_mpih_object.F90): each rank under-estimates its share"
         echo "                          by the node count, so extra nodes buy no capacity."
         fail=1
      fi
   fi
fi

if [[ $fail -eq 0 ]]; then
   echo "PASS [memory-scaling/$backend] capacity scaling matches the backend's memory model"
   exit 0
else
   echo "FAIL [memory-scaling/$backend] capacity scaling violates the backend's memory model"
   exit 1
fi
