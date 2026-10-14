#!/usr/bin/env bash
# Run the FNL TSC conserving-current tests of the esirkepov_Q<N> filter pairs sequentially, one after the other:
#   conserving_test_TSC_Q5    Q5/P10  (current_conserving_solver = esirkepov_Q5)
#   conserving_test_TSC_Q6    Q6/P11  (current_conserving_solver = esirkepov_Q6)
#   conserving_test_TSC_Q8    Q8/P13  (current_conserving_solver = esirkepov_Q8)
# Same case as conserving_test_TSC (esirkepov_Q10, the former esirkepov-modified): TSC, no filter, 1e-11 tail cleanup on.
#
# Usage (from WSL):
#   ./run_Q_filters.sh             # run the three tests
#   ./run_Q_filters.sh Q5 Q8       # run only the listed tests: directory name, or suffix of conserving_test_TSC_<...>
#   ./run_Q_filters.sh conserving_test_TSC   # any test directory of this folder (e.g. the Q10 reference)
#   NP=2 ./run_Q_filters.sh        # MPI ranks per test (default 1)
#
# Each test writes its stdout/stderr to <test_dir>/<test_name>.log; a failing test does not stop the
# sequence. A summary (with wall time per test) is appended to run_Q_filters_summary.log in this directory.

set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
NP="${NP:-1}"

# NVHPC 25.11 environment (equivalent of `module load nvhpc`, no Lmod needed)
NVROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/25.11
export PATH="$NVROOT/cuda/bin:$NVROOT/compilers/bin:$NVROOT/comm_libs/mpi/bin:$PATH"
export LD_LIBRARY_PATH="$NVROOT/cuda/lib64:$NVROOT/compilers/lib:$NVROOT/math_libs/lib64:${LD_LIBRARY_PATH:-}"

# WSL2 GPU+MPI workarounds (see CLAUDE.md "WSL2 GPU+MPI Caveats"); do NOT copy to cluster job scripts
export OMPI_MCA_pml="ucx"
export UCX_TLS="^cma"
export OMPI_MCA_coll_hcoll_enable="0"
export UCX_RNDV_THRESH="inf"
export UCX_MEMTYPE_CACHE="n"

TESTS=(conserving_test_TSC_Q5 conserving_test_TSC_Q6 conserving_test_TSC_Q8)
[[ $# -gt 0 ]] && TESTS=("$@")

SUMMARY="$HERE/run_Q_filters_summary.log"
echo "=== run_Q_filters started $(date '+%F %T') (np=$NP) ===" | tee -a "$SUMMARY"

for t in "${TESTS[@]}"; do
   name="$t" ; [[ -d "$HERE/$name" ]] || name="conserving_test_TSC_$t"
   dir="$HERE/$name"
   if [[ ! -f "$dir/input.ini" ]]; then
      echo "[SKIP] $name: no input.ini" | tee -a "$SUMMARY"
      continue
   fi
   echo ">> $(date '+%F %T') start $name" | tee -a "$SUMMARY"
   t0=$(date +%s)
   (cd "$dir" && mpirun -np "$NP" ./adam_prism_fnl > "$name.log" 2>&1)
   rc=$?
   dt=$(( $(date +%s) - t0 ))
   status=OK ; [[ $rc -ne 0 ]] && status="FAIL(rc=$rc)"
   printf '[%s] %s  wall=%02d:%02d:%02d\n' "$status" "$name" $((dt/3600)) $((dt%3600/60)) $((dt%60)) | tee -a "$SUMMARY"
done

echo "=== run_Q_filters finished $(date '+%F %T') ===" | tee -a "$SUMMARY"
