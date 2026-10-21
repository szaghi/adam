#!/usr/bin/env bash
# Run the 7 FNL conserving-current tests sequentially with the COHERENT gather (field_weighting_model = name of the
# conserving solver: E and B gathered with the transpose of its charge operator S^rho, the adjoint of the charge deposit),
# one per conserving solver, TSC charge shape and FD6. Same inputs as ../gather_non_coherent except the gather:
#   conserving_test_TSC_Q10           esirkepov_Q10, no filter, 1e-11 tail cleanup on
#   conserving_test_TSC_Q8|Q6|Q5      esirkepov_Q8/Q6/Q5, no filter, cleanup on
#   conserving_test_TSC_totorica      totorica, no filter, cleanup on
#   conserving_test_TSC_londrillo     londrillo, no filter, cleanup off (J not local)
#   conserving_test_TSC_londrillo_f   londrillo + binomial filter = variant B7 of Metodi_conserving.pdf, cleanup off
# londrillo and londrillo_f: S^rho is the (filtered) particle shape, so the coherent gather equals the 2D one.
#
# Usage (from WSL):
#   ./run_all.sh                     # run all tests
#   ./run_all.sh totorica londrillo  # run only the listed tests: directory name, or suffix of conserving_test_TSC_<...>
#   NP=2 ./run_all.sh            # MPI ranks per test (default 1)
#
# Each test writes its stdout/stderr to <test_dir>/<test_name>.log; a failing test does not stop the
# sequence. A summary (with wall time per test) is appended to run_all_summary.log in this directory.
# Every test directory holds adam_prism_fnl, a link to exe/adam_prism_fnl (build prism-fnl-nvf first).

set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
NP="${NP:-1}"

# NVHPC 26.9 environment, the toolchain of the local_nvf varset (equivalent of `module load nvhpc`, no Lmod needed)
NVROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/26.9
export PATH="$NVROOT/cuda/bin:$NVROOT/compilers/bin:$NVROOT/comm_libs/mpi/bin:$PATH"
export LD_LIBRARY_PATH="$NVROOT/cuda/lib64:$NVROOT/compilers/lib:$NVROOT/math_libs/lib64:${LD_LIBRARY_PATH:-}"

# WSL2 GPU+MPI workarounds (see CLAUDE.md "WSL2 GPU+MPI Caveats"); do NOT copy to cluster job scripts
export OMPI_MCA_pml="ucx"
export UCX_TLS="^cma"
export OMPI_MCA_coll_hcoll_enable="0"
export UCX_RNDV_THRESH="inf"
export UCX_MEMTYPE_CACHE="n"

TESTS=(conserving_test_TSC_Q10 conserving_test_TSC_Q8 conserving_test_TSC_Q6 conserving_test_TSC_Q5 \
       conserving_test_TSC_totorica conserving_test_TSC_londrillo conserving_test_TSC_londrillo_f)
[[ $# -gt 0 ]] && TESTS=("$@")

SUMMARY="$HERE/run_all_summary.log"
echo "=== run_all started $(date '+%F %T') (np=$NP) ===" | tee -a "$SUMMARY"

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

echo "=== run_all finished $(date '+%F %T') ===" | tee -a "$SUMMARY"
