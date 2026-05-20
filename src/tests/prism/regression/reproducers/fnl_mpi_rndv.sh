#!/usr/bin/env bash
# Reproducer for issue #12 — FNL SIGABRT in the UCX rendezvous send path.
#
# ROOT CAUSE (diagnosed 2026-05-20, verified by isolation):
#   `field_fnl%update_ghost_mpi_gpu` posts MPI_Isend/Irecv directly on the
#   device-resident ghost buffers (GPU-direct / CUDA-aware MPI). With >=2 ranks
#   UCX moves those (large) device buffers via its RENDEZVOUS protocol, whose
#   device-memory transports (cuda_copy / gdr) must retrieve the GPU primary
#   context. On WSL2 the `/dev/dxg` shim cannot satisfy that, so the rendezvous
#   send aborts:
#       ucp_proto_rndv_send_start() ... -> SIGABRT (signal 6)
#   The abort surfaces at the FIRST ghost exchange (prism_fnl initialize_forest),
#   NOT inside MPI_Waitall, and NOT in a Fortran request-handle bug — the request
#   is valid; the fault is in the UCX rendezvous path on a broken WSL stack.
#
# WORKAROUND: UCX_RNDV_THRESH=inf forces every message EAGER, so the rendezvous
#   protocol is never entered and no device pointer reaches the broken transport.
#   This is a blunt, WSL-only crutch — a pure performance regression on real
#   InfiniBand + GPUDirect RDMA, where rendezvous IS the fast path. It lives in
#   run-fnl-local.sh and must never propagate to a cluster job script or the app.
#   There is NO host-staging path in the code (the old ADAM_FNL_HOST_STAGED_MPI
#   toggle was removed); this env var is the only knob.
#
# This script demonstrates BOTH:
#   * default (rendezvous)         -> reproduces the SIGABRT on a WSL box
#   * UCX_RNDV_THRESH=inf (eager)  -> runs clean
#
# The crash is `-np 2`-only (cross-rank GPU-buffer exchange); `-np 1` never hits
# the MPI path and always passes. The trigger needs >=2 ranks on a CUDA-aware
# MPI build whose UCX device-memory rendezvous cannot serve the GPU primary
# context (WSL2 today; possibly a degraded cluster node).
#
# Usage:
#   ./fnl_mpi_rndv.sh [iterations]      # default 20
#
# Prereqs: prism-fnl-nvf already built (exe/adam_prism_fnl), nvhpc on Lmod.

set -uo pipefail

ITERS="${1:-20}"
REPRO_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$REPRO_DIR/../../../../.." && pwd)"
EXE="$REPO_ROOT/exe/adam_prism_fnl"
CASE_INI_DIR="$REPRO_DIR/../rmf"

# --- nvhpc environment (mirror run-fnl-local.sh) ---------------------------
for lmod_init in /usr/share/lmod/lmod/init/bash /etc/profile.d/lmod.sh; do
   [[ -r "$lmod_init" ]] && { source "$lmod_init"; break; }
done
command -v module >/dev/null 2>&1 && module load nvhpc 2>/dev/null
export OMPI_MCA_pml="ucx" UCX_TLS="^cma" OMPI_MCA_coll_hcoll_enable="0"

[[ -x "$EXE" ]] || { echo "ERROR: $EXE not found — build prism-fnl-nvf first" >&2; exit 2; }

run_burst() {
   # $1 = label, $2 = UCX_RNDV_THRESH value ("" leaves the UCX default = rendezvous)
   local label="$1" rndv="$2" scratch pass=0 fail=0 out rc
   scratch="$(mktemp -d /tmp/adam_repro12_XXXX)"
   cp "$CASE_INI_DIR"/*.ini "$scratch"/ 2>/dev/null
   echo "=== $label: $ITERS x (mpirun -np 2) ==="
   ( cd "$scratch" || exit 1
     for ((i=1; i<=ITERS; i++)); do
        if [[ -n "$rndv" ]]; then export UCX_RNDV_THRESH="$rndv"; else unset UCX_RNDV_THRESH; fi
        out="$(timeout 180 mpirun -np 2 "$EXE" 2>&1)"; rc=$?
        if [[ $rc -eq 0 ]] && ! grep -qE "signal 6|signal 11|Aborted|Fatal|rndv|Segmentation" <<<"$out"; then
           pass=$((pass+1)); printf "%d " "$i"
        else
           fail=$((fail+1)); echo; echo "iter $i FAIL (rc=$rc):"
           grep -E "signal|Aborted|Fatal|rndv|Segmentation" <<<"$out" | head -3
           break
        fi
     done
     echo; echo ">> $label: pass=$pass fail=$fail" )
   rm -rf "$scratch"
}

run_burst "Rendezvous (UCX default — expected to SIGABRT on WSL2)" ""
echo
run_burst "Eager (UCX_RNDV_THRESH=inf — expected clean)" "inf"
