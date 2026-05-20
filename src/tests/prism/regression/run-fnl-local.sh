#!/usr/bin/env bash
# Local FNL regression trigger — workstation-only convenience wrapper.
#
# GitHub-hosted runners have no GPU, so the FNL (OpenACC / nvfortran) backend
# cannot run in CI. This wrapper is the deliberate, manually-invoked gate:
# run it on a GPU workstation before an "important" push to confirm the FNL
# backend still matches its committed golden.
#
# Only two things distinguish the FNL run from the CPU run:
#   1. `module load nvhpc` — puts nvfortran and the bundled MPI on PATH.
#   2. fobis must build with `--varset local_nvf` — that varset defines the
#      nvf-built HDF5 prefix and NVF_CC. Those are fobos *variables*, not
#      shell environment variables; the wrapper does not export them.
# This wrapper does (1), sets MPI tuning (see below), then hands off to the
# backend-agnostic run.sh with `--varset local_nvf`.
#
# Usage:
#   ./run-fnl-local.sh                # build prism-fnl-nvf, run, diff FNL golden
#   ./run-fnl-local.sh --no-build     # skip the build (use existing exe/)
#
# Capturing / refreshing the FNL golden (only after an intentional change):
#   ./run-fnl-local.sh                # produces rmf/work-fnl/digest.txt
#   mkdir -p rmf/golden/fnl
#   cp rmf/work-fnl/digest.txt rmf/work-fnl/*-residuals.dat rmf/golden/fnl/
#
# Prerequisites on the workstation:
#   - Lmod with an `nvhpc` module (provides nvfortran + bundled mpirun).
#   - fobis 3.8+ on PATH.
#   - the `local_nvf` varset in the repo `fobos` resolves to an HDF5 build
#     that exists on disk (currently lib/hdf5/develop/nvf/26.1).

set -euo pipefail

REGRESSION_DIR="$(cd "$(dirname "$0")" && pwd)"

# ---------------------------------------------------------------------------
# NVHPC environment
# ---------------------------------------------------------------------------
# `module` is a shell function, not a binary — source Lmod's init so it is
# defined inside this non-interactive script.
for lmod_init in /usr/share/lmod/lmod/init/bash /etc/profile.d/lmod.sh; do
   if [[ -r "$lmod_init" ]]; then
      # shellcheck disable=SC1090
      source "$lmod_init"
      break
   fi
done
if ! command -v module >/dev/null 2>&1; then
   echo "ERROR: Lmod 'module' command not available — cannot load nvhpc" >&2
   exit 2
fi

echo ">> loading nvhpc module"
module load nvhpc

if ! command -v nvfortran >/dev/null 2>&1; then
   echo "ERROR: nvfortran not on PATH after 'module load nvhpc'" >&2
   exit 2
fi

# ---------------------------------------------------------------------------
# MPI / UCX tuning — WSL2-ONLY workarounds. Do NOT propagate to clusters.
# ---------------------------------------------------------------------------
# The nvhpc modulefile already exports all of these (and aliases mpirun), but
# bash aliases do NOT expand in non-interactive scripts and the module may not
# be loaded by whoever invokes run.sh — so re-export them here as a defensive
# fallback. Each is correct on this WSL2 dev box and a performance regression
# (or meaningless) on real InfiniBand + GPUDirect RDMA:
#
#   OMPI_MCA_pml=ucx               select the UCX point-to-point layer.
#   UCX_TLS=^cma                   exclude CMA: WSL's cross-memory-attach
#                                  same-node shmem transport is broken/restricted
#                                  (ptrace_scope), so do not let UCX pick it.
#   OMPI_MCA_coll_hcoll_enable=0   disable hcoll: Mellanox hardware-offloaded
#                                  collectives need IB HW this box lacks.
#   UCX_RNDV_THRESH=inf            issue #12: the FNL halo exchange posts MPI on
#                                  device-resident ghost buffers (GPU-direct).
#                                  With >=2 ranks UCX moves those via RENDEZVOUS,
#                                  whose device-memory transports (cuda_copy/gdr)
#                                  cannot get the GPU primary context through
#                                  WSL's /dev/dxg shim -> SIGABRT in
#                                  ucp_proto_rndv_send_start at the first
#                                  exchange. Forcing every message EAGER never
#                                  enters rendezvous, so no device pointer reaches
#                                  the broken path. On real IB+GPUDirect this is
#                                  a pure perf regression: rendezvous IS the fast
#                                  path there. There is no host-staging fallback
#                                  in the code; this env var is the only knob.
export OMPI_MCA_pml="ucx"
export UCX_TLS="^cma"
export OMPI_MCA_coll_hcoll_enable="0"
export UCX_RNDV_THRESH="inf"

# ---------------------------------------------------------------------------
# Hand off to the backend-agnostic harness.
# ---------------------------------------------------------------------------
# --varset local_nvf is the one piece of fobos config the FNL build needs;
# everything else (HDF5 prefix, NVF_CC) is resolved by fobos from that varset.
echo ">> nvfortran   = $(command -v nvfortran)"
echo ">> mpirun      = $(command -v mpirun)"
echo ">> handing off to run.sh fnl --varset local_nvf"
echo

exec "$REGRESSION_DIR/run.sh" fnl --varset local_nvf "$@"
