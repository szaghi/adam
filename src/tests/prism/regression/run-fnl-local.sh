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
# MPI tuning
# ---------------------------------------------------------------------------
# `module load nvhpc` aliases `mpirun` with these flags, but bash aliases do
# NOT expand in non-interactive scripts — run.sh would call the bare binary
# and lose the tuning. Express the alias as MCA / UCX environment variables
# instead, which the bare mpirun picks up automatically:
#   --mca pml ucx               -> OMPI_MCA_pml=ucx
#   -x UCX_TLS=^cma             -> UCX_TLS=^cma
#   --mca coll_hcoll_enable 0   -> OMPI_MCA_coll_hcoll_enable=0
export OMPI_MCA_pml="ucx"
export UCX_TLS="^cma"
export OMPI_MCA_coll_hcoll_enable="0"

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
