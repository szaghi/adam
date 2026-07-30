#!/usr/bin/env bash
# Local AMD regression trigger — Thera (AMD MI300 / gfx942) workstation wrapper.
#
# The GitHub-hosted CI runners are x86/NVIDIA only, so neither the amdflang CPU
# backend nor the ROCm OpenMP-offload GPU backend can run there. This wrapper is
# the deliberate, manually-invoked AMD gate: run it on an AMD box to confirm the
# amdflang toolchain builds and the produced fields still match the committed
# golden references (within the usual cross-compiler / cross-runtime tolerance).
#
# It mirrors run-fnl-local.sh, but instead of `module load nvhpc` it puts a
# self-consistent AMD toolchain on PATH:
#
#   * amdflang / amdclang   — AMD "AFAR" flang drop (LLVM flang 23.0.0).
#   * OpenMPI 5.0.10        — whose mpif90/mpicc wrappers drive that same
#                             amdflang, and against which the HDF5 below links.
#   * parallel HDF5 (Fortran ON) built with that exact amdflang + OpenMPI —
#                             resolved by fobos from `--varset thera_amd`
#                             ($HDF5_PREFIX). fobos *variables* are not shell
#                             env vars, so the wrapper does not export them.
#
# The three must agree: Fortran .mod files are compiler-specific, so the HDF5
# .mod files only load under the amdflang that built them, and HDF5 must link
# the same MPI the app links. The paths below are the toolchain the thera_amd
# HDF5 in `fobos` was built against — keep them in lock-step with that varset.
#
# Backend selector (first positional arg, default `amd`):
#   amd       CPU backend   -> fobis mode prism-cpu-amd    -> diffed vs golden/cpu
#   amd-omp   GPU offload   -> fobis mode prism-fnl-omp-amd -> diffed vs golden/fnl
#
# Usage:
#   ./run-amd-local.sh                 # CPU amdflang build, run, diff vs golden/cpu
#   ./run-amd-local.sh amd --no-build  # reuse exe/
#   ./run-amd-local.sh amd-omp         # ROCm OpenMP-offload build, run, diff vs golden/fnl
#
# GPU_ARCH: prism-fnl-omp-amd compiles for $GPU_ARCH (fobos [common-variables],
# default gfx942 = MI300, which matches this node's gfx94X drop). Override for a
# different GPU by adding the arch varset, e.g.
#   ./run-amd-local.sh amd-omp --varset "thera_amd gfx90a"

set -euo pipefail

REGRESSION_DIR="$(cd "$(dirname "$0")" && pwd)"

# ---------------------------------------------------------------------------
# AMD toolchain locations (must match the thera_amd HDF5 build in `fobos`)
# ---------------------------------------------------------------------------
# OpenMPI 5.0.10 built with amdflang; its mpif90 wrapper points at the AFAR
# amdflang drop below. This is the MPI the thera_amd HDF5 links against.
AMD_OMPI_PREFIX="${AMD_OMPI_PREFIX:-/home/mbycklin/code/spack/install/linux-zen4/openmpi-5.0.10-3rvtw7bdkzbm3tkapamurnizsr5v7hjn}"
# AFAR amdflang/amdclang drop (LLVM flang 23.0.0, gfx94X). Only its lib dir is
# needed on the runtime/link path; mpif90 already invokes its amdflang by path.
AMD_FLANG_PREFIX="${AMD_FLANG_PREFIX:-/home/mbycklin/amd/therock-23.2.1-gfx94X-7.13.0-7357b5084b}"
# ROCm runtime (provides libamdhip64 for the offload build's `-lamdhip64` and at
# exec time). rocm/6.4.1 on this box.
ROCM_PATH="${ROCM_PATH:-/opt/rocm-6.4.1}"

# ---------------------------------------------------------------------------
# Backend selection
# ---------------------------------------------------------------------------
BACKEND="amd"
if [[ $# -gt 0 && "$1" != --* ]]; then
   BACKEND="$1"
   shift
fi
case "$BACKEND" in
   amd | amd-omp) ;;
   *)
      echo "ERROR: run-amd-local.sh backend must be 'amd' or 'amd-omp' (got '$BACKEND')" >&2
      exit 2
      ;;
esac

# ---------------------------------------------------------------------------
# Put the AMD toolchain on PATH / library path
# ---------------------------------------------------------------------------
if [[ ! -x "$AMD_OMPI_PREFIX/bin/mpif90" ]]; then
   echo "ERROR: mpif90 not found at $AMD_OMPI_PREFIX/bin — fix AMD_OMPI_PREFIX" >&2
   exit 2
fi
export PATH="$AMD_OMPI_PREFIX/bin:$AMD_FLANG_PREFIX/bin:$ROCM_PATH/bin:$PATH"
export LD_LIBRARY_PATH="$AMD_OMPI_PREFIX/lib:$AMD_FLANG_PREFIX/lib:$ROCM_PATH/lib:${LD_LIBRARY_PATH:-}"

# The spack OpenMPI 5.0.10 has hard NEEDED deps on librdmacm.so.1 / libibverbs.so.1
# (rdma-core / OFED). Compute nodes ship these in the default loader path; login
# nodes on Thera do NOT, so the executable dies with "librdmacm.so.1: cannot open
# shared object file". Only when they are unresolvable do we fall back to the
# in-tree MLNX OFED copy (harmless no-op on a compute node, and same-node runs
# never touch the RDMA fast path anyway).
if ! ldconfig -p 2>/dev/null | grep -q 'librdmacm\.so\.1'; then
   OFED_LIB="${AMD_OFED_LIB:-/share/modules/mlnxofed/5.2-2.2.0.0/lib64}"
   if [[ -e "$OFED_LIB/librdmacm.so.1" ]]; then
      echo ">> librdmacm not in loader cache — adding OFED fallback: $OFED_LIB"
      export LD_LIBRARY_PATH="$OFED_LIB:$LD_LIBRARY_PATH"
   fi
fi

# Make the OpenMPI wrappers drive the AMD compilers explicitly (defensive: the
# wrapper defaults already point here, but a stray site OMPI_* in the env would
# otherwise win).
export OMPI_FC="$AMD_FLANG_PREFIX/bin/amdflang"
export OMPI_CC="$AMD_FLANG_PREFIX/bin/amdclang"
export OMPI_CXX="$AMD_FLANG_PREFIX/bin/amdclang++"

if ! command -v mpirun >/dev/null 2>&1; then
   echo "ERROR: mpirun not on PATH after AMD env setup" >&2
   exit 2
fi

echo ">> backend     = $BACKEND"
echo ">> amdflang    = $OMPI_FC"
echo ">> mpif90      = $(command -v mpif90)  (FC -> $("$AMD_OMPI_PREFIX/bin/mpif90" -show 2>/dev/null | awk '{print $1}'))"
echo ">> mpirun      = $(command -v mpirun)"
echo ">> handing off to run.sh $BACKEND --varset thera_amd"
echo

# ---------------------------------------------------------------------------
# Hand off to the backend-agnostic harness. --varset thera_amd supplies the
# amdflang-built HDF5 prefix (a fobos variable, resolved by fobos, not exported).
# ---------------------------------------------------------------------------
exec "$REGRESSION_DIR/run.sh" "$BACKEND" --varset thera_amd "$@"
