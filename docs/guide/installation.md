---
title: Installation
---

# Installation

ADAM is provided as source files and must be compiled to produce executables.

ADAM has two types of dependencies: source libraries fetched by FoBiS (into `src/third_party`) and external libraries that must be installed on the system.

## Supported Compilers

| Compiler | Status |
|----------|--------|
| NVIDIA HPC SDK, nvfortran **26.9** (HPC-X Open MPI bundled) | GPU (OpenACC, `*-fnl-nvf`) and CPU builds; see [Install NVIDIA HPC SDK](#install-nvidia-hpc-sdk) |
| GNU gfortran **16.2.0** + Open MPI **5.0.11** | CPU builds (`*-cpu-gnu`); CI uses Ubuntu 24.04's gfortran 14 (varset `ci_gnu`) |
| AMD flang (ROCm) | OpenMP offload builds (`*-omp-amd`) |
| Intel oneAPI (ifx) | OpenMP offload templates, experimental |

These are the versions the regression suites are run with. Fortran `.mod` files work only with the compiler that wrote
them, so HDF5 must be built with the same compiler and MPI as ADAM, and a build must start clean after a compiler
switch.

## Obtain ADAM

Clone the repository, then fetch source dependencies with FoBiS:

```bash
git clone https://github.com/szaghi/adam
cd adam
fobis fetch
```

`fobis fetch` reads `src/third_party/.deps_config.ini` and downloads all required source libraries into `src/third_party/`.

ADAM tracks the newest `main`/`master` of every dependency: no commit is pinned, and the `fobos.lock` that `fobis fetch`
writes stays local (it is gitignored). A fresh clone therefore gets the current heads, but an existing checkout does not
move by itself. Bring it to the current heads with:

```bash
fobis fetch --update
```

then rebuild clean, since the dependencies' module files change.

Alternatively, use the provided `install.sh` script which automates both steps:

```bash
# clone + fetch deps + build
scripts/install.sh --download git --build fobis --mode <mode>

# download a release tarball instead of cloning
scripts/install.sh --download wget --build fobis --mode <mode> --tag v1.0.0
```

Run `scripts/install.sh --help` for the full option list.

## Compile ADAM

The preferred build system is [FoBiS](https://github.com/szaghi/FoBiS) — the CLI binary is `fobis` (version 3.8+); flags use double-dash long form. Install it from PyPI:

```bash
pip install FoBiS.py
```

(`FoBiS.py` is the PyPI package name; it installs the `fobis` command.)

The repository root structure:

```
.
├── exe/
├── fobos
├── LICENSE.lgpl3.md
├── LICENSE.mit.md
├── README.md
├── scripts/
└── src/
```

List all available build modes:

```bash
fobis build --lmodes
```

Modes are named `<app>-<backend>-<compiler>[-debug]`, e.g. `prism-cpu-gnu`, `flume-fnl-nvf`, `flume-cpu-gnu-debug`.
Every mode needs a **varset** that points it at the HDF5 build of its toolchain (see
[Varsets](#varsets-one-per-machine-and-toolchain)); the `*-nvf` modes also take the GPU compute capability from it:

```bash
fobis build --mode prism-cpu-gnu                       # the default varset, local_gnu
fobis build --mode flume-fnl-nvf --varset local_nvf    # an NVF mode always needs an NVF varset
```

## Install NVIDIA HPC SDK

Download NVIDIA HPC SDK **26.9** from the [NVIDIA HPC SDK downloads page](https://developer.nvidia.com/hpc-sdk-downloads),
unpack it and run its `install` script; the default installation path is `/opt/nvidia/hpc_sdk/Linux_x86_64/26.9`.

Pick the CUDA bundle your GPU driver can run. A bundle whose CUDA is newer than the driver supports fails at the first
kernel launch: on a node limited to CUDA 12.9, use the `cuda_multi` bundle and build with `-gpu=cuda12.9`.

Then put the compilers and the bundled MPI on `PATH`, with the SDK's module files:

```bash
module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc/26.9
```

or directly:

```bash
NVHPC=/opt/nvidia/hpc_sdk/Linux_x86_64/26.9
export PATH=$NVHPC/compilers/bin:$NVHPC/comm_libs/mpi/bin:$PATH
export MANPATH=$NVHPC/compilers/man:$MANPATH
```

In a non-interactive shell (a batch job, a script) `module load` may not populate `PATH`; export it directly there.

### Known issues: NVHPC 26.9 with HPC-X 2.50

Both were reported on Ubuntu 24.04 under WSL2 with the SDK's bundled HPC-X 2.50 (Open MPI 5) in
[#96](https://github.com/szaghi/adam/issues/96).

- **`prterun` needs `libnl-route-3`.** Without it every `mpirun` stops with
  `libnl-route-3.so.200: cannot open shared object file`. Install it once:

  ```bash
  sudo apt-get install libnl-route-3-200
  ```

- **Debug modes die in `MPI_Init`.** The `*-fnl-nvf-debug` modes (`template-nvf-debug`) compile and link with
  `-Ktrap=fp`, and a floating-point exception raised inside UCX's `ucp_worker_create` during `MPI_Init` traps
  with signal 8, before any ADAM code runs. Until that changes, debug with a local copy of the mode without
  `-Ktrap=fp`; the bounds, pointer and stack checks (`-Mbounds -Mchkptr -Mchkstk`) still work.

## Install External Dependencies

ADAM needs parallel HDF5, with zlib and szip, built with the **same** compiler and MPI as ADAM itself.
`scripts/hdf5_build.sh` fetches the HDF5 `develop` branch with szip and zlib and builds all three with the `mpicc` and
`mpif90` found on `PATH`, so load the toolchain first. The repository's own builds use HDF5 2.3.0 (`develop` at
`ab5df552`).

```bash
# with the toolchain on PATH: -sdk NVF for nvfortran, -sdk GNU for gfortran
scripts/hdf5_build.sh -build -sdk NVF -abs-paths -lsrc /tmp/hdf5-src -hdf5 $PWD/lib/hdf5/develop/nvf/26.9
```

On Rocky Linux or RHEL, HDF5 2.3 installs into `lib64/` while ADAM expects `lib/`: add
`-DCMAKE_INSTALL_LIBDIR=lib` to the script's `cmake` line. An HDF5 built on a newer glibc than the machine that links
ADAM fails with `undefined reference to __isoc23_strtol`; build it on that machine.

### Varsets: one per machine and toolchain

The HDF5 path is the fobos variable `$HDF5_PREFIX`, defined per `[varset:*]` in the repository `fobos`. It is a fobos
variable, not a shell environment variable: exporting it has no effect. NVF varsets also set `$NVF_CC`, the
`-gpu=ccXX` compute capability. Selecting a varset selects the HDF5 build:

| Varset | Machine and toolchain | `$HDF5_PREFIX` | `$NVF_CC` |
|---|---|---|---|
| `local_gnu` (default) | workstation: gfortran 16.2.0 + Open MPI 5.0.11 | `lib/hdf5/develop/gnu/16.2.0` | — |
| `local_nvf` | workstation: nvfortran 26.9 | `lib/hdf5/develop/nvf/26.9` | `cc89` |
| `ci_gnu` | GitHub Actions: Ubuntu 24.04 gfortran 14 + Open MPI | `lib/hdf5/develop/gnu/14.2.0` | — |
| `iac_gnu` | medooza cluster: gfortran 16.2.0 + Open MPI 5.0.11 | `lib/hdf5/develop/iac/gnu/16.2.0/openmpi-5.0.11` | — |
| `iac_nvf` | medooza cluster: nvfortran 26.9, CUDA 12.9 | `lib/hdf5/develop/iac/nvf/26.9` | `cc80` |
| `leonardo` | Leonardo, Spack NVHPC 24.5 | Spack path | `cc80` |
| `spacehpc` | Cray PE | `/opt/cray/pe/hdf5/1.14.3.3/nvidia/23.3` | `cc90` |

The AMD builds use `local_amd` or `thera_amd` together with a GPU-architecture varset (`gfx90a`, `gfx942`, ...), e.g.
`--varset "thera_amd gfx942"`.

The toolchain must match the varset exactly. A mismatch stops the build at `use hdf5` ("created by a different version
of GNU Fortran") or, under nvfortran, with an internal compiler error on a module symbol. Bring up a new machine by
adding its own `[varset:*]` block; the [regression bring-up tutorial](/tests/prism-regression-tutorial) walks through
it.
