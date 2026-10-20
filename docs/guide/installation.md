---
title: Installation
---

# Installation

ADAM is provided as source files and must be compiled to produce executables.

ADAM has two types of dependencies: source libraries fetched by FoBiS (into `src/third_party`) and external libraries that must be installed on the system.

## Supported Compilers

| Compiler | Status |
|----------|--------|
| NVIDIA HPC SDK (nvfortran) ≥ 25.3 | Supported — see [Install NVIDIA HPC SDK](#install-nvidia-hpc-sdk) |
| GNU gfortran ≥ 9 | CPU builds only |
| Intel oneAPI (ifort/ifx) | GMP backend (experimental) |
| AMD flang | Planned |

## Obtain ADAM

Clone the repository, then fetch source dependencies with FoBiS:

```bash
git clone https://github.com/szaghi/adam
cd adam
fobis fetch --frozen
```

`fobis fetch` downloads the source libraries listed in the `fobos` `[dependencies]` section into `src/third_party/`.
With `--frozen` it checks each one out at the commit pinned in `src/third_party/fobos.lock`, which is tracked in the
repository: the dependency state ADAM was tested with. Without `--frozen` it takes each dependency's current head and
only warns where that differs from the lock.

To move the pins, run `fobis fetch --update`, which fetches the current heads and rewrites the lock. Rebuild, run the
regression suites, and commit the new `fobos.lock` with the change that needs it.

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

See [Architecture → Build Commands](/guide/architecture#build-commands) for a full list of build targets.

## Install NVIDIA HPC SDK

Download the latest NVIDIA HPC SDK:

```bash
wget https://developer.download.nvidia.com/hpc-sdk/25.3/nvhpc_2025_253_Linux_x86_64_cuda_12.8.tar.gz
```

Decompress and install:

```bash
tar xf nvhpc_2025_253_Linux_x86_64_cuda_12.8.tar.gz
cd nvhpc_2025_253_Linux_x86_64_cuda_12.8/
sudo ./install
```

The default installation path is `/opt/nvidia/hpc_sdk/Linux_x86_64/25.3`.

After installation, update your PATH. If the module system is available:

```bash
module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc/25.3
module load nvhpc/25.3
```

Alternatively, set environment variables directly:

```bash
export PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/25.3/compilers/bin:$PATH
export MANPATH=/opt/nvidia/hpc_sdk/Linux_x86_64/25.3/compilers/man:$MANPATH
export PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/25.3/comm_libs/mpi/bin:$PATH
```

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

### SZIP Library

```bash
git clone https://github.com/erdc/szip
cd szip
export CC=nvc
export CXX=nvc++
export FC=nvfortran
export CFLAGS='-O3'
export CXXFLAGS='-O3'
export FCFLAGS='-O3'
./configure --prefix=/opt/szip/bin/2.1.1/nvf/25.3
make
make check
make install
cd ..
```

::: tip
The installation directory must be created before `make install` and you must have proper permissions.
:::

### ZLIB Library

```bash
wget https://zlib.net/current/zlib.tar.gz
tar xf zlib.tar.gz
cd zlib
export CC=nvc
export CFLAGS='-O3 -fPIC'
./configure --prefix=/opt/zlib/bin/1.3.1/nvf/25.3
make
make check
make install
cd ..
```

### HDF5 Library

Once ZLIB and SZIP are installed:

```bash
git clone https://github.com/HDFGroup/hdf5
cd hdf5
CFLAGS="-fPIC" FCFLAGS="-fPIC" ./configure \
  --prefix=/opt/HDF5/bin/1.14.6/nvf/25.3 \
  --enable-shared --enable-parallel --enable-fortran \
  --disable-libtool-lock \
  --with-szlib=/opt/szip/bin/2.1.1/nvf/25.3 \
  --with-zlib=/opt/zlib/bin/1.3.1/nvf/25.3 \
  CC=mpicc FC=mpif90
make
make check
make install
cd ..
```
