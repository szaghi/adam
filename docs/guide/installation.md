---
title: Installation
---

# Installation

ADAM is provided as source files and must be compiled to produce executables.

ADAM has two types of dependencies: sources available as git submodules (in `src/third_party`) and external libraries that must be installed on the system.

## Supported Compilers

| Compiler | Status |
|----------|--------|
| NVIDIA HPC SDK (nvfortran) ≥ 25.3 | Supported — see [Install NVIDIA HPC SDK](#install-nvidia-hpc-sdk) |
| GNU gfortran ≥ 9 | CPU builds only |
| Intel oneAPI (ifort/ifx) | GMP backend (experimental) |
| AMD flang | Planned |

## Obtain ADAM

Clone the repository, then initialise the submodules in a separate step:

```bash
git clone https://github.com/szaghi/adam
cd adam
git submodule update --init
```

## Compile ADAM

The preferred build system is [FoBiS](https://github.com/szaghi/FoBiS). Install it with:

```bash
pip install FoBiS.py
```

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
FoBiS.py build -lmodes
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
