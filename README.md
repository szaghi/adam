<a name="top"></a>
# ADAM

> ADAM, Adaptive Mesh Refinement (AMR) with Immersed Boundary (IB) fluid dynamic solver tailored for High Performance GPU Computing

### Authors

+ Andrea di Mascio, [andrea.dimascio@univaq.it](andrea.dimascio@univaq.it)
+ Federico Negro, [federico.negro.01@gmail.com](federico.negro.01@gmail.com)
+ Giacomo Rossi, [giacomo.rossi@intel.com](giacomo.rossi@intel.com)
+ Francesco Salvadore, [f.salvadore@cineca.it](mailto:f.salvadore@cineca.it)
+ Stefano Zaghi, [stefano.zaghi@cnr.it](stefano.zaghi@cnr.it)

| [Main Features](#main-features) | [Copyrights](#copyrights) | [Install](#install) | [Test](#test) | [API Documentation](#api-documentation) |

Go to [Top](#top)

### Main Features

+ KISS, keep it simple and stupid;
+ simulate different types of flows;
+ use high order finite difference scheme;
+ exploit Adaptive Mesh Refinement (AMR) to accurate simulate complex geometries;
+ exploit Immersed Boundary (IB) techniques to easy handle complex moving geometries;
+ ready for very High Performance Computing (HPC) by means of GPU parallel paradigms exploitation;

# Copyrights

ADAM is currently a closed project:

> Copyright (C) Di Mascio/Negro/Rossi/Salvadore/Zaghi, Inc - All Rights Reserved.
>
> Unauthorized copying of these source files, via any medium is strictly prohibited, proprietary and confidential.
> Written by Andrea di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore and Stefano Zaghi, April 2025.

Future versions could be released with a more Free Open Source Software (FOSS) licence.

Go to [Top](#top)

# Install

The ADAM framework is provided as source files archive and it must be compiled in order to have the executable ready to be installed.

ADAM has two type of dependencies: the first one is provided by sources directly available as git submodules (in `src/third_party` directory) while
second are external libraries that must be installed in the used system and available in user space; the installation of this kind of external dependencies are
documented in [Install External Dependencies](#install-external-dependencies).

The compiler currently supported are:

[x] NVIDIA HPC SDK, in particular nvfortran compiler, see [Install NVIDIA HPC SDK](#install-nvidia-hpc-sdk);
[ ] GCC;
[ ] Intel One Api;
[ ] AMD flang;

## Obtain ADAM framework

The best way to obtain ADAM framework is cloning its github repository:

```bash
git clone --recurse-submodules https://github.com/szaghi/adam
```

Once cloned, it can be compiled. However, you need the right compiler and all external dependencies must be installed previously,
see [Install NVIDIA HPC SDK](#install-nvidia-hpc-sdk) and [Install External Dependencies](#install-external-dependencies).

## Compile ADAM framework

To compile ADAM framework the preferred method is to use [FoBiS](https://github.com/szaghi/FoBiS).

The root of ADAM framework should look like:

```bash
┌╼ stefano@enlil
├───╼ ~/fortran/adam
└──────╼ tree -L 1
.
├── exe
├── fobos
├── LICENSE
├── README.md
├── scripts
└── src
```

To be completed.

Go to [Top](#top)

### Install NVIDIA HPC SDK

Download the latest NVIDIA HPC SDK (at time of writing is 12.8):

```bash
wget https://developer.download.nvidia.com/hpc-sdk/25.3/nvhpc_2025_253_Linux_x86_64_cuda_12.8.tar.gz
```

Decompress the downloaded archive and install the SDK using the provided script:

```bash
tar xf nvhpc_2025_253_Linux_x86_64_cuda_12.8.tar.gz
cd nvhpc_2025_253_Linux_x86_64_cuda_12.8/
sudo ./install
```

The default installation path is in `/opt/nvidia/hpc_sdk/Linux_x86_64` with a subdirectory corresponding
to the version installed, e.g. `/opt/nvidia/hpc_sdk/Linux_x86_64/25.3`.

After the installation the user PATH must be updated. If module app is available in the system
it can be used to load the SDK on demand, i.e.:

```bash
module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc/25.3
module load nvhpc/25.3
```

Alternatively, the shell environment may be initialized to use the HPC SDK.
In bash, sh, or ksh, use these commands:

```bash
export PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/25.3/compilers/bin:$PATH
export MANPATH=/opt/nvidia/hpc_sdk/Linux_x86_64/25.3/compilers/man:$MANPATH
export PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/25.3/comm_libs/mpi/bin:$PATH
```

Go to [Top](#top)

### Install External Dependencies

In the following the instructions to download and compile external dependencies are provided. Keep in mind that
web link can be changed from the time of writing this notes. Moreover, the compilations instructions are for
NVIDIA SDK, for other compilers may differ.

#### Install HDF5 library

ADAM exploits the HDF5 for performing HPC IO operations, e.g. saving output files or loading restart files.
HDF5, in turn, depends (optionally) on ZLIB and SZIP for compressing output files, thus it is strongly suggested
to install also those libraries.

##### Install SZIP and ZLIB libraries

The SZIP library can be downloaded from a github hosted fork of the latest version:

```bash
git clone https://github.com/erdc/szip
```

Once cloned it can be compiled by:

```bash
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
Note that the directory of installation must be created before `make install` and the user must have the proper permissions.

The ZLIB library can be downloaded by the permalink of the latest version:

```bash
wget https://zlib.net/current/zlib.tar.gz
tar xf zlib.tar.gz
```

Once downloaded and decompressed it can be compiled by:

```bash
cd zlib
export CC=nvc
export CFLAGS='-O3 -fPIC'
./configure --prefix=/opt/zlib/bin/1.3.1/nvf/25.3
make
make check
make install
cd ..
```
Note that the directory of installation must be created before `make install` and the user must have the proper permissions.

##### Install HDF5 library

Once ZLIB and SZIP libraries are installed the HDF5 can be installed too.

To obtain the it you can use the github repository:

```bash
git clone https://github.com/HDFGroup/hdf5
```

Once cloned it can be compiled and installed by:

```bash
cd hdf5
CFLAGS="-fPIC" FCFLAGS="-fPIC" ./configure --prefix=/opt/HDF5/bin/1.14.6/nvf/25.3 --enable-shared --enable-parallel --enable-fortran --disable-libtool-lock --with-szlib=/opt/szip/bin/2.1.1/nvf/25.3 --with-zlib=/opt/zlib/bin/1.3.1/nvf/25.3 CC=mpicc FC=mpif90
make
make check
make install
cd ..
```
Note that the directory of installation must be created before `make install` and the user must have the proper permissions.

Go to [Top](#top)
