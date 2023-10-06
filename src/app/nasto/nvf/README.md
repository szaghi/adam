<a name="top"></a>

# NASTO NVF

> ADAM for Navier-Stokes equations, GPU-Nvidia CUDAFortran (NVF) backend.

NASTO is an application developed on top of ADAM framework to solve compressible Navier-Stokes conservation equations. These are the sources
of GPU accelerated backend, with MPI/Nvidia CUDAFortran GPU offloading parallelization.

The main NASTO NVF documentation is contained in the following sections:

| [Main Features](#main-features) | [Copyrights](#copyrights) | [Install](#install) | [Test](#test) | [API Documentation](#api-documentation) |

Go to [Top](#top)

# Main Features

NASTO GMP provides the main GPU accelerated backend (objects and kernels) with the executable program for performing a NASTO simulation on a
GPU accelerated architecture by means of MPI/Nvidia CUDAFortran GPU offloading parallel pardigms.

Go to [Top](#top)

# Copyrights

ADAM-NASTO is currently a closed project:

> Copyright (C) Di Mascio/Rossi/Salvadore/Zaghi, Inc - All Rights Reserved.
>
> Unauthorized copying of these source files, via any medium is strictly prohibited, proprietary and confidential.
> Written by Andrea di Mascio, Giacomo Rossi, Francesco Salvadore and Stefano Zaghi, September 2023.

Future versions could be released with a more Free Open Source Software (FOSS) licence.

Go to [Top](#top)

# Install

NASTO, like the ADAM framework, is provided as source files archive and it must be compiled in order to have the executable ready to be installed.

### Compile

To compile NASTO NVF application the preferred method is to use [FoBiS](https://github.com/szaghi/FoBiS).

The root of ADAM framework should look like:

```bash
┌╼ stefano@enlil
├───╼ ~/fortran/adam
└──────╼ tree -L 1
.
├── exe
├── fobos -> fobos_nasto
├── fobos_nasto
├── LICENSE
├── README.md
├── scripts
├── src
└── tests
```

The `fobos_nasto` contains the instructions to compile NASTO application by means of FoBiS tool.

The `src` subdirectory contains the sources of ADAM and NASTO, in particular it should look like:

```bash
┌╼ stefano@enlil
├───╼ ~/fortran/adam
└──────╼ tree src/ -L 1
src/
├── app
├── lib
├── tests
├── third_party
└── third_party_manual

┌╼ stefano@enlil
├───╼ ~/fortran/adam
└──────╼ tree src/app/nasto/
src/app/nasto/
├── common/
├── cpu/
├── gmp/
├── nvf/
├── adam_equation_nasto_gpu_object.F90
├── adam_nasto_gpu.F90
├── adam_nasto_sphere_shock.ini
└── README.md

┌╼ stefano@enlil
├───╼ ~/fortran/adam
└──────╼ tree src/app/nasto/nvf/
src/app/nasto/nvf/
├── adam_nasto.ini
├── adam_nasto_nvf.F90
├── adam_nasto_nvf_kernels.F90
├── adam_nasto_nvf_object.F90
└── README.md
```

The `adam_nastro.ini` contains an example of the `ini` input.

To compile NASTO use FoBiS tool as in the following examples.

#### Prepare some third-party libraries

```bash
FoBiS.py rule -ex makethirdpartymanual
```

This should produce something like the following:

```bash
┌╼ stefano@enlil
├───╼ ~/fortran/adam
└──────╼ FoBiS.py rule -ex makethirdpartymanual
Executing rule "makethirdpartymanual"
   Command => mkdir -p exe/mod exe/obj
   Command => cd src/third_party_manual/CGAL/ ; g++ -std=c++17 -frounding-math -O2 -I/opt/cgal/5.2.1/include/ -c cgal_c_wrappers.cpp ; cd -
   Command => cd src/third_party_manual/getmemory/ ; gcc -c getmemory.c ; cd -
   Command => mv src/third_party_manual/getmemory/*.o exe/obj/
   Command => mv src/third_party_manual/CGAL/*.o exe/obj/
```

#### Compile NASTO NVF

After the third-part libraries are compiled, you can compule NASTO by

```bash
FoBiS.py rule -ex makethirdpartymanual
```

This should produce something like the following:

```bash
┌╼ stefano@enlil
├───╼ ~/fortran/adam
└──────╼ FoBiS.py build -mode nasto-nvf-mpi-cuda
Builder options
  Directories
    Building directory: "exe"
    Compiled-objects .o   directory: "exe/obj"
    Compiled-objects .mod directory: "exe/mod"
  External libraries directories: /opt/HDF5/bin/1.12.1/openmpi/3.1.5/nvidia/22.3/lib
  Included paths: /opt/HDF5/bin/1.12.1/openmpi/3.1.5/nvidia/22.3/include src/third_party/VecFor/src/lib
  Linked libraries with full path: /opt/HDF5/bin/1.12.1/openmpi/3.1.5/nvidia/22.3/lib/libhdf5hl_fortran.a /opt/HDF5/bin/1.12.1/openmpi/3.1.5/nvidia/22.3/lib/libhdf5_hl.a /opt/HDF5/bin/1.12.1/openmpi/3.1.5/nvidia/22.3/lib/libhdf5_fortran.a /opt/HDF5/bin/1.12.1/openmpi/3.1.5/nvidia/22.3/lib/libhdf5.a -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/11 /usr/lib/gcc/x86_64-linux-gnu/11/libstdc++.so /opt/zlib/bin/1.2.11/lib/libz.a /opt/szip/bin/2.1.1/lib/libsz.a ./exe/obj/cgal_c_wrappers.o ./exe/obj/getmemory.o
  Linked libraries in path: dl m
  Compiler options
    Vendor: "nvfortran"
    Compiler command: "mpif90"
    Module directory switch: "-module"
    Compiling flags: "-cpp -c -Mcuda=cc75,cuda12.0,ptxinfo -O2 -D_NVF -D_MPI_"
    Linking flags: "-Mcuda=cc75,cuda12.0,ptxinfo -O2"
    Preprocessing flags: "-D_NVF -D_MPI_"
    Coverage: False
    Profile: False
  Preprocessor used: None
  Preprocessor output directory: None
  Preprocessor extensions processed: []

Building src/app/nasto/nvf/adam_nasto_nvf.F90
Compiling src/third_party/PENF/src/lib/penf_global_parameters_variables.F90 serially
ptxas info    : 128 bytes gmem

Compiling src/third_party/PENF/src/lib/penf_b_size.F90 serially
ptxas info    : 4 bytes gmem

Compiling src/third_party/PENF/src/lib/penf_stringify.F90 serially
ptxas info    : 4 bytes gmem
...
...
...
Linking exe/adam_nasto_gpu
Target src/app/nasto/nvf/adam_nasto_nvf.F90 has been successfully built
```

NASTO executable is now present into subdirectory `exe`

```bash
┌╼ stefano@enlil
├───╼ ~/fortran/adam
└──────╼ ls exe/
adam_nasto_nvf  build_adam_nasto_nvf.log  mod  obj
```

Note that the `exe` subdirectory contains also the compiled objects, i.e. files into the subrdirectories
`exe/mod` and `exe/obj`, as well as the compilation log file.

Go to [Top](#top)

# Test

To be written.

Go to [Top](#top)

# API Documentation

Currently, NASTO NVF objects are made by the following source files:

### Standalone objects

+ `adam_nasto_nvf_object.F90 `

### Modules

+ `adam_nasto_nvf_kernels.F90 `

### Main programs

+ `adam_nasto_nvf.F90 `

### Standalone objects

#### `adam_nasto_nvf_object.F90`

This contains the definition of NASTO NVF object that extends the NASTO common object with all the necessary data and methods for the NVF backend.
See [nasto NVF object API documentantion](https://szaghi.github.io/adam/type/nasto_nvf_object.html) for more details.

### Modules

#### `adam_nasto_nvf_kernels.F90`

This contains the definition of all NASTO NVF CUDAFortran kernels procedures, defined without the explicit knownledge of of any ADAM-NASTO derived types:
all arguments passed to NASTO NVF kernels must be primitive types.
See [nasto NVF kernels API documentantion](https://szaghi.github.io/adam/module/nasto_nvf_kernels.html) for more details.

### Main programs

#### `adam_nasto_nvf.F90`

This is only the main program that instantiates a `type(nasto_nvf_object)` object and invoke its `simulate` method.
See [nasto NVF program API documentantion](https://szaghi.github.io/adam/program/nasto_nvf.html) for more details.

Go to [Top](#top)
