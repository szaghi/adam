<a name="top"></a>

# NASTO CPU

> ADAM for Navier-Stokes equations, CPU backend.

NASTO is an application developed on top of ADAM framework to solve compressible Navier-Stokes conservation equations. These are the sources
of pure CPU backend, with MPI/OpenMP only CPU parallelization.

The main NASTO CPU documentation is contained in the following sections:

| [Main Features](#main-features) | [Copyrights](#copyrights) | [Install](#install) | [Test](#test) | [API Documentation](#api-documentation) |

Go to [Top](#top)

# Main Features

NASTO CPU provides the main CPU backend with the executable program for performing a NASTO simulation on a pure CPU architecture by means of MPI/OpenMP
parallel pardigms.

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

To compile NASTO CPU application the preferred method is to use [FoBiS](https://github.com/szaghi/FoBiS).

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
└──────╼ tree src/app/nasto/cpu/
src/app/nasto/cpu/
├── adam_nasto.ini
├── adam_nasto_cpu.F90
├── adam_nasto_cpu_object.F90
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

#### Compile NASTO CPU

After the third-part libraries are compiled, you can compule NASTO by

```bash
┌╼ stefano@enlil
├───╼ ~/fortran/adam
└──────╼ FoBiS.py build -mode nasto-cpu-mpi-gnu
Builder options
  Directories
    Building directory: "exe"
    Compiled-objects .o   directory: "exe/obj"
    Compiled-objects .mod directory: "exe/mod"
  External libraries directories: /opt/HDF5/bin/1.12.2/openmpi/4.1.4/gnu/11.2.0/lib
  Included paths: /opt/HDF5/bin/1.12.2/openmpi/4.1.4/gnu/11.2.0/include src/third_party/VecFor/src/lib
  Linked libraries with full path: /opt/HDF5/bin/1.12.2/openmpi/4.1.4/gnu/11.2.0/lib/libhdf5hl_fortran.a /opt/HDF5/bin/1.12.2/openmpi/4.1.4/gnu/11.2.0/lib/libhdf5_hl.a /opt/HDF5/bin/1.12.2/openmpi/4.1.4/gnu/11.2.0/lib/libhdf5_fortran.a /opt/HDF5/bin/1.12.2/openmpi/4.1.4/gnu/11.2.0/lib/libhdf5.a /opt/zlib/bin/1.2.11/lib/libz.a /opt/szip/bin/2.1.1/lib/libsz.a ./exe/obj/getmemory.o
  Linked libraries in path: dl m
  Compiler options
    Vendor: "gnu"
    Compiler command: "mpif90"
    Module directory switch: "-J"
    Compiling flags: "-cpp -c -O2 -D_MPI_"
    Linking flags: "-O2"
    Preprocessing flags: "-D_MPI_"
    Coverage: False
    Profile: False
  Preprocessor used: None
  Preprocessor output directory: None
  Preprocessor extensions processed: []

Building src/app/nasto/cpu/adam_nasto_cpu.F90
Compiling src/third_party/PENF/src/lib/penf_global_parameters_variables.F90 serially
Compiling src/third_party/PENF/src/lib/penf_b_size.F90 serially
Compiling src/third_party/PENF/src/lib/penf_stringify.F90 serially
Compiling src/third_party/PENF/src/lib/penf.F90 serially
...
...
...
Linking exe/adam_nasto_cpu
Target src/app/nasto/cpu/adam_nasto_cpu.F90 has been successfully built
```

NASTO executable is now present into subdirectory `exe`

```bash
┌╼ stefano@enlil
├───╼ ~/fortran/adam
└──────╼ ls exe/
adam_nasto_cpu  build_adam_nasto_cpu.log  mod  obj
```

Note that the `exe` subdirectory contains also the compiled objects, i.e. files into the subrdirectories
`exe/mod` and `exe/obj`, as well as the compilation log file.

Go to [Top](#top)

# Test

To be written.

Go to [Top](#top)

# API Documentation

Currently, NASTO CPU objects are made by the following source files:

### Standalone objects

+ `adam_nasto_cpu_object.F90 `

### Main programs

+ `adam_nasto_cpu.F90 `

### Standalone objects

#### `adam_nasto_cpu_object.F90`

This contains the definition of NASTO CPU object that extends the NASTO common object with all the necessary data and methods for the CPU backend.
See [nasto CPU object API documentantion](https://szaghi.github.io/adam/type/nasto_cpu_object.html) for more details.

### Main programs

#### `adam_nasto_cpu.F90`

This is only the main program that instantiates a `type(nasto_cpu_object)` object and invoke its `simulate` method.
See [nasto CPU program API documentantion](https://szaghi.github.io/adam/program/nasto_cpu.html) for more details.

Go to [Top](#top)
