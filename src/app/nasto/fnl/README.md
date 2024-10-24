<a name="top"></a>

# NASTO NFL

> ADAM for Navier-Stokes equations, FUNDAL-based (NFL) backend.

NASTO is an application developed on top of ADAM framework to solve compressible Navier-Stokes conservation equations. These are the sources
of device accelerated backend by means of FUNDAL library for device offloading parallelization.

The main NASTO NVF documentation is contained in the following sections:

| [Main Features](#main-features) | [Copyrights](#copyrights) | [Install](#install) | [Test](#test) | [API Documentation](#api-documentation) |

Go to [Top](#top)

# Main Features

NASTO FNL provides the main device accelerated backend (objects and kernels) with the executable program for performing a NASTO simulation on a
device accelerated architecture.

Go to [Top](#top)

# Copyrights

ADAM is currently a closed project:

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
├── fnl/
├── gmp/
├── nvf/
└── README.md

┌╼ stefano@enlil
├───╼ ~/fortran/adam
└──────╼ tree src/app/nasto/fnl/
src/app/nasto/fnl/
├── adam_nasto_fnl_cns_kernels.F90
├── adam_nasto_fnl_kernels.F90
├── adam_nasto_fnl_object.F90
├── adam_nasto_fnl.F90
└── README.md
```

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
└──────╼ FoBiS.py build -mode nasto-fnl-mpi
Builder options
  Directories
    Building directory: "exe"
    Compiled-objects .o   directory: "exe/obj"
    Compiled-objects .mod directory: "exe/mod"
...
Linking exe/adam_nasto_fnl
Target src/app/nasto/fnl/adam_nasto_fnl.F90 has been successfully built
```

NASTO executable is now present into subdirectory `exe`

```bash
┌╼ stefano@enlil
├───╼ ~/fortran/adam
└──────╼ ls exe/
adam_nasto_fnl  build_adam_nasto_fnl.log  mod  obj
```

Note that the `exe` subdirectory contains also the compiled objects, i.e. files into the subrdirectories
`exe/mod` and `exe/obj`, as well as the compilation log file.

Go to [Top](#top)

# Test

NASTO FNL tests are contained into the following subdirectory:

```bash
┌╼ stefano@enlil
├───╼ ~/fortran/adam/src/tests/nasto/fnl
└──────╼ tree
.
├── i-vortex
└── shock-sphere
```

+ `i-vortes` contains the isentropic vortex test, see [its own documentation](../../../tests/nasto/fnl/i-vortex/README.md) for more details;
+ `shock-sphere` contains the shock-sphere interaction test, see [its own documentation](../../../tests/nasto/fnl/shock-sphere/README.md) for more details.

Go to [Top](#top)

# API Documentation

Currently, NASTO FNL app is made by the following source files:

### Standalone objects

+ `adam_nasto_fnl_object.F90 `

### Modules

+ `adam_nasto_fnl_kernels.F90 `

### Main programs

+ `adam_nasto_fnl.F90 `

### Standalone objects

#### `adam_nasto_fnl_object.F90`

This contains the definition of NASTO FNL object class that extends the NASTO common object class with all the necessary data and methods for the FNL backend.
See [nasto FNL object API documentantion](https://szaghi.github.io/adam/type/nasto_fnl_object.html) for more details.

### Modules

#### `adam_nasto_fnl_kernels.F90`

This contains the definition of all NASTO FNL kernels procedures, defined without the explicit knownledge of of any ADAM-NASTO derived types:
all arguments passed to NASTO FNL kernels must be primitive types.
See [nasto FNL kernels API documentantion](https://szaghi.github.io/adam/module/nasto_fnl_kernels.html) for more details.

### Main programs

#### `adam_nasto_fnl.F90`

This is only the main program that instantiates a `type(nasto_fnl_object)` object and invoke its `simulate` method.
See [nasto FNL program API documentantion](https://szaghi.github.io/adam/program/nasto_fnl.html) for more details.

Go to [Top](#top)
