<a name="top"></a>

# NASTO

> ADAM for Navier-Stokes equations.

NASTO is an application developed on top of ADAM framework to solve compressible Navier-Stokes conservation equations.

The main NASTO documentation is contained in the following sections:

| [Main Features](#main-features) | [Copyrights](#copyrights) | [Install](#install) | [Test](#test) | [API Documentation](#api-documentation) |

Go to [Top](#top)

# Main Features

The main features of NASTO application are:

+ the mathematical model is the compressible, ideal gas 3D Navier-Stokes conservation laws;
+ finite difference numerical approximation is used for the spatial discretization:
    + high order Weighted Essentially Non Oscillatory (WENO) reconstructions (upwind schemes up to 5th order) are used for achieving high order accuracy;
    + accurate Immersed Boundary (IB) method is used to discretize solid body;
+ explicit Runge-Kutta multi-stages methods are used for the temporal discretization;
+ 4 different natural boundary conditions are currently implemented:
    + supersonic inflow;
    + pure extrapolation outflow;
    + solid wall;
+ 1 numerical boundary conditions is currently implemented:
    + periodic boundary;
+ 2 different initial conditions are currently implemented:
    + 2 regions initial conditions for shock-sphere interaction simulation;
    + 2 regions initial conditions for vortex advection simulation;
+ simulation setup by means of a single input `ini` file.

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

To compile NASTO application the preferred method is to use [FoBiS](https://github.com/szaghi/FoBiS).

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
└── README.md
```

Each subdirectory, contains a specific sources type:

+ `common` contains sources that are in common to all backends;
+ `cpu` contains sources specific for CPU only backend, namely MPI/OpneMP only;
+ `gmp` contains sources specific for GPU accelerated backend by means of OpenMP offloading and MPI communications;
+ `nvf` contains sources specific for GPU accelerated backend by means of Nvidia CUDAFortran offloading and MPI communications.

Typically, a specific backend uses its own sources and the common ones. Currently, NASTO has only the `nvf` backend completed. The
development of `gmp` backend is ongoing.

To compile NASTO use FoBiS tool as in the examples contained into each backend sources subdirectory.

Go to [Top](#top)

# Test

To be written.

Go to [Top](#top)

# API Documentation

Refer to the documentation contained into each backend sources subdirectory:

+ [`common`](./nasto/common/README.md)
+ [`cpu`](./nasto/cpu/README.md)
+ [`gmp`](./nasto/gmp/README.md)
+ [`nvf`](./nasto/nvf/README.md)

Go to [Top](#top)
