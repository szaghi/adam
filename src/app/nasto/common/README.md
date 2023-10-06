<a name="top"></a>

# NASTO (common)

> ADAM for Navier-Stokes equations, common data.

NASTO is an application developed on top of ADAM framework to solve compressible Navier-Stokes conservation equations. These are the sources common
to all backends.

The main NASTO common documentation is contained in the following sections:

| [Main Features](#main-features) | [Copyrights](#copyrights) | [Install](#install) | [API Documentation](#api-documentation) |

Go to [Top](#top)

# Main Features

NASTO common data provides the following features:

+ *standalone* objects for easy handling all NASTO simulation setups that are not peculiar of a specific backend, e.g. the setup of initial conditions,
  boundary conditions, simulation IO configurations, ecc...
    + each standalone object provides helpful methods (type bound procedures) for its own specific aims as well as for easy initialization by parsing
      input config (INI) file;
+ all standalone objects are *collected* by a main object, namely the `nasto_common_object` that is designed to be **extended** by all specific backend
  classes in order to obtain easily the handlers for easy simulation setup.
    + the main common object provides helpful methods for non-backend-specific aims, e.g. IO tasks.

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

The common sources does not provide a standalone program it being designed as library for specific backend, instead compile one of the backend provided,
e.g. the `nvf` backend, following its own documentation.

Go to [Top](#top)

# API Documentation

Currently, NASTO common objects are made by the following source files:

### standalone objects

+ `adam_nasto_bc_object.F90`
+ `adam_nasto_eos_object.F90`
+ `adam_nasto_ic_object.F90`
+ `adam_nasto_io_object.F90`
+ `adam_nasto_physics_object.F90`
+ `adam_nasto_schemes_object.F90`
+ `adam_nasto_time_object.F90`

### modules

+ `adam_nasto_parameters.F90`

### main common object

+ `adam_nasto_common_object.F90 `

### Standalone objects API

#### `adam_nasto_bc_object.F90`

This contains the definition of NASTO boundary conditions handler, e.g. the type of boundary conditions implemented, their parsing from config file,
the allocation of the neccessary data, ecc... See [nasto BC object API documentantion](https://szaghi.github.io/adam/type/nasto_bc_object.html)
for more details.

#### `adam_nasto_eos_object.F90`

This contains the definition of NASTO equations of state handler, e.g. the relation between fluids state variables, fluids constants, their parsing from
config file, the allocation of the neccessary data, ecc... It is designed to handle multi-fluids simulations.
See [nasto EOS object API documentantion](https://szaghi.github.io/adam/type/nasto_eos_object.html) for more details.

#### `adam_nasto_ic_object.F90`

This contains the definition of NASTO intial conditions handler, e.g. the type of intial conditions implemented, their parsing from config file,
the allocation of the neccessary data, ecc... See [nasto IC object API documentantion](https://szaghi.github.io/adam/type/nasto_ic_object.html)
for more details.

#### `adam_nasto_io_object.F90`

This contains the definition of NASTO Input/Output (IO) handler, e.g. the frequency and location of IO, their parsing from config file,
the allocation of the neccessary data, ecc... See [nasto IO object API documentantion](https://szaghi.github.io/adam/type/nasto_io_object.html)
for more details.

#### `adam_nasto_physics_object.F90`

This contains the definition of NASTO fluids physics, e.g. the fluid state variables transformations, their parsing from
config file, the allocation of the neccessary data, ecc... It is designed to handle multi-fluids simulations.
See [nasto PHYSICS object API documentantion](https://szaghi.github.io/adam/type/nasto_physics_object.html) for more details.

#### `adam_nasto_schemes_object.F90`

This contains the definition of NASTO numerical schemes handler, e.g. the numerical time and spatial operators, their parsing from
config file, the allocation of the neccessary data, ecc...
See [nasto SCHEMES object API documentantion](https://szaghi.github.io/adam/type/nasto_SCHEMES_object.html) for more details.

#### `adam_nasto_time_object.F90`

This contains the definition of NASTO time handler, e.g. the simulation time-related setup ()maximum number of iteration, CFL limit, ec...), their parsing from
config file, the allocation of the neccessary data, ecc...
See [nasto TIME object API documentantion](https://szaghi.github.io/adam/type/nasto_time_object.html) for more details.

### Modules API

#### `adam_nasto_parameters.F90`

This contains the main global parameters of NASTO.
See [nasto parameters module API documentantion](https://szaghi.github.io/adam/module/nasto_parameters.html) for more details.

### Main common object API

#### `adam_nasto_common_object.F90`

This contains the definition of NASTO common class, a class that is used (extended) by all NASTO backends.
See [nasto common object API documentantion](https://szaghi.github.io/adam/type/nasto_common_object.html) for more details.

Go to [Top](#top)
