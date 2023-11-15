# ADAM

> ADAM, Adaptive Mesh Refinement (AMR) with Immersed Boundary (IB) fluid dynamic solver tailored for High Performance GPU Computing

### Authors

+ Francesco Salvadore, [f.salvadore@cineca.it](mailto:f.salvadore@cineca.it)
+ Andrea di Mascio, [andrea.dimascio@univaq.it](andrea.dimascio@univaq.it)
+ Giacomo Rossi, [giacomo.rossi@intel.com](giacomo.rossi@intel.com)
+ Stefano Zaghi, [stefano.zaghi@cnr.it](stefano.zaghi@cnr.it)

### Features

+ KISS, keep it simple and stupid;
+ simulate different types of flows;
+ use high order finite difference scheme;
+ exploit Adaptive Mesh Refinement (AMR) to accurate simulate complex geometries;
+ exploit Immersed Boundary (IB) techniques to easy handle complex moving geometries;
+ ready for very High Performance Computing (HPC) by means of GPU parallel paradigm exploitation;

## Contributing

The project home is hosted on GitHub and the development is organized by means of multiple branches as described in the following.

### Branches Organization

Currently, ADAM development tree contains the following branches:

+ master
+ develop
+ fortran-cuda-stable
+ fortran-cuda-develop
+ openmp-gpu-develop
+ euler-cpu-develop

#### master

Soon or later a stable v1.0.0 will born and `master` branch will be its home, for now it is an empty branch.

#### develop

It has been the home for the Fortran-CUDA development until now. In 15th September 2023 it has been copied into branch `fortran-cuda-stable` and it is currently freeze.

#### fortran-cuda-stable

It is a **reference** branch, used for comparison reason, a sort of *named commit*.

#### fortran-cuda-develop

It contains the ongoing Fortran-CUDA development.

Currently, the main CUDA development is focused into the NASTO (CUDA) application. More details about NASTO app can be found in its [readme](src/app/nasto/README.md).

#### openmp-gpu-develop

It contains the ongoing OpenMP development, GPU offloading.

#### euler-cpu-develop

It contains the ongoing development of a very simplified application for Euler flows simulation with a simple CPU backend.
