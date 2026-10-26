---
title: About ADAM
---

<img src="/adam-logo-bga.png" alt="ADAM logo" width="480" class="light-only">
<img src="/adam-logo-bgai.png" alt="ADAM logo" width="480" class="dark-only">

# About ADAM

**ADAM** (Accelerated fluid Dynamics on Adaptive Mesh refinement grids) is a high-performance fluid dynamics framework written in modern Fortran for device(GPU)-accelerated HPC simulations and traditional CPU-based superpc. It provides physics-agnostic SDK library offering AMR, IB, WENO, RK, parellel I/O, ecc... on top of which is easy to develop CFD solvers. The maintained applications are PRISM, an electromagnetic and plasma solver, and FLUME, a compressible Euler and MHD solver.

The primary focus is on large-scale scientific computing: ADAM targets exascale deployments on device-accelerated clusters, where thousands of devices (GPU) nodes work in concert via MPI with node-level parallelism provided by OpenACC (NVIDIA GPUs) or OpenMP offloading (AMD GPUs).

## Key capabilities

- **AMR** — Morton-order linearized octree/quadtree; refinement at initialization by geometric, gradient and total-variation markers, and runtime regridding in FLUME (`[amr] frequency`)
- **Immersed Boundary** — complex and moving geometries without body-fitted meshes
- **High-order WENO** — upwind reconstructions of order 1 to 9 and centred of order 2 to 8
- **GPU backend** — one FNL codebase compiles for OpenACC (nvfortran) and OpenMP offloading (amdflang)
- **Elliptic solvers** — FLAIL: FFT and AMR multigrid Poisson solvers, on CPU and GPU
- **I/O** — per-rank HDF5 files with an XDMF index, stream restart files, MPI-IO slices

## Applications

| Application | Description |
|-------------|-------------|
| [FLUME](/applications/flume/) | Compressible Euler and MHD solver |
| [PRISM](/applications/prism/) | Maxwell equations solver for electromagnetic and plasma (PIC) simulations |

NASTO, CHASE, PATCH and ASCOT are in the source tree but unmaintained, and not documented on this site.

## Next

- [Install](/start/install): compilers, dependencies, the varset of your machine, the first build.
- [Concepts](/concepts/architecture): the realm, the forest and the library objects.
- [Citing ADAM](/project/citing), the authors and the [license](/project/license).
