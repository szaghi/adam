---
title: About ADAM
---

# About ADAM

**ADAM** (Accelerated fluid Dynamics on Adaptive Mesh refinement grids) is a high-performance fluid dynamics framework written in modern Fortran for device(GPU)-accelerated HPC simulations and traditional CPU-based superpc. It provides physics-agnostic SDK library offering AMR, IB, WENO, RK, parellel I/O, ecc... on top of which is easy to develop CFD solvers (currently a compressible Navier-Stokes solver named NASTO and an electromagnetic plasma solver named PRSIM are the principal applications developed with ADAM).

The primary focus is on large-scale scientific computing: ADAM targets exascale deployments on device-accelerated clusters, where thousands of devices (GPU) nodes work in concert via MPI with node-level parallelism provided by CUDA Fortran, OpenACC, or OpenMP offloading.

## Key capabilities

- **AMR** — Morton-order linearized octree/quadtree with automatic refinement and coarsening
- **Immersed Boundary** — complex and moving geometries without body-fitted meshes
- **High-order WENO** — finite difference reconstructions up to 5th order for shock-capturing flows
- **Multi-backend GPU** — a single codebase compiles for CUDA Fortran (NVF), OpenACC (FNL), and OpenMP offloading (GMP)
- **Parallel I/O** — HDF5 output and restart files via MPI-IO

## Applications

| Application | Description |
|-------------|-------------|
| [NASTO](/applications/nasto/) | Compressible 3D Navier-Stokes solver |
| [PRISM](/applications/prism) | Maxwell equations solver for electromagnetic simulations |
| [CHASE](/applications/chase) | CFD application |
| [PATCH](/applications/patch) | Patch-based application |
| [ASCOT](/applications/ascot) | Binary-to-ASCII output converter |

## Authors

- Andrea di Mascio — [andrea.dimascio@univaq.it](mailto:andrea.dimascio@univaq.it)
- Federico Negro — [federico.negro.01@gmail.com](mailto:federico.negro.01@gmail.com)
- Giacomo Rossi — [giacomo.rossi@amd.com](mailto:giacomo.rossi@amd.com)
- Francesco Salvadore — [f.salvadore@cineca.it](mailto:f.salvadore@cineca.it)
- Stefano Zaghi — [stefano.zaghi@cnr.it](mailto:stefano.zaghi@cnr.it)

## Copyrights

ADAM is dual-licensed under the [MIT License](/license#mit) and the [GNU Lesser General Public License v3.0](/license#lgplv3) (LGPLv3). You may choose either license.

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.

## Citing ADAM

If you use ADAM in work that leads to a scientific publication, please cite the following paper:

> S. Zaghi, F. Salvadore, A. Di Mascio, G. Rossi —
> **Efficient GPU parallelization of adaptive mesh refinement technique for high-order compressible solver with immersed boundary** —
> *Computers and Fluids*, 266 (2023) 106040.
> DOI: [10.1016/j.compfluid.2023.106040](https://doi.org/10.1016/j.compfluid.2023.106040)

The paper describes the ADAM framework architecture, the AMR/IB coupling strategy, the GPU parallelization approach (CUDA Fortran), and demonstrates strong scaling on a shock–sphere interaction benchmark. A preprint is available in [docs/papers/zaghi-2023-computer_fluids.pdf](../../papers/zaghi-2023-computer_fluids.pdf).

BibTeX entry:

```bibtex
@article{zaghi2023adam,
  author  = {Zaghi, S. and Salvadore, F. and {Di Mascio}, A. and Rossi, G.},
  title   = {Efficient {GPU} parallelization of adaptive mesh refinement technique
             for high-order compressible solver with immersed boundary},
  journal = {Computers \& Fluids},
  volume  = {266},
  pages   = {106040},
  year    = {2023},
  doi     = {10.1016/j.compfluid.2023.106040},
}
```
