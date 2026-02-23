---
title: Features
---

# Features

## GPU backends

ADAM supports three GPU-acceleration backends, all sharing the same numerical core:

| Backend | Macro | Compiler | Parallelism model |
|---------|-------|----------|-------------------|
| **NVF** | `_NVF` | NVIDIA nvfortran | CUDA Fortran device arrays and kernels |
| **FNL** | `_FNL` | NVIDIA nvfortran | OpenACC + FUNDAL memory management |
| **GMP** | `_GMP` | Intel ifx / GNU gfortran | OpenMP target offloading (experimental) |

Each application (`nasto`, `prism`, …) compiles to separate executables per backend, sharing a common source tree under `src/app/<name>/common/`.

---

## Adaptive Mesh Refinement

ADAM uses a **Morton-order linearized octree/quadtree** for AMR:

- Blocks are sorted by Z-order (Morton) curve for cache-friendly traversal
- Refinement and coarsening preserve conservation properties
- Ghost cell exchange via MPI before every stencil operation
- Configurable ghost cell width `ngc` (WENO5: 3, WENO7: 4)

---

## Immersed Boundary Method

Complex and moving geometries are handled without body-fitted meshes:

- Solid body positions are tracked on a background Cartesian grid
- IB forcing reconstructs boundary conditions at cut cells
- Eikonal-based signed distance field for geometry representation
- Compatible with AMR — IB markers trigger local refinement

---

## WENO Schemes

Weighted Essentially Non-Oscillatory finite difference reconstructions:

| Scheme | Order | Ghost cells |
|--------|-------|-------------|
| WENO-3 | 3rd | `ngc = 2` |
| WENO-5 | 5th | `ngc = 3` |
| WENO-7 | 7th | `ngc = 4` |
| WENO-9 | 9th | `ngc = 5` |
| WENO-11 | 11th | `ngc = 6` |

Reconstruction variables: conservative or characteristic. Order-of-refinement (ROR) switching near shocks is supported.

---

## Temporal Integration

Explicit Runge-Kutta multi-stage methods:

- SSP-RK3, SSP-RK4, SSP-RK5(4), and standard RK schemes
- CFL-based time step control
- Optional leapfrog and Adams-Bashforth integrators

---

## Field Memory Layout

5D field arrays `(nv, ni, nj, nk, nb)` with Fortran column-major storage:

- `nv` — number of variables (density, momentum, energy, …)
- `ni, nj, nk` — grid cells per block (including ghost cells)
- `nb` — number of blocks (AMR + domain decomposition)

The leftmost index `nv` is contiguous in memory and is the innermost loop for CPU vectorization and GPU coalesced access.

---

## Parallel I/O

- HDF5 parallel output with ZLIB/SZIP compression
- Restart files for checkpoint/resume
- Slice extraction for in-situ visualization
- XDMF metadata support for ParaView/VisIt

---

## Compiler Support

| Compiler | GPU backends | CPU |
|----------|-------------|-----|
| NVIDIA nvfortran ≥ 25.3 | NVF, FNL | Yes |
| GNU gfortran ≥ 9 | — | Yes |
| Intel ifort / ifx | GMP (experimental) | Yes |
| AMD flang | — | Planned |

---

## Design Principles

- **KISS** — keep it simple; prefer explicit Fortran over macro magic
- **Portable** — a single source tree, per-backend preprocessor guards only where unavoidable
- **Evidence-based** — profile first, optimize with benchmarks, document trade-offs
- **Standard-conforming** — Fortran 2003+ OOP; `implicit none` everywhere; explicit `intent` on all arguments
