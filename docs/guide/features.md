---
title: Features
---

# Features

## GPU backends

The maintained GPU backend is **FNL**, one source tree with two device paths:

| Modes | Macro | Compiler | Parallelism model |
|-------|-------|----------|-------------------|
| `*-fnl-nvf` | `_FNL` | NVIDIA nvfortran | OpenACC + FUNDAL memory management (NVIDIA GPUs) |
| `*-fnl-omp-amd` | `_FNL`, `DEV_OMP` | AMD amdflang | OpenMP target offloading + FUNDAL (AMD GPUs) |

The kernels carry both directive sets (`!$acc` and `!$omp` twins). The `src/lib/nvf` (CUDA Fortran) and `src/lib/gmp` (OpenMP offloading) libraries are in the tree, but no PRISM or FLUME build mode uses them.

Each application (`prism`, `flume`) compiles to separate executables per backend, sharing a common source tree under `src/app/<name>/common/`.

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
| `weno-u-1` | 1st (upwind) | `ngc = 1` |
| `weno-u-3` | 3rd (upwind) | `ngc = 2` |
| `weno-u-5` | 5th (upwind) | `ngc = 3` |
| `weno-u-7` | 7th (upwind) | `ngc = 4` |
| `weno-u-9` | 9th (upwind) | `ngc = 5` |
| `weno-c-2` … `weno-c-8` | 2nd–8th (centred) | `ngc = order/2` |

Nonlinear weights: Jiang–Shu (`js`, the default) or scale-invariant (`si`). Reconstruction variables: conservative or characteristic. ROR iterations lower the stencil where a variable check fails, and the order is reduced close to immersed solids.

---

## Temporal Integration

Explicit Runge-Kutta multi-stage methods:

- `runge-kutta-1`, `-2`, `-3`; SSP `runge-kutta-ssp-11`, `-22`, `-33`, `-54`; the 4-stage symplectic `runge-kutta-yoshida`
- CFL-based time step control
- Leapfrog, Blanes–Moan (`blanes-moan4`) and commutator-free Magnus (`CFM-4`) integrator objects

---

## Field Memory Layout

5D field arrays `(nv, ni, nj, nk, nb)` with Fortran column-major storage:

- `nv` — number of variables (density, momentum, energy, …)
- `ni, nj, nk` — grid cells per block (including ghost cells)
- `nb` — number of blocks (AMR + domain decomposition)

The leftmost index `nv` is contiguous in memory and is the innermost loop for CPU vectorization and GPU coalesced access.

---

## Parallel I/O

- Field output: one HDF5 file per MPI rank (`<basename>-procNNNNNN.h5`) plus one XDMF file (`<basename>.xdmf`) for ParaView/VisIt
- Restart: Fortran stream files, the blocks per rank (`-procNNNNNN.fbd`), the tree nodes (`.tnd`) and the time (`.time`)
- Slices: raw binary written with MPI-IO (`.mat`), converted by ASCOT or `slices_converter.py`

---

## Compiler Support

| Compiler | GPU backends | CPU |
|----------|-------------|-----|
| NVIDIA nvfortran 26.9 (tested) | FNL (OpenACC) | Yes |
| GNU gfortran 16.2.0 + Open MPI 5.0.11 (tested) | — | Yes |
| AMD amdflang | FNL (OpenMP offloading, `*-fnl-omp-amd`) | Yes (`prism-cpu-amd`) |

---

## Design Principles

- **KISS** — keep it simple; prefer explicit Fortran over macro magic
- **Portable** — a single source tree, per-backend preprocessor guards only where unavoidable
- **Evidence-based** — profile first, optimize with benchmarks, document trade-offs
- **Standard-conforming** — Fortran 2003+ OOP; `implicit none` everywhere; explicit `intent` on all arguments
