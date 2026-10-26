# FLUME — CPU Backend

> FLUME is part of the ADAM framework, released under the
> [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.en.html) (LGPLv3).
> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.

The CPU backend (`flume-cpu`) is the MPI (+ optional OpenMP) implementation of FLUME, extending the backend-independent
infrastructure of `common/`.

> **Status: development** (milestones M1 to M13 complete, the list in `docs/applications/flume/index.md`).

## Source Files

| File | Description |
|------|-------------|
| `adam_flume_cpu.F90` | Entry point (program) |
| `adam_flume_cpu_object.F90` | `flume_cpu_object`: ghost update and boundary conditions, WENO residual (with the immersed-boundary variant), Runge-Kutta integration, reflux accumulation and application, output |
| `adam_flume_cpu_<model>[_<solver>]_kernels.F90` | One module per physical model (`euler`, `euler_de`, `euler_ps`, `euler_de_ps`, `euler_iso`, `euler_baro`, `mhd`, `mhd_glm`, `mhd_eglm`, `mhd_iso`, `mhd_glm_iso`, `mhd_baro`, `mhd_glm_baro`) and per Riemann solver of `weno-riemann` (`llf`, `hll`, `hllc` for Euler, `hlld` for MHD): the face fluxes of that model, so no kernel branches on the model |
| `adam_flume_cpu_*_kernels_agnostic.INC` | The model-agnostic kernel bodies (face, Riemann face, ghost, auxiliary, positivity, MHD) the model modules include |
| `adam_flume_cpu_dissipation_kernels.F90` | Viscous, heat-conduction and Ohmic (MHD) face fluxes of order 2 and 4 (issue #65), model-agnostic, added to the inviscid face fluxes before the seam accumulation |

## Build modes

```bash
fobis build --mode flume-cpu-gnu          # MPI
fobis build --mode flume-cpu-gnu-omp      # MPI + OpenMP
fobis build --mode flume-cpu-gnu-debug    # MPI, runtime checks
```
