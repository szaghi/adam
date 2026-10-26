# FLUME — CPU Backend

The CPU backend (`flume-cpu`) is the MPI (+ optional OpenMP) implementation of FLUME, extending the backend-independent
infrastructure of `common/`.

**Documentation:** <https://szaghi.github.io/adam/applications/flume/cpu>

## Source Files

| File | Description |
|------|-------------|
| `adam_flume_cpu.F90` | Entry point (program) |
| `adam_flume_cpu_object.F90` | `flume_cpu_object`: ghost update and boundary conditions, WENO residual (with the immersed-boundary variant), Runge-Kutta integration, reflux accumulation and application, output |
| `adam_flume_cpu_<model>[_<solver>]_kernels.F90` | One module per physical model (`euler`, `euler_de`, `euler_ps`, `euler_de_ps`, `euler_iso`, `euler_baro`, `mhd`, `mhd_glm`, `mhd_eglm`, `mhd_iso`, `mhd_glm_iso`, `mhd_baro`, `mhd_glm_baro`) and per Riemann solver of `weno-riemann` (`llf`, `hll`, `hllc` for Euler, `hlld` for MHD): the face fluxes of that model, so no kernel branches on the model |
| `adam_flume_cpu_*_kernels_agnostic.INC` | The model-agnostic kernel bodies (face, Riemann face, ghost, auxiliary, positivity, MHD) the model modules include |
| `adam_flume_cpu_dissipation_kernels.F90` | Viscous, heat-conduction and Ohmic (MHD) face fluxes of order 2 and 4 (issue #65), model-agnostic, added to the inviscid face fluxes before the seam accumulation |

ADAM is dual-licensed under MIT or LGPLv3 (`LICENSE.mit.md`, `LICENSE.lgpl3.md`).
