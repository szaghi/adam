# FLUME — FNL (OpenACC) Backend

The FNL backend (`flume-fnl`) is the GPU implementation of FLUME on the **FNL library** (OpenACC directives, each with
its OpenMP offload twin), extending the backend-independent infrastructure of `common/`. The fields live on the device;
the host copies are refreshed for output, restart and the coarse-fine seam fluxes of the reflux register.

> **Status: development** (milestones M1 to M13 complete, the list in `docs/applications/flume/index.md`). The seam
> face fluxes of a rank go to the host register in one copy per stage (issue #90).

## Source Files

| File | Description |
|------|-------------|
| `adam_flume_fnl.F90` | Entry point (program) |
| `adam_flume_fnl_object.F90` | `flume_fnl_object`: device data, ghost update, residual, Runge-Kutta integration, reflux, output |
| `adam_flume_fnl_dissipation_kernels.F90` | Device twins of the CPU viscous, heat and Ohmic face-flux kernels (issue #65) |
| `adam_flume_fnl_kernels.F90` | Device kernels: face fluxes and flux difference (with the immersed-boundary variant), boundary conditions, time-step and conservation reductions, auxiliary fields, seam skin packing and reflux application, radiative cooling and the cooling time (issue #93) |
| `adam_flume_fnl_<model>[_<solver>]_kernels.F90` | Device twins of the CPU per-model and per-Riemann-solver kernel modules (same model and solver names) |
| `adam_flume_fnl_*_kernels_agnostic.INC` | The model-agnostic device kernel bodies the model modules include |
| `adam_flume_fnl_library.F90` | Barrel re-export of all FNL modules |

## Build modes

```bash
fobis build --mode flume-fnl-nvf --varset local_nvf          # NVIDIA, OpenACC (verified)
fobis build --mode flume-fnl-nvf-debug --varset local_nvf    # NVIDIA, OpenACC, runtime checks
fobis build --mode flume-fnl-omp-amd                          # AMD, OpenMP offload (not verified)
```
