# FLUME — FNL (OpenACC) Backend

> FLUME is part of the ADAM framework, released under the
> [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.en.html) (LGPLv3).
> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.

The FNL backend (`flume-fnl`) is the GPU implementation of FLUME on the **FNL library** (OpenACC directives, each with
its OpenMP offload twin), extending the backend-independent infrastructure of `common/`. The fields live on the device;
the host copies are refreshed for output, restart and the coarse-fine seam fluxes of the reflux register.

> **Status: development** (milestone M1, compressible Euler). Every verification test passes on both backends; the
> per-stage host copies of the seam faces dominate the run time of AMR cases.

## Source Files

| File | Description |
|------|-------------|
| `adam_flume_fnl.F90` | Entry point (program) |
| `adam_flume_fnl_object.F90` | `flume_fnl_object`: device data, ghost update, residual, Runge-Kutta integration, reflux, output |
| `adam_flume_fnl_kernels.F90` | Device kernels: face fluxes and flux difference (with the immersed-boundary variant), boundary conditions, time-step and conservation reductions, auxiliary fields, seam skin packing and reflux application |
| `adam_flume_fnl_library.F90` | Barrel re-export of all FNL modules |

## Build modes

```bash
fobis build --mode flume-fnl-nvf --varset local_nvf          # NVIDIA, OpenACC (verified)
fobis build --mode flume-fnl-nvf-debug --varset local_nvf    # NVIDIA, OpenACC, runtime checks
fobis build --mode flume-fnl-omp-amd                          # AMD, OpenMP offload (not verified)
```
