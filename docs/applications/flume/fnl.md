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
| `adam_flume_fnl_*_kernel_body.INC` | Per-stencil-width bodies of the face-flux and gravity edge-correction kernels, flux splitting and `weno-riemann` (`face_fluxes`, `gravity_edge`, `riemann_face_fluxes`, `riemann_gravity_edge`), included by the agnostic kernel files |
| `adam_flume_fnl_hybrid_kernels.F90` | Device twins of the CPU `hybrid-kep-weno` kernels (issue #119) |
| `adam_flume_fnl_forcing_kernels.F90` | Device twins of the CPU mass-flow forcing kernels (issue #119) |
| `adam_flume_fnl_library.F90` | Barrel re-export of all FNL modules |

## Build modes

```bash
fobis build --mode flume-fnl-nvf --varset local_nvf          # NVIDIA, OpenACC (verified)
fobis build --mode flume-fnl-nvf-debug --varset local_nvf    # NVIDIA, OpenACC, runtime checks
fobis build --mode flume-fnl-nvf-fft --varset local_nvf      # NVIDIA, OpenACC + cuFFT (FFT Poisson solver of self-gravity), exe adam_flume_fnl_fft
fobis build --mode flume-fnl-nvf-fft-debug --varset local_nvf  # the same, runtime checks
fobis build --mode flume-fnl-omp-amd                          # AMD, OpenMP offload
fobis build --mode flume-fnl-omp-amd-debug                    # AMD, OpenMP offload, runtime checks
```

The AMD build runs on an MI300X (amdflang): the `profile-couette3d` case, size `large`, WENO-5, 20 steps, on one and on
eight GPUs; commit `e4caa8ea` took its loop from 632.2 s to 50.7 s on one GPU (12.5x) ( [profile-couette3d README](https://github.com/szaghi/adam/tree/develop/src/tests/flume/speedup/profile-couette3d)).
Its results have not been compared with the FNL regression goldens.
