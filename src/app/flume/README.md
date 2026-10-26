# FLUME

**FLUME** — **F**luid **L**aboratory for **U**nsteady **M**ultiphysics at **E**xascale — is the compressible-flow
application of the ADAM framework: compressible Euler, Navier–Stokes and ideal MHD (GLM or EGLM cleaning, Ohmic
resistivity); external gravity, self-gravity and sink particles; ideal, isothermal and barotropic gases; optically thin
radiative cooling; dual energy and passive scalars; immersed boundaries; high-order WENO schemes with a positivity
limiter; static and runtime AMR and multi-realm forests; on CPUs (MPI, OpenMP) and GPUs (OpenACC, OpenMP offload).

**Documentation:** <https://szaghi.github.io/adam/applications/flume/>

## Source Layout

```
src/app/flume/
├── common/          # Shared across all backends (physics, numerics, BC, IC, I/O, ...)
├── cpu/             # CPU backend (MPI + OpenMP)
└── fnl/             # GPU backend (FNL library: OpenACC, or OpenMP offload on amdflang)
```

Module and file naming follows the ADAM convention: `adam_flume_<name>_object.F90` for types,
`adam_flume_<name>_library.F90` for procedure collections, `adam_flume_cpu.F90` / `adam_flume_fnl.F90` for the program
entry points. The per-directory READMEs list the modules.

| Backend | Entry point | Key type |
|---------|-------------|----------|
| CPU | `cpu/adam_flume_cpu.F90` | `flume_cpu_object` |
| FNL | `fnl/adam_flume_fnl.F90` | `flume_fnl_object` |

Build modes: `fobis build --lmodes`, or the [build-modes page](https://szaghi.github.io/adam/reference/build-modes).

## License

FLUME is part of ADAM, dual-licensed under the [MIT License](../../../LICENSE.mit.md) and the
[GNU LGPLv3](../../../LICENSE.lgpl3.md); you may choose either.
