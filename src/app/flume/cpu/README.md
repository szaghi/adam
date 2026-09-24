# FLUME — CPU Backend

> FLUME is part of the ADAM framework, released under the
> [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.en.html) (LGPLv3).
> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.

The CPU backend (`flume-cpu`) is the MPI (+ optional OpenMP) implementation of FLUME, extending the backend-independent
infrastructure of `common/`.

> **Status: development** (milestone M1, compressible Euler).

## Source Files

| File | Description |
|------|-------------|
| `adam_flume_cpu.F90` | Entry point (program) |
| `adam_flume_cpu_object.F90` | `flume_cpu_object`: ghost update and boundary conditions, WENO residual (with the immersed-boundary variant), Runge-Kutta integration, reflux accumulation and application, output |

## Build modes

```bash
fobis build --mode flume-cpu-gnu          # MPI
fobis build --mode flume-cpu-gnu-omp      # MPI + OpenMP
fobis build --mode flume-cpu-gnu-debug    # MPI, runtime checks
```
