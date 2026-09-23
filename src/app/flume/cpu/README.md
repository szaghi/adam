# FLUME — CPU Backend

> FLUME is part of the ADAM framework, released under the
> [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.en.html) (LGPLv3).
> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.

The CPU backend (`flume-cpu`) will provide a CPU/MPI implementation of the FLUME compressible MHD solver, inheriting the backend-independent infrastructure from `common/`.

> **Status: under development.** No sources exist yet.

## Planned Source Files

| File | Description |
|------|-------------|
| `adam_flume_cpu.F90` | Entry point (program) |
| `adam_flume_cpu_object.F90` | `flume_cpu_object` implementation |

## Build mode

Not defined yet; expected form:

```bash
fobis build --mode flume-cpu-gnu
```
