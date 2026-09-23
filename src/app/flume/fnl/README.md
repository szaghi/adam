# FLUME — FNL (OpenACC) Backend

> FLUME is part of the ADAM framework, released under the
> [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.en.html) (LGPLv3).
> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.

The FNL backend (`flume-fnl`) will be the GPU-accelerated implementation of FLUME using the
**FNL library** (OpenACC-based GPU framework), inheriting the backend-independent infrastructure from `common/`.

> **Status: under development.** No sources exist yet.

## Planned Source Files

| File | Description |
|------|-------------|
| `adam_flume_fnl.F90` | Entry point (program) |
| `adam_flume_fnl_object.F90` | `flume_fnl_object` implementation |
| `adam_flume_fnl_kernels.F90` | OpenACC GPU kernels |
| `adam_flume_fnl_library.F90` | Barrel re-export of all FNL modules |

## Build modes

Not defined yet; expected form:

```bash
fobis build --mode flume-fnl-nvf --varset local_nvf
```
