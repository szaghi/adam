# FLUME Common

The `common/` directory will contain all backend-independent modules shared by the CPU and FNL backends of FLUME. Following the PRISM design, every type defined here is expected to be aggregated into a `flume_common_object`, the base type for both backends.

> **Status: under development.** No modules exist yet.

## Planned Content

The expected set of modules, by analogy with PRISM (names and scope subject to change):

| Concern | Description |
|---------|-------------|
| Parameters | Physical constants and eigensystem of the ideal MHD equations |
| Physics | Physical model, variable counts, equation of state |
| Numerics | Temporal and spatial scheme selection |
| Common object | `flume_common_object` base aggregate type |
| Boundary conditions | Boundary conditions handler |
| Initial conditions | Initial conditions handler |
| I/O | Input/output handler |
| Time | Time integration handler |
| Riemann solvers | MHD Riemann solvers and flux routines |
| Library | Barrel re-export of all common modules |

## License

FLUME is part of the ADAM framework, released under the [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.html) (LGPLv3).

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.
