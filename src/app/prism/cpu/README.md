# PRISM — CPU Backend

The CPU backend (`prism-cpu`) provides a pure CPU/MPI implementation of the PRISM Maxwell equations
solver. It inherits the full PRISM common infrastructure (grid, field, physics, numerics, I/O, PIC,
coils, boundary conditions) and adds flux storage arrays, spatial-operator procedure pointers, and the
Runge-Kutta time integration.

**Documentation:** <https://szaghi.github.io/adam/applications/prism/backends>

## Source files

| File | Module |
|------|--------|
| `adam_prism_cpu.F90` | entry point — parses `--filename` CLI argument, delegates to `prism%simulate` |
| `adam_prism_cpu_object.F90` | `adam_prism_cpu_object` — type declaration, all procedures |

ADAM is dual-licensed under MIT or LGPLv3 (`LICENSE.mit.md`, `LICENSE.lgpl3.md`).
