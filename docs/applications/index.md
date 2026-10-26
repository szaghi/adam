# Applications

Two applications are maintained against the current library, **FLUME** and **PRISM**. Both are built on the same
library objects (`adam_grid_object`, `adam_field_object`, `adam_weno_object`, …) through `realm_object`, share its
infrastructure (AMR, ghost-cell exchange, I/O, immersed boundaries, FLAIL) and specialise only in the equations and the
methods they need.

| Application | Equations solved | Physics domain | Backends |
|-------------|------------------|----------------|----------|
| [FLUME](/applications/flume/) | Compressible Euler and ideal MHD | Compressible flow, magnetohydrodynamics, astrophysics | CPU · FNL |
| [PRISM](/applications/prism/) | Maxwell + PIC | Electromagnetics and plasmas | CPU · FNL |

NASTO, CHASE, PATCH and ASCOT are in the source tree but unmaintained, and not documented on this site: NASTO, CHASE and
PATCH predate the realm refactor (CHASE and PATCH do not build), and FLUME supersedes CHASE.

## Common design pattern

Both applications have the same directory layout:

```
src/app/<name>/
├── common/        # Shared physics, BC, IC, IO: backend-independent
├── cpu/           # CPU: MPI (+ OpenMP)
└── fnl/           # GPU: OpenACC (nvfortran) or OpenMP offloading (amdflang), via FUNDAL
```

Each `<name>_<backend>_object` extends `<name>_common_object`, specialising the compute kernels while inheriting
configuration, I/O and physics logic.

## Configuration

Every application is configured through a single INI file, parsed by [FiNeR](https://github.com/szaghi/FiNeR). The
sections the library reads are shared; each application adds its own. Common sections:

| Section | Purpose |
|---------|---------|
| `[IO]` | Output basename, save frequency, restart options |
| `[time]` | Maximum simulated time or iterations, CFL number |
| `[grid]` | Domain bounds, cell counts `(ni, nj, nk)`, ghost-cell width `ngc` |
| `[amr]`, `[amr_marker_N]` | Refinement levels and the refinement markers |
| `[weno]`, `[runge_kutta]` | The spatial reconstruction and the time integration |
| `[numerics]`, `[physics]` | The application's schemes and physical model |
| `[bc_*]` | Boundary conditions for each domain face |
| `[solids]` | Immersed boundary geometry (analytical sphere, circle, rectangle) |

The full key lists are in the [FLUME input reference](/applications/flume/input) and the PRISM
[common page](/applications/prism/common).
