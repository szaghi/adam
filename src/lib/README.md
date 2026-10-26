# ADAM Library — `src/lib`

The ADAM SDK: the physics-agnostic, hardware-portable core every application is built on. PRISM and FLUME build the
common layer and the FNL backend; no maintained application builds NVF or GMP.

**Documentation:** <https://szaghi.github.io/adam/concepts/library>

| Directory | Role | Parallelism | Entry point |
|-----------|------|-------------|-------------|
| `common/` | Portable core objects: grid, tree, field, numerics, FLAIL, I/O, realm and forest | CPU · MPI · OpenMP | `adam_common_library` |
| `fnl/` | GPU backend (FUNDAL-based), built by PRISM and FLUME | OpenACC or OpenMP target · MPI | `adam_fnl_library` |
| `nvf/` | CUDA Fortran backend; no build mode builds it | CUDA Fortran · MPI | `adam_nvf_library` |
| `gmp/` | OpenMP target offloading backend; no build mode builds it | OpenMP target · MPI | `adam_gmp_library` |

Each subdirectory has its own `README.md` with its module table.

ADAM is dual-licensed under MIT or LGPLv3 (`LICENSE.mit.md`, `LICENSE.lgpl3.md`).
