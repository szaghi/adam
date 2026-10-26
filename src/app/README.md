# Applications — `src/app`

Solver applications built on the ADAM library. **FLUME** and **PRISM** are maintained against the current library;
NASTO, CHASE and PATCH predate the realm refactor and are not maintained (CHASE and PATCH do not build), and ASCOT is a
standalone slice converter.

**Documentation:** <https://szaghi.github.io/adam/applications/>

| Directory | Application | Status |
|-----------|-------------|--------|
| `flume/` | Compressible Euler and ideal MHD (CPU · FNL) | Maintained |
| `prism/` | Maxwell + PIC (CPU · FNL) | Maintained |
| `nasto/` | Compressible Navier-Stokes | Not maintained |
| `chase/` | Euler; superseded by FLUME | Not maintained; does not build |
| `patch/` | Poisson | Not maintained; does not build |
| `ascot/` | Binary slices to ASCII (Tecplot) converter, no ADAM objects | Standalone utility |

Each maintained application has `common/`, `cpu/` and `fnl/` subdirectories; `<name>_<backend>_object` extends
`<name>_common_object`.

ADAM is dual-licensed under MIT or LGPLv3 (`LICENSE.mit.md`, `LICENSE.lgpl3.md`).
