# Unit tests

The unit tests pin the building blocks below the solver: the eigensystems, the fluxes and Riemann solvers, the
positivity limiter, the WENO tables and weights, and the library pieces FLUME relies on (the tree lookup, the quadtree
seam exchange, the Poisson solvers). Each device twin (`_fnl`) also asserts the device result equal to the host one.

## FLUME

Sources in `src/tests/flume/unit/`.

| Test | What it pins |
|---|---|
| `test_flume_euler_library` (+ `_fnl`) | Euler eigensystem, flux, Roe average, split consistency; RS(q, q) = f(q) for LLF/HLL/HLLC, HLLC exact on a contact, positive first-order updates; the `si` descaler of planar and static states within $10^8$ of the largest (host); device = host |
| `test_flume_mhd_library` (+ `_fnl`) | MHD, GLM and EGLM eigensystems (including degenerate states), fluxes, auxiliary variables, cyclic invariance; the `si` descaler of planar, field-free and static states within $10^8$ of the largest; $\alpha_f$, $\alpha_s$ to full relative accuracy at high and low β with $\mathbf{B}$ near the normal, against quadruple precision ([#79](https://github.com/szaghi/adam/issues/79)) (host); device = host |
| `test_flume_mhd_riemann` (+ `_fnl`) | MHD LLF/HLL/HLLD without cleaning, with GLM and with EGLM: consistency, HLLD exact on contact, tangential and rotational discontinuities, cyclic invariance bitwise, positive updates, EGLM = GLM bitwise at $\psi = 0$; device = host |
| `test_flume_positivity` | PV-0, the limiter on random admissible states with perturbed fluxes (Euler, MHD, EGLM): every limited update positive and above the relative floor, the limiter needed and acting; a NaN or infinite high-order flux replaced by the backbone flux; inadmissible face ghosts blended above the floors, every other value untouched; a cell with a 2:1 seam face (mean donor-state backbone) admissible for any seam factor up to its own, and not without the limiter |
| `test_flume_weno_interpolation` (+ `_fnl`) | WENO interpolation tables: exactness, convergence, device = host for both weights (`js`, `si`); the reconstruction tables unchanged; the `si` weights scale-covariant bitwise on the device |
| `test_flume_weno_weights` | NV-1 of issue [#49](https://github.com/szaghi/adam/issues/49): with `[weno] weights = si`, WENO($2^n v$) = $2^n$ WENO($v$) bitwise for every scheme, both tables, random stencils from $10^{-12}$ to $10^{12}$; the `js` weights fail it (negative control); quadratic data reproduced by both; zero data stay zero |

```bash
fobis build --mode test-flume-<name>-gnu                         # host
fobis build --mode test-flume-<name>-fnl-nvf --varset local_nvf  # device twin
```

## Library tests used by the FLUME cases

| Test | Sources | What it pins | Case |
|---|---|---|---|
| `test_tree_closest_block` | `src/tests/amr/` | the tree lookup of the seam peers on binary trees, quadtrees and octrees, including a single-block tree ([#54](https://github.com/szaghi/adam/issues/54)) | [multi-realm runs](./euler/multi-realm) |
| `test_quadtree_seam_ghost` (+ `-fnl-nvf`) | `src/tests/amr/` | the 2:1 seam exchange on a linear field: octree, quadtree `nk = 4` and `nk = 1`, 1 to 3 ranks; the device exchange equal to the CPU one, bit for bit ([#46](https://github.com/szaghi/adam/issues/46)) | [quadtree AMR](./mhd/quadtree-amr) |
| `test_poisson_fft` (+ `-fnl-nvf`) | `src/tests/poisson/` | the FFT Poisson solver, periodic and isolated (SV-2, SV-6) | [FFT Poisson solver](./self-gravity/fft-poisson), [isolated solver](./self-gravity/isolated-solver) |
| `test_poisson_mg` | `src/tests/poisson/` | the multigrid Poisson solver on AMR grids (SV-9) | [multigrid on AMR](./self-gravity/multigrid-amr) |

Build modes: `test-tree-closest-block-gnu`, `test-quadtree-seam-ghost-gnu` (`-fnl-nvf`), `test-poisson-fft-gnu`
(`-fnl-nvf`), `test-poisson-mg-gnu`.
