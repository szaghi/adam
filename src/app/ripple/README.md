# RIPPLE

RIPPLE solves the scalar linear advection `u_t + a · ∇u = 0` on a periodic box. It is the worked example of the
"build your own solver on the ADAM library" tutorial ([issue #59](https://github.com/szaghi/adam/issues/59)): an app
small enough that each tutorial chapter can quote whole files, and built and checked in CI so the tutorial cannot drift
from the code.

## What it shows

- **A realm on the library.** `ripple_common_object` extends `realm_object`, which owns the grid, tree, field, maps,
  WENO and Runge-Kutta objects. The app adds only the state `q`, the velocity, the time bookkeeping, the initial
  condition and the initial AMR.
- **The scheme** is Shu-Osher conservative finite difference on point values. The split fluxes `a+ u` and `a- u` are
  reconstructed at every face by the library's upwind WENO and advanced by the library's SSP Runge-Kutta.
- **The realm contract, two of its three layers:**
  - **core** (the forest's fast path, a uniform grid): `initialize_forest`, `compute_local_dt_forest`,
    `advance_one_step_forest`, `post_step_forest`, `is_done_forest`, `finalize_forest`;
  - **staged** (any 2:1 AMR face): `stages_per_step_forest`, `open_step_forest`, `begin_stage_forest` (sets
    `stage_active`), `end_stage_forest` (the seam skins into the flux register), `close_step_forest` (resets
    `stage_active`) and the end-of-step reflux.

  RIPPLE is single-realm, so it has no seams layer.
- **The FNL backend** as a sibling of the CPU one: the state lives on the device in transposed `(b, i, j, k, v)` order,
  every `dev_alloc` is checked, and the library's `field_fnl`, `rk_fnl` and `weno_fnl` do the device side.
- **A refinement criterion of your own** ([issue #110](https://github.com/szaghi/adam/issues/110)): a type extending the
  library's `amr_criterion_object`, registered by the program (`register_amr_criterion`) and selected by
  `[amr_marker_N] mode = 6, name = ...`; the example `crest` criterion and its driver are
  `src/tests/ripple/test_ripple_crest_criterion.F90`.
- **Refusals at initialization** for every input it cannot honour: a non-SSP scheme, a centred WENO, `weights = si`,
  `ngc` below the WENO stencil half-width, runtime AMR (`[amr] frequency > 0`), immersed solids, a marker other than
  the box and a registered criterion, an unregistered criterion name, a criterion that acts at runtime regrids only.

## Layout

| Path | Content |
|---|---|
| `common/adam_ripple_common_object.F90` | The common realm: input, initial condition, exact solution, initial AMR, history and output, the backend-independent contract TBPs |
| `cpu/adam_ripple_cpu_object.F90` | The CPU realm: residuals, Runge-Kutta stages, flux register and reflux |
| `cpu/adam_ripple_cpu.F90` | The driver: one realm handed to `forest%simulate` |
| `fnl/adam_ripple_fnl_object.F90` | The FNL realm: device state, device stages, device reflux |
| `fnl/adam_ripple_fnl_kernels.F90` | The OpenACC kernels |
| `fnl/adam_ripple_fnl.F90` | The FNL driver |

## Build and check

```bash
fobis build --mode ripple-cpu-gnu                       # CPU
fobis build --mode ripple-fnl-nvf --varset local_nvf    # FNL (OpenACC, nvfortran)
fobis build --mode test-ripple-crest-criterion-cpu-gnu  # the driver with a registered criterion (oracle r6)
./src/tests/ripple/check.sh                             # oracles r1-r3, r5, r6 (CPU); --legs r4 compares the FNL run
```

The oracles (`src/tests/ripple/check.sh`, `oracle.py`) compare with the exact translated profile:
- r1, uniform grid: L1 order 5;
- r2, a static 2:1 region: L1 order 2, the cap of the library's 2:1 transfers, and the total of `u dV` conserved to
  round-off;
- r3, the AMR run identical on 1 and 2 ranks;
- r4, CPU and FNL agree to 1e-12;
- r5, CRLF input bitwise the LF one;
- r6, a registered criterion (`crest`) refines exactly the blocks where the initial profile rises above 0.9, the
  total of `u dV` conserved, and an unregistered name refused.

CI runs r1–r3, r5 and r6.

**Backends:** CPU and FNL on NVIDIA (OpenACC). The kernels carry no `!$omp` twins and there is no AMD
(`*-fnl-omp-amd`) mode; FLUME shows that pattern.
