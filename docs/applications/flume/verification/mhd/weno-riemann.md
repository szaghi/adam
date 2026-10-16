# MHD on `weno-riemann`

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| RV-5, RV-6, RV-7 | Ideal MHD | M3 ([#47](https://github.com/szaghi/adam/issues/47), M3-P3c) | CPU, FNL | yes | `rj2a-hlld` |

## Case and scope

The MHD verifications of M2 were built on the WENO flux splitting. The `weno-riemann` scheme (WENO interpolation of
the face states, an LLF, HLL or HLLD Riemann solver, a high-order flux correction) reuses their inputs and oracles
through the shared `--numerics` option, which rewrites the `[numerics]` block of the generated inputs:

```
--numerics SOLVER[:RECON[:CORRECTION[:SENSOR]]]
```

with `SOLVER` = `llf | hll | hlld`, `RECON` = `primitive | characteristic` (default `characteristic`, the MHD default),
`CORRECTION` = `none | 4th | 6th` (default `6th`), `SENSOR` = `weno | none` (default `weno`). Seven scripts take it:
`linear-wave`, `cpaw`, `riemann`, `rj2a`, `orszag-tang`, `rotor` and `field-loop`.

What the case catches:
- **a Riemann-solver path that loses accuracy:** the specs `hlld:primitive` and `hlld:characteristic` (6th, weno)
  assert bounds measured on them in `rj2a`, `linear-wave`, `cpaw` and `field-loop`; HLLD must be at least as accurate as
  HLL on Ryu–Jones 2a (RV-6);
- **a scheme-dependent defect:** every spec runs the scheme-independent checks (order, bitwise x/y/z and left/right,
  symmetry, conservation, positivity); the other specs are measurements.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | ideal MHD, no cleaning, GLM or EGLM, as in each verification | [models](../../models#ideal-mhd) |
| Space | `scheme_space = weno-riemann`: WENO interpolation of the face states (primitive or characteristic variables) | `compute_residuals_riemann`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `mhd_face_interpolation_fields`, `mhd_face_states`, [`adam_flume_mhd_riemann_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_riemann_library.F90); [numerics](../../numerics#scheme-space-weno-riemann-interpolation-riemann-flux-and-correction) |
| Riemann solvers | LLF, HLL, HLLD (with a fallback to HLL, reported in the log); GLM and EGLM variants | `mhd_riemann_llf`, `mhd_riemann_hll`, `mhd_riemann_hlld`, `mhd_glm_riemann_hlld`, `mhd_eglm_riemann_hlld`, [`adam_flume_mhd_riemann_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_riemann_library.F90); [`adam_flume_cpu_mhd_hlld_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_mhd_hlld_kernels.F90) (and the `_llf`, `_hll`, `_glm_*`, `_eglm_*` twins) |
| Flux correction | `6th` with the WENO sensor (default) | [numerics](../../numerics#scheme-space-weno-riemann-interpolation-riemann-flux-and-correction) |
| Time, boundaries, grids | as in each verification | — |
| FNL backend | the same kernels on the device | `compute_residuals_riemann_dev`, [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

## Initial and boundary conditions

The inputs are those of the verifications, unchanged except for the `[numerics]` block:

| Problem | Setup |
|---|---|
| [Ryu–Jones 2a](./ryu-jones-2a) | 256 and 512 cells, $t = 0.2$; RV-6 adds an HLL run at 256 |
| [linear waves, CPAW](./order-of-accuracy) | the MV-5 and MV-6 ladders |
| [Brio–Wu, Ryu–Jones 4d](./brio-wu-rj4d) | 512 cells, floors armed at $10^{-10}$ |
| [field loop](./field-loop) | the uniform 64/128 ladder |
| [Orszag–Tang](./orszag-tang), [rotor](./rotor) | as in MV-12, MV-13 |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/mhd/numerics.sh` | `numerics_check`, `numerics_tag`, `numerics_apply`: validates the spec and rewrites `[numerics]` |
| `src/tests/flume/verification/make_doc_runs.sh <runs> rj2a mhd`, `make_doc_figures.py --only weno-riemann` | `rj2a` re-runs the splitting and HLLD-characteristic runs, `mhd` the HLLD-primitive run `rj2a-hlld-prim` |
| `src/tests/flume/verification/mhd/{rj2a,linear-wave,cpaw,riemann,orszag-tang,rotor,field-loop}/check.sh` | the seven scripts that take `--numerics` |

```bash
cd src/tests/flume/verification/mhd
rj2a/check.sh --numerics hlld                     # HLLD, characteristic interpolation
rotor/check.sh --numerics hlld:primitive          # HLLD, primitive interpolation
```

## Results

![Ryu-Jones 2a on the splitting scheme and on HLLD](/flume/weno-riemann.png)

*Ryu–Jones 2a at 256 cells, CPU, over the exact solution: the splitting scheme and `weno-riemann` with HLLD on characteristic and on primitive interpolation. The three agree at this scale; the primitive interpolation has the largest $L_1$ sum of the three (table below).*

Measured in issue #47, M3-P3c; CPU and FNL agree to the printed digits:

| Problem | Splitting | HLLD primitive | HLLD characteristic |
|---|---|---|---|
| Ryu–Jones 2a, $L_1$ sum at 256 / 512 cells | $3.84 / 2.12\cdot10^{-2}$ | $5.13 / 2.93\cdot10^{-2}$ | $3.76 / 2.12\cdot10^{-2}$ |
| linear waves, CPAW: observed order | 4.96–5.02 | 4.96–5.00 | 4.96–4.99 |
| Brio–Wu, Ryu–Jones 4d: floored cells, HLLD fallbacks | 0 | 0, 0 | 0, 0 |
| field loop, magnetic energy kept at 64 / 128 cells | 0.872 / 0.937 | 0.869 / 0.936 | 0.887 / 0.944 |

Bounds asserted under the two HLLD specs (CPU M3-P3c + 2%; energies − 0.2%):

| Check | HLLD primitive | HLLD characteristic |
|---|---|---|
| Ryu–Jones 2a, $L_1$ sum at 256 / 512 | $5.235 / 2.985\cdot10^{-2}$ | $3.840 / 2.166\cdot10^{-2}$ |
| CPAW $\epsilon$ at 32 / 64 / 128 | $1.124\cdot10^{-3}$ / $3.774\cdot10^{-5}$ / $1.176\cdot10^{-6}$ | $1.098\cdot10^{-3}$ / $3.726\cdot10^{-5}$ / $1.174\cdot10^{-6}$ |
| linear waves, CPAW: observed order | $\ge 4.8$ | $\ge 4.8$ |
| field loop $\langle\lvert B_z\rvert\rangle/A_0$ at 64 / 128 | $9.950\cdot10^{-4}$ / $5.546\cdot10^{-4}$ | $1.083\cdot10^{-3}$ / $5.923\cdot10^{-4}$ |
| field loop $E_B(T)/E_B(0)$ at 64 / 128 | $\ge 0.86757$ / $0.93417$ | $\ge 0.88559$ / $0.94189$ |
| Ryu–Jones 2a, GLM against no cleaning | $10^{-11}$ | $10^{-11}$ |
| rotor symmetry at step 100 / at the end | $10^{-8}$ / $10^{-4}$ | $10^{-8}$ / $10^{-4}$ |

The linear-wave bounds are per run, the larger of the two interpolations (`linear-wave/check.sh`).

Related cases: [Sod](../euler/sod), [Lax](../euler/lax) and [Shu–Osher](../euler/shu-osher) compare the Euler Riemann
solvers; [EGLM cleaning](./eglm) runs on HLLD characteristic; [positivity limiter](./positivity-limiter) runs HLLD and
HLLC with the limiter.

## Findings

- On Ryu–Jones 2a HLLD is more accurate than HLL ($5.98$ and $4.32\cdot10^{-2}$): RV-6 asserts the HLLD $L_1$ sum at
  most the HLL one.
- HLLD with characteristic interpolation matches the splitting scheme on Ryu–Jones 2a; primitive interpolation is
  about 35% less accurate there. On the smooth waves the two interpolations agree to 3 digits (linear WENO weights).
- The GLM run differs from the run without cleaning by round-off ($10^{-12}$; bitwise with the splitting), because the
  interpolated uniform $B_n$ seeds $\psi$.
- The rotor's symmetry defect is round-off amplified by the flow and grows as the dissipation decreases (final
  $5\cdot10^{-9}$ to $6\cdot10^{-6}$ depending on the interpolation and the backend), so under `--numerics` it is
  checked at step 100 within $10^{-8}$ and at the end within $10^{-4}$ ([rotor](./rotor)).
- The scripts `glm-pulse`, `divb-peak`, `vortex` and `rotated-shock-tube` do not take `--numerics`: MV-3, MV-7, MV-8 and
  MV-10 run on the splitting scheme only.
