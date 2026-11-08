# Absorbing layers on AMR

| Cases | Issues | Type | Backends | Figure | Regression |
|---|---|---|---|---|---|
| `rmf-fwl`, `rmf-fwl-amr`, `rmf-pml-amr` | [#31](https://github.com/szaghi/adam/issues/31), [#58](https://github.com/szaghi/adam/issues/58), [#95](https://github.com/szaghi/adam/issues/95) (P5.1) | C | CPU, FNL | no | `rmf-fwl` (digest, FNL only); `rmf-fwl-amr`, `rmf-pml-amr` (digest + residuals, CPU and FNL) |

## Case and scope

Three digest-goldened cases anchor the absorbing layers: the far-wave layer (`[fWLayer]`) and the perfectly matched
layer (`[PML]`) on the six faces of the four-coil `rmf` problem, on a uniform grid and after an initialization-time
refinement. They carry no oracle of their own: the digest pins the solution, and the defects they guard against moved
it far outside the tolerance.

What the cases catch:
- **the FNL host→device copy of the layer** (`rmf-fwl`, [#31](https://github.com/szaghi/adam/issues/31)): a
  buffered transposed copy through a staging buffer whose shape differed from the device destination made nvfortran
  pass a host temporary to the asynchronous copy (`CUDA_ERROR_INVALID_VALUE` at more than one rank);
- **a layer width that does not fit a block** ([#58](https://github.com/szaghi/adam/issues/58)): the layer is set by a
  physical `width`, and the cell count per face is derived, $C_\text{face} = \min(n_i, \lceil \text{width}/\Delta s\rceil)$;
- **per-block face data left on the pre-AMR blocks** (`rmf-fwl-amr`, `rmf-pml-amr`,
  [#95](https://github.com/szaghi/adam/issues/95) P5.1): the layer ranges, the PML support and its RK storage must be
  rebuilt on every change of the block set. Before the fix they stayed on the old blocks, 48 of the 96 face blocks per
  rank, and the initial condition was built with the stale layer.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | Maxwell, SI, 4 AC rectangular coils | [models](../models#models); [coils](../initial-conditions#coils) |
| Space, time | centred finite differences of order 6; SSP Runge–Kutta (5,4), CFL 0.4997 | [numerics](../numerics#time-integration) |
| Far-wave layer | six faces, `width = 0.045` m | [`adam_prism_fWLayer_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/prism/common/adam_prism_fWLayer_object.F90); [absorbing layers](../pml#far-wave-layer) |
| PML (`rmf-pml-amr`) | `CFS`, six faces, `width = 0.045` m, `gamma_max = 7.8e10`, `alfa_max = 2.5e8`, `k_max = 5` | [`adam_prism_pml_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/prism/common/adam_prism_pml_object.F90), [`adam_prism_rk_pml_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/prism/common/adam_prism_rk_pml_object.F90); [absorbing layers](../pml#perfectly-matched-layer) |
| Face data rebuild | after every change of the block set | `initialize_block_face_data`, [`adam_prism_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/prism/common/adam_prism_common_object.F90) |

## Initial and boundary conditions

| Item | `rmf-fwl` | `rmf-fwl-amr` | `rmf-pml-amr` |
|---|---|---|---|
| Domain | $[-0.24, 0.24]^3$ m | the same | the same |
| Grid | 64 blocks of $16^3$ (`iu_ref_levels = 2`), $\Delta x = 7.5\cdot10^{-3}$: the layer is 6 cells | plus a box marker fired at initialization (`amr_iterations = 1`) refining the slab $x < -0.13$ to level 3 | the same as `rmf-fwl-amr` |
| Layer | fWLayer, six faces | fWLayer, six faces | CFS PML, six faces |
| Initial condition | `vacuum`, 4 AC coils (`time_delay = 1e-8`) | the same | the same |
| Steps, ranks | 10, 2 | 5, 2 | 5, 2 |

## Input files

| File | Role |
|---|---|
| `src/tests/prism/regression/rmf-fwl/input.ini` | the uniform fWLayer case |
| `src/tests/prism/regression/rmf-fwl-amr/input.ini`, `rmf-pml-amr/input.ini` | the refined cases |
| `src/tests/prism/regression/<case>/golden/{cpu,fnl}/` | the digests `run.sh` compares |
| `src/tests/lint/check-removed-input-keys.sh` | rejects `[fWLayer] C` (the former cell count) in every committed input; `run.sh` runs it first |

```bash
REGRESSION_CASES="rmf-fwl rmf-fwl-amr rmf-pml-amr" src/tests/prism/regression/run.sh cpu --no-build
src/tests/prism/regression/run-fnl-local.sh --no-build                  # FNL, the whole sweep
```

## Results

| Case | Asserts | Bound | CPU | FNL |
|---|---|---|---|---|
| `rmf-fwl` | digest against the golden | `rtol = 1e-6`, `atol = 1e-3` | no golden (SKIP) | pass |
| `rmf-fwl-amr` | digest and residuals against the goldens | the same | pass | pass |
| `rmf-pml-amr` | digest and residuals against the goldens | the same | pass | pass |
| removed-keys lint | no `[fWLayer] C` in a committed input | none found | pass | — |

Related pages: [input validation and refusals](./refusals) for the inputs PRISM refuses; [hazards](../hazards#pml-constraints)
for the PML parameter sets. The PML reflection against the angle of incidence is an open analytic oracle
([#106](https://github.com/szaghi/adam/issues/106)).

## Findings

- The fWLayer field now holds no device state and `copy_cpu_gpu` copies only `q`, whose staging buffer has the shape
  of its device destination, so the #31 crash cannot recur. The rule stands for every buffered transposed copy: pass a
  staging buffer only when its shape equals the destination's.
- On a refined grid the layer and the PML follow the final block set: both refined cases reproduce their goldens on
  both backends.
- These cases verify that the layers are applied where they should be, not how well they absorb.
