# Static AMR and the FV reflux register

| Cases | Issues | Type | Backends | Figure | Regression |
|---|---|---|---|---|---|
| `rmf-amr` | [#13](https://github.com/szaghi/adam/issues/13) (M1–M3), [#23](https://github.com/szaghi/adam/issues/23), [#24](https://github.com/szaghi/adam/issues/24), [#28](https://github.com/szaghi/adam/issues/28) | C | CPU, FNL | no | `rmf-amr` (digest + residuals, CPU and FNL) |

## Case and scope

The first PRISM case with a 2:1 coarse–fine face. A geometric box marker refines the $x < 0$ half of the domain once
at initialization; the run then advances on the finite-volume path, where the Berger–Colella reflux replaces the
coarse flux on each 2:1 face by the sum of the fine fluxes.

The oracle is structural: it asserts that each stage of the machinery exists and is wired to the real tree, against a
control run in which the marker is off.

What the case catches:
- **a refine path that does nothing:** no 2:1 jump forms, or the uniform refinement alone forms one;
- **a registration pass that ignores the tree:** seam faces registered without a coarse–fine adjacency, or none with
  one;
- **a silent reflux:** no non-zero coarse/fine flux mismatch across the 2:1 face;
- **a decomposition-dependent register** ([#28](https://github.com/szaghi/adam/issues/28)): a face printed twice at
  2 ranks, a stuck zero entry, a shifted value.

The digest golden of the case was refreshed twice, each time adjudicated, not rubber-stamped: for the reflux fixes
(the face-side code of intra-realm seams, `e883ff93`, and the stage-weighted flux accumulation) and for the
scheme-aware coil stamping (`bad09954`).

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | Maxwell, SI, 4 AC rectangular coils | [models](../models#models); [coils](../initial-conditions#coils) |
| Space | centred finite volumes of order 6 (`fv_centered`) | [numerics](../numerics#centred-finite-volumes-fv-centered) |
| Time | SSP Runge–Kutta (5,4), CFL 0.5 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Refinement | AMR_GEO `primitive-box` marker, fired once at initialization (`amr_iterations = 1`) | [numerics](../numerics#adaptive-mesh-refinement) |
| Seam registration | intra-realm 2:1 faces registered from the tree adjacency | `register_intra_realm_amr_seams`, [`adam_forest_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_forest_object.F90) |
| Reflux | fine fluxes accumulated per RK stage with the stage weight, applied to the coarse cells | `apply_reflux_corrections`, [`adam_forest_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_forest_object.F90); `apply_reflux_to_stage_forest`, [`adam_prism_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/prism/cpu/adam_prism_cpu_object.F90); [`adam_flux_register_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_flux_register_object.F90) |
| FNL backend | the same reflux on the device | `apply_reflux_to_stage_forest`, [`adam_prism_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/prism/fnl/adam_prism_fnl_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | $[-0.06, 0.06]^3$ m |
| Grid | octree, 8 blocks of $16^3$ at level 1 ($32^3$ cells, $\Delta x = 3.75\cdot10^{-3}$); the blocks of $x < 0$ refined to level 2 ($\Delta x = 1.875\cdot10^{-3}$); `ngc = 3` |
| Initial condition | `vacuum`; the four coils carry an AC current |
| Boundaries | `Neumann` on all six faces |
| Steps | 5 |
| Control | the same input with `markers_number = 0` |
| Ranks | 1 (main, control); 1 and 2 (register parity) |

## Input files

| File | Role |
|---|---|
| `src/tests/prism/regression/rmf-amr/input.ini` | the case |
| `src/tests/prism/regression/rmf-amr/check.sh` | runs the main, control and parity runs (`--build`; `PRISM_EXE` selects the backend) |
| `src/tests/prism/regression/rmf-amr/golden/{cpu,fnl}/` | the digest `run.sh` compares |

```bash
cd src/tests/prism/regression/rmf-amr && ./check.sh                # CPU
PRISM_EXE=$PWD/../../../../../exe/adam_prism_fnl ./check.sh        # FNL
```

## Results

Each leg of `check.sh`, as coded:

| Leg | Asserts | CPU | FNL |
|---|---|---|---|
| main, health | no error/abort/NaN; the time loop reaches 100 % | pass | pass |
| main, blocks | `nodes number` > 1 | 36 | 36 |
| main, levels | both level 1 and level 2 present (a 2:1 jump) | 1, 2 | 1, 2 |
| main, registration | `registered intra-realm AMR seam faces` > 0 | 4 | 4 |
| main, reflux | last `reflux max\|F_coarse-F_fine_sum\|` > 0 | 1.510 | 1.510 |
| control, levels | no level 2 | pass | pass |
| control, registration | exactly 0 faces | 0 | 0 |
| control, reflux | no reflux line printed | pass | pass |
| parity, face count | np 1 = np 2, > 0 | 4 = 4 | 4 = 4 |
| parity, per face | every face printed once per step at np 2, the values string-equal to np 1 | pass, 4 faces × 5 steps | pass, 4 faces × 5 steps |

Related pages: [seam div(B)](./seam-divergence) runs the same grid on the FD path; [FV reflux on a capped last
step](./reflux-capped-step) measures the conservation the register buys.

## Findings

- The box marker, not the uniform refinement, creates the 2:1 face: the control run has one level and registers no
  face.
- Field digests cannot see the reflux correction at this scale, so the parity leg asserts the register content
  directly, face by face and step by step.
- At 2 ranks each face is printed once per step, by its coarse-owner rank, with the values of the 1-rank run. The
  block groups are rank-complete at this decomposition, so the cross-rank reduction adds exact zeros and string
  equality is the intended strength.
