# MHD plumbing

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| — (`mhd/plumbing`, M2-P1) | Ideal MHD | M2 ([#41](https://github.com/szaghi/adam/issues/41), P1) | CPU, FNL | no | — |

## Case and scope

M2-P1 wired the ideal-MHD model through FLUME end to end: the input predicate (state width and names), the MHD
conversions and fast magnetosonic speed, the per-model kernel instances, the initial- and boundary-condition keys and
the wall rule, the outputs and the restart. The leg checks that plumbing, for both variants
`divergence_control = none` ($n_v = 8$) and `glm` ($n_v = 9$, with $\psi$), on a seeded uniform state:

1. **time step:** one step of the unperturbed uniform state ($s = 0$). The step must equal
   $\mathrm{CFL} / \sum_d (|u_d| + c_{f,d})/\Delta x_d$, with the fast speed computed by the oracle independently of the
   Fortran; under GLM it is bounded by $\mathrm{CFL} / (c_h \sum_d 1/\Delta x_d)$, the active bound here (M2-P4);
2. **auxiliary fields:** the saved MHD auxiliaries equal the values recomputed from the saved conservative fields;
3. **restart:** 5 steps, a restart and 5 more steps of the seeded ($s = 0.05$) state equal the continuous 10-step run
   bit for bit, the residual and conservation histories byte for byte, and every model variable
   (`bx`, `by`, `bz` [, `psi`]) is saved;
4. **refusals** (variant `none`): an unknown `divergence_control` and `isentropic-vortex` with MHD stop with their
   error messages.

What the case catches: a wrong fast-speed formula, a wrong $\mathbf{B}$ index, a wrong auxiliary layout on the device
(the dt leg exercises the device auxiliary and fast-speed kernels on FNL), a kernel instance of the wrong model,
$c_h$ missing from the time step, a missing $\mathbf{B}$ or $\psi$ field in the output, a restart that loses $\psi$,
and a silently accepted invalid input.

::: details The retired zero-residual legs
At M2-P1 the MHD face fluxes did not exist yet, so the residual was exactly zero and the leg also required constant
integrals and 10 steps equal to 1 step. Since M2-P3 the face fluxes are live and those legs are retired; the exact
steady state is checked by MV-2 ([zero field](./zero-field), uniform state, bitwise). The oracle's module docstring
still describes the P1 state ("no face fluxes yet"); its `--static` mode is no longer called by `check.sh`.
:::

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | ideal MHD, $\gamma = c_p/c_v = 2.5/1.5 = 5/3$; `none` ($n_v = 8$) or mixed GLM ($n_v = 9$) | [models](../../models#ideal-mhd); `initialize`, [`adam_flume_physics_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_physics_object.F90) |
| GLM | $c_h = 3$, `glm_alpha = 0.18`, `glm_damping_length = 1`, `glm_ch_check = warning` | [models](../../models#divergence-control-mixed-glm); [`adam_flume_mhd_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_object.F90) |
| Space | WENO5 (`weno-u-5`) flux splitting on characteristic variables | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `mhd_face_split_fluxes`, `mhd_glm_face_split_fluxes`, [`adam_flume_mhd_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_library.F90) |
| Time step | $\mathrm{CFL} / \max \sum_d (|u_d| + c_{f,d})/\Delta x_d$, GLM bound $c_h \sum_d 1/\Delta x_d$ | `compute_local_dt_forest`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `compute_lambda_max`, [`adam_flume_cpu_mhd_kernels_agnostic.INC`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_mhd_kernels_agnostic.INC); `mhd_fast_speed`, [`adam_flume_mhd_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_library.F90); [models](../../models#time-step-bound) |
| Time | SSP Runge–Kutta, 3 stages, CFL 0.8 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Auxiliary fields | $\rho, u, v, w, p, T, H, a, B_x, B_y, B_z$ from the conservative state | `compute_q_aux`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `mhd_conservative_to_auxiliary`, [`adam_flume_mhd_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_library.F90); [models](../../models#derived-output-fields) |
| Output and restart | HDF5 checkpoints; restart files | `save_xh5f`, `save_restart_files`, `load_restart_files`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90) |
| Refusals | unknown `divergence_control`; `isentropic-vortex` requires `euler` | `load_from_file`, [`adam_flume_mhd_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_object.F90), [`adam_flume_ic_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_ic_object.F90) |
| Boundaries | `extrapolation`, `wall-inviscid`, `inflow` (with $\mathbf{B}$) | [`adam_flume_bc_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_bc_object.F90); [boundary conditions](../../boundary-conditions#wall-inviscid-slip-wall-by-mirroring-with-a-sign-vector) |
| FNL backend | the same kernels on the device | [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | $[0, 1]^3$ |
| Grid | octree, two uniform levels (`iu_ref_levels = 2`): 64 blocks of $8^3$ cells, $32^3$ cells, $\Delta x = 1/32$, `ngc = 3` (both ranks own blocks, so the restart leg is meaningful) |
| Ranks | 2 |
| Initial condition | [`uniform`](../../initial-conditions#uniform-seeded-uniform-state): $(\rho, u, v, w, p) = (1, 0.3, -0.2, 0.1, 1)$, $\mathbf{B} = (0.8, 0.5, -0.3)$, seeded density/pressure perturbation $s = 0.05$ ($s = 0$ in the dt leg) |
| Boundaries | $x$: `extrapolation`; $y$: `wall-inviscid`; $z_{\min}$: `inflow` with the region state; $z_{\max}$: `extrapolation` |
| Steps | dt leg 1; auxiliary leg 10; restart leg 5 + 5 against 10 |

| Leg | Variant | Change | Expected |
|---|---|---|---|
| dt | `none`, `glm` | $s = 0$, `it_max = 1` | step 1 equals the oracle's dt to $10^{-13}$ relative |
| aux | `none`, `glm` | — | auxiliaries equal the recomputed ones to $10^{-12}$ relative |
| restart | `none`, `glm` | `restart_save = 5`, then `restart = .true.` | bitwise fields, identical histories, model variables saved |
| refused | `none` | `divergence_control = bogus` | `unknown [mhd].(divergence_control) "bogus"` |
| refused | `none` | `type = isentropic-vortex` | `isentropic-vortex requires [physics].(physical_model) = euler` |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/mhd/plumbing/mhd-uniform-none.ini` | the case, `divergence_control = none` |
| `src/tests/flume/verification/mhd/plumbing/mhd-uniform-glm.ini` | the case, `divergence_control = glm` |
| `src/tests/flume/verification/mhd/plumbing/check.sh` | runs every leg (`--np N`, `--build`, `--legs LIST`: `none`, `glm`, `refused`; `FLUME_EXE` selects the backend) |
| `src/tests/flume/verification/mhd/plumbing/mhd_plumbing_oracle.py` | `--dt`, `--aux`, `--compare`; also `--zero` and `--steady` used by [MV-2](./zero-field) |

```bash
cd src/tests/flume/verification/mhd/plumbing && ./check.sh           # CPU
FLUME_EXE=$PWD/../../../../../../exe/adam_flume_fnl ./check.sh       # FNL
```

## Results

A contract leg: no figure. The work directories of both backends are kept
(`plumbing/work-adam_flume_{cpu,fnl}-np2-*`), with their logs and histories; the checkpoints are deleted, and the
oracle's printed lines are not kept. The table quotes what those files still show:

| Check | Bound | CPU | FNL |
|---|---|---|---|
| dt of step 1, `none` (history `*-dt`) | oracle $4.7620653358110355\cdot10^{-3}$, $10^{-13}$ relative | $4.7620653358110363\cdot10^{-3}$ ($1.7\cdot10^{-16}$) | the same digits |
| dt of step 1, `glm`: the $c_h$ bound $0.8/(3 \cdot 3 \cdot 32)$ is active ($c_h \cdot 3 = 9$ against $\sum_d (\lvert u_d\rvert + c_{f,d}) = 5.25$) | oracle $2.7777777777777778\cdot10^{-3}$ | $2.7777777777777779\cdot10^{-3}$ | the same digits |
| restart: residual and conservation histories against the continuous run (`*-A` vs `*-B`, `cmp`) | byte-identical | identical, both variants | identical, both variants |
| auxiliary fields against the recomputed ones | $10^{-12}$ relative | not kept | not kept |
| restart: interior fields, every variable, model variables present | bitwise | not kept | not kept |
| unknown `divergence_control` | refused | `unknown [mhd].(divergence_control) "bogus"; expected one of glm, none` | the same |
| `isentropic-vortex` with MHD | refused | `[initial_conditions].(type) = isentropic-vortex requires [physics].(physical_model) = euler` | the same |

The dt oracle values are recomputed from the oracle's formula and the input; the measured values are the step-1 rows
of the kept conservation histories.

Related cases: [zero field](./zero-field) (MV-2) took over the exact-steady-state legs; [GLM divergence
cleaning](./glm-cleaning) verifies the $c_h$ time-step bound and the damping physically.

## Findings

- The fast-speed time step agrees with the independent textbook formula to one ulp on both backends, and the GLM
  step is set by $c_h$, so a $c_h$ missing from the time step would fail the leg.
- The restart reproduces the continuous run's histories byte for byte on both backends and both variants.
- The work directories `*-none-C` and `*-glm-C` (dated before the others) belong to the retired "10 steps = 1 step"
  leg; `check.sh` no longer writes them.
- The auxiliary and bitwise-field verdicts rest on the `check.sh` run; the checkpoints they read are deleted after the
  oracle, so this page cannot re-derive them.
