# MHD with a zero field

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| MV-2 | Ideal MHD | M2 ([#41](https://github.com/szaghi/adam/issues/41), P3) | CPU, FNL | no | `uniform-amr-mhd` |

## Case and scope

The MHD face kernels (the $7 \times 7$ Roe–Balsara core at the face average, with the decoupled $B_n$ or the GLM pair)
replace the Euler ones whenever `physical_model = mhd-ideal`. With $\mathbf{B} = 0$ the MHD system is the Euler one
plus inert fields, which gives exact expectations before any MHD-specific test ([MV-4](./ryu-jones-2a)). Both variants,
`divergence_control = none` and `glm`, run two legs:

1. **Sod along $x$, $y$, $z$ through the MHD path:** $L_1(\rho)$ against the exact Riemann solution within the M1
   bound of [V1](../euler/sod); the $y$ and $z$ runs equal the $x$ run **bit for bit** after the permutation of the
   axes; $\mathbf{B}$ (and $\psi$) stay exactly zero;
2. **a uniform MHD state** ($\rho$, $\mathbf{u}$, $p$ and $\mathbf{B}$ all non-zero) on the periodic box of
   [V3](../euler/amr-conservation):
   - single level: the last checkpoint equals the step-0 state **bit for bit** (the MHD fluxes of equal states cancel);
   - with the 2:1 refined octant: equal to round-off, $10^{-14}$ relative.

What the case catches: an MHD kernel that does not reduce to Euler at $\mathbf{B} = 0$, a direction-dependent defect
of the MHD eigensystem, a spurious source of $\mathbf{B}$ or $\psi$, and MHD fluxes of equal states that do not
cancel, on a uniform grid and across 2:1 faces.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | ideal MHD, $\gamma = c_p/c_v = 1.4$ (the $c_p, c_v$ of the Euler base inputs, kept by `make_input.py`) | [models](../../models#ideal-mhd) |
| Space | WENO5 (`weno-u-5`) flux splitting on characteristic variables; the 7×7 Roe–Balsara core at the arithmetic face average | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `mhd_face_average`, `mhd_eigenvectors`, `mhd_face_split_fluxes`, `mhd_glm_face_split_fluxes`, [`adam_flume_mhd_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_library.F90); [numerics](../../numerics#scheme-space-weno-flux-splitting) |
| Cleaning | `none`, or GLM with $c_h = 3$, `glm_alpha = 0.18`, `glm_damping_length = 1`, `glm_ch_check = warning` | [models](../../models#divergence-control-mixed-glm); `add_glm_damping`, [`adam_flume_cpu_mhd_glm_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_mhd_glm_kernels.F90) |
| Time | SSP Runge–Kutta, 3 stages, CFL 0.5 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| AMR (uniform, AMR leg) | geometric box marker on the octant; 2:1 ghost restriction and interpolation; reflux | `mark_by_geometry`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90); `interp_seam_ghost`, [`adam_field_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_field_object.F90); [boundary conditions](../../boundary-conditions#_2-1-amr-coarse-fine-ghosts) |
| Boundaries | Sod: `extrapolation`; uniform: `periodic` | [boundary conditions](../../boundary-conditions#extrapolation-zeroth-order-constant), [periodic](../../boundary-conditions#periodic) |
| FNL backend | the same kernels on the device | [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

## Initial and boundary conditions

Every input is derived from a committed Euler input by `make_input.py`: `physical_model = mhd-ideal`, an `[mhd]`
section, $n_v = 8$ (`none`) or 9 (`glm`), and `bx, by, bz` added to every region and inflow state.

| Leg | Base input | Grid | State | Boundaries | Steps |
|---|---|---|---|---|---|
| Sod $x$, $y$, $z$ | `sod/sod-{x,y,z}.ini` | as [V1](../euler/sod): 64 blocks of $50 \times 4 \times 4$ cells along the active axis, one level | [`riemann-problem`](../../initial-conditions#riemann-problem-piecewise-constant-axis-aligned-regions): $(\rho, u, p) = (1, 0, 1)$ / $(0.125, 0, 0.1)$, $\mathbf{B} = 0$ | `extrapolation` | $t = 0.2$ (174 steps `none`, 241 steps `glm`: the $c_h$ bound sets $\Delta t$) |
| uniform, single | `conservation/amr-periodic.ini`, `markers_number = 0` | $[0, 1]^3$, octree level 1: 8 blocks of $8^3$ cells | [`uniform`](../../initial-conditions#uniform-seeded-uniform-state), $s = 0$: $(\rho, u, v, w, p) = (1, 0.3, -0.2, 0.1, 1)$, $\mathbf{B} = (0.8, 0.5, -0.3)$ | `periodic` | 10 |
| uniform, AMR | `conservation/amr-periodic.ini` | the octant $[0, 0.5]^3$ refined once more: 15 blocks, 2:1 faces on every side of the octant | the same | `periodic` | 10 |

Ranks: 2 in every leg.

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/mhd/zero-field/make_input.py` | derives the MHD input from an Euler one (`--divergence-control`, `--b`, `--set section.key=value`) |
| `src/tests/flume/verification/sod/sod-{x,y,z}.ini` | the Sod base inputs |
| `src/tests/flume/verification/conservation/amr-periodic.ini` | the V3 base input of the uniform legs |
| `src/tests/flume/verification/mhd/zero-field/check.sh` | runs every leg (`--np N`; `FLUME_EXE` selects the backend) |
| `src/tests/flume/verification/sod/sod_oracle.py` | $L_1(\rho)$ against the exact solution (`--l1-max 3.31e-03`), the $y$/$z$ against $x$ bitwise check |
| `src/tests/flume/verification/mhd/plumbing/mhd_plumbing_oracle.py` | `--zero` ($\mathbf{B}$, $\psi$ exactly zero), `--steady` (last checkpoint against step 0) |

```bash
cd src/tests/flume/verification/mhd/zero-field && ./check.sh         # CPU
FLUME_EXE=$PWD/../../../../../../exe/adam_flume_fnl ./check.sh       # FNL
```

## Results

A contract leg: no figure. Work directories of both backends are kept (`zero-field/work-adam_flume_{cpu,fnl}-np2-*`)
with logs and histories; the checkpoints are deleted after the oracles and the oracle lines are not kept, so the
cell-level verdicts below rest on the `check.sh` run. The values quoted are read from the kept histories:

| Check | Bound | CPU | FNL |
|---|---|---|---|
| Sod, $L_1(\rho)$, both variants | $3.31\cdot10^{-3}$ (V1) | not kept | not kept |
| Sod, $y$ and $z$ against $x$ | bitwise | not kept | not kept |
| Sod, $\mathbf{B}$ (and $\psi$) in the last checkpoint | exactly zero | not kept; $\int B_x, \int B_y, \int B_z$ (and $\int\psi$) are 0 in every history row, and the `glm` div(B) history is 0 at every step | the same |
| uniform, single level: last against step-0 checkpoint | bitwise | not kept; every volume integral constant to the printed digit, $\int\psi = 0$ | the same |
| uniform, 2:1 octant: last against step-0 checkpoint | $10^{-14}$ relative | not kept; integrals constant to the printed digit, $\max\lvert\int\psi\rvert = 1.2\cdot10^{-18}$ | integrals constant, $\max\lvert\int\psi\rvert = 2.0\cdot10^{-19}$ |

The regression case `uniform-amr-mhd` freezes the uniform AMR leg (GLM, 10 steps; `regression/run.sh`).

Related cases: [MHD plumbing](./plumbing) checks the time step, outputs and restart of the same model; [Sod shock
tube](../euler/sod) is the Euler reference; [Ryu–Jones 2a](./ryu-jones-2a) (MV-4) is the first test with a non-zero
field; [conservation across AMR seams](../euler/amr-conservation) (V3) is the base of the uniform legs and
[MHD conservation](./conservation) (MV-11) its MHD counterpart.

## Findings

- The Sod legs are **not** bitwise against the Euler path: at $\mathbf{B} = 0$ the MHD eigenspaces are degenerate and
  WENO acts in a different characteristic basis, so the expectation is the V1 accuracy bound, not equality.
- The uniform state across the 2:1 patch is exact only to round-off. The cause is the interface, not the MHD kernels:
  the coarse–fine ghost fill reproduces a constant only to round-off for general values, and the Euler path does the
  same (a uniform Euler state with $u, v, w = 0.8, 0.5, -0.3$ drifts by $1.3\cdot10^{-15}$ on this grid, while the
  committed all-ones state happens to be exact; `check.sh` header).
- The work directories `*-none-uniform` and `*-glm-uniform` (CPU, dated before the others) predate the single/AMR
  split; `check.sh` no longer writes them.
