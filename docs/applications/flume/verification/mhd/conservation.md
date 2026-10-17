# MHD conservation across AMR seams

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| MV-11 | Ideal MHD | M2 ([#41](https://github.com/szaghi/adam/issues/41), P7a) | CPU, FNL | yes | — |

## Case and scope

The V3 case of [conservation across AMR seams](../euler/amr-conservation) (a triply periodic box, one refined octant,
6 coarse–fine faces, a uniform flow with a 1% seeded density/pressure perturbation) run through the MHD path with the
uniform field $\mathbf{B} = (0.8, 0.5, -0.3)$. The perturbation drives the induction, so the flux register carries the
$\mathbf{B}$ and $\psi$ rows at every seam. Nothing crosses the periodic boundary, so every volume integral of the
conservative state must stay constant to round-off when the Berger–Colella reflux is on.

Four legs:

1. **GLM without damping** (`glm_alpha = 0`), reflux on: all 9 integrals ($\rho$, $\rho\mathbf{u}$, $E$, each
   component of $\mathbf{B}$, $\psi$) within $10^{-13}$; $\psi$ starts at zero, so its drift is relative to the
   largest initial integral;
2. **the same, reflux off:** the negative control must drift by at least $10^{-10}$, so the seams are exercised;
3. **GLM with damping** (`glm_alpha = 0.18`), reflux on: the 8 physical integrals within $10^{-13}$; $\int\psi$ is
   reported, not bounded;
4. **no divergence control** (`none`, $n_v = 8$), reflux on: the 8 integrals within $10^{-13}$.

Every run must register the 6 coarse–fine faces (`registered intra-realm AMR seam faces: +6` in the log).

What the case catches: a reflux register without the $\mathbf{B}$ or $\psi$ rows, a wrong seam flux of the MHD or GLM
kernels, a vacuous test (the negative control), and a seam layout that registers fewer faces than the six of the
octant.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | ideal MHD, $\gamma = c_p/c_v = 1.4$ (the V3 input's); `none` ($n_v = 8$) or mixed GLM ($n_v = 9$) | [models](../../models#ideal-mhd) |
| Space | WENO5 (`weno-u-5`) flux splitting on characteristic variables | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `mhd_face_split_fluxes`, `mhd_glm_face_split_fluxes`, [`adam_flume_mhd_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_library.F90) |
| Cleaning | GLM, $c_h = 3$, `glm_damping_length = 1`, `glm_alpha = 0` or $0.18$, `glm_ch_check = warning`; or `none` | [models](../../models#divergence-control-mixed-glm); `add_glm_damping`, [`adam_flume_cpu_mhd_glm_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_mhd_glm_kernels.F90) |
| Time | SSP Runge–Kutta, 3 stages (`runge-kutta-ssp-33`), CFL 0.5 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| AMR | geometric box marker, one more level on the octant; 2:1 ghost restriction and interpolation | `mark_by_geometry`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90); `interp_seam_ghost`, [`adam_field_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_field_object.F90); [boundary conditions](../../boundary-conditions#_2-1-amr-coarse-fine-ghosts) |
| Reflux | Berger–Colella, stage fluxes weighted by the SSP coefficients, every conservative row | `accumulate_seam_fluxes`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `accumulate_coarse_flux`, `accumulate_fine_flux`, [`adam_flux_register_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_flux_register_object.F90); `apply_reflux_corrections`, [`adam_forest_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_forest_object.F90); [numerics](../../numerics#adaptive-mesh-refinement) |
| Integrals | volume integrals of the conservative variables, every step | `compute_conservation`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Boundaries | periodic on all six faces | [boundary conditions](../../boundary-conditions#periodic) |
| FNL backend | the same kernels on the device | [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Base input | `conservation/amr-periodic.ini` (V3), turned into MHD by `mhd/zero-field/make_input.py` |
| Domain | $[0, 1]^3$ |
| Grid | octree, one uniform level (8 blocks of $8^3$ cells), the octant $[0, 0.5]^3$ refined once more: 15 blocks, 6 coarse–fine faces, `ngc = 3` |
| Ranks | 2 |
| Initial condition | [`uniform`](../../initial-conditions#uniform-seeded-uniform-state): $(\rho, u, v, w, p) = (1, 1, 1, 1, 1)$, $\mathbf{B} = (0.8, 0.5, -0.3)$, $\psi = 0$, seeded density/pressure perturbation $s = 0.01$ |
| Boundaries | `periodic` on all six faces |
| Final time | 20 steps ($t = 0.0347$) |

| Leg | `divergence_control` | `glm_alpha` | `reflux` | Bound |
|---|---|---|---|---|
| `glm-reflux` | `glm` | 0 | on | all 9 integrals $\le 10^{-13}$ |
| `glm-noreflux` | `glm` | 0 | off | max drift $\ge 10^{-10}$ (negative control) |
| `glm-damped-reflux` | `glm` | 0.18 | on | the 8 physical integrals $\le 10^{-13}$; $\psi$ reported |
| `none-reflux` | `none` | — | on | the 8 integrals $\le 10^{-13}$ |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/conservation/amr-periodic.ini` | the V3 base input |
| `src/tests/flume/verification/mhd/zero-field/make_input.py` | derives the MHD input (`--divergence-control`, `--b 0.8 0.5 -0.3`, `--set mhd.glm_alpha=…`, `--set numerics.reflux=…`) |
| `src/tests/flume/verification/mhd/conservation/check.sh` | runs the four legs and checks the six registered faces (`--np N`, `--legs LIST`: `mv11-glm` for legs 1 and 2, `mv11-damped`, `mv11-none`; `FLUME_EXE` selects the backend) |
| `src/tests/flume/verification/conservation/conservation_oracle.py` | the drifts (`--conserved`, `--max-drift`, `--leaky`, `--min-drift`, `--exclude psi`) |
| `src/tests/flume/verification/make_doc_figures.py` (`fig_mhd_conservation`) | the figure, from the kept CPU conservation histories of `mhd/conservation/` (`--only mhd-conservation`) |

```bash
cd src/tests/flume/verification/mhd/conservation && ./check.sh        # CPU
FLUME_EXE=$PWD/../../../../../../exe/adam_flume_fnl ./check.sh        # FNL
cd ../.. && python3 make_doc_figures.py --out ../../../../docs/public/flume --only mhd-conservation
```

## Results

![MV-11: the drift of the MHD integrals across 2:1 seams](/flume/mhd-conservation.png)

*One panel per CPU leg (GLM with reflux, damped GLM with reflux, GLM without reflux, no cleaning with reflux):
$\lvert I(t) - I(0)\rvert / \lvert I(0)\rvert$ of every volume integral against the step ($\psi$ absolute), floored
at $10^{-17}$.*

Maximum drift over the 20 steps, computed from the kept conservation histories
(`mhd/conservation/work-adam_flume_{cpu,fnl}-np2-<leg>/amr-periodic-conservation_history.dat`) with the oracle's
definition (relative to $\lvert I(0)\rvert$; $\psi$ relative to the largest initial integral, $\int E = 4.490$):

| Leg | Bound | CPU | FNL |
|---|---|---|---|
| `glm-reflux`, 8 physical integrals | $10^{-13}$ | 0 (constant to the printed digit) | $\le 3.7\cdot10^{-16}$ ($B_z$) |
| `glm-reflux`, $\psi$ | $10^{-13}$ | $3.8\cdot10^{-20}$ | $9.5\cdot10^{-20}$ |
| `glm-noreflux`, max over the integrals | $\ge 10^{-10}$ | $1.69\cdot10^{-5}$ ($\rho u$); $1.5\cdot10^{-6}$ to $1.7\cdot10^{-5}$ over the 9 | the same to 3 digits |
| `glm-damped-reflux`, 8 physical integrals | $10^{-13}$ | $\le 2.2\cdot10^{-16}$ | $\le 3.7\cdot10^{-16}$ |
| `glm-damped-reflux`, $\psi$ (reported) | — | $6.0\cdot10^{-11}$ ($\lvert\int\psi\rvert = 2.70\cdot10^{-10}$) | the same |
| `none-reflux`, 8 integrals | $10^{-13}$ | 0 | $\le 3.7\cdot10^{-16}$ |
| registered coarse–fine faces, every leg, both ranks | 6 | 6 | 6 |

Related cases: [conservation across AMR seams](../euler/amr-conservation) is the Euler V3 base; [zero
field](./zero-field) (MV-2) runs the same box with a uniform MHD state; [GLM divergence cleaning](./glm-cleaning) checks
the damping; the M6 integral-drift study near seams is [integral drift](../seams/integral-drift).

## Findings

- With reflux the MHD integrals are constant to round-off on both backends, $\psi$ included: the register carries the
  $\mathbf{B}$ and $\psi$ rows.
- With damping, $\int\psi$ is not conserved: the source $-k\psi$ couples the Runge–Kutta stages. The seam-flux
  mismatch of stage $s$ reaches $\psi^{n+1}$ through the later stages' damping, so it enters with weight
  $b_s(1 - O(k\,\Delta t))$ while the register corrects $b_s$; the remainder is $O(k\,\Delta t)$ of the uncorrected
  leak. Measured at M2-P7a (`check.sh` header): $\lvert\int\psi\rvert = 2.7\cdot10^{-11}$, $2.7\cdot10^{-10}$,
  $1.3\cdot10^{-9}$ for `glm_alpha` 0.018, 0.18, 0.9, linear in the rate; without the seams it stays at $10^{-21}$.
  The kept damped histories give $2.70\cdot10^{-10}$ at 0.18, as recorded.
- The `check.sh` header records the M2-P7a CPU drift with reflux as $\le 4.7\cdot10^{-14}$ ($B_z$ the largest) and
  $7\cdot10^{-6}$ to $1.7\cdot10^{-5}$ without. The kept CPU histories show exact constancy with reflux and
  $1.5\cdot10^{-6}$ ($B_x$) to $1.7\cdot10^{-5}$ without; the $7\cdot10^{-6}$ lower end matches $E$
  ($7.07\cdot10^{-6}$), not the $\mathbf{B}$ integrals.
- The bounds hold on FNL as well. The FNL drifts with reflux are one to two ulps, not zero, and the FNL step-0
  integrals already differ from the CPU ones in the last digit ($\int B_z = -0.29999999999999993$ against
  $-0.29999999999999999$): the two backends do not compute the integrals with the same arithmetic.
