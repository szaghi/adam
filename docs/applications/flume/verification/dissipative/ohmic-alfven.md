# Ohmic decay and the visco-resistive Alfvén wave

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| VV-5 | Dissipative terms (M4) | M4 ([#65](https://github.com/szaghi/adam/issues/65)) | CPU, FNL | yes | `alfven-resistive` |

## Case and scope

A transverse field $b_t = A\sin(\mathbf k\cdot\mathbf x)$ on a gas at rest decays by Ohmic diffusion at the rate
$\eta k^2$; with a normal field $B_n$ the same mode is an Alfvén wave damped by viscosity and resistivity together.
The oracle (`viscous/waves.py`) integrates the $2\times2$ linear system for $(v_t, b_t)$ exactly; a fourth leg reruns
the acoustic wave of [VV-2](./viscous-waves) on the MHD model with $\mathbf B = 0$.

What the case catches. The Ohmic flux has its own kernels, a curl form whose cross terms an axis-aligned wave never
reads, and a sign that the plan in #65 printed wrong (an anti-diffusion; see
[numerics](../../numerics#dissipative-fluxes-navier–stokes)). Each defect is a number the exact linear mode pins:
- a **sign error** grows the field;
- a **missing cross term** changes the diagonal rate;
- a **coefficient error** changes every rate;
- a **wrong temperature on MHD** (a wrong slot, or a wrong $T$ read from the MHD auxiliary state by the viscous and
  heat fluxes) shows on the acoustic leg on MHD and nowhere else: the constant laws of the Alfvén leg never read $T$.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | ideal MHD (`mhd-ideal`), $\gamma = 1.4$, with resistivity $\eta$ and, on two legs, viscosity and conductivity | [models](../../models#ideal-mhd) |
| Divergence control | none (`mhd-none`) on the Ohmic and acoustic legs; GLM, $c_h = 3$, $\alpha = 0$ (`mhd-glm`) on the Alfvén leg | [GLM](../../models#divergence-control-mixed-glm); [`adam_flume_cpu_mhd_glm_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_mhd_glm_kernels.F90) |
| Ohmic flux | curl form, `dissipative_order` 4 or 2 | `add_resistive_fluxes_o4`, `add_resistive_fluxes_o2`, [`adam_flume_cpu_dissipation_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_dissipation_kernels.F90); `compute_resistive_flux`, [`adam_flume_dissipation_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_dissipation_library.F90); [numerics](../../numerics#dissipative-fluxes-navier–stokes) |
| Viscous and heat fluxes | as in VV-1, the temperature from the MHD auxiliary state | `add_dissipative_fluxes_o4`, [`adam_flume_cpu_dissipation_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_dissipation_kernels.F90) |
| Inviscid fluxes | WENO5 flux splitting, characteristic variables | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Time | `runge-kutta-ssp-54`, CFL 0.2 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| FNL backend | the same fluxes as device kernels | `add_resistive_fluxes_o4_dev`, [`adam_flume_fnl_dissipation_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_dissipation_kernels.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | $[0, 1]^2$, $z$ null |
| Grid | quadtree, $4 \times 4$ blocks, $N/4$ cells per block along the active axes; an axis-aligned magnetic wave keeps the transverse axis active and thin (4 blocks of 6 cells), the acoustic wave on MHD a null transverse axis (1 cell) |
| Ranks | 2 |
| Initial condition | [`sine-wave`](../../initial-conditions#sine-wave-a-small-sine-wave-on-a-uniform-state-verification) in its `magnetic` mode ($b_t = A\sin(\mathbf k\cdot\mathbf x)$, $A = 10^{-5}$) or `acoustic` mode, $\rho_0 = p_0 = 1$, $u = 0$; wavelength 1 along an axis, $1/\sqrt2$ along the diagonal |
| Boundaries | `periodic` on the $x$ and $y$ faces, `extrapolation` on the $z$ faces |

| Leg | Model | Direction | `dissipative_order` | Coefficients | $\mathbf B$ background | $N$ | $t$ |
|---|---|---|---|---|---|---|---|
| Ohmic x | `mhd-none` | $x$ | 4 and 2 | $\eta = 0.01$ | 0 | 32 / 64 / 128 | 0.5 |
| Ohmic diagonal | `mhd-none` | 45° (the cross terms $\partial_i B_d$ of the curl form) | 4 | $\eta = 0.01$ | 0 | 24 / 48 / 96 | 0.1 |
| Alfvén x | `mhd-glm` | $x$ | 4 | $\mu = 0.01$, $\eta = 0.005$ | $B_n = 1$ | 32 / 64 / 128 | 0.5 |
| acoustic on MHD | `mhd-none` | $x$ | 4 | $\mu = 0.01$, $k = 0.02$; odd part of the $\pm A$ twins | 0 | 32 / 64 / 128 | 0.5 |

The rates the oracle holds: Ohmic $\eta k^2$; Alfvén damping $(\nu + \eta)k^2/2$.

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/viscous/waves.py` | `make --mode magnetic` writes each input; `oracle` integrates the linear $(v_t, b_t)$ system exactly |
| `src/tests/flume/verification/viscous/check.sh --leg vv5` | runs the five ladders (`--np N`; `FLUME_EXE` selects the backend) |
| `src/tests/flume/regression/alfven-resistive/input.ini` | the frozen visco-resistive Alfvén wave at 45°, MHD-GLM, $B_n = 1$, $\mu = 0.01$, $\eta = 0.005$, $32^2$, $t = 0.05$ |
| `make_doc_figures.py` (`fig_viscous`) | the figure, from the CPU work directories of `viscous/check.sh` |

```bash
cd src/tests/flume/verification/viscous && ./check.sh --leg vv5             # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh --leg vv5       # FNL
```

## Results

![VV-1 to VV-5: convergence, and Becker's shock against the exact profile](/flume/viscous.png)

*Middle panel: the five VV-5 ladders, each normalised by its coarsest error, against the order-2 and order-4 slopes.*

Oracle: as VV-1, the observed order of the finest pair $\ge 3.8$ at order 4, $\ge 1.9$ at order 2.

| Leg | N | RMS error (CPU) | Order (CPU) | FNL |
|---|---|---|---|---|
| Ohmic x, order 4 | 32 / 64 / 128 | $3.64\cdot10^{-11}$ / $2.29\cdot10^{-12}$ / $1.43\cdot10^{-13}$ | +3.99, +4.00 | the same |
| Ohmic x, order 2 | 32 / 64 / 128 | $3.68\cdot10^{-9}$ / $9.20\cdot10^{-10}$ / $2.30\cdot10^{-10}$ | +2.00, +2.00 | the same |
| Ohmic diagonal, order 4 | 24 / 48 / 96 | $1.14\cdot10^{-10}$ / $7.19\cdot10^{-12}$ / $4.50\cdot10^{-13}$ | +3.99, +4.00 | the same |
| Alfvén x, order 4 | 32 / 64 / 128 | $6.37\cdot10^{-11}$ / $1.10\cdot10^{-12}$ / $2.18\cdot10^{-14}$ | +5.85, +5.66 | the same |
| acoustic on MHD, order 4 | 32 / 64 / 128 | $8.20\cdot10^{-11}$ / $1.69\cdot10^{-12}$ / $2.19\cdot10^{-14}$ | +5.60, +6.27 | $2.21\cdot10^{-14}$ at 128 cells, +6.26 |

CPU and FNL on 2 ranks, identical to the printed digits except the round-off of the acoustic leg's odd part.

| Quantity | Measured | Reference |
|---|---|---|
| Ohmic decay | $0.3947841$ | $\eta k^2 = 0.3947842$ |
| Alfvén damping | $0.29610$ | $(\nu + \eta)k^2/2 = 0.29609$ |
| acoustic amplitude on MHD, 32 cells (CPU) | $8.6252132315\cdot10^{-6}$ | VV-2 on Euler, CPU: $8.6252132316\cdot10^{-6}$ (FNL: $8.6252132321\cdot10^{-6}$) |

Related cases: the same oracle on Euler is [VV-1 and VV-2](./viscous-waves); the energy budget and the divergence of
the same resistive flux with GLM and EGLM are in [VV-6](./ohmic-heating); the diagonal Ohmic wave crosses 2:1 seams in
[VV-8](./seam-accuracy).

## Findings

- The measured Ohmic rate is the derived one: the sign printed in the #65 plan was wrong, the code's is right.
- Along $x$ the Ohmic decay of $b_y$ is the same equation as the shear wave's $u_y$ in VV-1, and its errors are the
  VV-1 errors scaled by the amplitude, to three digits: two kernels, one stencil.
- The Alfvén wave's error is dominated by the inviscid WENO part, as on the acoustic wave, hence the order above four.
- The acoustic wave on MHD reproduces VV-2 on Euler: the measured amplitude agrees to ten digits at 32 cells; the RMS
  errors agree to three, the inviscid parts of the two models being different kernels.
- **The diagonal runs without divergence control.** With GLM or EGLM the same ladder measured +7.65 then +1.72, and
  +3.32 from 96 to 192 cells. The cause is the ideal scheme: run with $\eta = 0$, where the exact field is static,
  GLM's upwinding ($c_h = 3$) damps it at fifth order ($8.5\cdot10^{-12}$ at 48 cells, $2.7\cdot10^{-13}$ at 96), with
  the opposite sign to the Ohmic error and as large, so the two cancel near 48 cells and a 24/48/96 ladder is not
  asymptotic. Without divergence control the $\eta = 0$ field stays at round-off. A smaller time step (CFL 0.05)
  changed no printed digit, so neither error is temporal. GLM and EGLM with resistivity are covered by
  [VV-6](./ohmic-heating).
- The magnetic amplitude is $10^{-5}$: the wave's magnetic pressure drives an $O(A^2)$ flow that advects it at
  $O(A^3)$, odd in $A$ (no twin cancels it), relative $A^2 = 10^{-10}$ here.
- #85 asks for the decay against the exact rate as a figure; the published figure shows the error ladders only.
