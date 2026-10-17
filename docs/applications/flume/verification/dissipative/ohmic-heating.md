# Ohmic heating and the divergence constraint

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| VV-6 | Dissipative terms (M4) | M4 ([#65](https://github.com/szaghi/adam/issues/65)) | CPU, FNL | no | — |

## Case and scope

The diagonal Ohmic wave of [VV-5](./ohmic-alfven) with GLM and with EGLM, judged by its energy budget and by its
divergence instead of its error.

What the case catches:
- **a wrong Ohmic energy flux:** $G^E_d = \sum_i B_i\,G^{B_i}_d$ must return to the gas, as heat, exactly the energy
  the field loses. A wrong (or wrong-signed) energy flux leaves the field equations right and shows only in this
  budget; the total energy, conserved by the flux form whatever the fluxes, would not see it either;
- **a resistive flux that feeds $\nabla\cdot\mathbf B$:** the curl form must not raise the constraint GLM and EGLM are
  there to keep.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | ideal MHD (`mhd-ideal`), $\gamma = 1.4$, $\eta = 0.01$, $\mu = 0$ | [models](../../models#ideal-mhd) |
| Divergence control | GLM (`mhd-glm`) and EGLM (`mhd-eglm`), $c_h = 3$, $\alpha = 0$ | [GLM](../../models#divergence-control-mixed-glm), [EGLM](../../models#divergence-control-eglm); [`adam_flume_cpu_mhd_glm_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_mhd_glm_kernels.F90), [`adam_flume_cpu_mhd_eglm_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_mhd_eglm_kernels.F90) |
| Ohmic flux | curl form with its energy flux, `dissipative_order = 4` | `add_resistive_fluxes_o4`, [`adam_flume_cpu_dissipation_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_dissipation_kernels.F90); `compute_resistive_flux`, [`adam_flume_dissipation_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_dissipation_library.F90); [numerics](../../numerics#dissipative-fluxes-navier–stokes) |
| Inviscid fluxes | WENO5 flux splitting, characteristic variables | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Time | `runge-kutta-ssp-54`, CFL 0.2 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| FNL backend | the same fluxes as device kernels | `add_resistive_fluxes_o4_dev`, [`adam_flume_fnl_dissipation_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_dissipation_kernels.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | $[0, 1]^2$, $z$ null |
| Grid | quadtree, $4 \times 4$ blocks of $12^2$ cells ($48^2$), uniform |
| Ranks | 2 |
| Initial condition | [`sine-wave`](../../initial-conditions#sine-wave-a-small-sine-wave-on-a-uniform-state-verification), `magnetic` mode at 45°, $A = 10^{-4}$ (the budget compares energies of $O(A^2)$ with differences of the internal energy, so it needs the signal), $\rho_0 = p_0 = 1$, wavelength $1/\sqrt2$ |
| Boundaries | `periodic` on the $x$ and $y$ faces, `extrapolation` on the $z$ faces |
| Final time | $t = 0.2$ |
| Runs | `mhd-glm` and `mhd-eglm`, each with its ideal twin ($\eta = 0$) |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/viscous/waves.py` | `make` writes the inputs; `budget` judges the energy budget and $\nabla\cdot\mathbf B$ against the ideal twin |
| `src/tests/flume/verification/viscous/check.sh --legs vv6` | runs the four runs (`--np N`; `FLUME_EXE` selects the backend) |

```bash
cd src/tests/flume/verification/viscous && ./check.sh --legs vv6             # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh --legs vv6       # FNL
```

No figure (#85 lists none for VV-6).

## Results

Oracle (`waves.py budget`):
- heat gained over field and kinetic energy lost within $10^{-3}$ of 1;
- the drift of the total energy at most $10^{-13}$ of the total;
- $\max|\nabla\cdot\mathbf B|$ (central differences) at most 1.05 times the ideal twin's plus $10^{-14}$.

| Check | Bound | CPU | FNL |
|---|---|---|---|
| field energy lost (of $2.5\cdot10^{-9}$) | — | $6.770374\cdot10^{-10}$ | the same |
| internal energy gained | — | $6.770375\cdot10^{-10}$ | $6.770384\cdot10^{-10}$ |
| heat / (field + kinetic loss) | $1 \pm 10^{-3}$ | 1.0000002 | 1.0000015 |
| total energy drift, relative | $10^{-13}$ | $1.8\cdot10^{-16}$ | $1.8\cdot10^{-16}$ |
| $\max\lvert\nabla\cdot\mathbf B\rvert$, resistive / ideal (GLM; EGLM alike) | $1.05\times$ ideal $+ 10^{-14}$ | $2.1\cdot10^{-18}$ / $2.6\cdot10^{-18}$ | $1.6\cdot10^{-18}$ / $2.3\cdot10^{-18}$ (GLM) |

Both models give the same budget.

Related cases: the wave itself is the diagonal Ohmic leg of [VV-5](./ohmic-alfven); the divergence at 2:1 seams is
exercised by [VV-7](./seam-conservation) and [VV-8](./seam-accuracy).

## Findings

- The departures of the heat ratio from 1 are round-off: the internal energy is the total minus the kinetic and
  magnetic parts, and a difference of $6.8\cdot10^{-10}$ out of a total of 2.5 keeps six digits.
- **What the divergence leg does not show.** On a uniform periodic grid the diagonal field is discretely solenoidal
  from the start, and the curl-form flux keeps it so; both runs sit at round-off, so the leg checks only that the
  resistive flux adds no divergence there. The test with teeth is the 2:1 AMR seam, where the discrete operators no
  longer commute (PRISM's [#29](https://github.com/szaghi/adam/issues/29) floor): it belongs to M4 P4
  ([VV-7](./seam-conservation), [VV-8](./seam-accuracy)).
