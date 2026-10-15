# MHD rotor

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| MV-13 | Ideal MHD | M2 ([#41](https://github.com/szaghi/adam/issues/41)) | CPU, FNL | yes | `rotor` |

## Case and scope

A dense disc spins in a magnetised ambient gas (Balsara & Spicer 1999; the first rotor of Tóth 2000) and launches
torsional Alfvén waves that brake it, compresses the field and drives strong shocks: a stringent 2-D test of positivity
and of the multi-dimensional field evolution. There is no exact solution; the oracle checks:

- **symmetry:** the 180° rotation composed with $\mathbf{B} \to -\mathbf{B}$ (the ambient field is uniform) maps the
  solution to itself. A symmetry defect of the scheme is at truncation level ($10^{-3}$) from the first steps; a correct
  scheme leaves a round-off seed amplified by the flow;
- **positivity:** the floors are disabled, so a non-positive density or pressure stops the run.

There is no conservation check: the outflow boundaries do not conserve.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | ideal MHD, $\gamma = 1.4$ ($c_p = 3.5$, $c_v = 2.5$) | [models](../../models#ideal-mhd) |
| Space | WENO5 flux splitting on characteristic variables; `weno-riemann` under `--numerics` | `compute_residuals_weno`, `compute_residuals_riemann`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Divergence control | GLM, $c_h = 6$ (`glm_ch_check = error`), `glm_alpha = 0.18`, `glm_damping_length = 1.0` | `mhd_glm_face_split_fluxes`, [`adam_flume_mhd_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_library.F90); [models](../../models#divergence-control-mixed-glm) |
| Floors | disabled | [models](../../models#positivity-floors) |
| Time | SSP Runge–Kutta, 5 stages (`runge-kutta-ssp-54`), CFL 0.4 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Boundaries | zeroth-order extrapolation (outflow) | `set_boundary_conditions`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [boundary conditions](../../boundary-conditions#extrapolation-zeroth-order-constant) |
| FNL backend | the same kernels on the device | [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | $[0, 1]^2$, $z$ null |
| Grid | the vortex grid of [V2](../euler/isentropic-vortex) (quadtree, $4 \times 4$ blocks), $256^2$ cells ($64^2$ per block); a power of two, so the cell centres are exact binary fractions and the initial state is bitwise symmetric (Tóth uses $200^2$) |
| Ranks | 2 |
| Initial condition | [`mhd-rotor`](../../initial-conditions#mhd-rotor): ambient $\rho = 1$ at rest, $p = 1$, $\mathbf{B} = (5/\sqrt{4\pi}, 0, 0)$; disc $\rho = 10$ for $r < 0.1$ around $(0.5, 0.5)$, rigid rotation at speed 2 at $r = 0.1$, linear taper to $r = 0.115$ |
| Boundaries | `extrapolation` on $x$ and $y$ |
| Final time | $t = 0.15$ (1152 steps) |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/mhd/rotor/make_rotor.py` | writes the input from `vortex/vortex-n064.ini` (`--cells`, `--cfl`, `--time-max`, `--glm-ch`, `--it-save`) |
| `src/tests/flume/verification/mhd/rotor/check.sh` | runs the case (`--np N`, `--numerics SPEC`) |
| `src/tests/flume/verification/mhd/orszag-tang/orszag_tang_oracle.py` | the symmetry check, `--field-parity even` (the $\mathbf{B} \to -\mathbf{B}$ composition) |
| `src/tests/flume/verification/make_doc_figures.py` `fig_rotor` | draws `rotor.png` from `<runs>/rotor` (`--only rotor`) |

```bash
cd src/tests/flume/verification/mhd
rotor/check.sh                                                     # CPU (about 10 min)
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl rotor/check.sh    # FNL
rotor/check.sh --numerics hlld:primitive                           # weno-riemann, HLLD, primitive interpolation
```

`make_doc_runs.sh` has no group that writes `<runs>/rotor`: the figure run is not reproducible from the scripts alone
([#85](https://github.com/szaghi/adam/issues/85)).

## Results

![MHD rotor](/flume/rotor.png)

*Density and magnetic pressure $\lvert B\rvert^2/2$ of the rotor.*

| Check | Bound | CPU | FNL |
|---|---|---|---|
| symmetry (180° rotation with $\mathbf{B} \to -\mathbf{B}$), splitting, $t = 0.15$ | $10^{-8}$ | $1.8\cdot10^{-9}$ ($10^{-11}$ at step 100) | pass |
| positivity (floors disabled) | run completes | pass | pass |
| `--numerics`: symmetry at step 100 | $10^{-8}$ | $3.9\cdot10^{-12}$ HLLD primitive, $4.9\cdot10^{-10}$ HLLD characteristic | $2.1\cdot10^{-12}$, $7.9\cdot10^{-10}$ |
| `--numerics`: symmetry at the end | $10^{-4}$ (gross failures only) | $3.4\cdot10^{-11}$ HLL; $5.4\cdot10^{-9}$ HLLD primitive; $7.4\cdot10^{-8}$ HLLD primitive without the sensor; $1.1\cdot10^{-7}$ HLLD characteristic | $2.3\cdot10^{-7}$ and $6.1\cdot10^{-6}$ for the two HLLD |

$w$ and $B_z$ stay exactly 0 on the CPU.

Related cases: [Orszag–Tang](./orszag-tang) shares the symmetry oracle; [MHD on `weno-riemann`](./weno-riemann)
explains the two-point symmetry check; [multi-realm runs](../euler/multi-realm) split this rotor on a quadtree.

## Findings

- The final symmetry defect grows smoothly from round-off (the time series has no jump, so no switch flips) and orders
  by the dissipation of the solver and by the backend's round-off: the least dissipative interpolations amplify the
  seed most, and FNL contracts differently into FMAs.
- A single end-of-run tolerance cannot separate a round-off seed from a defect across solvers, so under `--numerics`
  the symmetry is asserted twice: at step 100, before the flow amplifies the seed, within $10^{-8}$, and at the end
  within $10^{-4}$.
