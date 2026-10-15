# EGLM cleaning

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| EV-2, EV-3, EV-4 | Ideal MHD | M3 ([#47](https://github.com/szaghi/adam/issues/47), M3-P4) | CPU, FNL | yes | `orszag-tang-eglm` |

## Case and scope

The energy-consistent cleaning (`divergence_control = eglm`, Derigs et al. 2018) puts $\psi^2/2$ in the total energy
and adds nonconservative sources,

$$-(\nabla\cdot\mathbf{B})\,(0, \mathbf{B}, \mathbf{u}\cdot\mathbf{B}, \mathbf{u}, 0) - (\mathbf{u}\cdot\nabla\psi)\,(0, 0, \psi, 0, 1)$$

in the order $(\rho, \rho\mathbf{u}, E, \mathbf{B}, \psi)$. On a periodic box only the density integral is then
conserved; the momentum, $\mathbf{B}$ and energy integrals move by the time integral of their sources. The legs run
through the `--divergence-control eglm` option of the glm-pulse, Orszag–Tang and linear-wave scripts, and
`rj2a/check.sh` compares EGLM with GLM:

- **EV-2:** the GLM pulse against the telegraph solution, which EGLM shares with $\psi$ in field units (the oracle
  compares $c_h\psi$). It catches a wrong EGLM cleaning block or damping;
- **EV-3:** Ryu–Jones 2a, EGLM against GLM. With $\psi = 0$ the EGLM energy coupling and sources vanish, so the runs
  must agree bitwise on the 8 shared variables: it catches a source or an energy term that does not vanish;
- **EV-4:** Orszag–Tang, the drifts bounded by their sources, and the 180° symmetry. It catches a source with a wrong
  magnitude or sign;
- MV-5 under EGLM: the linear-wave orders;
- the Balsara–Spicer blast without limiter and floors, the case EGLM was expected to help.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | ideal MHD with EGLM, $\psi$ in field units, $\psi^2/2$ in $E$ | [models](../../models#divergence-control-eglm); `mhd_eglm_conservative_to_auxiliary`, `mhd_eglm_eigenvectors`, `mhd_eglm_flux`, [`adam_flume_mhd_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_library.F90) |
| Space, splitting | WENO5 flux splitting on characteristic variables | `mhd_eglm_face_split_fluxes`, [`adam_flume_mhd_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_library.F90) |
| Space, `weno-riemann` | HLLD, characteristic interpolation (and primitive on the blast) | `mhd_eglm_riemann_hlld`, [`adam_flume_mhd_riemann_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_riemann_library.F90) |
| Sources and damping | the nonconservative sources above; the damped cleaning energy becomes heat | `add_eglm_sources`, `add_glm_damping`, [`adam_flume_cpu_mhd_eglm_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_mhd_eglm_kernels.F90) |
| Time, boundaries | as in each verification (SSP-33 for the pulse and Ryu–Jones 2a, SSP-54 for Orszag–Tang, the waves and the blast) | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| FNL backend | the same kernels on the device | [`adam_flume_fnl_mhd_eglm_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_mhd_eglm_kernels.F90), [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

## Initial and boundary conditions

| Leg | Setup |
|---|---|
| EV-2 | the MV-3 legs of [GLM cleaning](./glm-cleaning) (d'Alembert and damped at $N = 64, 128, 256$, `min-cell`, walls) at pulse amplitude $10^{-3}$ instead of 0.1 |
| EV-3 | [Ryu–Jones 2a](./ryu-jones-2a), $N = 256$, $x$, `glm_ch = 2.0`, `glm_alpha = 0.18` |
| EV-4 | [Orszag–Tang](./orszag-tang), $128^2$, $t = 0.5$, a checkpoint every 25 steps (the oracle integrates the $\psi$ source from them) |
| MV-5 under EGLM | the MV-5 ladders of [order of accuracy](./order-of-accuracy), $c_h = 1.1(\lvert u_n\rvert + c_f)$, no damping |
| Blast | Balsara–Spicer, $\beta = 2.5\cdot10^{-4}$, $128^2$ (the [positivity limiter](./positivity-limiter) PV-1 setup), no limiter, no floors |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/mhd/glm-pulse/check.sh --divergence-control eglm` | EV-2 |
| `src/tests/flume/verification/mhd/rj2a/check.sh` | EV-3 (the EGLM leg is part of the default run) |
| `src/tests/flume/verification/mhd/orszag-tang/check.sh --divergence-control eglm` | EV-4 and the symmetry |
| `src/tests/flume/verification/mhd/eglm_conservation_oracle.py` | the EV-4 drift bounds |
| `src/tests/flume/verification/mhd/linear-wave/check.sh --divergence-control eglm` | MV-5 under EGLM |
| `src/tests/flume/verification/mhd/positivity-probe/make_blast.py --eglm` | the blast input |
| `src/tests/flume/verification/make_doc_figures.py` `fig_eglm_energy` | draws `eglm-energy.png` from the conservation histories of the `orszag-tang` work directories (`--only eglm-energy`) |

```bash
cd src/tests/flume/verification/mhd
glm-pulse/check.sh --divergence-control eglm                                  # EV-2
orszag-tang/check.sh --divergence-control eglm [--numerics hlld]              # EV-4, symmetry
linear-wave/check.sh --divergence-control eglm [--numerics hlld]              # MV-5 under EGLM
rj2a/check.sh [--numerics hlld]                                               # EV-3 (EGLM leg)
```

## Results

![Orszag-Tang, drift of the total energy under GLM and EGLM](/flume/eglm-energy.png)

*The total energy of the periodic Orszag–Tang box, $128^2$: GLM keeps it to round-off; EGLM changes it through its
nonconservative sources, by $3.0\cdot10^{-4}$ (splitting) and $4.4\cdot10^{-4}$ (HLLD) at $t = 0.5$, within the EV-4
bound.*

Measured in issue #47, M3-P4; CPU and FNL agree to the printed digits unless a cell says otherwise:

| Check | Bound | Splitting | HLLD characteristic |
|---|---|---|---|
| EV-2, pulse against the telegraph solution (amplitude $10^{-3}$): $L_1$ at 256 cells, d'Alembert / damped | $6.484\cdot10^{-9}$ / $1.406\cdot10^{-9}$ (CPU + 2%) | $6.36\cdot10^{-9}$ / $1.38\cdot10^{-9}$ | — |
| EV-2, $\int B_x$ drift | $10^{-14}$ | $3.7\cdot10^{-16}$ ($N = 128$) | — |
| EV-3, Ryu–Jones 2a, EGLM against GLM | bitwise (splitting); $10^{-11}$ (`--numerics`) | bitwise, $\psi = 0$ | $6.3\cdot10^{-13}$ (FNL: bitwise) |
| EV-4, Orszag–Tang: $\int\rho$ drift; $\int E$ drift / its source bound | round-off; $\le 1$ | $3\cdot10^{-17}$; $2.8\cdot10^{-3}$ | 0; $7.0\cdot10^{-3}$ |
| Orszag–Tang symmetry at step 100 / at $t = 0.5$ | $10^{-8}$ / $10^{-4}$ | $5.4\cdot10^{-13}$ / $4.3\cdot10^{-7}$ | $1.2\cdot10^{-12}$ / $1.1\cdot10^{-11}$ |
| linear waves, observed order | $\ge 4.7$ | 4.96–4.98 | 4.96–4.97 |
| Balsara–Spicer blast ($\beta = 2.5\cdot10^{-4}$, $128^2$, no limiter, no floors) | — | non-positive pressure at step 3 | reaches $t = 0.01$ (396 steps) |

The EV-2 legs at 64 and 128 cells, `min-cell` and walls have their own bounds (CPU M3-P4c + 2%:
$6.686\cdot10^{-7}$, $5.789\cdot10^{-8}$; damped $1.487\cdot10^{-7}$, $1.266\cdot10^{-8}$; `min-cell`
$1.201\cdot10^{-9}$; walls $9.021\cdot10^{-8}$). The linear-wave bounds under EGLM are per run, for the splitting
scheme and for `hlld:characteristic`.

Related cases: [GLM cleaning](./glm-cleaning) is the GLM version of EV-2; [Orszag–Tang](./orszag-tang) the GLM version
of EV-4; [positivity limiter](./positivity-limiter) makes the splitting scheme run the blast;
[G3b](../seams-studies/eglm-reflux) handles the EGLM reflux at 2:1 seams.

## Findings

- EV-2 runs at amplitude $10^{-3}$ because EGLM couples the $(B_x, \psi)$ pair to the fluid: the damped cleaning
  energy becomes heat, the heat drives a flow, and the source $-(\nabla\cdot\mathbf{B})\mathbf{u}$ moves $B_x$. At the
  GLM amplitude 0.1 this floors the damped legs at a relative $L_1$ of $8\cdot10^{-5}$ (against $7\cdot10^{-6}$ for
  GLM); at $10^{-3}$ EGLM equals GLM to 3 digits (relative $6.6\cdot10^{-6}$ both, on the damped leg at 256), so the
  telegraph oracle is exact for EGLM there. The $\int B_x$ drift scales as the amplitude cubed.
- EV-4 asserts $\int\rho$ to round-off and bounds every other drift, at every step, by twice the time integral of its
  source magnitude: $\max\lvert\mathbf{B}\rvert\int\|\nabla\cdot\mathbf{B}\|_1$ (momentum),
  $\max\lvert\mathbf{u}\rvert\int\|\nabla\cdot\mathbf{B}\|_1$ ($\mathbf{B}$),
  $\max\lvert\mathbf{u}\cdot\mathbf{B}\rvert\int\|\nabla\cdot\mathbf{B}\|_1 + \int\|\psi\,\mathbf{u}\cdot\nabla\psi\|_1$
  (energy). It is a necessary bound, not a tight one: the signed sources largely cancel (energy drift $2.8\cdot10^{-3}$
  of its bound; $\rho\mathbf{u}$ and $\mathbf{B}$ $\sim10^{-11}$, where the 180° symmetry cancels their sources).
- The symmetry defect is round-off amplified by the flow, more by the EGLM sources on the splitting scheme
  ($4.3\cdot10^{-7}$ against $5.9\cdot10^{-10}$ for GLM at $t = 0.5$), so it is checked at step 100 within $10^{-8}$
  and at the end within $10^{-4}$, as for the [rotor](./rotor).
- On the blast HLLD with primitive interpolation also reaches $t = 0.01$; the splitting scheme needs the
  [positivity limiter](./positivity-limiter) (M3-P5).
- `--amr` is not combined with EGLM in `orszag-tang/check.sh`: the EGLM sources at 2:1 seams are out of scope there.
