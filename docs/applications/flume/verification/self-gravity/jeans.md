# Jeans modes under self-gravity

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| SV-3, SV-4, SV-5 | Self-gravity | M8 ([#81](https://github.com/szaghi/adam/issues/81)), P3 | CPU, FNL | yes | `selfgravity-jeans` |

## Case and scope

P3 couples the FFT solver to the scheme: `[gravity] self_gravity`, one Poisson solve per Runge–Kutta stage, applied
through the potential path of [SV-1](./potential). A linear Jeans mode on a uniform gas obeys
$\omega^2 = c^2k^2 - 4\pi G\rho_0$: it oscillates for wavelengths below the Jeans length and grows above it.

What the case catches:
- **SV-3, the stable wave:** a wrong factor $4\pi G$, a wrong sign of the force, the mean not removed, or a solve that
  does not use each stage's own density; any of them shifts the frequency, and the L1 error against the exact wave
  leaves second order;
- **SV-4, the unstable mode:** a wrong growth rate; a broken point symmetry (an exact symmetry of the equations);
  a leak in the conservation the design promises, exact in the semi-discrete form, so that only the time integration
  may drift: $E + W$ must drift at the Runge–Kutta order, and the momentum must stay at round-off;
- **SV-5, scaling covariance:** a wrong dimension of $G$ in the `[reference]` layer.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | compressible Euler, ideal gas, $\gamma = 5/3$, with self-gravity | [models](../../models#self-gravity) |
| Gravity | `self_gravity = .true.`, $\Phi$ solved every Runge–Kutta stage from that stage's density | `solve_potential`, `set_potential`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90); `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `solve_potential_gpu`, [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |
| Poisson solver | FFT, periodic ([SV-2](./fft-poisson)): FFTW on CPU, cuFFT on FNL | [`adam_poisson_fft_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_poisson_fft_object.F90), [`adam_fnl_poisson_fft_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/fnl/adam_fnl_poisson_fft_object.F90) |
| Potential sources | the face form of SV-1 | `add_potential_sources`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Space | WENO5 (`weno-u-5`) flux splitting on characteristic variables (`weno`) | [numerics](../../numerics#scheme-space-weno-flux-splitting) |
| Time | SSP Runge–Kutta, 3 stages (`runge-kutta-ssp-33`), CFL 0.4 (collapse: 0.4, 0.2, 0.1) | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Scaling (SV-5) | `[reference]` converts $G$ by $u_0^2/(\rho_0L_0^2)$ | [`adam_flume_gravity_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_gravity_object.F90); [models](../../models#units-and-scaling) |
| Executable | the FFT build: `exe/adam_flume_cpu_fft` (mode `flume-cpu-gnu-fft`), `exe/adam_flume_fnl_fft` (mode `flume-fnl-nvf-fft`) | `fobos.d/flume.fobos` |

## Initial and boundary conditions

| Item | Value |
|---|---|
| State | $\rho_0 = 1$, $p_0 = 3/5$, $\gamma = 5/3$ (sound speed $c = 1$), at rest |
| Mode | [`sine-wave`](../../initial-conditions#sine-wave-a-small-sine-wave-on-a-uniform-state-verification) with `wave_mode = jeans`: $\mathbf k = 2\pi(1, 1)$ (`wave_angle` 45°, `wavelength` $1/\sqrt2$), $k^2 = 8\pi^2$; $\rho = \rho_0(1 + A\sin\theta)$, $p = p_0(1 + \gamma A\sin\theta)$, the velocity of the mode |
| Grid | the unit square (2-D quadtree, 4×4 blocks, two uniform levels) or the unit cube (3-D octree, 2×2×2 blocks, one level), uniform, 4 ghosts; $4n$ (2-D) or $2n$ (3-D) cells per side for $n$ cells per block side |
| Boundaries | every active face `periodic` (2-D: `extrapolation` on the null $z$ faces) |
| Ranks | 2 |

| Leg | $G$ | Physics | Grid | Amplitude, end |
|---|---|---|---|---|
| SV-3, `jeans2d` | $\pi/2$: $4\pi G\rho_0 = 2\pi^2$ | stable, $\lambda = \lambda_J/2$, $\omega = \pi\sqrt6$ | $32^2$, $64^2$, $128^2$ | $A = 10^{-6}$; one period, $t = 2/\sqrt6$ |
| SV-3, `jeans3d` | $\pi/2$ | the same wave in the $x$-$y$ plane, $z$ uniform | $16^3$, $32^3$ | $10^{-6}$; $t = 2/\sqrt6$ |
| SV-4, `jeansu2d` (linear) | $9\pi/2$: $4\pi G\rho_0 = 18\pi^2$ | unstable, $\lambda = 1.5\,\lambda_J$, $\sigma = \pi\sqrt{10}$; the box shifted to $[1/8, 9/8]^2$, so that the density's centre of symmetry is its corner | $64^2$, a checkpoint every 5 steps | $10^{-6}$; $t = 0.5$ |
| SV-4, `jeansc2d` (collapse) | $9\pi/2$ | as `jeansu2d`, into the nonlinear collapse (the peak density reaches 1.73) | $64^2$, CFL 0.4, 0.2, 0.1 | $10^{-2}$; $t = 0.4$ |
| SV-5 | $\pi/2$, and its dimensional twin | `jeans2d` in $L_0 = 4$, $u_0 = 1/2$, $\rho_0 = 1/16$ with `[reference]`, against the code-unit run | $32^2$ | $10^{-6}$; $t = 0.2$ |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/selfgravity/make_jeans.py` | writes `jeans2d`, `jeans3d`, `jeansu2d`, `jeansc2d` (`--cells`, `--cfl`, `--time-max`, `--it-save`, `--amplitude`) |
| `src/tests/flume/verification/selfgravity/jeans.sh` | runs the legs `sv3-2d`, `sv3-3d`, `sv4`, `sv5` (`--np N`, `--legs LIST`, `--keep`; `FLUME_EXE` selects the backend) |
| `src/tests/flume/verification/selfgravity/jeans_oracle.py` | `wave` (SV-3), `growth` and `symmetry` (SV-4 linear), `conservation` (SV-4 collapse) |
| `src/tests/flume/verification/scaling/reference.sh` | SV-5: the dimensional twin against the code-unit run, bitwise (NV-5) |
| `src/tests/flume/regression/selfgravity-jeans/input.ini` | the frozen `jeansc2d` at $32^2$, 20 steps |
| `src/tests/flume/verification/make_doc_runs.sh`, group `jeans` | `sv3-2d-{8,16,32}`, `sv3-3d-{8,16}`, `sv4-linear`, `sv4-collapse-{0.4,0.2,0.1}` |
| `src/tests/flume/verification/make_doc_figures.py`, `fig_selfgravity_jeans` | draws the figure (`--only selfgravity-jeans`) |

```bash
cd src/tests/flume/verification/selfgravity && ./jeans.sh --legs "sv3-2d sv3-3d sv4 sv5"        # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl_fft ./jeans.sh --legs "sv3-2d sv3-3d sv4 sv5"  # FNL
```

## Results

![SV-3 and SV-4: the stable wave's errors against N, the unstable mode's growth, the collapse's energy drift and density](/flume/selfgravity-jeans.png)

*Top left: the stable wave's error falls as $\Delta x^2$ in 2-D and 3-D. Top right: the unstable mode grows as
$e^{\sigma t}$, the fit on the exact line. Bottom left: in the collapse $E$ and $W$ exchange 17 % (their curves
coincide) while $E + W$ drifts at the Runge–Kutta order, a factor 8 per halving of the CFL number. Bottom right: the
collapsing density at $t = 0.4$, peak 1.73.*

The exact linear solution is $\rho_0(1 + A\sin(\theta - \omega t))$ (stable) and $\rho_0(1 + Ae^{\sigma t}\sin\theta)$
(unstable). The growth rate is the least-squares fit of $\log a(t)$, $a = 2\langle(\rho/\rho_0 - 1)\sin\theta\rangle$,
over the checkpoints in $[0.1, 0.5]$. The symmetry is the point reflection about the box corner: $\rho$ even,
$\mathbf u$ odd.

| Leg | Asserted | Bound | CPU | FNL |
|---|---|---|---|---|
| SV-3, 2-D, $N$ = 32, 64, 128 | L1 error of $\rho$ over $\rho_0A$: $2.14\cdot10^{-3}$, $5.35\cdot10^{-4}$, $1.34\cdot10^{-4}$ (CPU = FNL to 4 digits); order of the finest pair | $\ge 1.8$ | 2.00, 2.00 | 2.00, 2.00 |
| SV-3, 3-D, $N$ = 16, 32 | $8.49\cdot10^{-3}$, $2.14\cdot10^{-3}$; order | $\ge 1.8$ | 1.99 | 1.99 |
| SV-4, growth ($N$ = 64, 27 checkpoints in $[0.1, 0.5]$) | the fitted $\sigma$ = 9.9275 against $\pi\sqrt{10}$ = 9.9346, relative | $10^{-2}$ | $7.2\cdot10^{-4}$ | $7.2\cdot10^{-4}$ |
| SV-4, mirror symmetry | defect over the mode amplitude at $t = 0.5$ | $10^{-9}$ | $\rho$ $4.8\cdot10^{-11}$, $\mathbf u$ $3.7\cdot10^{-11}$ | $\rho$ $9.8\cdot10^{-11}$, $\mathbf u$ $7.7\cdot10^{-11}$ |
| SV-4, collapse, CFL 0.4, 0.2, 0.1 | drift of $E + W$: $2.74\cdot10^{-7}$, $3.45\cdot10^{-8}$, $4.33\cdot10^{-9}$; order in $\Delta t$ | $\ge 2.5$ | 2.99, 2.99 | 2.99, 2.99 (the same drifts) |
| SV-4, collapse | the largest $\lvert\int\rho\mathbf u\rvert$ | $10^{-14}$ | $8.9\cdot10^{-18}$ | $5.1\cdot10^{-18}$ |
| SV-5, scaling covariance | the dimensional run (with $G$ converted) against the code-unit run | bitwise | bitwise (fields, histories, 14 conversions) | bitwise with nvfortran 26.1 (WSL GPU); nvfortran 25.5 (A30) stops reading the converted input, as NV-5 does: a defect older than P3, [#86](https://github.com/szaghi/adam/issues/86) |

Related cases: the same modes across 2:1 seams are [SV-10](./amr-jeans), under runtime AMR [SV-11](./runtime-amr); the
isolated analogue of the collapse is [SV-8](./polytrope); the `[reference]` layer is
[scaling covariance](../infrastructure/scaling-covariance) (NV-5).

## Findings

- **The coupling is second order,** as P0 designed it (the 7-point Laplacian and the face form): the WENO-5 fluxes
  cannot lift it. The error after one period is $1.3\cdot10^{-4}$ of the amplitude at $N = 128$.
- **Only the time integration drifts.** In the collapse $E$ and $W$ each change by 17 % of the total, yet $E + W$
  drifts by $3\cdot10^{-7}$ at CFL 0.4, and the drift falls as $\Delta t^3$, SSP-RK3's order: the semi-discrete
  scheme conserves $E + W$ exactly, as the design predicted. Momentum stays at round-off.
- **The symmetry holds to round-off, not bitwise:** the summation orders of WENO and of the FFT are not mirror
  invariant. The defect relative to the mode is $5\cdot10^{-11}$ at the end; it is largest ($7\cdot10^{-10}$) early,
  when the amplitude is $10^{-6}$ and round-off is a larger fraction of it.
- **Cost.** CPU (WSL, indicative): a Poisson solve of the $32^2$ grid takes 0.2 ms; with self-gravity the $32^2$ run
  takes 31 s against 29 s without (one solve per stage and one per history row). FNL (one A30 per rank): 1.3 ms per
  solve at $32^2$, 4.4 ms at $128^2$, the density and the potential crossing the host each stage. A first version
  copied the whole device potential every stage, sized to the block capacity (23 270 blocks, 16 in use): 0.73 s per
  solve, whatever the grid. Only the active blocks cross now; the results are bitwise those of the first version.
- **CPU and FNL agree** to the printed digits on every SV-3 and SV-4 measure; the symmetry defects differ (round-off
  of different summation orders), both within the bound by a factor 10.
- **SV-5 on nvfortran 25.5** fails to read the converted input, as NV-5 does ([#86](https://github.com/szaghi/adam/issues/86)).
- The regression case `selfgravity-jeans` freezes `jeansc2d` at $32^2$ (4×4 blocks of $8^2$) for 20 steps, not the
  verification's $64^2$ collapse.
