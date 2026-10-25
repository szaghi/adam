# Jeans modes on a static AMR grid

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| SV-10 | Self-gravity | M8 ([#81](https://github.com/szaghi/adam/issues/81)), P5b, P5c | CPU, FNL | yes | `selfgravity-amr-stage`; `amr-periodic-reflux-stage` (the per-stage reflux, without gravity) |

## Case and scope

On an AMR grid three things change in self-gravity:
- the potential comes from the multigrid solver of [SV-9](./multigrid-amr);
- the face gravity at the 2:1 faces is that solver's first-order mean-flux gradient;
- the energy exchange needs the seam gravity work carried by the reflux ([models](../../models#self-gravity)).

The Jeans modes of [SV-3 and SV-4](./jeans) check the coupling, now with a 2:1 seam the wave crosses on every side of
the refined box.

What the case catches:
- a multigrid coupling that differs from the FFT one on a grid where both apply;
- a loss of order across the seams (the stable wave);
- a broken growth rate or point symmetry of the unstable mode;
- a missing or wrong seam gravity work: the semi-discrete $E + W$ is no longer exact, and the drift stops converging in
  $\Delta t$;
- with the per-stage reflux (`reflux_cadence = stage`, P5c), a drift that does not recover the Runge–Kutta order;
- a solver selection that accepts what it cannot solve (isolated boundaries on a periodic grid with the multigrid
  solver, `poisson_solver = fft` on AMR).

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | compressible Euler, $\gamma = 5/3$, periodic self-gravity | [models](../../models#self-gravity) |
| Solver selection | `poisson_solver = auto`: the FFT solver on a uniform grid, the multigrid on an AMR grid; `fft` on AMR is fatal, and so is isolated on a periodic grid (the multigrid's own check; isolated boundaries on AMR are [SV-13](./multigrid-amr) and [SV-14](./polytrope), [#87](https://github.com/szaghi/adam/issues/87)) | `set_potential`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90); [`adam_flume_gravity_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_gravity_object.F90) |
| Poisson solver | multigrid (FCG + V-cycle) on the mimetic composite Laplacian, on the host (CPU) or the device (FNL, #84 E6a); its coarsest level on the backend's FFT solver (FFTW, or cuFFT on FNL) | [`adam_poisson_mg_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_poisson_mg_object.F90), [`adam_fnl_poisson_mg_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/fnl/adam_fnl_poisson_mg_object.F90), [`adam_fnl_poisson_fft_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/fnl/adam_fnl_poisson_fft_object.F90) |
| Seam gravity work | added per stage to the flux register's energy entry: the coarse side takes its own seam source back, the fine side adds the exact two-point work less its own source | `add_seam_gravity_work`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90); `accumulate_seam_fluxes`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90), [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |
| Reflux | `reflux = .true.`; `reflux_cadence = step` (default, end of step) or `stage` (P5c: the coarse seam cells corrected at every Runge–Kutta stage) | `reflux_stage`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [`adam_flume_numerics_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_numerics_object.F90); [numerics](../../numerics#adaptive-mesh-refinement) |
| Space | WENO5 (`weno-u-5`) flux splitting on characteristic variables (`weno`) | [numerics](../../numerics#scheme-space-weno-flux-splitting) |
| Time | SSP Runge–Kutta, 3 stages, CFL 0.4 (collapse: 0.4, 0.2, 0.1) | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |

## Initial and boundary conditions

`make_jeans.py --box` refines the central box (a quarter of the period from every face; in 3-D also along $z$) one
level more. The symmetry of the unstable mode survives, because the box is symmetric about the point the oracle
reflects through.

| Item | Value |
|---|---|
| Physics | as SV-3 / SV-4: $\rho_0 = 1$, $c = 1$, $\gamma = 5/3$, $\mathbf k = 2\pi(1, 1)$, every active face periodic; the stable wave $\lambda = \lambda_J/2$ ($G = \pi/2$), the unstable mode $\lambda = 1.5\lambda_J$ ($G = 9\pi/2$), [`sine-wave`](../../initial-conditions#sine-wave-a-small-sine-wave-on-a-uniform-state-verification) with `wave_mode = jeans` |
| Grid | base $4n$ (2-D, 4×4 blocks) or $2n$ (3-D, 2×2×2 blocks) for $n$ cells per block side, the box at twice the resolution; 4 ghosts |
| Ranks, backends | 2; CPU and FNL (the multigrid on the device since #84 E6a, its coarsest level on cuFFT) |

| Leg (`jeans.sh`) | Case | Base grid | Amplitude, end |
|---|---|---|---|
| `sv10`: multigrid against FFT | `jeans2d`, uniform, `poisson_solver = fft` and `multigrid` | $32^2$ | $10^{-6}$; one period |
| `sv10`: unstable mode | `jeansu2d --box`, a checkpoint every 5 steps | $64^2$ | $10^{-6}$; $t = 0.5$ |
| `sv10`: collapse | `jeansc2d --box`, CFL 0.4, 0.2, 0.1 | $64^2$ | $10^{-2}$; $t = 0.4$ |
| `sv10`: collapse, per-stage reflux | the same with `--reflux-cadence stage` | $64^2$ | $10^{-2}$; $t = 0.4$ |
| `sv10`: refusals | `jeans2d --solver multigrid --poisson-bc isolated`; `jeans2d --box --solver fft` | $32^2$ | — |
| `sv10-wave`: stable wave, 2-D | `jeans2d --box` | $32^2$, $64^2$, $128^2$ | $10^{-6}$; one period |
| `sv10-wave`: stable wave, 3-D | `jeans3d --box` | $16^3$, $32^3$ | $10^{-6}$; one period |
| Regression `selfgravity-amr-stage` | `jeansc2d --box --reflux-cadence stage` | $32^2$ | $10^{-2}$; 20 steps |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/selfgravity/make_jeans.py` | the Jeans cases; `--box`, `--solver`, `--poisson-bc`, `--reflux-cadence` |
| `src/tests/flume/verification/selfgravity/jeans.sh` | legs `sv10` and `sv10-wave` (`FLUME_EXE` selects the backend) |
| `src/tests/flume/verification/selfgravity/jeans_oracle.py` | `compare`, `wave`, `growth`, `symmetry`, `conservation` (with `--tol-mass` for the per-stage reflux) |
| `src/tests/flume/regression/selfgravity-amr-stage/input.ini` | the frozen per-stage collapse on the box-refined grid, $32^2$ base, 20 steps |
| `src/tests/flume/regression/amr-periodic-reflux-stage/input.ini` | the per-stage reflux on a triply periodic octree (no gravity) |
| `src/tests/flume/verification/make_doc_runs.sh`, group `amr` | `sv10-2d-{8,16,32}`, `sv10-uniform-{8,16,32}`, `sv10-collapse-{0.4,0.2,0.1}`, `sv10-stage-{0.4,0.2,0.1}` |
| `src/tests/flume/verification/make_doc_figures.py`, `fig_selfgravity_amr_jeans` | draws the figure (`--only selfgravity-amr-jeans`) |

```bash
cd src/tests/flume/verification/selfgravity && ./jeans.sh --legs "sv10 sv10-wave"        # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl_fft ./jeans.sh --legs "sv10 sv10-wave"  # FNL
```

## Results

![SV-10: the wave error on AMR against N, the drift of E + W in the collapse, the collapse density on the composite grid](/flume/selfgravity-amr-jeans.png)

*Left: the L1 error of the stable wave falls at second order across the seams. Middle: in the collapse $E$ and $W$
trade 17 % of $\lvert E + W\rvert$ while $E + W$ drifts by about $10^{-7}$; with the end-of-step reflux (solid) the
drift falls 2.8 and then 2.2 times per halving of the CFL number, with the per-stage reflux (dashed) 8 times, the
SSP-RK3 order. Right: the collapse density at $t = 0.4$ on the composite grid; the refined blocks are the central
box.*

| Asserted | Bound | CPU | FNL (WSL GPU) |
|---|---|---|---|
| on a uniform $32^2$ grid, `poisson_solver = multigrid` against `fft`, $\max\lvert\rho_{mg} - \rho_{fft}\rvert / (\rho_0A)$ after one period | $10^{-6}$ | $3.8\cdot10^{-9}$ | $2.7\cdot10^{-9}$ |
| stable wave, order of the L1 error, finest pair, 2-D | $\ge 1.8$ | 1.98 (1.97, 1.98) | 1.98 |
| the same, 3-D ($16^3 \to 32^3$) | $\ge 1.8$ | 2.00 | 2.00 |
| unstable mode (base $64^2$), growth rate, relative error | $10^{-2}$ | $5.8\cdot10^{-4}$ | $5.8\cdot10^{-4}$ |
| its mirror symmetry at $t = 0.5$ | $10^{-9}$ | $5.7\cdot10^{-11}$ | $4.7\cdot10^{-11}$ |
| collapse (amplitude $10^{-2}$, base $64^2$, $t = 0.4$), max $\lvert$momentum$\rvert$ | $10^{-14}$ | $3.5\cdot10^{-18}$ | $6.9\cdot10^{-18}$ |
| its $E + W$ drift at CFL 0.4, 0.2, 0.1, order of the finest pair | $\ge 0.8$ | $8.8\cdot10^{-8}$, $3.1\cdot10^{-8}$, $1.4\cdot10^{-8}$ (1.50, 1.16) | the same to 3 digits |
| the same with `reflux_cadence = stage` (P5c), order of the finest pair; mass drift | $\ge 2.5$; $10^{-14}$ | $3.5\cdot10^{-8}$, $4.3\cdot10^{-9}$, $5.4\cdot10^{-10}$ (2.99, 3.00); $2.2\cdot10^{-16}$ | the same to 3 digits; $3.3\cdot10^{-16}$ |
| isolated boundaries on a periodic grid with the multigrid solver; `poisson_solver = fft` on AMR | refused, exit 1 | refused | refused |

Related cases: the uniform-grid modes are [SV-3 to SV-5](./jeans); the solver is [SV-9](./multigrid-amr); the same
collapse under runtime AMR is [SV-11](./runtime-amr); conservation across seams without gravity is
[conservation across AMR seams](../euler/amr-conservation).

## Findings

- **The seam gravity work makes the energy exchange exact in space.** A run with the work switched off (a temporary
  ablation build, not a leg) drifts by $1.03\cdot10^{-5}$, $9.8\cdot10^{-6}$, $9.6\cdot10^{-6}$ at the three CFL
  numbers, orders 0.07 and 0.02. That is a spatial leak at the seams, which no $\Delta t$ removes. With the work the
  drift is 12 times smaller at CFL 0.4 and converges in $\Delta t$.
- **With the end-of-step reflux the drift is first order in $\Delta t$, not third.** The reflux corrects the coarse
  cells at the end of the step, and the stages in between advance them with their own seam fluxes. On a uniform grid
  (SV-4) the order is that of SSP-RK3. The per-stage reflux (`reflux_cadence = stage`, P5c) restores it: orders 2.99
  and 3.00, and the drift at CFL 0.1 falls from $1.4\cdot10^{-8}$ to $5.4\cdot10^{-10}$.
- **The multigrid's stopping rule needed a round-off floor.** With a $10^{-6}$ mode on a unit density the source is
  known to $10^{-16}/10^{-6} = 10^{-10}$ of itself. The $10^{-10}$ relative tolerance was unreachable, and every solve
  ran to the 200-iteration cap. The solver now stops at $\max(10^{-10},\ 100\,\epsilon\max\lvert4\pi G\rho\rvert/
  \max\lvert4\pi G(\rho - \bar\rho)\rvert)$ and warns if it hits the cap.
- **The AMR wave is more accurate than the uniform one at the same base N.** At base $32^2$ the L1 error is
  $1.9\cdot10^{-3}$ against $2.1\cdot10^{-3}$, because the box adds resolution. The first-order seam gravity shows in
  no norm the legs measure.
- Isolated boundaries on AMR came with [#87](https://github.com/szaghi/adam/issues/87): the multigrid's symmetric
  James step ([SV-13](./multigrid-amr)), the polytrope on AMR ([SV-14](./polytrope)).
