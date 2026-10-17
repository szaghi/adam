# Gravity from a cell potential

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| SV-1 | Self-gravity | M8 ([#81](https://github.com/szaghi/adam/issues/81)), P1 | CPU, FNL | yes | — |

## Case and scope

Self-gravity applies gravity through a cell potential $\Phi$. That path must be right before a Poisson solver feeds
it: `add_potential_sources` (face gravity $-\partial\Phi/\partial x$, cell gravity the mean of the cell's two faces,
momentum $\rho$ times the cell gravity, energy from the face mass fluxes times the face gravity), the face and cell
gravity in the positivity limiter, and the hydrostatic ghost from the potential drop
$-(\Phi_{\text{ghost}} - \Phi_1)$.

`[gravity] potential = linear` writes the constant $\mathbf g$ of M7 as $\Phi = -\mathbf g\cdot\mathbf x$ and applies it
through that path alone (the constant-g kernels get zero), so the M7 oracles apply unchanged.

What the case catches:
- a wrong sign, spacing or face average in the potential sources: the free fall and the hydrostatic ghosts miss their
  exact values;
- an energy source that is not the telescoping face form: the closed box loses $E + \rho\Phi$;
- a change of the constant-g path: that path must stay bitwise as before, which the regression suite holds.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | compressible Euler, ideal gas, $\gamma = 1.4$; `iso2d-mhd`, `box2d-mhd`: ideal MHD with EGLM cleaning (`glm_ch` 3, `glm_alpha` 0.18, `glm_damping_length` 1) | [models](../../models#external-gravity) |
| Gravity | `[gravity] potential = linear`: $\Phi = -\mathbf g\cdot\mathbf x$ at the cell centres, ghosts included | `set_potential`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90); `load_from_file`, [`adam_flume_gravity_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_gravity_object.F90) |
| Potential sources | face gravity, cell gravity, momentum, energy from the face mass fluxes (second order, no fourth-order correction) | `add_potential_sources`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `add_potential_sources_dev`, [`adam_flume_fnl_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_kernels.F90) |
| Space | WENO5 (`weno-u-5`) flux splitting on characteristic variables | [numerics](../../numerics#scheme-space-weno-flux-splitting) |
| Time | SSP Runge–Kutta, 3 stages (`runge-kutta-ssp-33`), CFL 0.4; the `wave2d` ladder scales it as GV-9, $0.4\,(8/n)^{2/3}$ for $n$ cells per block side | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Positivity limiter (`box2d-lim`, `box2d-mhd`) | `cell`, with the face and cell gravity from $\Phi$ | `limit_positivity`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [numerics](../../numerics#positivity-limiter) |
| Hydrostatic ghosts | the potential drop $\mathbf g\cdot\mathbf s = -(\Phi_{\text{ghost}} - \Phi_1)$ | `set_boundary_conditions`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `hydrostatic_ghost`, [`adam_flume_bc_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_bc_object.F90) |
| FNL backend | the same path on the device | [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

## Initial and boundary conditions

The cases are those of M7 (`gravity/make_atmosphere.py`), written with `--potential linear`. Every case: the unit
square, a quadtree of 4×4 blocks (two uniform levels), 4 ghost cells, 2 ranks.

| Case (M7 origin) | State, gravity | Boundaries | Cells per block side, end |
|---|---|---|---|
| `iso2d` (G0) | isothermal atmosphere, $\mathbf g = (0, -1)$ ([IC](../../initial-conditions#isothermal-atmosphere-polytropic-atmosphere-hydrostatic-equilibria-under-gravity)); ghosts poisoned | $x$ periodic, $y$ `hydrostatic` | 8; one step at CFL $10^{-30}$ |
| `iso2d-mhd` (G0) | `iso2d` with EGLM and $\mathbf B = (0.5, 0.7, 0.2)$ | as `iso2d` | 8; one step |
| `freefall2d` (G1) | uniform $\rho = 1$, $p = 1$ at rest ([IC](../../initial-conditions#uniform-seeded-uniform-state)), $\mathbf g = (0, -1)$ | periodic on every face | 8; 50 steps |
| `box2d` (G1) | as `freefall2d` | $x$ periodic, $y$ `wall-inviscid` | 8; 200 steps |
| `box2d-mhd` (G1) | `box2d` with EGLM, $\mathbf B = (0.5, 0, 0)$, the limiter on | as `box2d` | 8; 200 steps |
| `box2d-lim` (G1) | `box2d` with $p = 0.05$ and the positivity limiter (the gas falls away from the top wall toward vacuum) | as `box2d` | 8; 300 steps |
| `wave2d` (G2a, GV-9) | $\rho = 1 + 0.2\sin(2\pi(x + y))$ ([`sine-wave`](../../initial-conditions#sine-wave-a-small-sine-wave-on-a-uniform-state-verification), entropy mode), $\mathbf u_0 = (1, 0.5)$, $p = 1$, $\mathbf g = (-0.6, -0.8)$ | periodic on every face | 8, 16, 32 ($32^2$ to $128^2$); $t = 0.2$ |

The figure's closed box is `box2d` at 16 cells per block side ($N = 64$) to $t = 0.2$, checkpointed every 2 steps,
once through each path.

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/gravity/make_atmosphere.py` | writes the M7 cases; `--potential linear` selects the potential path |
| `src/tests/flume/verification/selfgravity/check.sh` | runs SV-1 (`--np N`, `--legs LIST`: `iso2d`, `iso2d-mhd`, `freefall2d`, `box2d`, `box2d-mhd`, `box2d-lim`, `wave2d` (the ladder), default all; `--keep`; `FLUME_EXE` selects the backend) |
| `src/tests/flume/verification/gravity/atmosphere_oracle.py` | the M7 oracles: `profile`, `freefall`, `energy`, `wave` |
| `src/tests/flume/verification/make_doc_runs.sh`, group `selfgravity` | `sv1-box2d-{uniform,linear}`, `sv1-wave-{uniform,linear}-{8,16,32}`: both paths |
| `src/tests/flume/verification/make_doc_figures.py`, `fig_selfgravity_potential` | draws the figure (`--only selfgravity-potential`) |

```bash
cd src/tests/flume/verification/selfgravity && ./check.sh                   # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh                   # FNL
cd .. && ./make_doc_runs.sh <dir> selfgravity                               # the figure's runs
python3 make_doc_figures.py --out ../../../../docs/public/flume --runs <dir> --only selfgravity-potential
```

## Results

![SV-1: energy budget of the closed box through both paths, and the falling-wave errors against N](/flume/selfgravity-potential.png)

*Left: the closed box at $N = 64$, each path's conserved energy at round-off while $E$ alone exchanges with $\rho\Phi$
by up to $6\cdot10^{-3}$. Right: the falling wave to $t = 0.2$, second order through the potential path (5.0 and 4.6
through the constant path).*

| Case | Asserted | Bound | CPU | FNL |
|---|---|---|---|---|
| `iso2d`, `iso2d-mhd` | hydrostatic ghosts against the profile, relative | $10^{-13}$ | $2.2\cdot10^{-16}$ | $2.2\cdot10^{-16}$ |
| `freefall2d` | $\rho u = \rho g t$, $E = E_0 + \rho g^2t^2/2$, relative | $10^{-12}$ | $5.6\cdot10^{-16}$ | $3.3\cdot10^{-16}$ |
| `box2d`, `box2d-mhd` | mass and $E + \rho\Phi$ in the closed box, relative | $10^{-12}$ | $1.3\cdot10^{-14}$, $7.5\cdot10^{-15}$ | $1.7\cdot10^{-14}$, $8.8\cdot10^{-15}$ |
| `box2d-lim` | the run completes (limiter acting); drift reported | — | $4.4\cdot10^{-7}$ | $4.4\cdot10^{-7}$ |
| `wave2d`, 3 resolutions | L1 order of $\rho$ and $p$ on the finest pair | $\ge 1.8$ | 2.12, 2.00 | 2.12, 2.00 |

Related cases: the constant-g path of the same cases is verified in [the force](../gravity/force),
[the order of the energy source](../gravity/energy-order) and [the stratified atmospheres](../gravity/atmospheres);
the Poisson solver that feeds this path is [SV-2](./fft-poisson).

## Findings

- In the closed box the potential path's energy budget closes on $E + \rho\Phi$ alone. The constant path closes on
  $W$, which adds the $O(\Delta x^2)$ wall term of its fourth-order source: without that term its $E + \rho\Phi$ drifts
  by $2\cdot10^{-6}$.
- With the limiter acting, the drift of the potential path ($4.4\cdot10^{-7}$) is of the size of the constant path's
  ($1.2\cdot10^{-6}$).
- The potential path is second order (2.12 in $\rho$, 2.00 in $p$), as its face form without a fourth-order correction
  must be; the constant path gives 5.0 and 4.6 (GV-9).
- The constant-g path is unchanged bitwise: the regression suite (35/35 on CPU and FNL at the time of P1) holds it.
