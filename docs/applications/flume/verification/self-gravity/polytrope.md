# The n = 1 polytrope: equilibrium and collapse

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| SV-7, SV-8 | Self-gravity | M8 ([#81](https://github.com/szaghi/adam/issues/81)), P4 | CPU, FNL | yes | `selfgravity-star` |

## Case and scope

An n = 1 polytrope star under the isolated self-gravity of [SV-6](./isolated-solver), centred in a closed box. The
n = 1 polytrope is the Lane–Emden equilibrium $p = K\rho^2$ with a closed-form density, so the discrete scheme holds
it up to its truncation error, and halving its pressure makes it collapse, bounce and ring.

What the case catches:
- **SV-7, the equilibrium:** a potential or a force inconsistent with the pressure gradient. The spurious motion of
  the discrete scheme around the continuum equilibrium (the drift of the central density, the Mach number inside the
  star) must shrink with the grid;
- **SV-8, the collapse:** a leak in the exact exchange of $E$ and $W$ through the isolated solver. Its symmetry (an
  even kernel) must give the same exactness as the periodic solver: $E + W$ drifts only at the Runge–Kutta order while
  $E$ and $W$ exchange, and the momentum stays at round-off.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | compressible Euler, $\gamma = 2$ (the adiabatic exponent of $p = K\rho^2$), isolated self-gravity, $G = 1$ | [models](../../models#self-gravity) |
| Poisson solver | FFT with isolated boundaries (`poisson_bc = isolated`): zero padding and the lattice Green's function | [`adam_poisson_fft_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_poisson_fft_object.F90), [`adam_fnl_poisson_fft_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/fnl/adam_fnl_poisson_fft_object.F90) |
| Coupling | one solve per Runge–Kutta stage, the potential path | `solve_potential`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90); `add_potential_sources`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Space | WENO5 (`weno-u-5`) flux splitting on characteristic variables (`weno`) | [numerics](../../numerics#scheme-space-weno-flux-splitting) |
| Time | SSP Runge–Kutta, 3 stages, CFL 0.4 (SV-8: 0.4, 0.2, 0.1) | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Positivity limiter | SV-7: `cell` (its gravity energy correction, M7-G1 option B); SV-8: `none` | `limit_positivity`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [numerics](../../numerics#positivity-limiter) |
| Boundaries | `wall-inviscid` on every face: a closed box, no mass or energy crosses it, so $E + W$ has no boundary term | [boundary conditions](../../boundary-conditions#wall-inviscid-slip-wall-by-mirroring-with-a-sign-vector); [`adam_flume_bc_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_bc_object.F90) |
| Initial condition | `polytrope`: the star and its atmosphere | `set_initial_conditions`, [`adam_flume_ic_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_ic_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Star | [`polytrope`](../../initial-conditions#polytrope-an-n-1-polytrope-star-self-gravity): n = 1, radius $a = 1/2$, $\rho_c = 1$, $G = 1$, $p = K\rho^2$ with $K = 2Ga^2/\pi$ (the Lane–Emden equilibrium), $\gamma = 2$; mass $4\rho_ca^3/\pi$, mean density $3\rho_c/\pi^2$, free-fall time of the mean density $t_{ff} = \sqrt{3\pi/(32G\bar\rho)} = 0.98$ |
| Atmosphere | isothermal at the central temperature ($p_1 = K\rho_c\rho_1$) and hydrostatic in the star's field, $\rho = \rho_1 e^{GM/T\,(1/r - 1/a)}$ outside; inside, $\rho_1$ and $p_1$ added to the star |
| Box | $[-1, 1]^3$, an octree, uniform; 8 cells per block side, 4 ghosts; `wall-inviscid` on every face |
| Ranks | 2 |
| SV-7, equilibrium (`star`) | `star_pressure_factor` 1, atmosphere $\rho_1 = 10^{-4}\rho_c$, positivity limiter `cell`; $32^3$ and $64^3$ (levels 2 and 3) to $t = 2$, CFL 0.4, a checkpoint every 10 steps |
| SV-8, collapse (`starc`) | `star_pressure_factor` 0.5, atmosphere $\rho_1 = 10^{-2}\rho_c$, no limiter; $32^3$ to $t = 1.5$ at CFL 0.4, 0.2, 0.1 |
| Regression `selfgravity-star` | `star --levels 1`: $16^3$ (2×2×2 blocks of $8^3$), atmosphere $10^{-2}\rho_c$, no limiter, 20 steps |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/selfgravity/make_polytrope.py` | writes `star` and `starc` (`--cells`, `--levels`, `--cfl`, `--time-max`, `--it-save`, `--atmosphere`, `--limiter`) |
| `src/tests/flume/verification/selfgravity/star.sh` | runs the legs `sv7`, `sv8` (`--np N`, `--legs LIST`, `--keep`; `FLUME_EXE` selects the backend) |
| `src/tests/flume/verification/selfgravity/star_oracle.py` | `convergence` (SV-7): the drifts on two grids and their reduction |
| `src/tests/flume/verification/selfgravity/jeans_oracle.py` | `conservation` (SV-8): $E + W$ against $\Delta t$, the momentum |
| `src/tests/flume/regression/selfgravity-star/input.ini` | the frozen equilibrium at $16^3$, 20 steps |
| `src/tests/flume/verification/make_doc_runs.sh`, group `isolated` | `sv7-l{2,3}`, `sv8-collapse-{0.4,0.2,0.1}` (a checkpoint every 10 steps); needs `exe/adam_flume_cpu_fft` and `exe/test_poisson_fft` |
| `src/tests/flume/verification/make_doc_figures.py`, `fig_selfgravity_star` | draws the figure (`--only selfgravity-star`) |

```bash
cd src/tests/flume/verification/selfgravity && ./star.sh                   # CPU (exe/adam_flume_cpu_fft)
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl_fft ./star.sh               # FNL
```

## Results

![SV-7 and SV-8: the equilibrium star's core and spurious motion at two resolutions; the collapse's energy drift and density](/flume/selfgravity-star.png)

*Top: at $32^3$ and $64^3$ the core first expands (the central density falls, by 10 % and 4.5 %) and then recompresses
past its initial value within two time units; on the finer grid the swing is smaller and the interior Mach number
about half. Bottom left: in the collapse $E$ and $W$ exchange 74 % of $\lvert E + W\rvert$ (their curves coincide)
while $E + W$ drifts at the Runge–Kutta order, a factor 8 per halving of the CFL number. Bottom right: the collapsing
star at its densest checkpoint, $t = 1.03$, its central density four times the initial one.*

The central density is the mean over the cells within one cell of the centre; the Mach numbers are taken inside
$0.8a$ and over the whole star (the cells denser than a tenth of the central density).

| Leg | Asserted | Bound | CPU | FNL |
|---|---|---|---|---|
| SV-7 | the largest central-density drift to $t = 2$: $9.91\cdot10^{-2}$ ($32^3$), $5.71\cdot10^{-2}$ ($64^3$); reduction | $\ge 1.5$ | 1.74 | 1.73 |
| SV-7 | the largest Mach number inside $0.8a$: 0.154, 0.075; reduction | $\ge 1.5$ | 2.05 | 2.06 |
| SV-7 | the largest Mach number in the whole star (set at its surface): 0.43, 0.29; reduction | reported | 1.50 | 1.49 |
| SV-8 | drift of $E + W$ at CFL 0.4, 0.2, 0.1: $5.97\cdot10^{-7}$, $8.29\cdot10^{-8}$, $9.71\cdot10^{-9}$; order in $\Delta t$ | $\ge 2.5$ | 2.85, 3.09 | 2.85, 3.09 (the same drifts) |
| SV-8 | the largest $\lvert\int\rho\mathbf u\rvert$ | $10^{-14}$ | $7.3\cdot10^{-18}$ | $6.9\cdot10^{-18}$ |

Related cases: the periodic analogue of the collapse is [SV-4](./jeans); the solver is [SV-6](./isolated-solver);
the limiter's gravity energy correction comes from [the force](../gravity/force) (M7-G1).

## Findings

- **Only the time integration drifts.** In the collapse $E$ and $W$ each change by 74 % of $\lvert E + W\rvert$, while
  $E + W$ drifts by $6\cdot10^{-7}$ at CFL 0.4 and falls as $\Delta t^3$: the isolated solver conserves as exactly as
  the periodic one (SV-4), through its symmetry.
- **The positivity limiter breaks the exchange where it acts.** Its gravity energy correction (M7-G1, option B) is not
  part of the flux form. A star in a $10^{-4}$ atmosphere needs it at its surface (without it the gas there reaches
  $\lvert\mathbf u\rvert \approx 0.4$ in ten steps and its pressure goes negative), and with it $E + W$ drifts by
  $2.8\cdot10^{-3}$ by $t = 4$ at $32^3$, against $3\cdot10^{-9}$ without it in a $10^{-2}$ atmosphere. SV-8 therefore
  runs in the heavier atmosphere, without the limiter.
- **A heavy atmosphere is not an equilibrium test.** The $10^{-2}$ atmosphere holds 22 % of the star's mass, and the
  star then pulsates with 13 % of its central density at $32^3$; SV-7 uses the light one.
- **The equilibrium converges at about first order, set by the star's surface.** The n = 1 surface has a kink in the
  density ($\rho_s \propto a - r$, $p \propto (a - r)^2$), the light atmosphere needs the limiter there (first order
  where it acts), and the surface is where the largest Mach number sits ($r = 0.454$, $a = 0.5$): it shrinks only
  1.5-fold from $32^3$ to $64^3$. The pulsation it drives sets the core's drift, which shrinks 1.7-fold, and the
  interior motion 2.1-fold. A star of radius 8 cells pulsates by 10 % of its central density, of 16 cells by 6 %.
- **The bound asserts convergence, not an order.** P0 asked for 100 free-fall times; the runs reach 2 ($64^3$ took
  37 min on the WSL CPU, 20 min on the A30).
- The regression case `selfgravity-star` freezes the equilibrium at $16^3$ in the $10^{-2}$ atmosphere without the
  limiter, not SV-7's light-atmosphere setup.
