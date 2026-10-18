# Radiative cooling

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| CO (co1 to co7) | Euler | M11 ([#93](https://github.com/szaghi/adam/issues/93)) | CPU, FNL | yes | `radiative-shock` |

## Case and scope

`[cooling]` adds optically thin cooling, $\partial_t(\rho e) = -n^2\Lambda(T)$, $n = \rho/m$, $T = p/(\rho R)$, to the
hydrodynamic models ([numerics](../../numerics#radiative-cooling)). It is applied once per step by Townsend's exact
integration of the isochoric cooling, so on a uniform state, where the flux divergence vanishes, every cell must follow
the cooling ODE to round-off whatever the step, and stop exactly at the floor $T_1$. On a flow the cooling meets the
dynamics only through the splitting; the radiative shock checks that coupling against a structure known in closed form.

What the case catches:
- **a cooling that is not exact:** a wrong segment, exponent or inverse, a lost digit near an exponent of 1, a step
  that depends on $\Delta t$;
- **a floor that leaks or overshoots;**
- **energy that is not accounted for:** $E$ plus the radiated energy must stay constant;
- **a dual-energy entropy left on the old adiabat:** the next sync would lift the cooled cells back and undo the
  cooling;
- **a backend that disagrees, a unit that is not converted, or a silent input;**
- **a cooling-dynamics coupling that is wrong:** the shock speed, the layer and the dense layer of the radiative shock.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Model | `euler` (and `euler-de`, co4) with `[cooling]`: $\Lambda$ the piecewise power law through the nodes, the first node the floor | [`adam_flume_cooling_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_cooling_object.F90) |
| Cooling | once per step after the dual-energy sync, Townsend's exact integration (D-M11-1), `dt_factor` limiter optional | [`adam_flume_cooling_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_cooling_library.F90), `apply_cooling` (CPU), `apply_cooling_dev` (FNL) |
| Space, time | as [Sod](./sod) (V1): WENO5 characteristic flux splitting, SSP-RK(3,3), CFL 0.5 | `adam_flume_cpu_euler_kernels.F90`, FNL `adam_flume_fnl_euler_kernels.F90` |

## Initial and boundary conditions

| Leg | Setup |
|---|---|
| co1 | the uniform box: the `sod-x` grid, 50 cells, $\rho = 1$ at rest, $\gamma = 5/3$, $R = 1$; the table $T = 0.1, 0.3, 1, 3$, $\Lambda = 0.2, 1, 0.5, 1.5$ (exponents 1.46, −0.58, 1, the last one continued above 3), $m = 1$. From $T = 5$ to $t = 6.5$ at the CFL step; from $T = 1$ to $t = 2$ with `dt_factor = 0.01` |
| co2 | the box from $T = 1$ to $t = 3$: the floor reached at $t = 2.2635$ |
| co3 | `[cooling]` with `mhd-ideal`; non-increasing temperatures; `lambdas` of another length; 65 nodes; no `particle_mass`; a negative `dt_factor` |
| co4 | the box moving at $u = 10$ ($e = 0.03E$, below the dual-energy switch) with `euler-de`, and with `euler`, to $t = 1.5$ |
| co5 | the co1 limited box to $t = 1$ written in the units $L_0 = 2$, $u_0 = 2$, $\rho_0 = 4$ through `[reference]` |
| co6 | the co1 CFL and co4 `euler-de` runs of CPU and FNL |
| co7 | the radiative shock: two streams $\rho = 1$, $u = \pm 1$ at $T = 0.1$ (Mach 2.45) collide at $x = 0.5$, 400 cells, $\Lambda = 0.18\,(T/0.1)^{1.5}$ (the floor at 0.1), to $t = 3$ |

**The radiative shock.** The shocked gas cools back to the floor and piles into a dense layer at rest. Taken whole,
shock plus cooling layer is an isothermal jump at $T_1$: in the frame of the structure, mass and momentum alone give
its speed and the dense layer,

$$V_s\,(u_0 + V_s) = R T_1, \qquad \rho_d = \rho_0\,\frac{u_0 + V_s}{V_s},$$

$V_s = 0.0916$, $\rho_d = 11.92$, with no energy equation. The cooling layer is steady in that frame: $\rho w = m$,
$p + \rho w^2 = P$ and $m\,dH/d\xi = -n^2\Lambda(T)$, $H = \gamma RT/(\gamma - 1) + w^2/2$, from the Rankine–Hugoniot
state behind the shock ($\rho = 2.82$, $T = 0.308$) to $T_1$, where $w$ reaches $V_s$: the end state of the ODE is the
dense layer of the jump, a consistency check of the two. The layer is 0.049 long, 20 cells. A planar radiative shock is
overstable for small cooling exponents (Chevalier & Imamura 1982): at $\alpha = 0.5$ the shock still oscillates at
$t = 1.5$ and its mean speed is 17 % below $V_s$; at $\alpha = 1.5$ the structure settles after $t \approx 1$.

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/cooling/make_co.py` | writes the inputs (`--case box|shock`, the state, the table, `dt_factor`) |
| `src/tests/flume/verification/cooling/check.sh` | runs the legs (`--legs co1,...,co7`, default all but co6; `FLUME_EXE` selects the backend; co6 compares the work directories of both) |
| `src/tests/flume/verification/cooling/co_oracle.py` | the isochoric ODE by Gauss–Legendre quadrature in $\ln T$ (independent of the closed form), the floor, the entropy, the CPU–FNL rule, the radiative shock's structure |
| `src/tests/flume/regression/radiative-shock/input.ini` | the frozen co7 case, 100 cells to $t = 1$ |
| `make_doc_runs.sh <runs> cooling`, `make_doc_figures.py --only radiative-cooling` | the figure (3 runs, 2 minutes) |

```bash
cd src/tests/flume/verification/cooling && ./check.sh                                  # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh                           # FNL
./check.sh --legs co6                                                                 # CPU against FNL
```

## Results

![Radiative cooling: the uniform box against the ODE and the radiative shock](/flume/radiative-cooling.png)

*Top: the uniform box from $T = 5$ follows the isochoric ODE through the extension above the table and its three
segments, to round-off at the CFL step and with the limiter. Bottom: the radiative shock at $t = 3$, density and
temperature against the steady cooling layer of the ODE aligned at the shock; the shock moves at $V_s$ once the
structure has settled.*

| Leg | Asserted | Bound | CPU | FNL |
|---|---|---|---|---|
| co1 | $T(t)$ against the ODE, from $T = 5$, CFL step | $10^{-12}$ | $5.8\cdot10^{-15}$ | $3.7\cdot10^{-14}$ |
| co1 | $T(t)$ against the ODE, `dt_factor = 0.01`; every step within $0.01\,t_\mathrm{cool}$, the limit reached | $10^{-12}$; 1 | $1.3\cdot10^{-14}$; reached on 47 % of the steps | $9.4\cdot10^{-15}$; the same |
| co1 | $E$ plus the radiated energy | $10^{-13}$ | $7\cdot10^{-16}$ | $1.9\cdot10^{-15}$ |
| co2 | $T = T_1$ from the step crossing the ODE's floor time on, $E$ constant | $10^{-14}$, bitwise | exact, bitwise | $2\cdot10^{-16}$, bitwise |
| co3 | the six invalid inputs | refused | refused | refused |
| co4 | `euler-de` moving box: $T(t)$; $\rho s$ on the cooled adiabat; the history of `euler` | $10^{-12}$; $10^{-13}$; bitwise | $1.7\cdot10^{-13}$; $4.9\cdot10^{-15}$; bitwise | the same |
| co5 | the box through `[reference]` against the code-units box | bitwise | bitwise | bitwise |
| co6 | FNL against CPU, every row: time, $T$, radiated | $10^{-12}$ | (the same comparison) | $5\cdot10^{-15}$, $2.4\cdot10^{-14}$, $1.0\cdot10^{-14}$ |
| co7 | the shock speed against $V_s$ (fitted from $t = 1.5$) | 1 % | $5.9\cdot10^{-3}$ | $5.9\cdot10^{-3}$ |
| co7 | the cooling layer's density against the ODE (L1); its $\rho = 6$ point | 5 %; one cell | $2.6\cdot10^{-2}$; 0.22 cells | $2.6\cdot10^{-2}$; 0.22 cells |
| co7 | the outer half of the dense layer: $\rho$ against $\rho_d$, $T$ against $T_1$ | 2 % | $+2.4\cdot10^{-3}$, $-4.0\cdot10^{-3}$ | $+2.4\cdot10^{-3}$, $-4.0\cdot10^{-3}$ |

Related cases: [dual energy](./dual-energy) holds the entropy co4 keeps consistent, and the cold Noh shock whose
collision co7 cools; [Sod](./sod) is the grid of the box.

## Findings

- **The cooling is exact at any step.** Townsend's integration turns the stiff ODE into a closed form: the box follows
  the quadrature of the ODE to $10^{-14}$ at the CFL step, with no cooling limit on $\Delta t$, and the floor is reached
  exactly. The `dt_factor` limiter is there for the coupling with the dynamics, not for the cooling.
- **The entropy of the dual energy must follow the radiation.** Cooling lowers the entropy for real; the cooled cells
  take the entropy of the internal energy the stored $E$ holds. From the unrounded value, one ulp of rounding in $E$ left
  every cold cooled cell below its adiabat and the next sync lifted it.
- **The radiative shock checks the splitting.** Its speed and its dense layer follow from the isothermal jump alone,
  its layer from a steady ODE whose end state is that jump; the run matches the three. The case needs a stable cooling
  exponent: at $\alpha = 0.5$ the overstable shock oscillates for longer than the run.
- **The centre of the dense layer keeps the start-up transient.** The gas of the first collision sits about 6 % below
  the floor temperature, at the dense layer's pressure; the floor stops the cooling, not the dissipation of the
  scheme, as in Noh's wall heating with the opposite sign.
