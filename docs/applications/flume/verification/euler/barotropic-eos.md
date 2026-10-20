# Barotropic equation of state

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| EO (eb1 to eb8) | Euler | M13 ([#100](https://github.com/szaghi/adam/issues/100)), P3 | CPU, FNL | yes | `baro-tube`, `baro-collide-hllc` |

## Case and scope

`[eos] type = barotropic` closes the Euler equations with a continuous piecewise polytrope, $P = K_k\rho^{\Gamma_k}$
([numerics](../../numerics#isothermal-and-barotropic-equations-of-state),
[input](../../input#eos-optional-equation-of-state)): the energy slot is inert, as for the
[isothermal gas](./isothermal-eos), and the sound speed $c^2 = dP/d\rho$ jumps where $\Gamma$ changes. On a convex
table, where $\Gamma$ only increases upward, the Riemann problem still has an exact solution built from mass and
momentum alone, and `eos_oracle.py` constructs it independently of the code; the default table of Mathew & Federrath
(2021) is not convex, so it is checked by self-convergence.

What the case catches:
- **a wrong table or segment lookup:** the plateau at a break density, which a rarefaction through a convex break
  leaves between its two fans, and the shocks that jump across segments;
- **a wrong Roe speed:** the HLL and HLLC speeds use the secant $\sqrt{[P]/[\rho]}$;
- **a wrong $c^2 = P'(\rho)$:** the acoustic wave's speed, hence its order;
- **a table that disagrees with the closed form:** the one-segment table $\Gamma = 1$ against the isothermal gas;
- **a scheme that does not converge** where no exact solution is available, across a non-convex break;
- **a backend that disagrees, or a table the device does not hold.**

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Model | `euler` with `[eos] type = barotropic`: `MODEL_EULER_BARO`, the energy slot inert (D-M13-1), the table registered in `adam_flume_eos_tables_library` | [`adam_flume_euler_baro_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_euler_baro_library.F90), [`adam_flume_eos_tables_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_eos_tables_library.F90) |
| Space | WENO5 characteristic flux splitting (eb1 to eb3, eb5, eb6); `weno-riemann` with `llf`, `hll`, `hllc` (eb7) | `adam_flume_cpu_euler_baro_kernels.F90`, `adam_flume_cpu_euler_baro_{llf,hll,hllc}_kernels.F90`; FNL `adam_flume_fnl_euler_baro_*` |
| Time | SSP-RK(3,3), CFL 0.5; SSP-RK(5,4) at CFL 0.1 for the order (eb3) | — |

## Initial and boundary conditions

The convex test table is `POLY` of `eos_oracle.py`: $c_0 = 1$, breaks $\rho = 2, 8$, exponents $1, 1.4, 5/3$.

| Leg | Setup |
|---|---|
| eb1 | `POLY` tube, $\rho = 10 \mid 0.5$, $v = 0 \mid 0.7$, at rest, 1600 cells, $y$ active, $t = 0.15$: $\rho_* = 2.602$; the left fan runs from $\rho = 10$ through the break 8, where it splits into two fans with a plateau at $\rho = 8$ |
| eb2 | `POLY` collision $(\rho, u) = (1, 3.4) \mid (1.5, -2.6)$, 400 cells, $t = 0.3$: $\rho_* = 8.056$, both shocks into the $\Gamma = 5/3$ segment, moving at $-0.466$ and $0.612$ |
| eb3 | the acoustic wave of amplitude $10^{-6}$ on $\rho_0 = 5$ ($\Gamma = 1.4$, $c^2 = P'(5) = 2.0198$), periodic $[0, 1]$, $N = 12, 24, 48$, one period |
| eb4 | decreasing `densities`, `exponents` of another count, an exponent below 1, 9 segments, a negative `reference_density`, the gas keys, `[reference]`, no `sound_speed` |
| eb5 | the one-segment table `exponents = 1` and the isothermal gas, both $c_0 = 1.3$, on the eo1 tube, $t = 0.15$ |
| eb6 | the default table (cgs, $c_0 = 0.2$ km/s), the periodic tube $\rho = 10^{-7} \mid 10^{-8}$ g cm$^{-3}$ at rest on 1 cm, $t = 3\cdot10^{-7}$ s: the left fan runs down through the non-convex break $3.84\cdot10^{-8}$; $N = 200, 400, 800$ and the reference $N = 1600$, `[weno] weights = si` (the data are far from order 1) |
| eb7 | eb2 with `weno-riemann` and each solver |
| eb8 | the eb1, eb2, eb3 ($N = 48$) and eb7 `hllc` runs of CPU and FNL |

**The frame of the collision.** eb2 is the collision $(1, 4) \mid (1.5, -2)$ seen from a frame moving at $0.6$. In the
original frame the left shock moves at $0.134$ relative to the grid, and a slowly moving shock sheds persistent
post-shock noise into the flow behind it (Arora & Roe 1997): under HLL and HLLC the layer's mean misses $\rho_*$ by
$1.3\cdot10^{-3}$, with a 4.6 % ripple. In the shifted frame every shock moves at least $0.36$, the ripple is 0.47 %
and the mean error $3.8\cdot10^{-5}$. A defect of the solvers would move with the problem; this one does not, so the
legs use the shifted frame, and the figure shows both.

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/eos/make_eos.py` | writes the inputs (`--eos barotropic --densities --exponents --rho-ref`, `--left`/`--right` states, `--rho0`) |
| `src/tests/flume/verification/eos/check.sh` | runs the legs (`--legs eb1,...,eb8`, default all but eb8 and eo7; `FLUME_EXE` selects the backend) |
| `src/tests/flume/verification/eos/eos_oracle.py` | the exact barotropic Riemann solver (convex tables), the linear wave, the comparisons: `riemann` (with the break plateaus), `wave`, `agree`, `selfconv` |
| `src/tests/flume/regression/baro-tube/input.ini`, `baro-collide-hllc/input.ini` | the frozen eb1 tube (WENO5 splitting) and eb2 collision (`weno-riemann` HLLC), 100 cells in four blocks |
| `make_doc_runs.sh <runs> eos`, `make_doc_figures.py --only barotropic-eos` | the figure (10 runs, 20 minutes on one CPU node) |

```bash
cd src/tests/flume/verification/eos && ./check.sh --legs eb1,eb2,eb3,eb4,eb5,eb6,eb7      # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh --legs eb1,eb2,eb3,eb7      # FNL
./check.sh --legs eb8                                                                 # CPU against FNL
```

## Results

![Barotropic equation of state: the tube across the breaks, the collision in two frames, the order and the default table](/flume/barotropic-eos.png)

*Top left: the tube across the breaks of `POLY`; the inset shows the plateau at the break $\rho = 8$ between the two
fans. Top right: the collision under HLL, the density behind the left shock in the slow-shock frame and in the
shifted one. Bottom left: the acoustic wave on the $\Gamma = 1.4$ segment, fifth order. Bottom right: the default
table, its fan through the non-convex break, converging to the finest run.*

| Leg | Asserted | Bound | CPU | FNL |
|---|---|---|---|---|
| eb1 | the plateau mean; the break plateau at $\rho = 8$; the right shock; the contact; `rE` | $10^{-3}$; $10^{-3}$; 1.5 cells; 2 cells; 0 | $6.0\cdot10^{-5}$; $6.6\cdot10^{-6}$; 0.37; 0.29; 0 | the same |
| eb2 | the layer mean against $\rho_* = 8.056$; both shocks | $10^{-3}$; 1.5 cells | $8.1\cdot10^{-6}$; 0.13, 0.20 | the same |
| eb3 | the L1 order of the finest pair | $> 4.5$ | 4.71 (4.94 coarser) | the same |
| eb4 | the 8 invalid tables and inputs | refused | refused | refused |
| eb5 | the table $\Gamma = 1$ against the isothermal gas, per field over its scale | $10^{-10}$ | $1.9\cdot10^{-14}$ | $2.9\cdot10^{-14}$ |
| eb6 | the L1 distance to the $N = 1600$ run, order per pair; $\rho > 0$; `rE`; the total mass | $\ge 0.5$; —; 0; $10^{-13}$ | 1.21, 1.58; $\min\rho = 10^{-8}$; 0; $1.2\cdot10^{-16}$ | the same orders; $2.4\cdot10^{-16}$ |
| eb7 | eb2 with `llf`, `hll`, `hllc` | as eb2 | layer $2.6\cdot10^{-5}$, $3.8\cdot10^{-5}$, $3.8\cdot10^{-5}$; shocks $\le 0.21$ cells | the same |
| eb8 | FNL against CPU, per field over its scale; `rE` | $10^{-10}$; 0 | (the same comparison) | $9.1\cdot10^{-13}$ (eb7 `hllc`, `ru`); 0 |

Related cases: the [isothermal gas](./isothermal-eos) is the one-segment table $\Gamma = 1$.

## Findings

- **A rarefaction through a convex break splits, and the code keeps the plateau.** The fan from $\rho = 10$ to
  $\rho_* = 2.6$ meets the break $\rho = 8$, where $c$ jumps up: the exact solution has two fans with a constant state
  at the break density between them, and the run holds it to $6.6\cdot10^{-6}$ over 24 of its 34 cells.
- **The slowly moving shock, not the solver.** A Galilean shift of $-0.6$ cuts the HLL ripple behind the left shock
  from 4.6 % to 0.47 % and the plateau error from $1.3\cdot10^{-3}$ to $3.8\cdot10^{-5}$: the noise is the
  frame-dependent slowly-moving-shock anomaly, worst for the least dissipative fluxes (LLF: 1.5 %).
- **The default table converges where no exact solution exists.** Across the non-convex break of Mathew & Federrath
  (2021) the runs approach the finest one at orders 1.21 and 1.58, with the density kept at $10^{-8}$ and the
  mass to $2.4\cdot10^{-16}$.
