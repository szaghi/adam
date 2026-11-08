# Plane channel and the wall model (T3)

| ID | Milestone | Type | Legs | Backends | Figure |
|---|---|---|---|---|---|
| WM | [#119](https://github.com/szaghi/adam/issues/119) · T3.1–T3.6 | P+C | wm1–wm11 | CPU, FNL | yes |

## Case and scope

T3 adds the channel infrastructure (mass-flow forcing, isothermal walls, T3.1), the `wall-modelled` boundary (T3.3) with
the equilibrium model (T3.4) and the inverse-velocity-transformation model (T3.5), its input filter, the eddy viscosity
below the matching height, and 2:1 faces along a modelled wall (T3.6). The methods are on the
[turbulence](../../turbulence#wall-model) page. The legs check each piece against a laminar exact solution, a
restart or a rank change. The campaign runs the compressible channel of Coleman, Kim & Moser (1995) at bulk Mach 1.5
and 3, against a wall-resolved reference (T3.2) and the DNS.

## Legs

`src/tests/flume/verification/turbulence/check.sh`, CPU and FNL.

| Leg | Phase | What it verifies | Measured |
|---|---|---|---|
| wm1 | T3.1 | laminar channel, isothermal walls, Re 100: force against $3/Re$; one rank against 2 | order 1.99, $4.8\cdot10^{-4}$ at 64 across; $4\cdot10^{-11}$ |
| wm2 | T3.3, T3.4 | `test_wall_model`; laminar channel on modelled walls | see [unit tests](../unit-tests); order 1.84, 1.9 % at 64 |
| wm3 | T3.6 | half the wall refined, 2:1: force between the coarse and the fine uniform runs | 0.03 % from their mean, mass $10^{-16}$ |
| wm4 | T3.3 | hot wall, modelled against resolved: energy change; one rank against 2 on modelled walls | 7.8 % (bound 10 %); $4\cdot10^{-11}$ |
| wm5 | T3.4 | matching height 0.3, laminar channel at Re 5 against the exact slip solution $3/Re/(1+3h^2/2)$; semi-local against classical damping | $6.4\cdot10^{-3}$, $2.2\cdot10^{-3}$, $3.0\cdot10^{-5}$ at 16, 32, 64; $4\cdot10^{-7}$ |
| wm6 | T3.5 | semi-local against classical damping on a hot wall, Re 100 | force 0.7 % apart (min 0.01 %) |
| wm7 | T3.5 | input filter: $T_f\to0$ against no filter; $T_f=0.05$ near the steady state; $T_f=1$ over a transient; one rank against 2 | 0 (bit for bit); $3.0\cdot10^{-4}$ (bound $10^{-3}$); $1.4\cdot10^{-2}$ (min $10^{-3}$); $4\cdot10^{-11}$ (CPU), $10^{-13}$ (FNL); CPU against FNL $5\cdot10^{-11}$ |
| wm8 | T3.5 | inverse-transformation model: laminar channel against the equilibrium model; hot wall; one rank against 2 with the filter | $1.0\cdot10^{-4}$, the heat flux of the Reynolds analogy moving the energy by $2.3\cdot10^{-3}$ (bound $10^{-4}$); $9.8\cdot10^{-3}$ (min $10^{-3}$); $4\cdot10^{-11}$ (CPU), 0 (FNL); CPU against FNL $4\cdot10^{-11}$ |
| wm9 | T3.5 | restart with the filter averages (`.wmf`, one per rank): 10 steps + restart to 20 against 20 straight, inverse-transformation model, $T_f=0.01$; control with the `.wmf` removed | fields and histories bit for bit (CPU and FNL); the control differs and warns |
| wm10 | T3.1 ([#49](https://github.com/szaghi/adam/issues/49) NV-5, NV-9) | dimensional units on a Vreman channel with the wall model and the forcing, gas `cp, cv`: input bit for bit (NV-5), output against `scaling.py dimensionalize` (NV-9), forcing and turbulence histories included | bit for bit (CPU and FNL); bit for bit, 4 histories |
| wm11 | Mach 3 channel | the wall model's eddy viscosity below the matching height (`wall_model_eddy_viscosity`), perturbed channel at Re 1000, $h = 0.3$: it acts; 2 against 4 blocks per direction; one rank against 2; the same with Vreman | not recorded (bounds: force change above $10^{-3}$; blocks within `channel_oracle agree`'s bounds, $2.4\cdot10^{-6}$ against $10^{-9}$ with the ghost columns left out) |

## Mach 1.5

| Setup | |
|---|---|
| case | Coleman, Kim & Moser (1995): bulk Mach 1.5, bulk Re 3000 on the half height and the wall viscosity |
| box | $4\pi\times2\times\frac43\pi$, periodic in $x$ and $z$, mass-flow forcing |
| gas, walls | viscosity $(T/T_w)^{0.7}$, Pr 0.7; isothermal walls |
| scheme | `hybrid-kep-weno` |
| reference | $Re_\tau$ about 220–222 for the DNS, as cited by Qi et al. (2022), Cao et al. (2024), Alberti et al. (2025b) |
| run | `turbulence/channel_campaign.sh` with `channel_cases_t35.txt`, 4 A30 |
| measures | $\tau_w$: the mean force per unit volume times $h$; $Re_\tau=3000\sqrt{\tau_w\rho_w}$ with $\rho_w = 1.30$ for every row, so the rows compare among themselves (the reference run's own wall density, extrapolated from its last checkpoint, is 1.25, which lowers every $Re_\tau$ by 2 %); $u_c$, $T_c/T_w$: the Favre plane means at the centreline of the last checkpoint |
| inverse model | edge state given: $u_e$, $T_e/T_w = 1.29$, the centreline of the $h=0.15$ semi-local runs |

![The plane channel](/flume/channel-wm.png)

*$\tau_w$ against time; dashed: $Re_\tau$ 220 and 222; dotted: the laminar value.*

| Case | Cells | $\tau_w$, $100<t\le150$ | $150<t\le200$ | $Re_\tau$ | $u_c$, $T_c/T_w$ |
|---|---|---|---|---|---|
| resolved walls | $256\times192\times128$, $\Delta y^+\approx2.3$ | 3.82e-3 | – | 211 ($50<t\le100$: 218) | 1.16, 1.26 |
| resolved walls | $64\times32\times32$ | 2.35e-3 | 2.72e-3 | 166, 178 | 1.22, 1.26 |
| equilibrium, 3rd cell | $64\times32\times32$ | 4.37e-3 | 4.49e-3 | 226, 229 | 1.21, 1.29 |
| equilibrium, 3rd cell | $96\times48\times48$ | 4.36e-3 | 4.52e-3 | 226, 230 | 1.19, 1.30 |
| equilibrium, 3rd cell | $128\times64\times64$ | 4.69e-3 | 4.66e-3 | 234, 234 | 1.18, 1.31 |
| equilibrium, 3rd cell, semi-local | $128\times64\times64$ | 4.33e-3 | 4.33e-3 | 225, 225 | 1.18, 1.30 |
| equilibrium, $h=0.15$ | $96\times48\times48$ | 4.71e-3 | 4.77e-3 | 235, 236 | 1.19, 1.30 |
| equilibrium, $h=0.15$ | $128\times64\times64$ | 4.89e-3 | 4.90e-3 | 239, 239 | 1.19, 1.30 |
| equilibrium, $h=0.15$, semi-local | $96\times48\times48$ | 4.15e-3 | 4.09e-3 | 220, 219 | 1.18, 1.29 |
| equilibrium, $h=0.15$, semi-local | $128\times64\times64$ | 4.13e-3 | 4.21e-3 | 220, 222 | 1.17, 1.29 |
| … and the filter $T_f=1$ | $96\times48\times48$ | 4.15e-3 | 4.03e-3 | 220, 217 | 1.18, 1.29 |
| … and the filter $T_f=1$ | $128\times64\times64$ | 4.15e-3 | 4.16e-3 | 220, 221 | 1.19, 1.29 |
| inverse, $h=0.15$, $u_e=1.18$ | $128\times64\times64$ | 4.13e-3 | 4.08e-3 | 220, 218 | 1.16, 1.21 |
| … and the filter $T_f=1$ | $128\times64\times64$ | 4.04e-3 | 4.11e-3 | 217, 219 | 1.17, 1.21 |
| inverse, $h=0.15$, $u_e=1.18$ | $96\times48\times48$ | 2.03e-3 | 1.62e-3 | 154, 138: relaminarises | 1.51, 1.11 |
| inverse, 3rd cell, $u_e=1.18$ | $96\times48\times48$ | 1.61e-3 | 1.40e-3 | 137, 128: relaminarises | 1.44, 1.12 |
| inverse, $h=0.15$, $u_e=1.10$ | $96\times48\times48$ | 3.61e-3 | 4.15e-3 | 206, 220 | 1.18, 1.21 |

Every run conserves mass to $10^{-16}$.

**What it shows.**
- **Matching height 0.15 with semi-local damping puts the equilibrium model in the DNS band on both grids** (219–222).
  Either change alone is worse: the 3rd-cell sample is 226–234 and rises with the grid, $h=0.15$ with the classical
  damping 235–239, the semi-local damping at the 3rd cell 225. The fixed height removes the grid trend; the semi-local
  damping removes the excess. This is one Mach and Reynolds number, single runs, and the band sits within the $\pm2\%$
  of the wall density.
- **The input filter does not move the means** (within 1.5 %); on the coarser grid it slows the transient.
- **The inverse model is in the band on the finer grid** with a cooler centreline ($T_c/T_w$ 1.21 against 1.29 for the
  equilibrium model, 1.26 for the resolved run): its heat flux differs.
- **On $96\times48\times48$ the inverse model relaminarises with $u_e = 1.18$**, at $h=0.15$ and at the 3rd cell, but
  not with $u_e=1.10$. The model is sensitive to the given edge velocity, which the flow does not set (the paper takes
  the edge from the flow; a face kernel cannot).
- **The wall-resolved reference is 4–5 % short of the DNS** (211). The modelled rows in the band agree with the DNS, not
  with it.
- **The coarse grid with resolved walls falls 20 % short**: without a model the grid cannot carry the wall layer.

## Mach 3

Coleman, Kim & Moser's second case: bulk Mach 3, bulk Re 4880, $Re_\tau$ about 451 for the DNS; the rest as at Mach 1.5,
the equilibrium model only.

**Result.** With the wall model's eddy viscosity below $h$ (`wall_model_eddy_viscosity`) and run to a thermal steady
state at $t=800$, $h=0.15$ with semi-local damping gives $Re_\tau = 443$–$445$ against about 451 for the DNS; the
classical damping overshoots by 13–21 %. Without the key the first cell drains, as the runs to $t = 200$ show.

### Without a flux path below $h$, to $t = 200$

| Setup | |
|---|---|
| run | `channel_cases_t3_mach3.txt`, to $t = 200$; the first-cell history from `channel_cases_t3_mach3_drain.txt` (the same runs, bitwise, a checkpoint every 500 steps) |
| $\rho_w$, $Re_\tau$ | the wall density changes with the run at this Mach number, so each row takes its own: $\rho_w$ is the first-cell mean pressure of the last checkpoint over $T_w$, $Re_\tau = 4880\sqrt{\tau_w\rho_w}$ |
| $T_1/T_w$ | the first cell (Favre, both walls) |
| $Q_w/W$ | wall heat flux over forcing work for $180<t\le200$, from the total-energy budget (`channel_wall_cells.py budget`; $Q_w$ is the residual): 1 at a thermal steady state |

![The plane channel at Mach 3](/flume/channel-m3.png)

*Left: $\tau_w$ against time (dotted: the laminar value). Middle and right: the first cell of the two $h=0.15$ runs on $128\times64\times64$ against
time, from a checkpoint every 500 steps.*

| Case | Cells | $\tau_w$, $100<t\le150$ | $150<t\le200$ | $\rho_w$ | $Re_\tau$ | $u_c$, $T_c/T_w$ | $T_1/T_w$ | $Q_w/W$ |
|---|---|---|---|---|---|---|---|---|
| equilibrium, 3rd cell | $96\times48\times48$ | 1.77e-3 | 1.43e-3 | 1.65 | 264, 238: relaminarises | 1.43, 1.63 | 0.86 | 0.71 |
| equilibrium, $h=0.15$, semi-local | $96\times48\times48$ | 1.89e-3 | 1.59e-3 | 1.54 | 263, 241: relaminarises | 1.50, 1.53 | 0.64 | 0.67 |
| equilibrium, 3rd cell | $128\times64\times64$ | 4.29e-3 | 4.44e-3 | 2.12 | 465, 473 | 1.25, 2.22 | 1.03 | 0.87 |
| equilibrium, 3rd cell, semi-local | $128\times64\times64$ | 3.20e-3 | 3.03e-3 | 2.02 | 392, 381 | 1.20, 2.04 | 1.28 | 0.85 |
| equilibrium, $h=0.15$, semi-local | $128\times64\times64$ | 3.58e-3 | 3.46e-3 | 1.91 | 404, 397 | 1.28, 2.01 | 0.81 | 0.83 |
| equilibrium, $h=0.15$ | $128\times64\times64$ | 6.46e-3 | 7.49e-3 | 2.12 | 571, 615: still rising | 1.48, 2.63 | 0.41 | 0.88 |

**What it shows.**
- **No run is at a thermal steady state at $t=200$**: the total energy still rises in every run ($Q_w/W$ 0.67–0.88), so
  these $Re_\tau$ are values on the way, not means.
- **The Mach 1.5 setting does not carry over**: $h=0.15$ with semi-local damping is 12 % below the DNS; the 3rd-cell
  sample is 5 % above. The coarse grid relaminarises with either setting.
- **The first cell is drained.** The wall face takes the modelled $\tau_w$ and $q_w$ from the sample at $h$; the first
  cell's own state never enters, and nothing below $h$ carries those fluxes (no subgrid model): the molecular plus
  resolved stress at the face between cells 1 and 2 is 0.54 $\tau_w$ with semi-local damping and 0.21–0.39 with the
  classical one (after $t=30$).
- **With $h=0.15$ and the classical damping the run does not settle.** The first cell falls to $0.26\,T_w$, reverses its
  flow and holds 17 % of the mass at $t\approx80$; $\tau_w$ follows it, down ($t$ 80–110) and up again.
- **With semi-local damping the drain is slower** and $\tau_w$ stays flat, but the first cell is still below $T_w$ from
  $t\approx125$. A first cell colder than an isothermal wall that takes heat out of it is not a boundary layer: the
  treatment below $h$ fails, not the model's law.

### With the eddy viscosity below $h$, to $t = 800$

| Setup | |
|---|---|
| key | `wall_model_eddy_viscosity = .true.`: $\mu_t = \max(\mu_{t,\mathrm{SGS}},\ \kappa y\sqrt{\rho\tau_w}\,D^2)$ in the cells whose centre is below $h$ ([turbulence](../../turbulence#eddy-viscosity-below-the-matching-height)) |
| runs | the three $128\times64\times64$ inputs above plus the key, to $t=200$ (`channel_cases_t3_mach3_ev.txt`), continued from the restart to $t=800$ (`channel_cases_t3_mach3_ev800.txt`) |
| measures | $\rho_w$ and $T_1/T_w$ at the last checkpoint as above; $Q_w/W$ over $700<t\le800$ |

![The plane channel at Mach 3 with the eddy viscosity below h](/flume/channel-m3-ev.png)

*Left: $Re_\tau$ against time, with the wall density of each checkpoint (dashed: the DNS). Middle: the first cell. Right: the wall heat flux over the
forcing work in 50-unit windows, 1 at a thermal steady state.*

| Case | $\tau_w$, $600<t\le700$ | $700<t\le800$ | $\rho_w$ | $Re_\tau$ | $u_c$, $T_c/T_w$ | $T_1/T_w$ | $Q_w/W$ |
|---|---|---|---|---|---|---|---|
| equilibrium, $h=0.15$, semi-local | 3.52e-3 | 3.48e-3 | 2.36 | **445, 443** | 1.17, 2.44 | 1.79 | 0.99 |
| equilibrium, 3rd cell | 4.51e-3 | 4.35e-3 | 2.51 | 519, 510 | 1.21, 2.65 | 1.57 | 0.99 |
| equilibrium, $h=0.15$ | 5.17e-3 | 5.15e-3 | 2.44 | 548, 547 | 1.21, 2.53 | 1.84 | 1.00 |

**What it shows.**
- **The drain is gone.** The first cell is $1.6$–$1.8\,T_w$ at every time after the transient, the two walls agree within
  2 %, and the two first-cell layers hold 4–5 % of the mass (3.1 % if uniform), against 16 % at $0.41\,T_w$ without the
  key.
- **The runs are thermally steady at $t=800$**: $Q_w/W$ is 0.99–1.00 and $dE/dt$ about 60 times smaller than at $t=200$.
- **$h=0.15$ with semi-local damping matches the DNS within 2 %.** With the classical damping $Re_\tau$ is 13 % (3rd
  cell) and 21 % ($h=0.15$) above.
- **$Re_\tau$ needs the thermal steady state at this Mach number.** The bulk heats until $t\approx600$, $\rho_w$ rises
  from about 2.1 to 2.4–2.5, and $Re_\tau$ rises with it: at $t=200$ the 3rd-cell run gave 438, 448, which looked like a
  match.
- **Cost.** From the same $t=800$ state of the semi-local run to $t=810$ (4 A30), the key on takes 1.69 times the wall
  time of the key off: 1.28 times the steps (the eddy viscosity enters the viscous time-step limit) and 1.32 times the
  time per step.

## The flow

![Instantaneous streamwise velocity in the two channels](/flume/channel-fields.png)

*The streamwise velocity over the bulk velocity of the last checkpoint: the $x$-$y$ plane at mid-span (left), the wall-parallel plane at the
matching height $h=0.15$ above the lower wall (middle, its own colour range) and the cross-flow plane (right). Top: bulk Mach 1.5,
$128\times64\times64$, $h=0.15$ with semi-local damping, $t=200$. Bottom: bulk Mach 3, the same grid and model with the eddy viscosity below $h$,
$t=800$. Planes from `turbulence/field_slices.py`; the figure is `make_doc_figures.py --only channel-fields`.*

- **At the matching height the flow is in streaks** along $x$, but outer-scaled ones. At Mach 3 the first minimum of the
  spanwise correlation of $u$ in that plane puts them $1.05\,h$ apart, about 470 wall units at $Re_\tau\approx445$,
  against the 100 of the near-wall cycle; at Mach 1.5 the correlation has no clear minimum within the span. The
  near-wall cycle lives below the matching height ($y^+\approx33$ at Mach 1.5), in the wall model, and the grid
  ($\Delta z^+\approx14$, $\Delta x^+\approx22$ at Mach 1.5) does not form it above.
- **At Mach 3 the streaks look longer and more regular** than at Mach 1.5. One snapshot each: a description, not a
  measurement.

## Not done

- Averaged profiles against the DNS: only $\tau_w$ and the centreline are compared, and the DNS profiles are not in this tree.
- The cause of the wall-resolved reference's 4–5 % shortfall: its resolution or the hybrid scheme's dissipation.
- The relaminarisation of the inverse model on $96\times48\times48$ with $u_e = 1.18$; an edge state taken from the flow.
- Why the classical damping overshoots at Mach 3; a cheaper flux path below $h$.
- Runtime regridding along a modelled wall.
