# Barotropic gases under gravity

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| EO (eg1 to eg10) | Euler | M13 ([#100](https://github.com/szaghi/adam/issues/100)), P4 | CPU, FNL | yes | `baro-jeans-box`, `iso-sinks-shu` |

## Case and scope

The [isothermal](./isothermal-eos) and [barotropic](./barotropic-eos) gases have no energy equation: their energy slot
is inert (D-M13-1). P4 couples them to gravity, self-gravity, the Jeans refinement marker and the sink particles
([numerics](../../numerics#isothermal-and-barotropic-equations-of-state),
[input](../../input#eos-optional-equation-of-state)). Gravity must then move the momentum only, and every quantity the
ideal gas reads from its internal energy (the sound speed of the Jeans length, the sinks' threshold density, the
thermal energy of their Jeans test) must come from the equation of state, $c_s^2 = dP/d\rho$ (D-M13-4).

What the legs catch:
- **gravity writing into the energy slot:** the free fall and the AMR collapse keep `rE` and its integral at exactly 0;
- **a sound speed taken as $P/\rho$:** the Jeans wave and the Jeans marker on a barotropic background where
  $dP/d\rho$ and $P/\rho$ differ by the factor $\Gamma$;
- **a broken self-gravity coupling:** the order of the stable Jeans wave, the growth rate of the unstable mode, the
  mirror symmetry, the conservation across the 2:1 seams;
- **sinks reading the energy slot:** with $c_s^2 = 0$ every cell would pass the threshold, with $E_\mathrm{th} = 0$
  every candidate the Jeans test;
- **the isothermal sphere:** its accretion rate against Shu's similarity solution, now with an exactly isothermal gas;
- **the stability of a Bonnor–Ebert sphere** on either side of its critical radius;
- **a limiter that extrapolates the velocity:** the Boss–Bodenheimer cloud's envelope crossing the corners of the
  refined cube, where seam ghosts and cells reach the density floor;
- **a backend that disagrees.**

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Model | `euler` with `[eos] type = isothermal` or `barotropic` (`MODEL_EULER_ISO`, `MODEL_EULER_BARO`) | [`adam_flume_physics_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_physics_object.F90) (`eos_pressure`, `eos_sound_speed2`) |
| Gravity | the momentum sources of the ideal gas, the energy residual cleared (`clear_inert_energy`), no fourth-order energy correction, no seam gravity work | `adam_flume_cpu_object.F90`, `adam_flume_fnl_object.F90`, `adam_flume_fnl_kernels.F90` |
| Self-gravity | FFT (periodic, isolated) and multigrid (AMR) Poisson solvers, the potential path | [self-gravity](../self-gravity/) |
| Sinks | threshold $\pi c_s^2/(4G r_\mathrm{acc}^2)$ with $c_s^2 = dP/d\rho$; Jeans test with $E_\mathrm{th} = \tfrac32\int P\,dV$ | `adam_flume_common_object.F90` (`accrete_cells`, `create_sinks`) |
| Positivity | the cell limiter on the density, plus a kinetic chord (`internal_energy` = $\rho$ and `kinetic_sound_speed2` in the models' kernels; the reflux limiter's `cs2_inert`); the seam fill clamps the density and the velocity of a fine ghost to an MP bound of the coarse cells (#101) | `adam_flume_cpu_euler_{iso,baro}_kernels.F90`, `adam_flume_common_object.F90` (`reflux_positivity_theta`), `adam_seam_interpolation_library.F90` (`seam_mp_bound`), `adam_field_object.F90` (`interp_seam_ghost`) |
| Space, time | WENO5 characteristic flux splitting, SSP-RK(3,3), CFL 0.4 | — |

## Initial and boundary conditions

The barotropic cases use `POLY` of `eos_oracle.py` ($c_0 = 1$, breaks $\rho = 2, 8$, exponents $1, 1.4, 5/3$).

| Leg | Setup |
|---|---|
| eg1 | a uniform periodic box at rest (2 blocks of 24 cells), $g_x = 0.5$, $t = 0.4$; isothermal $\rho = 1$, barotropic $\rho = 5$ |
| eg2 | the stable Jeans wave of `make_jeans.py` (4×4 blocks, $\mathbf k = 2\pi(1, 1)$, amplitude $10^{-6}$), $4\pi G\rho_0 = c^2k^2/4$, one period, 32², 64², 128²; isothermal $c_0 = 1$, $\rho_0 = 1$; barotropic $\rho_0 = 5$: $c^2 = P'(5) = 2.0198$, while $P/\rho = 1.4427$ |
| eg3 | the unstable mode, $4\pi G\rho_0 = \tfrac94 c^2k^2$ ($\sigma = \sqrt{5/4}\,ck$), 64², $t = 0.5/c$, both gases |
| eg4 | eg2 (isothermal) with the central box refined 2:1, the multigrid solver; the barotropic collapse (amplitude $10^{-2}$, 64² base, $t = 0.4/c$) on the same grid |
| eg5 | a static periodic box $[-1, 1]^3$, 32³, a Jeans marker of 5.4 cells per Jeans length, a regrid at step 4: isothermal at 5 cells per $\lambda_J$; barotropic $\rho_0 = 10$ (the $\Gamma = 5/3$ segment) at 6 cells by $dP/d\rho$, 4.65 by $P/\rho$ |
| eg6 | the sinks' box (32³, $G = 1$): a Jeans-stable periodic box at 1.2 times the threshold ($c_0^2 = 0.0259$, creation on, 20 steps); accretion by a sink of mass 0.1 from $\rho = 1$, $c_0^2 = 2\cdot10^{-3}$; the singular isothermal sphere $A = 29.3$, $c_0 = 0.5$, a seed sink, $t = 0.25$; a supercritical Bonnor–Ebert sphere ($\xi = 10$ at the walls, $c_0 = 0.2$, density doubled), creation on, $t = 3$ |
| eg7 | the closed cube $[-1, 1]^3$, 32³, isolated self-gravity, $c_0 = 1$, $G = 1$, the Bonnor–Ebert profile with $\xi = 3, 6$ (sub-critical) and $9$ (super-critical) at the walls, $t = 4$ |
| eg8 | `bonnor-ebert` with the barotropic gas or without self-gravity, `boss-bodenheimer` with the ideal gas, `singular-isothermal-sphere` with the barotropic gas |
| eg9 | the eg1, eg2 (64²) and eg4 runs and the eg6 sphere and Bonnor–Ebert runs of CPU and FNL |
| eg10 | the Boss–Bodenheimer cloud of the figure (32 cells per block, the Jeans marker to 2 levels more, sinks), 300 steps (iteration-driven), still advancing at the end |

**The Bonnor–Ebert sphere in a cube.** The Lane–Emden profile is the equilibrium of an isolated sphere; in a closed
cube the gas in the corners, outside the inscribed sphere, carries 46 % of the sphere's mass at $\xi = 3$, and its
cubic multipoles pull on the sphere by 3 % of its gravity at $r = 0.6$ and 10 % at $r = 0.9$ (outward along the axes,
inward along the diagonals; direct summation over the 48³ profile). The profile is therefore not a discrete
equilibrium, and refining does not change that (16³ and 32³ drift alike); an isothermal gas offers no confining
medium either ($P(\rho)$ is single-valued). The leg tests what survives: a sub-critical sphere rings about its
initial state, a super-critical one runs away.

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/eos/make_eos_gravity.py` | writes the cases from the gamma-law Jeans and sinks inputs, the gas of `[eos]` in place of $\gamma$ |
| `src/tests/flume/verification/eos/gravity.sh` | runs the legs (`--legs eg1,...,eg10`, default all but eg9; `FLUME_EXE`, `EG_NP` ranks) |
| `eos_oracle.py` (`freefall`, `refined`, `inert-history`, `advancing`), `selfgravity/jeans_oracle.py`, `selfgravity/star_oracle.py` (`unstable`), `sinks/si_oracle.py` | the checks; the Jeans oracle reads $c^2$ from `[eos]` |
| `src/tests/flume/regression/baro-jeans-box/input.ini`, `iso-sinks-shu/input.ini` | the frozen eg4 collapse (20 steps) and eg6 sphere (16³, 30 steps) |
| `make_doc_runs.sh <runs> eosgrav`, `make_doc_figures.py --only barotropic-gravity` | the figure |

```bash
cd src/tests/flume/verification/eos && ./gravity.sh                                              # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl_fft EG_NP=2 ./gravity.sh                        # FNL
./gravity.sh --legs eg9                                                                          # CPU against FNL
```

## Results

![Barotropic gases under gravity: the Jeans wave and mode, Bonnor-Ebert stability, the singular isothermal sphere and the Boss-Bodenheimer cloud](/flume/barotropic-gravity.png)

*Top left: the stable Jeans wave, the L1 error of $\rho$ after one period, isothermal, barotropic and across the 2:1
seams (slope 2 dashed). Top centre: the unstable mode against $e^{\sigma t}$. Top right: the central density of the
Bonnor–Ebert spheres, sub-critical $\xi = 3, 6$ ringing, super-critical $\xi = 9$ departing. Bottom left: the sink of
the singular isothermal sphere against Shu's rate. Bottom centre and right: the Boss–Bodenheimer cloud, the column
density along $z$, with the sinks.*

| Leg | Asserted | Bound | CPU | FNL |
|---|---|---|---|---|
| eg1 | $\rho$ uniform, $\rho u - \rho g t$; `rE` | $10^{-12}$; 0 | 0, $8.9\cdot10^{-17}$; 0 | 0, $4.4\cdot10^{-17}$; 0 |
| eg2 | the L1 order of the finest pair, isothermal; barotropic | $\ge 1.8$ | 2.00; 2.00 | 2.00; 2.00 |
| eg3 | the growth rate; the mirror symmetry | $10^{-2}$; $10^{-10}$ | $7.2\cdot10^{-4}$ both; $7.3\cdot10^{-11}$, $4.2\cdot10^{-11}$ | $7.2\cdot10^{-4}$ both; $1.4\cdot10^{-11}$, $1.5\cdot10^{-11}$ |
| eg4 | the order across the seams; the momentum, the mass; the energy integral | $\ge 1.8$; $10^{-14}$, $10^{-13}$; 0 | 1.98; $5.6\cdot10^{-17}$, $1.8\cdot10^{-16}$; 0 | 1.98; $2.8\cdot10^{-17}$, $1.8\cdot10^{-16}$; 0 |
| eg5 | isothermal refined, barotropic not | — | 8 → 64 blocks; 8 → 8 | 8 → 64 blocks; 8 → 8 |
| eg6 | no sink; accreted mass, momentum; Shu's rate; one central sink | 0; $10^{-13}$, 5 %; 6 %; 1 cell | 0; 0, 2.6 %; 4.1 %; $8.9\cdot10^{-15}$ | 0; $2.2\cdot10^{-16}$, 2.8 %; 4.1 %; $3.8\cdot10^{-15}$ |
| eg7 | sub-critical drift and Mach; super-critical departure | 0.25, 0.25; $\ge 0.5$ | 0.147, 0.084 ($\xi = 3$), 0.136, 0.147 ($\xi = 6$); 0.695 | the same to the printed digits |
| eg8 | the 4 invalid inputs | refused | refused | refused |
| eg9 | FNL against CPU | $10^{-10}$, $10^{-6}$ of the amplitude, $10^{-8}$ | (the same comparison) | free fall $6.7\cdot10^{-17}$; Jeans $2.0\cdot10^{-9}$, $1.4\cdot10^{-9}$, collapse $2.9\cdot10^{-13}$; sinks $8\cdot10^{-11}$, $4\cdot10^{-16}$ |
| eg10 | the time advanced over the last 50 steps, over 50 mean steps; sinks formed | $\ge 0.1$; $\ge 1$ | 0.372 ($t = 1.649$); 1 | 0.372 ($t = 1.649$); 1 |
| eg10 | the edge cells of the shell against the rest: largest speed; smallest density ([#101](https://github.com/szaghi/adam/issues/101)) | $\le 2$; $\ge 0.5$ | 1.42; 1.00 (8 ranks) | 1.42; 1.00 (2 ranks) |

Related cases: the [Jeans modes](../self-gravity/jeans) and the [sink particles](../self-gravity/sink-particles) of the
ideal gas.

## Findings

- **The sound speed decides the Jeans physics, not the pressure.** On the barotropic background $\rho_0 = 5$ the wave
  is checked against $\omega^2 = P'(\rho_0)k^2 - 4\pi G\rho_0$; with $P/\rho$ in place of $P'$ its phase after one
  period would be off by 1.3 rad, and the marker of eg5 would refine a box it leaves alone.
- **The sinks need the virial thermal energy.** The ideal gas's internal energy $\int P/(\gamma - 1)\,dV$ diverges as
  $\gamma \to 1$; $\tfrac32\int P\,dV$ is the virial term, equal at $\gamma = 5/3$. With it the Jeans-stable box forms
  no sink, the supercritical sphere one at its centre, and the singular isothermal sphere, now exactly isothermal,
  accretes at 0.959 of Shu's rate (95.9 %), on both backends.
- **Accretion needs the density-only limiter.** A sink draining a cold background to its threshold density leaves a jump
  the scheme undershoots: without the limiter the run fails (the ideal-gas control at the same sound speed fails too).
  With it, mass is conserved to round-off and the momentum to 2.6 % of the exchange (FNL 2.8 %). The limiter makes the
  cold accretion run differ between backends: equal to round-off for five steps, apart by $10^{-6}$ at step 15, as the
  other limited cases (so eg9 leaves it out).
- **A limiter must not extrapolate the velocity of an inert gas.** When the cloud's envelope, expanding into an ambient
  100 times thinner, reaches the corners of the refined cube (step 130), the tricubic seam ghosts of the fine-block
  corner cells come out with $\rho \le 0$. Two limiters then take states to the density floor along a segment whose far
  end has negative density, the seam ghost blend (to $\kappa\rho_\mathrm{in}$) and the cell limiter (to $\kappa\rho_b$):
  the velocity $\rho\mathbf{u}/\rho$ of the result is extrapolated, up to $1/\kappa = 10$ times the anchor's. For the
  ideal gas the internal-energy chord rejects that state; for the inert gases nothing did, the corner cells' velocity
  grew by a factor per step, and $dt \to 0$. The ghost blend alone stalled the run at $t = 1.3721$ (2940 steps,
  $1.24\,t_\mathrm{ff}$); with the ghosts taking their interior cell, the cell limiter stalled it at $1.4637$, the
  ambient beside the corners drained to $10^{-5}$ at $|u|$ up to 365. The cell and reflux limiters now add a kinetic
  chord on $\rho W - |\mathbf{m}|^2/(2\rho)$, $W = (4|\mathbf{u}_b|^2 + c^2)/2$ of the backbone
  ([numerics](../../numerics#isothermal-and-barotropic-equations-of-state)): the cloud runs 400 steps to $t = 1.83$
  ($1.65\,t_\mathrm{ff}$), still advancing (eg10): the ambient dips once to $1.4\cdot10^{-4}$ and ends above
  $6.6\cdot10^{-4}$, and the fastest gas at the end, $|u| = 4.4$, falls into the sink. The chord's factor 4 matters for cold gas: with $W = (|\mathbf{u}_b|^2 + c^2)/2$ the cold accretion of eg6
  drifted its momentum by 13 % of the exchange.
- **The seam fill must bound the velocity at the corners of a refined region** ([#101](https://github.com/szaghi/adam/issues/101)).
  The chord keeps the run going, it does not keep the corner ghosts sound. eg10 compares, every 20 steps, the cells of
  the shell $1.55 < r < 1.85$ near the edges of the refined cube (3680 cells) with the rest of the shell (42216): their
  largest speed (bound: twice the rest's) and their smallest density (bound: half the rest's). With the tricubic ghosts
  unbounded the edge cells reach $|u| = 22.7\,c_0$, 14.8 times the rest (step 200), and their density falls to 0.04
  times the rest's (step 240). Three fills were measured before the one adopted (400 steps, 8 ranks):
  a density guard (a ghost below 0.1 times its anchor's density takes the anchor's state) holds the density (0.90) but
  not the speed (2.95); clamping the density and the momenta to a monotonicity-preserving bound holds the density
  (1.00) but not the speed either (2.93): the velocity, their ratio, stays free; injection everywhere passes (1.07,
  1.00) at first order. The fill adopted clamps the density and the **velocity** to the MP bound of the coarse cells
  around the ghost ([numerics](../../numerics#adaptive-mesh-refinement)): **1.42 and 1.00**, the same on CPU (400 steps)
  and FNL (eg10, 300 steps, 2 GPUs), while the order across a 2:1 corner on smooth data is unchanged (AO `euler-l6`:
  L1 8.548e-8 at $N = 256$, orders 5.33 and 5.53, CPU and FNL to every printed digit).

![Corner cells against the rest of the shell, five seam fills](/flume/eos-corner-guard.png)
- **A Bonnor–Ebert sphere in a cube is no equilibrium, but its stability is testable.** See above; sub-critical
  spheres ring within 15 % of their central density (0.147 at $\xi = 3$, 0.136 at $\xi = 6$) below Mach 0.15,
  the super-critical one departs by 70 % (0.695; $-0.67$ at $t = 4$).
