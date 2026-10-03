# Verification gallery

Every FLUME capability is backed by a verification test in `src/tests/flume/verification/`: a `check.sh` that runs the
case (CPU by default, `FLUME_EXE=exe/adam_flume_fnl` for the GPU backend) and an oracle that asserts the physics (an
exact solution, a reference, a symmetry or a conservation property), not a stored golden. This page shows the main
cases, with the input that defines them and the measured results. The figures are drawn from the runs by
`src/tests/flume/verification/make_doc_figures.py`; all runs use 2 MPI ranks and WENO5 unless stated, and the CPU and
GPU (FNL) backends agree to the printed digits (bitwise wherever the table says so) unless a section says otherwise.

Run any case from its directory:

```bash
cd src/tests/flume/verification/sod && ./check.sh                  # CPU backend
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh        # FNL backend (GPU)
```

## Euler

### Sod shock tube

The classical shock tube (Sod 1978): $(\rho, u, p) = (1, 0, 1)$ for $x \le 0.5$, $(0.125, 0, 0.1)$ beyond, $\gamma = 1.4$,
200 cells, $t = 0.2$, CFL 0.5, SSP-RK3. The same problem is run along $x$, $y$ and $z$ (the other two directions null).

![Sod shock tube](/flume/sod.png)

| Check | Result |
|---|---|
| $L_1(\rho)$ against the exact solution, `weno` (split, characteristic) | $3.244\cdot10^{-3}$ |
| $L_1(\rho)$, `weno-riemann` LLF / HLL / HLLC | $3.426$ / $3.013$ / $2.942\cdot10^{-3}$ |
| $y$ and $z$ runs against $x$ | **bitwise** equal after the permutation of the axes |
| reflecting-wall variant (`sod-wall-*`) | mirror symmetric, mass and energy of the closed box constant |

`sod/check.sh` (V1) and `riemann-flux/check.sh --leg sod` (RV-3). Input: `sod/sod-x.ini`
([riemann-problem](./initial-conditions#riemann-problem-piecewise-constant-axis-aligned-regions) initial condition,
extrapolation boundaries); the complete file is on the [input reference](./input) page.

### Lax problem

$(\rho, u, p) = (0.445, 0.698, 3.528)$ on the left, $(0.5, 0, 0.571)$ on the right, $t = 0.13$ on $[0, 1]$, 200 cells: a
stronger contact than Sod, where the dissipation of the solver shows.

![Lax problem](/flume/lax.png)

| Scheme | $L_1(\rho)$ |
|---|---|
| `weno` (split) | $1.064\cdot10^{-2}$ |
| `weno-riemann` LLF | $1.214\cdot10^{-2}$ |
| `weno-riemann` HLL | $1.125\cdot10^{-2}$ |
| `weno-riemann` HLLC | $1.081\cdot10^{-2}$ |

The acceptance bound of `weno-riemann` is $1.05\times$ the split scheme's (`riemann-flux/check.sh --leg lax`).

### Shu–Osher shock / entropy-wave interaction

A Mach 3 shock runs into a sinusoidal density field (Shu & Osher 1989): left state $(3.857143, 2.629369, 10.33333)$ for
$x \le -4$, $\rho = 1 + 0.2\sin 5x$, $u = 0$, $p = 1$ beyond, $t = 1.8$ on $[-5, 5]$, 400 cells
([shu-osher](./initial-conditions#shu-osher-shock-and-density-wave-interaction) initial condition). The reference is the
split scheme on 6400 cells, averaged to the grid.

![Shu-Osher](/flume/shu-osher.png)

| Scheme | $L_1(\rho)$ | max $\rho$ (reference 4.672) |
|---|---|---|
| `weno` (split) | 0.2605 | 4.563 |
| `weno-riemann` LLF | 0.3009 | 4.550 |
| `weno-riemann` HLL | 0.2490 | 4.590 |
| `weno-riemann` HLLC | **0.2425** | 4.592 |

`riemann-flux/check.sh --leg shu-osher`; the x, y and z runs are bitwise equal.

### Isentropic vortex: order of accuracy

An exact smooth solution: a vortex convected diagonally through a doubly periodic box for one period
([isentropic-vortex](./initial-conditions#isentropic-vortex) initial condition, strength 5, radius 0.07). The error is
measured against the exact solution on 64², 128² and 256² cells.

![Isentropic vortex](/flume/vortex.png)

| Scheme | $L_1(\rho)$, 256² | Order 128 → 256 |
|---|---|---|
| `weno` (split) | $8.38\cdot10^{-8}$ | 5.79 |
| `weno-riemann` LLF | $1.42\cdot10^{-7}$ | 5.64 |
| `weno-riemann` HLL | $8.98\cdot10^{-8}$ | 6.26 |
| `weno-riemann` HLLC | $8.04\cdot10^{-8}$ | 6.49 |

`vortex/check.sh` (V2) and `riemann-flux/check.sh --leg vortex` (RV-2, order at least 4.8). The order from 64² to 128²
is lower (about 4): at 64² the vortex spans 4.5 cells and the WENO weights are still far from the linear ones.

### Shock over a cylinder

A Mach 2 shock (inflow state at $x_{\min}$) hits a cylinder of radius 0.13 imposed by the immersed boundary; the tree is
refined at initialisation on the solid surface (176 blocks, 96 crossed by the surface, 32 coarse–fine faces); slip walls
at $y = 0$ and $y = 1$. The oracle checks the refinement, the mirror symmetry about $y = 0.5$ and the positivity of the
fluid.

![Shock over a cylinder](/flume/shock-cylinder.png)

| Check | `weno` (split) | `weno-riemann` HLLC |
|---|---|---|
| relative asymmetry about $y = 0.5$ | $1.2\cdot10^{-11}$ (FNL $4.8\cdot10^{-12}$) | $4.9\cdot10^{-12}$ (FNL $1.9\cdot10^{-12}$) |
| min $\rho$, min $p$ in the fluid | positive | 1.365, 0.968 |

`shock-cylinder/check.sh` (V6), `riemann-flux/check.sh --leg cylinder`. Input: `shock-cylinder/shock-cylinder.ini`
(`[solids]`, `[solid_1] definition = analytical_circle`, an `[amr]` solid marker).

### Conservation across AMR seams

A triply periodic box with one refined octant (6 coarse–fine faces, periodic ones included) carries a uniform flow with a
1% seeded perturbation. Nothing crosses the boundary, so the five volume integrals must stay constant: with the
Berger–Colella reflux they do, to round-off; without it (the negative control) they drift.

![Conservation across AMR seams](/flume/conservation.png)

| Run | Largest relative drift of the integrals |
|---|---|
| reflux on, `weno` | $\le 6\cdot10^{-15}$ (CPU), $2.2\cdot10^{-16}$ (FNL) |
| reflux on, `weno-riemann` HLLC | 0 (CPU), $2.2\cdot10^{-16}$ (FNL) |
| reflux off | $1$–$2\cdot10^{-5}$ |

`conservation/check.sh` (V3), `riemann-flux/check.sh --leg conservation` (RV-4); the MHD version (MV-11) adds the three
components of $\mathbf{B}$ and $\psi$ (drift $\le 4.7\cdot10^{-14}$).

### Multi-realm runs

`multirealm/check.sh` splits Sod along $x$ into two realms at the diaphragm, glued by a mirror seam with stage-coincident
cadence, and compares the union with the single-realm run: **bitwise** equal on all cells at $t = 0.2$, CPU and FNL. The
same holds with one realm refined, the shock crossing a 2:1 face. Leg 3 splits Sod along $z$ instead, where the two
sides of the seam sit on different ranks: bitwise equal to the single-realm run on 1, 2, 3 and 4 ranks (the seam ghosts
and the reflux register cross ranks, [#40](https://github.com/szaghi/adam/issues/40)).

## Ideal MHD

### Ryu–Jones 2a: all seven waves

The 1-D Riemann problem of Ryu & Jones (1995), case 2a, produces the full MHD wave fan: two fast shocks, two rotational
discontinuities, two slow shocks and a contact. Its exact solution is known (7 states); the problem is run along $x$,
$y$ and $z$ with 256 cells and along $x$ with 512.

![Ryu-Jones 2a](/flume/rj2a.png)

The figure overlays the splitting scheme and `weno-riemann` with HLLD (characteristic interpolation), which are
indistinguishable at this scale. $L_1$ below the recorded bound and decreasing from 256 to 512 cells (first order, as
the discontinuities dominate); the $y$ and $z$ runs equal the $x$ run **bitwise** in the rotated frame.
`mhd/rj2a/check.sh` (MV-4).

### Brio–Wu and Ryu–Jones 4d

The Brio–Wu shock tube ($\gamma = 2$, $B_x = 0.75$, $B_y = \pm 1$) develops a compound wave; Ryu–Jones 4d a switch-on
structure. Both must reach the final time with zero floored cells (`mhd/riemann/check.sh`); a negative-pressure state is
floored when the floors are armed and stops the run when they are not.

![Brio-Wu](/flume/brio-wu.png)

![Ryu-Jones 4d](/flume/rj4d.png)

### GLM divergence cleaning

A Gaussian pulse in $B_x$ along $x$ is a pure divergence error. The pair $(B_x, \psi)$ then obeys a (damped) telegraph
equation with an exact solution: the error splits into two pulses travelling at $\pm c_h$ and is damped at the rate
$\alpha c_h / L$.

![GLM pulse](/flume/glm-pulse.png)

The figure is the undamped case (`glm_alpha = 0`, periodic, 128 cells, $c_h = 2$, $t = 0.3$): FLUME on the exact solution.

`mhd/glm-pulse/check.sh` (MV-3) checks the $L_1$ error and the observed order with and without damping, the `min-cell`
damping length, $c_h$ in the time step and the reflection of $\psi$ at walls. `mhd/divb-peak/check.sh` (MV-10) checks the
decay of $\|\nabla\cdot\mathbf{B}\|_1$ from an initial peak of $\nabla\cdot\mathbf{B}$, on a uniform grid and across AMR seams.

### Order of accuracy: linear waves, CPAW, magnetised vortex

| Test | Case | Measured |
|---|---|---|
| MV-5 `mhd/linear-wave/` | the fast, Alfvén, slow and entropy eigenmodes, amplitude $10^{-7}$ | about 5th order for every family |
| MV-6 `mhd/cpaw/` | circularly polarised Alfvén wave (exact nonlinear solution), oblique in 2-D | orders 4.89 and 4.99 on 32/64/128; left and right polarisations **bitwise** equal |
| MV-7 `mhd/vortex/` | magnetised vortex of Balsara (2004), GLM | $L_1$ order 5.01, $L_\infty$ 4.41 (pre-asymptotic); div B at truncation level |

![MHD linear waves, order of accuracy](/flume/linear-wave-order.png)

The 1-D linear waves on 16, 32 and 64 cells per wavelength, on the two schemes; the legend gives the orders of the two
pairs. From 16 to 32 cells every family is at 4.96–4.97; at 64 cells the error is $10^{-13}$, six decades below the
amplitude, and the round-off floor lowers the slow and entropy orders to 4.85 (which is why MV-5 asserts on the 16/32
pair).

### Orszag–Tang vortex

The standard 2-D test of MHD turbulence and shock interactions ([orszag-tang](./initial-conditions#orszag-tang) initial
condition), 128², $\gamma = 5/3$, GLM, $t = 0.5$.

![Orszag-Tang](/flume/orszag-tang.png)

The oracle checks the 180° rotational symmetry (relative defect $1.3\cdot10^{-9}$ at $t = 0.5$: a round-off seed amplified
by the flow), the conservation of $\rho$, $\rho\mathbf{u}$, $E$ and $\mathbf{B}$ ($1.2\cdot10^{-14}$) and the positivity. With
`--amr` the centre box is refined once (2:1 seams crossed by the shocks); symmetry and conservation hold there too, and the
div(B) beside the seams stays below the uniform-grid peak (`mhd/orszag-tang/check.sh`, MV-12).

![Orszag-Tang on AMR](/flume/orszag-tang-amr.png)

### MHD rotor

A dense disc spins in a magnetised ambient gas (Balsara & Spicer 1999; the first rotor of Tóth 2000) and launches
torsional Alfvén waves ([mhd-rotor](./initial-conditions#mhd-rotor)). The oracle checks the symmetry under the 180°
rotation composed with $\mathbf{B} \to -\mathbf{B}$ (`mhd/rotor/check.sh`, MV-13).

![MHD rotor](/flume/rotor.png)

### Field loop advection

A weak magnetic loop (plasma $\beta \approx 10^6$) is advected across the periodic box
([field-loop](./initial-conditions#field-loop-advected-field-loop)). With a velocity component along $z$,
$\partial B_z/\partial t = w\,\nabla\cdot\mathbf{B}$, so the generated $\langle|B_z|\rangle$ measures the divergence error; the oracle
requires it at the GLM level and decreasing under refinement, and bounds the decay of the magnetic energy
(`mhd/field-loop/check.sh`, MV-9). The `--amr` variant crosses 2:1 seams: the seam div(B) peaks at the first step (the
initial loop edge) and GLM removes it.

![Field loop](/flume/field-loop.png)

### Rotated shock tube

Ryu–Jones 1a rotated to 63.4° and 45° on a periodic square (`mhd/rotated-shock-tube/`, MV-8): the exact normal field is
uniform, so its deviation measures the divergence error. At 256 cells $\delta B_\parallel = 2.79\cdot10^{-3}$ (63.4°) and
$4.17\cdot10^{-3}$ (45°), first order; Tóth's second-order base scheme reports $3.7\cdot10^{-3}$ at 63.4°.

### MHD on `weno-riemann`

Every MHD script above also runs on the `weno-riemann` scheme ([numerics](./numerics)) through
`--numerics SOLVER[:RECON[:CORRECTION[:SENSOR]]]`, which rewrites the `[numerics]` block of the generated inputs
(defaults: `characteristic`, `6th`, `weno`):

```bash
cd src/tests/flume/verification/mhd
rj2a/check.sh --numerics hlld                     # HLLD, characteristic interpolation
rotor/check.sh --numerics hlld:primitive          # HLLD, primitive interpolation
```

The specs `hlld:primitive` and `hlld:characteristic` assert bounds measured on them; any other spec runs the
scheme-independent checks only (order, bitwise x/y/z and left/right, symmetry, conservation, positivity). Measured
(issue #47, M3-P3c; CPU and FNL agree to the printed digits):

| Problem | Splitting | HLLD primitive | HLLD characteristic |
|---|---|---|---|
| Ryu–Jones 2a, $L_1$ sum at 256 / 512 cells | $3.84 / 2.12\cdot10^{-2}$ | $5.13 / 2.93\cdot10^{-2}$ | $3.76 / 2.12\cdot10^{-2}$ |
| linear waves, CPAW: observed order | 4.96–5.02 | 4.96–5.00 | 4.96–4.99 |
| Brio–Wu, Ryu–Jones 4d: floored cells, HLLD fallbacks | 0 | 0, 0 | 0, 0 |
| field loop, magnetic energy kept at 64 / 128 cells | 0.872 / 0.937 | 0.869 / 0.936 | 0.887 / 0.944 |

On Ryu–Jones 2a HLLD is more accurate than HLL ($5.98$ and $4.32\cdot10^{-2}$). The GLM run differs from the run
without cleaning by round-off ($10^{-12}$; bitwise with the splitting), because the interpolated uniform $B_n$ seeds
$\psi$. The rotor's symmetry defect is round-off amplified by the flow and grows as the dissipation decreases (final
$5\cdot10^{-9}$ to $6\cdot10^{-6}$ depending on the interpolation and the backend), so under `--numerics` it is checked at
step 100 within $10^{-8}$ and at the end within $10^{-4}$.

### EGLM cleaning

The energy-consistent cleaning (`divergence_control = eglm`, [models](./models#divergence-control-eglm)) runs through
the `--divergence-control eglm` option of the glm-pulse, Orszag–Tang and linear-wave scripts, and `rj2a/check.sh`
compares it with GLM (issue #47, M3-P4; CPU and FNL agree to the printed digits):

```bash
cd src/tests/flume/verification/mhd
glm-pulse/check.sh --divergence-control eglm                                  # EV-2
orszag-tang/check.sh --divergence-control eglm [--numerics hlld]              # EV-4, symmetry
linear-wave/check.sh --divergence-control eglm [--numerics hlld]              # MV-5 under EGLM
rj2a/check.sh [--numerics hlld]                                               # EV-3 (EGLM leg)
```

| Check | Splitting | HLLD characteristic |
|---|---|---|
| EV-2, pulse against the telegraph solution (amplitude $10^{-3}$): $L_1$ at 256 cells, d'Alembert / damped | $6.36\cdot10^{-9}$ / $1.38\cdot10^{-9}$ | — |
| EV-3, Ryu–Jones 2a, EGLM against GLM | bitwise, $\psi = 0$ | $6.3\cdot10^{-13}$ (FNL: bitwise) |
| EV-4, Orszag–Tang: $\int\rho$ drift; $\int E$ drift / its source bound | $3\cdot10^{-17}$; $2.8\cdot10^{-3}$ | 0; $7.0\cdot10^{-3}$ |
| Orszag–Tang symmetry at step 100 / at $t = 0.5$ | $5.4\cdot10^{-13}$ / $4.3\cdot10^{-7}$ | $1.2\cdot10^{-12}$ / $1.1\cdot10^{-11}$ |
| linear waves, observed order | 4.96–4.98 | 4.96–4.97 |
| Balsara–Spicer blast ($\beta = 2.5\cdot10^{-4}$, $128^2$, no limiter, no floors) | non-positive pressure at step 3 | reaches $t = 0.01$ (396 steps) |

![Orszag-Tang, drift of the total energy under GLM and EGLM](/flume/eglm-energy.png)

The total energy of the periodic Orszag–Tang box: GLM keeps it to round-off, EGLM changes it through its
nonconservative sources, by $3.0\cdot10^{-4}$ (splitting) and $4.4\cdot10^{-4}$ (HLLD) at $t = 0.5$, within the EV-4
bound.

EV-2 runs at amplitude $10^{-3}$: EGLM couples the $(B_x, \psi)$ pair to the fluid (the damped cleaning energy becomes
heat, the heat drives a flow, the source $-(\nabla\cdot\mathbf{B})\mathbf{u}$ moves $B_x$), which at the GLM amplitude
0.1 floors the damped legs at a relative $L_1$ of $8\cdot10^{-5}$; at $10^{-3}$ EGLM equals GLM to 3 digits. EV-4
asserts $\int\rho$ to round-off and bounds every other drift, at every step, by twice the time integral of its source
magnitude: $\max|\mathbf{B}|\int\|\nabla\cdot\mathbf{B}\|_1$ (momentum), $\max|\mathbf{u}|\int\|\nabla\cdot\mathbf{B}\|_1$
($\mathbf{B}$), $\max|\mathbf{u}\cdot\mathbf{B}|\int\|\nabla\cdot\mathbf{B}\|_1 + \int\|\psi\,\mathbf{u}\cdot\nabla\psi\|_1$
(energy; `mhd/eglm_conservation_oracle.py`): a necessary bound, not a tight one, since the signed sources largely
cancel. The symmetry defect is round-off amplified by the flow, more by the EGLM sources on the splitting scheme, so it
is checked at step 100 within $10^{-8}$ and at the end within $10^{-4}$, as for the rotor. On the blast HLLD with
primitive interpolation also reaches $t = 0.01$; the splitting scheme needs the positivity limiter (M3-P5).

### Positivity limiter

`mhd/positivity/check.sh` runs the limiter (`positivity_limiter = cell`, [numerics](./numerics#positivity-limiter)) with
the floors disabled, on the splitting scheme and on `weno-riemann` (HLLD or HLLC, characteristic; issue #47, M3-P5). A
leg passes when the run reaches its final time and no stage meets an inadmissible backbone; the Euler legs, which have
no floors to stop a run, also assert a positive density and pressure in the last checkpoint.

```bash
cd src/tests/flume/verification/mhd/positivity
./check.sh                       # PV-1, PV-2, PV-3, PV-4 (about 70 min on the CPU, PV-3 takes 55)
./check.sh --legs pv1-3d         # the blast as a sphere on 64^3 (about 1 h per scheme on the CPU)
./check.sh --legs pv5            # the blast across 2:1 AMR seams (about 14 min per scheme on the CPU)
```

| Problem | Without the limiter | Splitting + limiter | Riemann + limiter |
|---|---|---|---|
| PV-1, Balsara–Spicer blast, $\beta = 2.5\cdot10^{-4}$, $128^2$, EGLM | splitting fails at step 3; HLLD runs | passes, up to 1540 faces limited per stage | passes, up to 280 |
| PV-1 as a sphere, $64^3$ (FNL) | not run | passes, up to 51566 | passes, up to 2976 |
| PV-2, Wu–Shu blast, $\beta = 2.51\cdot10^{-6}$, $128^2$, EGLM | splitting fails at step 1, HLLD at step 41 | passes, up to 10333 | passes, up to 352 |
| PV-3, LeBlanc shock tube, 400 cells | both run | passes, never limited | passes, up to 512 |
| PV-3, double rarefaction, 400 cells | splitting runs; HLLC fails at step 79 | passes, never limited | passes, up to 1024 |
| PV-3, planar Sedov blast, 800 cells | splitting runs; HLLC fails at step 69 | passes, never limited | passes, up to 768 |
| PV-4, isentropic vortex ($64^2$, $128^2$) and fast wave (16, 32 cells, EGLM) | — | bitwise equal to the run without the limiter | — |
| PV-5, PV-1 at $64^2$ with $[0.25, 0.75]^2$ refined 2:1 (issue #50) | not run (uniform: splitting fails at step 3) | passes; mass to round-off, seam flux mismatch $3.6\cdot10^{-12}$ | passes; the same |

![Balsara-Spicer blast with the limiter](/flume/blast.png)

![Wu-Shu blast with the limiter](/flume/blast-wushu.png)

Both schemes reach the final time, but not with the same quality: the density of the splitting scheme carries
grid-scale noise around the blast, mild at $\beta = 2.5\cdot10^{-4}$ and plain at $\beta = 2.51\cdot10^{-6}$, where HLLD
stays clean. The limiter guarantees positivity, not accuracy: on strongly magnetised blasts `weno-riemann` with HLLD is
the scheme to use.

![Euler near vacuum with the limiter](/flume/near-vacuum.png)

The near-vacuum problems on the two schemes, with the exact solution of LeBlanc and of the double rarefaction. On
LeBlanc both schemes put the shock ahead of its exact position $x = 7.96$: by 0.33 at 400 cells and by 0.19 at 800
(measured on the splitting scheme), an error that converges slowly under refinement, as reported for this problem in the
literature, and that the limiter does not cause (the splitting scheme never limits here). On the double rarefaction the
velocity is ill-defined where the density nearly vanishes, hence the kink at $x = 0$; on Sedov HLLC reaches a lower
centre density than the splitting scheme ($4\cdot10^{-3}$ against $7\cdot10^{-3}$) and the ambient pressure stays at
its initial $4\cdot10^{-13}$.

PV-2 is the second blast of Wu & Shu (2018, Example 4.4): $p = 10^4$ in the disc, $\mathbf{B} = (1000/\sqrt{4\pi}, 0, 0)$,
$t = 0.001$. PV-3 takes the near-vacuum problems of Zhang & Shu (2010): LeBlanc ($\gamma = 5/3$, $(\rho, e) = (1, 0.1)$
and $(10^{-3}, 10^{-7})$, $t = 6$), the double rarefaction ($\rho = 7$, $u = \mp1$, $p = 0.2$, $t = 0.6$) and the planar
Sedov blast ($p = 4\cdot10^{-13}$ but $2.56\cdot10^8$ in the centre cell, $t = 0.001$); the final states keep
$\rho \ge 9.9\cdot10^{-4}$ and $p \ge 4\cdot10^{-13}$, the ambient value. The Euler splitting scheme needs no limiting on
any of the three at these resolutions: the limiter matters for `weno-riemann`. The face counts include the copies along
the null directions (256 per 1-D face). The counts are those of the CPU; the device ones differ in the splitting
blasts (1468 and 10955), whose limited runs are not bitwise reproducible across backends, and agree elsewhere. The 3-D
blast was run with the relative floor on the device only.

Sedov on HLLC is the case that set the limiter's relative floor. With the absolute floor alone the run stopped at step
10, at every CFL from 0.2 to 0.5: one stage drained the centre cell to $\rho = 10^{-13}$ with its energy kept, and the
next stage, whose time step is that of the step's first state, had no admissible backbone. With the relative floor
($\kappa = 0.1$) the limiter acts in 15 stages and the run completes; the smooth cases of PV-4 stay bitwise. The
splitting scheme takes 614 steps on PV-2 against 284 for HLLD (779 with the absolute floor alone): its limited states
keep low-density cells with large Alfvén speeds, which cost time steps but not admissibility.

PV-5 runs the limited blast across 2:1 AMR seams that the outer shock crosses. Before issue #50 neither scheme
reached the end there: the splitting scheme stopped at step 184, when the end-of-step reflux drove coarse cells beside
the seam to $\rho e < 0$ (the reflux replaced the limited coarse flux with the restricted fine one, a correction no cell
factor bounds), and HLLD stopped at step 261, after the tricubic seam ghost fill had produced inadmissible fine ghosts.
With the ghost positivity blend and the per-stage seam flux synchronisation ([numerics](./numerics#positivity-limiter))
both schemes reach $t = 0.01$ on the CPU and on FNL, the mass stays constant to round-off, and the forest's
$\max|F_\text{coarse} - F_\text{fine}|$ at the end-of-step reflux is $3.6\cdot10^{-12}$: the coarse seam flux is the mean of
the fine ones at every stage. The leg asserts both (mass drift below $10^{-12}$, mismatch below $10^{-9}$).

## Unit tests

| Test | What it pins |
|---|---|
| `test_flume_euler_library` (+ `_fnl`) | Euler eigensystem, flux, Roe average, split consistency; RS(q, q) = f(q) for LLF/HLL/HLLC, HLLC exact on a contact, positive first-order updates; device = host |
| `test_flume_mhd_library` (+ `_fnl`) | MHD, GLM and EGLM eigensystems (including degenerate states), fluxes, auxiliary variables, cyclic invariance; device = host |
| `test_flume_mhd_riemann` (+ `_fnl`) | MHD LLF/HLL/HLLD without cleaning, with GLM and with EGLM: consistency, HLLD exact on contact, tangential and rotational discontinuities, cyclic invariance bitwise, positive updates, EGLM = GLM bitwise at $\psi = 0$; device = host |
| `test_flume_positivity` | PV-0, the limiter on random admissible states with perturbed fluxes (Euler, MHD, EGLM): every limited update positive and above the relative floor, the limiter needed and acting; a NaN or infinite high-order flux replaced by the backbone flux; inadmissible face ghosts blended above the floors, every other value untouched; a cell with a 2:1 seam face (mean donor-state backbone) admissible for any seam factor up to its own, and not without the limiter |
| `test_flume_weno_interpolation` (+ `_fnl`) | WENO interpolation tables: exactness, convergence, device = host; the reconstruction tables unchanged |

Build and run with `fobis build --mode test-flume-<name>-gnu` (`-fnl-nvf --varset local_nvf` for the device twin).

## Regression suite

`src/tests/flume/regression/` holds 26 goldened cases (Sod along x/y/z, AMR, multi-realm, immersed boundary, the MHD
cases RJ2a, Brio–Wu, GLM pulse, Orszag–Tang, rotor, field loop, rotated shock tube, uniform AMR, and the M3 cases: Sod
with HLLC and RJ2a with HLLD on `weno-riemann`, Orszag–Tang with EGLM, the blast with the positivity limiter in 2-D and
in 3-D, the only 3-D flow case of the suite, and, issue #50, the blast with the limiter across two 2:1 seams through
its centre): `run.sh cpu` runs in
CI, `run-fnl-local.sh` on a GPU workstation, and `run-omp-bitwise.sh` checks that the OpenMP build reproduces the serial
one bit for bit.
