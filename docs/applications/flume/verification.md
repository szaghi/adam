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
and the reflux register cross ranks, [#40](https://github.com/szaghi/adam/issues/40)). Leg 4 splits the refined Sod at
its 2:1 face into a coarse and a fine realm glued by `coupling = refined`
([#52](https://github.com/szaghi/adam/issues/52)): the seam ghosts are interpolated (coarse to fine) and restricted (fine
to coarse) with the formulas of the intra-realm AMR seams, and the 2:1 reflux register spans the two realms. The union
is bitwise equal to the single-realm refined run (135168 cells, 216 steps) on 2, 3 and 4 ranks, CPU and FNL, and so is
the same split along $z$, where every 2:1 seam row crosses ranks. The pair declared `mirror` is refused at
initialization, and so is the positivity limiter on any multi-realm run. The MHD twin (`mhd/multirealm/check.sh` leg 4,
RJ2a with GLM on the same cells) is bitwise on all nine fields over 349 steps; in 1-D $\psi$ stays at round-off, so it
tests $\mathbf{B}$, not the cleaning, across the 2:1 seam. Leg 5 ([#51](https://github.com/szaghi/adam/issues/51))
splits Sod, along $x$ and along $z$, into realms whose seam blocks do not line up (one realm on blocks twice as large,
the same cells): bitwise equal to the single-realm run on 2 and 4 ranks. Sod is 1-D, so the MHD rotor (MV-14 leg 5, on an
octree, split both ways) checks where the scattered fluxes land: bitwise on all nine fields, and a deliberately
misplaced overlap drives the density negative at step 2. MV-14 leg 6 ([#54](https://github.com/szaghi/adam/issues/54))
repeats the split on the rotor's own quadtree, blocks lined up and not: bitwise on all nine fields. Before #54 the
tree lookup of the seam peers used the 3-D Morton code on every tree, wrong on a quadtree from level 2 on, and the
lined-up split stopped at initialization; the library unit test `test_tree_closest_block` pins the lookup on binary
trees, quadtrees and octrees, including a single-block tree, whose root the lookup now returns instead of aborting.

### Forward-facing step: a three-realm forest on quadtrees

The Mach 3 wind tunnel with a step of Woodward & Colella (1984, *J. Comput. Phys.* 54, 115–173, §IV b): the tunnel
$[0, 3] \times [0, 1]$ with a step of height 0.2 at $x = 0.6$, $\gamma = 1.4$, the gas entering from the left with
$\rho = 1.4$, $u = 3$, $p = 1$ (sound speed 1, Mach 3), leaving on the right; every other boundary, the step faces
included, is a reflecting wall. The tunnel is L-shaped, so FLUME builds it as **three rectangular realms** glued by two
1:1 mirror seams (`step/make_step.py`), and the step is body-fitted: its two faces are physical walls of the realms, with
no immersed boundary.

```
 y = 1   +--------+-----------------------------+   wall
         |   B    |              C              |
  inflow |        | seam                        | outflow
 y = 0.2 +--------+   +------ wall (step top) --+
         |   A    | wall (step front)
 y = 0   +--------+
       x = 0    x = 0.6                       x = 3
```

| Realm | Extent | Faces |
|---|---|---|
| A | $[0, 0.6] \times [0, 0.2]$ | inflow, wall (bottom), wall (step front), seam with B (top) |
| B | $[0, 0.6] \times [0.2, 1]$ | inflow, seam with A (bottom), seam with C (right), wall (top) |
| C | $[0.6, 3] \times [0.2, 1]$ | seam with B (left), wall (step top), wall (top), outflow |

The realms share the base cell $1/N$. B and C have $N/10$ blocks per axis, narrow enough to refine close to the
seams. A has $N/40$, so that its blocks keep 8 cells along y: a refined block needs at least $2\,n_{gc} = 6$ cells along
every active axis (see below). Both seams are misaligned 1:1 mirror seams
([#51](https://github.com/szaghi/adam/issues/51)). One more level refines a box in each realm, on **quadtrees**
(`ratio = 4`, `nk = 1`, [#46](https://github.com/szaghi/adam/issues/46)).

::: warning The refinement is static
The boxes are placed by hand on the features of the reference solution at $t = 4$, and the grid is frozen after
initialisation: FLUME has no runtime AMR (milestone M5). A solution-driven marker could not place them either, since
the tunnel starts uniform. The boxes therefore follow the flow only once it has settled, and the case verifies the
seam machinery, not adaptivity.
:::

At $N = 80$ the boxes are:
- **A**, $y \le 0.1$: the bow shock standing on the lower wall.
- **B**, $x \in [0.3, 0.525]$, $y \in [0.3, 0.7]$: the bow shock bending up towards the Mach stem. The column
  $x > 0.525$ touches the B–C seam.
- **C**, $x \ge 0.9$: all of the realm except the block column on the seam. That covers the reflected shock reaching
  the step top, the slip line and the re-reflection off the top wall.

A refined block may not touch an inter-realm seam (C's seam is only its face $x = 0.6$; its face $y = 0.2$ is the
step top, a wall). The strongest features, the Mach stem and the triple point at $x \approx 0.6$, therefore sit on the
B–C seam and stay on the base grid. Refining across that seam would take the 2:1 inter-realm coupling
(`coupling = refined`, [#52](https://github.com/szaghi/adam/issues/52)), which is available on octrees only. The case
carries both kinds of seam: two inter-realm mirror seams on different axes, and the intra-realm 2:1 seams of the boxes.

WENO5 characteristic flux splitting, SSP-54, CFL 0.5, reflux on. The corner is left untreated (Woodward and Colella
reset the entropy near it).

**Trees** (the default leg of `step/check.sh`). The same forest at $N = 40$, run to $t = 0.5$ (440 steps) on quadtrees
and on octrees with a null z axis (`nk = 4`, 16 in A: the layout of every AMR verification before #46), compared column
by column with `mhd/quadtree/quadtree_oracle.py` (momentum scaled as one vector). At this resolution A is a single
block that touches its seam and stays unrefined; B is refined on $x \in [0.3, 0.45]$, $y \in [0.4, 0.8]$, C from
$x = 1.2$.

| | Octree, null z | Quadtree (`nk = 1`) |
|---|---|---|
| z spread of a column | $6.1\cdot10^{-13}$ (FNL $3.1\cdot10^{-13}$) | 0 (FNL 0) |
| against the octree, per $(x, y)$ column | — | $6.5\cdot10^{-13}$ (FNL $2.7\cdot10^{-13}$) |
| min $\rho$, min $p$ (CPU and FNL) | 0.469, 0.715 | 0.469, 0.715 |

The trees order the restriction sums differently (8 fine cells over two identical z layers, against 4), so they agree to
round-off, not bit for bit; the octree is z-invariant to round-off only, its tricubic seam weights depending on the z
sub-position. The difference map carries no structure at the seams, inter-realm or 2:1. The leg takes about 10 min on
the CPU, 9 of them for the octree (its blocks split along z too).

::: details What this leg caught: blocks thinner than 2 ngc ([#66](https://github.com/szaghi/adam/issues/66))
A first layout kept A's blocks at 4 cells along y and refined them. The quadtree then differed from the octree by
$2.4\cdot10^{-5}$, with no difference in either tree's z invariance. Bisection put it on the quadtree alone, and only
across ranks:

| Comparison | Max difference |
|---|---|
| Quadtree, np 1, against the octree, np 1 or np 2 | $2\cdot10^{-14}$ |
| Quadtree, np 1, against np 2 | $4.6\cdot10^{-6}$ after one step |

The difference was independent of `seam_ghost_fill` and of `reflux`. At the bisected cell, the coarse ghost density
next to the step wall differed by 0.20.

The coarse ghost layers beside a finer block are restricted from $2\,n_{gc}$ fine cells. A 4-cell block holds fewer,
so the outer layer was restricted from the fine block's own ghost cells, which the same exchange is still filling. They
are fresh or a stage stale depending on the exchange order: the local order happened to refresh them first, the MPI
path did not. The octree passed only because its rank layout kept the two fine blocks together.

The library now refuses blocks thinner than $2\,n_{gc}$ along a refined, non-null axis, the constraint block-structured
AMR codes such as PARAMESH impose. No committed case is affected: every one below 6 cells is so only along a null axis,
where the field is invariant.
:::

![Step forest at t = 0.5 on quadtrees and octrees](/flume/step-trees.png)

**Woodward and Colella's run** (the `full` leg, not default: about 41 min on the CPU, plus 13 min for the uniform
control). $N = 80$, their coarse grid, the boxes at $1/160$, to $t = 4$ (7545 steps): density and pressure stay
positive (min 0.344, 0.370). The picture has the features of their reference:
- the bow shock standing at $x \approx 0.3$ on the lower wall;
- its Mach reflection off the top wall at $x \approx 0.6$, with the triple point near $(0.6, 0.78)$ and the slip line
  leaving it along $y \approx 0.8$;
- the reflected shock striking the step top at $x \approx 1.2$ in a Mach reflection;
- its re-reflection reaching the top wall at $x \approx 2.4$.

The low-Mach layer along the step top is the numerical boundary layer of the untreated corner, the artifact Woodward
and Colella remove with their entropy fix.

The zoom compares the refined run with the same run on a uniform $1/80$ grid (no box). The Mach stem and the triple
point sit on the B–C seam, on the base grid in both runs, and agree. C's 2:1 seam at $x = 0.9$ is crossed by the
reflected shock and the slip line without a kink. Downstream of it, the box resolves the shock pattern near the step top
($x \approx 1.2$) more sharply than the uniform grid, which smears it into a single blob.

![Woodward-Colella step at t = 4](/flume/step.png)

`step/check.sh` (EV-step: `--legs trees`, the default, and `full`; `--keep` keeps the checkpoints for
`make_doc_figures.py --only step step-trees`). Input: `step/make_step.py <dir> --tree quad|oct --cells N --refine
--time-max T`, which writes the manifest `step.ini` and the realms `step-A.ini`, `step-B.ini`, `step-C.ini`.

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

### Quadtree AMR

A quadtree (`ratio = 4`) refines x and y only, so its 2:1 seams are 1:1 in z. Until
[#46](https://github.com/szaghi/adam/issues/46) the seam machinery assumed an octree, and the coarse cells beside a
quadtree seam picked up a spurious z dependence ($10^{-2}$ at the first step, negative pressure in Orszag–Tang at step
26), so quadtree AMR was refused and every AMR verification ran on an octree with a null z axis and `nk = 4`. MV-15
(`mhd/quadtree/check.sh`, `quadtree_oracle.py`) runs the same 2-D problem on that octree, on a quadtree with `nk = 1`
and on a quadtree with an active z axis (`nk = 4`), keys the cells by $(x, y)$ and compares them per variable:

| Leg | Quadtree `nk = 1` against the octree | Octree z spread | Quadtree `nk = 4` z spread |
|---|---|---|---|
| Orszag–Tang, $32^2$ + $[0.25, 0.75]^2$ refined, $t = 0.2$ | $4.1\cdot10^{-12}$ ($\psi$; $\le 3\cdot10^{-13}$ elsewhere) | $5.7\cdot10^{-12}$ | 0 |
| Balsara–Spicer blast, limiter on, $t = 0.006$ | $1.9\cdot10^{-5}$ ($\psi$), $3.7\cdot10^{-6}$ ($\rho$) | $2.8\cdot10^{-5}$, $9.3\cdot10^{-6}$ | — |

The runs agree to round-off, not bit for bit: an octree restriction averages 8 fine cells (two identical z layers), a
quadtree 4, and even the octree's identical layers differ at round-off, its tricubic weights depending on the z
sub-position. In the blast the limiter's switches amplify that round-off: the octree's own layers drift apart by
$2.8\cdot10^{-5}$, and the quadtree differs from the octree by no more than that spread (ratio per variable 0.40 to
0.80 on the CPU, 0.87 to 1.02 on FNL; the same 579 limited stages on both); the leg bounds the agreement by twice the
reference spread. FNL gives the same picture on Orszag–Tang ($6.4\cdot10^{-12}$ against the octree, z spread 0 with
`nk = 4`). The `nk = 4` quadtree computes
the z fluxes and keeps the solution z-invariant bit for bit; its z cells enter the time step, so it is compared with
itself only. The unit test `test_quadtree_seam_ghost` pins the seam exchange on a linear field (octree, quadtree
`nk = 4` and `nk = 1`, 1 to 3 ranks; its FNL build also asserts the device exchange equal to the CPU one, bit for bit).

## Ghost cells

### GP: every face and edge ghost on a linear field

**Why.** The directional WENO stencils read only the face ghosts: the cells outside a block along one axis. The
dissipative terms of M4 ([#65](https://github.com/szaghi/adam/issues/65)) take cross derivatives, such as the
tangential derivatives of the viscous stress or of the current at a face. Those read the edge ghosts too: the cells
outside a block along two axes, in the rows next to a block edge. Before M4 nothing checked what an edge ghost held.
GP checks every one, on every path that fills it: the same-level copy and the MPI exchange, the 2:1 restriction and
tricubic fill, the inter-realm seams (1:1 mirror, misaligned, 2:1 refined), and the physical boundary conditions.

**Oracle.** Each case starts from the [`linear`](./initial-conditions#linear-a-linear-field-verification) initial
condition, $\mathbf q = \mathbf q_1 (1 + \mathbf g \cdot \mathbf x)$, with every velocity component non-zero. The
copy, the restriction and the tricubic fill reproduce such a field exactly. Each case takes one step with
`CFL = 1e-30`, which moves the field by ~$10^{-30}$ relative. The post-step write then fills the ghosts in the order
of a Runge–Kutta stage (seams, then the intra-realm exchange and the boundary conditions) and writes them.
`ghosts/ghost_probe.py` compares every ghost with the value a stencil must read:
- the field at the ghost's own centre, inside the realm or across a seam;
- along an axis where the ghost lies beyond a physical face, that face's condition applied to the value of its
  donor: the coordinate mirrored with the normal momentum negated for an inviscid wall (and the normal field on MHD);
  mirrored with the velocity reflected about the wall velocity for a no-slip wall (and the temperature $T_w^2/T$ for
  an isothermal one); the first interior cell for extrapolation; the inflow state for inflow; the wrapped coordinate
  for periodic. Beyond two faces the conditions compose in the backends' order (inflow first, else the first face in
  axis order acts on a donor valued the same way).

The no-slip walls are not linear in the state; the probe evaluates them with an independent implementation of the rule
of [Boundary conditions](./boundary-conditions#wall-noslip-and-wall-isothermal-no-slip-walls), so `walls3d` checks the
plumbing (which donor, which face, which wall parameters, in which order), not the physics of the wall. The physics
comes with the viscous fluxes (VV-3, Couette).

Faces and edges must hold that value to $10^{-12}$ of the field scale. Corners are reported, not asserted, since no M4
stencil reads them, and neither are the ghosts that lie outside every realm (`solid`, inside the step).

| Case (`ghosts/make_probe.py`) | Paths |
|---|---|
| `box3d` (Euler and MHD) | octree, the block at the origin refined: 2:1 seams meeting an inflow face and two walls in 3-D |
| `channel2d` | quadtree, periodic x and a wall in y, the 2:1 seam crossing the periodic boundary at the wall |
| `mirror3d` | two realms and a 1:1 mirror seam, walls on the other faces: seam edges in every plane |
| `refined3d` | as `mirror3d`, the second realm one level finer: a `coupling = refined` (2:1) seam |
| `walls3d` (Euler and MHD with EGLM) | the no-slip walls of #65 P1 around a refined corner: a moving `wall-noslip`, a resting one, a `wall-isothermal` with a wall velocity, an inviscid wall, inflow and extrapolation |
| `step` | the [Woodward–Colella forest](#forward-facing-step-a-three-realm-forest-on-quadtrees) at $N = 80$ with its boxes: misaligned mirror seams, a re-entrant corner, inflow, outflow and walls |

**The poison.** The negligible step leaves a stale ghost equal to a fresh one, so the oracle alone cannot see a ghost
read before its donor is written. Every case therefore runs with `[diagnostics] ghost_poison = .true.`: each ghost is
set to NaN before the fill, and a ghost that reads an unfilled donor is written as NaN, which the probe reports as an
infinite error. Its negative control is the defect it found (below).

**Results.** Every case runs on 1, 2 and 3 ranks on the CPU, and on 1 and 2 ranks on FNL. Every face and edge group
holds the linear field to round-off: at most $1.3\cdot10^{-15}$ of the field scale on the CPU, the same on FNL.

**What it found.** Before P0, a ghost beyond two realm faces held a copy of its inward diagonal cell. That covers wall
plus wall at the step corner, wall plus seam where the step face meets the A–B seam, inflow plus wall, and the like.
On the step forest those ghosts were off by up to 0.78 of the field scale: the wrong position, and the wrong sign of
the normal momentum. They now compose the conditions of the faces, as
[Boundary conditions](./boundary-conditions#edges-and-corners-fec-6) describes. The face ghosts, and the edges with
another block beyond one of their faces, were already right.

**What the poison found** (P1, after P0 had been committed). An extrapolation ghost copied the previous ghost of its
chain along the normal. When the ghost also lies beyond a block interface along another axis, deeper than along the
normal, that previous ghost is in the same crown, so in the same device launch: on FNL the rows raced. Two identical FNL
runs of `sod-x` wrote different edge and corner ghosts (the interior identical), and with the poison GP failed the
`block+extrapolation` edges on every case (NaN read). The CPU walks each row outward and read its donors fresh, by
the row order alone. The P0 probe could not see it: without the poison the stale value equalled the fresh one. An
extrapolation ghost now copies the first interior cell along its normal directly (the same value, and a donor the
exchange has filled), on both backends; with the poison every case passes and two FNL runs are bit for bit equal.

![Ghost probe on the step forest, before and after #65 P0](/flume/ghosts.png)

Each point is a face or edge ghost of the step forest, coloured by its error. Top: before P0, the realm edges at the
step corner, the domain corners and the inflow corners carry $O(10^{-1})$ errors; every other ghost is at round-off.
Bottom: after.

The probe also caught one gap in its own instrument: the step-0 file is written while each realm initialises, before
the forest connects the realms, so its seam ghosts are unfilled. GP therefore reads the file written after the
negligible step, whose write refills the seams first.

```bash
cd src/tests/flume/verification/ghosts && ./check.sh                    # CPU, np 1 2 3, every case
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh --np "1 2"  # FNL
```

## Dissipative terms

### DC: the input contract

**Why.** M4 ([#65](https://github.com/szaghi/adam/issues/65)) reads each dissipative coefficient either as a
coefficient or as its dimensionless number (#49 N4), refuses the positivity limiter with any of them (D-M4-5) and adds
the no-slip walls (D-M4-7). An invalid combination must stop the run with a message that names it; a valid one must
imply exactly the coefficient it states; nothing may be dropped silently. The keys are listed in
[input](./input#dissipative-terms-issue-65-m4).

**Cases** (`dissipation/make_contract.py`, on the regression inputs `sod-x` for Euler and `uniform-amr-mhd` for MHD):

| Kind | Cases | Oracle |
|---|---|---|
| refuse | both keys of a term (3), `prandtl` without a viscosity, resistivity on Euler, a negative viscosity, a zero Reynolds number, the power law without a viscosity, its keys without the law, a missing reference temperature, an unknown law, `dissipative_order = 3`, the limiter with a viscosity, a normal wall velocity, an isothermal wall without (or with a negative) temperature, `dissipative_order = 4` with `ngc = 2`, a negative power-law exponent | the run stops with the expected message |
| log | `reynolds` + `prandtl`; `viscosity` + `conductivity`; `magnetic_reynolds` + `reynolds` (MHD) | the logged $\mu$, $k$, $\eta$ equal $1/Re$, $\mu c_p/Pr$, $1/Rm$ exactly; every run takes ten steps and logs its diffusive time-step limit |
| ideal | zero coefficients, `dissipative_order = 2` | the run equals the base run bit for bit |
| convert | the Euler case with viscosity, conductivity, the power law and an isothermal moving wall; the MHD case with resistivity; dimensionalised with $L_0 = 2$, $u_0 = 1/2$, $\rho_0 = 4$ | every dimensional key logged converted back exactly (`scaling.py check-log`), the conductivity and the temperatures with the gas constant |
| convert-refuse | `lundquist` under a `[reference]` without the Alfvénic preset | refused |

**Results.** All 25 cases pass, on the CPU (2 ranks) and on FNL; the conversion check covers 34 keys on the Euler
case and 30 on the MHD one.

```bash
cd src/tests/flume/verification/dissipation && ./check.sh
```

### VV-1 and VV-2: viscous waves

**Why.** A dissipative flux can be consistent and still wrong in its order (a missing correction term), its coefficient
(a missing 4/3, a conductivity without $R$) or its direction (a cross derivative read with the wrong stencil). A small
sine wave on a gas at rest is one Fourier mode of the linearised compressible Navier–Stokes equations; the oracle
(`viscous/waves.py`) integrates that mode exactly, $\exp(Mt)$ of the $3\times3$ system for $(\rho', u_n, T')$ plus the
decoupled transverse velocity, from the Fourier coefficients of the initial state, so it holds every mode the start
excites, not only a decay rate. The error is the RMS over the cells of the point values.

**Cases** ([`sine-wave`](./initial-conditions#sine-wave-a-small-sine-wave-on-a-uniform-state-verification), $\gamma = 1.4$,
$\rho_0 = p_0 = 1$, $\mu = 0.01$, SSP-RK(5,4), CFL 0.2, WENO5 characteristic):

- VV-1, shear wave $A = 10^{-4}$, $t = 0.5$: along x and along y at `dissipative_order` 4, along x at order 2, and along
  the diagonal (both directions active, so the cross terms $\partial_i u_d$ and the tangential derivatives are
  exercised; $t = 0.1$). The velocity axis is kept active (4 blocks of 6 cells across it): FLUME freezes the momentum
  normal to a null direction, which would freeze the shear.
- VV-2, the right-running acoustic wave with $k = 0.02$ as well, so $\mu$ and $k$ act together. A sound wave steepens at
  $O(A^2)$: at $A = 10^{-4}$ that residue ($10^{-8}$) flattened the error at every resolution. The leg therefore runs
  the twin with amplitude $-A$ and judges the odd part $(f(+A) - f(-A))/2$, whose residue is $O(A^3)$, at
  $A = 10^{-5}$.

**Oracle.** The observed order of the finest pair $\ge 3.8$ at order 4, $\ge 1.9$ at order 2.

**Results** (CPU, 2 ranks):

| Leg | N | RMS error | Order |
|---|---|---|---|
| shear x, order 4 | 32 / 64 / 128 | $3.64\cdot10^{-10}$ / $2.29\cdot10^{-11}$ / $1.43\cdot10^{-12}$ | +3.99, +4.00 |
| shear y, order 4 | 32 / 64 / 128 | identical to x, bit for bit | +3.99, +4.00 |
| shear x, order 2 | 32 / 64 / 128 | $3.68\cdot10^{-8}$ / $9.20\cdot10^{-9}$ / $2.30\cdot10^{-9}$ | +2.00, +2.00 |
| shear diagonal, order 4 | 24 / 48 / 96 | $6.72\cdot10^{-10}$ / $1.16\cdot10^{-11}$ / $2.43\cdot10^{-13}$ | +5.86, +5.57 |
| acoustic, order 4 | 32 / 64 / 128 | $8.20\cdot10^{-11}$ / $1.69\cdot10^{-12}$ / $2.18\cdot10^{-14}$ | +5.60, +6.28 |

The measured decay rate of the shear wave is $\nu k^2$ to $2\cdot10^{-8}$ relative at 128 cells. The diagonal wave
converges faster than its design order: on the diagonal the leading truncation errors of the x and y fluxes appear to
cancel. On the acoustic wave the error is dominated by the inviscid WENO part, not by the viscous one, so VV-2 checks
that $\mu$ and $k$ enter with the right coefficients (the measured amplitude equals the exact one to $6\cdot10^{-10}$
relative), and VV-1 pins the order of the dissipative flux. The measured acoustic rate, $0.29576$, differs by 4% from the
first-order Kirchhoff–Stokes rate $k^2/(2\rho_0)\,(\tfrac43\mu + (\gamma-1)k/c_p) = 0.30831$: here
$\nu k/a = 0.053$, and the exact linear oracle carries that correction.

### VV-3: compressible Couette flow

**Why.** The first test of the no-slip walls with the viscous fluxes, and of a nonlinear energy flux (the viscous
heating $u\,\tau$).

**Case** ([`couette`](./initial-conditions#couette-exact-compressible-couette-flow-verification), `viscous/profiles.py`):
a fixed `wall-isothermal` at $y = 0$ ($T_w = 1$), a `wall-noslip` moving at $U = 1$ at $y = 1$, $\mu = 0.01$,
$Pr = 0.75$, $p = 1$. The exact profile is an exact steady solution, so the run starts from it and the departure at
$t = 0.2$ is the error of the scheme, with no relaxation transient.

**Oracle.** Order $\ge 1.9$ on $u/U$ and $T/T_w$: the moving-wall mirror is exact for the linear $u$ and the symmetric
$T$, the isothermal geometric mirror $T_w^2/T$ is second-order, so the run is second-order overall.

**Results** (CPU and FNL, 2 ranks): RMS of $u/U$ $1.71\cdot10^{-7}$ / $3.84\cdot10^{-8}$ / $8.72\cdot10^{-9}$, of $T/T_w$
$3.01\cdot10^{-6}$ / $6.62\cdot10^{-7}$ / $1.49\cdot10^{-7}$ at 32 / 64 / 128 cells: orders +2.14 and +2.15. (With the
linear isothermal mirror of P1 the errors were 18% lower, at the same orders: the geometric mirror, adopted because the
linear one made negative ghost densities beside hot gas, differs from it by $O(\Delta x^2)$.)

The fields show where that second order comes from. Velocity and temperature are uniform along the periodic $x$, and
the profiles lie on the exact ones (the linear $u$, the viscous-heating parabola of $T$). The error enters at the
isothermal wall, whose geometric mirror is second order, and diffuses into the channel: at $t = 0.2$ its front has
reached about $y = 0.5$, and beyond it the solution is exact to round-off at 128 cells. The moving adiabatic wall adds
nothing, its mirror being exact for the linear $u$ and the symmetric $T$.

![VV-3: the Couette channel, its fields, the profiles against the exact ones and where the error lives](/flume/couette.png)

### VV-4: Becker's viscous shock

**Why.** The strongest nonlinear viscous flux of the suite: a shock resolved by the viscosity, with the energy flux
carrying both the viscous work and the heat flux.

**Case** ([`becker-shock`](./initial-conditions#becker-shock-becker-s-exact-viscous-shock-verification)): Becker's
exact profile in the shock frame (constant $\mu$, $Pr = 3/4$), supersonic `inflow` upstream, `extrapolation`
downstream, $x \in [-0.25, 0.25]$, $t = 0.05$; Mach 2 with $\mu = 0.005$ and Mach 3 with $\mu = 0.01$, which put about 8
cells across the shock at 256 cells (at Mach 3 with $\mu = 0.005$ the shock spans 4 cells, WENO's nonlinear weights still
act, and the measured order was +2.25). Like Couette, the run starts from the exact solution.

**Oracle.** Order $\ge 3$ on the RMS and the maximum of $(u - u_\mathrm{exact})/(u_1 - u_2)$.

**Results** (CPU, 2 ranks):

| Mach | 64 | 128 | 256 | Orders (RMS) |
|---|---|---|---|---|
| 2 | RMS $1.05\cdot10^{-2}$, max $6.2\cdot10^{-2}$ | $1.92\cdot10^{-3}$, $1.15\cdot10^{-2}$ | $1.30\cdot10^{-4}$, $8.3\cdot10^{-4}$ | +2.45, +3.88 |
| 3 | RMS $1.47\cdot10^{-2}$, max $8.9\cdot10^{-2}$ | $3.09\cdot10^{-3}$, $1.89\cdot10^{-2}$ | $2.77\cdot10^{-4}$, $1.7\cdot10^{-3}$ | +2.25, +3.48 |

**FNL.** Every leg passes on FNL (2 ranks) with the same numbers to the printed digits, except where the round-off of
the run itself shows: the acoustic odd part, whose $\pm A$ difference amplifies the round-off of a $10^{-5}$ wave
($8.196853\cdot10^{-11}$ against $8.196891\cdot10^{-11}$ at 32 cells), and the last digit of the Couette $u$ error at
128 cells.

### VV-5: Ohmic decay and the visco-resistive Alfvén wave

**Why.** The Ohmic flux has its own kernels, a curl form whose cross terms an axis-aligned wave never reads, and a sign
that the plan in #65 printed wrong (an anti-diffusion; see [numerics](./numerics#dissipative-fluxes-navier-stokes)). A sign error
grows the field, a missing cross term changes the diagonal rate, a coefficient error changes every rate: each is a
number the exact linear mode pins.

**Cases** (`sine-wave` in its `magnetic` mode, a transverse field $b_t = A\sin(\mathbf k\cdot\mathbf x)$ on a gas at
rest, $A = 10^{-5}$, $\eta = 0.01$, the same numerics as VV-1; oracle `viscous/waves.py`, the $2\times2$ system for
$(v_t, b_t)$ integrated exactly):

- Ohmic decay along x at `dissipative_order` 4 and 2 (`mhd-none`, no background field): the rate is $\eta k^2$.
- Ohmic decay along the diagonal, $t = 0.1$ (`mhd-none`): the cross terms $\partial_i B_d$ of the curl form.
- The visco-resistive Alfvén wave along x, $B_n = 1$, $\mu = 0.01$, $\eta = 0.005$ (`mhd-glm`): $\mu$ and $\eta$
  together, damping $(\nu + \eta)k^2/2$.
- VV-2's acoustic wave on `mhd-none` with $\mathbf B = 0$, $\mu = 0.01$, $k = 0.02$, odd part of the $\pm A$ twins: the
  viscous and heat fluxes on MHD, which read the temperature from the MHD auxiliary state, so a wrong slot or a wrong
  $T$ shows here and nowhere else (the constant laws of the Alfvén leg never read $T$).

**Oracle.** As VV-1: the observed order of the finest pair $\ge 3.8$ at order 4, $\ge 1.9$ at order 2.

**Results** (CPU and FNL, 2 ranks, identical to the printed digits except as noted):

| Leg | N | RMS error | Order |
|---|---|---|---|
| Ohmic x, order 4 | 32 / 64 / 128 | $3.64\cdot10^{-11}$ / $2.29\cdot10^{-12}$ / $1.43\cdot10^{-13}$ | +3.99, +4.00 |
| Ohmic x, order 2 | 32 / 64 / 128 | $3.68\cdot10^{-9}$ / $9.20\cdot10^{-10}$ / $2.30\cdot10^{-10}$ | +2.00, +2.00 |
| Ohmic diagonal, order 4 | 24 / 48 / 96 | $1.14\cdot10^{-10}$ / $7.19\cdot10^{-12}$ / $4.50\cdot10^{-13}$ | +3.99, +4.00 |
| Alfvén x, order 4 | 32 / 64 / 128 | $6.37\cdot10^{-11}$ / $1.10\cdot10^{-12}$ / $2.18\cdot10^{-14}$ | +5.85, +5.66 |
| acoustic on MHD, order 4 | 32 / 64 / 128 | $8.20\cdot10^{-11}$ / $1.69\cdot10^{-12}$ / $2.19\cdot10^{-14}$ | +5.60, +6.27 |

The measured Ohmic rate is $0.3947841$ against $\eta k^2 = 0.3947842$: the sign is the derived one. Along x the Ohmic
decay of $b_y$ is the same equation as the shear wave's $u_y$ in VV-1, and its errors are the VV-1 errors scaled by the
amplitude, to three digits: two kernels, one stencil. The Alfvén wave decays at $0.29610$ against
$(\nu + \eta)k^2/2 = 0.29609$; as on the acoustic wave, its error is dominated by the inviscid WENO part, hence the order
above four. The acoustic wave on MHD reproduces VV-2 on Euler: the measured amplitude agrees to ten digits at 32 cells
($8.6252132315\cdot10^{-6}$ against $8.6252132321\cdot10^{-6}$); the RMS errors agree to three, the inviscid parts of
the two models being different kernels. FNL differs from the CPU only in the round-off of that leg's odd part
($2.21\cdot10^{-14}$ against $2.19\cdot10^{-14}$ at 128 cells, +6.26).

**The diagonal runs without divergence control.** With GLM or EGLM the same ladder measured +7.65 then +1.72, and +3.32
from 96 to 192 cells. The cause is the ideal scheme: run with $\eta = 0$, where the exact field is static, GLM's
upwinding ($c_h = 3$) damps it at fifth order ($8.5\cdot10^{-12}$ at 48 cells, $2.7\cdot10^{-13}$ at 96), with the
opposite sign to the Ohmic error and as large, so the two cancel near 48 cells and a 24/48/96 ladder is not
asymptotic. Without divergence control the $\eta = 0$ field stays at round-off. A smaller time step (CFL 0.05) changed
no printed digit, so neither error is temporal. GLM and EGLM with resistivity are covered by VV-6.

### VV-6: Ohmic heating and the divergence constraint

**Why.** The Ohmic energy flux $G^E_d = \sum_i B_i\,G^{B_i}_d$ must return to the gas, as heat, exactly the energy the
field loses: a wrong energy flux leaves the field equations right and shows only in this budget. And the curl form
must not feed $\nabla\cdot\mathbf B$, the constraint GLM and EGLM are there to keep.

**Case** (`waves.py budget`): the diagonal Ohmic wave, 48 cells, $A = 10^{-4}$ (the budget compares energies of
$O(A^2)$ with differences of the internal energy, so it needs the signal), $\eta = 0.01$, $t = 0.2$, with `mhd-glm` and
with `mhd-eglm`, each against its ideal twin ($\eta = 0$).

**Oracle.** Heat gained over field and kinetic energy lost within $10^{-3}$ of 1; the drift of the total energy at most
$10^{-13}$ of the total; $\max|\nabla\cdot\mathbf B|$ (central differences) at most 1.05 times the ideal twin's plus
$10^{-14}$.

**Results.** Both models: the field loses $6.770374\cdot10^{-10}$ of its $2.5\cdot10^{-9}$, the gas gains
$6.770375\cdot10^{-10}$, a ratio of 1.0000002 on the CPU and 1.0000015 on FNL. Those departures are round-off: the
internal energy is the total minus the kinetic and magnetic parts, and a difference of $6.8\cdot10^{-10}$ out of a total
of 2.5 keeps six digits. The total energy drifts by $1.8\cdot10^{-16}$ of itself. $\max|\nabla\cdot\mathbf B|$ is
$2.1\cdot10^{-18}$ with resistivity and $2.6\cdot10^{-18}$ without (GLM; EGLM alike).

**What the divergence leg does not show.** On a uniform periodic grid the diagonal field is discretely solenoidal from
the start, and the curl-form flux keeps it so; both runs sit at round-off, so the leg checks only that the resistive
flux adds no divergence there. The test with teeth is the 2:1 AMR seam, where the discrete operators no longer commute
(PRISM's [#29](https://github.com/szaghi/adam/issues/29) floor): it belongs to M4 P4.

![VV-1 to VV-5: convergence, and Becker's shock against the exact profile](/flume/viscous.png)

```bash
cd src/tests/flume/verification/viscous && ./check.sh                  # CPU, every leg
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh             # FNL
```

### VV-7: conservation across 2:1 seams and forest seams

**Why.** The dissipative fluxes are added to the inviscid ones before the seam accumulation, so the Berger–Colella
reflux and the inter-realm register must carry them: at a 2:1 face the coarse and the fine side compute different
viscous, heat and Ohmic fluxes, and without the register the grid leaks. The stencils also reach further than WENO's:
the tangential derivatives read the edge and corner ghosts, which no inviscid stencil touches.

**Cases.** A diagonal wave of amplitude $0.01$ (large, so the seams carry a real flux) in the periodic box with its
centre $[0.25, 0.75]^2$ refined 2:1 (a quadtree, 8 coarse-fine faces), $t = 0.05$, 48 cells per side at the base
level: the acoustic wave on Euler ($\mu = 0.01$, $k = 0.02$: density, normal velocity and temperature all vary) and the
Alfvén wave on MHD-GLM ($B_n = 1$, $\mu = 0.01$, $\eta = 0.01$). Each runs three times: with reflux, without, and
without reflux and without the dissipative coefficients (the ideal twin). An octree with null $z$ gave the same drifts to
the printed digits at four times the cells, so the leg keeps the quadtree; the refined forest below runs on an octree.

**Oracle** (`waves.py conserve`). The drift of each conserved integral is divided by $\int|q|\,dV$, the scale of its
round-off (the momentum of a sine wave integrates to zero, so its own value is no scale; $\psi$, sourced only by the
discrete $\nabla\cdot\mathbf B$, is measured against $c_h\int|\mathbf B|\,dV$):

- with reflux, every drift $\le 10^{-13}$;
- without reflux, the largest drift $\ge 10^{-10}$ (the seams are exercised);
- the leak without reflux differs from the ideal twin's by $\ge 10^{-10}$: the dissipative fluxes cross the seams
  unmatched, and the first item shows the register removes that too.

**Results** (CPU and FNL, 2 ranks, identical to the printed digits except the round-off of the conserved runs, at most
$1.1\cdot10^{-16}$ on FNL):

| Case | With reflux | Without | Dissipative minus ideal leak |
|---|---|---|---|
| Euler | $\le 1.1\cdot10^{-17}$ | $3.0\cdot10^{-5}$ | $2.0\cdot10^{-6}$ |
| MHD-GLM | $\le 1.5\cdot10^{-16}$ | $5.9\cdot10^{-5}$ | $5.1\cdot10^{-6}$ |

**Forests.** The same waves with outflow ($x$ extrapolation) faces, split at $x = 0.5$ into two realms (`waves.py
split`): `mirror` (1:1, quadtree) and `refined` (the single run refines $x > 0.5$; the coarse and the fine realm are
glued 2:1, on an octree, because the `refined` coupling needs ratio 2 along every axis, null $z$ included). The union
must reproduce the single realm: fields within $10^{-10}$ relative. Not bitwise, unlike MV-14's piecewise-constant
states: each realm evaluates the sine from its own block origins, one ulp apart at step 0 ($3\cdot10^{-15}$); the
largest difference after the run is $8.7\cdot10^{-12}$ on the CPU and $4.4\cdot10^{-12}$ on FNL, at the outer $x$
faces, not at the seam.

This leg found a library defect. A seam slab's edge and corner ghosts past a tangential boundary of the realm lie
outside the peer domain, and were left to the boundary condition; along a periodic axis no condition fills them, and
they kept their initial values. WENO never reads them, the dissipative tangential derivatives do: the mirror split went
NaN at the first step (on a single realm of the same half domain the run is clean). The seam enumeration now wraps
those ghost centres by one period when the axis is periodic on both realms
([forest guide](/guide/forest#seams-across-ranks-issue-40)).

### VV-8: accuracy across 2:1 seams

**Why.** Conservation says nothing about accuracy: a seam can conserve and still degrade the order.

**Cases.** The diagonal shear wave ($\mu = 0.01$, Euler) and the diagonal Ohmic wave ($\eta = 0.01$, `mhd-none`) of
VV-1 and VV-5, on the quadtree with its centre refined 2:1, at 24 / 48 / 96 cells at the base level (48 / 96 / 192 in
the refined box), against the exact linear solution, the RMS weighted by the cell area.

**Oracle.** Order $\ge 1.8$ on the finest pair; the Ohmic ladder without reflux, the negative control, $\le 1.3$.

**Results** (CPU and FNL identical, 2 ranks; N is the finest cell count; the rows without coefficients are probes,
not legs):

| Leg | N | RMS error | Order |
|---|---|---|---|
| shear, $\mu$ | 48 / 96 / 192 | $4.26\cdot10^{-8}$ / $1.09\cdot10^{-8}$ / $2.78\cdot10^{-9}$ | +1.96, +1.98 |
| shear, $\mu = 0$ | 48 / 96 / 192 | $5.39\cdot10^{-8}$ / $1.36\cdot10^{-8}$ / $3.42\cdot10^{-9}$ | +1.99, +1.99 |
| Ohmic, $\eta$ | 48 / 96 / 192 | $1.15\cdot10^{-9}$ / $2.94\cdot10^{-10}$ / $8.12\cdot10^{-11}$ | +1.97, +1.86 |
| Ohmic, $\eta = 0$ | 48 / 96 | $7.9\cdot10^{-15}$ / $2.6\cdot10^{-15}$ | round-off |
| Ohmic, no reflux | 48 / 96 / 192 | $2.37\cdot10^{-9}$ / $1.18\cdot10^{-9}$ / $5.82\cdot10^{-10}$ | +1.00, +1.02 |

The composite grid is second order, against fourth on a uniform grid. The rows without coefficients say where that
comes from. The shear wave is second order without viscosity too: its seam error is the inviscid one. With the mean
restriction (a coarse ghost is the mean of the fine point values under it) and the Berger–Colella reflux of
point-value fluxes (a conservative but O(H²)-inconsistent coarse flux), the composite error is second order
([issue #68](https://github.com/szaghi/adam/issues/68) F1: each of the two caps the order at 2 on its own; the #21 bound
on restriction-compatible prolongations does not apply to the tricubic fill). A static field has no inviscid flux at
first order, so the Ohmic composite at $\eta = 0$ stays at round-off, and the Ohmic composite error is the dissipative
flux crossing the seam: second order with reflux, first without (the gradient of mean-restricted ghosts is an O(1)
error, which the reflux replaces by the fine flux). The plan in #65 asked for the seam error "within a stated factor of the
uniform fine run"; between a second-order seam and a fourth-order interior that factor grows with $N$, so the leg
asserts the order instead.

The error maps show where the seam error lives. On the uniform grid the error of the shear wave is smooth and spread
over the domain, near $10^{-13}$. On the refined grid the Ohmic error is a ring on the 2:1 faces, two orders above the
rest of the domain: the dissipative flux at the face. The shear error fills the refined box and its surroundings
instead, being the inviscid seam error that the flow carries and the viscosity spreads.

![VV-1 and VV-8 in 2-D: the shear wave and its error on a uniform grid; the shear and Ohmic errors with the centre refined 2:1](/flume/viscous-fields.png)

### CPU and FNL agree run by run

**Why.** The legs above pass on both backends, but a pass bounds an error norm, and two backends can pass with
different fields. Every dissipative kernel exists twice, as a CPU loop and an FNL device kernel (OpenACC, or OpenMP
offload), so the same run on both must give the same state up to the round-off of the operation order.

**Check** (`viscous/check.sh --leg agree`, after every leg has run on both executables; `waves.py agree`): each run on
the CPU against the same run on FNL, cell by cell on the last checkpoint, blocks matched by their origin. Each
difference is divided by the physical scale of its variable, not by its own maximum: density and energy by their
maximum, the momentum by $\rho_0 a_0$ (an Ohmic wave at rest carries momentum of $O(A^2)$, far below the round-off of
the state), the field by its maximum, $\psi$ by $c_h \max|\mathbf B|$. Bound: $10^{-10}$.

**Results** (one clean build of each backend, 2 ranks; 72 runs, every resolution of every ladder, the $-A$ twins, the
ideal twins and both sides of each forest):

| Leg | Runs | Largest difference |
|---|---|---|
| VV-1 shear | 12 | $2.7\cdot10^{-14}$ |
| VV-2 acoustic | 6 | $4.9\cdot10^{-15}$ |
| VV-3 Couette | 3 | $1.4\cdot10^{-13}$ |
| VV-4 Becker | 6 | $4.3\cdot10^{-13}$ |
| VV-5 Ohmic, Alfvén, acoustic on MHD | 18 | $4.7\cdot10^{-14}$ |
| VV-6 Ohmic budget | 4 | $5.3\cdot10^{-15}$ |
| VV-7 conservation | 6 | $9.3\cdot10^{-15}$ |
| VV-7 forests | 8 | $3.4\cdot10^{-14}$ |
| VV-8 seams | 9 | $8.4\cdot10^{-14}$ |

The largest differences sit in the Becker shock, the steepest gradients of the suite, at $4\cdot10^{-13}$ of the
momentum scale: no run differs by more than a few hundred units of round-off. The DC contract passes on both backends as
well.

```bash
cd src/tests/flume/verification/viscous
./check.sh && FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh && ./check.sh --leg agree
```

## Runtime regridding (M5)

Milestone M5 ([#74](https://github.com/szaghi/adam/issues/74)) regrids during the time loop. `[amr] frequency = 0`, the
value of every input so far, keeps the AMR of the initial condition only; `n > 0` regrids every `n` steps. The library already holds a complete regrid step (`adam_object%amr_update`: the tree adapts with
2:1 balance, the interior data is prolonged or restricted, the blocks are redistributed, the maps rebuilt); P0 tests it
on its own and puts the hooks in place. P1 adds the conservative prolongation FLUME regrids with
(`[amr] regrid_prolongation`, default `conservative`). P2 regrids FLUME on the CPU: the markers combined, the Löhner
estimator, the solids' distance function and the flux register rebuilt, restarts across regrids
([numerics](./numerics#adaptive-mesh-refinement)). P3 regrids the FNL backend by a host round trip.

### RG: the regrid round trip and the input contract

**Why.** Until M5 nothing ran `amr_update` after initialisation and no test ever derefined. A regrid must keep the
bookkeeping consistent over the ranks (the redistribution moves blocks) and transfer the data exactly where it can.

**The conservative prolongation.** Each parent cell is split along the refined axes into children set, variable by
variable, to

$$
q_c = q + \phi \sum_d \frac{\sigma_d\, s_d}{4}, \qquad \sigma_d = \pm 1,
$$

with $s_d$ the monotonized-central slope of the parent along axis $d$ (in parent-cell units, from its two face
neighbours) and $\phi \in [0, 1]$ the largest factor that keeps every child in the range of the parent and its face
neighbours (Barth and Jespersen 1989). The offsets cancel over the children, so their mean is the parent: the mean
restriction undoes the prolongation exactly and $\sum q\,\Delta V$ does not change. The MC limit bounds each axis on
its own, not the sum over the axes a corner child takes; $\phi$ closes that gap. Linear data keep $\phi = 1$ and are
reproduced exactly.

**Legs** (`regrid/check.sh`):

- rg0, `tests/amr/test_amr_regrid_roundtrip` on 1, 2 and 3 ranks: a realm refined uniformly to level 1 holds a field
  in every cell, ghosts included (a perfect ghost fill, so the transfers are tested alone); the block at the origin is
  refined, then its children derefined, with each prolongation. Fields: linear, quadratic, and a steep front from 1 to
  0.01 (a tanh half a coarse cell wide, across $x+y+z=0.45$ on the octree, $x+y=0.3$ on the quadtree). Asserted: the
  leaves and the blocks summed over the ranks agree (octree 8 → 15 → 8, quadtree 4 → 7 → 4); a linear field exact
  after the refine; with `conservative`, the round trip returns the initial state and $\sum q\,\Delta V$ is unchanged
  after the refine and after the round trip, on every field, and the front stays in $(0.01, 1)$; with `linear`, the
  quadratic field drifts (the negative control).
- rg1, the `frequency` contract on the sod-x input: 0 runs, a negative value is fatal, `n > 0` without markers is
  fatal, each with its message.
- rg2, the `regrid_prolongation` contract on sod-x: absent resolves to `conservative`, `linear` and `conservative` are
  taken as given, any other value is fatal.

**Results** (CPU): every leg passes on 1, 2 and 3 ranks, with the same numbers.

| Field | Tree | `linear`: drift of $\sum q\,\Delta V$ | `linear`: round trip | `conservative`: drift | `conservative`: round trip |
|---|---|---|---|---|---|
| quadratic | octree | $1.3\cdot10^{-4}$ | $3.4\cdot10^{-4}$ | 0 | $4\cdot10^{-17}$ |
| quadratic | quadtree | $2.7\cdot10^{-4}$ | $3.4\cdot10^{-4}$ | 0 | $2\cdot10^{-17}$ |
| front | octree | $2.6\cdot10^{-2}$ | $1.8\cdot10^{-1}$ | $4.8\cdot10^{-15}$ | $1\cdot10^{-16}$ |
| front | quadtree | $1.8\cdot10^{-2}$ | $1.3\cdot10^{-1}$ | $2.6\cdot10^{-16}$ | $1\cdot10^{-16}$ |

Drift is relative to the initial integral, after the refine; round trip is the largest $|q - f|$ over the field scale
after refine then derefine. The `linear` prolongation is second order (a quadratic field's children miss by
$2.6\cdot10^{-4}$ of the scale, the conservative one's by $1.1\cdot10^{-4}$) but changes the integral at every
refine. The bound scaling is load-bearing: with $\phi = 1$ forced (a mutation run), an octree child of the front
reaches $-0.094$, below zero, and the leg fails; on the quadtree the MC limit alone suffices for this front.

### RG: runtime regrids in FLUME

**Why.** A regrid during the run must not change what the scheme conserves, must leave the run reproducible across a
restart, and must keep every grid-dependent piece (the flux register, the solids' distance function) in step with the
grid. These legs regrid for real, on cases where the grid keeps changing.

**Legs** (`regrid/check.sh`, inputs from `regrid/make_regrid.py`):

- rg3, conservation through regrids: the isentropic vortex of V2 travelling across a periodic quadtree (base level 1,
  `max_level` 3, 8×8 cells per block), a Löhner marker on the density (`refine_tol` 0.3, `derefine_tol` 0.1, buffer 2),
  a regrid every 5 steps, reflux on, 200 steps. Asserted: the grid refines and coarsens during the run; with the
  conservative prolongation the five volume integrals stay constant within $10^{-13}$; with the linear one they drift
  above $10^{-10}$ (the negative control). Nothing crosses the periodic boundary, so a drift can only come from a
  regrid or a seam: the leg covers the prolongation, the restriction, the redistribution and the flux register rebuilt
  at every regrid, with seams that appear and disappear.
- rg4, restart across regrids: 30 steps against 20 steps, a restart and 10 more, the restart saved at step 20, a regrid
  step that changes the grid (asserted). The last fields bitwise and the histories identical. The forest regrids before
  the per-step output, so the restart holds the grid the next step runs on; saved before the regrid, it would not.
- rg5, immersed solid: the shock over the cylinder (V6) with the solid marker and a Löhner marker (`refine_tol` 0.6),
  a regrid every 5 steps, 60 steps: the run completes with a regrid that changes the grid, each regrid checking the
  new state (density and pressure positive). Its accuracy against the initial-AMR and uniform runs is #74 P4 (AV-7).

- rg6, CPU against FNL: the rg3 conservative runs of the two backends regrid at the same steps into the same grids (the
  regrid log lines equal) and end with the same fields within $10^{-10}$, relative to each variable's largest
  magnitude. The FNL backend regrids on the host, so the grid decisions are the same code on both; this leg checks
  that the device state reaches the host and comes back intact.

rg3, rg4 and rg5 run on both backends.

**Results** (np 2, the same on CPU and FNL except where noted):

| Leg | Regrids that changed the grid | Blocks refined / families coarsened | Measured |
|---|---|---|---|
| rg3 `conservative` | 10 | 20 / 4 | drift $2.2\cdot10^{-16}$ (ρ), $1.3\cdot10^{-16}$ (ρE) |
| rg3 `linear` | 7 | 14 / 6 | drift $2.6\cdot10^{-5}$ (ρ), $2.2\cdot10^{-4}$ (ρu) |
| rg4 | 3 (steps 5, 20, 25) | — | 40 blocks bitwise, histories identical |
| rg5 | 1 (176 → 344 blocks) | 24 / 0 | completes, admissible |
| rg6 | 10, the same on both | 20 / 4 | fields within $7.0\cdot10^{-14}$; FNL round trips 2.3 to 3.1 s each |

On FNL the conservative rg3 run drifts by at most $2.6\cdot10^{-16}$, the same round-off as the CPU. Its round trips cost
2.3–3.1 s each on a grid of at most 52 blocks: the copies move the whole device state, sized by the block capacity
(17195 blocks per rank here), not the blocks in use ([#75](https://github.com/szaghi/adam/issues/75)).

The two rg3 runs regrid differently: the linear prolongation changes the data the Löhner estimator reads, so the grids
part after the first regrid. The integrals still separate the two by eleven orders of magnitude.

P0 also fixes the capacity check of a regrid: the new blocks take indices on their parent's rank before the
redistribution, so the bound is per rank (`nb`), not the former `procs_number · nb` over every node of the replicated
tree, which let one rank overflow its arrays undetected.

### AV: accuracy of runtime AMR

**Why.** RG proves the regrid machinery. AV asks what a user of runtime AMR wants: is a run whose grid follows the
solution about as accurate as the uniform run at its finest level, with fewer cells, and does a regrid keep what the
scheme guarantees? Each case compares a tracked run (a regrid every 5 steps) with uniform runs, against bounds fixed in
the #74 P4 plan before any run. Five bounds turned out to rest on false premises, all pre-existing properties of the
scheme, of a marker or of the backends that runtime AMR merely exposed; they were reported and replaced, each with its
issue, never moved to fit a number.

**Legs** (`regrid/accuracy.sh`, inputs from `regrid/make_regrid.py`, checks by `regrid/av_oracle.py` on top of the V1,
V2, V6 and MV-9 oracles):

| Leg | Case | Asserted |
|---|---|---|
| av2 | isentropic vortex (V2), N = 32 tracked to 128, gradient marker on ρ (tol 0.05), t = 0.2 | L1(ρ) ≤ 1.5× uniform-128 and ≤ 1/8 of uniform-32 |
| av3a | Sod (V1), quadtree with y null, 48 cells tracked to 192 (gradient OR Löhner on ρ), on [0, 1.2] | L1(ρ) ≤ 1.3× uniform-192; contact and shock on finest cells; null copies within 1e-12 |
| av3b | Balsara–Spicer MHD blast (EGLM, positivity limiter), N = 64 tracked to 128 (gradient OR Löhner on p) | outer shock radius on the four half-axes within 2 finest cells of uniform-128; ρ conserved to 1e-13; ρE and B drift ≤ 1.1× uniform-128's; every regrid admissible |
| av4 | field loop (MV-9), GLM and EGLM, N = 64 tracked to 128 (gradient on B_x OR Löhner on B_x, B_y, floor 1e-4), t = 1 | E_B(T)/E_B(0) ≥ uniform-128's − 0.5 %; ⟨\|B_z\|⟩/A0 ≤ 1.2× (GLM; EGLM reported); final max\|div B\| ≤ 2×; div B at each regrid reported |
| av6 | the av2 tracked vortex on 1, 2, 3 ranks (1, 2 on FNL) | the same regrids; fields within 1e-13 |
| av7 | shock over the cylinder (V6) to t = 0.25, solid OR Löhner on ρ | bow-shock stand-off within 2 finest cells of uniform at V6's finest level; mirror symmetry 1e-10; positivity (the init-AMR V6 run's stand-off reported) |
| av8 | CPU against FNL on the tracked runs of av2, av3a, av3b, av4, av7 | the same regrids; per field, the tracked CPU–FNL difference ≤ 2× the uniform-fine runs' (or ≤ 1e-10) |

**Results** (np 2; CPU numbers, FNL the same within the av8 agreement):

| Leg | Tracked run | Reference | Measured | Regrids | Cost (tracked / uniform-fine) |
|---|---|---|---|---|---|
| av2 | L1(ρ) 1.126e-5 | uniform-128 8.43e-6; uniform-32 1.91e-3 | 1.34× fine, 0.006× coarse | 4 | 28 blocks / 64; 58 s / 65 s |
| av3a | L1(ρ) 4.1031e-3 | uniform-192 4.1029e-3 | 1.00003×; contact and shock on finest cells | 7 | 156 cells / 192; 28 s / 28 s |
| av3b | shock radius 0.4219–0.4375 | uniform-128 0.4297 | 1 finest cell on each half-axis; ρ drift 5.6e-16, ρE 5.33e-5 against 5.35e-5 | 3 | ends fully refined; 115 s / 122 s |
| av4 GLM | E_B 0.943937, ⟨\|B_z\|⟩ 5.75e-4, div B 2.56e-4 | 0.943953, 6.86e-4, 2.76e-4 | ⟨\|B_z\|⟩ 0.84×, div B 0.93×; at the regrids max \|div B\| changes by ×0.94–1.05 | 78 | 40–46 blocks / 64; 237 s / 244 s |
| av4 EGLM | E_B 0.943769, ⟨\|B_z\|⟩ 2.37e-5, div B 2.58e-4 | 0.943774, 2.72e-6, 2.72e-4 | ⟨\|B_z\|⟩ 8.7× (reported, #78), div B 0.95× | 82 | 242 s / 250 s |
| av6 | np 2 and 3 against np 1 (FNL: np 2, two GPUs) | — | the same 4 regrids, fields bitwise on both backends | 4 | 97 / 55 / 48 s |
| av7 | stand-off 0.12000 | uniform at V6's finest level 0.12000 (init-AMR V6: 0.12195) | 0 finest cells; mirror symmetry 6.5e-12 | 7 | ends fully refined; 315 s / 389 s |
| av8 | FNL against CPU, the six tracked runs | the same pair of uniform-fine runs | the same regrids in every case; Euler: CPU and FNL agree to round-off on both runs (tracked 5.5e-14 to 6.9e-12); MHD: the tracked difference is 0.5–1.7× the uniform one, per field (see below) | — | FNL tracked: av4 556 s, av7 1055 s (round trips, #75) |

The savings are modest on these cases, and the reason is visible: the tracked vortex and loop hold half the blocks of
the uniform runs, but the time step is set by the finest cells either way (no subcycling, D-M5-3) and the regrids, the
2:1 seams and the reflux cost their share, so the wall time drops by 3 to 19 % (av7 the most); the blast and the
cylinder end fully refined, their waves filling the box by the final time. The accuracy is within the stated factor of
the uniform-fine run in every case, and equal to it where the refined region covers the whole feature (av3a, av7).

**What runtime AMR exposed** (each reported before its bound changed):

- **A jump on a block face is invisible to the gradient marker** at initialisation (it reads interior cells only):
  Sod's jump at 0.5 is a face at every level on [0, 1], the initial grid stayed coarse, and the error of the first
  coarse steps persisted, 2.57× uniform-fine whatever the Löhner threshold. On [0, 1.2] the jump is off every face and
  the tracked run matches uniform-fine ([#76](https://github.com/szaghi/adam/issues/76)).
- **The Löhner estimator needs an absolute floor** on a variable that vanishes in quiet regions: outside the field loop
  B is GLM residue (~10⁻⁶), the relative filter `epsilon` filters nothing, and the estimator read E ≈ 0.98 everywhere,
  refining the whole box. The `floor` key (default 0, the FLASH/PLUTO form) fixes it
  ([numerics](./numerics#adaptive-mesh-refinement)).
- **The MHD blast breaks point symmetry on a uniform grid** (0.2 relative at N = 128 from an exactly symmetric initial
  condition), and does not conserve ρE and B (it drifts 5e-5 on the uniform grid too): the planned symmetry bound was
  dropped and the conservation bound limited to ρ, with ρE and B bounded against the uniform drift
  ([#77](https://github.com/szaghi/adam/issues/77)).
- **EGLM loses its B_z advantage across 2:1 seams**: ⟨|B_z|⟩/A0 is 2.7e-6 on the uniform grid, 2.4e-5 tracked and
  1.9e-4 on a static AMR grid with no regrid at all, so the seams, not the regrid, raise it; the tracked run, whose
  seams follow the loop, is 8× better than the static one. The EGLM bound on it is reported only
  ([#78](https://github.com/szaghi/adam/issues/78)).
- **On MHD the two backends do not agree to round-off**, with no regrid involved: on the uniform field loop the fluid
  variables agree to 1e-13 but B only to ~1e-6 ([#79](https://github.com/szaghi/adam/issues/79)), and the uniform
  blast differs at O(0.1) ([#77](https://github.com/szaghi/adam/issues/77)). The planned av8 bound (1e-10 on every
  case) tested that, not the regrid; av8 now bounds the tracked CPU–FNL difference by the uniform one, case by case.
- Across 2:1 seams the copies along a null direction differ by ~1e-13, on a static grid too: round-off of the seam
  ghost fill; av3a bounds it at 1e-12.

![Runtime regridding: tracked vortex and field loop with their blocks](/flume/regrid.png)

## Scaling covariance

Ideal Euler and MHD in FLUME's units carry no dimensionless number, so an input rescaled by powers of two (lengths
$\times 2^j$, velocities $\times 2^k$, density $\times 4^m$) is the same problem in other units, and a scheme without
absolute constants reproduces it bit for bit once the outputs are scaled back (issue #49).
`scaling/check.sh` runs Sod, the isentropic vortex, RJ2a and the limited Balsara–Spicer blast at their base scale and
at four rescalings (length only, velocity only, density only, all together; default $j = 2$, $k = -1$, $m = -2$) and
compares the last checkpoint of each with the base one scaled exactly (`scaling.py`).

| Weights | Result (CPU and FNL) |
|---|---|
| `js` (default) | length-only bitwise; every velocity or density rescaling differs by $10^{-4}$ to $10^{-1}$ (the absolute $\varepsilon$) |
| `si` (`--weights si --expect-bitwise`) | 16/16 bitwise; and the blast at density $\times 4^{-24}$, where the limiter floor acts, bitwise |

The second row needs the relative-only positivity floor and the unguarded immersed-boundary cut spacing (issue #49,
N1b). `scaling/weights-exe.sh` runs any verification script with `[weno] weights = si` (`FLUME_EXE` pointed at it,
`WEIGHTS_EXE` at the executable), so the orders and the shock errors of the new weights are judged by the same oracles
and bounds as the default ones. One tolerance differs: RJ2a's GLM-against-no-cleaning and EGLM-against-GLM legs are
bitwise with `js` but within $10^{-11}$ with `si` (measured $\le 2.5\cdot10^{-13}$), the round-off of a uniform field
reconstructed with face-varying weights ([numerics](./numerics#weno-reconstruction)).

### Reference layer (NV-5)

A dimensional input with a [`[reference]`](./input#reference-optional-dimensional-input) section must run as its
hand-normalised twin. `scaling/reference.sh` writes each case in the units $L_0 = 2^j$, $u_0 = 2^k$, $\rho_0 = 4^m$
(`scaling.py dimensionalize`, default $j = 2$, $k = -1$, $m = -2$) with the section that converts it back; the
references are powers of two, so the solver sees the base numbers exactly and the run must equal the base run bit for
bit (every conservative field, the history and slice files byte for byte). `scaling.py` classifies every option
independently of the Fortran layer, and `check-log` compares every conversion the layer logs with the base value (exact
equality, every dimensional option converted), which also covers the options whose conversion does not show in the
solution.

| Leg | Cases | Result (CPU and FNL) |
|---|---|---|
| regression cases | the 21 single-realm ones (orszag-tang and the 5 multi-realm ones are refused by the layer) | 21/21 bitwise, conversions exact and complete |
| generated inputs | the MHD linear waves (fast GLM, Alfvén EGLM, slow), CPAW, the magnetised vortex, the div(B) peak, the EGLM pulse, Shu–Osher, the isentropic vortex, Sod with a slice, Sod with two gradient AMR markers | 11/11 bitwise, conversions exact and complete |

NV-8 (`reference.sh --restart`): on sod-x (fast path), shock-cylinder-ib (immersed boundary, AMR),
amr-periodic-reflux (staged path, reflux) and blast-amr-limiter (MHD EGLM, AMR, positivity limiter), dimensionalised as
above, a run of 20 steps equals a run of 10 steps plus a restart to 20 bit for bit (conservative fields, residual and
conservation histories); restarting with another density reference, or without the `.reference` record, is refused.

NV-9 (`reference.sh --output`): on sod-x, shock-cylinder-ib, field-loop (GLM), blast-amr-limiter (EGLM, AMR) and
uniform-amr-mhd (3-D, AMR), each base input with residual and auxiliary fields and a slice along $x$, the
dimensionalised twin with `output_units = dimensional` computes the base numbers and writes them times the references:
every field, the block spacing, the XDMF times, every history, the slices and the `.units` record equal the base output
times its exact power of two, bit for bit (`scaling.py compare --output`). The temperature is the exception: the base
computes $p/(\rho R)$ with the regression inputs' $c_p$, $c_v$, the twin $(p/\rho)/R$ with $R = 1$ in the solver, a few
roundings apart ($2.3 \cdot 10^{-16}$ on sod-x; tolerance $2 \cdot 10^{-15}$). A power-of-two scaling is exact except
where it underflows: the FNL build flushes subnormal results to zero (`-fast`), so a written 0 matches an expected
subnormal (one `rv` of $3 \cdot 10^{-308}$ on the shock-cylinder slice, scaled by $2^{-5}$).

NV-6 (`scaling/interstellar.sh`): the V1 Sod problem in SI interstellar units, $L_0 = 1$ pc, $u_0 = 10$ km/s,
$\rho_0 = 10^{-21}$ kg/m³ (`scaling.py physical`; not powers of two), on both flux paths (`weno`, and `weno-riemann`
HLLC with the 6th-order correction and the `weno` sensor). With the classic weights the absolute `zeps = 1e-6`
dominates smoothness indicators of order $10^{-42}$: the weights collapse to the linear ones and the scheme stops
limiting, as issue #48 (A3) predicted from the formula. Measured: the total variation of $\rho/\rho_0$ grows from 0.877
to 1.162 (`weno`) and from 0.884 to 2.455 (`weno-riemann`, whose blind sensor leaves the unlimited correction on: L1 ×20,
an 8% undershoot). Both cures restore the code-units solution within round-off ($10^{-13}$ in $\rho/\rho_0$): the
scale-invariant weights (`[weno] weights = si`) on the raw SI input, and the reference layer with the classic weights.
The raw-SI classic run is kept as a negative control: its total variation must exceed the code-units one by 10%, or it
must stop on a non-finite state. The unlimited `weno-riemann` run undershoots towards vacuum and, round-off deciding,
survives on the CPU and stops at step 90 on the FNL.

The first run found two options the scaling tool had misclassified since N0: `rho_amplitude` taken as dimensionless
(it is a density) and `loop_amplitude` as a vector potential (it is the field of the loop). The two independent tables
disagreed, and the field-loop case was not bitwise.

## Unit tests

| Test | What it pins |
|---|---|
| `test_flume_euler_library` (+ `_fnl`) | Euler eigensystem, flux, Roe average, split consistency; RS(q, q) = f(q) for LLF/HLL/HLLC, HLLC exact on a contact, positive first-order updates; the `si` descaler of planar and static states within $10^8$ of the largest (host); device = host |
| `test_flume_mhd_library` (+ `_fnl`) | MHD, GLM and EGLM eigensystems (including degenerate states), fluxes, auxiliary variables, cyclic invariance; the `si` descaler of planar, field-free and static states within $10^8$ of the largest (host); device = host |
| `test_flume_mhd_riemann` (+ `_fnl`) | MHD LLF/HLL/HLLD without cleaning, with GLM and with EGLM: consistency, HLLD exact on contact, tangential and rotational discontinuities, cyclic invariance bitwise, positive updates, EGLM = GLM bitwise at $\psi = 0$; device = host |
| `test_flume_positivity` | PV-0, the limiter on random admissible states with perturbed fluxes (Euler, MHD, EGLM): every limited update positive and above the relative floor, the limiter needed and acting; a NaN or infinite high-order flux replaced by the backbone flux; inadmissible face ghosts blended above the floors, every other value untouched; a cell with a 2:1 seam face (mean donor-state backbone) admissible for any seam factor up to its own, and not without the limiter |
| `test_flume_weno_interpolation` (+ `_fnl`) | WENO interpolation tables: exactness, convergence, device = host for both weights (`js`, `si`); the reconstruction tables unchanged; the `si` weights scale-covariant bitwise on the device |
| `test_flume_weno_weights` | NV-1 of issue #49: with `[weno] weights = si`, WENO($2^n v$) = $2^n$ WENO($v$) bitwise for every scheme, both tables, random stencils from $10^{-12}$ to $10^{12}$; the `js` weights fail it (negative control); quadratic data reproduced by both; zero data stay zero |

Build and run with `fobis build --mode test-flume-<name>-gnu` (`-fnl-nvf --varset local_nvf` for the device twin).

## Regression suite

`src/tests/flume/regression/` holds 35 goldened cases (Sod along x/y/z, AMR, multi-realm (mirror seams, lined up or not, and 2:1 refined seams), immersed boundary, the MHD
cases RJ2a, Brio–Wu, GLM pulse, Orszag–Tang, rotor, field loop, rotated shock tube, uniform AMR, and the M3 cases: Sod
with HLLC and RJ2a with HLLD on `weno-riemann`, Orszag–Tang with EGLM, the blast with the positivity limiter in 2-D and
in 3-D, the only 3-D flow case of the suite, and, issue #50, the blast with the limiter across two 2:1 seams through
its centre; and the M4 dissipative cases, issue #65: the diagonal shear wave, Couette flow, Becker's shock, the
visco-resistive Alfvén wave and the same Alfvén wave across 2:1 seams; and the M5 runtime AMR cases, issue #74: the
periodic vortex and the field loop regridding every 5 steps, six and ten regrids that refine and coarsen): `run.sh cpu`
runs in
CI, `run-fnl-local.sh` on a GPU workstation, and `run-omp-bitwise.sh` checks that the OpenMP build reproduces the serial
one bit for bit.
