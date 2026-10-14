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
initialization, and so is the positivity limiter on any multi-realm run. Leg 6
([#68](https://github.com/szaghi/adam/issues/68) C2) repeats leg 4, along $x$ and $z$, with `[amr] seam_restriction =
lagrange6` in every input (the null axes widened to 8 cells, which the 6-point footprint needs): the fine owner packs
the 6-point codes in the restrict rows, and the unions stay bitwise equal to the single-realm runs, CPU and FNL. The
MHD twin (`mhd/multirealm/check.sh` leg 4,
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
decay of $\|\nabla\cdot\mathbf{B}\|_1$ from an initial peak of $\nabla\cdot\mathbf{B}$, on a uniform grid and across AMR seams,
and the `divb_tol` monitor on the normalised D_flux (warning, stop, and the refusal without `divb_flux`; #68 G1).

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
initial loop edge) and GLM removes it. `--agree` compares the CPU and FNL runs of the ladder field by field, B relative
to the loop amplitude, within $10^{-10}$: the backends agree to round-off on B as on the fluid
([#79](https://github.com/szaghi/adam/issues/79): they differed by ~10⁻⁶ until the eigenvector normalisation was made
cancellation-free, see [numerics](./numerics#scheme-space-weno-flux-splitting)).

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

**The polynomial probe** ([#68](https://github.com/szaghi/adam/issues/68) D1). `check.sh --degree N` starts from the
[`polynomial`](./initial-conditions#polynomial-a-polynomial-field-verification) field
$\mathbf q_1 (1 + \sum_d g_d x_d^N)$ instead. The copies and the boundary conditions still reproduce it, but each seam
operator only up to its degree (the mean restriction 1, `tricubic` 3, `lagrange6` and `quintic` 5), so at $N > 1$ the
probe checks the high-order operators in situ: their footprints, anchor shifts and every exchange path. With
`--degree 5 --seam-restriction lagrange6 --seam-ghost-fill quintic` on `box3d`, `channel2d` and `refined3d` every face
and edge group holds the field to $5.8\cdot10^{-15}$ (CPU 1–3 ranks, FNL 1–2). Negative controls, both required to
fail: `tricubic` at degree 5 misses by $\sim10^{-5}$ and the mean restriction at degree 2 by $\sim4\cdot10^{-5}$ of the
field scale, in the seam groups only.

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
| rg6 | 10, the same on both | 20 / 4 | fields within $7.0\cdot10^{-14}$; FNL round trips 0.015 to 0.029 s each (2.3 to 3.1 s before #75) |

On FNL the conservative rg3 run drifts by at most $2.6\cdot10^{-16}$, the same round-off as the CPU. Its round trips cost
0.015–0.031 s each on a grid of at most 52 blocks: the copies move the blocks in use, packed on the device. They moved
the whole device state, sized by the block capacity (17195 blocks per rank here), and cost 2.3–3.1 s each until
[#75](https://github.com/szaghi/adam/issues/75).

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
| av3a | Sod (V1), quadtree with y null, 48 cells tracked to 192 (gradient OR Löhner on ρ), on [0, 1], the jump on a block face at every level | L1(ρ) ≤ 1.3× uniform-192; contact and shock on finest cells; null copies within 1e-12 |
| av3b | Balsara–Spicer MHD blast (EGLM, positivity limiter), N = 64 tracked to 128 (gradient OR Löhner on p) | outer shock radius on the four half-axes within 2 finest cells of uniform-128; ρ conserved to 1e-13; ρE and B drift ≤ 1.1× uniform-128's; every regrid admissible |
| av3c | the av3b blast on `weno-riemann` (HLLD, 6th-order correction), uniform N = 128 and tracked from 64 | point symmetry within 1e-10 on both ([#77](https://github.com/szaghi/adam/issues/77)) |
| av4 | field loop (MV-9), GLM and EGLM, N = 64 tracked to 128 (gradient on B_x OR Löhner on B_x, B_y, floor 1e-4), t = 1 | E_B(T)/E_B(0) ≥ uniform-128's − 0.5 %; ⟨\|B_z\|⟩/A0 ≤ 1.2× (GLM; EGLM reported); final max\|div B\| ≤ 2×; div B at each regrid reported |
| av6 | the av2 tracked vortex on 1, 2, 3 ranks (1, 2 on FNL) | the same regrids; fields within 1e-13 |
| av7 | shock over the cylinder (V6) to t = 0.25, solid OR Löhner on ρ | bow-shock stand-off within 2 finest cells of uniform at V6's finest level; mirror symmetry 1e-10; positivity (the init-AMR V6 run's stand-off reported) |
| av8 | CPU against FNL on the tracked runs of av2, av3a, av3b, av4, av7 | the same regrids; per field, the tracked CPU–FNL difference ≤ 2× the uniform-fine runs' (or ≤ 1e-10) |

**Results** (np 2; CPU numbers, FNL the same within the av8 agreement):

| Leg | Tracked run | Reference | Measured | Regrids | Cost (tracked / uniform-fine) |
|---|---|---|---|---|---|
| av2 | L1(ρ) 1.039e-5 | uniform-128 8.43e-6; uniform-32 1.91e-3 | 1.23× fine, 0.005× coarse | 4 | 28 blocks / 64; 58 s / 67 s |
| av3a | L1(ρ) 3.5392e-3 | uniform-192 3.5420e-3 | 0.999× (2.56× before #76); contact and shock on finest cells | 7 | 144 cells / 192; 32 s / 33 s |
| av3b | shock radius 0.4219–0.4375 | uniform-128 0.4297 | 1 finest cell on each half-axis; ρ drift 5.6e-16, ρE 5.33e-5 against 5.35e-5 | 3 | ends fully refined; 115 s / 122 s |
| av3c | symmetry 3.3e-12 (FNL 2.6e-12) | uniform-128 6.1e-12 (FNL 2.6e-12) | both within 1e-10; shock radius 0.4297 on every half-axis, 0 finest cells from uniform | 3 | 128 s / 155 s |
| av4 GLM | E_B 0.943937, ⟨\|B_z\|⟩ 5.75e-4, div B 2.56e-4 | 0.943953, 6.86e-4, 2.76e-4 | ⟨\|B_z\|⟩ 0.84×, div B 0.93×; at the regrids max \|div B\| changes by ×0.94–1.05 | 78 | 40–46 blocks / 64; 237 s / 244 s |
| av4 EGLM | E_B 0.943769, ⟨\|B_z\|⟩ 2.37e-5, div B 2.58e-4 | 0.943774, 2.72e-6, 2.72e-4 | ⟨\|B_z\|⟩ 8.7× (reported, #78), div B 0.95× | 82 | 242 s / 250 s |
| av6 | np 2 and 3 against np 1 (FNL: np 2, two GPUs) | — | the same 4 regrids, fields bitwise on both backends | 4 | 97 / 55 / 48 s |
| av7 | stand-off 0.12000 | uniform at V6's finest level 0.12000 (init-AMR V6: 0.12195) | 0 finest cells; mirror symmetry 6.5e-12 | 7 | ends fully refined; 315 s / 389 s |
| av8 | FNL against CPU, the six tracked runs | the same pair of uniform-fine runs | the same regrids in every case; Euler: round-off on both runs (tracked 5.7e-14 to 6.9e-12); field loop: 1.4e-11 GLM, 3.0e-10 EGLM (relative to each field's scale, B_z ~1e-6; ~1e-6 before #79); blast: O(1) as on the uniform run (#77); every tracked difference at most 1× the uniform one but av3a's (1.8e-13, round-off) | — | FNL tracked: av4 60 s (550 s before the bounded copies of #75), av7 973 s |

The savings are modest on these cases, and the reason is visible: the tracked vortex and loop hold half the blocks of
the uniform runs, but the time step is set by the finest cells either way (no subcycling, D-M5-3) and the regrids, the
2:1 seams and the reflux cost their share, so the wall time drops by 3 to 19 % (av7 the most); the blast and the
cylinder end fully refined, their waves filling the box by the final time. The accuracy is within the stated factor of
the uniform-fine run in every case, and equal to it where the refined region covers the whole feature (av3a, av7).

**What runtime AMR exposed** (each reported before its bound changed):

- **A jump on a block face was invisible to the gradient marker** at initialisation (it read interior cells only):
  Sod's jump at 0.5 is a face at every level on [0, 1], the initial grid stayed coarse, and the error of the first
  coarse steps persisted, 2.56× uniform-fine whatever the Löhner threshold. The marker now differences across the
  faces between blocks, reading one ghost layer (the initial AMR refreshes the host ghosts first), and av3a runs on
  [0, 1] again: 0.999× ([#76](https://github.com/szaghi/adam/issues/76)).
- **The Löhner estimator needs an absolute floor** on a variable that vanishes in quiet regions: outside the field loop
  B is GLM residue (~10⁻⁶), the relative filter `epsilon` filters nothing, and the estimator read E ≈ 0.98 everywhere,
  refining the whole box. The `floor` key (default 0, the FLASH/PLUTO form) fixes it
  ([numerics](./numerics#adaptive-mesh-refinement)).
- **The MHD blast breaks point symmetry on a uniform grid** (0.2 relative at N = 128 from an exactly symmetric initial
  condition), and does not conserve ρE and B (it drifts 5e-5 on the uniform grid too): the planned symmetry bound was
  dropped and the conservation bound limited to ρ, with ρE and B bounded against the uniform drift
  ([#77](https://github.com/szaghi/adam/issues/77)): a grid-scale instability of the flux splitting at low β
  ([below](#the-mhd-blast-loses-its-point-symmetry-on-the-flux-splitting-path-77)); av3c asserts the symmetry on the
  Riemann-flux path.
- **EGLM loses its B_z advantage across 2:1 seams**: ⟨|B_z|⟩/A0 is 2.7e-6 on the uniform grid, 2.4e-5 tracked and
  1.9e-4 on a static AMR grid with no regrid at all, so the seams, not the regrid, raise it; the tracked run, whose
  seams follow the loop, is 8× better than the static one. The EGLM bound on it is reported only
  ([#78](https://github.com/szaghi/adam/issues/78)).
- **On MHD the two backends did not agree to round-off**, with no regrid involved: on the uniform field loop the fluid
  variables agreed to 1e-13 but B only to ~1e-6, a cancellation in the eigenvector normalisation since fixed
  ([#79](https://github.com/szaghi/adam/issues/79), MV-9 `--agree`), and the uniform blast differs at O(0.1)
  ([#77](https://github.com/szaghi/adam/issues/77)). The planned av8 bound (1e-10 on every case) tested that, not the
  regrid; av8 bounds the tracked CPU–FNL difference by the uniform one, case by case.
- Across 2:1 seams the copies along a null direction differ by ~1e-13, on a static grid too: round-off of the seam
  ghost fill; av3a bounds it at 1e-12.

![Runtime regridding: tracked vortex and field loop with their blocks](/flume/regrid.png)

### The MHD blast loses its point symmetry on the flux-splitting path ([#77](https://github.com/szaghi/adam/issues/77))

**Why.** The Balsara–Spicer blast of `blast-limiter` and av3b starts exactly point-symmetric. Ideal MHD maps it onto
itself under the rotation by π about the centre composed with B → −B: ρ, ρE, B even, ρu, ψ odd. On a uniform grid
the run broke the symmetry at O(0.1) and the two backends disagreed at O(0.1), with no AMR involved. #77 asks where
it comes from, and either fixes it or explains and bounds it.

**Setup.** The uniform blast of av3b: N = 128, β = 2.5·10⁻⁴ outside the hot disc, EGLM, `weno` flux splitting,
characteristic reconstruction, positivity limiter, t = 0.01 (428 steps); np 2, a checkpoint every 10 or 20 steps. The
asymmetry of a field is max |f(x) − s f(1 − x)| relative to max |f|.

| Probe | Result |
|---|---|
| HEAD, after the #79 fix (CPU) | ρ 0.16 pointwise; in L1 only 2.3·10⁻²: the pointwise figure is a shock shifted by a few cells. CPU against FNL 0.16 |
| first steps | 1.2·10⁻¹⁰ after one step, at the front cells in the low-β gas; exponential growth from step 60 (e-folding in 10–15 steps), saturation at O(0.1) by step 220 |
| N = 64 / 128 / 256 | the asymmetry reaches 10⁻² at t = 0.0094 / 0.0047 / 0.0024: the growth rate scales as 1/Δx, a grid-scale instability, not a physical one (whose rate would converge) |
| CFL 0.4 → 0.2 | 10⁻² at t = 0.0047 / 0.0041: the spatial operator, not the time integrator |
| β × 10 / × 100 (B / √10, B / 10) | 6.8·10⁻⁸ at t = 0.01, still growing / 2·10⁻¹⁰, flat: a low-β effect |
| conservative reconstruction | the same growth, about 100 steps later |
| limiter off (floors 0.01), EGLM or GLM | NaN at step 8: the limiter keeps the run alive |
| `weno-riemann`: HLLD with and without the 6th-order correction, LLF, HLLD with primitive interpolation | symmetric to 10⁻¹² – 10⁻¹⁰ through t = 0.01 |
| positivity limiter activity | flux splitting: faces limited at every step (428 of 428; β × 100: 46 of 384); `weno-riemann`: only in the first 17–27 steps. The backbone flux the limiter falls back to is the same on both paths |

**Diagnosis.** At β ≈ 2.5·10⁻⁴ the `weno` flux splitting keeps producing nearly inadmissible states at the blast
front, the positivity limiter clips them at every step, and the round-off asymmetry of the arithmetic grows at the
grid scale from 10⁻¹⁰ to O(0.1). The Riemann-flux path does not produce those states past the first steps and stays
symmetric. The limiter's part rests on correlation and elimination (it acts at every step exactly where the scheme
is unstable, and the scheme cannot run without it), not on a direct measurement of its amplification.

**A fix tried and not kept.** One Lax–Friedrichs speed for the seven MHD waves (the stencil's fast bound) instead of one
per wave, the cleaning pair unchanged, delays the growth (N = 128: 5.8·10⁻⁵ at t = 0.01; N = 256 saturates at
4·10⁻³ instead of 0.17) but does not stop it, keeps the limiter acting at every step, and raises the MV-5 errors 2 to 4
times on smooth waves (orders kept, every fixed bound failed).

**Resolution (accept and bound).** On the flux-splitting path at β ≲ 10⁻³ the point symmetry, and the agreement of
the two backends, hold only to the grid-scale instability above: the shock positions agree within a few cells, the
pointwise fields do not. That path keeps its golden (`blast-limiter`: each backend is compared with its own golden,
which a deterministic run reproduces; symmetry is not asserted). For low-β problems where symmetry or backend agreement matter, use
`scheme_space = weno-riemann` with `riemann_solver = hlld`. **av3c** bounds it: the av3b blast on that path,
point-symmetric within 10⁻¹⁰ on the uniform and on the tracked run (measured: uniform 6.1·10⁻¹² CPU, 2.6·10⁻¹²
FNL; tracked 3.3·10⁻¹² CPU, 2.6·10⁻¹² FNL, 3 regrids; the shock radius on the four half-axes the same on both runs).

## Seam accuracy and divergence (M6)

Milestone M6 ([issue #68](https://github.com/szaghi/adam/issues/68)) measures, and then lifts, the accuracy cap of the
static 2:1 seams. Phase A adds the instruments; the transfer operators and the reflux policy change only behind new
keys, whose defaults keep today's behaviour.

### AO: order across static 2:1 seams

**Why.** Every smooth order test (V2, MV-5 to MV-7) runs on a uniform grid; VV-8 measures a seam with dissipative
waves. No test measured the inviscid order of a smooth feature crossing a seam, where the #68 audit (F1) predicts that
the mean restriction of the coarse ghosts and the Berger–Colella reflux of point-value fluxes each cap the order at 2.

**Cases** (`verification/amr-order/`). A vortex starts centred on the corner of a refined quadrant, so its core
straddles an x seam, a y seam and their corner, and the free stream $(1, 1)$ carries it into the fine quadrant:
- `euler`: the V2 isentropic vortex on $[0, 1]^2$, quadrant $[0.5, 1]^2$ refined, $t = 0.1$, CFL 0.4 (V2's 0.1 gives
  the same errors to five digits);
- `mhd`: the MV-7 magnetised vortex on $[-7, 7]^2$, GLM, quadrant $[0, 7]^2$ refined, $t = 0.5$;
- both on a quadtree, base $N$ = 64 / 128 / 256, each with a uniform control (no refined quadrant). An octree leg,
  `euler-octree`, runs on request: it refines along the null z too, at about 16 times the cost.

**Oracle** (`amr_order_oracle.py`). Each cell is compared with the exact convected vortex (the density for Euler; the
eight conservative variables for MHD, combined as MV-7 does). The error is split into the **seam band**, the cells
within `ngc` of their own spacing from the boundary of the refined quadrant (periodic images included), and the
**interior**. Bounds: the finest refined errors at most the CPU baseline plus 2 %, so that an improvement passes and a
degradation fails; the uniform control at L1 order $\ge 4.7$.

**Results** (CPU and FNL identical to every printed digit, 2 ranks; orders of the finest pair, N = 128 → 256):

| Leg | Refined L1 at N = 256 | Uniform L1 at N = 256 | L1 order: all / band / interior | L∞ order | Uniform L1 order |
|---|---|---|---|---|---|
| euler | $4.05\cdot10^{-6}$ | $8.38\cdot10^{-8}$ | 2.01 / 1.71 / 1.91 | 1.46 | 5.79 |
| mhd | $1.83\cdot10^{-6}$ | $2.62\cdot10^{-8}$ | 2.07 / 1.85 / 1.92 | 1.89 | 5.34 |

Three findings:
- **The seams cap the whole solution at second order, not only the band.** The error made at a seam is carried into
  the interior by the flow, so every region a feature reaches after crossing a seam is second order. The plan's
  expectation of "design order in the interior" holds only for regions the flow does not reach from a seam.
- **On this smooth flow the refinement makes the solution worse.** At N = 256 the refined composite error is 48 times
  (Euler) and 70 times (MHD) the error of the uniform run at the coarse spacing, which reaches design order.
- **Euler's L∞ order is 1.46,** below the 2 of the 1-D prototype in #68 (Appendix A). The vortex crosses the seam
  corner, where the two seam directions meet; the 2-D prototype of Phase B will check whether the corner explains it.

Phases C (high-order restriction), D (quintic prolongation) and E (reflux policy) are measured against these
baselines.

**C1: the 6-point restriction** (`[amr] seam_restriction = lagrange6`, [#68](https://github.com/szaghi/adam/issues/68)
C1, C2). The coarse seam ghosts take the tensor 6-point Lagrange interpolation of the fine cells to the coarse centre
(exact 1D weights $k/256$, the footprint shifted inside the donor block at its edge, as B3 recommends) instead of the
mean. The leg `euler-l6` runs the euler ladder with it and the reflux off; CPU and FNL (the device twin of C2) agree to
every printed digit, and on the same run to $3\cdot10^{-14}$:

| N = 256 | L1 all | L∞ all | Orders 128 → 256, L1 / L∞ |
|---|---|---|---|
| euler (mean, reflux on, today) | $4.05\cdot10^{-6}$ | – | 2.01 / 1.46 |
| euler-l6 (lagrange6, reflux off) | $8.55\cdot10^{-8}$ | $6.36\cdot10^{-6}$ | **5.33 / 5.53** |
| uniform at the coarse spacing | $8.38\cdot10^{-8}$ | $5.34\cdot10^{-6}$ | 5.79 / 7.38 |

With the 6-point restriction and without the reflux the refined composite error equals the uniform coarse run's, 47
times below today's: at these resolutions the restriction was the seam error, and the tricubic prolongation's $O(h^4)$
term has not shown yet. With the reflux on, the restriction changes nothing (band L∞ order 1.43 against 1.46 at N =
64 → 128): the reflux caps the order at 2 by itself, as B1 predicts. No restriction row fell back to the mean (every
fine pair lies inside its donor block). Reaching this with the reflux kept at the shocks is Phase E's job. On FNL the
6-point restriction costs about 16 % of the wall time of this run (N = 256, reflux off: 56 s against 48 s with the mean,
WSL, indicative only): the restriction rows read 36 cells instead of 4 on this quadtree, in the same kernels as the
copies.

### AD: the divergence the cleaning acts on, across seams

**Why.** Every div(B) column above is the centred divergence of the cell values, which is not the divergence GLM and
EGLM act on (#68 F4). `[diagnostics] divb_flux` adds D_flux, the divergence of the interface normal field read from
the scheme's own $\psi$ face fluxes, the coarse 2:1 faces taking the fine fluxes; its signed integral over a periodic
composite domain is round-off (measured: $\le 10^{-17}$ every step on the refined MHD vortex, against
$4\cdot10^{-7}$ without the fine fluxes), and CPU and FNL agree to $10^{-11}$ relative.

**Cases** (`verification/amr-divb/`). `cpaw`: the source-free CPAW of MV-6 under GLM, inclined, its two middle block
columns refined (two x seams the wave crosses all the time), 10 periods, refined and uniform N = 64 / 128 (and 256 on
request). `loop`: the field loop of MV-9 under GLM and EGLM, refined N = 64 against uniform N = 128, $t = 1$
([#78](https://github.com/szaghi/adam/issues/78)).

**Oracle** (`amr_divb_oracle.py`). `cpaw`: the seam D_flux floor (median over the second half), its growth (the mean
over the last quarter over that of the second quarter) at most 1.5, and the floor's order at least 1. `loop`: a report.

**Results** (FNL, 2 ranks):

| N | Seam D_flux floor | Growth | Uniform floor | Seam / uniform |
|---|---|---|---|---|
| 64 | $5.39\cdot10^{-4}$ | 0.95 | $9.97\cdot10^{-6}$ | 54 |
| 128 | $5.43\cdot10^{-5}$ | 0.72 | $4.57\cdot10^{-7}$ | 119 |
| 256 | $2.12\cdot10^{-5}$ | 0.80 | $6.74\cdot10^{-9}$ | 3146 |

| Field loop, $t = 1$ | D_flux (seam, refined) | D_flux (max, uniform) | ⟨\|B_z\|⟩/A0 refined | uniform |
|---|---|---|---|---|
| GLM | $2.35\cdot10^{-4}$ | $2.03\cdot10^{-4}$ | $1.10\cdot10^{-3}$ | $6.86\cdot10^{-4}$ |
| EGLM | $2.45\cdot10^{-4}$ | $2.05\cdot10^{-4}$ | $1.91\cdot10^{-4}$ | $2.72\cdot10^{-6}$ |

Findings:
- **GLM keeps the seam divergence bounded.** Over 10 periods of a source-free flow the seam D_flux has no secular growth
  at any resolution: the upwinded $(B_n, \psi)$ pair of FLUME bounds what PRISM's centred GLM lets run away
  ([#29](https://github.com/szaghi/adam/issues/29); #68 F3).
- **The seam floor converges at first order** (orders +3.31 then +1.36: the first pair is pre-asymptotic), the O(H)
  divergence of the O(H²) flux inconsistency of a 2:1 face (#68 F1), while the uniform floor converges at 4 to 6: the
  seams are the divergence error of a refined run, 3146 times the uniform level at N = 256.
- **#78 is not a divergence difference.** GLM and EGLM leave the same D_flux, refined and uniform; EGLM's 250-fold
  smaller B_z on a uniform grid comes from its source term, whose centred divergence cancels the B_z the flux generates.
  At a seam that centred divergence reads the second-order restricted ghosts, the cancellation fails, and B_z rises 70
  times (reproduced here to the printed digits of #78). The fix is in Phase G3: a source divergence consistent with the
  face fluxes at the seams. **Corrected by [G3a](#g3a-78-is-the-reflux-not-the-restricted-ghosts):** the ghosts are
  not the cause; the reflux is.

### AS: a shock crossing a static 2:1 face

**Why.** The committed shock-across-a-seam checks (the multirealm legs, the `sod-amr` goldens) compare FLUME runs with
each other, so an error common to all of them passes. AS compares with the exact solution.

**Cases** (`verification/amr-shock/`). Sod (V1) and RJ2a (MV-4, no divergence control) on the V1 grid, four blocks
along x, with one 2:1 face at $x = 0.75$; the right shock (RJ2a: the fast shock) crosses it at $t = 0.143$ (0.110) and
is past it at $t = 0.2$. Two sides: `c2f` refines $x > 0.75$ (the shock goes from coarse to fine), `f2c` refines
$x < 0.75$ (fine to coarse). Base $N$ = 200 (Sod) and 256 (RJ2a), with uniform controls at $N$ and $2N$; `--long`
adds the refined runs at $2N$. The octree refines the null y and z too, so the refined RJ2a runs cost 5 to 14 minutes
on FNL: run that leg there.

**Oracle** (`amr_shock_oracle.py`). At $t = 0.2$ the face sits in the constant state behind the shock, the plateau,
so every deviation there is error with a known reference:
- mass: the relative defect of the mass integral against the initial mass plus $t$ times the net inflow at the ends;
- shock: the mass-equivalent shock position (the step holding the mass of the cells within $W = 0.03$ of the exact
  shock) minus the exact one, in local cells;
- plateau: between the left wave (Sod: the contact; RJ2a: the rotational discontinuity) and the shock, $W$ away from
  both: $\max|q - q_2| / |q_1 - q_2|$ for $\rho$, $u_n$ and $p$, and $\max|s/s_2 - 1|$ of $s = p/\rho^\gamma$;
- band: the same over the cells within `ngc` of their own spacing from the face.

$W$ is a length, the same for every run: a margin in cells of the run would let the contact, wider in cells on the
fine controls, into their plateau. Bounds: mass $\le 10^{-8}$ (on RJ2a at N = 256 the uniform run itself leaves
$3\cdot10^{-10}$: numerical precursors of the left fast rarefaction reach the extrapolated end), $|$shock$| \le 0.5$
cells, every plateau and band metric at most the baseline plus 10 %.

**Results** (2 ranks; Sod: CPU and FNL agree to every printed digit; RJ2a: FNL). Plateau maxima of the refined run
and, in brackets, of its uniform control at the spacing of the side that holds the shock:

| Case | Shock (cells) | $\rho$ | $u_n$ | $p$ | $s$ |
|---|---|---|---|---|---|
| Sod c2f, N = 200 | 0.16 (0.17) | $2.5\cdot10^{-2}$ ($1.3\cdot10^{-4}$) | $4.5\cdot10^{-3}$ ($1.7\cdot10^{-5}$) | $6.9\cdot10^{-3}$ ($2.4\cdot10^{-5}$) | $1.8\cdot10^{-2}$ ($9.4\cdot10^{-5}$) |
| Sod c2f, N = 400 | 0.16 (0.17) | $2.4\cdot10^{-2}$ ($1.3\cdot10^{-4}$) | $4.7\cdot10^{-3}$ ($1.2\cdot10^{-6}$) | $7.2\cdot10^{-3}$ ($1.5\cdot10^{-6}$) | $1.7\cdot10^{-2}$ ($9.4\cdot10^{-5}$) |
| Sod f2c, N = 200 | 0.18 (0.17) | $1.6\cdot10^{-2}$ ($9.5\cdot10^{-4}$) | $6.7\cdot10^{-3}$ ($1.1\cdot10^{-4}$) | $1.0\cdot10^{-2}$ ($1.2\cdot10^{-4}$) | $1.2\cdot10^{-2}$ ($7.1\cdot10^{-4}$) |
| Sod f2c, N = 400 | 0.18 (0.17) | $1.5\cdot10^{-2}$ ($1.3\cdot10^{-4}$) | $6.6\cdot10^{-3}$ ($1.7\cdot10^{-5}$) | $1.0\cdot10^{-2}$ ($2.4\cdot10^{-5}$) | $1.1\cdot10^{-2}$ ($9.4\cdot10^{-5}$) |
| RJ2a c2f, N = 256 | 0.11 (0.11) | $1.4\cdot10^{-2}$ ($1.9\cdot10^{-4}$) | $8.1\cdot10^{-4}$ ($8.7\cdot10^{-5}$) | $6.2\cdot10^{-3}$ ($7.4\cdot10^{-5}$) | $5.9\cdot10^{-3}$ ($5.2\cdot10^{-5}$) |
| RJ2a f2c, N = 256 | 0.11 (0.11) | $1.9\cdot10^{-2}$ ($2.1\cdot10^{-4}$) | $9.1\cdot10^{-4}$ ($2.1\cdot10^{-4}$) | $6.5\cdot10^{-3}$ ($2.0\cdot10^{-4}$) | $7.6\cdot10^{-3}$ ($5.8\cdot10^{-5}$) |

With the reflux turned off (`make_amr_shock.py --no-reflux`, a measurement only): mass defect $1.2$ to
$2.2\cdot10^{-4}$ (Sod) and $0.9$ to $1.5\cdot10^{-4}$ (RJ2a), against at most $3\cdot10^{-10}$ with it; the shock
off by 0.2 to 0.35 cells; the plateau maxima at most 30 % larger than with the reflux.

Findings:
- **The reflux does what it is for.** Mass is exact (to the floor the uniform RJ2a run shares), and the shock lands
  where the uniform run puts it. Without the reflux the mass drifts by $10^{-4}$ and the shock moves by a fraction of a
  cell.
- **A shock crossing a 2:1 face leaves an O(1) error behind it.** The plateau carries a density and entropy defect of
  0.6 to 2.5 % of the shock jump and acoustic waves of 0.08 to 1 % in velocity and pressure: 4 to 4800 times the
  uniform level. This error **does not converge**: from N = 200 to 400 the maxima stay the same while the uniform
  level falls up to 14 times; only the width of the defect shrinks (the mean entropy error halves).
- **It is not the reflux.** Without the reflux the plateau defect changes by at most 30 %. Its source is the
  coarse-fine transfer of a discontinuity, the mean restriction and the prolongation of the ghosts, which Phases C
  and D replace; AS is their shock baseline (D2's "no worse than tricubic").
- **The entropy defect leaves the face.** It moves with the flow, $u_2 t$ downstream of the face by $t = 0.2$, so the
  band entropy maxima are 100 to 1900 times below the plateau ones; the acoustic maxima of Sod c2f still sit at the
  face.

**C3: the 6-point restriction at shocks** ([#68](https://github.com/szaghi/adam/issues/68) C3). `lagrange6` is linear
and unlimited, so it could overshoot where a shock lies in its footprint. `make_amr_shock.py --seam-restriction R`
reruns the AS cases under either restriction; the octree refines the null $y$ and $z$ too, so their 4-cell blocks are
widened to 8, which the 6-point footprint needs (the `mean` run on the widened grid reproduces the AS baseline to every
printed digit). FNL, 2 ranks:

| Case | Plateau $\rho$ | Plateau $p$ | Plateau $s$ | Band $\rho$ | Band $p$ | Shock (cells) |
|---|---|---|---|---|---|---|
| Sod c2f, mean / lagrange6 | $2.463$ / $2.457\cdot10^{-2}$ | $6.936$ / $6.911\cdot10^{-3}$ | $1.792$ / $1.787\cdot10^{-2}$ | $6.352$ / $6.330\cdot10^{-3}$ | $6.936$ / $6.911\cdot10^{-3}$ | 0.158 / 0.158 |
| Sod f2c, mean / lagrange6 | $1.576$ / $1.575\cdot10^{-2}$ | $1.018$ / $1.015\cdot10^{-2}$ | $1.180$ / $1.179\cdot10^{-2}$ | $3.481$ / $3.430\cdot10^{-3}$ | $3.718$ / $3.663\cdot10^{-3}$ | 0.177 / 0.177 |
| RJ2a c2f, mean / lagrange6 | $1.413$ / $1.410\cdot10^{-2}$ | $6.220$ / $6.218\cdot10^{-3}$ | $5.906$ / $5.895\cdot10^{-3}$ | $1.005$ / $1.025\cdot10^{-4}$ | $1.402$ / $1.300\cdot10^{-4}$ | 0.109 / 0.109 |

The PV-5 blast (the limited Balsara–Spicer blast across 2:1 seams, `nk` widened to 8) passes under both restrictions
and both schemes: $t = 0.01$ reached, no inadmissible stage, no non-finite flux, mass and seam flux mismatch unchanged;
the largest number of faces limited in a stage is 46336 against 48512 with the mean (split) and 17920 for both (HLLD).

Findings:
- **The 6-point restriction does not worsen a shock crossing.** Every plateau and band defect is within 2 % of the
  mean's (the RJ2a band density 2 % above, the others equal or below), and the shock lands in the same place with the
  mass exact. The defect behind a crossing shock (above) does not come from the restriction's order.
- **It does not worsen positivity.** The limited blast passes as with the mean, with no more faces limited.
- **On FNL it is slow on an octree.** The PV-5 runs take 2.4 to 2.5 times the wall time of the mean (4900 s against
  1955 s, split; 4395 s against 1841 s, HLLD, with the same 433 steps; WSL, indicative only), against 16 % on the
  quadtree AO case: a restricted ghost reads 216 fine cells, one after the other, in the loop that copies the other
  ghosts. The nonlinear `weno6` follows as an option, not as a fix.
- **The cost was the ghost kernel, for every regime (C3b).** Profiled (Nsight Systems, 40 PV-5 steps, one rank), the
  local ghost update was 98 % of the FNL GPU time even with the mean: 0.70 s per call. nvfortran parallelised its loop
  over the variables only (9 threads) and ran the map rows in sequence (a "loop carried dependence" through the
  indirect writes): the schedule [#27](https://github.com/szaghi/adam/issues/27) specifies the fix for, and that
  [#53](https://github.com/szaghi/adam/issues/53) measured as the dominant FNL cost. With the two loops collapsed, as
  the other map kernels are, the kernel takes 0.84 s over the 444 calls instead of 310 s (mean) and 1.39 s instead of
  815 s (lagrange6), bitwise identical; the full PV-5 split run drops from 1955 s to 276 s (mean) and from 4900 s to
  288 s (lagrange6, whose premium becomes 4 %, with 5 % more steps). WSL timings: #27 asks for the benchmark on a real
  node.

**C3c: the WENO restriction** (`[amr] seam_restriction = weno6`, [#68](https://github.com/szaghi/adam/issues/68) C3).
Per direction, the WENO-Z combination (Borges et al. 2008, $p = 2$) of the cubic interpolants of the 4-point
sub-stencils of the 6 nodes that contain the coarse centre, with Jiang–Shu indicators and $\tau$ the highest undivided
difference of the candidates' union. An interval whose difference exceeds twice its neighbours' is a jump; a sub-stencil
holding one is dropped, and with none left (a jump inside the fine pair) the restriction is the mean of the pair, the
troubled-cell fallback. No sub-stencil extrapolates (across a jump it would import the far side's state), so the design
order is 6 centred, 5 one cell from the donor block's edge and 4 at it. Applied direction by direction (x, y, z) on CPU,
FNL and the inter-realm rows; the regime travels in the row's meta (1..216 Lagrange, 217..432 WENO), so the
inter-realm packers need no new argument. Unit test (`test_seam_interpolation`): exact on cubics at any scale and
offset, design orders on smooth data, the fine pair's own state across a jump elsewhere and the mean of a jump inside
the pair, mirror-symmetric on random data.

| Check (FNL, 2 ranks) | mean | lagrange6 | weno6 |
|---|---|---|---|
| AO euler, reflux off, N = 256: L1 all / band, Linf | — | $8.55\cdot10^{-8}$ / $2.61\cdot10^{-7}$, $6.36\cdot10^{-6}$ | $8.43\cdot10^{-8}$ / $2.26\cdot10^{-7}$, $6.39\cdot10^{-6}$ |
| AO orders L1 / Linf (128 → 256) | — | 5.33 / 5.53 | 5.34 / 5.53 |
| Sod c2f plateau $\rho$, band $\rho$ | $2.463\cdot10^{-2}$, $6.352\cdot10^{-3}$ | $2.457\cdot10^{-2}$, $6.330\cdot10^{-3}$ | $2.464\cdot10^{-2}$, $6.365\cdot10^{-3}$ |
| Sod f2c plateau $\rho$, band $\rho$ | $1.576\cdot10^{-2}$, $3.481\cdot10^{-3}$ | $1.575\cdot10^{-2}$, $3.430\cdot10^{-3}$ | $1.576\cdot10^{-2}$, $3.516\cdot10^{-3}$ |
| RJ2a c2f plateau $\rho$, band $p$ | $1.413\cdot10^{-2}$, $1.402\cdot10^{-4}$ | $1.410\cdot10^{-2}$, $1.300\cdot10^{-4}$ | $1.414\cdot10^{-2}$, $1.360\cdot10^{-4}$ |
| PV-5 split / HLLD: most faces limited in a stage | 48512 / 17920 | 46336 / 17920 | 48064 / 17920 |
| PV-5 split wall (WSL, after C3b) | 276 s | 288 s | 278 s |

Findings:
- **On the smooth vortex `weno6` equals `lagrange6`** (the band L1 13 % lower): no jump is detected, and the
  fourth-order rows at the donor block's edge do not show at these resolutions. AO leg `euler-w6` asserts it.
- **At shocks it is as neutral as `lagrange6`**: shock positions identical, mass exact, every defect within 1.5 % of
  the mean's; PV-5 passes under both schemes with no inadmissible stage. Neither 6-point restriction makes a shock
  crossing worse, so the default stays `mean` until Phase D/E evidence; `weno6` is the option for data with jumps
  inside the restriction footprint that these cases do not produce.
- Multirealm leg 6 runs both regimes (bitwise unions), GP probes `weno6` on box3d, channel2d and refined3d.

### D1: the quintic coarse-to-fine fill

**Method** (`[amr] seam_ghost_fill = quintic`, [#68](https://github.com/szaghi/adam/issues/68) D1). The tensor 6-point
Lagrange interpolation at the fine ghost's centre (1-D weights $k/8192$, exact in binary, `SEAM_W_QUINTIC`), the
footprint shifted inside the coarse donor block at its edge as the tricubic one is (blocks of at least 6 cells along
the refined axes); along an unrefined axis only the anchor is read. CPU, FNL (summed straight from the device field,
no per-thread footprint) and the inter-realm rows; the anchor shifts travel in the fill meta (bits 18–26). The default
stays `tricubic`.

| Check (FNL, 2 ranks) | lagrange6 + tricubic | lagrange6 + quintic |
|---|---|---|
| AO euler, reflux off, N = 256: L1 all / band, Linf | $8.55\cdot10^{-8}$ / $2.61\cdot10^{-7}$, $6.36\cdot10^{-6}$ | $6.03\cdot10^{-8}$ / $1.74\cdot10^{-7}$, $4.78\cdot10^{-6}$ |
| AO orders L1 / Linf (128 → 256) | 5.33 / 5.53 | 6.11 / 6.35 |
| Sod c2f plateau $\rho$ / $p$, band $\rho$ | $2.457$ / $0.691$, $0.633\cdot10^{-2}$ | $2.631$ / $0.858$, $0.795\cdot10^{-2}$ |
| Sod f2c plateau $\rho$ / $p$, band $\rho$ | $1.575$ / $1.015$, $0.343\cdot10^{-2}$ | $1.573$ / $0.985$, $0.338\cdot10^{-2}$ |
| RJ2a c2f plateau $\rho$ / $p$, band $\rho$ | $1.410$ / $0.622$, $1.03\cdot10^{-4}$ | $2.297$ / $0.986$, $0.68\cdot10^{-4}$ |
| PV-5 split / HLLD: most faces limited, inadmissible stages | 46336 / 17920, 0 | 45644 / 17920, 0 |

(Plateau values in units of $10^{-2}$; band values absolute.) Findings:
- **The smooth order reaches the scheme's.** With the reflux off, the 6-point restriction and the quintic fill, AO
  euler is sixth order (6.11 / 6.35) and its error at N = 256 is below the uniform coarse run's ($8.38\cdot10^{-8}$):
  the Phase D exit (L1 $\ge 4.7$, Linf $\ge 4.5$) holds, and the 2:1 seams no longer cost accuracy on this flow.
- **Unlimited, it worsens a shock entering a finer region.** On the c2f crossings the plateau defect behind the
  shock grows (Sod +7 % in $\rho$, +24 % in $p$; RJ2a +63 % / +59 %): the shock's coarse values reach the fine ghosts
  through a 6-point stencil with larger overshoots than the 4-point one. f2c is unchanged, shock positions are the
  same, mass is exact, and PV-5 passes with no inadmissible stage. The Phase D exit requires A3 no worse than
  tricubic, which is D2's limited variants' job; until then `quintic` is an opt-in for smooth flows.
- AO leg `euler-q6` asserts the order; GP checks the fill on a degree-5 field (above); multirealm leg 6 runs it across
  a refined inter-realm seam, bitwise; np 1/2/3 bitwise on the CPU, np 1/2 on FNL, CPU vs FNL $2.7\cdot10^{-13}$.

### D2: limited coarse-to-fine fills, and what the shock-crossing defect is not

**Methods** ([#68](https://github.com/szaghi/adam/issues/68) D2; both opt-in, the default stays `tricubic`):
- `seam_ghost_fill = weno6`: per direction the WENO-Z combination of the cubics of the 4-point sub-stencils that
  contain the ghost's interval (the C3 restriction's kernel, `seam_weno6_core`, now shared), linear weights reproducing
  the Lagrange interpolant of their union (all positive), sub-stencils holding a jump dropped, the clamped linear
  interpolation of the interval when none is left. Order 6 at the centred anchor positions, 5 one cell from the
  donor block's edge, 4 at it.
- `seam_ghost_fill = quintic-mp`: the D1 quintic clamped per direction to a Suresh–Huynh monotonicity-preserving bound
  (the interval's values and the curvature-corrected linear value $L + \theta(\theta-1)/2\, d^{M4}$, widened by the
  curvature term; at the block edge, where the target lies a quarter cell outside the footprint, the extrapolation is
  allowed only when the edge-side intervals hold no jump). Order 6 at every target.

**The diagnostic.** Neither limited fill brought the shock crossings of A3 to `tricubic`'s level, so every fill was run
on the c2f crossings (FNL, 2 ranks, `lagrange6` restriction; plateau and band defects as in A3):

| Fill (order) | Sod plateau $\rho$ | Sod band $\rho$ | Sod shock (cells) | RJ2a plateau $\rho$ | RJ2a band $\rho$ |
|---|---|---|---|---|---|
| `injection` (1) | $2.98\cdot10^{-2}$ | $7.91\cdot10^{-3}$ | 0.173 | $2.49\cdot10^{-2}$ | $2.16\cdot10^{-4}$ |
| `restriction-compatible` (2) | $2.32\cdot10^{-2}$ | $5.31\cdot10^{-3}$ | 0.164 | $1.98\cdot10^{-2}$ | $1.28\cdot10^{-4}$ |
| `tricubic` (4) | $2.46\cdot10^{-2}$ | $6.33\cdot10^{-3}$ | 0.158 | $1.41\cdot10^{-2}$ | $1.03\cdot10^{-4}$ |
| `quintic` (6) | $2.63\cdot10^{-2}$ | $7.95\cdot10^{-3}$ | 0.159 | $2.30\cdot10^{-2}$ | $0.68\cdot10^{-4}$ |
| `weno6` (6, limited) | $3.01\cdot10^{-2}$ | $7.75\cdot10^{-3}$ | 0.178 | $1.48\cdot10^{-2}$ | $1.52\cdot10^{-4}$ |
| `quintic-mp` (6, limited) | $3.14\cdot10^{-2}$ | $8.82\cdot10^{-3}$ | 0.173 | $2.47\cdot10^{-2}$ | $6.25\cdot10^{-4}$ |

Findings:
- **The defect has no trend with the fill's order, and limiting makes it worse.** The best fill differs by case and
  by metric (`restriction-compatible` on Sod, `tricubic` on the RJ2a plateau, `quintic` on its band), every fill
  sits within roughly $\pm 30\,\%$ of `tricubic` but the limited ones, and the two limited fills are the worst on
  Sod, their shock furthest off: the plateau defect behind a crossing shock is not an overshoot of the fill, so no
  fill can be held to "no worse than `tricubic`". The best Sod result belongs to the one fill with exact $R \circ P =
  I$: the defect points to the consistency of the fill with the restriction at the crossing, and to the reflux (A3:
  reflux on/off moves it by up to 30 %), which is Phase E's ground, not the fill's order.
- **On smooth flows both limited fills keep their order.** AO `euler-ww` (`weno6` restriction and fill, reflux off):
  L1 orders 5.32 / 5.50, L1 $8.60\cdot10^{-8}$ at N = 256; `euler-mp` (`lagrange6` + `quintic-mp`): 6.10 / 6.33,
  $6.04\cdot10^{-8}$, `euler-q6`'s (the bound never clips the vortex). PV-5 passes with no inadmissible stage under
  both; GP holds a cubic field to round-off with either fill.

### E0: the shock-crossing defect is a capturing error

**Why.** A3 found an O(1), non-converging defect behind a shock that crosses a 2:1 seam; D2 showed that no
coarse-to-fine fill reduces it and that limiting the fill makes it worse. E0
([#68](https://github.com/szaghi/adam/issues/68)) asks whether it comes from the consistency of the fill with the
restriction (only `restriction-compatible` with the `mean` restriction has $R \circ P = I$ exactly), from the reflux,
or from neither.

**Runs.** The AS c2f crossings (FNL, 2 ranks), every combination of fill (`restriction-compatible`, `tricubic`),
restriction (`mean`, `lagrange6`) and reflux (on, off); Sod also at N = 400. Plateau density defect behind the shock:

| Fill + restriction | Sod, reflux on | Sod, reflux off | RJ2a, reflux on | RJ2a, reflux off |
|---|---|---|---|---|
| compatible + mean ($R \circ P = I$) | $2.32\cdot10^{-2}$ | $1.82\cdot10^{-2}$ | $1.97\cdot10^{-2}$ | $1.83\cdot10^{-2}$ |
| compatible + lagrange6 | $2.32\cdot10^{-2}$ | $1.85\cdot10^{-2}$ | $1.98\cdot10^{-2}$ | $1.84\cdot10^{-2}$ |
| tricubic + mean | $2.46\cdot10^{-2}$ | $3.10\cdot10^{-2}$ | $1.41\cdot10^{-2}$ | $1.68\cdot10^{-2}$ |
| tricubic + lagrange6 | $2.46\cdot10^{-2}$ | $3.17\cdot10^{-2}$ | $1.41\cdot10^{-2}$ | $1.68\cdot10^{-2}$ |

Sod with the reflux on at N = 400: $2.35\cdot10^{-2}$ (tricubic) and $2.14\cdot10^{-2}$ (compatible), while the band
defect next to the seam converges ($6.35\cdot10^{-3} \to 1.71\cdot10^{-3}$, order 1.9). With the reflux off the mass
drifts ($10^{-5}$ to $10^{-4}$) and the shock lands 0.06 to 0.22 cells off.

Findings:
- **The restriction is irrelevant** (`mean` and `lagrange6` within 0.5 % in every case).
- **Consistency does not control it**: the $R \circ P = I$ pair is best on Sod and worse than `tricubic` on RJ2a.
- **The reflux does not cause it**: off, it moves the defect by up to 30 % either way and costs conservation.
- **It does not converge**: N = 200 to 400 lowers it by 5 to 8 %.

So it is not a seam-operator defect. The evidence points to the start-up error of a captured shock whose numerical
structure, a few cells wide, is squeezed into finer cells (an interpretation, not measured directly): the transition emits a spurious wave (the plateau entropy defect is about 1.8 %)
whose amplitude scales with the shock's strength, not with the spacing. The remedy is to keep shocks off the seams,
refining ahead of them (a gradient or Löhner marker with a buffer and a short regrid frequency, M5), which
[Numerics](./numerics#adaptive-mesh-refinement) now recommends; Phase E (the reflux policy) proceeds on its own merit.

### E2: the reflux policy `troubled` against `always`

**Why.** The Berger–Colella reflux is what keeps a 2:1 seam conservative, and what caps the composite at second order
(A1, B1): it replaces the coarse face flux by the restricted fine one everywhere. `[numerics] reflux_policy = troubled`
([#68](https://github.com/szaghi/adam/issues/68) E1) applies it only at the seam face cells flagged in a stage of the
step: those the positivity limiter blended, and those where the WENO-5 sensor ($\min_k w_k/d_k$ of density, normal
velocity or pressure across the face) is below 0.2. E2 measures it against `always` on the four checks of the plan
(FNL, 2 ranks; the FNL flags equal the CPU's to the cell).

| Check | `always` | `troubled` |
|---|---|---|
| AO euler, default operators, N = 256: L1, order | $4.05\cdot10^{-6}$, 2.01 | $1.80\cdot10^{-6}$, 2.34 |
| AO euler, `lagrange6` + `quintic`, N = 256: L1, order | $4.05\cdot10^{-6}$, 2.01 | $6.03\cdot10^{-8}$, **6.11** |
| AO mhd, default operators, N = 256: L1, order | $1.83\cdot10^{-6}$, 2.07 | $1.10\cdot10^{-6}$, 2.04 |
| AO face cells flagged, N = 64 / 128 / 256 | — | 3–5 % / 0–0.6 % / 0 % |
| A3 shock position (Sod, RJ2a, c2f and f2c) | reference | equal within $2\cdot10^{-4}$ cells |
| A3 plateau and band defects | reference | within 1 % (RJ2a c2f band +44 %) |
| A3 mass drift | $\le 3\cdot10^{-10}$ | $10^{-11}$ to $1.5\cdot10^{-7}$ |
| A4 drift, vortex N = 64 to $t = 1$ (r, ru, rv, rE) | $\le 3\cdot10^{-16}$ | $1.8$–$4.4\cdot10^{-5}$, growth 1.3–1.6 (rv 8.8) |
| A4 drift, Sod c2f | $3\cdot10^{-16}$ | $1.9\cdot10^{-7}$ |
| CPAW D_flux seam floor, N = 64 / 128 (growth) | $5.39\cdot10^{-4}$ / $5.43\cdot10^{-5}$ (0.95 / 0.72) | $5.52\cdot10^{-4}$ / $3.80\cdot10^{-5}$ (0.76 / 0.68) |

Findings:
- **The reflux was what capped the high-order seam operators.** With `lagrange6` and `quintic` the reflux on
  (`always`) holds the vortex at second order, $4.05\cdot10^{-6}$; under `troubled` no smooth face is flagged once
  resolved and the run is sixth order, $6.03\cdot10^{-8}$, the error of the uniform coarse run (67 times smaller).
  AO leg `euler-tq` asserts it.
- **Shocks are unchanged**: positions equal to $2\cdot10^{-4}$ cells, defects within 1 % but one band value; the
  flags cover 35 to 42 % of the face-cell steps of a crossing.
- **Conservation is the price**: an unflagged face keeps its flux mismatch, so the integrals drift at the reflux-off
  rate on smooth flows ($10^{-5}$ per period at N = 64, bounded but one component accumulating) and by $10^{-7}$
  across a shock.
- **MHD**: the D_flux seam floor is bounded and 1.4 times lower at N = 128 (83 against 119 times the uniform floor).

**Decision (D-M6-2, owner):** the default stays `always`, exactly conservative; `troubled` is the opt-in that
makes the high-order seam operators worth their cost with the reflux on.

### G3a: #78 is the reflux, not the restricted ghosts

**Why.** A2 traced EGLM's 70-fold $\langle|B_z|\rangle$ on a static 2:1 grid
([#78](https://github.com/szaghi/adam/issues/78)) to the centred divergence of its Derigs source reading the
second-order restricted ghosts. Phases C to E gave the seams sixth-order ghosts and a reflux policy, so G3a
([#68](https://github.com/szaghi/adam/issues/68)) asks first whether #78 survives them.

**Runs.** The A2 field loop (refined N = 64 against uniform N = 128, $t = 1$, FNL, 2 ranks), GLM and EGLM, under
the default seam operators (`mean` + `tricubic`), under `lagrange6` + `quintic`, each with the reflux on (`always`),
with `troubled` (which flags no face of this smooth loop, so it equals the reflux off bit for bit) and off:

| Seam operators, reflux | EGLM $\langle\lvert B_z\rvert\rangle/A_0$ (× uniform) | GLM (× uniform) |
|---|---|---|
| uniform N = 128 | $2.72\cdot10^{-6}$ (1) | $6.86\cdot10^{-4}$ (1) |
| default, `always` | $1.91\cdot10^{-4}$ (70) | $1.10\cdot10^{-3}$ (1.6) |
| default, off (= `troubled`) | $9.31\cdot10^{-6}$ (**3.4**) | $1.18\cdot10^{-3}$ (1.7) |
| `lagrange6` + `quintic`, `always` | $2.89\cdot10^{-4}$ (106) | $1.12\cdot10^{-3}$ (1.6) |
| `lagrange6` + `quintic`, off (= `troubled`) | $3.14\cdot10^{-5}$ (11.5) | $1.56\cdot10^{-3}$ (2.3, run with `troubled`) |

Findings:
- **The reflux causes #78.** With the default operators and no reflux, EGLM's $B_z$ is 3.4 times its uniform level
  instead of 70; higher-order ghosts do not help (106 with the reflux, 11.5 without). GLM does not care (1.6 to 2.3
  times in every row). The A2 attribution to the restricted ghosts is wrong.
- **Mechanism (inferred, not measured directly):** EGLM's Derigs source $-\nabla\cdot\mathbf B\,(0, \mathbf B,
  \mathbf u\cdot\mathbf B, \mathbf u, \ldots)$ cancels, on a uniform grid, the $B_z$ that the flux divergence
  generates. The reflux corrects a coarse skin cell's state by the flux mismatch $F_c - \sum F_f$ after the step,
  but the source that step used saw the coarse fluxes' divergence: the correction changes the discrete divergence the
  fluxes imply without the matching source term, so the cancellation fails exactly at the refluxed cells.
- **Workaround before G3b:** EGLM with `reflux_policy = troubled` keeps #78 at 3.4 times on smooth fields (the reflux is
  applied only where a shock or the limiter flags a face), at the conservation price measured in E2.
- **The fix:** G3b, below: the EGLM reflux leaves B alone (3.3 times).

### G2: the GLM damping length on a source-free seam test (study)

**Why.** The GLM damping rate is $c_h^2/c_p^2 = \alpha\, c_h / L$ with `[mhd] glm_damping_length` $L$ a fixed
length or `min-cell` (the finest spacing; with runtime AMR it is frozen at the spacing of `max_level`, D-M5-9). A
`per-level` length, each level damping at its own spacing, does not exist; G2
([#68](https://github.com/szaghi/adam/issues/68)) asks whether it would matter, bracketing it between the coarse and
the finest spacings on A2's source-free CPAW (GLM, two $x$ seams, 10 periods, FNL, 2 ranks):

| `glm_damping_length` | rate $c_h^2/c_p^2$ (N = 128) | seam D_flux floor N = 64 / 128 | × uniform (N = 128) | growth N = 128 |
|---|---|---|---|---|
| 1.0 (A2) | 0.36 | $5.39\cdot10^{-4}$ / $5.43\cdot10^{-5}$ | 119 | 0.72 |
| coarse spacing $L_x/N$ | 20.6 | $6.96\cdot10^{-4}$ / $5.32\cdot10^{-5}$ | 116 | 1.02 |
| `min-cell` (finest spacing) | 41.2 | $8.30\cdot10^{-4}$ / $1.07\cdot10^{-4}$ | 233 | 1.01 |

Findings:
- **Every length keeps the seam divergence bounded** (growth 0.72 to 1.02 over 10 periods).
- **Strong damping raises the seam floor**: `min-cell` doubles it at N = 128 against the coarse spacing or the
  domain-scale length; a per-level length (between the coarse and the finest spacings) would sit between them.
- **The D-M5-9 freeze stays** (the plan's choice): its cost on this test is at most that factor 2 on the floor,
  with the damping rate fixed through regrids. No `per-level` option is added.

### G3b: the EGLM reflux leaves B alone

**Why.** G3a showed that [#78](https://github.com/szaghi/adam/issues/78), EGLM's $\langle|B_z|\rangle$ 70 times its
uniform level on a static 2:1 grid, is caused by the reflux. G3b ([#68](https://github.com/szaghi/adam/issues/68))
looked for the fix.

**A first attempt that failed: the source divergence from the face fluxes.** Computing the Derigs source's $\nabla\cdot
\mathbf B$ from the $\psi$ face fluxes ($\sum_d (F_\psi(i+\tfrac12) - F_\psi(i-\tfrac12)) / (c_h \Delta x_d)$, the
divergence of the interface field, with the matching source correction at the reflux) keeps the MV-5 orders (4.96 to
4.97 on every wave) but raises EGLM's uniform $\langle|B_z|\rangle$ from $2.72\cdot10^{-6}$ to $1.28\cdot10^{-4}$ (47
times; refined N = 64 $2.64\cdot10^{-4}$): EGLM's cancellation of the $B_z$ the flux splitting generates needs the
centred derivative of the cell values, not the face divergence. The code was not kept.

**The fix: no reflux of B under EGLM.** The field loop (refined N = 64, default seam operators, reflux on) with the
reflux of B removed, the rest refluxed:

| EGLM, reflux | $\langle\lvert B_z\rvert\rangle/A_0$ | × uniform |
|---|---|---|
| all components (before) | $1.91\cdot10^{-4}$ | 70 |
| all but B (G3b) | $9.11\cdot10^{-6}$ | **3.3** |
| all but B and $\psi$ | $9.31\cdot10^{-6}$ | 3.4 |
| none | $9.31\cdot10^{-6}$ | 3.4 |

Refluxing B is the whole of #78. Under EGLM the induction equation carries the Derigs source $-\mathbf u\,\nabla\cdot
\mathbf B$, so $\int \mathbf B$ is not a conserved integral: refluxing B buys a conservation the model does not have
and breaks the source/flux cancellation it is built on. Since G3b (owner decision, no key) the EGLM reflux corrects
$\rho$, $\rho\mathbf u$, $E$ and $\psi$ and leaves B alone, on CPU and FNL; GLM is unchanged (it refluxes B, and its
$\langle|B_z|\rangle$ does not depend on it, G3a).

**Gate.** The FLUME regressions pass unchanged on CPU and FNL (35 of 35 each). The one case with EGLM across 2:1 seams,
`blast-amr-limiter`, moves but stays inside its digest tolerance (rtol $10^{-6}$), so its golden is kept:
- step 0 is identical, and $B_z$ (zero in 2D) does not move;
- at step 100 seven of 26 digest rows differ (`r`, `ru`, `rv`, `rE`, `bx`, `by`, `psi`), the largest by
  $2.4\cdot10^{-7}$ relative on CPU and $9.8\cdot10^{-8}$ on FNL;
- the integrals drift exactly as in the golden: $\rho$ to round-off ($3\cdot10^{-16}$), the rest at the golden's
  levels ($\rho\mathbf u$ $2.5\cdot10^{-4}$, $E$ $1.9\cdot10^{-2}$, B $4.5\cdot10^{-4}$, $\psi$ $1.3\cdot10^{-4}$, all
  set by the blast's boundaries and EGLM's non-conservative sources, not by the reflux).

The field loop on FNL confirms the table above: EGLM $9.11\cdot10^{-6}$, GLM $1.10\cdot10^{-3}$ (unchanged). MV-5
under EGLM keeps its orders (4.96 to 4.97 on every wave), and the multirealm suite passes.


### AL: integral drift and the WENO sensor near seams

**Why.** The reflux policy of Phase E (`reflux = troubled`) would turn the reflux off at the seam faces that a WENO
sensor judges smooth. Two quantities decide whether that is safe and where its threshold sits, and no check logged
either: how the conserved integrals drift when the reflux is off, and what the sensor reads on the seam faces of a
smooth flow (false positives) and of a shock crossing (detections).

**Cases** (`verification/amr-log/`, logs only). `vortex`: the AO Euler vortex, quadtree, quadrant $[0.5, 1]^2$
refined, N = 64, run to $t = 1$, one period of the free stream $(1, 1)$, so the vortex crosses every seam; the domain
is periodic and every integral is conserved. `sod`: the AS Sod c2f, N = 200, $t = 0.2$; mass and energy are conserved.
Each leg runs with the reflux on and off, with a checkpoint every 10 steps.

**Logs** (`amr_log.py`):
- drift: the relative drift of the integrals of the conservation history, its maximum, and its growth (the mean
  $|$drift$|$ over the last quarter of the run over that of the second quarter: about 1 when bounded, $7/3$ for a
  linear drift);
- sensor: on every checkpoint, $r = \min_k w_k/d_k$ of the WENO-5 weights of every face, the minimum over $\rho$,
  $u_n$, $p$ and over the two reconstructions beside the face (the weno-riemann sensor, read on primitive variables),
  with the weights of `adam_weno_object` (Jiang–Shu, $\epsilon = 10^{-6}$, exponent 3), recomputed from the
  checkpoints, whose ghosts are filled right before the save as the stencils read them. Per region (the seam band, the
  interior): the minimum $r$ and the fraction of faces below $\tau = 0.2$.

The default scheme weights its characteristic split fluxes, so the map is the sensor a seam policy would read, not a
dump of the scheme's own weights. Asserted: with the reflux on, every integral drifts at most $10^{-13}$.

**Results** (2 ranks):

| Leg | Reflux on: max drift | Reflux off: max drift (mass, energy) | Reflux off: growth |
|---|---|---|---|
| vortex | $\le 4.5\cdot10^{-16}$ | $2.2\cdot10^{-5}$, $2.5\cdot10^{-5}$ (momentum $3.9$, $2.7\cdot10^{-5}$) | 1.4 to 7.7 |
| sod | $\le 3.3\cdot10^{-16}$ | $1.2\cdot10^{-4}$, $2.4\cdot10^{-4}$ | a jump at the crossing |

| Leg | Seam band: $r < \tau$, largest (mean) fraction | Interior: largest (mean) fraction |
|---|---|---|
| vortex | 0.110 (0.053) | 0.051 (0.020) |
| sod | 1.000 (0.172) | 0.128 (0.078) |

Findings:
- **Without the reflux, a smooth flow's drift accumulates.** Every seam crossing adds to it: by one period the vortex
  has lost about $2\cdot10^{-5}$ of its mass and energy, and the growth (1.4 to 7.7) says the drift is not bounded. A
  shock adds $10^{-4}$ in one crossing. The reflux keeps both at round-off.
- **The sensor detects the crossing shock.** In the Sod run the seam band reads $r = 1$ until the shock reaches it,
  67 to 100 % of its faces are below $\tau$ while the shock crosses ($t = 0.126$ to 0.152), and 0 to 13 % after,
  where the A3 defect sits.
- **It also flags smooth seam faces.** On the vortex, up to 11 % of the band faces (mean 5 %) read $r < \tau$, 2.6
  times the interior's mean fraction: the JS weights drop at the extrema of a smooth profile, and more so where the
  seam's second-order error adds kinks. For `reflux = troubled` a false positive keeps the reflux, the conservative
  choice: the cost is the order on those faces, which Phase B2 measures. The threshold is not tuned here.

### B1: the 2:1 caps in a 2-D prototype

**Why.** #68 Appendix A measured in 1-D that the mean restriction and the reflux of point-value fluxes each cap the
order at 2, and that a 6-point restriction with a quintic prolongation and no reflux recovers order 5. A 1-D model has
no tangential direction: no tangential averaging in the register, no tangential footprint of the prolongation, no seam
corners. Phase B repeats the measurement in 2-D, in NumPy, without ADAM runs.

**Model** (`verification/amr-proto/proto2d.py b1`). $u_t + u_x + u_y = 0$ on the periodic unit square, initial
$\sin 2\pi x \sin 2\pi y$, one period. Conservative finite differences of point values with the 5th-order linear upwind
flux (the WENO5 optimal weights), RK4 with $\Delta t \propto h^{5/4}$ so that the spatial error dominates. The box
$[0.25, 0.75]^2$ is refined 2:1; the coarse cells it covers hold the restriction of the fine ones, as FLUME's coarse
ghosts do. Arms:
- restriction: `mean` (the 2×2 mean, FLUME's) or `lag6` (6×6 tensor Lagrange, shifted inside the box at its edge);
- prolongation: `cubic` (FLUME's tricubic) or `quintic`, from coarse cells outside the box only (shifted away from the
  box along the normal axis, as FLUME's donor-block footprint);
- reflux: on (a coarse box face takes the mean of its two fine fluxes) or off.

The maximum error is also split by region: within 3 coarse cells of a box corner, of a box edge, or neither.

**Results** (orders of N = 128 → 256; mass drift after one period at N = 256):

| Arm | L1 | L∞ | L∞ corner / edge / interior | Mass drift |
|---|---|---|---|---|
| uniform (coarse) | 5.01 | 5.01 | – / – / 5.01 | $5\cdot10^{-18}$ |
| mean, cubic, on (today) | 2.85 | 1.93 | 1.93 / 1.97 / 1.97 | $1\cdot10^{-18}$ |
| mean, cubic, off | 2.01 | 1.95 | 1.95 / 1.97 / 1.97 | $5\cdot10^{-8}$ |
| lag6, cubic, on | 2.85 | 1.93 | 1.93 / 1.97 / 1.95 | $2\cdot10^{-18}$ |
| lag6, cubic, off | 3.98 | 4.01 | 4.01 / 3.96 / 3.95 | $1.3\cdot10^{-9}$ |
| lag6, quintic, on | 2.85 | 1.93 | 1.93 / 1.95 / 1.96 | $6\cdot10^{-19}$ |
| lag6, quintic, off | **4.99** | **4.99** | 4.99 / 4.98 / 4.98 | $9.1\cdot10^{-11}$ |
| mean, quintic, off | 2.02 | 1.93 | 1.93 / 1.98 / 1.97 | $7\cdot10^{-8}$ |

Findings:
- **The 1-D findings hold in 2-D.** The mean restriction caps the order at 2 whatever the prolongation; the reflux caps
  it at 2 whatever the transfers (its tangential mean of point fluxes included). The 6-point restriction with the
  quintic prolongation and no reflux recovers order 5 in L1 and L∞, with the tangential footprint and the corners
  present; with the cubic prolongation, order 4.
- **Without the reflux, the mass drift converges at the scheme's order** when the transfers are high-order
  ($9\cdot10^{-11}$ at N = 256, against $5\cdot10^{-8}$ with the mean restriction).
- **The corners hold the maximum error, not a lower order.** In every arm the maximum sits within 3 cells of a box
  corner, but the corner, edge and interior maxima converge at the same rate. In this linear model the corners do not
  explain the L∞ order 1.46 of A1's Euler vortex; what the model lacks is the nonlinear WENO weights and the
  characteristic projection (#68 F6), which B2's nonlinear prototype brings in.

`check.sh` (1 minute; `--long` for N = 256) asserts lower bounds only: the uniform control at 4.8, lag6-quintic-off at
4.7 (L1 and L∞), lag6-cubic-off at 3.7 (L∞).

### B2: a sensor-gated reflux, the go/no-go for Phase E

**Why.** B1 recovers order 5 only without the reflux, and A3 showed that a shock needs it. Phase E would apply the
reflux only where a WENO sensor flags the face (`reflux = troubled`); D-M6-4 builds it only if a prototype shows the
order recovered with the shock positions intact.

**Model** (`proto2d.py b2`). The B1 composite grid with WENO5-JS (Lax–Friedrichs splitting, the weights of
`adam_weno_object`: $\epsilon = 10^{-6}$, exponent 3) and the policies `always` (today), `never` and `troubled`: a
coarse box face takes the mean of its fine fluxes only where the sensor $r = \min_k w_k/d_k$ of the coarse face or of
either fine face is below $\tau = 0.2$. Arms: today's mean/cubic/always, and lag6/quintic under each policy. Legs:
- `smooth`: the B1 advection; orders, and the share of seam-face evaluations the sensor flags;
- `c2f`, `f2c`: Burgers pulses, $u = 1$ on $[0, 0.15]$ and on $[0.45, 0.6]$, whose shocks cross $x = 0.25$ (coarse to
  fine, $t = 0.2$) and $x = 0.75$ (fine to coarse, $t = 0.3$), run to $t = 0.5$ and 0.6, the rarefaction catching the
  shock on the way;
- `weak`: a $10^{-3}$ pulse riding on $u = 1$, crossing $x = 0.25$. Its smoothness indicators sit near $\epsilon$, so
  the JS sensor misses it: the missed-detection control.

The pulses start from exact cell averages, so the discrete mass is the exact one (sampling the step at the centres
offsets it by up to a cell, which made the first version of the shock metric erratic). Against the exact solution: the
mass defect over the pulse's own mass, and the mass-equivalent shock position of every row, in local cells.

**Results.** Smooth leg, orders of N = 64 → 128:

| Arm | L1 | L∞ | L∞ corner / edge | Faces flagged |
|---|---|---|---|---|
| uniform | 5.01 | 4.95 | – | – |
| mean, cubic, always (today) | 3.03 | 2.22 | 2.22 / 2.05 | – |
| lag6, quintic, always | 2.96 | 2.10 | 2.10 / 1.96 | – |
| lag6, quintic, troubled | **5.01** | **4.95** | 5.10 / 5.02 | 0 % |
| lag6, quintic, never | 5.01 | 4.95 | 5.10 / 5.02 | – |

Shock legs (lag6/quintic; mean shock position in local cells, mass defect over the pulse's mass, the uniform control
in brackets):

| Case, N | always | troubled | never | Faces flagged |
|---|---|---|---|---|
| c2f, 64 | −0.53, $0$ | −0.53, $5.0\cdot10^{-5}$ | −1.12, $1.8\cdot10^{-2}$ | 15 % |
| c2f, 128 (−0.30) | −0.50, $2\cdot10^{-16}$ | −0.50, $3.6\cdot10^{-5}$ | −1.07, $8.6\cdot10^{-3}$ | 8 % |
| f2c, 64 | −0.25, $2\cdot10^{-16}$ | −0.25, $6.5\cdot10^{-5}$ | +0.05, $1.7\cdot10^{-2}$ | 15 % |
| f2c, 128 (−0.37) | −0.25, $4\cdot10^{-16}$ | −0.25, $6.9\cdot10^{-6}$ | +0.06, $8.4\cdot10^{-3}$ | 7 % |
| weak, 64 (+0.02) | −0.003, $10^{-12}$ | −0.014, $3.9\cdot10^{-4}$ | −0.014, $3.9\cdot10^{-4}$ | 0 % |
| weak, 128 (+0.01) | −0.015, $10^{-12}$ | −0.022, $3.2\cdot10^{-5}$ | −0.022, $3.2\cdot10^{-5}$ | 0 % |

`check.sh` asserts: the smooth L1 and L∞ orders of `troubled` at least 4.7; on c2f and f2c, `troubled` within 0.05
cells of `always` and its mass defect at most $10^{-3}$, while `never`'s is at least $10^{-3}$ (the negative control:
the metrics see a missing reflux).

Findings:
- **`troubled` recovers the order.** On the smooth leg the sensor flags no seam face, so `troubled` is `never`: order
  5 in L1 and L∞, as the uniform control, against about 2 with the reflux always on.
- **`troubled` keeps the shock where `always` puts it.** The two differ by at most 0.003 cells; `never` moves the shock
  by 0.3 to 0.6 cells and loses $10^{-2}$ of the pulse's mass. The reflux at the flagged faces (7 to 15 % of the
  seam-face evaluations) is what holds the shock.
- **`troubled` is near-conservative, not conservative.** Its mass defect, $7\cdot10^{-6}$ to $6.5\cdot10^{-5}$ of the
  pulse, comes from the unflagged smooth faces (the rarefaction crossing the seam) and from the steps before the shock
  reaches the sensor's stencil; it is two to three orders below `never`'s.
- **A missed detection is a small error here.** The weak shock is never flagged, so `troubled` runs as `never`: its
  shock moves by 0.007 to 0.011 cells against `always` and it loses $3\cdot10^{-5}$ to $4\cdot10^{-4}$ of its own
  mass. What the JS sensor misses is what the reflux matters least for, but not nothing.
- **The corners are again not the L∞ cause.** Under the nonlinear weights corner and edge maxima still converge at the
  same rate; A1's Euler L∞ order 1.46 is not reproduced by a scalar model, and stays open for the system (the
  characteristic projection).

**Verdict (D-M6-4): go for Phase E**, with three conditions on record. The order is recovered only together with the
high-order transfers of Phases C and D: with today's mean restriction the order stays at 2 whatever the reflux does
(B1). The policy is
near-conservative ($10^{-5}$ of the shocked mass), not exact. And A4 measured 5 % false positives on the Euler vortex
with the sensor read on primitive variables, where this scalar model flags none; E must measure its own false-positive
rate, since each flagged smooth face brings the order cap back.

### B3: the footprint of a high-order restriction

**Why.** Phase C replaces the mean restriction with a 6-point one. Its stencil reaches 2 to 3 fine cells beyond the two
cells of the coarse cell, so near the edge of a fine block it needs data the block does not hold. Either the stencil
shifts inside the donor block (shift-inward: one-sided, at least 6 fine cells per block), or the fine block's ghosts
are filled first and the stencil stays centred (two-pass: any block size, at the price of an ordering constraint).

**Model** (`proto2d.py b3`). The B1 problem with the reflux off, where the restriction's order shows. The fine box is
split into blocks of B fine cells per axis: `lag6` shifts inside its donor block, `lag6c` reads the fine ghosts (exact
fine neighbours inside the box, the prolongation at its edge). B = 4, 6, 8, 16 and the whole box; both prolongations;
N = 48 / 96 / 192 (so that 6, 8 and 16 divide the fine cells).

**Results** (N = 192; orders of 96 → 192):

| Footprint, prolongation | Block size | L1 | L∞ | Order L1 / L∞ |
|---|---|---|---|---|
| shift-inward, cubic | 6 / 8 / 16 / box | $9.2$ to $9.4\cdot10^{-8}$ | $3.7$ to $3.9\cdot10^{-7}$ | 3.98 / 4.00 |
| two-pass, cubic | any | $2.4\cdot10^{-8}$ | $1.9\cdot10^{-7}$ | 4.14 / 4.03 |
| shift-inward, quintic | 6 / 8 / 16 / box | $2.38\cdot10^{-9}$ | $7.83\cdot10^{-9}$ | 4.98 / 4.99 |
| two-pass, quintic | any | $2.42\cdot10^{-9}$ | $7.85\cdot10^{-9}$ | 5.00 / 5.00 |
| shift-inward | 4 | infeasible | | |

Findings:
- **The block size does not matter for shift-inward.** From 6 fine cells per block to the whole box the error changes
  by at most 2.5 % in L1 and 5 % in L∞, at the same order.
- **With the quintic prolongation the two footprints are equivalent** (within 1.5 %), both at order 5.
- **With today's cubic prolongation, two-pass is 4 times more accurate in L1** (2 times in L∞), at the same order 4:
  the one-sided shift has a larger error constant than the centred stencil reading 4th-order ghosts.
- **What each costs in FLUME.** `update_ghost_local` fills copies, coarse-to-fine interpolations and fine-to-coarse
  restrictions in one pass, in map order, the restriction reading only the donor block's interior cells. Shift-inward
  keeps that single pass and needs 6 fine cells per refined axis, which the #66 guard already requires
  ($n \ge 2\,n_{gc} = 6$). Two-pass needs the fine ghosts complete before any restriction: a second local pass and, on
  several ranks, a second exchange round, because the sender must update its fine ghosts before it packs.

**Recommendation for Phase C: shift-inward.** It fits the current ghost update and block guard, and with the quintic
prolongation of Phase D it is as accurate as two-pass. Two-pass would buy a constant factor only while the cubic
prolongation stays, at the cost of a second synchronisation per ghost update.

`check.sh` asserts the L∞ orders of shift-inward (B = 6) and two-pass: 3.7 with the cubic prolongation, 4.7 with the
quintic (N = 24 / 48 / 96, about 20 s).

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

## Gravity (M7)

External gravity lands in phases ([#73](https://github.com/szaghi/adam/issues/73), plan in the issue): G0 the inputs,
G1 the force, G2 a well-balanced mode, G3 Rayleigh–Taylor. Self-gravity is M8
([#81](https://github.com/szaghi/adam/issues/81)). The figures of this section are drawn from re-runs of the
`check.sh` inputs (CPU, 2 ranks), whose checkpoints the checks delete: `make_doc_runs.sh <dir> gravity`, then
`make_doc_figures.py --out docs/public/flume --runs <dir> --only gravity-atmospheres gravity-force gravity-order
gravity-balance`.

### GV (G0): the stratified atmospheres and the `hydrostatic` boundary

**Why.** Every later gravity test starts from a hydrostatic atmosphere closed by `hydrostatic` faces: G0 checks that
the initial state and its ghosts are what the formulas say before any force acts on them.

**Cases** (`gravity/check.sh`, inputs from `gravity/make_atmosphere.py`, checks by `gravity/atmosphere_oracle.py`):
each runs one negligible step (CFL = 10⁻³⁰) with the ghosts poisoned, so the checkpoint holds the initial atmosphere
and its ghosts filled in the order of a Runge–Kutta stage.

| Case | Atmosphere | Gravity | Faces | Grid |
|---|---|---|---|---|
| iso2d | isothermal | (0, −1) | x periodic, y hydrostatic | quadtree, 4×4 blocks of 8² |
| poly2d | polytropic, $G = 1.2$ | (0, −1) | as iso2d | as iso2d |
| tilted2d | isothermal | (−0.6, −0.8) | every face hydrostatic | as iso2d, one block refined (2:1 seams) |
| tilted3d | isothermal | (−0.48, −0.64, −0.6) | every face hydrostatic | octree, 4×4×4 blocks of 8³ |
| iso2d-mhd | isothermal, EGLM, B = (0.5, 0.7, 0.2) | (0, −1) | as iso2d | as iso2d |
| froude2d | isothermal | `froude = 1`, `axis = -y` | as iso2d | as iso2d |

**Asserted** (relative, 10⁻¹³):
- the interior cells against the analytic atmosphere;
- every ghost beyond a `hydrostatic` face against the boundary rule. For the isothermal atmosphere the rule is the
  analytic profile, so the realm edges and corners are checked too; for the polytropic one, the ghosts beyond one face.
  The velocity in those ghosts must be zero and the field uniform;
- the periodic ghosts against the analytic profile at the wrapped position.

Boundary ghosts beside their own block (their donor along the normal is an interface or seam ghost, not exact on an
exponential) are counted, not checked; the exchange is the ghost probe's job (GP).

A **negative control** runs iso2d with `extrapolation` on the y faces and judges it as if they were `hydrostatic`: the
oracle must reject it. **NV-5** runs `scaling/reference.sh --inputs` on iso2d, iso2d-mhd and froude2d: the gravity
components convert as an acceleration ($u_0^2/L_0$), `froude` does not, and the dimensional run equals the base one bit
for bit.

Two uniform levels spread the blocks over the ranks (with one level every block stayed on rank 0, so an np 2 run
decomposed nothing).

**Results** (CPU and FNL, np 1 and 2, every case PASS):

| Check | Cells (per run) | Max relative error |
|---|---|---|
| interior ρ, p | 1 024 (2-D), 32 768 (3-D) | 2.2·10⁻¹⁶ CPU, 6.7·10⁻¹⁶ FNL |
| hydrostatic ghosts ρ, p | 192 (iso2d) to 22 104 (tilted3d, edges and corners) | 8.9·10⁻¹⁶ |
| hydrostatic ghosts u, B | as above | velocity 0 to 10⁻³², field uniform |
| negative control | 192 | 9.8·10⁻² (rejected, as it must be) |
| NV-5 | iso2d, iso2d-mhd, froude2d | bitwise; conversions exact and complete |

![GV (G0): the iso2d and poly2d atmospheres with their hydrostatic ghosts, and their relative error](/flume/gravity-atmospheres.png)

Left, the two atmospheres along y, every interior cell (dots) and every ghost beyond the y faces (crosses, shaded),
against the analytic profiles. Beyond a polytropic face the ghosts leave the polytropic curve on purpose: the
`hydrostatic` rule extends the boundary cell isothermally, and they are checked against that extension. Right, the
relative error of every one of those cells: zero or one unit in the last place, ghosts included.

**Found on the way** (pre-existing, not gravity): on FNL an MHD run died when a rank owns no blocks
(`dev_assign_to_device: device memory not allocated`), with or without gravity and with the executable before G0; the
Euler model and the CPU backend ran. The G0 cases give every rank blocks. Fixed and asserted since, see
[ranks without blocks](#ranks-without-blocks-dv-1) ([#82](https://github.com/szaghi/adam/issues/82)).

### GV (G1): the force

**Why.** G1 adds the gravity sources, momentum $\rho\mathbf g$ and the energy from the face mass fluxes. Three
properties must hold by construction and are checked to round-off: free fall is exact, a closed box conserves its mass
and $E + \rho\Phi$, and the scaling covariance of #49 survives. The positivity limiter must keep the run alive where the
gas falls toward vacuum.

**Cases** (`gravity/check.sh`, np 2), a uniform $\rho = 1$, $p = 1$ at rest under $\mathbf g = (0, -1)$, CFL 0.4:

| Leg | Case | Asserted |
|---|---|---|
| GV-1 | `freefall2d`: every face periodic, 50 steps | $\rho\mathbf u = \rho\mathbf g t$, $E = E_0 + \rho g^2 t^2/2$, $\rho$ unchanged, 10⁻¹² |
| GV-3 | `box2d`: x periodic, y `wall-inviscid` (the gas falls and sloshes), 200 steps; `box2d-mhd` (EGLM, $B_x = 0.5$, limiter on) | mass and the total energy constant at every checkpoint, 10⁻¹²: $E + \rho\Phi$ ($\Phi = -\mathbf g\cdot\mathbf x$) in G1, $W$ (its wall-row term added) since G2a; $E$ alone reported |
| limiter | `box2d-lim` ($p = 0.05$, 300 steps: the gas falls away from the lid toward vacuum); `box2d-stiff` ($p = 10^{-3}$, $g\Delta t \approx 4.5\,c$) | the runs complete, the mass is conserved; $E + \rho\Phi$ reported (exact only where the limiter does not act) |
| GV-2 | `iso2d` with the force, t = 0.5, N = 32, 64, 128 | reported: the spurious $\max|\mathbf u|/c$ of the unbalanced scheme |
| GV-6 | `box2d` with `weights = si` against its input rescaled by $2^2$, $2^{-1}$, $4^{-2}$ | bitwise |
| GV-7 | `box2d`, `box2d-mhd` on FNL against CPU | $10^{-10}$ |

**Results** (CPU; FNL the same, GV-7):

| Leg | Measured |
|---|---|
| GV-1 | $\rho$ exact, $\rho\mathbf u$ 2.2·10⁻¹⁶, $E$ 5.6·10⁻¹⁶ |
| GV-3 | `box2d`: mass 1.1·10⁻¹⁴, $E + \rho\Phi$ 1.3·10⁻¹⁴, while $E$ alone moves 3.8·10⁻²; `box2d-mhd`: 8.4·10⁻¹⁵, 7.5·10⁻¹⁵ ($E$ 1.8·10⁻²) |
| limiter | `box2d-lim`: 74 304 faces limited, no inadmissible backbone, mass 1.0·10⁻¹⁴, $E + \rho\Phi$ 4.4·10⁻⁷; `box2d-stiff`: 324 424 faces limited, 32 inadmissible backbones, mass 1.2·10⁻¹⁴, $E + \rho\Phi$ 1.6·10⁻² |
| GV-2 | $\max|\mathbf u|/c$ = 2.9·10⁻¹², 1.1·10⁻¹², 2.5·10⁻¹⁴ (scale height 1 on the unit box: a mild stratification) |
| GV-6 | bitwise |
| GV-7 | FNL against CPU: `box2d` ≤ 9.3·10⁻¹⁴, `box2d-mhd` ≤ 9.3·10⁻¹⁴ per field (a vector component scaled by its vector's largest: $\rho u$ is zero by symmetry, round-off only). FNL repeats every CPU number above to the last digit shown, but GV-2 at N = 128, 2.6·10⁻¹⁴ |

![GV (G1): free fall, the energy budget of the closed boxes and the settled box2d](/flume/gravity-force.png)

Drawn with the current scheme (G2a's fourth-order energy source). Left, GV-1: the momentum of the periodic box follows
$\rho g t$, with a pointwise error of one or two units in the last place. Centre, GV-3 over the whole run: the mass
and $W$ stay at round-off, while the plain $E + \rho\Phi$ drifts to 3·10⁻⁵ (`box2d`) through the $O(\Delta x^2)$
wall-row term of G2a; $E$ alone (not drawn) exchanges 3.8 per cent with $\rho\Phi$. Right, `box2d` after 200 steps:
the uniform gas has settled into a stratification, $\rho$ from 1.58 at the floor to 0.42 under the lid.

**The limiter and gravity (owner decision: option B).** The first form put the face-form energy into the backbone. A
cell beside a wall then gets half of the face form (the wall face carries no mass) while its momentum gains the whole
$\rho\mathbf g\,\Delta t$: the stage takes about $\tfrac12\rho|\mathbf u||\mathbf g|\Delta t$ out of its internal energy,
which near vacuum is negative. `box2d-lim` died at step 128 and `box2d-stiff` at step 1 (every backbone inadmissible),
with and without the limiter. Option B gives the backbone the time-centred energy $\rho\mathbf g\cdot(\mathbf u +
\mathbf g\Delta t/2)$, which leaves its internal energy unchanged, and limits the face-form excess with $\Lambda$; both
cases then run (the stiff one's 32 inadmissible backbones come from the fluxes, not from gravity).

### GV (G2a): the order of the energy source

**Why.** The fluxes of a finite-difference scheme are point values of the primitive function $h$ of $\rho\mathbf u$
($\rho u = \frac1{\Delta x}\int h$ over the cell), not of $\rho\mathbf u$. G1's energy source
$\tfrac12 g(F^\rho_{i-1/2} + F^\rho_{i+1/2})$ is then $g\rho u_i + O(\Delta x^2)$, and it took the whole scheme to
second order on any smooth flow under gravity: none of the G1 legs measured an order (free fall is uniform, the
conservation legs are exact at any order). G2a adds the correction
$\tfrac1{24} g(-F^\rho_{i-3/2} + F^\rho_{i-1/2} + F^\rho_{i+1/2} - F^\rho_{i+3/2})$, the cell value of $h$ to fourth
order; each face's weights still sum to one over the cells, so $E + \rho\Phi$ stays conserved exactly with periodic
boundaries and across block interfaces. At a reflecting wall the correction of the wall cell reads the face beyond
it: the mirrored ghosts give $F_{-1} = -F_1$ and $F_0 = 0$, and the wall row's mass balance $F_1 = -\Delta x\,
d\rho_1/dt$ turns the left-over $g\Delta x(F_1 - F_{-1})/24$ into $d/dt(-g\Delta x\rho_1 V/12)$. The exactly conserved
total energy is then $W = \sum(E + \rho\Phi)V + \sum_{\text{lower wall}} g_n\Delta x_n\rho V/12 - \sum_{\text{upper
wall}} g_n\Delta x_n\rho V/12$, an $O(\Delta x^2)$ wall term away from the plain sum (owner decision: W is the
invariant; dropping the correction flux at walls would keep the plain sum but cost the wall cells an $O(\Delta x)$
energy source error). GV-3 asserts W and reports the plain sum. The faces
$i \pm 3/2$ of the block-end cells lie beyond the block: two edge kernels (one per face family, `weno` and
`weno-riemann`) compute them from the ghosts, which takes one more ghost layer (`ngc >= S + 1`, fatal otherwise). With
the positivity limiter the correction is a high-order cell source, scaled by $\Lambda$ with the rest of the face form.

**Case** `wave2d` (`gravity/check.sh`, np 2): every face periodic, $\rho = 1 + 0.2\sin(2\pi(x + y))$ at
$\mathbf u_0 = (1, 0.5)$, $p = 1$, under $\mathbf g = (-0.6, -0.8)$, the `sine-wave` initial condition's `entropy` mode.
The density wave falls freely and exactly: advected by $\mathbf u_0 t + \mathbf g t^2/2$, $\mathbf u = \mathbf u_0 +
\mathbf g t$, $p$ uniform; the momentum source is exact on it, the energy source is not (its $(\rho u)''$ is not zero).
N = 32, 64, 128 to $t = 0.2$, CFL $0.4\,(32/N)^{2/3}$ so the RK3 error stays below the spatial one. Asserted: the L1
orders of $\rho$ and $p$ of the finest pair at least 3.5 (GV-9), on `weno` and on `weno-riemann` (HLLC, sixth-order
correction).

**Results** (CPU and FNL agree to every digit shown):

| Scheme | L1 $p/p_1$ at N = 32, 64, 128 | Orders ($\rho$, $p$), finest pair |
|---|---|---|
| G1 face form alone (`weno`, negative control) | 1.1·10⁻⁵, 2.8·10⁻⁶, 7.0·10⁻⁷ | 2.12, 2.00 |
| `weno` | 3.3·10⁻⁷, 1.2·10⁻⁸, 5.2·10⁻¹⁰ | 5.01, 4.57 |
| `weno-riemann` (HLLC) | 3.4·10⁻⁷, 1.3·10⁻⁸, 5.3·10⁻¹⁰ | 5.00, 4.59 |

At N = 128 the pressure error falls 1300-fold. The negative control is a build with the correction switched off,
not a run-time option; it is not part of `check.sh`.

![GV-9: L1 errors of the falling density wave against N, weno and weno-riemann](/flume/gravity-order.png)

The two schemes lie on top of each other, $\rho$ at slope 5 and $p$ at 4.6 on the finest pair. The pressure error is
the energy source's: on this wave the momentum source is exact and $p$ is uniform, so any pressure error comes from the
energy equation. The control is not drawn (it needs the special build); its errors are in the table.

The other GV legs with the correction (CPU and FNL): GV-1 unchanged (uniform: the correction is zero); GV-3 `box2d`
mass 7.7·10⁻¹⁵, $W$ 7.7·10⁻¹⁵ (the plain $E + \rho\Phi$ 3.3·10⁻⁵, the wall term), `box2d-mhd` 1.2·10⁻¹⁴, 6.1·10⁻¹⁵
(1.6·10⁻⁵); `box2d-lim` 107 984 faces limited (74 304 in G1: the correction is a high-order cell source inside
$\Lambda$), no inadmissible backbone, mass 9.4·10⁻¹⁵, $W$ 1.2·10⁻⁶ (limited cells only); `box2d-stiff` runs, mass
1.2·10⁻¹⁴; GV-2 unchanged; GV-6 and NV-5 bitwise; GV-7 FNL against CPU ≤ 7.9·10⁻¹⁴ per field.

### GV (G2b): the well-balanced mode

**Why.** Without balancing, the truncation errors of the flux difference and of the gravity source do not cancel on a
stratified atmosphere at rest: it drives a stationary flow (GV-2: Mach $6\cdot10^{-4}$ on ten scale heights at N = 32,
$10^{-2}$ across a 2:1 seam) that swamps any perturbation smaller than it. G2b keeps a known equilibrium as an exact
discrete steady state with residual subtraction, the deviation method (Edelmann et al. 2021, §4.3; the alternative of
Xing & Shu 2013 and Li & Xing 2018, a modified Lax–Friedrichs split with the source reconstructed by the same WENO
operator, was prototyped in NumPy and matched it to three digits, but covers the flux-splitting path only).

**The mode.** `[gravity] equilibrium = atmosphere`: the initial atmosphere at rest $q_e$ (zero velocity, region 1's
field, no pulse). Its residual $R(q_e)$ is computed with the run's own space operator (ghosts, boundary conditions,
seam fill, limiter, every source) once per grid, at the first step, after a restart and after every regrid, and
subtracted from every stage residual. At $q = q_e$ the two are the same operator on the same state, at 2:1 seams too.
The end-of-step reflux would still add the equilibrium's coarse–fine flux mismatch, so that mismatch is stored per
register face and the reflux subtracts it. Folding it into $R(q_e)$ instead was tried and fails: the coarse seam cells
then carry a stage residual, the stages leave the equilibrium and the fine ghosts follow (Mach $6\cdot10^{-5}$ on the
AMR box; exact with the reflux off). The residual and the reflux each conserve for the state and for the
equilibrium, so the difference conserves.

**Cases** (`gravity/check.sh`, np 2, N = 32, CFL 0.4): `strat2d`, isothermal, $\mathbf g = (0, -10)$ on the unit box,
ten scale heights (density ratio $2.2\cdot10^4$); `polystrat2d`, polytropic $G = 1.2$, $\mathbf g = (0, -4)$ (ratio
243); `strat2d-mhdy`, `strat2d-mhdx`, EGLM with a uniform field along and across $\mathbf g$ ($|\mathbf B| = 0.05$,
$\beta \approx 7$ at the top); `strat2d-amr`, a refined box $[0.25, 0.75]^2$; x periodic, y `hydrostatic`.

| Leg | Asserted | Unbalanced | Balanced (CPU = FNL) |
|---|---|---|---|
| GV-2b `strat2d`, t = 0.5 | $\max|\mathbf u|/c < 10^{-12}$ | 5.9·10⁻⁴ | 0 |
| GV-2b `polystrat2d`, t = 0.5 | as above | 1.5·10⁻³ | 0 |
| GV-5 `strat2d-mhdy`, `weno` / `weno-riemann` HLLD | as above | 5.9·10⁻⁴ | 0 / 0 |
| GV-5 `strat2d-mhdx`, `weno` / `weno-riemann` HLLD | as above | 2.3·10⁻⁴ | 0 / 0 |
| GV-8 `strat2d-amr`, t = 0.5 | as above | 1.2·10⁻² | 2.8·10⁻¹⁵ (FNL), 3.3·10⁻¹⁵ (CPU) |
| GV-10 pulses $10^{-6}$, $10^{-4}$ on `strat2d`, t = 0.1 | responses per unit amplitude agree to $10^{-2}$ | differ by 97 (reported) | 5.7·10⁻⁵ |
| GV-7 FNL against CPU, the $10^{-4}$ pulse and `strat2d-amr` | per field $\le 10^{-10}$ | — | ≤ 7.8·10⁻¹⁶ |

(The unbalanced MHD values are at t = 0.5 from the trial; `check.sh` runs the MHD legs to t = 0.2.) The balanced rest
is exactly zero where the state equals the stored equilibrium bit for bit; on the AMR box the reflux arithmetic leaves
round-off. GV-10 shows what the mode is for: the balanced response is linear in the amplitude to its second order
(5.7·10⁻⁵ against $A = 10^{-4}$), while the unbalanced one is dominated by its own stationary flow.

![GV (G2b): spurious flows unbalanced and balanced, the refined box, the pulse responses](/flume/gravity-balance.png)

Top left, $\max|\mathbf u|/c$ over the run: the unbalanced atmospheres reach their stationary flow within a few steps
(`strat2d` levels at 6·10⁻⁴, `polystrat2d` at 1.5·10⁻³), the refined box keeps growing to 1.2·10⁻²; balanced, `strat2d`
and `polystrat2d` stay exactly at rest (drawn at 10⁻¹⁸) and the refined box at round-off. Top centre, the unbalanced
refined box at t = 0.5: the flow peaks at (0.52, 0.80), in the second coarse row above the top seam ($y = 0.75$), where
the coarse and fine truncation errors meet. Top right, the same run balanced: round-off, largest in the fine cells
(3.3·10⁻¹⁵), exactly zero in parts of the coarse rows. Bottom, the pulse at t = 0.1 along $x = 0.484$: balanced, the
responses to $10^{-6}$ and $10^{-4}$ coincide (peak $(p/p_e - 1)/A = 0.40$); unbalanced, the curve is the stationary
truncation flow divided by $A$ (peak 333 for $10^{-6}$, 3.4 for $10^{-4}$) and hides the pulse entirely. Bottom right,
the balanced $10^{-6}$ response: the acoustic front, stronger upward, where the density falls. GV-7 scales the momentum by at least
$\max\rho c$: a balanced run's momentum is the perturbation, round-off at rest, and the backends' round-off of the
$O(1)$ residuals that the subtraction cancels (≤ 6·10⁻¹⁵ $\rho c$) would read $O(1)$ against its own maximum.

**Limits.** The subtracted residual is outside the positivity limiter's factor: at the equilibrium the factor is 1 and
the balance exact, but the limiter's guarantee does not cover that term (the equilibrium's truncation error, small
where the atmosphere is resolved). Refused: multi-realm runs (the peer realms' equilibrium would be needed at the
seams) and `reflux_policy = troubled` (its reflux cells change per step). One field and one residual evaluation per
grid.

**Found on the way** (FNL, nvfortran 26.1 with FiNeR): a character `get` right after the failed `get` of an absent
numeric key returns an empty value with no error; `equilibrium = atmosphere` read as "" after the absent `froude`. The
gravity section now reads its string keys first.

### GV (G3): Rayleigh-Taylor

**Why.** Rayleigh–Taylor is the gravity test of the literature, and there it is compared by eye (Liska & Wendroff
2003; Stone et al. 2008, §8.3; Almgren et al. 2010, §8.3): its nonlinear stage has no converged solution. Its linear
stage has one, the growth rate of the seeded mode, and GV-4 asserts it against the linearised compressible Euler
equations of the exact equilibrium profile (owner decision: B, with the nonlinear stage as figures, C). The
incompressible sharp-interface formula $\sigma^2 = Agk$ is not a reference for that: the finite interface width and
compressibility move $\sigma$ by several per cent.

**The equilibrium** (`rayleigh-taylor`, [initial conditions](./initial-conditions#rayleigh-taylor-two-isothermal-layers-heavy-over-light)):
two isothermal layers, cold and dense over hot and light, joined by a tanh in temperature at equal pressure. Isothermal
layers are stably stratified, so the interface is the only instability, and the `hydrostatic` boundary is exact on
them. (Two layers of uniform density would not do: under gravity a uniform-density layer has entropy falling with
height and is convectively unstable, a second instability competing with the first.)

**The oracle** (`gravity/rt_oracle.py`). With perturbations $\propto e^{ikx + \sigma t}$, eliminating $u$ and
$\rho_1$, the adiabatic linear equations reduce to two first-order ODEs in $(v, p_1)$ along the height, with
$v = 0$ at both faces. RK4 shooting from the bottom with bisection on $\sigma$ gives the fundamental mode (no node in
$v$). It is independent of FLUME and checked three ways (`rt_oracle.py selftest`):

| Check | Result |
|---|---|
| closed-form pressure against quadrature of $dp/dh = -\rho g$ | 9.3·10⁻¹³ |
| incompressible ($p_1 = 10^6$), sharp ($\delta$ = 0.004, 0.002) against the two-layer $\sigma^2 = gk\,\Delta\rho/(\rho_h\coth kH + \rho_l\coth kH)$ | 2.3·10⁻², 1.1·10⁻²: $O(k\delta)$, halving with $\delta$ |
| light over heavy | no growing mode |
| 2000 against 4000 RK4 steps | 2.4·10⁻¹¹ |

**Case** `rt2d` (`gravity/check.sh`, np 2): $\mathbf g = (0, -1)$, layers 2:1 ($\rho_1 = 1$, $p_1 = 1$ at the
interface height 0.5), interface width 0.05, one wavelength 0.5 on $[0, 0.5]\times[0.25, 0.75]$ (the mode is down to
$e^{-\pi}$ at the y faces, which the oracle includes as walls), x periodic, y `hydrostatic`, CFL 0.4,
$\sigma = 1.523265$. The seed (amplitude $10^{-5}$) is divergence-free: a seed in $v$ alone starts sound waves that
beat against the mode, and with one the fitted rate was 0.28 % to 0.92 % off depending on the window
(g = 0.1, N = 64).

| Leg | Asserted | Measured (CPU) |
|---|---|---|
| GV-4a unseeded, one negligible step | the layers and their hydrostatic ghosts against the closed form, 10⁻¹³ | interior 1.1·10⁻¹⁶, ghosts 2.2·10⁻¹⁶ |
| GV-4b unseeded to t = 1 | balanced $\max|\mathbf u|/c < 10^{-12}$; unbalanced reported | balanced 0; unbalanced 2.2·10⁻⁶ |
| GV-4c seeded, N = 32 and 64 to t = 4 | the slope of $\ln a$ over $t \in [2, 4]$ against $\sigma$, 10⁻³ | 1.522543 (4.7·10⁻⁴), 1.522724 (3.6·10⁻⁴); the two agree to 1.2·10⁻⁴ |
| GV-7 FNL against CPU | `rt-8`, 10⁻¹⁰ | ≤ 1.3·10⁻¹² ($\rho$): round-off amplified by the instability, $e^{\sigma t} \approx 400$ over the run; FNL repeats every GV-4 number to the digits shown |

$a$ is the amplitude of the seeded mode in the vertical velocity, $(\sum_y |\tfrac2{N_x}\sum_x v\,e^{-ikx}|^2
\Delta y)^{1/2}$. The unbalanced scheme's stationary flow (GV-4b) is as large as the seed but uniform in x: its
content at the seeded wavenumber is 5.8·10⁻²³, so GV-4c runs unbalanced. The two resolutions agree more closely
than either agrees with the theory, and the error falls as the window starts later (1.1·10⁻³, 4.7·10⁻⁴, 2.6·10⁻⁴
for $t_0$ = 1.5, 2, 2.5 at N = 32): the residual is the seed's transient decaying in the window, not the
discretisation, which is below 1.2·10⁻⁴ with three cells across the interface width.

![GV-4: the seeded mode, its local growth rate and its vertical profile against the linear theory](/flume/gravity-rt-linear.png)

Left, the mode grows from $3.6\cdot10^{-6}$ to $7.9\cdot10^{-4}$, at the theory's slope once the transient has passed.
Centre, the local rate between checkpoints: it settles by $t \approx 1.5$ ($\sigma t \approx 2.3$) and then
oscillates around the theory by a few $10^{-3}$, the beat of the seed's stable components. Right, the vertical
profile $|v_k(y)|$ at t = 4 against the oracle's eigenfunction: 1.6·10⁻³ (N = 32) and 1.3·10⁻³ (N = 64) apart at
most. FLUME was never given that shape: the seed is a Gaussian.

**The nonlinear stage (figures, no oracle).** `rt2d-nl`: the same layers with a sharp interface (width 0.01) and a
$10^{-2}$ seed, two wavelengths across the unit square, N = 128 to t = 5 ($\sigma = 1.869$ by the oracle); `rt2d-mhd`:
the same with EGLM and a horizontal field $B_x = 0.1$, half the critical field of the mode
($\sqrt{\Delta\rho\,g/(2k)} = 0.20$ for a sharp incompressible interface; $\beta \approx 200$). Runs by
`make_doc_runs.sh <dir> rt`.

![Rayleigh-Taylor single mode, hydro and MHD, density at t = 2, 3, 4, 5](/flume/gravity-rt.png)

Hydro (top): the spikes roll up through secondary Kelvin–Helmholtz instabilities by t = 3 and reach the floor by
t = 5. The local growth rate of the seeded mode peaks at 1.77 at t ≈ 1.5, 5 % under the linear 1.869, before the
mode saturates. This is not a measurement: the 0.01 interface spans 1.3 cells and the $10^{-2}$ seed is nonlinear
almost from the start, which is why GV-4 runs its own resolved case. MHD (bottom; each row has its own colour scale):
the field slows the mode, with a peak rate of 1.49, 0.80 of the hydro rate against the sharp-interface estimate
$\sqrt{1 - (B/B_c)^2} = 0.87$, and suppresses the roll-ups, as in Stone et al. (2008), §8.4. Both runs are periodic
with two identical wavelengths, each mirror-symmetric about its own axis. Hydro keeps both symmetries to round-off,
which grows with the flow to 6·10⁻¹⁰ at t = 5. MHD does not. Its mirror asymmetry stays at round-off to t ≈ 2,
then grows at a rate of 6 to 10, against the mode's 1.5: 2.3·10⁻¹⁰ at t = 3, 1.5·10⁻⁷ at t = 4, 2.8·10⁻³ at t = 5
(invisible in the figure). `weno-riemann` with HLLD delays it by about one time unit, 9.8·10⁻¹² at t = 3 and
2.1·10⁻⁵ at t = 5, but it grows as fast. The peak growth rates of the mode agree (1.4931 and 1.4914). It is therefore
not the low-β flux-splitting instability of [#77](https://github.com/szaghi/adam/issues/77): it appears under both
flux families and at $\beta \approx 200$, and only with the field on. The next paragraph traces it to the divergence
cleaning.

**The mirror symmetry of the MHD run ([#83](https://github.com/szaghi/adam/issues/83)).** Why: a symmetry that the
equations keep and the scheme loses is either physics (an unstable mode of the flow, fed by round-off) or numerics,
and a user who reads the asymmetry needs to know which. The measure is
$a(t) = \max |\rho(x, y) - \rho(\lambda - x, y)|$, which the exact solution keeps at zero
(`atmosphere_oracle.py symmetry`, with its mean growth rate over $t \in [2, 4]$). The runs are `rt2d-mhd` at N = 128,
with one cleaning parameter changed at a time (`make_doc_runs.sh <dir> rt-sym`, the FNL build, about 2 min each on two
A30 GPUs). The figures are drawn by `make_doc_figures.py --only gravity-rt-symmetry`.

| Run (`rt2d-mhd`, N = 128) | rate over t ∈ [2, 4] | $a$ at t = 4 | $a$ at t = 5 | max \|div B\| at t = 4 |
|---|---|---|---|---|
| hydro (`rt2d-nl`, CPU) | 1.78 | 1.2·10⁻¹⁰ | 6.3·10⁻¹⁰ | — |
| `divergence_control = none` | 2.04 | 1.1·10⁻¹⁰ | 7.0·10⁻⁶ | **8.6** |
| `glm`, `glm_ch` = 3 | 3.27 | 3.5·10⁻¹⁰ | 4.2·10⁻⁶ | 0.87 |
| `eglm`, `glm_ch` = 1.5 | 4.38 | 9.9·10⁻¹⁰ | 7.0·10⁻⁶ | 1.00 |
| `eglm`, `glm_ch` = 3 (the case's default) | 6.25 | 1.4·10⁻⁷ | 2.5·10⁻³ | 0.84 |
| `eglm`, `glm_ch` = 6 | 10.91 | 2.7·10⁻² | 5.2·10⁻² | 0.81 |
| `eglm`, `glm_ch` = 3, `glm_alpha` = 0 | 6.87 | 5.1·10⁻⁸ | 2.8·10⁻³ | 0.86 |
| `eglm`, `glm_ch` = 3, `glm_alpha` = 1 | 2.91 | 1.4·10⁻⁹ | 2.8·10⁻³ | 0.87 |
| `eglm`, `glm_ch` = 3, `glm_alpha` = 5 | 4.85 | 4.0·10⁻⁹ | 2.8·10⁻³ | 0.91 |
| `eglm`, `glm_ch` = 3, `glm_damping_length = min-cell` | 4.09 | 5.4·10⁻¹⁰ | 3.0·10⁻³ | 0.76 |

![Magnetised Rayleigh-Taylor: mirror asymmetry per cleaning and damping, rate against c_h](/flume/gravity-rt-symmetry.png)

- **The cleaning drives it.** Without cleaning, the MHD run keeps its symmetry about as well as hydro does (rate 2.04
  against 1.78). GLM raises the rate to 3.27, and EGLM to 4.4, 6.3 and 10.9 at `glm_ch` = 1.5, 3 and 6. The excess
  over the uncleaned run grows by about 1.5 `glm_ch`. It is first seen on the shear flanks between bubble and spike.
- **It is not a time-stepping artefact.** At CFL 0.2, which takes 2.5 times the steps, the EGLM rate at `glm_ch` = 3
  is 5.31, against 6.25 at CFL 0.4. A fixed gain per step would have doubled the excess.
- **It is not in the cleaned equations.** `gravity/glm_linear.py` linearises the GLM and EGLM systems about an exact,
  divergence-free double shear layer ($V_0 = 0.5$, width 0.05, field 0.3 along the flow) and solves the spectrum by
  Fourier collocation. At $k = 2\pi$ the largest growth rate is 1.008 (N = 128) for both systems and for every
  $c_h$ from 0.75 to 6 (right panel). Damping ($\alpha c_h / L = c_h$) and a sheared field leave it independent of
  $c_h$. At $k = 4\pi$ it is 0.52 at N = 128 and 0.40 at N = 256, independent of $c_h$ at both: the critical layer is
  under-resolved there, so the value itself is not converged. The $c_h$-proportional growth is therefore a property of
  the discretised cleaning. It acts where the discrete div B is not small: max $|\nabla\cdot\mathbf{B}| \approx 0.8$,
  so $\Delta x\,|\nabla\cdot\mathbf{B}| / |\mathbf{B}| \approx 6\,\%$ per cell. The GLM system without
  cleaning ($c_h = 0$, conservative MHD without the Powell terms) is the negative control: it is not hyperbolic, and
  its rate rises with $k$ (1.6 at $4\pi$, 3.2 at $8\pi$).
- **The damping does not cure it.** `glm_alpha` from 0 to 5, or $L = \Delta x$, moves $a$ at t = 4 within a factor of
  about 40, with no monotone trend. Every `glm_ch` = 3 run ends with the same jump between t = 4.5 and 5, to
  $a \approx 3 \cdot 10^{-3}$. The `glm_ch` = 1.5, GLM and uncleaned runs show a smaller late jump instead, to about
  $10^{-5}$: a current sheet on the spike's stem that narrows to 2 cells at N = 64. It is a grid-scale effect, common to
  every run.
- **The cleaning still does its job on div B:** 0.8–1.0 at t = 4 under every cleaning, against 8.6 without.

**Guidance.** For a magnetised run whose symmetry matters (mirror-symmetric instabilities, convergence studies of the
small scales), set `glm_ch` only just above the fastest signal speed. `glm_ch_check = error` stops the run if the
flow outruns it. Here 1.5 is enough: the check never fired in the `glm_ch` = 1.5 run. Prefer `glm` to `eglm`: at the
same speed EGLM's excess is about three times GLM's (4.2 against 1.2).
`none` keeps the symmetry best but lets div B grow tenfold, so use it only as a diagnostic. With EGLM at
`glm_ch` = 3, take the asymmetry above about $10^{-6}$ after $t \approx 4.3$ as numerical, not physical. Run at a
second `glm_ch` to check: a physical asymmetry does not depend on it.

## Self-gravity (M8)

Self-gravity ([#81](https://github.com/szaghi/adam/issues/81)) solves $\nabla^2\Phi = 4\pi G(\rho - \bar\rho)$ every
Runge–Kutta stage and applies gravity through the cell potential. The design (the P0 comment on #81) uses the spatial
discretisation of Mullen et al. (2021), with the 7-point Laplacian:
- **face gravity:** $g_{d,i+1/2} = -(\Phi_{i+1} - \Phi_i)/\Delta x_d$;
- **cell gravity:** the mean of the cell's two faces;
- **momentum:** $\rho$ times the cell gravity;
- **energy:** $\tfrac12\sum_d (g_{d,i+1/2} F^\rho_{d,i+1/2} + g_{d,i-1/2} F^\rho_{d,i-1/2})$, from the mass fluxes the
  continuity equation uses.

Summed over the cells, the energy source is $\sum_i \Phi_i (F_{i+1/2} - F_{i-1/2})/\Delta x$. So $E + \rho\Phi$ is
conserved exactly with a static potential, and $E + \tfrac12\rho\Phi$ with a self-consistent one: in the semi-discrete
form, given a Poisson solution at round-off. The scheme is second order.

### SV-1 (P1): gravity from a cell potential

**Why.** The potential path (`add_potential_sources`, the face and cell gravity in the positivity limiter, the
hydrostatic ghost from the potential drop $-(\Phi_{\text{ghost}} - \Phi_1)$) must be right before a Poisson solver feeds
it. `[gravity] potential = linear` writes the constant g as $\Phi = -\mathbf g\cdot\mathbf x$ and applies it through that
path alone (the constant-g kernels get zero), so the M7 oracles apply unchanged.

**The leg** (`selfgravity/check.sh`, CPU and FNL, 2 ranks): the M7 cases with `--potential linear`.

| Case | Asserted | CPU | FNL |
|---|---|---|---|
| iso2d, iso2d-mhd | hydrostatic ghosts against the profile, 1e-13 | 2.2·10⁻¹⁶ | 2.2·10⁻¹⁶ |
| freefall2d | ρu = ρgt, E = E₀ + ρg²t²/2, 1e-12 | 5.6·10⁻¹⁶ | 3.3·10⁻¹⁶ |
| box2d, box2d-mhd | mass and E + ρΦ in the closed box, 1e-12 | 1.3·10⁻¹⁴, 7.5·10⁻¹⁵ | 1.7·10⁻¹⁴, 8.8·10⁻¹⁵ |
| box2d-lim | the run completes (limiter acting); drift reported | 4.4·10⁻⁷ | 4.4·10⁻⁷ |
| wave2d, 3 resolutions | L1 order of ρ and p ≥ 1.8 on the finest pair | 2.12, 2.00 | 2.12, 2.00 |

In the closed box the potential path's energy budget closes on $E + \rho\Phi$ alone. The constant path closes on W,
which adds the $O(\Delta x^2)$ wall term of its fourth-order source: without that term its $E + \rho\Phi$ drifts by
$2\cdot10^{-6}$. With the limiter acting, the drift (4.4·10⁻⁷) is that of the constant path (1.2·10⁻⁶) in size. The
constant-g path is unchanged bitwise: the regression suite (35/35 on CPU and FNL) holds it.

![SV-1: energy budget of the closed box through both paths, and the falling-wave errors against N](/flume/selfgravity-potential.png)

Left, the closed box at N = 64: the drift of each path's conserved energy is at round-off, while E alone exchanges with
$\rho\Phi$ by up to $6\cdot10^{-3}$. Right, the falling wave (GV-9's exact solution) to t = 0.2. The potential path is
second order (2.12 in ρ, 2.00 in p), as its face form without a fourth-order correction must be; the constant path
gives 5.0 and 4.6. Run with `make_doc_runs.sh <dir> selfgravity` and `make_doc_figures.py --only selfgravity-potential`.

## Ranks without blocks (DV-1)

**Why.** The tree never splits a family of sibling blocks across ranks. `tree%mpi_redistribute` cuts the Morton-ordered
blocks into equal chunks, and when a cut falls inside a family it moves the whole family to one side. A grid with
fewer families than ranks therefore leaves ranks with no block: four blocks at level 1 are one family, all on rank 0
at np 2. Every code path that sizes an array by the local `blocks_number` meets a zero there. On FNL that is a
zero-byte device allocation: `acc_malloc(0)` returns a null pointer (measured, nvhpc 25.5), and FUNDAL reports a null
pointer as a failed allocation. Every MHD run therefore died at initialisation
([#82](https://github.com/szaghi/adam/issues/82)): the div(B) history copied its seam-face flags,
`divb_seam(blocks_number, 6)`, to the device. The array is now sized by the block capacity `nb`, as every other
block-indexed device array is. Small multi-rank tests rarely leave a rank empty, so this leg makes the empty rank the
test.

**The leg** (`ranks/check.sh`, `ranks_oracle.py`):
- **Grid:** four blocks of 16² on $[0,1]^2$ (one family), 20 steps at CFL 0.4.
- **Cases:** nine, which together exercise every model and option that allocates block-indexed state:
  - `repro`, the issue's input: a uniform EGLM state, y walls, no gravity;
  - the `strat2d` atmospheres with a 10 % pressure pulse: Euler; MHD under `none`, `glm` and `eglm`; EGLM with the
    positivity limiter, with `equilibrium = atmosphere`, and on `weno-riemann` with HLLD;
  - EGLM with runtime regrids every 5 steps to level 2. A Löhner marker on the density acts at the regrids only, so
    rank 1 starts empty and the first regrid repopulates it.
- **Runs:** each case on 1 and on N = 2 ranks.
- **Asserted:**
  - **ownership:** in the np 2 run every rank but 0 owns no block at step 0, so the leg cannot become vacuous if the
    distribution changes; in the AMR case, every rank owns blocks at the end;
  - **identity:** np 2 = np 1 **bitwise**, every field and every cell at the last checkpoint. With every block on rank
    0 the arithmetic is the same, and the reductions add exact zeros.

| Case (FNL, 2 A30 GPUs, nvhpc 25.5) | blocks per rank, step 0 → 20 | before the fix, np 2 | after: np 2 vs np 1 |
|---|---|---|---|
| `repro` | 4 / 0 → 4 / 0 | **aborts** at initialisation | bitwise |
| Euler | 4 / 0 → 4 / 0 | runs | bitwise |
| MHD, `none` | 4 / 0 → 4 / 0 | **aborts** | bitwise |
| MHD, `glm` | 4 / 0 → 4 / 0 | **aborts** | bitwise |
| MHD, `eglm` | 4 / 0 → 4 / 0 | **aborts** | bitwise |
| EGLM, positivity limiter | 4 / 0 → 4 / 0 | **aborts** | bitwise |
| EGLM, `equilibrium = atmosphere` | 4 / 0 → 4 / 0 | **aborts** | bitwise |
| EGLM, `weno-riemann` + HLLD | 4 / 0 → 4 / 0 | **aborts** | bitwise |
| EGLM, runtime AMR | 4 / 0 → 8 / 8 | **aborts** | bitwise (4096 cells) |

The leg takes 3.6 min on FNL, and the pre-fix executable fails it (negative control: 8 of 9 cases abort; Euler, which
never computes the div(B) history, runs). The CPU backend passes the leg as well, all nine cases bitwise in 2.8 min.
The issue reported it unaffected, since a zero-size host array is legal. `divergence_control = none` failed too,
because the div(B) history is skipped only for Euler.

![DV-1: blocks per rank before and after the first regrid, and the np 2 - np 1 difference](/flume/ranks.png)

Left, the speed at step 4, before the first regrid, with the block outlines coloured by owner: the four blocks are all
rank 0's. Centre, step 20: the first regrid (step 5) refined the family into four, and the distribution gave each rank
two families, 8 blocks. Right, the blocks of each rank against the step. Over all 21 checkpoints the largest difference
between the np 2 and np 1 states, over every field and every cell, is 0.

Run with `make_doc_runs.sh <dir> ranks` (the AMR case, a checkpoint every step) and
`make_doc_figures.py --only ranks`.

## Unit tests

| Test | What it pins |
|---|---|
| `test_flume_euler_library` (+ `_fnl`) | Euler eigensystem, flux, Roe average, split consistency; RS(q, q) = f(q) for LLF/HLL/HLLC, HLLC exact on a contact, positive first-order updates; the `si` descaler of planar and static states within $10^8$ of the largest (host); device = host |
| `test_flume_mhd_library` (+ `_fnl`) | MHD, GLM and EGLM eigensystems (including degenerate states), fluxes, auxiliary variables, cyclic invariance; the `si` descaler of planar, field-free and static states within $10^8$ of the largest; $\alpha_f$, $\alpha_s$ to full relative accuracy at high and low β with $\mathbf{B}$ near the normal, against quadruple precision ([#79](https://github.com/szaghi/adam/issues/79)) (host); device = host |
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
