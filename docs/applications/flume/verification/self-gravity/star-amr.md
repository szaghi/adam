# The polytrope on AMR: static boxes, runtime regrids, Jeans refinement

| ID | Group | Milestone | Type | Backends | Figure | Regression |
|---|---|---|---|---|---|---|
| SV-14 to SV-17, rg7 | Self-gravity | [#87](https://github.com/szaghi/adam/issues/87) (follow-up) · P3, P4; [#88](https://github.com/szaghi/adam/issues/88) (follow-up); M9 ([#91](https://github.com/szaghi/adam/issues/91)) | P+C | CPU, FNL | yes | `selfgravity-star-amr` |

## Case and scope

The n = 1 polytrope of [SV-7 and SV-8](./polytrope) on AMR grids, under isolated self-gravity by the multigrid's
symmetric James step ([SV-13](./multigrid-amr)). The star is smaller than on the uniform grid (radius 0.3 instead of
0.5), so that the 2:1 seams lie in the atmosphere and not on the stellar surface.

What the case catches:
- **SV-14, the star on a static AMR box** ([#87](https://github.com/szaghi/adam/issues/87) P3): the equilibrium with the
  star refined must have less spurious motion than on the uniform base; the collapse must keep the exact exchange of
  $E$ and $W$ (the per-stage reflux carries the seam gravity work); a regrid that refines a boundary block must stop the
  run;
- **SV-15, the collapse under runtime AMR** ([#88](https://github.com/szaghi/adam/issues/88)): a prolongation that
  builds cold, fast cells mid-collapse drives the run non-finite; `conservative-state` must run it to the end;
- **SV-16, Jeans refinement** ([#91](https://github.com/szaghi/adam/issues/91)): a Jeans marker must refine the star
  before the collapse, so that no regrid happens mid-collapse;
- **SV-17, dual energy** ([#91](https://github.com/szaghi/adam/issues/91)): `euler-de` must keep the cold infall from
  draining with the default prolongation;
- **rg7, the Jeans-marker contract:** a Jeans marker without self-gravity, or with an inverted hysteresis band, must be
  refused at initialization.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations, scheme, boundaries | as [SV-7, SV-8](./polytrope): Euler, $\gamma = 2$, isolated self-gravity with $G = 1$, WENO5 flux splitting, SSP-RK3, `wall-inviscid` on every face | [models](../../models#self-gravity) |
| Poisson solver | the multigrid's James step for isolated boundaries, its coarsest level on the FFT solver | [`adam_poisson_mg_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_poisson_mg_object.F90), [`adam_poisson_fft_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_poisson_fft_object.F90) |
| Reflux | `reflux_cadence = stage`: the seam fluxes, and the seam gravity work, reconciled every Runge–Kutta stage | [numerics](../../numerics#adaptive-mesh-refinement) |
| Boundary refinement | a regrid that refines a block on the domain boundary is fatal (the James step needs every leaf on the boundary at the coarsest level) | [`adam_poisson_mg_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_poisson_mg_object.F90) |
| `[amr] regrid_prolongation` | `conservative` (default) or `conservative-state` (below) | `refine_conservative`, `refine_conservative_state`, [`adam_field_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_field_object.F90) |
| Jeans marker | `[amr_marker_N] mode = 5`: refine where the cell spacing exceeds $\lambda_J$ / `jeans_cells`; `jeans_cells_derefine` > `jeans_cells` sets the hysteresis band; needs self-gravity | [numerics](../../numerics#adaptive-mesh-refinement); [input](../../input#setting-up-runtime-amr) |
| Dual energy | `physical_model = euler-de`, switch 0.05 | [DE](../euler/dual-energy); [models](../../models#dual-energy) |

### The conservative-state prolongation

The opt-in `conservative-state` ([#88](https://github.com/szaghi/adam/issues/88)) sets the $2^d$ children of a parent
cell $c$ to

$$
q_k = q_c + \theta\,\big(\delta_k - \langle\delta\rangle\big), \qquad
\delta_k = \sum_{\mathbf a \ne \mathbf 0} w_{\mathbf a}\,\big(q_{c + \mathbf a \odot \boldsymbol\sigma_k} - q_c\big),
\qquad w_{\mathbf a} = \prod_{\text{refined axes}} \big(\tfrac34, \tfrac14\big)_{a_i},
$$

with $\boldsymbol\sigma_k$ the side of child $k$ per axis ($\pm 1$), $\mathbf a \in \{0, 1\}^d$ and
$\langle\delta\rangle$ the mean of the $\delta_k$ over the children. $q_c + \delta_k$ is the `linear` prolongation, a
mixture of whole coarse states with the same weights for every variable; removing the mean makes the children average
to the parent. $\theta \in [0, 1]$, one factor for every variable of the cell, is the largest that keeps every child in
the range of the parent and its $3^d$ neighbours. Linear data have $\langle\delta\rangle = 0$ and keep $\theta = 1$.

The default `conservative` prolongation limits each variable by its own slope
([regrid round trip](../runtime-amr/regrid-round-trip#the-conservative-prolongation)), so its children are not
mixtures of whole states; in a collapse that is what builds the cold, fast cells of SV-15. A bound on the children's
specific internal energy or entropy, minmod slopes, and one $\theta$ per variable each failed a case of the sweep over
the star's radius and regrid cadence that the single $\theta$ runs to the end.

The single $\theta$ is also its price: an extremum of any one variable in the neighbourhood (a vortex core) flattens
the children of every variable. The tracked vortex of [AV](../runtime-amr/accuracy) (av2) then gives an L1 error 1.68
times the uniform fine run's, above the bound 1.5 (1.23 with `conservative`). Hence an option for collapsing flows, not
the default. The round-trip identities of `conservative-state` (exact integral, range) are asserted by rg0 in the
[regrid round trip](../runtime-amr/regrid-round-trip).

## Initial and boundary conditions

| Item | Value |
|---|---|
| Star | the [SV-7](./polytrope) polytrope with radius $a = 0.3$: the 2:1 seam of the box $[-1/2, 1/2]^3$ lies in the atmosphere, 0.2 off the surface |
| Box, base | $[-1, 1]^3$, octree, base $32^3$ (blocks of $8^3$, 4 ghosts), `wall-inviscid` on every face; 2 ranks |
| SV-14, grid | the box $[-1/2, 1/2]^3$ refined once over the $32^3$ base (`--box`: 120 blocks, the star at the $64^3$ spacing, 9.6 cells per radius) |
| SV-14, equilibrium | as SV-7 (atmosphere $10^{-4}\rho_c$, limiter `cell`, to $t = 2$), on the AMR grid and on the uniform $32^3$ grid |
| SV-14, collapse | as SV-8 (atmosphere $10^{-2}\rho_c$, no limiter, to $t = 1.5$, CFL 0.4, 0.2, 0.1) on the AMR grid with `reflux_cadence = stage` |
| SV-14, refusal | `starc --regrid 2 --boundary-box`: a runtime marker on $[1/2, 1]^3$; the first regrid refines a boundary block |
| SV-15, runtime-AMR collapse | `starc --radius 0.3 --reflux-cadence stage --regrid 5`: the uniform $32^3$ base (4.8 cells per radius), atmosphere $10^{-2}\rho_c$, no limiter, to $t = 1.5$; a regrid every 5 steps by a Löhner marker on the density (one level more, the atmosphere's density as its floor), the first one refining the 8 inner blocks; `regrid_prolongation = conservative-state`, and the default `conservative` reported |
| SV-16, Jeans refinement | `starc --radius 0.22 --reflux-cadence stage --regrid 5 --jeans 8`: 3.5 cells per radius on the base; its Jeans length is $\sqrt 2\,a$, 5 base cells, since $\lambda_J^2 = \pi\gamma p/(G\rho^2) = \pi K/G = 2a^2$ for the $\gamma = 2$ polytrope, uniform inside it and unchanged by an adiabatic collapse; the SV-15 Löhner marker and a Jeans marker of 8 cells per Jeans length, the initial AMR on (`amr_iterations = 1`), the default prolongation; the same star without the Jeans marker reported |
| SV-17, dual energy | the SV-15 case with the default prolongation, no Jeans marker, and `--physical-model euler-de` |
| rg7 | the regression `sod-x` input ([Sod](../euler/sod)) to $t = 0.002$, 2 ranks, with one Jeans marker (`mode = 5`, `jeans_cells = 8`): as is (no `[gravity]`), and with `jeans_cells_derefine = 4` |
| Regression `selfgravity-star-amr` | `starc --radius 0.3 --box --reflux-cadence stage`: SV-14's collapse, 120 blocks of $8^3$, 20 steps |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/selfgravity/make_polytrope.py` | writes `star` and `starc`; on AMR `--radius`, `--box`, `--regrid`, `--prolongation`, `--boundary-box`, `--reflux-cadence`, `--jeans`, `--physical-model` |
| `src/tests/flume/verification/selfgravity/star.sh` | runs the legs `sv14`, `sv15`, `sv16`, `sv17` (`--np N`, `--legs LIST`, `--keep`; `FLUME_EXE` selects the backend) |
| `src/tests/flume/verification/selfgravity/star_oracle.py`, `jeans_oracle.py` | the drift reductions (equilibrium) and $E + W$ against $\Delta t$ (collapse) |
| `src/tests/flume/verification/regrid/check.sh --legs rg7` | rg7, the Jeans-marker refusals |
| `src/tests/flume/regression/selfgravity-star-amr/input.ini` | the frozen SV-14 collapse on AMR, 20 steps ([#87](https://github.com/szaghi/adam/issues/87) P4) |
| `src/tests/flume/verification/make_doc_runs.sh`, group `isolated` | `sv14-{uniform,box}`, `sv14-collapse-{0.4,0.2,0.1}` |
| `src/tests/flume/verification/make_doc_runs.sh`, group `collapse-regrid` | `sv15-{conservative-state,conservative}`, `sv16-{jeans,control}` and `sv17-de`, a checkpoint every 2 steps (the `conservative` and `control` runs stop by design) |
| `src/tests/flume/verification/make_doc_figures.py`, `fig_selfgravity_star_amr`, `fig_selfgravity_star_regrid` | draws the figures (`--only selfgravity-star-amr selfgravity-star-regrid`) |

```bash
cd src/tests/flume/verification/selfgravity && ./star.sh --legs sv14,sv15,sv16,sv17       # CPU (exe/adam_flume_cpu_fft)
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl_fft ./star.sh --legs sv14,sv15,sv16,sv17   # FNL
cd ../regrid && ./check.sh --legs rg7                                                    # the Jeans-marker contract
```

## Results

**SV-14, the star on a static AMR box:**

![SV-14: the smaller star on AMR, its core and spurious motion against the uniform grid; the collapse's energy drift; the star on its composite grid](/flume/selfgravity-star-amr.png)

*Top: the star of radius 0.3 on the uniform $32^3$ grid (4.8 cells per radius) and with the box refined (9.6): on
the AMR grid the core swings half as far and the interior Mach number stays at half. Bottom left: the collapse on the
AMR grid, $E + W$ at the Runge–Kutta order while $E$ and $W$ exchange 81 % of $\lvert E + W\rvert$ (their curves
coincide). Bottom right: the equilibrium star at $t = 2$ on its composite grid, the radius dashed; the 2:1 seam lies
in the atmosphere.*

| Leg | Asserted | Bound | CPU | FNL |
|---|---|---|---|---|
| SV-14 | the largest central-density drift to $t = 2$: $1.65\cdot10^{-1}$ (uniform $32^3$), $8.24\cdot10^{-2}$ (AMR); reduction | $\ge 1.5$ | 2.00 | 2.00 |
| SV-14 | the largest Mach number inside $0.8a$: 0.414, 0.189; reduction | $\ge 1.5$ | 2.19 | 2.19 |
| SV-14 | the largest Mach number in the whole star: 15.2 (uniform: the surface spans 5 cells), 0.36; reduction | reported | 42 | 42 |
| SV-14 | collapse on AMR, drift of $E + W$ at CFL 0.4, 0.2, 0.1: $4.76\cdot10^{-7}$, $6.07\cdot10^{-8}$, $7.72\cdot10^{-9}$; order in $\Delta t$ | $\ge 2.5$ | 2.97, 2.97 | 2.97, 2.97 (the same drifts) |
| SV-14 | the largest $\lvert\int\rho\mathbf u\rvert$; mass drift | $10^{-14}$; $10^{-14}$ | $2.2\cdot10^{-19}$; $2.4\cdot10^{-16}$ | $4.3\cdot10^{-19}$; $2.4\cdot10^{-16}$ |
| SV-14 | a regrid that refines a boundary block | refused, exit 1 | refused | refused |

**SV-15 to SV-17, the runtime-AMR collapse:**

![SV-15 and SV-16: the coldest cell of the runtime-AMR collapse with the two conservative prolongations and with the Jeans marker; the drift of the mass and of E + W](/flume/selfgravity-star-regrid.png)

*Top row, SV-15: the SV-14 star collapsing on the uniform $32^3$ base, a regrid every 5 steps. Left: the regrid at step 5 refines the
8 inner blocks and lowers the smallest $e/\rho$ from $7.0\cdot10^{-3}$ to $5.0\cdot10^{-3}$ with either prolongation.
With the default `conservative` the coldest cells then keep cooling while an infall reaches Mach 10: $3.5\cdot10^{-4}$
by step 62, non-finite at step 64 ($t = 0.70$). With `conservative-state` the smallest $e/\rho$ never falls below
$5.07\cdot10^{-3}$, the Mach number stays below 1.8 and the run reaches $t = 1.5$. Right: the mass at round-off in both;
$E + W$ jumps at the regrid (the potential is solved on the new grid), by $1.5\cdot10^{-3}$ with `conservative` and
$3.5\cdot10^{-5}$ with `conservative-state`, and drifts at the time integration's error between regrids. Bottom row,
SV-16: the smaller star refined only by the Löhner marker at step 5 cools to non-finite at step 27; with the Jeans
marker the initial AMR refines it before the collapse (120 blocks, 64 on the base), no regrid changes the grid, the
coldest cell recovers and the run reaches $t = 1.5$ with the default prolongation, $E + W$ drifting by
$7.6\cdot10^{-7}$ with no regrid jump. Bottom row, SV-17: the SV-15 star with `euler` (the red curve of the top row)
and with `euler-de`: the floor lifts the coldest cells to their adiabat in 6 steps, adding $4.1\cdot10^{-7}$ to the
total energy, and the run reaches $t = 1.5$ with the default prolongation.*

| Leg | Asserted | Bound | CPU | FNL |
|---|---|---|---|---|
| SV-15 | `conservative-state`: reaches $t = 1.5$, with a regrid that changes the grid | completes | $t = 1.5$, 1 regrid changes the grid | $t = 1.5$, 1 regrid changes the grid |
| SV-15 | `conservative-state`: mass and momentum drifts | $10^{-13}$ | mass $1.2\cdot10^{-16}$, momentum $3.8\cdot10^{-18}$ | mass $1.2\cdot10^{-16}$, momentum $7.6\cdot10^{-18}$ |
| SV-15 | `conservative`, the default | reported | non-finite at step 64 | non-finite at step 65 |
| SV-16 | Jeans marker: the initial grid refined (blocks at step 0 above the 64 of the base); reaches $t = 1.5$ | completes | 120 blocks at step 0; $t = 1.5$ | 120 blocks at step 0; $t = 1.5$ |
| SV-16 | Jeans marker: mass and momentum drifts | $10^{-13}$ | mass $2.1\cdot10^{-16}$, momentum $1.7\cdot10^{-18}$ | mass $2.1\cdot10^{-16}$, momentum $3.3\cdot10^{-18}$ |
| SV-16 | without the Jeans marker | reported | non-finite at step 27 | non-finite at step 33 |
| SV-17 | `euler-de`: reaches $t = 1.5$ | completes | $t = 1.5$ | $t = 1.5$ |
| SV-17 | `euler-de`: mass and momentum drifts | $10^{-13}$ | mass $1.2\cdot10^{-16}$, momentum $4.8\cdot10^{-19}$ | mass $1.2\cdot10^{-16}$, momentum $9.6\cdot10^{-19}$ |
| SV-17 | the energy the entropy floor added | reported | $4.1\cdot10^{-7}$ in 6 steps | $4.1\cdot10^{-7}$ in 6 steps |

**rg7, the Jeans-marker contract:**

| Input | Expected | CPU | FNL |
|---|---|---|---|
| Jeans marker, no self-gravity | fatal: "needs [gravity].(self_gravity) = .true." | pass | pass |
| `jeans_cells_derefine = 4` | fatal: "needs 0 < jeans_cells < jeans_cells_derefine" | pass | pass |

Related cases: the uniform-grid star is [SV-7, SV-8](./polytrope); the solver on AMR is [SV-13](./multigrid-amr);
self-gravity through regrids is [SV-11](./runtime-amr); the dual-energy model is [DE](../euler/dual-energy); the
library regrid is the [regrid round trip](../runtime-amr/regrid-round-trip).

## Findings

- **Isolated self-gravity on AMR holds the star and keeps the exchange exact** (SV-14). Refining the star halves its
  spurious motion against the uniform base, as a uniform grid of the fine spacing would; the collapse drifts at third
  order in $\Delta t$ with the momentum at round-off, as on the uniform grid (SV-8): the James step is symmetric and
  the per-stage reflux carries the seam gravity work.
- **Keep refinement boundaries off the stellar surface.** With the SV-7 star ($a = 1/2$) the box's faces lie on the
  surface; in the light atmosphere the end-of-step reflux run goes non-finite at $t = 1.41$, the per-stage one
  reaches $t = 2$ ([#88](https://github.com/szaghi/adam/issues/88)). Production codes keep coarse–fine boundaries out
  of dense gas for the same reason (Castro refines everything far below the stellar density, Katz et al. 2016;
  Berger and Colella 1989 on the spurious waves of a resolution jump at a steep gradient).
- **A runtime-AMR star changes its grid once, and needs `conservative-state`** (SV-15). One level above the $32^3$
  base can only refine the 8 inner blocks (0.5 wide); two levels would carry the 2:1 balance to the boundary blocks,
  which is refused; 4-cell blocks are below $2\,n_{gc}$. The first regrid refines the 8 inner blocks and later regrids
  keep that grid. With the default prolongation the collapse then drains a cell's internal energy to negative;
  `conservative-state`, one range limit for the whole state, runs it to the end. On the CPU it also runs the star of
  radius 0.25 regridded every 5 steps and of radius 0.3 every 8, on which the default goes non-finite (steps 53 and
  48); it fails on the radius-0.22 star regridded every 5 steps and the radius-0.25 star every 8 (3.5 and 4.0 cells
  per radius on the base), which the non-conservative `linear` prolongation runs to the end.
- **Refine before the collapse: the Jeans criterion** (SV-16). A star prolonged mid-collapse is the fragile case; a
  Jeans marker refines it at the start instead, and the collapse never regrids. With 8 cells per Jeans length the
  radius-0.22 star regridded every 5 steps, which no conservative prolongation survives, reaches $t = 1.5$ with the
  default one, and so do the radius-0.25 star regridded every 8 and the SV-15 star (CPU). For this polytrope the Jeans
  length does not shrink in the collapse ($\lambda_J = \sqrt2\,a$ at $\gamma = 2$), so one refinement at the start
  holds; in an adiabatic collapse $\lambda_J \propto \rho^{(\gamma - 2)/2}$, so for any $\gamma < 2$ it shortens as
  the density rises, and the runtime regrids then follow it.
- **Or carry the entropy: dual energy** (SV-17). The coldest cells of the SV-15 collapse sit at $e/E \approx 10^{-2}$
  two steps before the default run goes non-finite; `euler-de` lifts them to their adiabat when they drop below it, and
  the run ends with the default prolongation. On the CPU it also runs the radius-0.25 star regridded every 5 and every
  8 steps; the radius-0.22 star (3.5 cells per radius) still needs the Jeans marker.
- **Three remedies, one per cause:** `conservative-state` keeps a mid-collapse regrid from building cold fast cells,
  the Jeans marker avoids the mid-collapse regrid, dual energy keeps the cold infall from draining.
