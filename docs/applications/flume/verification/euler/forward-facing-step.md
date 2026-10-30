# Forward-facing step: a three-realm forest on quadtrees

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| EV-step | Euler | — ([#46](https://github.com/szaghi/adam/issues/46), quadtree AMR; [#66](https://github.com/szaghi/adam/issues/66)) | CPU, FNL (`trees` leg); CPU only (`full` leg) | yes | — |

## Case and scope

The Mach 3 wind tunnel with a step of Woodward & Colella (1984, *J. Comput. Phys.* 54, 115–173, §IV b): the gas enters
from the left at Mach 3, hits a step of height 0.2, and forms a bow shock, a Mach reflection off the top wall and a
re-reflection down the tunnel. The tunnel is L-shaped, so FLUME builds it as **three rectangular realms** glued by two
1:1 mirror seams, and the step is body-fitted: its two faces are physical walls of the realms, with no immersed
boundary.

```
 y = 1   +--------+-----------------------------+   wall
         |   B    |              C              |
  inflow |        | seam                        | outflow
 y = 0.2 +--------+   +------ wall (step top) --+
         |   A    | wall (step front)
 y = 0   +--------+
       x = 0    x = 0.6                       x = 3
```

The case puts together inter-realm mirror seams (two, on different axes, misaligned), intra-realm quadtree 2:1 seams
(refined boxes off the inter-realm seams) and the walls, inflow and outflow of a supersonic tunnel. What it catches:
- **a quadtree that differs from the octree** ([#46](https://github.com/szaghi/adam/issues/46)): the same forest on
  quadtrees (`nk = 1`) and on octrees with a null $z$ axis must agree column by column to round-off;
- **a $z$ dependence** in a 2-D flow: every run must be $z$-invariant;
- **loss of positivity** at the seams, the corner and the reflected shocks;
- **seam defects** that show as structure in the quadtree–octree difference map at the inter-realm or 2:1 seams.

::: warning The refinement is static
The boxes are placed by hand on the features of the reference solution at $t = 4$, and the grid is frozen after
initialisation: FLUME had no runtime AMR when this case was written (milestone M5). A solution-driven marker could not
place them either, since the tunnel starts uniform. The boxes therefore follow the flow only once it has settled, and
the case verifies the seam machinery, not adaptivity.

The realms of a forest now regrid at runtime ([#110](https://github.com/szaghi/adam/issues/110)), but this forest's
seams are `mirror` seams between blocks that do not line up, and a mirror seam keeps one level per seam side. A
Löhner marker tried on the octree at $N = 80$ (2026-10-29) refined whole seam sides at its first regrid, and every
later sweep balanced to no mark, so the grid stayed fixed. A `coupling = auto` seam varies its levels face by face but
needs the blocks of a 2:1 face to nest, which realm A (0.2 high) cannot offer realm B on an octree at this size.
:::

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | compressible Euler, ideal gas, $\gamma = c_p/c_v = 3.5/2.5 = 1.4$ | [models](../../models#compressible-euler) |
| Space | WENO5 (`weno-u-5`) flux splitting on characteristic variables; no positivity limiter (refused on multi-realm runs) | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Time | SSP Runge–Kutta, 5 stages (`runge-kutta-ssp-54`), CFL 0.5 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Forest | three realms, two 1:1 mirror seams (A–B at $y = 0.2$, B–C at $x = 0.6$), β cadence, misaligned seam blocks ([#51](https://github.com/szaghi/adam/issues/51)) | `populate_inter_realm_topology`, `build_seam_rows`, [`adam_forest_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_forest_object.F90); [numerics](../../numerics#multi-realm-runs) |
| AMR | one more level on a static box per realm (geometric box marker, block centroids), quadtree `ratio = 4`, reflux on | `mark_by_geometry`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90); [numerics](../../numerics#adaptive-mesh-refinement) |
| Thin-block guard | blocks thinner than $2\,n_{gc}$ along a refined, non-null axis refused ([#66](https://github.com/szaghi/adam/issues/66)) | `check_amr_block_cells`, [`adam_realm_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_realm_object.F90) |
| Boundaries | inflow, extrapolation (outflow), inviscid walls; realm edges compose the conditions of their faces | `set_boundary_conditions`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `realm_edge_face`, `realm_edge_donor`, [`adam_flume_bc_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_bc_object.F90); [boundary conditions](../../boundary-conditions#edges-and-corners-fec-6) |
| Corner | untreated (Woodward and Colella reset the entropy near it) | — |

## Initial and boundary conditions

| Realm | Extent | Faces |
|---|---|---|
| A | $[0, 0.6] \times [0, 0.2]$ | inflow, wall (bottom), wall (step front), seam with B (top) |
| B | $[0, 0.6] \times [0.2, 1]$ | inflow, seam with A (bottom), seam with C (right), wall (top) |
| C | $[0.6, 3] \times [0.2, 1]$ | seam with B (left), wall (step top), wall (top), outflow (extrapolation) |

The realms share the base cell $1/N$. B and C have $N/10$ blocks per axis, narrow enough to refine close to the seams;
A has $N/40$, so that its blocks keep 8 cells along $y$: a refined block needs at least $2\,n_{gc} = 6$ cells along
every active axis.

| Item | `trees` leg | `full` leg |
|---|---|---|
| Base cell | $1/40$ | $1/80$ (the coarse grid of Woodward and Colella) |
| Blocks, quadtree | A: 1 block of $24 \times 8$; B: $4 \times 4$ of $6 \times 8$; C: $4 \times 4$ of $24 \times 8$ | A: $2 \times 2$ of $24 \times 8$; B: $8 \times 8$ of $6 \times 8$; C: $8 \times 8$ of $24 \times 8$ |
| Trees | quadtree (`nk = 1`) and octree with a null $z$ axis (`nk = 4`, 16 in A) | quadtree |
| Refined boxes (one level, block centroids) | A: none (its single block touches its seam); B: $x \in [0.3, 0.45]$, $y \in [0.4, 0.8]$; C: $x \ge 1.2$ | A: $y \le 0.1$; B: $x \in [0.3, 0.525]$, $y \in [0.3, 0.7]$; C: $x \ge 0.9$ (cells $1/160$) |
| Initial state | the inflow state everywhere: $(\rho, u, v, p) = (1.4, 3, 0, 1)$ (sound speed 1, Mach 3) | the same |
| Final time | $t = 0.5$ (440 steps) | $t = 4$ (7545 steps); plus the uniform $1/80$ control (no box) |
| Ranks | 2 | 2 |

The boxes of the `full` leg, at $N = 80$:
- **A**, $y \le 0.1$: the bow shock standing on the lower wall.
- **B**, $x \in [0.3, 0.525]$, $y \in [0.3, 0.7]$: the bow shock bending up towards the Mach stem. The column
  $x > 0.525$ touches the B–C seam.
- **C**, $x \ge 0.9$: all of the realm except the block column on the seam. That covers the reflected shock reaching
  the step top, the slip line and the re-reflection off the top wall.

A refined block may not touch an inter-realm seam (C's seam is only its face $x = 0.6$; its face $y = 0.2$ is the step
top, a wall). The strongest features, the Mach stem and the triple point at $x \approx 0.6$, therefore sit on the B–C
seam and stay on the base grid. Refining across that seam would take the 2:1 inter-realm coupling
(`coupling = refined`, [#52](https://github.com/szaghi/adam/issues/52)), available on quadtrees too since
[#71](https://github.com/szaghi/adam/issues/71) but not used by this case.

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/step/make_step.py <dir> --tree quad\|oct --cells N --refine --time-max T` | writes the manifest `step.ini` and the realms `step-A.ini`, `step-B.ini`, `step-C.ini` |
| `src/tests/flume/verification/step/check.sh` | `--legs trees` (default) and `full`; `--keep` keeps the checkpoints for the figures; checks positivity of every run |
| `src/tests/flume/verification/mhd/quadtree/quadtree_oracle.py` | the $z$ invariance and the quadtree against the octree per $(x, y)$ column, momentum scaled as one vector |
| `src/tests/flume/verification/make_doc_figures.py` (`fig_step`, `fig_step_trees`) | the figures, from the kept work directories of `step/` |

```bash
cd src/tests/flume/verification/step && ./check.sh                    # trees leg, CPU (about 10 min)
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh           # trees leg, FNL
./check.sh --legs trees,full --keep                                   # both legs, checkpoints kept
python3 ../make_doc_figures.py --out ../../../../../docs/public/flume --only step step-trees
```

## Results

![Step forest at t = 0.5 on quadtrees and octrees](/flume/step-trees.png)

*`trees` leg, $t = 0.5$, base $1/40$: density on the quadtree and on the octree with the blocks outlined, and the
quadtree–octree difference per $(x, y)$ column.*

| `trees` leg | Bound | Octree, null $z$ (CPU / FNL) | Quadtree, `nk = 1` (CPU / FNL) |
|---|---|---|---|
| $z$ spread of a column | $10^{-10}$ | $6.1\cdot10^{-13}$ / $3.1\cdot10^{-13}$ | 0 / 0 |
| against the octree, per $(x, y)$ column | $10^{-10}$ | — | $6.5\cdot10^{-13}$ / $2.7\cdot10^{-13}$ |
| min $\rho$, min $p$ | positive, finite | 0.469, 0.715 (both) | 0.469, 0.715 (both) |

![Woodward-Colella step at t = 4](/flume/step.png)

*`full` leg, $t = 4$: density with the blocks outlined (static boxes, inter-realm seams dashed), Mach number, and the
zoom on the Mach stem and C's 2:1 seam at $x = 0.9$, refined run against the uniform $1/80$ control.*

| `full` leg (CPU) | Bound | Result |
|---|---|---|
| min $\rho$, min $p$, refined run | positive, finite | 0.344, 0.370 |
| uniform control | positive, finite | pass |

The picture has the features of the reference of Woodward and Colella:
- the bow shock standing at $x \approx 0.3$ on the lower wall;
- its Mach reflection off the top wall at $x \approx 0.6$, with the triple point near $(0.6, 0.78)$ and the slip line
  leaving it along $y \approx 0.8$;
- the reflected shock striking the step top at $x \approx 1.2$ in a Mach reflection;
- its re-reflection reaching the top wall at $x \approx 2.4$.

Related cases: [Multi-realm runs](./multi-realm) (the seam legs on Sod), [Quadtree AMR](../mhd/quadtree-amr) (the
quadtree on MHD), [Ghost probe](../infrastructure/ghost-probe) (the step forest's edge ghosts).

## Findings

- The trees order the restriction sums differently (8 fine cells over two identical $z$ layers, against 4), so they
  agree to round-off, not bit for bit; the octree is $z$-invariant to round-off only, its tricubic seam weights
  depending on the $z$ sub-position. The difference map carries no structure at the seams, inter-realm or 2:1.
- The `trees` leg takes about 10 min on the CPU, 9 of them for the octree (its blocks split along $z$ too). The `full`
  leg is not default: about 41 min on the CPU, plus 13 min for the uniform control.
- The low-Mach layer along the step top is the numerical boundary layer of the untreated corner, the artifact
  Woodward and Colella remove with their entropy fix.
- The zoom compares the refined run with the uniform $1/80$ run. The Mach stem and the triple point sit on the B–C
  seam, on the base grid in both runs, and agree. C's 2:1 seam at $x = 0.9$ is crossed by the reflected shock and the
  slip line without a kink. Downstream of it, the box resolves the shock pattern near the step top ($x \approx 1.2$)
  more sharply than the uniform grid, which smears it into a single blob.

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
