# Forest seams: the inter-realm contract

| ID | Group | Milestone | Type | Backends | Figure | Regression |
|---|---|---|---|---|---|---|
| multirealm legs 1–6, 15; MV-14 | Forest | [#37](https://github.com/szaghi/adam/issues/37) (follow-up) · legs 1–2; [#40](https://github.com/szaghi/adam/issues/40), [#51](https://github.com/szaghi/adam/issues/51), [#52](https://github.com/szaghi/adam/issues/52), [#54](https://github.com/szaghi/adam/issues/54) (follow-up) · legs 3–5; M2 ([#41](https://github.com/szaghi/adam/issues/41)) · P7c (MV-14); M6 ([#68](https://github.com/szaghi/adam/issues/68)) · leg 6; [#71](https://github.com/szaghi/adam/issues/71) (follow-up) · leg 15 | C | CPU, FNL | yes | `sod-2realm`, `sod-2realm-blocks`, `sod-amr`, `sod-amr-2realm`, `sod-amr-refined` |

## Case and scope

A forest manifest glues realms through inter-realm seams. The case splits the [Sod](./sod) tube into realms whose union
is the cell set of a single-realm run, and requires the union to reproduce that run **bit for bit**: interior cells
compared by their centres with tolerance 0. With the seam ghosts filled at every Runge–Kutta stage (β,
`coupling_cadence = stage_coincident`) a 1:1 mirror seam is a block interface like any other, and a 2:1 refined seam
must use the formulas of the intra-realm AMR seams; any difference in the seam fill, the seam exchange across ranks or
the inter-realm reflux register breaks the bitwise identity. The split at the diaphragm makes the seam carry the
Riemann problem from the first step.

The contract is built issue by issue: the 1:1 mirror seam ([#37](https://github.com/szaghi/adam/issues/37), legs 1–2),
the seam across ranks ([#40](https://github.com/szaghi/adam/issues/40), leg 3), the 2:1 inter-realm seam
([#52](https://github.com/szaghi/adam/issues/52), leg 4), misaligned seam blocks ([#51](https://github.com/szaghi/adam/issues/51),
leg 5), the high-order seam operators of M6 (leg 6), the seam peer lookup on quadtrees
([#54](https://github.com/szaghi/adam/issues/54), MV-14 leg 6) and the 2:1 seams of quadtrees
([#71](https://github.com/szaghi/adam/issues/71), leg 15). MV-14 carries the same contract for the nine MHD fields.

What the legs catch:
- **seam fill and inter-realm reflux** at a mirror seam (leg 1): a boundary routine that does not skip the forest's
  `BC_SEAM` crown rows, a fine side of the inter-realm register restricted like an intra-realm AMR seam, seam and
  boundary maps missing on the device;
- **both kinds of register face in one realm**, the intra-realm 2:1 face and the fine side of the inter-realm seam
  (leg 2): the forest must compose the two registrations;
- **rank-local seam machinery** ([#40](https://github.com/szaghi/adam/issues/40)): the seam ghosts and the reflux
  register must cross ranks (leg 3); a rank-local register ends normally with a wrong solution;
- **a 2:1 inter-realm seam** ([#52](https://github.com/szaghi/adam/issues/52)) that interpolates (coarse to fine) or
  restricts (fine to coarse) differently from an intra-realm 2:1 face (leg 4), and seam operators of higher order
  ([#68](https://github.com/szaghi/adam/issues/68) C2, C3, D1, D2) that differ across the seam (leg 6);
- **seam blocks that do not line up** ([#51](https://github.com/szaghi/adam/issues/51)): the coverage of the register
  faces (leg 5), and where the overlaps land along a face (MV-14 leg 5);
- **a seam peer lookup that assumes an octree** ([#54](https://github.com/szaghi/adam/issues/54)): the closest leaf of
  the peer's tree on binary trees, quadtrees and octrees (MV-14 leg 6, `test_tree_closest_block`);
- **the MHD state across a seam and a restart** (MV-14): every forest and restart channel widened to the nine fields;
- **invalid couplings**: a 2:1 pair declared `mirror`, and the positivity limiter on a multi-realm run, must be refused
  at initialization.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations, scheme | as [Sod](./sod): Euler, $\gamma = 1.4$, WENO5 flux splitting on characteristic variables, SSP-33, CFL 0.5 | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Forest | manifest, inter-realm topology, coupling checks (a mirror pair needs equal cells, a refined pair ratio 2) | `populate_inter_realm_topology`, `check_couplings`, `register_inter_realm_seams`, [`adam_forest_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_forest_object.F90); [numerics](../../numerics#multi-realm-runs) |
| Seam fill | β cadence, every Runge–Kutta stage; cross-rank rows | `fill_seam_from_peer_forest`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `seam_fill`, [`adam_seam_exchange.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_seam_exchange.F90); [boundary conditions](../../boundary-conditions#inter-realm-seam) |
| 2:1 seam operators | tricubic fill and mean restriction (default); `seam_restriction = lagrange6, weno6`, `seam_ghost_fill = quintic, weno6, quintic-mp` (leg 6) | `seam_interpolate_tricubic`, `seam_interpolate_quintic`, `seam_restriction_weno6`, `seam_interpolate_weno6`, `seam_interpolate_quintic_mp`, [`adam_seam_interpolation_library.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_seam_interpolation_library.F90) |
| Inter-realm reflux | the register spans the realms; misaligned faces through overlaps | `share_flux_register`, `apply_reflux_corrections`, [`adam_forest_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_forest_object.F90); `accumulate_fine_overlaps`, [`adam_flux_register_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_flux_register_object.F90) |
| Seam peer lookup | the closest leaf of the peer's tree ([#54](https://github.com/szaghi/adam/issues/54)) | `get_closest_block`, [`adam_tree_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_tree_object.F90) |
| Refusals | limiter, `[reference]`, gravity equilibrium, self-gravity on multi-realm runs | `initialize_flume`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90), [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

## Initial and boundary conditions

Every leg is the Sod problem of [Sod](./sod): $(\rho, u, p) = (1, 0, 1)$ for $x \le 0.5$, $(0.125, 0, 0.1)$ beyond,
`extrapolation` on the physical faces, $t = 0.2$; the seam faces are declared in the manifest
(`[forest.topology.face_1]`, `coupling_cadence = stage_coincident`).

| Leg | Split | Realms | Coupling | Reference | Ranks |
|---|---|---|---|---|---|
| 1 | $x = 0.5$ | each 64 blocks of $25 \times 4 \times 4$ (the `sod-x` cells, $\Delta x = 0.005$) | `mirror` | `sod-x`, 51200 cells, 174 steps | 2 |
| 2 | $x = 0.5$, realm 2 refined on $x > 0.75$ | `sod-amr` cells: blocks of $24 \times 4 \times 4$ ($\Delta x = 1/192$), the 2:1 face at $x = 0.75$ crossed by the shock at $t \approx 0.14$; the seam blocks unrefined | `mirror` | `sod-amr` (`ni = 48`, $x > 0.75$ refined), 135168 cells, 216 steps | 2 |
| 3 | $z = 0.5$ | the leg-1 split along $z$; on $N \ge 2$ ranks the two sides of the seam sit on different ranks (np 2: 3468 seam rows each way, none local) | `mirror` | `sod-z` | 1 and $N$ |
| 4 | $x = 0.75$, the 2:1 face of `sod-amr`; also along $z$, where every 2:1 seam row crosses ranks | coarse realm $x < 0.75$: 4 blocks per axis, `ni = 36` ($\Delta x = 1/192$); fine realm $x > 0.75$: 8 blocks per axis, `ni = 12` ($\Delta x = 1/384$) | `refined` | `sod-amr` (135168 cells, 216 steps), `sod-amr-z` | 2, 3, 4 |
| 4, refused | as leg 4 | the pair declared `mirror`; the pair with `positivity_limiter = cell` | — | must stop at initialization | 2 |
| 5 | $x = 0.5$ and $z = 0.5$ | the same cells, seam blocks not lined up: one realm on 2 blocks per axis of twice the cells ($50 \times 8 \times 8$; along $z$ the register side, $8 \times 8 \times 50$), each of its seam blocks facing $2 \times 2$ blocks of the other | `mirror` | `sod-x`, `sod-z` | 2, 4 |
| 6 | leg 4, along $x$ and $z$ | every input with `[amr] seam_restriction = lagrange6` or `weno6`, or `seam_ghost_fill = quintic`, `weno6` or `quintic-mp`; the null axes widened to 8 cells (the 6-point footprint needs them) | `refined` | `sod-amr`, `sod-amr-z` with the same option | 2 |
| 15 | leg 4 on quadtrees (`[amr] ratio = 4`, $z$ not refined), along $x$, then turned to $y$ ([#71](https://github.com/szaghi/adam/issues/71)) | leg 4's realms: along $z$ the 2:1 seam ghosts copy the coincident cell, a coarse ghost is the mean of $2 \times 2$ fine cells, a register face holds two fine blocks | `refined` | `sod-amr` on a quadtree (21504 cells, 216 steps), and its $y$ twin | 2 |
| 15, $z$ | leg 4's $z$ split on quadtrees; leg 3's $z$ split on quadtrees (`nk = 25`, odd) | a seam normal to the axis the quadtree does not refine | `refined`; `auto` | must stop at initialization; resolves to `mirror`, against `sod-z` on a quadtree | 2 |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/multirealm/sod-2realm{,-r1,-r2}.ini` | leg 1: the manifest and the two realms |
| `src/tests/flume/verification/multirealm/sod-amr.ini`, `sod-amr-2realm{,-r1,-r2}.ini` | leg 2: the single-realm reference and the split |
| `src/tests/flume/verification/multirealm/sod-2realm-z{,-r1,-r2}.ini` | leg 3 |
| `src/tests/flume/verification/multirealm/sod-amr-refined{,-r1,-r2}.ini`, `sod-amr-z.ini`, `sod-amr-refined-z{,-r1,-r2}.ini` | legs 4 and 6 |
| `src/tests/flume/verification/multirealm/sod-2realm-blocks{,-r1,-r2}.ini`, `sod-2realm-blocks-z{,-r1,-r2}.ini` | leg 5 |
| `src/tests/flume/verification/sod/sod-{x,z}.ini` | the single-realm references of legs 1, 3 and 5 |
| `src/tests/flume/verification/multirealm/check.sh` | runs these legs and the refusals (`--np N`, `--legs LIST`: `leg1` .. `leg6` and `leg15`, the refusals with `leg4`; leg4 needs leg2, leg5 needs leg1 and leg3; leg 15 writes its quadtree inputs from leg 4's and leg 3's; `FLUME_EXE` selects the backend) |
| `src/tests/flume/verification/make_doc_runs.sh <runs> multirealm`, `make_doc_figures.py --only multirealm` | re-runs legs 1 and 4 and their single-realm twins with checkpoints kept (4 runs, about 4 min) |
| `src/tests/flume/verification/multirealm/multirealm_oracle.py` | the union against the single-realm run, cell by cell, tolerance 0, and the summed conservation histories |
| `src/tests/flume/verification/mhd/multirealm/check.sh` | the MHD twin (MV-14; `--legs LIST`: `mv14-1` .. `mv14-6`, mv14-2 needs mv14-1) |
| `src/tests/amr/test_tree_closest_block.F90` | the unit test of the tree lookup ([#54](https://github.com/szaghi/adam/issues/54)) |

```bash
cd src/tests/flume/verification/multirealm && ./check.sh                 # CPU, np 2 (leg 3 also on np 1)
./check.sh --np 4                                                         # the same legs on 4 ranks
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh               # FNL
```

## Results

![Multi-realm Sod unions against the single realm](/flume/multirealm.png)

*Sod at $t = 0.2$, CPU, 2 ranks: leg 1 (two realms glued by a 1:1 mirror seam at the diaphragm) over `sod-x`, and leg 4 (a coarse and a fine realm glued by the 2:1 `refined` coupling at $x = 0.75$) over the single-realm `sod-amr`; the legend gives the largest cell difference, 0 in both (bitwise). The rotor split on the quadtree (MV-14 leg 6) has no figure: it is bitwise too, and a map of zero carries nothing.*

| Leg | Check | Bound | CPU | FNL |
|---|---|---|---|---|
| 1 | union against `sod-x`, 51200 cells | bitwise | bitwise | bitwise |
| 1 | summed conservation histories | round-off | within $2\cdot10^{-13}$ | |
| 2 | union against `sod-amr`, 135168 cells, 216 steps | bitwise | bitwise | bitwise |
| 3 | union against `sod-z` | bitwise | bitwise on 1, 2, 3, 4 ranks | bitwise on 1, 2 ranks |
| 4 | union against `sod-amr` and `sod-amr-z` | bitwise | bitwise on 2, 3, 4 ranks | bitwise on 2, 3, 4 ranks |
| 4 | `refined` pair declared `mirror`; limiter on | refused | refused | refused |
| 5 | unions against `sod-x` and `sod-z` | bitwise | bitwise on 2, 4 ranks | bitwise on 2 ranks |
| 6 | unions, every restriction and fill option | bitwise | bitwise | bitwise |
| 15 | quadtree unions against the single quadtree, $x$ and $y$ splits, 21504 cells | bitwise | bitwise | bitwise |
| 15 | summed conservation histories | round-off | within $2\cdot10^{-16}$ | within $4\cdot10^{-16}$ |
| 15, $z$ | `refined` across the null axis; `auto` across it | refused; bitwise (12800 cells) | refused; bitwise | refused; bitwise |

**MV-14, the MHD twin** (`mhd/multirealm/check.sh`, M2 P7c; np 2, CPU):

| Leg | Check | Result |
|---|---|---|
| 1 | RJ2a (MV-4, GLM, $N = 256$) split at the diaphragm by a mirror seam, against the single realm | bitwise on all nine fields over 297 steps; summed conservation within $1.7\cdot10^{-16}$ |
| 2 | restart round trip of the two-realm run (148 steps, restart, to $t = 0.2$) | bitwise on all nine fields; residuals, conservation and div(B) histories byte-identical |
| 3 | restart round trip of the MV-11 box (3-D, 2:1 AMR, reflux, damped GLM) | bitwise on all nine fields, with $\max\lvert\psi\rvert = 1.7\cdot10^{-3}$ |
| 4 | RJ2a with GLM on the leg-4 cells, `coupling = refined` ([#52](https://github.com/szaghi/adam/issues/52)) | bitwise on all nine fields over 349 steps |
| 5 | the MHD rotor on an octree, split at $x = 0.5$ with misaligned seam blocks, both ways, 50 steps ([#51](https://github.com/szaghi/adam/issues/51)) | bitwise on all nine fields; a deliberately misplaced overlap drives the density negative at step 2 |
| 6 | the rotor on its own quadtree, blocks lined up and not ([#54](https://github.com/szaghi/adam/issues/54)) | bitwise on all nine fields |

Related cases: [Sod](./sod) (the reference runs), [Forward-facing step](./forward-facing-step) (three realms on
quadtrees), [MHD rotor](../mhd/rotor) (the MV-14 split), [Ghost probe](../infrastructure/ghost-probe) (the seam edge
ghosts).

## Findings

- **A 1:1 mirror seam with β cadence is a block interface like any other:** the union reproduces the single realm bit
  for bit, across ranks and with misaligned seam blocks. The seam fill must run at every Runge–Kutta stage
  (`coupling_cadence = stage_coincident`).
- **A `refined` seam reproduces an intra-realm 2:1 face bit for bit,** with the default operators and with every
  6-point restriction and fill of M6 (leg 6): the fine owner packs the 6-point restriction codes in the restrict rows,
  the coarse owner the anchor shifts of the fill; the shock crosses the seam, so a different operator on either side
  would show.
- Sod is 1-D, so legs 4 and 5 test the coverage of the register faces, not where the overlaps land along a face; the
  rotor of MV-14 leg 5 does. In 1-D $\psi$ stays at round-off, so MV-14 legs 2 and 4 test $\mathbf{B}$, not the
  cleaning; MV-14 leg 3 carries a non-zero $\psi$ through a restart.
- The seam peer lookup takes the closest leaf of the peer's tree by the Morton code of the tree's own dimension; on a
  single-block tree it returns the root. `test_tree_closest_block` pins it on binary trees, quadtrees and octrees.
- **Refused at initialization:** a 2:1 pair declared `mirror`; a `refined` seam normal to the axis a quadtree does not
  refine; the positivity limiter, `[reference]`, the gravity equilibrium and self-gravity on a multi-realm run. The
  inter-realm seam faces carry no limiting factor (the 2:1 seam flux synchronisation of
  [#50](https://github.com/szaghi/adam/issues/50) is intra-realm).
- Runtime AMR across the seams of a forest is [runtime AMR on a forest](../runtime-amr/multi-realm)
  ([#110](https://github.com/szaghi/adam/issues/110), legs 7–13).
