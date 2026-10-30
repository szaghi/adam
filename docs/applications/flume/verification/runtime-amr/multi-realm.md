# Runtime AMR on a forest

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| multirealm legs 7–13 | Euler | [#110](https://github.com/szaghi/adam/issues/110) P1–P5 | CPU, FNL | yes | — |

## Case and scope

The realms of a forest regrid during the time loop, each by its own `[amr]` section. At a regrid step:
- every realm marks its blocks;
- the forest balances the marks across the seams;
- every realm adapts;
- the inter-realm topology (register faces, seam ghost rows, boundary crown, device maps) is rebuilt from the new
  trees.

The [forest cookbook, R12](/cookbook/forest#r12-runtime-amr-on-a-forest) describes the mechanism; the
[forest concepts](/concepts/forest) give the rules.

A `mirror` or `refined` seam keeps one level per seam side. A `coupling = auto` seam is decided from the cells at every
topology build, and its levels may differ along it, face by face (equal cells or 2:1).

The legs extend the [multi-realm runs](../euler/multi-realm) of legs 1–6. Where the blocks of the realms can match a
single-realm run, the union must reproduce that run bit for bit. Where they cannot (2-D realms whose roots are not
square), the legs assert invariants:
- conservation in a closed box;
- the same run on 1 and 2 ranks;
- the same run with the topology rebuilt after every step (`ADAM_FOREST_REBUILD_TOPOLOGY=1`, a test hook);
- a restart across regrids.

What the legs catch:
- **marking across a seam:** the forest fills the seams before the realms mark, so a ghost-reading marker sees the
  peer realm (leg 7);
- **the cross-realm balance:** a seam turning 2:1 and back as Sod's waves cross it (leg 8). The balance compares cell
  sizes, not levels, so equal cells may sit at different levels (leg 9);
- **refinement criteria the program registers** (leg 10);
- **mixed levels along a seam:** 2:1 faces in both directions on one seam (leg 11), and a refined patch moving along a
  seam (leg 12). The balance counts the neighbours that touch at an edge or a corner, since the seam slab's edge and
  corner ghosts read them; leg 12's first runs found that gap;
- **restart across regrids** (leg 13).

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Regrid of a forest | the realm phases run by the forest: mark, balance, adapt, rebuild, sweep after sweep | `regrid_realms`, `balance_seam_marks`, `balance_mixed_pair`, [`adam_forest_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_forest_object.F90) |
| Mixed seam | per-face register faces and checks | `check_mixed_seam`, `leaf_peers`, `seam_register_faces`, [`adam_forest_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_forest_object.F90) |
| Markers | Löhner on the density (legs 8, 9, 12), primitive boxes (legs 7, 10, 11), a registered criterion (leg 10) | `mark_runtime`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Seam cadence | β (`stage_coincident`) in every leg | [cookbook R4](/cookbook/forest#r4-seam-cadence-α-or-β) |

## Initial and boundary conditions

| Leg | Problem | Forest | Reference | Ranks |
|---|---|---|---|---|
| 7 | Sod, leg 2's inputs, `[amr] frequency = 4` in realm 2 (a box marker: the regrids keep the grid) | `mirror` | `sod-amr` | 2 |
| 8 | Sod split at the diaphragm, Löhner (refine 0.3, derefine 0.1), a regrid every 4 steps, levels 1–2 in the realms (2–3 in the single run) | `auto` | the single Sod with the same marker; `seam_ghost_fill = injection` and tricubic | 2 |
| 9 | leg 5's forest (misaligned blocks), equal cells at levels 2–3 and 1–2, Löhner | `mirror` | conservation | 2 |
| 10 | a registered criterion (`mode = 6`) voting as the box marker, on one realm and on leg 2's forest | `mirror` | the built-in box marker | 2 |
| 11 | leg 8's grids, a box next to the seam in each realm over half of it (2:1 faces in both directions), a regrid every 4 steps | `auto` | the single realm with the same boxes | 2 |
| 12 | the isentropic vortex (radius 0.07) at (0.47, 0.3), free stream (0, 1) along the seam, $[0, 1]^2 \times [0, 0.5]$ split at $x = 0.5$, octree $8 \times 16 \times 8$ cells per block, levels 1–3, Löhner, a regrid every 5 steps, $t = 0.3$; y periodic, x inviscid walls, z extrapolated with $w = 0$ (a closed box) | `auto` | conservation, ranks, hook | 1, 2 |
| 13 | leg 12's case, 195 steps against 155 + restart + 40 (the restart at a regrid that coarsens and leaves 2:1 faces on the seam; a regrid into a new grid at step 190) | `auto` | the uninterrupted run | 2 (and 1: refused) |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/multirealm/check.sh` | runs the legs (`--legs leg7 ... leg13`; `FLUME_EXE` selects the backend; leg 10 needs `exe/test_flume_user_criterion`) |
| `src/tests/flume/verification/multirealm/make_vortex_forest.py` | the forest of legs 12 and 13 (`--it-save N` for the figure) |
| `src/tests/flume/verification/multirealm/multirealm_oracle.py` | the union against the single-realm run cell by cell (tolerance 0 unless stated), and the summed conservation histories |
| `src/tests/flume/verification/make_doc_runs.sh <runs> multirealm-amr`, `make_doc_figures.py --only multirealm-amr` | re-runs leg 12's case with a checkpoint every 20 steps and draws the figure (about 6 min) |

```bash
cd src/tests/flume/verification/multirealm
./check.sh --legs leg2,leg7,leg8,leg9,leg10,leg11,leg12,leg13          # CPU (leg 7 compares against leg 2's run)
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh --legs leg8,leg9,leg11,leg12,leg13    # FNL
```

## Results

![Runtime AMR on a forest](/flume/multirealm-amr.png)

*Leg 12's case, CPU, 2 ranks, a checkpoint every 20 steps. Top: the density and the blocks of both realms at three
steps; the seam is the dashed line at $x = 0.5$, and the refined patch moves along it with the vortex. Bottom left:
the mass of each realm and the drift of their sum. Bottom right: the blocks per level in each realm.*

Measured on the development box (WSL), 2026-10-29:

| Leg | Check | Bound | CPU | FNL |
|---|---|---|---|---|
| 7 | union against `sod-amr` with a regridding realm | bitwise | bitwise | — |
| 8 | union against the single run, `seam_ghost_fill = injection` | bitwise | bitwise | bitwise |
| 8 | the same, tricubic fill | $10^{-12}$ | $1.7\cdot10^{-13}$ | $1.3\cdot10^{-13}$ |
| 9 | mass of the two realms, the seam stays a mirror | $10^{-12}$ | $2.0\cdot10^{-16}$ | $4.0\cdot10^{-16}$ |
| 10 | registered criterion against the built-in marker, one realm and the forest | bitwise | bitwise | — (CPU driver) |
| 11 | union against the single run, injection; and with the topology rebuilt every step | bitwise | bitwise | bitwise |
| 11 | the same, tricubic fill | $10^{-12}$ | $1.5\cdot10^{-13}$ | $9.3\cdot10^{-14}$ |
| 12 | mass of the two realms over 519 steps | $10^{-12}$ | $1.1\cdot10^{-16}$ | $1.1\cdot10^{-16}$ |
| 12 | 1 rank and the rebuilt topology against 2 ranks | bitwise | bitwise | bitwise |
| 12 | the seam mixed, the grid refined and coarsened | at least one each | the seam mixed at 25 topology builds; 15 realm regrids at 13 steps, 7 refining and 9 coarsening (one both) | the same |
| 13 | restarted against uninterrupted: fields, histories | bitwise, identical | bitwise, identical | bitwise, identical |
| 13 | restart on fewer ranks than saved | refused | refused | refused |

The tricubic legs differ from the single realm only next to the 2:1 faces. One level lower, the realms' coarse
blocks are larger in y and z, which moves the 4-point footprints of their 2:1 faces within the block; the same
constant-in-y data is then summed with other weights.

## Findings

- **Edge and corner neighbours (leg 12).** Leg 12's first octree runs stopped at step 10. A tangential edge ghost of a
  realm-1 seam block at level 3 fell on a realm-2 block at level 1, four times coarser. The balance compared only
  leaves whose faces overlap, while the seam slab's edge and corner ghosts also read the diagonal neighbours. It now
  counts neighbours that touch at an edge or a corner, including across a tangential axis periodic in both realms.
- **A closed box for a conservation bound.** Leg 12 was first designed with extrapolated x faces and drifted by
  $2.7\cdot10^{-6}$: the vortex tail carried mass through the open boundaries. The net boundary flux of the final
  cells matched the measured mass change. A drift bound needs a closed box; inviscid walls give one.
- **Restart and cadence (leg 13).** Every leg here uses β seams. When this page was written an `end_of_step` (α) seam
  was filled into the committed state only, never into the Runge–Kutta stage buffers, so a restarted α forest
  differed at the seam from its first step, and α and β differed by 18–45 % next to the seam after four steps. Fixed
  by [#113](https://github.com/szaghi/adam/issues/113): `multirealm/check.sh` leg 14 measures α against β and
  restarts an α forest bitwise ([cookbook R4](/cookbook/forest#r4-seam-cadence-α-or-β)).
- **Inter-realm 2:1 faces need the same tree kind in both realms and a seam normal to a refined axis.** Quadtrees
  and binary trees are accepted since [#71](https://github.com/szaghi/adam/issues/71) (static 2:1 seams checked by
  `multirealm/check.sh` leg 15; this page's runtime legs run on octrees).
