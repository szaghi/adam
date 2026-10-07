# Forest cookbook

This page is the task-oriented companion of the [forest reference](/guide/forest). Each recipe is a configuration ADAM
runs today:
- a picture of what it builds;
- the manifest and realm keys it needs;
- what the machinery does with them;
- what it refuses;
- how to tell from the log that it did the right thing;
- the test that keeps it working.

The figures are drawn to scale by `scripts/make_forest_figures.py` (standard library only; run it to regenerate
`docs/public/forest/*.svg`).

::: tip What "verified" means here
Every recipe names its tests. Unless stated otherwise the oracle is **bitwise**: the forest run must reproduce the
single-realm run on the same cells to the last bit, on the CPU and the FNL (OpenACC) backends. Where a test holds only
to round-off or to a tolerance, the recipe says so and why.
:::

## Vocabulary

![A two-realm forest](/forest/anatomy.svg)

| Term | Meaning |
|---|---|
| **realm** | One simulation domain with its own INI (grid, numerics, physics), its own block tree and its own partition over the MPI ranks. |
| **forest** | A set of realms advanced together with one global time step. A plain INI is a one-realm forest. |
| **manifest** | The small INI with a `[forest]` section that lists the realms and the seams between them. |
| **seam** | A face where two realms touch (an *inter-realm* seam, declared in the manifest), or a 2:1 face inside one realm's tree (an *intra-realm* AMR seam, created by refinement). |
| **face pair** | One `[forest.topology.face_N]` entry: `realm_a`/`face_a` glued to `realm_b`/`face_b`. |
| **ghost slab** | The `ngc` layers of ghost cells beyond a realm's seam face, filled from the peer realm. |
| **row** | One seam ghost to fill: a **copy** (same cell size), an **interpolate** (a fine ghost from coarse cells) or a **restrict** (a coarse ghost from fine cells). |
| **flux register** | The per-face store of the coarse-side flux and the sum of the fine-side fluxes; after the step the coarse cells get the Berger-Colella correction. |
| **cadence** | When the seam ghosts are refreshed: before every RK stage (β, `stage_coincident`) or once per step (α, `end_of_step`). |

## Which recipe?

| You want | Recipe |
|---|---|
| one domain, uniform grid | [R1](#r1-one-realm-uniform) |
| one domain, finer cells in a region | [R2](#r2-one-realm-with-amr-intra-realm-2-1-faces) |
| two domains, same cells, glued | [R3](#r3-two-realms-mirror-seam-lined-up-blocks) |
| choose how the seam is synchronised, or mix RK schemes | [R4](#r4-seam-cadence-α-or-β) |
| two domains, same cells, different block sizes | [R5](#r5-mirror-seam-between-blocks-that-do-not-line-up) |
| two domains, one twice as fine as the other | [R6](#r6-refined-seam-a-2-1-jump-between-realms) |
| a seam **and** refinement inside a realm | [R7](#r7-both-families-a-mirror-seam-and-amr-inside-a-realm) |
| any of the above on many ranks | [R8](#r8-many-ranks) |
| restart a forest | [R9](#r9-restart) |
| build a split from an existing single-realm input | [R10](#r10-splitting-an-existing-input) |

---

## R1. One realm, uniform

A plain INI, no manifest: the forest has one realm and takes the fused fast path (`advance_one_step_forest`); no
seam machinery runs.

```ini
[amr]
max_level      = 2
ratio          = 8        ; octree
iu_ref_levels  = 2        ; 4 x 4 x 4 blocks
markers_number = 0
```

- **Blocks.** `iu_ref_levels = L` gives `2^L` blocks per refined axis (x, y, z on an octree; x, y on a quadtree,
  `ratio = 4`), each with `ni x nj x nk` cells.
- **Verified by** every single-realm regression case (FLUME `sod-x`, PRISM `rmf`, ...).

## R2. One realm with AMR: intra-realm 2:1 faces

![Intra-realm AMR](/forest/intra-amr.svg)

A refinement marker creates blocks one level finer; where a level-`l` block meets a level-`l+1` block the tree has a
**2:1 face**. No manifest is involved.

What happens at a 2:1 face:
- **Ghosts.** The fine block's ghosts that overlap the coarse block are **interpolated** with the
  `[amr] seam_ghost_fill` regime; the coarse block's ghosts that overlap fine cells take the **mean of the 2x2x2**
  fine cells.
- **Reflux.** The face is registered (`SEAM_KIND_INTRA_REALM_AMR`); after the step the coarse cells get the
  Berger-Colella correction (FLUME: `[numerics] reflux = .true.`).

| `seam_ghost_fill` | Order | Footprint | Use |
|---|---|---|---|
| `tricubic` (**default**) | q = 4 | 4x4x4 coarse cells | production |
| `restriction-compatible` | q = 2 | 3x3x3 | exact restriction of the interpolant (`R∘P = I`); diagnostics |
| `injection` | 0 | the anchor cell | legacy, comparisons |

**Requirements.** The first two are checked at initialisation:
1. On an **octree** (`ratio = 8`) with a null z axis, `nk >= 4` (the tricubic footprint is 4 cells along every refined
   axis). A **quadtree** (`ratio = 4`) never refines z, so any `nk` works, `nk = 1` included
   ([#46](https://github.com/szaghi/adam/issues/46); the NVF and GMP backends still refuse quadtree AMR).
2. **Even** `ni`/`nj`/`nk` along every refined, non-null axis when 2:1 refinement is possible (markers and
   `max_level > iu_ref_levels`): an odd count is refused ([#39](https://github.com/szaghi/adam/issues/39); it used to
   give NaN silently).
3. At least 4 cells per block on every axis, null axes included, for the coarse-fine interpolation footprint.

**Log.** `forest: registered intra-realm AMR seam faces: +N`.

**PRISM.** A 2:1 face injects an O(h^p) div(B) source that refinement reduces but time does not
([#29](https://github.com/szaghi/adam/issues/29), accepted). The run-time monitor `[IO] seam_divB_tol` turns it into a
warning or a stop; see [the reference](/guide/forest#seam-div-b-guard-rail-io).

**Verified by** FLUME `sod-amr` (a 2:1 face crossed by a shock), `amr-periodic-reflux` (V3: conservation to round-off,
the negative control without reflux, and the odd-cell refusal), MHD MV-11; PRISM `rmf-amr`, `rmf-amr-fd`,
`rmf-amr-fd-pulse`.

## R3. Two realms, mirror seam, lined-up blocks

![Mirror seam](/forest/mirror-aligned.svg)

The manifest (any file name; the apps recognise it by its `[forest]` section):

```ini
[forest]
realms_number = 2

[realm.1]
ini = sod-2realm-r1.ini        ; relative to the manifest's directory

[realm.2]
ini = sod-2realm-r2.ini

[forest.topology]
inter_realm_faces_number = 1

[forest.topology.face_1]
realm_a          = 1
face_a           = +x          ; +x | -x | +y | -y | +z | -z
realm_b          = 2
face_b           = -x
coupling         = mirror      ; the default
coupling_cadence = stage_coincident
```

The realm INIs are complete inputs. For a split of a single-realm domain, each realm keeps the cell size and the
block layout and takes its part of the extent:

```ini
; realm 1, x in [0, 0.5]                     ; realm 2, x in [0.5, 1]
[grid]                                         [grid]
ni     = 25          ; 50 / 2                  ni     = 25
emin_x = 0.0                                   emin_x = 0.5
emax_x = 0.5                                   emax_x = 1.0
[IO]                                           [IO]
output_basename  = sod-2realm-r1               output_basename  = sod-2realm-r2
restart_basename = sod-2realm-r1-restart       restart_basename = sod-2realm-r2-restart
```

What happens:
- **Ghosts.** Every seam ghost is a **copy** of the peer cell with the same centre. The forest builds the rows at
  initialisation from the replicated trees; it fills the realm 1 slab from realm 2 and the realm 2 slab from realm 1,
  edge and corner ghosts included.
- **Register.** One register face per realm_a seam block. On a mirror seam both sides compute the same flux, so the
  correction is exactly zero: the forest prints `forest: reflux max|F_coarse-F_fine_sum|` only when it is not.
- **Physical faces.** A realm's faces that are not on a seam keep their `[bc_*]` boundary conditions.

**Requirements.**
1. The same cell size on both sides, with cell centres that coincide across the seam.
2. The two seam faces cover each other.
3. One cell size along each realm's seam face (no refinement touching the seam face).

**Log.** `forest: realm R seam rows: L local, R received, S sent` (per rank).

**Verified by** FLUME `multirealm/check.sh` leg 1 (Sod) and the regression cases `sod-2realm` (`equivalent_to sod-x`)
and `rj2a-2realm` (MHD, `equivalent_to rj2a-x`); PRISM `rmf-2realm*` and `rmf-2realm-fd-pulse` (div(B) = div(D) = 0
across the seam).

## R4. Seam cadence: α or β

![Seam cadence](/forest/cadence.svg)

| `coupling_cadence` | Ghosts refreshed | Requires | Gives |
|---|---|---|---|
| `stage_coincident` (β) | before every RK stage's residual, from the peer's stage buffer | the same `scheme_time`, RK scheme, stage count and `nv` on both sides (checked) | the seam has the time order of the interior; the split reproduces the single realm (bitwise in every FLUME split) |
| `end_of_step` (α, **default**) | once, after the step, from the peer's committed state | nothing | first-order seam coupling in time; the realms may use different RK schemes (asymmetric K) |

**Choose β** whenever both realms use the same time integrator. It is the cadence under which a split reproduces the
single-realm run, and the only one that keeps a PRISM 1:1 seam divergence-free
([#31](https://github.com/szaghi/adam/issues/31): under α the substage ghosts are stale and div(B) leaks).

**Choose α** for different RK schemes per realm (for example SSP-33 next to SSP-54): the realm with fewer stages idles
the trailing ones.

**Rules.**
1. β on a pair whose realms differ in scheme, stages or `nv` stops at initialisation, naming the field; it is never
   silently downgraded.
2. All face pairs joining the same two realms must declare the same cadence.

**Verified by** PRISM `rmf-2realm` (α), `rmf-2realm-asymK` (α, K = 3 and 5), `rmf-2realm-stagesync` (β, the
single-realm golden within the digest tolerances), `rmf-2realm-fd-pulse` (β div-free, with α as the leaking negative
control). Every FLUME forest test uses β.

## R5. Mirror seam between blocks that do not line up

![Misaligned blocks](/forest/mirror-misaligned.svg)

The two realms keep the **same cells** but use different block sizes along the seam: for example realm 2 with one
block level fewer and twice the cells per block. Nothing changes in the manifest
([#51](https://github.com/szaghi/adam/issues/51)).

```ini
; realm 1: 4 blocks per axis                ; realm 2: 2 blocks per axis, the same cells
ni = 25   nj = 4   nk = 4                      ni = 50   nj = 8   nk = 8
iu_ref_levels = 2                              iu_ref_levels = 1
```

What happens:
- **Ghosts.** They are rows per cell, so block boundaries do not matter.
- **Register.** At initialisation the forest intersects each register face with the realm_b seam leaves, in cells.
  Each realm_b block face gets one **overlap row** per register face it touches (`maps%seam_overlap`); at every
  accumulation its skin is scattered into those faces (`flux_register%accumulate_fine_overlaps`). A lined-up block has
  one row covering its whole skin, which is exactly R3.

![Register cases](/forest/register.svg)

**Requirements.** The overlaps must be cell-aligned and cover each register face exactly; otherwise the forest stops
("the cells of realm A and realm B do not line up across the seam", or "the seam block ... faces N of its M seam
cells").

**Verified by:**
- FLUME `multirealm/check.sh` leg 5: Sod split along x and along z with misaligned blocks, bitwise, 2 and 4 ranks.
- MV-14 leg 5: the MHD rotor on an octree, misaligned both ways, bitwise on all 9 fields. This is the test with teeth:
  with the overlap offsets forced to zero the density goes negative at step 2.
- Regression `sod-2realm-blocks` (`equivalent_to sod-x`).
- PRISM `rmf-2realm-fv-pulse-blocks`: on the FV path the register closes exactly (mismatch 0). The fields agree to
  1e-13 only, because PRISM FV changes at round-off with the block layout alone.

## R6. Refined seam: a 2:1 jump between realms

![Refined seam](/forest/refined.svg)

One realm is exactly twice as fine as the other along every axis ([#52](https://github.com/szaghi/adam/issues/52)):

```ini
[forest.topology.face_1]
realm_a          = 1
face_a           = +x
realm_b          = 2
face_b           = -x
coupling         = refined           ; the coarse side is found from the cell sizes
coupling_cadence = stage_coincident
```

```ini
; realm 1 (coarse), x in [0, 0.75]          ; realm 2 (fine), x in [0.75, 1]
ni = 36   nj = 4   nk = 4                      ni = 12   nj = 4   nk = 4
iu_ref_levels = 2      ; dx = 1/192            iu_ref_levels = 3      ; dx = 1/384
seam_ghost_fill = tricubic ; [amr]             seam_ghost_fill = tricubic ; [amr]
```

What happens:
- **Rows.** The seam reuses the intra-realm 2:1 formulas. A fine-side ghost is an **interpolate** row: the owner of
  the coarse cells evaluates the `seam_ghost_fill` interpolant around the coarse cell containing the ghost. A
  coarse-side ghost is a **restrict** row: the mean of the 2x2x2 fine cells under it. Rows whose cells live on the
  same rank travel as messages to self.
- **Register.** One register face per coarse seam block. Its four fine blocks are each restricted 2:1 into their
  quadrant (`SEAM_KIND_INTER_REALM_REFINED`).
- **Time.** No subcycling: both realms take the global time step.

**Requirements.** Each is checked, with a named error:
1. A cell-size ratio of exactly 2 along every axis.
2. The same block cell counts along the seam (the tangential `n`).
3. The same `[amr] seam_ghost_fill` in both realms.
4. Even `ni`, `nj`, `nk` in the fine realm.
5. Nested blocks: every fine seam block covers one quadrant of one coarse seam block face.
6. Declaring such a pair `mirror` is refused ("coupling = mirror joins cells of the same size").

**Log.** `forest: face_pair 1 is a 2:1 seam, realm 1 coarse, realm 2 fine (issue #52)`.

**PRISM.** The 2:1 seam carries the same O(h^p) div(B) source as an intra-realm 2:1 face
([#29](https://github.com/szaghi/adam/issues/29)), and the `[IO] seam_divB_tol` monitor covers it.

**Verified by:**
- FLUME `multirealm/check.sh` leg 4: `sod-amr` split at its 2:1 face, along x and along z, bitwise, CPU np 2-4 and FNL.
- MV-14 leg 4: RJ2a on the same cells, all 9 fields bitwise.
- Regression `sod-amr-refined` (`equivalent_to sod-amr`).
- PRISM `rmf-2realm-fd-pulse-refined`: bitwise on all 27 fields, the #29 source unchanged, the monitor fires.

## R7. Both families: a mirror seam and AMR inside a realm

![Mixed forest](/forest/mixed.svg)

A realm of a forest may carry its own refinement, as in R2, away from its seam face. The register then holds both
kinds of face, the intra-realm 2:1 faces first, and one reduction completes them.

**Requirement.** Refinement must not reach the seam face: each realm's seam face keeps one cell size ("realm R has
seam cells of different sizes (refined blocks along the seam)"). To put a resolution jump **on** the seam, use R6.

**Verified by** FLUME `multirealm/check.sh` leg 2 and the regression case `sod-amr-2realm` (`equivalent_to sod-amr`).
Before the register composed the two kinds, the 2:1 face got no reflux.

## R8. Many ranks

![Seams across ranks](/forest/ranks.svg)

Nothing to configure. Each realm is partitioned over all ranks on its own (Morton order), so the two sides of a seam
may live on different ranks ([#40](https://github.com/szaghi/adam/issues/40)):
- **Rows.** Every rank enumerates every seam row from the replicated trees, in one canonical order. It keeps local rows
  (a direct copy, device to device on FNL), receive rows and send rows. The sends and receives match one to one with no
  index exchange (point-to-point MPI, tag 4040).
- **Register.** It is replicated and completed by a reduction per face.

**Check.** The `seam rows` log line: on an x split the rows are local, while on a z split at 2 ranks they are all
received and sent.

**Cost.** Set `ADAM_SEAM_TIMING=1` (`mpirun -x ADAM_SEAM_TIMING`) to print, at the end, the time per phase of the
seam fill and of the reflux ([#53](https://github.com/szaghi/adam/issues/53)). On FNL the cross-rank rows go through
host buffers, measured at 2.4-5.5% of the step on the test cases (WSL, not a benchmark).

**Verified by** FLUME `multirealm/check.sh` legs 1-5 at np 2, 3 and 4 (leg 3: the z split equals `sod-z` bitwise on
1 and N ranks).

A rank may own no block of a realm. Every collective in the step must then still be entered by that rank, so a
code change that skips work on empty ranks must keep the collectives. The PRISM FNL residuals once hung this way;
`rmf-2realm-fd-pulse-refined` at np 2 covers it.

## R9. Restart

Each realm restarts from its own files, with the usual keys in its own INI:

```ini
[IO]
restart_save     = 100                      ; write every 100 steps
restart          = .false.                  ; .true. to resume
restart_basename = sod-2realm-r1-restart    ; DIFFERENT in every realm
```

Give every realm its own `restart_basename` (and `output_basename`). Two realms writing to the same name overwrite each
other's files.

**Verified by** MV-14 legs 2 and 3: a 2-realm run stopped at half its steps and restarted equals the straight run
bitwise, with byte-identical histories; leg 3 does the same on an AMR box with a non-zero psi.

## R10. Splitting an existing input

`src/tests/flume/verification/mhd/multirealm/make_split.py` splits a single-realm input in two at the centre of x,
halving `ni` and writing the manifest:

```bash
make_split.py single.ini out-dir name                       # R3: lined-up blocks
make_split.py single.ini out-dir name --coarse-blocks 2     # R5: realm 2 on blocks twice as large
make_split.py single.ini out-dir name --set time.it_max=50  # override a key in both realms
```

For a refined split (R6), start from a single-realm run with a box marker. Give the coarse realm the coarse part,
unrefined; give the fine realm the refined part, built uniform one level finer (`iu_ref_levels + 1`, cell counts
chosen so that `dx` halves). Drop the markers in both. `verification/multirealm/sod-amr-refined*.ini` and
`src/tests/prism/regression/rmf-2realm-fd-pulse-refined/` are worked examples.

---

## Diagnostics

| Log line | Meaning |
|---|---|
| `forest: realm R seam rows: L local, R received, S sent` | rows this rank fills locally, receives and sends for realm R |
| `forest: face_pair N is a 2:1 seam, realm C coarse, realm F fine` | a `refined` pair was accepted |
| `forest: registered intra-realm AMR seam faces: +N` | 2:1 faces inside the realms (R2, R7) |
| `forest: reflux max\|F_coarse-F_fine_sum\| = ...` | printed only when non-zero: expected on 2:1 faces, never on a mirror seam |
| `forest timing (max over ranks, s): ...`, `forest seam timing realm R ...` | with `ADAM_SEAM_TIMING=1` |
| `WARNING: seam div(B) ... exceeds [IO].seam_divB_tol` | PRISM, with the monitor armed (R2, R6) |

## Refusals and fixes

| Message (abridged) | Cause | Fix |
|---|---|---|
| `only coupling = mirror or refined is implemented` | `periodic` / `interpolate` | reserved names; before #52 they silently ran as `mirror` |
| `coupling = mirror joins cells of the same size, the seam cells differ` | different cell sizes on a mirror pair | equal cells, or `coupling = refined` for a 2:1 jump |
| `coupling = refined needs ...` (a ratio of exactly 2, the same block cell counts along the seam, the same `seam_ghost_fill`, even counts, nested blocks) | an R6 requirement | see [R6](#r6-refined-seam-a-2-1-jump-between-realms) |
| `realm R has seam cells of different sizes` | refinement reaches the seam face | keep refinement away from the seam (R7), or put the jump on the seam (R6) |
| `the cells of realm A and realm B do not line up` / `faces N of its M seam cells` | misaligned cells, or faces that do not cover each other | align the extents to the cell size; make the seam faces coincide |
| `the seam ghost (i,j,k) of block b of realm R does not match the cells of realm P ...; ghost centre ..., peer cell ...` | a peer cell of the wrong size or position | as above |
| `stage_coincident requires equal scheme_time / rk_scheme / physics nv / K` | β on unlike realms | `end_of_step`, or align the integrators |
| `conflicting coupling_cadence between realm A and realm B` | two pairs between the same realms with different cadences | one cadence per realm pair |
| `[amr] ratio = 4 (quadtree) with AMR markers ... on the NVF and GMP backends` | #46: their ghost kernels are octree-only | `ratio = 8`, or the CPU/FNL backend |
| `[grid].(ni)=+7 is odd: 2:1 refinement ...` | #39 | even counts along refined axes |
| `(positivity_limiter)=cell is not supported on multi-realm runs` | the FLUME limiter | no limiter in a forest |

## Current limits

- **Topologies.** The schema takes any number of realms and face pairs, but the verified forests have two realms and
  one face pair.
- **Time stepping.** There is no subcycling: one global time step, set by the most restrictive realm.
- **FLUME limiter.** The positivity limiter is refused in a forest.
- **FNL cost.** On FNL the cross-rank seam rows go through the host on every fill
  ([#53](https://github.com/szaghi/adam/issues/53)).
- **Divergence at 2:1 seams.** A 2:1 seam, inter- or intra-realm, is not divergence-free for collocated B in PRISM
  ([#29](https://github.com/szaghi/adam/issues/29), accepted, monitored).
- **Reflux by app.** FLUME accumulates the seam fluxes of every stage weighted by its SSP coefficient. PRISM
  accumulates its final stage (α.r1), and only on the FV path: the FD path does not use the register.

## Verification map

| Scenario | FLUME | PRISM | Oracle |
|---|---|---|---|
| R2 intra-realm 2:1 | `sod-amr`, `amr-periodic-reflux` (V3), MV-11 | `rmf-amr*` | goldens; conservation to round-off |
| R3 mirror, lined up | leg 1, `sod-2realm`, `rj2a-2realm` | `rmf-2realm`, `rmf-2realm-fd-pulse` | bitwise; div-free |
| R4 α / β / asymmetric K | (all β) | `rmf-2realm`, `-asymK`, `-stagesync`, `-fd-pulse` | goldens; β matches the single realm |
| R5 misaligned blocks | leg 5, MV-14 leg 5, `sod-2realm-blocks` | `rmf-2realm-fv-pulse-blocks` | bitwise; register exact |
| R6 refined 2:1 | leg 4, MV-14 leg 4, `sod-amr-refined` | `rmf-2realm-fd-pulse-refined` | bitwise |
| R7 both families | leg 2, `sod-amr-2realm` | — | bitwise |
| R8 many ranks | legs 1-5 at np 2-4, leg 3 (z split) | np 2 in every case | bitwise |
| R9 restart | MV-14 legs 2-3 | — | bitwise, histories identical |

"leg N" is `src/tests/flume/verification/multirealm/check.sh`; "MV-14 leg N" is
`src/tests/flume/verification/mhd/multirealm/check.sh`. FLUME case names are regression cases under
`src/tests/flume/regression/`, PRISM ones under `src/tests/prism/regression/`.
