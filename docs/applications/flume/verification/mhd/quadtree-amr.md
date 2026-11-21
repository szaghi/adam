# Quadtree AMR

| ID | Group | Milestone | Type | Backends | Figure | Regression |
|---|---|---|---|---|---|---|
| MV-15 | Ideal MHD | [#46](https://github.com/szaghi/adam/issues/46) (follow-up) | C | CPU, FNL | yes | — |

## Case and scope

A quadtree (`ratio = 4`) refines $x$ and $y$ only, so its 2:1 seams are 2:1 in $x$ and $y$ and 1:1 in $z$. Until
[#46](https://github.com/szaghi/adam/issues/46) the seam machinery (ghost restriction and interpolation, reflux, the
limiter's seam flux synchronisation) assumed an octree: the coarse cells beside a quadtree seam picked up a spurious
$z$ dependence ($10^{-2}$ at the first step, negative pressure in Orszag–Tang at step 26). Quadtree AMR was therefore
refused, and every AMR verification ran on an octree with a null $z$ axis and `nk = 4`.

MV-15 runs the same 2-D problem on three trees, keys the cells by $(x, y)$ and compares them per variable:
- the octree with null $z$, `nk = 4`: the reference, the layout of every AMR verification before #46;
- a quadtree with `nk = 1`: true 2-D, 4 times fewer cells per block;
- a quadtree with an active $z$ axis, `nk = 4` (Orszag–Tang only): the $z$ fluxes are computed, and the solution must
  stay $z$-invariant.

What the case catches: a seam operator that mixes cells of different $z$ on a quadtree (the #46 defect), as a $z$
dependence of a $z$-invariant solution or as a difference from the octree reference above round-off.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | ideal MHD; GLM (Orszag–Tang), EGLM (blast) | [models](../../models#ideal-mhd) |
| Space | WENO5 flux splitting on characteristic variables; the blast with the positivity limiter | `compute_residuals_weno`, `limit_positivity`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Time | SSP Runge–Kutta, 5 stages, CFL 0.4 | [`adam_stepper_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_stepper_object.F90) (SELAH) |
| Seam ghosts | fine→coarse restriction (mean of 8 fine cells on an octree, 4 on a quadtree), coarse→fine tricubic fill | `interp_seam_ghost`, `update_ghost_local`, `update_ghost_mpi`, [`adam_field_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_field_object.F90); `seam_interpolate_tricubic`, [`adam_seam_interpolation_library.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_seam_interpolation_library.F90); [boundary conditions](../../boundary-conditions#_2-1-amr-coarse-fine-ghosts) |
| Reflux | fine face flux restricted per axis (2×1 fine faces per coarse face cell on a quadtree) | `restrict_fine_face_to_quadrant`, `face_tangential_ratios`, [`adam_flux_register_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_flux_register_object.F90) |
| Limiter seam synchronisation (blast) | per stage | `seam_sync_fluxes`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `seam_fine_to_coarse`, [`adam_flume_seam_sync_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_seam_sync_object.F90) |
| Boundaries | periodic in $x$, $y$; extrapolated $z$ boundaries keep a $z$-invariant start $z$-invariant | [boundary conditions](../../boundary-conditions#periodic) |
| FNL backend | the same kernels on the device | [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

## Initial and boundary conditions

| Leg | Problem | Grid | Trees | Final time |
|---|---|---|---|---|
| `ot` | [Orszag–Tang](./orszag-tang) ([`orszag-tang`](../../initial-conditions#orszag-tang)), GLM, $c_h = 4$ | $[0, 1]^2$, $32^2$ base, the blocks of $[0.25, 0.75]^2$ refined once (2:1) | octree `nk = 4`; quadtree `nk = 1`; quadtree `nk = 4` | 0.2 |
| `blast` | Balsara–Spicer blast ([positivity limiter](./positivity-limiter) PV-1), EGLM, limiter on | the same | octree `nk = 4`; quadtree `nk = 1` | 0.006 (the outer shock crosses the seams) |

Both legs: 2 ranks.

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/mhd/orszag-tang/make_orszag_tang.py` | writes the `ot` inputs (`--cells 32 --time-max 0.2 --refine-box 0.25 0.25 0.75 0.75 --ratio 8\|4 [--nk 4]`) |
| `src/tests/flume/verification/mhd/positivity-probe/make_blast.py` | writes the `blast` inputs (`--eglm --cells 32 --time-max 0.006 --refine-box ... --ratio 8\|4`) |
| `src/tests/flume/verification/mhd/amr_box.py` | the one-level box refinement on either tree |
| `src/tests/flume/verification/mhd/quadtree/check.sh` | runs the legs (`--legs ot,blast`, `--np N`, `--keep`) |
| `src/tests/flume/verification/make_doc_runs.sh <runs> quadtree`, `make_doc_figures.py --only quadtree` | re-runs the `ot` leg on the octree and the quadtree with checkpoints kept (2 runs, about 3 min) |
| `src/tests/flume/verification/mhd/quadtree/quadtree_oracle.py` | the per-variable $z$ invariance and the agreement with the reference at the same $(x, y)$ |
| `src/tests/amr/test_quadtree_seam_ghost.F90` | the unit test of the seam exchange on a linear field |

```bash
cd src/tests/flume/verification/mhd
quadtree/check.sh                                                     # CPU (about 6 min)
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl quadtree/check.sh    # FNL
```

## Results

![Orszag-Tang on the quadtree against the octree](/flume/quadtree.png)

*Orszag–Tang, $32^2$ with $[0.25, 0.75]^2$ refined 2:1, $t = 0.2$, CPU: the density on the quadtree ($nk = 1$) with its blocks; the per-cell density difference against the octree with a null $z$ axis (its first $z$ layer), largest along the refined box and with no structure; and the largest difference per variable relative to the variable's largest value ($\psi$: $8\cdot10^{-12}$, the others $\le 6\cdot10^{-13}$): round-off.*

| Leg | Quadtree `nk = 1` against the octree | Octree z spread | Quadtree `nk = 4` z spread |
|---|---|---|---|
| Orszag–Tang, $32^2$ + $[0.25, 0.75]^2$ refined, $t = 0.2$ | $4.1\cdot10^{-12}$ ($\psi$; $\le 3\cdot10^{-13}$ elsewhere) | $5.7\cdot10^{-12}$ | 0 |
| Balsara–Spicer blast, limiter on, $t = 0.006$ | $3.5\cdot10^{-4}$ ($\psi$), $1.0\cdot10^{-4}$ ($\rho$) | $1.9\cdot10^{-3}$, $3.0\cdot10^{-4}$ | — |

| Check | Bound | CPU | FNL |
|---|---|---|---|
| `ot`: quadtree `nk = 1` against the octree, relative | $10^{-10}$ | $4.1\cdot10^{-12}$ | $6.4\cdot10^{-12}$ |
| `ot`: octree $z$ spread | $10^{-10}$ | $5.7\cdot10^{-12}$ | pass |
| `ot`: quadtree `nk = 4` $z$ spread | 0 (bitwise) | 0 | 0 |
| `blast`: octree $z$ spread | $10^{-4}$ | $1.9\cdot10^{-3}$ ($\psi$): fails since #134 D4 ([findings](#findings)) | $2.5\cdot10^{-3}$: the same |
| `blast`: quadtree against the octree, per variable | $\max(10^{-10}, 2\times$ the octree $z$ spread$)$ | ratio 0.09 to 0.33 | ratio 0.56 to 1.06 |
| `blast`: limited stages | reported | 578 | 578 |

Related cases: [Orszag–Tang](./orszag-tang) is the `--amr` problem of MV-12; [positivity limiter](./positivity-limiter)
PV-5 runs the limited blast across octree seams; [forward-facing step](../euler/forward-facing-step) runs a forest on
quadtrees; [ghost probe](../infrastructure/ghost-probe) checks every seam ghost on a linear field.

## Findings

- The runs agree to round-off, not bit for bit: an octree restriction averages 8 fine cells (two identical $z$ layers),
  a quadtree 4, and even the octree's identical layers differ at round-off, its tricubic weights depending on the $z$
  sub-position.
- In the blast the octree's own layers drift apart: the tricubic seam ghosts start them apart at round-off, the
  flux-splitting scheme's low-β grid-scale instability (the one of [the MHD blast symmetry
  study](../runtime-amr/mhd-blast-symmetry)) grows the difference by an e-fold about every 10 steps, and when the outer
  shock reaches the 2:1 seam the limiter's switches amplify it erratically, by an order of magnitude from one step to
  the next. At $t = 0.006$ the spread is $1.9\cdot10^{-3}$ ($\psi$) since [#134](https://github.com/szaghi/adam/issues/134)
  D4, $2.3\cdot10^{-5}$ before it, and $1.8\cdot10^{-4}$ on SSP-33, which D4 leaves bitwise: D4 changes the sample at
  the seam crossing, not the growth before it (a factor $\le 2$ at step 88). On `weno-riemann` HLLD the spread stays
  at $10^{-10}$. The quadtree differs from the octree by less than that spread (ratio 0.09 to 0.33 on the CPU, 0.56 to 1.06 on
  FNL, whose octree spread is $2.5\cdot10^{-3}$; the same 578 limited stages): the leg bounds the agreement by twice the reference spread, not by a fixed tolerance. The
  fixed bound on the octree spread itself ($10^{-4}$, 3.5 times the single pre-D4 measurement) lies inside the
  spread's variability and fails since D4; its value is under review.
- FNL gives the same picture on Orszag–Tang ($6.4\cdot10^{-12}$ against the octree, $z$ spread 0 with `nk = 4`).
- The `nk = 4` quadtree computes the $z$ fluxes and keeps the solution $z$-invariant bit for bit. Its $z$ cells enter
  the time step, so its step sequence differs and it is compared with itself only.
- The unit test `test_quadtree_seam_ghost` pins the seam exchange on the linear field $f = 1 + 0.3x + 0.7y + 1.1z$
  (octree, quadtree `nk = 4` and `nk = 1`, 1 to 3 ranks): every face, edge and corner ghost inside the domain must hold
  $f$ at its own centre. Its FNL build also asserts the device exchange equal to the CPU one, bit for bit.
