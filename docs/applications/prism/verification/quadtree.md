# Quadtree 2:1 seams

| Cases | Issues | Type | Backends | Figure | Regression |
|---|---|---|---|---|---|
| `rmf-amr-fd-pulse`, quadtree leg | [#46](https://github.com/szaghi/adam/issues/46) | C | CPU, FNL | yes | — |

## Case and scope

A quadtree (`[amr] ratio = 4`) refines $x$ and $y$ only: its blocks span the domain along $z$ at every level, so a 2:1
seam is 2:1 in $x$ and $y$ and 1:1 in $z$. The seam machinery (ghost restriction, coarse–fine interpolation, reflux)
handles both trees ([#46](https://github.com/szaghi/adam/issues/46)). Before #46 the coarse cells beside a quadtree
seam picked up a spurious $z$ dependence, and quadtree AMR was refused at initialization.

The leg runs the source-free pulse of [seam div(B)](./seam-divergence) on the octree of the regression case and on a
quadtree, and compares them column by column. The pulse is uniform in $x$ and $z$, so the seam is the only thing that
can perturb it: it generates $B_y$ and the seam $\nabla\cdot\mathbf B$ source of
[#29](https://github.com/szaghi/adam/issues/29) ([hazards](../hazards#seam-div-b)).

What the case catches: a seam operator that mixes cells of different $z$ on a quadtree, as a $z$ dependence of a
$z$-invariant solution or as a difference from the octree above round-off.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | Maxwell, SI, no source | [models](../models#models) |
| Space | centred finite differences of order 6 | [`adam_fdv_operators_library.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_fdv_operators_library.F90) |
| Time | SSP Runge–Kutta (5,4), CFL 0.5; $\Delta t = \mathrm{CFL}\,\min(\Delta x, \Delta y, \Delta z)/(2c)$, set on both trees by the finest $x$ and $y$ cells, so the two runs take the same steps | [numerics](../numerics#time-integration) |
| Seam ghosts | fine→coarse restriction (mean of 8 fine cells on an octree, 4 on a quadtree), coarse→fine tricubic fill | `seam_interpolate_tricubic`, [`adam_seam_interpolation_library.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_seam_interpolation_library.F90); [`adam_field_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_field_object.F90) |
| Reflux | fine face flux restricted per axis | `restrict_fine_face_to_quadrant`, [`adam_flux_register_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_flux_register_object.F90) |

## Initial and boundary conditions

| Item | Octree (reference) | Quadtree |
|---|---|---|
| Domain | $[-0.06, 0.06]^3$ m | the same |
| Tree | `ratio = 8` | `ratio = 4` |
| Grid | 8 blocks of $16^3$ at level 1, $x < 0$ refined to level 2; finest $\Delta z = 1.875\cdot10^{-3}$ | one $z$ layer of blocks, `nk = 16`: $\Delta z = 7.5\cdot10^{-3}$ |
| 2:1 seam | $x = 0$, 2:1 in $x$, $y$, $z$ | $x = 0$, 2:1 in $x$, $y$, 1:1 in $z$ |
| Initial condition | [Gaussian pulse](../initial-conditions#gaussian-pulse) along $+y$, `sigma = 0.01`, `B0 = 10` | the same |
| Boundaries | `Neumann` | the same |
| Steps, ranks | 5, 1 | 5, 1 |

## Input files

| File | Role |
|---|---|
| `src/tests/prism/regression/rmf-amr-fd-pulse/input.ini` | the octree case; the leg rewrites `ratio` for the quadtree |
| `src/tests/prism/regression/rmf-amr-fd-pulse/check.sh` | the quadtree leg is default-on; `--keep` keeps its checkpoints |
| `src/tests/flume/verification/mhd/quadtree/quadtree_oracle.py` | shared with FLUME's MV-15: keys the cells by their $(x, y)$ centre |
| `src/tests/prism/regression/make_doc_figures.py` | draws the figure from the `--keep` work directories |

```bash
cd src/tests/prism/regression/rmf-amr-fd-pulse && ./check.sh --keep   # CPU
PRISM_EXE=$PWD/../../../../../exe/adam_prism_fnl ./check.sh           # FNL
```

The oracle asserts, per field:
- the $z$ spread of every $(x, y)$ column of each run at most $10^{-13}$;
- the quadtree against the octree at the same $(x, y)$ at most $10^{-13}$;
- the quadtree's max $\lvert\nabla\cdot\mathbf D\rvert$ at most $10^{-13}$.

Each field is measured against the largest magnitude of its vector ($\mathbf D$, $\mathbf B$, $\mathbf J$). The pulse
carries $D_z$ and $B_x$, and the seam adds $B_y$. The other components are round-off on the octree ($D_x$, $D_y$
$\sim 3\cdot10^{-18}$, $B_z \sim 5\cdot10^{-16}$ against $\lvert B_x\rvert = 10$) and exactly zero on the quadtree:
scaled by their own magnitude they would compare round-off with round-off.

## Results

![PRISM pulse across a 2:1 seam, octree against quadtree](/prism/quadtree-pulse.png)

*Top row: $B_x$ on $z = 0$ (quadtree, blocks outlined: the refined half $x < 0$ and the seam at $x = 0$), and $B_x$ on
$y = 0$ for both trees: the octree refines its $x < 0$ blocks in $z$ as well, the quadtree keeps one $z$ layer of
blocks. Bottom row: the quadtree against the octree per $(x, y)$ column (relative to each vector's magnitude); $B_y$,
the field the seam generates; and the divergence histories of the two runs, superposed.*

| Check | Bound | Octree | Quadtree |
|---|---|---|---|
| $z$ spread of a column, max over $\mathbf D$ and $\mathbf B$ | $10^{-13}$ | $2.6\cdot10^{-16}$ (FNL $3.6\cdot10^{-16}$) | 0 |
| against the octree, per $(x, y)$ column, max over $\mathbf D$ and $\mathbf B$ | $10^{-13}$ | — | $1.8\cdot10^{-16}$ (FNL $3.6\cdot10^{-16}$) |
| seam max $\lvert\nabla\cdot\mathbf B\rvert$ at step 5 | reported | 12.65004 | 12.65004 |
| max $\lvert\nabla\cdot\mathbf D\rvert$ | $10^{-13}$ | $3.1\cdot10^{-15}$ | 0 |

CPU values; FNL in parentheses where it differs.

Related pages: FLUME's [quadtree AMR](/applications/flume/verification/mhd/quadtree-amr) (MV-15) runs the same
oracle on ideal MHD.

## Findings

- The quadtree reproduces the octree to round-off, including the seam's $\nabla\cdot\mathbf B$ source to every printed
  digit. On a $z$-invariant field the per-axis seam formulas (restriction over $2\times2$ fine cells, coincident
  interpolation rows along $z$) give the same values as the octree's $2\times2\times2$ ones on two identical $z$ layers.
- The octree's own $z$ spread and its non-zero $\nabla\cdot\mathbf D$ are round-off of its $z$ interpolation: its
  tricubic weights depend on the $z$ sub-position. The quadtree has neither.
