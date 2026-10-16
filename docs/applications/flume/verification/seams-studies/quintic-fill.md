# Quintic coarse-to-fine fill

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| D1 | Seam studies | M6 ([#68](https://github.com/szaghi/adam/issues/68)) | CPU, FNL (table: FNL) | yes | — |

## Case and scope

With the 6-point restriction ([C1](../seams/order-across-seams)) and the reflux off, the smooth AO vortex reaches the
uniform coarse error, but its order (5.33) is that of the tricubic coarse-to-fine fill's $O(h^4)$ term waiting to show.
D1 asks whether a 6-point (quintic) fill brings the composite to the scheme's order, and what it does to a shock
crossing a 2:1 face (the Phase D exit: AO L1 order $\ge 4.7$ and L∞ $\ge 4.5$, and A3 no worse than `tricubic`).

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Method | `[amr] seam_ghost_fill = quintic`: the tensor 6-point Lagrange interpolation at the fine ghost's centre, 1-D weights $k/8192$ (exact in binary); the footprint shifted inside the coarse donor block at its edge as the tricubic one is (blocks of at least 6 cells along the refined axes); along an unrefined axis only the anchor is read | `seam_interpolate_quintic`, `SEAM_W_QUINTIC`, `seam_shift_anchor_pos`, [`adam_seam_interpolation_library.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_seam_interpolation_library.F90); `interp_seam_ghost`, [`adam_field_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_field_object.F90) |
| Anchor shifts | carried in the fill meta of every map row, bits 18–26 | `seam_meta_pack`, [`adam_seam_interpolation_library.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_seam_interpolation_library.F90); `make_local_map_ghost_cell`, [`adam_maps_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_maps_object.F90) |
| FNL | summed straight from the device field, no per-thread footprint | `update_ghost_local_gpu_dev`, [`adam_fnl_field_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/fnl/adam_fnl_field_kernels.F90) |
| Inter-realm | the same fill on the inter-realm rows | [`adam_forest_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_forest_object.F90) |
| Restriction | `lagrange6` in every run below | [C1](../seams/order-across-seams) |
| Key | read in `initialize`; the default stays `tricubic` | [`adam_adam_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_adam_object.F90) |
| Scheme, time, reflux | as in [AO](../seams/order-across-seams) and [AS](../seams/shock-crossing): WENO5 flux splitting; reflux off on AO, on for the crossings | [numerics](../../numerics#adaptive-mesh-refinement) |

## Initial and boundary conditions

| Run | Setup |
|---|---|
| AO `euler-q6` | the [AO](../seams/order-across-seams) Euler vortex, quadtree, $N$ = 64 / 128 / 256, `lagrange6` + `quintic`, reflux off |
| Shock crossings | the [AS](../seams/shock-crossing) Sod ($N = 200$) and RJ2a ($N = 256$) `c2f` and `f2c` crossings, `lagrange6` restriction, reflux on, the null axes widened to 8 cells |
| PV-5 | the limited Balsara–Spicer blast across 2:1 seams of the [positivity limiter](../mhd/positivity-limiter), split and HLLD |
| GP | the [ghost probe](../infrastructure/ghost-probe) on a degree-5 field |
| Multirealm leg 6 | a refined inter-realm seam, bitwise against the single-realm run ([multi-realm runs](../euler/multi-realm)) |

All runs on 2 ranks unless stated.

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/amr-order/check.sh --leg euler-q6` | the AO ladder: orders $\ge 5.5$ and the errors at most the baseline plus 2 % |
| `src/tests/flume/verification/amr-shock/make_amr_shock.py` | the crossings (`--seam-restriction`, `--seam-ghost-fill`) |
| `src/tests/flume/verification/ghosts/check.sh --seam-restriction lagrange6 --seam-ghost-fill quintic --degree 5` | GP on a degree-5 field |
| `src/tests/flume/verification/make_doc_figures.py --only quintic-fill` | the figure, from the probe runs of `ghosts/check.sh --degree 5 --keep` (default operators, then lagrange6 + quintic) and the kept AO work directories (CPU) |
| `src/tests/flume/verification/multirealm/check.sh` | leg 6, every seam operator across a refined inter-realm seam |
| `src/tests/flume/verification/mhd/positivity/check.sh` | PV-5 |

```bash
cd src/tests/flume/verification/amr-order && ./check.sh --leg euler-q6
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh --leg euler-q6
```

## Results

![D1: the quintic fill in situ and on the AO vortex](/flume/quintic-fill.png)

*CPU. Left and centre: the ghost probe on a degree-5 field, box3d (octree) and channel2d (quadtree), 2 ranks, the largest error per ghost group: mean + tricubic misses the field by up to $5\cdot10^{-5}$ in every group a seam operator fills (the expected failure: exact to degree 1 and 3), lagrange6 + quintic matches it to $10^{-14}$. Right: the AO ladder without reflux, lagrange6 with the tricubic and the quintic fill.*

FNL, 2 ranks (plateau values in units of $10^{-2}$; band values absolute):

| Check | lagrange6 + tricubic | lagrange6 + quintic |
|---|---|---|
| AO euler, reflux off, N = 256: L1 all / band, Linf | $8.55\cdot10^{-8}$ / $2.61\cdot10^{-7}$, $6.36\cdot10^{-6}$ | $6.03\cdot10^{-8}$ / $1.74\cdot10^{-7}$, $4.78\cdot10^{-6}$ |
| AO orders L1 / Linf (128 → 256) | 5.33 / 5.53 | 6.11 / 6.35 |
| Sod c2f plateau $\rho$ / $p$, band $\rho$ | $2.457$ / $0.691$, $0.633\cdot10^{-2}$ | $2.631$ / $0.858$, $0.795\cdot10^{-2}$ |
| Sod f2c plateau $\rho$ / $p$, band $\rho$ | $1.575$ / $1.015$, $0.343\cdot10^{-2}$ | $1.573$ / $0.985$, $0.338\cdot10^{-2}$ |
| RJ2a c2f plateau $\rho$ / $p$, band $\rho$ | $1.410$ / $0.622$, $1.03\cdot10^{-4}$ | $2.297$ / $0.986$, $0.68\cdot10^{-4}$ |
| PV-5 split / HLLD: most faces limited, inadmissible stages | 46336 / 17920, 0 | 45644 / 17920, 0 |

Backend agreement: np 1 / 2 / 3 bitwise on the CPU, np 1 / 2 on FNL, CPU against FNL $2.7\cdot10^{-13}$; the AO
`euler-q6` ladder agrees on CPU and FNL to every printed digit.

Related cases: [D2](./limited-fills) limits this fill; [E2](./reflux-policy) keeps its order with the reflux on; the
baselines are [AO](../seams/order-across-seams) and [AS](../seams/shock-crossing).

## Findings

- **The smooth order reaches the scheme's.** With the reflux off, the 6-point restriction and the quintic fill, AO
  euler is sixth order (6.11 / 6.35) and its error at $N = 256$ is below the uniform coarse run's
  ($8.38\cdot10^{-8}$): the Phase D exit (L1 $\ge 4.7$, L∞ $\ge 4.5$) holds, and the 2:1 seams no longer cost accuracy
  on this flow.
- **Unlimited, it worsens a shock entering a finer region.** On the `c2f` crossings the plateau defect behind the
  shock grows (Sod +7 % in $\rho$, +24 % in $p$; RJ2a +63 % / +59 %): the shock's coarse values reach the fine ghosts
  through a 6-point stencil with larger overshoots than the 4-point one. `f2c` is unchanged, shock positions are the
  same, mass is exact, and PV-5 passes with no inadmissible stage.
- The Phase D exit requires A3 no worse than `tricubic`, which is the job of [D2](./limited-fills)'s limited variants;
  until then `quintic` is an opt-in for smooth flows.
- Coverage: AO leg `euler-q6` asserts the order; GP checks the fill on a degree-5 field; multirealm leg 6 runs it
  across a refined inter-realm seam, bitwise.
