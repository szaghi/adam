# Shock-crossing defect is a capturing error

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| E0 | Seam studies | M6 ([#68](https://github.com/szaghi/adam/issues/68)) | FNL only | no ([#85](https://github.com/szaghi/adam/issues/85)) | — |

## Case and scope

[A3](../seams/shock-crossing) found an $O(1)$, non-converging defect behind a shock that crosses a 2:1 seam;
[D2](./limited-fills) showed that no coarse-to-fine fill reduces it and that limiting the fill makes it worse. E0 asks
whether the defect comes from the consistency of the fill with the restriction (only `restriction-compatible` with the
`mean` restriction has $R \circ P = I$ exactly), from the reflux, or from neither. The answer decides whether a seam
operator can fix it.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | Sod: compressible Euler, $\gamma = 1.4$; RJ2a: ideal MHD, $\gamma = 5/3$, no divergence control | [models](../../models#compressible-euler) |
| Space, time | WENO5 flux splitting on characteristic variables; SSP Runge–Kutta 3 stages, CFL 0.5 | as in [AS](../seams/shock-crossing) |
| Fill | `restriction-compatible` (order 2) or `tricubic` (order 4) | `seam_interpolate_compatible`, `seam_interpolate_tricubic`, [`adam_seam_interpolation_library.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_seam_interpolation_library.F90) |
| Restriction | `mean` or `lagrange6` | `restrict_seam_high_order`, [`adam_field_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_field_object.F90) |
| Reflux | on or off | `apply_reflux_to_stage_forest`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) and [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Cases | the [AS](../seams/shock-crossing) `c2f` crossings: Sod $N = 200$ (and 400), RJ2a $N = 256$; octree, $x > 0.75$ refined, $t = 0.2$, `extrapolation` on all faces |
| Arms | every combination of fill (`restriction-compatible`, `tricubic`) × restriction (`mean`, `lagrange6`) × reflux (on, off) |
| Ranks | 2, FNL |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/amr-shock/make_amr_shock.py` | the crossings (`--seam-restriction`, `--no-reflux`; the fill key is added to the input by hand: the script has no `--seam-ghost-fill` option) |
| `src/tests/flume/verification/amr-shock/amr_shock_oracle.py` | the plateau and band defects, the mass and the shock position |

No `check.sh` leg runs E0.

## Results

No figure yet ([#85](https://github.com/szaghi/adam/issues/85): the plateau density defect against fill × restriction
× reflux, and Sod at $N$ = 200 and 400).

Plateau density defect behind the shock:

| Fill + restriction | Sod, reflux on | Sod, reflux off | RJ2a, reflux on | RJ2a, reflux off |
|---|---|---|---|---|
| compatible + mean ($R \circ P = I$) | $2.32\cdot10^{-2}$ | $1.82\cdot10^{-2}$ | $1.97\cdot10^{-2}$ | $1.83\cdot10^{-2}$ |
| compatible + lagrange6 | $2.32\cdot10^{-2}$ | $1.85\cdot10^{-2}$ | $1.98\cdot10^{-2}$ | $1.84\cdot10^{-2}$ |
| tricubic + mean | $2.46\cdot10^{-2}$ | $3.10\cdot10^{-2}$ | $1.41\cdot10^{-2}$ | $1.68\cdot10^{-2}$ |
| tricubic + lagrange6 | $2.46\cdot10^{-2}$ | $3.17\cdot10^{-2}$ | $1.41\cdot10^{-2}$ | $1.68\cdot10^{-2}$ |

Sod with the reflux on at $N = 400$: $2.35\cdot10^{-2}$ (tricubic) and $2.14\cdot10^{-2}$ (compatible), while the band
defect next to the seam converges ($6.35\cdot10^{-3} \to 1.71\cdot10^{-3}$, order 1.9). With the reflux off the mass
drifts ($10^{-5}$ to $10^{-4}$) and the shock lands 0.06 to 0.22 cells off.

Related cases: [AS](../seams/shock-crossing) found the defect; [D2](./limited-fills) ruled out the fill's order;
[E2](./reflux-policy) proceeds with the reflux policy on its own merit.

## Findings

- **The restriction is irrelevant** (`mean` and `lagrange6` within 0.5 % in every case).
- **Consistency does not control it:** the $R \circ P = I$ pair is best on Sod and worse than `tricubic` on RJ2a.
- **The reflux does not cause it:** off, it moves the defect by up to 30 % either way and costs conservation.
- **It does not converge:** $N$ = 200 to 400 lowers it by 5 to 8 %.
- **So it is not a seam-operator defect.** The evidence points to the start-up error of a captured shock whose
  numerical structure, a few cells wide, is squeezed into finer cells (an interpretation, not measured directly): the
  transition emits a spurious wave (the plateau entropy defect is about 1.8 %) whose amplitude scales with the shock's
  strength, not with the spacing.
- **Remedy:** keep shocks off the seams by refining ahead of them (a gradient or Löhner marker with a buffer and a
  short regrid frequency, [M5](../runtime-amr/)), which [numerics](../../numerics#adaptive-mesh-refinement) now
  recommends. Phase E (the reflux policy) proceeds on its own merit.
