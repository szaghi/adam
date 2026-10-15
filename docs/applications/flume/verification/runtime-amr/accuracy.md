# Accuracy of runtime AMR

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| AV (av2 to av8) | Runtime AMR (M5) | M5 ([#74](https://github.com/szaghi/adam/issues/74)) | CPU, FNL | yes | `field-loop-regrid` |

## Case and scope

[RG](./runtime-regrids) proves the regrid machinery. AV asks what a user of runtime AMR wants: is a run whose grid
follows the solution about as accurate as the uniform run at its finest level, with fewer cells, and does a regrid
keep what the scheme guarantees? Each case compares a tracked run (a regrid every 5 steps) with uniform runs, against
bounds fixed in the #74 P4 plan before any run.

What the case catches: a regrid that costs accuracy (a feature that leaves its refined region, a marker that misses a
feature, a transfer error), a regrid that breaks conservation, positivity, symmetry or the divergence control, and a
regrid that makes the result depend on the number of ranks or on the backend.

Five bounds turned out to rest on false premises, all pre-existing properties of the scheme, of a marker or of the
backends that runtime AMR merely exposed; they were reported and replaced, each with its issue, never moved to fit a
number (see [Findings](#findings)).

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Euler cases (av2, av3a, av6, av7) | WENO5 flux splitting, characteristic variables; av2 `runge-kutta-ssp-54`, CFL 0.4; av3a and av7 `runge-kutta-ssp-33`, CFL 0.5 | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [numerics](../../numerics#scheme-space-weno-flux-splitting) |
| MHD blast (av3b) | ideal MHD, $\gamma = 1.4$, EGLM ($c_h = 60$, $\alpha = 0.18$), `weno` flux splitting, characteristic, positivity limiter (`cell`), floors disabled; `runge-kutta-ssp-54`, CFL 0.4 | `mhd_eglm_face_split_fluxes`, [`adam_flume_mhd_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_library.F90); `limit_positivity`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [EGLM](../../models#divergence-control-eglm), [positivity limiter](../../numerics#positivity-limiter) |
| MHD blast on `weno-riemann` (av3c) | HLLD, 6th-order correction with the WENO sensor, characteristic interpolation | `compute_residuals_riemann`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [`adam_flume_cpu_mhd_eglm_hlld_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_mhd_eglm_hlld_kernels.F90); [numerics](../../numerics#scheme-space-weno-riemann-interpolation-riemann-flux-and-correction) |
| Field loop (av4) | ideal MHD, $\gamma = 5/3$, GLM ($c_h = 4$, $\alpha = 0.18$) or EGLM, WENO5 characteristic, `runge-kutta-ssp-54`, CFL 0.4 | [GLM](../../models#divergence-control-mixed-glm), [EGLM](../../models#divergence-control-eglm) |
| Gradient marker | centred differences, across the faces between blocks too (#76); spacing `delta_fine` where $\max\lvert\nabla v\rvert >$ `tol`, else `delta_coarse`; acts at initialisation and at the regrids | `mark_by_gradient`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90) |
| Löhner marker | the FLASH/PLUTO dimension-sum estimator with an absolute `floor` (#74 P4); acts at the regrids | `mark_by_lohner`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90) |
| Solid marker (av7) | blocks crossed by the cylinder surface | `mark_by_solid`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90) |
| Regrid | every 5 steps, `conservative` prolongation, mean restriction, reflux on; no subcycling (D-M5-3) | `regrid_forest`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90), [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90); [numerics](../../numerics#adaptive-mesh-refinement) |

At a runtime regrid the markers vote per block (`mark_runtime`): a block is refined when any marker asks for it,
derefined only when every marker agrees, and a block at the base level is never derefined.

## Initial and boundary conditions

| Leg | Case | Tracked run | Uniform runs | Final time |
|---|---|---|---|---|
| av2 | [isentropic vortex](../../initial-conditions#isentropic-vortex) (V2), periodic quadtree on $[0, 1]^2$, $16^2$ cells per block | base level 1 ($N = 32$) to level 3 ($N = 128$); gradient marker on $\rho$, `tol` 0.05 | $N = 32$, $N = 128$ | $t = 0.2$ |
| av3a | [Sod](../euler/sod) (V1) on $[0, 1]$, quadtree with $y$ null, $12 \times 4$ cells per block; the jump at 0.5 on a block face at every level | base level 2 (48 cells) to level 4 (192); gradient on $\rho$ (`tol` 1.0) OR Löhner on $\rho$ (`refine_tol` 0.5, `derefine_tol` 0.2) | 192 cells | $t = 0.2$ |
| av3b | Balsara–Spicer MHD blast of the `blast-limiter` golden: $\rho = 1$, $p = 1000$ for $r < 0.1$ around $(0.5, 0.5)$ (sampled at the cell centres of $N = 64$), $p = 0.1$ outside, $\mathbf B = 100/\sqrt{4\pi}$ at 45° (ambient $\beta = 2.5\cdot10^{-4}$); periodic quadtree on $[0, 1]^2$, $16^2$ cells per block | base level 2 ($N = 64$) to level 3 ($N = 128$); gradient on $p$ (`tol` 100) OR Löhner on $p$ (0.5, 0.2) | $N = 128$ | $t = 0.01$ |
| av3c | the av3b blast on `weno-riemann` (HLLD) | as av3b | $N = 128$ | $t = 0.01$ |
| av4 | [field loop](../../initial-conditions#field-loop-advected-field-loop) (MV-9): $A_0 = 10^{-3}$, $R = 0.3$, $\rho = p = 1$, $\mathbf v = (2, 1, 1)$, on $[-1, 1] \times [-0.5, 0.5]$, periodic quadtree, $4 \times 4$ blocks of $N/4 \times N/8$ cells; GLM and EGLM | $N = 64$ to $N = 128$; gradient on $B_x$ (`tol` 0.01) OR Löhner on $B_x$ and $B_y$ (0.3, 0.1, `floor` $10^{-4}$) | $N = 128$ | $t = 1$ |
| av6 | the av2 tracked vortex | on 1, 2 and 3 ranks (1, 2 on FNL: the WSL box has two GPUs, and three ranks sharing them run about 170 times slower) | — | $t = 0.2$ |
| av7 | [shock over the cylinder](../euler/shock-cylinder) (V6), octree with $z$ null, $16 \times 16 \times 4$ cells per block | V6's initial AMR (base level 2, solid marker to level 3) and a regrid every 5 steps: solid OR Löhner on $\rho$ (0.6, 0.2) | uniform at V6's finest level (level 3); the init-AMR V6 run (no regrid) | $t = 0.25$ |
| av8 | the tracked runs of av2, av3a, av3b, av4 (GLM, EGLM), av7 on both backends | — | the uniform-fine pair of each | — |

Every run on 2 ranks unless stated (av6).

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/regrid/make_regrid.py` | `vortex-av`, `sod-av`, `blast-av`, `loop-av`, `cylinder-av` write the inputs |
| `src/tests/flume/verification/regrid/accuracy.sh` | runs the legs (`--leg av2` … `--leg av8`; `FLUME_EXE` selects the backend; av8 compares the work directories of both) |
| `src/tests/flume/verification/regrid/av_oracle.py` | the bounds, on top of the V1, V2, V6 and MV-9 oracles |
| `src/tests/flume/verification/conservation/conservation_oracle.py` | the av6 rank comparison |
| `src/tests/flume/regression/field-loop-regrid/input.ini` | the frozen av4 GLM tracked run, 150 steps, 10 regrids |
| `make_doc_figures.py` (`fig_regrid`) | the figure, from the CPU work directories `av2-tracked` and `av4-glm-tracked` |

```bash
cd src/tests/flume/verification/regrid && ./accuracy.sh                     # CPU, av2 to av7
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./accuracy.sh              # FNL
./accuracy.sh --leg av8                                                     # CPU against FNL
```

## Results

![Runtime regridding: tracked vortex and field loop with their blocks](/flume/regrid.png)

*The tracked vortex at $t = 0.2$ with its blocks and its pointwise density error, and the tracked GLM field loop at
$t = 1$ with its blocks: the refined blocks follow the vortex core and the loop edge.*

The bounds:

| Leg | Asserted |
|---|---|
| av2 | L1(ρ) ≤ 1.5× uniform-128 and ≤ 1/8 of uniform-32 |
| av3a | L1(ρ) ≤ 1.3× uniform-192; contact and shock on finest cells; null copies within 1e-12 |
| av3b | outer shock radius on the four half-axes within 2 finest cells of uniform-128; ρ conserved to 1e-13; ρE and B drift ≤ 1.1× uniform-128's; every regrid admissible |
| av3c | point symmetry within 1e-10 on the uniform and the tracked run ([#77](https://github.com/szaghi/adam/issues/77)) |
| av4 | E_B(T)/E_B(0) ≥ uniform-128's − 0.5 %; ⟨\|B_z\|⟩/A0 ≤ 1.2× (GLM; EGLM reported); final max\|div B\| ≤ 2×; div B at each regrid reported |
| av6 | the same regrids; fields within 1e-13 |
| av7 | bow-shock stand-off within 2 finest cells of uniform at V6's finest level; mirror symmetry 1e-10; positivity (the init-AMR V6 run's stand-off reported) |
| av8 | the same regrids; per field, the tracked CPU–FNL difference ≤ 2× the uniform-fine runs' (or ≤ 1e-10) |

The measured values (np 2; CPU numbers, FNL the same within the av8 agreement):

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

Related cases: the regrid machinery is [RG](./runtime-regrids); the point-symmetry loss behind av3b and av3c is
[the MHD blast on the flux-splitting path](./mhd-blast-symmetry); the uniform cases are the
[isentropic vortex](../euler/isentropic-vortex), [Sod](../euler/sod), the [shock over a cylinder](../euler/shock-cylinder),
the [field loop](../mhd/field-loop) and the [positivity limiter](../mhd/positivity-limiter) blast; the EGLM seam
level of #78 is followed up in [EGLM reflux](../seams-studies/eglm-reflux).

## Findings

- **Accuracy:** within the stated factor of the uniform-fine run in every case, and equal to it where the refined
  region covers the whole feature (av3a, av7).
- **Cost:** the savings are modest on these cases, and the reason is visible. The tracked vortex and loop hold half
  the blocks of the uniform runs, but the time step is set by the finest cells either way (no subcycling, D-M5-3) and
  the regrids, the 2:1 seams and the reflux cost their share, so the wall time drops by 3 to 19 % (av7 the most); the
  blast and the cylinder end fully refined, their waves filling the box by the final time.

What runtime AMR exposed (each reported before its bound changed):

- **A jump on a block face was invisible to the gradient marker** at initialisation (it read interior cells only):
  Sod's jump at 0.5 is a face at every level on $[0, 1]$, the initial grid stayed coarse, and the error of the first
  coarse steps persisted, 2.56× uniform-fine whatever the Löhner threshold. The marker now differences across the
  faces between blocks, reading one ghost layer (the initial AMR refreshes the host ghosts first), and av3a runs on
  $[0, 1]$ again: 0.999× ([#76](https://github.com/szaghi/adam/issues/76)).
- **The Löhner estimator needs an absolute floor** on a variable that vanishes in quiet regions: outside the field
  loop B is GLM residue (~10⁻⁶), the relative filter `epsilon` filters nothing, and the estimator read E ≈ 0.98
  everywhere, refining the whole box. The `floor` key (default 0, the FLASH/PLUTO form) fixes it
  ([numerics](../../numerics#adaptive-mesh-refinement)).
- **The MHD blast breaks point symmetry on a uniform grid** (0.2 relative at N = 128 from an exactly symmetric initial
  condition), and does not conserve ρE and B (it drifts 5e-5 on the uniform grid too): the planned symmetry bound was
  dropped and the conservation bound limited to ρ, with ρE and B bounded against the uniform drift
  ([#77](https://github.com/szaghi/adam/issues/77)): a grid-scale instability of the flux splitting at low β
  ([the MHD blast on the flux-splitting path](./mhd-blast-symmetry)); av3c asserts the symmetry on the Riemann-flux
  path.
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
