# Inter-realm 1:1 seams: α and β cadence

| Cases | Issues | Type | Backends | Figure | Regression |
|---|---|---|---|---|---|
| `rmf-2realm`, `rmf-2realm-asymK`, `rmf-2realm-stagesync`, `rmf-2realm-fd-pulse` | [#13](https://github.com/szaghi/adam/issues/13) (Phase D), [#16](https://github.com/szaghi/adam/issues/16), [#18](https://github.com/szaghi/adam/issues/18), [#31](https://github.com/szaghi/adam/issues/31) | C | CPU, FNL | no | `rmf-2realm`, `rmf-2realm-asymK`, `rmf-2realm-stagesync` (digest + residuals, CPU and FNL) |

## Case and scope

A forest manifest splits the domain at $x = 0$ into two realms of the same resolution, glued by a `mirror` seam: the
seam ghosts of each realm are a copy of the peer's interior. The union of the two realms' cells is the cell set of the
single-realm case. The seam is filled at one of two cadences:

- **α**, `coupling_cadence = end_of_step`: once per step, after every realm has closed its step. The RK stages read the
  ghosts of the previous step (the Berger–Oliger convention of AMR subcycling);
- **β**, `coupling_cadence = stage_coincident`: at every RK stage. Both realms must agree on the time scheme, the RK
  scheme, `nv` and the stage count (`check_beta_admissibility`).

Four cases:
- `rmf-2realm` (α) and `rmf-2realm-asymK` (α, SSP (3,3) in realm 1 against SSP (5,4) in realm 2, the same global
  $\Delta t$): digest goldens of their own. Under α the seam is first order in time, so the split is not expected to
  equal the single realm;
- `rmf-2realm-stagesync` (β): its digest must also match the single-realm `rmf` golden
  ([#18](https://github.com/szaghi/adam/issues/18));
- `rmf-2realm-fd-pulse` (β, source-free): the divergence oracle of [#31](https://github.com/szaghi/adam/issues/31). A
  1:1 mirror seam is numerically a 1:1 block face, so both divergences must stay at round-off.

What the cases catch:
- a seam fill that differs from the single-realm block exchange (β digest against the `rmf` golden);
- RK stages reading unfilled seam ghosts: under α the stages after the first read the stage buffer without seam
  ghosts, and the curl evolves the seam skin from garbage (#31);
- a divergence diagnostic that reads stale seam ghosts: `post_step_forest` refills the seam from the peers before the
  divergence (#31);
- a K-equality guard that refuses realms with different stage counts under α
  ([#16](https://github.com/szaghi/adam/issues/16)).

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | Maxwell, SI; 4 AC coils (`rmf-2realm*`) or no source (`rmf-2realm-fd-pulse`) | [models](../models#models) |
| Space | centred finite differences of order 6 | [`adam_fdv_operators_library.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_fdv_operators_library.F90) |
| Time | SSP Runge–Kutta (5,4), CFL 0.5 (realm 1 of `asymK`: SSP (3,3)); one global $\Delta t$ (minimum over the realms) | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Seam fill | peer-interior copy, per step (α) or per stage (β) | `fill_seam_from_peer_forest`, [`adam_prism_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/prism/cpu/adam_prism_cpu_object.F90), [`adam_prism_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/prism/fnl/adam_prism_fnl_object.F90); `evolve_one_step`, `check_beta_admissibility`, [`adam_forest_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_forest_object.F90) |
| Seam face BC | the manifest overrides each realm's BC on the seam face | [numerics](../numerics#multi-realm-runs) |

## Initial and boundary conditions

| Item | `rmf-2realm`, `-asymK`, `-stagesync` | `rmf-2realm-fd-pulse` |
|---|---|---|
| Realms | $x \in [-0.06, 0]$ and $[0, 0.06]$, $y, z \in [-0.06, 0.06]$ m | the same |
| Grid per realm | 64 blocks of $8\times16\times16$ (`iu_ref_levels = 2`), $\Delta x = 1.875\cdot10^{-3}$: the cells of `rmf` | 8 blocks of $8\times16\times16$ (`iu_ref_levels = 1`), $\Delta x = 3.75\cdot10^{-3}$: the cells of `rmf-amr-fd-pulse` without its marker |
| Initial condition | `vacuum`, 4 AC coils | [Gaussian pulse](../initial-conditions#gaussian-pulse) along $+y$, `sigma = 0.01`, `B0 = 10` |
| Boundaries | `Neumann` on the outer faces | the same |
| Cadence | α, α, β | β; α in the negative-control leg |
| Steps, ranks | 10, 2 (`run.sh`) | 5, 1 (the seam is a same-rank peer copy) |

## Input files

| File | Role |
|---|---|
| `src/tests/prism/regression/rmf-2realm*/input.ini` | the forest manifests |
| `src/tests/prism/regression/rmf-2realm*/realm_{1,2}.ini` | the realm inputs |
| `src/tests/prism/regression/rmf-2realm/README.md` | the α semantics and the seam-localisation measurement |
| `src/tests/prism/regression/rmf-2realm-fd-pulse/check.sh` | the divergence oracle (`--build`; `PRISM_EXE`) |
| `src/tests/prism/regression/run.sh` | the digests; for `rmf-2realm-stagesync` also the comparison with `rmf/golden/<backend>/digest.txt` |

```bash
REGRESSION_CASES="rmf-2realm rmf-2realm-asymK rmf-2realm-stagesync" src/tests/prism/regression/run.sh cpu --no-build
cd src/tests/prism/regression/rmf-2realm-fd-pulse && ./check.sh
```

## Results

| Case, leg | Asserts | Bound | CPU | FNL |
|---|---|---|---|---|
| `rmf-2realm`, `-asymK` | digest and residuals against their α goldens | `rtol = 1e-6`, `atol = 1e-3` | pass | pass |
| `rmf-2realm-stagesync` | digest against its β golden and against the single-realm `rmf` golden | `rtol = 1e-6`, `atol = 1e-3` | pass | pass |
| fd-pulse, β | max $\lvert\nabla\cdot\mathbf D\rvert$, max $\lvert\nabla\cdot\mathbf B\rvert$, both realms | $\le 10^{-13}$ | 0, 0 (both realms) | 0, 0 (both realms) |
| fd-pulse, health | both realms reach 100 %, both divergence histories exist | — | pass | pass |
| fd-pulse, α (negative control) | max $\lvert\nabla\cdot\mathbf B\rvert$ above the bound on at least one realm | $> 10^{-13}$ | $4.56\cdot10^{2}$ (both realms) | $4.56\cdot10^{2}$ (both realms) |

Seam localisation under α (`rmf-2realm`, measured before the α golden was ratified), max $\lvert B_z\rvert$ by
distance from the seam:

| Distance from $x = 0$ | max $\lvert B_z\rvert$ | $\lvert B_z\rvert/\lvert B_y\rvert$ |
|---|---|---|
| 0–0.015 (seam-adjacent block) | $1.0\cdot10^{-7}$ | 6.4 % |
| 0.015–0.030 | $5.4\cdot10^{-10}$ | 0.04 % |
| 0.030–0.045 | $4.3\cdot10^{-10}$ | 0.03 % |
| 0.045–0.060 | $4.3\cdot10^{-10}$ | 0.03 % |

Related pages: [misaligned blocks](./inter-realm-misaligned) and [inter-realm 2:1 seam](./inter-realm-refined) extend
the seam to other block layouts and to a resolution jump; FLUME's [forest seams](/applications/flume/verification/euler/multi-realm)
carry the same contract.

## Findings

- Under β the 1:1 seam is divergence-free to the last bit, as a 1:1 block face is: the defect of #31 was machinery, not
  a fundamental obstruction.
- Under α the same seam leaks ($4.56\cdot10^{2}$ in 5 steps). α's end-of-step lag suits subcycled AMR, not a 1:1 seam
  with one $\Delta t$; such a seam needs β. The α leg asserts only that α leaks: it gives the β leg teeth, because a
  seam that stopped being exercised would also read zero.
- The α coupling error stays local: it drops by more than 200 times within one block (8 cells) of the seam.
