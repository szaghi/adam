# Inter-realm 2:1 seam

| Cases | Issues | Type | Backends | Figure | Regression |
|---|---|---|---|---|---|
| `rmf-2realm-fd-pulse-refined` | [#52](https://github.com/szaghi/adam/issues/52) | C | CPU, FNL | no | — |

## Case and scope

`rmf-amr-fd-pulse` ([seam div(B)](./seam-divergence)) split at its 2:1 face $x = 0$ into two realms: a fine realm
($x < 0$, built uniform one level finer, no marker) and a coarse realm ($x > 0$), glued by `coupling = refined`. The
union of the realms is the single-realm cell set. The inter-realm 2:1 seam reuses the intra-realm 2:1 formulas
(coarse→fine ghost interpolation, the fine→coarse $2\times2\times2$ mean, the 2:1 flux register), so the split must
reproduce the single realm **bit for bit**, its [#29](https://github.com/szaghi/adam/issues/29) seam
$\nabla\cdot\mathbf B$ source included: the same source, not a smaller or a larger one.

What the case catches:
- any difference between the inter-realm and the intra-realm 2:1 seam operators;
- a seam $\nabla\cdot\mathbf B$ monitor blind to inter-realm seams;
- a collective ghost update that hangs or misbehaves on an empty rank: at 2 ranks one rank holds no block of one
  realm.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | Maxwell, SI, no source | [models](../models#models) |
| Space | centred finite differences of order 6 | [`adam_fdv_operators_library.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_fdv_operators_library.F90) |
| Time | SSP Runge–Kutta (5,4), CFL 0.5, β cadence | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Seam | `coupling = refined`: the intra-realm 2:1 interpolation, restriction and register across the realm boundary | `fill_seam_from_peer_forest`, [`adam_prism_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/prism/cpu/adam_prism_cpu_object.F90); [`adam_seam_exchange.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_seam_exchange.F90); [`adam_seam_interpolation_library.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_seam_interpolation_library.F90) |
| Monitor | `[IO] seam_divB_tol`, warn-only | `save_divergence_history`, [`adam_prism_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/prism/common/adam_prism_common_object.F90); [hazards](../hazards#seam-div-b) |

## Initial and boundary conditions

| Item | Fine realm | Coarse realm | Single realm (reference) |
|---|---|---|---|
| Extent | $x \in [-0.06, 0]$ | $x \in [0, 0.06]$ | $[-0.06, 0.06]^3$ |
| Grid | 64 blocks of $8\times16\times16$ (`iu_ref_levels = 2`), $\Delta x = 1.875\cdot10^{-3}$ | 8 blocks of $8\times16\times16$ (`iu_ref_levels = 1`), $\Delta x = 3.75\cdot10^{-3}$ | `rmf-amr-fd-pulse`: $x < 0$ refined by the box marker |
| Cells | 131072 | 16384 | 147456 |

Common: [Gaussian pulse](../initial-conditions#gaussian-pulse) along $+y$, `sigma = 0.01`, `B0 = 10`; `Neumann` outer
faces; 5 steps; 2 ranks (`--np N`), so every seam row crosses ranks. The script arms the monitor with
`seam_divB_tol = 1` in both realms.

## Input files

| File | Role |
|---|---|
| `src/tests/prism/regression/rmf-2realm-fd-pulse-refined/input.ini`, `realm_{1,2}.ini` | the manifest and the realms |
| `src/tests/prism/regression/rmf-amr-fd-pulse/input.ini` | the single-realm reference |
| `src/tests/prism/regression/rmf-2realm-fd-pulse-refined/check.sh` | runs both and compares (`--build`, `--np N`; `PRISM_EXE`) |
| `src/tests/flume/verification/multirealm/multirealm_oracle.py` | the union of the realms against the single realm, `--tol 0` |

```bash
cd src/tests/prism/regression/rmf-2realm-fd-pulse-refined && ./check.sh
PRISM_EXE=$PWD/../../../../../exe/adam_prism_fnl ./check.sh        # FNL
```

## Results

| Leg | Asserts | Bound | CPU | FNL |
|---|---|---|---|---|
| health | both runs reach 100 % with no error/abort/NaN | — | pass | pass |
| fields | every checkpoint field ($\mathbf B$, $\mathbf D$, $\mathbf J$, `div_*`, `res_*`) of the union equals the single realm on every cell | bitwise | 147456 cells, 27 fields | 147456 cells, 27 fields |
| divergence histories | at every step the single-realm max $\lvert\nabla\cdot\mathbf D\rvert$ and max $\lvert\nabla\cdot\mathbf B\rvert$ equal the larger of the two realms' maxima | printed digits | pass | pass |
| seam max $\lvert\nabla\cdot\mathbf B\rvert$, steps 1 → 5 | reported | — | 6.86 → 12.65 | 6.86 → 12.65 |
| monitor | `WARNING: seam div(B) ... exceeds [IO].seam_divB_tol` in the split's log | fires | 15 warnings | 15 warnings |

Related pages: [seam div(B)](./seam-divergence) for the single-realm source; [inter-realm 1:1 seams](./inter-realm-mirror).

## Findings

- The split equals the single realm bit for bit on every field, on both backends: an inter-realm 2:1 seam is an
  intra-realm 2:1 seam across a realm boundary.
- The split carries the single-realm seam source exactly. The coarse realm holds the maximum at steps 1–2, the fine
  realm at steps 3–5.
- The seam $\nabla\cdot\mathbf B$ monitor covers inter-realm 2:1 seams.
