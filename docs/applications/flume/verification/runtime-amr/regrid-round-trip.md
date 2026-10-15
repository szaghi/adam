# Regrid round trip and input contract

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| RG (rg0, rg1, rg2) | Runtime AMR (M5) | M5 ([#74](https://github.com/szaghi/adam/issues/74)) | CPU (rg0, a host program); rg1 and rg2 CPU, FNL | no | — |

## Case and scope

Until M5 nothing ran `amr_update` after initialisation and no test ever derefined. rg0 exercises the library regrid
step on its own: a block is refined, then its children derefined, on a field known in every cell; rg1 and rg2 pin the
two input keys of runtime AMR.

What the case catches:
- **inconsistent bookkeeping over the ranks:** the redistribution moves blocks, and the leaves and the blocks summed
  over the ranks must agree after every step;
- **an inexact transfer where an exact one is possible:** a linear field must be reproduced by the refine; with the
  conservative prolongation the round trip must return the initial state and $\sum q\,\Delta V$ must not change;
- **a prolongation that leaves the range of the data:** a steep front must stay in $(0.01, 1)$, so a positive
  variable stays positive;
- **a silent input:** a negative `frequency`, a `frequency > 0` without markers, an unknown `regrid_prolongation`.

### The conservative prolongation

Each parent cell is split along the refined axes into children set, variable by variable, to

$$
q_c = q + \phi \sum_d \frac{\sigma_d\, s_d}{4}, \qquad \sigma_d = \pm 1,
$$

with $s_d$ the monotonized-central slope of the parent along axis $d$ (in parent-cell units, from its two face
neighbours) and $\phi \in [0, 1]$ the largest factor that keeps every child in the range of the parent and its face
neighbours (Barth and Jespersen 1989). The offsets cancel over the children, so their mean is the parent: the mean
restriction undoes the prolongation exactly and $\sum q\,\Delta V$ does not change. The MC limit bounds each axis on
its own, not the sum over the axes a corner child takes; $\phi$ closes that gap. Linear data keep $\phi = 1$ and are
reproduced exactly.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Regrid step | tree adapted with 2:1 balance, interior data prolonged or restricted, blocks redistributed, maps rebuilt | `amr_update`, `check_blocks_number`, [`adam_adam_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_adam_object.F90); `adapt`, [`adam_field_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_field_object.F90) |
| `linear` prolongation | tensor linear (weights 1/4, 3/4); the library default | `refine2D`, `refine3D`, [`adam_field_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_field_object.F90) |
| `conservative` prolongation | limited MC slopes with the Barth–Jespersen factor $\phi$ (#74 D-M5-2); the FLUME default | `refine_conservative`, [`adam_field_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_field_object.F90) |
| Restriction | the mean of the children | `derefine2D`, `derefine3D`, [`adam_field_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_field_object.F90) |
| `[amr] regrid_prolongation` | read by the library; FLUME resolves an absent key to `conservative` and logs it | `load_from_file`, `prolongation_name`, [`adam_amr_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_amr_object.F90); `check_runtime_amr`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90) |
| `[amr] frequency` | 0: the AMR of the initial condition only; $n > 0$: a regrid every $n$ steps, which needs markers | `check_runtime_amr`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90); [input](../../input#setting-up-runtime-amr) |

The regrid algorithms are described in [numerics](../../numerics#adaptive-mesh-refinement).

## Initial and boundary conditions

| Leg | Setup |
|---|---|
| rg0 | `tests/amr/test_amr_regrid_roundtrip` on 1, 2 and 3 ranks. Domain $[0, 1]^3$; a realm refined uniformly to level 1: octree (`ratio = 8`, $8^3$ cells per block, 8 blocks) and quadtree (`ratio = 4`, $8 \times 8 \times 1$ cells, $z$ null, 4 blocks), `ngc = 3`. The field is set in every cell, ghosts included (a perfect ghost fill, so the transfers are tested alone). The block at the origin is refined, then its children derefined, with each prolongation |
| rg0 fields | linear $f = 1 + 0.3x + 0.7y + 0.5z$; quadratic, $f + 2x^2 + 1.5y^2 + xy$; a steep front from 1 to 0.01, $0.505 + 0.495\tanh((c - s)/0.03)$, a tanh half a coarse cell wide across $s = x+y+z = c = 0.45$ on the octree, $s = x+y = c = 0.3$ on the quadtree |
| rg1 | the regression `sod-x` input ([Sod](../euler/sod)) to $t = 0.002$, 2 ranks, with `[amr] frequency` = 0, −1, 5 (and no markers) |
| rg2 | the same input with `[amr] regrid_prolongation` absent, `linear`, `conservative`, `bogus` |

## Input files

| File | Role |
|---|---|
| `src/tests/amr/test_amr_regrid_roundtrip.F90` | the rg0 unit test (`fobis build --mode test-amr-regrid-roundtrip-gnu`) |
| `src/tests/flume/verification/regrid/check.sh --leg rg0`, `--leg rg1`, `--leg rg2` | runs the legs (`FLUME_EXE` selects the backend of rg1 and rg2) |
| `src/tests/flume/regression/sod-x/input.ini` | the base input of rg1 and rg2 |

```bash
cd src/tests/flume/verification/regrid && ./check.sh --leg rg0 --leg rg1 --leg rg2         # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh --leg rg1 --leg rg2           # FNL
```

No figure: the legs are a transfer identity and an input contract (#85).

## Results

rg0 asserts:
- the leaves and the blocks summed over the ranks agree: octree 8 → 15 → 8, quadtree 4 → 7 → 4;
- a linear field exact after the refine, with both prolongations;
- with `conservative`, the round trip returns the initial state and $\sum q\,\Delta V$ is unchanged after the refine
  and after the round trip, on every field, to round-off ($10^{-13}$ of the field scale); the front stays in
  $(0.01, 1)$;
- with `linear`, the quadratic field drifts (the negative control).

Every leg passes on 1, 2 and 3 ranks, with the same numbers (CPU):

| Field | Tree | `linear`: drift of $\sum q\,\Delta V$ | `linear`: round trip | `conservative`: drift | `conservative`: round trip |
|---|---|---|---|---|---|
| quadratic | octree | $1.3\cdot10^{-4}$ | $3.4\cdot10^{-4}$ | 0 | $4\cdot10^{-17}$ |
| quadratic | quadtree | $2.7\cdot10^{-4}$ | $3.4\cdot10^{-4}$ | 0 | $2\cdot10^{-17}$ |
| front | octree | $2.6\cdot10^{-2}$ | $1.8\cdot10^{-1}$ | $4.8\cdot10^{-15}$ | $1\cdot10^{-16}$ |
| front | quadtree | $1.8\cdot10^{-2}$ | $1.3\cdot10^{-1}$ | $2.6\cdot10^{-16}$ | $1\cdot10^{-16}$ |

Drift is relative to the initial integral, after the refine; round trip is the largest $|q - f|$ over the field scale
after refine then derefine.

| Leg | Input | Expected | CPU | FNL |
|---|---|---|---|---|
| rg1 | `frequency = 0` | runs | pass | pass |
| rg1 | `frequency = -1` | fatal: "must be 0 (no runtime regridding) or positive" | pass | pass |
| rg1 | `frequency = 5`, no markers | fatal: "nothing would mark a block" | pass | pass |
| rg2 | `regrid_prolongation` absent | runs, resolved to `conservative` | pass | pass |
| rg2 | `linear`, `conservative` | runs, taken as given | pass | pass |
| rg2 | `bogus` | fatal: `unknown [amr].(regrid_prolongation) "bogus"` | pass | pass |

The FNL column of rg1 and rg2 is read from the logs of the `work-adam_flume_fnl-*` directories of `regrid/check.sh`.

Related cases: the regrids during a FLUME run are [runtime regrids](./runtime-regrids); the inviscid conservation of
static AMR is [Euler AMR conservation](../euler/amr-conservation).

## Findings

- The `linear` prolongation is second order (a quadratic field's children miss by $2.6\cdot10^{-4}$ of the scale,
  the conservative one's by $1.1\cdot10^{-4}$) but changes the integral at every refine.
- The bound scaling is load-bearing: with $\phi = 1$ forced (a mutation run), an octree child of the front reaches
  $-0.094$, below zero, and the leg fails; on the quadtree the MC limit alone suffices for this front.
