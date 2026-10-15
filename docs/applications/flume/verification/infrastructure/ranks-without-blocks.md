# Ranks without blocks

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| DV-1 | Infrastructure | — ([#82](https://github.com/szaghi/adam/issues/82)) | CPU, FNL | yes | — |

## Case and scope

The tree never splits a family of sibling blocks across ranks. `tree%mpi_redistribute` cuts the Morton-ordered blocks
into equal chunks, and when a cut falls inside a family it moves the whole family to one side. A grid with fewer
families than ranks therefore leaves ranks with no block: four blocks at level 1 are one family, all on rank 0 at np 2.
Every code path that sizes an array by the local `blocks_number` meets a zero there.

On FNL that is a zero-byte device allocation: `acc_malloc(0)` returns a null pointer (measured, nvhpc 25.5), and FUNDAL
reports a null pointer as a failed allocation. Every MHD run therefore died at initialisation
([#82](https://github.com/szaghi/adam/issues/82)): the div(B) history copied its seam-face flags,
`divb_seam(blocks_number, 6)`, to the device. The array is now sized by the block capacity `nb`, as every other
block-indexed device array is. Small multi-rank tests rarely leave a rank empty, so this leg makes the empty rank the
test.

What the leg catches: any block-indexed allocation, copy or reduction that fails, or changes the result, on a rank that
owns no block, in every model and option that allocates block-indexed state; and a regrid that does not repopulate an
empty rank.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Models | Euler; ideal MHD under `none`, `glm` and `eglm`; external gravity $g = (0, -10)$ | [models](../../models#external-gravity) |
| Space | WENO5 flux splitting on characteristic variables; `weno-riemann` HLLD (6th-order correction, `weno` sensor) in one case | `compute_residuals_weno`, `compute_residuals_riemann`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Time | SSP Runge–Kutta, 3 stages (`runge-kutta-ssp-33`), CFL 0.4 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Distribution | sibling families kept on one rank | `mpi_redistribute`, [`adam_tree_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_tree_object.F90) |
| div(B) history | seam-face flags sized by the block capacity `nb` | `set_divb_seam`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90); `compute_divb_history`, [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |
| Options | positivity limiter, `equilibrium = atmosphere`, runtime AMR with a Löhner marker | `limit_positivity`, `regrid_forest`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `mark_by_lohner`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | $[0, 1]^2$, $z$ null |
| Grid | quadtree, four blocks of $16^2$ cells at level 1 (one sibling family), `ngc = 4` |
| Ranks | 1 and $N = 2$ |
| Steps | 20 at CFL 0.4 |
| Boundaries | $x$ periodic; $y$ `hydrostatic` (the atmospheres), `wall-inviscid` (`repro`) |

| Case | Input (`gravity/make_atmosphere.py`) |
|---|---|
| `repro` | the issue's input: `iso2d-mhd` turned into a uniform EGLM state with $\mathbf B = (0.5, 0.7, 0.2)$, $y$ walls, no gravity |
| Euler | `strat2d`: the isothermal atmosphere under $g = (0, -10)$ with a 10% pressure pulse at the centre |
| MHD, `none` / `glm` / `eglm` | `strat2d-mhdx`: `strat2d` with a field across gravity, $\mathbf B = (0.05, 0, 0)$, under each divergence control |
| EGLM, positivity limiter | `strat2d-mhdx` with `positivity_limiter = cell` |
| EGLM, `equilibrium = atmosphere` | `strat2d-mhdx --balanced` |
| EGLM, `weno-riemann` + HLLD | `strat2d-mhdx --riemann hlld` |
| EGLM, runtime AMR | `strat2d-mhdx` with runtime regrids every 5 steps to level 2: a Löhner marker on the density (`refine_tol = 0.1`, `derefine_tol = 0.02`, `buffer = 2`), acting at the regrids only, so rank 1 starts empty and the first regrid repopulates it |

Asserted (`ranks_oracle.py`):
- **ownership:** in the np 2 run every rank but 0 owns no block at step 0, so the leg cannot become vacuous if the
  distribution changes; in the AMR case, every rank owns blocks at the end;
- **identity:** np 2 = np 1 **bitwise**, every field and every cell at the last checkpoint. With every block on rank 0
  the arithmetic is the same, and the reductions add exact zeros.

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/gravity/make_atmosphere.py` | writes the inputs |
| `src/tests/flume/verification/ranks/check.sh` | rewrites them to one family, runs each case on 1 and $N$ ranks (`--np N`, `--cases`, `--keep`, `--it-save N`) |
| `src/tests/flume/verification/ranks/ranks_oracle.py` | `ownership` (`--empty-at-start`, `--filled-at-end`) and `same` (bitwise) |
| `src/tests/flume/verification/make_doc_runs.sh` (group `ranks`) | the AMR case on 1 and 2 ranks, a checkpoint every step |
| `src/tests/flume/verification/make_doc_figures.py` (`fig_ranks`) | the figure |

```bash
cd src/tests/flume/verification/ranks && ./check.sh                    # CPU, np 1 and 2
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh            # FNL
cd .. && ./make_doc_runs.sh <dir> ranks && python3 make_doc_figures.py --out ../../../../docs/public/flume --runs <dir> --only ranks
```

## Results

![DV-1: blocks per rank before and after the first regrid, and the np 2 - np 1 difference](/flume/ranks.png)

*Left: the speed at step 4, before the first regrid, the block outlines coloured by owner: the four blocks are all rank
0's. Centre: step 20; the first regrid (step 5) refined the family into four, and the distribution gave each rank two
families, 8 blocks. Right: the blocks of each rank against the step. Over all 21 checkpoints the largest difference
between the np 2 and np 1 states, over every field and every cell, is 0.*

| Case (FNL, 2 A30 GPUs, nvhpc 25.5) | blocks per rank, step 0 → 20 | before the fix, np 2 | after: np 2 vs np 1 |
|---|---|---|---|
| `repro` | 4 / 0 → 4 / 0 | **aborts** at initialisation | bitwise |
| Euler | 4 / 0 → 4 / 0 | runs | bitwise |
| MHD, `none` | 4 / 0 → 4 / 0 | **aborts** | bitwise |
| MHD, `glm` | 4 / 0 → 4 / 0 | **aborts** | bitwise |
| MHD, `eglm` | 4 / 0 → 4 / 0 | **aborts** | bitwise |
| EGLM, positivity limiter | 4 / 0 → 4 / 0 | **aborts** | bitwise |
| EGLM, `equilibrium = atmosphere` | 4 / 0 → 4 / 0 | **aborts** | bitwise |
| EGLM, `weno-riemann` + HLLD | 4 / 0 → 4 / 0 | **aborts** | bitwise |
| EGLM, runtime AMR | 4 / 0 → 8 / 8 | **aborts** | bitwise (4096 cells) |

| Backend | Bound | Result |
|---|---|---|
| FNL, after the fix | 9/9 bitwise, ownership as asserted | pass, 3.6 min |
| FNL, the pre-fix executable (negative control) | must fail | 8 of 9 cases abort; Euler, which never computes the div(B) history, runs |
| CPU | 9/9 bitwise | pass, 2.8 min |

Related cases: [Scaling covariance](./scaling-covariance) (the other identity leg of this group); the atmospheres are
the [external-gravity](../gravity/atmospheres) inputs.

## Findings

- The issue reported the CPU backend unaffected, since a zero-size host array is legal; it passes the leg as well.
- `divergence_control = none` failed too, because the div(B) history is skipped only for Euler.
- The failure was specific to FNL and invisible to every small multi-rank test that happened to give each rank a block.
