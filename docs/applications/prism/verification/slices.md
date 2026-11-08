# Slices output

| Cases | Issues | Type | Backends | Figure | Regression |
|---|---|---|---|---|---|
| `slices` | — | C | CPU, FNL | no | — |

## Case and scope

The `[slices]` output writes the fields interpolated on lines, planes and boxes into MPI-IO `.mat` files. No run reads
them back, so a wrong value goes unnoticed. The case recomputes every written value from the checkpoint of the slice
step and checks the decomposition independence, the non-dimensional scaling, the accuracy against an analytic wave and
the input validation of the `[slice_N]` sections.

What the case catches:
- a slice value that is not the trilinear (or inverse-distance) combination of the right eight cells, a point written
  by no rank, a point in the wrong place;
- stale auxiliary ghosts: the curl and the divergence are computed on the interior cells only, and `save_slices` must
  refresh their intra-realm ghosts before interpolating;
- an output that depends on the rank count (each point is written by the one rank owning it);
- slice extents left dimensional under a non-dimensional model: the points were then in the wrong place, or outside
  the domain and the run aborted at the first save;
- a missing or invalid `[slice_N]` option that left the slice undefined (`n_save = 0` crashed on `mod(it, 0)`);
- a stale, longer slice file not truncated by the save.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Case | the source-free pulse of `rmf-amr-fd-pulse` (2:1 seam at $x = 0$, 36 blocks), curl and divergence saved | [seam div(B)](./seam-divergence) |
| Interpolation | `trilinear` or `inverse_distance` over the eight cells around the point, ghost cells included | `interpolate_at_point`, `save_slice`, [`adam_adam_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_adam_object.F90); [`adam_slices_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_slices_object.F90) |
| Ghost refresh | auxiliary fields refreshed before the save | `save_slices`, [`adam_prism_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/prism/common/adam_prism_common_object.F90) |
| `.mat` layout | 4 `int32` ($n_i$, $n_j$, $n_k$, $3 + n_q$), then per point 3 coordinates and $n_q$ values in `float64`, $i$ fastest | — |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain, grid | $[-0.06, 0.06]^3$ m; 8 blocks of $16^3$ at level 1, $x < 0$ refined to level 2 |
| Initial condition | [Gaussian pulse](../initial-conditions#gaussian-pulse) along $+y$, `sigma = 0.01`, `B0 = 10` |
| Steps | 5; slices saved every 5 steps |
| `slice_1` | `trilinear` line, 97 points along $x \in [-0.06, 0.06]$ at $y = z = 0.0013$ (through the seam) |
| `slice_2` | `inverse_distance` plane, $33\times29$, $x \in [-0.05, 0.05]$, $y \in [-0.04, 0.04]$, $z = 0.0007$ |
| `slice_3` | `trilinear` box, $7^3$, $[-0.004, 0.004]^3$, around the block corner at the origin |
| `slice_4` | `inverse_distance` plane, $16\times32$, $x \in [0, 0.06]$, $y \in [-0.06, 0.06]$, $z = 0.001875$: exactly on cell centres |
| Ranks | 1, 2, 4 (`--np`) |
| Non-dimensional leg | `AdimEM`, `L0 = 0.01`, `B0 = 10`, a plane wave along $y$ with $\lambda = 0.04$; the last rank count |
| Convergence leg | the same plane wave, no AMR, 16 and 32 cells per block side, 1 step, one line slice off the cell centres |

## Input files

| File | Role |
|---|---|
| `src/tests/prism/regression/slices/slices.ini` | the case and its four slices |
| `src/tests/prism/regression/slices/check.sh` | the runs and the legs (`--build`, `--np "N ..."`, `--convergence`; `PRISM_EXE`) |
| `src/tests/prism/regression/slices/slices_oracle.py` | recomputes every slice value; checks the face ghosts |
| `src/tests/prism/regression/slices/slices_convergence.py` | the slice error against the analytic plane wave |

```bash
cd src/tests/prism/regression/slices && ./check.sh --convergence
PRISM_EXE=$PWD/../../../../../exe/adam_prism_fnl ./check.sh --np "1 2"   # FNL
```

## Results

| Leg | Asserts | Bound | CPU | FNL |
|---|---|---|---|---|
| oracle, every slice and field set (q, curl, div), per rank count | coordinates against the `[slice_N]` box; values against the recomputed interpolation, relative to the largest of the eight stencil values | coordinates $10^{-14}$; values $10^{-12}$ | ⟨run⟩ | ⟨run⟩ |
| oracle, face ghosts | ghost slab equal to the same-level neighbour's interior slab, every field set (faces across the 2:1 seam counted, not checked) | bitwise | ⟨run⟩ | ⟨run⟩ |
| decomposition | every `.mat` file at 2 and 4 ranks identical to the 1-rank file | bitwise | ⟨run⟩ | ⟨run⟩ |
| truncation | the pre-placed 100000-byte file replaced by a file of the right size | file size | ⟨run⟩ | ⟨run⟩ |
| non-dimensional | the oracle, with the extents divided by `L0` | as the oracle | ⟨run⟩ | ⟨run⟩ |
| `--convergence` | observed order of the slice error, 16 → 32 cells | $\ge 1.8$ | ⟨run⟩ | ⟨run⟩ |

Input validation (1 rank): each leg must stop with the message and write no slice.

| Leg | Change | Message |
|---|---|---|
| `missing-n_save` | `n_save` deleted | `failed to load [slice_1].(n_save)` |
| `n_save-zero` | `n_save = 0` | `[slice_1].(n_save) = 0 must be >= 1` |
| `itype-unknown` | `itype = bilinear` | `unknown [slice_1].(itype) "bilinear"` |
| `nj-zero` | `nj = 0` | `[slice_1].(ni, nj, nk)` |
| `emin-above-max` | `emin_y = 0.002` | `[slice_1] emin` |
| `number-negative` | `slices_number = -1` | `[slices].(slices_number) = -1 must be >= 0` |
| `outside-domain` | `emax_x = 0.07` | `outside the realm domain` |
| `itype-crlf` | `itype = trilinear` with a trailing CR | accepted: the run completes and saves the slice |

Related pages: [refusals](./refusals) for the input validation of the other sections.

## Findings

- The slices are recomputed point by point from the checkpoint, so the oracle covers the interpolation weights, the
  stencil choice, the owner of each point and the ghost refresh at once.
- A point on a block face (on the seam, $x = 0$) is owned by one side or the other by its last bits (gfortran writes
  0, nvfortran $-2.7\cdot10^{-18}$): the oracle recomputes each value at the written coordinates.
- After one step the solver error is far below the interpolation error of a line slice off the cell centres, so the
  convergence leg measures the trilinear error, $O(h^2)$.
