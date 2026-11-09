# Restart round trip and output channels

| ID | Group | Milestone | Type | Backends | Figure | Regression |
|---|---|---|---|---|---|---|
| V7 | Euler | M1 ([#35](https://github.com/szaghi/adam/issues/35)) · P7; [#42](https://github.com/szaghi/adam/issues/42), [#109](https://github.com/szaghi/adam/issues/109) (follow-up) | C | CPU, FNL | no | — |

## Case and scope

Restart files, slices and auxiliary fields are outputs that no run reads back in normal operation, so a defect in them
is silent. The case runs each case twice: run A goes $N$ steps; run B goes $N/2$ steps saving a restart, restarts and
completes $N$. The two final states must be **bitwise** identical on the interior cells, and the residual and
conservation histories byte-identical: the restarted run appends to them and does not re-save its starting step.

Four cases cover the three integration paths of FLUME and an empty rank:
- `sod-x`, the single-realm fast path;
- `amr-periodic`, the staged path with static 2:1 AMR and the reflux;
- `shock-cylinder`, the immersed boundary with AMR;
- `uniform`, the skeleton case whose 8 blocks all sit on rank 0 at np 2, so rank 1 owns no block when the restart is
  saved.

What the case catches:
- **a restart channel that misses part of the state:** any field, time, step counter or history line not saved or not
  read back breaks the bitwise identity;
- **an empty rank** that writes no blocks file, so the restart aborts ([#42](https://github.com/szaghi/adam/issues/42));
- **a restart snapshot placed off its cells:** the XH5F snapshot written with every restart holds the interior cells
  only, and its block origin must be the first stored cell, not the ghost corner
  ([#109](https://github.com/szaghi/adam/issues/109));
- **a restart on fewer ranks than saved,** which must stop instead of dropping the missing ranks' blocks
  ([#109](https://github.com/szaghi/adam/issues/109));
- **wrong derived outputs:** the auxiliary fields against the conservative ones, a slice against the cell values.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Restart files | `<restart_basename>.time`, `.tnd`, `-procNNNNNN.fbd` per rank, and the XH5F snapshot `<restart_basename>.xdmf` + `.h5` (interior cells) | `save_restart_files`, `load_restart_files`, [`adam_adam_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_adam_object.F90); [input](../../input#io) |
| Rank check | a restart whose blocks file was saved on more ranks than the run has stops | [`adam_tree_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_tree_object.F90) |
| Auxiliary fields | `[IO] save_auxiliary_fields`: $u, v, w, p, H, a$ in every checkpoint | [input](../../input#io) |
| Slices | `[slices]`, one trilinear slice | library `slices_object`, `save_slice`, [`adam_adam_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_adam_object.F90) |
| Cases | the numerics of [Sod](./sod) (V1), [conservation across AMR seams](./amr-conservation) (V3) and [shock over a cylinder](./shock-cylinder) (V6) | — |

## Initial and boundary conditions

| Leg | Case | Blocks | Run A | Run B | Checks |
|---|---|---|---|---|---|
| `v7` | `sod/sod-x.ini` | 64 | 40 steps | 20 + restart + 20 | interior fields bitwise, histories byte-identical; the snapshot of step 20 against the regular save of step 20; a restart on np − 1 ranks |
| `v7` | `conservation/amr-periodic.ini` | 15 | 20 steps | 10 + restart + 10 | interior fields bitwise, histories byte-identical |
| `v7` | `shock-cylinder/shock-cylinder.ini` | 176 | 40 steps | 20 + restart + 20 | interior fields bitwise, histories byte-identical |
| `v7` | `uniform/input.ini` | 8, all on rank 0 | 10 steps | 5 + restart + 5 | interior fields bitwise, histories byte-identical (the state is uniform: the leg checks that an empty rank saves and restarts) |
| `slices` | `sod/sod-x.ini`, `save_auxiliary_fields = .true.`, one trilinear slice of $200 \times 1 \times 1$ points on the $x$ cell centres | 64 | 40 steps | — | auxiliary fields against the conservative ones; slice against the cell values |

Every run on 2 ranks (`--np N`). Interior cells only: the edge and corner ghost cells that no map fills keep stale values
the directional stencils never read.

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/io/check.sh` | runs the legs (`--legs v7,slices`, default both; `--np N`; `--build` builds `flume-cpu-gnu`; `FLUME_EXE` selects the backend) |
| `src/tests/flume/verification/conservation/conservation_oracle.py --compare A B --tol 0 --ngc 3` | runs A and B cell by cell, tolerance 0 |
| `src/tests/flume/verification/io/snapshot_oracle.py` | the snapshot block origins against the regular save of the same step |
| `src/tests/flume/verification/io/io_oracle.py` | the auxiliary fields and the slice (relative tolerance $10^{-12}$) |

```bash
cd src/tests/flume/verification/io && ./check.sh                       # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh            # FNL
```

No figure: the legs are a bitwise identity and an output contract.

## Results

| Leg | Check | Bound | CPU | FNL |
|---|---|---|---|---|
| `v7` | `sod-x`: interior fields, A against B | bitwise | bitwise | bitwise |
| `v7` | `amr-periodic`: interior fields | bitwise | bitwise | bitwise |
| `v7` | `shock-cylinder`: interior fields | bitwise | bitwise | bitwise |
| `v7` | `uniform`, an empty rank at the save: interior fields | bitwise | bitwise | bitwise |
| `v7` | residual and conservation histories, every case | byte-identical | identical | identical |
| `v7` | restart snapshot origin against the first stored cell (64 blocks) | 0 cells ($< 10^{-9}$) | 0 | — |
| `v7` | restart on 1 of 2 ranks | stops with the rank message | stops | — |
| `slices` | auxiliary fields $u, v, w, p, H, a$ against the recomputed ones | $\le 10^{-12}$ relative | 0 | $1.9\cdot10^{-16}$ |
| `slices` | slice points on the cell centres; slice against the cell values | 0 off; $\le 10^{-12}$ | 0; $5.6\cdot10^{-16}$ | 0; $5.6\cdot10^{-16}$ |

The FNL column of the two [#109](https://github.com/szaghi/adam/issues/109) rows has not been run since the fix.

Related cases: restart across runtime regrids is [runtime regrids](../runtime-amr/runtime-regrids) (rg3–rg6); the MHD
restart round trips are MV-14 legs 2 and 3 on [forest seams](./multi-realm); the forest restart is
[runtime AMR on a forest](../runtime-amr/multi-realm).

## Findings

- A restart is invisible: on every integration path, and with a rank that owns no block, the restarted run ends bitwise
  on the straight one, and its histories continue byte for byte.
- The restart snapshot is a regular XH5F file of the interior cells: a reader places it on the cells of the field
  saves of the same step.
- A restart needs the rank count it was saved with; fewer ranks stop the run.
- The auxiliary fields and the slices reproduce the state they are derived from to round-off.
