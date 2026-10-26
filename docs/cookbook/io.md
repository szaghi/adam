# Input and output

Recipes for the files a run writes: the field checkpoints, the restart files, the histories and the slices, and how to
read each of them. The keys are in the [library input reference](/reference/library-input#io); the application pages
list their own histories ([FLUME](/applications/flume/input#output-formats-no-format-key-exists),
[PRISM](/applications/prism/#outputs)).

## IO1. What a run writes

**Problem.** Know which files a run produces and which key controls each.

**Input** (`src/tests/flume/regression/sod-x/input.ini`):

```ini
[IO]
output_basename        = sod-x
it_save                = 100000
restart                = .false.
restart_basename       = sod-x-restart
restart_save           = 0
residuals_save         = 1
save_memory_status     = .false.
save_residual_fields   = .false.
save_auxiliary_fields  = .false.
save_curl_fields       = .false.
save_divergence_fields = .false.
save_gradient_fields   = .false.
save_laplacian_fields  = .false.
```

**What you get.**

| Files | Content | Written |
|---|---|---|
| `<output_basename>-<it:9>.xdmf` + `<output_basename>-<it:9>-procNNNNNN.h5` | The fields of every block, **ghost cells included**, one HDF5 file per rank and one XDMF index per save | step 0, every `it_save` steps and the last step (`it_save > 0`; a restarted FLUME run skips its first) |
| `<restart_basename>.time`, `.tnd`, `-procNNNNNN.fbd` | Restart: the time, the tree nodes, the blocks of each rank | every `restart_save` steps (FLUME: also the last step) |
| `<restart_basename>.xdmf` + `<restart_basename>-procNNNNNN.h5` | A snapshot of the restart state, interior cells only, viewable as a field file (origin at the first interior cell, [#109](https://github.com/szaghi/adam/issues/109)) | with every restart |
| FLUME: units and sinks restart files; PRISM: PIC restart files | The app's extra restart state | with every restart |
| `<output_basename>-residuals.dat` | L2 norm of every residual, one row per step | every `residuals_save` steps |
| application histories | FLUME: `-conservation_history.dat`, `-divb_history.dat` (MHD); PRISM: `-energy_history.dat`, `-divergence_history.dat`, PIC and probe histories | see the application pages |
| `<output_basename>-slice_NN-<it:9>.mat` | Slices (recipe IO5) | `[slice_N] n_save` |

**Reference.** [`[IO]`](/reference/library-input#io); the
[FLUME output table](/applications/flume/input#output-formats-no-format-key-exists); the
[PRISM outputs](/applications/prism/#outputs).

## IO2. Open the fields in ParaView

**Problem.** Look at a checkpoint.

**Do.** Open `<output_basename>-<it:9>.xdmf` in ParaView. It indexes the
`-procNNNNNN.h5` files of every rank; each block is a uniform Cartesian grid, and the ghost cells are part of it, so
neighbouring blocks overlap by `2 ngc` cells. To hide the overlap, apply *Extract Subset* or colour by cell data and
ignore the block edges.

::: warning Windows ParaView reading files on the WSL file system
A Windows ParaView opening an XDMF under `\\wsl$\<distro>\...` (or `\\wsl.localhost\...`, or a drive mapped to them)
loads the file **with no warning and empty data**: HDF5 1.10+ locks every file it opens, and the lock fails on the
WSL file server. Set the Windows **user** environment variable `HDF5_USE_FILE_LOCKING=FALSE` and restart ParaView. It is
safe for read-only viewing and needs no copy to NTFS. The same applies to any Windows HDF5 reader (VisIt, h5py on
Windows).
:::

**Reference.** IO1 for the file names.

## IO3. Read a checkpoint in Python

**Problem.** Post-process the fields with NumPy.

**Layout.** In every `-procNNNNNN.h5` each block `b` of the rank has four datasets, named
`block_<b:9>-proc<rank:6>-<name>`:

| Dataset | Shape | Content |
|---|---|---|
| `...-origin` | `(3,)` | the block's minimum corner **in (z, y, x) order**, at the corner of the ghost cells: `emin - ngc dx` |
| `...-dxdydz` | `(3,)` | the cell spacing, (z, y, x) order |
| `...-time_iteration` | `(1,)` | the step |
| `...-<variable>` | `(nk + 2 ngc, nj + 2 ngc, ni + 2 ngc)` | one array per variable, (z, y, x) order |

**Code** (the reader of `src/tests/ripple/oracle.py`, which holds for FLUME and PRISM field files):

```python
import h5py

with h5py.File("sod-x-000000000-proc000000.h5", "r") as h5:     # step 0; sod-x also saves its last step
    for blk in sorted({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")}):
        origin = h5[f"{blk}-origin"][()][::-1]        # (x, y, z), ghost corner
        dxyz = h5[f"{blk}-dxdydz"][()][::-1]          # (x, y, z)
        r = h5[f"{blk}-r"][()].transpose(2, 1, 0)      # density r[i, j, k], ghosts included
```

The interior of a block is `r[ngc:-ngc, ngc:-ngc, ngc:-ngc]`. The variable names are the application's: FLUME writes
`r`, `ru`, `rv`, `rw`, `rE` (plus the MHD and scalar slots, and the [derived output
fields](/applications/flume/input#derived-output-fields) with `save_auxiliary_fields`); PRISM writes `Dx`, `Dy`, `Dz`,
`Bx`, `By`, `Bz`, then `phi`, `psi`, `Jx`, `Jy`, `Jz` and `rho` as its model has them.

**Reference.** `open_block_xh5f` in `src/lib/common/adam_realm_object.F90`.

## IO4. Restart a run

**Problem.** Stop a run and continue it later.

**Input.** Write restarts during the first run, then resume with `restart = .true.` and the same basename:

```ini
[IO]
restart          = .true.
restart_basename = sod-x-restart
restart_save     = 100
```

**Rules** (`load_restart_files`, `src/lib/common/adam_adam_object.F90`; `load_blocks`,
`src/lib/common/adam_field_object.F90`):
- **Same rank count.** Each rank reads its own `<restart_basename>-procNNNNNN.fbd`. A missing file stops the run
  (`file ... does not exist`, error code -104), so a restart on **more** ranks stops. A restart on **fewer** ranks
  stops too: the tree file names the rank of every block, and a block of a rank this run does not have stops the run
  ("saved on more ranks than ... of this run", `tree%load_nodes`, [#109](https://github.com/szaghi/adam/issues/109)).
  A missing `.tnd` stops the run.
- **Same block shape.** `nv`, `ni`, `nj`, `nk` and `ngc` must equal the saved ones; a mismatch stops the run (-103).
- **Same backend budget.** A rank that saved more blocks than the restarted run can allocate stops (-102).
- The time and the step come from `<restart_basename>.time`, the tree from `.tnd`. The residuals history is appended,
  so a run split by restarts has the history of the uninterrupted one.
- In a forest each realm restarts from its own files: give every realm its own `restart_basename`
  ([forest R9](/cookbook/forest#r9-restart)).

**What you get.** The run continues from the saved step; FLUME's multi-realm restart is bitwise identical to the
straight run (MV-14, [forest R9](/cookbook/forest#r9-restart)).

**Reference.** [`[IO] restart`, `restart_basename`, `restart_save`](/reference/library-input#io).

## IO5. Sample the solution on a line or a plane: slices

**Problem.** Record the state on a lattice of points at a cadence, without writing full checkpoints.

**Input** (`src/tests/prism/regression/slices/slices.ini`):

```ini
[slices]
slices_number = 4

[slice_1]
itype  = trilinear
n_save = 5
ni     = 97
nj     = 1
nk     = 1
emin_x = -0.06
emin_y = 0.0013
emin_z = 0.0013
emax_x = 0.06
emax_y = 0.0013
emax_z = 0.0013
```

**What you get.** `<output_basename>-slice_01-<it:9>.mat` every `n_save` steps and at the last step: one file for
all ranks, written with MPI-IO. The format (`save_slice`, `src/lib/common/adam_adam_object.F90`):
- a header of four 32-bit integers, `ni, nj, nk, nc`, with `nc` the columns per point;
- then, for every point with `i` fastest, `nc = 3 + nv` 64-bit reals: `x, y, z` and the state variables (`save_slice` can
  also write the level set `phi`, but no caller passes it).

Reading it in Python (`src/tests/prism/regression/slices/slices_convergence.py`):

```python
import numpy as np

raw = open("slices-slice_01-000000005.mat", "rb").read()
ni, nj, nk, nc = np.frombuffer(raw[:16], dtype=np.int32)
data = np.frombuffer(raw[16:], dtype=np.float64).reshape(nk, nj, ni, nc)
x, dz = data[0, 0, :, 0], data[0, 0, :, 5]          # x and Dz (columns x, y, z, Dx, Dy, Dz, ...) along the line
```

To Tecplot ASCII: ASCOT, `src/app/ascot/ascot.F90` (`ascot -i <file>.mat -o <file>.dat`, `-v` for the variable names), which
`scripts/slices_converter.py mat2dat` calls for a whole series, followed by `dat2png` and `png2mp4` for a movie.
ASCOT's modes are `ascot-nvf-cuda` and `ascot-omp-amd` ([build modes](/reference/build-modes)); the program itself uses
only PENF.

**Reference.** [`[slices]`, `[slice_N]`](/reference/library-input#slices-and-slice-n); the slice oracle
`src/tests/prism/regression/slices/check.sh`.
