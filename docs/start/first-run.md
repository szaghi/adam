# First run

Build FLUME, run the Sod shock tube, look at what it writes and check it against the exact solution. The toolchain and
HDF5 are set up on the [install page](./install).

## Build

```bash
git clone https://github.com/szaghi/adam && cd adam
fobis fetch                              # the dependencies, into src/third_party/
fobis build --mode flume-cpu-gnu         # exe/adam_flume_cpu (gfortran + Open MPI)
```

## Run

Run in a directory of your own: a run writes its outputs next to its input.

```bash
mkdir -p run/sod && cp src/tests/flume/regression/sod-x/input.ini run/sod/
cd run/sod && mpirun -np 2 ../../exe/adam_flume_cpu input.ini > log.txt
```

The input is the Sod shock tube along $x$: $(\rho, u, p) = (1, 0, 1)$ for $x < 0.5$ and $(0.125, 0, 0.1)$ beyond, 200
cells (64 blocks of $50\times4\times4$ cells, four along $x$; $y$ and $z$ are null axes), WENO5 on characteristic
variables, SSP Runge–Kutta 3, CFL 0.5, to $t = 0.2$. The run ends with

```text
[mpi-00000]save HDF5 files t: 174, time: 0.20000000000000001E+000
[mpi-00000]ADAM run complete: all realms finalized, entering MPI_Finalize
```

## What it writes

```text
sod-x-000000000.xdmf  sod-x-000000000-proc000000.h5  sod-x-000000000-proc000001.h5
sod-x-000000174.xdmf  sod-x-000000174-proc000000.h5  sod-x-000000174-proc000001.h5
sod-x-residuals.dat   sod-x-conservation_history.dat
```

- **Field files**, at step 0 and at the last step (`[IO] it_save = 100000` saves nothing in between): one HDF5 file
  per rank and one XDMF index per save, the conservative variables `r`, `ru`, `rv`, `rw`, `rE` of every block, ghost
  cells included.
- **`sod-x-residuals.dat`**: the L2 norm of every residual, one row per step.
- **`sod-x-conservation_history.dat`**: the domain integrals of the conservative variables, one row per step. With
  extrapolation at both ends mass is conserved up to the flux through the ends, which the waves have not reached:

```text
         0  0.0000000000000000E+000  5.6250000000000000E-001  0.0000000000000000E+000  ...  1.3750000000000002E+000
       174  2.0000000000000001E-001  5.6249999999999989E-001  1.7999999999999999E-001  ...  1.3750000000000002E+000
```

## Look at it

Open `sod-x-000000174.xdmf` in ParaView, apply *Cell Data to Point Data* if you want a line plot, and plot `r` over a
line along $x$ (*Plot Over Line*). The 64 blocks come in as one multiblock dataset; the ghost cells overlap at the block
faces. On Windows with the files on the WSL share, set the Windows user variable `HDF5_USE_FILE_LOCKING=FALSE` first,
or ParaView shows empty data ([IO cookbook](/cookbook/io)).

## Check it

FLUME's Sod verification compares the density with the exact Riemann solution:

```bash
python3 -m venv ../../exe/.regression-venv && ../../exe/.regression-venv/bin/pip install h5py numpy   # once
../../exe/.regression-venv/bin/python ../../src/tests/flume/verification/sod/sod_oracle.py input.ini . --l1-max 3.31e-3
```

```text
/sod-x: axis x, it 174, cells 200, L1(rho) = 3.244140e-03  PASS
```

$3.244\cdot10^{-3}$ is the value recorded on the [Sod verification page](/applications/flume/verification/euler/sod).

::: tip What you learned
- A run is one executable and one INI file; it writes next to the input.
- Field files are per rank, indexed by one XDMF per save; histories are plain text.
- Every FLUME verification case has an oracle you can run on your own output.

Next: the [FLUME tutorial](/tutorials/flume/sod), or the other [tutorials](/tutorials/).
:::
