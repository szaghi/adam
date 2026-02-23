# Sod-Z Test

> 1D Riemann Problem of Sod along Z axis — NASTO NVF backend.

This test solves a 1D Riemann Problem of Sod along the Z axis using ADAM NASTO compiled with the NVF backend. The test is configured to run quickly on a single NVIDIA GPU.

## Prerequisites

Compile NASTO NVF app — see [NVF Backend documentation](/applications/nasto/nvf).

## Running the Test

The test directory contains:

```
.
├── adam_nasto_nvf -> exe/adam_nasto_nvf   # symlink to compiled executable
├── adam-nasto-sod-x.ini                   # simulation configuration
├── check-sod-x.pvsm                      # ParaView state file
├── clean.sh                               # cleanup script
└── run_test.sh                            # test runner
```

### Clean

```bash
./clean.sh
```

### Run

```bash
./run_test.sh
```

Output is written to stdout and saved to `run_test.log`.

### Debug mode

```bash
./run_test.sh -debug
```

Runs via `compute-sanitizer` for additional diagnostics.

## Expected Results

After a successful run, HDF5 output files are produced. Visualize results with ParaView using the provided state file.
