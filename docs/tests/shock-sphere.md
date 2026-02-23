# Shock-Sphere Test

> Shock-sphere interaction test — NASTO NVF backend.

This test performs a shock-sphere interaction simulation using ADAM NASTO compiled with the NVF backend. The test is configured to run quickly on a single NVIDIA GPU.

## Prerequisites

Compile NASTO NVF app — see [NVF Backend documentation](/applications/nasto/nvf).

## Running the Test

The test directory contains:

```
.
├── adam_nasto_nvf -> exe/adam_nasto_nvf      # symlink to compiled executable
├── adam-nasto-shock-sphere.ini               # simulation configuration
├── check-schlieren.png                       # reference image for validation
├── check-schlieren.pvsm                      # ParaView state file
├── clean.sh                                  # cleanup script
└── run_test.sh                               # test runner
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

After a successful run, HDF5 output files are produced at multiple time steps. Visualize results with ParaView using the provided `check-schlieren.pvsm` state file. The results should show the schlieren visualization of the shock-sphere interaction.
