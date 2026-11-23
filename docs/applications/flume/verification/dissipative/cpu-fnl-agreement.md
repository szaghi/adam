# CPU and FNL agree run by run

| ID | Group | Milestone | Type | Backends | Figure | Regression |
|---|---|---|---|---|---|---|
| VV `agree` | Dissipative terms (M4) | M4 ([#65](https://github.com/szaghi/adam/issues/65)) | C | CPU, FNL | no | — |

## Case and scope

Every run of the VV legs on the CPU against the same run on FNL, cell by cell on the last checkpoint.

What the case catches: the legs pass on both backends, but a pass bounds an error norm, and two backends can pass with
different fields. Every dissipative kernel exists twice, as a CPU loop and an FNL device kernel (OpenACC, the
`flume-fnl-nvf` build), so the same run on both must give the same state up to the round-off of the operation order. A device
kernel that differs from its CPU twin in a term below the oracle's resolution shows here.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| CPU dissipative kernels | Navier–Stokes and Ohmic fluxes, orders 2 and 4 | `add_dissipative_fluxes_o2`, `add_dissipative_fluxes_o4`, `add_resistive_fluxes_o2`, `add_resistive_fluxes_o4`, [`adam_flume_cpu_dissipation_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_dissipation_kernels.F90) |
| FNL dissipative kernels | the device twins | `add_dissipative_fluxes_o2_dev`, `add_dissipative_fluxes_o4_dev`, `add_resistive_fluxes_o2_dev`, `add_resistive_fluxes_o4_dev`, [`adam_flume_fnl_dissipation_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_dissipation_kernels.F90); `add_dissipative_fluxes_dev`, [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |
| Comparison | blocks matched by their origin (rank independent); each difference divided by the physical scale of its variable, not by its own maximum | `agree`, `src/tests/flume/verification/viscous/waves.py` |

The scales: density and energy by their maximum; the momentum by the larger of its maximum and $\rho_0 a_0$ (an Ohmic
wave at rest carries momentum of $O(A^2)$, far below the round-off of the state); the field by its maximum (it is the
signal); $\psi$ by $c_h \max|\mathbf B|$.

## Initial and boundary conditions

The runs of the legs, unchanged: 72 runs, every resolution of every ladder, the $-A$ twins, the ideal twins and both
sides of each forest, on 2 ranks, one clean build of each backend. The setups are on the case pages:

| Leg | Runs | Setup |
|---|---|---|
| VV-1 shear | 12 | [viscous waves](./viscous-waves) |
| VV-2 acoustic | 6 | [viscous waves](./viscous-waves) |
| VV-3 Couette | 3 | [Couette](./couette) |
| VV-4 Becker | 6 | [Becker](./becker-shock) |
| VV-5 Ohmic, Alfvén, acoustic on MHD | 18 | [Ohmic decay](./ohmic-alfven) |
| VV-6 Ohmic budget | 4 | [Ohmic heating](./ohmic-heating) |
| VV-7 conservation | 6 | [seam conservation](./seam-conservation) |
| VV-7 forests | 8 | [seam conservation](./seam-conservation) |
| VV-8 seams | 9 | [seam accuracy](./seam-accuracy) |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/viscous/check.sh --legs agree` | runs nothing: pairs every `work-adam_flume_cpu-np2-*` directory with its `work-adam_flume_fnl-np2-*` twin, after every leg has run on both executables |
| `src/tests/flume/verification/viscous/waves.py agree` | the cell-by-cell comparison |

```bash
cd src/tests/flume/verification/viscous
./check.sh && FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh && ./check.sh --legs agree
```

No figure: a run-by-run identity has nothing a plot would show (#85).

## Results

Bound: $10^{-10}$ on every run.

| Leg | Runs | Largest difference |
|---|---|---|
| VV-1 shear | 12 | $5.8\cdot10^{-15}$ |
| VV-2 acoustic | 6 | $1.1\cdot10^{-15}$ |
| VV-3 Couette | 3 | $9.4\cdot10^{-16}$ |
| VV-4 Becker | 6 | $2.2\cdot10^{-13}$ |
| VV-5 Ohmic, Alfvén, acoustic on MHD | 18 | $4.8\cdot10^{-15}$ |
| VV-6 Ohmic budget | 4 | $2.6\cdot10^{-15}$ |
| VV-7 conservation | 6 | $2.7\cdot10^{-15}$ |
| VV-7 forests | 8 | $3.5\cdot10^{-14}$ |
| VV-8 seams | 9 | $2.0\cdot10^{-14}$ |

*Measured: ADAM `0bcdf870` + the night-2026-11-22 tree (the paper commit, tag pending), CPU gfortran 16.2.0 + Open
MPI 5.0.11 (`iac_gnu`) on medooza (gnode22), FNL nvfortran 26.9, CUDA 12.9 (`iac_nvf`, the OpenACC build
`flume-fnl-nvf`) on NVIDIA A30, 2026-11-23 (`viscous/check.sh` on both executables, then `viscous/check.sh --legs
agree`, jobs 1369376 CPU, 1369378 FNL, 1369379 agree).*

The [DC contract](./input-contract) passes on both backends as well.

### The inviscid path

The same comparison for the inviscid runs of the order and MHD studies, on the last checkpoint of each run:
`make_doc_runs.sh` groups `vortex`, `order`, `mhd` and `shock-cylinder` run on each executable, then
`src/tests/regression_ab.py <cpu run> <fnl run> --ngc 3`. The measure differs from the table above: the largest
difference of any dataset (one variable of one block) divided by that dataset's own maximum, a stricter scale than the
physical scale of `waves.py agree`. Ghost cells are left out.

| Runs | Count | Largest difference, last checkpoint |
|---|---|---|
| isentropic vortex, splitting and `weno-riemann` LLF / HLL / HLLC, $N = 64, 128, 256$ | 12 | $2.1\cdot10^{-14}$ |
| MHD linear waves, splitting and HLLD, four families, $N = 16, 32, 64$ | 24 | $2.7\cdot10^{-14}$ (23 runs; see below) |
| field loop | 1 | $1.4\cdot10^{-14}$ |
| shock on a cylinder (IB), splitting / HLLC | 2 | $3.5\cdot10^{-12}$ / $1.9\cdot10^{-12}$ |
| Orszag–Tang uniform / `--amr` | 2 | $\psi$: $2.4\cdot10^{-7}$ / $1.5\cdot10^{-11}$ ($1.1\cdot10^{-8}$ / $8\cdot10^{-13}$ absolute) |
| rotor | 1 | $\psi$: $5.6\cdot10^{-8}$ ($5.9\cdot10^{-9}$ absolute) |

*Measured: as above, the runs of jobs 1369375 (CPU) and 1369377 (FNL), compared in job 1369379, 2026-11-23.*

- The smooth runs (vortex, waves, field loop) agree to about a hundred units of round-off ($2.2\cdot10^{-16}$): the
  two backends integrate the same discrete operator.
- The splitting Alfvén wave at $N = 16$ is the exception in the count: CPU takes 321 steps and FNL 320, a round-off
  difference in the last time step's length near `time_max`, so the final checkpoints are at different steps and are
  not compared (the initial ones agree to $1.4\cdot10^{-16}$).
- The shocked runs differ by $10^{-12}$–$10^{-11}$: round-off amplified where a WENO stencil or the IB switches.
- On Orszag–Tang and the rotor the largest difference is in the GLM scalar $\psi$, measured against its own block
  maximum, which is small ($\sim 0.05$ on Orszag–Tang): $10^{-8}$ in absolute value. `regression_ab.py` reports only
  the largest dataset of each file, so the difference of the conserved variables there is below that and was not
  printed separately. These are also the runs whose symmetry defect grows from round-off to $10^{-9}$ (see
  [Orszag–Tang](../mhd/orszag-tang) and [rotor](../mhd/rotor)).

Related cases: the Sod tube shows the same CPU–FNL round-off agreement on a Riemann problem ([Sod](../euler/sod)).

## Findings

- The largest differences sit in the Becker shock, the steepest gradients of the suite, at $2\cdot10^{-13}$ of the
  momentum scale: no run differs by more than about a thousand units of round-off.
