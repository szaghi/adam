# CPU and FNL agree run by run

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| VV `agree` | Dissipative terms (M4) | M4 ([#65](https://github.com/szaghi/adam/issues/65)) | CPU, FNL | no | — |

## Case and scope

Every run of the VV legs on the CPU against the same run on FNL, cell by cell on the last checkpoint.

What the case catches: the legs pass on both backends, but a pass bounds an error norm, and two backends can pass with
different fields. Every dissipative kernel exists twice, as a CPU loop and an FNL device kernel (OpenACC, or OpenMP
offload), so the same run on both must give the same state up to the round-off of the operation order. A device
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
| `src/tests/flume/verification/viscous/check.sh --leg agree` | runs nothing: pairs every `work-adam_flume_cpu-np2-*` directory with its `work-adam_flume_fnl-np2-*` twin, after every leg has run on both executables |
| `src/tests/flume/verification/viscous/waves.py agree` | the cell-by-cell comparison |

```bash
cd src/tests/flume/verification/viscous
./check.sh && FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh && ./check.sh --leg agree
```

No figure: a run-by-run identity has nothing a plot would show (#85).

## Results

Bound: $10^{-10}$ on every run.

| Leg | Runs | Largest difference |
|---|---|---|
| VV-1 shear | 12 | $2.7\cdot10^{-14}$ |
| VV-2 acoustic | 6 | $4.9\cdot10^{-15}$ |
| VV-3 Couette | 3 | $1.4\cdot10^{-13}$ |
| VV-4 Becker | 6 | $4.3\cdot10^{-13}$ |
| VV-5 Ohmic, Alfvén, acoustic on MHD | 18 | $4.7\cdot10^{-14}$ |
| VV-6 Ohmic budget | 4 | $5.3\cdot10^{-15}$ |
| VV-7 conservation | 6 | $9.3\cdot10^{-15}$ |
| VV-7 forests | 8 | $3.4\cdot10^{-14}$ |
| VV-8 seams | 9 | $8.4\cdot10^{-14}$ |

The [DC contract](./input-contract) passes on both backends as well.

Related cases: the Sod tube shows the same CPU–FNL round-off agreement on the inviscid path ([Sod](../euler/sod)).

## Findings

- The largest differences sit in the Becker shock, the steepest gradients of the suite, at $4\cdot10^{-13}$ of the
  momentum scale: no run differs by more than a few hundred units of round-off.
