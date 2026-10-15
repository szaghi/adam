# Input contract of the dissipative terms

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| DC | Dissipative terms (M4) | M4 ([#65](https://github.com/szaghi/adam/issues/65)) | CPU, FNL | no | — |

## Case and scope

M4 reads each dissipative coefficient either as a coefficient or as its dimensionless number (#49 N4), refuses the
positivity limiter with any of them (D-M4-5) and adds the no-slip walls (D-M4-7). The case is a contract test, not a
flow: a list of 25 inputs, each a base input with some keys set, each with an expected outcome.

What the case catches:
- **an invalid combination that runs:** it must stop the run with a message that names it;
- **a wrong implied coefficient:** a valid combination must imply exactly the coefficient it states ($\mu = 1/Re$,
  $k = \mu c_p/Pr$, $\eta = 1/Rm$);
- **a key dropped silently:** a coefficient that is read but not used, or a dimensional key not converted back;
- **a side effect of the dissipative path on an ideal run:** zero coefficients must leave the run bit for bit equal to
  the base run.

The keys are listed in [input](../../input#dissipative-terms-issue-65-m4).

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Coefficients and numbers | `viscosity` or `reynolds`, `conductivity` or `prandtl`, `resistivity` or `magnetic_reynolds` or `lundquist`; `viscosity_law = power-law` with `viscosity_exponent`, `reference_temperature` | `load_from_file`, `laws`, [`adam_flume_dissipation_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_dissipation_object.F90) |
| Order of the dissipative fluxes | `[numerics] dissipative_order`, 2 or 4; 4 needs `ngc >= 3` | `load_from_file`, [`adam_flume_numerics_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_numerics_object.F90); `check_dissipation`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90) |
| Positivity limiter | refused with any dissipative coefficient (D-M4-5) | `check_positivity_limiter`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90) |
| No-slip walls | `wall-noslip` (tangential `wall_u/v/w`, a normal component is fatal), `wall-isothermal` (`wall_temperature > 0`) | `load_from_file`, [`adam_flume_bc_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_bc_object.F90); [boundary conditions](../../boundary-conditions#wall-noslip-and-wall-isothermal-no-slip-walls) |
| Diffusive time-step limit | logged at every run (`dissipative dt limit (issue #65)`) | `log_dt_limit`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90); [time-step bound](../../models#time-step-bound) |
| Dimensional input | `[reference]` section; `lundquist` needs the Alfvénic velocity preset | `classify`, [`adam_flume_reference_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_reference_object.F90); [units and scaling](../../models#units-and-scaling) |

## Initial and boundary conditions

The cases modify two regression inputs; the flow itself is not judged.

| Item | Value |
|---|---|
| Euler base | `regression/sod-x`: the [Sod](../euler/sod) tube along $x$, 64 blocks of $50 \times 4 \times 4$ cells, `extrapolation` on every face |
| MHD base | `regression/uniform-amr-mhd`: a [uniform](../../initial-conditions#uniform-seeded-uniform-state) MHD state on a periodic octree, $8^3$ cells per block, level 1 with the box $[0, 0.5]^3$ refined to level 2, GLM, 10 steps |
| Ranks | 2 |
| Steps | the `log` and `convert` cases stop after 10 steps (`it_max = 10`); the `ideal` case runs the base input to its end |
| Dimensionalisation (`convert`) | `scaling.py dimensionalize --j 1 --k -1 --m 1`: $L_0 = 2$, $u_0 = 1/2$, $\rho_0 = 4$ |

The 25 cases (`make_contract.py --list`):

| Kind | Cases | Oracle |
|---|---|---|
| refuse (18) | both keys of a term (3: `viscosity` + `reynolds`, `conductivity` + `prandtl`, `magnetic_reynolds` + `lundquist`), `prandtl` without a viscosity, resistivity on Euler, a negative viscosity, a zero Reynolds number, the power law without a viscosity, its keys without the law, a missing reference temperature, an unknown law, `dissipative_order = 3`, the limiter with a viscosity, a normal wall velocity, an isothermal wall without (or with a negative) temperature, `dissipative_order = 4` with `ngc = 2`, a negative power-law exponent | the run stops with the expected message |
| log (3) | `reynolds` + `prandtl` (64, 0.75); `viscosity` + `conductivity` (0.015625, 0.03125); `magnetic_reynolds` + `reynolds` (128, 32, MHD) | the logged $\mu$, $k$, $\eta$ equal $1/Re$, $\mu c_p/Pr$, $1/Rm$ exactly; every run takes ten steps and logs its diffusive time-step limit |
| ideal (1) | zero coefficients, `dissipative_order = 2` | the run equals the base run bit for bit |
| convert (2) | the Euler case with viscosity, conductivity, the power law and an isothermal moving wall; the MHD case with resistivity; dimensionalised with $L_0 = 2$, $u_0 = 1/2$, $\rho_0 = 4$ | every dimensional key logged converted back exactly (`scaling.py check-log`), the conductivity and the temperatures with the gas constant |
| convert-refuse (1) | `lundquist` under a `[reference]` without the Alfvénic preset | refused |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/dissipation/make_contract.py` | the 25 cases (`--list`), writes each input, judges the logged coefficients (`--check-log`) |
| `src/tests/flume/verification/dissipation/check.sh` | runs every case and applies its oracle (`--np N`; `FLUME_EXE` selects the backend) |
| `src/tests/flume/verification/scaling/scaling.py` | `dimensionalize` and `check-log` of the `convert` cases |
| `src/tests/flume/regression/sod-x/input.ini`, `src/tests/flume/regression/uniform-amr-mhd/input.ini` | the base inputs |

```bash
cd src/tests/flume/verification/dissipation && ./check.sh             # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh           # FNL
```

No figure: a contract leg has nothing a plot would show (#85).

## Results

| Check | Bound | CPU | FNL |
|---|---|---|---|
| 18 refusals | expected message | pass | pass |
| 3 logged coefficients | exact | pass | pass |
| ideal run against the base run | bitwise | pass | pass |
| 2 conversions | exact | pass (34 keys on the Euler case, 30 on the MHD one) | pass |
| 1 conversion refusal | expected message | pass | pass |

All 25 cases pass, on the CPU (2 ranks) and on FNL.

Related cases: the coefficients are exercised by [viscous waves](./viscous-waves) and [Ohmic decay](./ohmic-alfven);
the dimensional input by [scaling covariance](../infrastructure/scaling-covariance).

## Findings

- Every invalid combination of D-M4-5, D-M4-6 and D-M4-7 stops the run with a message that names the keys; none runs
  with a coefficient dropped.
- The coefficient a number implies is logged bit-exactly, so the conversion $1/Re$, $\mu c_p/Pr$, $1/Rm$ has no hidden
  rounding step.
- The dissipative path is inert at zero coefficients: the run is bitwise the base run.
