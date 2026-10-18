# Passive scalars

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| PS (ps1 to ps4) | Euler | M10 ([#92](https://github.com/szaghi/adam/issues/92)) | CPU, FNL | yes | `sod-ps` |

## Case and scope

`[physics] passive_scalars = N` carries up to four scalars $\rho s_k$ advected with the flow (composition, tracers,
metallicity), on `euler` and on `euler-de` ([numerics](../../numerics#passive-scalars)). A passive scalar must not act
on the flow, must be conserved, and must follow the flow exactly where its specific value is uniform.

What the case catches:
- **a scalar that acts on the flow:** the flow fields of a run with scalars must be those of the same run without them,
  bitwise;
- **a scalar flux inconsistent with the mass flux:** a uniform specific value must stay uniform, through every regrid
  and reflux;
- **a slot that leaks:** the unused scalar slots must stay empty;
- **a backend that disagrees, or a silent input.**

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Model | `euler` or `euler-de` with `passive_scalars = N`: four more conserved variables (nv 9 or 10), the first N used | `adam_flume_euler_ps_library`, `adam_flume_euler_de_ps_library` (the body `adam_flume_euler_library_body.INC`) |
| Space, time | as [Sod](./sod) (V1) and [runtime regrids](../runtime-amr/runtime-regrids): WENO5 characteristic flux splitting; each scalar a characteristic field of eigenvalue $u_n$ | `adam_flume_cpu_euler_ps_kernels.F90`, FNL `adam_flume_fnl_euler_ps_kernels.F90` |
| Initial values | `s1`..`sN` per region (the first region whose box holds the cell), and per inflow face | `set_initial_conditions`, [`adam_flume_ic_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_ic_object.F90) |

## Initial and boundary conditions

| Leg | Setup |
|---|---|
| ps1 | the `sod-x` input with two scalars, $s = 1$ left of the diaphragm and 0 right of it, on `euler` and on `euler-de`, each against the same run without scalars |
| ps2 | `passive_scalars = 5`; `passive_scalars = 1` with `mhd-ideal` |
| ps3 | the ps1 runs of CPU and FNL |
| ps4 | the `vortex-regrid` input (the isentropic vortex on a periodic quadtree, levels 1–3, a regrid every 5 steps, reflux on, 100 steps) with one scalar, $s_1 = 1$ everywhere |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/scalars/make_ps.py` | writes the inputs (`--case sod|vortex`, `--model`, `--scalars`, `--left`, `--right`) |
| `src/tests/flume/verification/scalars/check.sh` | runs the legs (`--legs ps1,ps2,ps3,ps4`, default all but ps3; `FLUME_EXE` selects the backend; ps3 compares the work directories of both) |
| `src/tests/flume/verification/scalars/ps_oracle.py` | the bitwise comparison and the drift bound (ps1), the CPU–FNL rule (ps3), the uniform scalar (ps4) |
| `src/tests/flume/regression/sod-ps/input.ini` | the frozen ps1 `euler` case |
| `make_doc_runs.sh <runs> scalars`, `make_doc_figures.py --only passive-scalars` | the figure (2 runs, 2 minutes) |

```bash
cd src/tests/flume/verification/scalars && ./check.sh                                  # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh                           # FNL
./check.sh --legs ps3                                                                 # CPU against FNL
```

## Results

![Passive scalars: the Sod contact and the uniform scalar through runtime AMR](/flume/passive-scalars.png)

*Left: Sod at $t = 0.2$, 50 cells. The scalar rides the contact (exact at $x = 0.6855$), smeared over a few cells like
the density's contact; WENO keeps it within $[-2.5\cdot10^{-4}, 1.0005]$, not strictly in $[0, 1]$. Right: the
regridded vortex with $s_1 = 1$ everywhere: the largest $|s_1 - 1|$ stays at $10^{-15}$ and the integral at round-off
through the six regrids.*

| Leg | Asserted | Bound | CPU | FNL |
|---|---|---|---|---|
| ps1 | the flow fields (and the entropy of `euler-de`) against the run without scalars | bitwise | bitwise equal | bitwise equal |
| ps1 | each scalar integral | $10^{-13}$ | $1.1\cdot10^{-16}$ | 0 |
| ps1 | the unused slots 3 and 4 | zero | zero | zero |
| ps2 | `passive_scalars = 5`; with `mhd-ideal` | refused | refused | refused |
| ps3 | per field, FNL against CPU, relative to the runs without scalars | $\le 2\times$ theirs, or $10^{-10}$ | the flow fields as the base ($\le 6.3\cdot10^{-14}$), the scalars $2.2\cdot10^{-14}$ | (the same comparison) |
| ps4 | the largest $\lvert s_1 - 1\rvert$; the integral drift | $10^{-12}$; $10^{-13}$ | $1.6\cdot10^{-15}$; $2.2\cdot10^{-16}$ | $1.0\cdot10^{-15}$; $1.1\cdot10^{-16}$ |

Related cases: [Sod](./sod) is the flow of ps1; [runtime regrids](../runtime-amr/runtime-regrids) the grid of ps4; the
entropy of `euler-de` uses the same machinery ([dual energy](./dual-energy)).

## Findings

- **A passive scalar is invisible to the flow.** Each scalar is a characteristic field of its own, and the Euler waves'
  left eigenvectors have no component on it, so the flow fields of a run with scalars are bitwise those without.
- **A uniform scalar stays uniform to round-off, through AMR.** Its characteristic field $\rho s - s\rho$ vanishes, so
  the scalar flux is $s$ times the mass flux whatever the reconstruction; the prolongation, the restriction and the
  reflux treat it like the density.
- **WENO bounds a scalar only to its own accuracy.** A mass fraction across a contact over- and undershoots by
  $5\cdot10^{-4}$; a model that needs $s \in [0, 1]$ strictly (chemistry) would need a bound-preserving limiter on the
  scalars.
- **The fixed capacity costs width.** The kernels carry four slots whatever $N$ (D-M10-1): a scalar run projects the
  wider state, the unused slots stay exactly zero.
