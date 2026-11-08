# t0 elliptic solves on FLAIL

| Cases | Issues | Type | Backends | Figure | Regression |
|---|---|---|---|---|---|
| `rmf-pec-t0`, `rmf-open-t0`, `pic-pec-t0`, `pic-open-t0` | [#84](https://github.com/szaghi/adam/issues/84) (E3b, E3c; [#105](https://github.com/szaghi/adam/issues/105) open) | P | CPU, FNL | no | the four cases (digest + residuals, CPU and FNL) |

## Case and scope

At $t = 0$ PRISM builds the electric field of the initial charge and the magnetic field of the initial current by
elliptic solves on FLAIL, with PRISM's own discrete operators ([numerics](../numerics#t0-elliptic-solves)):

- electric: $D(G)\,\phi = -\rho/\varepsilon$, then $\mathbf D \leftarrow \mathbf D - \varepsilon\,G\phi$;
- magnetic: $D(G)\,A_d = -\mu J_d$ for each component, then $\mathbf B \leftarrow \mathbf B + \nabla_h\times\mathbf A$.

These four cases are the first that enter the t0 solves; before them none of the 18 regression cases did. Two face
kinds are covered:

- **PEC** (E3b): the potential is odd about the faces, $A_d$ even about the faces normal to $d$ and odd about the
  others;
- **exact/open**, PRISM's `Neumann` (E3c): every ghost layer holds the continuum free-space potential of the source,
  $-\sum f\,dV/(4\pi r)$ (`FLAIL_BC_FREE_GHOSTS`), summed on rank 0 over the gathered source.

What the cases catch:
- a t0 field that does not satisfy the constraint it was solved for (Gauss's law, $\nabla\cdot\mathbf B = 0$);
- a fall-back to the legacy smoother: the log must report every FLAIL solve;
- a t0 field that depends on the rank count: the legacy smoother summed the free-space potential over each rank's own
  blocks.

The PEC electric solve has one exact null mode. PRISM's centred divergence and gradient of order $p$ annihilate the
checkerboard $\chi = (-1)^{i+j+k}$, so with PEC faces $D(G)$ is singular on $\chi$ and the component of $\rho$ along
it is out of reach of any $\mathbf D$: $\nabla\cdot\mathbf D - \rho = -c\,\chi$, $c = \langle\rho,\chi\rangle/N$.
The oracle measures the residual after removing it.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | Maxwell, SI (`rmf-*`, 4 rectangular coils carrying their current at t0, `time_delay = 0`); PIC (`pic-*`, a plasma box) | [models](../models#models); [the fields at t0](../pic#the-fields-at-t0) |
| Space | centred finite differences of order 6, `ngc = 6` | [`adam_fdv_operators_library.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_fdv_operators_library.F90) |
| t0 solve | FLAIL; PEC → `FLAIL_BC_DIRICHLET` ($\phi$, tangential $A$) or `FLAIL_BC_NEUMANN` (normal $A$); exact/open → `FLAIL_BC_FREE_GHOSTS` | `impose_pic_fields_time_zero`, `solve_time_zero`, [`adam_prism_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/prism/common/adam_prism_common_object.F90); [`adam_flail_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_flail_object.F90); [hazards](../hazards#elliptic-boundary-map) |
| PIC deposit | TSC, non-conserving current (the conserving CPU deposit is single-block, single-rank) | [particle-in-cell](../pic#deposition-and-gather) |
| Time | SSP Runge–Kutta (5,4), CFL 0.5; the oracles read the step-0 checkpoint | — |

## Initial and boundary conditions

| Case | Domain (m) | Grid | Faces | Source | Steps |
|---|---|---|---|---|---|
| `rmf-pec-t0` | $[-0.06, 0.06]^3$ | 64 blocks of $16^3$ ($64^3$ cells) | `PEC` | 4 AC coils at t0 | 10 |
| `rmf-open-t0` | $[-0.06, 0.06]^3$ | 64 blocks of $8^3$ ($32^3$ cells) | `Neumann` | the same coils, reaching the faces | 10 |
| `pic-pec-t0` | $[-0.02, 0.02]^3$ | 8 blocks of $16^3$ ($32^3$ cells) | `PEC` | plasma box, density $2\cdot10^{9}$ | 4 |
| `pic-open-t0` | $[-0.02, 0.02]^3$ | 64 blocks of $8^3$ ($32^3$ cells) | `Neumann` | the same plasma | 4 |

Each `check.sh` runs its case at 1 and 2 ranks. The open cases use 64 blocks so that 2 ranks split them 32/32;
`pic-pec-t0` has 8 blocks, all on rank 0 at t0, so its rank-count leg checks the empty rank only.

## Input files

| File | Role |
|---|---|
| `src/tests/prism/regression/{rmf,pic}-{pec,open}-t0/input.ini` | the cases |
| `src/tests/prism/regression/{rmf,pic}-{pec,open}-t0/check.sh` | the runs at 1 and 2 ranks and the oracle calls (`--build`; `PRISM_EXE`, `PRISM_PY`) |
| `src/tests/prism/regression/t0_elliptic.py` | the oracles: `gauss`, `open`, `divb`, `same` |
| `src/tests/poisson/test_flail_bc.F90` (E3a, mode `test-flail-bc-gnu`) | PRISM's operators applied to FLAIL's solution reproduce the source to $10^{-15}$ |

```bash
cd src/tests/prism/regression/pic-pec-t0 && ./check.sh             # CPU
PRISM_EXE=$PWD/../../../../../exe/adam_prism_fnl ./check.sh        # FNL
```

## Results

Each case asserts, at 1 and 2 ranks: the run exits 0; the log reports the expected number of `FLAIL t0 solve` lines
(3 for the `rmf` cases, $A$ only; 4 for the `pic` cases, $\phi$ then $A$); the oracle below on the step-0 checkpoint;
then the step-0 fields of the two rank counts agree.

| Case | Oracle | Bound | CPU | FNL |
|---|---|---|---|---|
| `rmf-pec-t0` | `divb`: max $\lvert\nabla\cdot\mathbf B\rvert\,h/\max\lvert\mathbf B\rvert$ | $10^{-12}$ | ⟨run⟩ | ⟨run⟩ |
| `rmf-pec-t0` | `same`: $B$ at np 1 against np 2, relative | $10^{-12}$ | ⟨run⟩ | ⟨run⟩ |
| `rmf-open-t0` | `divb`, beyond 3 cells of the faces | $10^{-12}$ | ⟨run⟩ | ⟨run⟩ |
| `rmf-open-t0` | `same`: $B$ | $10^{-12}$ | ⟨run⟩ | ⟨run⟩ |
| `pic-pec-t0` | `gauss`: max $\lvert\nabla\cdot\mathbf D - \rho + c\chi\rvert/\max\lvert\rho\rvert$ | $10^{-12}$ | $\approx 10^{-16}$ | ⟨run⟩ |
| `pic-pec-t0` | checkerboard share $c/\max\lvert\rho\rvert$ (reported) | — | $3\cdot10^{-5}$ | ⟨run⟩ |
| `pic-pec-t0` | `same`: $D$, $B$ | $10^{-12}$ | ⟨run⟩ | ⟨run⟩ |
| `pic-open-t0` | `open`: max $\lvert\nabla\cdot\mathbf D - \rho\rvert/\max\lvert\rho\rvert$ beyond 3 cells of the faces | $10^{-11}$ (the CG tolerance) | $6.6\cdot10^{-13}$ | ⟨run⟩ |
| `pic-open-t0` | `same`: $D$, $B$ | $10^{-12}$ | ⟨run⟩ | ⟨run⟩ |

Against the legacy smoother FLAIL replaced (measured in #84):

| Case | Legacy smoother | FLAIL |
|---|---|---|
| `pic-pec-t0` | $3.4\cdot10^{-5}\max\lvert\rho\rvert$ of Gauss-law error beyond $\chi$ | $\approx 10^{-16}$ |
| `rmf-pec-t0` | stopped after 11 iterations with max $\lvert\mathbf B\rvert$ 2.5 times smaller; its $\mathbf B$ differs from FLAIL's by 1.6 (relative max norm) after 11 iterations, 0.48 after 40, 0.034 after 160 | reference |
| `rmf-open-t0` | $\mathbf B$ at 2 ranks differs from 1 rank by 3.5 (relative) | rank-independent |
| `pic-open-t0` | $\mathbf D$ and $\mathbf B$ differ between rank counts by $2\cdot10^{-2}$; Gauss law $4.6\cdot10^{-16}$ after 10000 sweeps | rank-independent |

Related pages: the [PIC cases](./pic) run the conserving deposit with the same t0 solves; [refusals](./refusals)
holds the refusal of a t0 solve on a refined grid.

## Findings

- With PEC faces Gauss's law holds to round-off up to the checkerboard mode, the one component of $\rho$ no
  $\mathbf D$ can carry with PRISM's centred operators ($3\cdot10^{-5}$ of max $\lvert\rho\rvert$ on `pic-pec-t0`).
- $\nabla\cdot\mathbf B = 0$ at t0 follows from $\nabla_h\cdot\nabla_h\times = 0$ whatever the solve's accuracy; the
  rank-count leg is the one that tests the solve.
- With open faces the operator is nonsingular, and Gauss's law holds to the CG tolerance on the cells that read no
  field-BC ghost. Nearer the faces the divergence reads ghosts that PRISM's own `Neumann` field BC rewrote, so the
  oracles skip 3 cells.
- The free-space ghosts are summed over the whole domain on rank 0: the t0 field no longer depends on the rank count.
