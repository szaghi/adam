# Verification

PRISM's regression suite (`src/tests/prism/regression/`) pins every case to a golden digest
([regression reference](/reference/regression#prism-suite)). The pages below document the cases that also carry an
**oracle**: a `check.sh` that asserts a physical or structural property of the run (an invariant, an agreement
between two discretisations of the same problem, a refusal) instead of a stored number.

The pages are listed in the order the capabilities were implemented, so the list reads as the path the solver took.
Each page has the same sections: case and scope, models and numerics with references to the code, initial and
boundary conditions, input files, results, findings.

**Type** marks what a page verifies:

- **P**, physics: a run against an exact or reference solution, an order of accuracy or a conservation law;
- **C**, contract: a bitwise identity, rank or decomposition invariance, a structural invariant, an input refusal.

**Figure** "no" marks a page whose results are a table of bounds and measured values only.

## Timeline

| # | Page | Issues | Cases / legs | Type | Figure |
|---|---|---|---|---|---|
| 1 | [Uniform grid: energy and divergence](./uniform-grid) | — (base) | `periodic-pulse-blocks` (digest) | P | no |
| 2 | [Static AMR and the FV reflux register](./amr-reflux) | [#13](https://github.com/szaghi/adam/issues/13) (M1–M3), [#23](https://github.com/szaghi/adam/issues/23), [#24](https://github.com/szaghi/adam/issues/24), [#28](https://github.com/szaghi/adam/issues/28) | `rmf-amr/check.sh`: structural, control, register parity | C | no |
| 3 | [2:1 seam div(B): truncation order](./seam-divergence) | [#21](https://github.com/szaghi/adam/issues/21), [#25](https://github.com/szaghi/adam/issues/25), [#26](https://github.com/szaghi/adam/issues/26), [#29](https://github.com/szaghi/adam/issues/29) | `rmf-amr-fd/check.sh`; `rmf-amr-fd-pulse/check.sh` (control, uniform, seam, RK legs, `--convergence`) | P | no |
| 4 | [Quadtree 2:1 seams](./quadtree) | [#46](https://github.com/szaghi/adam/issues/46) | `rmf-amr-fd-pulse/check.sh`, quadtree leg | C | yes |
| 5 | [Inter-realm 1:1 seams: α and β cadence](./inter-realm-mirror) | [#16](https://github.com/szaghi/adam/issues/16), [#18](https://github.com/szaghi/adam/issues/18), [#31](https://github.com/szaghi/adam/issues/31) | `rmf-2realm`, `rmf-2realm-asymK`, `rmf-2realm-stagesync` (digests); `rmf-2realm-fd-pulse/check.sh` | C | no |
| 6 | [Inter-realm seams between misaligned blocks](./inter-realm-misaligned) | [#51](https://github.com/szaghi/adam/issues/51), [#54](https://github.com/szaghi/adam/issues/54) | `rmf-2realm-fv-pulse-blocks/check.sh`: `small1`, `small2`, `one1`, `one2` | C | no |
| 7 | [Inter-realm 2:1 seam](./inter-realm-refined) | [#52](https://github.com/szaghi/adam/issues/52) | `rmf-2realm-fd-pulse-refined/check.sh` | C | no |
| 8 | [FV reflux on a capped last step](./reflux-capped-step) | [#38](https://github.com/szaghi/adam/issues/38) | `rmf-amr-fv-pulse-timed/check.sh` | P | no |
| 9 | [Absorbing layers on AMR](./absorbing-layers-amr) | [#31](https://github.com/szaghi/adam/issues/31), [#58](https://github.com/szaghi/adam/issues/58), [#95](https://github.com/szaghi/adam/issues/95) (P5.1) | `rmf-fwl`, `rmf-fwl-amr`, `rmf-pml-amr` (digests) | C | no |
| 10 | [Input validation and refusals](./refusals) | [#11](https://github.com/szaghi/adam/issues/11), [#55](https://github.com/szaghi/adam/issues/55), [#56](https://github.com/szaghi/adam/issues/56), [#57](https://github.com/szaghi/adam/issues/57), [#61](https://github.com/szaghi/adam/issues/61), [#84](https://github.com/szaghi/adam/issues/84), [#111](https://github.com/szaghi/adam/issues/111) | `refusals.sh` (24 legs); `rmf-hyperbolic` (digest, `nv = 11`) | C | no |
| 11 | [t0 elliptic solves on FLAIL](./t0-elliptic) | [#84](https://github.com/szaghi/adam/issues/84) (E3b, E3c) | `rmf-pec-t0`, `rmf-open-t0`, `pic-pec-t0`, `pic-open-t0` (`check.sh` each) | P | no |
| 12 | [PIC: periodic faces and the Fresnel slab](./pic) | [#111](https://github.com/szaghi/adam/issues/111), [#112](https://github.com/szaghi/adam/issues/112) | `pic-periodic`, `pic-slab-pulse` (`pic_check.sh` legs) | P | no |
| 13 | [Slices output](./slices) | — | `slices/check.sh` (oracle, decomposition, adimensional, `--convergence`, input) | C | no |

## Running a case

Every `check.sh` runs the CPU executable `exe/adam_prism_cpu` by default; `PRISM_EXE` selects another one. The PIC
cases take the backend as their first argument. `mpirun` must be on `PATH`, and the scripts read the Python venv
`exe/.regression-venv` (h5py, numpy) that `run.sh` creates on its first run.

```bash
cd src/tests/prism/regression/rmf-amr-fd-pulse && ./check.sh       # CPU backend
PRISM_EXE=$PWD/../../../../../exe/adam_prism_fnl ./check.sh        # FNL backend (GPU)
./check.sh --convergence                                           # adds the opt-in refinement leg
cd ../pic-periodic && ./check.sh cpu                               # PIC cases: the backend is an argument
cd .. && PRISM_EXE=$PWD/../../../../exe/adam_prism_fnl ./refusals.sh   # the init refusals, FNL
```

| Script | Ranks | Options |
|---|---|---|
| `rmf-amr/check.sh` | 1 and 2 | `--build` |
| `rmf-amr-fd/check.sh` | 1 | `--build` |
| `rmf-amr-fd-pulse/check.sh` | 1 | `--build`, `--convergence`, `--keep` (keeps the quadtree checkpoints for the figure) |
| `rmf-2realm-fd-pulse/check.sh` | 1 | `--build` |
| `rmf-2realm-fv-pulse-blocks/check.sh`, `rmf-2realm-fd-pulse-refined/check.sh` | 2 | `--build`, `--np N` |
| `rmf-amr-fv-pulse-timed/check.sh` | 2 | `--np N` |
| `{rmf,pic}-{pec,open}-t0/check.sh` | 1 and 2 | `--build`; `PRISM_PY` selects the Python |
| `pic-{periodic,slab-pulse}/check.sh <cpu\|fnl>` | CPU 1; FNL 1 and 2 | `--capture`; `FNL_NP`, `PIC_CHECK_LAUNCH`, `PRISM_PY` |
| `slices/check.sh` | 1, 2, 4 | `--build`, `--np "N ..."`, `--convergence` |
| `refusals.sh` | 1 | none |

`--build` always builds `prism-cpu-gnu`, never the `PRISM_EXE` override. The FNL runs need the NVHPC toolchain on
`PATH`; on the WSL box they also need the UCX settings of `run-fnl-local.sh`
([WSL caveats](/cookbook/gpu#gpu5-run-on-a-development-box-with-wsl2-not-on-a-cluster)). The goldened cases are checked by `run.sh`
([running the suite](/reference/regression#running)); the oracles are not part of that sweep and are run by hand.

`run_verification_legs.sh` runs every oracle of this section and logs each case to `<outdir>/<case>.log`:

```bash
src/tests/prism/regression/run_verification_legs.sh $PWD $PWD/exe/adam_prism_cpu /tmp/prism-verif          # CPU
src/tests/prism/regression/run_verification_legs.sh --fnl $PWD $PWD/exe/adam_prism_fnl /tmp/prism-verif    # FNL
```

The values in the result tables of this section come from these logs (2026-11-08, code at `e3da0473`): the CPU legs on
a medooza node (gfortran 16.2.0, Open MPI 5.0.11), the FNL legs on A30 GPUs (nvfortran 26.9, CUDA 12.9), and
`slices` on the WSL2 workstation (the same compilers, its two GPUs for FNL).

## Not done

Analytic oracles for the physics are tracked and not yet implemented:

- [#106](https://github.com/szaghi/adam/issues/106): plane-wave and Gaussian-pulse convergence, PEC image, cavity
  mode, PML reflection against the angle of incidence, the infinite wire against Biot–Savart, the Larmor orbit, the
  plasma frequency;
- [#105](https://github.com/szaghi/adam/issues/105): the Fresnel slab against its exact spectral solution. The
  plane-mean agreement with the exact cold slab is already measured on `pic-slab-pulse` ([PIC page](./pic)).
