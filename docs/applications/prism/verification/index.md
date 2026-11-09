# Verification

PRISM's regression suite (`src/tests/prism/regression/`) pins every case to a golden digest
([regression reference](/reference/regression#prism-suite)). The pages below document the cases that also carry an
**oracle**: a `check.sh` that asserts a physical or structural property of the run (an invariant, an agreement
between two discretisations of the same problem, a refusal) instead of a stored number.

Each page has the same sections: case and scope, models and numerics with references to the code, initial and
boundary conditions, input files, results, findings.

The tag on each card says what its page verifies:

- **physics**: a run against an exact or reference solution, an order of accuracy or a conservation law;
- **contract**: a bitwise identity, rank or decomposition invariance, a structural invariant, an input refusal.

A page without a figure reports its results as a table of bounds and measured values.

## Fields on uniform and refined grids

The source-free Maxwell system and its discrete invariants on a uniform grid, then the 2:1 faces of static AMR on octrees and quadtrees: the reflux on the finite-volume path and the divergence source the faces inject.

<div class="adam-gallery adam-catalog">

[<strong>Uniform grid: energy and divergence</strong><span>a Gaussian pulse once around a periodic box: the energy kept up to the time integrator, $\nabla_h\cdot\mathbf D$ and $\nabla_h\cdot\mathbf B$ exactly</span><em>physics</em>](./uniform-grid)

[<strong>Static AMR and the FV reflux register</strong><span>the 2:1 coarse–fine face on the finite-volume path: each stage of the Berger–Colella reflux present and wired to the tree, against a control without the face</span><em>contract</em>](./amr-reflux)

[<strong>2:1 seam div(B): truncation order</strong><span>the $\nabla\cdot\mathbf B$ source a 2:1 face injects: truncation order, convergent under refinement, growing in time at fixed $h$</span><em>physics</em>](./seam-divergence)

[![Quadtree 2:1 seams](/prism/quadtree-pulse.png)<strong>Quadtree 2:1 seams</strong><span>2:1 seams of a quadtree, refined in $x$ and $y$ only: the pulse against the octree, without a spurious $z$ dependence</span><em>contract</em>](./quadtree)

[<strong>FV reflux on a capped last step</strong><span>the reflux scaled by the capped last step of a time-driven run: the domain integrals conserved on the last step too</span><em>physics</em>](./reflux-capped-step)

</div>

## Forests of realms

Realms glued by inter-realm seams run as one domain: 1:1 seams, seams between misaligned blocks and a 2:1 seam, each against the single-realm run.

<div class="adam-gallery adam-catalog">

[<strong>Inter-realm 1:1 seams: α and β cadence</strong><span>two realms glued by a 1:1 mirror seam against the single realm: the divergences held with the stage-coincident cadence</span><em>contract</em>](./inter-realm-mirror)

[<strong>Inter-realm seams between misaligned blocks</strong><span>1:1 seams between blocks that do not line up: every register face closes on the finite-volume path</span><em>contract</em>](./inter-realm-misaligned)

[<strong>Inter-realm 2:1 seam</strong><span>a 2:1 seam between two realms through the intra-realm 2:1 formulas, against the single-realm 2:1 run</span><em>contract</em>](./inter-realm-refined)

</div>

## Absorbing layers

The far-wave layer and the perfectly matched layer, on uniform and refined grids.

<div class="adam-gallery adam-catalog">

[<strong>Absorbing layers on AMR</strong><span>the far-wave layer and the perfectly matched layer on the six faces, on a uniform grid and after refinement, pinned by their digests</span><em>contract</em>](./absorbing-layers-amr)

</div>

## Elliptic solves and particles

The $t = 0$ fields of the initial charge and current by FLAIL solves, and the particle-in-cell path of the production runs.

<div class="adam-gallery adam-catalog">

[<strong>t0 elliptic solves on FLAIL</strong><span>the fields of the initial charge and current at $t = 0$ by elliptic solves with PRISM's own operators, PEC and open boundaries</span><em>physics</em>](./t0-elliptic)

[<strong>PIC: periodic faces and the Fresnel slab</strong><span>the production PIC path on a periodic plasma and a Fresnel slab: the periodic current fold, charge conservation, the slab against the exact cold slab</span><em>physics</em>](./pic)

</div>

## Input and output

The inputs PRISM refuses, and the slices it writes.

<div class="adam-gallery adam-catalog">

[<strong>Input validation and refusals</strong><span>inputs PRISM cannot run stop at initialization with a message that names the option and the value</span><em>contract</em>](./refusals)

[<strong>Slices output</strong><span>every value of the slice files recomputed from the checkpoint: decomposition-independent, scaled in non-dimensional runs, accurate on an analytic wave</span><em>contract</em>](./slices)

</div>

## How it was built

The pages in the order the capabilities were implemented, with the issues that built them.

| Built under | Page | Cases / legs |
|---|---|---|
| — (base) | [Uniform grid: energy and divergence](./uniform-grid) | `periodic-pulse-blocks` (digest) |
| [#13](https://github.com/szaghi/adam/issues/13) (M1–M3), [#23](https://github.com/szaghi/adam/issues/23), [#24](https://github.com/szaghi/adam/issues/24), [#28](https://github.com/szaghi/adam/issues/28) | [Static AMR and the FV reflux register](./amr-reflux) | `rmf-amr/check.sh`: structural, control, register parity |
| [#21](https://github.com/szaghi/adam/issues/21), [#25](https://github.com/szaghi/adam/issues/25), [#26](https://github.com/szaghi/adam/issues/26), [#29](https://github.com/szaghi/adam/issues/29) | [2:1 seam div(B): truncation order](./seam-divergence) | `rmf-amr-fd/check.sh`; `rmf-amr-fd-pulse/check.sh` (control, uniform, seam, RK legs, `--convergence`) |
| [#46](https://github.com/szaghi/adam/issues/46) | [Quadtree 2:1 seams](./quadtree) | `rmf-amr-fd-pulse/check.sh`, quadtree leg |
| [#16](https://github.com/szaghi/adam/issues/16), [#18](https://github.com/szaghi/adam/issues/18), [#31](https://github.com/szaghi/adam/issues/31) | [Inter-realm 1:1 seams: α and β cadence](./inter-realm-mirror) | `rmf-2realm`, `rmf-2realm-asymK`, `rmf-2realm-stagesync` (digests); `rmf-2realm-fd-pulse/check.sh` |
| [#51](https://github.com/szaghi/adam/issues/51), [#54](https://github.com/szaghi/adam/issues/54) | [Inter-realm seams between misaligned blocks](./inter-realm-misaligned) | `rmf-2realm-fv-pulse-blocks/check.sh`: `small1`, `small2`, `one1`, `one2` |
| [#52](https://github.com/szaghi/adam/issues/52) | [Inter-realm 2:1 seam](./inter-realm-refined) | `rmf-2realm-fd-pulse-refined/check.sh` |
| [#38](https://github.com/szaghi/adam/issues/38) | [FV reflux on a capped last step](./reflux-capped-step) | `rmf-amr-fv-pulse-timed/check.sh` |
| [#31](https://github.com/szaghi/adam/issues/31), [#58](https://github.com/szaghi/adam/issues/58), [#95](https://github.com/szaghi/adam/issues/95) (P5.1) | [Absorbing layers on AMR](./absorbing-layers-amr) | `rmf-fwl`, `rmf-fwl-amr`, `rmf-pml-amr` (digests) |
| [#11](https://github.com/szaghi/adam/issues/11), [#55](https://github.com/szaghi/adam/issues/55), [#56](https://github.com/szaghi/adam/issues/56), [#57](https://github.com/szaghi/adam/issues/57), [#61](https://github.com/szaghi/adam/issues/61), [#84](https://github.com/szaghi/adam/issues/84), [#111](https://github.com/szaghi/adam/issues/111) | [Input validation and refusals](./refusals) | `refusals.sh` (24 legs); `rmf-hyperbolic` (digest, `nv = 11`) |
| [#84](https://github.com/szaghi/adam/issues/84) (E3b, E3c) | [t0 elliptic solves on FLAIL](./t0-elliptic) | `rmf-pec-t0`, `rmf-open-t0`, `pic-pec-t0`, `pic-open-t0` (`check.sh` each) |
| [#111](https://github.com/szaghi/adam/issues/111), [#112](https://github.com/szaghi/adam/issues/112) | [PIC: periodic faces and the Fresnel slab](./pic) | `pic-periodic`, `pic-slab-pulse` (`pic_check.sh` legs) |
| — | [Slices output](./slices) | `slices/check.sh` (oracle, decomposition, adimensional, `--convergence`, input) |

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
