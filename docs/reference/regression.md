# Regression suites

ADAM is tested in two layers. The **regression suites** run on every change and detect any change of a committed
result, wanted or not. The **verification legs** assert the physics (an exact solution, a reference, a symmetry, a
conservation property) and run for the part a change touches. This page is the reference for the suites; the tasks
(run a suite, add a case, refresh a golden, bring up a new machine) are in the [regression cookbook](/cookbook/regression).

| Suite | Directory | Cases | Backends | Built modes |
|---|---|---|---|---|
| FLUME | `src/tests/flume/regression/` | 55 | CPU (CI), FNL (workstation) | `flume-cpu-gnu` + `flume-cpu-gnu-fft`, `flume-fnl-nvf` + `flume-fnl-nvf-fft` |
| PRISM | `src/tests/prism/regression/` | 24 | CPU (CI), FNL (workstation) | `prism-cpu-gnu`, `prism-fnl-nvf` |

## The harness

Both suites share one harness pattern, `run.sh {cpu|fnl} [--no-build] [--varset <name>]`:

1. It runs the lints: the GPU race-shape check (`src/tests/lint/check-gpu-race-shapes.sh`), then the removed-keys check
   (`src/tests/lint/check-removed-input-keys.sh`).
2. It creates a Python venv at `exe/.regression-venv/` on the first run and installs `h5py` and `numpy` for
   `digest.py`. The venv lives under `exe/`, not `src/`, because `fobis` and the API-doc generator scan `src/`.
3. It builds the backend (skipped with `--no-build`).
4. It discovers the cases: every immediate subdirectory with an `input.ini`.
5. It runs each case under `mpirun -np 2` in `<case>/work-<backend>/` and compares the outputs with `<case>/golden/<backend>/`.

| Control | Effect |
|---|---|
| `--no-build` | reuse the executables in `exe/`; FLUME refuses an executable older than any of its sources |
| `--varset <name>` | build with that varset ([install](/start/install#varsets-one-per-machine-and-toolchain)) |
| `REGRESSION_CASES="a b"` | run only the named cases; in PRISM a name that is not a case directory stops the harness with exit 2 before the build |
| `REGRESSION_RUN_GOLDENLESS=1` | run the cases that have no golden for the backend (golden capture) |
| `REGRESSION_ALLOW_MPIRUN_FAILURE=1` | turn the failure of a run that exits non-zero into a warning, on an environment whose crash is known |

A case **fails** when its outputs differ from its goldens, and also when its run exits non-zero (a crash or the hard
timeout), even if the outputs match: a crash after the last checkpoint is still a crash. The outputs are compared
anyway, for the report. The one exemption is the harness's own termination of a run that hung in `MPI_Finalize` after
its end-of-run marker. A case **skips** when it has no golden for the backend. The summary is `PASS: n FAIL: m SKIP: k`,
and the exit status is non-zero on any failure.

Every case runs on 2 ranks, so a golden is a two-rank golden: a one-rank rerun does not match it.

### Goldens and tolerance

The raw HDF5 checkpoints (about 100 MB per rank for PRISM) are not committed. Each case commits compact references,
compared tolerance-aware:

| Golden | Content | Comparison |
|---|---|---|
| `digest.txt` | for every checkpoint and every field variable, `count, min, max, sum, sum_sq` over all blocks, ranks and (for forests) realms (`src/tests/prism/regression/digest.py`, shared by both suites) | `digest.py compare`, `rtol = 1e-6`, `atol = 1e-3` |
| `*-residuals.dat` | per-iteration residuals | header verbatim, integer columns exact, float columns within the same tolerances |
| `*-conservation_history.dat` (FLUME) | per-iteration volume integrals of the conservative variables (and $W$, $E + W$ under self-gravity) | as the residuals |

`rtol = 1e-6` absorbs FMA, libm and MPI-reduction-order differences between compilers and machines (about 9
significant digits are still checked); `atol = 1e-3` is the floor for antisymmetric fields whose `sum` is pure
cancellation noise. The `div_*` diagnostics are excluded from the digest (constraint residues at round-off). Per-block
metadata that differs only by block partitioning (`dxdydz`, `origin`, `time_iteration`) is downgraded to
`SKIP_METADATA`, which lets a multi-realm run be compared with a single-realm golden.

Bitwise identity is checked separately where it is the contract: FLUME's OpenMP build against the serial one
(`run-omp-bitwise.sh`, in CI), and RIPPLE's rank and CRLF legs.

**Golden authority.** The CPU goldens are captured in CI (Ubuntu 24.04, gfortran 14, varset `ci_gnu`); a golden captured
on a workstation with another compiler can fail in CI with no real change. The FNL goldens are captured on the GPU
workstation; the CI runners have no GPU.

## Running

```bash
src/tests/flume/regression/run.sh cpu                      # builds flume-cpu-gnu and flume-cpu-gnu-fft; the CI job
src/tests/flume/regression/run.sh fnl --varset local_nvf   # builds flume-fnl-nvf and flume-fnl-nvf-fft
src/tests/prism/regression/run.sh cpu                      # builds prism-cpu-gnu; the CI job
src/tests/prism/regression/run.sh fnl --varset local_nvf   # builds prism-fnl-nvf
REGRESSION_CASES="sod-x brio-wu" src/tests/flume/regression/run.sh cpu --no-build
```

On the WSL2 workstation the FNL suites run through `run-fnl-local.sh [--no-build]` in each suite directory: it loads
the `nvhpc` module, exports the WSL-only MPI settings of [GPU5](/cookbook/gpu#gpu5-run-on-a-development-box-with-wsl2-not-on-a-cluster)
and calls `run.sh fnl --varset local_nvf`. It needs Lmod; on any other machine call `run.sh` directly.

On a one-GPU box the two ranks share the device. Each FNL rank counts the node-local ranks bound to its device
(through `CUDA_VISIBLE_DEVICES`, `ROCR_VISIBLE_DEVICES`, `HIP_VISIBLE_DEVICES` and the device it picks) and sizes its
block pool for its share; the startup log prints `device shared by N rank(s)`.

## FLUME suite

Each case directory `src/tests/flume/regression/<case>/` holds a frozen `input.ini` (a multi-realm case: a forest
manifest plus one INI per realm), most of them copies of a [verification](/applications/flume/verification/) input cut to
a few dozen steps, and the goldens `golden/cpu/` and `golden/fnl/`.

| Script | Role |
|---|---|
| `run.sh cpu` | builds `flume-cpu-gnu` and `flume-cpu-gnu-fft`, runs every case, compares; runs in CI |
| `run-fnl-local.sh` | the same on the GPU (`flume-fnl-nvf` and `flume-fnl-nvf-fft`, `--varset local_nvf`); a workstation gate |
| `run-omp-bitwise.sh` | checks that the OpenMP build reproduces the serial one bit for bit (CI) |

Two executables per backend: every case with `self_gravity = .true.` runs on the `-fft` executable (FFTW on the CPU,
cuFFT on FNL), since every self-gravity solve needs an FFT; the CI job installs `libfftw3-dev`.

### FLUME cases

| Case | Freezes | Issue |
|---|---|---|
| `sod-x`, `sod-y`, `sod-z` | [Sod shock tube](/applications/flume/verification/euler/sod), along each axis | [#35](https://github.com/szaghi/adam/issues/35) |
| `sod-hllc` | [Sod](/applications/flume/verification/euler/sod) on `weno-riemann`, HLLC | [#47](https://github.com/szaghi/adam/issues/47) |
| `vortex-periodic` | [isentropic vortex](/applications/flume/verification/euler/isentropic-vortex), 100 steps | [#35](https://github.com/szaghi/adam/issues/35) |
| `vortex-rk3` | `vortex-periodic` on the low-storage `runge-kutta-3` at CFL 0.8, final state digested: the only case on the low-storage stage update (every other case runs an SSP scheme, and a low-storage scheme is refused on AMR seams and multi-realm runs). At CFL 0.1 the two schemes agree to $10^{-7}$, inside the digest tolerances; at 0.8 the digest rejects the SSP-54 run | [#95](https://github.com/szaghi/adam/issues/95) |
| `amr-periodic-reflux` | [conservation across AMR seams](/applications/flume/verification/euler/amr-conservation): init-time AMR, reflux | [#35](https://github.com/szaghi/adam/issues/35) |
| `amr-periodic-reflux-stage` | the same with `reflux_cadence = stage` | [#81](https://github.com/szaghi/adam/issues/81) |
| `shock-cylinder-ib` | [shock over a cylinder](/applications/flume/verification/euler/shock-cylinder), immersed boundary, 120 steps | [#35](https://github.com/szaghi/adam/issues/35) |
| `sod-2realm`, `sod-2realm-blocks`, `sod-amr`, `sod-amr-2realm`, `sod-amr-refined` | [multi-realm runs](/applications/flume/verification/euler/multi-realm): mirror seams lined up or not, a 2:1 face, the refined coupling | [#37](https://github.com/szaghi/adam/issues/37), [#51](https://github.com/szaghi/adam/issues/51), [#52](https://github.com/szaghi/adam/issues/52) |
| `rj2a-x`, `rj2a-y`, `rj2a-z` | [Ryu–Jones 2a](/applications/flume/verification/mhd/ryu-jones-2a) | [#41](https://github.com/szaghi/adam/issues/41) |
| `rj2a-hlld` | Ryu–Jones 2a on `weno-riemann`, HLLD | [#47](https://github.com/szaghi/adam/issues/47) |
| `rj2a-2realm` | Ryu–Jones 2a in two realms | [#41](https://github.com/szaghi/adam/issues/41) |
| `brio-wu` | [Brio–Wu](/applications/flume/verification/mhd/brio-wu-rj4d) | [#41](https://github.com/szaghi/adam/issues/41) |
| `glm-pulse` | [GLM cleaning](/applications/flume/verification/mhd/glm-cleaning) | [#41](https://github.com/szaghi/adam/issues/41) |
| `orszag-tang` | [Orszag–Tang](/applications/flume/verification/mhd/orszag-tang) | [#41](https://github.com/szaghi/adam/issues/41) |
| `orszag-tang-eglm` | Orszag–Tang with [EGLM](/applications/flume/verification/mhd/eglm) | [#47](https://github.com/szaghi/adam/issues/47) |
| `rotor` | [MHD rotor](/applications/flume/verification/mhd/rotor) | [#41](https://github.com/szaghi/adam/issues/41) |
| `field-loop` | [field loop](/applications/flume/verification/mhd/field-loop) | [#41](https://github.com/szaghi/adam/issues/41) |
| `rotated-shock-tube` | [rotated shock tube](/applications/flume/verification/mhd/rotated-shock-tube) | [#41](https://github.com/szaghi/adam/issues/41) |
| `uniform-amr-mhd` | MHD on a uniform state across 2:1 faces | [#41](https://github.com/szaghi/adam/issues/41) |
| `blast-limiter`, `blast3d-limiter` | the Balsara–Spicer blast with the [positivity limiter](/applications/flume/verification/mhd/positivity-limiter), 2-D and 3-D | [#47](https://github.com/szaghi/adam/issues/47) |
| `blast-amr-limiter` | the limited blast across two 2:1 seams | [#50](https://github.com/szaghi/adam/issues/50) |
| `shear-viscous` | [viscous waves](/applications/flume/verification/dissipative/viscous-waves) | [#65](https://github.com/szaghi/adam/issues/65) |
| `couette` | [Couette flow](/applications/flume/verification/dissipative/couette) | [#65](https://github.com/szaghi/adam/issues/65) |
| `becker-shock` | [Becker's viscous shock](/applications/flume/verification/dissipative/becker-shock) | [#65](https://github.com/szaghi/adam/issues/65) |
| `alfven-resistive` | [the visco-resistive Alfvén wave](/applications/flume/verification/dissipative/ohmic-alfven) | [#65](https://github.com/szaghi/adam/issues/65) |
| `seam-viscous-amr` | [the dissipative fluxes across 2:1 seams](/applications/flume/verification/dissipative/seam-conservation) | [#65](https://github.com/szaghi/adam/issues/65) |
| `vortex-regrid` | [runtime regrids](/applications/flume/verification/runtime-amr/runtime-regrids): the vortex, 6 regrids | [#74](https://github.com/szaghi/adam/issues/74) |
| `field-loop-regrid` | [runtime regrids](/applications/flume/verification/runtime-amr/runtime-regrids) on MHD: the field loop, 10 regrids | [#74](https://github.com/szaghi/adam/issues/74) |
| `selfgravity-jeans` | [Jeans modes](/applications/flume/verification/self-gravity/jeans): the unstable mode at amplitude $10^{-2}$, FFT solver, 20 steps | [#81](https://github.com/szaghi/adam/issues/81) |
| `selfgravity-amr-stage` | [self-gravity on a static AMR grid](/applications/flume/verification/self-gravity/amr-jeans): multigrid, the seam gravity work, per-stage reflux | [#81](https://github.com/szaghi/adam/issues/81) |
| `selfgravity-star` | [the n = 1 polytrope](/applications/flume/verification/self-gravity/polytrope), isolated boundaries, $16^3$ | [#81](https://github.com/szaghi/adam/issues/81) |
| `selfgravity-star-amr` | [the n = 1 polytrope on AMR](/applications/flume/verification/self-gravity/star-amr) (SV-14): isolated boundaries by the multigrid's James step, the star of radius 0.3 collapsing, the box refined over $32^3$, per-stage reflux, 20 steps | [#87](https://github.com/szaghi/adam/issues/87) |
| `sod-ps` | [passive scalars](/applications/flume/verification/euler/passive-scalars) (PS ps1): the `sod-x` shock tube with two scalars, 1 left and 0 right of the diaphragm | [#92](https://github.com/szaghi/adam/issues/92) |
| `radiative-shock` | [radiative cooling](/applications/flume/verification/euler/radiative-cooling) (CO co7): two streams at the cooling floor $T = 0.1$ collide, 100 cells to $t = 1$, $\Lambda = 0.18\,(T/0.1)^{1.5}$: the shock heating, the cooling layer, the dense layer | [#93](https://github.com/szaghi/adam/issues/93) |
| `noh-de` | [dual energy](/applications/flume/verification/euler/dual-energy) (DE de2): the cold planar Noh shock, Mach 845, with `physical_model = euler-de`, 50 cells, $t = 0.6$; every cell starts cold, so the entropy floor acts every step | [#91](https://github.com/szaghi/adam/issues/91) |
| `iso-tube` | [isothermal equation of state](/applications/flume/verification/euler/isothermal-eos) (EO eo1): the isothermal shock tube, $\rho = 1 \mid 0.125$, $v = 0.3 \mid -0.2$, 100 cells in four blocks, $y$ active, WENO5 flux splitting, $t = 0.2$: the isothermal library, the inert energy slot | [#100](https://github.com/szaghi/adam/issues/100) |
| `periodic-1block` | the isothermal acoustic wave of [EO](/applications/flume/verification/euler/isothermal-eos) on one block not refined at init, every face periodic, 48 cells, $t = 1$: the tree and the maps rebuilt after the BC types are known ([#103](https://github.com/szaghi/adam/issues/103)); equal to the 8-block run bitwise on CPU and FNL; the BC crown map with no row | [#103](https://github.com/szaghi/adam/issues/103) |
| `iso-collide-hllc` | [isothermal equation of state](/applications/flume/verification/euler/isothermal-eos) (EO eo2): two Mach-5 isothermal streams colliding, 100 cells, `weno-riemann` with the isothermal HLLC, $t = 0.5$: the isothermal Riemann kernels, the dense layer at $M_1^2$ | [#100](https://github.com/szaghi/adam/issues/100) |
| `baro-tube` | [barotropic equation of state](/applications/flume/verification/euler/barotropic-eos) (EO eb1): the tube $\rho = 10 \mid 0.5$ on the convex table `POLY`, its left fan crossing the break $\rho = 8$, 100 cells in four blocks, WENO5 flux splitting, $t = 0.15$: the table registry, the barotropic library | [#100](https://github.com/szaghi/adam/issues/100) |
| `baro-collide-hllc` | [barotropic equation of state](/applications/flume/verification/euler/barotropic-eos) (EO eb2): the collision $(1, 3.4) \mid (1.5, -2.6)$ on `POLY`, both shocks into the $\Gamma = 5/3$ segment, 100 cells, `weno-riemann` HLLC, $t = 0.3$: the barotropic Riemann kernels, the Roe secant speed | [#100](https://github.com/szaghi/adam/issues/100) |
| `baro-jeans-box` | [barotropic gases under gravity](/applications/flume/verification/euler/barotropic-gravity) (EO eg4): the barotropic Jeans collapse (`POLY` on $\rho_0 = 5$, amplitude $10^{-2}$), 32² base with the central box refined 2:1, the multigrid Poisson solver, the seam reflux, 20 steps: gravity on the momentum only, the energy row cleared | [#100](https://github.com/szaghi/adam/issues/100) |
| `iso-sinks-shu` | [barotropic gases under gravity](/applications/flume/verification/euler/barotropic-gravity) (EO eg6): the overdense singular isothermal sphere ($A = 29.3$, $c_0 = 0.5$) with a seed sink, isothermal, 16³, isolated self-gravity, the density-only positivity limiter, 30 steps: the sinks' $c_s^2 = c_0^2$; per-backend goldens | [#100](https://github.com/szaghi/adam/issues/100) |
| `iso-mhd-tube2` | [isothermal and barotropic MHD](/applications/flume/verification/mhd/barotropic-mhd) (EO em4): Mignone's isothermal tube 2 (rotational discontinuities), 200 cells in four blocks of one transverse cell, WENO5 flux splitting on Kim's eigenvectors, $t = 0.2$: the isothermal MHD library, the inert energy slot; per-backend goldens | [#100](https://github.com/szaghi/adam/issues/100) |
| `baro-mhd-glm-hlld` | [isothermal and barotropic MHD](/applications/flume/verification/mhd/barotropic-mhd) (EO em4): tube 2 with the densities times 5 on `POLY` (the $\Gamma = 1.4$ segment), mixed GLM, `weno-riemann` with the isothermal HLLD on $P(\rho)$, 200 cells, $t = 0.2$: the barotropic MHD library, the HLLD and GLM kernels; per-backend goldens | [#100](https://github.com/szaghi/adam/issues/100) |
| `sinks-shu` | [sink particles](/applications/flume/verification/sinks/sink-particles) (SI si11): the overdense singular isothermal sphere ($A = 29.3$, $\gamma = 1.001$) with a seed sink at the centre, $16^3$, 30 steps: the deposit into the Poisson source, the CIC gather, threshold accretion, the kick-drift-kick, the sink columns of the history | [#94](https://github.com/szaghi/adam/issues/94) |
| `tgv-small` | [the Taylor–Green vortex](/applications/flume/verification/turbulence/iles-baseline) (T0.2): Re 1600, Mach 0.4, $16^3$ cells on 8 blocks, 20 steps, with the turbulence history on | [#119](https://github.com/szaghi/adam/issues/119) |

The self-gravity goldens agree across the backends to $9\cdot10^{-14}$ relative in the volume integrals and $E + W$,
and to $1.8\cdot10^{-18}$ absolute in the momentum, which is round-off.


## PRISM suite

The suite is the structural-change baseline of PRISM. Each case runs in about a minute per backend and is checked
against committed goldens or, for the seam cases, against a physics oracle. The operational details (golden capture,
tolerance rationale) are also in the suite's
[`README.md`](https://github.com/szaghi/adam/blob/develop/src/tests/prism/regression/README.md).

::: tip Coverage
The suite has 24 case directories (an `input.ini` each). 15 carry a digest golden for the CPU backend and 16 for FNL
(`rmf-fwl` is FNL-goldened only), so a green sweep reports **15 PASS, 11 SKIP** on CPU and **16 PASS, 10 SKIP** on FNL.
A SKIP is a case with no golden for that backend, or one of the two non-case directories `reproducers/` and `slices/`.
The other cases are checked by their own `check.sh`, `scaling.sh` or `coils.sh`.
:::

### Design

- **Structural regressions, not physics.** The cases are CI-sized: `16³` or `8×16×16` cells per realm, `it_max` of
  5 to 10, per-step residual logging. A refactor that changes a reduction, a kernel chain or a seam fill shows up here.
- **Two backends, one baseline.** Every case runs on `prism-cpu-gnu` and `prism-fnl-nvf`; a case is fully covered only
  when it has goldens for both.
- **Fixed two-rank runs.** Every case runs at `mpirun -np 2`, so every run exercises the MPI ghost exchange and the
  decomposition.
- **Deterministic across compilers.** Runtime refinement is off (`amr_iterations = 0`, `[amr] frequency = 999999`), so
  no refinement decision depends on a floating-point comparison mid-run; the initial uniform refinement
  (`iu_ref_levels = 2`) gives 64 leaf blocks over the two ranks.

### PRISM cases

All cases use `fdv_order = 6`. The `rmf*`, `periodic-pulse-blocks` and `memory-scaling*` cases are `physical_model = EM`; the four `pic-*` cases are `physical_model = PIC`. The `pic-*`, `rmf-open-t0` and `rmf-pec-t0` cases use `ngc = 6`, every other case `ngc = 3`. All but `rmf-hyperbolic` run without divergence correction (`constrained_transport = no`, `divergence_correction = no`); `rmf-hyperbolic` is the only case with hyperbolic (GLM) cleaning, on D and B (`nv = 11`).

#### Digest/residuals-goldened cases

| Case | Realm | Cadence | Scheme | Exercises |
|---|---|---|---|---|
| `rmf`                  | single | — | fd | baseline EM case (4 AC coils); the digest anchor and the β cross-config reference target |
| `rmf-2realm`           | inter (x-split) | α | fd | single-realm `rmf` split at x=0; aggregated digest vs the `rmf` golden ([#13](https://github.com/szaghi/adam/issues/13) Phase D) |
| `rmf-2realm-asymK`     | inter | α | fd | the asymmetric-K validation: realm_1 SSP-33 (K=3) ∥ realm_2 SSP-54 (K=5), the α K-equality-guard removal ([#16](https://github.com/szaghi/adam/issues/16)) |
| `rmf-2realm-stagesync` | inter | β | fd | the load-bearing β oracle — digest must match the single-realm `rmf` golden bit-for-bit ([#18](https://github.com/szaghi/adam/issues/18)) |
| `rmf-fwl`              | single | — | fd | the only fWLayer case (six faces) + 4 coils — guards the FNL fWLayer host→device transfer fix ([#31](https://github.com/szaghi/adam/issues/31)); **FNL-goldened, CPU golden pending CI** |
| `rmf-fwl-amr`          | single | — | fd | `rmf-fwl` with an init-time box marker refining the x-minus slab to level 3: the fWLayer ranges rebuilt on the refined blocks ([#95](https://github.com/szaghi/adam/issues/95) P5.1) |
| `rmf-pml-amr`          | single | — | fd | the same grid with a CFS PML instead of the layer: the PML support and its RK storage rebuilt on the refined blocks ([#95](https://github.com/szaghi/adam/issues/95) P5.1) |
| `rmf-hyperbolic`       | single | — | fd | `rmf` with hyperbolic cleaning on D and B (`chi = 1`, damping `alpha = 1e9` 1/s): the GLM state layout, φ and ψ in slots 7 and 8, J in 9..11 ([#61](https://github.com/szaghi/adam/issues/61)). Its two init refusals (`hyperbolic` without constrained transport, `poisson`) are legs of `refusals.sh` |
| `periodic-pulse-blocks` | single | — | fd | a Gaussian pulse along +x in a box periodic on every face, 48³ cells as 64 blocks split 32/32 over the two ranks, 192 steps = one period: library periodicity across blocks and ranks |
| `rmf-pec-t0`           | single | — | fd | the four `rmf` coils carrying their current at t0 inside PEC walls: FLAIL's t0 solves of A give a divergence-free B ([#84](https://github.com/szaghi/adam/issues/84) E3b); its `check.sh` asserts `max\|div B\| h / max\|B\| ≤ 1e-12` and np 1 = np 2 to 1e-12 |
| `rmf-open-t0`          | single | — | fd | the same coils inside exact/open faces (`Neumann`): the free-space potential in the ghosts, summed on rank 0 ([#84](https://github.com/szaghi/adam/issues/84) E3c); same `check.sh` acceptance, away from the faces |
| `pic-pec-t0`           | single | — | fd | the t0 electric field of a PIC plasma box between PEC walls: `div D = rho` up to the checkerboard mode ([#84](https://github.com/szaghi/adam/issues/84) E3b); `check.sh` bound `1e-12 max\|rho\|`, np 1 = np 2 |
| `pic-open-t0`          | single | — | fd | the same plasma inside exact/open faces ([#84](https://github.com/szaghi/adam/issues/84) E3c); `check.sh` bound `1e-11 max\|rho\|` beyond 3 cells of the faces, np 1 = np 2 |

#### PIC cases (`check.sh`-driven, goldens in `golden-check/`)

Both run the production PIC path (TSC, conserving current `esirkepov_Q10` with the coherent gather, the t0 FLAIL solves,
FD6, SSP(5,4), fixed seed) through the shared driver `pic_check.sh`, as legs: CPU np 1 on one block for the first 4
steps, FNL np 1, and FNL np 2 on the same grid as 64 blocks. Every leg also asserts Gauss's law at t0 and charge
conservation to round-off. Their goldens sit in `golden-check/`, not `golden/`, because `run.sh` runs every goldened case
at `-np 2`, which the single-rank CPU leg refuses; `run.sh` reports them as SKIP.

| Case | Exercises |
|---|---|
| `pic-periodic`   | uniform plasma periodic on every face, plane wave of one box length, 32×24×24 cells, 40 steps: periodic current fold and particle wrap on three axes |
| `pic-slab-pulse` | a plasma slab 3 c/ω_pe thick, PEC along x and periodic in y, z, `gaussian_derivative_pulse` through the slab, 384×24×24 cells, 1152 steps |

#### Capacity probes (no golden)

`memory-scaling` runs a minimal EM case (`it_max = 1`) at `-np 1` and `-np 2` and asserts, through `scaling.sh`, how the
per-rank block capacity `nb` scales: it halves on the CPU (host memory shared by the ranks of a node) and stays constant
on FNL (one device per rank). `memory-scaling-coils` is the same case with the four `rmf` coils; `coils.sh` asserts that
the capacity budget counts 3 block-sized fields per coil.

#### Inter-realm 1:1 seam anchor

`rmf-2realm-fd-pulse` — a source-free Gaussian pulse split at x=0 into two same-resolution realms glued by a **β** mirror seam. It isolates the seam *mechanism* (inter-realm peer copy) from any resolution *jump*: a 1:1 seam must be div-free exactly like a 1:1 intra-block interface, and it holds `div(B) = div(D) = 0` — **but only under β**. Under α the seam ghosts go unfilled during RK substages and div(B) leaks (this is the [#31](https://github.com/szaghi/adam/issues/31) diagnosis). For a 1:1 same-`dt` seam, β is correct and required.

It carries a **`check.sh` oracle** (three legs, `-np 1`, CPU by default, `PRISM_EXE=` to override):

1. **β cadence, as committed** — `max|div(D)|` and `max|div(B)|` must both sit at round-off (`1.0E-13`) on **both** realms. Measured: **exactly `0.0`** for all four quantities.
2. **Health** — both realms reach 100% with no error/abort/NaN and both divergence histories exist.
3. **α cadence, negative control** — rewriting `coupling_cadence` to `end_of_step` must push `div(B)` *above* tolerance. Measured: **`4.56E+02`**, identically on both realms. Without this leg the first two could pass vacuously (a seam that stopped being exercised would also read zero), so it is what gives the oracle teeth. It asserts only *that* α leaks, not a magnitude — α's drift is not a quantity this suite has reason to hold fixed.

::: warning Not run by `run.sh`
Like the AMR-seam oracles, this `check.sh` is invoked **by hand**, not by the sweep — and the case has no digest golden, so `run.sh` reports it as a SKIP. Run it directly after touching `fill_seam_from_peer_forest`, `post_step_forest`, or any forest seam path.
:::

#### FV reflux on a capped last step

`rmf-amr-fv-pulse-timed` ([#38](https://github.com/szaghi/adam/issues/38)) — the source-free pulse of `rmf-amr-fd-pulse` on the FV path, a band x ∈ [-0.03, 0.03] refined 2:1 (two seam planes, both x boundaries coarse) and **time-driven**: `time_max` is 4.5 steps, so PRISM caps the fifth step at half a step. With J = 0 the FV update with reflux keeps the domain integrals of D and B constant, and the Berger-Colella correction must be scaled by the step the update used, not by the forest's uncapped `dt`. Its `check.sh` asserts that the last step is capped, that the reflux diagnostic is printed, and that every integral drifts by at most `1e-9` (relative to the largest D, or B, integral). The y extent is twice that of `rmf-amr-fd-pulse` (12σ each side): with the narrower box the Gaussian tails leak through the y faces and even the grid without seams drifts by 1e-8 per step. The reflux-sensitive integral is ∫Dz (its x-flux is Hy); ∫Bx has no x-flux and cannot see the defect. Measured on CPU and FNL: ∫Dz drifts 1.0e-11 over the four full steps (the seams' disturbance reaching the x boundaries) and 2.2e-11 after the capped one. **Negative control** (the fix reverted): the same drift up to step 4, then **1.4e-7** on the capped step, 6500 times larger; the tolerance sits 45× above the fixed run and 140× below the reverted one.

#### Inter-realm misaligned-block anchor

`rmf-2realm-fv-pulse-blocks` ([#51](https://github.com/szaghi/adam/issues/51)) — the 1:1 mirror pulse on the FV path (`fv_centered`, 20 steps) with the seam blocks of the two realms not lined up: one realm on 4 blocks per axis (the same cells), first realm 1 (the register side: each register face holds 2x2 realm-2 faces), then realm 2 (each realm-2 face scattered into 4 register faces). Its `check.sh` asserts that the register closes exactly, `max|F_coarse-F_fine_sum|` = 0 on every face and step (320 and 80 lines), and that the fields match the single-realm FV pulse within `1e-13`. They are not compared bitwise: PRISM FV is not invariant under a change of block layout at round-off (the single realm alone on 4 blocks per axis instead of 2 differs by 6.6e-16). Measured on CPU and FNL: mismatch 0, fields within 3.7e-15.

#### Inter-realm 2:1 seam anchor

`rmf-2realm-fd-pulse-refined` ([#52](https://github.com/szaghi/adam/issues/52)) — `rmf-amr-fd-pulse` split at its 2:1 face x=0 into a fine realm (x<0, one level finer, no marker) and a coarse realm (x>0) glued by `coupling = refined`, β cadence. The union is the single-realm cell set, so the split must reproduce `rmf-amr-fd-pulse` **bit for bit**, including its [#29](https://github.com/szaghi/adam/issues/29) seam div(B) source: the inter-realm 2:1 seam reuses the intra-realm formulas, so the source must be neither smaller nor larger. Its `check.sh` (`-np 2` by default, `--np N`, `PRISM_EXE=` to override) runs both and asserts:

1. **Health** — both runs reach 100% with no error/abort/NaN.
2. **Fields** — every checkpoint field (B, D, J, `div_*`, `res_*`) bitwise on every cell (`multirealm_oracle.py --fields-only --tol 0`). Measured: 147456 cells, 27 fields, CPU and FNL.
3. **Divergence histories** — at every step the single-realm `max|div(D)|` and `max|div(B)|` equal the larger of the two realms' values. Measured: `max|div(B)|` 6.86 → 12.65 in both, the coarse realm holding the maximum at steps 1–2 and the fine one at steps 3–5.
4. **Seam div(B) monitor** — `[IO].seam_divB_tol = 1` (warn-only) fires on the split.

At `-np 2` the refined split leaves one rank with no block of a realm, so the case also exercises the collective ghost update on an empty rank. Like the 1:1 anchor, the case has no digest golden and is run by hand.

#### AMR-seam cases (digest-goldened **and** `check.sh`-driven)

Three single-realm cases carry a static intra-realm **2:1 AMR jump** (`markers_number = 1`, an AMR_GEO box marker covering the x<0 half). They are doubly covered: digest + residuals goldens on **both** backends ([#24](https://github.com/szaghi/adam/issues/24)), which `run.sh` checks like any other case, **plus** a bespoke `check.sh` asserting a physics oracle with a marker-off control leg for contrast. The oracles are *not* invoked by `run.sh` — run them directly (they are the only thing that checks the seam `div(B)` behaviour).

| Case | Scheme | `check.sh` oracle |
|---|---|---|
| `rmf-amr`          | `fv_centered` | refinement happens → 2:1 jump exists; the AMR seam registers (`SEAM_KIND_INTRA_REALM_AMR`); Berger-Colella reflux fires non-zero; np1 ≡ np2 register parity ([#28](https://github.com/szaghi/adam/issues/28)) |
| `rmf-amr-fd`       | `fd_centered` | control `max\|div(B)\| ≤ 1e-13` (interior identity intact); seam `max\|div(B)\|` within ±5% of a pinned baseline; div(J) truthfulness within a factor-3 band ([#26](https://github.com/szaghi/adam/issues/26)) |
| `rmf-amr-fd-pulse` | `fd_centered` | source-free two-invariant oracle (`div(D)` structural-zero, `div(B)` seam baseline ±5%); a `--convergence` `p_obs ≥ 0.8` refinement ladder; and the RK-contract legs — leap-splitting refusal + SSP-11 instrument ([#25](https://github.com/szaghi/adam/issues/25)) |

The pulse case is the valid convergence instrument: with `J = 0` the Maxwell system conserves both discrete divergences, so any seam `div(B)` is attributable to the 2:1 jump alone (no coil `h⁻²` contamination). Its `--convergence` leg is what pins the seam `div(B)` source as O(h^p) truncation-order (see [`div(B)` at 2:1 seams](/concepts/forest#the-2-1-seam-is-not-div-free-accept-truncation)).

### Reproducers

`src/tests/prism/regression/reproducers/` holds standalone minimal reproducers (not cases): `fnl_mpi_rndv.sh` for the FNL UCX-rendezvous SIGABRT on WSL ([#12](https://github.com/szaghi/adam/issues/12), the `UCX_RNDV_THRESH=inf` workaround) and a gfortran-16 debug-mode `-fcheck=bounds` false-positive reproducer.


### Research cases

Research and development inputs live under `src/tests/prism/<backend>/{fd,fv}/<family>/`: long integrations, full AMR,
sized for physics validation. They are not regression anchors, and `run.sh` never descends into them.

## Unit tests

The library and application unit tests are separate build modes (`test-*` in [build modes](/reference/build-modes)),
run as `mpirun -np <N> exe/<test_executable>`. The families are the FDV operators on trigonometric and step fields
(`test-fdv-*`, CPU and FNL), AMR (`test-amr-*`, `test-tree-closest-block`, `test-periodic-ghost`), seams
(`test-seam-interpolation`, `test-quadtree-seam-ghost`), the Poisson solvers and FLAIL (`test-poisson-fft`,
`test-poisson-mg`, `test-flail-bc`), the RK tableaux (`test-rk-ssp-tableau`), PRISM reflux (`test-prism-reflux-apply`),
turbulence (`test-sgs`, `test-kep`, `test-wm`, `test-df`), and FLUME's libraries and user criteria (`test-flume-*`,
`test-ripple-crest-criterion`). The tests the FLUME cases rely on are described in
[FLUME unit tests](/applications/flume/verification/unit-tests).

Changes to core numerics or parallel logic are gated on both maintained toolchains (versions in
[install](/start/install)): gfortran for the CPU backend (the `*-debug` modes add bounds checking and FPE traps) and
nvfortran for FNL.

## See also

- [Regression cookbook](/cookbook/regression): run a suite, add a case, refresh a golden, bring up a new machine.
- [Forest (multi-realm)](/concepts/forest): the seam machinery the PRISM cases exercise.
- [`kernel_inventory.md`](https://github.com/szaghi/adam/blob/develop/src/tests/prism/regression/kernel_inventory.md):
  the PRISM-FNL OpenACC kernel-class audit.
