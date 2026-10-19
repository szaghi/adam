# Regression suite

The verification cases assert the physics: an exact solution, a reference, a symmetry or a conservation property. The
regression suite does a different job. It freezes 39 inputs, most of them copies of a verification input cut to a few
dozen steps, and detects **any** change of their outputs, wanted or not.

## How it works

Each case directory `src/tests/flume/regression/<case>/` holds a frozen `input.ini` (a multi-realm case: a forest
manifest plus one INI per realm) and the goldens `golden/cpu/` and `golden/fnl/`:

| Golden | Content | Comparison |
|---|---|---|
| `digest.txt` | per-variable digest of every checkpoint field (`src/tests/prism/regression/digest.py`) | within the digest tolerances (default rtol $10^{-6}$) |
| `*-residuals.dat` | per-iteration residuals | header verbatim, integer columns exact, float columns within the tolerances |
| `*-conservation_history.dat` | per-iteration volume integrals of the conservative variables (and $W$, $E + W$ under self-gravity) | as the residuals |

The raw HDF5 checkpoints are not committed. Every case runs on 2 MPI ranks.

| Script | Role |
|---|---|
| `src/tests/flume/regression/run.sh cpu` | builds `flume-cpu-gnu` and `flume-cpu-gnu-fft`, runs every case, compares; runs in CI |
| `src/tests/flume/regression/run-fnl-local.sh` | the same on the GPU (`flume-fnl-nvf` and `flume-fnl-nvf-fft`, `--varset local_nvf`); a workstation gate, since the CI runners have no GPU |
| `src/tests/flume/regression/run-omp-bitwise.sh` | checks that the OpenMP build reproduces the serial one bit for bit (CI) |

Two executables per backend: every case with `self_gravity = .true.` runs on the `-fft` executable (FFTW on CPU, cuFFT
on FNL), since every self-gravity solve needs an FFT; the CI job installs `libfftw3-dev`. `--no-build` skips the
build; `REGRESSION_CASES="a b"` restricts the run; `REGRESSION_RUN_GOLDENLESS=1` runs a case that has no golden yet (the
capture of a new case). Before running, the harness refuses an executable older than any of its sources.

A case also fails when its run exits non-zero, a crash or the hard timeout, even if its outputs match the goldens: a
crash after the last checkpoint is still a crash (hcoll segfaulting in `MPI_Finalize` on a cluster passed unseen
before). The outputs are compared anyway, for the report. The one exemption is the harness's own termination of a run
that hung in `MPI_Finalize` after its end-of-run marker (the WSL GPU stack). `REGRESSION_ALLOW_MPIRUN_FAILURE=1`
turns the failure into a warning on an environment whose crash is known.

```bash
src/tests/flume/regression/run.sh cpu                     # build, run, compare (CPU)
src/tests/flume/regression/run-fnl-local.sh               # the same on the GPU
REGRESSION_CASES="sod-x rotor" src/tests/flume/regression/run.sh cpu --no-build
```

## Cases

| Case | Freezes | Issue |
|---|---|---|
| `sod-x`, `sod-y`, `sod-z` | [Sod shock tube](./euler/sod), along each axis | [#35](https://github.com/szaghi/adam/issues/35) |
| `sod-hllc` | [Sod](./euler/sod) on `weno-riemann`, HLLC | [#47](https://github.com/szaghi/adam/issues/47) |
| `vortex-periodic` | [isentropic vortex](./euler/isentropic-vortex), 100 steps | [#35](https://github.com/szaghi/adam/issues/35) |
| `amr-periodic-reflux` | [conservation across AMR seams](./euler/amr-conservation): init-time AMR, reflux | [#35](https://github.com/szaghi/adam/issues/35) |
| `amr-periodic-reflux-stage` | the same with `reflux_cadence = stage` | [#81](https://github.com/szaghi/adam/issues/81) |
| `shock-cylinder-ib` | [shock over a cylinder](./euler/shock-cylinder), immersed boundary, 120 steps | [#35](https://github.com/szaghi/adam/issues/35) |
| `sod-2realm`, `sod-2realm-blocks`, `sod-amr`, `sod-amr-2realm`, `sod-amr-refined` | [multi-realm runs](./euler/multi-realm): mirror seams lined up or not, a 2:1 face, the refined coupling | [#37](https://github.com/szaghi/adam/issues/37), [#51](https://github.com/szaghi/adam/issues/51), [#52](https://github.com/szaghi/adam/issues/52) |
| `rj2a-x`, `rj2a-y`, `rj2a-z` | [Ryu–Jones 2a](./mhd/ryu-jones-2a) | [#41](https://github.com/szaghi/adam/issues/41) |
| `rj2a-hlld` | Ryu–Jones 2a on `weno-riemann`, HLLD | [#47](https://github.com/szaghi/adam/issues/47) |
| `rj2a-2realm` | Ryu–Jones 2a in two realms | [#41](https://github.com/szaghi/adam/issues/41) |
| `brio-wu` | [Brio–Wu](./mhd/brio-wu-rj4d) | [#41](https://github.com/szaghi/adam/issues/41) |
| `glm-pulse` | [GLM cleaning](./mhd/glm-cleaning) | [#41](https://github.com/szaghi/adam/issues/41) |
| `orszag-tang` | [Orszag–Tang](./mhd/orszag-tang) | [#41](https://github.com/szaghi/adam/issues/41) |
| `orszag-tang-eglm` | Orszag–Tang with [EGLM](./mhd/eglm) | [#47](https://github.com/szaghi/adam/issues/47) |
| `rotor` | [MHD rotor](./mhd/rotor) | [#41](https://github.com/szaghi/adam/issues/41) |
| `field-loop` | [field loop](./mhd/field-loop) | [#41](https://github.com/szaghi/adam/issues/41) |
| `rotated-shock-tube` | [rotated shock tube](./mhd/rotated-shock-tube) | [#41](https://github.com/szaghi/adam/issues/41) |
| `uniform-amr-mhd` | MHD on a uniform state across 2:1 faces | [#41](https://github.com/szaghi/adam/issues/41) |
| `blast-limiter`, `blast3d-limiter` | the Balsara–Spicer blast with the [positivity limiter](./mhd/positivity-limiter), 2-D and 3-D | [#47](https://github.com/szaghi/adam/issues/47) |
| `blast-amr-limiter` | the limited blast across two 2:1 seams | [#50](https://github.com/szaghi/adam/issues/50) |
| `shear-viscous` | [viscous waves](./dissipative/viscous-waves) | [#65](https://github.com/szaghi/adam/issues/65) |
| `couette` | [Couette flow](./dissipative/couette) | [#65](https://github.com/szaghi/adam/issues/65) |
| `becker-shock` | [Becker's viscous shock](./dissipative/becker-shock) | [#65](https://github.com/szaghi/adam/issues/65) |
| `alfven-resistive` | [the visco-resistive Alfvén wave](./dissipative/ohmic-alfven) | [#65](https://github.com/szaghi/adam/issues/65) |
| `seam-viscous-amr` | [the dissipative fluxes across 2:1 seams](./dissipative/seam-conservation) | [#65](https://github.com/szaghi/adam/issues/65) |
| `vortex-regrid` | [runtime regrids](./runtime-amr/runtime-regrids): the vortex, 6 regrids | [#74](https://github.com/szaghi/adam/issues/74) |
| `field-loop-regrid` | [runtime regrids](./runtime-amr/runtime-regrids) on MHD: the field loop, 10 regrids | [#74](https://github.com/szaghi/adam/issues/74) |
| `selfgravity-jeans` | [Jeans modes](./self-gravity/jeans): the unstable mode at amplitude $10^{-2}$, FFT solver, 20 steps | [#81](https://github.com/szaghi/adam/issues/81) |
| `selfgravity-amr-stage` | [self-gravity on a static AMR grid](./self-gravity/amr-jeans): multigrid, the seam gravity work, per-stage reflux | [#81](https://github.com/szaghi/adam/issues/81) |
| `selfgravity-star` | [the n = 1 polytrope](./self-gravity/polytrope), isolated boundaries, $16^3$ | [#81](https://github.com/szaghi/adam/issues/81) |
| `selfgravity-star-amr` | [the n = 1 polytrope on AMR](./self-gravity/polytrope) (SV-14): isolated boundaries by the multigrid's James step, the star of radius 0.3 collapsing, the box refined over $32^3$, per-stage reflux, 20 steps | [#87](https://github.com/szaghi/adam/issues/87) |
| `sod-ps` | [passive scalars](./euler/passive-scalars) (PS ps1): the `sod-x` shock tube with two scalars, 1 left and 0 right of the diaphragm | [#92](https://github.com/szaghi/adam/issues/92) |
| `radiative-shock` | [radiative cooling](./euler/radiative-cooling) (CO co7): two streams at the cooling floor $T = 0.1$ collide, 100 cells to $t = 1$, $\Lambda = 0.18\,(T/0.1)^{1.5}$: the shock heating, the cooling layer, the dense layer | [#93](https://github.com/szaghi/adam/issues/93) |
| `noh-de` | [dual energy](./euler/dual-energy) (DE de2): the cold planar Noh shock, Mach 845, with `physical_model = euler-de`, 50 cells, $t = 0.6$; every cell starts cold, so the entropy floor acts every step | [#91](https://github.com/szaghi/adam/issues/91) |
| `sinks-shu` | [sink particles](./self-gravity/sink-particles) (SI si11): the overdense singular isothermal sphere ($A = 29.3$, $\gamma = 1.001$) with a seed sink at the centre, $16^3$, 30 steps: the deposit into the Poisson source, the CIC gather, threshold accretion, the kick-drift-kick, the sink columns of the history | [#94](https://github.com/szaghi/adam/issues/94) |

The self-gravity goldens agree across the backends to $9\cdot10^{-14}$ relative in the volume integrals and $E + W$,
and to $1.8\cdot10^{-18}$ absolute in the momentum, which is round-off.

## Refreshing a golden

Only after an intentional change, adjudicated case by case (which fields move, when, by how much, and why). The run
fails against the old golden but leaves its outputs in `work-cpu/` (`work-fnl/` on the GPU):

```bash
REGRESSION_CASES=<case> src/tests/flume/regression/run.sh cpu --no-build
mkdir -p src/tests/flume/regression/<case>/golden/cpu
cp src/tests/flume/regression/<case>/work-cpu/{digest.txt,*-residuals.dat,*-conservation_history.dat} \
   src/tests/flume/regression/<case>/golden/cpu/
```
