# Verification

Every FLUME capability is backed by a verification case in `src/tests/flume/verification/`: a `check.sh` that runs
the case (CPU by default; `FLUME_EXE=exe/adam_flume_fnl` for the GPU backend) and an oracle that asserts the physics
(an exact solution, a reference, a symmetry or a conservation property), not a stored golden. The stored goldens are
the job of the [regression suite](./regression); the building blocks below the solver are pinned by the
[unit tests](./unit-tests).

Each case has its own page with the same sections: the case and its scope, the models and numerics with references
to the code, the initial and boundary conditions, the input files, the results on CPU and FNL, and the findings.
All runs use 2 MPI ranks and WENO5 unless a page says otherwise. The figures are drawn by
`src/tests/flume/verification/make_doc_figures.py`, from the work directories of the `check.sh` scripts or from the
re-runs of `make_doc_runs.sh`.

Run any case from its directory:

```bash
cd src/tests/flume/verification/sod && ./check.sh                  # CPU backend
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh        # FNL backend (GPU)
```

In the catalog, **Figure** "no" marks a contract or bitwise leg (an input refusal, CPU = FNL, a restart, a
round trip): a setup table documents it better than a plot of a difference that is zero.


## [Euler](./euler/)

| Case | ID | Milestone | Backends | Figure |
|---|---|---|---|---|
| [Sod shock tube](./euler/sod) | V1, RV-3 | M1 ([#35](https://github.com/szaghi/adam/issues/35)), M3 ([#47](https://github.com/szaghi/adam/issues/47)) | CPU, FNL | yes |
| [Lax problem](./euler/lax) | RV-3 (`lax` leg) | M3 ([#47](https://github.com/szaghi/adam/issues/47)) | CPU (all schemes), FNL (`weno-riemann` HLLC) | yes |
| [Shu–Osher shock / entropy-wave interaction](./euler/shu-osher) | RV-3 (`shu-osher` leg) | M3 ([#47](https://github.com/szaghi/adam/issues/47)) | CPU (all schemes), FNL (`weno-riemann` HLLC) | yes |
| [Isentropic vortex: order of accuracy](./euler/isentropic-vortex) | V2, RV-2 | M1 ([#35](https://github.com/szaghi/adam/issues/35)), M3 ([#47](https://github.com/szaghi/adam/issues/47)) | CPU (all schemes), FNL (`weno`, `weno-riemann` LLF and HLLC) | yes |
| [Shock over a cylinder](./euler/shock-cylinder) | V6 | M1 ([#35](https://github.com/szaghi/adam/issues/35)), M3 ([#47](https://github.com/szaghi/adam/issues/47)) | CPU, FNL | yes |
| [Conservation across AMR seams](./euler/amr-conservation) | V3, RV-4 | M1 ([#35](https://github.com/szaghi/adam/issues/35)), M3 ([#47](https://github.com/szaghi/adam/issues/47)) | CPU, FNL | yes |
| [Multi-realm runs](./euler/multi-realm) | multirealm legs 1–6 | M2 entry ([#37](https://github.com/szaghi/adam/issues/37)); [#40](https://github.com/szaghi/adam/issues/40), [#51](https://github.com/szaghi/adam/issues/51), [#52](https://github.com/szaghi/adam/issues/52), [#54](https://github.com/szaghi/adam/issues/54); M6 ([#68](https://github.com/szaghi/adam/issues/68)) | CPU, FNL | yes |
| [Forward-facing step: a three-realm forest on quadtrees](./euler/forward-facing-step) | EV-step | — ([#46](https://github.com/szaghi/adam/issues/46), quadtree AMR; [#66](https://github.com/szaghi/adam/issues/66)) | CPU, FNL (`trees` leg); CPU only (`full` leg) | yes |

## [Ideal MHD](./mhd/)

| Case | ID | Milestone | Backends | Figure |
|---|---|---|---|---|
| [MHD on `weno-riemann`](./mhd/weno-riemann) | RV-5, RV-6, RV-7 | M3 ([#47](https://github.com/szaghi/adam/issues/47), M3-P3c) | CPU, FNL | yes |
| [Ryu–Jones 2a: all seven waves](./mhd/ryu-jones-2a) | MV-4 | M2 ([#41](https://github.com/szaghi/adam/issues/41)) | CPU, FNL | yes |
| [Brio–Wu and Ryu–Jones 4d](./mhd/brio-wu-rj4d) | — (`mhd/riemann`, M2-P3 exit) | M2 ([#41](https://github.com/szaghi/adam/issues/41)) | CPU, FNL | yes |
| [MHD plumbing](./mhd/plumbing) | — (`mhd/plumbing`, M2-P1) | M2 ([#41](https://github.com/szaghi/adam/issues/41), P1) | CPU, FNL | no |
| [MHD with a zero field](./mhd/zero-field) | MV-2 | M2 ([#41](https://github.com/szaghi/adam/issues/41), P3) | CPU, FNL | no |
| [GLM divergence cleaning](./mhd/glm-cleaning) | MV-3, MV-10 | M2 ([#41](https://github.com/szaghi/adam/issues/41)); the MV-10 monitor M6 ([#68](https://github.com/szaghi/adam/issues/68), G1) | CPU, FNL | yes |
| [Order of accuracy: linear waves, CPAW, magnetised vortex](./mhd/order-of-accuracy) | MV-5, MV-6, MV-7 | M2 ([#41](https://github.com/szaghi/adam/issues/41)) | CPU, FNL | yes |
| [MHD conservation across AMR seams](./mhd/conservation) | MV-11 | M2 ([#41](https://github.com/szaghi/adam/issues/41), P7a) | CPU, FNL | yes |
| [Orszag–Tang vortex](./mhd/orszag-tang) | MV-12 | M2 ([#41](https://github.com/szaghi/adam/issues/41)) | CPU, FNL | yes |
| [MHD rotor](./mhd/rotor) | MV-13 | M2 ([#41](https://github.com/szaghi/adam/issues/41)) | CPU, FNL | yes |
| [Field loop advection](./mhd/field-loop) | MV-9 | M2 ([#41](https://github.com/szaghi/adam/issues/41)); `--agree` [#79](https://github.com/szaghi/adam/issues/79) | CPU, FNL | yes |
| [Rotated shock tube](./mhd/rotated-shock-tube) | MV-8 | M2 ([#41](https://github.com/szaghi/adam/issues/41)) | CPU, FNL | yes |
| [EGLM cleaning](./mhd/eglm) | EV-2, EV-3, EV-4 | M3 ([#47](https://github.com/szaghi/adam/issues/47), M3-P4) | CPU, FNL | yes |
| [Positivity limiter](./mhd/positivity-limiter) | PV-1 to PV-5 | M3 ([#47](https://github.com/szaghi/adam/issues/47), M3-P5); PV-5 [#50](https://github.com/szaghi/adam/issues/50) | CPU, FNL (PV-1 as a sphere: FNL only) | yes |
| [Quadtree AMR](./mhd/quadtree-amr) | MV-15 | [#46](https://github.com/szaghi/adam/issues/46) (no milestone, after M3) | CPU, FNL | yes |

## [Infrastructure](./infrastructure/)

| Case | ID | Milestone | Backends | Figure |
|---|---|---|---|---|
| [Ghost probe: every face and edge ghost on a linear field](./infrastructure/ghost-probe) | GP | M4 ([#65](https://github.com/szaghi/adam/issues/65) P0, P1), M6 ([#68](https://github.com/szaghi/adam/issues/68) C2, D1) | CPU, FNL | yes |
| [Scaling covariance and the reference layer](./infrastructure/scaling-covariance) | NV-0, NV-4, NV-5, NV-6, NV-8, NV-9 | — ([#49](https://github.com/szaghi/adam/issues/49), non-dimensional form, N0–N2) | CPU, FNL | no |
| [Ranks without blocks](./infrastructure/ranks-without-blocks) | DV-1 | — ([#82](https://github.com/szaghi/adam/issues/82)) | CPU, FNL | yes |
| [Non-finite state guard](./infrastructure/nonfinite) | — (`nonfinite`) | — ([#45](https://github.com/szaghi/adam/issues/45)) | CPU, FNL | no |

## [Dissipative terms (M4)](./dissipative/)

| Case | ID | Milestone | Backends | Figure |
|---|---|---|---|---|
| [Input contract of the dissipative terms](./dissipative/input-contract) | DC | M4 ([#65](https://github.com/szaghi/adam/issues/65)) | CPU, FNL | no |
| [Viscous waves](./dissipative/viscous-waves) | VV-1, VV-2 | M4 ([#65](https://github.com/szaghi/adam/issues/65)) | CPU, FNL | yes |
| [Compressible Couette flow](./dissipative/couette) | VV-3 | M4 ([#65](https://github.com/szaghi/adam/issues/65)) | CPU, FNL | yes |
| [Becker's viscous shock](./dissipative/becker-shock) | VV-4 | M4 ([#65](https://github.com/szaghi/adam/issues/65)) | CPU, FNL | yes |
| [Ohmic decay and the visco-resistive Alfvén wave](./dissipative/ohmic-alfven) | VV-5 | M4 ([#65](https://github.com/szaghi/adam/issues/65)) | CPU, FNL | yes |
| [Ohmic heating and the divergence constraint](./dissipative/ohmic-heating) | VV-6 | M4 ([#65](https://github.com/szaghi/adam/issues/65)) | CPU, FNL | no |
| [Conservation across 2:1 seams and forest seams](./dissipative/seam-conservation) | VV-7 | M4 ([#65](https://github.com/szaghi/adam/issues/65)) | CPU, FNL | yes |
| [Accuracy across 2:1 seams](./dissipative/seam-accuracy) | VV-8 | M4 ([#65](https://github.com/szaghi/adam/issues/65)) | CPU, FNL | yes |
| [CPU and FNL agree run by run](./dissipative/cpu-fnl-agreement) | VV `agree` | M4 ([#65](https://github.com/szaghi/adam/issues/65)) | CPU, FNL | no |

## [Runtime AMR (M5)](./runtime-amr/)

| Case | ID | Milestone | Backends | Figure |
|---|---|---|---|---|
| [Regrid round trip and input contract](./runtime-amr/regrid-round-trip) | RG (rg0, rg1, rg2) | M5 ([#74](https://github.com/szaghi/adam/issues/74)) | CPU (rg0, a host program); rg1 and rg2 CPU, FNL | no |
| [Runtime regrids in FLUME](./runtime-amr/runtime-regrids) | RG (rg3 to rg6) | M5 ([#74](https://github.com/szaghi/adam/issues/74)) | CPU, FNL | yes |
| [Accuracy of runtime AMR](./runtime-amr/accuracy) | AV (av2 to av8) | M5 ([#74](https://github.com/szaghi/adam/issues/74)) | CPU, FNL | yes |
| [MHD blast symmetry on the flux-splitting path](./runtime-amr/mhd-blast-symmetry) | #77 study, bounded by av3c | M5 ([#74](https://github.com/szaghi/adam/issues/74)), [#77](https://github.com/szaghi/adam/issues/77) | CPU, FNL | yes |

## [Seam accuracy and divergence (M6)](./seams/)

| Case | ID | Milestone | Backends | Figure |
|---|---|---|---|---|
| [Order across static 2:1 seams](./seams/order-across-seams) | AO (A1, C1) | M6 ([#68](https://github.com/szaghi/adam/issues/68)) | CPU, FNL | yes |
| [Divergence across seams](./seams/divergence-across-seams) | AD (A2) | M6 ([#68](https://github.com/szaghi/adam/issues/68)), [#78](https://github.com/szaghi/adam/issues/78) | FNL (tables); CPU for the loop leg and the D_flux cross-check | yes |
| [Shock crossing a static 2:1 face](./seams/shock-crossing) | AS (A3, C3, C3b, C3c) | M6 ([#68](https://github.com/szaghi/adam/issues/68)) | CPU, FNL (Sod); FNL only (RJ2a, C3, C3c) | yes |
| [Integral drift and the WENO sensor near seams](./seams/integral-drift) | AL (A4) | M6 ([#68](https://github.com/szaghi/adam/issues/68)) | CPU, FNL (the two backends log the same values) | yes |

## [Seam studies (M6)](./seams-studies/)

| Case | ID | Milestone | Backends | Figure |
|---|---|---|---|---|
| [Quintic coarse-to-fine fill](./seams-studies/quintic-fill) | D1 | M6 ([#68](https://github.com/szaghi/adam/issues/68)) | CPU, FNL (table: FNL) | yes |
| [Limited coarse-to-fine fills](./seams-studies/limited-fills) | D2 | M6 ([#68](https://github.com/szaghi/adam/issues/68)) | FNL (crossings); CPU, FNL (AO legs) | yes |
| [Shock-crossing defect is a capturing error](./seams-studies/capturing-error) | E0 | M6 ([#68](https://github.com/szaghi/adam/issues/68)) | FNL; CPU re-runs agree to the printed digits | yes |
| [Reflux policy: troubled against always](./seams-studies/reflux-policy) | E2 (E1) | M6 ([#68](https://github.com/szaghi/adam/issues/68)) | FNL (flags checked equal to the CPU's); AO `euler-tq`: CPU, FNL | yes |
| [#78 is the reflux, not the restricted ghosts](./seams-studies/reflux-not-ghosts) | G3a | M6 ([#68](https://github.com/szaghi/adam/issues/68)), [#78](https://github.com/szaghi/adam/issues/78) | FNL only | yes |
| [GLM damping length on a seam](./seams-studies/glm-damping) | G2 | M6 ([#68](https://github.com/szaghi/adam/issues/68)) | FNL only | yes |
| [EGLM reflux leaves B alone](./seams-studies/eglm-reflux) | G3b | M6 ([#68](https://github.com/szaghi/adam/issues/68)), [#78](https://github.com/szaghi/adam/issues/78) | CPU, FNL | yes |
| [2:1 caps in a 2-D prototype](./seams-studies/caps-prototype) | B1 | M6 ([#68](https://github.com/szaghi/adam/issues/68)) | NumPy prototype (no ADAM run) | yes |
| [Sensor-gated reflux prototype](./seams-studies/sensor-gated-reflux) | B2 | M6 ([#68](https://github.com/szaghi/adam/issues/68)) | NumPy prototype (no ADAM run) | yes |
| [Footprint of a high-order restriction](./seams-studies/high-order-restriction) | B3 | M6 ([#68](https://github.com/szaghi/adam/issues/68)) | NumPy prototype (no ADAM run) | yes |

## [External gravity (M7)](./gravity/)

| Case | ID | Milestone | Backends | Figure |
|---|---|---|---|---|
| [Stratified atmospheres and the `hydrostatic` boundary](./gravity/atmospheres) | GV (G0) | M7 ([#73](https://github.com/szaghi/adam/issues/73)) | CPU, FNL | yes |
| [The gravity force](./gravity/force) | GV-1, GV-2, GV-3, GV-6, GV-7 (G1) | M7 ([#73](https://github.com/szaghi/adam/issues/73)) | CPU, FNL | yes |
| [The order of the energy source](./gravity/energy-order) | GV-9 (G2a) | M7 ([#73](https://github.com/szaghi/adam/issues/73)) | CPU, FNL | yes |
| [The well-balanced mode](./gravity/well-balanced) | GV-2b, GV-5, GV-8, GV-10, GV-7 (G2b) | M7 ([#73](https://github.com/szaghi/adam/issues/73)) | CPU, FNL | yes |
| [Rayleigh–Taylor](./gravity/rayleigh-taylor) | GV-4 (G3) | M7 ([#73](https://github.com/szaghi/adam/issues/73)); symmetry study [#83](https://github.com/szaghi/adam/issues/83) | CPU, FNL | yes |

## [Self-gravity (M8)](./self-gravity/)

| Case | ID | Milestone | Backends | Figure |
|---|---|---|---|---|
| [Gravity from a cell potential](./self-gravity/potential) | SV-1 | M8 ([#81](https://github.com/szaghi/adam/issues/81)), P1 | CPU, FNL | yes |
| [The FFT Poisson solver](./self-gravity/fft-poisson) | SV-2 | M8 ([#81](https://github.com/szaghi/adam/issues/81)), P2 | CPU, FNL | yes |
| [Jeans modes under self-gravity](./self-gravity/jeans) | SV-3, SV-4, SV-5 | M8 ([#81](https://github.com/szaghi/adam/issues/81)), P3 | CPU, FNL | yes |
| [The isolated Poisson solver](./self-gravity/isolated-solver) | SV-6 | M8 ([#81](https://github.com/szaghi/adam/issues/81)), P4 | CPU, FNL | yes |
| [The n = 1 polytrope: equilibrium and collapse](./self-gravity/polytrope) | SV-7, SV-8 | M8 ([#81](https://github.com/szaghi/adam/issues/81)), P4 | CPU, FNL | yes |
| [The multigrid Poisson solver on AMR grids](./self-gravity/multigrid-amr) | SV-9 | M8 ([#81](https://github.com/szaghi/adam/issues/81)), P5a | CPU only | yes |
| [Jeans modes on a static AMR grid](./self-gravity/amr-jeans) | SV-10 | M8 ([#81](https://github.com/szaghi/adam/issues/81)), P5b, P5c | CPU, FNL | yes |
| [Self-gravity under runtime AMR](./self-gravity/runtime-amr) | SV-11 | M8 ([#81](https://github.com/szaghi/adam/issues/81)), P5c | CPU, FNL | yes |
