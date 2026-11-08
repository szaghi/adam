# Verification

Every FLUME capability is backed by a verification case in `src/tests/flume/verification/`: a `check.sh` that runs the
case and an oracle that asserts the physics (an exact solution, a reference, a symmetry or a conservation property),
not a stored golden. The goldens are the job of the [regression suite](/reference/regression#flume-suite); the building
blocks below the solver are pinned by the [unit tests](./unit-tests).

The cases are listed below in the order they were implemented, milestone by milestone, so the list reads as the path
the solver took. Each case page has the same sections: case and scope, models and numerics with references to the
code, initial and boundary conditions, input files, results on CPU and FNL, findings. All runs use 2 MPI ranks and
WENO5 unless a page says otherwise.

**Type** marks what a page verifies:

- **P**, physics: a run against an exact or reference solution, an order of accuracy or a symmetry;
- **C**, contract: a bitwise identity, rank or backend invariance (CPU = FNL), a restart, conservation to round-off,
  an input refusal.

**Figure** "no" marks a contract page: a setup table documents it better than a plot of a difference that is zero.

## Running a case

```bash
cd src/tests/flume/verification/sod && ./check.sh                  # CPU backend
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh        # FNL backend (GPU)
./check.sh --legs vv1,vv3                                          # the listed legs (commas or blanks)
./check.sh --legs all                                              # every leg
```

The cases are complete simulations: run the legs a change touches, not the whole catalog; the regression suite is the
layer that runs on every change. Every script selects its legs through the shared `legs.sh`, and lists its legs and
its default selection in its `Usage:` header. An unknown leg, or a selection that names no leg, stops the script with
exit status 2 before anything runs; each log opens with the resolved selection (`>> legs: ...`). The figures are drawn
by `make_doc_figures.py`, from the work directories of the `check.sh` scripts or from the re-runs of
`make_doc_runs.sh`; each group page lists its commands.

## M1 · Euler core (#35)

Issue [#35](https://github.com/szaghi/adam/issues/35). The Euler solver on one realm: Riemann problems, the design
order, the immersed boundary, the reflux at 2:1 faces and the restart. [Group page](./euler/).

| Page | IDs / legs | Type | Figure |
|---|---|---|---|
| [Sod shock tube](./euler/sod) | V1, RV-3 | P | yes |
| [Isentropic vortex: order of accuracy](./euler/isentropic-vortex) | V2, RV-2 | P | yes |
| [Shock over a cylinder](./euler/shock-cylinder) | V6 | P | yes |
| [Conservation across AMR seams](./euler/amr-conservation) | V3, RV-4 | C | yes |
| [Restart round trip and output channels](./euler/restart) | V7 | C | no |

## #37 · Forest entry

Issue [#37](https://github.com/szaghi/adam/issues/37). Several realms glued by inter-realm seams run as one forest,
bitwise equal to a single realm. The same page carries the legs of the later seam issues and the MHD twin MV-14 of M2.

| Page | IDs / legs | Type | Figure |
|---|---|---|---|
| [Forest seams: the inter-realm contract](./euler/multi-realm) | multirealm legs 1–6, 15; MV-14 | C | yes |

## M2 · Ideal MHD (#41)

Issue [#41](https://github.com/szaghi/adam/issues/41). The conservative MHD system with cell-centred $\mathbf{B}$ and
mixed GLM cleaning, on the flux-splitting scheme. [Group page](./mhd/).

| Page | IDs / legs | Type | Figure |
|---|---|---|---|
| [MHD plumbing](./mhd/plumbing) | — (`mhd/plumbing`) | C | no |
| [MHD with a zero field](./mhd/zero-field) | MV-2 | C | no |
| [Ryu–Jones 2a: all seven waves](./mhd/ryu-jones-2a) | MV-4 | P | yes |
| [Brio–Wu and Ryu–Jones 4d](./mhd/brio-wu-rj4d) | — (`mhd/riemann`) | P | yes |
| [GLM divergence cleaning](./mhd/glm-cleaning) | MV-3; MV-10 (M6) | P | yes |
| [Order of accuracy: linear waves, CPAW, magnetised vortex](./mhd/order-of-accuracy) | MV-5, MV-6, MV-7 | P | yes |
| [Orszag–Tang vortex](./mhd/orszag-tang) | MV-12 | P | yes |
| [MHD rotor](./mhd/rotor) | MV-13 | P | yes |
| [Field loop advection](./mhd/field-loop) | MV-9 | P | yes |
| [Rotated shock tube](./mhd/rotated-shock-tube) | MV-8 | P | yes |
| [MHD conservation across AMR seams](./mhd/conservation) | MV-11 | C | yes |

## #45 · Non-finite state guard

Issue [#45](https://github.com/szaghi/adam/issues/45). A run whose state turns non-finite stops at once, on every rank.
[Group page](./infrastructure/).

| Page | IDs / legs | Type | Figure |
|---|---|---|---|
| [Non-finite state guard](./infrastructure/nonfinite) | — (`nonfinite`) | C | no |

## M3 · Riemann fluxes, EGLM, positivity limiter (#47)

Issue [#47](https://github.com/szaghi/adam/issues/47). The `weno-riemann` scheme (LLF, HLL, HLLC, HLLD face fluxes with
a flux correction), the energy-consistent EGLM cleaning and the positivity limiter.

| Page | IDs / legs | Type | Figure |
|---|---|---|---|
| [Lax problem](./euler/lax) | RV-3 (`lax`) | P | yes |
| [Shu–Osher shock / entropy-wave interaction](./euler/shu-osher) | RV-3 (`shu-osher`) | P | yes |
| [MHD on `weno-riemann`](./mhd/weno-riemann) | RV-5, RV-6, RV-7 | P | yes |
| [EGLM cleaning](./mhd/eglm) | EV-2, EV-3, EV-4 | P | yes |
| [Positivity limiter](./mhd/positivity-limiter) | PV-1 to PV-5 | P | yes |

## #40, #51, #52, #54 · Inter-realm seams

Issues [#40](https://github.com/szaghi/adam/issues/40), [#51](https://github.com/szaghi/adam/issues/51),
[#52](https://github.com/szaghi/adam/issues/52), [#54](https://github.com/szaghi/adam/issues/54). Inter-realm seams
across ranks, between misaligned blocks and with a 2:1 jump; the tree lookup on quadtrees and single-block roots.
Verified on the forest page of #37.

| Page | IDs / legs | Type | Figure |
|---|---|---|---|
| [Forest seams: the inter-realm contract](./euler/multi-realm) | multirealm legs 3–5 | C | yes |

## #49 · Non-dimensional form

Issue [#49](https://github.com/szaghi/adam/issues/49). The scheme carries no absolute constant: rescaled inputs
reproduce the run bit for bit, and the `[reference]` layer reads dimensional input.

| Page | IDs / legs | Type | Figure |
|---|---|---|---|
| [Scaling covariance and the reference layer](./infrastructure/scaling-covariance) | NV-0, NV-4, NV-5, NV-6, NV-8, NV-9 | C | no |

## #46 · Quadtrees

Issue [#46](https://github.com/szaghi/adam/issues/46). 2:1 seams on quadtrees, per axis, against the octree reference;
blocks thinner than two ghost layers are refused ([#66](https://github.com/szaghi/adam/issues/66)).

| Page | IDs / legs | Type | Figure |
|---|---|---|---|
| [Quadtree AMR](./mhd/quadtree-amr) | MV-15 | C | yes |
| [Forward-facing step: a three-realm forest on quadtrees](./euler/forward-facing-step) | EV-step | P+C | yes |

## M4 · Dissipative terms (#65)

Issue [#65](https://github.com/szaghi/adam/issues/65). Viscosity, heat conduction and resistivity on Euler and MHD, and
the no-slip walls. [Group page](./dissipative/).

| Page | IDs / legs | Type | Figure |
|---|---|---|---|
| [Ghost probe: every face and edge ghost on a linear field](./infrastructure/ghost-probe) | GP | C | yes |
| [Input contract of the dissipative terms](./dissipative/input-contract) | DC | C | no |
| [Viscous waves](./dissipative/viscous-waves) | VV-1, VV-2 | P | yes |
| [Compressible Couette flow](./dissipative/couette) | VV-3 | P | yes |
| [Becker's viscous shock](./dissipative/becker-shock) | VV-4 | P | yes |
| [Ohmic decay and the visco-resistive Alfvén wave](./dissipative/ohmic-alfven) | VV-5 | P | yes |
| [Ohmic heating and the divergence constraint](./dissipative/ohmic-heating) | VV-6 | C | no |
| [Conservation across 2:1 seams and forest seams](./dissipative/seam-conservation) | VV-7 | C | yes |
| [Accuracy across 2:1 seams](./dissipative/seam-accuracy) | VV-8 | P | yes |
| [CPU and FNL agree run by run](./dissipative/cpu-fnl-agreement) | VV `agree` | C | no |

## M5 · Runtime AMR (#74)

Issue [#74](https://github.com/szaghi/adam/issues/74). Regrids during the time loop, on both backends, and their
accuracy against uniform-fine runs. [Group page](./runtime-amr/).

| Page | IDs / legs | Type | Figure |
|---|---|---|---|
| [Regrid round trip and input contract](./runtime-amr/regrid-round-trip) | RG (rg0 to rg2) | C | no |
| [Runtime regrids in FLUME](./runtime-amr/runtime-regrids) | RG (rg3 to rg6) | C | yes |
| [Accuracy of runtime AMR](./runtime-amr/accuracy) | AV (av2 to av8) | P | yes |
| [MHD blast symmetry on the flux-splitting path](./runtime-amr/mhd-blast-symmetry) | [#77](https://github.com/szaghi/adam/issues/77), av3c | P | yes |

## M6 · Seam accuracy and divergence (#68)

Issue [#68](https://github.com/szaghi/adam/issues/68). M6 measures the accuracy cap of static 2:1 seams, then lifts it
behind new keys: high-order restriction and fills, a sensor-gated reflux, an EGLM reflux that leaves $\mathbf{B}$ alone.
The forest page adds leg 6, the GLM page the MV-10 monitor. [Group page](./seams/).

| Page | IDs / legs | Type | Figure |
|---|---|---|---|
| [Order across static 2:1 seams](./seams/order-across-seams) | AO (A1, C1, D1, D2) | P | yes |
| [Divergence across seams](./seams/divergence-across-seams) | AD (A2) | P | yes |
| [Shock crossing a static 2:1 face](./seams/shock-crossing) | AS (A3, C3, C3c) | P | yes |
| [Integral drift and the WENO sensor near seams](./seams/integral-drift) | AL (A4) | P+C | yes |
| [Reflux policy: troubled against always](./seams/reflux-policy) | E1, E2 | P | yes |
| [EGLM reflux leaves B alone](./seams/eglm-reflux) | G3a, G3b | P | yes |
| [Seam design record (appendix)](./seams/design-record) | B1, B2, B3, E0, G2 | P | yes |

## M7 · External gravity (#73)

Issue [#73](https://github.com/szaghi/adam/issues/73). A constant gravitational acceleration: initial atmospheres, the
force, a fourth-order energy source, a well-balanced mode, Rayleigh–Taylor. [Group page](./gravity/).

| Page | IDs / legs | Type | Figure |
|---|---|---|---|
| [Stratified atmospheres and the `hydrostatic` boundary](./gravity/atmospheres) | GV (G0) | C | yes |
| [The gravity force](./gravity/force) | GV-1, GV-2, GV-3, GV-6, GV-7 | P+C | yes |
| [The order of the energy source](./gravity/energy-order) | GV-9 | P | yes |
| [The well-balanced mode](./gravity/well-balanced) | GV-2b, GV-5, GV-7, GV-8, GV-10 | P | yes |
| [Rayleigh–Taylor](./gravity/rayleigh-taylor) | GV-4 | P | yes |

## #82 · Ranks without blocks

Issue [#82](https://github.com/szaghi/adam/issues/82). A rank that owns no block runs every model and option, and the
result does not change.

| Page | IDs / legs | Type | Figure |
|---|---|---|---|
| [Ranks without blocks](./infrastructure/ranks-without-blocks) | DV-1 | C | yes |

## M8 · Self-gravity (#81)

Issue [#81](https://github.com/szaghi/adam/issues/81). The Poisson equation solved every Runge–Kutta stage, gravity
applied through the cell potential: periodic and isolated FFT solvers, multigrid on AMR grids. [Group
page](./self-gravity/).

| Page | IDs / legs | Type | Figure |
|---|---|---|---|
| [Gravity from a cell potential](./self-gravity/potential) | SV-1 | P | yes |
| [The FFT Poisson solver](./self-gravity/fft-poisson) | SV-2 | P | yes |
| [Jeans modes under self-gravity](./self-gravity/jeans) | SV-3, SV-4, SV-5 | P+C | yes |
| [The isolated Poisson solver](./self-gravity/isolated-solver) | SV-6 | P | yes |
| [The n = 1 polytrope: equilibrium and collapse](./self-gravity/polytrope) | SV-7, SV-8 | P | yes |
| [The multigrid Poisson solver on AMR grids](./self-gravity/multigrid-amr) | SV-9, SV-12, SV-13 | P | yes |
| [Jeans modes on a static AMR grid](./self-gravity/amr-jeans) | SV-10 | P | yes |
| [Self-gravity under runtime AMR](./self-gravity/runtime-amr) | SV-11 | C | yes |

## #87, #88 · Isolated self-gravity on AMR

Issues [#87](https://github.com/szaghi/adam/issues/87), [#88](https://github.com/szaghi/adam/issues/88). The star on
static AMR boxes and its collapse under runtime AMR: the conservative-state prolongation, the Jeans marker and, with M9,
dual energy.

| Page | IDs / legs | Type | Figure |
|---|---|---|---|
| [The polytrope on AMR: static boxes, runtime regrids, Jeans refinement](./self-gravity/star-amr) | SV-14 to SV-17, rg7 | P+C | yes |

## M9 · Dual energy (#91)

Issue [#91](https://github.com/szaghi/adam/issues/91). An advected entropy that floors the internal energy of cold cells
(`euler-de`).

| Page | IDs / legs | Type | Figure |
|---|---|---|---|
| [Dual energy: the advected entropy of euler-de](./euler/dual-energy) | DE (de1 to de4) | P+C | yes |

## M10 · Passive scalars (#92)

Issue [#92](https://github.com/szaghi/adam/issues/92). Up to four advected scalars that do not act on the flow.

| Page | IDs / legs | Type | Figure |
|---|---|---|---|
| [Passive scalars](./euler/passive-scalars) | PS (ps1 to ps4) | P+C | yes |

## M11 · Radiative cooling (#93)

Issue [#93](https://github.com/szaghi/adam/issues/93). Optically thin cooling by Townsend's exact integration.

| Page | IDs / legs | Type | Figure |
|---|---|---|---|
| [Radiative cooling](./euler/radiative-cooling) | CO (co1 to co7) | P | yes |

## M12 · Sink particles (#94)

Issue [#94](https://github.com/szaghi/adam/issues/94). Point masses that act through self-gravity: creation, accretion,
mergers. [Group page](./sinks/).

| Page | IDs / legs | Type | Figure |
|---|---|---|---|
| [Sink particles](./sinks/sink-particles) | SI (si1 to si11) | P+C | yes |

## M13 · Isothermal and barotropic equations of state (#100)

Issue [#100](https://github.com/szaghi/adam/issues/100). Gases without an energy equation, $P = P(\rho)$, on Euler and
MHD, with gravity, self-gravity and sinks.

| Page | IDs / legs | Type | Figure |
|---|---|---|---|
| [Isothermal equation of state](./euler/isothermal-eos) | EO (eo1 to eo7) | P | yes |
| [Barotropic equation of state](./euler/barotropic-eos) | EO (eb1 to eb8) | P | yes |
| [Barotropic gases under gravity](./euler/barotropic-gravity) | EO (eg1 to eg10) | P | yes |
| [Isothermal and barotropic MHD](./mhd/barotropic-mhd) | EO (em1 to em10) | P | yes |

## #110, #71 · Runtime AMR on a forest

Issues [#110](https://github.com/szaghi/adam/issues/110), [#71](https://github.com/szaghi/adam/issues/71). The realms of
a forest regrid across their seams. The inter-realm 2:1 seams on quadtrees of #71 are leg 15 of the
[forest page](./euler/multi-realm).

| Page | IDs / legs | Type | Figure |
|---|---|---|---|
| [Multi-realm runtime AMR](./runtime-amr/multi-realm) | multirealm legs 7–13 | C | yes |

## #119 · Turbulence

Issue [#119](https://github.com/szaghi/adam/issues/119). What the scheme does to turbulence, then the closures in order:
subgrid-scale models, a hybrid central/WENO scheme, wall models and a turbulent inflow. [Group page](./turbulence/).

| Page | IDs / legs | Type | Figure |
|---|---|---|---|
| [The ILES baseline (T0)](./turbulence/iles-baseline) | TT, HT: tt1–tt6, ht1–ht2 | P+C | yes |
| [Subgrid-scale models (T1)](./turbulence/sgs-models) | SG: sg1–sg4 | P+C | yes |
| [The hybrid central/WENO scheme (T2)](./turbulence/hybrid-scheme) | HY: hy1–hy13 | P+C | yes |
| [Plane channel and the wall model (T3)](./turbulence/channel) | WM: wm1–wm11 | P+C | yes |
| [Compressible boundary layer (T3.7)](./turbulence/boundary-layer) | TB: tb1–tb2 | P+C | yes |

## #121 · STREAmS roadmap

Issue [#121](https://github.com/szaghi/adam/issues/121). The inflow and boundary toolchain of a shock/boundary-layer
interaction, validated against the STREAmS DNS and the IUSTI experiment. [Group page](./streams/).

| Page | IDs / legs | Type | Figure |
|---|---|---|---|
| [Non-reflecting boundary (S2)](./streams/non-reflecting) | NR: nr1–nr4 | P | yes |
| [Shock generator (S3, S4a)](./streams/shock-generator) | SG: sg1–sg4 (boundary) | P | yes |
| [Shock/boundary-layer interaction (S4b, S4c)](./streams/sbli) | campaign, no leg | P | yes |
