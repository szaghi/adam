# Turbulence (#119)

[#119](https://github.com/szaghi/adam/issues/119) measures what the scheme does to turbulence before anything is built on
it, then adds the closures in order: subgrid-scale models, a hybrid central/WENO scheme, and wall models with a turbulent
inflow. The methods and keys are on the [turbulence](../../turbulence) page and in the [input](../../input) reference.
The shock/boundary-layer interaction that continues T3.7 is in the [STREAmS roadmap](../streams/) (#121).

| Page | Phases | Legs | Type |
|---|---|---|---|
| [The ILES baseline](./iles-baseline) | T0.1–T0.6 | tt1–tt6, ht1–ht2 | P+C |
| [Subgrid-scale models](./sgs-models) | T1 | sg1–sg4 | P+C |
| [The hybrid central/WENO scheme](./hybrid-scheme) | T2.1–T2.7 | hy1–hy13 | P+C |
| [Plane channel and the wall model](./channel) | T3.1–T3.6 | wm1–wm11 | P+C |
| [Compressible boundary layer](./boundary-layer) | T3.7 | tb1–tb2 | P+C |

The contract legs (C) are exact integrals, rank and backend invariance, restarts and laminar exact solutions; the
physics (P) is the comparison against DNS. All legs are in `src/tests/flume/verification/turbulence/check.sh`
(`--legs`), CPU and FNL.

| Phase | What | Legs | Page |
|---|---|---|---|
| T0.1 | Diagnostics: volume integrals of the kinetic-energy budget, spectra | tt1, tt2, tt4, tt5, tt6 | [ILES baseline](./iles-baseline) |
| T0.2 | `taylor-green` initial condition | tt1 | [ILES baseline](./iles-baseline) |
| T0.3 | Taylor-Green vortex, Re 1600, M 0.1, against the 512³ DNS | tt3 | [ILES baseline](./iles-baseline#taylor-green-vortex-t0-3-t0-4) |
| T0.4 | Taylor-Green vortex, $M_0 = 1.25$ | — | [ILES baseline](./iles-baseline#taylor-green-vortex-t0-3-t0-4) |
| T0.5 | Decaying compressible isotropic turbulence | ht1, ht2 | [ILES baseline](./iles-baseline#decaying-compressible-isotropic-turbulence-t0-5) |
| T0.6 | A static 2:1 box on the Taylor-Green vortex | — | [ILES baseline](./iles-baseline#a-static-2-1-box-t0-6) |
| T1 | Algebraic subgrid-scale models (Vreman, WALE, Sigma) and the go/no-go | sg1–sg4 | [Subgrid-scale models](./sgs-models) |
| T2.1 | Ducros-type sensor | hy4 | [Hybrid scheme](./hybrid-scheme) |
| T2.2 | 6th-order kinetic-energy-preserving central flux | hy3 | [Hybrid scheme](./hybrid-scheme) |
| T2.3 | Blending, one flux per face | hy2 | [Hybrid scheme](./hybrid-scheme) |
| T2.4 | Skin at 2:1 seams and physical faces | hy5–hy7, hy10–hy12 | [Hybrid scheme](./hybrid-scheme#on-a-2-1-box-and-what-the-skin-costs) |
| T2.5 | Pressure equilibrium | hy8 | [Hybrid scheme](./hybrid-scheme) |
| T2.6 | The model only where the scheme is central | hy9 | [Hybrid scheme](./hybrid-scheme) |
| T2.7 | Validation: inviscid Taylor-Green, T0.3–T0.6 rerun, shock robustness | hy1, hy13 | [Hybrid scheme](./hybrid-scheme#against-the-dns) |
| T3.1 | Channel infrastructure: mass-flow forcing, isothermal walls | wm1, wm10 | [Channel](./channel) |
| T3.2 | Wall-resolved reference, Coleman Mach 1.5 | — | [Channel](./channel#mach-1-5) |
| T3.3 | `wall-modelled` boundary | wm2, wm4 | [Channel](./channel) |
| T3.4 | Equilibrium wall model | wm2, wm5 | [Channel](./channel) |
| T3.5 | Compressibility: semi-local damping, inverse-transformation model, input filter | wm6–wm9 | [Channel](./channel#mach-1-5) |
| T3.6 | 2:1 faces along a modelled wall | wm3 | [Channel](./channel) |
| T3.7 | Coleman channel at Mach 1.5 and 3; eddy viscosity below the matching height | wm11 | [Channel](./channel#mach-3) |
| T3.7 | Compressible boundary layer, Mach 2.5 and 5.86, digital-filter inflow | tb1, tb2 | [Boundary layer](./boundary-layer) |

## Results in one line each

- **T0.** Numerical dissipation dominates the ILES baseline: 77 % and 58 % of the energy lost up to the dissipation peak
  at $64^3$ and $128^3$ with WENO-5 at M 0.1.
- **T1.** No model beats no model on WENO-5: each lowers the enstrophy peak further from the DNS.
- **T2.** The hybrid scheme at $128^3$ is closer to the DNS than WENO-5 at $256^3$ on every L2 error; a model makes it worse.
- **T3.** With $h = 0.15$ and semi-local damping the wall model gives $Re_\tau$ 219–222 at Mach 1.5 (DNS 220–222) and
  443–445 at Mach 3 (DNS 451, with the eddy viscosity below $h$).
- **T3.7.** $C_f$ 0.2 % below the DNS at Mach 2.5 ($\delta_i/16$), 3.8 % below at Mach 5.86 with $B_q$ as the DNS;
  $H$ the DNS's against the edge state.

Each page lists what it leaves open under "Not done".
