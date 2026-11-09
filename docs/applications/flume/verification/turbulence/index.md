# Turbulence

The turbulence verification of FLUME measures what the scheme does to turbulence before anything is built on it, then
adds the closures in order: subgrid-scale models, a hybrid central/WENO scheme, and wall models with a turbulent inflow.
The contract legs (C) are exact integrals, rank and backend invariance, restarts and laminar exact solutions; the
physics (P) is the comparison against DNS. All legs are in `src/tests/flume/verification/turbulence/check.sh`
(`--legs`), CPU and FNL. The methods and keys are on the [turbulence](../../turbulence) page and in the
[input](../../input) reference. The shock/boundary-layer interaction that continues the
[compressible boundary layer](./boundary-layer) is in [supersonic boundary layers](../streams/).

| Case | Type | What it verifies |
|---|---|---|
| [The ILES baseline](./iles-baseline) | P+C | what the scheme alone does to the Taylor-Green vortex (Re 1600, M 0.1, against the 512³ DNS; $M_0 = 1.25$), on a static 2:1 box and on decaying compressible isotropic turbulence: numerical dissipation dominates, 77 % and 58 % of the energy lost up to the dissipation peak at $64^3$ and $128^3$ with WENO-5 at M 0.1 |
| [Subgrid-scale models](./sgs-models) | P+C | the algebraic models (Vreman, WALE, Sigma) on WENO-5: no model beats no model, each lowers the enstrophy peak further from the DNS |
| [The hybrid central/WENO scheme](./hybrid-scheme) | P+C | a Ducros-type sensor, a 6th-order kinetic-energy-preserving central flux, one flux per face, the skin at 2:1 seams and physical faces, pressure equilibrium: at $128^3$ closer to the DNS than WENO-5 at $256^3$ on every L2 error; a model makes it worse |
| [Plane channel and the wall model](./channel) | P+C | mass-flow forcing, isothermal and `wall-modelled` walls, the equilibrium wall model with semi-local damping, against the Coleman channel: with $h = 0.15$, $Re_\tau$ 219–222 at Mach 1.5 (DNS 220–222) and 443–445 at Mach 3 (DNS 451, with the eddy viscosity below $h$) |
| [Compressible boundary layer](./boundary-layer) | P+C | a digital-filter inflow at Mach 2.5 and 5.86 against the DNS: $C_f$ 0.2 % below at Mach 2.5 ($\delta_i/16$), 3.8 % below at Mach 5.86 with $B_q$ as the DNS; $H$ the DNS's against the edge state |

Each page lists what it leaves open under "Not done".

## How it was built

The turbulence front is [#119](https://github.com/szaghi/adam/issues/119), in phases T0 (the ILES baseline), T1 (the
subgrid-scale models), T2 (the hybrid scheme) and T3 (the wall model and the boundary layer); the shock/boundary-layer
interaction continues T3.7 under [#121](https://github.com/szaghi/adam/issues/121).

| Page | Phases | Legs | Type |
|---|---|---|---|
| [The ILES baseline](./iles-baseline) | T0.1–T0.6 | tt1–tt6, ht1–ht2 | P+C |
| [Subgrid-scale models](./sgs-models) | T1 | sg1–sg4 | P+C |
| [The hybrid central/WENO scheme](./hybrid-scheme) | T2.1–T2.7 | hy1–hy13 | P+C |
| [Plane channel and the wall model](./channel) | T3.1–T3.6 | wm1–wm11 | P+C |
| [Compressible boundary layer](./boundary-layer) | T3.7 | tb1–tb2 | P+C |

| Phase | What | Legs | Page |
|---|---|---|---|
| T0.1 | Diagnostics: volume integrals of the kinetic-energy budget, spectra | tt1, tt2, tt4, tt5, tt6 | [ILES baseline](./iles-baseline) |
| T0.2 | `taylor-green` initial condition | tt1 | [ILES baseline](./iles-baseline) |
| T0.3 | Taylor-Green vortex, Re 1600, M 0.1, against the 512³ DNS | tt3 | [ILES baseline](./iles-baseline#taylor-green-vortex) |
| T0.4 | Taylor-Green vortex, $M_0 = 1.25$ | — | [ILES baseline](./iles-baseline#taylor-green-vortex) |
| T0.5 | Decaying compressible isotropic turbulence | ht1, ht2 | [ILES baseline](./iles-baseline#decaying-compressible-isotropic-turbulence) |
| T0.6 | A static 2:1 box on the Taylor-Green vortex | — | [ILES baseline](./iles-baseline#a-static-2-1-box) |
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
