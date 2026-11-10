# Verification

Every FLUME capability is backed by a verification case in `src/tests/flume/verification/`: a `check.sh` that runs the
case and an oracle that asserts the physics (an exact solution, a reference, a symmetry or a conservation property),
not a stored golden. The goldens are the job of the [regression suite](/reference/regression#flume-suite); the building
blocks below the solver are pinned by the [unit tests](./unit-tests).

The catalog below is grouped by capability. Each card names what its page verifies; each page has the same sections:
case and scope, models and numerics with references to the code, initial and boundary conditions, input files, results
on CPU and FNL, findings. All runs use 2 MPI ranks and WENO5 unless a page says otherwise. The tag of a card says what
it verifies:

- **physics**: a run against an exact or reference solution, an order of accuracy or a symmetry;
- **contract**: a bitwise identity, rank or backend invariance (CPU = FNL), a restart, conservation to round-off, an
  input refusal. A contract page documents its setup in a table; most have no figure, since a plot of a difference that
  is zero shows nothing.

## Compressible Euler

Riemann problems against their exact solutions, the design order on smooth flows, the immersed boundary, the equations of state (ideal, isothermal, barotropic), dual energy, passive scalars and radiative cooling; conservation across 2:1 faces and forest seams, restarts. [Group page](./euler/).

<div class="adam-gallery adam-catalog">

[![Sod shock tube](/flume/sod.png)<strong>Sod shock tube</strong><span>$L_1(\rho)$ against the exact solution; $x$, $y$, $z$ bitwise equal; reflecting walls</span><em>physics</em>](./euler/sod)

[![Isentropic vortex: order of accuracy](/flume/vortex.png)<strong>Isentropic vortex: order of accuracy</strong><span>the design order on a smooth exact solution, periodic blocks across ranks</span><em>physics</em>](./euler/isentropic-vortex)

[![Shock over a cylinder](/flume/shock-cylinder.png)<strong>Shock over a cylinder</strong><span>the immersed boundary with solid AMR: refinement, mirror symmetry, positivity</span><em>physics</em>](./euler/shock-cylinder)

[![Conservation across AMR seams](/flume/conservation.png)<strong>Conservation across AMR seams</strong><span>the reflux at 2:1 faces: the volume integrals constant to round-off</span><em>contract</em>](./euler/amr-conservation)

[<strong>Restart round trip and output channels</strong><span>a restarted run equals the uninterrupted one bit for bit on four integration paths, an empty rank included; slices and auxiliary fields</span><em>contract</em>](./euler/restart)

[![Forest seams: the inter-realm contract](/flume/multirealm.png)<strong>Forest seams: the inter-realm contract</strong><span>mirror, misaligned and 2:1 inter-realm seams, across ranks: bitwise equal to one realm</span><em>contract</em>](./euler/multi-realm)

[![Lax problem](/flume/lax.png)<strong>Lax problem</strong><span>the Riemann solvers on a strong contact, bound $1.05\times$ the splitting scheme</span><em>physics</em>](./euler/lax)

[![Shu–Osher shock / entropy-wave interaction](/flume/shu-osher.png)<strong>Shu–Osher shock / entropy-wave interaction</strong><span>the entropy waves behind a Mach 3 shock against a 16× finer reference</span><em>physics</em>](./euler/shu-osher)

[![Forward-facing step](/flume/step-trees.png)<strong>Forward-facing step</strong><span>a three-realm forest on quadtrees against octrees; the Woodward–Colella step at $t = 4$</span><em>physics + contract</em>](./euler/forward-facing-step)

[![Dual energy: the advected entropy of euler-de](/flume/dual-energy-noh.png)<strong>Dual energy: the advected entropy of euler-de</strong><span>`euler-de` changes nothing where no cell is cold (Sod bitwise) and keeps a cold Noh shock heated</span><em>physics + contract</em>](./euler/dual-energy)

[![Passive scalars](/flume/passive-scalars.png)<strong>Passive scalars</strong><span>up to four scalars that do not act on the flow, are conserved, and stay uniform through runtime AMR</span><em>physics + contract</em>](./euler/passive-scalars)

[![Radiative cooling](/flume/radiative-cooling.png)<strong>Radiative cooling</strong><span>optically thin cooling by Townsend's exact integration: the isochoric ODE to round-off, the floor, a radiative shock against its steady structure</span><em>physics</em>](./euler/radiative-cooling)

[![Isothermal equation of state](/flume/isothermal-eos.png)<strong>Isothermal equation of state</strong><span>$P = c_0^2\rho$: the shock tube and a Mach-5 collision against the exact solution, the fifth order on an acoustic wave, the inert energy slot, the contract</span><em>physics</em>](./euler/isothermal-eos)

[![Barotropic equation of state](/flume/barotropic-eos.png)<strong>Barotropic equation of state</strong><span>a piecewise polytrope: a fan split at a convex break, shocks across segments, the order with $c^2 = P'(\rho_0)$, the Mathew & Federrath table by self-convergence</span><em>physics</em>](./euler/barotropic-eos)

[![Barotropic gases under gravity](/flume/barotropic-gravity.png)<strong>Barotropic gases under gravity</strong><span>gravity on the momentum only, Jeans waves with $c^2 = P'(\rho_0)$, the Jeans marker, sinks with $c_s^2 = dP/d\rho$ and the virial thermal energy, Bonnor–Ebert stability, Boss–Bodenheimer</span><em>physics</em>](./euler/barotropic-gravity)

</div>

## Ideal MHD

The seven MHD waves and the classic MHD problems against exact or reference solutions, the order of accuracy on linear waves, the divergence control (GLM, EGLM), the Riemann-flux path and the positivity limiter, conservation on AMR and quadtrees. [Group page](./mhd/).

<div class="adam-gallery adam-catalog">

[<strong>MHD plumbing</strong><span>the MHD model end to end: the fast-speed and $c_h$ time step, the auxiliary fields, a bitwise restart, the refused inputs</span><em>contract</em>](./mhd/plumbing)

[<strong>MHD with a zero field</strong><span>MHD with $\mathbf{B} = 0$ reduces to Euler (Sod along x/y/z, bitwise across axes); a uniform MHD state is steady, bitwise on one level and to round-off across a 2:1 patch</span><em>contract</em>](./mhd/zero-field)

[![Ryu–Jones 2a: all seven waves](/flume/rj2a.png)<strong>Ryu–Jones 2a: all seven waves</strong><span>the whole MHD path against an exact 7-state solution; x/y/z bitwise; GLM inert in 1-D</span><em>physics</em>](./mhd/ryu-jones-2a)

[![Brio–Wu and Ryu–Jones 4d](/flume/brio-wu.png)<strong>Brio–Wu and Ryu–Jones 4d</strong><span>clean runs near the eigensystem degeneracies; the positivity floors engage and stop as they should</span><em>physics</em>](./mhd/brio-wu-rj4d)

[![GLM divergence cleaning](/flume/glm-pulse.png)<strong>GLM divergence cleaning</strong><span>the GLM telegraph system against its exact solution; div(B) decay of the Dedner peak, uniform and across seams; the `divb_tol` monitor</span><em>physics</em>](./mhd/glm-cleaning)

[![Order of accuracy](/flume/linear-wave-order.png)<strong>Order of accuracy</strong><span>fifth order on linear waves, the circularly polarised Alfvén wave and the magnetised vortex</span><em>physics</em>](./mhd/order-of-accuracy)

[![Orszag–Tang vortex](/flume/orszag-tang.png)<strong>Orszag–Tang vortex</strong><span>180° symmetry, conservation and positivity, uniform and across 2:1 seams</span><em>physics</em>](./mhd/orszag-tang)

[![MHD rotor](/flume/rotor.png)<strong>MHD rotor</strong><span>symmetry under rotation with $\mathbf{B} \to -\mathbf{B}$, positivity</span><em>physics</em>](./mhd/rotor)

[![Field loop advection](/flume/field-loop.png)<strong>Field loop advection</strong><span>the divergence error GLM leaves, the loop dissipation, CPU–FNL agreement on $\mathbf{B}$</span><em>physics</em>](./mhd/field-loop)

[![Rotated shock tube](/flume/rotated-shock-tube.png)<strong>Rotated shock tube</strong><span>$B_\parallel$ preservation and agreement with the 1-D problem at 63.4° and 45°</span><em>physics</em>](./mhd/rotated-shock-tube)

[![MHD conservation across AMR seams](/flume/mhd-conservation.png)<strong>MHD conservation across AMR seams</strong><span>every MHD integral, $\mathbf{B}$ and $\psi$ included, constant to round-off across 2:1 faces with reflux; the no-reflux negative control drifts</span><em>contract</em>](./mhd/conservation)

[![MHD on `weno-riemann`](/flume/weno-riemann.png)<strong>MHD on `weno-riemann`</strong><span>the HLLD scheme on the M2 verifications</span><em>physics</em>](./mhd/weno-riemann)

[![EGLM cleaning](/flume/eglm-energy.png)<strong>EGLM cleaning</strong><span>the energy-consistent cleaning: telegraph solution, EGLM = GLM at $\psi = 0$, drifts bounded by the sources</span><em>physics</em>](./mhd/eglm)

[![Positivity limiter](/flume/blast.png)<strong>Positivity limiter</strong><span>positivity on magnetised blasts and near vacuum, inactivity on smooth flow, 2:1 seams</span><em>physics</em>](./mhd/positivity-limiter)

[![Quadtree AMR](/flume/quadtree.png)<strong>Quadtree AMR</strong><span>2:1 seams on quadtrees against the octree reference</span><em>contract</em>](./mhd/quadtree-amr)

[![Isothermal and barotropic MHD](/flume/barotropic-mhd.png)<strong>Isothermal and barotropic MHD</strong><span>ideal MHD without an energy equation: Kim's eigensystem and Mignone's isothermal HLLD on $P(\rho)$, fifth-order linear waves, Mignone's tubes, the inert energy slot</span><em>physics</em>](./mhd/barotropic-mhd)

</div>

## Viscosity, heat conduction and resistivity

Viscous and resistive waves against their exact decay, Becker's shock and Couette flow against exact profiles, Ohmic heating and the divergence constraint, accuracy and conservation across 2:1 seams, and the input contract. [Group page](./dissipative/).

<div class="adam-gallery adam-catalog">

[<strong>Input contract</strong><span>34 input cases: invalid combinations refused with a message, implied coefficients logged exactly, zero coefficients bitwise inert, dimensional keys converted back exactly</span><em>contract</em>](./dissipative/input-contract)

[![Viscous waves](/flume/viscous.png)<strong>Viscous waves</strong><span>order 4 and 2 of the viscous flux on a shear wave (along $x$, $y$ and the diagonal), $\mu$ and $k$ together on an acoustic wave, against the exact linear mode</span><em>physics</em>](./dissipative/viscous-waves)

[![Compressible Couette flow](/flume/couette.png)<strong>Compressible Couette flow</strong><span>the no-slip walls (isothermal and moving) and the viscous heating, against the exact steady profile</span><em>physics</em>](./dissipative/couette)

[![Becker's viscous shock](/flume/viscous.png)<strong>Becker's viscous shock</strong><span>the nonlinear viscous and heat fluxes through a resolved shock, Mach 2 and 3</span><em>physics</em>](./dissipative/becker-shock)

[![Ohmic decay and the visco-resistive Alfvén wave](/flume/viscous.png)<strong>Ohmic decay and the visco-resistive Alfvén wave</strong><span>the Ohmic flux (sign, cross terms, coefficient), $\mu$ and $\eta$ together, the viscous fluxes on MHD</span><em>physics</em>](./dissipative/ohmic-alfven)

[<strong>Ohmic heating and the divergence constraint</strong><span>the Ohmic energy returned to the gas as heat; no divergence added by the curl form (GLM, EGLM)</span><em>contract</em>](./dissipative/ohmic-heating)

[![Conservation across 2:1 seams and forest seams](/flume/seam-conservation.png)<strong>Conservation across 2:1 seams and forest seams</strong><span>the dissipative fluxes carried by the reflux and the inter-realm register; forest union equal to the single realm</span><em>contract</em>](./dissipative/seam-conservation)

[![Accuracy across 2:1 seams](/flume/viscous-fields.png)<strong>Accuracy across 2:1 seams</strong><span>the order of the composite grid (second), and the reflux as the carrier of the dissipative seam flux</span><em>physics</em>](./dissipative/seam-accuracy)

[<strong>CPU and FNL agree run by run</strong><span>every one of the 72 VV runs equal on the two backends within $10^{-10}$ of the physical scales</span><em>contract</em>](./dissipative/cpu-fnl-agreement)

</div>

## External gravity

The gravity force and its energy source to design order, stratified atmospheres held at rest by the well-balanced mode and the hydrostatic boundary, the Rayleigh–Taylor growth rate against linear theory. [Group page](./gravity/).

<div class="adam-gallery adam-catalog">

[![Stratified atmospheres and the `hydrostatic` boundary](/flume/gravity-atmospheres.png)<strong>Stratified atmospheres and the `hydrostatic` boundary</strong><span>the initial atmospheres and the ghosts of the `hydrostatic` faces against the closed form, before any force acts</span><em>contract</em>](./gravity/atmospheres)

[![The gravity force](/flume/gravity-force.png)<strong>The gravity force</strong><span>exact free fall, mass and total energy of a closed box to round-off, the limiter near vacuum, scaling covariance</span><em>physics + contract</em>](./gravity/force)

[![The order of the energy source](/flume/gravity-order.png)<strong>The order of the energy source</strong><span>the fourth-order energy source restores the order of the scheme on a falling density wave</span><em>physics</em>](./gravity/energy-order)

[![The well-balanced mode](/flume/gravity-balance.png)<strong>The well-balanced mode</strong><span>an atmosphere at rest stays at rest to round-off, across 2:1 seams too, and small perturbations respond linearly</span><em>physics</em>](./gravity/well-balanced)

[![Rayleigh–Taylor](/flume/gravity-rt-linear.png)<strong>Rayleigh–Taylor</strong><span>the linear growth rate against the compressible linear theory; the nonlinear stage and the mirror symmetry of the MHD run</span><em>physics</em>](./gravity/rayleigh-taylor)

</div>

## Self-gravity

The Poisson solvers (FFT, isolated, multigrid on AMR) against exact potentials, Jeans modes against the dispersion relation, the exact exchange of kinetic and potential energy, and the $n=1$ polytrope in equilibrium and collapse on uniform and adaptive grids. [Group page](./self-gravity/).

<div class="adam-gallery adam-catalog">

[![Gravity from a cell potential](/flume/selfgravity-potential.png)<strong>Gravity from a cell potential</strong><span>the potential path reproduces the constant-g oracles of M7</span><em>physics</em>](./self-gravity/potential)

[![The FFT Poisson solver](/flume/selfgravity-poisson.png)<strong>The FFT Poisson solver</strong><span>the periodic solver is exact for $L_h$, second order against the continuum</span><em>physics</em>](./self-gravity/fft-poisson)

[![Jeans modes](/flume/selfgravity-jeans.png)<strong>Jeans modes</strong><span>the coupling: dispersion relation, growth rate, symmetry, exact $E + W$ exchange, scaling covariance</span><em>physics + contract</em>](./self-gravity/jeans)

[![The isolated solver](/flume/selfgravity-isolated.png)<strong>The isolated solver</strong><span>zero padding with the lattice Green's function: free-space potential of the mass on the grid</span><em>physics</em>](./self-gravity/isolated-solver)

[![The n = 1 polytrope](/flume/selfgravity-star.png)<strong>The n = 1 polytrope</strong><span>a star in equilibrium and a collapsing one under isolated self-gravity, uniform grid</span><em>physics</em>](./self-gravity/polytrope)

[![The multigrid solver on AMR grids](/flume/selfgravity-amr.png)<strong>The multigrid solver on AMR grids</strong><span>the mimetic 2:1 coupling (symmetric and consistent), the FCG solver, Dirichlet faces, isolated boundaries by the symmetric James step</span><em>physics</em>](./self-gravity/multigrid-amr)

[![Jeans modes on a static AMR grid](/flume/selfgravity-amr-jeans.png)<strong>Jeans modes on a static AMR grid</strong><span>the multigrid coupling, the seam gravity work, the per-stage reflux</span><em>physics</em>](./self-gravity/amr-jeans)

[![Self-gravity under runtime AMR](/flume/selfgravity-amr-regrid.png)<strong>Self-gravity under runtime AMR</strong><span>the solver and $\Phi$ follow every regrid; $E + W$ exact between regrids</span><em>contract</em>](./self-gravity/runtime-amr)

[![The polytrope on AMR](/flume/selfgravity-star-amr.png)<strong>The polytrope on AMR</strong><span>the star on a static AMR box; its collapse under runtime AMR with the conservative-state prolongation, a Jeans marker and dual energy; the Jeans-marker contract</span><em>physics + contract</em>](./self-gravity/star-amr)

</div>

## Sink particles

Point masses coupled through the Poisson source: rest and Kepler orbits, creation, accretion at Shu's rate, mergers, restarts, CPU against FNL and runtime AMR. [Group page](./sinks/).

<div class="adam-gallery adam-catalog">

[![Sink particles](/flume/sink-particles.png)<strong>Sink particles</strong><span>point masses through the Poisson source: rest, Kepler binary, the input contract, restart, creation, accretion at Shu's rate, mergers, CPU against FNL, runtime AMR</span><em>physics + contract</em>](./sinks/sink-particles)

</div>

## Accuracy across 2:1 seams

What a static 2:1 face costs: the order of accuracy across it, the divergence the cleaning sees, a shock crossing it against the exact solution, the drift of the integrals, and the operators that lift the cap. [Group page](./seams/).

<div class="adam-gallery adam-catalog">

[![Order across static 2:1 seams](/flume/order-across-seams.png)<strong>Order across static 2:1 seams</strong><span>the order of a smooth vortex crossing a refined quadrant, per region, against the uniform control; the 6-point restriction; the quintic and limited coarse-to-fine fills</span><em>physics</em>](./seams/order-across-seams)

[![Divergence across seams](/flume/divergence-across-seams.png)<strong>Divergence across seams</strong><span>D_flux, the divergence the cleaning acts on: bounded and first-order convergent at the seams; #78 is not a divergence difference</span><em>physics</em>](./seams/divergence-across-seams)

[![Shock crossing a static 2:1 face](/flume/shock-crossing.png)<strong>Shock crossing a static 2:1 face</strong><span>mass, shock position and the plateau behind a shock crossing a 2:1 face, against the exact solution; the 6-point and WENO restrictions at shocks</span><em>physics</em>](./seams/shock-crossing)

[![Integral drift and the WENO sensor near seams](/flume/integral-drift.png)<strong>Integral drift and the WENO sensor near seams</strong><span>the drift of the integrals without the reflux, and what the WENO sensor reads at the seams</span><em>physics + contract</em>](./seams/integral-drift)

[![Reflux policy: troubled against always](/flume/reflux-policy.png)<strong>Reflux policy: troubled against always</strong><span>the sensor-gated reflux recovers the sixth order with shocks unchanged, at a conservation cost</span><em>physics</em>](./seams/reflux-policy)

[![EGLM reflux leaves B alone](/flume/field-loop-seams.png)<strong>EGLM reflux leaves B alone</strong><span>#78 is caused by the reflux of $\mathbf{B}$; without it EGLM's $B_z$ drops from 70 to 3.3 times the uniform level</span><em>physics</em>](./seams/eglm-reflux)

[![Seam design record](/flume/seam-proto-b1.png)<strong>Seam design record</strong><span>appendix: the 2-D prototypes, the restriction footprint, the capturing error behind a crossing shock, the GLM damping length</span><em>physics</em>](./seams/design-record)

</div>

## Runtime AMR

Regridding during the run: the regrid step and its conservative prolongation, conservation and bitwise restarts through regrids, tracked runs against uniform-fine controls, and forests whose realms regrid across their seams. [Group page](./runtime-amr/).

<div class="adam-gallery adam-catalog">

[<strong>Regrid round trip and input contract</strong><span>the library regrid in isolation: bookkeeping over 1, 2, 3 ranks, the conservative prolongation exact on the integral and in range, the `frequency` and `regrid_prolongation` contracts</span><em>contract</em>](./runtime-amr/regrid-round-trip)

[![Runtime regrids in FLUME](/flume/regrid-runtime.png)<strong>Runtime regrids in FLUME</strong><span>conservation through regrids, bitwise restart across regrids, the immersed solid, CPU against FNL</span><em>contract</em>](./runtime-amr/runtime-regrids)

[![Accuracy of runtime AMR](/flume/regrid.png)<strong>Accuracy of runtime AMR</strong><span>tracked runs against uniform-fine runs on the vortex, Sod, the MHD blast, the field loop and the cylinder; rank and backend invariance</span><em>physics</em>](./runtime-amr/accuracy)

[![MHD blast symmetry on the flux-splitting path](/flume/mhd-blast-symmetry.png)<strong>MHD blast symmetry on the flux-splitting path</strong><span>why the low-β blast loses its point symmetry on `weno` flux splitting, and the Riemann-flux path that keeps it</span><em>physics</em>](./runtime-amr/mhd-blast-symmetry)

[![Runtime AMR on a forest](/flume/multirealm-amr.png)<strong>Runtime AMR on a forest</strong><span>the realms of a forest regridding across their seams: the cross-realm balance, mixed levels along a seam, restart across regrids</span><em>contract</em>](./runtime-amr/multi-realm)

</div>

## Turbulence

What the scheme does to turbulence, then the closures: subgrid-scale models, a hybrid central/WENO scheme, and wall models with a turbulent inflow, against DNS of the Taylor–Green vortex, isotropic turbulence, channels and boundary layers. [Group page](./turbulence/).

<div class="adam-gallery adam-catalog">

[![The ILES baseline](/flume/tgv-iles.png)<strong>The ILES baseline</strong><span>what WENO-5 does to turbulence before any model: 77 % and 58 % of the energy lost to numerical dissipation up to the peak at $64^3$ and $128^3$; Taylor–Green and isotropic decay against DNS</span><em>physics + contract</em>](./turbulence/iles-baseline)

[![Subgrid-scale models](/flume/tgv-models.png)<strong>Subgrid-scale models</strong><span>Vreman, WALE and Sigma on WENO-5 against the DNS: no model beats no model, each lowers the enstrophy peak further from it</span><em>physics + contract</em>](./turbulence/sgs-models)

[![The hybrid central/WENO scheme](/flume/tgv-hybrid.png)<strong>The hybrid central/WENO scheme</strong><span>the central/WENO hybrid at $128^3$ closer to the DNS than WENO-5 at $256^3$ on every L2 error; a subgrid model makes it worse</span><em>physics + contract</em>](./turbulence/hybrid-scheme)

[![Plane channel and the wall model](/flume/channel-wm.png)<strong>Plane channel and the wall model</strong><span>the equilibrium wall model in a compressible channel: $Re_\tau$ 219–222 at Mach 1.5 (DNS 220–222) and 443–445 at Mach 3 (DNS 451)</span><em>physics + contract</em>](./turbulence/channel)

[![Compressible boundary layer](/flume/tbl-m2p5.png)<strong>Compressible boundary layer</strong><span>Mach 2.5 and 5.86 layers from a digital-filter inflow: $C_f$ 0.2 % below the DNS at $\delta_i/16$, 3.8 % below at Mach 5.86 with $B_q$ as the DNS</span><em>physics + contract</em>](./turbulence/boundary-layer)

</div>

## Supersonic boundary layers: boundaries and shock interaction

The boundaries a shock/boundary-layer interaction needs, a non-reflecting outflow and an oblique-shock generator, verified against exact relations, and the wall-modelled interaction against DNS and experiment. [Group page](./streams/).

<div class="adam-gallery adam-catalog">

[![Non-reflecting boundary](/flume/nrbc.png)<strong>Non-reflecting boundary</strong><span>the characteristic outflow and far field: a leaving packet reflects $2.9\cdot10^{-3}$ of its amplitude, a 2 % overpressure drains to 0.099 (`extrapolation` keeps 1.000)</span><em>physics</em>](./streams/non-reflecting)

[![Shock generator](/flume/shockgen.png)<strong>Shock generator</strong><span>an oblique shock through a $y$ face against Rankine–Hugoniot ($p$, $\rho$ within $3\cdot10^{-3}$, angle within $0.05^\circ$) and its regular reflection off a slip wall</span><em>physics</em>](./streams/shock-generator)

[![Shock/boundary-layer interaction](/flume/sbli.png)<strong>Shock/boundary-layer interaction</strong><span>the Mach 2.28, $8^\circ$ interaction, wall-modelled, against the STREAmS DNS and the IUSTI experiment: $L_{int}/\delta_0$ 3.15 against 3.30; the DNS's separation ($2.2\,\delta_0$) missing, the compression half an interaction length late</span><em>physics</em>](./streams/sbli)

</div>

## Infrastructure

Contracts below the physics: every ghost cell on a linear field, the non-finite state guard, ranks that own no block, and the covariance of the non-dimensional form. [Group page](./infrastructure/).

<div class="adam-gallery adam-catalog">

[<strong>Non-finite state guard</strong><span>a NaN planted in the Sod state stops the run at step 1 on the guard's message, on every rank, with the first non-finite value located</span><em>contract</em>](./infrastructure/nonfinite)

[<strong>Scaling covariance and the reference layer</strong><span>rescaled inputs reproduce the base run bit for bit; dimensional input, restart and output through `[reference]`</span><em>contract</em>](./infrastructure/scaling-covariance)

[![Ghost probe](/flume/ghosts.png)<strong>Ghost probe</strong><span>every face and edge ghost, on every fill path (copy, MPI, 2:1, inter-realm seams, physical boundaries), holds the value a stencil must read; a poison exposes reads before the fill</span><em>contract</em>](./infrastructure/ghost-probe)

[![Ranks without blocks](/flume/ranks.png)<strong>Ranks without blocks</strong><span>a rank that owns no block runs every model and option, and np 2 equals np 1 bit for bit</span><em>contract</em>](./infrastructure/ranks-without-blocks)

</div>

## How it was built

The capabilities were added one milestone or issue at a time; this table maps each to the pages that verify it, in the
order they were implemented.

| Built under | Capability | Pages |
|---|---|---|
| M1, [#35](https://github.com/szaghi/adam/issues/35) | Euler core | [Sod shock tube](./euler/sod), [Isentropic vortex: order of accuracy](./euler/isentropic-vortex), [Shock over a cylinder](./euler/shock-cylinder), [Conservation across AMR seams](./euler/amr-conservation), [Restart round trip and output channels](./euler/restart) |
| [#37](https://github.com/szaghi/adam/issues/37) | Forest entry | [Forest seams: the inter-realm contract](./euler/multi-realm) |
| M2, [#41](https://github.com/szaghi/adam/issues/41) | Ideal MHD | [MHD plumbing](./mhd/plumbing), [MHD with a zero field](./mhd/zero-field), [Ryu–Jones 2a: all seven waves](./mhd/ryu-jones-2a), [Brio–Wu and Ryu–Jones 4d](./mhd/brio-wu-rj4d), [GLM divergence cleaning](./mhd/glm-cleaning), [Order of accuracy: linear waves, CPAW, magnetised vortex](./mhd/order-of-accuracy), [Orszag–Tang vortex](./mhd/orszag-tang), [MHD rotor](./mhd/rotor), [Field loop advection](./mhd/field-loop), [Rotated shock tube](./mhd/rotated-shock-tube), [MHD conservation across AMR seams](./mhd/conservation) |
| [#45](https://github.com/szaghi/adam/issues/45) | Non-finite state guard | [Non-finite state guard](./infrastructure/nonfinite) |
| M3, [#47](https://github.com/szaghi/adam/issues/47) | Riemann fluxes, EGLM, positivity limiter | [Lax problem](./euler/lax), [Shu–Osher shock / entropy-wave interaction](./euler/shu-osher), [MHD on `weno-riemann`](./mhd/weno-riemann), [EGLM cleaning](./mhd/eglm), [Positivity limiter](./mhd/positivity-limiter) |
| [#40](https://github.com/szaghi/adam/issues/40), [#51](https://github.com/szaghi/adam/issues/51), [#52](https://github.com/szaghi/adam/issues/52), [#54](https://github.com/szaghi/adam/issues/54) | Inter-realm seams | [Forest seams: the inter-realm contract](./euler/multi-realm) |
| [#49](https://github.com/szaghi/adam/issues/49) | Non-dimensional form | [Scaling covariance and the reference layer](./infrastructure/scaling-covariance) |
| [#46](https://github.com/szaghi/adam/issues/46) | Quadtrees | [Quadtree AMR](./mhd/quadtree-amr), [Forward-facing step: a three-realm forest on quadtrees](./euler/forward-facing-step) |
| M4, [#65](https://github.com/szaghi/adam/issues/65) | Dissipative terms | [Ghost probe: every face and edge ghost on a linear field](./infrastructure/ghost-probe), [Input contract of the dissipative terms](./dissipative/input-contract), [Viscous waves](./dissipative/viscous-waves), [Compressible Couette flow](./dissipative/couette), [Becker's viscous shock](./dissipative/becker-shock), [Ohmic decay and the visco-resistive Alfvén wave](./dissipative/ohmic-alfven), [Ohmic heating and the divergence constraint](./dissipative/ohmic-heating), [Conservation across 2:1 seams and forest seams](./dissipative/seam-conservation), [Accuracy across 2:1 seams](./dissipative/seam-accuracy), [CPU and FNL agree run by run](./dissipative/cpu-fnl-agreement) |
| M5, [#74](https://github.com/szaghi/adam/issues/74) | Runtime AMR | [Regrid round trip and input contract](./runtime-amr/regrid-round-trip), [Runtime regrids in FLUME](./runtime-amr/runtime-regrids), [Accuracy of runtime AMR](./runtime-amr/accuracy), [MHD blast symmetry on the flux-splitting path](./runtime-amr/mhd-blast-symmetry) |
| M6, [#68](https://github.com/szaghi/adam/issues/68) | Seam accuracy and divergence | [Order across static 2:1 seams](./seams/order-across-seams), [Divergence across seams](./seams/divergence-across-seams), [Shock crossing a static 2:1 face](./seams/shock-crossing), [Integral drift and the WENO sensor near seams](./seams/integral-drift), [Reflux policy: troubled against always](./seams/reflux-policy), [EGLM reflux leaves B alone](./seams/eglm-reflux), [Seam design record (appendix)](./seams/design-record) |
| M7, [#73](https://github.com/szaghi/adam/issues/73) | External gravity | [Stratified atmospheres and the `hydrostatic` boundary](./gravity/atmospheres), [The gravity force](./gravity/force), [The order of the energy source](./gravity/energy-order), [The well-balanced mode](./gravity/well-balanced), [Rayleigh–Taylor](./gravity/rayleigh-taylor) |
| [#82](https://github.com/szaghi/adam/issues/82) | Ranks without blocks | [Ranks without blocks](./infrastructure/ranks-without-blocks) |
| M8, [#81](https://github.com/szaghi/adam/issues/81) | Self-gravity | [Gravity from a cell potential](./self-gravity/potential), [The FFT Poisson solver](./self-gravity/fft-poisson), [Jeans modes under self-gravity](./self-gravity/jeans), [The isolated Poisson solver](./self-gravity/isolated-solver), [The n = 1 polytrope: equilibrium and collapse](./self-gravity/polytrope), [The multigrid Poisson solver on AMR grids](./self-gravity/multigrid-amr), [Jeans modes on a static AMR grid](./self-gravity/amr-jeans), [Self-gravity under runtime AMR](./self-gravity/runtime-amr) |
| [#87](https://github.com/szaghi/adam/issues/87), [#88](https://github.com/szaghi/adam/issues/88) | Isolated self-gravity on AMR | [The polytrope on AMR: static boxes, runtime regrids, Jeans refinement](./self-gravity/star-amr) |
| M9, [#91](https://github.com/szaghi/adam/issues/91) | Dual energy | [Dual energy: the advected entropy of euler-de](./euler/dual-energy) |
| M10, [#92](https://github.com/szaghi/adam/issues/92) | Passive scalars | [Passive scalars](./euler/passive-scalars) |
| M11, [#93](https://github.com/szaghi/adam/issues/93) | Radiative cooling | [Radiative cooling](./euler/radiative-cooling) |
| M12, [#94](https://github.com/szaghi/adam/issues/94) | Sink particles | [Sink particles](./sinks/sink-particles) |
| M13, [#100](https://github.com/szaghi/adam/issues/100) | Isothermal and barotropic equations of state | [Isothermal equation of state](./euler/isothermal-eos), [Barotropic equation of state](./euler/barotropic-eos), [Barotropic gases under gravity](./euler/barotropic-gravity), [Isothermal and barotropic MHD](./mhd/barotropic-mhd) |
| [#110](https://github.com/szaghi/adam/issues/110), [#71](https://github.com/szaghi/adam/issues/71) | Runtime AMR on a forest | [Runtime AMR on a forest](./runtime-amr/multi-realm) |
| [#119](https://github.com/szaghi/adam/issues/119) | Turbulence | [The ILES baseline (T0)](./turbulence/iles-baseline), [Subgrid-scale models (T1)](./turbulence/sgs-models), [The hybrid central/WENO scheme (T2)](./turbulence/hybrid-scheme), [Plane channel and the wall model (T3)](./turbulence/channel), [Compressible boundary layer (T3.7)](./turbulence/boundary-layer) |
| [#121](https://github.com/szaghi/adam/issues/121) | Boundaries and shock interaction for supersonic boundary layers (the STREAmS roadmap) | [Non-reflecting boundary](./streams/non-reflecting), [Shock generator and regular reflection](./streams/shock-generator), [Shock/boundary-layer interaction](./streams/sbli) |

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
