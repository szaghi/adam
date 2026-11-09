---
layout: home

hero:
  name: FLUME
  text: Fluid Laboratory for Unsteady Multiphysics at Exascale
  tagline: Compressible Euler, Navier–Stokes and ideal MHD on ADAM's adaptive grids, with gravity, self-gravity, sinks and turbulence models, on CPUs and GPUs.
  image:
    light: /flume/logo.svg
    dark: /flume/logo-dark.svg
    alt: FLUME, a brushed black hole, a tilted disk, its lensed arcs and the photon ring, over an adaptive mesh
  actions:
    - theme: brand
      text: Overview
      link: /applications/flume/overview
    - theme: alt
      text: Tutorial
      link: /tutorials/flume/sod
    - theme: alt
      text: Input reference
      link: /applications/flume/input
    - theme: alt
      text: Verification
      link: /applications/flume/verification/

features:
  - icon: 💨
    title: Euler and Navier–Stokes
    details: Compressible Euler with optional dual energy and up to 4 passive scalars; viscosity (constant, power law or Sutherland), heat conduction, no-slip and isothermal walls.
    link: /applications/flume/models#compressible-euler
    linkText: Physical models
  - icon: 🧲
    title: Ideal MHD
    details: Cell-centred B with mixed GLM or energy-consistent EGLM cleaning, Ohmic resistivity, HLLD fluxes, and a cell-based positivity limiter for very low plasma β.
    link: /applications/flume/models#ideal-mhd
    linkText: MHD and cleaning
  - icon: 📐
    title: High-order space schemes
    details: WENO flux splitting (orders 3–9), WENO interpolation with LLF, HLL, HLLC or HLLD fluxes and a high-order correction, and a hybrid central/WENO scheme for turbulence.
    link: /applications/flume/numerics
    linkText: Numerical methods
  - icon: 🌡️
    title: Gases and cooling
    details: Ideal, isothermal and barotropic equations of state for Euler and MHD; optically thin radiative cooling.
    link: /applications/flume/models#isothermal-and-barotropic-gases
    linkText: Equations of state
  - icon: 🌌
    title: Gravity, self-gravity, sinks
    details: External gravity with a well-balanced mode; self-gravity, periodic or isolated, uniform or AMR, solved by FLAIL every stage; sink particles on the finest level.
    link: /applications/flume/models#external-gravity
    linkText: Gravity
  - icon: 🌀
    title: Turbulence and wall models
    details: Vreman, WALE and Sigma subgrid models, equilibrium and inverse-transformation wall models, digital-filter inflow, mass-flow forcing, non-reflecting faces and a shock generator.
    link: /applications/flume/turbulence
    linkText: Turbulence
  - icon: 🌲
    title: Static and runtime AMR
    details: Refinement at initialization and during the run (box, gradient, Löhner, Jeans and solid markers), 2:1 faces with conservative reflux, restarts across regrids.
    link: /applications/flume/numerics#adaptive-mesh-refinement
    linkText: AMR
  - icon: 🧩
    title: Forests and immersed solids
    details: Several realms glued by 1:1 or 2:1 seams in one run; static solids by an immersed boundary (Euler).
    link: /applications/flume/numerics#multi-realm-runs
    linkText: Multi-realm runs
  - icon: ⚡
    title: CPU and GPU
    details: The same INI input on the CPU backend (MPI + OpenMP) and on the FNL backend (MPI + OpenACC, or OpenMP offload); every milestone is verified on both.
    link: /developer/flume-source
    linkText: Source layout
---

## Start here

The [overview](./overview) lists what is implemented, with a [quick start](./overview#quick-start) and the
[known limitations](./overview#known-limitations). The [first run](/start/first-run) builds FLUME and checks a Sod
shock tube against its exact solution; the [tutorial](/tutorials/flume/sod) continues from there.

## See it run

<div class="adam-casts">
<figure class="adam-cast" style="--ar: 0.97">
<video class="light-only" autoplay loop muted playsinline poster="/gallery/orszag-tang-light.png" aria-label="Orszag–Tang vortex, density"><source src="/gallery/orszag-tang-light.webm" type="video/webm"><source src="/gallery/orszag-tang-light.mp4" type="video/mp4"></video>
<video class="dark-only" autoplay loop muted playsinline poster="/gallery/orszag-tang-dark.png" aria-label="Orszag–Tang vortex, density"><source src="/gallery/orszag-tang-dark.webm" type="video/webm"><source src="/gallery/orszag-tang-dark.mp4" type="video/mp4"></video>
<figcaption>Orszag–Tang vortex, ideal MHD with GLM cleaning: density to t = 0.5 on 128², the blocks outlined. Case MV-12.</figcaption>
</figure>
<figure class="adam-cast" style="--ar: 1.89">
<video class="light-only" autoplay loop muted playsinline poster="/gallery/rayleigh-taylor-light.png" aria-label="Rayleigh–Taylor instability, hydro and MHD, density"><source src="/gallery/rayleigh-taylor-light.webm" type="video/webm"><source src="/gallery/rayleigh-taylor-light.mp4" type="video/mp4"></video>
<video class="dark-only" autoplay loop muted playsinline poster="/gallery/rayleigh-taylor-dark.png" aria-label="Rayleigh–Taylor instability, hydro and MHD, density"><source src="/gallery/rayleigh-taylor-dark.webm" type="video/webm"><source src="/gallery/rayleigh-taylor-dark.mp4" type="video/mp4"></video>
<figcaption>Rayleigh–Taylor single mode under constant gravity, density to t = 5 on 128²: hydrodynamics (left) and MHD with a horizontal field at half the critical value of the mode (right). Case GV-4.</figcaption>
</figure>
<figure class="adam-cast adam-cast-wide" style="--ar: 2.75">
<video class="light-only" autoplay loop muted playsinline poster="/gallery/forward-facing-step-light.png" aria-label="Mach 3 forward-facing step on a three-realm forest, density"><source src="/gallery/forward-facing-step-light.webm" type="video/webm"><source src="/gallery/forward-facing-step-light.mp4" type="video/mp4"></video>
<video class="dark-only" autoplay loop muted playsinline poster="/gallery/forward-facing-step-dark.png" aria-label="Mach 3 forward-facing step on a three-realm forest, density"><source src="/gallery/forward-facing-step-dark.webm" type="video/webm"><source src="/gallery/forward-facing-step-dark.mp4" type="video/mp4"></video>
<figcaption>Woodward–Colella Mach 3 step, density to t = 4 on a forest of three quadtree realms, base spacing 1/80 with static refined boxes; the inter-realm seams dashed. Case EV-step.</figcaption>
</figure>
</div>

<div class="adam-gallery">

[![Shock/boundary-layer interaction at Mach 2.28](/flume/sbli-fields.png)<strong>Shock/boundary-layer interaction</strong><span>Mach 2.28, wall-modelled LES: schlieren, pressure and velocity planes</span>](./verification/streams/sbli)

[![Mach 2.5 turbulent boundary layer and the interaction in 3-D](/flume/tbl-3d.png)<strong>Turbulent boundary layer</strong><span>Mach 2.5 from a digital-filter inflow: skin friction 0.2 % below the DNS at δ_i/16</span>](./verification/turbulence/boundary-layer)

[![Vorticity of the Taylor-Green vortex and of isotropic turbulence](/flume/tgv-fields.png)<strong>Taylor–Green and isotropic turbulence</strong><span>vorticity, hybrid-kep-weno at 128³</span>](./verification/turbulence/hybrid-scheme)

[![Streamwise velocity in two wall-modelled channels](/flume/channel-fields.png)<strong>Wall-modelled channel</strong><span>Mach 1.5 and 3: Re_τ 219–222 and 443–445, DNS 220–222 and 451</span>](./verification/turbulence/channel)

[![Balsara-Spicer blast with the positivity limiter](/flume/blast.png)<strong>Balsara–Spicer blast</strong><span>β = 2.5·10⁻⁴ with the positivity limiter, both space schemes</span>](./verification/mhd/positivity-limiter)

[![MHD rotor, density and magnetic pressure](/flume/rotor.png)<strong>MHD rotor</strong><span>density and magnetic pressure, symmetric to 1.8·10⁻⁹</span>](./verification/mhd/rotor)

[![Shock over a cylinder imposed by the immersed boundary](/flume/shock-cylinder.png)<strong>Shock over a cylinder</strong><span>Mach 2 shock on an immersed-boundary cylinder, AMR</span>](./verification/euler/shock-cylinder)

[![A self-gravitating star on a static AMR box](/flume/selfgravity-star-amr.png)<strong>Self-gravitating star on AMR</strong><span>the refined box halves the core's spurious swing</span>](./verification/self-gravity/star-amr)

[![Sink particles: binary, accretion rate, collapse on runtime AMR](/flume/sink-particles.png)<strong>Sink particles</strong><span>a binary orbit, the Shu accretion rate, a collapse on runtime AMR</span>](./verification/sinks/sink-particles)

[![Woodward-Colella step at t = 4](/flume/step.png)<strong>Forward-facing step</strong><span>three quadtree realms, refined run against the uniform 1/80 control</span>](./verification/euler/forward-facing-step)

</div>

The animations are made by `docs/scripts/make_gallery_runs.sh` (committed verification inputs, only the checkpoint
period changed) and `docs/scripts/make_gallery.py`; the figures by `src/tests/flume/verification/make_doc_figures.py`.

## What is verified

Each capability links to its section of the [verification catalog](./verification/), where every case has an oracle
run on both backends.

| Capability | Content |
|---|---|
| [Compressible Euler](./verification/#compressible-euler) | Riemann problems, design order, immersed boundary, reflux, restart; `weno-riemann`; isothermal and barotropic equations of state, dual energy (`euler-de`), passive scalars, optically thin radiative cooling (M1, M3, M9–M11, M13) |
| [Ideal MHD](./verification/#ideal-mhd) | the MHD waves and classic problems, mixed GLM and EGLM cleaning, `weno-riemann`, the positivity limiter, quadtrees (M2, M3, M13) |
| [Viscosity, heat conduction and resistivity](./verification/#viscosity-heat-conduction-and-resistivity) | viscosity, heat conduction, Ohmic resistivity, no-slip walls, conservative across seams (M4) |
| [External gravity](./verification/#external-gravity) | the force and its energy source, stratified atmospheres, the well-balanced mode, Rayleigh–Taylor (M7) |
| [Self-gravity](./verification/#self-gravity) | periodic and isolated, uniform and AMR, the Jeans refinement marker (M8, M9) |
| [Sink particles](./verification/#sink-particles) | creation, accretion, mergers, restarts, runtime AMR (M12) |
| [Accuracy across 2:1 seams](./verification/#accuracy-across-2-1-seams) | order, divergence and shocks at 2:1 faces, the operators that lift the cap (M6) |
| [Runtime AMR](./verification/#runtime-amr) | regrids during the run, restarts across regrids, forests that regrid across their seams (M5) |
| [Turbulence](./verification/#turbulence) (open) | subgrid-scale models, `hybrid-kep-weno`, wall models, digital-filter inflow, mass-flow forcing ([#119](https://github.com/szaghi/adam/issues/119)) |
| [Supersonic boundary layers](./verification/#supersonic-boundary-layers-boundaries-and-shock-interaction) (open) | Sutherland viscosity, non-reflecting faces, shock generator, wall-modelled shock/boundary-layer interaction ([#121](https://github.com/szaghi/adam/issues/121)) |
| [Infrastructure](./verification/#infrastructure) | ghost cells, the non-finite guard, ranks without blocks, the non-dimensional form |
