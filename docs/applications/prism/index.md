---
layout: home

hero:
  name: PRISM
  text: Plasma Research usIng Simulation Methods
  tagline: Maxwell's equations and particle-in-cell plasmas on ADAM's adaptive grids, on CPUs and GPUs.
  image:
    light: /prism/logo.svg
    dark: /prism/logo-dark.svg
    alt: PRISM, a brushed helicon plasma thruster, four coils, the field lines and the nozzle, over an adaptive mesh
  actions:
    - theme: brand
      text: Overview
      link: /applications/prism/overview
    - theme: alt
      text: Tutorial
      link: /tutorials/prism/pulse
    - theme: alt
      text: Input reference
      link: /applications/prism/input
    - theme: alt
      text: Verification
      link: /applications/prism/verification/

features:
  - icon: 🔌
    title: Maxwell's equations
    details: The electric displacement D and the magnetic induction B, in SI or non-dimensional form, with hyperbolic (GLM) cleaning of D, B or both.
    link: /applications/prism/models
    linkText: Physical models
  - icon: 📐
    title: Centred operators
    details: Centred finite differences or finite volumes of order 2–10; Runge–Kutta in time with a CFL step on the speed of light.
    link: /applications/prism/numerics#space-discretization
    linkText: Numerical methods
  - icon: 🌲
    title: AMR and forests
    details: Octrees and quadtrees refined at initialization (box and current-variation markers), 2:1 faces with FV reflux, several realms glued with α or β seam cadence.
    link: /applications/prism/numerics#adaptive-mesh-refinement
    linkText: AMR and seams
  - icon: 🧲
    title: Coils and external fields
    details: Rectangular, circular, solenoid and helicon coils, DC or AC with a smooth ramp; external rotating or uniform fields acting on the particles.
    link: /applications/prism/initial-conditions
    linkText: Sources
  - icon: 🧱
    title: Boundaries and absorbing layers
    details: PEC, first-order Silver–Müller, periodic, extrapolation, Neumann and zero Dirichlet faces; four PML formulations and a far-wave damping layer.
    link: /applications/prism/pml
    linkText: Absorbing layers
  - icon: ✨
    title: Particle-in-cell
    details: Plasmas or a single test particle; deposition from NGP to sixth order and Gaussian; charge-conserving currents with periodic faces; Runge–Kutta push.
    link: /applications/prism/pic
    linkText: Particle-in-cell
  - icon: 🎯
    title: t0 fields by FLAIL
    details: The electric field of the initial charge and the magnetic field of the initial current, solved with PRISM's own operator on PEC, open and periodic faces.
    link: /applications/prism/numerics#t0-elliptic-solves
    linkText: t0 elliptic solves
  - icon: ⚡
    title: CPU and GPU
    details: The CPU backend (MPI + OpenMP) and the FNL backend (MPI + OpenACC); regression digests on both.
    link: /applications/prism/backends
    linkText: Backends
  - icon: 🛑
    title: Refusals at initialization
    details: One INI file per realm, every value checked; unsupported combinations stop the run before the first step.
    link: /applications/prism/hazards
    linkText: Hazards
---

## Start here

The [overview](./overview) lists what is implemented, with a [quick start](./overview#quick-start), the
[outputs](./overview#outputs) and the [known limitations](./overview#known-limitations). The
[tutorial](/tutorials/prism/pulse) runs a Gaussian pulse in a periodic box, then against reflecting and absorbing faces.

## See it run

<div class="adam-casts adam-casts-wide">
<figure class="adam-cast">
<video class="light-only" autoplay loop muted playsinline poster="/gallery/prism-faces-light.png" aria-label="A pulse leaving through PEC, Silver–Müller and PML faces"><source src="/gallery/prism-faces-light.webm" type="video/webm"><source src="/gallery/prism-faces-light.mp4" type="video/mp4"></video>
<video class="dark-only" autoplay loop muted playsinline poster="/gallery/prism-faces-dark.png" aria-label="A pulse leaving through PEC, Silver–Müller and PML faces"><source src="/gallery/prism-faces-dark.webm" type="video/webm"><source src="/gallery/prism-faces-dark.mp4" type="video/mp4"></video>
<figcaption>A Gaussian pulse in vacuum meets the x faces of the box. Left to right: PEC, first-order Silver–Müller, Silver–Müller with a classic PML (outlined). log₁₀ of √(|E|² + c²|B|²) relative to the initial peak, plane z ≈ 0, to t = 6.67·10⁻¹⁰ s; the three runs share one time step. Inputs of <a href="/adam/tutorials/prism/faces.html">tutorial 2</a> and <a href="/adam/tutorials/prism/pml.html">tutorial 3</a>.</figcaption>
</figure>
</div>

<div class="adam-gallery">

[![PRISM pulse across a 2:1 seam, octree against quadtree](/prism/quadtree-pulse.png)<strong>Quadtree 2:1 seams</strong><span>a pulse across a 2:1 seam, quadtree against octree, per column</span>](./verification/quadtree)

</div>

The animation is made by `docs/scripts/make_gallery_runs.sh` (the inputs of tutorials 2 and 3, run at CFL 0.18 with a
shorter checkpoint period) and `docs/scripts/make_gallery.py`; the figure by `src/tests/prism/regression/make_doc_figures.py`.

## Verification

The [verification timeline](./verification/) lists the oracle checks in the order the capabilities were implemented:
uniform-grid identities, static AMR and the FV reflux, 2:1 seam divergence, quadtrees, inter-realm seams, absorbing
layers on AMR, input refusals, the t0 solves on FLAIL, PIC with periodic faces and the Fresnel slab, slices. Analytic
oracles for the field physics are open ([#106](https://github.com/szaghi/adam/issues/106)).
