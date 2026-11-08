---
layout: home

hero:
  name: ADAM
  text: Accelerated fluid Dynamics on Adaptive Mesh refinement grids
  tagline: A physics-agnostic Fortran SDK and its solvers, from a laptop to GPU clusters.
  image:
    light: /logo.svg
    dark: /logo-dark.svg
    alt: ADAM, a brushed bow shock over an adaptive mesh
  actions:
    - theme: brand
      text: Get started
      link: /start/install
    - theme: alt
      text: FLUME
      link: /applications/flume/
    - theme: alt
      text: PRISM
      link: /applications/prism/
    - theme: alt
      text: Developer
      link: /developer/
    - theme: alt
      text: View on GitHub
      link: https://github.com/szaghi/adam

features:
  - icon: 🌲
    title: Adaptive mesh refinement
    details: Morton-ordered octrees and quadtrees of Cartesian blocks, 2:1 balanced, with conservative reflux at coarse–fine faces. FLUME refines at initialization and during the run (box, gradient, Löhner, Jeans-length, solid and user-registered markers); PRISM refines at initialization.
    link: /developer/library#adam-amr-object
    linkText: The AMR objects
  - icon: 🧩
    title: Realms and forests
    details: An application is a realm on the library's objects. A forest couples several realms through 1:1 and 2:1 seams in one time loop, with the seam cadence chosen per seam.
    link: /concepts/forest
    linkText: The forest
  - icon: ⚡
    title: One GPU code base
    details: The FNL backend compiles as OpenACC with nvfortran or as OpenMP offloading with amdflang, with GPU-aware MPI. The CPU backend runs the same inputs with MPI and OpenMP.
    link: /developer/fnl
    linkText: The FNL library
  - icon: 🎯
    title: FLAIL elliptic solvers
    details: One entry point for every elliptic solve. A direct FFT solver on uniform grids and a multigrid-preconditioned conjugate gradient on 2:1 AMR grids, each with a device twin.
    link: /developer/flail
    linkText: FLAIL
  - icon: 🌊
    title: FLUME, Euler and MHD
    details: Compressible Euler, Navier–Stokes and ideal MHD with GLM or EGLM cleaning, WENO flux splitting or WENO with Riemann fluxes, a positivity limiter, ideal, isothermal and barotropic gases.
    link: /applications/flume/models
    linkText: Physical models
  - icon: 🌌
    title: FLUME, gravity and sinks
    details: External gravity with a well-balanced mode, self-gravity (periodic or isolated, uniform or AMR), radiative cooling, dual energy and sink particles that form on the finest level.
    link: /applications/flume/models#external-gravity
    linkText: Gravity and self-gravity
  - icon: 🌀
    title: FLUME, turbulence and walls
    details: Vreman, WALE and Sigma subgrid models, a kinetic-energy-preserving hybrid scheme, equilibrium and inverse-transformation wall models, digital-filter inflow, mass-flow forcing.
    link: /applications/flume/turbulence
    linkText: Turbulence
  - icon: 🔌
    title: PRISM, Maxwell and PIC
    details: Maxwell's equations with centred finite differences or finite volumes, coils, PEC, Silver–Müller and periodic faces, perfectly matched layers, particle-in-cell plasmas and t0 fields from FLAIL.
    link: /applications/prism/
    linkText: PRISM
  - icon: 🔧
    title: Build your own solver
    details: RIPPLE, a scalar advection realm written for the tutorial, shows the contract layer by layer, from one uniform realm to AMR and the GPU.
    link: /tutorials/build-your-own/core
    linkText: The tutorial
---

## In a nutshell

Build FLUME, run the Sod shock tube on two ranks and check it against the exact solution (the
[first run](/start/first-run) explains each step):

```bash
fobis fetch && fobis build --mode flume-cpu-gnu
mkdir -p run/sod && cp src/tests/flume/regression/sod-x/input.ini run/sod/
cd run/sod && mpirun -np 2 ../../exe/adam_flume_cpu input.ini > log.txt
grep '^\[mpi-00000\]' log.txt
```

```text
…
[mpi-00000]t:            1
[mpi-00000]nodes number: 64
[mpi-00000]time step:    0.21128856368212916E-002
[mpi-00000]time:         0.21128856368212916E-002
[mpi-00000]progress:     1%
…
[mpi-00000]t:            174
[mpi-00000]nodes number: 64
[mpi-00000]time step:    0.40494782553279762E-003
[mpi-00000]time:         0.20000000000000001E+000
[mpi-00000]progress:     100%
[mpi-00000]
[mpi-00000]save HDF5 files t: 174, time: 0.20000000000000001E+000
[mpi-00000]ADAM run complete: all realms finalized, entering MPI_Finalize
```

```bash
../../exe/.regression-venv/bin/python ../../src/tests/flume/verification/sod/sod_oracle.py input.ini . --l1-max 3.31e-3
```

```text
/sod-x: axis x, it 174, cells 200, L1(rho) = 3.244140e-03  PASS
```

The run writes one HDF5 file per rank and one XDMF index per save, for ParaView. The oracle compares the density with
the exact Riemann solution. Next: the [tutorials](/tutorials/).

## See it run

<div class="adam-casts">
<figure class="adam-cast">
<video class="light-only" autoplay loop muted playsinline poster="/gallery/orszag-tang-light.png" aria-label="Orszag–Tang vortex, density"><source src="/gallery/orszag-tang-light.webm" type="video/webm"><source src="/gallery/orszag-tang-light.mp4" type="video/mp4"></video>
<video class="dark-only" autoplay loop muted playsinline poster="/gallery/orszag-tang-dark.png" aria-label="Orszag–Tang vortex, density"><source src="/gallery/orszag-tang-dark.webm" type="video/webm"><source src="/gallery/orszag-tang-dark.mp4" type="video/mp4"></video>
<figcaption>Orszag–Tang vortex, ideal MHD with GLM cleaning: density to t = 0.5 on 128², the blocks outlined. FLUME case MV-12.</figcaption>
</figure>
<figure class="adam-cast">
<video class="light-only" autoplay loop muted playsinline poster="/gallery/rayleigh-taylor-light.png" aria-label="Rayleigh–Taylor instability, hydro and MHD, density"><source src="/gallery/rayleigh-taylor-light.webm" type="video/webm"><source src="/gallery/rayleigh-taylor-light.mp4" type="video/mp4"></video>
<video class="dark-only" autoplay loop muted playsinline poster="/gallery/rayleigh-taylor-dark.png" aria-label="Rayleigh–Taylor instability, hydro and MHD, density"><source src="/gallery/rayleigh-taylor-dark.webm" type="video/webm"><source src="/gallery/rayleigh-taylor-dark.mp4" type="video/mp4"></video>
<figcaption>Rayleigh–Taylor single mode under constant gravity, density to t = 5 on 128²: hydrodynamics (left) and MHD with a horizontal field at half the critical value of the mode (right). FLUME case GV-4.</figcaption>
</figure>
<figure class="adam-cast">
<video class="light-only" autoplay loop muted playsinline poster="/gallery/forward-facing-step-light.png" aria-label="Mach 3 forward-facing step on a three-realm forest, density"><source src="/gallery/forward-facing-step-light.webm" type="video/webm"><source src="/gallery/forward-facing-step-light.mp4" type="video/mp4"></video>
<video class="dark-only" autoplay loop muted playsinline poster="/gallery/forward-facing-step-dark.png" aria-label="Mach 3 forward-facing step on a three-realm forest, density"><source src="/gallery/forward-facing-step-dark.webm" type="video/webm"><source src="/gallery/forward-facing-step-dark.mp4" type="video/mp4"></video>
<figcaption>Woodward–Colella Mach 3 step, density to t = 4 on a forest of three quadtree realms, base spacing 1/80 with static refined boxes; the inter-realm seams dashed. FLUME case EV-step.</figcaption>
</figure>
<figure class="adam-cast">
<video class="light-only" autoplay loop muted playsinline poster="/gallery/prism-faces-light.png" aria-label="A pulse leaving through PEC, Silver–Müller and PML faces"><source src="/gallery/prism-faces-light.webm" type="video/webm"><source src="/gallery/prism-faces-light.mp4" type="video/mp4"></video>
<video class="dark-only" autoplay loop muted playsinline poster="/gallery/prism-faces-dark.png" aria-label="A pulse leaving through PEC, Silver–Müller and PML faces"><source src="/gallery/prism-faces-dark.webm" type="video/webm"><source src="/gallery/prism-faces-dark.mp4" type="video/mp4"></video>
<figcaption>A Gaussian pulse in vacuum meets the x faces of the box. Left to right: PEC, first-order Silver–Müller, Silver–Müller with a classic PML (outlined). log₁₀ of √(|E|² + c²|B|²) relative to the initial peak, plane z ≈ 0. PRISM tutorials 2 and 3.</figcaption>
</figure>
</div>

Verification figures, each linked to its case in the [FLUME](/applications/flume/verification/) and
[PRISM](/applications/prism/verification/) timelines:

<div class="adam-gallery">

[![Shock/boundary-layer interaction at Mach 2.28](/flume/sbli-fields.png)<strong>Shock/boundary-layer interaction</strong><span>Mach 2.28, wall-modelled LES: schlieren, pressure and velocity planes</span>](/applications/flume/verification/streams/sbli)

[![Mach 2.5 turbulent boundary layer and the interaction in 3-D](/flume/tbl-3d.png)<strong>Turbulent boundary layer</strong><span>Mach 2.5 from a digital-filter inflow: skin friction 0.2 % below the DNS at δ_i/16</span>](/applications/flume/verification/turbulence/boundary-layer)

[![Vorticity of the Taylor-Green vortex and of isotropic turbulence](/flume/tgv-fields.png)<strong>Taylor–Green and isotropic turbulence</strong><span>vorticity, hybrid-kep-weno at 128³</span>](/applications/flume/verification/turbulence/hybrid-scheme)

[![Streamwise velocity in two wall-modelled channels](/flume/channel-fields.png)<strong>Wall-modelled channel</strong><span>Mach 1.5 and 3: Re_τ 219–222 and 443–445, DNS 220–222 and 451</span>](/applications/flume/verification/turbulence/channel)

[![Woodward-Colella step at t = 4](/flume/step.png)<strong>Forward-facing step</strong><span>three quadtree realms, refined run against the uniform 1/80 control</span>](/applications/flume/verification/euler/forward-facing-step)

[![Rayleigh-Taylor single mode, hydro and MHD](/flume/gravity-rt.png)<strong>Rayleigh–Taylor</strong><span>hydro and MHD, density at t = 2, 3, 4, 5</span>](/applications/flume/verification/gravity/rayleigh-taylor)

[![Balsara-Spicer blast with the positivity limiter](/flume/blast.png)<strong>Balsara–Spicer blast</strong><span>β = 2.5·10⁻⁴ with the positivity limiter, both space schemes</span>](/applications/flume/verification/mhd/positivity-limiter)

[![A self-gravitating star on a static AMR box](/flume/selfgravity-star-amr.png)<strong>Self-gravitating star on AMR</strong><span>the refined box halves the core's spurious swing</span>](/applications/flume/verification/self-gravity/star-amr)

[![Sink particles: binary, accretion rate, collapse on runtime AMR](/flume/sink-particles.png)<strong>Sink particles</strong><span>a binary orbit, the Shu accretion rate, a collapse on runtime AMR</span>](/applications/flume/verification/sinks/sink-particles)

[![PRISM pulse across a 2:1 seam, octree against quadtree](/prism/quadtree-pulse.png)<strong>PRISM, quadtree 2:1 seams</strong><span>a pulse across a 2:1 seam, quadtree against octree</span>](/applications/prism/verification/quadtree)

</div>

The animations are made by `docs/scripts/make_gallery_runs.sh`, which runs committed verification and tutorial inputs
with only the checkpoint period changed, and `docs/scripts/make_gallery.py`, which renders each one on both grounds.
The verification figures are drawn by `make_doc_figures.py` in `src/tests/flume/verification/` and
`src/tests/prism/regression/`. The nutshell output above was captured from `exe/adam_flume_cpu` (gfortran 16.2.0,
Open MPI 5.0.11, 2 ranks), not typed.

## Applications

| Application | Description |
|-------------|-------------|
| [FLUME](/applications/flume/) | Compressible Euler and MHD solver |
| [PRISM](/applications/prism/) | Maxwell equations solver for electromagnetic and plasma (PIC) simulations |
| RIPPLE (`src/app/ripple`) | Scalar advection, the application of the [build-your-own tutorial](/tutorials/build-your-own/core) |

NASTO, CHASE and PATCH are in the source tree but unmaintained, and not documented on this site. ASCOT is the
slice-file converter ([cookbook IO](/cookbook/io)). The [concepts](/concepts/architecture) explain the realm, the
forest and the library objects.

## Authors

Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi ([contacts](/project/citing#authors)).

If you use ADAM in published work, please [cite it](/project/citing). ADAM is dual-licensed under
[MIT or LGPLv3](/project/license).
