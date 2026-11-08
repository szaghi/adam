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
      text: Install
      link: /start/install
    - theme: alt
      text: FLUME
      link: /applications/flume/
    - theme: alt
      text: PRISM
      link: /applications/prism/
    - theme: alt
      text: View on GitHub
      link: https://github.com/szaghi/adam

features:
  - icon: 🌲
    title: Adaptive Mesh Refinement
    details: Morton-ordered octrees and quadtrees with 2:1 balance and conservative refluxing. FLUME refines by geometric, gradient, Löhner, Jeans-length and user-registered markers and also regrids at runtime; PRISM refines at initialization by geometric and total-variation markers.
    link: /developer/library#adam-amr-object
    linkText: AMR markers
  - icon: ⚡
    title: One GPU code base
    details: The FNL backend compiles as OpenACC with nvfortran or as OpenMP offloading with amdflang, with GPU-aware MPI; the CPU backend runs the same inputs with MPI and OpenMP.
    link: /developer/fnl
    linkText: The FNL library
  - icon: 🧩
    title: Realms and forests
    details: An application is a realm on the library's objects; a forest couples several realms through mirror and 2:1 seams in one synchronized time loop.
    link: /concepts/forest
    linkText: The forest
  - icon: 🌌
    title: FLUME
    details: Compressible Euler and ideal MHD with dissipative terms, ideal, isothermal and barotropic gases, gravity, self-gravity, cooling, sink particles, and turbulence (SGS models for LES, wall models, digital-filter inflow).
    link: /applications/flume/
    linkText: FLUME
  - icon: 🔌
    title: PRISM
    details: Maxwell equations with particle-in-cell plasmas, coils, PML and perfect-conductor boundaries, and t0 fields from the FLAIL elliptic solvers.
    link: /applications/prism/
    linkText: PRISM
  - icon: 🔧
    title: INI-driven configuration
    details: A run is configured by one INI file, and an unrecognised option value stops it at initialization (the input reference marks the keys left unchecked); no recompilation to change physics, numerics or grid.
    link: /applications/flume/input
    linkText: Input reference
---

## What ADAM is

ADAM (Accelerated fluid Dynamics on Adaptive Mesh refinement grids) is a modern Fortran library for block-structured
adaptive grids: AMR, immersed boundaries, WENO reconstruction, Runge-Kutta integration, elliptic solvers (FLAIL) and
parallel I/O. It runs on CPU clusters with MPI and OpenMP, and on GPU clusters with MPI between ranks and OpenACC
(NVIDIA GPUs) or OpenMP offloading (AMD GPUs) on each node. Solvers are built on it by supplying only their equations.

| Application | Description |
|-------------|-------------|
| [FLUME](/applications/flume/) | Compressible Euler and MHD solver |
| [PRISM](/applications/prism/) | Maxwell equations solver for electromagnetic and plasma (PIC) simulations |
| RIPPLE (`src/app/ripple`) | Scalar advection, the application of the [build-your-own tutorial](/tutorials/build-your-own/core) |

NASTO, CHASE and PATCH are in the source tree but unmaintained, and not documented on this site. ASCOT is the
slice-file converter ([cookbook IO](/cookbook/io)).

Start with [Install](/start/install), then the [first run](/start/first-run) and the [tutorials](/tutorials/); the
[concepts](/concepts/architecture) explain the realm, the forest and the library objects.

If you use ADAM in published work, please [cite it](/project/citing). ADAM is dual-licensed under
[MIT or LGPLv3](/project/license).
