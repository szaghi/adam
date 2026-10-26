---
layout: home

hero:
  name: ADAM
  text: Accelerated fluid Dynamics on Adaptive Mesh refinement grids
  tagline: A physics-agnostic Fortran SDK and its solvers, from a laptop to GPU clusters.
  image:
    light: /adam-logo-bga.png
    dark:  /adam-logo-bgai.png
    alt:   ADAM logo
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
    details: Morton-ordered octrees and quadtrees with 2:1 balance and conservative refluxing. Refinement at initialization by geometric, gradient and total-variation markers; FLUME also regrids at runtime.
    link: /concepts/library
    linkText: AMR markers
  - icon: ⚡
    title: One GPU code base
    details: The FNL backend compiles as OpenACC with nvfortran or as OpenMP offloading with amdflang, with GPU-aware MPI; the CPU backend runs the same inputs with MPI and OpenMP.
    link: /concepts/library-fnl
    linkText: The FNL library
  - icon: 🧩
    title: Realms and forests
    details: An application is a realm on the library's objects; a forest couples several realms through mirror and 2:1 seams in one synchronized time loop.
    link: /concepts/forest
    linkText: The forest
  - icon: 🌌
    title: FLUME
    details: Compressible Euler and ideal MHD with dissipative terms, ideal, isothermal and barotropic gases, gravity, self-gravity, cooling and sink particles.
    link: /applications/flume/
    linkText: FLUME
  - icon: 🔌
    title: PRISM
    details: Maxwell equations with particle-in-cell plasmas, coils, PML and perfect-conductor boundaries, and t0 fields from the FLAIL elliptic solvers.
    link: /applications/prism/
    linkText: PRISM
  - icon: 🔧
    title: INI-driven configuration
    details: A run is configured by one INI file, and a missing key stops it at initialization; no recompilation to change physics, numerics or grid.
    link: /applications/flume/input
    linkText: Input reference
---

## What ADAM is

ADAM is a modern Fortran library for block-structured adaptive grids: AMR, immersed boundaries, WENO reconstruction,
Runge-Kutta integration, elliptic solvers (FLAIL) and parallel I/O, on CPUs and GPUs. Solvers are built on it by
supplying only their equations; the two maintained ones are [FLUME](/applications/flume/) and
[PRISM](/applications/prism/). Start with [About](/start/) and [Install](/start/install).

If you use ADAM in published work, please [cite it](/project/citing). ADAM is dual-licensed under
[MIT or LGPLv3](/project/license).
