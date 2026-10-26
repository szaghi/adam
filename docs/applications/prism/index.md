# PRISM

**PRISM** (**P**lasma **R**esearch us**I**ng **S**imulation **M**ethods) is the electromagnetics and plasma application
of ADAM. It solves Maxwell's equations for the electric displacement $\mathbf D$ and the magnetic induction
$\mathbf B$ with centred finite differences or finite volumes on block-structured adaptive grids, optionally coupled to
charged macro-particles by a particle-in-cell (PIC) method, on CPUs (MPI + OpenMP) and on GPUs (MPI + OpenACC).

## What is implemented

| Area | Capabilities | Details |
|---|---|---|
| Physics | Maxwell's equations in SI or non-dimensional form, with a prescribed current density (coils) or the current and charge of PIC macro-particles; hyperbolic (GLM) divergence cleaning of $\mathbf D$, $\mathbf B$ or both | [Physical models](./models) |
| Space | Centred finite differences or finite volumes of order 2–10 (`[fdv]`); WENO flux splitting (CPU only) | [Numerical methods](./numerics#space-discretization) |
| Time | Runge–Kutta (SSP 11/22/33/54, low-storage 1/2/3, Yoshida on the CPU), CFL time step on the speed of light | [Numerical methods](./numerics#time-integration) |
| Grid | Octree/quadtree of Cartesian blocks, AMR at initialization (box and current-variation markers), 2:1 coarse–fine faces, FV reflux | [Numerical methods](./numerics#adaptive-mesh-refinement) |
| Coupling | Multi-realm runs glued by a forest manifest, α or β seam cadence | [Numerical methods](./numerics#multi-realm-runs) |
| Sources | Rectangular, circular, solenoid and helicon coils, DC or AC, with a smooth ramp; external rotating or uniform fields (acting on the particles) | [Initial conditions and sources](./initial-conditions) |
| Boundaries | Extrapolation, Neumann, Dirichlet (zero), first-order Silver–Müller, PEC, periodic; perfectly matched layers (four formulations) and a far-wave damping layer | [Boundary conditions](./boundary-conditions), [Absorbing layers](./pml) |
| Particles | Plasma (uniform, cylinder, slab, box) or a single test particle; Maxwellian or drift velocities; charge and current deposition from NGP to sixth order and Gaussian; charge-conserving currents with periodic faces; Runge–Kutta push | [Particle-in-cell](./pic) |
| t0 fields | The electric field of the initial charge and the magnetic field of the initial current, solved by FLAIL with PRISM's own operator | [Numerical methods](./numerics#t0-elliptic-solves) |
| Input | One INI file per realm, every value checked, refusals at initialization | [Input reference](./input) |
| Output | XDMF + HDF5 field files (ParaView), restart files, slices, probes, residual, energy, divergence, RMS-field and charge-conservation histories | [Outputs](#outputs) |
| Verification | Regression digests on both backends; oracle checks on AMR seams | [Verification](./verification) |

## Quick start

Build (FoBiS, from the repository root):

```bash
fobis build --mode prism-cpu-gnu                         # CPU backend
fobis build --mode prism-fnl-nvf --varset local_nvf      # GPU backend (OpenACC, NVIDIA HPC SDK)
```

Run the rotating-magnetic-field regression case:

```bash
cd src/tests/prism/regression/rmf
mpirun -np 2 ../../../../../exe/adam_prism_cpu input.ini
```

Without an argument the executable reads `input.ini` in the working directory. A file with a `[forest]` section is a
[forest manifest](/reference/library-input#forest-manifest-multi-realm) and runs several realms.

## Outputs

| File | Content | Enabled by |
|---|---|---|
| `<basename>-<it>-procNNNNNN.h5` + `<basename>-<it>.xdmf` | Field files: the state (`[IO] save_fields = all`) or D and B only (`em`), plus residuals, curls, divergences and the coils' current densities on request | `[IO] it_save` and the `save_*_fields` flags ([input](./input#io)) |
| restart files | State and time, same rank count and grid to restart | `[IO] restart_save` |
| `<basename>-residuals.dat` | L2 norm of every residual | `[IO] residuals_save` |
| `<basename>-energy_history.dat` | the electric and magnetic energies $\tfrac12\sum\lvert\mathbf D\rvert^2/\varepsilon_0\,\Delta V$ and $\tfrac12\sum\lvert\mathbf B\rvert^2/\mu_0\,\Delta V$, the coil power and the Poynting flux | every 10 steps (fixed) |
| `<basename>-divergence_history.dat` | max $\lvert\nabla\cdot\mathbf D\rvert$ ($\lvert\nabla\cdot\mathbf D - \rho\rvert$ with PIC), $\lvert\nabla\cdot\mathbf B\rvert$, $\lvert\nabla\cdot\mathbf J\rvert$ outside the absorbing layers, and the residual of the current solver | every step |
| `<basename>-charge_conservation_history.dat` | PIC: particle and deposited charge, continuity and Gauss residuals | PIC models |
| `<basename>-grms_history.dat`, `<basename>-magnetic_field_at_center_domain.dat` | RMS of the field over a region; $\mathbf B$ at the domain centre | `[grms]`, `[magnetic_field_at_center_domain]` |
| `<basename>-probe_NN.dat` | The state at a fixed point, trilinear | `[probes]` |
| `rectangular_current_coil_N.dat`, `helicon_current_coil_N.dat` | The current of rectangular and helicon coils | every step, when such a coil exists |
| `single_particle_output.dat` | Position, velocity and fields of the test particle | `problem_type = single_particle` |

## Pages

- [Physical models](./models): equations, state vector, constants, non-dimensional form, GLM cleaning.
- [Numerical methods](./numerics): the centred operators, WENO, time integration, AMR and seams, multi-realm runs, the
  t0 elliptic solves.
- [Initial conditions and sources](./initial-conditions): field initial conditions, coils, external fields.
- [Boundary conditions](./boundary-conditions) and [absorbing layers](./pml): the face kinds, the PML, the far-wave layer.
- [Particle-in-cell](./pic): particles, deposition, gather, push, particle boundaries, charge conservation.
- [Input reference](./input): every section and key, with the refused combinations.
- [Backends](./backends): CPU and FNL, what each implements, the step structure, device data.
- [Verification](./verification): the oracle checks.

## Known limitations

- **Seam div(B).** At a 2:1 AMR face the centred operators no longer commute, and the seam injects an $O(h^p)$
  $\nabla\cdot\mathbf B$ source that converges under refinement but grows in time at fixed $h$
  ([#29](https://github.com/szaghi/adam/issues/29), accepted as truncation; an upwinded GLM coupling is
  [#99](https://github.com/szaghi/adam/issues/99)). `[IO] seam_divB_tol` makes it visible
  ([numerics](./numerics#the-2-1-seam-and-div-b)).
- **Refinement at initialization only.** PRISM has no runtime regrid; the t0 elliptic solves refuse a refined grid.
- **Time integration.** Only Runge–Kutta runs: leapfrog, Blanes–Moan and commutator-free Magnus are refused, for the
  fields and for the particles ([#55](https://github.com/szaghi/adam/issues/55)).
- **The GPU backend** has no WENO and no Yoshida; its PML needs the FD operators and an SSP scheme
  ([backends](./backends#what-each-backend-implements)).
- **The CPU charge-conserving deposit and the CPU coherent matched gather** run on one block and one rank only.
- **Coils** (defects, [#108](https://github.com/szaghi/adam/issues/108)): a circular coil injects no current, and
  the coils after the first take the first coil's amplitude, not their own `Amplitude`.
- **The PML and far-wave-layer profiles** are computed per rank, so a layer that spans the blocks of several ranks
  gives results that depend on the decomposition ([#108](https://github.com/szaghi/adam/issues/108)).
- **Verification.** Analytic oracles (plane wave, cavity, PML reflection, Biot–Savart, Larmor orbit, plasma frequency)
  are [#106](https://github.com/szaghi/adam/issues/106).
