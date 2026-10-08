# FLUME

**FLUME** (**F**luid **L**orentz-coupled **U**nsteady **M**agnetohydrodynamic **E**quations) is the compressible-flow
and magnetohydrodynamics application of ADAM. It solves the compressible Euler equations and the ideal MHD equations with
high-order WENO finite differences on block-structured adaptive grids, on CPUs (MPI + OpenMP) and on GPUs (MPI +
OpenACC), and it supersedes the deprecated CHASE.

![Orszag-Tang vortex](/flume/orszag-tang.png)

## What is implemented

| Area | Capabilities | Details |
|---|---|---|
| Physics | Compressible Euler; ideal MHD with mixed GLM divergence cleaning (or without control, for diagnostics); ideal gas | [Physical models](./models) |
| Space | WENO flux splitting (`weno`, orders 3–9, characteristic or conservative); WENO interpolation + Riemann flux + high-order correction (`weno-riemann`, LLF/HLL/HLLC for Euler) | [Numerical methods](./numerics) |
| Time | SSP and low-storage Runge–Kutta, CFL time step (with the GLM cleaning speed) | [Numerical methods](./numerics#time-integration) |
| Grid | Octree/quadtree of Cartesian blocks, initialisation-time AMR (box, gradient and solid markers), 2:1 coarse–fine faces with conservative reflux | [Numerical methods](./numerics#adaptive-mesh-refinement) |
| Geometry | Immersed boundary for static solids (Euler) | [Numerical methods](./numerics#immersed-boundary-euler) |
| Coupling | Multi-realm runs glued by a forest manifest | [Numerical methods](./numerics#multi-realm-runs) |
| Setup | 17 initial conditions (Riemann problems, vortices, Shu–Osher, Orszag–Tang, rotor, field loop, linear waves, exact viscous states, ...), 6 boundary conditions (no-slip walls included) | [Initial conditions](./initial-conditions), [Boundary conditions](./boundary-conditions) |
| Input | One INI file per realm, every value checked | [Input reference](./input) |
| Output | XDMF + HDF5 checkpoints (ParaView), restart, slices, residual, conservation and div(B) histories | [Input reference](./input#io) |
| Verification | Exact solutions, convergence orders, symmetry and conservation oracles on both backends | [Verification gallery](./verification) |

**Status.** Milestones M1 (Euler, [#35](https://github.com/szaghi/adam/issues/35)) and M2 (ideal MHD with GLM,
[#41](https://github.com/szaghi/adam/issues/41)) are complete; M3 ([#47](https://github.com/szaghi/adam/issues/47)) adds
the Riemann-solver scheme (`weno-riemann`: Euler done, MHD in progress), the energy-consistent EGLM cleaning and a
positivity-preserving limiter. M4 ([#65](https://github.com/szaghi/adam/issues/65), in progress) adds the dissipative
terms: compressible Navier–Stokes (viscosity and heat conduction, 4th-order conservative fluxes, no-slip walls) for
every model and Ohmic resistivity for MHD, on both backends, conservative across AMR seams and realms (second order
at a 2:1 seam, as the inviscid fluxes).
Runtime AMR (M5) is planned.

## Quick start

Build (FoBiS, from the repository root):

```bash
fobis build --mode flume-cpu-gnu                         # CPU backend
fobis build --mode flume-fnl-nvf --varset local_nvf      # GPU backend (OpenACC, NVIDIA HPC SDK)
```

Run the Sod shock tube of the verification suite:

```bash
cd src/tests/flume/verification/sod
mpirun -np 2 ../../../../../exe/adam_flume_cpu sod-x.ini
```

The run writes `sod-x-<iteration>.xdmf` (open it with ParaView) and the residual and conservation histories. The
[input reference](./input) explains every key of `sod-x.ini`. To run the same problem with the Riemann-solver scheme,
change its `[numerics]` section to

```ini
[numerics]
scheme_space             = weno-riemann
reconstruction_variables = characteristic
riemann_solver           = hllc
flux_correction          = 6th
flux_correction_sensor   = weno
reflux                   = .true.
```

`./check.sh` in the same directory runs the three directions and compares them with the exact solution.

## Pages

- [Physical models](./models): equations, thermodynamics, eigensystems, GLM cleaning, time-step bound.
- [Numerical methods](./numerics): finite differences, WENO, flux splitting, Riemann solvers and correction, time
  integration, AMR and reflux, immersed boundary, multi-realm.
- [Initial conditions](./initial-conditions) and [boundary conditions](./boundary-conditions): formulas, keys, examples.
- [Input reference](./input): every section and key, with two complete inputs.
- [Verification gallery](./verification): the test cases with figures and measured results.
- Source layout: [common](./common), [CPU backend](./cpu), [FNL backend](./fnl).

## Known limitations

- **Seam div(B).** A cell-centred $\mathbf{B}$ cannot be kept divergence-free across a 2:1 coarse–fine face (Tóth & Roe
  2002): the seams inject a truncation-level $\nabla\cdot\mathbf{B}$, which GLM keeps bounded; `seam_max_divb` in the div(B)
  history monitors it.
- **Positivity at very low plasma β.** The floors are a heuristic under the mixed GLM; the Balsara–Spicer strong blast
  ($\beta = 2.5\cdot10^{-4}$) fails within a few steps. EGLM and a positivity limiter are planned in M3.
- **GLM damping and reflux.** With damping, $\int\psi$ is not conserved across 2:1 faces; the 8 physical integrals are.
- The GPU backend copies the coarse–fine seam faces to the host at every stage, which dominates its run time on AMR
  cases; inter-realm seams must join blocks of the same resolution.
