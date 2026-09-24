# FLUME

**FLUME** — **F**luid **L**orentz-coupled **U**nsteady **M**agnetohydrodynamic **E**quations — is a compressible magnetohydrodynamics (MHD) application built on the ADAM framework. It solves the equations of an electrically conducting, compressible fluid interacting with its own magnetic field through the Lorentz force.

> **Status: development.** Milestone M1 of [issue #35](https://github.com/szaghi/adam/issues/35) is delivered: the
> **compressible Euler** equations on both backends, the fluid core on which the MHD model is built. Ideal MHD (the
> magnetic field, its wave families and the divergence control) is the next milestone and is not implemented yet.
> FLUME supersedes CHASE, which is deprecated.

## Physical Models

| Model | `[physics] physical_model` | Status |
|-------|----------------------------|--------|
| Compressible Euler | `euler` | Implemented (M1): inviscid, ideal gas (`cp`, `cv`), $\mathbf{q} = (\rho, \rho u, \rho v, \rho w, E)^\top$ |
| Ideal MHD | — | Target: inviscid, compressible, perfectly conducting single fluid (ideal Ohm's law $\mathbf{E} + \mathbf{u} \times \mathbf{B} = 0$) |

Dissipative effects (viscosity, thermal conduction, resistivity) are outside the initial scope; the application name and
structure do not preclude adding them later.

## Implemented Capabilities (M1)

| Area | What is available | INI |
|------|-------------------|-----|
| Space | WENO flux splitting, per-face Roe eigenvectors, per-wave local Lax-Friedrichs; reconstruction in characteristic or conservative variables; orders `weno-u-3` to `weno-u-9` | `[numerics] scheme_space = weno`, `reconstruction_variables`; `[weno] scheme` |
| Time | Library Runge-Kutta schemes (SSP and low-storage), CFL time step | `[runge_kutta] scheme`; `[time] CFL, it_max, time_max` |
| Boundary conditions | `extrapolation`, `inflow` (primitive state `r, u, v, w, p`), `wall-inviscid`, `periodic` (both faces of an axis or neither) | `[bc_{x,y,z}_{min,max}] type` |
| Initial conditions | `uniform` (optionally with a seeded perturbation), `isentropic-vortex`, `riemann-problem` (piecewise-constant regions) | `[initial_conditions] type` |
| AMR | Init-time refinement (`amr_iterations` passes) by geometric box, variable gradient or immersed-solid surface; 2:1 coarse-fine faces with stage-weighted conservative reflux | `[amr]`, `[initial_conditions] amr_iterations`, `[numerics] reflux` |
| Immersed boundary | Static solids, Euler (inviscid) wall: distance function, eikonal extrapolation into the solid, cut-cell spacing, solid masks in the Runge-Kutta stages | `[solids]` |
| Output | XH5F checkpoints (optionally with the auxiliary fields `u, v, w, p, H, a`), slices, residuals and conservation histories, restart | `[IO]`, `[slices]` |

Every option value is matched against a fixed list; an unknown value stops the run with a message naming the accepted
spellings.

## Governing Equations

FLUME targets the conservative hyperbolic form of the ideal MHD equations:

$$\frac{\partial \mathbf{q}}{\partial t} + \nabla \cdot \mathbf{F}(\mathbf{q}) = \mathbf{S}$$

In units where the magnetic permeability is absorbed into $\mathbf{B}$ ($\mathbf{B} \to \mathbf{B}/\sqrt{\mu_0}$), the system reads:

$$\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{u}) = 0$$

$$\frac{\partial (\rho \mathbf{u})}{\partial t} + \nabla \cdot \left[\rho \mathbf{u} \otimes \mathbf{u} + \left(p + \tfrac{1}{2}|\mathbf{B}|^2\right)\mathbf{I} - \mathbf{B} \otimes \mathbf{B}\right] = 0$$

$$\frac{\partial E}{\partial t} + \nabla \cdot \left[\left(E + p + \tfrac{1}{2}|\mathbf{B}|^2\right)\mathbf{u} - (\mathbf{u} \cdot \mathbf{B})\,\mathbf{B}\right] = 0$$

$$\frac{\partial \mathbf{B}}{\partial t} + \nabla \cdot (\mathbf{u} \otimes \mathbf{B} - \mathbf{B} \otimes \mathbf{u}) = 0$$

with total energy $E = \dfrac{p}{\gamma - 1} + \tfrac{1}{2}\rho|\mathbf{u}|^2 + \tfrac{1}{2}|\mathbf{B}|^2$ (ideal gas closure) and the solenoidal constraint

$$\nabla \cdot \mathbf{B} = 0 .$$

The conservative state vector is $\mathbf{q} = (\rho, \rho u, \rho v, \rho w, E, B_x, B_y, B_z)^\top$. The system is hyperbolic, with seven wave families per direction: two fast magnetosonic, two Alfvén, two slow magnetosonic and one entropy wave.

## Divergence Control

The induction equation preserves $\nabla \cdot \mathbf{B} = 0$ analytically, but not, in general, discretely. Unlike the Maxwell case, a divergence error in MHD feeds back into the momentum and energy equations as a spurious force parallel to $\mathbf{B}$. The divergence control strategy of FLUME is **still to be defined**; the ADAM framework already provides machinery used by PRISM (hyperbolic cleaning, constrained transport) that is a natural starting point.

## Source Layout

The layout mirrors PRISM: backend-independent modules in `common/`, one directory per backend.

```
src/app/flume/
├── common/          # Shared across all backends (physics, numerics, BC, IC, I/O, ...)
├── cpu/             # CPU backend (MPI + OpenMP)
└── fnl/             # OpenACC GPU backend (FNL library)
```

Module and file naming follows the ADAM convention: `adam_flume_<name>_object.F90` for types, `adam_flume_<name>_library.F90` for procedure collections, `adam_flume_cpu.F90` / `adam_flume_fnl.F90` for the program entry points.

## Backends

| Backend | Entry point | Key type | Accelerator |
|---------|-------------|----------|-------------|
| CPU | `adam_flume_cpu.F90` | `flume_cpu_object` | MPI + OpenMP |
| FNL | `adam_flume_fnl.F90` | `flume_fnl_object` | MPI + OpenACC (NVIDIA GPU) |

## Building

```bash
fobis build --mode flume-cpu-gnu                         # CPU backend (GNU compiler)
fobis build --mode flume-cpu-gnu-omp                     # CPU backend with OpenMP threads
fobis build --mode flume-fnl-nvf --varset local_nvf      # FNL (OpenACC) backend
```

A run takes the INI file as its argument: `mpirun -np 2 exe/adam_flume_cpu input.ini`.

## Verification and Regression

`src/tests/flume/verification/` holds one `check.sh` per test, each asserting the physics against an oracle:

| Test | Case | Pass criterion |
|------|------|----------------|
| V0 | `unit/`: pointwise Euler library on random states | eigenvector, Jacobian, round-trip and flux-split identities |
| V1 | `sod/`: Sod along x, y, z; reflecting-wall variant | L1 error against the exact solution; the three directions bitwise identical |
| V2 | `vortex/`: isentropic vortex, 64/128/256 cells | observed order of accuracy |
| V3 | `conservation/`: periodic AMR box with reflux | volume integrals constant to round-off; drift without reflux (negative control) |
| V6 | `shock-cylinder/`: Mach 2 shock over a cylinder, IB + solid AMR | refined surface blocks, mirror symmetry, positivity |
| V7 | `io/`: restart round trip, slices, auxiliary fields | bitwise restart; slice and auxiliary values exact |
| — | `multirealm/`: sod-x split in two realms at the diaphragm, mirror seam, beta cadence; the same with one realm refined (2:1 face crossed by the shock) | union bitwise equal to the single-realm run (issue #37) |

`src/tests/flume/regression/` is the goldened regression suite (a copy of the PRISM harness): `run.sh cpu` runs in CI,
`run-fnl-local.sh` on a GPU workstation. `run-omp-bitwise.sh` (also in CI) requires the OpenMP CPU build to reproduce
the serial one bit for bit on the immersed-boundary case: one unexplained single-ulp divergence is on record (issue #35),
never reproduced since.

Known limitations: the FNL backend copies the coarse-fine seam faces to the host at every stage, which dominates its run
time on AMR cases; inter-realm seams must join blocks of the same resolution (mirror coupling), and the realms are
verified with the same block partition on both sides of the seam (the seam fluxes are matched rank-locally).

## License

FLUME is part of the ADAM framework, released under the [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.html) (LGPLv3).

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.
