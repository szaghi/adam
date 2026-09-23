# FLUME

**FLUME** — **F**luid **L**orentz-coupled **U**nsteady **M**agnetohydrodynamic **E**quations — is a compressible magnetohydrodynamics (MHD) application built on the ADAM framework. It solves the equations of an electrically conducting, compressible fluid interacting with its own magnetic field through the Lorentz force.

> **Status: under development.** This document describes the target of the application. The solver is not implemented yet: the directory layout, the physical model and the backends below are the design intent, and every implementation detail is subject to change.

## Physical Models

The first physical model targeted by FLUME is the **inviscid (ideal) compressible MHD** system. Dissipative effects (viscosity, thermal conduction, resistivity) are outside the initial scope; the application name and structure do not preclude adding them later.

| Model | Description |
|-------|-------------|
| Ideal MHD | Inviscid, compressible, perfectly conducting single fluid (ideal Ohm's law $\mathbf{E} + \mathbf{u} \times \mathbf{B} = 0$) |

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

Build modes are not defined yet. They will follow the ADAM naming scheme, e.g.:

```bash
fobis build --mode flume-cpu-gnu                         # CPU backend (GNU compiler)
fobis build --mode flume-fnl-nvf --varset local_nvf      # FNL (OpenACC) backend
```

## License

FLUME is part of the ADAM framework, released under the [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.html) (LGPLv3).

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.
