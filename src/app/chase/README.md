# CHASE

> CFD-HPC enabled, Adaptive mesh, Simulation code for Euler equations.

CHASE solves the compressible **inviscid** Euler equations on block-structured
adaptive meshes. It is built on the ADAM SDK and shares the same
`common/` + backend pattern as NASTO, but targets the Euler regime (no
viscosity or heat diffusion). Only a CPU backend is provided.

## Governing equations

The compressible Euler equations in conservation form:

$$
\frac{\partial \mathbf{q}}{\partial t}
+ \nabla \cdot \mathbf{F}(\mathbf{q}) = 0
$$

with conservative variables $\mathbf{q} = (\rho,\, \rho u,\, \rho v,\, \rho w,\, \rho E)^\top$
and ideal-gas closure $p = \rho R T$, $E = c_v T + \tfrac{1}{2}(u^2+v^2+w^2)$.

The Euler flux in the $x$-direction is:

$$
\mathbf{F}_x =
\begin{pmatrix}
\rho u \\ \rho u^2 + p \\ \rho v u \\ \rho w u \\ \rho u H
\end{pmatrix},\quad
H = E + \frac{p}{\rho}
$$

with analogous expressions for $\mathbf{F}_y$ and $\mathbf{F}_z$.

---

## Source layout

```
src/app/chase/
├── common/                              # Backend-independent modules
│   ├── adam_chase_common_object.F90     # Base object — ADAM SDK + CHASE handlers
│   ├── adam_chase_physics_object.F90    # Fluid physics — EOS, variable layout, flux functions
│   ├── adam_chase_parameters.F90        # Global constants and eigenvector matrices
│   ├── adam_chase_bc_object.F90         # Boundary conditions handler
│   ├── adam_chase_ic_object.F90         # Initial conditions handler
│   ├── adam_chase_io_object.F90         # IO handler (XH5F output, restart, residuals)
│   ├── adam_chase_time_object.F90       # Time integration parameters and state
│   ├── adam_chase_riemann_library.F90   # Riemann solver library (LLF, HLL)
│   └── adam_chase_common_library.F90    # Barrel re-export of all common modules
└── cpu/                                 # CPU backend
    ├── adam_chase_cpu.F90               # Entry point — parses CLI, calls simulate
    └── adam_chase_cpu_object.F90        # Backend object — working arrays, all methods
```

---

## Variable layout

### Conservative variables — `q(nv=5, ni, nj, nk, nb)`

| Index | Symbol | Description |
|-------|--------|-------------|
| 1 | $\rho$ | Density |
| 2 | $\rho u$ | x-momentum |
| 3 | $\rho v$ | y-momentum |
| 4 | $\rho w$ | z-momentum |
| 5 | $\rho E$ | Total energy |

### Auxiliary variables — `q_aux(nv_aux=11, ni, nj, nk, nb)`

| Index | Symbol | Description |
|-------|--------|-------------|
| 1 | $\rho$ | Density |
| 2–4 | $u, v, w$ | Velocity components |
| 5 | $p$ | Pressure |
| 6 | $T$ | Temperature |
| 7 | $E$ | Total specific energy |
| 8 | $H$ | Total specific enthalpy |
| 9 | $a$ | Speed of sound |
| 10 | $R_f = c_p - c_v$ | Fluid constant |
| 11 | $\gamma = c_p/c_v$ | Specific heats ratio |

---

## `adam_chase_common_object`

Base (non-extending) type that holds all ADAM SDK objects and CHASE handlers.
Both backends compose the CPU object by extending this type.

### Data members

| Member | Type | Description |
|--------|------|-------------|
| `mpih` | `mpih_object` | MPI handler |
| `adam` | `adam_object` | ADAM SDK — grid, tree, field, AMR, IO |
| `field` | `field_object` (pointer) | Field variable storage and ghost-cell maps |
| `grid` | `grid_object` (pointer) | Structured block grid |
| `amr` | `amr_object` | AMR marker handler |
| `ib` | `ib_object` | Immersed Boundary handler |
| `slices` | `slices_object` | Slice sampling for output |
| `rk` | `rk_object` | Runge-Kutta integrator |
| `weno` | `weno_object` | WENO reconstructor |
| `io` | `chase_io_object` | IO handler |
| `physics` | `chase_physics_object` | Physics — EOS, WENO basis |
| `ic` | `chase_ic_object` | Initial conditions |
| `bc` | `chase_bc_object` | Boundary conditions |
| `time` | `chase_time_object` | Time state and parameters |
| `q` | `real(R8P)(nv,…,nb)` | Conservative variables |
| `q_aux` | `real(R8P)(nv_aux,…,nb)` | Auxiliary variables |

### `initialize_common(field, filename, …)`

Parses the INI file; initialises all handlers; builds the ADAM grid, tree, and
field; performs uniform refinement and pruning; binds all pointer members.
Output variable names registered with ADAM IO: `r_a u_a v_a w_a p_a a_a T_a E_a H_a Rfa g_a`.

---

## `adam_chase_physics_object`

Encapsulates the fluid EOS and WENO reconstruction basis selection.

### Key data members

| Member | Description |
|--------|-------------|
| `nv = 5` | Conservative variable count |
| `nv_aux = 11` | Auxiliary variable count |
| `eos` | `eos_ic_object` — ideal-gas state functions (R, cv, g, …) |
| `weno_rec_var` | `'CONSERVATIVE'` or `'CHARACTERISTICS'` |
| `erw`, `elw` | Right/left eigenvector matrices for WENO projection (pointer) |

**WENO reconstruction basis** (set via `[physics]` INI section, key `weno_rec_var`):

| Value | Eigenvectors used | Description |
|-------|------------------|-------------|
| `CONSERVATIVE` | Identity (`IERL`) | Reconstruct conservative variables directly |
| `CHARACTERISTICS` | `ER` / `EL` | Project to characteristics before reconstruction |

### Module-level procedures

| Procedure | Description |
|-----------|-------------|
| `compute_auxiliary(cv, Rf, g, conservative, auxiliary)` | Conservative → auxiliary for a single cell |
| `compute_conservative(auxiliary, conservative)` | Auxiliary → conservative for a single cell |
| `compute_fluxes_conservative(sir, auxiliary, fluxes)` | Euler flux vector in direction `sir` |

---

## `adam_chase_bc_object`

Reads and stores boundary conditions for all 6 domain faces from INI sections
`[bc_x_min]`, `[bc_x_max]`, `[bc_y_min]`, `[bc_y_max]`, `[bc_z_min]`, `[bc_z_max]`.

| Constant | Value | INI keyword | Description |
|----------|-------|-------------|-------------|
| `BC_EXTRAPOLATION` | 1 | `extrapolation` | Zero-order extrapolation into ghost cells |
| `BC_INFLOW` | 2 | `inflow` | Supersonic inflow — prescribed $(\rho, u, v, w, p)$ |
| `BC_WALL_INVISCID` | 3 | `wall-inviscid` | Slip wall (reflective) |

Inflow state is read from INI options `r`, `u`, `v`, `w`, `p` (and `s` for a
perturbation seed) under the relevant face section.

---

## `adam_chase_ic_object`

Sets initial conditions at simulation start and after AMR grid adaptation.

| IC type | INI keyword | Description |
|---------|-------------|-------------|
| Uniform | `uniform` | Constant state with optional random velocity perturbation scaled by `q(6)` |
| Isentropic vortex | `isentropic-vortex` | Analytical isentropic vortex centred at `(emin_x, emin_y)`, radius from `q(6)` |
| Riemann problem | `riemann-problem` | Piecewise constant regions defined by bounding boxes `[emin, emax]` |

INI section: `[initial_conditions]`. Multi-region state per box in
`[initial_conditions_region_N]` (options `r, u, v, w, p, s, emin_x/y/z, emax_x/y/z`).

---

## `adam_chase_time_object`

| INI option (`[time]`) | Default | Description |
|----------------------|---------|-------------|
| `it_max` | −1 | Maximum time steps (−1 = time-based termination) |
| `time_max` | 1.0 | Maximum simulation time |
| `CFL` | 0.3 | CFL stability coefficient |

---

## `adam_chase_riemann_library`

Riemann solver library providing LLF and HLL solvers and Maxwell flux functions.

| Procedure | Description |
|-----------|-------------|
| `compute_riemann_maxwell_llf(sir, nv, q1, q4, f, lmax)` | LLF (Rusanov) solver; wave speed $c_0 = 1/\sqrt{\mu_0 \varepsilon_0}$ |
| `compute_riemann_maxwell_hll(sir, nv, q1, q4, f, lmax, lmin)` | HLL solver |
| `compute_convective_fluxes_maxwell(sir, q, f)` | Maxwell flux vector |
| `compute_convective_fluxes_maxwell_div_d(sir, q, f)` | Maxwell flux with $\nabla\cdot\mathbf{D}$ correction |
| `compute_convective_fluxes_maxwell_div_d_b(sir, q, f)` | Maxwell flux with $\nabla\cdot\mathbf{D}$ and $\nabla\cdot\mathbf{B}$ corrections |

---

## CPU backend — `chase_cpu_object`

Extends `chase_common_object` with working arrays and the full suite of compute
methods.

### Inheritance

```
chase_common_object          (common/)
    └── chase_cpu_object     (cpu/)   ← backend object
```

### Additional data members

| Member | Shape | Description |
|--------|-------|-------------|
| `dq` | `(nv, 1-ngc:ni+ngc, …, nb)` | Residual (right-hand side) |
| `flx` | `(nv, 1-ngc:ni+ngc, …, nb)` | Convective face fluxes in x |
| `fly` | `(nv, 1-ngc:ni+ngc, …, nb)` | Convective face fluxes in y |
| `flz` | `(nv, 1-ngc:ni+ngc, …, nb)` | Convective face fluxes in z |

### Methods

**Lifecycle**

| Procedure | Purpose |
|-----------|---------|
| `initialize(filename)` | Calls `initialize_common`; then `allocate_cpu` for working arrays |
| `allocate_cpu()` | Allocates `dq`, `flx`, `fly`, `flz`; zero-initialises |
| `simulate(filename)` | Full simulation driver: init → IC/restart → time loop → finalise |

**AMR**

| Procedure | Purpose |
|-----------|---------|
| `amr_update()` | Iterates AMR markers; calls `adam%amr_update`; recomputes IB `phi` after each pass |
| `mark_by_geo(delta_fine, delta_coarse)` | Marks blocks by proximity to IB surface (`phi` sign change) |
| `mark_by_grad_var(grad_tol, delta_type, …)` | Marks blocks by max gradient of a field variable |
| `refine_uniform(levels)` | Uniformly refines all blocks |
| `compute_phi()` | Recomputes signed-distance IB field after grid changes |
| `move_phi(velocity, s)` | Translates IB solid `s` by `velocity` and updates `phi` |

**Immersed boundary**

| Procedure | Purpose |
|-----------|---------|
| `integrate_eikonal(q)` | Runs `n_eikonal` Godunov sweeps; applies `invert_eikonal` |

**I/O**

| Procedure | Purpose |
|-----------|---------|
| `save_xh5f([output_basename])` | Writes `q` (`r,ru,rv,rw,rE`) and `q_aux` to XH5F with Morton cell ordering |
| `save_restart_files()` | Saves ADAM binary restart + an XH5F snapshot |
| `load_restart_files(t, time)` | Reads restart binary; rebuilds MPI communication maps |
| `save_residuals()` | Computes per-variable L2 norm of `dq` via `MPI_ALLREDUCE(SUM)` |
| `save_simulation_data()` | Schedules output: XH5F snapshots, restart checkpoints, slice sampling |

Slice data is written in `.mat` format with variable names `rho rhu rhv rhw rhe`.

**Initial and boundary conditions**

| Procedure | Purpose |
|-----------|---------|
| `set_initial_conditions()` | Delegates to `ic%set_initial_conditions` |
| `set_boundary_conditions(q)` | Iterates `local_map_bc_crown`; applies extrapolation or inflow BCs |
| `update_ghost(q [,step])` | Local inter-block copy → MPI exchange → BC application |

**Numerical methods**

| Procedure | Purpose |
|-----------|---------|
| `compute_dt()` | CFL stability criterion (see below) |
| `compute_q_auxiliary(q, q_aux)` | Full-field conservative → auxiliary conversion |
| `compute_residuals(q, q_aux, dq)` | Ghost sync → eikonal → aux vars → fluxes → divergence |
| `integrate()` | Runge-Kutta time advance; same scheme dispatch as NASTO CPU |
| `simulate(filename)` | Full time loop driver |

### Time-step control

`compute_dt` uses a convective-only CFL bound (no viscous Fourier term):

$$
\Delta t = \text{CFL} \cdot \left[
  \max_{b,i,j,k} \left(
    \frac{|u|+a}{\Delta x/2} + \frac{|v|+a}{\Delta y/2} + \frac{|w|+a}{\Delta z/2}
  \right)
\right]^{-1}
$$

A global minimum is then taken via `MPI_ALLREDUCE(MPI_MIN)`.

### Residual computation pipeline

```
1. update_ghost(q)                   ! sync ghost cells + apply BC
2. integrate_eikonal(q)              ! update IB distance field
3. compute_q_auxiliary(q, q_aux)     ! q → primitive/auxiliary vars
4. compute_fluxes_convective(x,y,z)  ! WENO + LLF Euler fluxes per direction
5. compute_fluxes_difference(…)      ! flux divergence → dq; IB masking
```

---

## `adam_chase_common_library`

Barrel re-export module. A single `use adam_chase_common_library` brings all
common modules into scope:

```fortran
use adam_chase_bc_object
use adam_chase_common_object
use adam_chase_ic_object
use adam_chase_io_object
use adam_chase_parameters
use adam_chase_physics_object
use adam_chase_time_object
```

---

## Build

```bash
# GNU compiler (production)
FoBiS.py build -mode chase-gnu

# GNU compiler (debug)
FoBiS.py build -mode chase-gnu-debug
```

The executable `adam_chase_cpu` is placed in `exe/`.

## Running

```bash
# Single rank, default INI filename (input.ini)
mpirun -np 1 exe/adam_chase_cpu

# Multiple ranks, custom INI
mpirun -np 8 exe/adam_chase_cpu my_case.ini
```

An example INI file is provided at `src/app/chase/cpu/adam_nasto.ini`.

---

## Copyrights

CHASE is part of the ADAM framework, released under the
[GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.html)
(LGPLv3).

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.
