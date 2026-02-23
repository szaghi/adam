# NASTO — CPU Backend

> Pure CPU backend: MPI distributed memory + OpenMP shared memory.

The `cpu/` directory provides the reference CPU implementation of the NASTO
solver. It extends the backend-independent `nasto_common_object` with all
compute kernels running exclusively on the host, using MPI for distributed
memory and OpenMP `collapse` loops for shared-memory parallelism within each
rank.

## Source layout

| Source file | Module / Program | Role |
|-------------|-----------------|------|
| `adam_nasto_cpu.F90` | program `adam_nasto_cpu` | Entry point — parses CLI, calls `simulate` |
| `adam_nasto_cpu_object.F90` | module `adam_nasto_cpu_object` | Backend object — AMR, IB, IO, BC, numerics |
| `adam_nasto_cpu_cns.F90` | module `adam_nasto_cpu_cns` | Pure CNS kernel library — fluxes, Riemann, aux vars |

---

## `adam_nasto_cpu` — main program

The entry point is minimal by design:

```fortran
type(nasto_cpu_object) :: nasto
! Read optional CLI argument for INI filename; default "input.ini"
call nasto%simulate(filename=...)
```

The INI filename is taken from the first command-line argument if present,
otherwise it defaults to `input.ini`.

---

## `nasto_cpu_object`

Extends `nasto_common_object` (from `common/`) with CPU-specific working arrays
and the full suite of compute methods.

### Inheritance

```
nasto_common_object          (common/)
    └── nasto_cpu_object     (cpu/)   ← backend object
```

### Additional data members

| Member | Shape | Description |
|--------|-------|-------------|
| `dq` | `(nv, 1-ngc:ni+ngc, …, nb)` | Residual (right-hand side) of the CNS equations |
| `flx` | `(nv, 1-ngc:ni+ngc, …, nb)` | Convective + diffusive face fluxes in x |
| `fly` | `(nv, 1-ngc:ni+ngc, …, nb)` | Convective + diffusive face fluxes in y |
| `flz` | `(nv, 1-ngc:ni+ngc, …, nb)` | Convective + diffusive face fluxes in z |

All arrays follow the standard 5D layout `(nv, ni, nj, nk, nb)` with
ghost-cell extensions.

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
| `amr_update()` | Iterates AMR markers; calls `adam%amr_update`; recomputes IB phi after each pass; exits early when grid is stable |
| `mark_by_geo(delta_fine, delta_coarse)` | Marks blocks by proximity to IB surface (`phi` sign change) |
| `mark_by_grad_var(grad_tol, delta_type, …)` | Marks blocks by the max gradient of a selected field variable |
| `refine_uniform(levels)` | Uniformly refines all blocks by the requested number of levels |
| `compute_phi()` | Recomputes the signed-distance IB field after grid changes |
| `move_phi(velocity, s)` | Translates IB solid `s` by `velocity` and updates `phi` |

**Immersed boundary**

| Procedure | Purpose |
|-----------|---------|
| `integrate_eikonal(q)` | Runs `n_eikonal` Godunov sweeps of the eikonal equation; applies ghost-cell sync and `invert_eikonal` to establish the signed-distance field inside solids |

**I/O**

| Procedure | Purpose |
|-----------|---------|
| `save_hdf5([output_basename])` | Writes `q` (conservative) and `q_aux` (primitive + derived) to HDF5; optionally appends `phi` for IB runs |
| `save_restart_files()` | Writes ADAM binary restart + an HDF5 snapshot with the restart basename |
| `load_restart_files(t, time)` | Reads restart binary; rebuilds MPI communication maps |
| `save_residuals()` | Computes global L2 norm of `dq` via `MPI_ALLREDUCE`; appends to residuals file |
| `save_simulation_data()` | Orchestrates scheduled output: HDF5 snapshots, restart checkpoints, and slice sampling |

HDF5 variable names written to file:

| Array | Variable names |
|-------|----------------|
| `q` | `rho rhu rhv rhw rhe` |
| `q_aux` | `rhob u v w ya tem pres ental csp` |

**Initial and boundary conditions**

| Procedure | Purpose |
|-----------|---------|
| `set_initial_conditions()` | Delegates to `ic%set_initial_conditions` |
| `set_boundary_conditions(q)` | Iterates `local_map_bc_crown` (precomputed ghost-cell index map); applies extrapolation or inflow from `bc%q`; wall BC applied via MPI-exchanged interior values |
| `update_ghost(q [,step])` | Local inter-block copy → MPI exchange → BC application; optional `step` parameter enables async overlap |

**Numerical methods**

| Procedure | Purpose |
|-----------|---------|
| `compute_dt()` | CFL + viscous Fourier stability criterion (see below) |
| `compute_q_auxiliary(q, q_aux)` | Converts conservative `q` → primitive `q_aux` via `compute_q_aux` from CNS library |
| `compute_residuals(q, dq)` | Ghost sync → eikonal → aux vars → convective fluxes → diffusive correction → flux divergence into `dq` |
| `integrate()` | Runge-Kutta time advance; dispatches between low-storage (RK1/2/3) and SSP (SSP-RK2/3, SSP-RK54) schemes |

### Time-step control

`compute_dt` uses a combined explicit stability bound across all active blocks:

$$
\Delta t = \text{CFL} \cdot \left[
  \max_{b,i,j,k} \left(
    \frac{|u|+c}{\Delta x/2} + \frac{|v|+c}{\Delta y/2} + \frac{|w|+c}{\Delta z/2}
    + \frac{2\mu}{\rho(\Delta x/2)^2}
    + \frac{2\mu}{\rho(\Delta y/2)^2}
    + \frac{2\mu}{\rho(\Delta z/2)^2}
  \right)
\right]^{-1}
$$

The convective (CFL) and viscous (Fourier) contributions are summed before
taking the global minimum. An `MPI_ALLREDUCE` with `MPI_MIN` reduces across
all ranks.

### Residual computation pipeline

`compute_residuals` executes the following sequence for each time step:

```
1. update_ghost(q)                   ! sync ghost cells + apply BC
2. integrate_eikonal(q)              ! update IB distance field
3. compute_q_auxiliary(q, q_aux)     ! q → primitive vars
4. compute_fluxes_convective(x,y,z)  ! LLF flux splitting per direction
5. compute_fluxes_diffusive(…)       ! central 2nd-order viscous correction (if μ > 0)
6. compute_fluxes_difference(…)      ! flux divergence → dq; IB masking if solids present
```

### Integration loop

`integrate` supports two Runge-Kutta families:

| Scheme constant | Stages | Type | Description |
|----------------|--------|------|-------------|
| `RK_1` | 1 | Low-storage | Forward Euler |
| `RK_2` | 2 | Low-storage | 2-stage low-storage RK |
| `RK_3` | 3 | Low-storage | 3-stage low-storage RK |
| `RK_SSP_22` | 2 | SSP | Shu-Osher SSP-RK(2,2) |
| `RK_SSP_33` | 3 | SSP | Shu-Osher SSP-RK(3,3) |
| `RK_SSP_54` | 4 | SSP | Gottlieb SSP-RK(5,4) |

Low-storage schemes update `q` in-place at each stage.
SSP schemes store all intermediate stages in `rk%q_rk(:,:,:,:,:,s)` then
combine via the Butcher tableau coefficients.

### Simulation driver

`simulate` orchestrates the full time loop:

```
initialize(filename)
  ├── initialize_common()      ! parse INI, build grid, AMR, field, WENO, IB
  └── allocate_cpu()           ! allocate dq, flx, fly, flz

if restart: load_restart_files()
else:
  loop ic%amr_iterations:
    set_initial_conditions() → compute_phi() → amr_update()
  set_initial_conditions()     ! final IC on converged grid

save_simulation_data()         ! t=0 snapshot

time loop:
  if (mod(it, amr%frequency)==0): amr_update()
  compute_dt()
  integrate()                  ! one RK step
  time += dt; it += 1
  save_simulation_data()       ! HDF5 / restart / slices (scheduled)
  if termination criterion: exit

save_simulation_data()         ! final snapshot
close residuals file
```

---

## `adam_nasto_cpu_cns` — CNS kernel library

Pure, side-effect-free procedures implementing the compressible Navier-Stokes
building blocks. Used by `nasto_cpu_object` but designed as a standalone library
(no type dependencies).

### Procedures

**Conservative / auxiliary variable conversions**

| Procedure | Description |
|-----------|-------------|
| `compute_conservatives(b,i,j,k,ngc,q_aux,q)` | `q_aux → q` for a single cell: `ρ = ρ`, `ρu = ρ·u`, `ρE = ρ·H − p` |
| `compute_conservatives_scalar(…)` | Scalar variant |
| `compute_q_aux(ni,nj,nk,ngc,nb,R,cv,g,q,q_aux)` | Full-field conservative → auxiliary bulk conversion |

**Convective fluxes**

| Procedure | Description |
|-----------|-------------|
| `compute_conservative_fluxes(sir,b,i,j,k,ngc,q_aux,f)` | Euler flux vector in direction `sir` for a single cell |
| `compute_conservative_fluxes_scalar(…)` | Scalar variant |
| `compute_max_eigenvalues(si,sir,weno_s,…,evmax)` | Local Lax-Friedrichs wave-speed estimate over the WENO stencil |
| `compute_eigenvectors(si,sir,…,el,er)` | Left/right eigenvector matrices of the flux Jacobian |
| `decompose_fluxes_convective_llf(sir,nv,q_aux,evmax,fmp)` | LLF flux splitting: `f±` decomposition |

**Riemann solvers** (all elemental / pure)

| Procedure | Description |
|-----------|-------------|
| `compute_riemann_exact` | Exact iterative Riemann solver |
| `compute_riemann_exact_2`, `compute_riemann_exact_3` | Exact solver variants |
| `compute_riemann_hllc` | HLLC approximate Riemann solver |
| `compute_riemann_llf` | Local Lax-Friedrichs (Rusanov) solver |
| `compute_riemann_ts` | Two-shock approximate Riemann solver |

---

## Build

```bash
# GNU compiler (production)
FoBiS.py build -mode nasto-cpu-gnu

# GNU compiler (debug — bounds checking + traceback)
FoBiS.py build -mode nasto-cpu-gnu-debug
```

The executable `adam_nasto_cpu` is placed in `exe/`.

## Running

```bash
# Single rank, default INI filename (input.ini)
mpirun -np 1 exe/adam_nasto_cpu

# Multiple ranks, custom INI
mpirun -np 8 exe/adam_nasto_cpu my_case.ini
```

An example INI file is provided at `src/app/nasto/cpu/adam_nasto.ini`.

---

## Copyrights

NASTO is part of the ADAM framework, released under the
[GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.html)
(LGPLv3).

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.
