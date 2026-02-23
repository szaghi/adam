# PRISM — CPU Backend

> PRISM is part of the ADAM framework, released under the
> [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.en.html) (LGPLv3).
> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.

The CPU backend (`prism-cpu`) provides a pure CPU/MPI implementation of the PRISM Maxwell equations
solver. It inherits the full PRISM common infrastructure (grid, field, physics, numerics, I/O, PIC,
coils, boundary conditions) and adds flux storage arrays, spatial-operator procedure pointers, and a
comprehensive set of time integration schemes.

## Source files

| File | Module |
|------|--------|
| `adam_prism_cpu.F90` | entry point — parses `--filename` CLI argument, delegates to `prism%simulate` |
| `adam_prism_cpu_object.F90` | `adam_prism_cpu_object` — type declaration, all procedures |

## Build mode

```bash
FoBiS.py build -mode prism-gnu        # CPU build with GNU
FoBiS.py build -mode prism-gnu-debug  # debug build
```

## Type: `prism_cpu_object`

`prism_cpu_object` extends `prism_common_object` (from `adam_prism_common_object`) and adds:

### Additional flux storage

| Member | Shape | Description |
|--------|-------|-------------|
| `flxyz_c` | `(nv, 2, 3, ni, nj, nk, nb)` | Cell-center fluxes with ±-decomposition for all 3 directions |
| `flx_f` | `(nv, 0:ni, nj, nk, nb)` | X-face fluxes (index 0 = left boundary) |
| `fly_f` | `(nv, ni, 0:nj, nk, nb)` | Y-face fluxes |
| `flz_f` | `(nv, ni, nj, 0:nk, nb)` | Z-face fluxes |

### Procedure pointers

Nine type-bound procedure pointers select the active spatial/temporal operators at run time:

| Pointer | Signature |
|---------|-----------|
| `compute_curl` | `(ivar, q, curl)` |
| `compute_derivative1` | `(dir, ivar, q, dq_ds)` |
| `compute_derivative2` | `(dir, ivar, q, d2q_ds2)` |
| `compute_derivative4` | `(dir, ivar, q, d4q_ds4)` |
| `compute_divergence` | `(ivar, q, divergence)` |
| `compute_gradient` | `(ivar, q, gradient)` |
| `compute_laplacian` | `(ivar, q, laplacian)` |
| `compute_residuals` | `(q, dq [, s])` |
| `integrate` | `()` |

A module-level procedure pointer `compute_fluxes_Maxwell` selects the Maxwell flux function based
on the divergence-cleaning configuration.

## Dispatch table (`initialize`)

After calling `initialize_common`, `allocate_cpu` allocates the flux arrays, then all pointers are
set based on INI parameters:

### Spatial operator pointers

| `scheme_space` | `compute_curl/divergence/gradient/laplacian/derivative1/2` | `compute_residuals` |
|----------------|-------------------------------------------------------------|---------------------|
| `WENO` | `*_fv` variants | `compute_residuals_weno` |
| `FD_CENTERED` | `*_fd` variants | `compute_residuals_fd_centered` |
| `FV_CENTERED` | `*_fv` variants | `compute_residuals_fv_centered` |

### `compute_derivative4` pointer

Always bound to `compute_derivative4_fd` (Kreiss-Oliger dissipation, FD only).

### Time integration pointer

| Physical model | `scheme_time` | RK variant | `integrate` |
|----------------|--------------|------------|-------------|
| Maxwell (non-PIC) | `LEAPFROG` | — | `integrate_leapfrog` |
| Maxwell (non-PIC) | `RK` | `LS` | `integrate_rk_ls` |
| Maxwell (non-PIC) | `RK` | `SSP` | `integrate_rk_ssp` |
| Maxwell (non-PIC) | `RK` | `YOSHIDA` | `integrate_rk_yoshida` |
| Maxwell (non-PIC) | `BLANES_MOAN` | — | `integrate_blanesmoan` |
| Maxwell (non-PIC) | `CFM` | — | `integrate_cfm` |
| PIC | `LEAPFROG` | — | `integrate_leapfrog_pic` |
| PIC | `RK` | `SSP` | `integrate_rk_ssp_pic` |

### `compute_fluxes_Maxwell` pointer

| `div_corr_var` | `constrained_transport_B` | `constrained_transport_D` | Pointer |
|----------------|--------------------------|--------------------------|---------|
| `NONE` | false | false | `compute_fluxes_Maxwell_6v` |
| `HYPER` | false | false | `compute_fluxes_Maxwell_hyper_D` |
| `HYPER` | true | false | `compute_fluxes_Maxwell_hyper_B` |
| `HYPER` | true | true | `compute_fluxes_Maxwell_hyper_DB` |

(All variants from `adam_prism_common_library`.)

## Methods

### Lifecycle and I/O

#### `allocate_cpu`
Allocates the four flux arrays. Face flux arrays use one-sided extended bounds to store both
left- and right-boundary face values:

```fortran
allocate(flxyz_c(nv, 2, 3, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, nb))
allocate(flx_f  (nv, 0:ni, nj, nk, nb))
allocate(fly_f  (nv, ni, 0:nj, nk, nb))
allocate(flz_f  (nv, ni, nj, 0:nk, nb))
```

#### `save_residuals`
Computes per-block L2 residual norms, reduces via `MPI_ALLREDUCE(SUM)`, normalizes by
`sqrt(ni × nj × nk)`, and writes from rank 0 via `io%write_residuals`.

#### `save_simulation_data`
Triggered by `it_save`, `restart_save`, or active slice monitors. Calls `update_ghost` and
`compute_auxiliary_fields`, then saves XH5F field files, optional restart files, and optional
slices. In `single_particle` mode also writes `single_particle_output.dat`.

### Simulation loop (`simulate`)

```
initialize
if restart: load_restart_files
else:        (AMR iterations) → set_initial_conditions
             set_initial_conditions (final)

compute initial divergence of D, B, J
save_simulation_data
compute_energy
open residuals file

if LEAPFROG: bootstrap with explicit Euler step
if PIC+LEAPFROG: bootstrap PIC with explicit Euler step

do
  compute_dt (CFL, MPI_ALLREDUCE MIN)
  clip dt to time_max if needed
  integrate          ← dispatched pointer
  time += dt
  save_simulation_data
  compute_energy
  if converged: exit
enddo

compute_energy_error
save final data
print initial/final energies and RMS errors
MPI finalize
```

Exit conditions: `time >= time_max` (for `it_max <= 0`) or `it >= it_max` (for `it_max > 0`).

### Ghost cell update (`update_ghost`)

Two-step asynchronous or synchronous ghost update:

| Step | Action |
|------|--------|
| 1 | `field%update_ghost_local` — block-local halo fill |
| MPI | `field%update_ghost_mpi` — MPI nearest-neighbor exchange |
| 3 | `apply_fWL_correction`, `set_boundary_conditions` |

When `step` is not present all three are performed synchronously. When `step=1` only the local
update runs; when `step=3` only the BC part runs (allowing overlap with MPI communication).

### Boundary conditions (`set_boundary_conditions`)

Iterates `local_map_bc_crown` over all `ngc` ghost crowns. For each ghost cell the `bc_type`
selects:

| BC type | Ghost cell value |
|---------|-----------------|
| `EXTRAPOLATION` | Copy from interior: `q(ghost) = q(interior)` |
| `NEUMANN` | Mirror: `q(ghost) = q(mirror)` |
| `Silver_Müller` | 1st-order absorbing (Barbas notation): face-direction-dependent `α,β,γ` indices; `q_D_α = s₁ C₀ q_B_β ε₀`, `q_D_β = −s₁ C₀ q_B_α ε₀`, `q_D_γ = q_D_γ`, `q_B_α = −s₁/C₀ q_D_β/ε₀`, `q_B_β = s₁/C₀ q_D_α/ε₀`, `q_B_γ = q_B_γ`; particle variables extrapolated |
| `DIRICHLET` | Set all variables to zero |
| `PERIODIC` | Wrap-around: `q(ghost) = q(interior + N)` per face direction |
| `BC_radiative` | RK-BC absorbing: `q(ghost) = 2 q_bc_rk(stage) − q(interior)` |

#### Radiative BC residuals (`compute_residuals_BC`)
First-order upwind along the outward normal:

$$\frac{\partial q}{\partial t} = -C_0 \frac{q_\text{bc} - q_\text{interior}}{\Delta s}$$

with an additional damping term `−χ C₀ / Δs` applied to the divergence-cleaning variables.

### Spatial discretization

All FDV operators are dispatched through procedure pointers. Concrete implementations:

#### `compute_residuals_fd_centered`
Centered finite-difference curl Maxwell RHS:

$$\frac{\partial \mathbf{D}}{\partial t} = \nabla \times \frac{\mathbf{B}}{\mu_0} - \mathbf{J}, \qquad
\frac{\partial \mathbf{B}}{\partial t} = -\nabla \times \frac{\mathbf{D}}{\varepsilon_0}$$

Calls `compute_curl_fd_centered` at each cell using stencil half-width `s1 = fdv_half_stencils(1)`.

#### `compute_residuals_fv_centered`
Finite-volume flux-based RHS:
1. Evaluate Maxwell fluxes at all cell centers (including ghosts) via `compute_fluxes_Maxwell`.
2. Reconstruct face values using `compute_reconstruction_r_fv_centered`.
3. Compute flux divergence:

$$\frac{\partial q_v}{\partial t} = -\frac{F^x_{i+1/2} - F^x_{i-1/2}}{\Delta x}
                                   -\frac{F^y_{j+1/2} - F^y_{j-1/2}}{\Delta y}
                                   -\frac{F^z_{k+1/2} - F^z_{k-1/2}}{\Delta z}$$

4. Subtract J source on Dx, Dy, Dz components.

#### `compute_residuals_weno`
High-order WENO upwind reconstruction:
1. `compute_fluxes_convective_weno` in each direction: characteristic decomposition
   (`decompose_fluxes_convective`) → WENO reconstruction (`weno_reconstruct_upwind`) → back
   projection into conservative variables (`erw` eigenvectors).
2. `compute_fluxes_difference`: flux divergence + J source subtraction.

Characteristic decomposition uses global Lax-Friedrichs splitting:

$$f^{\pm}_v = \frac{1}{2}(g_v \pm \lambda_{\max} w_v)$$

with $\lambda_{\max} = \texttt{evmax}$ (speed of light), $w = \mathbf{L} q$,
$g = \mathbf{L} f$, $\mathbf{L}$ the left eigenvector matrix.

#### FDV operator library

| Operator | FD | FV | Stencil half-width |
|----------|----|----|-------------------|
| `compute_curl` | `compute_curl_fd_centered` | `compute_curl_fv_centered` | `fdv_half_stencils(1)` |
| `compute_derivative1` | `compute_derivative1_fd_centered` | `compute_derivative1_fv_centered` | `fdv_half_stencils(1)` |
| `compute_derivative2` | `compute_derivative2_fd_centered` | *(stub)* | `fdv_half_stencils(2)` |
| `compute_derivative4` | `compute_derivative4_fd_centered` | — | `fdv_half_stencil` |
| `compute_divergence` | `compute_divergence_fd_centered` | `compute_divergence_fv_centered` | `fdv_half_stencils(1)` |
| `compute_gradient` | `compute_gradient_fd_centered` | `compute_gradient_fv_centered` | `fdv_half_stencils(1)` |
| `compute_laplacian` | `compute_laplacian_fd_centered` | *(stub)* | `fdv_half_stencils(2)` |

All operators iterate `(b, k, j, i)` loops and call the corresponding scalar/vector kernel from
`adam_fdv_operators_library`.

### Time integration schemes

#### Leapfrog (`integrate_leapfrog`)
```
compute_coils_current
compute_residuals → dq
save_residuals
leapfrog%integrate(dt, q, dq)
impose_div_free
```

#### Leapfrog-PIC (`integrate_leapfrog_pic`)
```
pic%particle_cartesian_grid_index
pic%current_weighting              ← J from particles
compute_coils_current
compute_residuals → dq
save_residuals
pic%field_weighting                ← fields at particle positions
leapfrog%integrate(dt, q, dq)
leapfrog_pic%integrate(dt, q_pic, pic_fields)
impose_div_free
```

#### Low-storage RK (`integrate_rk_ls`)
Classic low-storage RK (Williamson form). Stages accumulate in place on `q` using `rk%compute_stage_ls`.

#### SSP RK (`integrate_rk_ssp`)
Strong Stability Preserving RK. Each stage `s`:
1. `rk%compute_stage` — advance `q_rk(:,:,:,:,:,s)` from `q_rk(:,:,:,:,:,s-1)`.
2. `compute_coils_current(q=q_rk(:,:,:,:,:,s), gamma=rk%gamm(s))` — interpolated coil current.
3. `compute_residuals` — Maxwell RHS.
4. `rk%assign_stage` — store stage residuals.

After all stages: `rk%update_q`, `update_q_BC` (radiative BC), `compute_coils_current`,
`impose_div_free`. Supports external fields via `sub_external_fields` / `add_external_fields`.

#### SSP RK-PIC (`integrate_rk_ssp_pic`)
Same as SSP RK but co-integrates field (`rk`) and particle (`rk_pic`) stages:
- `pic%current_weighting` at each stage for J source terms.
- `pic%field_weighting` at each stage for Lorentz forces on particles.
- Final `rk_pic%update_q_pic` advances particle positions/velocities.

#### Yoshida symplectic (`integrate_rk_yoshida`)
4th-order symplectic using Yoshida coefficients `ssa`/`ssb`. Alternates B-half-kick and
D-full-kick in `nrk−1` sub-steps, then a final B-half-kick:

```
for s=1..nrk-1:
  compute_residuals → dq
  B += ssa(s)*dt*dq_B
  compute_residuals → dq
  D += ssb(s)*dt*dq_D
compute_residuals → dq
B += ssa(nrk)*dt*dq_B
impose_div_free
```

#### Blanes-Moan symplectic (`integrate_blanesmoan`)
6th-order symplectic. Same alternating B/D structure with Blanes-Moan `b(s)` (for B) and `a(s)`
(for D) coefficients over `nc` stages.

#### Commutator-Free Magnus (`integrate_cfm`)
Exponential time integration. Evaluates `n_stages` residuals at intermediate times, then applies
exponential updates with coefficients `s_coeffs` / `e_coeffs` via `cfm%compute_exponential_update`.

### Divergence cleaning (`impose_div_free`)

Called at the end of every time step. Dispatches based on `div_corr_var` and CT flags:

| Condition | Action |
|-----------|--------|
| `constrained_transport_D` and `div_corr_var == POISS` | `impose_ct_correction(ivar=1)` (D field) |
| `constrained_transport_B` and `div_corr_var == POISS` | `impose_ct_correction(ivar=4)` (B field) |

#### `impose_ct_correction(ivar)`
Poisson-based divergence correction via FLAIL multigrid:
1. Compute `divF = ∇·q(ivar:ivar+2)`.
2. Solve `−∇²φ = divF` using `compute_smoothing_multigrid` (with tolerance `flail%tolerance`,
   up to `flail%iterations` outer iterations, `iterations_init/fine/coarse` V-cycle sweeps).
3. Correct: `q(ivar+v-1) += (∇φ)_v` for `v = 1, 2, 3`.

The `divergence` array (size `(7,:,:,:,:)`) is used as scratch buffer.

### Coil current initialization

Rectangular coils are initialized by `set_rectangular_coil_x/y/z` (called once in
`set_initial_conditions`). Each routine:

1. Constructs a smooth vector potential **A** for the four wire segments using error-function
   transitions (`erf_function`) and tangential window functions (`tangential_window`).
2. Takes the discrete curl `J = ∇×A` via the dispatched `compute_curl` pointer.
3. Calls `compute_coil_current_density_flux` to calibrate amplitude against the requested current.
   For the z-normal coil this also applies a Poisson correction to suppress residual divergence.

The amplitude calibration integrates `J_z` (or `J_y`) over the coil cross-section and applies a
multiplicative correction factor.

At runtime, `compute_coils_current` injects the current density into the field `q(VAR_JX/Y/Z)` at
each time step using a smooth start transient:

$$g(t) = 10\!\left(\frac{t}{t_d}\right)^3 - 15\!\left(\frac{t}{t_d}\right)^4 + 6\!\left(\frac{t}{t_d}\right)^5, \quad t < t_d$$

followed by a step to steady-state or AC amplitude using branchless arithmetic on the activation
weight `w_c_`.

### Time step computation (`compute_dt`)

CFL condition with light-speed wave propagation:

$$\Delta t = \text{CFL} \cdot \frac{\Delta x_{\min}}{c}$$

where $c = \texttt{evmax}$ and $\Delta x_{\min}$ is the minimum cell width across all blocks
(OpenMP reduction). The result is reduced to the global minimum via `MPI_ALLREDUCE(MIN)`.

### Energy diagnostics

#### `compute_energy`
Accumulates field energy at each time step:

$$E_\mathbf{D} = \tfrac{1}{2} \sum_{b,i,j,k} |\mathbf{D}|^2, \qquad
E_\mathbf{B} = \tfrac{1}{2} \sum_{b,i,j,k} |\mathbf{B}|^2$$

Both sums are reduced via `MPI_ALLREDUCE(SUM)` and appended to `energy_D(it)` / `energy_B(it)`.

#### `compute_energy_error`
At the end of the simulation computes normalized RMS energy drift:

$$\text{RMS}_E = \sqrt{\frac{1}{N} \sum_{n=1}^{N} \left(\frac{|E_n - E_0|}{|E_0|}\right)^2}$$

separately for D and B fields. Printed to stdout at the end of `simulate`.

### Auxiliary fields (`compute_auxiliary_fields`)

Conditionally computed (flags from `io` object) before each save:

| Flag | Computed quantities |
|------|---------------------|
| `save_divergence_fields` | `∇·D`, `∇·B`, `∇·J` → `divergence(1:3,:,:,:,:)` |
| `save_curl_fields` | `∇×D`, `∇×B`, `∇×J` → `curl(1:9,:,:,:,:)` |
