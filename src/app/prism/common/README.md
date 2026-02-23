# PRISM Common

The `common/` directory contains all backend-independent modules shared by the CPU and FNL backends of PRISM. Every type described here is aggregated into `prism_common_object`, which serves as the base type for both backends.

## Source Layout

| File | Module | Description |
|------|--------|-------------|
| `adam_prism_parameters.F90` | `adam_prism_parameters` | Physical constants, eigenvalue and eigenvector matrices |
| `adam_prism_physics_object.F90` | `adam_prism_physics_object` | Physical model selection, variable counts, eigenvector dispatch |
| `adam_prism_numerics_object.F90` | `adam_prism_numerics_object` | Temporal/spatial scheme selection, divergence correction config |
| `adam_prism_common_object.F90` | `adam_prism_common_object` | Base aggregate type for all backends |
| `adam_prism_bc_object.F90` | `adam_prism_bc_object` | Boundary conditions handler |
| `adam_prism_ic_object.F90` | `adam_prism_ic_object` | Initial conditions handler |
| `adam_prism_io_object.F90` | `adam_prism_io_object` | I/O configuration and output routines |
| `adam_prism_time_object.F90` | `adam_prism_time_object` | Time integration control |
| `adam_prism_coil_object.F90` | `adam_prism_coil_object` | Current-carrying coil current density sources |
| `adam_prism_fWLayer_object.F90` | `adam_prism_fWLayer_object` | Far-wave absorbing boundary layer |
| `adam_prism_external_fields_object.F90` | `adam_prism_external_fields_object` | Externally prescribed electromagnetic fields |
| `adam_prism_rk_bc_object.F90` | `adam_prism_rk_bc_object` | RK sub-step integrator for boundary conditions |
| `adam_prism_pic_object.F90` | `adam_prism_pic_object` | Particle-in-Cell configuration and particle-grid weighting |
| `adam_prism_particle_injection_object.F90` | `adam_prism_particle_injection_object` | Initial particle space and velocity distribution |
| `adam_prism_leapfrog_pic_object.F90` | `adam_prism_leapfrog_pic_object` | Leapfrog (Boris) PIC particle time integrator |
| `adam_prism_rk_pic_object.F90` | `adam_prism_rk_pic_object` | Runge-Kutta PIC particle time integrator |
| `adam_prism_riemann_library.F90` | `adam_prism_riemann_library` | Maxwell Riemann solvers (LLF, HLL) and convective flux routines |
| `adam_prism_common_library.F90` | `adam_prism_common_library` | Barrel re-export of all common modules |

---

## `adam_prism_parameters` — Physical Constants

Module-level parameters (all `save`):

| Constant | Value | Description |
|----------|-------|-------------|
| `NV_MAX` | `11` | Maximum supported variable count (static array sizing) |
| `MU0` | `1.256637e-6` | Vacuum magnetic permeability (H/m) |
| `EPS0` | `8.854187e-12` | Vacuum permittivity (F/m) |
| `C0` | `1/√(MU0·EPS0)` | Vacuum speed of light (m/s) |
| `E_CHARGE` | `-1.6e-19` | Electron charge (C) |
| `E_MASS` | `9.11e-31` | Electron mass (kg) |
| `K_B` | `1.380649e-23` | Boltzmann constant (J/K) |
| `PI` | `3.141592653589793` | π |

Derived constants: `MU0_SQ = √MU0`, `MU0_SQ_I2 = 1/(2√MU0)`, `EPS0_SQ = √EPS0`, `EPS0_SQ_I2 = 1/(2√EPS0)`.

### Eigensystem Arrays

The Maxwell system has wave speeds $\pm c_0$ and two degenerate zero speeds. The module provides target arrays used by `prism_physics_object` to set eigenvector pointers:

| Symbol | Shape | Description |
|--------|-------|-------------|
| `EV(6)` | `[0, 0, c0, c0, −c0, −c0]` | Eigenvalues of the 6-variable system |
| `ER(6,6,3)` | `real, target` | Right eigenvectors for x, y, z directions |
| `EL(6,6,3)` | `real, target` | Left eigenvectors for x, y, z directions |
| `IEV(6)` | all ones | Identity eigenvalues (conservative reconstruction) |
| `IERL(6,6,3)` | identity | Identity right/left eigenvectors (conservative) |
| `IEV_D(7)`, `IERL_D(7,7,3)` | 7-variable | With D divergence-cleaning scalar φ |
| `IEV_B(7)`, `IERL_B(7,7,3)` | 7-variable | With B divergence-cleaning scalar ψ |
| `IEV_D_B(8)`, `IERL_D_B(8,8,3)` | 8-variable | With both φ and ψ |

---

## `adam_prism_physics_object` — Physical Model

The `prism_physics_object` reads the physical model from the `[physics]` INI section and assembles the variable layout and eigensystem pointers.

### INI Keys (`[physics]`)

| Key | Values | Description |
|-----|--------|-------------|
| `physical_model` | `electromagnetic` / `PIC` | Select EM or Particle-in-Cell model |
| `chi` | real > 0 | Divergence-cleaning propagation speed coefficient (hyperbolic only); `evmax = chi · c₀` |

Accepted aliases for `electromagnetic`: `Maxwell`, `EM`, `em`, `maxwell`. For `PIC`: `pic`, `ParticleInCell`.

### Data Members

| Member | Default | Description |
|--------|---------|-------------|
| `physical_model` | — | Selected model string |
| `nv_c` | `6` | Conservative variable count (D, B — grows by `nv_cl`) |
| `nv_s` | `3` | Source variable count (J, always zero flux) |
| `nv_cl` | `0` | Divergence-cleaning variable count (0, 1, or 2) |
| `nv_pic` | `0` | PIC charge density variable count (0 or 1) |
| `nv` | `nv_c + nv_s + nv_cl + nv_pic` | Total variable count |
| `var_Jx/Jy/Jz` | 7/8/9 (no cleaning) | Q-vector indices of J components; shift with cleaning |
| `chi` | `0` | Cleaning speed coefficient |
| `evmax` | `c₀` or `chi·c₀` | Maximum eigenvalue (maximum signal speed) |
| `erw`, `elw` | pointer | Right/left eigenvectors for WENO characteristic reconstruction |

### Variable Index Constants

```fortran
VAR_DX = 1,  VAR_DY = 2,  VAR_DZ = 3
VAR_BX = 4,  VAR_BY = 5,  VAR_BZ = 6
! J indices shift depending on cleaning: var_Jx = 7, 8, or 9
```

### Eigenvector Dispatch

At `initialize`, the pointers `erw` and `elw` are set based on the combination of `reconstruction_vars` and `div_corr_var`:

| `reconstruction_vars` | Div correction | `erw` / `elw` |
|-----------------------|----------------|----------------|
| `CONSERVATIVE` | any | `IERL` / `IERL_D` / `IERL_B` / `IERL_D_B` (identity) |
| `CHARACTERISTICS` | none / POISSON | `ER(6,6,3)` / `EL(6,6,3)` |
| `CHARACTERISTICS` | HYPERBOLIC, D only | `ER_D(7,7,3)` / `EL_D(7,7,3)` |
| `CHARACTERISTICS` | HYPERBOLIC, B only | `ER_B(7,7,3)` / `EL_B(7,7,3)` |
| `CHARACTERISTICS` | HYPERBOLIC, D+B | `ER_D_B(8,8,3)` / `EL_D_B(8,8,3)` |

For hyperbolic divergence cleaning the enlarged eigensystem uses the cleaning speed `ch = chi·c₀`: eigenvalues include `0, ±ch, ±c₀`.

---

## `adam_prism_numerics_object` — Numerical Schemes

Configured from the `[numerics]` INI section.

### INI Keys (`[numerics]`)

| Key | Values | Description |
|-----|--------|-------------|
| `scheme_time` | `RUNGE_KUTTA` / `LEAPFROG` / `BLANES_MOAN` / `COMMUTATOR_FREE_MAGNUS` | Temporal integration scheme |
| `scheme_space` | `WENO` / `FD_CENTERED` / `FV_CENTERED` | Spatial discretisation scheme |
| `fdv_order` | integer | Order of centered FD/FV schemes |
| `reconstruction_variables` | `CONSERVATIVE` / `CHARACTERISTICS` | WENO reconstruction basis |
| `constrained_transport` | `NO` / `D` / `B` / `DB` | Enable Constrained Transport on D and/or B |
| `divergence_correction` | `POISSON` / `HYPERBOLIC` | Divergence error control strategy (absent = none) |

The half-stencil width is derived as `fdv_half_stencil = fdv_order / 2`. Per-derivative stencil widths in `fdv_half_stencils(6)` are incremented progressively for higher-order derivatives.

---

## `adam_prism_common_object` — Base Aggregate Type

`prism_common_object` aggregates every framework and PRISM sub-object and owns all simulation field arrays. Both backends (`prism_cpu_object`, `prism_fnl_object`) extend this type.

### Data Members

**Framework objects:**

| Member | Type | Description |
|--------|------|-------------|
| `mpih` | `mpih_object` | MPI handler |
| `adam` | `adam_object` | Top-level ADAM framework handle |
| `field` | `field_object` pointer | Field storage and metrics |
| `grid` | `grid_object` pointer | Grid geometry |
| `amr` | `amr_object` | AMR marker handler |
| `ib` | `ib_object` | Immersed Boundary handler |
| `slices` | `slices_object` | Slice output handler |
| `leapfrog` | `leapfrog_object` | Leapfrog field integrator |
| `blanesmoan` | `blanesmoan_object` | Blanes-Moan field integrator |
| `cfm` | `cfm_object` | Commutator-Free Magnus integrator |
| `rk` | `rk_object` | Runge-Kutta field integrator |
| `weno` | `weno_object` | WENO reconstructor |
| `flail` | `flail_object` | Linear algebra methods (Poisson solver) |

**PRISM sub-objects:**

| Member | Type | Description |
|--------|------|-------------|
| `io` | `prism_io_object` | I/O configuration |
| `numerics` | `prism_numerics_object` | Scheme selection |
| `physics` | `prism_physics_object` | Physical model and variable layout |
| `ic` | `prism_ic_object` | Initial conditions |
| `bc` | `prism_bc_object` | Boundary conditions |
| `rk_bc` | `prism_rk_bc_object` | RK sub-step BC integrator |
| `time` | `prism_time_object` | Time control |
| `fWLayer` | `prism_fWLayer_object` | Far-wave absorbing layer |
| `coil` | `prism_coil_object` | Coil current density sources |
| `external_fields` | `prism_external_fields_object` | Prescribed external fields |
| `pic` | `prism_pic_object` | PIC configuration and weighting |
| `particle_injection` | `prism_particle_injection_object` | Initial particle distribution |
| `leapfrog_pic` | `prism_leapfrog_pic_object` | Leapfrog PIC integrator |
| `rk_pic` | `prism_rk_pic_object` | RK PIC integrator |

**Convenience pointers** (bound to `field`/`physics` members at init):

`ngc`, `ni`, `nj`, `nk`, `nb`, `blocks_number`, `nv`, `nv_c`, `nv_s`, `nv_cl`

**Field arrays** (allocated in `allocate_common`):

| Array | Shape | Description |
|-------|-------|-------------|
| `q` | `(nv, 1-ngc:ni+ngc, …, nb)` | Conservative state vector (D, B, J, ±φ/ψ, ±ρ) |
| `dq` | same as `q` | Residual right-hand side |
| `curl` | same as `q` | Curl of field components (optional output) |
| `divergence` | same as `q` | Divergence of field components (optional output) |
| `q_pic` | `(8, particle_number)` | PIC particle state (PIC model only) |
| `pic_fields` | `(8, particle_number)` | Field values at particle locations (PIC model only) |

**Energy diagnostics:**

| Member | Description |
|--------|-------------|
| `energy_D(:)` | Time history of electric field energy |
| `energy_B(:)` | Time history of magnetic field energy |
| `rms_energy_error_D/B` | RMS relative energy error |

### `initialize_common` Sequence

The 34-step initialization sequence called by both backends:

1. `mpih%initialize` — MPI setup
2. `io%initialize(filename)` — load and parse INI file
3. `bc%initialize` — read boundary condition types
4. `numerics%initialize` — read scheme selections and divergence correction
5. `physics%initialize` — configure physical model and eigensystem (uses `numerics` outputs)
6. `adam%grid%initialize` — initialize grid from INI (with BC type array)
7. `adam%compute_blocks_number` — estimate blocks from available memory
8. `adam%initialize` — allocate tree, maps, and field storage for `nv` variables
9. Associate convenience pointers (`field`, `grid`, `ngc`, `ni`, …, `nv_cl`)
10. `pic%initialize` — PIC configuration and particle count (PIC model only)
11. `particle_injection%initialize` — particle distribution config (PIC model only)
12. `adam%refine_uniform` — uniform mesh refinement to configured levels
13. `adam%prune` — prune tree regions defined by `ijkl_prune`
14. `amr%initialize` — AMR marker configuration
15. `field%compute_metrics` — cell geometry and metrics
16. `time%initialize` — read `it_max`, `time_max`, `CFL`
17. `ic%initialize` — read initial conditions type and parameters
18. `fWLayer%initialize` — compute far-wave layer damping function
19. `coil%initialize` — compute coil current density fields
20. `ib%initialize` — immersed boundary setup
21. `slices%initialize` — slice output configuration
22. `external_fields%initialize` — external field type and parameters
23. `rk%initialize` — RK coefficients (if `RUNGE_KUTTA`)
24. `rk_bc%initialize` — RK-BC stage arrays (if `RUNGE_KUTTA`)
25. `leapfrog%initialize` — leapfrog weights (if `LEAPFROG`)
26. `blanesmoan%initialize` — Blanes-Moan weights (if `BLANES_MOAN`)
27. `cfm%initialize` — CFM weights (if `COMMUTATOR_FREE_MAGNUS`)
28. `weno%initialize` — WENO stencil configuration (if `WENO`)
29. `leapfrog_pic%initialize` — Boris leapfrog PIC (if PIC + LEAPFROG)
30. `rk_pic%initialize` — RK PIC (if PIC + RUNGE_KUTTA)
31. `flail%initialize` — linear algebra methods
32. `check_ngc_number` — verify `ngc ≥ weno%S` and `ngc ≥ fdv_half_stencil`; stop on failure
33. `allocate_common` — allocate `dq`, `divergence`, `curl`, and PIC arrays
34. `io_initialize` — register auxiliary output fields with the I/O layer

**`io_initialize`** registers auxiliary fields conditionally:

| Condition | Registered field | Names |
|-----------|-----------------|-------|
| `io%save_residual_fields` | `dq` | `res_Dx` … `res_Jz` [± `res_ph`/`res_ps`/`res_rh`] |
| `io%save_divergence_fields` | `divergence` | `div_D`, `div_B`, `div_J`, `fWL_x/y/z`, … |
| `coil%total_coils_number > 0` | `coil%j_vec`, `coil%coil_flag` | `j_vec_1/2/3`, `f_Gauss`, `coil_flag` |
| `io%save_curl_fields` | `curl` | `curlD_x/y/z`, `curlB_x/y/z`, `curlJ_x/y/z` |

### Methods

| Method | Description |
|--------|-------------|
| `allocate_common` | Allocate `dq`, `divergence`, `curl`, and PIC arrays |
| `initialize_common` | Full 34-step initialization (called by backend `initialize`) |
| `load_restart_files` | Load `q` and time from HDF5 restart files |
| `save_energy_error` | Write energy error history if save frequency triggered |
| `save_restart_files` | Save HDF5 restart files and an XH5F snapshot |
| `save_xh5f` | Save simulation data in XDMF+HDF5 format with variable names |

---

## `adam_prism_bc_object` — Boundary Conditions

Reads one BC type per face from six INI sections.

### INI Sections

`[bc_x_min]`, `[bc_x_max]`, `[bc_y_min]`, `[bc_y_max]`, `[bc_z_min]`, `[bc_z_max]`

Each section requires the key `type`:

| `type` value | Integer constant | Description |
|--------------|-----------------|-------------|
| `extrapolation` | `BC_EXTRAPOLATION = 1` | Zero-gradient extrapolation |
| `Neumann` | `BC_NEUMANN = 2` | Neumann (zero normal derivative) |
| `Dirichlet` | `BC_DIRICHLET = 3` | Dirichlet (prescribed field values) |
| `Silver_Muller` | `BC_Silver_Muller = 4` | Silver-Müller absorbing BC |
| `periodic` | `BC_PERIOD = 5` | Periodic |
| `radiative` | `BC_radiative = 6` | Radiative (first-order absorbing) |

The object stores `bc_type(6)` (integer codes) and `q(9,6)` (field values for Dirichlet BCs, one set per face).

---

## `adam_prism_ic_object` — Initial Conditions

Configured from the `[initial_conditions]` INI section.

### IC Types

| `type` value | Description |
|-------------|-------------|
| `vacuum` | All fields set to zero |
| `riemann-problem` | Piecewise constant regions; each region in `[initial_conditions_region_N]` |
| `plane_wave` | Sinusoidal electromagnetic plane wave |
| `rmf_field` | Rotating Magnetic Field (reads from `[external_fields]`) |
| `uniform_field` | Uniform background B and/or D |
| `magnetic_nozzle` | Magnetic nozzle configuration (parameters TBD) |
| `rmf_magnetic_nozzle` | Combined RMF + magnetic nozzle |

### INI Keys (`[initial_conditions]`)

| Key | Description |
|-----|-------------|
| `amr_iterations` | Number of AMR iterations during IC application |
| `type` | IC type string (see table above) |
| `regions_number` | Number of regions (Riemann problem only) |
| `kx`, `ky`, `kz` | Wave number direction cosines (plane wave) |
| `lambda` | Wavelength (plane wave) |
| `B0` | Background magnetic field amplitude (plane wave) |
| `B_x`, `B_y`, `B_z` | Uniform background B components |

**Region sections** (`riemann-problem`): `[initial_conditions_region_N]` with keys `Dx`, `Dy`, `Dz`, `Bx`, `By`, `Bz`, `Jx`, `Jy`, `Jz`, `emin_x/y/z`, `emax_x/y/z`.

---

## `adam_prism_io_object` — I/O Handler

Configured from the `[IO]` INI section.

### INI Keys (`[IO]`)

| Key | Default | Description |
|-----|---------|-------------|
| `output_basename` | — | Prefix for all output files |
| `it_save` | `100` | Field output save frequency (iterations) |
| `restart` | `.false.` | Enable restart from existing files |
| `restart_basename` | — | Prefix for restart files |
| `restart_save` | `100` | Restart save frequency |
| `residuals_save` | `10` | Residuals log save frequency |
| `save_memory_status` | `.false.` | Log memory usage during allocations |
| `save_residual_fields` | `.false.` | Include `dq` in field output |
| `save_curl_fields` | `.false.` | Include curl field in output |
| `save_divergence_fields` | `.false.` | Include divergence field in output |
| `save_gradient_fields` | `.false.` | Include gradient field in output |
| `save_laplacian_fields` | `.false.` | Include Laplacian field in output |

### Output Files

| File | Columns | Trigger |
|------|---------|---------|
| `<basename>-energy_error.dat` | `it  blocks  time  error_D  error_B  rms_D  rms_B` | Every `energy_error_save` iterations |
| `<basename>-residuals.dat` | `it  time  blocks  rq1 … rqN` | Every `residuals_save` iterations |

---

## `adam_prism_time_object` — Time Handler

Configured from the `[time]` INI section.

### INI Keys (`[time]`)

| Key | Default | Description |
|-----|---------|-------------|
| `it_max` | `-1` | Maximum iteration count (`−1` = terminate by time) |
| `time_max` | `1.0` | Maximum physical time |
| `CFL` | `0.3` | CFL number for time-step control |

### Termination Logic (`is_to_save`)

| Condition | Meaning |
|-----------|---------|
| `mod(it, it_save) == 0` | Periodic save |
| `it == it_max` | Final iteration (when `it_max > 0`) |
| `it_max ≤ 0 .and. time ≥ time_max` | Physical time limit reached |
| `it_max > 0 .and. it ≥ it_max` | Iteration limit reached |

`print_progress` prints `it`, time-step `dt`, current `time`, and a percentage completion based on whichever termination criterion is active.

---

## `adam_prism_coil_object` — Coil Sources

Models current-carrying coils as volumetric current density **J** sources.

### INI Section (`[coils_input]`)

| Key | Values | Description |
|-----|--------|-------------|
| `rectangular_coils_number` | integer | Number of rectangular coils |
| `circular_coils_number` | integer | Number of circular coils |

Per-coil configuration (section `[coil_N]`):

| Key | Values | Description |
|-----|--------|-------------|
| `coil_type` | `rectangular` / `circular` | Loop geometry |
| `current_type` | `DC_current` / `AC_current` | Current waveform |
| `normal` | `+x` / `-x` / `+y` / `-y` / `+z` / `-z` | Coil plane normal |
| `A` | real | Current amplitude (A) |
| `f` | real | Frequency (Hz) — AC only |
| `phase` | real | Initial phase (rad) — AC only |
| `x_center`, `y_center`, `z_center` | real | Coil centre coordinates |
| `lx`, `ly` | real | Rectangular loop dimensions |
| `r_coil` | real | Circular loop radius |
| `sigma` | real | Gaussian current distribution width |

### Internal Arrays

| Array | Shape | Description |
|-------|-------|-------------|
| `J_vec` | `(total_coils_number, 4, ni, nj, nk, nb)` | Per-coil current density vector and Gaussian weight |
| `coil_flag` | `(ni, nj, nk, nb)` | Integer flag: which coil passes through each cell |

`J_vec(:,1:3,…)` holds current density direction components; `J_vec(:,4,…)` holds the Gaussian distribution value `f_Gauss`.

---

## `adam_prism_fWLayer_object` — Far-Wave Absorbing Layer

Implements a sponge-like damping layer at domain boundaries to absorb outgoing electromagnetic waves, following Barbas notation.

### INI Section (`[fWLayer]`)

| Key | Description |
|-----|-------------|
| `C` | Layer width in cells (0 = layer disabled) |
| `layer_xm`, `layer_xp`, … | Boolean flags for each of the six sides (−x, +x, −y, +y, −z, +z) |

### Damping Function

The damping function `f(3, ni, nj, nk, nb)` is initialized once at startup. For a layer of width `C` cells:

$$f_i = \begin{cases} \frac{1}{150}(-7C^2 + 255C + 250) & C < 40 \\ 25 & C \geq 40 \end{cases}$$

The values decay smoothly from the interior (where `f = 1`) toward the boundary. The three components of `f` correspond to the x, y, z sides of the layer, enabling correct treatment of corner and edge regions.

The `divergence` array is initialized with `fWLayer%f` components at indices 4–6, making the layer damping coefficients available as output fields (`fWL_x`, `fWL_y`, `fWL_z`).

---

## `adam_prism_external_fields_object` — External Fields

Prescribed electromagnetic fields added to (and subtracted from) the state vector around each time step.

### INI Section (`[external_fields]`)

| Key | Values | Description |
|-----|--------|-------------|
| `ef_type` | `RMF` / `magnetic_nozzle` / `RMF_and_magnetic_nozzle` / `none` | External field type |
| `RMF_frequency` | real | Rotating magnetic field frequency (Hz) |
| `RMF_B_amplitude` | real | Rotating magnetic field amplitude (T) |
| `RMF_rotation_axis` | `X` / `Y` / `Z` | Rotation axis |

At initialization, the `add_external_fields` and `sub_external_fields` procedure pointers are dispatched to the corresponding implementation (`add_external_fields_rmf`, etc.) based on `ef_type`.

---

## `adam_prism_rk_bc_object` — RK-BC Integrator

Stores the per-stage boundary condition state during Runge-Kutta sub-stepping.

### Data Members

| Member | Description |
|--------|-------------|
| `q_bc_rk(:,:,:,:,:,:)` | Per-stage BC field values |
| `dq_bc_rk(:,:,:,:,:)` | Per-stage BC residuals |
| `ark`, `brk`, `crk` | Low-storage RK coefficients |
| `alph`, `beta`, `gamm` | SSP RK coefficients |
| `ssa`, `ssb` | Symplectic-splitting RK coefficients |
| `nrk` | Number of RK stages |

### Methods

| Method | Description |
|--------|-------------|
| `initialize` | Mirror the parent `rk_object` coefficient set |
| `initialize_stages` | Zero the stage arrays |
| `assign_stage` | Copy current q to the stage buffer |
| `compute_stage` | Advance the BC state for one RK sub-step |

---

## `adam_prism_pic_object` — Particle-in-Cell

Manages macro-particle configuration, grid-particle mapping, and weighting.

### INI Section (`[PIC]`)

| Key | Values | Description |
|-----|--------|-------------|
| `problem_type` | `plasma` / `single_particle` | Simulation scenario |
| `plasma_density` | real | Initial plasma number density (m⁻³) — plasma only |
| `neutral_fraction` | real | Fraction of neutral particles — plasma only |
| `particle_weighting_model` | `CIC` / `NGP` / `TSC` | Charge deposition scheme |
| `current_weighting_model` | `CIC` / `NGP` / `TSC` | Current deposition scheme |
| `field_weighting_model` | `0D` / `1D` | Field interpolation order at particle positions |
| `scheme_time` | `LEAPFROG` / `RUNGE_KUTTA` | Particle time integrator |

### Particle Count

For `plasma`: `particle_number = plasma_density × domain_volume`. The population is split into ions, electrons, and neutrals according to `neutral_fraction`. For `single_particle`: `particle_number = 1`.

### Weighting Models

| Scheme | Method | Description |
|--------|--------|-------------|
| `NGP` | Nearest Grid Point | Particle charge/current deposited to the containing cell only |
| `CIC` | Cloud-in-Cell | Trilinear interpolation to 3×3×3 neighbourhood |
| `TSC` | Triangular Shaped Cloud | Quadratic weight function over 3×3×3 neighbourhood |

Field interpolation models:

| Model | Description |
|-------|-------------|
| `0D` | Zero-order: field value at the nearest grid point |
| `1D` | First-order: piecewise-linear field interpolation at particle position |

At initialization, the three procedure pointers `particle_weighting`, `current_weighting`, and `field_weighting` are dispatched to the selected method implementations.

### `neighbour_list(4, particle_number)`

Stores the current grid indices `(block, i, j, k)` for each particle, updated each time step by `particle_cartesian_grid_index`.

---

## `adam_prism_particle_injection_object` — Particle Injection

Sets the initial spatial and velocity distributions for PIC particles.

### INI Section (`[particle_injection]`)

| Key | Values | Description |
|-----|--------|-------------|
| `space_distribution` | `Uniform_domain_space_distribution` / `Uniform_boxes_space_distribution` / `Uniform_cell_space_distribution` / `Space_random_number_generator` / `Space_layered_number_generator` | Initial spatial distribution |
| `velocity_distribution` | `Uniform_Maxwellian` / `Non_uniform_Maxwellian` / `Velocity_random_number_generator` / `Velocity_layered_number_generator` | Initial velocity distribution |
| `T_i`, `T_e`, `T_n` | real | Ion/electron/neutral temperature (uniform Maxwellian) |
| `T_i_x/y/z`, `T_e_x/y/z`, `T_n_x/y/z` | real | Per-axis temperatures (non-uniform Maxwellian) |
| `v_drift_x/y/z` | real | Drift velocity components |
| `box_number` | real | Number of boxes for charge-neutral injection |
| `space_pairing` | logical | Pair ions and electrons in space |
| `velocity_pairing` | logical | Pair ion/electron velocities |
| `v_av_correction` | logical | Correct mean velocity after sampling |
| `x_position`, `y_position`, `z_position` | real | Initial position (single particle) |
| `charge`, `mass` | real | Charge and mass (single particle) |

At initialization, the procedure pointers `particle_space_injection` and `particle_velocity_injection` are dispatched to the selected distribution subroutines.

---

## `adam_prism_leapfrog_pic_object` — Leapfrog PIC Integrator

Implements the Boris leapfrog particle pusher, integrating particle equations of motion staggered in time relative to the field update. The Boris method is energy-conserving and second-order accurate for the Lorentz force equation.

---

## `adam_prism_rk_pic_object` — Runge-Kutta PIC Integrator

Implements explicit Runge-Kutta time integration for PIC particle trajectories, sharing RK coefficients with the parent field integrator `rk_object`.

---

## `adam_prism_riemann_library` — Riemann Solvers

Provides the Maxwell convective flux routines and Riemann solvers used in all spatial schemes.

### Public Routines

| Routine | Description |
|---------|-------------|
| `compute_riemann_maxwell_llf` | LLF (Rusanov) solver |
| `compute_riemann_maxwell_hll` | HLL (Harten–Lax–van Leer) solver (internal) |
| `compute_convective_fluxes_maxwell` | Base flux for the 6-variable system |
| `compute_convective_fluxes_maxwell_div_d` | Flux extended with D-cleaning scalar φ |
| `compute_convective_fluxes_maxwell_div_b` | Flux extended with B-cleaning scalar ψ |
| `compute_convective_fluxes_maxwell_div_d_b` | Flux with both φ and ψ |
| `compute_eigenvalues_vector` | Maximum wave speed for time-step control |

### LLF Solver

$$F_{\text{LLF}}(q_L, q_R) = \tfrac{1}{2}\bigl[F(q_L) + F(q_R) - c_0(q_R - q_L)\bigr]$$

applied to the 6 conservative variables (D, B). The J source variables have zero flux and are not passed through the Riemann solver.

### HLL Solver

$$F_{\text{HLL}} = \frac{\lambda^+ F_L - \lambda^- F_R + \lambda^+ \lambda^- (q_R - q_L)}{\lambda^+ - \lambda^-}, \quad \lambda^+ = c_0,\; \lambda^- = -c_0$$

### `compute_convective_fluxes_interface`

Abstract interface for flux routines, allowing procedure pointer dispatch in backend implementations:

```fortran
subroutine compute_convective_fluxes_interface(sir, q, f, chi)
  real(R8P), intent(in)    :: sir(3)  ! directional unit increment (x, y, or z)
  real(R8P), intent(in)    :: q(1:)   ! state vector
  real(R8P), intent(inout) :: f(1:)   ! computed fluxes
  real(R8P), intent(in)    :: chi     ! divergence-cleaning speed coefficient
end subroutine
```

---

## `adam_prism_common_library` — Barrel Re-export

A single convenience module that re-exports all 17 common modules. Both backends use only this one `use` statement:

```fortran
use adam_prism_common_library
```

Modules re-exported: `adam_prism_parameters`, `adam_prism_physics_object`, `adam_prism_numerics_object`, `adam_prism_common_object`, `adam_prism_bc_object`, `adam_prism_ic_object`, `adam_prism_io_object`, `adam_prism_time_object`, `adam_prism_coil_object`, `adam_prism_fWLayer_object`, `adam_prism_external_fields_object`, `adam_prism_rk_bc_object`, `adam_prism_pic_object`, `adam_prism_particle_injection_object`, `adam_prism_leapfrog_pic_object`, `adam_prism_rk_pic_object`, `adam_prism_riemann_library`.

---

## License

PRISM is part of the ADAM framework, released under the [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.html) (LGPLv3).

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.
