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
| `adam_prism_pic_bc_object.F90` | `adam_prism_pic_bc_object` | Particle closure at the domain faces (periodic or legacy) |
| `adam_prism_pml_object.F90` | `adam_prism_pml_object` | Perfectly matched layer (`[PML]`) |
| `adam_prism_rk_pml_object.F90` | `adam_prism_rk_pml_object` | RK support for the PML auxiliary variables |
| `adam_prism_absorbing_layer_geometry.F90` | `adam_prism_absorbing_layer_geometry` | Layer geometry shared by fWLayer and PML |
| `adam_prism_probes_object.F90` | `adam_prism_probes_object` | Field probes (`[probes]`): time histories at fixed points |
| `adam_prism_grms_object.F90` | `adam_prism_grms_object` | Grms diagnostic (`[grms]`) |
| `adam_prism_magnetic_field_at_center_domain_object.F90` | `adam_prism_magnetic_field_at_center_domain_object` | Field-at-domain-centre diagnostic |
| `adam_prism_riemann_library.F90` | `adam_prism_riemann_library` | Maxwell Riemann solvers (LLF, HLL) and convective flux routines |
| `adam_prism_common_library.F90` | `adam_prism_common_library` | Barrel re-export of all common modules |

---

## `adam_prism_parameters` — Physical Constants

Module-level parameters (all `save`):

| Constant | Value | Description |
|----------|-------|-------------|
| `NV_MAX` | `12` | Maximum supported variable count (static array sizing); a larger `nv` is fatal at initialization |
| `MU0` | `1.256637e-6` | Vacuum magnetic permeability (H/m) |
| `EPS0` | `8.854187e-12` | Vacuum permittivity (F/m) |
| `C0` | `1/√(MU0·EPS0)` | Vacuum speed of light (m/s) |
| `E_CHARGE` | `-1.6e-19` | Electron charge (C) |
| `E_MASS` | `9.11e-31` | Electron mass (kg) |
| `ION_MASS` | `1.67262192595e-27` | Ion mass, H⁺ (proton, CODATA 2022) (kg): mass of the injected plasma ions |
| `NEUTRAL_MASS` | `ION_MASS + E_MASS` | Mass of the injected plasma neutrals (kg) |
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
| `physical_model` | `electromagnetic` / `AdimEM` / `PIC` / `AdimPIC` | Select the EM or Particle-in-Cell model, dimensional or non-dimensional |
| `L0`, `B0` | real | Reference length and field of the non-dimensional models (`AdimEM`, `AdimPIC` only) |
| `chi` | real > 0 | Divergence-cleaning propagation speed coefficient (hyperbolic only, required); `evmax = max(1, chi) · c₀` (`max(1, chi)` non-dimensional) |
| `alpha` | real ≥ 0 | Damping rate of the cleaning scalars, 1/s (hyperbolic only, optional, default 0) |

Accepted aliases for `electromagnetic`: `Electromagnetic`, `Maxwell`, `EM`, `em`, `maxwell`. For `PIC`: `pic`,
`ParticleInCell`, `particleincell`. For `AdimEM`: `adimem`, `AdimensionalizedElectromagnetic`,
`adimensionalelectromagnetic`. For `AdimPIC`: `adimpic`, `adim_PIC`, `adim_pic`. Any other value is fatal.

### Data Members

| Member | Default | Description |
|--------|---------|-------------|
| `physical_model` | — | Selected model string |
| `nv_c` | `6` | Conservative variable count (D, B — grows by `nv_cl`) |
| `nv_s` | `3` | Source variable count (J, always zero flux) |
| `nv_cl` | `0` | Divergence-cleaning variable count (0, 1, or 2) |
| `nv_pic` | `0` | PIC charge density variable count (0 or 1) |
| `nv` | `nv_c + nv_s + nv_cl + nv_pic` | Total variable count |
| `var_phi`, `var_psi` | `0` | Q-vector indices of the GLM scalars φ (D) and ψ (B), 0 when absent. The **single source** of the φ/ψ slots: `io_initialize` names the variables from them and stops if the names do not number `nv` (issue #61) |
| `var_Jx/Jy/Jz` | 7/8/9 (no cleaning) | Q-vector indices of J components; shift with cleaning |
| `chi` | `0` | Cleaning speed coefficient |
| `alpha` | `0` | Cleaning damping rate |
| `evmax` | `c₀` or `max(1, chi)·c₀` | Maximum eigenvalue (maximum signal speed) |
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
| `CHARACTERISTICS` | none | `ER(6,6,3)` / `EL(6,6,3)` |
| `CHARACTERISTICS` | HYPERBOLIC, D only | `ER_D(7,7,3)` / `EL_D(7,7,3)` |
| `CHARACTERISTICS` | HYPERBOLIC, B only | `ER_B(7,7,3)` / `EL_B(7,7,3)` |
| `CHARACTERISTICS` | HYPERBOLIC, D+B | `ER_D_B(8,8,3)` / `EL_D_B(8,8,3)` |
| any | HYPERBOLIC, no CT | none: `error_stop` (issue #61; `numerics` refuses the input first) |

For hyperbolic divergence cleaning the enlarged eigensystem uses the cleaning speed `ch = chi·c₀`: eigenvalues include `0, ±ch, ±c₀`.

---

## `adam_prism_numerics_object` — Numerical Schemes

Configured from the `[numerics]` INI section.

### INI Keys (`[numerics]`)

| Key | Values | Description |
|-----|--------|-------------|
| `scheme_time` | `RUNGE_KUTTA` / `LEAPFROG` / `BLANES_MOAN` / `COMMUTATOR_FREE_MAGNUS` | Temporal integration scheme. Only `RUNGE_KUTTA` works: the other three are parsed, then refused by both backends at init (issue #55) |
| `scheme_space` | `WENO` / `FD_CENTERED` / `FV_CENTERED` | Spatial discretisation scheme |
| `reconstruction_variables` | `CONSERVATIVE` / `CHARACTERISTICS` | WENO reconstruction basis |
| `constrained_transport` | `NO` / `D` / `B` / `DB` | Fields the hyperbolic cleaning acts on (ignored, with a warning, without it) |
| `divergence_correction` | `NO` / `HYPERBOLIC` | Divergence error control strategy. `HYPERBOLIC` requires `constrained_transport` ≠ `NO`; `POISSON` is refused (not implemented). Issue #61 |

Every key is required, and an unknown value stops at load naming the option and the accepted values (issue #57);
before, `scheme_time`, `scheme_space`, `reconstruction_variables`, `constrained_transport` and `divergence_correction`
warned and fell back to a default. The same holds for `[external_fields] external_fields_applied` (`none`, `RMF`,
`uniform_field`) and the `[fWLayer]` face flags (`yes`, `no`), which fell back silently.

The order of the centred FD/FV schemes is not a `[numerics]` key: `[fdv] fdv_scheme` and `fdv_order` are read by the
library's `realm_object` (`adam_realm_object.F90`). The half-stencil width is derived as `fdv_half_stencil = fdv_order / 2`. Per-derivative stencil widths in `fdv_half_stencils(6)` are incremented progressively for higher-order derivatives.

---

## `adam_prism_common_object` — Base Aggregate Type

`prism_common_object` extends the library's `realm_object` and owns the PRISM sub-objects and field arrays. Both
backends (`prism_cpu_object`, `prism_fnl_object`) extend this type.

### Data Members

**Framework objects** (components of `realm_object`): `io`, `amr`, `slices`, `blanesmoan`, `cfm`, `leapfrog`, `flail`
(the elliptic solver library of the t0 solves, issue #84), `weno`, `ib`, `rk`, and `adam` (`grid`, `tree`, `field`,
`maps`). MPI goes through the program-scope `mpih` singleton (`adam_mpih_global`).

**PRISM sub-objects** (`adam_prism_common_object.F90`):

| Member | Type | Description |
|--------|------|-------------|
| `bc` | `prism_bc_object` | Boundary conditions |
| `coil` | `prism_coil_object` | Coil current density sources |
| `external_fields` | `prism_external_fields_object` | Prescribed external fields |
| `fWLayer` | `prism_fWLayer_object` | Far-wave absorbing layer |
| `grms` | `prism_grms_object` | Grms diagnostic |
| `ic` | `prism_ic_object` | Initial conditions |
| `leapfrog_pic` | `prism_leapfrog_pic_object` | Leapfrog PIC integrator |
| `magnetic_field_at_center_domain` | `prism_magnetic_field_at_center_domain_object` | Field-at-domain-centre diagnostic |
| `numerics` | `prism_numerics_object` | Scheme selection |
| `particle_injection` | `prism_particle_injection_object` | Initial particle distribution |
| `physics` | `prism_physics_object` | Physical model and variable layout |
| `pic` | `prism_pic_object` | PIC configuration and weighting |
| `pml` | `prism_pml_object` | Perfectly matched layer |
| `probes` | `prism_probes_object` | Field probes |
| `rk_pic` | `prism_rk_pic_object` | RK PIC integrator |
| `rk_pml` | `prism_rk_pml_object` | RK PML auxiliary handler |
| `time` | `prism_time_object` | Time control |

**Field arrays** (allocated in `allocate_common`):

| Array | Shape | Description |
|-------|-------|-------------|
| `q` | `(nv, 1-ngc:ni+ngc, …, nb)` | Conservative state vector (D, B, ±φ/ψ, J, ±ρ) |
| `dq` | same as `q` | Residual right-hand side |
| `curl` | same as `q` | Curl of field components (optional output) |
| `divergence` | same as `q` | Divergence of field components (optional output) |
| `q_pic` | `(8, particle_number)` | PIC particle state (PIC model only) |
| `pic_fields` | `(6, particle_number)` | D and B at the particle locations (PIC model only) |

**Energy diagnostics:**

| Member | Description |
|--------|-------------|
| `energy_D(:)` | Time history of electric field energy |
| `energy_B(:)` | Time history of magnetic field energy |
| `rms_energy_error_D/B` | RMS relative energy error |

### `initialize` sequence

`prism_common_object%initialize` (`adam_prism_common_object.F90`), called by both backends:

1. `mpih%initialize`; `io%initialize(filename)` (load the INI file); `[IO] save_fields`
2. `numerics%initialize`, then `physics%initialize` (uses `numerics`); `nv > NV_MAX` is fatal
3. `fields_number` (given, or computed from the input)
4. `realm_object%initialize`: `adam` (grid, tree, maps, field), AMR markers, IB, slices, RK, WENO, the `[fdv]` keys
   (`L0` for the non-dimensional models); FLAIL is set up later, by each t0 solve
5. `bc%initialize`, `grid%set_bc_type`, then `tree%rebuild_neighborhood` and the ghost/BC maps (periodic faces need the
   BC types first)
6. PIC only: `pic%initialize`, `particle_injection%initialize`
7. `time`, `ic`, `coil`, `grms`, `magnetic_field_at_center_domain`, `probes`
8. PIC only: `leapfrog_pic%initialize` or `rk_pic%initialize`, by `[PIC] scheme_time`
9. `check_ngc_number` (`ngc ≥ weno%S` with WENO, `ngc ≥ fdv_half_stencil`), `allocate_common`
10. uniform refinement to `iu_ref_levels`; PIC Gaussian support cells; block face data
11. `io_initialize` — the variable names

The field integrators `leapfrog`, `blanesmoan` and `cfm` are never initialized (`realm_object%initialize` leaves them
commented out), which is why both backends refuse those schemes (issue #55).

**`io_initialize`** registers auxiliary fields conditionally:

| Condition | Registered field | Names |
|-----------|-----------------|-------|
| always | `q` | `Dx Dy Dz Bx By Bz` [`phi`] [`psi`] `Jx Jy Jz` [`rho`] |
| `io%save_residual_fields` | `dq` | `res_Dx` … `res_Bz` [`res_phi`] [`res_psi`] `res_Jx` … `res_Jz` [`res_rho`] |
| `io%save_divergence_fields` | `divergence` | `div_D`, `div_B`, `div_J`, `totvar`, `div_05`, `div_06` [`div_07`] [`div_08`] `div_09`…`div_11` [`div_12`] |
| `io%save_curl_fields` | `curl` | `curl_Dx` … `curl_Bz` [`curl_07`] [`curl_08`] `curl_Jx` … `curl_Jz` [`curl_12`] |
| coils present | `coil%j_vec` | `coil_NN_j_vec_1/2/3` per coil (`adam_prism_coil_object.F90`) |

The `q` names must number `physics%nv`, otherwise `io_initialize` stops (issue #61).

### Methods

| Method | Description |
|--------|-------------|
| `allocate_common` | Allocate `dq`, `divergence`, `curl`, and PIC arrays |
| `initialize` | The initialization sequence above (called by the backends) |
| `load_restart_files` | Load `q` and time from HDF5 restart files |
| `save_energy_error` | Write energy error history if save frequency triggered |
| `save_restart_files` | Save HDF5 restart files and an XH5F snapshot |
| `save_xh5f` | Save simulation data in XDMF+HDF5 format with variable names |

---

## `adam_prism_bc_object` — Boundary Conditions

Reads one BC type per face from six INI sections.

### INI Sections

`[bc_x_min]`, `[bc_x_max]`, `[bc_y_min]`, `[bc_y_max]`, `[bc_z_min]`, `[bc_z_max]`

Each section requires the key `type`; matching is case-sensitive and an unknown value is fatal. The two faces of an
axis must be both periodic or neither:

| `type` value | Integer constant | Description |
|--------------|-----------------|-------------|
| `extrapolation` | `BC_EXTRAPOLATION = 1` | Zero-gradient extrapolation |
| `Neumann` | `BC_NEUMANN = 2` | Neumann (zero normal derivative) |
| `Dirichlet` | `BC_DIRICHLET = 3` | Dirichlet (prescribed field values) |
| `Silver_Muller` | `BC_Silver_Muller = 4` | Silver-Müller absorbing BC |
| `periodic` | `BC_PERIODIC = -1` (library) | Periodic: library periodicity, the tree builds periodic neighbours |
| `PEC` (or `pec`) | `BC_PEC = 7` | Perfect electric conductor |
| `radiative` | — | **Refused at load**: not implemented on any backend (issue #57) |

The object stores `bc_type(6)` (integer codes) and `q(9,6)` (field values for Dirichlet BCs, one set per face).

---

## `adam_prism_ic_object` — Initial Conditions

Configured from the `[initial_conditions]` INI section.

### IC Types

| `type` value | Description |
|-------------|-------------|
| `vacuum` | All fields set to zero |
| `plane_wave` | Sinusoidal electromagnetic plane wave |
| `gaussian_pulse` | Gaussian-modulated plane wave packet (`B0`, `k0`, `r0`, `sigma`, `pulse_direction` = `±x`, `±y`, `±z`) |
| `gaussian_derivative_pulse` | Single-cycle pulse, derivative of a Gaussian vector potential: zero area, no ω = 0 content (`B0` = peak \|B\|, `r0`, `sigma` of the potential, `pulse_direction`; no `k0`) |
| `rmf_field` | Rotating Magnetic Field (reads `RMF_frequency`, `RMF_B_amplitude` and `RMF_rotation_axis` = `X`, `Y` or `Z` from `[external_fields]`; any other axis stops at load) |
| `uniform_field` | Uniform background B (`B_x`, `B_y`, `B_z`); D is set to zero |

Any other `type` stops at load, naming the accepted values ([issue #56](https://github.com/szaghi/adam/issues/56)). That
includes three names the loader still knows but `set_initial_conditions` never applied, so they used to start the run
from the zero field silently:
- `riemann-problem`: the regions were read and dropped (its apply code was commented out, then removed in 2025);
- `magnetic_nozzle` and `rmf_magnetic_nozzle`: for a nozzle field use `[external_fields] external_fields_applied`.

### INI Keys (`[initial_conditions]`)

| Key | Description |
|-----|-------------|
| `amr_iterations` | Number of AMR iterations during IC application |
| `type` | IC type string (see table above) |
| `kx`, `ky`, `kz` | Wave number direction cosines (plane wave) |
| `lambda` | Wavelength (plane wave) |
| `B0` | Background magnetic field amplitude (plane wave) |
| `B_x`, `B_y`, `B_z` | Uniform background B components |

The loader still reads the `riemann-problem` region sections (`[initial_conditions_region_N]`, keys `Dx..Jz`,
`emin_x/y/z`, `emax_x/y/z`; `emin_x` was never read before #56), so restoring that type needs its apply case and the removal of
its refusal.

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
| `solenoid_coils_number` | integer | Number of solenoids |
| `helicon_coils_number` | integer | Number of helicon antennas |
| `time_delay` | real | Ramp-up time (required when any coil exists): the current is multiplied by g(s) = 10s³ − 15s⁴ + 6s⁵, s = t/`time_delay` clamped to [0, 1], so it grows smoothly from zero at t0 to full at `time_delay`. With `time_delay ≤ 0` the coils carry full current at t0, and the t0 magnetic field is solved for (the t0 elliptic solve) |

Per-coil configuration, section `[coils_input_coil_N]` (N = 1..total coils), every key required:

| Key | Values | Description |
|-----|--------|-------------|
| `coil_type` | `rectangular` / `circular` / `solenoid` / `helicon` | Geometry |
| `current_type` | `DC_current` / `AC_current` | Current waveform |
| `sigma` | real | Gaussian current distribution width |
| `x_center`, `y_center`, `z_center` | real | Coil centre coordinates |
| `normal` | `+x` / `-x` / `+y` / `-y` / `+z` / `-z` | Coil plane normal |
| `lx`, `ly` | real | Rectangular loop dimensions (`rectangular`) |
| `r_coil` | real | Radius (`circular`, `solenoid`, `helicon`) |
| `l_solenoid`, `windings` | real, integer | Length and windings (`solenoid`) |
| `N_points`, `x_point_K`, `y_point_K`, `z_point_K` | integer (≤ 50), real | The antenna's points (`helicon`) |
| `Amplitude` | real | Current amplitude (A) |
| `Frequency`, `Phase` | real | Frequency (Hz) and initial phase (rad), `AC_current` only |

### Internal Arrays

| Array | Shape and description |
|-------|-------------|
| `J_vec` | `(3, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, nb, total_coils_number)`: per-coil current density, written to the field output as `coil_NN_j_vec_1/2/3` |

---

## `adam_prism_fWLayer_object` — Far-Wave Absorbing Layer

Implements a sponge-like damping layer at domain boundaries to absorb outgoing electromagnetic waves, following Barbas notation.

### INI Section (`[fWLayer]`)

| Key | Description |
|-----|-------------|
| `width` | Physical layer width (real, required even when no face is active); the per-face cell count is `min(ni, ceiling(width/ds))` |
| `x_minus_layer`, `x_plus_layer`, `y_minus_layer`, `y_plus_layer`, `z_minus_layer`, `z_plus_layer` | `yes` / `no` per face; any other value is fatal |

The former cell-count key `C` is removed: the lint `src/tests/lint/check-removed-input-keys.sh` rejects it in committed inputs.

### Damping Function

The damping factor is **not** stored: `compute_fwl_factor` evaluates it per cell, on demand, from the
cell's distance to the boundary and the face's profile geometry. It is declared `!$acc routine seq` /
`!$omp declare target`, so one function serves both the CPU and the device backend.

For a face spanning `profile_cells` cells over a physical `profile_extent`, a cell at `center_distance`
from the boundary-cell center gets the profile-*shape* coefficient

$$f_i = \begin{cases} \frac{1}{150}(-7C_r^2 + 255C_r + 250) & C_r < 40 \\ 25 & C_r \geq 40 \end{cases}$$

with `C_r = profile_cells`, and from it the damping factor

$$f = \frac{1}{f_i}\log_{10}\!\left[r\left(10^{f_i} - 1\right) + 1\right], \qquad
r = \mathrm{clamp}\!\left(\frac{\texttt{center\_distance}}{\texttt{profile\_extent}},\, 0,\, 1\right)$$

A face with `profile_cells <= 0` or `profile_extent <= 0` returns `f = 1` (no damping). `profile_cells`
and `profile_extent` are derived per face at initialization from the requested physical `width` (see
`compute_absorbing_face_range`); on a uniform mesh `profile_extent` equals `C*ds`, which is what the
former cell-count `C` input expressed directly.

`7e39318b` removed the `f` array from `prism_fWLayer_object` together with its HDF5 save block, so
`fWL_x`, `fWL_y` and `fWL_z` are no longer emitted as output fields.

---

## `adam_prism_external_fields_object` — External Fields

Prescribed electromagnetic fields added to (and subtracted from) the state vector around each time step.

### INI Section (`[external_fields]`)

| Key | Values | Description |
|-----|--------|-------------|
| `external_fields_applied` | `none` / `RMF` / `uniform_field` (also `magnetic_nozzle`, `RMF_and_magnetic_nozzle`: parsed, then refused at initialization, not implemented) | External field type; an unknown value is fatal |
| `RMF_frequency` | real | Rotating magnetic field frequency (Hz), `RMF` |
| `RMF_B_amplitude` | real | Rotating magnetic field amplitude (T), `RMF` |
| `RMF_rotation_axis` | `X` / `Y` / `Z` | Rotation axis, `RMF` |
| `Uniform_field_D_amplitude`, `Uniform_field_B_amplitude` | real | Uniform field amplitudes, `uniform_field` (aliases `Uniform_D_amplitude`, `Uniform_B_amplitude`) |
| `Uniform_field_direction` | `x` / `y` / `z` | Uniform field direction, `uniform_field` (alias `Uniform_direction`); any other value is fatal |

At initialization, the `add_external_fields` and `sub_external_fields` procedure pointers are dispatched to the corresponding implementation (`add_external_fields_rmf`, `add_external_fields_uniform`) based on `external_fields_applied`.

---

## `adam_prism_rk_bc_object` — RK-BC Integrator

Stores the per-stage boundary condition state during Runge-Kutta sub-stepping. **Not used**: the `rk_bc` component
and its initialization are commented out in `adam_prism_common_object.F90`.

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
| `particle_weighting_model` | `NGP` / `CIC` / `TSC` / `cubic` / `quartic` / `quintic` / `sextic` / `Gaussian` | Charge deposition scheme |
| `current_weighting_model` | the same, or `conserving` (charge-conserving, `current_conserving_solver`) | Current deposition scheme |
| `field_weighting_model` | `0D` / `1D` / `2D` / `3D` / `4D` / `5D` / `6D` / `Gaussian` | Field interpolation at particle positions |
| `scheme_time` | `LEAPFROG` / `RUNGE_KUTTA` | Particle time integrator; only `RUNGE_KUTTA` with SSP 22/33/54 runs, `LEAPFROG` stops at initialization on both backends |

The section also reads `plasma_domain` (with the `cilinder_*`, `box_*` and `slab_*` geometry keys), `sigma` and
`cutoff_sigma` (Gaussian kinds), `gaussian_support_cells`, `filter_deposition`, `current_conserving_solver`,
`esirkepov_tail_cleanup`, `initialization`, and the t0 solve switches `electric_field_elliptic_initialization` and
`magnetic_field_elliptic_initialization`. `elliptic_correction` is removed (issue #84 E3d): `.true.` stops the run.

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

Plasma problems (`[PIC] problem_type = plasma`); the spatial distribution follows `[PIC] plasma_domain`
(`space_distribution` is no longer read):

| Key | Values | Description |
|-----|--------|-------------|
| `random_seed` | integer | Optional: a fixed seed makes the injected plasma reproducible (default: the processor-dependent seed) |
| `space_random_number_generator` | `random` / `layered` | Position sampling |
| `space_pairing` | `yes` / `no` | Electron i at the position of ion i (ions sampled as usual): zero charge density at t0 |
| `ions_velocity_distribution`, `electrons_velocity_distribution`, `neutrals_velocity_distribution` | `uniform_Maxwellian` (also `Maxwellian`, `uniform`) / `non_uniform_Maxwellian` / `drift` | Velocity distribution per species |
| `ionic_temperature`, `electronic_temperature`, `neutrals_temperature` | real | Temperatures (uniform Maxwellian) |
| `ionic_temperature_x/y/z`, `electronic_temperature_x/y/z`, `neutrals_temperature_x/y/z` | real | Per-axis temperatures (non-uniform Maxwellian) |
| `velocity_random_number_generator` | `random` / `layered` | Velocity sampling, read when a Maxwellian is selected |
| `velocity_pairing` | `yes` / `no` | Electron i with the velocity of ion i (ions sampled as usual); needs the same distribution, temperatures and drifts for ions and electrons, otherwise fatal; read only when a Maxwellian is selected |
| `ions_drift_x/y/z`, `electrons_drift_x/y/z`, `neutrals_drift_x/y/z` | real | Per-species drift velocity components (all required) |
| `v_av_correction` | `yes` / `no` | Correct the mean velocity after sampling |

Single-particle problems: `x_position`, `y_position`, `z_position`, `v_drift_x/y/z`, `charge`, `mass`.

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
