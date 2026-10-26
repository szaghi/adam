# Input reference

A PRISM run reads one INI file: `mpirun -np N exe/adam_prism_cpu input.ini` (or `adam_prism_fnl`). This page lists
every section and key PRISM reads, with its type, whether it is required, its default, the accepted values and what
happens on an invalid one. The sections the library reads for every application (`[IO]`, `[grid]`, `[amr]`,
`[amr_marker_N]`, `[runge_kutta]`, `[weno]`, `[fdv]`, `[solids]`, `[slices]`, the forest manifest) are in the
[library input reference](/reference/library-input), with PRISM's use of each; PRISM's additional rules for them are
[below](#library-sections-in-prism). The [library page](/reference/library-input#parsing-rules) also gives the parsing
rules and the legend of the **Req.** column.

**Read the [hazards](#hazards) first**: they collect the combinations PRISM refuses, and the ones it accepts but
handles differently than the key names suggest.

## Program entry and read order

`adam_prism_cpu [file]` / `adam_prism_fnl [file]` (default `input.ini`). A file with a `[forest]` section and
`realms_number >= 1` is a [forest manifest](/reference/library-input#forest-manifest-multi-realm); otherwise it is
itself a single-realm input.

A realm input is read in this order: `[IO]` → `[numerics]` → `[physics]` → `[runge_kutta] scheme` and the coil counts
(the memory budget) → the library sections (`[grid]`, `[amr]`, markers, `[weno]`, `[runge_kutta]`, `[fdv]`, `[solids]`,
`[slices]`) → `[bc_*]` → `[PIC]` and `[particle_injection]` (PIC models) → `[time]` → `[initial_conditions]` →
`[coils_input]` → `[grms]` → `[magnetic_field_at_center_domain]` → `[probes]` → `[leapfrog]` (PIC leapfrog) →
`[external_fields]` → `[fWLayer]` → `[PML]`. The backend then checks the combinations (time and space schemes, PML, PIC
current) and stops on the ones it does not implement.

**Strings and CRLF.** Save inputs with LF line endings: on nvfortran (the FNL build) a CRLF input stops at its first
numeric key, whatever the string keys do (see the [library parsing rules](/reference/library-input#parsing-rules)). On
gfortran a CRLF input runs as an LF one. The string keys of `[numerics]`, `[bc_*]`,
`[external_fields] external_fields_applied`, `[fWLayer]` and `[PML] PML_type` go through `strip_control`, and
`[initial_conditions] type` and `[IO] save_fields` blank their control characters.

**Units.** SI, except under the non-dimensional models (`AdimEM`, `AdimPIC`), where lengths are divided by `L0`, fields
by `B0` and `D0 = B0 c0 ε0`, times by `T0 = L0/c0`; the input stays in SI and PRISM converts it.

---

## [IO]

The library keys are in the [library reference](/reference/library-input#io); PRISM uses `save_memory_status`,
`save_curl_fields`, `save_divergence_fields` and the seam div(B) monitor `seam_divB_tol` / `seam_divB_error`, which FLUME
does not. PRISM adds one key.

| Key | Type | Req. | Default | Accepted / invalid | Meaning |
|-----|------|------|---------|--------------------|---------|
| `save_fields` | string | no | `all` | `all`, `em`; other values are fatal | `all`: every field in the field files. `em`: D and B only (production runs); the `save_residual/curl/divergence_fields` flags still add their arrays. |

A cadence `<= 0` disables that output, end of run included ([#97](https://github.com/szaghi/adam/issues/97)).

---

## [numerics]

Every key is required; an unknown value is fatal.

| Key | Accepted values | Meaning |
|-----|-----------------|---------|
| `scheme_time` | `runge_kutta`, `leapfrog`, `blanes_moan`, `commutator_free_magnus` (also UPPER and Capitalised) | Time integration. **Only `runge_kutta` runs**: the other three are refused at initialization on both backends ([#55](https://github.com/szaghi/adam/issues/55)). |
| `scheme_space` | `fd_centered`, `fv_centered`, `weno` (also UPPER and Capitalised) | Spatial scheme. `weno` runs on the CPU backend only (the FNL backend refuses it). |
| `reconstruction_variables` | `conservative`, `characteristics` | Variables the WENO path reconstructs. |
| `constrained_transport` | `no`, `D`, `B`, `DB` (case variants accepted) | Which GLM scalars to carry: φ for D, ψ for B. See the [divergence-correction matrix](#divergence-correction). |
| `divergence_correction` | `no`, `hyperbolic`, `poisson` | `hyperbolic`: the GLM cleaning of the scalars selected by `constrained_transport` (required: `D`, `B` or `DB`). `poisson`: **refused** at load. `no`: no correction; a `constrained_transport` other than `no` is then ignored with a warning. |

---

## [physics]

| Key | Type | Req. | Default | Accepted / invalid | Meaning |
|-----|------|------|---------|--------------------|---------|
| `physical_model` | string | yes | — | `EM` (`em`, `Electromagnetic`, `electromagnetic`, `Maxwell`, `maxwell`), `PIC` (`pic`, `ParticleInCell`, `particleincell`), `AdimEM` (`adimem`, `AdimensionalizedElectromagnetic`, `adimensionalelectromagnetic`), `AdimPIC` (`adimpic`, `adim_PIC`, `adim_pic`); other values are fatal | Maxwell's equations, with particles (`PIC`), in SI or non-dimensional form. A PIC model adds the charge density `rho` to the state. |
| `L0` | real | cond. (`AdimEM`, `AdimPIC`) | — | `> 0`, otherwise fatal | Reference length (m). |
| `B0` | real | cond. (`AdimEM`, `AdimPIC`) | — | `> 0`, otherwise fatal | Reference magnetic field (T). |
| `chi` | real | cond. (`divergence_correction = hyperbolic`) | — | `> 0`, otherwise fatal | GLM propagation speed, in units of c: the fastest wave is `max(1, chi) c0`. |
| `alpha` | real | no | 0 | `>= 0`, otherwise fatal | GLM damping rate (1/s, or 1/T0 non-dimensional); `> 0` also caps the time step at `1/alpha`. Read only with `hyperbolic`. |

---

## Library sections in PRISM

| Section | PRISM rule |
|---------|------------|
| `[grid]` | `ngc >= fdv_order/2` (fatal); with `scheme_space = weno` also `ngc >= [weno] S` (fatal). Extents are divided by `L0` in a non-dimensional model. |
| `[amr]` | Refinement at initialization only, during the `[initial_conditions] amr_iterations` passes; `frequency` is read and unused, `regrid_prolongation` is unused. |
| `[amr_marker_N]` | Mode 1 `primitive-box` refines; mode 3 (`AMR_TV`) refines on the total variation of the coils' current density (`tol`, `delta_type`, `delta_fine`, `delta_coarse`). `stl` is fatal. **Mode 1 `solid`, mode 2, mode 4 and mode 5 do nothing, silently.** |
| `[runge_kutta]` | EM on CPU: every scheme, Yoshida included. EM on FNL: Yoshida is refused. PIC: `runge-kutta-ssp-22/33/54` only. The PML, the charge-conserving PIC current and the staged path (2:1 faces, forests) need an SSP scheme (the PIC current one of SSP-22/33/54). |
| `[weno]` | Used only with `scheme_space = weno` (CPU). `weights` (`js`, `si`) is used; ROR is not implemented. |
| `[fdv]` | `fdv_scheme` (`FD`, `FV`) selects PRISM's operators and the coil stamping; `fdv_order` the order of the centred operators, which the t0 elliptic solves also use (they accept 2, 4, 6, 8). |
| `[solids]` | Loaded and passed to the Runge-Kutta stages; no committed PRISM input uses solids. |

---

## [bc_x_min], [bc_x_max], [bc_y_min], [bc_y_max], [bc_z_min], [bc_z_max]

One key per section, required.

| Key | Accepted values | Meaning |
|-----|-----------------|---------|
| `type` | `extrapolation`, `Neumann`, `Dirichlet`, `Silver_Muller`, `periodic`, `PEC` (or `pec`); case-sensitive, other values are fatal. `radiative` is **refused** ([#57](https://github.com/szaghi/adam/issues/57)). | Face boundary condition. |

- **`periodic` is the library's periodicity** (`BC_PERIODIC = -1`): the tree builds periodic neighbours and the ghost
  exchange fills the periodic ghosts, across blocks and ranks. The two faces of an axis must be both periodic or
  neither (fatal).
- **`Dirichlet` writes zero** into the ghosts; there is no key for a prescribed value.
- Edge and corner ghost rows receive no boundary condition.
- `Silver_Muller` is spelled with an underscore.
- The [elliptic map](#elliptic-boundary-map) gives the boundary the t0 solves see.

---

## [time]

Every key is required.

| Key | Type | Meaning |
|-----|------|---------|
| `it_max` | int | `> 0`: stop after `it_max` steps. `<= 0`: run to `time_max`. |
| `time_max` | real | End time (s; divided by `T0` in a non-dimensional model). |
| `CFL` | real | CFL number on the fastest wave, `max(1, chi) c0` (upper-case key). The PML (`C_gamma`) and the GLM damping (`alpha`) may cap the step further. |

---

## [initial_conditions]

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `amr_iterations` | int | yes | negative values are clamped to 0 | Initial-condition + AMR passes: each sets the initial condition and runs the markers. |
| `type` | string | yes | see below; other values are fatal (CRLF-safe) | Initial field. |

| `type` | Extra keys (all required, real unless stated) | Field |
|--------|-----------------------------------------------|-------|
| `vacuum` | — | Every variable zero. |
| `plane_wave` | `kx`, `ky`, `kz` (direction components), `lambda` (wavelength), `B0` (amplitude, T) | A plane wave of wavelength `lambda` along `(kx, ky, kz)`, D = B0 c0 ε0. |
| `gaussian_pulse` | `B0`, `k0` (carrier wavenumber), `r0` (centre along the direction), `sigma` (`> 0`, width), `pulse_direction` (string) | A Gaussian-enveloped wave packet. |
| `gaussian_derivative_pulse` | `B0`, `r0`, `sigma` (`> 0`), `pulse_direction` (string) | A single-cycle pulse with no zero-frequency content. |
| `rmf_field` | `[external_fields]` `RMF_frequency`, `RMF_B_amplitude`, `RMF_rotation_axis` (`X`, `Y`, `Z`, either case; other values are fatal) | A rotating magnetic field. |
| `uniform_field` | `B_x`, `B_y`, `B_z` | A uniform magnetic field. |
| `riemann-problem` | — | **Refused**: its regions were read but never applied ([#56](https://github.com/szaghi/adam/issues/56)). |
| `magnetic_nozzle`, `rmf_magnetic_nozzle` | — | **Refused**: use `[external_fields] external_fields_applied` for a nozzle field (also refused, see there). |

`pulse_direction` accepts `+x`, `-x`, `+y`, `-y`, `+z`, `-z`; another value is fatal when the initial condition is set.
Under the non-dimensional models the pulse lengths are given in SI and converted.

---

## [coils_input] and [coils_input_coil_N]

Coil current sources, stamped into the current density J as a Gaussian distribution around the coil path.

### [coils_input]

| Key | Type | Req. | Meaning |
|-----|------|------|---------|
| `rectangular_coils_number` | int | yes | Rectangular loops. |
| `circular_coils_number` | int | yes | Circular loops. |
| `solenoid_coils_number` | int | yes | Solenoids. |
| `helicon_coils_number` | int | yes | Helicon antennas (polylines). |
| `time_delay` | real | cond. (any coil) | Ramp-up time: the current is multiplied by `g(s) = 10s³ − 15s⁴ + 6s⁵`, `s = t/time_delay` clamped to [0, 1], so it grows smoothly from zero at t0 to full at `time_delay`. With `time_delay <= 0` the coils carry full current at t0, and the t0 magnetic field is solved for (see `[PIC] magnetic_field_elliptic_initialization`). |

The four counts are always required, even when they are all zero.

### [coils_input_coil_N] (N = 1..total coils)

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `coil_type` | string | yes | `rectangular`, `circular`, `solenoid`, `helicon`; **another value is accepted silently and stamps no coil** | Geometry. |
| `current_type` | string | yes | `DC_current`, `AC_current`; **another value is accepted silently and reads no amplitude** | Waveform. |
| `sigma` | real | yes | — | Width of the Gaussian current distribution. |
| `x_center`, `y_center`, `z_center` | real | yes (every kind) | — | Coil centre. |
| `normal` | string | yes (every kind) | `+x`, `-x`, `+y`, `-y`, `+z`, `-z`; **another value stamps no coil, silently** | Coil plane normal (the axis of a solenoid). |
| `lx`, `ly` | real | `rectangular` | — | Loop sides. |
| `r_coil` | real | `circular`, `solenoid`, `helicon` | — | Radius. |
| `l_solenoid` | real | `solenoid` | — | Length. |
| `windings` | int | `solenoid` | — | Turns. |
| `N_points` | int | `helicon` | `<= 50`, otherwise fatal | Points of the antenna path. |
| `x_point_K`, `y_point_K`, `z_point_K` | real | `helicon`, K = 1..`N_points` | — | Antenna path points. |
| `Amplitude` | real | yes (`DC_current`, `AC_current`) | — | Current amplitude (A). |
| `Frequency`, `Phase` | real | `AC_current` | — | Frequency (Hz) and initial phase (**degrees**). |

---

## [external_fields]

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `external_fields_applied` | string | yes | `none` (`no`), `RMF`, `uniform_field` (`UniformField`), and the case variants; `magnetic_nozzle` and `rmf_and_magnetic_nozzle` are parsed and **refused** at initialization (no backend implements them); other values are fatal | External field added to the evolved one. |
| `RMF_frequency`, `RMF_B_amplitude` | real | `RMF` | — | Rotating field frequency (Hz) and amplitude (T). |
| `RMF_rotation_axis` | string | `RMF` | `X`, `Y`, `Z` (either case); **another value leaves the rotation indices undefined, silently** | Rotation axis. |
| `Uniform_field_D_amplitude` (or `Uniform_D_amplitude`) | real | `uniform_field` | — | Uniform D amplitude. |
| `Uniform_field_B_amplitude` (or `Uniform_B_amplitude`) | real | `uniform_field` | — | Uniform B amplitude. |
| `Uniform_field_direction` (or `Uniform_direction`) | string | `uniform_field` | `x`, `y`, `z` (either case); other values are fatal | Direction of the uniform field. |

For each pair, the first name is read first and the second is the fallback.

---

## [fWLayer]

A far-wave damping layer along chosen domain faces.

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `width` | real | yes (even when no layer is used) | `>= 0`, otherwise fatal | Physical layer width; the cells per face are `min(ni, ceiling(width/ds))`. `0` disables the layer. |
| `x_minus_layer`, `x_plus_layer`, `y_minus_layer`, `y_plus_layer`, `z_minus_layer`, `z_plus_layer` | string | yes | `yes`, `no` (case variants); other values are fatal | Faces carrying the layer. |

The removed key `C` (cells per face) is rejected by the removed-keys lint.

---

## [PML]

A perfectly matched layer, stored on the boundary blocks only. The section is optional: without `PML_type` the PML is
off.

| Key | Type | Req. | Default | Accepted / invalid | Meaning |
|-----|------|------|---------|--------------------|---------|
| `PML_type` | string | no | none | `no`/`none`, `CLASSIC` (`PML`, `STANDARD`), `CLASSIC_DIRECT` (`DIRECT`, `CLASSIC-DIRECT`), `CFS` (`CPML`, `CFS_PML`), `BERMUDEZ`, with case variants; other values are fatal | PML formulation. |
| `width` | real | yes (PML on) | 0 | `> 0`, otherwise fatal | Physical layer width. |
| `C_gamma` (or `c_gamma`) | real | no | 0.1 | `> 0`, otherwise fatal | Time-step cap `dt <= C_gamma / gamma_max`. |
| `x_minus_layer` … `z_plus_layer` | logical | yes (PML on) | — | at least one `.true.`, otherwise fatal | Faces carrying the PML. |
| `gamma_max` | real | `CLASSIC`, `CLASSIC_DIRECT`, `CFS` | — | `>= 0`, otherwise fatal | Maximum damping. |
| `gamma_exponent` | real | `CLASSIC`, `CLASSIC_DIRECT`, `BERMUDEZ` | — | `> 0`, otherwise fatal | Exponent of the damping profile. |
| `alfa_max` | real | `CFS` | — | `>= 0`, otherwise fatal | CFS frequency shift (spelled `alfa`). |
| `k_max` | real | `CFS` | 1 | `> 0`, otherwise fatal | CFS stretching. |
| `beta` | real | `BERMUDEZ` | — | `>= 0`, otherwise fatal | Bermúdez profile coefficient. |

See the [PML hazards](#pml).

---

## [PIC] (PIC models only)

| Key | Type | Req. | Default | Accepted / invalid | Meaning |
|-----|------|------|---------|--------------------|---------|
| `problem_type` | string | yes | — | `plasma`, `single_particle`; **another value is accepted silently and injects no particle** | Problem kind. |
| `plasma_domain` | string | `plasma` | — | `uniform` (`all`, `full`, `uniform_domain`), `cilinder` (`uniform_cilinder`), `slab` (`uniform_slab`), `box` (`uniform_box`), with case variants; other values are fatal | Region the plasma fills. Note the spelling **`cilinder`**. |
| `cilinder_radius`, `cilinder_length`, `cilinder_x_center`, `cilinder_y_center`, `cilinder_z_center` | real | `cilinder` | — | — | Cylinder. |
| `cilinder_axis` | char | `cilinder` | — | `x`, `y`, `z`; **another value leaves the particle positions undefined, silently** | Cylinder axis. |
| `slab_normal` | string | `slab` | — | `x`, `y`, `z` (either case); other values are fatal | Slab normal. |
| `slab_thickness`, `slab_center` | real | `slab` | — | thickness `> 0`, otherwise fatal | Slab along its normal. |
| `box_length_x`, `box_length_y`, `box_length_z`, `box_x_center`, `box_y_center`, `box_z_center` | real | `box` | — | — | Box. |
| `plasma_density` | real | `plasma` | — | — | Macro-particles per unit volume: the particle count is `nint(plasma_density × volume)` of the plasma region. |
| `neutral_fraction` | real | `plasma` | — | — | Fraction of neutrals; the rest is split evenly into ions and electrons. |
| `initialization` | string | no | `coherent` | `coherent` only; other values are fatal | Field initialization. |
| `electric_field_elliptic_initialization` | logical | no | `.true.` | logical | Solve D(G) φ = −ρ/ε at t0 and correct D, so div D = ρ to round-off. |
| `magnetic_field_elliptic_initialization` | logical | no | `.true.` | logical | Solve D(G) A = −μ J at t0 and add curl A to B, when J ≠ 0 at t0. |
| `elliptic_correction` | logical | — | — | `.true.` is fatal ([#84](https://github.com/szaghi/adam/issues/84) E3d); the lint rejects the key | Removed. |
| `particle_weighting_model` | string | yes | — | `NGP`, `CIC`, `TSC`, `cubic` (`3rd`), `quartic` (`4th`), `quintic` (`5th`), `sextic` (`6th`), `gaussian`, with case variants; other values are fatal | Charge deposition. |
| `current_weighting_model` | string | yes | — | the same list, plus `conserving` (`charge_conserving`); other values are fatal | Current deposition. `conserving` needs `fd_centered`, `runge_kutta` and SSP-22/33/54 (fatal otherwise). |
| `field_weighting_model` | string | yes | — | `0D` … `6D`, `gaussian`, `londrillo`, `totorica` (`esirkepov_Q2`), `esirkepov_Q5`, `esirkepov_Q6`, `esirkepov_Q8`, `esirkepov_Q10`, with case variants; other values are fatal | Field gather to the particles. |
| `current_conserving_solver` | string | no | `londrillo` | `londrillo`, `totorica` (`esirkepov_Q2`), `esirkepov_Q5/Q6/Q8/Q10`; the old names `esirkepov`, `direct`/`line`, `esirkepov-modified` stop with the new name; other values are fatal. Ignored unless `current_weighting_model = conserving`. | Matched current of the conserving deposit. With periodic faces `londrillo` stops at the first step (CPU). |
| `esirkepov_tail_cleanup` | logical | no | `.false.` | logical | Tail cleanup of the Esirkepov current; stops at the first step with periodic faces (CPU). |
| `filter_deposition` | logical | no | `.false.` | logical | Filter the deposited source. |
| `sigma`, `cutoff_sigma` | real | cond. (any Gaussian model) | — | `> 0`, otherwise fatal | Gaussian width and cut-off (in widths). |
| `gaussian_support_cells` | int | no | automatic | — | Gaussian support in cells; automatic from the final cell size when absent. |
| `scheme_time` | string | yes | — | `runge_kutta`, `leapfrog` (case variants); **another value warns and falls back to `runge_kutta`** | Particle push. `leapfrog` is refused at initialization. |

The particle closure at the domain faces is not an input key: it follows the field `[bc_*] type` (periodic faces wrap
the particles, every other kind keeps the legacy closure, where a particle leaving the domain stops the run).

---

## [particle_injection] (PIC models only)

For `problem_type = plasma`:

| Key | Type | Req. | Default | Accepted / invalid | Meaning |
|-----|------|------|---------|--------------------|---------|
| `random_seed` | int | no | processor-dependent | any | Fixed seed, for a reproducible plasma. |
| `space_random_number_generator` | string | yes | — | `random`, `layered` (case variants); other values are fatal | Position sampling. |
| `space_pairing` | string | yes | — | `yes`, `no`; other values are fatal | Electrons placed on the ions' positions. |
| `ions_velocity_distribution`, `electrons_velocity_distribution`, `neutrals_velocity_distribution` | string | yes | — | `uniform_maxwellian` (`maxwellian`, `uniform`), `non_uniform_maxwellian`, `drift` (`pure_drift`), with case variants; other values are fatal | Velocity distribution per species. |
| `ionic_temperature`, `electronic_temperature`, `neutrals_temperature` | real | cond. (that species `uniform_maxwellian`) | — | — | Temperature (K). |
| `ionic_temperature_x/_y/_z`, `electronic_temperature_x/_y/_z`, `neutrals_temperature_x/_y/_z` | real | cond. (that species `non_uniform_maxwellian`) | — | — | Temperature per axis (K). |
| `velocity_random_number_generator` | string | cond. (any Maxwellian species) | — | `random`, `layered`; other values are fatal | Velocity sampling. |
| `velocity_pairing` | string | cond. (any Maxwellian species) | — | `yes`, `no`; other values are fatal. `yes` needs the same distribution, temperatures and drifts for ions and electrons (fatal otherwise). | Electrons take the ions' velocities. |
| `ions_drift_x/_y/_z`, `electrons_drift_x/_y/_z`, `neutrals_drift_x/_y/_z` | real | yes | — | — | Drift velocity per species (m/s). |
| `v_av_correction` | string | yes | — | `yes`, `no`; other values are fatal | Remove the sampled mean velocity. |

For `problem_type = single_particle`: `x_position`, `y_position`, `z_position`, `charge`, `mass`, `v_drift_x`,
`v_drift_y`, `v_drift_z`, all real and required.

---

## [leapfrog] (PIC leapfrog only)

`is_filtered` (logical), `nu`, `alpha` (real), all required, read only with `[PIC] scheme_time = leapfrog`, which is then
refused at initialization: the section has no working use.

---

## Diagnostics

### [grms]

| Key | Type | Req. | Default | Accepted / invalid | Meaning |
|-----|------|------|---------|--------------------|---------|
| `save_history` | logical | no | `.false.` | logical | Write the RMS-field history. |
| `history_save` | int | no | 10 | — | Its cadence in steps. |
| `region_type` | string | no | `all` | `all`, `cylinder`; other values are fatal | Region of the RMS. |
| `center_x`, `center_y`, `center_z`, `axis_x`, `axis_y`, `axis_z`, `radius`, `length` | real | `cylinder` | — | axis non-zero, radius and length `> 0`, otherwise fatal | Cylinder (the axis is normalised). |

### [magnetic_field_at_center_domain]

| Key | Type | Req. | Default | Meaning |
|-----|------|------|---------|---------|
| `save_history` | logical | no | `.false.` | Write the history of B at the domain centre. |
| `history_save` | int | no | 10 | Its cadence in steps. |

### [probes] and [probe_N]

Optional: without a `[probes]` section there are no probes. Each probe samples the whole state at a point by trilinear
interpolation; one file per probe, `<output_basename>-probe_NN.dat`.

| Section | Key | Type | Req. | Accepted / invalid | Meaning |
|---------|-----|------|------|--------------------|---------|
| `[probes]` | `probes_number` | int | yes (section present) | `>= 0`, otherwise fatal; `0`: no probes | Number of probes. |
| | `save_every` | int | yes (probes on) | `> 0`, otherwise fatal | Cadence in steps; the last step is always saved. |
| `[probe_N]` | `x`, `y`, `z` | real | yes (the section too) | inside the domain, otherwise fatal | Probe point, in the grid's units (divided by `L0` non-dimensional). |

---

## Hazards

### Divergence correction

The GLM scalars φ (for D) and ψ (for B) exist only under `hyperbolic`, and their slots follow `physics`
([#11](https://github.com/szaghi/adam/issues/11), [#61](https://github.com/szaghi/adam/issues/61)):

| `constrained_transport` | `divergence_correction` | Result |
|---|---|---|
| `no` | `no` | `nv = 9` |
| `D`, `B`, `DB` | `no` | constrained transport ignored with a warning, `nv = 9` |
| `D` or `B` | `hyperbolic` | `nv = 10` (φ or ψ in slot 7) |
| `DB` | `hyperbolic` | `nv = 11` (φ in 7, ψ in 8) |
| `no` | `hyperbolic` | **refused** at initialization |
| any | `poisson` | **refused** at initialization |

PIC adds `rho`, one more slot. At a 2:1 AMR seam div(B) is not kept at round-off: the seam injects an O(h^p) source
that converges under refinement but grows in time at fixed h ([#29](https://github.com/szaghi/adam/issues/29)); `[IO]
seam_divB_tol` makes the growth visible.

### PML

- **CPU:** `CLASSIC`, `CLASSIC_DIRECT`, `BERMUDEZ`, `CFS`, on `fd_centered`, with no divergence correction and no
  constrained transport, with an SSP Runge-Kutta scheme; anything else is fatal.
- **FNL:** `fd_centered`, `runge_kutta` with an SSP scheme, no hyperbolic correction, no constrained transport; anything
  else is fatal.
- **`CLASSIC_DIRECT` is not a Runge-Kutta PML:** it skips the auxiliary-variable integration and applies direct
  per-face damping. It is a testing instrument, not the production path.
- Enabling the section with every face flag `.false.` is fatal; each `PML_type` has its own required parameters, and a
  missing one is fatal, not defaulted.

### Elliptic boundary map

The t0 solves (`[PIC] electric_field_elliptic_initialization`, `magnetic_field_elliptic_initialization`) run on FLAIL's
FFT solver with PRISM's own operator D(G) of order `fdv_order`. Each face kind maps to an elliptic kind, then to a FLAIL
face kind:

| `[bc_*] type` | Elliptic kind | FLAIL face kind |
|---|---|---|
| `periodic` | periodic | `FLAIL_BC_PERIODIC` |
| `Dirichlet` | Dirichlet | `FLAIL_BC_ZERO_GHOSTS` (every ghost layer zero) |
| `PEC` | PEC | `FLAIL_BC_DIRICHLET` (φ, tangential A), `FLAIL_BC_NEUMANN` (normal A) |
| `Neumann`, `Silver_Muller` | exact/open | `FLAIL_BC_FREE_GHOSTS` (the free-space potential of the source in the ghosts, summed over the whole domain) |
| `extrapolation` | — | **fatal**, at the first t0 solve |

The map is total and stops on what it cannot map, but **only when a solve first touches the face**: an input with an
`extrapolation` face passes initialization and stops at t0 if a PIC field or a t0 current needs a solve. A null axis
takes the periodic kind.

### Refused combinations

| Input | Where it stops |
|---|---|
| `[numerics] scheme_time = leapfrog`, `blanes_moan`, `commutator_free_magnus` | initialization, both backends ([#55](https://github.com/szaghi/adam/issues/55)) |
| `[PIC] scheme_time = leapfrog` | initialization, both backends |
| `[numerics] divergence_correction = poisson` | load |
| `hyperbolic` without constrained transport | load |
| `[bc_*] type = radiative` | load ([#57](https://github.com/szaghi/adam/issues/57)) |
| `[initial_conditions] type = riemann-problem`, `magnetic_nozzle`, `rmf_magnetic_nozzle` | load ([#56](https://github.com/szaghi/adam/issues/56)) |
| `[external_fields] external_fields_applied = magnetic_nozzle`, `rmf_and_magnetic_nozzle` | initialization |
| `[numerics] scheme_space = weno` on FNL | initialization |
| `runge-kutta-yoshida` on FNL; a non-SSP scheme with PIC, the PML or the staged path | initialization |
| `[PIC] current_weighting_model = conserving` on the CPU backend with more than one block or one rank | the first step (the CPU charge-conserving deposit is single-block, single-rank) |
| A refined grid that enters a t0 elliptic solve | the solve ([#84](https://github.com/szaghi/adam/issues/84) E7: the wide operator across 2:1 faces is not defined) |
| `fdv_order = 10` with a t0 elliptic solve | the solve (FLAIL accepts 2, 4, 6, 8) |
