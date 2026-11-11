# Input reference

A FLUME run reads one INI file: `mpirun -np N exe/adam_flume_cpu input.ini` (or `adam_flume_fnl`). This page lists
every section and key the application reads, with its type, whether it is required, its default, the accepted values
and what happens on an invalid one. The two complete inputs at the end are verification cases; the
[verification gallery](./verification/) shows their results.

## Program entry and read order

`adam_flume_cpu [file]` / `adam_flume_fnl [file]` (default `input.ini`). If the file has a `[forest]` section with
`realms_number >= 1`, it is a **forest manifest** and every `[realm.N]` names a per-realm INI of the form documented
below; otherwise the file is itself a single-realm input.

A realm input is read in this order: `[IO]` → `[reference]` (if present it converts the file to code units, then
`[IO]` is read again) → `[numerics]` → `[physics]`
(+ `[eos]`, `[mhd]`, `[gravity]`) → `[runge_kutta].scheme` (memory budget) → `[grid]` → `[amr]` (tree, seam fill, markers) → `[solids]`
→ `[slices]` → `[runge_kutta]` → `[weno]` → `[fdv]` → `[bc_*]` → `[time]` →
`[initial_conditions]` → `[diagnostics]` → `[IO].save_auxiliary_fields`. After reading, FLUME applies
cross-section checks (slices, WENO scheme, Riemann solver per model, ghost cells, even block cells for AMR).

## Units

Any consistent unit system works: the equations carry no dimensionless number
([units and scaling](./models#units-and-scaling)), so SI, cgs and units of order 1 give the same run. Magnetic fields
are always $\mathbf{B}_{SI}/\sqrt{\mu_0}$ (a Gaussian field divided by $\sqrt{4\pi}$). With the classic WENO
weights the data should be of order 1: either write the input in such units, add a
[`[reference]`](#reference-optional-dimensional-input) section that converts a physical input, or set
`[weno] weights = si`.

## Parsing rules

The [library input reference](/reference/library-input#parsing-rules) gives the rules of the INI parser (case,
comments, missing keys, numbers, logicals, CRLF, the removed-keys lint). A Windows (CRLF) input reads
as an LF one on both backends.

## Legend for the tables

- **Req.**: `yes` means a missing key causes `error_stop`. `cond.` means it is required only under the stated condition.
  `no` means optional (the default is applied). `no check` means the error is ignored; see the notes for what happens.
- **Used**: `yes` means FLUME's behaviour depends on the value. `no` means the key is read but FLUME never uses it;
  the key usually must still be present, because the library reader requires it.

---

## [IO]

The library reads `[IO]` ([library reference](/reference/library-input#io)); FLUME adds one key. `[IO]` is loaded
twice (FLUME common, then the realm), and both loads are identical.

| Key | Type | Req. | Default | Accepted / invalid | Used | Meaning |
|-----|------|------|---------|--------------------|------|---------|
| `save_auxiliary_fields` | logical | yes | — | logical | yes | Also write the auxiliary fields (`rho,u,v,w,p,T,H,a`, plus `s` with dual energy, `s1..s4` with passive scalars, `Bx,By,Bz` for MHD). For MHD it also writes the derived fields `pt, beta, bmag, divb`. |

Of the library keys FLUME uses `output_basename`, `it_save`, `restart`, `restart_basename`, `restart_save`,
`residuals_save` and `save_residual_fields`. The field save writes the first step (`it = 0`, not on restart) and the
last step whenever `it_save > 0`; the last step is also a restart save. The other `[IO]` keys are read and unused by
FLUME (its div(B) history cadence is `[diagnostics]`, its div(B) monitor `[mhd] divb_tol`).

### Output formats (no format key exists)

FLUME has no output-format selector. It always writes these files:

| Output | File name | Controlled by |
|--------|-----------|---------------|
| Fields: XDMF + HDF5 (XH5F), ghost cells included | `<output_basename>-<it:9 digits>.xdmf` + `...-procNNNNNN.h5` | `it_save`, `save_residual_fields`, `save_auxiliary_fields` |
| Restart | `<restart_basename>.time`, `.tnd`, `-procNNNNNN.fbd` + an XH5F snapshot `<restart_basename>.xdmf/.h5` | `restart_save` |
| Residuals history (ASCII, Tecplot header) | `<output_basename>-residuals.dat` | `residuals_save` |
| Conservation history (ASCII) | `<output_basename>-conservation_history.dat` | `[diagnostics].conservation_history_save` |
| Turbulence history (ASCII, optional) | `<output_basename>-turbulence_history.dat` | `[diagnostics].turbulence_history_save` |
| Wall-model history (ASCII, optional) | `<output_basename>-wall_<face>_history.dat` | `[diagnostics].wall_history_save` |
| div(B) history (MHD only) | `<output_basename>-divb_history.dat` | `[diagnostics].conservation_history_save` |
| Slices (MPI-IO binary) | `<output_basename>-slice_NN-<it:9>.mat` | `[slices]`, `[slice_N]` |

---

## Library sections in FLUME

`[grid]`, `[amr]`, `[amr_marker_N]`, `[field]`, `[runge_kutta]`, `[weno]`, `[fdv]`, `[solids]`, `[solid_N]`, `[slices]`
and `[slice_N]` are read by the library and documented in the [library input reference](/reference/library-input).
FLUME adds these rules:

| Section | FLUME rule |
|---------|------------|
| `[grid]` | `ngc >= [weno] S` (fatal); `>= 2` with `scheme_space = weno-riemann`; with gravity on a gas with an energy equation, one more (`S + 1`, or `max(S, 2) + 1` with `weno-riemann`), the face beyond each block end ([#73](https://github.com/szaghi/adam/issues/73)); MHD: `fdv_order/2 <= ngc` (CPU) or `<= min(ngc, 3)` (FNL), checked at the first div(B) diagnostic. With `[initial_conditions] amr_iterations > 0`, `ni`/`nj`/`nk` must be even along every non-null axis (fatal). With markers and `max_level > iu_ref_levels`, every refined non-null axis needs an even cell count ([#39](https://github.com/szaghi/adam/issues/39)) and at least `2 ngc` cells ([#66](https://github.com/szaghi/adam/issues/66)); both fatal at initialisation. For Euler a null direction also freezes the momentum along it. |
| `[amr]` | `frequency > 0` regrids at runtime ([#74](https://github.com/szaghi/adam/issues/74)): negative is fatal, and so is runtime regridding without markers. On a multi-realm forest each realm regrids by its own `[amr]` section and the forest balances the marks of the realms across each seam ([#110](https://github.com/szaghi/adam/issues/110); see [Setting up runtime AMR](#setting-up-runtime-amr)). `regrid_prolongation` defaults to `conservative` in FLUME (the library default is `linear`); `conservative-state` ([#88](https://github.com/szaghi/adam/issues/88)) is for self-gravitating collapse ([numerics](./numerics#runtime-regridding)). Derefinement is always the mean of the children. A static grid (`frequency = 0`) is balanced by block count after the initial AMR, `balance_static` ([#123](https://github.com/szaghi/adam/issues/123), [library input](/reference/library-input)). |
| `[amr_marker_N]` | Modes 1 (`solid`, `primitive-box`), 2, 4, 5 and 6 (a registered criterion); mode 3, `stl` and an unknown mode are fatal at initialisation. A Löhner marker (mode 4) acts at the runtime regrids only, with `buffer <= ngc - 1` (fatal otherwise). A Jeans marker (mode 5) reads `jeans_cells` (default 8) and `jeans_cells_derefine` (default twice `jeans_cells`) and needs `[gravity] self_gravity = .true.` (fatal otherwise) ([numerics](./numerics#runtime-regridding)). At a runtime regrid the markers vote per block: refine when any asks, coarsen only when all agree, keep otherwise, never below `iu_ref_levels`; the box marker votes to coarsen outside its box. The gradient is centred, across the faces between blocks too ([#76](https://github.com/szaghi/adam/issues/76)), one-sided at the domain boundary. |
| `[runge_kutta]` | `runge-kutta-yoshida` is fatal (no FLUME integrator); the staged path (2:1 faces, forests) and the positivity limiter need an SSP scheme. |
| `[weno]` | `weno-c-*` is fatal. `weights = si` makes the scheme independent of the units of the data ([numerics](./numerics#weno-reconstruction)); `si` with `weno-riemann` and `reconstruction_variables = primitive` is fatal. With `scheme_space = weno-riemann` the same `scheme` sets the order of the WENO interpolation tables. ROR is not implemented: its keys are read and unused. |
| `[fdv]` | `fdv_scheme` is read and unused (the div(B) diagnostic is always centred FD); `fdv_order` sets the order of that diagnostic, MHD only. |
| `[solids]` | Euler only: `number > 0` with `mhd-ideal` is fatal. |
| `[slices]` | Only the conservative variables are sliced. |

Variable indices for the gradient marker:

| Index | Conservative (`field = 1`) | Auxiliary (`field = 2`) |
|-------|----------------------------|-------------------------|
| 1 | `r` | `rho` |
| 2 | `ru` | `u` |
| 3 | `rv` | `v` |
| 4 | `rw` | `w` |
| 5 | `rE` | `p` |
| 6 | `bx` (MHD) | `T` |
| 7 | `by` (MHD) | `H` |
| 8 | `bz` (MHD) | `a` |
| 9 | `psi` (MHD + GLM or EGLM) | `Bx` (MHD) |
| 10 | — | `By` (MHD) |
| 11 | — | `Bz` (MHD) |

### Setting up runtime AMR

Runtime AMR ([#74](https://github.com/szaghi/adam/issues/74)) regrids during the run, on both backends. A minimal setup
tracks a feature from a base level to a finest level (the spacings below are for a unit box with 16 cells per block):

```ini
[amr]
max_level      = 3          ; finest level
iu_ref_levels  = 1          ; base level: no block is ever coarsened below it
frequency      = 5          ; regrid after every 5th step
iters          = 2          ; up to 2 sweeps per regrid (a sweep moves a block by one level)
markers_number = 2
; regrid_prolongation = conservative   (the FLUME default; the refined children average to the parent)

[amr_marker_1]               ; refines the initial grid too (Loehner does not act at initialisation)
mode         = 2             ; gradient
delta_type   = max
delta_fine   = 0.0079        ; ~1 % above the finest spacing: the blocks where the gradient fires reach it
delta_coarse = 0.0316        ; ~1 % above the base spacing
field        = 1
var          = 1             ; density
tol          = 0.05

[amr_marker_2]
mode         = 4             ; Loehner: scale-free, ~0 where smooth, ~1 at a jump
delta_type   = max           ; read but unused by mode 4
delta_fine   = 0.0
delta_coarse = 0.0
field        = 1
var          = 1
refine_tol   = 0.5           ; refine above
derefine_tol = 0.2           ; coarsen below (keep between: the hysteresis band)
buffer       = 2             ; ghost layers read: a feature at a neighbour's edge marks the block
; floor      = 1e-4          ; needed when the variable vanishes in quiet regions (e.g. B outside a loop)
```

How the pieces combine, and what to expect:

- **Markers vote.** A block is refined when any marker asks for it and coarsened only when all agree; the box marker
  votes to coarsen outside its box. Use a gradient (or solid, or box) marker to shape the initial grid and Loehner to
  follow the solution: without an initial marker the first `frequency` steps run on the base grid, and their error
  persists.
- **Thresholds.** Loehner's `refine_tol` around 0.5 follows shocks and contacts; smooth structures (a vortex, a
  rarefaction) score low and may stay coarse, so mark them with a gradient marker. On a variable that is zero in quiet
  regions set `floor` to a fraction of the feature's jump, or round-off noise marks the whole grid.
- **Cost.** The time step is global (no subcycling), so a run regridding to level L advances at the step of level L
  everywhere; the savings come from the cells, and on the verification cases the wall time dropped by 3–19 %
  ([verification](./verification/runtime-amr/accuracy)). On the FNL backend each regrid is a host round trip
  whose copies scale with the blocks in use ([#75](https://github.com/szaghi/adam/issues/75)).
- **Limits.** On a multi-realm forest ([#110](https://github.com/szaghi/adam/issues/110)) the blocks along each side
  of a `mirror` or `refined` seam keep one level: a refine vote on one refines the whole side. A `coupling = auto` seam
  may hold different levels along it, face by face (equal cells or 2:1), with per-axis checks
  ([#71](https://github.com/szaghi/adam/issues/71)); mixed levels along a seam are verified on octrees only. With EGLM the 2:1 seams raise the divergence error B_z well above its uniform-grid level
  ([#78](https://github.com/szaghi/adam/issues/78)).
- **Restart.** A run restarted from any step regrids on the same steps and continues bitwise (the regrid happens before
  the step's output and restart files). On a forest every realm restarts from its own files and the inter-realm seams
  are rebuilt from the restored trees, mixed levels along a seam included, bitwise under either `coupling_cadence`
  ([#113](https://github.com/szaghi/adam/issues/113)).
  Restart on the rank count the run was saved with: fewer ranks stop the run.
- **Log.** Each regrid prints `flume: regrid at step N: blocks a -> b, R refined, C coarsened, S sweeps` (plus its wall
  time on FNL); a state with non-positive density or pressure after a regrid stops the run with its location.
- **Your own criterion** ([#110](https://github.com/szaghi/adam/issues/110)). Extend `amr_criterion_object`
  (`adam_amr_criterion_object`) with a `vote(realm, b)` returning `TO_BE_REFINED`, `TO_NOT_TOUCH` or `TO_BE_DEREFINED`
  for block `b`, register an instance on every realm in a copy of the driver (`call realm(i)%register_amr_criterion(c)`,
  before `simulate`; `c` a TARGET variable of the main program, since the realm keeps a pointer to it; an unregistered
  name stops the run at initialisation), and select it with `[amr_marker_N] mode = 6` and `name = <its name>`; its votes combine with the
  other markers' as theirs do, at runtime regrids. The example `src/tests/flume/verification/regrid/test_flume_user_criterion.F90`
  is the CPU driver plus a criterion that reproduces the built-in box marker (bitwise the same run).

---

## [reference] (optional: dimensional input)

FLUME's equations carry no dimensionless number (`B` is `B_SI/sqrt(mu0)`), so a dimensional input runs once every
value is divided by the reference of its dimension ([#49](https://github.com/szaghi/adam/issues/49)). With this section present, FLUME does that once, on
the loaded file, before any other section is parsed; without it nothing changes. References: density `rho0`, length
`L0`, velocity `u0`; derived: time `L0/u0`, pressure and energy density `rho0 u0^2`, field `u0 sqrt(rho0)`, the GLM
`psi` `u0^2 sqrt(rho0)` (the EGLM one as the field). `cp`, `cv` are replaced by `gamma = cp/cv` (the gas constant
becomes 1, the temperature unit `u0^2/R`); the physics derives `cv = R/(gamma-1)`, `cp = gamma cv` on both paths, so the
conversion is bit for bit. Restart files and logs stay in code units; the other outputs follow
`output_units`. Every restart save also writes `<restart_basename>.reference` (the three references, 1 without this section); a restart under different
references is refused, and restart files without that record (written before it existed, so in code units) restart only
when the references are all 1.

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `density` | real | no (1) | `> 0`, otherwise fatal | `rho0`. |
| `length` | real | no (1) | `> 0`, otherwise fatal | `L0`. |
| `velocity` | real or string | no (1) | `> 0`; `acoustic` (`sqrt(gamma pressure / density)`, needs `pressure`); `alfvenic` (`field / sqrt(density)`, needs `field` and `mhd-ideal`). Other → fatal. | `u0`. |
| `pressure` | real | cond. | `> 0` | Reference pressure of the `acoustic` preset. |
| `field` | real | cond. | `> 0` | Reference field of the `alfvenic` preset. |
| `output_units` | string | no (`code`) | `code`, `dimensional`; other → fatal | Units of the fields, grid, time, slices and histories. |

With `output_units = dimensional` every written value is multiplied, at write time only, by the reference of its
dimension: the fields (conservative, residuals, auxiliary, MHD derived: `beta` is dimensionless, `divb` a field per
length), the block origins and spacings, the time, the slice points and values, the conservation integrals (times
`L0^3`, the cell volume includes the null directions), the div(B) and residual histories. The temperature is written
as `p / (rho R)` with the gas constant of the input (`cp - cv`, or 1 with `gamma`). `<output_basename>.units` (written
whenever the section is present) records the output units, the references, the derived ones and the factor of every
written variable. The state, the restart files and the solver are unaffected.

Any other key is fatal. With the section active **every** section and key of the file must be known to the layer:
an unknown one is fatal, so no dimensional value can pass unconverted. Keys whose dimension depends on the context are
resolved from it: `[amr_marker_N].tol` of a gradient marker takes the dimension of the marked variable per length,
`floor` of a Löhner marker that of the marked variable (a marker on the temperature is refused; the Löhner
`refine_tol`, `derefine_tol` and `epsilon` are dimensionless), `[initial_conditions].wave_amplitude` that of the
eigenvector it multiplies.
`type = orszag-tang` (a state hard-coded in code units) and multi-realm runs are refused. The log lists every
converted key with its old and new value.

The dissipative coefficients and the wall keys convert as follows; the numbers are
dimensionless and pass unchanged. $R$ is the gas constant of the input ($c_p - c_v$, or 1 with `gamma`), since the layer
replaces `cp`, `cv` by `gamma` and the code temperature becomes $p/\rho$:

| Key | Divided by |
|---|---|
| `viscosity` | $\rho_0 u_0 L_0$ |
| `conductivity` | $\rho_0 u_0 L_0 R$ |
| `resistivity` | $u_0 L_0$ |
| `[bc_*]` `wall_model_height` | $L_0$ |
| `[bc_*]` `wall_model_filter` | $L_0 / u_0$ |
| `[bc_*]` `wall_model_pressure_gradient_filter` | $L_0 / u_0$ |
| `[bc_*]` `wall_model_edge_velocity` | $u_0$ |
| `[bc_*]` `wall_model_edge_temperature` | $u_0^2 / R$ |
| `[forcing]` `bulk_momentum` | $\rho_0 u_0$ |
| `[initial_conditions]` `ch_rho` | $\rho_0$ |
| `[initial_conditions]` `ch_ub` | $u_0$ |
| `[initial_conditions]` `ch_h`, `ch_lx`, `ch_lz` | $L_0$ |
| `[initial_conditions]` `ch_mach`, `ch_pert`; every `[turbulence]` key; `[forcing]` `type`, `axis` | — (dimensionless) |
| `wall_u`, `wall_v`, `wall_w` | $u_0$ |
| `wall_temperature`, `reference_temperature`, `sutherland_temperature` | $u_0^2 / R$ |
| `temperature_unit_kelvin` | its inverse: kelvin per input temperature unit, times $u_0^2 / R$ |
| `reynolds`, `prandtl`, `magnetic_reynolds`, `viscosity_exponent` | — (dimensionless) |
| `lundquist` | — (refused unless `velocity = alfvenic`) |
| `[gravity]` `g_x`, `g_y`, `g_z` | $u_0^2 / L_0$ |
| `[gravity]` `G` | $u_0^2 / (\rho_0 L_0^2)$ |
| `[gravity]` `froude`, `G_star`; `[initial_conditions]` `polytropic_gamma` | — (dimensionless) |
| `[cooling]` `temperatures` (each value) | $u_0^2 / R$ |
| `[cooling]` `lambdas` (each value) | $L_0^5 u_0^3 \rho_0$ |
| `[cooling]` `particle_mass` | $\rho_0 L_0^3$ |
| `[cooling]` `dt_factor` | — (dimensionless) |
| `[sinks]` `x`, `y`, `z` (each value) | $L_0$ |
| `[sinks]` `u`, `v`, `w` (each value) | $u_0$ |
| `[sinks]` `mass` (each value) | $\rho_0 L_0^3$ |

The conversion divides by the power-of-two reference first and applies $R$ last, a single rounding, so the
verification DC can check every converted value exactly.

---

## [physics]

`physical_model` is required, with either `gamma` or both `cp` and `cv` (neither with
[`[eos]`](#eos-optional-equation-of-state) `isothermal` or `barotropic`).

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `physical_model` | string | yes | `euler` (nv = 5), `euler-de` (nv = 6) or `mhd-ideal` (nv = 8, or 9 with GLM or EGLM). Other → fatal. | Equation set. `euler-de` is Euler with dual energy, the advected entropy `rs` ([#91](https://github.com/szaghi/adam/issues/91), [models](./models#dual-energy); CPU and FNL). `mhd-ideal` also loads `[mhd]`. |
| `passive_scalars` | int | no (`euler`, `euler-de`) | `0..4`, default 0; above 4 or with `mhd-ideal` → fatal | Passive scalars ([#92](https://github.com/szaghi/adam/issues/92), [models](./models#passive-scalars)): `N > 0` selects the scalar variant of the model, 4 more conserved variables `rs1`..`rs4` (`nv` 9, or 10 with `euler-de`), the first `N` set by the `s1`..`sN` keys of the initial-condition regions and of the inflows, the others zero. CPU and FNL. |
| `dual_energy_eta` | real | no (`euler-de`) | `(0, 1)`, default 0.05; otherwise fatal | Dual-energy switch: a cell whose internal energy is below `dual_energy_eta` times its total energy takes the entropy as the floor of its internal energy. |
| `gamma` | real | cond. | `> 1`; with `cp` or `cv` → fatal | Specific heats ratio; the gas constant is 1 (code units). |
| `cp` | real | cond. | `cp > cv > 0`, otherwise fatal; required without `gamma` | Specific heat at constant pressure. |
| `cv` | real | cond. | as above | Specific heat at constant volume. `gamma = cp/cv`, `R = cp - cv`; the physics then uses `cv = R/(gamma-1)`, `cp = gamma cv` (an input `cp` may move by an ulp). |

<a id="dissipative-terms-issue-65-m4"></a>

### Dissipative terms

Each term is given either as its coefficient or as the dimensionless number it stands for. In code units the
references are 1, so a number is the reciprocal coefficient. Giving two keys of a term is fatal; giving none leaves the
term off, so an input without these keys stays ideal and bitwise unchanged. The viscosity and the conductivity drive
the Navier–Stokes fluxes of every model, the resistivity the Ohmic ones of `mhd-ideal`
([models](./models#dissipative-terms), [numerics](./numerics#dissipative-fluxes-navier–stokes)).

| Key | Type | Accepted / invalid | Meaning |
|-----|------|--------------------|---------|
| `viscosity` | real | `>= 0`; with `reynolds` → fatal | Dynamic viscosity $\mu$ (code units, or dimensional with [`[reference]`](#reference-optional-dimensional-input)). |
| `reynolds` | real | `> 0` | $Re$: $\mu = 1/Re$. |
| `conductivity` | real | `>= 0`; with `prandtl` → fatal | Thermal conductivity $k$; the heat flux is $-k\nabla T$ with $T = p/(\rho R)$. |
| `prandtl` | real | `> 0`; needs `viscosity` or `reynolds` | $Pr$: $k = \mu c_p / Pr$ (follows $\mu(T)$ under the power and Sutherland laws). |
| `resistivity` | real | `>= 0`; MHD only; with a number → fatal | Magnetic diffusivity $\eta$ (Ohmic resistivity). |
| `magnetic_reynolds` | real | `> 0`; MHD only | $Rm$: $\eta = 1/Rm$. |
| `lundquist` | real | `> 0`; MHD only; with `[reference]` it needs `velocity = alfvenic` | $S$: $\eta = 1/S$, the magnetic Reynolds number at the Alfvén speed. |
| `viscosity_law` | string | `constant` (default), `power-law`, `sutherland`; other → fatal; `power-law` and `sutherland` need a viscosity | Temperature law of the viscosity. |
| `viscosity_exponent` | real | `>= 0`; required by `power-law`, fatal without it | $\omega$ in $\mu(T) = \mu\,(T/T_\mathrm{ref})^\omega$ (non-negative: the diffusive time step bound evaluates the laws at the largest temperature a face can see). |
| `reference_temperature` | real | `> 0`; required by `power-law` and `sutherland`, fatal without it | $T_\mathrm{ref}$ (code temperature $p/(\rho R)$), where $\mu$ is the given viscosity. |
| `temperature_unit_kelvin` | real | `> 0`; optional | The kelvin value of one temperature unit of the input (`p/(\rho R)`): the absolute-temperature anchor of a non-dimensional input, read by what needs kelvin (Sutherland's default $S$; the thermally perfect gas to come). With a dimensional `[reference]` input in kelvin it is 1, and the layer converts it to kelvin per code unit. |
| `sutherland_temperature` | real | `> 0`; required by `sutherland` unless `temperature_unit_kelvin` is given (then $S$ = 110.4 K, air), fatal with another law; `viscosity_exponent` is fatal with `sutherland` | $S$ in Sutherland's law $\mu(T) = \mu\,(T/T_\mathrm{ref})^{3/2}\,(T_\mathrm{ref}+S)/(T+S)$, in the units of `reference_temperature` (110.4 K for air: with a dimensional `[reference]` input both are given in kelvin; in code units $S = 110.4\,T_\mathrm{ref}/T_\mathrm{ref,K}$). |

Fatal with any non-zero coefficient: `[numerics] positivity_limiter = cell`
([why](./numerics#dissipative-fluxes-navier–stokes)) and immersed solids. `dissipative_order = 4` with `[grid] ngc < 3`
is fatal too.

---

## [eos] (optional: equation of state)

Without the section the gas is the ideal gas of `[physics]` (`gamma`, or `cp` and `cv`). `type = isothermal` and
`type = barotropic` ([#100](https://github.com/szaghi/adam/issues/100)) select a pressure $P(\rho)$ and drop the energy
equation; the energy slot `rE` stays 0 and is still written ([models](./models#isothermal-and-barotropic-gases)).

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `type` | string | yes | `ideal`, `isothermal`, `barotropic`; other → fatal | The equation of state. `ideal` is the default gas of `[physics]`. |
| `sound_speed` | real | with `isothermal`, `barotropic` | `> 0`; missing or non-positive → fatal | $c_0$: the sound speed of the isothermal gas; for a table, the isothermal speed $\sqrt{P/\rho}$ at `reference_density`, so $K_1 = c_0^2\rho_\mathrm{ref}^{1-\Gamma_1}$ (the sound speed there when $\Gamma_1 = 1$). |
| `reference_density` | real | no (1) | `> 0`, otherwise fatal | $\rho_\mathrm{ref}$ of `barotropic`. |
| `densities` | list of reals | no | all `> 0`, strictly increasing; otherwise fatal | The breaks $\rho_1 < \dots < \rho_{n-1}$ of `barotropic`; a break belongs to the segment above it. |
| `exponents` | list of reals | no | one per segment (one more than `densities`), all `>= 1`, at most 8 segments; otherwise fatal | $\Gamma_1, \dots, \Gamma_n$ of `barotropic`. $K_{k+1} = K_k \rho_k^{\Gamma_k - \Gamma_{k+1}}$ keeps $P$ continuous. |

Without `densities` and `exponents`, `barotropic` takes the table of Mathew & Federrath (2021, MNRAS 507, 2448),
$\Gamma = 1, 1.1, 1.4, 1.1, 5/3$ with breaks at $2.5\cdot10^{-16}$, $3.84\cdot10^{-13}$, $3.84\cdot10^{-8}$ and
$3.84\cdot10^{-3}$ g cm$^{-3}$: a cgs input (their $c_0$ is 0.2 km/s). A break where $\Gamma$ decreases upward is
accepted and reported at start-up.

With `type = isothermal` or `barotropic`:

- `physical_model` must be `euler` or `mhd-ideal`, without `passive_scalars`; `gamma`, `cp` and `cv` are fatal (they
  would be ignored);
- the regions of `[initial_conditions]` and the `inflow` faces need no `p` (it is ignored);
- the accepted initial conditions are `uniform`, `riemann-problem`, `rotated-riemann`, `sine-wave` (`wave_mode =
  acoustic` or `jeans`), `linear`, `polynomial`, `shu-osher` and `boss-bodenheimer`, with `isothermal` also
  `singular-isothermal-sphere` and `bonnor-ebert`; the others are fatal;
- `[gravity]` acts on the momentum only, so the `ngc >= S + 1` requirement of gravity does not apply;
- `[numerics] positivity_limiter = cell` is the density-only limiter ([numerics](./numerics#isothermal-and-barotropic-equations-of-state));
  the sinks need it;
- fatal: `[cooling]`, the dissipative terms, `[reference]`, immersed solids, and the `wall-noslip`, `wall-isothermal`
  and `hydrostatic` faces;
- the schemes are `weno` and `weno-riemann`, all reconstructions, with `llf`, `hll`, `hllc` (Euler) or `llf`, `hll`,
  `hlld` (MHD), on CPU and FNL;
- with `physical_model = mhd-ideal` ([models](./models#isothermal-and-barotropic-mhd)): `[mhd] divergence_control`
  must be `none` or `glm` (`eglm` is fatal); `rho_floor` and `p_floor` must be 0; the accepted initial conditions are
  those above plus `glm-pulse`, `divb-peak`, `mhd-linear-wave` (`wave = entropy` is fatal), `mhd-cpaw`, `mhd-rotor`
  and `field-loop`; `positivity_limiter = cell` is fatal with `glm`.

---

## [mhd] (only with `physical_model = mhd-ideal`)

The section is not read for `euler`.

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `divergence_control` | string | yes | `glm` (nv = 9, psi added), `eglm` (nv = 9, psi also in the energy; the `glm_*` keys apply) or `none` (nv = 8). Other → fatal. | div(B) control ([models](./models#divergence-control-mixed-glm)). Under `eglm` the AMR reflux leaves B alone ([numerics](./numerics#reflux)). |
| `glm_ch` | real | `glm`, `eglm` | `> 0`, otherwise fatal | Constant cleaning speed `c_h`. It also bounds the time step (`c_h * sum 1/dx`). Keep it just above the fastest signal speed: a larger one speeds up the growth of small asymmetries in sheared, magnetised flows ([#83](https://github.com/szaghi/adam/issues/83)). |
| `glm_alpha` | real | `glm`, `eglm` | `>= 0`, otherwise fatal | Damping parameter: `c_h^2/c_p^2 = glm_alpha * c_h / L`. |
| `glm_damping_length` | real or `min-cell` | `glm`, `eglm` | A positive real, or the literal `min-cell` (the minimum cell spacing over the active directions, MPI-reduced; with runtime regridding, `[amr] frequency > 0`, the spacing of the finest level `max_level` allows, fixed for the run). A non-number or a value `<= 0` is fatal. | Damping length `L`. |
| `glm_ch_check` | string | `glm`, `eglm` | `warning` or `error`. Other → fatal. Without GLM it is forced to `warning` and not read. | What happens when `max(\|u\|+c_f) > glm_ch`: warn (logged each time a new maximum appears) or stop. |
| `divb_tol` | real | yes | `>= 0`, otherwise fatal. `0` disables the monitor; `> 0` needs `[diagnostics] divb_flux = .true.` (hence GLM or EGLM), fatal otherwise. | Monitor on the normalised D_flux ([numerics](./numerics#output-and-diagnostics)): `max h\|D_flux\| / max\|B\| > divb_tol` warns or stops. Dimensionless, not converted by `[reference]`. |
| `divb_error` | logical | yes | logical | `.true.`: exceeding `divb_tol` is fatal. |
| `rho_floor` | real | yes | `>= 0`, otherwise fatal. `0` disables it. | Density positivity floor. |
| `p_floor` | real | yes | `>= 0`, otherwise fatal. `0` disables it. | Pressure positivity floor. With **both** floors at 0, any non-positive density or pressure stops the run. |

---

## [gravity] (optional: external gravity and self-gravity)

A constant acceleration $\mathbf g$ ([#73](https://github.com/szaghi/adam/issues/73)) and the self-gravity of the gas
([#81](https://github.com/szaghi/adam/issues/81)); the sources are on the [models](./models#external-gravity) page.
Without the section there is no gravity and nothing changes.

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `g_x`, `g_y`, `g_z` | real | one of the two forms | any; each optional (default 0) | The acceleration. Dimensional (length / time²) with [`[reference]`](#reference-optional-dimensional-input). |
| `froude` | real | the other form, with `axis` | `> 0`, otherwise fatal | The number form: $g = 1/\mathrm{Fr}^2$ along `axis` ($g L_0/u_0^2$ in the reference units; not converted by `[reference]`). |
| `axis` | string | with `froude` | `+x`, `-x`, `+y`, `-y`, `+z`, `-z`; other → fatal | The signed direction of $\mathbf g$. |
| `equilibrium` | string | no (default `none`) | `none` or `atmosphere`; other → fatal. `atmosphere` needs an `isothermal-atmosphere`, `polytropic-atmosphere` or `rayleigh-taylor` initial condition, and is fatal with `[numerics] reflux_policy = troubled` and on multi-realm runs | The well-balanced mode: the initial atmosphere at rest is kept to round-off ([models](./models#external-gravity)). |
| `potential` | string | no (default `uniform`) | `uniform` or `linear`; other → fatal. `linear` is fatal with `equilibrium = atmosphere` and with `self_gravity` | `uniform`: the constant acceleration in every kernel, with the fourth-order energy source. `linear`: the same field as the cell potential Φ = −g·x, applied through the potential path of self-gravity ([models](./models#external-gravity)). |
| `self_gravity` | logical | no (default `.false.`) | with it: `equilibrium = atmosphere`, `potential = linear`, multi-realm runs, immersed solids and a build without an FFT library (the `*-fft` modes) are fatal, and so are the boundaries `poisson_bc` cannot hold (below) | Self-gravity, solved every Runge–Kutta stage ([models](./models#self-gravity)). |
| `G` | real | one of `G`, `G_star` with `self_gravity` | `> 0`, otherwise fatal; fatal without `self_gravity`, and together with `G_star` | The gravitational constant. Dimensional with [`[reference]`](#reference-optional-dimensional-input). |
| `G_star` | real | the other form | as `G` | The number $G\rho_0L_0^2/u_0^2$; not converted by `[reference]`. |
| `poisson_solver` | string | no (default `auto`) | `auto`, `fft` or `multigrid`; other → fatal. `fft` is fatal on an AMR grid | `auto`: FFT on a uniform grid, multigrid on an AMR grid (leaves at more than one level). `fft`: the library FFT solver only. `multigrid`: the multigrid solver on any grid ([models](./models#self-gravity)). |
| `poisson_bc` | string | no (default `periodic`) | `periodic` or `isolated`; other → fatal. `periodic` is fatal with a non-periodic active direction; `isolated` is fatal in 2-D, with a periodic direction, with non-cubic cells, and on AMR grids with a leaf above the coarsest level on the domain boundary (at start or at the regrid that makes one, [#87](https://github.com/szaghi/adam/issues/87)) | `periodic`: every active direction periodic, the mean density removed. `isolated`: the free-space potential of the mass on the grid; the fluid boundaries are the `[bc_*]` ones, e.g. `wall-inviscid` for a closed box. |

Giving `g_*` together with `froude` is fatal, as is a section with neither (unless `self_gravity`, where $\mathbf g$
defaults to 0).

Constraints:

- external gravity on a gas with an energy equation needs `[grid] ngc >= S + 1` ($S$ the WENO half-stencil, 3 for
  WENO-5; with `weno-riemann` $\max(S, 2) + 1$), fatal otherwise;
- gravity along a null direction and gravity with immersed solids are fatal;
- `self_gravity` needs a `*-fft` build: `flume-cpu-gnu-fft` (FFTW 3, executable `adam_flume_cpu_fft`) or
  `flume-fnl-nvf-fft` (cuFFT, `adam_flume_fnl_fft`); the default modes refuse it at the first solve.

With `self_gravity` the conservation history gains the columns `W` $= \tfrac12\sum\rho\Phi\,dV$ and `E_plus_W`, the
energy the scheme conserves.

## [cooling] (optional: optically thin radiative cooling)

The section switches on optically thin radiative cooling ([#93](https://github.com/szaghi/adam/issues/93);
[models](./models#radiative-cooling), [numerics](./numerics#radiative-cooling)). Without the section nothing changes.
The hydrodynamic models only (`euler`, `euler-de`, with or without passive scalars): with `mhd-ideal` it is fatal.

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `temperatures` | list of reals | yes | 2 to 64 values (`COOLING_NODES_MAX`: the device kernels take fixed-size tables), all `> 0`, strictly increasing; otherwise fatal | The nodes $T_1 < \dots < T_N$ of the cooling function. $T_1$ is the floor: below it nothing cools, and a cooling cell stops exactly at it. |
| `lambdas` | list of reals | yes | as many values as `temperatures`, all `> 0`; otherwise fatal | $\Lambda_k$ at the nodes. $\Lambda$ is the piecewise power law through them, $\Lambda_k (T/T_k)^{\alpha_k}$ on $[T_k, T_{k+1}]$; above $T_N$ the last segment continues. |
| `particle_mass` | real | yes | `> 0`, otherwise fatal | The mass $m$ per particle of $n = \rho/m$ (for $n_\mathrm{H}$: $m_\mathrm{H}/X$). |
| `dt_factor` | real | no (0) | `>= 0`, otherwise fatal | `0`: off. Otherwise the step is also limited to `dt_factor` times the shortest cooling time $t_\mathrm{cool} = \rho e/(n^2\Lambda)$, for runs that need the coupling of cooling and dynamics resolved (the cooling itself is exact for any step). |

The values are in the units of the rest of the input: $T$ uses the gas constant of `[physics]` ($c_p - c_v$, or 1 with
`gamma`), so a table in kelvin needs `cp`, `cv` in J/(kg K) and $\Lambda$, $m$ in the same unit system (with
$\Lambda$ in erg cm³ s⁻¹ and SI elsewhere, multiply $\Lambda$ by $10^{-13}$). The `[reference]` layer converts the
three keys (above). The conservation history gains two columns: `radiated`, the energy radiated since the start of the
run (of this run segment: it restarts from 0 on a restart), and `E_plus_radiated`, which the scheme conserves.

---

## [sinks] (optional: sink particles)

The section switches on sink particles ([#94](https://github.com/szaghi/adam/issues/94);
[models](./models#sink-particles), [numerics](./numerics#sink-particles)). It needs `[gravity] self_gravity = .true.`
and a hydrodynamic model (`euler`, `euler-de`, with or without passive scalars); otherwise it is fatal. The keys list
the initial sinks, one value each; their number is that of `mass` (none: the section with no initial sink).

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `mass` | list of reals | no (no sink) | all `> 0`; otherwise fatal | The sinks' masses. |
| `x`, `y`, `z` | lists of reals | with `mass` | as many values as `mass`, inside the domain; otherwise fatal | The positions. |
| `u`, `v`, `w` | lists of reals | no (0) | as many values as `mass`; otherwise fatal | The velocities. |
| `creation` | logical | no (`.false.`) | | Sinks form where the gas passes the checks of Federrath et al. (2010) ([models](./models#sink-particles)). |
| `accretion` | logical | no (`.true.`) | | Sinks accrete the gas above the threshold within their accretion radius. |
| `merging` | logical | no (`.true.`) | | Two bound, converging sinks within an accretion radius merge. |
| `accretion_radius_cells` | real | no (2.5) | in `[1, ngc]`; otherwise fatal | The accretion radius $r_\mathrm{acc}$ in cells of the sink's level (with runtime regridding, of the finest level allowed). It sets the threshold density $\rho_\mathrm{thr} = \pi c_s^2/(4 G r_\mathrm{acc}^2)$. |

The conservation history gains `M_sinks`, `Px_sinks`, `Py_sinks`, `Pz_sinks` and `K_sinks` (`W` then includes the
sinks' mass), and `<output_basename>-sinks_history.dat` gets one row per sink per save (`it time id m x y z u v w`).
Restart files carry `<restart_basename>.sinks`. A sink that leaves a closed domain is fatal.

---

## [turbulence] (optional: algebraic subgrid-scale model)

A library section ([#119](https://github.com/szaghi/adam/issues/119); [turbulence](./turbulence#subgrid-scale-models)). Without it, or with `sgs_model = none`, nothing changes.

| Key | Type | Required | Values | Meaning |
|---|---|---|---|---|
| `sgs_model` | string | no (`none`) | `none`, `vreman`, `wale`, `sigma`. Other → fatal. | Eddy viscosity of the model. |
| `sgs_constant` | real | no | `> 0`, otherwise fatal; defaults 0.07 (Vreman: $c = 2.5 C_s^2$), 0.5 (WALE), 1.35 (Sigma) | Model constant. Taken from other discretizations, not calibrated on WENO. |
| `prandtl_t` | real | no (0.9) | `> 0`, otherwise fatal | Turbulent Prandtl number. |

A model needs `physical_model = euler` or `euler-de` (fatal otherwise) and `ngc >= 3` (fatal otherwise). It runs without a viscosity (`viscosity = 0`). It counts as a dissipative term, so the positivity limiter and immersed solids are refused with it. The turbulence history gains a last column, `eps_sgs`.

## [forcing] (optional: mass-flow forcing of a periodic channel)

| Key | Type | Required | Accepted | Meaning |
|---|---|---|---|---|
| `type` | string | no (`none`) | `none`, `mass-flow`; other → fatal; Euler only | Constant mass flow along `axis` ([turbulence](./turbulence#mass-flow-forcing)). |
| `axis` | string | with `mass-flow` | `x`, `y`, `z` | Direction of the flow. |
| `bulk_momentum` | real | with `mass-flow` | | Target mean of the momentum along the axis. |

The force per unit volume goes to `<output_basename>-forcing_history.dat` on the conservation-history cadence; a restart appends to it.

## [numerics]

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `scheme_space` | string | yes | `weno`, `weno-riemann` or `hybrid-kep-weno` (Euler, `ngc >= 3`, no limiter, no immersed solid, all fatal otherwise). Other → fatal. | Spatial operator ([numerics](./numerics), [turbulence](./turbulence#hybrid-central-weno)). MHD at very low β (≲ 1e-3): use `weno-riemann` with `hlld` where symmetry or backend agreement matter ([#77](https://github.com/szaghi/adam/issues/77), [verification](./verification/runtime-amr/mhd-blast-symmetry)). |
| `hybrid_threshold` | real | no (0.2), only with `hybrid-kep-weno` | `> 0`, otherwise fatal | Ducros sensor above which a cell keeps WENO ([turbulence](./turbulence#hybrid-central-weno)). |
| `hybrid_floor` | real | no ($10^{-8}$), only with `hybrid-kep-weno` | `> 0`, otherwise fatal | Floor of the sensor denominator over $(a/\Delta x_\min)^2$. Not calibrated. |
| `hybrid_jump` | real | no (0.1), only with `hybrid-kep-weno` | `> 0`, otherwise fatal | Largest relative density jump of a central cell. Not calibrated. |
| `hybrid_skin` | string | no (`weno`), only with `hybrid-kep-weno` | `weno`, `central-seams`, `central-bcs`, `central`, otherwise fatal | Which block faces keep WENO on their skin: all that are not same-level; all but resolution jumps; all but boundaries and forest seams; none. The relaxed values are for measurement ([verification](./verification/turbulence/hybrid-scheme#on-a-2-1-box-and-what-the-skin-costs)). With `central-bcs` the sensor ghost at a forest seam is computed locally: not tested. |
| `reconstruction_variables` | string | yes | With `weno`: `characteristic` or `conservative`. With `weno-riemann`: `characteristic` or `primitive`. Other → fatal. | Variables reconstructed (or interpolated) at the faces. |
| `riemann_solver` | string | `weno-riemann` | `llf`, `hll`, `hllc`, `hlld`. Other → fatal. By model: `euler` allows `llf`, `hll`, `hllc`; `mhd-ideal` allows `llf`, `hll`, `hlld` (fatal otherwise). Not read with `weno`. | Face Riemann solver. |
| `flux_correction` | string | `weno-riemann` | `6th`, `4th`, `none` (2nd order). Other → fatal. | Face-flux correction `F^ = c1 F + c2 (f_i+f_i+1) + c3 (f_i-1+f_i+2)`. |
| `flux_correction_sensor` | string | `weno-riemann` | `weno` or `none`. Other → fatal. | `weno`: the correction is switched off at faces where `min_k w_k/d_k < 0.2` (fixed threshold `FLUX_CORRECTION_SENSOR_TAU`). `none`: always on. |
| `reflux` | logical | yes | logical | Berger–Colella reflux at AMR coarse-fine faces (and inter-realm seams). `.false.` is a diagnostic only: the run is then not conservative across 2:1 faces. |
| `reflux_policy` | string | no | `always`, `troubled` (default `always`; any other value is fatal) | With `reflux = .true.`, which seam face cells the reflux corrects: `always` every one, conservative to round-off; `troubled` only those flagged by the positivity limiter or the WENO-5 sensor ([numerics](./numerics#reflux), [#68](https://github.com/szaghi/adam/issues/68)). The log reports the flagged face cells per step. |
| `reflux_cadence` | string | no | `step`, `stage` (default `step`; any other value is fatal; `stage` is fatal on multi-realm runs) | With `reflux = .true.`, when the coarse cells are corrected: once at the end of the step, or at every Runge–Kutta stage ([numerics](./numerics#reflux)). |
| `positivity_limiter` | string | no (default `none`) | `none` or `cell`. Other → fatal. `cell` is fatal with `[mhd] divergence_control = glm`, a non-SSP `[runge_kutta] scheme`, immersed solids, multi-realm runs and any dissipative coefficient. | `cell`: the cell-based positivity limiter ([numerics](./numerics#positivity-limiter)); the limited faces are logged per stage. |
| `dissipative_order` | int | no (default 4) | `2` or `4`. Other → fatal; `4` needs `[grid] ngc >= 3` (fatal otherwise). | Order of the dissipative face fluxes ([numerics](./numerics#dissipative-fluxes-navier–stokes)). An ideal run ignores it. |

---

## [bc_x_min], [bc_x_max], [bc_y_min], [bc_y_max], [bc_z_min], [bc_z_max]

All six sections are required, even along null directions.

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `type` | string | yes (every face) | `extrapolation`, `inflow`, `inflow-digital-filter`, `non-reflecting`, `shock-generator`, `wall-inviscid`, `wall-noslip`, `wall-isothermal`, `wall-modelled`, `periodic`, `hydrostatic`. Other → fatal. The two faces of an axis must be **both** periodic or **neither** (fatal). | See the list below. |
| `r`, `u`, `v`, `w`, `p` | real | `type = inflow` | no range check | Prescribed primitive inflow state. |
| `bx`, `by`, `bz` | real | `type = inflow` **and** `mhd-ideal` | no range check | Prescribed inflow magnetic field. `psi` is set to 0. |
| `s1` .. `sN` | real | no (`type = inflow` with `passive_scalars = N`) | default 0, no range check | Specific values (mass fractions) of the passive scalars of the inflow. |
| `wall_u`, `wall_v`, `wall_w` | real | no (default 0) with `wall-noslip`, `wall-isothermal` | the component normal to the face must be 0 (fatal) | Velocity of the wall (tangential): a moving lid, a Couette wall. |
| `wall_temperature` | real | `type = wall-isothermal` | `> 0`, otherwise fatal | Wall temperature, $T = p/(\rho R)$ in code units. |
| `wall_temperature` with `wall-modelled` | real | yes | `> 0`, otherwise fatal | Wall temperature of a modelled wall ([turbulence](./turbulence#wall-model)). Euler with a viscosity and a `prandtl` only, 3 cells per block along the normal, `[grid] ngc >= 3` (fatal otherwise). |
| `wall_model` with `wall-modelled` | string | no (`equilibrium`) | `equilibrium`, `inverse-transformation`; other → fatal | The equilibrium thin-layer model, or the inverse-velocity-transformation model of Griffin, Fu & Moin (2023), which needs the two edge keys below. |
| `wall_model_height` with `wall-modelled` | real | no (0) | `>= 0` | Matching height, a length; 0 samples the 3rd cell. Otherwise interpolated between the cells around it, within the block (beyond it, its last cell). |
| `wall_model_damping` with `wall-modelled` | string | no (`classical`) | `classical`, `semi-local`; other → fatal; given with `wall_model = inverse-transformation` → fatal | Van Driest damping of the equilibrium model in wall units, or with the local density and viscosity. |
| `wall_model_filter` with `wall-modelled` | real | no (0) | `>= 0` | Time scale $T_f$ of the exponential filter of the sample; 0 for none ([turbulence](./turbulence#input-filter)). |
| `wall_model_eddy_viscosity` with `wall-modelled` | logical | no (`.false.`) | | The model's eddy viscosity in the cells below the matching height ([turbulence](./turbulence#eddy-viscosity-below-the-matching-height)). |
| `wall_model_pressure_gradient` with `wall-modelled` | logical | no (`.false.`) | | The pressure gradient along the sample velocity in the equilibrium model's momentum balance; refused with `wall_model = inverse-transformation` ([turbulence](./turbulence#pressure-gradient)). |
| `wall_model_pressure_gradient_filter` with `wall_model_pressure_gradient` | real | no (0) | `>= 0` | Time scale $T_p$ of the exponential filter of the pressure gradient; 0 for none; ignored without `wall_model_pressure_gradient` ([turbulence](./turbulence#pressure-gradient)). |
| `wall_model_pressure_gradient_ypg` with `wall_model_pressure_gradient` | real | no (0) | `>= 0` | $y^*_{pg}$ of the Bernoulli-consistent variant (Hickel et al. 2012: 4), the height in viscous units of the gradient above which convection balances it; 0 for none; ignored without `wall_model_pressure_gradient` ([turbulence](./turbulence#bernoulli-consistent-variant)). |
| `r`, `u`, `p` with `inflow-digital-filter` | real | yes | `> 0`, otherwise fatal | Free-stream density, speed along the inward normal, pressure. **Euler** only, one such face per realm. |
| `df_profile` | string | with `inflow-digital-filter` | file must exist, 7 numbers per row, $y$ increasing, $T > 0$ (fatal) | The layer table: `y / delta`, `U / U_inf`, `T / T_inf`, `R_uu`, `R_vv`, `R_ww`, `R_uv` in units of $U_\infty^2$ (v wall-normal, w spanwise), `#` comments; the free stream beyond its last row. `tbl_oracle.py profile` writes it from a DNS station file. |
| `df_delta`, `df_spacing` | real | with `inflow-digital-filter` | `> 0` | Layer thickness (the table's unit of $y$); spacing of the filter plane, about the grid's. |
| `df_wall_axis` | string | with `inflow-digital-filter` | `x`, `y`, `z`, not the face normal (fatal) | Wall-normal axis, the wall on its minimum side; the span is the third axis (wrapped when periodic). |
| `df_length_streamwise`, `df_length_normal`, `df_length_spanwise` | real | with `inflow-digital-filter` | `> 0` | Integral lengths in units of `df_delta`; the streamwise one becomes the time scale $L_x / U_\infty$. |
| `df_seed` | integer | no (1) | | Seed of the random stream (the same plane on every rank). |
| `r`, `u`, `v`, `w`, `p` with `non-reflecting` | real | yes | `r, p > 0` (fatal) | The target state: its pressure is the one the face relaxes to, its entropy and tangential velocity those of flow entering through it. **Euler** only (no dual energy, no passive scalars). |
| `nr_length` | real | with `non-reflecting` | `> 0` (fatal) | $L$ of the relaxation $K = \sigma(1 - M^2)c/L$, a length of the domain. |
| `nr_sigma` | real | no (0.28) | `>= 0` (fatal) | $\sigma$: 0 is extrapolation at an outflow; larger holds the pressure faster and reflects more. |
| `r`, `u`, `v`, `w`, `p` with `shock-generator` | real | yes | `r, p > 0`; supersonic in the x-y plane (fatal) | The state upstream of the shock. **Euler** only, `bc_y_min` and `bc_y_max` only. |
| `shock_deflection`, `shock_angle` | real | exactly one, with `shock-generator` | degrees; $0 < \theta <$ detachment, $\mu < \beta \le 90$ (fatal) | The flow turning (the weak shock) or the shock angle from the upstream velocity. |
| `shock_x`, `shock_x_end` | real | `shock_x` yes, `shock_x_end` no (no end) | `shock_x_end > shock_x` (fatal) | Where the shock enters: the ghosts with `shock_x` $\le x <$ `shock_x_end` hold the post-shock state, the others are extrapolated. |
| `wall_model_edge_velocity`, `wall_model_edge_temperature` | real | with `wall_model = inverse-transformation` | `> 0`, otherwise fatal | Speed and temperature at the edge of the layer (the free stream; the centreline of a channel), constants of the run. |

What each `type` does:

- `extrapolation`: zeroth-order extrapolation.
- `inflow`: prescribed state.
- `non-reflecting` ([#121](https://github.com/szaghi/adam/issues/121)): characteristic outflow and far field, the incoming acoustic wave relaxed toward the target pressure ([boundary conditions](./boundary-conditions#non-reflecting-a-characteristic-outflow-and-far-field)).
- `shock-generator` ([#121](https://github.com/szaghi/adam/issues/121)): the post-shock state of an oblique shock in a range of x on a y face, extrapolation elsewhere ([boundary conditions](./boundary-conditions#shock-generator-an-oblique-shock-through-a-y-face)).
- `inflow-digital-filter` ([#119](https://github.com/szaghi/adam/issues/119)): a turbulent boundary layer through the face ([turbulence](./turbulence#digital-filter-inflow)).
- `wall-inviscid`: slip wall (mirror state with the normal momentum negated; for MHD the normal B is negated too).
- `wall-noslip`: adiabatic no-slip wall; `wall-isothermal`: isothermal no-slip wall
  ([boundary conditions](./boundary-conditions#wall-noslip-and-wall-isothermal-no-slip-walls)).
- `wall-modelled` ([#119](https://github.com/szaghi/adam/issues/119)): isothermal wall with a wall-model flux ([turbulence](./turbulence#wall-model)).
- `periodic`: library periodicity, so true periodic neighbours across blocks and ranks.
- `hydrostatic` ([#73](https://github.com/szaghi/adam/issues/73)): closes a stratified atmosphere under `[gravity]`
  ([boundary conditions](./boundary-conditions#hydrostatic-an-atmosphere-under-gravity)).

In a forest, a face glued to another realm by the manifest still needs a `type` here. The seam overrides the ghost
fill.

---

## [initial_conditions]

`amr_iterations` and `type` are always required. The other keys depend on `type`, and **every key the
selected type uses is required** (fatal "failed to load" otherwise).

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `amr_iterations` | int | yes | Negative values are clamped to 0. | Number of initial-condition + AMR passes (each pass: set IC, run the markers). `> 0` requires even `ni/nj/nk`. |
| `type` | string | yes | See the table below. Other → fatal. | Initial condition kind. |
| `regions_number` | int | `riemann-problem` | `>= 1`, otherwise fatal. Ignored by the other types (they fix the count). | Number of `[initial_conditions_region_N]` sections. |

### Accepted `type` values and their extra keys

| `type` | Model | Regions read | Extra `[initial_conditions]` keys (all real unless stated) and checks |
|--------|-------|--------------|-----------------------------------------------------------------------|
| `uniform` | any | 1 | `s`: signed relative amplitude of a deterministic hashed perturbation of `rho` and `p`. |
| `boundary-layer` | **euler** | 0 | `bl_profile` (the table of `df_profile`), `bl_delta`, `bl_rho`, `bl_u`, `bl_p` (free stream), `bl_wall` (default 0) ([initial conditions](./initial-conditions#boundary-layer-mean-profile-of-a-boundary-layer)). |
| `channel` | **euler** | 0 | `ch_rho`, `ch_ub` (bulk velocity), `ch_mach` (on the wall sound speed), `ch_h` (half height; the walls at $y=\pm h$), `ch_lx`, `ch_lz`, `ch_pert` (amplitude over `ch_ub`) ([initial conditions](./initial-conditions#channel-plane-channel-laminar-profile-with-a-perturbation)). |
| `isentropic-vortex` | **euler** (fatal otherwise) | 1 (free stream) | `x0`, `y0`, `radius` (`> 0`), `strength`; optional `shock_x`: past it the post-shock state of a stationary shock, region 1 its supersonic upstream state with $u > 0$. |
| `isotropic-turbulence` | **euler, euler-de** (fatal otherwise) | 0 | `hit_density`, `hit_urms`, `hit_mach`, `hit_k0` (all `> 0`), `hit_seed` (default 1). The box $[0, 2\pi]^3$, its periodicity and the resolution of $3k_0$ are not checked ([initial conditions](./initial-conditions#isotropic-turbulence-decaying-compressible-isotropic-turbulence)). |
| `taylor-green` | **euler, euler-de** (fatal otherwise) | 0 | `tg_density`, `tg_velocity`, `tg_mach` (all `> 0`), `tg_length` (default 1). The box $[0, 2\pi L]^3$ and its periodicity are not checked ([initial conditions](./initial-conditions#taylor-green-compressible-taylor–green-vortex)). |
| `riemann-problem` | any | `regions_number` | `regions_number`. Each region also needs extents (see the next section). Every cell must be covered by a region (fatal otherwise, at set time), and the **first** matching region wins. |
| `glm-pulse` | **mhd-ideal** | 1 | `pulse_axis` (string `x`/`y`/`z`, fatal otherwise), `pulse_center`, `pulse_width` (`> 0`), `pulse_amplitude`. |
| `divb-peak` | **mhd-ideal** | 1 | `peak_x0`, `peak_y0`, `peak_radius` (`> 0`), `peak_amplitude`. |
| `mhd-linear-wave` | **mhd-ideal** | 1 | `wave` (string `fast`/`alfven`/`slow`/`entropy`, fatal otherwise), `wave_angle` (degrees, x-y plane), `wave_amplitude`, `wavelength` (`> 0`). |
| `mhd-cpaw` | **mhd-ideal** | 1 (u, v, w, bx, by, bz must be 0, fatal otherwise) | `polarisation` (string `right`/`left`, fatal otherwise), `wave_angle`, `wave_amplitude`, `wavelength` (`> 0`), `b_par` (`/= 0`). |
| `mhd-vortex` | **mhd-ideal** | 1 (bx, by must be 0, fatal otherwise) | `x0`, `y0`, `radius` (`> 0`), `kappa`, `mu`. |
| `orszag-tang` | **mhd-ideal** | **none** | No extra keys. The unit-period Orszag–Tang state is hard-coded (`gamma` comes from `[physics]`). |
| `mhd-rotor` | **mhd-ideal** | 1 (ambient) | `x0`, `y0`, `r0`, `r1` (`0 < r0 < r1`), `rho_in` (`> 0`), `v0`. |
| `field-loop` | **mhd-ideal** | 1 | `x0`, `y0`, `loop_radius` (`> 0`), `loop_amplitude`. |
| `rotated-riemann` | any | 2 (states given in the **normal frame**: `u`, `bx` normal; `v`, `by` tangential) | `normal_x`, `normal_y` (not both 0), `interface_1`, `interface_2`, `period` (`> 0`), `interface_2_width` (`>= 0`). Requires `interface_1 < interface_2` and `interface_2 + interface_2_width < interface_1 + period`, fatal otherwise. |
| `shu-osher` | **euler** | 2 | `axis` (string `x`/`y`/`z`), `interface`, `rho_amplitude`, `rho_wavenumber`. Requires `rho_amplitude < region_2 r`. Region 1 applies where `s <= interface`; region 2 applies elsewhere, with `r + rho_amplitude * sin(rho_wavenumber * s)`. |
| `linear` | any | 1 | `gradient_x`, `gradient_y`, `gradient_z` (inverse length). Every conservative variable is `q_1 (1 + g . x)`, `q_1` the state of region 1: the verification field of the ghost probe ([#65](https://github.com/szaghi/adam/issues/65)), not a flow. |
| `isothermal-atmosphere` | any | 1 (the state at the origin; its velocity and field uniform) | Needs `[gravity]` (fatal otherwise). Optional, both atmospheres: a relative pressure pulse, `pulse_amplitude` (default 0), `pulse_x`, `pulse_y`, `pulse_z`, `pulse_width` (`> 0` when the amplitude is not 0, fatal otherwise) ([initial conditions](./initial-conditions#isothermal-atmosphere-polytropic-atmosphere-hydrostatic-equilibria-under-gravity)). |
| `polytropic-atmosphere` | any | 1 (as above) | `polytropic_gamma` $G > 1$ (fatal otherwise); needs `[gravity]`. $p = p_1(\rho/\rho_1)^G$, $\rho = \rho_1 s^{1/(G-1)}$, $s = 1 + \frac{G-1}{G}\frac{\rho_1}{p_1}\mathbf g\cdot\mathbf x$; fatal where $s \le 0$ (the top of the atmosphere inside the domain). |
| `rayleigh-taylor` | any | 1 (the light layer at the interface pressure; its field uniform) | Needs `[gravity]` along one axis (fatal otherwise). `rt_height`, `rt_width`, `rt_density_ratio` (heavy over light; below 1 is the stable control), `rt_wavelength`, `rt_amplitude`, `rt_seed_width`; every key required, width, ratio, wavelength and seed width `> 0` (fatal otherwise) ([initial conditions](./initial-conditions#rayleigh-taylor-two-isothermal-layers-heavy-over-light)). |
| `polynomial` | any | 1 | `gradient_x`, `gradient_y`, `gradient_z`, `degree` (an integer $n \ge 1$). Every conservative variable is `q_1 (1 + sum_d g_d x_d^n)`: the ghost probe's field for the high-order seam operators ([#68](https://github.com/szaghi/adam/issues/68)), not a flow. Fatal under a non-dimensional reference (`g_d` has dimension $L^{-n}$). |
| `bonnor-ebert` | **euler**, `[eos] type = isothermal` | 1 (its velocity) | Needs `[gravity] self_gravity`. `be_density` (`> 0`, required), `be_x0`, `be_y0`, `be_z0` (default 0), `be_overdensity` (`> 0`, default 1) ([initial conditions](./initial-conditions#bonnor-ebert-the-isothermal-lane–emden-sphere-self-gravity)). |
| `boss-bodenheimer` | **euler**, `[eos]` isothermal or barotropic | 1 (the ambient medium) | Needs `[gravity] self_gravity`. `bb_radius`, `bb_density` (`> 0`, required), `bb_omega` (required), `bb_amplitude` (default 0.1), `bb_x0`, `bb_y0`, `bb_z0` (default 0) ([initial conditions](./initial-conditions#boss-bodenheimer-the-rotating-cloud-with-an-m-2-perturbation-self-gravity)). |

Five more types are accepted and documented, with their keys, in [initial conditions](./initial-conditions): `sine-wave`, `couette`, `becker-shock`, `polytrope` and `singular-isothermal-sphere`. A model mismatch (for example `glm-pulse` with `euler`) is fatal. The complete formulas are in the module header.
With [`[eos] type = isothermal` or `barotropic`](#eos-optional-equation-of-state) only `uniform`, `riemann-problem`,
`rotated-riemann`, `sine-wave`, `linear`, `polynomial`, `shu-osher` and `boss-bodenheimer` are accepted (with
`isothermal` also `singular-isothermal-sphere` and `bonnor-ebert`; with `mhd-ideal` also `glm-pulse`, `divb-peak`,
`mhd-linear-wave` without the entropy wave, `mhd-cpaw`, `mhd-rotor` and `field-loop`), and the regions need no `p`.

### [initial_conditions_region_N]

N runs over the regions the selected type reads (see the table above).

| Key | Type | Req. | Meaning |
|-----|------|------|---------|
| `r`, `u`, `v`, `w`, `p` | real | yes | Primitive state (density, velocity, pressure). No positivity check. |
| `bx`, `by`, `bz` | real | `mhd-ideal` | Magnetic field. `psi` starts at 0. |
| `s1` .. `sN` | real | no (`passive_scalars = N`) | Specific values (mass fractions) of the passive scalars, default 0. Each cell takes those of the first region whose box holds its centre; the closed-form types (vortex, Shu–Osher, ...) take region 1's everywhere. |
| `emin_x`, `emin_y`, `emin_z`, `emax_x`, `emax_y`, `emax_z` | real | `type = riemann-problem` only | Region box. A cell belongs to it iff `emin < centre <= emax` on every axis. |

---

## [time]

Every key is required, and none has a range check.

| Key | Type | Req. | Meaning |
|-----|------|------|---------|
| `it_max` | int | yes | `> 0`: stop after `it_max` steps (iteration-driven). `<= 0`: time-driven run. |
| `time_max` | real | yes | End time (time-driven runs). The last step is clipped to hit it. |
| `CFL` | real | yes | CFL number (note the upper-case key). |

---

## [diagnostics]

| Key | Type | Req. | Meaning |
|-----|------|------|---------|
| `conservation_history_save` | int | yes | Cadence (steps) of `<output_basename>-conservation_history.dat` (volume integrals of every conservative variable). For MHD it is also the cadence of `-divb_history.dat` and of the `divb_tol` monitor. `<= 0` disables them. |
| `ghost_poison` | logical | no (`.false.`) | Verification instrument (issue #65): before each field write every ghost cell is set to NaN and refilled in the order of a stage (seams, exchange, boundary conditions), so a ghost the fill misses, or reads before its donor is written, is written as NaN instead of a stale value. Used by the ghost probe GP; it costs one pass over the ghosts per write. |
| `divb_flux` | logical | no (`.false.`) | MHD with GLM or EGLM ([#68](https://github.com/szaghi/adam/issues/68)): adds four columns to `-divb_history.dat` (`max_dflux`, `l1_dflux`, `seam_max_dflux`, `max_h_dflux_over_bmax`), the divergence the cleaning acts on ([numerics](./numerics#output-and-diagnostics)). Fatal without GLM/EGLM or with a one-stage scheme. The step-0 row has no reflux substitution and no seam column. |
| `turbulence_history_save` | int | no (0, disabled) | Cadence (steps) of `<output_basename>-turbulence_history.dat`: the volume integrals `K`, `Z`, `eps_s`, `eps_d`, `eps`, `pd` of the kinetic-energy budget ([turbulence](./turbulence#kinetic-energy-budget)), in the output units of [`[reference]`](#reference-optional-dimensional-input) (`K` an energy $\rho_0u_0^2L_0^3$, `Z` $\rho_0u_0^2L_0$, the others a power $\rho_0u_0^3L_0^2$). Needs `ngc >= 2`. |
| `wall_history_save` | int | no (0, disabled) | Cadence (steps) of `<output_basename>-wall_<face>_history.dat`, one file per `wall-modelled` face: wall shear stress and heat flux in `wall_history_bins` bins along `wall_history_axis`, code units, appended on restart ([turbulence](./turbulence#wall-histories)). |
| `wall_history_bins` | int | with `wall_history_save` | Bins along the axis, `>= 1` (fatal otherwise). |
| `wall_history_axis` | string | no (`x`) | `x`, `y` or `z`, not the normal of a modelled face (fatal). |

---

## Forest manifest (multi-realm)

A multi-realm FLUME run reads a forest manifest, documented in the
[library input reference](/reference/library-input#forest-manifest-multi-realm); every `[realm.N]` points at a full
FLUME input of the form documented here. Each realm regrids by its own `[amr]` section, balanced across the seams by
the forest (see `[amr]` above and [the forest](/concepts/forest)); a forest restarts on the rank count it was saved
with.

---

## Complete example 1: Euler, Sod shock tube (verbatim)

The file is `src/tests/flume/verification/sod/sod-x.ini`, copied without changes.

```ini
; FLUME verification V1 (issue #35, section 11): Sod shock tube along x, 1-D (the other directions are null).
; 64 blocks, four along x (200 cells); the 8 sibling groups split 32/32 over two ranks (the library never splits a
; sibling group, so one refinement level would leave a rank empty). WENO-5 characteristic, SSP-33, CFL 0.5, t = 0.2.
; check.sh compares density with the exact Riemann solution (L1 bound) and the three directions with each other (bitwise).
[IO]
output_basename        = sod-x
it_save                = 100000
restart                = .false.
restart_basename       = sod-x-restart
restart_save           = 0
residuals_save         = 1
save_memory_status     = .false.
save_residual_fields   = .false.
save_auxiliary_fields  = .false.
save_curl_fields       = .false.
save_divergence_fields = .false.
save_gradient_fields   = .false.
save_laplacian_fields  = .false.

[grid]
ni     = 50
nj     = 4
nk     = 4
ngc    = 3
emin_x = 0.0
emin_y = 0.0
emin_z = 0.0
emax_x = 1.0
emax_y = 1.0
emax_z = 1.0
null_x = .false.
null_y = .true.
null_z = .true.

[amr]
max_level      = 2
ratio          = 8
iu_ref_levels  = 2
i_prune        = 0
j_prune        = 0
k_prune        = 0
l_prune        = -1
frequency      = 0
iters          = 1
markers_number = 0

[field]
nv = 5

[runge_kutta]
scheme = runge-kutta-ssp-33

[weno]
scheme              = weno-u-5
ror_number          = 0
ror_threshold       = 0.9
ror_vars_number     = 0
enable_ror_stats    = .false.
ib_reduction_extent = 0
ib_reduced_order    = 2

[fdv]
fdv_scheme = fd
fdv_order  = 2

[solids]
number    = 0
n_eikonal = 0

[slices]
slices_number = 0

[physics]
physical_model = euler
cp             = 1040.004
cv             = 742.86

[numerics]
scheme_space             = weno
reconstruction_variables = characteristic
reflux                   = .true.

[bc_x_min]
type = extrapolation
[bc_x_max]
type = extrapolation
[bc_y_min]
type = extrapolation
[bc_y_max]
type = extrapolation
[bc_z_min]
type = extrapolation
[bc_z_max]
type = extrapolation

[initial_conditions]
type           = riemann-problem
regions_number = 2
amr_iterations = 0

[initial_conditions_region_1]
r      = 1.0
u      = 0.0
v      = 0.0
w      = 0.0
p      = 1.0
emin_x = 0.0
emin_y = 0.0
emin_z = 0.0
emax_x = 0.5
emax_y = 1.0
emax_z = 1.0

[initial_conditions_region_2]
r      = 0.125
u      = 0.0
v      = 0.0
w      = 0.0
p      = 0.1
emin_x = 0.5
emin_y = 0.0
emin_z = 0.0
emax_x = 1.0
emax_y = 1.0
emax_z = 1.0

[time]
it_max   = -1
time_max = 0.2
CFL      = 0.5

[diagnostics]
conservation_history_save = 1
```

### Annotations for the Sod input

- **Grid and block count.** With `ratio = 8` and `iu_ref_levels = 2` the domain is split into `8^2 = 64` blocks of
  `50 x 4 x 4` cells. The run is 1-D: `null_y` and `null_z` make y and z inactive, so there are 4 blocks along x and
  200 cells in total along x.
- **`it_save = 100000`.** This is effectively "first and last step only", because the last step is always saved.
- **`restart_save = 0` and `residuals_save = 1`.** Restarts are disabled, and the residual history is written every
  step.
- **Keys FLUME ignores.** `[field] nv`, `fdv_scheme` and the `ror_*`/`ib_*` keys of `[weno]` are ignored by FLUME.
  All of them except `[field]` must still be present.
- **Riemann-problem regions.** The regions tile x at 0.5: region 1 owns `0 < x <= 0.5` and region 2 owns
  `0.5 < x <= 1`. The domain minimum corner is covered because `emin` is the domain minimum and cell centres are
  strictly larger.

---

## Complete example 2: ideal MHD with GLM, Orszag–Tang vortex with a 2:1 AMR box

There is **no committed INI** under `src/tests/flume/verification/mhd/orszag-tang/`. `check.sh` generates the input
with `make_orszag_tang.py` from `verification/vortex/vortex-n064.ini`. The example below is what the generator writes
for the AMR leg of `check.sh`:

```bash
cd src/tests/flume/verification
python3 mhd/orszag-tang/make_orszag_tang.py vortex/vortex-n064.ini orszag-tang.ini --cells 32 --refine-box 0.25 0.25 0.75 0.75
```

It is written by Python `configparser`, which explains the section order and the spacing.

```ini
[IO]
output_basename = orszag-tang
it_save = 1000000
restart = .false.
restart_basename = vortex-n064-restart
restart_save = 0
residuals_save = 1
save_memory_status = .false.
save_residual_fields = .false.
save_auxiliary_fields = .false.
save_curl_fields = .false.
save_divergence_fields = .false.
save_gradient_fields = .false.
save_laplacian_fields = .false.

[grid]
ni = 8
nj = 8
nk = 4
ngc = 3
emin_x = 0.0
emin_y = 0.0
emin_z = 0.0
emax_x = 1.0
emax_y = 1.0
emax_z = 1.0
null_x = .false.
null_y = .false.
null_z = .true.

[amr]
max_level = 3
ratio = 8
iu_ref_levels = 2
i_prune = 0
j_prune = 0
k_prune = 0
l_prune = -1
frequency = 0
iters = 1
markers_number = 1

[field]
nv = 9

[runge_kutta]
scheme = runge-kutta-ssp-54

[weno]
scheme = weno-u-5
ror_number = 0
ror_threshold = 0.9
ror_vars_number = 0
enable_ror_stats = .false.
ib_reduction_extent = 0
ib_reduced_order = 2

[fdv]
fdv_scheme = fd
fdv_order = 2

[solids]
number = 0
n_eikonal = 0

[slices]
slices_number = 0

[physics]
physical_model = mhd-ideal
cp = 2.5
cv = 1.5

[numerics]
scheme_space = weno
reconstruction_variables = characteristic
reflux = .true.

[bc_x_min]
type = periodic

[bc_x_max]
type = periodic

[bc_y_min]
type = periodic

[bc_y_max]
type = periodic

[bc_z_min]
type = extrapolation

[bc_z_max]
type = extrapolation

[initial_conditions]
type = orszag-tang
amr_iterations = 1

[time]
it_max = -1
time_max = 0.5
CFL = 0.4

[diagnostics]
conservation_history_save = 1

[mhd]
divergence_control = glm
divb_tol = 0.0
divb_error = .false.
rho_floor = 0.0
p_floor = 0.0
glm_ch = 4.0
glm_alpha = 0.18
glm_damping_length = 1.0
glm_ch_check = error

[amr_marker_1]
mode = 1
geo_type = primitive-box
delta_type = max
delta_fine = 0.0
delta_coarse = 0.0
box_xmin = 0.25
box_ymin = 0.25
box_zmin = -1e+30
box_xmax = 0.75
box_ymax = 0.75
box_zmax = 1e+30
target_level = 3
```

### Annotations for the Orszag–Tang input

- **Physics.** `cp/cv = 5/3`. `divergence_control = glm` gives `nv = 9` (`psi`); `[field] nv = 9` is ignored.
- **Time step and cleaning speed.** `glm_ch = 4` enters the time step. `glm_ch_check = error` stops the run if a fast
  wave ever outruns `c_h`.
- **Monitors and floors.** `divb_tol = 0` disables the div(B) monitor (on the normalised D_flux, which needs
  `[diagnostics] divb_flux`), but `-divb_history.dat` is still written every step. Both floors are 0, so a non-positive density or pressure is fatal.
- **Grid.** The base grid is `4 x 4` blocks (`iu_ref_levels = 2` on a quasi-2-D octree, with z null) of `8 x 8` cells,
  i.e. 32² cells. `ni` and `nj` are even, as `amr_iterations = 1` requires.
- **Refinement.** `amr_iterations = 1` runs one pass of the single marker. The marker is `mode = 1` +
  `geo_type = primitive-box`: it refines the blocks whose centroid lies in `[0.25, 0.75]^2` to level 3, one level
  above `iu_ref_levels`, capped by `max_level = 3`. `delta_type`, `delta_fine` and `delta_coarse` are required but
  unused for a box.
- **Boundaries and IC.** x and y are periodic (both faces of each axis). The IC has no region sections.
- **Restart basename.** `restart_basename = vortex-n064-restart` is inherited from the base input. It is harmless here
  because `restart_save = 0`.

## Derived output fields

Both groups are written to the XH5F output only when `[IO] save_auxiliary_fields = .true.`. The key is **required**. The auxiliary variables are recomputed on the host from the saved `q`,
including ghost cells. The derived MHD group is written only for the MHD models.

**Field files include the ghost cells**, filled in the order of a Runge–Kutta stage, so the written ghosts are those
the stencils read. In a forest the step-0 file holds unfilled inter-realm seam ghosts.

### Auxiliary variables

| name | index | Euler | MHD |
|---|---|---|---|
| `rho` | `IA_R` | $\rho$ | $\rho$ |
| `u`, `v`, `w` | `IA_U..IA_W` | $(\rho u)/\rho$, … | same |
| `p` | `IA_P` | $(\gamma-1)\big(E-\tfrac12\rho\vert \mathbf u\vert ^2\big)$ | $(\gamma-1)\big(E-\tfrac12\rho\vert \mathbf u\vert ^2-\tfrac12\vert \mathbf B\vert ^2\big)$ |
| `T` | `IA_T` | $p/(\rho R)$, with $R=c_p-c_v$ | same |
| `H` | `IA_H` | $(E+p)/\rho$ | $(E+p+\tfrac12\vert \mathbf B\vert ^2)/\rho$, which includes the magnetic pressure |
| `a` | `IA_A` | $\sqrt{\gamma p/\rho}$ | $\sqrt{\gamma p/\rho}$ (the **sound** speed, not the fast speed) |
| `Bx`, `By`, `Bz` | `IA_BX..IA_BZ` | — | copies of $B_x, B_y, B_z$ (capitalised to avoid clashing with the conservative `bx, by, bz`) |
| `s` | `IA_S` | `euler-de`: $(\rho s)/\rho$, the advected entropy per unit mass | — |
| `s1` .. `s4` | after `IA_S` | passive scalars: $(\rho s_k)/\rho$, all four written whatever `passive_scalars` | — |

With `[eos] type = isothermal` or `barotropic` (no energy equation): `p` $= P(\rho)$, `T` $= P/(\rho R)$, `a` $= \sqrt{dP/d\rho}$ and `H` $= 0$.

The conservative variables are named `r, ru, rv, rw, rE [, rs] [, rs1 .. rs4]` (Euler) or `r, ru, rv, rw, rE, bx, by, bz [, psi]` (MHD), and the residuals `dq_<name>`. With an inert energy (isothermal, barotropic) `rE` is carried and stays 0.

### MHD derived fields

| name | definition | where |
|---|---|---|
| `pt` | $p + \tfrac12\vert \mathbf B\vert ^2$ (total pressure, rationalised units) | every cell, ghosts included |
| `beta` | $2p/\vert \mathbf B\vert ^2$, or `huge(1._R8P)` $\approx 1.8\times10^{308}$ where $\vert \mathbf B\vert ^2 = 0$ | every cell |
| `bmag` | $\vert \mathbf B\vert  = \sqrt{B_x^2+B_y^2+B_z^2}$ | every cell |
| `divb` | $\sum_d w_d\,\delta_d B_d$ (signed; see below) | interior only, **zero in the ghost cells** |

**Discretisation of `divb`.** Each $\delta_d$ is the library centred first derivative,
$\delta B = \frac{1}{\Delta}\sum_{m=1}^{s}c_{m,s}\,(B_{i+m}-B_{i-m})$.
The coefficients are:

| $s$ | order | $c_{m,s}$ |
|---|---|---|
| 1 | 2 | $\tfrac12$ |
| 2 | 4 | $\tfrac{8}{12}, -\tfrac{1}{12}$ |
| 3 | 6 | $\tfrac{45}{60}, -\tfrac{9}{60}, \tfrac1{60}$ |
. The half stencil is $s$ = `fdv_half_stencils(1)` $=$ `[fdv] fdv_order / 2`,
so the formal order equals `fdv_order`. Every FLUME test input sets `fdv_order = 2`, which gives
$(B_{i+1}-B_{i-1})/(2\Delta x)$.

Two further rules apply:
- Null directions are weighted $w_d = 0$.
- The scheme is **always FD-centred**, even when `[fdv] fdv_scheme = fv`.

The same operator gives the div(B) history ($\max|\nabla\!\cdot\!\mathbf B|$, $\sum|\nabla\!\cdot\!\mathbf B|\,\Delta V$
and a seam-band maximum), which needs $s \le \min(n_{gc}, 3)$.
