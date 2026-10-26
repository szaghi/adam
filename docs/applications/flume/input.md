# Input reference

A FLUME run reads one INI file: `mpirun -np N exe/adam_flume_cpu input.ini` (or `adam_flume_fnl`). This page lists
every section and key the application reads, with its type, whether it is required, its default, the accepted values
and what happens on an invalid one. The two complete inputs at the end are verification cases, copied verbatim; the
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
| `[amr]` | `frequency > 0` regrids at runtime ([#74](https://github.com/szaghi/adam/issues/74)): negative is fatal, and so are runtime regridding without markers and a multi-realm forest with `frequency > 0`. The FNL backend regrids by a host round trip (each regrid logs its wall time). `regrid_prolongation` defaults to `conservative` in FLUME (the library default is `linear`); `conservative-state` ([#88](https://github.com/szaghi/adam/issues/88)) keeps a runtime-AMR stellar collapse admissible ([SV-15](./verification/self-gravity/polytrope)) but is less accurate on smooth flow ([AV](./verification/runtime-amr/accuracy)). Derefinement is always the mean of the children. |
| `[amr_marker_N]` | Modes 1 (`solid`, `primitive-box`), 2, 4 and 5; mode 3 and an unknown mode are fatal at the first AMR pass, `stl` fatal. A Löhner marker (mode 4) acts at the runtime regrids only, with `buffer <= ngc - 1` (fatal otherwise); a Jeans marker (mode 5) needs `[gravity] self_gravity = .true.` (fatal otherwise). At a runtime regrid the markers vote per block: refine when any asks, coarsen only when all agree, keep otherwise, never below `iu_ref_levels`; the box marker votes to coarsen outside its box. The gradient is centred, across the faces between blocks too ([#76](https://github.com/szaghi/adam/issues/76)), one-sided at the domain boundary. |
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
  whose copies scale with the blocks in use, about 0.02 s per regrid on a 50-block grid on the development box
  ([#75](https://github.com/szaghi/adam/issues/75)).
- **Limits.** Single-realm runs only (a multi-realm forest with `frequency > 0` is refused); with EGLM the 2:1 seams
  raise the divergence error B_z well above its uniform-grid level ([#78](https://github.com/szaghi/adam/issues/78)).
- **Restart.** A run restarted from any step regrids on the same steps and continues bitwise (the regrid happens before
  the step's output and restart files).
- **Log.** Each regrid prints `flume: regrid at step N: blocks a -> b, R refined, C coarsened, S sweeps` (plus its wall
  time on FNL); a state with non-positive density or pressure after a regrid stops the run with its location.

---

## [reference] (optional: dimensional input)

FLUME's equations carry no dimensionless number (`B` is `B_SI/sqrt(mu0)`), so a dimensional input runs once every
value is divided by the reference of its dimension (issue #49). With this section present, FLUME does that once, on
the loaded file, before any other section is parsed; without it nothing changes. References: density `rho0`, length
`L0`, velocity `u0`; derived: time `L0/u0`, pressure and energy density `rho0 u0^2`, field `u0 sqrt(rho0)`, the GLM
`psi` `u0^2 sqrt(rho0)` (the EGLM one as the field). `cp`, `cv` are replaced by `gamma = cp/cv` (the gas constant
becomes 1, the temperature unit `u0^2/R`). Restart files and logs stay in code units; the other outputs follow
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

The dissipative coefficients and the wall keys (issue #49 N4, issue #65) convert as follows; the numbers are
dimensionless and pass unchanged. $R$ is the gas constant of the input ($c_p - c_v$, or 1 with `gamma`), since the layer
replaces `cp`, `cv` by `gamma` and the code temperature becomes $p/\rho$:

| Key | Divided by |
|---|---|
| `viscosity` | $\rho_0 u_0 L_0$ |
| `conductivity` | $\rho_0 u_0 L_0 R$ |
| `resistivity` | $u_0 L_0$ |
| `wall_u`, `wall_v`, `wall_w` | $u_0$ |
| `wall_temperature`, `reference_temperature` | $u_0^2 / R$ |
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
| `physical_model` | string | yes | `euler` (nv = 5), `euler-de` (nv = 6) or `mhd-ideal` (nv = 8, or 9 with GLM or EGLM). Other → fatal. | Equation set. `euler-de` is Euler with dual energy, the advected entropy `rs` ([#91](https://github.com/szaghi/adam/issues/91), [numerics](./numerics#dual-energy); CPU and FNL). `mhd-ideal` also loads `[mhd]`. |
| `passive_scalars` | int | no (`euler`, `euler-de`) | `0..4`, default 0; above 4 or with `mhd-ideal` → fatal | Passive scalars ([#92](https://github.com/szaghi/adam/issues/92), [numerics](./numerics#passive-scalars)): `N > 0` selects the scalar variant of the model, 4 more conserved variables `rs1`..`rs4` (`nv` 9, or 10 with `euler-de`), the first `N` set by the `s1`..`sN` keys of the initial-condition regions and of the inflows, the others zero. CPU and FNL. |
| `dual_energy_eta` | real | no (`euler-de`) | `(0, 1)`, default 0.05; otherwise fatal | Dual-energy switch: a cell whose internal energy is below `dual_energy_eta` times its total energy takes the entropy as the floor of its internal energy. |
| `gamma` | real | cond. | `> 1`; with `cp` or `cv` → fatal | Specific heats ratio; the gas constant is 1 (code units). |
| `cp` | real | cond. | `cp > cv > 0`, otherwise fatal; required without `gamma` | Specific heat at constant pressure. |
| `cv` | real | cond. | as above | Specific heat at constant volume. `gamma = cp/cv`, `R = cp - cv`. |

### Dissipative terms (issue #65, M4)

Each term is given either as its coefficient or as the dimensionless number it stands for. In code units the
references are 1, so a number is the reciprocal coefficient. Giving two keys of a term is fatal; giving none leaves the
term off, so an input without these keys stays ideal and bitwise unchanged. The viscosity and the conductivity drive
the Navier–Stokes fluxes of every model, the resistivity the Ohmic ones of `mhd-ideal`
([numerics](./numerics#dissipative-fluxes-navier–stokes)).

| Key | Type | Accepted / invalid | Meaning |
|-----|------|--------------------|---------|
| `viscosity` | real | `>= 0`; with `reynolds` → fatal | Dynamic viscosity $\mu$ (code units, or dimensional with [`[reference]`](#reference-optional-dimensional-input)). |
| `reynolds` | real | `> 0` | $Re$: $\mu = 1/Re$. |
| `conductivity` | real | `>= 0`; with `prandtl` → fatal | Thermal conductivity $k$; the heat flux is $-k\nabla T$ with $T = p/(\rho R)$. |
| `prandtl` | real | `> 0`; needs `viscosity` or `reynolds` | $Pr$: $k = \mu c_p / Pr$ (follows $\mu(T)$ under the power law). |
| `resistivity` | real | `>= 0`; MHD only; with a number → fatal | Magnetic diffusivity $\eta$ (Ohmic resistivity). |
| `magnetic_reynolds` | real | `> 0`; MHD only | $Rm$: $\eta = 1/Rm$. |
| `lundquist` | real | `> 0`; MHD only; with `[reference]` it needs `velocity = alfvenic` | $S$: $\eta = 1/S$, the magnetic Reynolds number at the Alfvén speed. |
| `viscosity_law` | string | `constant` (default), `power-law`; other → fatal; `power-law` needs a viscosity | Temperature law of the viscosity. |
| `viscosity_exponent` | real | `>= 0`; required by `power-law`, fatal without it | $\omega$ in $\mu(T) = \mu\,(T/T_\mathrm{ref})^\omega$ (non-negative: the diffusive time step bound evaluates the laws at the largest temperature a face can see). |
| `reference_temperature` | real | `> 0`; required by `power-law`, fatal without it | $T_\mathrm{ref}$ (code temperature $p/(\rho R)$). |

`[numerics] positivity_limiter = cell` is refused with any non-zero coefficient: the limiter's first-order backbone
with a central dissipative flux is not admissible without the extension of Zhang (2017, *J. Comput. Phys.* 328), which
is not implemented (D-M4-5). A coefficient is also fatal with immersed solids (their walls are inviscid), and
`dissipative_order = 4` with `[grid] ngc < 3`.

---

## [eos] (optional: equation of state)

Without the section the gas is the ideal gas of `[physics]` (`gamma`, or `cp` and `cv`). `type = isothermal`
([#100](https://github.com/szaghi/adam/issues/100), M13 P1) selects the isothermal gas $P = c_0^2\rho$, and
`type = barotropic` (P3) a piecewise polytrope $P = K_k\rho^{\Gamma_k}$, continuous across its density breaks
([numerics](./numerics#isothermal-and-barotropic-equations-of-state)). Both drop the energy equation: the state keeps
its five slots, the energy slot `rE` is inert (zero flux, zero source) and stays 0, and the checkpoints still write
it. The auxiliary fields report $p = P(\rho)$, $T = P/\rho$ (with $R = 1$), $a = \sqrt{dP/d\rho}$, $H = 0$.

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `type` | string | yes | `ideal`, `isothermal`, `barotropic`; other → fatal | The equation of state. `ideal` is the default gas of `[physics]`. |
| `sound_speed` | real | with `isothermal`, `barotropic` | `> 0`; missing or non-positive → fatal | $c_0$: the sound speed of the isothermal gas; for a table, the isothermal speed $\sqrt{P/\rho}$ at `reference_density`, so $K_1 = c_0^2\rho_\mathrm{ref}^{1-\Gamma_1}$ (the sound speed there when $\Gamma_1 = 1$). |
| `reference_density` | real | no (1) | `> 0`, otherwise fatal | $\rho_\mathrm{ref}$ of `barotropic`. |
| `densities` | list of reals | no | all `> 0`, strictly increasing; otherwise fatal | The breaks $\rho_1 < \dots < \rho_{n-1}$ of `barotropic`; a break belongs to the segment above it. |
| `exponents` | list of reals | no | one per segment (one more than `densities`), all `>= 1`, at most 8 segments; otherwise fatal | $\Gamma_1, \dots, \Gamma_n$ of `barotropic`. $K_{k+1} = K_k \rho_k^{\Gamma_k - \Gamma_{k+1}}$ keeps $P$ continuous. |

Without `densities` and `exponents`, `barotropic` takes the table of Mathew & Federrath (2021, MNRAS 507, 2448),
$\Gamma = 1, 1.1, 1.4, 1.1, 5/3$ with breaks at $2.5\cdot10^{-16}$, $3.84\cdot10^{-13}$, $3.84\cdot10^{-8}$ and
$3.84\cdot10^{-3}$ g cm$^{-3}$: a cgs input (their $c_0$ is 0.2 km/s). A break where $\Gamma$ decreases upward (that
table's $1.4 \to 1.1$) is not genuinely nonlinear: its Riemann solutions contain composite waves. It is accepted (the
schemes do not depend on convexity) and reported at start-up.

With `type = isothermal` or `barotropic`:
- `[physics]` must be `physical_model = euler` or `mhd-ideal` (P5, below) without `passive_scalars`, and must not set
  `gamma`, `cp` or `cv` (fatal: they would be ignored);
- the regions of `[initial_conditions]` and the `inflow` faces need no `p` (it is ignored: $P = P(\rho)$);
- the accepted initial conditions are `uniform`, `riemann-problem`, `rotated-riemann`, `sine-wave` (with
  `wave_mode = acoustic` the wave travels at $\sqrt{P'(\rho_0)}$, with `wave_mode = jeans` it obeys
  $\omega^2 = P'(\rho_0) k^2 - 4\pi G\rho_0$), `linear`, `polynomial`, `shu-osher` and `boss-bodenheimer`; with
  `isothermal` also `singular-isothermal-sphere` ($a^2 = c_0^2$) and `bonnor-ebert`
  ([initial conditions](./initial-conditions)); the others, whose density needs $\gamma$, are fatal;
- gravity (P4): `[gravity]`, external and self-gravity, acts on the momentum only, the energy slot stays 0 (no energy
  source, so no fourth-order energy correction and no `ngc >= 4` requirement, no seam gravity work in the AMR reflux);
  the Jeans marker (`[amr_marker_*] mode = 5`) and the sinks take $c_s^2 = dP/d\rho$ of the cell, and the sinks'
  Jeans test (`|E_grav| > 2 E_th`) the virial thermal energy $E_\mathrm{th} = \tfrac{3}{2}\int P\,dV$ (equal to the
  internal energy at $\gamma = 5/3$);
- `[numerics] positivity_limiter = cell` (P4) is the density-only limiter: the models' admissible quantity is the
  density, kept above `POSITIVITY_LIMITER_KAPPA` times its first-order backbone value (the stages and the AMR reflux),
  with a kinetic chord that keeps the velocity of a limited cell near its backbone's
  ([numerics](./numerics#isothermal-and-barotropic-equations-of-state)); the sinks need it, since accretion drains
  the cells within the accretion radius to the threshold density;
- refused, each with its own message: `[cooling]` (nothing to radiate), the dissipative terms, `[reference]` (it
  converts no barotropic law), immersed solids, and the `wall-noslip`, `wall-isothermal` and `hydrostatic` faces;
- CPU and FNL run the `weno` and `weno-riemann` schemes, all reconstructions, and the `llf`, `hll`, `hllc` solvers;
  the FNL kernels are the device twins of the CPU ones, on the same library;
- with `physical_model = mhd-ideal` (P5) the gas closes ideal MHD without an energy equation, the total pressure
  $P(\rho) + |\mathbf B|^2/2$ ([numerics](./numerics#isothermal-and-barotropic-mhd)):
  - `[mhd] divergence_control` must be `none` or `glm`. `eglm` is fatal, because its sources are built on the energy
    equation.
  - `rho_floor` and `p_floor` must be 0 (fatal otherwise: there is no pressure to floor).
  - The accepted initial conditions are those above plus `glm-pulse`, `divb-peak`, `mhd-linear-wave` (the eigenvectors
    of Kim et al. 1999; `wave = entropy` is fatal, since these gases have no entropy mode), `mhd-cpaw`, `mhd-rotor` and
    `field-loop`. `orszag-tang` and `mhd-vortex` are fatal: their states need $\gamma$.
  - The schemes are `weno` and `weno-riemann`, with the `llf`, `hll` and `hlld` solvers (Mignone's isothermal HLLD,
    on any $P(\rho)$).
  - `positivity_limiter = cell` works as for the Euler gases, and is fatal with `glm`, as for the ideal gas.
  - CPU and FNL run them (P5a, P5b); the FNL kernels are the device twins of the CPU ones, on the same library.

---

## [mhd] (only with `physical_model = mhd-ideal`)

The section is not read for `euler`.

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `divergence_control` | string | yes | `glm` (nv = 9, psi added), `eglm` (nv = 9, psi also in the energy; the `glm_*` keys apply) or `none` (nv = 8). Other → fatal. | div(B) control. Under `eglm` the AMR reflux leaves B alone (ρ, ρu, E and psi are refluxed): with the Derigs sources B is not a conserved integral, and refluxing it broke the sources' cancellation of the flux-generated $B_z$ at 2:1 seams ([#68](https://github.com/szaghi/adam/issues/68) G3b, #78). |
| `glm_ch` | real | `glm`, `eglm` | `> 0`, otherwise fatal | Constant cleaning speed `c_h`. It also bounds the time step (`c_h * sum 1/dx`). Keep it just above the fastest signal speed: a larger one speeds up the growth of small asymmetries in sheared, magnetised flows ([#83](https://github.com/szaghi/adam/issues/83)). |
| `glm_alpha` | real | `glm`, `eglm` | `>= 0`, otherwise fatal | Damping parameter: `c_h^2/c_p^2 = glm_alpha * c_h / L`. |
| `glm_damping_length` | real or `min-cell` | `glm`, `eglm` | A positive real, or the literal `min-cell` (the minimum cell spacing over the active directions, MPI-reduced; with runtime regridding, `[amr] frequency > 0`, the spacing of the finest level `max_level` allows, fixed for the run). A non-number or a value `<= 0` is fatal. | Damping length `L`. |
| `glm_ch_check` | string | `glm`, `eglm` | `warning` or `error`. Other → fatal. Without GLM it is forced to `warning` and not read. | What happens when `max(\|u\|+c_f) > glm_ch`: warn (logged each time a new maximum appears) or stop. |
| `divb_tol` | real | yes | `>= 0`, otherwise fatal. `0` disables the monitor; `> 0` needs `[diagnostics] divb_flux = .true.` (hence GLM or EGLM), fatal otherwise. | Monitor on the normalised D_flux ([#68](https://github.com/szaghi/adam/issues/68) G1): `max h\|D_flux\| / max\|B\| > divb_tol` warns or stops. Dimensionless (no `[reference]` conversion): the divergence of the interface field the cleaning acts on, scaled by the spacing and the field, bounded on a resolved run with an O(H) seam floor (A2). |
| `divb_error` | logical | yes | logical | `.true.`: exceeding `divb_tol` is fatal. |
| `rho_floor` | real | yes | `>= 0`, otherwise fatal. `0` disables it. | Density positivity floor. |
| `p_floor` | real | yes | `>= 0`, otherwise fatal. `0` disables it. | Pressure positivity floor. With **both** floors at 0, any non-positive density or pressure stops the run. |

---

## [gravity] (optional: external gravity and self-gravity)

A constant acceleration $\mathbf g$ ([#73](https://github.com/szaghi/adam/issues/73), milestone M7) and the
self-gravity of the gas ([#81](https://github.com/szaghi/adam/issues/81), milestone M8). Without the section there is
no gravity and nothing changes.

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `g_x`, `g_y`, `g_z` | real | one of the two forms | any; each optional (default 0) | The acceleration. Dimensional (length / time²) with [`[reference]`](#reference-optional-dimensional-input). |
| `froude` | real | the other form, with `axis` | `> 0`, otherwise fatal | The number form: $g = 1/\mathrm{Fr}^2$ along `axis` ($g L_0/u_0^2$ in the reference units; not converted by `[reference]`). |
| `axis` | string | with `froude` | `+x`, `-x`, `+y`, `-y`, `+z`, `-z`; other → fatal | The signed direction of $\mathbf g$. |
| `equilibrium` | string | no (default `none`) | `none` or `atmosphere`; other → fatal. `atmosphere` needs an `isothermal-atmosphere`, `polytropic-atmosphere` or `rayleigh-taylor` initial condition, and is fatal with `[numerics] reflux_policy = troubled` and on multi-realm runs | The well-balanced mode (M7-G2b): the initial atmosphere at rest (zero velocity, region 1's field, no pulse or seed) is kept to round-off. See below. |
| `potential` | string | no (default `uniform`) | `uniform` or `linear`; other → fatal. `linear` is fatal with `equilibrium = atmosphere` and with `self_gravity` | How g reaches the sources. `uniform`: the constant acceleration in every kernel, with the fourth-order energy source. `linear` ([#81](https://github.com/szaghi/adam/issues/81), M8-P1): the same field written as the cell potential Φ = −g·x and applied through the potential path of self-gravity (face gravity −ΔΦ/Δx, second order); a test of that path ([SV-1](./verification/self-gravity/potential)). |
| `self_gravity` | logical | no (default `.false.`) | with it: `equilibrium = atmosphere`, `potential = linear`, multi-realm runs, immersed solids and a build without an FFT library (the `*-fft` modes) are fatal, and so are the boundaries `poisson_bc` cannot hold (below) | Self-gravity (M8-P3): $L_h\Phi = 4\pi G(\rho - \bar\rho)$ (periodic) or $L_h\Phi = 4\pi G\rho$, $\Phi \to 0$ at infinity (isolated, M8-P4), solved every Runge–Kutta stage, applied through the potential path. See below. |
| `G` | real | one of `G`, `G_star` with `self_gravity` | `> 0`, otherwise fatal; fatal without `self_gravity`, and together with `G_star` | The gravitational constant. Dimensional with [`[reference]`](#reference-optional-dimensional-input). |
| `G_star` | real | the other form | as `G` | The number $G\rho_0L_0^2/u_0^2$; not converted by `[reference]`. |
| `poisson_solver` | string | no (default `auto`) | `auto`, `fft` or `multigrid`; other → fatal. `fft` is fatal on an AMR grid | `auto`: the FFT solver on a uniform grid, the multigrid solver on an AMR grid (leaves at more than one level, M8-P5b). `fft`: the library FFT solver only (FFTW on the CPU, cuFFT on FNL). `multigrid`: the multigrid solver on any grid (periodic, or isolated by the James step of [#87](https://github.com/szaghi/adam/issues/87)), a check against `fft` on a uniform one. |
| `poisson_bc` | string | no (default `periodic`) | `periodic` or `isolated`; other → fatal. `periodic` is fatal with a non-periodic active direction; `isolated` is fatal in 2-D, with a periodic direction, with non-cubic cells, and on AMR grids with a leaf above the coarsest level on the domain boundary (at start or at the regrid that makes one, [#87](https://github.com/szaghi/adam/issues/87)) | `periodic`: every active direction periodic, the mean density removed. `isolated` (M8-P4): the free-space potential of the mass on the grid (Hockney–Eastwood zero padding with the lattice Green's function of $L_h$; on AMR grids the multigrid's symmetric James step, [#87](https://github.com/szaghi/adam/issues/87)); the fluid boundaries are the `[bc_*]` ones, e.g. `wall-inviscid` for a closed box. |

Giving `g_*` together with `froude` is fatal, as is a section with neither (unless `self_gravity`, where $\mathbf g$
defaults to 0).

$\mathbf g$ sets the `isothermal-atmosphere` and `polytropic-atmosphere` initial conditions and the `hydrostatic`
boundary condition (M7-G0), and drives the sources (M7-G1):

- momentum $\rho\mathbf g$; energy from the face mass fluxes the continuity equation uses,
  $\tfrac1{24}\sum_d g_d\,(-F^\rho_{d,i-3/2} + 13F^\rho_{d,i-1/2} + 13F^\rho_{d,i+1/2} - F^\rho_{d,i+3/2})$, so that with
  this static potential the total energy $E + \rho\Phi$ ($\Phi = -\mathbf g\cdot\mathbf x$) is conserved exactly
  (Mikami et al. 2008, Mullen et al. 2021: each face's weights sum to one over the cells) with periodic boundaries and
  across block interfaces; at a `wall-inviscid` wall along $\mathbf g$ the source reads the mirrored face beyond the
  wall, and the exactly conserved total is $\sum(E + \rho\Phi)V \pm \sum g_n\Delta x_n\rho V/12$ over the cells
  beside the wall (+ at the lower wall, − at the upper), which differs from the plain sum by $O(\Delta x^2)$. The fluxes of a
  finite-difference scheme are point values of the primitive function of $\rho\mathbf u$, not of $\rho\mathbf u$
  itself, so the two-face average $\tfrac12(F^\rho_{i-1/2} + F^\rho_{i+1/2})$ alone is a second-order cell value and
  takes the whole scheme to second order; the four-face form is fourth order (M7-G2a,
  [verification](./verification/gravity/energy-order)). It reads the face beyond each block end, so
  **gravity needs `[grid] ngc >= S + 1`** ($S$ the WENO half-stencil, 3 for WENO-5; with `weno-riemann`
  $\max(S, 2) + 1$), fatal otherwise;
- with `positivity_limiter = cell`, the limiter's first-order backbone takes the time-centred energy
  $\rho\mathbf g\cdot(\mathbf u + \mathbf g\,\Delta t/2)$, which leaves its internal energy unchanged, and limits the rest
  of the face form with the cell factor $\Lambda$: gravity never makes the backbone inadmissible (near vacuum, or with
  $g\Delta t$ above the sound speed), and $E + \rho\Phi$ stays exact wherever the limiter does not act;
- with `equilibrium = atmosphere` (M7-G2b, the deviation method of Edelmann et al. 2021, §4.3), the residual of the
  atmosphere at rest is computed once per grid (at the first step, after a restart, after every regrid) with the
  run's own space operator and subtracted from every residual: the atmosphere is then an exact discrete steady state
  (round-off), on `weno` and `weno-riemann`, Euler and MHD, across AMR 2:1 seams (the equilibrium's reflux is folded
  into its residual), and a perturbation of it evolves without the stationary truncation flow of the unbalanced
  scheme ([verification](./verification/gravity/well-balanced)). It costs one more field (the stored
  residual) and one residual evaluation per grid. The subtracted residual is outside the positivity limiter's factor:
  it is the truncation error of the equilibrium, small against the state wherever the atmosphere is resolved;
- no time-step limit: the first step logs $|\mathbf g|\,\Delta t^2/(\mathrm{CFL}\,h_{\min})$, at least the velocity
  gravity adds in a step over the fastest wave speed (near 1 or above, gravity is stiff);
- refused (fatal) along a null direction and with immersed solids.

With `self_gravity` ([models](./models#self-gravity)):

- every Runge–Kutta stage solves $\Phi$ from that stage's density after the ghost update, with the library FFT solver:
  distributed over the z-planes of the transform array ([#84](https://github.com/szaghi/adam/issues/84) E5d; bitwise the same on every rank count),
  $\Phi$ returned with every ghost layer filled. Periodic: it divides each Fourier
  mode by the eigenvalue of the 7-point Laplacian $L_h$, so the discrete residual is round-off. Isolated: the density
  sits in a zero-padded array of $2(N + n_{gc})$ cells per direction and is convolved with the lattice Green's function
  of $L_h$, whose residual is about $10^{-13}$ per offset (`scripts/lattice_green.py`); the residual of $\Phi$ is that
  summed over the source, $2\cdot10^{-11}$ of $4\pi G\rho$ at $64^3$, and the ghosts beyond the domain hold the
  free-space potential. The solve is serial on rank 0 and, on FNL, the density and $\Phi$ cross the host every stage:
  correct, not scalable;
- on an AMR grid (M8-P5b) the multigrid solver (`adam_poisson_mg_object`) solves the composite Laplacian, whose 2:1
  faces carry the mimetic mean flux (symmetric and consistent; [SV-9](./verification/self-gravity/multigrid-amr)):
  flexible CG preconditioned by a V-cycle, its coarsest level solved by the FFT solver of the backend; on FNL it runs
  on the host. The face gravity at a 2:1 face is the mean-flux one on both sides (first order there, the potential
  second order), and the gravity work across each seam rides on the reflux (below), so $E + W$ stays exact in the
  semi-discrete form;
- across a 2:1 face the reflux, which corrects the coarse cells once at the end of the step, also hands the coarse cell
  the seam's gravity work, so that summed over the cells the energy sources equal $-\dot W$ exactly. Its stages,
  though, advance the coarse cells with their own seam fluxes: on AMR the drift of $E + W$ is first order in
  $\Delta t$, not of the Runge–Kutta order (SV-10); with `[numerics] reflux_cadence = stage` the correction is made
  every stage and the drift is of the Runge–Kutta order again. `reflux = .false.` or `reflux_policy = troubled` leave
  the work uncorrected where they leave the flux;
- under runtime AMR (`[amr] frequency > 0`, P5c) every regrid sets the solver up again on the new grid (the FFT solver
  when the grid is uniform, the multigrid otherwise) and solves $\Phi$ from the regridded state. The regrid conserves
  $\rho$, $\rho\mathbf u$ and $E$, but $W$ is the grid's own discrete energy: $E + W$ jumps at every regrid that
  changes the grid, by a truncation-sized amount ($2\cdot10^{-4}$ of $E + W$ on the SV-11 collapse at base $32^2$);
- the conservation history gains the columns `W` $= \tfrac12\sum\rho\Phi\,dV$ and `E_plus_W`, the energy the scheme
  conserves (each row solves $\Phi$ from the committed state);
- the first step logs $G$, the free-fall time $t_{ff} = \sqrt{3\pi/(32G\rho_{\max})}$ of the densest cell, $\Delta t/t_{ff}$
  (a warning above 0.1) and the wall time of one solve;
- a constant $\mathbf g$ may be given too; it adds through its own kernels;
- with `positivity_limiter = cell`, the limiter's gravity energy correction (above) acts on the potential path too:
  $E + W$ stays exact only where the limiter does not act. A star in a light atmosphere ($10^{-4}$ of its central
  density) needs the limiter at its surface and drifts by $2\cdot10^{-3}$ of $E + W$ over a free-fall time at $32^3$,
  against $10^{-9}$ without it (the SV-7 and SV-8 cases use a $10^{-2}$ atmosphere, which runs without it);
- **build:** `flume-cpu-gnu-fft` (FFTW 3, executable `adam_flume_cpu_fft`) or `flume-fnl-nvf-fft` (cuFFT,
  `adam_flume_fnl_fft`); the default modes refuse `self_gravity` at the first solve.

## [cooling] (optional: optically thin radiative cooling)

The section switches on radiative cooling ([#93](https://github.com/szaghi/adam/issues/93), M11):
$\partial_t(\rho e) = -n^2\Lambda(T)$ at fixed density and momentum, with $n = \rho/m$ and $T = p/(\rho R)$. It is
applied once per step, after the step, by Townsend's exact integration ([numerics](./numerics#radiative-cooling)).
Without the section nothing changes. The hydrodynamic models only (`euler`, `euler-de`, with or without passive
scalars): with `mhd-ideal` it is fatal. CPU and FNL run the same rule and the same arithmetic.

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

The section switches on sink particles ([#94](https://github.com/szaghi/adam/issues/94), M12): point masses that act on
the gas through the self-gravity solve only ([numerics](./numerics#sink-particles)). It needs `[gravity] self_gravity =
.true.` and a hydrodynamic model (`euler`, `euler-de`, with or without passive scalars); otherwise it is fatal. CPU
and FNL (on FNL, creation copies the state to the host every step: its checks run on the host). The keys list the
initial sinks, one value each; their
number is that of `mass` (none: the section with no initial sink).

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `mass` | list of reals | no (no sink) | all `> 0`; otherwise fatal | The sinks' masses. |
| `x`, `y`, `z` | lists of reals | with `mass` | as many values as `mass`, inside the domain; otherwise fatal | The positions. |
| `u`, `v`, `w` | lists of reals | no (0) | as many values as `mass`; otherwise fatal | The velocities. |
| `creation` | logical | no (`.false.`) | | Sinks form where the gas passes the checks of Federrath et al. (2010) ([numerics](./numerics#sink-particles)). |
| `accretion` | logical | no (`.true.`) | | Sinks accrete the gas above the threshold within their accretion radius. |
| `merging` | logical | no (`.true.`) | | Two bound, converging sinks within an accretion radius merge. |
| `accretion_radius_cells` | real | no (2.5) | in `[1, ngc]`; otherwise fatal | The accretion radius $r_\mathrm{acc}$ in cells of the sink's level (with runtime regridding, `[amr] frequency > 0`, of the finest level allowed, where the sinks form and stay). The threshold density is $\rho_\mathrm{thr} = \pi c_s^2/(4 G r_\mathrm{acc}^2)$, a Jeans length of $2 r_\mathrm{acc}$, with the cell's own $c_s^2 = \gamma p/\rho$. |

The conservation history gains `M_sinks`, `Px_sinks`, `Py_sinks`, `Pz_sinks` and `K_sinks` (`W` then includes the
sinks' mass), and `<output_basename>-sinks_history.dat` gets one row per sink per save (`it time id m x y z u v w`).
Restart files carry `<restart_basename>.sinks`. A sink that leaves a closed domain is fatal.

---

## [numerics]

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `scheme_space` | string | yes | `weno` (WENO flux splitting) or `weno-riemann` (WENO interpolation + Riemann flux + high-order correction). Other → fatal. | Spatial operator. MHD at very low β (≲ 1e-3, e.g. the Balsara–Spicer blast): `weno` needs the positivity limiter at every step and breaks point symmetry at the grid scale; use `weno-riemann` with `hlld` where symmetry or backend agreement matter ([#77](https://github.com/szaghi/adam/issues/77), [verification](./verification/runtime-amr/mhd-blast-symmetry)). |
| `reconstruction_variables` | string | yes | With `weno`: `characteristic` or `conservative`. With `weno-riemann`: `characteristic` or `primitive`. Other → fatal. | Variables reconstructed (or interpolated) at the faces. |
| `riemann_solver` | string | `weno-riemann` | `llf`, `hll`, `hllc`, `hlld`. Other → fatal. By model: `euler` allows `llf`, `hll`, `hllc`; `mhd-ideal` allows `llf`, `hll`, `hlld` (fatal otherwise). Not read with `weno`. | Face Riemann solver. |
| `flux_correction` | string | `weno-riemann` | `6th`, `4th`, `none` (2nd order). Other → fatal. | Face-flux correction `F^ = c1 F + c2 (f_i+f_i+1) + c3 (f_i-1+f_i+2)`. |
| `flux_correction_sensor` | string | `weno-riemann` | `weno` or `none`. Other → fatal. | `weno`: the correction is switched off at faces where `min_k w_k/d_k < 0.2` (fixed threshold `FLUX_CORRECTION_SENSOR_TAU`). `none`: always on. |
| `reflux` | logical | yes | logical | Berger–Colella reflux at AMR coarse-fine faces (and inter-realm seams). `.false.` is a diagnostic only: the run is then not conservative across 2:1 faces. |
| `reflux_policy` | string | no | `always`, `troubled` (default `always`; any other value is fatal) | With `reflux = .true.`, which seam face cells the reflux corrects ([#68](https://github.com/szaghi/adam/issues/68) E1). `always`: every one, conservative to round-off. `troubled`: only those flagged in some stage of the step, by the positivity limiter (a blended seam face, issue #50) or by the WENO-5 sensor (the minimum weight ratio $\min_k w_k/d_k$ of density, normal velocity or pressure across the face below 0.2, the weno-riemann flux-correction sensor); elsewhere the coarse flux stays, so smooth regions keep the accuracy of the high-order seam operators at the cost of conservation there (Sod c2f: mass drift $1.5\cdot10^{-7}$, shock unchanged). The log reports the flagged face cells per step. CPU and FNL (E1b: the FNL backend gathers the seam cells' auxiliary variables to the host and flags them there). The default stays `always` (owner decision after E2: exact conservation); `troubled` is the switch that lets `seam_restriction = lagrange6` with `seam_ghost_fill = quintic` keep sixth order across 2:1 seams with the reflux on (AO `euler-tq`), at the price of a bounded conservation drift (about $10^{-5}$ per period on the smooth vortex at N = 64, the reflux-off rate). |
| `reflux_cadence` | string | no | `step`, `stage` (default `step`; any other value is fatal; `stage` is fatal on multi-realm runs) | With `reflux = .true.`, when the coarse cells are corrected ([#81](https://github.com/szaghi/adam/issues/81) P5c). `step`: once, at the end of the step, from the step's fluxes (Berger–Colella, the AMReX cadence). `stage`: at every Runge–Kutta stage the coarse seam cells' residual is corrected from that stage's fluxes, so every stage advances them with the restricted fine fluxes and the scheme is the conservative composite one at every stage. Both conserve to round-off; `stage` keeps the time integration's order where `step` leaves an $O(\Delta t)$ inconsistency, e.g. the $E + W$ of self-gravity on AMR (first order with `step`, third with `stage`: [SV-10](./verification/self-gravity/amr-jeans)). Intra-realm 2:1 faces only; the exceptions of the step cadence (the equilibrium's mismatch, B under EGLM, `troubled`) apply per stage. |
| `positivity_limiter` | string | no (default `none`) | `none` or `cell`. Other → fatal. `cell` is fatal with `[mhd] divergence_control = glm`, a non-SSP `[runge_kutta] scheme`, immersed solids, multi-realm runs and any dissipative coefficient (issue #65, D-M4-5). | `cell`: the cell-based positivity limiter ([numerics](./numerics#positivity-limiter)): every face flux blended with the first-order Lax–Friedrichs backbone so that each stage keeps the density and the pressure positive; the limited faces are logged per stage. |

| `dissipative_order` | int | no (default 4) | `2` or `4`. Other → fatal; `4` needs `[grid] ngc >= 3` (fatal otherwise). | Order of the dissipative face fluxes (issue #65, D-M4-2): `4` is the conservative 4th-order flux with the Shu–Osher correction, `2` the compact central one ([numerics](./numerics#dissipative-fluxes-navier–stokes)). An ideal run ignores it. |

---

## [bc_x_min], [bc_x_max], [bc_y_min], [bc_y_max], [bc_z_min], [bc_z_max]

All six sections are required, even along null directions.

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `type` | string | yes (every face) | `extrapolation`, `inflow`, `wall-inviscid`, `wall-noslip`, `wall-isothermal`, `periodic`, `hydrostatic`. Other → fatal. The two faces of an axis must be **both** periodic or **neither** (fatal). | See the list below. |
| `r`, `u`, `v`, `w`, `p` | real | `type = inflow` | no range check | Prescribed primitive inflow state. |
| `bx`, `by`, `bz` | real | `type = inflow` **and** `mhd-ideal` | no range check | Prescribed inflow magnetic field. `psi` is set to 0. |
| `s1` .. `sN` | real | no (`type = inflow` with `passive_scalars = N`) | default 0, no range check | Specific values (mass fractions) of the passive scalars of the inflow. |
| `wall_u`, `wall_v`, `wall_w` | real | no (default 0) with `wall-noslip`, `wall-isothermal` | the component normal to the face must be 0 (fatal) | Velocity of the wall (tangential): a moving lid, a Couette wall. |
| `wall_temperature` | real | `type = wall-isothermal` | `> 0`, otherwise fatal | Wall temperature, $T = p/(\rho R)$ in code units. |

What each `type` does:

- `extrapolation`: zeroth-order extrapolation.
- `inflow`: prescribed state.
- `wall-inviscid`: slip wall (mirror state with the normal momentum negated; for MHD the normal B is negated too).
- `wall-noslip`: adiabatic no-slip wall (the velocity reflected about the wall velocity, the temperature mirrored).
- `wall-isothermal`: isothermal no-slip wall (the velocity as above, the ghost temperature $T_w^2/T$, the geometric mirror). See
  [Boundary conditions](./boundary-conditions#wall-noslip-and-wall-isothermal-no-slip-walls).
- `periodic`: library periodicity, so true periodic neighbours across blocks and ranks.
- `hydrostatic` ([#73](https://github.com/szaghi/adam/issues/73)): closes a stratified atmosphere under `[gravity]`.
  The ghost density and pressure continue the isothermal hydrostatic equilibrium from the first interior cell along
  the normal, $\rho_g = \rho_1 \exp(g_d s / T_1)$, $p_g = \rho_g T_1$ ($T_1 = p_1/\rho_1$, $s$ the signed distance
  from that cell to the ghost), so an isothermal atmosphere has exact ghosts. The velocity is mirrored with its normal
  component negated (no mass flux). The field and psi are mirrored **even**: field lines pass through the boundary, so
  a uniform field along or across gravity stays uniform (a reflecting conductor is `wall-inviscid`). Without gravity
  along the normal it is a slip wall that copies the thermodynamic state of the first interior cell.

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
| `isentropic-vortex` | **euler** (fatal otherwise) | 1 (free stream) | `x0`, `y0`, `radius` (`> 0`), `strength`. |
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
| `linear` | any | 1 | `gradient_x`, `gradient_y`, `gradient_z` (inverse length). Every conservative variable is `q_1 (1 + g . x)`, `q_1` the state of region 1: the verification field of the ghost probe (issue #65 P0), not a flow. |
| `isothermal-atmosphere` | any | 1 (the state at the origin; its velocity and field uniform) | Needs `[gravity]` (fatal otherwise). Optional, both atmospheres: a relative pressure pulse $p(1 + A e^{-|\mathbf x - \mathbf x_0|^2/w^2})$, `pulse_amplitude` $A$ (default 0), `pulse_x`, `pulse_y`, `pulse_z` ($\mathbf x_0$), `pulse_width` $w$ (`> 0` when $A \ne 0$, fatal otherwise; M7-G2b). $\rho = \rho_1 \exp(\mathbf g\cdot\mathbf x / T_1)$, $p = \rho T_1$, $T_1 = p_1/\rho_1$: the hydrostatic equilibrium of the constant gravity ([#73](https://github.com/szaghi/adam/issues/73)). |
| `polytropic-atmosphere` | any | 1 (as above) | `polytropic_gamma` $G > 1$ (fatal otherwise); needs `[gravity]`. $p = p_1(\rho/\rho_1)^G$, $\rho = \rho_1 s^{1/(G-1)}$, $s = 1 + \frac{G-1}{G}\frac{\rho_1}{p_1}\mathbf g\cdot\mathbf x$; fatal where $s \le 0$ (the top of the atmosphere inside the domain). |
| `rayleigh-taylor` | any | 1 (the light layer: $T_l = p_1/\rho_1$, $p_1$ the pressure at the interface; its field uniform) | Needs `[gravity]` along one axis (fatal otherwise). `rt_height` $h_0$ (along $-\mathbf g$, $h = -\mathbf g\cdot\mathbf x/|\mathbf g|$), `rt_width` $\delta$, `rt_density_ratio` $R$ (heavy over light at equal pressure; $R < 1$ is the stable control), `rt_wavelength` $\lambda$, `rt_amplitude` $A$, `rt_seed_width` $w$; every key required, $\delta, R, \lambda, w > 0$ (fatal otherwise). Two isothermal layers joined by $T = a + b\tanh((h - h_0)/\delta)$, $a = (T_l + T_h)/2$, $b = (T_h - T_l)/2$, $T_h = T_l/R$, hydrostatic in closed form, at rest; seeded with the solenoidal single mode $v_h = A\cos(k x_t) f(h)$, $u_t = -\frac{A}{k}\sin(k x_t) f'(h)$, $f = e^{-((h - h_0)/w)^2}$, $k = 2\pi/\lambda$, $x_t$ the first axis other than $\mathbf g$'s (M7-G3; its linear growth rate: `verification/gravity/rt_oracle.py`). |
| `polynomial` | any | 1 | `gradient_x`, `gradient_y`, `gradient_z`, `degree` (an integer $n \ge 1$). Every conservative variable is `q_1 (1 + sum_d g_d x_d^n)`: the ghost probe's field for the high-order seam operators ([#68](https://github.com/szaghi/adam/issues/68) D1), not a flow. Fatal under a non-dimensional reference (`g_d` has dimension $L^{-n}$). |
| `bonnor-ebert` | **euler**, `[eos] type = isothermal` | 1 (its velocity) | Needs `[gravity] self_gravity`. `be_density` $\rho_c$ (`> 0`, required), `be_x0`, `be_y0`, `be_z0` (default 0), `be_overdensity` $f$ (`> 0`, default 1). $\rho = f\rho_c e^{-\psi(r/r_0)}$, $r_0 = c_0/\sqrt{4\pi G\rho_c}$, $\psi$ the isothermal Lane–Emden function ([#100](https://github.com/szaghi/adam/issues/100) P4). |
| `boss-bodenheimer` | **euler**, `[eos]` isothermal or barotropic | 1 (the ambient medium) | Needs `[gravity] self_gravity`. `bb_radius` $R$, `bb_density` $\rho_0$ (`> 0`, required), `bb_omega` $\Omega$ (required), `bb_amplitude` $A$ (default 0.1), `bb_x0`, `bb_y0`, `bb_z0` (default 0). Inside $R$: $\rho = \rho_0(1 + A\cos 2\phi)$, solid-body rotation $\Omega\,\mathbf e_z \times (\mathbf x - \mathbf x_0)$ (P4). |

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
| `divb_flux` | logical | no (`.false.`) | MHD with GLM or EGLM ([#68](https://github.com/szaghi/adam/issues/68) A2): adds four columns to `-divb_history.dat` (`max_dflux`, `l1_dflux`, `seam_max_dflux`, `max_h_dflux_over_bmax`). D_flux is the divergence the cleaning acts on: that of the interface normal field $B_n = F_\psi/s$ read from the scheme's own $\psi$ face fluxes ($s = c_h^2$ GLM, $c_h$ EGLM), the coarse 2:1 faces taking the fine fluxes (the reflux), so its integral over a composite periodic domain is round-off. The last column is the scale-free $\max(h\,|D_\text{flux}|)/\max|\mathbf B|$. Each history row costs one residual evaluation of the committed state, on the idle RK stage buffers: fatal without GLM/EGLM or with a one-stage scheme. An inter-realm seam keeps each side's own flux; the step-0 row, written before the forest builds its register, has no reflux substitution and no seam column (as the centred one). |

---

## Forest manifest (multi-realm)

A multi-realm FLUME run reads a forest manifest, documented in the
[library input reference](/reference/library-input#forest-manifest-multi-realm); every `[realm.N]` points at a full
FLUME input of the form documented here. FLUME refuses runtime regridding (`[amr] frequency > 0`) in a multi-realm
forest.

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

## Complete example 2: ideal MHD with GLM, Orszag–Tang vortex with a 2:1 AMR box (verbatim)

There is **no committed INI** under `src/tests/flume/verification/mhd/orszag-tang/`. `check.sh` generates the input
with `make_orszag_tang.py` from `verification/vortex/vortex-n064.ini`. The file below is copied without changes from
the generated (git-ignored) `src/tests/flume/verification/mhd/orszag-tang/work-adam_flume_cpu-np2-amr/orszag-tang.ini`
(written 2026-09-29, `--cells 32 --refine-box 0.25 0.25 0.75 0.75`). It is written by Python `configparser`, which
explains the section order and the spacing.

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

`save_simulation_data` calls `save_xh5f(with_ghost=.true.)`, so **field files include the ghost cells**. In a forest
the inter-realm seam ghosts are refilled before the write, then the intra-realm ghosts and the boundary conditions
(the order of a Runge-Kutta stage), so the written ghosts are those the stencils read; the step-0 file, written while
each realm initialises and before the forest connects them, holds unfilled seam ghosts.

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

Here $|\mathbf B|^2$ is summed in plain $x,y,z$ order, not through `mhd_sum3`.

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

The operator is the same one used by the div(B) history kernels, `compute_divb_norms_dev`. Those kernels report $\max|\nabla\!\cdot\!\mathbf B|$,
$\sum|\nabla\!\cdot\!\mathbf B|\,\Delta V$ and a seam-band maximum, and cap $s$ at `HS_MAX = 3`. The CPU history also
requires $s\le n_{gc}$. The output path has no such cap or check.
