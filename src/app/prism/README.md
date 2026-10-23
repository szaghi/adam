# PRISM

**PRISM** — **P**lasma **R**esearch us**I**ng **S**imulation **M**ethods — is an electromagnetic and plasma simulation application built on the ADAM framework. It solves the Maxwell equations for computational electromagnetics and, optionally, couples them to charged-particle dynamics via a Particle-in-Cell (PIC) method.

## Physical Models

PRISM supports two physical models, selected in the INI file under `[physics] physical_model`:

| Model | INI key | Description |
|-------|---------|-------------|
| Electromagnetic | `electromagnetic` | Solves the Maxwell equations for **D** and **B** with a prescribed current source **J** |
| Particle-in-Cell | `PIC` | Couples Maxwell equations to individual charged-particle trajectories; current **J** and charge density ρ are deposited from particle positions |

## Governing Equations

PRISM solves the first-order hyperbolic form of the Maxwell equations:

$$\frac{\partial \mathbf{q}}{\partial t} + \nabla \cdot \mathbf{F}(\mathbf{q}) = \mathbf{S}$$

The conservative state vector is $\mathbf{q} = (D_x, D_y, D_z, B_x, B_y, B_z, J_x, J_y, J_z)^\top$ and the flux components are:

$$F_x = \begin{pmatrix} 0 \\ B_z/\mu \\ -B_y/\mu \\ 0 \\ -D_z/\varepsilon \\ D_y/\varepsilon \\ 0 \\ 0 \\ 0 \end{pmatrix} \quad F_y = \begin{pmatrix} -B_z/\mu \\ 0 \\ B_x/\mu \\ D_z/\varepsilon \\ 0 \\ -D_x/\varepsilon \\ 0 \\ 0 \\ 0 \end{pmatrix} \quad F_z = \begin{pmatrix} B_y/\mu \\ -B_x/\mu \\ 0 \\ -D_y/\varepsilon \\ D_x/\varepsilon \\ 0 \\ 0 \\ 0 \end{pmatrix}$$

The wave speeds are $\pm c_0 = \pm 1/\sqrt{\varepsilon_0 \mu_0}$; the current components $J_{x,y,z}$ have zero flux (source-only variables).

## Source Layout

```
src/app/prism/
├── common/                                      # Shared across all backends
│   ├── adam_prism_parameters.F90                # Physical constants and eigenvector matrices
│   ├── adam_prism_physics_object.F90            # Physical model, variable counts, eigenvectors
│   ├── adam_prism_numerics_object.F90           # Temporal and spatial scheme selection
│   ├── adam_prism_common_object.F90             # prism_common_object base type
│   ├── adam_prism_bc_object.F90                 # Boundary conditions handler
│   ├── adam_prism_ic_object.F90                 # Initial conditions handler
│   ├── adam_prism_rk_bc_object.F90              # RK sub-stepping for boundary conditions
│   ├── adam_prism_io_object.F90                 # I/O handler
│   ├── adam_prism_time_object.F90               # Time integration handler
│   ├── adam_prism_coil_object.F90               # Current-carrying coil sources
│   ├── adam_prism_fWLayer_object.F90            # Far-wave layer (absorbing BC layer)
│   ├── adam_prism_external_fields_object.F90    # Externally prescribed fields
│   ├── adam_prism_pic_object.F90                # Particle-in-Cell configuration
│   ├── adam_prism_particle_injection_object.F90 # Particle injection handler
│   ├── adam_prism_leapfrog_pic_object.F90       # Leapfrog PIC time integrator
│   ├── adam_prism_rk_pic_object.F90             # Runge-Kutta PIC time integrator
│   ├── adam_prism_riemann_library.F90           # Maxwell Riemann solvers (LLF, HLL) and flux routines
│   └── adam_prism_common_library.F90            # Barrel re-export of all common modules
├── cpu/                                         # CPU backend
│   ├── adam_prism_cpu.F90                       # Entry point (program)
│   └── adam_prism_cpu_object.F90                # prism_cpu_object implementation
└── fnl/                                         # OpenACC GPU backend
    ├── adam_prism_fnl.F90                       # Entry point (program)
    ├── adam_prism_fnl_object.F90                # prism_fnl_object implementation
    ├── adam_prism_fnl_kernels.F90               # OpenACC GPU kernels
    ├── adam_prism_fnl_external_fields_kernels.F90 # OpenACC external-fields kernels
    ├── adam_prism_fnl_coil_object.F90           # GPU coil handler
    ├── adam_prism_fnl_fWLayer_object.F90        # GPU far-wave layer handler
    └── adam_prism_fnl_library.F90               # Barrel re-export of all FNL modules
```

## Field Variables

Variable counts are determined at runtime by the physical model and divergence correction settings:

| Group | Variables | Count | Notes |
|-------|-----------|-------|-------|
| Conservative | $D_x, D_y, D_z, B_x, B_y, B_z$ | `nv_c = 6` | Always present |
| Source | $J_x, J_y, J_z$ | `nv_s = 3` | Always present; zero flux |
| Divergence cleaning | $\varphi$ (D) and/or $\psi$ (B) | `nv_cl = 0, 1, 2` | Optional hyperbolic cleaning scalars |
| PIC charge density | $\rho$ | `nv_pic = 0, 1` | Present only in PIC model |

Total: `nv = nv_c + nv_s + nv_cl + nv_pic`

| Configuration | nv |
|--------------|-----|
| EM, no cleaning | 9 |
| EM + D or B cleaning | 10 |
| EM + D and B cleaning | 11 |
| PIC, no cleaning | 10 |
| PIC + D or B cleaning | 11 |
| PIC + D and B cleaning | 12 |

**Maximum supported**: `NV_MAX = 11`.

### Saved fields

The XH5F checkpoints hold, by default (`[IO] save_fields = all`), the whole state vector (its width depends on the
model above) and the current density of every coil, plus the arrays switched on by `save_residual_fields`,
`save_curl_fields` and `save_divergence_fields`. `save_fields = em` writes only $D_x, D_y, D_z, B_x, B_y, B_z$ (the
fields as stored, no conversion to $E$) for production runs; the three `save_*_fields` flags still add their arrays, so
set them to `.false.` for an output with the electromagnetic fields alone. Any other value is fatal.

## Numerical Methods

### Temporal Schemes

| Scheme | INI key | Description |
|--------|---------|-------------|
| Runge-Kutta | `RUNGE_KUTTA` | Explicit RK, scheme from `[runge_kutta] scheme`. SSP 11/22/33/54 and low-storage 1/2/3 on both backends; Yoshida on CPU only (FNL stops at init) |
| Leapfrog | `LEAPFROG` | **Not functional**: stops at init on both backends |
| Blanes-Moan | `BLANES_MOAN` | **Not functional**: stops at init on both backends ([issue #55](https://github.com/szaghi/adam/issues/55)) |
| Commutator-Free Magnus | `COMMUTATOR_FREE_MAGNUS` | **Not functional**: stops at init on both backends (issue #55) |

The three non-RK field integrators exist as library objects, but `realm_object%initialize` never initializes them, so
PRISM refuses them instead of running on unallocated coefficients (CPU) or empty device stubs (FNL). PIC's own particle
leapfrog (`[PIC] scheme_time`) is a separate path.

### Spatial Schemes

| Scheme | INI key | Description |
|--------|---------|-------------|
| WENO | `WENO` | High-order weighted essentially non-oscillatory reconstruction |
| FD centered | `FD_CENTERED` | Finite-difference centered stencil |
| FV centered | `FV_CENTERED` | Finite-volume centered stencil |

WENO reconstruction operates on either conservative variables (`CONSERVATIVE`) or characteristic variables (`CHARACTERISTICS`), selected via `[numerics] reconstruction_variables`.

### Riemann Solvers (`adam_prism_riemann_library`)

| Solver | Description |
|--------|-------------|
| LLF (Rusanov) | Local Lax-Friedrichs: $F = \tfrac{1}{2}(F_L + F_R - c_0 (q_R - q_L))$ |
| HLL | Harten–Lax–van Leer two-wave solver |

## Divergence Control

Maxwell's equations require $\nabla \cdot \mathbf{D} = \rho$ and $\nabla \cdot \mathbf{B} = 0$. `[numerics] divergence_correction` selects the strategy, and `[numerics] constrained_transport = D | B | DB | NO` the fields it acts on:

| Strategy | INI value | Description |
|----------|-----------|-------------|
| None | `no` | No divergence control; `constrained_transport` is ignored (a warning says so) |
| Hyperbolic | `hyperbolic` | GLM cleaning: one scalar per field named by `constrained_transport` (φ for D, ψ for B), carrying divergence errors away at speed $\chi c_0$ (`[physics] chi`, required, > 0) and damped at the rate `[physics] alpha` (1/s, optional, default 0). **Requires `constrained_transport = D`, `B` or `DB`**: without it there is no scalar to carry, and the run stops at initialization |
| Poisson | `poisson` | Not implemented: the run stops at initialization |

With hyperbolic cleaning the state grows by one slot per scalar, placed after B: `nv = 10` for `D` or `B`, `nv = 11` for `DB` (φ in 7, ψ in 8), J follows. The regression case `rmf-hyperbolic` anchors the `DB` layout ([issue #61](https://github.com/szaghi/adam/issues/61)).

## Source Terms: Coils

The `prism_coil_object` models current-carrying coils as volumetric current density sources $\mathbf{J}$:

| Coil type | Description |
|-----------|-------------|
| `rectangular` | Rectangular loop coil with user-defined dimensions and normal direction |
| `circular` | Circular loop coil |

| Current type | Description |
|--------------|-------------|
| `DC_current` | Constant (DC) current |
| `AC_current` | Sinusoidal (AC) current with amplitude, frequency, and initial phase |

Coil normals are specified as `+x`, `-x`, `+y`, `-y`, `+z`, or `-z`. Coil $\mathbf{J}$ contributions are computed each time step and registered for optional output as auxiliary fields (`j_vec_1/2/3`, `f_Gauss`).

## Particle-in-Cell (PIC)

When `physical_model = PIC`, PRISM tracks charged macro-particles alongside the Maxwell field. Each particle carries 8 state variables stored in `q_pic(8, particle_number)`:

| Index | Variable | Description |
|-------|----------|-------------|
| 1 | x | Position x |
| 2 | y | Position y |
| 3 | z | Position z |
| 4 | vx | Velocity x |
| 5 | vy | Velocity y |
| 6 | vz | Velocity z |
| 7 | charge | Particle charge |
| 8 | mass | Particle mass |

Field values at particle locations are stored in `pic_fields(8, particle_number)`. Particle time integration uses either Leapfrog (`LEAPFROG`) or Runge-Kutta (`RUNGE_KUTTA`), configured independently from the field scheme.

Two problem types are supported:

| Type | INI key | Description |
|------|---------|-------------|
| Plasma | `plasma` | Ensemble of charged particles |
| Single particle | `single_particle` | Single test-particle trajectory |

## `prism_common_object` Overview

The base type aggregates all framework and PRISM sub-objects:

**ADAM framework objects**: `mpih_object`, `adam_object`, `field_object`, `grid_object`, `amr_object`, `ib_object`, `slices_object`, `leapfrog_object`, `blanesmoan_object`, `cfm_object`, `rk_object`, `weno_object`, `flail_object`

**PRISM sub-objects**: `prism_physics_object`, `prism_numerics_object`, `prism_io_object`, `prism_bc_object`, `prism_ic_object`, `prism_rk_bc_object`, `prism_time_object`, `prism_coil_object`, `prism_fWLayer_object`, `prism_external_fields_object`, `prism_pic_object`, `prism_particle_injection_object`, `prism_leapfrog_pic_object`, `prism_rk_pic_object`

**Field arrays** (all shape `(nv, ni, nj, nk, nb)`):

| Array | Description |
|-------|-------------|
| `q` | Conservative field: D, B, J (+ optional cleaning scalars and ρ) |
| `dq` | Residuals RHS |
| `curl` | Curl of field components (optional output) |
| `divergence` | Divergence of field components (optional output) |

**PIC arrays**: `q_pic(8, particle_number)`, `pic_fields(8, particle_number)`

**Energy diagnostics**: `energy_D(:)`, `energy_B(:)` (time histories), `rms_energy_error_D/B`

## Multi-realm support

PRISM realms participate in multi-realm forests as `prism_cpu_object` / `prism_fnl_object` arrays. Both backends override the `_forest`-suffixed TBP family (`begin_stage_forest`, `end_stage_forest`, `close_step_forest`, `fill_seam_from_peer_forest`, `apply_reflux_to_stage_forest`, `coupling_descriptor_forest`, ...) inherited from `realm_object` via `prism_common_object`.

- **Worked manifest example**: `src/tests/prism/regression/rmf-2realm/input.ini` and the sibling cases (`rmf-2realm-asymK`, `rmf-2realm-stagesync`).
- **Conceptual overview**: [Forest (multi-realm)](/guide/forest) — the project documentation site covers the manifest schema, the α/β seam-coupling cadence trade-offs, and the orchestrator's phase cycle.
- **Library contract reference**: `src/lib/common/README.md` → "Forest orchestration" — the per-realm and per-forest TBP surface, manifest data flow, and load-bearing invariants.

## Backends

| Backend | Entry point | Key type | Accelerator |
|---------|-------------|----------|-------------|
| CPU | `adam_prism_cpu.F90` | `prism_cpu_object` | MPI + OpenMP |
| FNL | `adam_prism_fnl.F90` | `prism_fnl_object` | MPI + OpenACC (NVIDIA GPU) |

The FNL backend adds GPU companion objects (`field_fnl_object`, `rk_fnl_object`, `weno_fnl_object`, etc.) and device-side pointer arrays (`q_gpu`, `dq_gpu`, `flx_f_gpu`, `curl_gpu`, …). The CPU backend implements all differential operators (curl, divergence, gradient, Laplacian, derivatives 1–4) via procedure pointers, dispatched at initialisation based on `scheme_space`.

## Building

```bash
# CPU backend (GNU compiler)
fobis build --mode prism-cpu-gnu

# CPU debug build
fobis build --mode prism-cpu-gnu-debug

# FNL (OpenACC) backend
fobis build --mode prism-fnl-nvf --varset local_nvf
```

Executables are written to `exe/adam_prism_cpu` and `exe/adam_prism_fnl`.

## Running

```bash
mpirun -np <N> exe/adam_prism_cpu [input.ini]
mpirun -np <N> exe/adam_prism_fnl [input.ini]
```

If no argument is given, `input.ini` in the working directory is used.

## License

PRISM is part of the ADAM framework, released under the [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.html) (LGPLv3).

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.
