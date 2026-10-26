# Backends

PRISM builds for two backends from one common code base:

| Backend | Program | Type | Parallelism | Build modes |
|---|---|---|---|---|
| CPU | `adam_prism_cpu` | `prism_cpu_object` | MPI + OpenMP | `prism-cpu-gnu`, `prism-cpu-gnu-debug`; `prism-cpu-amd` (amdflang) |
| FNL | `adam_prism_fnl` | `prism_fnl_object` | MPI + OpenACC (NVIDIA), or OpenMP offload (AMD) | `prism-fnl-nvf`, `-debug`, `-spacehpc` (with a `--varset`); `prism-fnl-omp-amd`, `-debug` |

The full list is on the [build-modes page](/reference/build-modes). Both types extend `prism_common_object`
(`src/app/prism/common`), which extends the library's `realm_object`; the FNL device helpers are components of the FNL
type ([library on FNL](/concepts/library-fnl)).

## What each backend implements

| Feature | CPU | FNL |
|---|---|---|
| `fd_centered`, `fv_centered` | yes | yes |
| `weno` | yes | refused |
| SSP and low-storage Runge–Kutta | yes | yes |
| Yoshida | yes (EM) | refused |
| leapfrog, Blanes–Moan, commutator-free Magnus | refused | refused |
| GLM cleaning (`hyperbolic`) | yes | yes |
| PML | `CLASSIC`, `CLASSIC_DIRECT`, `BERMUDEZ`, `CFS` with `fd_centered`, SSP, no correction | the same kinds with `fd_centered`, SSP, no hyperbolic correction, no constrained transport |
| far-wave layer | yes | yes |
| PIC, standard deposits | yes | yes |
| PIC, charge-conserving current | one block, one rank | yes (on a virtual block) |
| coherent matched gather without the conserving current | one block, one rank (warns) | refused |
| boundary kinds, periodic faces | all | all |
| t0 elliptic solves | FFTW | cuFFT |
| static AMR, 2:1 faces, FV reflux | yes | yes |
| multi-realm forests | yes | yes |
| probes, `[grms]`, field-at-centre histories | yes | yes |

Every refusal stops the run at initialization with a message naming the combination; on FNL a few come after the device
allocation (`weno`, the conserving-current refusals, the t0 refined-grid refusal); the full list
is in the [input hazards](./input#refused-combinations).

## Lifecycle

Both programs hand the realm (one, or the manifest's realms) to the library's forest, `forest%simulate` or
`forest%simulate_from_manifest`, which drives the realm contract ([forest](/concepts/forest)):

1. **`initialize_forest`.** The common initialization (below), the backend's dispatch and checks, then the initial
   conditions: `[initial_conditions] amr_iterations` passes of initial condition and AMR markers, the coils, the t0
   elliptic solves, the far-wave correction; or the restart files. The initial divergences are printed and the first
   outputs written.
2. **Every step.** The forest reduces the time step over the realms (`compute_local_dt_forest`), then advances each
   realm (`advance_one_step_forest`: one call of the dispatched integrator) or, with 2:1 faces or several realms, runs
   the staged path stage by stage (`open_step`, `begin_stage`, `end_stage`, `close_step`, the reflux).
3. **`post_step_forest`.** Field and restart files on their cadences, the energy, coil-current, RMS, probe and
   divergence histories (the inter-realm seams refilled first), the PIC charge diagnostic.
4. **`is_done_forest`** (`it_max` or `time_max`), then **`finalize_forest`**.

### Common initialization

`prism_common_object%initialize`, called by both backends:

1. `[IO]` and `[IO] save_fields`.
2. `numerics`, then `physics` (which uses it): the state layout; `nv > NV_MAX` is fatal.
3. `realm_object%initialize`: grid, tree, maps, field, AMR markers, IB, slices, Runge–Kutta, WENO, `[fdv]` (lengths in
   $L_0$ for the non-dimensional models), on the memory budget below.
4. The face kinds, `grid%set_bc_type`, then `tree%rebuild_neighborhood` and the ghost and boundary maps (periodic faces
   need the kinds before the maps).
5. PIC only: `pic`, `particle_injection`.
6. `time`, `ic`, `coil`, `grms`, `magnetic_field_at_center_domain`, `probes`; PIC: the particle integrator.
7. The ghost-width checks (`ngc >= fdv_order/2`, and `>= [weno] S` with WENO), the field arrays, the uniform
   refinement, the variable names of the field files (they must number `nv`).

The FNL backend first initializes `mpih_fnl` (device selection, once per process), then runs the PIC, PML and
time-scheme checks before any device allocation, then builds its device objects: `field_fnl`, `ib_fnl`, `rk_fnl`,
`weno_fnl`, the PML, the device arrays, `coil_fnl`, the PIC objects, FLAIL on cuFFT, the external-field kernels.

### Memory budget

The number of blocks a rank can hold is computed from a memory budget and the number of block-sized fields (the state,
the residual, the fluxes, the Runge–Kutta stages). On the CPU the budget is the rank's share of the host memory
(`mpih%memory_avail`); on FNL it is the device memory divided by the realms and by the ranks that share the device:
$M = M_{\rm device} / (N_{\rm realms}\,N_{\rm sharers})$. An allocation that does not fit stops the run with the name of
the array (every `dev_alloc` checks its error code).

## Device data (FNL)

The device arrays are **transposed**, `(nb, ni, nj, nk, nv)` against the host `(nv, ni, nj, nk, nb)`: the kernels
collapse the loops with the block index `b` innermost, and adjacent threads then read adjacent blocks, a coalesced
access along the stride-1 index. The state
crosses between host and device only for the initial conditions, the restart, and the outputs, through FUNDAL's
buffered transposed copies, whose staging buffer has the destination's shape
([#31](https://github.com/szaghi/adam/issues/31)).

| Array | Content |
|---|---|
| `q_gpu`, `dq_gpu` | state and residual |
| `flxyz_c_gpu`, `flx_f_gpu`, `fly_f_gpu`, `flz_f_gpu` | cell and face fluxes (FV path) |
| `curl_gpu`, `divergence_gpu` | diagnostics |
| `coil_fnl`, `pml_fnl`, `pic_fnl`, `conserving_fnl`, `probes_fnl` | device state of the coils, the PML, the particles, the conserving current, the probes |

The ghost exchange is GPU-direct MPI on the device buffers. On the WSL2 development box this needs
`UCX_RNDV_THRESH=inf` and other WSL-only settings, which must not be copied to a cluster
([library on FNL](/concepts/library-fnl#gpu-direct-ghost-exchange-and-the-wsl2-rendezvous-abort-issue-12)); WSL runs
are for correctness only, never for timings.

## Spatial-operator dispatch

Both backends bind their operators once, at initialization, by `scheme_space`: the curl, divergence, gradient,
Laplacian and first, second and fourth derivatives to the FD or FV variants of `adam_fdv_operators_library` (their device
twins on FNL), and `compute_residuals` to the FD, FV or WENO residual; the Maxwell flux function by the divergence
correction (`compute_convective_fluxes_maxwell`, `_div_d`, `_div_b`, `_div_d_b`, and their `_adim` twins,
`adam_prism_riemann_library`); the
integrator by the model and `[runge_kutta] scheme`.

## Source layout

| Directory | Content |
|---|---|
| `src/app/prism/common` | the physics, numerics, I/O, time, boundary and initial-condition objects; coils, external fields, far-wave layer, PML, PIC (configuration, boundaries, injection, integrators), probes and the other diagnostics; `adam_prism_common_object`, the realm shared by the backends |
| `src/app/prism/cpu` | `adam_prism_cpu_object` (operators, integrators, boundary kernels, PIC deposits on the host) and the program |
| `src/app/prism/fnl` | `adam_prism_fnl_object` and its device objects (coils, far-wave layer, external-field kernels, PML, PIC, conserving current, particle integrators, probes) and the program |

The module tables are in the READMEs of the three directories.
