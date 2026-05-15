# PRISM kernel inventory

**Step 3 of the forest-of-trees migration (issue #10).**

This file is the read-only audit deliverable: every kernel reachable from
`prism-cpu-gnu` and `prism-fnl-nvf`, classified by chain depth. The
classes drive the Step 4 / Step 5 refactor work (which kernels still need
to be made `self`-free).

## Class definitions

Per the plan (issue #10, Step 3.2):

- **A — clean.** Kernel reads only local dummies and scalars passed in by
  the caller. No `self%`, no use-imported module-scope state read inside
  the device region. **No work needed.**
- **B — module-singleton.** Kernel body reads through a use-imported
  global singleton (e.g. `field%q(...)`, `weno%coefficients(...)`).
  Mechanical conversion: lift the read to the caller, pass the array as
  an explicit dummy.
- **C — chain-through-self.** Kernel body reads `self%X(...)` or
  `self%X%Y(...)`. Same mechanical lift as B but the call site lives on a
  type-bound procedure.
- **D — chain-through-pointer-DT.** Kernel reads through a `pointer ::`
  derived-type component. **Highest risk;** miscompiles silently on
  nvfortran with `-fast` in some versions. Any D found here is a
  pre-existing bug.

## Scope

Only OpenACC kernels (`!$acc parallel`, `!$acc kernels`) are inventoried.
PRISM has no NVF backend and no GMP backend, so CUDA Fortran
(`attributes(global)`, `!$cuf kernel`) and OpenMP target (`!$omp
target`) regions are out of scope here — they exist in `src/lib/nvf/`,
`src/lib/gmp/`, and `src/app/nasto/{nvf,gmp}/`, all NASTO-only.

Counts at audit time:

| Region | Total | PRISM-reachable | Out of scope |
|--------|-------|-----------------|--------------|
| `!$acc parallel`  | 73 | 49 | 24 (NASTO FNL, READMEs) |
| `!$cuf kernel`    | 44 | 0  | 44 (NVF, NASTO-only) |
| `attributes(global)` | 7 | 0  | 7  (NVF, NASTO-only) |
| `!$omp target`    | 86 | 0  | 86 (GMP, NASTO-only) |

So the Step-4/5 inventory below covers 49 OpenACC regions across six
PRISM-FNL-reachable files. NVF/GMP kernels will be inventoried separately
when those backends come back into scope (out of scope per issue #10).

## Per-kernel classification

### `src/lib/fnl/adam_fnl_field_kernels.F90`

| Line | Subroutine | Class | Notes |
|------|------------|-------|-------|
|  39 | `compute_gradient_max_dev` | A | `q_gpu` dummy, `gradient` reduction. |
|  69 | `compute_norm_dev` | A | `dq_gpu` dummy, `norm_gpu` reduction. |
| 100 | `pack_send_ghost_buffer_dev` | A | `comm_map_send_ghost_cell_gpu`, `send_buffer_ghost_gpu`, `q_gpu` dummies. |
| 139 | `unpack_recv_ghost_buffer_dev` | A | `comm_map_recv_ghost_cell_gpu`, `recv_buffer_ghost_gpu`, `q_gpu` dummies. |
| 173 | `update_local_ghost_dev` | A | `l_map_ghost_cell_gpu`, `q_gpu` dummies. |

### `src/lib/fnl/adam_fnl_field_object.F90`

| Line | Subroutine | Class | Notes |
|------|------------|-------|-------|
| 151 | `set_q_to_q_t_dev` (in `field_fnl_object` method) | A | `q_gpu`, `q_t_gpu` dummies. Methods on field_fnl_object that wrap a kernel call have already adopted the dummy-arg pattern internally. |

### `src/lib/fnl/adam_fnl_rk_kernels.F90`

| Line | Subroutine | Class | Notes |
|------|------------|-------|-------|
|  39 | `rk_assign_stage_dev` (phi branch) | A | All dummies; phi optional. |
|  55 | `rk_assign_stage_dev` (no-phi branch) | A | All dummies. |
|  90 | `rk_compute_stage_dev` (phi branch) | A | All dummies; phi optional. |
| 111 | `rk_compute_stage_dev` (no-phi branch) | A | All dummies. |
| 151 | `rk_initialize_stages_dev` (phi branch) | A | All dummies. |
| 167 | `rk_initialize_stages_dev` (no-phi branch) | A | All dummies. |
| 195 | `rk_compute_stage_ls_dev` | A | All dummies; LS = low-storage RK. |
| 234 | `rk_update_q_dev` (phi branch) | A | All dummies. |
| 257 | `rk_update_q_dev` (no-phi branch) | A | All dummies. |

### `src/lib/fnl/adam_fnl_ib_kernels.F90`

| Line | Subroutine | Class | Notes |
|------|------------|-------|-------|
|  40 | `compute_eikonal_dq_phi_dev` | A | `dx_gpu/dy_gpu/dz_gpu/phi_gpu/q_gpu/dq_gpu` all dummies. |
| 104 | `init_phi_dev` | A | `phi_gpu` dummy. |
| 143 | `compute_phi_distance_dev` | A | `x_cell_gpu/y_cell_gpu/z_cell_gpu/phi_gpu` dummies. |
| 172 | `apply_eikonal_dev` | A | `phi_gpu/dq_gpu/q_gpu` dummies. |
| 208 | `nullify_solid_q_dev` (variant 1) | A | `phi_gpu/q_gpu` dummies. |
| 224 | `nullify_solid_q_dev` (variant 2) | A | `phi_gpu/q_gpu` dummies. |
| 275 | `compute_dphi_dev` (variant 1) | A | `phi_gpu/dphi_gpu` dummies. |
| 316 | `compute_dphi_dev` (variant 2) | A | `phi_gpu/dphi_gpu` dummies. |
| 347 | `classify_cell_scheme_dev` (variant 1) | A | `phi_gpu/cell_scheme_gpu` dummies. |
| 365 | `classify_cell_scheme_dev` (variant 2) | A | `phi_gpu/cell_scheme_gpu` dummies. |
| 383 | `classify_cell_scheme_dev` (variant 3) | A | `phi_gpu/cell_scheme_gpu` dummies. |

### `src/app/prism/fnl/adam_prism_fnl_external_fields_kernels.F90`

| Line | Subroutine | Class | Notes |
|------|------------|-------|-------|
| 135 | `apply_external_E_kernel` | A | All arrays as dummies including `e_ext_gpu, q_gpu`. |
| 223 | `apply_external_B_kernel` | A | All arrays as dummies including `b_ext_gpu, q_gpu`. |

### `src/app/prism/fnl/adam_prism_fnl_fWLayer_object.F90`

| Line | Subroutine | Class | Notes |
|------|------------|-------|-------|
| 117 | `apply_fwlayer_kernel` | A | `q_gpu`, `mask_gpu` dummies. |

### `src/app/prism/fnl/adam_prism_fnl_object.F90`

20 kernel directives, distributed across method bodies and contained
`*_dev_kernel` subroutines. Each directive line is mapped below to the
**actual** enclosing subroutine (verified by reading the file structure,
not by trusting `grep` line proximity).

| Line | Enclosing subroutine | Class | Notes |
|------|---------------------|-------|-------|
|  539 | `nullify_j_vec_vars_kernel` (contained in `compute_coils_current`) | A | All dummies. |
|  565 | `apply_j_vec_kernel` (contained in `compute_coils_current`) | A | All dummies. |
|  617 | `set_boundary_conditions_kernel` (contained in `set_boundary_conditions`) | A | All dummies. |
|  759 | `compute_curl_fd_dev_kernel` (contained in `compute_curl_fd_dev`) | A | All dummies. |
|  922 | `compute_divergence_fd_dev_kernel` (contained in `compute_divergence_fd_dev`) | A | All dummies. |
| 1096 | `compute_residuals_fd_centered_dev_kernel` branch 1 | A | All dummies; branch guarded by host `if (numerics%...)`. |
| 1168 | `compute_residuals_fd_centered_dev_kernel` branch 2 | A | Same. |
| 1236 | `compute_residuals_fd_centered_dev_kernel` branch 3 | A | Same. |
| 1317 | `compute_residuals_fd_centered_dev_kernel` branch 4 | A | Same. |
| 1663 | `compute_dxyz_min_kernel` (contained in `compute_dt`) | A | `dxyz_gpu` dummy; reduction kernel. |
| 1729 | `compute_e_dev_kernel` (contained in `compute_energy`) | A | All dummies; reduction. |
| 1754 | `compute_coil_power_dev_kernel` (contained in `compute_energy`) | A | All dummies; reduction. |
| 1783 | `compute_poynting_flux_dev_kernel` face -x | A | All dummies; six faces, six directives, one wrapper. |
| 1801 | `compute_poynting_flux_dev_kernel` face +x | A | Same. |
| 1818 | `compute_poynting_flux_dev_kernel` face -y | A | Same. |
| 1835 | `compute_poynting_flux_dev_kernel` face +y | A | Same. |
| 1852 | `compute_poynting_flux_dev_kernel` face -z | A | Same. |
| 1869 | `compute_poynting_flux_dev_kernel` face +z | A | Same. |
| 1951 | `compute_max_divergence_dev_kernel` (contained in `compute_max_divergence`) | A | All dummies. |
| 2026 | `impose_ct_correction` method body (no kernel wrapper) | **B/C** | `q_gpu, buffer` renamed via `associate` from `self%q_gpu, self%divergence_gpu`. The only inline kernel in this file. |

#### Hot file note

`adam_prism_fnl_object.F90` is **already almost fully refactored** to the
dummy-argument pattern. 19 of 20 kernel directives sit inside `*_dev_kernel`
contained subroutines that take explicit dummies — exactly the Step 4
target shape.

Exactly **one** kernel still sits inline in a method body: the
constrained-transport-correction update at line 2026, inside
`impose_ct_correction`. The surrounding `associate` renames
`self%q_gpu` and `self%divergence_gpu` to local `q_gpu` and `buffer`,
so the kernel body never spells `self%` — but per R3 of the plan,
`associate`-renamed reaches through `self%` ARE the pattern we are
eliminating. Step 4's PRISM-FNL work is to extract this one kernel
into an `impose_ct_correction_kernel` contained subroutine matching
the established template.

No class-D found. The `q_gpu`, `dq_gpu`, `divergence_gpu`, etc. members
of `prism_fnl_object` are `pointer :: ... => null()` to allocatable GPU
storage — that *is* the D-class pattern at the type level, but the
class-A wrappers pass these as plain dummies and never dereference the
type chain inside a device region. The dispatch isn't D in practice.

#### Audit-error note

An earlier version of this file (commit 57e17d5c) classified 14 of these
20 kernels as B/C. That was wrong: the lines I called "inline class B/C"
are in fact inside `*_dev_kernel` contained subroutines. The error came
from trusting `grep`-line proximity to subroutine boundaries instead of
verifying each kernel against its enclosing scope. The corrected table
above reflects what is actually in the source.

## Summary

| Class | Count | Where |
|-------|-------|-------|
| A — clean | 48 | All 26 in `lib/fnl/`, all 3 in `app/prism/fnl/` aux files, 19 of 20 in `adam_prism_fnl_object.F90`. |
| B/C — associate-renamed chain through self | 1 | `adam_prism_fnl_object.F90:2026`, kernel inside `impose_ct_correction`. |
| D — pointer-DT chain | 0 | none. |

So Step 4 work is **one kernel**: extract `impose_ct_correction`'s inline
kernel (line 2026) into an `impose_ct_correction_kernel` contained
subroutine taking explicit dummies (`ni, nj, nk, ngc, blocks_number,
ivar, q_gpu, buffer`).

Pre-existing examples in this file (e.g. `set_boundary_conditions_kernel`
at 603, `nullify_j_vec_vars_kernel` at 532, `apply_j_vec_kernel` at 555,
and the larger `compute_residuals_fd_centered_dev_kernel` at 1058) are
the template to follow. The PRISM-FNL backend is otherwise already at
the Step 4 end-state.

## Step 5 (CPU) scope — skipped

The plan calls for the same audit on the CPU backend. `prism-cpu-gnu`
has no OpenACC kernels by definition — the equivalent inner loops are
just `do`-loops without device directives. The CPU pass is therefore a
*style* refactor (inner-loop block → contained subroutine) for
testability per R4, not a portability fix.

**A CPU audit was performed and Step 5 was deliberately skipped.**
Findings:

- ~10 host methods in `adam_prism_cpu_object.F90` carry inline loop
  nests that read `self%` — directly or through `associate` aliasing.
  Candidates would have been `impose_ct_correction`,
  `compute_residuals_fd_centered` (4 branches), `compute_residuals_fv_centered`
  (5 branches), `compute_e`/`compute_coil_power`/`compute_Poynting_flux`
  (contained inside `compute_energy`; they read `self%q` via host
  association even though they ARE contained subroutines),
  `mark_by_j_vec_total_variation`, `update_q_BC`,
  `impose_div_coil_correction`.

- The plan explicitly classifies the CPU pass as "strictly speaking
  unnecessary for compiler portability (GNU and Intel don't have the
  chain bug)" — its only payoff is R4 testability and pattern
  consistency with FNL.

- No subsequent step (6, 7) depends on Step 5. Step 6 (`forest_object`
  cutover) modifies ownership and addressing, not kernel shape.

- The CPU side's pattern is internally consistent today: every method
  uses `associate` + inline loop. Refactoring some-but-not-all would
  introduce a *new* inconsistency. Refactoring all would be ~10 commits
  with no concrete payoff and real loop-ordering-changes-FP-noise risk.

Decision (2026-05-15): skip Step 5 entirely. Migration proceeds Step 4
→ Step 6. If R4 testability becomes the limiting factor later (e.g. for
a unit-test campaign on individual kernels), Step 5 can be revisited as
a self-contained refactor sweep with no migration dependencies blocking
it.
