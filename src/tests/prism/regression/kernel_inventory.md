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

20 kernels — the densest single file. Where the directive is inside an
explicit `contains`ed kernel subroutine the class is unambiguously A.
Where the directive sits in a top-level method body without a kernel-
subroutine wrapper, the class depends on whether the surrounding
`associate` is reaching through `self%`.

| Line | Caller method (host) | Wrapper subroutine | Class | Notes |
|------|---------------------|--------------------|-------|-------|
|  539 | `compute_coils_current` | `nullify_j_vec_vars_kernel` | A | All dummies. |
|  565 | `compute_coils_current` | `apply_j_vec_kernel` | A | All dummies. |
|  617 | `set_boundary_conditions` | `set_boundary_conditions_kernel` | A | All dummies. |
|  759 | `compute_curl` | `compute_curl_kernel` | A | All dummies. |
|  922 | `compute_divergence` | `compute_divergence_kernel` | A | All dummies. |
| 1096 | `compute_residuals_fd_centered` | (inline, branch 1) | **B/C** | Kernel uses `dxyz_gpu, q_gpu, dq_gpu` from enclosing `associate` reaching into `self%`. See "Hot file note" below. |
| 1168 | `compute_residuals_fd_centered` | (inline, branch 2) | **B/C** | Same as 1096. |
| 1236 | `compute_residuals_fd_centered` | (inline, branch 3) | **B/C** | Same as 1096. |
| 1317 | `compute_residuals_fd_centered` | (inline, branch 4) | **B/C** | Same as 1096. |
| 1663 | `compute_min_dxyz` | (inline) | **B/C** | `dxyz_gpu` from `associate(... =>self%...)`. Reduction kernel. |
| 1729 | `nullify_q` | (inline) | **B/C** | `q_gpu` from `associate(... =>self%q_gpu)`. |
| 1754 | `nullify_dq` | (inline) | **B/C** | `dq_gpu` from `associate(... =>self%dq_gpu)`. |
| 1783 | `compute_poynting_flux` (face -x) | (inline) | **B/C** | `q_gpu, dxyz_gpu` from `associate`. |
| 1801 | `compute_poynting_flux` (face +x) | (inline) | **B/C** | Same. |
| 1818 | `compute_poynting_flux` (face -y) | (inline) | **B/C** | Same. |
| 1835 | `compute_poynting_flux` (face +y) | (inline) | **B/C** | Same. |
| 1852 | `compute_poynting_flux` (face -z) | (inline) | **B/C** | Same. |
| 1869 | `compute_poynting_flux` (face +z) | (inline) | **B/C** | Same. |
| 1951 | `compute_max_divergence` (in `compute_max_divergence_dev_kernel`) | A | All dummies. |
| 2026 | `impose_ct_correction` | (inline) | **B/C** | `q_gpu, buffer` from `associate` reaching into `self%q_gpu, self%divergence_gpu`. |

#### Hot file note

`adam_prism_fnl_object.F90` mixes two patterns deliberately:

1. **Already-class-A:** kernels with `contains`ed wrapper subroutines
   taking explicit dummies (lines 539, 565, 617, 759, 922, 1951). These
   match the plan's target end-state. Six of twenty.
2. **Class-B/C disguised by `associate`:** kernels that sit inline in a
   host method body. The `associate(ni=>self%ni, q_gpu=>self%q_gpu, ...)`
   on entry to the method gives the kernel a flat set of local names, so
   the kernel *body* never writes `self%`. But the names are renamed
   references into derived-type components — exactly what R3 of the plan
   forbids in device regions ("`associate` blocks are not a substitute
   for dummy arguments in device code").

   Fourteen of twenty fall in this bucket. The classifier is B-vs-C in
   the plan's spelling, but the practical class is "associate-renamed
   chain through self." That is what Step 4 must eliminate.

No class-D found. The `q_gpu`, `dq_gpu`, `divergence_gpu`, etc. members
of `prism_fnl_object` are `pointer :: ... => null()` to allocatable GPU
storage — that *is* the D-class pattern at the type level, but the
kernels go through `associate` rather than direct chain walks, so the
device code never spells out `self%X`. The dispatch isn't D in practice.
That said: the `associate` and the pointer-to-DT both depend on the
program-local `prism` instance staying alive and bound for the duration
of the kernel — which is the same invariant the rest of the codebase
already depends on, so not a new risk surface.

## Summary

| Class | Count | Where |
|-------|-------|-------|
| A — clean | 35 | All of `lib/fnl/`, plus 6 in `app/prism/fnl/adam_prism_fnl_object.F90` (the kernels with `contains`ed wrappers), plus all 3 in the `app/prism/fnl/` aux files. |
| B/C — associate-renamed chain through self | 14 | All inside `adam_prism_fnl_object.F90`. |
| D — pointer-DT chain | 0 | none. |

So Step 4 work is concentrated in **one file**: `adam_prism_fnl_object.F90`.
The 14 inline kernels (lines 1096–1317, 1663–1869, 2026) each need to be
extracted into a `contains`ed subroutine taking explicit dummies, the
same shape the other 6 kernels in the same file already use.

Pre-existing examples in this file (e.g. `set_boundary_conditions_kernel`
at 603, `nullify_j_vec_vars_kernel` at 532, `apply_j_vec_kernel` at 555)
are the template to follow.

## Step 5 (CPU) scope

The plan calls for the same audit on the CPU backend. `prism-cpu-gnu`
has no OpenACC kernels by definition — the equivalent inner loops are
just `do`-loops without device directives. The CPU pass is therefore a
*style* refactor (kernel = inner-loop block → contained subroutine) for
testability per R4, not a portability fix. It applies to the same loop
bodies that are class-B/C above (the host-side loop nests inside
`compute_residuals`, `compute_poynting_flux`, `impose_ct_correction`,
`compute_min_dxyz`, `nullify_q`, `nullify_dq` — in `adam_prism_cpu_object.F90`).

A full CPU audit will be added here as a follow-up at Step 5 entry.
