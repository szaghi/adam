---
title: adam_nasto_fnl_kernels
---

# adam_nasto_fnl_kernels

> ADAM, NASTO FNL application kernels.

**Source**: `src/app/nasto/fnl/adam_nasto_fnl_kernels.F90`

**Dependencies**

```mermaid
graph LR
  adam_nasto_fnl_kernels["adam_nasto_fnl_kernels"] --> adam_fnl_weno_kernels["adam_fnl_weno_kernels"]
  adam_nasto_fnl_kernels["adam_nasto_fnl_kernels"] --> adam_nasto_fnl_cns_kernels["adam_nasto_fnl_cns_kernels"]
  adam_nasto_fnl_kernels["adam_nasto_fnl_kernels"] --> fundal["fundal"]
  adam_nasto_fnl_kernels["adam_nasto_fnl_kernels"] --> penf["penf"]
```

## Contents

- [compute_fluxes_convective_dev](#compute-fluxes-convective-dev)
- [compute_fluxes_difference_dev](#compute-fluxes-difference-dev)
- [compute_fluxes_diffusive_dev](#compute-fluxes-diffusive-dev)
- [compute_umax_dev](#compute-umax-dev)
- [set_bc_q_gpu_dev](#set-bc-q-gpu-dev)
- [compute_fluxes_convective_ri_dev](#compute-fluxes-convective-ri-dev)
- [decompose_fluxes_convective_device](#decompose-fluxes-convective-device)

## Subroutines

### compute_fluxes_convective_dev

Compute convective fluxes along x direction.

```fortran
subroutine compute_fluxes_convective_dev(dir, blocks_number, ni, nj, nk, ngc, nv, weno_s, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, ror_number, ror_schemes_gpu, ror_threshold, ror_ivar_gpu, ror_stats_gpu, g, q_aux_gpu, fluxes_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `dir` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Direction, 1=X, 2=Y, 3=Z. |
| `blocks_number` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of blocks. |
| `ni` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in I direction. |
| `nj` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in J direction. |
| `nk` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in K direction. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Ghost cells number. |
| `nv` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of conservative varibales. |
| `weno_s` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Weno stencils number/dimension. |
| `weno_a_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Optimal weights. |
| `weno_p_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Polinomials coefficients. |
| `weno_d_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Smoothness indicators coefficients. |
| `weno_zeps` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Parameter to avoid division by zero. |
| `ror_number` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of ROR iterations. |
| `ror_schemes_gpu` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Scheme (S value) for each ROR step. |
| `ror_threshold` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | ROR threshold triggering. |
| `ror_ivar_gpu` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Index variables to check in ROR. |
| `ror_stats_gpu` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout |  | Scheme (S value) for each ROR step. |
| `g` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Specific heats ratio. |
| `q_aux_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Auxiliary variables. |
| `fluxes_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout |  | Fluxes. |

**Call graph**

```mermaid
flowchart TD
  compute_residuals["compute_residuals"] --> compute_fluxes_convective_dev["compute_fluxes_convective_dev"]
  compute_fluxes_convective_dev["compute_fluxes_convective_dev"] --> compute_fluxes_convective_ri_dev["compute_fluxes_convective_ri_dev"]
  style compute_fluxes_convective_dev fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### compute_fluxes_difference_dev

Compute fluxes difference.

```fortran
subroutine compute_fluxes_difference_dev(null_x, null_y, null_z, blocks_number, ni, nj, nk, ngc, nv, ib_eps, dx_gpu, dy_gpu, dz_gpu, flx_gpu, fly_gpu, flz_gpu, phi_gpu, dq_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `null_x` | logical | in |  | Nullified directions tags. |
| `null_y` | logical | in |  | Nullified directions tags. |
| `null_z` | logical | in |  | Nullified directions tags. |
| `blocks_number` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of blocks. |
| `ni` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in I direction. |
| `nj` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in J direction. |
| `nk` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in K direction. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Ghost cells number. |
| `nv` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of conservative varibales. |
| `ib_eps` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Tolerance IB delta ratio. |
| `dx_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Space steps. |
| `dy_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Space steps. |
| `dz_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Space steps. |
| `flx_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | X direction fluxes. |
| `fly_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Y direction fluxes. |
| `flz_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Z direction fluxes. |
| `phi_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | optional | IB distance function. |
| `dq_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout |  | Fluxes differences. |

**Call graph**

```mermaid
flowchart TD
  compute_residuals["compute_residuals"] --> compute_fluxes_difference_dev["compute_fluxes_difference_dev"]
  style compute_fluxes_difference_dev fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### compute_fluxes_diffusive_dev

Compute diffusive fluxes.

```fortran
subroutine compute_fluxes_diffusive_dev(null_x, null_y, null_z, blocks_number, ni, nj, nk, ngc, nv, mu, kd, q_aux_gpu, dx_gpu, dy_gpu, dz_gpu, flx_gpu, fly_gpu, flz_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `null_x` | logical | in |  | Nullified directions tags. |
| `null_y` | logical | in |  | Nullified directions tags. |
| `null_z` | logical | in |  | Nullified directions tags. |
| `blocks_number` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Blocks number. |
| `ni` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid dimensionns. |
| `nj` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid dimensionns. |
| `nk` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid dimensionns. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of ghost cells. |
| `nv` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of conservative variables. |
| `mu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Viscosity. |
| `kd` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Thermal diffusivity. |
| `q_aux_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Auxiliary varibales |
| `dx_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Space steps. |
| `dy_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Space steps. |
| `dz_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Space steps. |
| `flx_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout |  | Fluxes along x. |
| `fly_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout |  | Fluxes along y. |
| `flz_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout |  | Fluxes along z. |

**Call graph**

```mermaid
flowchart TD
  compute_residuals["compute_residuals"] --> compute_fluxes_diffusive_dev["compute_fluxes_diffusive_dev"]
  style compute_fluxes_diffusive_dev fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### compute_umax_dev

Compute maximum speed.

```fortran
subroutine compute_umax_dev(ni, nj, nk, ngc, blocks_number, mu, dx_gpu, dy_gpu, dz_gpu, q_aux_gpu, umax)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `ni` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in I direction. |
| `nj` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in J direction. |
| `nk` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in K direction. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Ghost cells number. |
| `blocks_number` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of blocks. |
| `mu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Dynamic viscosity. |
| `dx_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | X space step. |
| `dy_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Y space step. |
| `dz_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Z space step. |
| `q_aux_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Auxiliary varibales. |
| `umax` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | out |  | Maximum speed. |

**Call graph**

```mermaid
flowchart TD
  compute_dt["compute_dt"] --> compute_umax_dev["compute_umax_dev"]
  style compute_umax_dev fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### set_bc_q_gpu_dev

Set BC over q.

```fortran
subroutine set_bc_q_gpu_dev(BC_EXTRAPOLATION, BC_INFLOW, nv, ngc, cv, R, l_map_bc_gpu, fec_1_6_array_gpu, q_bc_vars_gpu, q_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `BC_EXTRAPOLATION` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Extrapolation BC parameter. |
| `BC_INFLOW` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Inflow BC parameter. |
| `nv` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of variables. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Ghost cells number. |
| `cv` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Constant volume specific heat. |
| `R` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Gas constant. |
| `l_map_bc_gpu` | integer(kind=[I8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Local map for BC ghost cells. |
| `fec_1_6_array_gpu` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Local map for BC ghost cells. |
| `q_bc_vars_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Boundary variables. |
| `q_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout |  | Conservative variables. |

**Call graph**

```mermaid
flowchart TD
  set_boundary_conditions["set_boundary_conditions"] --> set_bc_q_gpu_dev["set_bc_q_gpu_dev"]
  style set_bc_q_gpu_dev fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### compute_fluxes_convective_ri_dev

Compute convective fluxes at right interface of b,i,j,k.

```fortran
subroutine compute_fluxes_convective_ri_dev(dir, si, sir, b, i, j, k, ngc, nv, weno_s, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, ror_number, ror_schemes_gpu, ror_threshold, ror_ivar_gpu, ror_stats_gpu, g, q_aux_gpu, fluxes_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `dir` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Direction, 1=X, 2=Y, 3=Z. |
| `si` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Stencil increment. |
| `sir` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Stencil increment, real cast. |
| `b` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Counter. |
| `i` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Counter. |
| `j` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Counter. |
| `k` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Counter. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Ghost cells number. |
| `nv` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of conservative varibales. |
| `weno_s` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Weno stencils number/dimension. |
| `weno_a_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Optimal weights. |
| `weno_p_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Polinomials coefficients. |
| `weno_d_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Smoothness indicators coefficients. |
| `weno_zeps` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Parameter to avoid division by zero. |
| `ror_number` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of ROR iterations. |
| `ror_schemes_gpu` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Scheme (S value) for each ROR step. |
| `ror_threshold` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | ROR threshold triggering. |
| `ror_ivar_gpu` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Index variables to check in ROR. |
| `ror_stats_gpu` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout |  | Scheme (S value) for each ROR step. |
| `g` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Specific heats ratio. |
| `q_aux_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Auxiliary variables. |
| `fluxes_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout |  | Fluxes. |

**Call graph**

```mermaid
flowchart TD
  compute_fluxes_convective_dev["compute_fluxes_convective_dev"] --> compute_fluxes_convective_ri_dev["compute_fluxes_convective_ri_dev"]
  compute_fluxes_convective_dev["compute_fluxes_convective_dev"] --> compute_fluxes_convective_ri_dev["compute_fluxes_convective_ri_dev"]
  compute_fluxes_convective_ri_dev["compute_fluxes_convective_ri_dev"] --> compute_eigenvectors_dev["compute_eigenvectors_dev"]
  compute_fluxes_convective_ri_dev["compute_fluxes_convective_ri_dev"] --> decompose_fluxes_convective_device["decompose_fluxes_convective_device"]
  compute_fluxes_convective_ri_dev["compute_fluxes_convective_ri_dev"] --> weno_reconstruct_upwind_dev["weno_reconstruct_upwind_dev"]
  style compute_fluxes_convective_ri_dev fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### decompose_fluxes_convective_device

Decompose convective fluxes.
 Flux vector splitting by local-Lax-Friedrics (Rusanov) with projection in pseudo-characteristics psace.

```fortran
subroutine decompose_fluxes_convective_device(si, sir, el, weno_s, b, i, j, k, ngc, nv, g, q_aux_gpu, fmpc)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `si` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Stencil increment. |
| `sir` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Stencil increment, real cast. |
| `el` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Left eigeinvectors. |
| `weno_s` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Weno stencils number/dimension. |
| `b` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Counter. |
| `i` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Counter. |
| `j` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Counter. |
| `k` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Counter. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Ghost cells number. |
| `nv` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of conservative varibales. |
| `g` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Specific heats ratio. |
| `q_aux_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Auxiliary variables. |
| `fmpc` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout |  | Fluxes -+ decomposition in characteristics space. |

**Call graph**

```mermaid
flowchart TD
  compute_fluxes_convective_device["compute_fluxes_convective_device"] --> decompose_fluxes_convective_device["decompose_fluxes_convective_device"]
  compute_fluxes_convective_device["compute_fluxes_convective_device"] --> decompose_fluxes_convective_device["decompose_fluxes_convective_device"]
  compute_fluxes_convective_ri_dev["compute_fluxes_convective_ri_dev"] --> decompose_fluxes_convective_device["decompose_fluxes_convective_device"]
  decompose_fluxes_convective_device["decompose_fluxes_convective_device"] --> compute_conservative_fluxes_dev["compute_conservative_fluxes_dev"]
  decompose_fluxes_convective_device["decompose_fluxes_convective_device"] --> compute_conservatives_dev["compute_conservatives_dev"]
  decompose_fluxes_convective_device["decompose_fluxes_convective_device"] --> compute_max_eigenvalues_dev["compute_max_eigenvalues_dev"]
  style decompose_fluxes_convective_device fill:#3e63dd,stroke:#99b,stroke-width:2px
```
