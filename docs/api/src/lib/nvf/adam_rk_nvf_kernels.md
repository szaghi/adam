---
title: adam_rk_nvf_kernels
---

# adam_rk_nvf_kernels

> ADAM, RK NVF kernels (NVF backend of [rk_nvf_object](/api/src/lib/nvf/adam_rk_nvf_object#rk-nvf-object)).

**Source**: `src/lib/nvf/adam_rk_nvf_kernels.F90`

**Dependencies**

```mermaid
graph LR
  adam_rk_nvf_kernels["adam_rk_nvf_kernels"] --> penf["penf"]
```

## Contents

- [rk_assign_stage_cuf](#rk-assign-stage-cuf)
- [rk_compute_stage_cuf](#rk-compute-stage-cuf)
- [rk_compute_stage_ls_cuf](#rk-compute-stage-ls-cuf)
- [rk_initialize_stages_cuf](#rk-initialize-stages-cuf)
- [rk_update_q_cuf](#rk-update-q-cuf)

## Subroutines

### rk_assign_stage_cuf

Assign q to RK stage.

```fortran
subroutine rk_assign_stage_cuf(ni, nj, nk, ngc, nv, blocks_number, s, phi_gpu, q_gpu, q_rk_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `ni` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Grid cells number in I direction. |
| `nj` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Grid cells number in J direction. |
| `nk` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Grid cells number in K direction. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Ghost cells number. |
| `nv` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Number of conservative varibales. |
| `blocks_number` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Number of blocks. |
| `s` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Stage index. |
| `phi_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device, optional | IB distance. |
| `q_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device | Conservative field. |
| `q_rk_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout | device | Conservative field stage. |

**Call graph**

```mermaid
flowchart TD
  assign_stage["assign_stage"] --> rk_assign_stage_cuf["rk_assign_stage_cuf"]
  style rk_assign_stage_cuf fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### rk_compute_stage_cuf

Sum RK stages up to s.

```fortran
subroutine rk_compute_stage_cuf(ni, nj, nk, ngc, nv, blocks_number, s, dt, alph, phi_gpu, q_rk_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `ni` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Grid cells number in I direction. |
| `nj` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Grid cells number in J direction. |
| `nk` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Grid cells number in K direction. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Ghost cells number. |
| `nv` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Number of conservative varibales. |
| `blocks_number` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Number of blocks. |
| `s` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Current stage. |
| `dt` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Time step. |
| `alph` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device | RK alpha coefficients. |
| `phi_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device, optional | IB distance. |
| `q_rk_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout | device | Conservative field stages. |

**Call graph**

```mermaid
flowchart TD
  compute_stage["compute_stage"] --> rk_compute_stage_cuf["rk_compute_stage_cuf"]
  style rk_compute_stage_cuf fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### rk_compute_stage_ls_cuf

Compute RK stage, low storage scheme using only q(n) and q.

```fortran
subroutine rk_compute_stage_ls_cuf(ni, nj, nk, ngc, nv, blocks_number, dt, ark, brk, crk, phi_gpu, q_n_gpu, dq_gpu, q_rk_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `ni` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Grid cells number in I direction. |
| `nj` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Grid cells number in J direction. |
| `nk` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Grid cells number in K direction. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Ghost cells number. |
| `nv` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Number of conservative varibales. |
| `blocks_number` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Number of blocks. |
| `dt` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Time step. |
| `ark` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | RK coefficients. |
| `brk` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | RK coefficients. |
| `crk` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | RK coefficients. |
| `phi_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device, optional | IB distance. |
| `q_n_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device | RK stage 0, Q^n. |
| `dq_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device | Residuals. |
| `q_rk_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout | device | Conservative field stage. |

**Call graph**

```mermaid
flowchart TD
  compute_stage_ls["compute_stage_ls"] --> rk_compute_stage_ls_cuf["rk_compute_stage_ls_cuf"]
  style rk_compute_stage_ls_cuf fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### rk_initialize_stages_cuf

Initialize RK stages.

```fortran
subroutine rk_initialize_stages_cuf(ni, nj, nk, ngc, nv, blocks_number, q_gpu, q_rk_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `ni` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in I direction. |
| `nj` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in J direction. |
| `nk` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in K direction. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Ghost cells number. |
| `nv` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of conservative varibales. |
| `blocks_number` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of blocks. |
| `q_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device | Conservative field. |
| `q_rk_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout | device | RK stage. |

**Call graph**

```mermaid
flowchart TD
  initialize_stages["initialize_stages"] --> rk_initialize_stages_cuf["rk_initialize_stages_cuf"]
  style rk_initialize_stages_cuf fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### rk_update_q_cuf

Update RK q.

```fortran
subroutine rk_update_q_cuf(ni, nj, nk, ngc, nv, blocks_number, nrk, dt, beta, phi_gpu, q_rk_gpu, q_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `ni` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Grid cells number in I direction. |
| `nj` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Grid cells number in J direction. |
| `nk` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Grid cells number in K direction. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Ghost cells number. |
| `nv` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Number of conservative varibales. |
| `blocks_number` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Number of blocks. |
| `nrk` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Number of RK stages. |
| `dt` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Time step. |
| `beta` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device | RK betaa coefficients. |
| `phi_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device, optional | IB distance. |
| `q_rk_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device | Conservative field stage. |
| `q_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout | device | Conservative field. |

**Call graph**

```mermaid
flowchart TD
  update_q["update_q"] --> rk_update_q_cuf["rk_update_q_cuf"]
  style rk_update_q_cuf fill:#3e63dd,stroke:#99b,stroke-width:2px
```
