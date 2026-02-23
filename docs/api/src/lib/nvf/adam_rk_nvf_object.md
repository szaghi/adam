---
title: adam_rk_nvf_object
---

# adam_rk_nvf_object

> ADAM, RK class NVF (NVF backend of [rk_object](/api/src/lib/common/adam_rk_object#rk-object)).

**Source**: `src/lib/nvf/adam_rk_nvf_object.F90`

**Dependencies**

```mermaid
graph LR
  adam_rk_nvf_object["adam_rk_nvf_object"] --> adam_memory_nvf_library["adam_memory_nvf_library"]
  adam_rk_nvf_object["adam_rk_nvf_object"] --> adam_mpih_nvf_object["adam_mpih_nvf_object"]
  adam_rk_nvf_object["adam_rk_nvf_object"] --> adam_rk_nvf_kernels["adam_rk_nvf_kernels"]
  adam_rk_nvf_object["adam_rk_nvf_object"] --> adam_rk_object["adam_rk_object"]
  adam_rk_nvf_object["adam_rk_nvf_object"] --> penf["penf"]
```

## Contents

- [rk_nvf_object](#rk-nvf-object)
- [assign_stage](#assign-stage)
- [compute_stage](#compute-stage)
- [compute_stage_ls](#compute-stage-ls)
- [initialize](#initialize)
- [initialize_stages](#initialize-stages)
- [update_q](#update-q)

## Derived Types

### rk_nvf_object

RK NVF class definition.

#### Components

| Name | Type | Attributes | Description |
|------|------|------------|-------------|
| `rk` | type([rk_object](/api/src/lib/common/adam_rk_object#rk-object)) | pointer | RK common handler. |
| `mpih` | type([mpih_nvf_object](/api/src/lib/nvf/adam_mpih_nvf_object#mpih-nvf-object)) |  | MPI handler. |
| `alph_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | allocatable, device | RK alpha coefficients. |
| `beta_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | allocatable, device | RK beta coefficients. |
| `gamm_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | allocatable, device | RK gamma coefficients. |
| `q_rk_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | allocatable, device | Field cell centered variables, RK stages. |

#### Type-Bound Procedures

| Name | Attributes | Description |
|------|------------|-------------|
| `assign_stage` | pass(self) | Assign q to RK stage. |
| `compute_stage` | pass(self) | Compute RK stage. |
| `compute_stage_ls` | pass(self) | Compute RK stage, low storage scheme. |
| `initialize` | pass(self) | Initialize class. |
| `initialize_stages` | pass(self) | Initialize RK stages. |
| `update_q` | pass(self) | Update RK q. |

## Subroutines

### assign_stage

Assign q to RK stage.

```fortran
subroutine assign_stage(self, s, q_gpu, phi_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([rk_nvf_object](/api/src/lib/nvf/adam_rk_nvf_object#rk-nvf-object)) | inout |  | RK object. |
| `s` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Current stage number. |
| `q_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device | Conservative variables. |
| `phi_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device, optional | IB distance. |

**Call graph**

```mermaid
flowchart TD
  integrate["integrate"] --> assign_stage["assign_stage"]
  integrate["integrate"] --> assign_stage["assign_stage"]
  integrate["integrate"] --> assign_stage["assign_stage"]
  integrate["integrate"] --> assign_stage["assign_stage"]
  integrate["integrate"] --> assign_stage["assign_stage"]
  integrate_rk_ssp["integrate_rk_ssp"] --> assign_stage["assign_stage"]
  assign_stage["assign_stage"] --> check_cuda_error["check_cuda_error"]
  assign_stage["assign_stage"] --> rk_assign_stage_cuf["rk_assign_stage_cuf"]
  style assign_stage fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### compute_stage

Compute RK stage.

```fortran
subroutine compute_stage(self, s, dt, phi_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([rk_nvf_object](/api/src/lib/nvf/adam_rk_nvf_object#rk-nvf-object)) | inout |  | RK object. |
| `s` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Current stage number. |
| `dt` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Current time step. |
| `phi_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device, optional | IB distance. |

**Call graph**

```mermaid
flowchart TD
  integrate["integrate"] --> compute_stage["compute_stage"]
  integrate["integrate"] --> compute_stage["compute_stage"]
  integrate["integrate"] --> compute_stage["compute_stage"]
  integrate["integrate"] --> compute_stage["compute_stage"]
  integrate["integrate"] --> compute_stage["compute_stage"]
  integrate_rk_ssp["integrate_rk_ssp"] --> compute_stage["compute_stage"]
  compute_stage["compute_stage"] --> check_cuda_error["check_cuda_error"]
  compute_stage["compute_stage"] --> rk_compute_stage_cuf["rk_compute_stage_cuf"]
  style compute_stage fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### compute_stage_ls

Compute RK stage, low storage scheme.
 The first (only) stage is assumed to be the previous time step q solution.

```fortran
subroutine compute_stage_ls(self, s, dt, phi_gpu, dq_gpu, q_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([rk_nvf_object](/api/src/lib/nvf/adam_rk_nvf_object#rk-nvf-object)) | in |  | RK object. |
| `s` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Current RK stage. |
| `dt` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Current time step. |
| `phi_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device, optional | IB distance. |
| `dq_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device | Conservative variables residuals. |
| `q_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout | device | Conservative variables stage. |

**Call graph**

```mermaid
flowchart TD
  integrate["integrate"] --> compute_stage_ls["compute_stage_ls"]
  integrate["integrate"] --> compute_stage_ls["compute_stage_ls"]
  integrate["integrate"] --> compute_stage_ls["compute_stage_ls"]
  integrate["integrate"] --> compute_stage_ls["compute_stage_ls"]
  integrate["integrate"] --> compute_stage_ls["compute_stage_ls"]
  compute_stage_ls["compute_stage_ls"] --> check_cuda_error["check_cuda_error"]
  compute_stage_ls["compute_stage_ls"] --> rk_compute_stage_ls_cuf["rk_compute_stage_ls_cuf"]
  style compute_stage_ls fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### initialize

Initialize class.

```fortran
subroutine initialize(self, rk, nb, ngc, ni, nj, nk, nv)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([rk_nvf_object](/api/src/lib/nvf/adam_rk_nvf_object#rk-nvf-object)) | inout |  | RK NVF object. |
| `rk` | type([rk_object](/api/src/lib/common/adam_rk_object#rk-object)) | in | target | RK object. |
| `nb` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Total blocks number for MPI. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of ghost cells. |
| `ni` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of cells in i direction. |
| `nj` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of cells in j direction. |
| `nk` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of cells in k direction. |
| `nv` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of conservative variables. |

**Call graph**

```mermaid
flowchart TD
  add_node["add_node"] --> initialize["initialize"]
  analize["analize"] --> initialize["initialize"]
  build_connectivity["build_connectivity"] --> initialize["initialize"]
  distribute_facets["distribute_facets"] --> initialize["initialize"]
  distribute_facets_tree["distribute_facets_tree"] --> initialize["initialize"]
  export_vtk_file["export_vtk_file"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> initialize["initialize"]
  initialize_common["initialize_common"] --> initialize["initialize"]
  initialize_common["initialize_common"] --> initialize["initialize"]
  initialize_common["initialize_common"] --> initialize["initialize"]
  initialize_common["initialize_common"] --> initialize["initialize"]
  load_from_file["load_from_file"] --> initialize["initialize"]
  load_from_file["load_from_file"] --> initialize["initialize"]
  open_file["open_file"] --> initialize["initialize"]
  open_file["open_file"] --> initialize["initialize"]
  open_file["open_file"] --> initialize["initialize"]
  resize["resize"] --> initialize["initialize"]
  save_into_file["save_into_file"] --> initialize["initialize"]
  save_vtk["save_vtk"] --> initialize["initialize"]
  simulate["simulate"] --> initialize["initialize"]
  simulate["simulate"] --> initialize["initialize"]
  simulate["simulate"] --> initialize["initialize"]
  simulate["simulate"] --> initialize["initialize"]
  simulate["simulate"] --> initialize["initialize"]
  simulate["simulate"] --> initialize["initialize"]
  simulate["simulate"] --> initialize["initialize"]
  test_stress["test_stress"] --> initialize["initialize"]
  initialize["initialize"] --> alloc_var_gpu["alloc_var_gpu"]
  initialize["initialize"] --> assign_allocatable_gpu["assign_allocatable_gpu"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> print_message["print_message"]
  style initialize fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### initialize_stages

Initialize RK stages.

```fortran
subroutine initialize_stages(self, q_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([rk_nvf_object](/api/src/lib/nvf/adam_rk_nvf_object#rk-nvf-object)) | inout |  | RK object. |
| `q_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device | Conservative variables. |

**Call graph**

```mermaid
flowchart TD
  integrate["integrate"] --> initialize_stages["initialize_stages"]
  integrate["integrate"] --> initialize_stages["initialize_stages"]
  integrate["integrate"] --> initialize_stages["initialize_stages"]
  integrate["integrate"] --> initialize_stages["initialize_stages"]
  integrate["integrate"] --> initialize_stages["initialize_stages"]
  integrate_rk_ssp["integrate_rk_ssp"] --> initialize_stages["initialize_stages"]
  initialize_stages["initialize_stages"] --> check_cuda_error["check_cuda_error"]
  initialize_stages["initialize_stages"] --> rk_initialize_stages_cuf["rk_initialize_stages_cuf"]
  style initialize_stages fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### update_q

Update RK q.

```fortran
subroutine update_q(self, dt, phi_gpu, q_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([rk_nvf_object](/api/src/lib/nvf/adam_rk_nvf_object#rk-nvf-object)) | in |  | RK object. |
| `dt` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Current time step. |
| `phi_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device, optional | IB distance. |
| `q_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout | device | Conservative variables. |

**Call graph**

```mermaid
flowchart TD
  integrate["integrate"] --> update_q["update_q"]
  integrate["integrate"] --> update_q["update_q"]
  integrate["integrate"] --> update_q["update_q"]
  integrate["integrate"] --> update_q["update_q"]
  integrate["integrate"] --> update_q["update_q"]
  integrate_rk_ssp["integrate_rk_ssp"] --> update_q["update_q"]
  update_q["update_q"] --> check_cuda_error["check_cuda_error"]
  update_q["update_q"] --> rk_update_q_cuf["rk_update_q_cuf"]
  style update_q fill:#3e63dd,stroke:#99b,stroke-width:2px
```
