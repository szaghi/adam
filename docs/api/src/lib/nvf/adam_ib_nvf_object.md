---
title: adam_ib_nvf_object
---

# adam_ib_nvf_object

> ADAM, IB class NVF (NVF backend of [ib_object](/api/src/lib/common/adam_ib_object#ib-object)).

**Source**: `src/lib/nvf/adam_ib_nvf_object.F90`

**Dependencies**

```mermaid
graph LR
  adam_ib_nvf_object["adam_ib_nvf_object"] --> adam_field_nvf_object["adam_field_nvf_object"]
  adam_ib_nvf_object["adam_ib_nvf_object"] --> adam_ib_nvf_kernels["adam_ib_nvf_kernels"]
  adam_ib_nvf_object["adam_ib_nvf_object"] --> adam_ib_object["adam_ib_object"]
  adam_ib_nvf_object["adam_ib_nvf_object"] --> adam_memory_nvf_library["adam_memory_nvf_library"]
  adam_ib_nvf_object["adam_ib_nvf_object"] --> adam_mpih_nvf_object["adam_mpih_nvf_object"]
  adam_ib_nvf_object["adam_ib_nvf_object"] --> penf["penf"]
```

## Contents

- [ib_nvf_object](#ib-nvf-object)
- [evolve_eikonal](#evolve-eikonal)
- [initialize](#initialize)
- [invert_eikonal](#invert-eikonal)

## Derived Types

### ib_nvf_object

IB NVF class definition.

#### Components

| Name | Type | Attributes | Description |
|------|------|------------|-------------|
| `ib` | type([ib_object](/api/src/lib/common/adam_ib_object#ib-object)) | pointer | IB common handler. |
| `mpih` | type([mpih_nvf_object](/api/src/lib/nvf/adam_mpih_nvf_object#mpih-nvf-object)) |  | MPI handler. |
| `field_gpu` | type([field_nvf_object](/api/src/lib/nvf/adam_field_nvf_object#field-nvf-object)) | pointer | Field NVF handler. |
| `q_bcs_vars_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | allocatable, device | Variables array for immersed boundary on GPU. |
| `phi_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | allocatable, device | Distance function on GPU. |
| `blocks_number` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Actual blocks number. |
| `nb` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Total blocks number for MPI. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Number of ghost cells. |
| `ni` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Number of cells in i direction. |
| `nj` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Number of cells in j direction. |
| `nk` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Number of cells in k direction. |
| `nv` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Number of conservative variables. |

#### Type-Bound Procedures

| Name | Attributes | Description |
|------|------------|-------------|
| `evolve_eikonal` | pass(self) | Evolve eikonal equation. |
| `initialize` | pass(self) | Initialize class. |
| `invert_eikonal` | pass(self) | Invert momentum eikonal equation. |

## Subroutines

### evolve_eikonal

Evolve eikonal equation.

```fortran
subroutine evolve_eikonal(self, dq_gpu, q_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([ib_nvf_object](/api/src/lib/nvf/adam_ib_nvf_object#ib-nvf-object)) | in |  | IB. |
| `dq_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout | device | State variables variations. |
| `q_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout | device | Conservative variables. |

**Call graph**

```mermaid
flowchart TD
  integrate_eikonal["integrate_eikonal"] --> evolve_eikonal["evolve_eikonal"]
  integrate_eikonal["integrate_eikonal"] --> evolve_eikonal["evolve_eikonal"]
  integrate_eikonal["integrate_eikonal"] --> evolve_eikonal["evolve_eikonal"]
  integrate_eikonal["integrate_eikonal"] --> evolve_eikonal["evolve_eikonal"]
  integrate_eikonal["integrate_eikonal"] --> evolve_eikonal["evolve_eikonal"]
  integrate_eikonal["integrate_eikonal"] --> evolve_eikonal["evolve_eikonal"]
  evolve_eikonal["evolve_eikonal"] --> check_cuda_error["check_cuda_error"]
  evolve_eikonal["evolve_eikonal"] --> compute_eikonal_dq_phi_cuf["compute_eikonal_dq_phi_cuf"]
  evolve_eikonal["evolve_eikonal"] --> evolve_eikonal_q_phi_cuf["evolve_eikonal_q_phi_cuf"]
  style evolve_eikonal fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### initialize

Initialize class.

```fortran
subroutine initialize(self, ib, field_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([ib_nvf_object](/api/src/lib/nvf/adam_ib_nvf_object#ib-nvf-object)) | inout |  | IB NVF object. |
| `ib` | type([ib_object](/api/src/lib/common/adam_ib_object#ib-object)) | in | target | IB object. |
| `field_gpu` | type([field_nvf_object](/api/src/lib/nvf/adam_field_nvf_object#field-nvf-object)) | in | target | The field. |

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
  initialize["initialize"] --> associate_adam_data["associate_adam_data"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> print_message["print_message"]
  style initialize fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### invert_eikonal

Invert momentum eikonal equation.

```fortran
subroutine invert_eikonal(self, q_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([ib_nvf_object](/api/src/lib/nvf/adam_ib_nvf_object#ib-nvf-object)) | in |  | IB. |
| `q_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout | device | Conservative variables. |

**Call graph**

```mermaid
flowchart TD
  integrate_eikonal["integrate_eikonal"] --> invert_eikonal["invert_eikonal"]
  integrate_eikonal["integrate_eikonal"] --> invert_eikonal["invert_eikonal"]
  integrate_eikonal["integrate_eikonal"] --> invert_eikonal["invert_eikonal"]
  integrate_eikonal["integrate_eikonal"] --> invert_eikonal["invert_eikonal"]
  integrate_eikonal["integrate_eikonal"] --> invert_eikonal["invert_eikonal"]
  integrate_eikonal["integrate_eikonal"] --> invert_eikonal["invert_eikonal"]
  invert_eikonal["invert_eikonal"] --> check_cuda_error["check_cuda_error"]
  invert_eikonal["invert_eikonal"] --> invert_eikonal_q_phi_cuf["invert_eikonal_q_phi_cuf"]
  style invert_eikonal fill:#3e63dd,stroke:#99b,stroke-width:2px
```
