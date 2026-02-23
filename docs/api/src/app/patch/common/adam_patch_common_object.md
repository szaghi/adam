---
title: adam_patch_common_object
---

# adam_patch_common_object

> ADAM, PATCH Poisson solver with Adpative mesh Refinement for HPC computing class definition, CPU backend.

**Source**: `src/app/patch/common/adam_patch_common_object.F90`

**Dependencies**

```mermaid
graph LR
  adam_patch_common_object["adam_patch_common_object"] --> adam_adam_object["adam_adam_object"]
  adam_patch_common_object["adam_patch_common_object"] --> adam_amr_object["adam_amr_object"]
  adam_patch_common_object["adam_patch_common_object"] --> adam_field_object["adam_field_object"]
  adam_patch_common_object["adam_patch_common_object"] --> adam_flail_object["adam_flail_object"]
  adam_patch_common_object["adam_patch_common_object"] --> adam_grid_object["adam_grid_object"]
  adam_patch_common_object["adam_patch_common_object"] --> adam_ib_object["adam_ib_object"]
  adam_patch_common_object["adam_patch_common_object"] --> adam_mpih_object["adam_mpih_object"]
  adam_patch_common_object["adam_patch_common_object"] --> adam_patch_bc_object["adam_patch_bc_object"]
  adam_patch_common_object["adam_patch_common_object"] --> adam_patch_ic_object["adam_patch_ic_object"]
  adam_patch_common_object["adam_patch_common_object"] --> adam_patch_io_object["adam_patch_io_object"]
  adam_patch_common_object["adam_patch_common_object"] --> adam_patch_time_object["adam_patch_time_object"]
  adam_patch_common_object["adam_patch_common_object"] --> penf["penf"]
```

## Contents

- [patch_common_object](#patch-common-object)
- [allocate_common](#allocate-common)
- [initialize_common](#initialize-common)

## Derived Types

### patch_common_object

Maxwell equations system class definition, common data to all backends.

**Inheritance**

```mermaid
classDiagram
  patch_common_object <|-- patch_cpu_object
```

#### Components

| Name | Type | Attributes | Description |
|------|------|------------|-------------|
| `mpih` | type([mpih_object](/api/src/lib/common/adam_mpih_object#mpih-object)) |  | MPI handler. |
| `adam` | type([adam_object](/api/src/lib/common/adam_adam_object#adam-object)) |  | ADAM. |
| `field` | type([field_object](/api/src/lib/common/adam_field_object#field-object)) | pointer | The field. |
| `grid` | type([grid_object](/api/src/lib/common/adam_grid_object#grid-object)) | pointer | The grid. |
| `amr` | type([amr_object](/api/src/lib/common/adam_amr_object#amr-object)) |  | AMR marker handler. |
| `ib` | type([ib_object](/api/src/lib/common/adam_ib_object#ib-object)) |  | Immersed Boundary (IB) handler. |
| `flail` | type([flail_object](/api/src/lib/common/adam_flail_object#flail-object)) |  | Linear algebra methods handler. |
| `io` | type([patch_io_object](/api/src/app/patch/common/adam_patch_io_object#patch-io-object)) |  | IO handler. |
| `ic` | type([patch_ic_object](/api/src/app/patch/common/adam_patch_ic_object#patch-ic-object)) |  | Initial Conditions (IC) handler. |
| `bc` | type([patch_bc_object](/api/src/app/patch/common/adam_patch_bc_object#patch-bc-object)) |  | Boundary Conditions (BC) handler. |
| `time` | type([patch_time_object](/api/src/app/patch/common/adam_patch_time_object#patch-time-object)) |  | Time handler. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Number of ghost cells. |
| `ni` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Number of cells in i direction. |
| `nj` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Number of cells in j direction. |
| `nk` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Number of cells in k direction. |
| `nb` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Total blocks number for MPI. |
| `blocks_number` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Actual blocks number. |
| `nv` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Number of conservative/primitive variables. |
| `q` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | allocatable | Potential field cell centered variable. |
| `r` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | allocatable | Rho function. |
| `q_name` | character(len=3) |  | Potential field name. |

#### Type-Bound Procedures

| Name | Attributes | Description |
|------|------------|-------------|
| `allocate_common` | pass(self) | Allocate common data. |
| `initialize_common` | pass(self) | Initialize the equation common data. |

## Subroutines

### allocate_common

Allocate common data.

```fortran
subroutine allocate_common(self)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([patch_common_object](/api/src/app/patch/common/adam_patch_common_object#patch-common-object)) | inout |  | The equation. |

**Call graph**

```mermaid
flowchart TD
  initialize_common["initialize_common"] --> allocate_common["allocate_common"]
  initialize_common["initialize_common"] --> allocate_common["allocate_common"]
  initialize_common["initialize_common"] --> allocate_common["allocate_common"]
  initialize_common["initialize_common"] --> allocate_common["allocate_common"]
  allocate_common["allocate_common"] --> allocate_variable["allocate_variable"]
  style allocate_common fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### initialize_common

Initialize the equation common data.

```fortran
subroutine initialize_common(self, field, filename, memory_avail, do_mpi_init, verbose)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([patch_common_object](/api/src/app/patch/common/adam_patch_common_object#patch-common-object)) | inout | target | The equation. |
| `field` | type([field_object](/api/src/lib/common/adam_field_object#field-object)) | inout |  | The field. |
| `filename` | character(len=*) | in |  | Input file name. |
| `memory_avail` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | value | Memory available for single MPI process. |
| `do_mpi_init` | logical | in | optional | Flag to activate MPI init call. |
| `verbose` | logical | in | optional | Trigger verbose output. |

**Call graph**

```mermaid
flowchart TD
  initialize["initialize"] --> initialize_common["initialize_common"]
  initialize["initialize"] --> initialize_common["initialize_common"]
  initialize["initialize"] --> initialize_common["initialize_common"]
  initialize["initialize"] --> initialize_common["initialize_common"]
  initialize["initialize"] --> initialize_common["initialize_common"]
  initialize["initialize"] --> initialize_common["initialize_common"]
  initialize["initialize"] --> initialize_common["initialize_common"]
  initialize_common["initialize_common"] --> allocate_common["allocate_common"]
  initialize_common["initialize_common"] --> associate_adam_data["associate_adam_data"]
  initialize_common["initialize_common"] --> compute_blocks_number["compute_blocks_number"]
  initialize_common["initialize_common"] --> initialize["initialize"]
  initialize_common["initialize_common"] --> print_message["print_message"]
  initialize_common["initialize_common"] --> prune["prune"]
  initialize_common["initialize_common"] --> refine_uniform["refine_uniform"]
  style initialize_common fill:#3e63dd,stroke:#99b,stroke-width:2px
```
