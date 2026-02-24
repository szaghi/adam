---
title: adam_prism_external_fields_object
---

# adam_prism_external_fields_object

> ADAM, PRISM external fields definition, CPU backend.

**Source**: `src/app/prism/common/adam_prism_external_fields_object.F90`

**Dependencies**

```mermaid
graph LR
  adam_prism_external_fields_object["adam_prism_external_fields_object"] --> adam_field_object["adam_field_object"]
  adam_prism_external_fields_object["adam_prism_external_fields_object"] --> adam_mpih_object["adam_mpih_object"]
  adam_prism_external_fields_object["adam_prism_external_fields_object"] --> adam_prism_parameters["adam_prism_parameters"]
  adam_prism_external_fields_object["adam_prism_external_fields_object"] --> finer["finer"]
  adam_prism_external_fields_object["adam_prism_external_fields_object"] --> penf["penf"]
```

## Contents

- [prism_external_fields_object](#prism-external-fields-object)
- [add_external_fields_interface](#add-external-fields-interface)
- [sub_external_fields_interface](#sub-external-fields-interface)
- [initialize](#initialize)
- [load_from_file](#load-from-file)
- [add_external_fields_rmf](#add-external-fields-rmf)
- [sub_external_fields_rmf](#sub-external-fields-rmf)
- [description](#description)

## Variables

| Name | Type | Attributes | Description |
|------|------|------------|-------------|
| `INI_SECTION_NAME` | character(len=15) | parameter | INI (config) file section name. |
| `EF_TYPE_MAGNETIC_NOZZLE` | character(len=15) | parameter | Magnetic Nozzle. |
| `EF_TYPE_NONE` | character(len=4) | parameter | Disable external field. |
| `EF_TYPE_RMF` | character(len=3) | parameter | Rotating Magnetic Field. |
| `EF_TYPE_RMF_AND_MAGNETIC_NOZZLE` | character(len=23) | parameter | Rotating Magnetic Field/Magnetic Nozzle. |

## Derived Types

### prism_external_fields_object

PRISM external fields object.

#### Components

| Name | Type | Attributes | Description |
|------|------|------------|-------------|
| `mpih` | type([mpih_object](/api/src/third_party/FUNDAL/src/lib/fundal_mpih_object#mpih-object)) |  | MPI handler. |
| `ef_type` | character(len=99) |  | Field type. |
| `RMF_frequency` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) |  | Rotating magnetic field frequency. |
| `RMF_B_amplitude` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) |  | Rotating magnetic field amplitude. |
| `RMF_rotation_axis` | character(len=99) |  | Rotating magnetic field rotation axis (X, Y, Z). |
| `alpha` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) |  | RMF rotation axis coordinate 1 |
| `beta` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) |  | RMF rotation axis coordinate 2 |
| `gamm` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) |  | RMF rotation axis coordinate 3 |
| `add_external_fields` | procedure(add_external_fields_interface) | pass(self), pointer | Add external fields. |
| `sub_external_fields` | procedure(sub_external_fields_interface) | pass(self), pointer | Subtract external fields. |

#### Type-Bound Procedures

| Name | Attributes | Description |
|------|------------|-------------|
| `description` | pass(self) | Return pretty-printed object description. |
| `initialize` | pass(self) | Initialize IC. |
| `load_from_file` | pass(self) | Load config from file. |
| `add_external_fields_rmf` | pass(self) | Add rotating magnetic field to the field. |
| `sub_external_fields_rmf` | pass(self) | Add rotating magnetic field to the field. |

## Interfaces

### add_external_fields_interface

### sub_external_fields_interface

## Subroutines

### initialize

Initialize external fields.

```fortran
subroutine initialize(self, file_parameters)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([prism_external_fields_object](/api/src/app/prism/common/adam_prism_external_fields_object#prism-external-fields-object)) | inout |  | External fields. |
| `file_parameters` | type([file_ini](/api/src/third_party/FiNeR/src/lib/finer_file_ini_t#file-ini)) | in |  | Simulation parameters ini file handler. |

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
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  style initialize fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### load_from_file

Load config from file.

```fortran
subroutine load_from_file(self, file_parameters, go_on_fail)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([prism_external_fields_object](/api/src/app/prism/common/adam_prism_external_fields_object#prism-external-fields-object)) | inout |  | External fields. |
| `file_parameters` | type([file_ini](/api/src/third_party/FiNeR/src/lib/finer_file_ini_t#file-ini)) | in |  | Simulation parameters ini file handler. |
| `go_on_fail` | logical | in | optional | Go on if load fails. |

**Call graph**

```mermaid
flowchart TD
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  initialize["initialize"] --> load_from_file["load_from_file"]
  load_from_file["load_from_file"] --> load_from_file["load_from_file"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> get["get"]
  load_from_file["load_from_file"] --> print_message["print_message"]
  style load_from_file fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### add_external_fields_rmf

Add rotating magnetic field to the field.

```fortran
subroutine add_external_fields_rmf(self, field, time, dt, gamm, q)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([prism_external_fields_object](/api/src/app/prism/common/adam_prism_external_fields_object#prism-external-fields-object)) | inout |  | External fields. |
| `field` | type([field_object](/api/src/lib/common/adam_field_object#field-object)) | inout |  | The field. |
| `time` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Current simulation time. |
| `dt` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | optional | Time step. |
| `gamm` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | optional | Gamma values of RK SSP |
| `q` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout |  | Primitive variables. |

### sub_external_fields_rmf

Add rotating magnetic field to the field.

```fortran
subroutine sub_external_fields_rmf(self, field, time, dt, gamm, q)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([prism_external_fields_object](/api/src/app/prism/common/adam_prism_external_fields_object#prism-external-fields-object)) | inout |  | External fields. |
| `field` | type([field_object](/api/src/lib/common/adam_field_object#field-object)) | inout |  | The field. |
| `time` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Current simulation time. |
| `dt` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | optional | Time step. |
| `gamm` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | optional | Gamma values of RK SSP |
| `q` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout |  | Primitive variables. |

## Functions

### description

Return a pretty-formatted object description.

**Attributes**: pure

**Returns**: `character(len=:)`

```fortran
function description(self) result(desc)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([prism_external_fields_object](/api/src/app/prism/common/adam_prism_external_fields_object#prism-external-fields-object)) | in |  | External fields. |

**Call graph**

```mermaid
flowchart TD
  description["description"] --> description["description"]
  description["description"] --> description["description"]
  description["description"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  initialize["initialize"] --> description["description"]
  description["description"] --> str["str"]
  style description fill:#3e63dd,stroke:#99b,stroke-width:2px
```
