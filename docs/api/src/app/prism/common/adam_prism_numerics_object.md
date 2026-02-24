---
title: adam_prism_numerics_object
---

# adam_prism_numerics_object

> ADAM, PRISM (Plasma Research usIng Simulation Methods) numerics class definition, common backend.

**Source**: `src/app/prism/common/adam_prism_numerics_object.F90`

**Dependencies**

```mermaid
graph LR
  adam_prism_numerics_object["adam_prism_numerics_object"] --> adam_mpih_object["adam_mpih_object"]
  adam_prism_numerics_object["adam_prism_numerics_object"] --> adam_prism_parameters["adam_prism_parameters"]
  adam_prism_numerics_object["adam_prism_numerics_object"] --> finer["finer"]
  adam_prism_numerics_object["adam_prism_numerics_object"] --> penf["penf"]
```

## Contents

- [prism_numerics_object](#prism-numerics-object)
- [initialize](#initialize)
- [load_from_file](#load-from-file)
- [description](#description)

## Variables

| Name | Type | Attributes | Description |
|------|------|------------|-------------|
| `INI_SECTION_NAME` | character(len=8) | parameter | INI file section name containing fluid numerics. |
| `NUM_SCHEME_TIME_BLANES_MOAN` | character(len=11) | parameter | Blanes-Moan numerical scheme for time operator. |
| `NUM_SCHEME_TIME_CFM` | character(len=22) | parameter | Commutator-Free Magnus scheme for time operator. |
| `NUM_SCHEME_TIME_LEAPFROG` | character(len=8) | parameter | Leapfrog numerical scheme for time operator. |
| `NUM_SCHEME_TIME_RUNGE_KUTTA` | character(len=11) | parameter | Runge-Kutta numerical scheme for time operator. |
| `NUM_SCHEME_SPACE_FD_CENTERED` | character(len=11) | parameter | FD centered scheme for space operator. |
| `NUM_SCHEME_SPACE_FV_CENTERED` | character(len=11) | parameter | FV centered scheme for space operator. |
| `NUM_SCHEME_SPACE_WENO` | character(len=4) | parameter | WENO numerical scheme for space operator. |
| `RECONSTRUCTION_VARS_CONS` | character(len=12) | parameter | High-order reconstruction on conservative vars. |
| `RECONSTRUCTION_VARS_CHAR` | character(len=15) | parameter | High-order reconstruction on charact. vars. |
| `DIV_CORR_VAR_POISS` | character(len=7) | parameter | Poisson divergence correction. |
| `DIV_CORR_VAR_HYPER` | character(len=10) | parameter | Hyperbolic divergence correction. |

## Derived Types

### prism_numerics_object

PRISM numerics class definition.

#### Components

| Name | Type | Attributes | Description |
|------|------|------------|-------------|
| `mpih` | type([mpih_object](/api/src/third_party/FUNDAL/src/lib/fundal_mpih_object#mpih-object)) |  | MPI handler. |
| `scheme_time` | character(len=:) | allocatable | Numerical scheme for time operator [runge_kutta, leapfrog,...]. |
| `scheme_space` | character(len=:) | allocatable | Numerical scheme for space operator [weno, centered]. |
| `fdv_order` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) |  | Order of finite difference/volume schemes, general order. |
| `fdv_half_stencil` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) |  | Half stencil length of finite difference/volume schemes. |
| `fdv_half_stencils` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) |  | Half stencil length of fdv schemes for each derivative up to 6. |
| `reconstruction_vars` | character(len=:) | allocatable | Type of WENO reconstruction variables (cons., charct.,...). |
| `constrained_transport_D` | logical |  | Enable Constrained Transport Correction on D. |
| `constrained_transport_B` | logical |  | Enable Constrained Transport Correction on D. |
| `div_corr_var` | character(len=:) | allocatable | Type of divergence correction variables (poisson, hyperbolic,...). |

#### Type-Bound Procedures

| Name | Attributes | Description |
|------|------------|-------------|
| `description` | pass(self) | Return pretty-printed object description. |
| `initialize` | pass(self) | Initialize numerics. |
| `load_from_file` | pass(self) | Load config from file. |

## Subroutines

### initialize

Initialize the equation.

```fortran
subroutine initialize(self, file_parameters)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([prism_numerics_object](/api/src/app/prism/common/adam_prism_numerics_object#prism-numerics-object)) | inout |  | numerics. |
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
| `self` | class([prism_numerics_object](/api/src/app/prism/common/adam_prism_numerics_object#prism-numerics-object)) | inout |  | numerics. |
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
| `self` | class([prism_numerics_object](/api/src/app/prism/common/adam_prism_numerics_object#prism-numerics-object)) | in |  | numerics. |

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
