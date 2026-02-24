---
title: adam_fnl_maps_object
---

# adam_fnl_maps_object

> ADAM, maps class definition, FNL backend.

**Source**: `src/lib/fnl/adam_fnl_maps_object.F90`

**Dependencies**

```mermaid
graph LR
  adam_fnl_maps_object["adam_fnl_maps_object"] --> adam_fnl_mpih_object["adam_fnl_mpih_object"]
  adam_fnl_maps_object["adam_fnl_maps_object"] --> adam_maps_object["adam_maps_object"]
  adam_fnl_maps_object["adam_fnl_maps_object"] --> fundal["fundal"]
  adam_fnl_maps_object["adam_fnl_maps_object"] --> mpi["mpi"]
  adam_fnl_maps_object["adam_fnl_maps_object"] --> penf["penf"]
```

## Contents

- [maps_fnl_object](#maps-fnl-object)
- [copy_cpu_gpu](#copy-cpu-gpu)
- [initialize](#initialize)

## Derived Types

### maps_fnl_object

Maps class, FNL backend.

#### Components

| Name | Type | Attributes | Description |
|------|------|------------|-------------|
| `mpih` | type([mpih_object](/api/src/third_party/FUNDAL/src/lib/fundal_mpih_object#mpih-object)) |  | MPI handler, FNL backend. |
| `maps` | type([maps_object](/api/src/lib/common/adam_maps_object#maps-object)) | pointer | The maps. |
| `local_map_ghost_cell_gpu` | integer(kind=[I8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Local map for ghost cells updating, cells order. |
| `comm_map_recv_ghost_cell_gpu` | integer(kind=[I8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Communication map, `fec` information, cell order. |
| `comm_map_send_ghost_cell_gpu` | integer(kind=[I8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Communication map, `fec` information, cell order. |
| `send_buffer_ghost_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Send buffer of ghost cells. |
| `recv_buffer_ghost_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Receive buffer of ghost cells. |
| `local_map_bc_crown_gpu` | integer(kind=[I8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Local map for face BC ghost cells, "crown" order. |

#### Type-Bound Procedures

| Name | Attributes | Description |
|------|------------|-------------|
| `copy_cpu_gpu` | pass(self) | Copy data from (maps_object) CPU to (maps_fnl_object) GPU. |
| `initialize` | pass(self) | Initialize MPI handler data. |

## Subroutines

### copy_cpu_gpu

Copy data from (maps_object) CPU to (maps_fnl_object) GPU.

```fortran
subroutine copy_cpu_gpu(self, verbose)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([maps_fnl_object](/api/src/lib/fnl/adam_fnl_maps_object#maps-fnl-object)) | inout |  | The maps. |
| `verbose` | logical | in | optional | Flag to activate verbose mode. |

**Call graph**

```mermaid
flowchart TD
  amr_update["amr_update"] --> copy_cpu_gpu["copy_cpu_gpu"]
  amr_update["amr_update"] --> copy_cpu_gpu["copy_cpu_gpu"]
  amr_update["amr_update"] --> copy_cpu_gpu["copy_cpu_gpu"]
  copy_cpu_gpu["copy_cpu_gpu"] --> copy_cpu_gpu["copy_cpu_gpu"]
  copy_cpu_gpu["copy_cpu_gpu"] --> copy_cpu_gpu["copy_cpu_gpu"]
  copy_cpu_gpu["copy_cpu_gpu"] --> copy_cpu_gpu["copy_cpu_gpu"]
  copy_cpu_gpu["copy_cpu_gpu"] --> copy_cpu_gpu["copy_cpu_gpu"]
  copy_cpu_gpu["copy_cpu_gpu"] --> copy_cpu_gpu["copy_cpu_gpu"]
  copy_cpu_gpu["copy_cpu_gpu"] --> copy_cpu_gpu["copy_cpu_gpu"]
  copy_cpu_gpu["copy_cpu_gpu"] --> copy_cpu_gpu["copy_cpu_gpu"]
  initialize["initialize"] --> copy_cpu_gpu["copy_cpu_gpu"]
  initialize["initialize"] --> copy_cpu_gpu["copy_cpu_gpu"]
  initialize["initialize"] --> copy_cpu_gpu["copy_cpu_gpu"]
  initialize["initialize"] --> copy_cpu_gpu["copy_cpu_gpu"]
  initialize["initialize"] --> copy_cpu_gpu["copy_cpu_gpu"]
  initialize["initialize"] --> copy_cpu_gpu["copy_cpu_gpu"]
  initialize["initialize"] --> copy_cpu_gpu["copy_cpu_gpu"]
  load_restart_files["load_restart_files"] --> copy_cpu_gpu["copy_cpu_gpu"]
  load_restart_files["load_restart_files"] --> copy_cpu_gpu["copy_cpu_gpu"]
  load_restart_files["load_restart_files"] --> copy_cpu_gpu["copy_cpu_gpu"]
  load_restart_files["load_restart_files"] --> copy_cpu_gpu["copy_cpu_gpu"]
  set_initial_conditions["set_initial_conditions"] --> copy_cpu_gpu["copy_cpu_gpu"]
  set_initial_conditions["set_initial_conditions"] --> copy_cpu_gpu["copy_cpu_gpu"]
  set_initial_conditions["set_initial_conditions"] --> copy_cpu_gpu["copy_cpu_gpu"]
  set_initial_conditions["set_initial_conditions"] --> copy_cpu_gpu["copy_cpu_gpu"]
  copy_cpu_gpu["copy_cpu_gpu"] --> dev_assign_to_device["dev_assign_to_device"]
  style copy_cpu_gpu fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### initialize

Initialize maps.

```fortran
subroutine initialize(self, maps)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([maps_fnl_object](/api/src/lib/fnl/adam_fnl_maps_object#maps-fnl-object)) | inout |  | The maps, FNL backend. |
| `maps` | type([maps_object](/api/src/lib/common/adam_maps_object#maps-object)) | in | target | The maps. |

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
  initialize["initialize"] --> copy_cpu_gpu["copy_cpu_gpu"]
  style initialize fill:#3e63dd,stroke:#99b,stroke-width:2px
```
