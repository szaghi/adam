---
title: adam_maps_gmp_object
---

# adam_maps_gmp_object

> ADAM, maps class definition, GMP backend.

**Source**: `src/lib/gmp/adam_maps_gmp_object.F90`

**Dependencies**

```mermaid
graph LR
  adam_maps_gmp_object["adam_maps_gmp_object"] --> adam_maps_object["adam_maps_object"]
  adam_maps_gmp_object["adam_maps_gmp_object"] --> adam_memory_gmp_library["adam_memory_gmp_library"]
  adam_maps_gmp_object["adam_maps_gmp_object"] --> adam_mpih_gmp_object["adam_mpih_gmp_object"]
  adam_maps_gmp_object["adam_maps_gmp_object"] --> mpi["mpi"]
  adam_maps_gmp_object["adam_maps_gmp_object"] --> penf["penf"]
```

## Contents

- [maps_gmp_object](#maps-gmp-object)
- [copy_cpu_gpu](#copy-cpu-gpu)
- [initialize](#initialize)

## Derived Types

### maps_gmp_object

Maps class, GMP backend.

#### Components

| Name | Type | Attributes | Description |
|------|------|------------|-------------|
| `mpih` | type([mpih_gmp_object](/api/src/lib/gmp/adam_mpih_gmp_object#mpih-gmp-object)) | pointer | MPI handler, GMP backend. |
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
| `copy_cpu_gpu` | pass(self) | Copy data from (maps_object) CPU to (maps_gmp_object) GPU. |
| `initialize` | pass(self) | Initialize MPI handler data. |

## Subroutines

### copy_cpu_gpu

Copy data from (maps_object) CPU to (maps_gmp_object) GPU.

```fortran
subroutine copy_cpu_gpu(self, verbose)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([maps_gmp_object](/api/src/lib/gmp/adam_maps_gmp_object#maps-gmp-object)) | inout |  | The maps. |
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
  copy_cpu_gpu["copy_cpu_gpu"] --> assign_allocatable_gpu["assign_allocatable_gpu"]
  copy_cpu_gpu["copy_cpu_gpu"] --> print_message["print_message"]
  style copy_cpu_gpu fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### initialize

Initialize maps.

```fortran
subroutine initialize(self, mpih, maps)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([maps_gmp_object](/api/src/lib/gmp/adam_maps_gmp_object#maps-gmp-object)) | inout |  | The maps, gmp backend. |
| `mpih` | type([mpih_gmp_object](/api/src/lib/gmp/adam_mpih_gmp_object#mpih-gmp-object)) | in | target | MPI handler, GMP backend. |
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
  initialize["initialize"] --> print_message["print_message"]
  style initialize fill:#3e63dd,stroke:#99b,stroke-width:2px
```
