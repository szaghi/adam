---
title: adam_mpih_object
---

# adam_mpih_object

> ADAM, MPI handler class definition.

**Source**: `src/lib/common/adam_mpih_object.F90`

**Dependencies**

```mermaid
graph LR
  adam_mpih_object["adam_mpih_object"] --> iso_c_binding["iso_c_binding"]
  adam_mpih_object["adam_mpih_object"] --> iso_fortran_env["iso_fortran_env"]
  adam_mpih_object["adam_mpih_object"] --> mpi["mpi"]
  adam_mpih_object["adam_mpih_object"] --> penf["penf"]
```

## Contents

- [mpih_object](#mpih-object)
- [abort](#abort)
- [barrier](#barrier)
- [error_stop](#error-stop)
- [finalize](#finalize)
- [initialize](#initialize)
- [print_message](#print-message)
- [tic](#tic)
- [description](#description)
- [tictoc_timing](#tictoc-timing)
- [toc](#toc)

## Derived Types

### mpih_object

MPI handler class.

**Inheritance**

```mermaid
classDiagram
  mpih_object <|-- mpih_gmp_object
  mpih_object <|-- mpih_nvf_object
```

#### Components

| Name | Type | Attributes | Description |
|------|------|------------|-------------|
| `myrank` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) |  | MPI rank process. |
| `myrankstr` | character(len=:) | allocatable | MPI rank process stringified. |
| `procs_number` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) |  | Number of MPI processes. |
| `memory_avail` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) |  | CPU memory available (GB) for each process. |
| `error` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) |  | Error traping flag. |
| `timing` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) |  | Tic toc timing. |
| `tictoc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) |  | Next is tic or toc? |
| `req_send_recv` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | allocatable | MPI request receive flags. |

#### Type-Bound Procedures

| Name | Attributes | Description |
|------|------------|-------------|
| `abort` | pass(self) | Handy MPI abort wrapper. |
| `barrier` | pass(self) | Handy MPI barrier wrapper. |
| `description` | pass(self) | Return pretty-printed object description. |
| `error_stop` | pass(self) | Stop run with error output. |
| `finalize` | pass(self) | Handy MPI finalize wrapper. |
| `initialize` | pass(self) | Initialize MPI handler data. |
| `print_message` | pass(self) | Print a message on stdout with rank prefix. |
| `tictoc_timing` | pass(self) | Return the last tic toc timing. |
| `tic` | pass(self) | Start a tic toc timing. |
| `toc` | pass(self) | Stop  a tic toc timing. |

## Subroutines

### abort

Handy MPI abort wrapper.

```fortran
subroutine abort(self, error_code, msg)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([mpih_object](/api/src/lib/common/adam_mpih_object#mpih-object)) | inout |  | MPI handler. |
| `error_code` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | optional | Abort error code. |
| `msg` | character(len=*) | in | optional | Error message. |

### barrier

Handy MPI barrier wrapper.

```fortran
subroutine barrier(self, tictoc, timing, single)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([mpih_object](/api/src/lib/common/adam_mpih_object#mpih-object)) | inout |  | MPI handler. |
| `tictoc` | logical | in | optional | Activate tic toc timing between 2 barrier calls. |
| `timing` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | out | optional | Current timing. |
| `single` | logical | in | optional | Single tictoc for one-shot timing. |

**Call graph**

```mermaid
flowchart TD
  integrate["integrate"] --> barrier["barrier"]
  integrate_eikonal["integrate_eikonal"] --> barrier["barrier"]
  integrate_eikonal["integrate_eikonal"] --> barrier["barrier"]
  integrate_eikonal["integrate_eikonal"] --> barrier["barrier"]
  integrate_eikonal["integrate_eikonal"] --> barrier["barrier"]
  integrate_eikonal["integrate_eikonal"] --> barrier["barrier"]
  integrate_eikonal["integrate_eikonal"] --> barrier["barrier"]
  save_hdf5["save_hdf5"] --> barrier["barrier"]
  save_hdf5["save_hdf5"] --> barrier["barrier"]
  save_hdf5["save_hdf5"] --> barrier["barrier"]
  save_hdf5["save_hdf5"] --> barrier["barrier"]
  save_mat["save_mat"] --> barrier["barrier"]
  save_restart_files["save_restart_files"] --> barrier["barrier"]
  save_restart_files["save_restart_files"] --> barrier["barrier"]
  save_restart_files["save_restart_files"] --> barrier["barrier"]
  save_restart_files["save_restart_files"] --> barrier["barrier"]
  save_restart_files["save_restart_files"] --> barrier["barrier"]
  save_restart_files["save_restart_files"] --> barrier["barrier"]
  save_restart_files["save_restart_files"] --> barrier["barrier"]
  save_xh5f["save_xh5f"] --> barrier["barrier"]
  save_xh5f["save_xh5f"] --> barrier["barrier"]
  save_xh5f["save_xh5f"] --> barrier["barrier"]
  simulate["simulate"] --> barrier["barrier"]
  simulate["simulate"] --> barrier["barrier"]
  simulate["simulate"] --> barrier["barrier"]
  simulate["simulate"] --> barrier["barrier"]
  simulate["simulate"] --> barrier["barrier"]
  simulate["simulate"] --> barrier["barrier"]
  style barrier fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### error_stop

Stop run with error output.

```fortran
subroutine error_stop(self, msg)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([mpih_object](/api/src/lib/common/adam_mpih_object#mpih-object)) | inout |  | MPI handler. |
| `msg` | character(len=*) | in | optional | Error message. |

**Call graph**

```mermaid
flowchart TD
  initialize["initialize"] --> error_stop["error_stop"]
  initialize["initialize"] --> error_stop["error_stop"]
  initialize["initialize"] --> error_stop["error_stop"]
  initialize["initialize"] --> error_stop["error_stop"]
  initialize["initialize"] --> error_stop["error_stop"]
  initialize["initialize"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_file["load_from_file"] --> error_stop["error_stop"]
  load_from_ini_file["load_from_ini_file"] --> error_stop["error_stop"]
  error_stop["error_stop"] --> finalize["finalize"]
  style error_stop fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### finalize

Handy MPI finalize wrapper.

```fortran
subroutine finalize(self)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([mpih_object](/api/src/lib/common/adam_mpih_object#mpih-object)) | inout |  | MPI handler. |

**Call graph**

```mermaid
flowchart TD
  error_stop["error_stop"] --> finalize["finalize"]
  error_stop["error_stop"] --> finalize["finalize"]
  export_vtk_file["export_vtk_file"] --> finalize["finalize"]
  finalize["finalize"] --> finalize["finalize"]
  finalize["finalize"] --> finalize["finalize"]
  finalize["finalize"] --> finalize["finalize"]
  finalize["finalize"] --> finalize["finalize"]
  initialize["initialize"] --> finalize["finalize"]
  save_vtk["save_vtk"] --> finalize["finalize"]
  simulate["simulate"] --> finalize["finalize"]
  simulate["simulate"] --> finalize["finalize"]
  simulate["simulate"] --> finalize["finalize"]
  style finalize fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### initialize

Initialize MPI handler data.

```fortran
subroutine initialize(self, do_mpi_init, do_device_init, myrankstr_char_length, verbose)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([mpih_object](/api/src/lib/common/adam_mpih_object#mpih-object)) | inout |  | MPI handler. |
| `do_mpi_init` | logical | in | optional | Flag to activate MPI init call. |
| `do_device_init` | logical | in | optional | Flag to activate device init call (used by backends). |
| `myrankstr_char_length` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | optional | MPI ID string length. |
| `verbose` | logical | in | optional | Trigger verbose output. |

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
  initialize["initialize"] --> get_memory_info["get_memory_info"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> strz["strz"]
  style initialize fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### print_message

Print a message on stdout with rank prefix.

```fortran
subroutine print_message(self, msg)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([mpih_object](/api/src/lib/common/adam_mpih_object#mpih-object)) | in |  | MPI handler. |
| `msg` | character(len=*) | in |  | Message to print. |

**Call graph**

```mermaid
flowchart TD
  allocate_cpu["allocate_cpu"] --> print_message["print_message"]
  allocate_cpu["allocate_cpu"] --> print_message["print_message"]
  allocate_cpu["allocate_cpu"] --> print_message["print_message"]
  allocate_gpu["allocate_gpu"] --> print_message["print_message"]
  allocate_gpu["allocate_gpu"] --> print_message["print_message"]
  allocate_gpu["allocate_gpu"] --> print_message["print_message"]
  allocate_gpu["allocate_gpu"] --> print_message["print_message"]
  amr_update["amr_update"] --> print_message["print_message"]
  amr_update["amr_update"] --> print_message["print_message"]
  amr_update["amr_update"] --> print_message["print_message"]
  amr_update["amr_update"] --> print_message["print_message"]
  amr_update["amr_update"] --> print_message["print_message"]
  amr_update["amr_update"] --> print_message["print_message"]
  amr_update["amr_update"] --> print_message["print_message"]
  check_blocks_number["check_blocks_number"] --> print_message["print_message"]
  compute_phi["compute_phi"] --> print_message["print_message"]
  compute_phi["compute_phi"] --> print_message["print_message"]
  compute_phi["compute_phi"] --> print_message["print_message"]
  compute_phi["compute_phi"] --> print_message["print_message"]
  compute_phi["compute_phi"] --> print_message["print_message"]
  compute_phi["compute_phi"] --> print_message["print_message"]
  compute_phi["compute_phi"] --> print_message["print_message"]
  compute_phi_all_solids["compute_phi_all_solids"] --> print_message["print_message"]
  copy_cpu_gpu["copy_cpu_gpu"] --> print_message["print_message"]
  copy_cpu_gpu["copy_cpu_gpu"] --> print_message["print_message"]
  copy_cpu_gpu["copy_cpu_gpu"] --> print_message["print_message"]
  copy_cpu_gpu["copy_cpu_gpu"] --> print_message["print_message"]
  copy_cpu_gpu["copy_cpu_gpu"] --> print_message["print_message"]
  copy_cpu_gpu["copy_cpu_gpu"] --> print_message["print_message"]
  copy_cpu_gpu["copy_cpu_gpu"] --> print_message["print_message"]
  copy_gpu_cpu["copy_gpu_cpu"] --> print_message["print_message"]
  impose_ct_correction["impose_ct_correction"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize["initialize"] --> print_message["print_message"]
  initialize_common["initialize_common"] --> print_message["print_message"]
  initialize_common["initialize_common"] --> print_message["print_message"]
  initialize_common["initialize_common"] --> print_message["print_message"]
  initialize_common["initialize_common"] --> print_message["print_message"]
  integrate["integrate"] --> print_message["print_message"]
  load_blocks["load_blocks"] --> print_message["print_message"]
  load_from_file["load_from_file"] --> print_message["print_message"]
  load_from_file["load_from_file"] --> print_message["print_message"]
  load_from_file["load_from_file"] --> print_message["print_message"]
  load_from_file["load_from_file"] --> print_message["print_message"]
  load_nodes["load_nodes"] --> print_message["print_message"]
  make_comm_local_maps["make_comm_local_maps"] --> print_message["print_message"]
  make_comm_local_maps_ghost["make_comm_local_maps_ghost"] --> print_message["print_message"]
  make_comm_local_maps_ghost_bc["make_comm_local_maps_ghost_bc"] --> print_message["print_message"]
  make_local_maps_bc["make_local_maps_bc"] --> print_message["print_message"]
  move_phi["move_phi"] --> print_message["print_message"]
  move_phi["move_phi"] --> print_message["print_message"]
  move_phi["move_phi"] --> print_message["print_message"]
  move_phi["move_phi"] --> print_message["print_message"]
  move_phi["move_phi"] --> print_message["print_message"]
  move_phi["move_phi"] --> print_message["print_message"]
  print_code_topology["print_code_topology"] --> print_message["print_message"]
  prune["prune"] --> print_message["print_message"]
  refine_uniform["refine_uniform"] --> print_message["print_message"]
  save_hdf5["save_hdf5"] --> print_message["print_message"]
  save_hdf5["save_hdf5"] --> print_message["print_message"]
  save_hdf5["save_hdf5"] --> print_message["print_message"]
  save_hdf5["save_hdf5"] --> print_message["print_message"]
  save_mat["save_mat"] --> print_message["print_message"]
  save_nodes["save_nodes"] --> print_message["print_message"]
  save_restart_files["save_restart_files"] --> print_message["print_message"]
  save_restart_files["save_restart_files"] --> print_message["print_message"]
  save_restart_files["save_restart_files"] --> print_message["print_message"]
  save_restart_files["save_restart_files"] --> print_message["print_message"]
  save_restart_files["save_restart_files"] --> print_message["print_message"]
  save_restart_files["save_restart_files"] --> print_message["print_message"]
  save_restart_files["save_restart_files"] --> print_message["print_message"]
  save_xh5f["save_xh5f"] --> print_message["print_message"]
  save_xh5f["save_xh5f"] --> print_message["print_message"]
  save_xh5f["save_xh5f"] --> print_message["print_message"]
  simulate["simulate"] --> print_message["print_message"]
  simulate["simulate"] --> print_message["print_message"]
  simulate["simulate"] --> print_message["print_message"]
  simulate["simulate"] --> print_message["print_message"]
  simulate["simulate"] --> print_message["print_message"]
  simulate["simulate"] --> print_message["print_message"]
  style print_message fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### tic

Start a tic toc timing.

```fortran
subroutine tic(self)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([mpih_object](/api/src/lib/common/adam_mpih_object#mpih-object)) | inout |  | MPI handler. |

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
| `self` | class([mpih_object](/api/src/lib/common/adam_mpih_object#mpih-object)) | in |  | MPI handler. |

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

### tictoc_timing

Return the last tic toc timing.

**Returns**: real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables))

```fortran
function tictoc_timing(self) result(timing)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([mpih_object](/api/src/lib/common/adam_mpih_object#mpih-object)) | in |  | MPI handler. |

**Call graph**

```mermaid
flowchart TD
  toc["toc"] --> tictoc_timing["tictoc_timing"]
  toc["toc"] --> tictoc_timing["tictoc_timing"]
  style tictoc_timing fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### toc

Stop a tic toc timing.

**Returns**: real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables))

```fortran
function toc(self) result(timing)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([mpih_object](/api/src/lib/common/adam_mpih_object#mpih-object)) | inout |  | MPI handler. |

**Call graph**

```mermaid
flowchart TD
  toc["toc"] --> tictoc_timing["tictoc_timing"]
  style toc fill:#3e63dd,stroke:#99b,stroke-width:2px
```
