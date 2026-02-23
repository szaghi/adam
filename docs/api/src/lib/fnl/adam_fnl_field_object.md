---
title: adam_fnl_field_object
---

# adam_fnl_field_object

> ADAM, field class definition, FNL backend.

**Source**: `src/lib/fnl/adam_fnl_field_object.F90`

**Dependencies**

```mermaid
graph LR
  adam_fnl_field_object["adam_fnl_field_object"] --> adam_common_library["adam_common_library"]
  adam_fnl_field_object["adam_fnl_field_object"] --> adam_fnl_field_kernels["adam_fnl_field_kernels"]
  adam_fnl_field_object["adam_fnl_field_object"] --> adam_fnl_maps_object["adam_fnl_maps_object"]
  adam_fnl_field_object["adam_fnl_field_object"] --> adam_fnl_mpih_object["adam_fnl_mpih_object"]
  adam_fnl_field_object["adam_fnl_field_object"] --> fundal["fundal"]
  adam_fnl_field_object["adam_fnl_field_object"] --> mpi["mpi"]
  adam_fnl_field_object["adam_fnl_field_object"] --> penf["penf"]
```

## Contents

- [field_fnl_object](#field-fnl-object)
- [compute_q_gradient](#compute-q-gradient)
- [copy_cpu_gpu](#copy-cpu-gpu)
- [copy_transpose_cpu_gpu](#copy-transpose-cpu-gpu)
- [copy_transpose_gpu_cpu](#copy-transpose-gpu-cpu)
- [initialize](#initialize)
- [update_ghost_local_gpu](#update-ghost-local-gpu)
- [update_ghost_mpi_gpu](#update-ghost-mpi-gpu)

## Derived Types

### field_fnl_object

Field class, FNL backend.

#### Components

| Name | Type | Attributes | Description |
|------|------|------------|-------------|
| `mpih` | type([mpih_object](/api/src/lib/common/adam_mpih_object#mpih-object)) |  | MPI handler. |
| `maps` | type([maps_fnl_object](/api/src/lib/fnl/adam_fnl_maps_object#maps-fnl-object)) |  | Maps handler. |
| `field` | type([field_object](/api/src/lib/common/adam_field_object#field-object)) | pointer | The field. |
| `q_t` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | allocatable | Transposed cell centered variables on CPU. |
| `q_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Field cell centered variables. |
| `q_t_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Transposed cell centered variables on GPU. |
| `fec_1_6_array_gpu` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Mapping fec1-26 to fec1-6 for boundaries (GPU). |
| `x_cell_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Cells x coordinates on GPU [nb,1-ngc:ni+ngc]. |
| `y_cell_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Cells y coordinates on GPU [nb,1-ngc:nj+ngc]. |
| `z_cell_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Cells z coordinates on GPU [nb,1-ngc:nk+ngc]. |
| `dxyz_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Delta cells GPU [nb,3]. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Number of ghost cells. |
| `ni` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Number of cells in i direction. |
| `nj` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Number of cells in j direction. |
| `nk` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Number of cells in k direction. |
| `nb` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Total blocks number for MPI. |
| `blocks_number` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Actual blocks number. |
| `nv` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Number of variables in q vector. |

#### Type-Bound Procedures

| Name | Attributes | Description |
|------|------------|-------------|
| `compute_q_gradient` | pass(self) | Compute maximum gradient module of q element of a block. |
| `copy_cpu_gpu` | pass(self) | Copy data from (field_object) CPU to (field_fnl_object) GPU. |
| `copy_transpose_cpu_gpu` | pass(self) | Transpose data from GPU to CPU. |
| `copy_transpose_gpu_cpu` | pass(self) | Transpose data from GPU to CPU. |
| `initialize` | pass(self) | Initialize field. |
| `update_ghost_local_gpu` | pass(self) | Update ghosts locally. |
| `update_ghost_mpi_gpu` | pass(self) | Update ghosts MPI. |

## Subroutines

### compute_q_gradient

Compute gradient (module) over q elements.

```fortran
subroutine compute_q_gradient(self, b, ivar, q_gpu, gradient)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([field_fnl_object](/api/src/lib/fnl/adam_fnl_field_object#field-fnl-object)) | in |  | The field. |
| `b` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Block index. |
| `ivar` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Index of q variable. |
| `q_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Field component to which apply gradient. |
| `gradient` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | out |  | Maximum gradient of q(ivar). |

**Call graph**

```mermaid
flowchart TD
  mark_by_grad_var["mark_by_grad_var"] --> compute_q_gradient["compute_q_gradient"]
  mark_by_grad_var["mark_by_grad_var"] --> compute_q_gradient["compute_q_gradient"]
  mark_by_grad_var["mark_by_grad_var"] --> compute_q_gradient["compute_q_gradient"]
  mark_by_grad_var["mark_by_grad_var"] --> compute_q_gradient["compute_q_gradient"]
  mark_by_grad_var["mark_by_grad_var"] --> compute_q_gradient["compute_q_gradient"]
  compute_q_gradient["compute_q_gradient"] --> compute_q_gradient_dev["compute_q_gradient_dev"]
  style compute_q_gradient fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### copy_cpu_gpu

Copy data from (field_object) CPU to (field_fnl_object) GPU.

```fortran
subroutine copy_cpu_gpu(self, verbose)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([field_fnl_object](/api/src/lib/fnl/adam_fnl_field_object#field-fnl-object)) | inout |  | The field. |
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
  copy_cpu_gpu["copy_cpu_gpu"] --> copy_cpu_gpu["copy_cpu_gpu"]
  copy_cpu_gpu["copy_cpu_gpu"] --> dev_assign_to_device["dev_assign_to_device"]
  copy_cpu_gpu["copy_cpu_gpu"] --> print_message["print_message"]
  style copy_cpu_gpu fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### copy_transpose_cpu_gpu

Copy transposed data from CPU to GPU.
 This routine is called by equation typically passing either q_gpu or q_aux_gpu.

```fortran
subroutine copy_transpose_cpu_gpu(self, nv, q_cpu, q_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([field_fnl_object](/api/src/lib/fnl/adam_fnl_field_object#field-fnl-object)) | inout |  | The field. |
| `nv` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of varibales. |
| `q_cpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Conservative variables on CPU. |
| `q_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | out |  | Conservative variables on GPU. |

**Call graph**

```mermaid
flowchart TD
  copy_cpu_gpu["copy_cpu_gpu"] --> copy_transpose_cpu_gpu["copy_transpose_cpu_gpu"]
  copy_cpu_gpu["copy_cpu_gpu"] --> copy_transpose_cpu_gpu["copy_transpose_cpu_gpu"]
  copy_cpu_gpu["copy_cpu_gpu"] --> copy_transpose_cpu_gpu["copy_transpose_cpu_gpu"]
  copy_cpu_gpu["copy_cpu_gpu"] --> copy_transpose_cpu_gpu["copy_transpose_cpu_gpu"]
  copy_transpose_cpu_gpu["copy_transpose_cpu_gpu"] --> dev_memcpy_to_device["dev_memcpy_to_device"]
  style copy_transpose_cpu_gpu fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### copy_transpose_gpu_cpu

Copy transposed data from GPU to CPU.
 This routine is called by equation typically passing either q_gpu or q_aux_gpu.

```fortran
subroutine copy_transpose_gpu_cpu(self, nv, q_gpu, q_cpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([field_fnl_object](/api/src/lib/fnl/adam_fnl_field_object#field-fnl-object)) | inout |  | The equation. |
| `nv` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of varibales. |
| `q_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Conservative variables on GPU. |
| `q_cpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout |  | Conservative variables on CPU. |

**Call graph**

```mermaid
flowchart TD
  copy_gpu_cpu["copy_gpu_cpu"] --> copy_transpose_gpu_cpu["copy_transpose_gpu_cpu"]
  copy_gpu_cpu["copy_gpu_cpu"] --> copy_transpose_gpu_cpu["copy_transpose_gpu_cpu"]
  copy_gpu_cpu["copy_gpu_cpu"] --> copy_transpose_gpu_cpu["copy_transpose_gpu_cpu"]
  copy_gpu_cpu["copy_gpu_cpu"] --> copy_transpose_gpu_cpu["copy_transpose_gpu_cpu"]
  copy_transpose_gpu_cpu["copy_transpose_gpu_cpu"] --> copy_transpose_gpu_cpu_dev["copy_transpose_gpu_cpu_dev"]
  copy_transpose_gpu_cpu["copy_transpose_gpu_cpu"] --> dev_memcpy_from_device["dev_memcpy_from_device"]
  style copy_transpose_gpu_cpu fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### initialize

Initialize field.

```fortran
subroutine initialize(self, field, nv_aux, verbose)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([field_fnl_object](/api/src/lib/fnl/adam_fnl_field_object#field-fnl-object)) | inout |  | The field. |
| `field` | type([field_object](/api/src/lib/common/adam_field_object#field-object)) | in | target | Field variable array. |
| `nv_aux` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | optional | Number of auxiliary variables. |
| `verbose` | logical | in | optional | Flag to activate verbose mode. |

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
  initialize["initialize"] --> allocate_variable["allocate_variable"]
  initialize["initialize"] --> copy_cpu_gpu["copy_cpu_gpu"]
  initialize["initialize"] --> dev_alloc["dev_alloc"]
  initialize["initialize"] --> dev_assign_to_device["dev_assign_to_device"]
  initialize["initialize"] --> initialize["initialize"]
  initialize["initialize"] --> print_message["print_message"]
  style initialize fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### update_ghost_local_gpu

Update (local) ghost cells.

```fortran
subroutine update_ghost_local_gpu(self, q_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([field_fnl_object](/api/src/lib/fnl/adam_fnl_field_object#field-fnl-object)) | in |  | The field. |
| `q_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout |  | Field component to be updated. |

**Call graph**

```mermaid
flowchart TD
  update_ghost["update_ghost"] --> update_ghost_local_gpu["update_ghost_local_gpu"]
  update_ghost["update_ghost"] --> update_ghost_local_gpu["update_ghost_local_gpu"]
  update_ghost["update_ghost"] --> update_ghost_local_gpu["update_ghost_local_gpu"]
  update_ghost_local_gpu["update_ghost_local_gpu"] --> update_ghost_local_gpu_dev["update_ghost_local_gpu_dev"]
  style update_ghost_local_gpu fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### update_ghost_mpi_gpu

Update ghost cells within other processes.

```fortran
subroutine update_ghost_mpi_gpu(self, q_gpu, step)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([field_fnl_object](/api/src/lib/fnl/adam_fnl_field_object#field-fnl-object)) | inout |  | The field. |
| `q_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout |  | Field component to be updated. |
| `step` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | optional | Step to be perfordmed in asyncronous comp. |

**Call graph**

```mermaid
flowchart TD
  update_ghost["update_ghost"] --> update_ghost_mpi_gpu["update_ghost_mpi_gpu"]
  update_ghost["update_ghost"] --> update_ghost_mpi_gpu["update_ghost_mpi_gpu"]
  update_ghost["update_ghost"] --> update_ghost_mpi_gpu["update_ghost_mpi_gpu"]
  update_ghost_mpi_gpu["update_ghost_mpi_gpu"] --> populate_send_buffer_ghost_gpu_dev["populate_send_buffer_ghost_gpu_dev"]
  update_ghost_mpi_gpu["update_ghost_mpi_gpu"] --> receive_recv_buffer_ghost_gpu_dev["receive_recv_buffer_ghost_gpu_dev"]
  style update_ghost_mpi_gpu fill:#3e63dd,stroke:#99b,stroke-width:2px
```
