---
title: adam_prism_fnl_coil_object
---

# adam_prism_fnl_coil_object

> ADAM, PRISM coil source definition, FNL backend.

**Source**: `src/app/prism/fnl/adam_prism_fnl_coil_object.F90`

**Dependencies**

```mermaid
graph LR
  adam_prism_fnl_coil_object["adam_prism_fnl_coil_object"] --> adam_fnl_mpih_object["adam_fnl_mpih_object"]
  adam_prism_fnl_coil_object["adam_prism_fnl_coil_object"] --> adam_prism_coil_object["adam_prism_coil_object"]
  adam_prism_fnl_coil_object["adam_prism_fnl_coil_object"] --> adam_prism_parameters["adam_prism_parameters"]
  adam_prism_fnl_coil_object["adam_prism_fnl_coil_object"] --> fundal["fundal"]
  adam_prism_fnl_coil_object["adam_prism_fnl_coil_object"] --> penf["penf"]
```

## Contents

- [prism_fnl_coil_object](#prism-fnl-coil-object)
- [copy_cpu_gpu](#copy-cpu-gpu)
- [copy_gpu_cpu](#copy-gpu-cpu)
- [initialize](#initialize)
- [compute_coils_current_dev](#compute-coils-current-dev)

## Derived Types

### prism_fnl_coil_object

ADAM, PRISM coil source definition, FNL (GPU) backend.

#### Components

| Name | Type | Attributes | Description |
|------|------|------------|-------------|
| `mpih` | type([mpih_object](/api/src/third_party/FUNDAL/src/lib/fundal_mpih_object#mpih-object)) |  | MPI handler. |
| `coil` | type([prism_coil_object](/api/src/app/prism/common/adam_prism_coil_object#prism-coil-object)) | pointer | Coil common handler. |
| `A_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Current amplitude (A) |
| `f_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Current frequency, if AC (Hz) |
| `phase_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Current initial phase, if AC |
| `d_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Coil wire diameter |
| `J_vec_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Matrice contenente versori corrente spire (se assente = 0) |
| `coil_flag_gpu` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Matrice contenente informazioni su quale spira pass. |
| `blocks_number` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Actual blocks number. |
| `nb` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Total blocks number for MPI. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Number of ghost cells. |
| `ni` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Number of cells in i direction. |
| `nj` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Number of cells in j direction. |
| `nk` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | pointer | Number of cells in k direction. |

#### Type-Bound Procedures

| Name | Attributes | Description |
|------|------------|-------------|
| `copy_cpu_gpu` | pass(self) | Copy data from CPU to GPU. |
| `copy_gpu_cpu` | pass(self) | Copy data from GPU to CPU. |
| `initialize` | pass(self) | Initialize class. |

## Subroutines

### copy_cpu_gpu

Copy data from CPU to GPU.

```fortran
subroutine copy_cpu_gpu(self, verbose)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([prism_fnl_coil_object](/api/src/app/prism/fnl/adam_prism_fnl_coil_object#prism-fnl-coil-object)) | inout |  | The field. |
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
  copy_cpu_gpu["copy_cpu_gpu"] --> dev_memcpy_to_device["dev_memcpy_to_device"]
  style copy_cpu_gpu fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### copy_gpu_cpu

Copy data from GPU to CPU.

```fortran
subroutine copy_gpu_cpu(self, verbose)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([prism_fnl_coil_object](/api/src/app/prism/fnl/adam_prism_fnl_coil_object#prism-fnl-coil-object)) | inout |  | The field. |
| `verbose` | logical | in | optional | Flag to activate verbose mode. |

**Call graph**

```mermaid
flowchart TD
  amr_update["amr_update"] --> copy_gpu_cpu["copy_gpu_cpu"]
  amr_update["amr_update"] --> copy_gpu_cpu["copy_gpu_cpu"]
  amr_update["amr_update"] --> copy_gpu_cpu["copy_gpu_cpu"]
  copy_gpu_cpu["copy_gpu_cpu"] --> copy_gpu_cpu["copy_gpu_cpu"]
  save_simulation_data["save_simulation_data"] --> copy_gpu_cpu["copy_gpu_cpu"]
  save_simulation_data["save_simulation_data"] --> copy_gpu_cpu["copy_gpu_cpu"]
  save_simulation_data["save_simulation_data"] --> copy_gpu_cpu["copy_gpu_cpu"]
  save_simulation_data["save_simulation_data"] --> copy_gpu_cpu["copy_gpu_cpu"]
  copy_gpu_cpu["copy_gpu_cpu"] --> dev_memcpy_from_device["dev_memcpy_from_device"]
  style copy_gpu_cpu fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### initialize

Initialize class.

```fortran
subroutine initialize(self, coil, blocks_number, nb, ngc, ni, nj, nk)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `self` | class([prism_fnl_coil_object](/api/src/app/prism/fnl/adam_prism_fnl_coil_object#prism-fnl-coil-object)) | inout |  | Coils. |
| `coil` | class([prism_coil_object](/api/src/app/prism/common/adam_prism_coil_object#prism-coil-object)) | in | target | Coils on host. |
| `blocks_number` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | target | Actual blocks number. |
| `nb` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | target | Maximum blocks number. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | target | Number of ghost cells. |
| `ni` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | target | Number of cells in i direction. |
| `nj` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | target | Number of cells in j direction. |
| `nk` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | target | Number of cells in k direction. |

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
  initialize["initialize"] --> dev_alloc["dev_alloc"]
  style initialize fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### compute_coils_current_dev

Compute current coils sources, device kernel.

```fortran
subroutine compute_coils_current_dev(ni, nj, nk, ngc, blocks_number, time_s, td, A_gpu, f_gpu, phase_gpu, coil_flag_gpu, J_vec_gpu, var_Jx, var_Jy, var_Jz, q_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `ni` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in I direction. |
| `nj` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in J direction. |
| `nk` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in K direction. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Ghost grid number. |
| `blocks_number` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of blocks. |
| `time_s` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Local time. |
| `td` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Delay coil start. |
| `A_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Current amplitude (A) |
| `f_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Current frequency, if AC (Hz) |
| `phase_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Current initial phase, if AC |
| `coil_flag_gpu` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Matrice contenente informazioni su quale spira pass. |
| `J_vec_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Matrice contenente versori corrente spire. |
| `var_Jx` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Indices of current density components in q vector. |
| `var_Jy` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Indices of current density components in q vector. |
| `var_Jz` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Indices of current density components in q vector. |
| `q_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout |  | Conservative variables on GPU. |

**Call graph**

```mermaid
flowchart TD
  compute_coils_current["compute_coils_current"] --> compute_coils_current_dev["compute_coils_current_dev"]
  style compute_coils_current_dev fill:#3e63dd,stroke:#99b,stroke-width:2px
```
