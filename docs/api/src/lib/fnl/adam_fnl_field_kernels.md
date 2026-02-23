---
title: adam_fnl_field_kernels
---

# adam_fnl_field_kernels

> ADAM, field class FNL kernels (FNL backend of [field_object](/api/src/lib/common/adam_field_object#field-object)).

**Source**: `src/lib/fnl/adam_fnl_field_kernels.F90`

**Dependencies**

```mermaid
graph LR
  adam_fnl_field_kernels["adam_fnl_field_kernels"] --> fundal["fundal"]
  adam_fnl_field_kernels["adam_fnl_field_kernels"] --> penf["penf"]
```

## Contents

- [compute_q_gradient_dev](#compute-q-gradient-dev)
- [compute_normL2_residuals_dev](#compute-norml2-residuals-dev)
- [copy_transpose_gpu_cpu_dev](#copy-transpose-gpu-cpu-dev)
- [populate_send_buffer_ghost_gpu_dev](#populate-send-buffer-ghost-gpu-dev)
- [receive_recv_buffer_ghost_gpu_dev](#receive-recv-buffer-ghost-gpu-dev)
- [update_ghost_local_gpu_dev](#update-ghost-local-gpu-dev)

## Subroutines

### compute_q_gradient_dev

Compute gradient of q(ivar).

```fortran
subroutine compute_q_gradient_dev(b, ni, nj, nk, ngc, dx, dy, dz, q_gpu, ivar, gradient)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `b` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Block index. |
| `ni` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in I direction. |
| `nj` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in J direction. |
| `nk` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in K direction. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Ghost cells number. |
| `dx` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | X space step. |
| `dy` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Y space step. |
| `dz` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Z space step. |
| `q_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Field component to which apply gradient. |
| `ivar` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Ghost cells number. |
| `gradient` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | out |  | Maximum gradient of q. |

**Call graph**

```mermaid
flowchart TD
  compute_q_gradient["compute_q_gradient"] --> compute_q_gradient_dev["compute_q_gradient_dev"]
  style compute_q_gradient_dev fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### compute_normL2_residuals_dev

Compute L2 norm of residuals.

```fortran
subroutine compute_normL2_residuals_dev(ni, nj, nk, ngc, nv, blocks_number, dq_gpu, norm)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `ni` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in I direction. |
| `nj` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in J direction. |
| `nk` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in K direction. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Ghost grid number. |
| `nv` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of states variables. |
| `blocks_number` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of blocks. |
| `dq_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Residuals. |
| `norm` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout |  | Residuals norm. |

**Call graph**

```mermaid
flowchart TD
  save_residuals["save_residuals"] --> compute_normL2_residuals_dev["compute_normL2_residuals_dev"]
  save_residuals["save_residuals"] --> compute_normL2_residuals_dev["compute_normL2_residuals_dev"]
  style compute_normL2_residuals_dev fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### copy_transpose_gpu_cpu_dev

Copy transposed data from GPU to CPU by CUF threads.

```fortran
subroutine copy_transpose_gpu_cpu_dev(ni, nj, nk, ngc, nv, blocks_number, q_gpu, q_t_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `ni` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in I direction. |
| `nj` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in J direction. |
| `nk` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in K direction. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Ghost cells number. |
| `nv` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of conservative varibales. |
| `blocks_number` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of blocks. |
| `q_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Conservative variables on GPU. |
| `q_t_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout |  | Conservative (transposed) variables on GPU. |

**Call graph**

```mermaid
flowchart TD
  copy_transpose_gpu_cpu["copy_transpose_gpu_cpu"] --> copy_transpose_gpu_cpu_dev["copy_transpose_gpu_cpu_dev"]
  style copy_transpose_gpu_cpu_dev fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### populate_send_buffer_ghost_gpu_dev

Polulate send buffer ghost GPU.

```fortran
subroutine populate_send_buffer_ghost_gpu_dev(ngc, comm_map_send_ghost_cell_gpu, send_buffer_ghost_gpu, q_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Ghost cells number. |
| `comm_map_send_ghost_cell_gpu` | integer(kind=[I8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | pointer | Comm map, cell information. |
| `send_buffer_ghost_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout | pointer | Send buffer of ghost cells. |
| `q_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout |  | Field component to be updated. |

**Call graph**

```mermaid
flowchart TD
  update_ghost_mpi_gpu["update_ghost_mpi_gpu"] --> populate_send_buffer_ghost_gpu_dev["populate_send_buffer_ghost_gpu_dev"]
  style populate_send_buffer_ghost_gpu_dev fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### receive_recv_buffer_ghost_gpu_dev

Receive recv buffer ghost GPU.

```fortran
subroutine receive_recv_buffer_ghost_gpu_dev(ngc, comm_map_recv_ghost_cell_gpu, recv_buffer_ghost_gpu, q_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Ghost cells number. |
| `comm_map_recv_ghost_cell_gpu` | integer(kind=[I8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | pointer | Comm map, cell information. |
| `recv_buffer_ghost_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout | pointer | Receive buffer of ghost cells. |
| `q_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout |  | Field component to be updated. |

**Call graph**

```mermaid
flowchart TD
  update_ghost_mpi_gpu["update_ghost_mpi_gpu"] --> receive_recv_buffer_ghost_gpu_dev["receive_recv_buffer_ghost_gpu_dev"]
  style receive_recv_buffer_ghost_gpu_dev fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### update_ghost_local_gpu_dev

Update (local) ghost cells.

```fortran
subroutine update_ghost_local_gpu_dev(ngc, l_map_ghost_cell_gpu, q_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Ghost cells number. |
| `l_map_ghost_cell_gpu` | integer(kind=[I8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | pointer | Local map of ghost cells. |
| `q_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout |  | Field component to be updated. |

**Call graph**

```mermaid
flowchart TD
  update_ghost_local_gpu["update_ghost_local_gpu"] --> update_ghost_local_gpu_dev["update_ghost_local_gpu_dev"]
  style update_ghost_local_gpu_dev fill:#3e63dd,stroke:#99b,stroke-width:2px
```
