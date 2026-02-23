---
title: adam_field_gmp_kernels
---

# adam_field_gmp_kernels

> ADAM, field class GMP kernels (GMP backend of [field_object](/api/src/lib/common/adam_field_object#field-object)).

**Source**: `src/lib/gmp/adam_field_gmp_kernels.F90`

**Dependencies**

```mermaid
graph LR
  adam_field_gmp_kernels["adam_field_gmp_kernels"] --> penf["penf"]
```

## Contents

- [compute_q_gradient_gmp](#compute-q-gradient-gmp)
- [compute_normL2_residuals_gmp](#compute-norml2-residuals-gmp)
- [copy_transpose_gpu_cpu_gmp](#copy-transpose-gpu-cpu-gmp)
- [populate_send_buffer_ghost_gmp](#populate-send-buffer-ghost-gmp)
- [receive_recv_buffer_ghost_gmp](#receive-recv-buffer-ghost-gmp)
- [update_ghost_local_gmp](#update-ghost-local-gmp)

## Subroutines

### compute_q_gradient_gmp

Compute gradient of q(ivar).

```fortran
subroutine compute_q_gradient_gmp(b, ni, nj, nk, ngc, dx, dy, dz, q_gpu, ivar, gradient)
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
  compute_q_gradient["compute_q_gradient"] --> compute_q_gradient_gmp["compute_q_gradient_gmp"]
  style compute_q_gradient_gmp fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### compute_normL2_residuals_gmp

Compute L2 norm of residuals.

```fortran
subroutine compute_normL2_residuals_gmp(ni, nj, nk, ngc, nv, blocks_number, dq_gpu, norm)
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
  save_residuals["save_residuals"] --> compute_normL2_residuals_gmp["compute_normL2_residuals_gmp"]
  style compute_normL2_residuals_gmp fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### copy_transpose_gpu_cpu_gmp

Copy transposed data from GPU to CPU by CUF threads.

```fortran
subroutine copy_transpose_gpu_cpu_gmp(ni, nj, nk, ngc, nv, blocks_number, q_gpu, q_t_gpu, q_cpu)
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
| `q_cpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | out |  | Conservative variables on CPU. |

**Call graph**

```mermaid
flowchart TD
  copy_transpose_gpu_cpu["copy_transpose_gpu_cpu"] --> copy_transpose_gpu_cpu_gmp["copy_transpose_gpu_cpu_gmp"]
  style copy_transpose_gpu_cpu_gmp fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### populate_send_buffer_ghost_gmp

Polulate send buffer ghost GPU.

```fortran
subroutine populate_send_buffer_ghost_gmp(ngc, comm_map_send_ghost_cell_gpu, send_buffer_ghost_gpu, q_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Ghost cells number. |
| `comm_map_send_ghost_cell_gpu` | integer(kind=[I8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | pointer | Comm map, cell information. |
| `send_buffer_ghost_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout | pointer | Send buffer of ghost cells. |
| `q_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Field component to be updated. |

**Call graph**

```mermaid
flowchart TD
  update_ghost_mpi["update_ghost_mpi"] --> populate_send_buffer_ghost_gmp["populate_send_buffer_ghost_gmp"]
  style populate_send_buffer_ghost_gmp fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### receive_recv_buffer_ghost_gmp

Receive recv buffer ghost GPU.

```fortran
subroutine receive_recv_buffer_ghost_gmp(ngc, comm_map_recv_ghost_cell_gpu, recv_buffer_ghost_gpu, q_gpu)
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
  update_ghost_mpi["update_ghost_mpi"] --> receive_recv_buffer_ghost_gmp["receive_recv_buffer_ghost_gmp"]
  style receive_recv_buffer_ghost_gmp fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### update_ghost_local_gmp

Update (local) ghost cells.

```fortran
subroutine update_ghost_local_gmp(ngc, local_map_ghost_cell_gpu, q_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Ghost cells number. |
| `local_map_ghost_cell_gpu` | integer(kind=[I8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | pointer | Local map of ghost cells. |
| `q_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout |  | Field component to be updated. |

**Call graph**

```mermaid
flowchart TD
  update_ghost_local["update_ghost_local"] --> update_ghost_local_gmp["update_ghost_local_gmp"]
  style update_ghost_local_gmp fill:#3e63dd,stroke:#99b,stroke-width:2px
```
