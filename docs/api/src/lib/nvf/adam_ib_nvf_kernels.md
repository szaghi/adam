---
title: adam_ib_nvf_kernels
---

# adam_ib_nvf_kernels

> ADAM, IB class NVF kernels (NVF backend of [ib_object](/api/src/lib/common/adam_ib_object#ib-object)).

**Source**: `src/lib/nvf/adam_ib_nvf_kernels.F90`

**Dependencies**

```mermaid
graph LR
  adam_ib_nvf_kernels["adam_ib_nvf_kernels"] --> penf["penf"]
```

## Contents

- [compute_eikonal_dq_phi_cuf](#compute-eikonal-dq-phi-cuf)
- [compute_phi_all_solids_cuf](#compute-phi-all-solids-cuf)
- [compute_phi_analytical_sphere_cuf](#compute-phi-analytical-sphere-cuf)
- [evolve_eikonal_q_phi_cuf](#evolve-eikonal-q-phi-cuf)
- [invert_eikonal_q_phi_cuf](#invert-eikonal-q-phi-cuf)
- [move_phi_cuf](#move-phi-cuf)
- [reduce_cell_order_phi_cuf](#reduce-cell-order-phi-cuf)

## Subroutines

### compute_eikonal_dq_phi_cuf

Compute eikonal dq inside IB.

```fortran
subroutine compute_eikonal_dq_phi_cuf(ib, ni, nj, nk, ngc, nv, blocks_number, dx_gpu, dy_gpu, dz_gpu, phi_gpu, q_gpu, dq_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `ib` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | IB solid index. |
| `ni` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in I direction. |
| `nj` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in J direction. |
| `nk` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in K direction. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Ghost cells number. |
| `nv` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of conservative varibales. |
| `blocks_number` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of blocks. |
| `dx_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device | X space steps. |
| `dy_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device | Y space steps. |
| `dz_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device | Z space steps. |
| `phi_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device | Distance function. |
| `q_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device | State variables. |
| `dq_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout | device | State variables variations. |

**Call graph**

```mermaid
flowchart TD
  evolve_eikonal["evolve_eikonal"] --> compute_eikonal_dq_phi_cuf["compute_eikonal_dq_phi_cuf"]
  style compute_eikonal_dq_phi_cuf fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### compute_phi_all_solids_cuf

Compute last phi index, all solids summary.

```fortran
subroutine compute_phi_all_solids_cuf(ni, nj, nk, ngc, blocks_number, phi_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `ni` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in I direction. |
| `nj` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in J direction. |
| `nk` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in K direction. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Ghost grid number. |
| `blocks_number` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of blocks. |
| `phi_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout | device | Distance function. |

**Call graph**

```mermaid
flowchart TD
  compute_phi["compute_phi"] --> compute_phi_all_solids_cuf["compute_phi_all_solids_cuf"]
  style compute_phi_all_solids_cuf fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### compute_phi_analytical_sphere_cuf

Compute phi, distance from ib solid, analytical sphere solid.

```fortran
subroutine compute_phi_analytical_sphere_cuf(ib, ni, nj, nk, ngc, blocks_number, sphere, x_cell_gpu, y_cell_gpu, z_cell_gpu, phi_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `ib` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | IB solid index. |
| `ni` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in I direction. |
| `nj` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in J direction. |
| `nk` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in K direction. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Ghost grid number. |
| `blocks_number` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of blocks. |
| `sphere` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Sphere center [1:3] and radius [4]. |
| `x_cell_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device | Cells x coordinates on GPU. |
| `y_cell_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device | Cells y coordinates on GPU. |
| `z_cell_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device | Cells z coordinates on GPU. |
| `phi_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout | device | Distance function. |

**Call graph**

```mermaid
flowchart TD
  compute_phi["compute_phi"] --> compute_phi_analytical_sphere_cuf["compute_phi_analytical_sphere_cuf"]
  style compute_phi_analytical_sphere_cuf fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### evolve_eikonal_q_phi_cuf

Evolve eikonal equation over q inside IB.

```fortran
subroutine evolve_eikonal_q_phi_cuf(ib, ni, nj, nk, ngc, nv, blocks_number, phi_gpu, dq_gpu, q_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `ib` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | IB solid index. |
| `ni` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in I direction. |
| `nj` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in J direction. |
| `nk` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in K direction. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Ghost cells number. |
| `nv` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of conservative varibales. |
| `blocks_number` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of blocks. |
| `phi_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device | Distance function. |
| `dq_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device | State variables variation. |
| `q_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout | device | State variables. |

**Call graph**

```mermaid
flowchart TD
  evolve_eikonal["evolve_eikonal"] --> evolve_eikonal_q_phi_cuf["evolve_eikonal_q_phi_cuf"]
  style evolve_eikonal_q_phi_cuf fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### invert_eikonal_q_phi_cuf

Invert eikonal equation over q inside IB.

```fortran
subroutine invert_eikonal_q_phi_cuf(BCS_VISCOUS, BCS_EULER, ib, ni, nj, nk, ngc, nv, blocks_number, bcs_type, phi_gpu, q_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `BCS_VISCOUS` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Viscous wall BCS parameter. |
| `BCS_EULER` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Euler wall BCS parameter. |
| `ib` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | IB solid index. |
| `ni` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in I direction. |
| `nj` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in J direction. |
| `nk` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in K direction. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Ghost cells number. |
| `nv` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of conservative varibales. |
| `blocks_number` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of blocks. |
| `bcs_type` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Immersed boundary type. |
| `phi_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device | Distance field. |
| `q_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout | device | Conservative field. |

**Call graph**

```mermaid
flowchart TD
  invert_eikonal["invert_eikonal"] --> invert_eikonal_q_phi_cuf["invert_eikonal_q_phi_cuf"]
  style invert_eikonal_q_phi_cuf fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### move_phi_cuf

Move phi and the actual ptree representation.

```fortran
subroutine move_phi_cuf(ni, nj, nk, ngc, blocks_number, velocity, phi_gpu, dphi_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `ni` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in I direction. |
| `nj` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in J direction. |
| `nk` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in K direction. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Ghost grid number. |
| `blocks_number` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of blocks. |
| `velocity` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Velocity of the movement. |
| `phi_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout | device | Distance function. |
| `dphi_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout | device | Distance function gradient. |

**Call graph**

```mermaid
flowchart TD
  move_phi["move_phi"] --> move_phi_cuf["move_phi_cuf"]
  style move_phi_cuf fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### reduce_cell_order_phi_cuf

Reduce local-cell order of spatial operator close to solids.

```fortran
subroutine reduce_cell_order_phi_cuf(ib, ni, nj, nk, ngc, blocks_number, iweno, ib_reduced_order, ib_reduction_extent, phi_gpu, cell_scheme_gpu)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `ib` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | IB solid index. |
| `ni` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in I direction. |
| `nj` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in J direction. |
| `nk` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Grid cells number in K direction. |
| `ngc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Ghost grid number. |
| `blocks_number` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Number of blocks. |
| `iweno` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | WENO (half) stencil lenght. |
| `ib_reduced_order` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Reduced order close to IB solids. |
| `ib_reduction_extent` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Extent of reduction close to IBi solids. |
| `phi_gpu` | real(kind=[R8P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | device | Distance function. |
| `cell_scheme_gpu` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout | device | Modified order close to solids. |

**Call graph**

```mermaid
flowchart TD
  compute_phi["compute_phi"] --> reduce_cell_order_phi_cuf["reduce_cell_order_phi_cuf"]
  style reduce_cell_order_phi_cuf fill:#3e63dd,stroke:#99b,stroke-width:2px
```
