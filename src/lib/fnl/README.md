<a name="top"></a>

# ADAM library FNL sources

> This subdirectory contains the ADAM library kernels of FUNDAL-based GPU-offloading backend.

All the FNL classes and libraries as well ass the *common* ones are exposed by the main **entry module** `adam_fnl_library.F90`, i.e.:

```fortran
! objects (classes)
public :: field_fnl_object
public :: ib_fnl_object
public :: maps_fnl_object
public :: mpih_fnl_object
public :: rk_fnl_object
public :: weno_fnl_object

! device kernels
public :: compute_q_gradient_dev
public :: compute_normL2_residuals_dev
public :: copy_transpose_gpu_cpu_dev
public :: populate_send_buffer_ghost_gpu_dev
public :: receive_recv_buffer_ghost_gpu_dev
public :: update_ghost_local_gpu_dev
public :: compute_eikonal_dq_phi_dev
public :: compute_phi_all_solids_dev
public :: compute_phi_analytical_sphere_dev
public :: evolve_eikonal_q_phi_dev
public :: invert_eikonal_q_phi_dev
public :: move_phi_dev
public :: reduce_cell_order_phi_dev
public :: rk_assign_stage_dev
public :: rk_compute_stage_dev
public :: rk_compute_stage_ls_dev
public :: rk_initialize_stages_dev
public :: rk_update_q_dev
public :: weno_reconstruct_upwind_dev

! global parameters
public :: S_max
public :: S_max_m1
```

The procedures with `_dev` suffix contain device calculations either by means of device-parallel loops or by device-offloading the entire procedure (pseudo kernel).

The main documentation of these sources is contained in the following sections:

| [Copyrights](#copyrights) | [API Documentation](#api-documentation) |

Go to [Top](#top)

# Copyrights

ADAM is currently a closed project:

> Copyright (C) Di Mascio/Rossi/Salvadore/Zaghi, Inc - All Rights Reserved.
>
> Unauthorized copying of these source files, via any medium is strictly prohibited, proprietary and confidential.
> Written by Andrea di Mascio, Giacomo Rossi, Francesco Salvadore and Stefano Zaghi, September 2023.

Future versions could be released with a more Free Open Source Software (FOSS) licence.

Go to [Top](#top)

# API Documentation

Currently, the following sources compose the subdirectory:

+ `adam_fnl_field_kernels.F90` is the **field kernels** library for FNL acceleration, see [field FNL kernels API](https://szaghi.github.io/adam/module/adam_fnl_field_kernels.html) for more details;
+ `adam_fnl_field_object.F90` is the **field** class for FNL acceleration, see [field FNL object API](https://szaghi.github.io/adam/type/adam_fnl_field_object.html) for more details;
+ `adam_fnl_library.F90` is the **FNL entry module** exposing all modules for FNL acceleration, see [FNL library API](https://szaghi.github.io/adam/module/adam_fnl_library.html) for more details;
+ `adam_fnl_ib_kernels.F90` is the **IB kernels** library for FNL acceleration, see [IB FNL kernels API](https://szaghi.github.io/adam/module/adam_fnl_ib_kernels.html) for more details;
+ `adam_fnl_ib_object.F90` is the **IB** class for FNL acceleration, see [IB FNL object API](https://szaghi.github.io/adam/type/adam_fnl_ib_object.html) for more details;
+ `adam_fnl_maps_object.F90` is the **maps** class for FNL acceleration, see [maps FNL object API](https://szaghi.github.io/adam/type/adam_fnl_maps_object.html) for more details;
+ `adam_fnl_mpih_object.F90` is the **MPI handler** class for FNL acceleration, see [MPIH FNL object API](https://szaghi.github.io/adam/type/adam_fnl_mpih_object.html) for more details;
+ `adam_fnl_rk_kernels.F90` is the **RK kernels** library for FNL acceleration, see [RK FNL kernels API](https://szaghi.github.io/adam/module/adam_fnl_rk_kernels.html) for more details;
+ `adam_fnl_rk_object.F90` is the **RK** class for FNL acceleration, see [RK FNL object API](https://szaghi.github.io/adam/type/adam_fnl_rk_object.html) for more details;
+ `adam_fnl_weno_kernels.F90` is the **WENO kernels** library for FNL acceleration, see [WENO FNL kernels API](https://szaghi.github.io/adam/module/adam_fnl_weno_kernels.html) for more details;
+ `adam_fnl_weno_object.F90` is the **WENO** class for FNL acceleration, see [WENO FNL object API](https://szaghi.github.io/adam/type/adam_fnl_weno_object.html) for more details;

Go to [Top](#top)

