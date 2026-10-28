   ! Model-agnostic FNL face-flux kernels (flux splitting): included in the `contains` part of every per-model kernel
   ! module `adam_flume_fnl_<model>_kernels` (issue #41, section 4), which provides, by use association:
   ! `NV_K`, `NV_AUX_K` (the model's conservative and auxiliary variables numbers, compile-time constants that size
   ! every private array: no model branching and no maximum-model buffers in the kernels), `S_MAX`, the model's
   ! split adapter `face_split_fluxes(gamma, ch, ...)` (one argument list for every model: the GLM cleaning speed `ch` is
   ! ignored by the models without cleaning), `compute_face_flux_back_projection`, and `weno_reconstruct_upwind_dev`.
   ! The host selects the model's module by `select case`, outside the kernels.
   !
   ! One kernel per stencil width. Each public routine is a host dispatcher: it selects, by the WENO half-width S, a
   ! kernel compiled with S a named constant (`_s1_dev` ... `_s5_dev`, WENO-1 to WENO-9, S_MAX = 5), or the generic
   ! kernel with S a runtime argument (`_sn_dev`, for a width without a specialized kernel). Every kernel compiles the
   ! same body (`adam_flume_fnl_face_fluxes_kernel_body.INC`, `adam_flume_fnl_gravity_edge_kernel_body.INC`). With S a
   ! constant the stencil loops have constant bounds, and the device link-time optimisation (-gpu=lto) inlines the
   ! device routines with that constant: on one A30 the face kernel is 33% (WENO-9) to 65% (WENO-3) faster than the
   ! generic one, 56 ms instead of 103 ms per launch at WENO-5 (src/tests/flume/speedup/profile-couette3d/README.md,
   ! "One face-flux kernel per stencil width").

   subroutine compute_face_fluxes_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, S, gamma, ch, is_characteristic, &
                                      weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, weno_sigma, q_gpu, q_aux_gpu,    &
                                      fl_gpu)
   !< Compute the WENO face fluxes of direction `d`: face `(i,j,k)` lies between cells `(i,j,k)` and
   !< `(i+di,j+dj,k+dk)`, so the face array starts at index 0 along `d`.
   !<
   !< Per face (device twin of the CPU `compute_face_fluxes`): gather the stencil `m = 1-S ... S` into constant-bound
   !< privates, project and split it, reconstruct every field (the private `v` is the packed stencil of the WENO
   !< primitive), back-project. Host dispatcher: the kernel of the width S (see the top of this file).
   integer(I4P), intent(in)    :: d                                       !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)    :: di, dj, dk                              !< Unit step along `d`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                         !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                           !< Actual blocks number.
   integer(I4P), intent(in)    :: S                                       !< WENO stencil half-width.
   real(R8P),    intent(in)    :: gamma                                   !< Specific heats ratio.
   real(R8P),    intent(in)    :: ch                                      !< GLM cleaning speed (models with GLM only).
   logical,      intent(in)    :: is_characteristic                       !< Characteristic (or conservative) variables.
   real(R8P),    intent(in)    :: weno_a_gpu(1:,0:,1:)                    !< WENO optimal weights.
   real(R8P),    intent(in)    :: weno_p_gpu(1:,0:,0:,1:)                 !< WENO polynomials coefficients.
   real(R8P),    intent(in)    :: weno_d_gpu(0:,0:,0:,1:)                 !< WENO smoothness indicators coefficients.
   real(R8P),    intent(in)    :: weno_zeps                               !< WENO parameter avoiding division by zero.
   real(R8P),    intent(in)    :: weno_sigma                              !< WENO descaler switch (0 Jiang-Shu, 1 SI).
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)       !< Conservative variables.
   real(R8P),    intent(in)    :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Auxiliary variables.
   real(R8P),    intent(inout) :: fl_gpu(1:,1-di:,1-dj:,1-dk:,1:)         !< Face fluxes of direction `d`.

   select case(S)
   case(1) ; call compute_face_fluxes_s1_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch,                &
                                             is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps,        &
                                             weno_sigma, q_gpu, q_aux_gpu, fl_gpu)
   case(2) ; call compute_face_fluxes_s2_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch,                &
                                             is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps,        &
                                             weno_sigma, q_gpu, q_aux_gpu, fl_gpu)
   case(3) ; call compute_face_fluxes_s3_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch,                &
                                             is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps,        &
                                             weno_sigma, q_gpu, q_aux_gpu, fl_gpu)
   case(4) ; call compute_face_fluxes_s4_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch,                &
                                             is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps,        &
                                             weno_sigma, q_gpu, q_aux_gpu, fl_gpu)
   case(5) ; call compute_face_fluxes_s5_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch,                &
                                             is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps,        &
                                             weno_sigma, q_gpu, q_aux_gpu, fl_gpu)
   case default ; call compute_face_fluxes_sn_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, S, gamma, ch,        &
                                                  is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps,   &
                                                  weno_sigma, q_gpu, q_aux_gpu, fl_gpu)
   endselect
   endsubroutine compute_face_fluxes_dev

   subroutine compute_face_fluxes_s1_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                         weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, weno_sigma, q_gpu, q_aux_gpu, &
                                         fl_gpu)
   !< `compute_face_fluxes_dev` for WENO-1: the stencil half-width S = 1 is a compile-time constant.
   use :: adam_fnl_weno_width_kernels, only : weno_reconstruct_upwind_s1_dev
   integer(I4P), parameter :: S = 1_I4P !< WENO stencil half-width.
#define FNL_WENO_WIDTH weno_reconstruct_upwind_s1_dev
#include "adam_flume_fnl_face_fluxes_kernel_body.INC"
#undef FNL_WENO_WIDTH
   endsubroutine compute_face_fluxes_s1_dev

   subroutine compute_face_fluxes_s2_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                         weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, weno_sigma, q_gpu, q_aux_gpu, &
                                         fl_gpu)
   !< `compute_face_fluxes_dev` for WENO-3: the stencil half-width S = 2 is a compile-time constant.
   use :: adam_fnl_weno_width_kernels, only : weno_reconstruct_upwind_s2_dev
   integer(I4P), parameter :: S = 2_I4P !< WENO stencil half-width.
#define FNL_WENO_WIDTH weno_reconstruct_upwind_s2_dev
#include "adam_flume_fnl_face_fluxes_kernel_body.INC"
#undef FNL_WENO_WIDTH
   endsubroutine compute_face_fluxes_s2_dev

   subroutine compute_face_fluxes_s3_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                         weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, weno_sigma, q_gpu, q_aux_gpu, &
                                         fl_gpu)
   !< `compute_face_fluxes_dev` for WENO-5: the stencil half-width S = 3 is a compile-time constant.
   use :: adam_fnl_weno_width_kernels, only : weno_reconstruct_upwind_s3_dev
   integer(I4P), parameter :: S = 3_I4P !< WENO stencil half-width.
#define FNL_WENO_WIDTH weno_reconstruct_upwind_s3_dev
#include "adam_flume_fnl_face_fluxes_kernel_body.INC"
#undef FNL_WENO_WIDTH
   endsubroutine compute_face_fluxes_s3_dev

   subroutine compute_face_fluxes_s4_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                         weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, weno_sigma, q_gpu, q_aux_gpu, &
                                         fl_gpu)
   !< `compute_face_fluxes_dev` for WENO-7: the stencil half-width S = 4 is a compile-time constant.
   use :: adam_fnl_weno_width_kernels, only : weno_reconstruct_upwind_s4_dev
   integer(I4P), parameter :: S = 4_I4P !< WENO stencil half-width.
#define FNL_WENO_WIDTH weno_reconstruct_upwind_s4_dev
#include "adam_flume_fnl_face_fluxes_kernel_body.INC"
#undef FNL_WENO_WIDTH
   endsubroutine compute_face_fluxes_s4_dev

   subroutine compute_face_fluxes_s5_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                         weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, weno_sigma, q_gpu, q_aux_gpu, &
                                         fl_gpu)
   !< `compute_face_fluxes_dev` for WENO-9: the stencil half-width S = 5 is a compile-time constant.
   use :: adam_fnl_weno_width_kernels, only : weno_reconstruct_upwind_s5_dev
   integer(I4P), parameter :: S = 5_I4P !< WENO stencil half-width.
#define FNL_WENO_WIDTH weno_reconstruct_upwind_s5_dev
#include "adam_flume_fnl_face_fluxes_kernel_body.INC"
#undef FNL_WENO_WIDTH
   endsubroutine compute_face_fluxes_s5_dev

   subroutine compute_face_fluxes_sn_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, S, gamma, ch,                 &
                                         is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps,            &
                                         weno_sigma, q_gpu, q_aux_gpu, fl_gpu)
   !< `compute_face_fluxes_dev` for a width without a specialized kernel: S a runtime argument.
   integer(I4P), intent(in)    :: S                                       !< WENO stencil half-width.
#include "adam_flume_fnl_face_fluxes_kernel_body.INC"
   endsubroutine compute_face_fluxes_sn_dev

   subroutine add_gravity_edge_correction_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, S, gamma, ch,            &
                                              is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps,       &
                                              weno_sigma, ir, gd, q_gpu, q_aux_gpu, gc_gpu)
   !< Add the faces beyond the block to the gravity energy correction of direction `d` (issue #73, M7-G2a; device twin
   !< of the CPU `add_gravity_edge_correction`): cell 1 takes `-g_d F^rho / 24` of face -1, cell `n_d` that of face
   !< `n_d + 1`, each face computed as `compute_face_fluxes_dev` does (the same body, its mass component kept). Host
   !< dispatcher: the kernel of the width S (see the top of this file).
   integer(I4P), intent(in)    :: d                                       !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)    :: di, dj, dk                              !< Unit step along `d`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                         !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                           !< Actual blocks number.
   integer(I4P), intent(in)    :: S                                       !< WENO stencil half-width.
   real(R8P),    intent(in)    :: gamma                                   !< Specific heats ratio.
   real(R8P),    intent(in)    :: ch                                      !< GLM cleaning speed (models with GLM only).
   logical,      intent(in)    :: is_characteristic                       !< Characteristic (or conservative) variables.
   real(R8P),    intent(in)    :: weno_a_gpu(1:,0:,1:)                    !< WENO optimal weights.
   real(R8P),    intent(in)    :: weno_p_gpu(1:,0:,0:,1:)                 !< WENO polynomials coefficients.
   real(R8P),    intent(in)    :: weno_d_gpu(0:,0:,0:,1:)                 !< WENO smoothness indicators coefficients.
   real(R8P),    intent(in)    :: weno_zeps                               !< WENO parameter avoiding division by zero.
   real(R8P),    intent(in)    :: weno_sigma                              !< WENO descaler switch (0 Jiang-Shu, 1 SI).
   integer(I4P), intent(in)    :: ir                                      !< Index of the density (mass) variable.
   real(R8P),    intent(in)    :: gd                                      !< Gravity along `d`.
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)       !< Conservative variables.
   real(R8P),    intent(in)    :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Auxiliary variables.
   real(R8P),    intent(inout) :: gc_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Gravity energy correction [..., 1].

   select case(S)
   case(1) ; call add_gravity_edge_correction_s1_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch,        &
                                                     is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu,           &
                                                     weno_zeps, weno_sigma, ir, gd, q_gpu, q_aux_gpu, gc_gpu)
   case(2) ; call add_gravity_edge_correction_s2_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch,        &
                                                     is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu,           &
                                                     weno_zeps, weno_sigma, ir, gd, q_gpu, q_aux_gpu, gc_gpu)
   case(3) ; call add_gravity_edge_correction_s3_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch,        &
                                                     is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu,           &
                                                     weno_zeps, weno_sigma, ir, gd, q_gpu, q_aux_gpu, gc_gpu)
   case(4) ; call add_gravity_edge_correction_s4_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch,        &
                                                     is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu,           &
                                                     weno_zeps, weno_sigma, ir, gd, q_gpu, q_aux_gpu, gc_gpu)
   case(5) ; call add_gravity_edge_correction_s5_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch,        &
                                                     is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu,           &
                                                     weno_zeps, weno_sigma, ir, gd, q_gpu, q_aux_gpu, gc_gpu)
   case default ; call add_gravity_edge_correction_sn_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, S, gamma,    &
                                                          ch, is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu,  &
                                                          weno_zeps, weno_sigma, ir, gd, q_gpu, q_aux_gpu, gc_gpu)
   endselect
   endsubroutine add_gravity_edge_correction_dev

   subroutine add_gravity_edge_correction_s1_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch,            &
                                                 is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps,    &
                                                 weno_sigma, ir, gd, q_gpu, q_aux_gpu, gc_gpu)
   !< `add_gravity_edge_correction_dev` for WENO-1: the stencil half-width S = 1 is a compile-time constant.
   use :: adam_fnl_weno_width_kernels, only : weno_reconstruct_upwind_s1_dev
   integer(I4P), parameter :: S = 1_I4P !< WENO stencil half-width.
#define FNL_WENO_WIDTH weno_reconstruct_upwind_s1_dev
#include "adam_flume_fnl_gravity_edge_kernel_body.INC"
#undef FNL_WENO_WIDTH
   endsubroutine add_gravity_edge_correction_s1_dev

   subroutine add_gravity_edge_correction_s2_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch,            &
                                                 is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps,    &
                                                 weno_sigma, ir, gd, q_gpu, q_aux_gpu, gc_gpu)
   !< `add_gravity_edge_correction_dev` for WENO-3: the stencil half-width S = 2 is a compile-time constant.
   use :: adam_fnl_weno_width_kernels, only : weno_reconstruct_upwind_s2_dev
   integer(I4P), parameter :: S = 2_I4P !< WENO stencil half-width.
#define FNL_WENO_WIDTH weno_reconstruct_upwind_s2_dev
#include "adam_flume_fnl_gravity_edge_kernel_body.INC"
#undef FNL_WENO_WIDTH
   endsubroutine add_gravity_edge_correction_s2_dev

   subroutine add_gravity_edge_correction_s3_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch,            &
                                                 is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps,    &
                                                 weno_sigma, ir, gd, q_gpu, q_aux_gpu, gc_gpu)
   !< `add_gravity_edge_correction_dev` for WENO-5: the stencil half-width S = 3 is a compile-time constant.
   use :: adam_fnl_weno_width_kernels, only : weno_reconstruct_upwind_s3_dev
   integer(I4P), parameter :: S = 3_I4P !< WENO stencil half-width.
#define FNL_WENO_WIDTH weno_reconstruct_upwind_s3_dev
#include "adam_flume_fnl_gravity_edge_kernel_body.INC"
#undef FNL_WENO_WIDTH
   endsubroutine add_gravity_edge_correction_s3_dev

   subroutine add_gravity_edge_correction_s4_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch,            &
                                                 is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps,    &
                                                 weno_sigma, ir, gd, q_gpu, q_aux_gpu, gc_gpu)
   !< `add_gravity_edge_correction_dev` for WENO-7: the stencil half-width S = 4 is a compile-time constant.
   use :: adam_fnl_weno_width_kernels, only : weno_reconstruct_upwind_s4_dev
   integer(I4P), parameter :: S = 4_I4P !< WENO stencil half-width.
#define FNL_WENO_WIDTH weno_reconstruct_upwind_s4_dev
#include "adam_flume_fnl_gravity_edge_kernel_body.INC"
#undef FNL_WENO_WIDTH
   endsubroutine add_gravity_edge_correction_s4_dev

   subroutine add_gravity_edge_correction_s5_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch,            &
                                                 is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps,    &
                                                 weno_sigma, ir, gd, q_gpu, q_aux_gpu, gc_gpu)
   !< `add_gravity_edge_correction_dev` for WENO-9: the stencil half-width S = 5 is a compile-time constant.
   use :: adam_fnl_weno_width_kernels, only : weno_reconstruct_upwind_s5_dev
   integer(I4P), parameter :: S = 5_I4P !< WENO stencil half-width.
#define FNL_WENO_WIDTH weno_reconstruct_upwind_s5_dev
#include "adam_flume_fnl_gravity_edge_kernel_body.INC"
#undef FNL_WENO_WIDTH
   endsubroutine add_gravity_edge_correction_s5_dev

   subroutine add_gravity_edge_correction_sn_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, S, gamma, ch,         &
                                                 is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps,    &
                                                 weno_sigma, ir, gd, q_gpu, q_aux_gpu, gc_gpu)
   !< `add_gravity_edge_correction_dev` for a width without a specialized kernel: S a runtime argument.
   integer(I4P), intent(in)    :: S                                       !< WENO stencil half-width.
#include "adam_flume_fnl_gravity_edge_kernel_body.INC"
   endsubroutine add_gravity_edge_correction_sn_dev
