   ! Model- and solver-agnostic FNL face-flux kernel of the scheme `weno-riemann` (issue #47): the body shared by every
   ! per-(model, solver) device kernel module, device twin of `adam_flume_cpu_riemann_face_kernels_agnostic.INC`.
   !
   ! Included in the `contains` part of `adam_flume_fnl_<model>_<solver>_kernels`, which provides, by use association:
   ! `NV_K`, `NV_AUX_K`, `S_MAX`, `weno_reconstruct_upwind_wratio_dev`, and the four device adapters
   ! `face_interpolation_fields`, `face_states`, `riemann_flux`, `cell_flux` (see the CPU twin). The host selects the
   ! module by `select case` on the model and the solver, outside the kernels.
   !
   ! One kernel per stencil width, as for the flux-splitting kernels (`adam_flume_fnl_face_kernels_agnostic.INC`): each
   ! public routine is a host dispatcher that selects, by the WENO half-width S, the kernel compiled with S a named
   ! constant (`_s1_dev` ... `_s5_dev`, WENO-1 to WENO-9, S_MAX = 5), or the generic kernel with S a run-time argument
   ! (`_sn_dev`, for a width without a specialized kernel). Every kernel compiles the same body
   ! (`adam_flume_fnl_riemann_face_fluxes_kernel_body.INC`, `adam_flume_fnl_riemann_gravity_edge_kernel_body.INC`) and the
   ! kernels of a width call the WENO reconstruction of that width (`adam_fnl_weno_width_kernels`), the stencil loops and the
   ! tables' stencil index being constants in it.

   subroutine compute_riemann_face_fluxes_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, S, gamma, ch,        &
                                              is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, &
                                              weno_sigma, cc, tau, q_gpu, q_aux_gpu, fl_gpu, fallbacks)
   !< Compute the `weno-riemann` face fluxes of direction `d`: face `(i,j,k)` lies between cells `(i,j,k)` and
   !< `(i+di,j+dj,k+dk)`, so the face array starts at index 0 along `d`.
   !<
   !< Per face (device twin of the CPU `compute_riemann_face_fluxes`): gather the stencil `m = min(1-S,-1) ... max(S,2)`
   !< into constant-bound privates, WENO-interpolate its fields to the two face states (`weno_a_gpu`, `weno_p_gpu` are
   !< the INTERPOLATION tables), take the Riemann flux `F` and add the correction
   !< `F^ = F + s (cc(1) F + cc(2) (f_0 + f_1) + cc(3) (f_-1 + f_2) - F)`, `s = 1` where every WENO weights ratio is
   !< `>= tau`, else 0; `fallbacks` counts the faces whose Riemann solver fell back (HLLD to HLL).
   integer(I4P), intent(in)    :: d                                       !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)    :: di, dj, dk                              !< Unit step along `d`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                         !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                           !< Actual blocks number.
   integer(I4P), intent(in)    :: S                                       !< WENO stencil half-width.
   real(R8P),    intent(in)    :: gamma                                   !< Specific heats ratio.
   real(R8P),    intent(in)    :: ch                                      !< GLM cleaning speed (models with GLM only).
   logical,      intent(in)    :: is_characteristic                       !< Characteristic (or primitive) variables.
   real(R8P),    intent(in)    :: weno_a_gpu(1:,0:,1:)                    !< WENO interpolation optimal weights.
   real(R8P),    intent(in)    :: weno_p_gpu(1:,0:,0:,1:)                 !< WENO interpolation polynomials coefficients.
   real(R8P),    intent(in)    :: weno_d_gpu(0:,0:,0:,1:)                 !< WENO smoothness indicators coefficients.
   real(R8P),    intent(in)    :: weno_zeps                               !< WENO parameter avoiding division by zero.
   real(R8P),    intent(in)    :: weno_sigma                              !< WENO descaler switch (0 Jiang-Shu, 1 SI).
   real(R8P),    intent(in)    :: cc(3)                                   !< Correction coefficients.
   real(R8P),    intent(in)    :: tau                                     !< Sensor threshold.
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)       !< Conservative variables.
   real(R8P),    intent(in)    :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Auxiliary variables.
   real(R8P),    intent(inout) :: fl_gpu(1:,1-di:,1-dj:,1-dk:,1:)         !< Face fluxes of direction `d`.
   integer(I4P), intent(out)   :: fallbacks                               !< Faces whose Riemann solver fell back.

   select case(S)
   case(1) ; call compute_riemann_face_fluxes_s1_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                                    is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, &
                                                    weno_sigma, cc, tau, q_gpu, q_aux_gpu, fl_gpu, fallbacks)
   case(2) ; call compute_riemann_face_fluxes_s2_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                                    is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, &
                                                    weno_sigma, cc, tau, q_gpu, q_aux_gpu, fl_gpu, fallbacks)
   case(3) ; call compute_riemann_face_fluxes_s3_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                                    is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, &
                                                    weno_sigma, cc, tau, q_gpu, q_aux_gpu, fl_gpu, fallbacks)
   case(4) ; call compute_riemann_face_fluxes_s4_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                                    is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, &
                                                    weno_sigma, cc, tau, q_gpu, q_aux_gpu, fl_gpu, fallbacks)
   case(5) ; call compute_riemann_face_fluxes_s5_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                                    is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, &
                                                    weno_sigma, cc, tau, q_gpu, q_aux_gpu, fl_gpu, fallbacks)
   case default ; call compute_riemann_face_fluxes_sn_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, S, gamma, &
                                                    ch, is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, &
                                                    weno_zeps, weno_sigma, cc, tau, q_gpu, q_aux_gpu, fl_gpu, &
                                                    fallbacks)
   endselect
   endsubroutine compute_riemann_face_fluxes_dev

   subroutine compute_riemann_face_fluxes_s1_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                         is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, &
                                         weno_sigma, cc, tau, q_gpu, q_aux_gpu, fl_gpu, fallbacks)
   !< `compute_riemann_face_fluxes_dev` for WENO-1: the stencil half-width S = 1 is a compile-time constant.
   use :: adam_fnl_weno_width_kernels, only : weno_width => weno_reconstruct_upwind_wratio_s1_dev
   integer(I4P), parameter :: S = 1_I4P !< WENO stencil half-width.
#define FNL_WENO_WRATIO_WIDTH
#include "adam_flume_fnl_riemann_face_fluxes_kernel_body.INC"
#undef FNL_WENO_WRATIO_WIDTH
   endsubroutine compute_riemann_face_fluxes_s1_dev

   subroutine compute_riemann_face_fluxes_s2_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                         is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, &
                                         weno_sigma, cc, tau, q_gpu, q_aux_gpu, fl_gpu, fallbacks)
   !< `compute_riemann_face_fluxes_dev` for WENO-3: the stencil half-width S = 2 is a compile-time constant.
   use :: adam_fnl_weno_width_kernels, only : weno_width => weno_reconstruct_upwind_wratio_s2_dev
   integer(I4P), parameter :: S = 2_I4P !< WENO stencil half-width.
#define FNL_WENO_WRATIO_WIDTH
#include "adam_flume_fnl_riemann_face_fluxes_kernel_body.INC"
#undef FNL_WENO_WRATIO_WIDTH
   endsubroutine compute_riemann_face_fluxes_s2_dev

   subroutine compute_riemann_face_fluxes_s3_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                         is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, &
                                         weno_sigma, cc, tau, q_gpu, q_aux_gpu, fl_gpu, fallbacks)
   !< `compute_riemann_face_fluxes_dev` for WENO-5: the stencil half-width S = 3 is a compile-time constant.
   use :: adam_fnl_weno_width_kernels, only : weno_width => weno_reconstruct_upwind_wratio_s3_dev
   integer(I4P), parameter :: S = 3_I4P !< WENO stencil half-width.
#define FNL_WENO_WRATIO_WIDTH
#include "adam_flume_fnl_riemann_face_fluxes_kernel_body.INC"
#undef FNL_WENO_WRATIO_WIDTH
   endsubroutine compute_riemann_face_fluxes_s3_dev

   subroutine compute_riemann_face_fluxes_s4_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                         is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, &
                                         weno_sigma, cc, tau, q_gpu, q_aux_gpu, fl_gpu, fallbacks)
   !< `compute_riemann_face_fluxes_dev` for WENO-7: the stencil half-width S = 4 is a compile-time constant.
   use :: adam_fnl_weno_width_kernels, only : weno_width => weno_reconstruct_upwind_wratio_s4_dev
   integer(I4P), parameter :: S = 4_I4P !< WENO stencil half-width.
#define FNL_WENO_WRATIO_WIDTH
#include "adam_flume_fnl_riemann_face_fluxes_kernel_body.INC"
#undef FNL_WENO_WRATIO_WIDTH
   endsubroutine compute_riemann_face_fluxes_s4_dev

   subroutine compute_riemann_face_fluxes_s5_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                         is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, &
                                         weno_sigma, cc, tau, q_gpu, q_aux_gpu, fl_gpu, fallbacks)
   !< `compute_riemann_face_fluxes_dev` for WENO-9: the stencil half-width S = 5 is a compile-time constant.
   use :: adam_fnl_weno_width_kernels, only : weno_width => weno_reconstruct_upwind_wratio_s5_dev
   integer(I4P), parameter :: S = 5_I4P !< WENO stencil half-width.
#define FNL_WENO_WRATIO_WIDTH
#include "adam_flume_fnl_riemann_face_fluxes_kernel_body.INC"
#undef FNL_WENO_WRATIO_WIDTH
   endsubroutine compute_riemann_face_fluxes_s5_dev

   subroutine compute_riemann_face_fluxes_sn_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, S, gamma, ch, &
                                         is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, &
                                         weno_sigma, cc, tau, q_gpu, q_aux_gpu, fl_gpu, fallbacks)
   !< `compute_riemann_face_fluxes_dev` for any width: the stencil half-width S is a run-time argument.
   integer(I4P), intent(in) :: S !< WENO stencil half-width.
#include "adam_flume_fnl_riemann_face_fluxes_kernel_body.INC"
   endsubroutine compute_riemann_face_fluxes_sn_dev

   subroutine add_gravity_edge_correction_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, S, gamma, ch,        &
                                              is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, &
                                              weno_sigma, cc, tau, ir, gd, q_gpu, q_aux_gpu, gc_gpu)
   !< Add the faces beyond the block to the gravity energy correction of direction `d` (issue #73, M7-G2a; device twin
   !< of the CPU `add_gravity_edge_correction` of `weno-riemann`): each face computed as
   !< `compute_riemann_face_fluxes_dev` does (the same body, its mass component kept).
   integer(I4P), intent(in)    :: d                                       !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)    :: di, dj, dk                              !< Unit step along `d`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                         !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                           !< Actual blocks number.
   integer(I4P), intent(in)    :: S                                       !< WENO stencil half-width.
   real(R8P),    intent(in)    :: gamma                                   !< Specific heats ratio.
   real(R8P),    intent(in)    :: ch                                      !< GLM cleaning speed (models with GLM only).
   logical,      intent(in)    :: is_characteristic                       !< Characteristic (or primitive) variables.
   real(R8P),    intent(in)    :: weno_a_gpu(1:,0:,1:)                    !< WENO interpolation optimal weights.
   real(R8P),    intent(in)    :: weno_p_gpu(1:,0:,0:,1:)                 !< WENO interpolation polynomials coefficients.
   real(R8P),    intent(in)    :: weno_d_gpu(0:,0:,0:,1:)                 !< WENO smoothness indicators coefficients.
   real(R8P),    intent(in)    :: weno_zeps                               !< WENO parameter avoiding division by zero.
   real(R8P),    intent(in)    :: weno_sigma                              !< WENO descaler switch (0 Jiang-Shu, 1 SI).
   real(R8P),    intent(in)    :: cc(3)                                   !< Correction coefficients.
   real(R8P),    intent(in)    :: tau                                     !< Sensor threshold.
   integer(I4P), intent(in)    :: ir                                      !< Index of the density (mass) variable.
   real(R8P),    intent(in)    :: gd                                      !< Gravity along `d`.
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)       !< Conservative variables.
   real(R8P),    intent(in)    :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Auxiliary variables.
   real(R8P),    intent(inout) :: gc_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Gravity energy correction [..., 1].

   select case(S)
   case(1) ; call add_gravity_edge_correction_s1_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                                    is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, &
                                                    weno_sigma, cc, tau, ir, gd, q_gpu, q_aux_gpu, gc_gpu)
   case(2) ; call add_gravity_edge_correction_s2_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                                    is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, &
                                                    weno_sigma, cc, tau, ir, gd, q_gpu, q_aux_gpu, gc_gpu)
   case(3) ; call add_gravity_edge_correction_s3_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                                    is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, &
                                                    weno_sigma, cc, tau, ir, gd, q_gpu, q_aux_gpu, gc_gpu)
   case(4) ; call add_gravity_edge_correction_s4_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                                    is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, &
                                                    weno_sigma, cc, tau, ir, gd, q_gpu, q_aux_gpu, gc_gpu)
   case(5) ; call add_gravity_edge_correction_s5_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                                    is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, &
                                                    weno_sigma, cc, tau, ir, gd, q_gpu, q_aux_gpu, gc_gpu)
   case default ; call add_gravity_edge_correction_sn_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, S, gamma, &
                                                    ch, is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, &
                                                    weno_zeps, weno_sigma, cc, tau, ir, gd, q_gpu, q_aux_gpu, gc_gpu)
   endselect
   endsubroutine add_gravity_edge_correction_dev

   subroutine add_gravity_edge_correction_s1_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                         is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, &
                                         weno_sigma, cc, tau, ir, gd, q_gpu, q_aux_gpu, gc_gpu)
   !< `add_gravity_edge_correction_dev` for WENO-1: the stencil half-width S = 1 is a compile-time constant.
   use :: adam_fnl_weno_width_kernels, only : weno_width => weno_reconstruct_upwind_wratio_s1_dev
   integer(I4P), parameter :: S = 1_I4P !< WENO stencil half-width.
#define FNL_WENO_WRATIO_WIDTH
#include "adam_flume_fnl_riemann_gravity_edge_kernel_body.INC"
#undef FNL_WENO_WRATIO_WIDTH
   endsubroutine add_gravity_edge_correction_s1_dev

   subroutine add_gravity_edge_correction_s2_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                         is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, &
                                         weno_sigma, cc, tau, ir, gd, q_gpu, q_aux_gpu, gc_gpu)
   !< `add_gravity_edge_correction_dev` for WENO-3: the stencil half-width S = 2 is a compile-time constant.
   use :: adam_fnl_weno_width_kernels, only : weno_width => weno_reconstruct_upwind_wratio_s2_dev
   integer(I4P), parameter :: S = 2_I4P !< WENO stencil half-width.
#define FNL_WENO_WRATIO_WIDTH
#include "adam_flume_fnl_riemann_gravity_edge_kernel_body.INC"
#undef FNL_WENO_WRATIO_WIDTH
   endsubroutine add_gravity_edge_correction_s2_dev

   subroutine add_gravity_edge_correction_s3_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                         is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, &
                                         weno_sigma, cc, tau, ir, gd, q_gpu, q_aux_gpu, gc_gpu)
   !< `add_gravity_edge_correction_dev` for WENO-5: the stencil half-width S = 3 is a compile-time constant.
   use :: adam_fnl_weno_width_kernels, only : weno_width => weno_reconstruct_upwind_wratio_s3_dev
   integer(I4P), parameter :: S = 3_I4P !< WENO stencil half-width.
#define FNL_WENO_WRATIO_WIDTH
#include "adam_flume_fnl_riemann_gravity_edge_kernel_body.INC"
#undef FNL_WENO_WRATIO_WIDTH
   endsubroutine add_gravity_edge_correction_s3_dev

   subroutine add_gravity_edge_correction_s4_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                         is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, &
                                         weno_sigma, cc, tau, ir, gd, q_gpu, q_aux_gpu, gc_gpu)
   !< `add_gravity_edge_correction_dev` for WENO-7: the stencil half-width S = 4 is a compile-time constant.
   use :: adam_fnl_weno_width_kernels, only : weno_width => weno_reconstruct_upwind_wratio_s4_dev
   integer(I4P), parameter :: S = 4_I4P !< WENO stencil half-width.
#define FNL_WENO_WRATIO_WIDTH
#include "adam_flume_fnl_riemann_gravity_edge_kernel_body.INC"
#undef FNL_WENO_WRATIO_WIDTH
   endsubroutine add_gravity_edge_correction_s4_dev

   subroutine add_gravity_edge_correction_s5_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                         is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, &
                                         weno_sigma, cc, tau, ir, gd, q_gpu, q_aux_gpu, gc_gpu)
   !< `add_gravity_edge_correction_dev` for WENO-9: the stencil half-width S = 5 is a compile-time constant.
   use :: adam_fnl_weno_width_kernels, only : weno_width => weno_reconstruct_upwind_wratio_s5_dev
   integer(I4P), parameter :: S = 5_I4P !< WENO stencil half-width.
#define FNL_WENO_WRATIO_WIDTH
#include "adam_flume_fnl_riemann_gravity_edge_kernel_body.INC"
#undef FNL_WENO_WRATIO_WIDTH
   endsubroutine add_gravity_edge_correction_s5_dev

   subroutine add_gravity_edge_correction_sn_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, S, gamma, ch, &
                                         is_characteristic, weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, &
                                         weno_sigma, cc, tau, ir, gd, q_gpu, q_aux_gpu, gc_gpu)
   !< `add_gravity_edge_correction_dev` for any width: the stencil half-width S is a run-time argument.
   integer(I4P), intent(in) :: S !< WENO stencil half-width.
#include "adam_flume_fnl_riemann_gravity_edge_kernel_body.INC"
   endsubroutine add_gravity_edge_correction_sn_dev
