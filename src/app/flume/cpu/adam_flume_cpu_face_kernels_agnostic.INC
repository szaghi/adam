   ! Model-agnostic CPU face-flux loop: the body shared by every per-model kernel module (issue #41, section 4).
   !
   ! Included in the `contains` part of `adam_flume_cpu_<model>_kernels`, which provides, by use association:
   ! `NV_K`, `NV_AUX_K` (the model's conservative and auxiliary variables numbers, compile-time constants that size
   ! every local array: no model branching and no maximum-model buffers in the loops), `S_MAX`, the model's
   ! split adapter `face_split_fluxes(gamma, ch, ...)` (one argument list for every model: the GLM cleaning speed `ch` is
   ! ignored by the models without cleaning), `compute_face_flux_back_projection`, and `weno_object` (the WENO
   ! configuration). The host selects the model's module by `select case`, outside the loops.
   !
   ! One kernel per stencil width, as on FNL (issue #139, S3a, decision D3): each public routine is a host dispatcher that
   ! selects, by the WENO half-width S, a kernel compiled with S a named constant (`_s1` ... `_s5`, WENO-1 to WENO-9);
   ! every kernel compiles the same body (`adam_flume_cpu_face_fluxes_body.INC`, `adam_flume_cpu_gravity_edge_body.INC`)
   ! and calls SETH's WENO-JS kernel of its width, renamed `weno_width` (`seth_weno<2S-1>_js`).

   subroutine compute_face_fluxes(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, weno, &
                                  q, q_aux, fl)
   !< Compute the WENO face fluxes of direction `d`: face `(i,j,k)` lies between cells `(i,j,k)` and
   !< `(i+di,j+dj,k+dk)`, so the face array starts at index 0 along `d`.
   !<
   !< Per face: gather the stencil `m = 1-S ... S` into constant-bound locals, project and split it, reconstruct every
   !< field (the local `v` is the packed stencil of the WENO primitive), back-project.
   !<
   !< Host dispatcher: the kernel of the WENO half-width S, compiled with S a named constant (`_s1` ... `_s5`).
   integer(I4P),      intent(in)    :: d                                     !< Direction, 1=x, 2=y, 3=z.
   integer(I4P),      intent(in)    :: di, dj, dk                            !< Unit step along `d`.
   integer(I4P),      intent(in)    :: ni, nj, nk, ngc                       !< Grid dimensions.
   integer(I4P),      intent(in)    :: blocks_number                         !< Actual blocks number.
   real(R8P),         intent(in)    :: gamma                                 !< Specific heats ratio.
   real(R8P),         intent(in)    :: ch                                    !< GLM cleaning speed (models with GLM only).
   logical,           intent(in)    :: is_characteristic                     !< Characteristic (or conservative) variables.
   type(weno_object), intent(in)    :: weno                                  !< WENO coefficients.
   real(R8P),         intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)         !< Conservative variables.
   real(R8P),         intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Auxiliary variables.
   real(R8P),         intent(inout) :: fl(1:,1-di:,1-dj:,1-dk:,1:)           !< Face fluxes of direction `d`.

   select case(weno%S)
   case(1) ; call compute_face_fluxes_s1(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                         weno, q, q_aux, fl)
   case(2) ; call compute_face_fluxes_s2(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                         weno, q, q_aux, fl)
   case(3) ; call compute_face_fluxes_s3(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                         weno, q, q_aux, fl)
   case(4) ; call compute_face_fluxes_s4(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                         weno, q, q_aux, fl)
   case(5) ; call compute_face_fluxes_s5(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                         weno, q, q_aux, fl)
   case default ; error stop 'compute_face_fluxes: the WENO half-width S must be 1 to 5'
   endselect
   endsubroutine compute_face_fluxes

   subroutine compute_face_fluxes_s1(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, weno, &
                                     q, q_aux, fl)
   !< `compute_face_fluxes` for WENO-1: the stencil half-width S = 1 is a compile-time constant.
   use :: seth, only : weno_width => seth_weno1_js
   integer(I4P), parameter :: S = 1_I4P !< WENO stencil half-width.
#include "adam_flume_cpu_face_fluxes_body.INC"
   endsubroutine compute_face_fluxes_s1

   subroutine compute_face_fluxes_s2(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, weno, &
                                     q, q_aux, fl)
   !< `compute_face_fluxes` for WENO-3: the stencil half-width S = 2 is a compile-time constant.
   use :: seth, only : weno_width => seth_weno3_js
   integer(I4P), parameter :: S = 2_I4P !< WENO stencil half-width.
#include "adam_flume_cpu_face_fluxes_body.INC"
   endsubroutine compute_face_fluxes_s2

   subroutine compute_face_fluxes_s3(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, weno, &
                                     q, q_aux, fl)
   !< `compute_face_fluxes` for WENO-5: the stencil half-width S = 3 is a compile-time constant.
   use :: seth, only : weno_width => seth_weno5_js
   integer(I4P), parameter :: S = 3_I4P !< WENO stencil half-width.
#include "adam_flume_cpu_face_fluxes_body.INC"
   endsubroutine compute_face_fluxes_s3

   subroutine compute_face_fluxes_s4(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, weno, &
                                     q, q_aux, fl)
   !< `compute_face_fluxes` for WENO-7: the stencil half-width S = 4 is a compile-time constant.
   use :: seth, only : weno_width => seth_weno7_js
   integer(I4P), parameter :: S = 4_I4P !< WENO stencil half-width.
#include "adam_flume_cpu_face_fluxes_body.INC"
   endsubroutine compute_face_fluxes_s4

   subroutine compute_face_fluxes_s5(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, weno, &
                                     q, q_aux, fl)
   !< `compute_face_fluxes` for WENO-9: the stencil half-width S = 5 is a compile-time constant.
   use :: seth, only : weno_width => seth_weno9_js
   integer(I4P), parameter :: S = 5_I4P !< WENO stencil half-width.
#include "adam_flume_cpu_face_fluxes_body.INC"
   endsubroutine compute_face_fluxes_s5

   subroutine add_gravity_edge_correction(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                          weno, ir, gd, q, q_aux, gc)
   !< Add the faces beyond the block to the gravity energy correction of direction `d` (issue #73, M7-G2a): cell 1 takes
   !< `-g_d F^rho / 24` of face -1 (between cells -1 and 0), cell `n_d` that of face `n_d + 1` (between cells `n_d + 1`
   !< and `n_d + 2`), the two terms of `g_d (-F_{i-3/2} + F_{i-1/2} + F_{i+1/2} - F_{i+3/2}) / 24` whose faces the
   !< block does not hold. Each face is computed as `compute_face_fluxes` does (the same body, its mass component kept),
   !< from the ghosts: its stencil reaches `S + 1` ghost cells, which the initialisation requires with gravity.
   !<
   !< Host dispatcher: the kernel of the WENO half-width S, compiled with S a named constant (`_s1` ... `_s5`).
   integer(I4P),      intent(in)    :: d                                     !< Direction, 1=x, 2=y, 3=z.
   integer(I4P),      intent(in)    :: di, dj, dk                            !< Unit step along `d`.
   integer(I4P),      intent(in)    :: ni, nj, nk, ngc                       !< Grid dimensions.
   integer(I4P),      intent(in)    :: blocks_number                         !< Actual blocks number.
   real(R8P),         intent(in)    :: gamma                                 !< Specific heats ratio.
   real(R8P),         intent(in)    :: ch                                    !< GLM cleaning speed (models with GLM only).
   logical,           intent(in)    :: is_characteristic                     !< Characteristic (or conservative) variables.
   type(weno_object), intent(in)    :: weno                                  !< WENO coefficients.
   integer(I4P),      intent(in)    :: ir                                    !< Index of the density (mass) variable.
   real(R8P),         intent(in)    :: gd                                    !< Gravity along `d`.
   real(R8P),         intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)         !< Conservative variables.
   real(R8P),         intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Auxiliary variables.
   real(R8P),         intent(inout) :: gc(1:,1-ngc:,1-ngc:,1-ngc:,1:)        !< Gravity energy correction [1, ...].

   select case(weno%S)
   case(1) ; call add_gravity_edge_correction_s1(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                                 is_characteristic, weno, ir, gd, q, q_aux, gc)
   case(2) ; call add_gravity_edge_correction_s2(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                                 is_characteristic, weno, ir, gd, q, q_aux, gc)
   case(3) ; call add_gravity_edge_correction_s3(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                                 is_characteristic, weno, ir, gd, q, q_aux, gc)
   case(4) ; call add_gravity_edge_correction_s4(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                                 is_characteristic, weno, ir, gd, q, q_aux, gc)
   case(5) ; call add_gravity_edge_correction_s5(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                                 is_characteristic, weno, ir, gd, q, q_aux, gc)
   case default ; error stop 'add_gravity_edge_correction: the WENO half-width S must be 1 to 5'
   endselect
   endsubroutine add_gravity_edge_correction

   subroutine add_gravity_edge_correction_s1(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                             weno, ir, gd, q, q_aux, gc)
   !< `add_gravity_edge_correction` for WENO-1: the stencil half-width S = 1 is a compile-time constant.
   use :: seth, only : weno_width => seth_weno1_js
   integer(I4P), parameter :: S = 1_I4P !< WENO stencil half-width.
#include "adam_flume_cpu_gravity_edge_body.INC"
   endsubroutine add_gravity_edge_correction_s1

   subroutine add_gravity_edge_correction_s2(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                             weno, ir, gd, q, q_aux, gc)
   !< `add_gravity_edge_correction` for WENO-3: the stencil half-width S = 2 is a compile-time constant.
   use :: seth, only : weno_width => seth_weno3_js
   integer(I4P), parameter :: S = 2_I4P !< WENO stencil half-width.
#include "adam_flume_cpu_gravity_edge_body.INC"
   endsubroutine add_gravity_edge_correction_s2

   subroutine add_gravity_edge_correction_s3(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                             weno, ir, gd, q, q_aux, gc)
   !< `add_gravity_edge_correction` for WENO-5: the stencil half-width S = 3 is a compile-time constant.
   use :: seth, only : weno_width => seth_weno5_js
   integer(I4P), parameter :: S = 3_I4P !< WENO stencil half-width.
#include "adam_flume_cpu_gravity_edge_body.INC"
   endsubroutine add_gravity_edge_correction_s3

   subroutine add_gravity_edge_correction_s4(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                             weno, ir, gd, q, q_aux, gc)
   !< `add_gravity_edge_correction` for WENO-7: the stencil half-width S = 4 is a compile-time constant.
   use :: seth, only : weno_width => seth_weno7_js
   integer(I4P), parameter :: S = 4_I4P !< WENO stencil half-width.
#include "adam_flume_cpu_gravity_edge_body.INC"
   endsubroutine add_gravity_edge_correction_s4

   subroutine add_gravity_edge_correction_s5(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                             weno, ir, gd, q, q_aux, gc)
   !< `add_gravity_edge_correction` for WENO-9: the stencil half-width S = 5 is a compile-time constant.
   use :: seth, only : weno_width => seth_weno9_js
   integer(I4P), parameter :: S = 5_I4P !< WENO stencil half-width.
#include "adam_flume_cpu_gravity_edge_body.INC"
   endsubroutine add_gravity_edge_correction_s5
