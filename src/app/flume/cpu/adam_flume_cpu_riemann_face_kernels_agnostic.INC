   ! Model- and solver-agnostic CPU face-flux loop of the scheme `weno-riemann` (issue #47): the body shared by every
   ! per-(model, solver) kernel module.
   !
   ! Included in the `contains` part of `adam_flume_cpu_<model>_<solver>_kernels`, which provides, by use association:
   ! `NV_K`, `NV_AUX_K` (compile-time constants sizing every local array), `S_MAX`, `weno_object` (the WENO
   ! configuration), and four adapters with one argument list for every model (the GLM cleaning speed
   ! `ch` is ignored by the models without cleaning): `face_interpolation_fields` (the stencil fields to interpolate, in
   ! the WENO upwind layout), `face_states` (the conservative face states, with the first-order fallback),
   ! `riemann_flux` (the solver, with a flag set when it fell back to a simpler solver, e.g. HLLD to HLL) and `cell_flux`
! (the physical flux of a cell). The host selects the module by
   ! `select case` on the model and the solver, outside the loops.
   !
   ! One kernel per stencil width, as on FNL (issue #139, S3a, decision D3): each public routine is a host dispatcher that
   ! selects, by the WENO half-width S, a kernel compiled with S a named constant (`_s1` ... `_s5`, WENO-1 to WENO-9);
   ! every kernel compiles the same body (`adam_flume_cpu_riemann_face_fluxes_body.INC`,
   ! `adam_flume_cpu_riemann_gravity_edge_body.INC`) and calls SETH's WENO-JS interpolation kernel of its width with the
   ! weights ratios, renamed `weno_width` (`seth_weno<2S-1>_js_interp_wr`).

   subroutine compute_riemann_face_fluxes(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                          weno, cc, tau, q, q_aux, fl, fallbacks)
   !< Compute the `weno-riemann` face fluxes of direction `d`: face `(i,j,k)` lies between cells `(i,j,k)` and
   !< `(i+di,j+dj,k+dk)`, so the face array starts at index 0 along `d`.
   !<
   !< Per face: gather the stencil `m = min(1-S,-1) ... max(S,2)`, WENO-interpolate its fields to the two face states
   !< (interpolation tables `weno%a_interp`, `weno%p_interp`), take the Riemann flux `F` of the states and add the
   !< correction `F^ = F + s (cc(1) F + cc(2) (f_0 + f_1) + cc(3) (f_-1 + f_2) - F)`, with `f_m` the physical flux of
   !< cell `m` and the sensor `s = 1` where the WENO weights ratio of every field on both sides is `>= tau`, else 0;
   !< `fallbacks` counts the faces whose Riemann solver fell back (HLLD to HLL; always 0 for the other solvers).
   !<
   !< Host dispatcher: the kernel of the WENO half-width S, compiled with S a named constant (`_s1` ... `_s5`).
   integer(I4P),      intent(in)    :: d                                 !< Direction, 1=x, 2=y, 3=z.
   integer(I4P),      intent(in)    :: di, dj, dk                        !< Unit step along `d`.
   integer(I4P),      intent(in)    :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P),      intent(in)    :: blocks_number                     !< Actual blocks number.
   real(R8P),         intent(in)    :: gamma                             !< Specific heats ratio.
   real(R8P),         intent(in)    :: ch                                !< GLM cleaning speed (models with GLM only).
   logical,           intent(in)    :: is_characteristic                 !< Characteristic (or primitive) variables.
   type(weno_object), intent(in)    :: weno                              !< WENO coefficients (with interpolation tables).
   real(R8P),         intent(in)    :: cc(3)                             !< Correction coefficients.
   real(R8P),         intent(in)    :: tau                               !< Sensor threshold.
   real(R8P),         intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Conservative variables.
   real(R8P),         intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),         intent(inout) :: fl(1:,1-di:,1-dj:,1-dk:,1:)       !< Face fluxes of direction `d`.
   integer(I4P),      intent(out)   :: fallbacks                         !< Faces whose Riemann solver fell back.

   select case(weno%S)
   case(1) ; call compute_riemann_face_fluxes_s1(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                                 is_characteristic, weno, cc, tau, q, q_aux, fl, fallbacks)
   case(2) ; call compute_riemann_face_fluxes_s2(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                                 is_characteristic, weno, cc, tau, q, q_aux, fl, fallbacks)
   case(3) ; call compute_riemann_face_fluxes_s3(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                                 is_characteristic, weno, cc, tau, q, q_aux, fl, fallbacks)
   case(4) ; call compute_riemann_face_fluxes_s4(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                                 is_characteristic, weno, cc, tau, q, q_aux, fl, fallbacks)
   case(5) ; call compute_riemann_face_fluxes_s5(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                                 is_characteristic, weno, cc, tau, q, q_aux, fl, fallbacks)
   case default ; error stop 'compute_riemann_face_fluxes: the WENO half-width S must be 1 to 5'
   endselect
   endsubroutine compute_riemann_face_fluxes

   subroutine compute_riemann_face_fluxes_s1(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                             weno, cc, tau, q, q_aux, fl, fallbacks)
   !< `compute_riemann_face_fluxes` for WENO-1: the stencil half-width S = 1 is a compile-time constant.
   use :: seth, only : weno_width => seth_weno1_js_interp_wr
   integer(I4P), parameter :: S = 1_I4P !< WENO stencil half-width.
#include "adam_flume_cpu_riemann_face_fluxes_body.INC"
   endsubroutine compute_riemann_face_fluxes_s1

   subroutine compute_riemann_face_fluxes_s2(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                             weno, cc, tau, q, q_aux, fl, fallbacks)
   !< `compute_riemann_face_fluxes` for WENO-3: the stencil half-width S = 2 is a compile-time constant.
   use :: seth, only : weno_width => seth_weno3_js_interp_wr
   integer(I4P), parameter :: S = 2_I4P !< WENO stencil half-width.
#include "adam_flume_cpu_riemann_face_fluxes_body.INC"
   endsubroutine compute_riemann_face_fluxes_s2

   subroutine compute_riemann_face_fluxes_s3(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                             weno, cc, tau, q, q_aux, fl, fallbacks)
   !< `compute_riemann_face_fluxes` for WENO-5: the stencil half-width S = 3 is a compile-time constant.
   use :: seth, only : weno_width => seth_weno5_js_interp_wr
   integer(I4P), parameter :: S = 3_I4P !< WENO stencil half-width.
#include "adam_flume_cpu_riemann_face_fluxes_body.INC"
   endsubroutine compute_riemann_face_fluxes_s3

   subroutine compute_riemann_face_fluxes_s4(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                             weno, cc, tau, q, q_aux, fl, fallbacks)
   !< `compute_riemann_face_fluxes` for WENO-7: the stencil half-width S = 4 is a compile-time constant.
   use :: seth, only : weno_width => seth_weno7_js_interp_wr
   integer(I4P), parameter :: S = 4_I4P !< WENO stencil half-width.
#include "adam_flume_cpu_riemann_face_fluxes_body.INC"
   endsubroutine compute_riemann_face_fluxes_s4

   subroutine compute_riemann_face_fluxes_s5(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                             weno, cc, tau, q, q_aux, fl, fallbacks)
   !< `compute_riemann_face_fluxes` for WENO-9: the stencil half-width S = 5 is a compile-time constant.
   use :: seth, only : weno_width => seth_weno9_js_interp_wr
   integer(I4P), parameter :: S = 5_I4P !< WENO stencil half-width.
#include "adam_flume_cpu_riemann_face_fluxes_body.INC"
   endsubroutine compute_riemann_face_fluxes_s5

   subroutine add_gravity_edge_correction(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                          weno, cc, tau, ir, gd, q, q_aux, gc)
   !< Add the faces beyond the block to the gravity energy correction of direction `d` (issue #73, M7-G2a), as the
   !< split kernel's `add_gravity_edge_correction`, with the `weno-riemann` face of `compute_riemann_face_fluxes` (the
   !< same body, its mass component kept): its stencil reaches `max(S, 2) + 1` ghost cells.
   !<
   !< Host dispatcher: the kernel of the WENO half-width S, compiled with S a named constant (`_s1` ... `_s5`).
   integer(I4P),      intent(in)    :: d                                 !< Direction, 1=x, 2=y, 3=z.
   integer(I4P),      intent(in)    :: di, dj, dk                        !< Unit step along `d`.
   integer(I4P),      intent(in)    :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P),      intent(in)    :: blocks_number                     !< Actual blocks number.
   real(R8P),         intent(in)    :: gamma                             !< Specific heats ratio.
   real(R8P),         intent(in)    :: ch                                !< GLM cleaning speed (models with GLM only).
   logical,           intent(in)    :: is_characteristic                 !< Characteristic (or primitive) variables.
   type(weno_object), intent(in)    :: weno                              !< WENO coefficients (with interpolation tables).
   real(R8P),         intent(in)    :: cc(3)                             !< Correction coefficients.
   real(R8P),         intent(in)    :: tau                               !< Sensor threshold.
   integer(I4P),      intent(in)    :: ir                                !< Index of the density (mass) variable.
   real(R8P),         intent(in)    :: gd                                !< Gravity along `d`.
   real(R8P),         intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Conservative variables.
   real(R8P),         intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),         intent(inout) :: gc(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Gravity energy correction [1, ...].

   select case(weno%S)
   case(1) ; call add_gravity_edge_correction_s1(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                                 is_characteristic, weno, cc, tau, ir, gd, q, q_aux, gc)
   case(2) ; call add_gravity_edge_correction_s2(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                                 is_characteristic, weno, cc, tau, ir, gd, q, q_aux, gc)
   case(3) ; call add_gravity_edge_correction_s3(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                                 is_characteristic, weno, cc, tau, ir, gd, q, q_aux, gc)
   case(4) ; call add_gravity_edge_correction_s4(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                                 is_characteristic, weno, cc, tau, ir, gd, q, q_aux, gc)
   case(5) ; call add_gravity_edge_correction_s5(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, &
                                                 is_characteristic, weno, cc, tau, ir, gd, q, q_aux, gc)
   case default ; error stop 'add_gravity_edge_correction: the WENO half-width S must be 1 to 5'
   endselect
   endsubroutine add_gravity_edge_correction

   subroutine add_gravity_edge_correction_s1(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                             weno, cc, tau, ir, gd, q, q_aux, gc)
   !< `add_gravity_edge_correction` for WENO-1: the stencil half-width S = 1 is a compile-time constant.
   use :: seth, only : weno_width => seth_weno1_js_interp_wr
   integer(I4P), parameter :: S = 1_I4P !< WENO stencil half-width.
#include "adam_flume_cpu_riemann_gravity_edge_body.INC"
   endsubroutine add_gravity_edge_correction_s1

   subroutine add_gravity_edge_correction_s2(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                             weno, cc, tau, ir, gd, q, q_aux, gc)
   !< `add_gravity_edge_correction` for WENO-3: the stencil half-width S = 2 is a compile-time constant.
   use :: seth, only : weno_width => seth_weno3_js_interp_wr
   integer(I4P), parameter :: S = 2_I4P !< WENO stencil half-width.
#include "adam_flume_cpu_riemann_gravity_edge_body.INC"
   endsubroutine add_gravity_edge_correction_s2

   subroutine add_gravity_edge_correction_s3(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                             weno, cc, tau, ir, gd, q, q_aux, gc)
   !< `add_gravity_edge_correction` for WENO-5: the stencil half-width S = 3 is a compile-time constant.
   use :: seth, only : weno_width => seth_weno5_js_interp_wr
   integer(I4P), parameter :: S = 3_I4P !< WENO stencil half-width.
#include "adam_flume_cpu_riemann_gravity_edge_body.INC"
   endsubroutine add_gravity_edge_correction_s3

   subroutine add_gravity_edge_correction_s4(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                             weno, cc, tau, ir, gd, q, q_aux, gc)
   !< `add_gravity_edge_correction` for WENO-7: the stencil half-width S = 4 is a compile-time constant.
   use :: seth, only : weno_width => seth_weno7_js_interp_wr
   integer(I4P), parameter :: S = 4_I4P !< WENO stencil half-width.
#include "adam_flume_cpu_riemann_gravity_edge_body.INC"
   endsubroutine add_gravity_edge_correction_s4

   subroutine add_gravity_edge_correction_s5(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                             weno, cc, tau, ir, gd, q, q_aux, gc)
   !< `add_gravity_edge_correction` for WENO-9: the stencil half-width S = 5 is a compile-time constant.
   use :: seth, only : weno_width => seth_weno9_js_interp_wr
   integer(I4P), parameter :: S = 5_I4P !< WENO stencil half-width.
#include "adam_flume_cpu_riemann_gravity_edge_body.INC"
   endsubroutine add_gravity_edge_correction_s5
