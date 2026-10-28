!< ADAM, WENO class FNL kernels, one reconstruction routine per stencil width (S a named constant).

module adam_fnl_weno_width_kernels
!< ADAM, WENO class FNL kernels, one reconstruction routine per stencil width (S a named constant).
!<
!< `weno_reconstruct_upwind_s<N>_dev` (`weno_reconstruct_upwind_wratio_s<N>_dev`) is `weno_reconstruct_upwind_dev`
!< (`weno_reconstruct_upwind_wratio_dev`) of [[adam_fnl_weno_kernels]] for the half-width S = N (WENO-1 ... WENO-9,
!< S_MAX = 5) without the argument S: the face kernels of the same width (`..._s<N>_dev` of the FLUME FNL backend) call
!< them. The generic routines stay for the callers whose width is a run-time value (the generic `_sn_dev` kernels, NASTO's
!< reduced-order fallback, RIPPLE). Both are bitwise the generic ones: the body is one include,
!< `adam_fnl_weno_width_body.INC`, and keeps the order of every operation.

! ADAM modules
use :: adam_weno_object, only : weno_weights_exponent
! third party modules
use :: penf, only : I4P, R8P

implicit none
private
public :: weno_reconstruct_upwind_s1_dev
public :: weno_reconstruct_upwind_s2_dev
public :: weno_reconstruct_upwind_s3_dev
public :: weno_reconstruct_upwind_s4_dev
public :: weno_reconstruct_upwind_s5_dev
public :: weno_reconstruct_upwind_wratio_s1_dev
public :: weno_reconstruct_upwind_wratio_s2_dev
public :: weno_reconstruct_upwind_wratio_s3_dev
public :: weno_reconstruct_upwind_wratio_s4_dev
public :: weno_reconstruct_upwind_wratio_s5_dev

contains
   subroutine weno_reconstruct_upwind_s1_dev(weno_a, weno_p, weno_d, weno_zeps, weno_rmu, V, VR)
   !< Reconstruct by WENO upwind method of 1 order, S = 1 a compile-time constant, non TBP.
   integer(I4P), parameter :: S = 1_I4P !< Number of stencils used.
#include "adam_fnl_weno_width_body.INC"
   endsubroutine weno_reconstruct_upwind_s1_dev

   subroutine weno_reconstruct_upwind_s2_dev(weno_a, weno_p, weno_d, weno_zeps, weno_rmu, V, VR)
   !< Reconstruct by WENO upwind method of 3 order, S = 2 a compile-time constant, non TBP.
   integer(I4P), parameter :: S = 2_I4P !< Number of stencils used.
#include "adam_fnl_weno_width_body.INC"
   endsubroutine weno_reconstruct_upwind_s2_dev

   subroutine weno_reconstruct_upwind_s3_dev(weno_a, weno_p, weno_d, weno_zeps, weno_rmu, V, VR)
   !< Reconstruct by WENO upwind method of 5 order, S = 3 a compile-time constant, non TBP.
   integer(I4P), parameter :: S = 3_I4P !< Number of stencils used.
#include "adam_fnl_weno_width_body.INC"
   endsubroutine weno_reconstruct_upwind_s3_dev

   subroutine weno_reconstruct_upwind_s4_dev(weno_a, weno_p, weno_d, weno_zeps, weno_rmu, V, VR)
   !< Reconstruct by WENO upwind method of 7 order, S = 4 a compile-time constant, non TBP.
   integer(I4P), parameter :: S = 4_I4P !< Number of stencils used.
#include "adam_fnl_weno_width_body.INC"
   endsubroutine weno_reconstruct_upwind_s4_dev

   subroutine weno_reconstruct_upwind_s5_dev(weno_a, weno_p, weno_d, weno_zeps, weno_rmu, V, VR)
   !< Reconstruct by WENO upwind method of 9 order, S = 5 a compile-time constant, non TBP.
   integer(I4P), parameter :: S = 5_I4P !< Number of stencils used.
#include "adam_fnl_weno_width_body.INC"
   endsubroutine weno_reconstruct_upwind_s5_dev

   subroutine weno_reconstruct_upwind_wratio_s1_dev(weno_a, weno_p, weno_d, weno_zeps, weno_rmu, V, VR, WRATIO)
   !< Reconstruct by WENO upwind method of 1 order, returning also the smoothness ratio of each interface, S = 1 a compile-time constant, non TBP.
   integer(I4P), parameter :: S = 1_I4P !< Number of stencils used.
#define FNL_WENO_WRATIO
#include "adam_fnl_weno_width_body.INC"
#undef FNL_WENO_WRATIO
   endsubroutine weno_reconstruct_upwind_wratio_s1_dev

   subroutine weno_reconstruct_upwind_wratio_s2_dev(weno_a, weno_p, weno_d, weno_zeps, weno_rmu, V, VR, WRATIO)
   !< Reconstruct by WENO upwind method of 3 order, returning also the smoothness ratio of each interface, S = 2 a compile-time constant, non TBP.
   integer(I4P), parameter :: S = 2_I4P !< Number of stencils used.
#define FNL_WENO_WRATIO
#include "adam_fnl_weno_width_body.INC"
#undef FNL_WENO_WRATIO
   endsubroutine weno_reconstruct_upwind_wratio_s2_dev

   subroutine weno_reconstruct_upwind_wratio_s3_dev(weno_a, weno_p, weno_d, weno_zeps, weno_rmu, V, VR, WRATIO)
   !< Reconstruct by WENO upwind method of 5 order, returning also the smoothness ratio of each interface, S = 3 a compile-time constant, non TBP.
   integer(I4P), parameter :: S = 3_I4P !< Number of stencils used.
#define FNL_WENO_WRATIO
#include "adam_fnl_weno_width_body.INC"
#undef FNL_WENO_WRATIO
   endsubroutine weno_reconstruct_upwind_wratio_s3_dev

   subroutine weno_reconstruct_upwind_wratio_s4_dev(weno_a, weno_p, weno_d, weno_zeps, weno_rmu, V, VR, WRATIO)
   !< Reconstruct by WENO upwind method of 7 order, returning also the smoothness ratio of each interface, S = 4 a compile-time constant, non TBP.
   integer(I4P), parameter :: S = 4_I4P !< Number of stencils used.
#define FNL_WENO_WRATIO
#include "adam_fnl_weno_width_body.INC"
#undef FNL_WENO_WRATIO
   endsubroutine weno_reconstruct_upwind_wratio_s4_dev

   subroutine weno_reconstruct_upwind_wratio_s5_dev(weno_a, weno_p, weno_d, weno_zeps, weno_rmu, V, VR, WRATIO)
   !< Reconstruct by WENO upwind method of 9 order, returning also the smoothness ratio of each interface, S = 5 a compile-time constant, non TBP.
   integer(I4P), parameter :: S = 5_I4P !< Number of stencils used.
#define FNL_WENO_WRATIO
#include "adam_fnl_weno_width_body.INC"
#undef FNL_WENO_WRATIO
   endsubroutine weno_reconstruct_upwind_wratio_s5_dev

endmodule adam_fnl_weno_width_kernels
