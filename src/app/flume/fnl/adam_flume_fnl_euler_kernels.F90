!< ADAM, FLUME FNL device kernels of the Euler model.

#include "fundal.H"

module adam_flume_fnl_euler_kernels
!< ADAM, FLUME FNL device kernels of the Euler model.
!<
!< The model-agnostic kernel bodies (`adam_flume_fnl_model_kernels_agnostic.INC`) instantiated on the Euler physics: the
!< private arrays are sized by the Euler constants `NV_K = NV_EULER`, `NV_AUX_K = NV_AUX` (issue #41, section 4), plus
!< the Euler kernels whose body depends on the variables (conservation integrals, signal speed). Same kernel rules as
!< `adam_flume_fnl_kernels` (issue #35, D-11/D-12).

! ADAM FNL classes, libraries
use :: adam_fnl_weno_kernels,    only : weno_reconstruct_upwind_dev
! FLUME modules
use :: adam_flume_euler_library, only : compute_face_flux_back_projection, compute_face_split_fluxes,                   &
                                        conservative_to_auxiliary
use :: adam_flume_parameters,    only : IA_A, IA_U, IA_V, IA_W, NV_AUX_K=>NV_AUX, NV_K=>NV_EULER, S_MAX
! third party modules
use :: penf,                     only : I4P, R8P

implicit none
private
public :: compute_conservation_dev
public :: compute_face_fluxes_dev
public :: compute_lambda_max_dev
public :: compute_q_aux_dev

contains
   ! public procedures
#include "adam_flume_fnl_model_kernels_agnostic.INC"

   subroutine compute_conservation_dev(ni, nj, nk, ngc, blocks_number, dxyz_gpu, q_gpu, integrals)
   !< Compute the volume integrals of the conservative variables (interior cells). The cell volume includes the null
   !< directions: the tree splits them too, so a refined block's cells are smaller along them (issue #37).
   integer(I4P), intent(in)  :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)  :: blocks_number                     !< Actual blocks number.
   real(R8P),    intent(in)  :: dxyz_gpu(1:,1:)                   !< Blocks space steps [nb, 3].
   real(R8P),    intent(in)  :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
   real(R8P),    intent(out) :: integrals(NV_K)                   !< Volume integrals.
   real(R8P)                 :: volume                            !< Cell volume.
   real(R8P)                 :: s1, s2, s3, s4, s5                !< Reduction accumulators.
   integer(I4P)              :: b, i, j, k                        !< Counters.

   s1 = 0._R8P ; s2 = 0._R8P ; s3 = 0._R8P ; s4 = 0._R8P ; s5 = 0._R8P
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu) &
   !$acc& firstprivate(ni,nj,nk,blocks_number) private(volume)                   &
   !$acc& reduction(+:s1,s2,s3,s4,s5)
   !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,q_gpu) &
   !$omp& firstprivate(ni,nj,nk,blocks_number) private(volume) &
   !$omp& reduction(+:s1,s2,s3,s4,s5)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      volume = dxyz_gpu(b,1) * dxyz_gpu(b,2) * dxyz_gpu(b,3)
      s1 = s1 + q_gpu(b,i,j,k,1) * volume
      s2 = s2 + q_gpu(b,i,j,k,2) * volume
      s3 = s3 + q_gpu(b,i,j,k,3) * volume
      s4 = s4 + q_gpu(b,i,j,k,4) * volume
      s5 = s5 + q_gpu(b,i,j,k,5) * volume
   enddo
   enddo
   enddo
   enddo
   integrals = [s1, s2, s3, s4, s5]
   endsubroutine compute_conservation_dev

   subroutine compute_lambda_max_dev(ni, nj, nk, ngc, blocks_number, gamma, R, dxyz_gpu, is_null, q_gpu, lambda_max)
   !< Compute `max(sum_d (|u_d| + a) / dx_d)` over the interior cells (null directions excluded).
   integer(I4P), intent(in)  :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)  :: blocks_number                     !< Actual blocks number.
   real(R8P),    intent(in)  :: gamma                             !< Specific heats ratio.
   real(R8P),    intent(in)  :: R                                 !< Gas constant.
   real(R8P),    intent(in)  :: dxyz_gpu(1:,1:)                   !< Blocks space steps [nb, 3].
   logical,      intent(in)  :: is_null(3)                        !< Null directions.
   real(R8P),    intent(in)  :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
   real(R8P),    intent(out) :: lambda_max                        !< Maximum of sum_d (|u_d| + a) / dx_d.
   real(R8P)                 :: wx, wy, wz                        !< Direction weights: 1 active, 0 null.
   real(R8P)                 :: q_(NV_K)                          !< Private conservative variables of one cell.
   real(R8P)                 :: qa_(NV_AUX_K)                     !< Private auxiliary variables of one cell.
   real(R8P)                 :: lambda                            !< Cell value.
   integer(I4P)              :: b, i, j, k, v                     !< Counters.

   wx = merge(0._R8P, 1._R8P, is_null(1))
   wy = merge(0._R8P, 1._R8P, is_null(2))
   wz = merge(0._R8P, 1._R8P, is_null(3))
   lambda_max = 0._R8P
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu) &
   !$acc& firstprivate(ni,nj,nk,blocks_number,gamma,R,wx,wy,wz) private(q_,qa_,lambda) &
   !$acc& reduction(max:lambda_max)
   !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,q_gpu) &
   !$omp& firstprivate(ni,nj,nk,blocks_number,gamma,R,wx,wy,wz) private(q_,qa_,lambda) &
   !$omp& reduction(max:lambda_max)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      !$acc loop seq
      do v=1, NV_K
         q_(v) = q_gpu(b,i,j,k,v)
      enddo
      call conservative_to_auxiliary(gamma=gamma, R=R, q=q_, qa=qa_)
      lambda = wx * (abs(qa_(IA_U)) + qa_(IA_A)) / dxyz_gpu(b,1) + &
               wy * (abs(qa_(IA_V)) + qa_(IA_A)) / dxyz_gpu(b,2) + &
               wz * (abs(qa_(IA_W)) + qa_(IA_A)) / dxyz_gpu(b,3)
      lambda_max = max(lambda_max, lambda)
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_lambda_max_dev
endmodule adam_flume_fnl_euler_kernels
