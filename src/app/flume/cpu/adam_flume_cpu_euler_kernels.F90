!< ADAM, FLUME CPU kernels of the Euler model.

module adam_flume_cpu_euler_kernels
!< ADAM, FLUME CPU kernels of the Euler model.
!<
!< The model-agnostic loop bodies (`adam_flume_cpu_{face,aux}_kernels_agnostic.INC`) instantiated on the Euler physics: the
!< local arrays are sized by the Euler constants `NV_K = NV_EULER`, `NV_AUX_K = NV_AUX` (issue #41, section 4), plus
!< the Euler signal-speed loop.

! ADAM classes, libraries, parameters
use :: adam_weno_object,         only : weno_object, weno_reconstruct_upwind
! FLUME modules
use :: adam_flume_euler_library, only : compute_face_flux_back_projection, compute_face_split_fluxes,                   &
                                        conservative_to_auxiliary
use :: adam_flume_parameters,    only : IA_A, IA_U, NV_AUX_K=>NV_AUX, NV_K=>NV_EULER, S_MAX
! third party modules
use :: penf,                     only : I4P, R8P

implicit none
private
public :: compute_face_fluxes
public :: compute_lambda_max
public :: compute_q_aux

contains
   ! public procedures
#include "adam_flume_cpu_face_kernels_agnostic.INC"

#include "adam_flume_cpu_aux_kernels_agnostic.INC"

   subroutine compute_lambda_max(ni, nj, nk, ngc, blocks_number, gamma, R, dxyz, is_null, q, lambda_max)
   !< Compute `max(sum_d (|u_d| + a) / dx_d)` over the interior cells (null directions excluded).
   integer(I4P), intent(in)  :: ni, nj, nk, ngc               !< Grid dimensions.
   integer(I4P), intent(in)  :: blocks_number                 !< Actual blocks number.
   real(R8P),    intent(in)  :: gamma                         !< Specific heats ratio.
   real(R8P),    intent(in)  :: R                             !< Gas constant.
   real(R8P),    intent(in)  :: dxyz(1:,1:)                   !< Blocks space steps [3, nb].
   logical,      intent(in)  :: is_null(3)                    !< Null directions.
   real(R8P),    intent(in)  :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
   real(R8P),    intent(out) :: lambda_max                    !< Maximum of sum_d (|u_d| + a) / dx_d.
   real(R8P)                 :: qa(NV_AUX_K)                  !< Auxiliary variables of one cell.
   integer(I4P)              :: b, i, j, k                    !< Counters.
   integer(I4P)              :: d                             !< Direction counter.

   lambda_max = 0._R8P
   !$omp parallel do collapse(4) default(firstprivate) shared(dxyz, q) reduction(max:lambda_max)
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               call conservative_to_auxiliary(gamma=gamma, R=R, q=q(:,i,j,k,b), qa=qa)
               lambda_max = max(lambda_max, sum([((abs(qa(IA_U+d-1)) + qa(IA_A)) / dxyz(d,b), d=1, 3)], mask=.not.is_null))
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine compute_lambda_max
endmodule adam_flume_cpu_euler_kernels
