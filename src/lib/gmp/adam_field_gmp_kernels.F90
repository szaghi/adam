!< ADAM, field class GMP kernels (GMP backend of [[field_object]]).
module adam_field_gmp_kernels
!< ADAM, field class NVF kernels (NVF backend of [[field_object]]).

use penf, only : I4P, R8P

implicit none
private
public :: compute_normL2_residuals_gmp

contains
   subroutine compute_normL2_residuals_gmp(ni, nj, nk, ngc, nv, blocks_number, dq_gpu, norm)
   !< Compute L2 norm of residuals.
   integer(I4P), intent(in)    :: ni                                 !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                                 !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                                 !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                                !< Ghost grid number.
   integer(I4P), intent(in)    :: nv                                 !< Number of states variables.
   integer(I4P), intent(in)    :: blocks_number                      !< Number of blocks.
   real(R8P),    intent(in)    :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Residuals.
   real(R8P),    intent(inout) :: norm(1:)                           !< Residuals norm.
   real(R8P)                   :: norm_gpu                           !< Residuals norm, local scalar buffer for reduction.
   integer(I4P)                :: b, i, j, k, v                      !< Counter.

   do v=1, nv
      norm_gpu = 0._R8P
      !$omp target teams distribute parallel do collapse(4) has_device_addr(dq_gpu) reduction(+:norm_gpu)
      do k=1, nk
         do j=1, nj
            do i=1, ni
               do b=1, blocks_number
                  norm_gpu = norm_gpu + dq_gpu(b, i, j, k, v)**2
               enddo
            enddo
         enddo
      enddo
      norm(v) = norm_gpu
   enddo
   endsubroutine compute_normL2_residuals_gmp
endmodule adam_field_gmp_kernels
