!< ADAM, PRISM (Plasma Research usIng Simulation Methods) fWLayer class definition, FNL backend.

#include "fundal.H"

module adam_prism_fnl_fWLayer_object
!< ADAM, PRISM (Plasma Research usIng Simulation Methods) fWLayer class definition, FNL backend.

! PRISM common classes, libraries, parameters
use :: adam_prism_common_library
! third party modules
use :: fundal
use :: penf

implicit none
private
public :: apply_fwl_correction_dev_kernel

contains
   subroutine apply_fwl_correction_dev_kernel(block_idx, ngc, ni, nj, nk, ni1, ni2, nj1, nj2, nk1, nk2, face, n, s2,  &
                                              alfa_D, beta_D, alfa_B, beta_B, domain_emin_n, domain_emax_n,           &
                                              profile_extent, profile_cells, x_cell_gpu, y_cell_gpu, z_cell_gpu, dxyz_gpu, q_gpu)
   !< Applay FWL correction, direction agnostic, device kernel.
   integer(I4P), intent(in)    :: block_idx                         !< Block index.
   integer(I4P), intent(in)    :: ngc                               !< Number of ghost cells.
   integer(I4P), intent(in)    :: ni, nj, nk                        !< Cells number in i, j and k directions.
   integer(I4P), intent(in)    :: ni1,ni2,nj1,nj2,nk1,nk2           !< Dimensions of FWL domain.
   integer(I4P), intent(in)    :: face                              !< Boundary face identifier.
   integer(I4P), intent(in)    :: n                                 !< f component.
   real(R8P),    intent(in)    :: s2                                !< Side coefficient.
   integer(I4P), intent(in)    :: alfa_D, beta_D                    !< Corrected var index of D (Barbas' notation).
   integer(I4P), intent(in)    :: alfa_B, beta_B                    !< Corrected var index of D (Barbas' notation).
   real(R8P),    intent(in)    :: domain_emin_n                     !< Domain minimum coordinate along the layer normal.
   real(R8P),    intent(in)    :: domain_emax_n                     !< Domain maximum coordinate along the layer normal.
   real(R8P),    intent(in)    :: profile_extent                    !< Discrete face-wise profile extent.
   integer(I4P), intent(in)    :: profile_cells                     !< Effective face-wise layer thickness in cells.
   real(R8P),    intent(in)    :: x_cell_gpu(1:,1-ngc:)             !< Cells x coordinates on GPU.
   real(R8P),    intent(in)    :: y_cell_gpu(1:,1-ngc:)             !< Cells y coordinates on GPU.
   real(R8P),    intent(in)    :: z_cell_gpu(1:,1-ngc:)             !< Cells z coordinates on GPU.
   real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)                   !< Block mesh spacing [nb,3].
   real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field variables.
   real(R8P)                   :: center_distance                   !< Physical distance from the boundary-cell center.
   real(R8P)                   :: f_value                           !< Local fWLayer factor.
   real(R8P)                   :: D_alfa, D_beta                    !< components of tangential fields before correction
   real(R8P)                   :: B_alfa, B_beta                    !< components of tangential fields before correction
   real(R8P)                   :: fm1, fp1                          !< fWLayer function values in -+ cell.
   integer(I4P)                :: i,j,k                             !< Counter.

   !$acc parallel loop independent gang vector collapse(3) &
   !$acc& DEVICEVAR(x_cell_gpu,y_cell_gpu,z_cell_gpu,dxyz_gpu,q_gpu) private(center_distance,f_value,fm1,fp1,D_alfa,D_beta,B_alfa,B_beta) &
   !$acc& firstprivate(block_idx,ni,nj,nk,ni1,ni2,nj1,nj2,nk1,nk2,face,n,s2,alfa_D,beta_D,alfa_B,beta_B,domain_emin_n,domain_emax_n, &
   !$acc& profile_extent, profile_cells)
   !$omp OMPLOOP collapse(3) &
   !$omp& DEVICEPTR(x_cell_gpu,y_cell_gpu,z_cell_gpu,dxyz_gpu,q_gpu) private(center_distance,f_value,fm1,fp1,D_alfa,D_beta,B_alfa,B_beta) &
   !$omp& firstprivate(block_idx,ni,nj,nk,ni1,ni2,nj1,nj2,nk1,nk2,face,n,s2,alfa_D,beta_D,alfa_B,beta_B,domain_emin_n,domain_emax_n, &
   !$omp& profile_extent, profile_cells)
   do k=nk1, nk2
   do j=nj1, nj2
   do i=ni1, ni2
      select case(face)
      case(1)
         center_distance = x_cell_gpu(block_idx,i) - (domain_emin_n + 0.5_R8P * dxyz_gpu(block_idx,1))
      case(2)
         center_distance = (domain_emax_n - 0.5_R8P * dxyz_gpu(block_idx,1)) - x_cell_gpu(block_idx,i)
      case(3)
         center_distance = y_cell_gpu(block_idx,j) - (domain_emin_n + 0.5_R8P * dxyz_gpu(block_idx,2))
      case(4)
         center_distance = (domain_emax_n - 0.5_R8P * dxyz_gpu(block_idx,2)) - y_cell_gpu(block_idx,j)
      case(5)
         center_distance = z_cell_gpu(block_idx,k) - (domain_emin_n + 0.5_R8P * dxyz_gpu(block_idx,3))
      case default
         center_distance = (domain_emax_n - 0.5_R8P * dxyz_gpu(block_idx,3)) - z_cell_gpu(block_idx,k)
      endselect
      f_value = compute_fwl_factor(center_distance=center_distance, profile_extent=profile_extent, profile_cells=profile_cells)
      fm1 = f_value - 1._R8P
      fp1 = f_value + 1._R8P
      D_alfa = q_gpu(block_idx,i,j,k,alfa_D)
      D_beta = q_gpu(block_idx,i,j,k,beta_D)
      B_alfa = q_gpu(block_idx,i,j,k,alfa_B)
      B_beta = q_gpu(block_idx,i,j,k,beta_B)
      q_gpu(block_idx,i,j,k,alfa_D) = MU0_SQ_I2  * ( s2*fm1*B_beta*EPS0_SQ +    fp1*D_alfa*MU0_SQ)
      q_gpu(block_idx,i,j,k,beta_D) = MU0_SQ_I2  * (-s2*fm1*B_alfa*EPS0_SQ +    fp1*D_beta*MU0_SQ)
      q_gpu(block_idx,i,j,k,alfa_B) = EPS0_SQ_I2 * (    fp1*B_alfa*EPS0_SQ - s2*fm1*D_beta*MU0_SQ)
      q_gpu(block_idx,i,j,k,beta_B) = EPS0_SQ_I2 * (    fp1*B_beta*EPS0_SQ + s2*fm1*D_alfa*MU0_SQ)
   enddo
   enddo
   enddo
   endsubroutine apply_fwl_correction_dev_kernel
endmodule adam_prism_fnl_fWLayer_object
