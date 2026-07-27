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
   subroutine apply_fwl_correction_dev_kernel(block_idx,ngc,ni1,ni2,nj1,nj2,nk1,nk2,n,s2,alfa_D,beta_D,alfa_B,beta_B,&
                                              C_face,dxyz_gpu,q_gpu)
   !< Applay FWL correction, direction agnostic, device kernel.
   integer(I4P), intent(in)    :: block_idx                         !< Block index.
   integer(I4P), intent(in)    :: ngc                               !< Number of ghost cells.
   integer(I4P), intent(in)    :: ni1,ni2,nj1,nj2,nk1,nk2           !< Dimensions of FWL domain.
   integer(I4P), intent(in)    :: n                                 !< f component.
   real(R8P),    intent(in)    :: s2                                !< Side coefficient.
   integer(I4P), intent(in)    :: alfa_D, beta_D                    !< Corrected var index of D (Barbas' notation).
   integer(I4P), intent(in)    :: alfa_B, beta_B                    !< Corrected var index of D (Barbas' notation).
   integer(I4P), intent(in)    :: C_face                            !< Layer thickness in cells on this block/face.
   real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)                   !< Block mesh spacing [nb,3].
   real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field variables.
   real(R8P)                   :: ds_b                              !< Block mesh spacing along the layer-normal direction.
   real(R8P)                   :: f_value                           !< Local fWLayer factor.
   real(R8P)                   :: D_alfa, D_beta                    !< components of tangential fields before correction
   real(R8P)                   :: B_alfa, B_beta                    !< components of tangential fields before correction
   real(R8P)                   :: fm1, fp1                          !< fWLayer function values in -+ cell.
   integer(I4P)                :: i,j,k                             !< Counter.
   integer(I4P)                :: offset                            !< Cell offset from the active boundary.

   !$acc parallel loop independent gang vector collapse(3) &
   !$acc& DEVICEVAR(dxyz_gpu,q_gpu) private(f_value,fm1,fp1,D_alfa,D_beta,B_alfa,B_beta,ds_b,offset) &
   !$acc& firstprivate(block_idx,ni1,ni2,nj1,nj2,nk1,nk2,n,s2,alfa_D,beta_D,alfa_B,beta_B,C_face)
   do k=nk1, nk2
   do j=nj1, nj2
   do i=ni1, ni2
      ds_b = dxyz_gpu(block_idx,n)
      select case(n)
      case(1)
         if (s2 > 0._R8P) then
            offset = i - ni1
         else
            offset = ni2 - i
         endif
      case(2)
         if (s2 > 0._R8P) then
            offset = j - nj1
         else
            offset = nj2 - j
         endif
      case default
         if (s2 > 0._R8P) then
            offset = k - nk1
         else
            offset = nk2 - k
         endif
      endselect
      f_value = compute_fwl_factor(offset=offset, cells_number=C_face, ds=ds_b)
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
