!< ADAM, PRISM (Plasma Research usIng Simulation Methods) fWLayer class definition, FNL backend.

#include "fundal.H"

module adam_prism_fnl_fWLayer_object
!< ADAM, PRISM (Plasma Research usIng Simulation Methods) fWLayer class definition, FNL backend.

! ADAM modules
use :: adam_field_object, only : field_object
use :: adam_fnl_library
! PRISM modules
use :: adam_prism_fWLayer_object
use :: adam_prism_parameters
! third party modules
use :: fundal
use :: penf

implicit none
private
public :: prism_fnl_fwlayer_object
public :: apply_fwl_correction_dev_kernel

type :: prism_fnl_fwlayer_object
   !< PRISM fWLayer class definition.
   ! ADAM library objects
   type(mpih_fnl_object) :: mpih_gpu !< MPI handler, FNL backend.
   ! PRISM library objects
   type(prism_fwlayer_object), pointer :: fwlayer=>null() !< Fwlayer common handler.
   ! device data
   real(R8P), pointer :: f_gpu(:,:,:,:,:)=>null() !< fWLayer function.
   contains
      ! public methods
      procedure, pass(self) :: copy_cpu_gpu !< Copy data from CPU to GPU.
      procedure, pass(self) :: copy_gpu_cpu !< Copy data from GPU to CPU.
      procedure, pass(self) :: initialize   !< Initialize object.
endtype prism_fnl_fwlayer_object

contains
   ! public methods
   subroutine copy_cpu_gpu(self, verbose)
   !< Copy data from CPU to GPU.
   class(prism_fnl_fwlayer_object), intent(inout)        :: self     !< The field.
   logical,                         intent(in), optional :: verbose  !< Flag to activate verbose mode.
   logical                                               :: verbose_ !< Flag to activate verbose mode, local var.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (verbose_) call self%mpih_gpu%print_message('prism_fnl_fwlayer_object%copy_cpu_gpu start')
   call dev_assign_to_device(src=self%fwlayer%f,dst=self%f_gpu,ij=[1,5])
   if (verbose_) call self%mpih_gpu%print_message('prism_fnl_fwlayer_object%copy_cpu_gpu finish')
   endsubroutine copy_cpu_gpu

   subroutine copy_gpu_cpu(self, verbose)
   !< Copy data from GPU to CPU.
   class(prism_fnl_fwlayer_object), intent(inout)        :: self     !< The field.
   logical,                         intent(in), optional :: verbose  !< Flag to activate verbose mode.
   logical                                               :: verbose_ !< Flag to activate verbose mode, local var.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (verbose_) call self%mpih_gpu%print_message('prism_fnl_fwlayer_object%copy_gpu_cpu start')
   call dev_assign_from_device(src=self%f_gpu,dst=self%fwlayer%f,ij=[1,5])
   if (verbose_) call self%mpih_gpu%print_message('prism_fnl_fwlayer_object%copy_gpu_cpu finish')
   endsubroutine copy_gpu_cpu

   subroutine initialize(self, field, fwlayer)
   !< Initialize the fWLayer.
   class(prism_fnl_fwlayer_object), intent(inout)      :: self    !< fWLayer.
   type(field_object),              intent(in)         :: field   !< Field.
   type(prism_fwlayer_object),      intent(in), target :: fwlayer !< Fwlayer common handler.
   integer(I4P)                                        :: ierr    !< Error status.

   call self%mpih_gpu%initialize(verbose=.true.)
   associate(ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, ngc=>field%grid%ngc, nb=>field%nb)
   print '(A)', self%mpih_gpu%myrankstr//'prism_fnl_fwlayer_object%initialize start'
   self%fwlayer => fwlayer
   call dev_alloc(fptr_dev=self%f_gpu, &
                  ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,size(fwlayer%f,dim=1)], &
                  lbounds=[1,1-ngc,1-ngc,1-ngc,1], init_value=0._R8P, ierr=ierr)
   print '(A)', self%mpih_gpu%myrankstr//'prism_fnl_fwlayer_object%initialize finish'
   endassociate
   endsubroutine initialize

   ! non TBP
   subroutine apply_fwl_correction_dev_kernel(blocks_number,ngc,ni1,ni2,nj1,nj2,nk1,nk2,n,s2,alfa_D,beta_D,alfa_B,beta_B,&
                                              f_gpu,q_gpu)
   !< Applay FWL correction, direction agnostic, device kernel.
   integer(I4P), intent(in)    :: blocks_number                     !< Blocks number.
   integer(I4P), intent(in)    :: ngc                               !< Number of ghost cells.
   integer(I4P), intent(in)    :: ni1,ni2,nj1,nj2,nk1,nk2           !< Dimensions of FWL domain.
   integer(I4P), intent(in)    :: n                                 !< f component.
   real(R8P),    intent(in)    :: s2                                !< Side coefficient.
   integer(I4P), intent(in)    :: alfa_D, beta_D                    !< Corrected var index of D (Barbas' notation).
   integer(I4P), intent(in)    :: alfa_B, beta_B                    !< Corrected var index of D (Barbas' notation).
   real(R8P),    intent(in)    :: f_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< fWLayer function values.
   real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field variables.
   real(R8P)                   :: fm1, fp1                          !< fWLayer function values in -+ cell.
   integer(I4P)                :: b,i,j,k                           !< Counter.

   !$acc parallel loop independent gang vector collapse(4) &
   !$acc& DEVICEVAR(f_gpu,q_gpu) firstprivate(n,s2,alfa_D,beta_D,alfa_B,beta_B)
   do b=1,blocks_number
   do k=nk1, nk2
   do j=nj1, nj2
   do i=ni1, ni2
      fm1 = f_gpu(b,i,j,k,n) - 1._R8P
      fp1 = f_gpu(b,i,j,k,n) + 1._R8P
      q_gpu(b,i,j,k,alfa_D) = MU0_SQ_I2  * ( s2*fm1*q_gpu(b,i,j,k,beta_B)*EPS0_SQ +    fp1*q_gpu(b,i,j,k,alfa_D)*MU0_SQ)
      q_gpu(b,i,j,k,beta_D) = MU0_SQ_I2  * (-s2*fm1*q_gpu(b,i,j,k,alfa_B)*EPS0_SQ +    fp1*q_gpu(b,i,j,k,beta_D)*MU0_SQ)
      q_gpu(b,i,j,k,alfa_B) = EPS0_SQ_I2 * (    fp1*q_gpu(b,i,j,k,alfa_B)*EPS0_SQ - s2*fm1*q_gpu(b,i,j,k,beta_D)*MU0_SQ)
      q_gpu(b,i,j,k,beta_B) = EPS0_SQ_I2 * (    fp1*q_gpu(b,i,j,k,beta_B)*EPS0_SQ + s2*fm1*q_gpu(b,i,j,k,alfa_D)*MU0_SQ)
   enddo
   enddo
   enddo
   enddo
   endsubroutine apply_fwl_correction_dev_kernel
endmodule adam_prism_fnl_fWLayer_object
