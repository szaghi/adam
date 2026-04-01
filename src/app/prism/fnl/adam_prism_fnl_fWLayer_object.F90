!< ADAM, PRISM (Plasma Research usIng Simulation Methods) fWLayer class definition, FNL backend.

#include "fundal.H"

module adam_prism_fnl_fWLayer_object
!< ADAM, PRISM (Plasma Research usIng Simulation Methods) fWLayer class definition, FNL backend.

! ADAM classes, libraries, parameters
use :: adam_common_library
! ADAM FNL classes, libraries, parameters
use :: adam_fnl_library
! PRISM common classes, libraries, parameters
use :: adam_prism_common_library
! third party modules
use :: fundal
use :: penf

implicit none
private
public :: prism_fnl_fwlayer_object
public :: apply_fwl_correction_dev_kernel

type :: prism_fnl_fwlayer_object
   !< PRISM fWLayer class definition.
   real(R8P), pointer :: f_gpu(:,:,:,:,:)=>null() !< fWLayer function.
   contains
      ! public methods
      procedure, pass(self) :: copy_cpu_gpu !< Copy data from CPU to GPU.
      procedure, pass(self) :: copy_gpu_cpu !< Copy data from GPU to CPU.
      procedure, pass(self) :: initialize   !< Initialize object from global singletons.
endtype prism_fnl_fwlayer_object

contains
   ! public methods
   subroutine copy_cpu_gpu(self, buffer, verbose)
   !< Copy data from CPU to GPU.
   class(prism_fnl_fwlayer_object), intent(inout)           :: self       !< The field.
   real(R8P),                       intent(inout), optional :: buffer(1:,                                 &
                                                                      1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,&
                                                                      1:) !< Buffer (host memory, device shape).
   logical,                         intent(in),    optional :: verbose    !< Flag to activate verbose mode.
   logical                                                  :: verbose_   !< Flag to activate verbose mode, local var.
   integer(I4P)                                             :: db(2,5)    !< Device data bounds.
   integer(I4P)                                             :: hb(2,5)    !< Host   data bounds.

   if (fwlayer%C ==0) return
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (verbose_) call mpih_fnl%print_message('prism_fnl_fwlayer_object%copy_cpu_gpu start')
   if (present(buffer)) then
      db(1,:) = lbound(self%f_gpu) ; db(2,:) = ubound(self%f_gpu)
      hb(1,:) = lbound(fwlayer%f ) ; hb(2,:) = ubound(fwlayer%f )
      call dev_memcpy_to_device(bb=db,ij=[1,5],tb=hb,dst=self%f_gpu,src=fwlayer%f,buf=buffer)
   else
      call dev_assign_to_device(src=fwlayer%f,dst=self%f_gpu,ij=[1,5])
   endif
   if (verbose_) call mpih_fnl%print_message('prism_fnl_fwlayer_object%copy_cpu_gpu finish')
   endsubroutine copy_cpu_gpu

   subroutine copy_gpu_cpu(self, buffer, verbose)
   !< Copy data from GPU to CPU.
   class(prism_fnl_fwlayer_object), intent(inout)           :: self       !< The field.
   real(R8P),                       intent(inout), optional :: buffer(1:,                                 &
                                                                      1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,&
                                                                      1:) !< Buffer (host memory, device shape).
   logical,                         intent(in), optional    :: verbose    !< Flag to activate verbose mode.
   logical                                                  :: verbose_   !< Flag to activate verbose mode, local var.
   integer(I4P)                                             :: db(2,5)    !< Device data bounds.
   integer(I4P)                                             :: hb(2,5)    !< Host   data bounds.

   if (fwlayer%C ==0) return
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (verbose_) call mpih_fnl%print_message('prism_fnl_fwlayer_object%copy_gpu_cpu start')
   if (present(buffer)) then
      db(1,:) = lbound(self%f_gpu) ; db(2,:) = ubound(self%f_gpu)
      hb(1,:) = lbound(fwlayer%f ) ; hb(2,:) = ubound(fwlayer%f )
      call dev_memcpy_from_device(bb=db,ij=[1,5],tb=hb,src=self%f_gpu,dst=fwlayer%f,buf=buffer)
   else
      call dev_assign_from_device(src=self%f_gpu,dst=fwlayer%f,ij=[1,5])
   endif
   if (verbose_) call mpih_fnl%print_message('prism_fnl_fwlayer_object%copy_gpu_cpu finish')
   endsubroutine copy_gpu_cpu

   subroutine initialize(self, fwlayer)
   !< Initialize the fWLayer from program-scope `field` (adam_field_global) and `grid` (adam_grid_global) singletons.
   !< Requires `mpih_fnl` (adam_fnl_mpih_global), `field` and `grid` singletons to be ready.
   class(prism_fnl_fwlayer_object), intent(inout)      :: self    !< fWLayer.
   type(prism_fwlayer_object),      intent(in), target :: fwlayer !< Fwlayer common handler.
   integer(I4P)                                        :: ierr    !< Error status.

   if (fwlayer%C ==0) return
   associate(ni=>grid%ni, nj=>grid%nj, nk=>grid%nk, ngc=>grid%ngc, nb=>field%nb)
   print '(A)', mpih_fnl%myrankstr//'prism_fnl_fwlayer_object%initialize start'
   call dev_alloc(fptr_dev=self%f_gpu, &
                  ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,size(fwlayer%f,dim=1)], &
                  lbounds=[1,1-ngc,1-ngc,1-ngc,1], init_value=0._R8P, ierr=ierr)
   print '(A)', mpih_fnl%myrankstr//'prism_fnl_fwlayer_object%initialize finish'
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

   if (fwlayer%C ==0) return
   !$acc parallel loop independent gang vector collapse(4) &
   !$acc& DEVICEVAR(f_gpu,q_gpu) private(fm1,fp1)          &
   !$acc& firstprivate(ni1,ni2,nj1,nj2,nk1,nk2,blocks_number,n,s2,alfa_D,beta_D,alfa_B,beta_B)
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
