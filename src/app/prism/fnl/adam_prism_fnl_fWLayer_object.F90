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
   contains
      ! public methods
      procedure, pass(self) :: copy_cpu_gpu !< Copy data from CPU to GPU.
      procedure, pass(self) :: copy_gpu_cpu !< Copy data from GPU to CPU.
      procedure, pass(self) :: initialize   !< Initialize object from global singletons.
endtype prism_fnl_fwlayer_object

contains
   ! public methods
   subroutine copy_cpu_gpu(self, fwlayer, grid, buffer, verbose)
   !< Copy data from CPU to GPU.
   class(prism_fnl_fwlayer_object), intent(inout)           :: self       !< The field.
   class(prism_fwlayer_object),     intent(in)              :: fwlayer    !< Fwlayer common handler (host).
   type(grid_object),               intent(in)              :: grid       !< Grid (sibling realm component, threaded in).
   real(R8P),                       intent(inout), optional :: buffer(1:,                                 &
                                                                      1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,&
                                                                      1:) !< Buffer (host memory, device shape).
   logical,                         intent(in),    optional :: verbose    !< Flag to activate verbose mode.
   logical                                                  :: verbose_   !< Flag to activate verbose mode, local var.
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (verbose_) call mpih_fnl%print_message('prism_fnl_fwlayer_object%copy_cpu_gpu start')
   if (verbose_) call mpih_fnl%print_message('prism_fnl_fwlayer_object%copy_cpu_gpu finish')
   endsubroutine copy_cpu_gpu

   subroutine copy_gpu_cpu(self, fwlayer, grid, buffer, verbose)
   !< Copy data from GPU to CPU.
   class(prism_fnl_fwlayer_object), intent(inout)           :: self       !< The field.
   class(prism_fwlayer_object),     intent(inout)           :: fwlayer    !< Fwlayer common handler (host).
   type(grid_object),               intent(in)              :: grid       !< Grid (sibling realm component, threaded in).
   real(R8P),                       intent(inout), optional :: buffer(1:,                                 &
                                                                      1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,&
                                                                      1:) !< Buffer (host memory, device shape).
   logical,                         intent(in), optional    :: verbose    !< Flag to activate verbose mode.
   logical                                                  :: verbose_   !< Flag to activate verbose mode, local var.
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (verbose_) call mpih_fnl%print_message('prism_fnl_fwlayer_object%copy_gpu_cpu start')
   if (verbose_) call mpih_fnl%print_message('prism_fnl_fwlayer_object%copy_gpu_cpu finish')
   endsubroutine copy_gpu_cpu

   subroutine initialize(self, fwlayer, field, grid)
   !< Initialize the fWLayer from program-scope `field` (adam_field_global) and `grid` (adam_grid_global) singletons.
   !< Requires `mpih_fnl` (adam_fnl_mpih_global), `field` and `grid` singletons to be ready.
   class(prism_fnl_fwlayer_object), intent(inout)      :: self    !< fWLayer.
   type(prism_fwlayer_object),      intent(in), target :: fwlayer !< Fwlayer common handler.
   type(field_object),         intent(in)         :: field !< Field (sibling realm component, threaded in).
   type(grid_object),          intent(in)         :: grid !< Grid (sibling realm component, threaded in).
   ! No device-side fWLayer state is needed: factors are computed on the fly during correction.
   print '(A)', mpih_fnl%myrankstr//'prism_fnl_fwlayer_object%initialize start'
   print '(A)', mpih_fnl%myrankstr//'prism_fnl_fwlayer_object%initialize finish'
   endsubroutine initialize

   ! non TBP
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
   !$omp OMPLOOP collapse(3) &
   !$omp& DEVICEPTR(dxyz_gpu,q_gpu) private(f_value,fm1,fp1,D_alfa,D_beta,B_alfa,B_beta,ds_b,offset) &
   !$omp& firstprivate(block_idx,ni1,ni2,nj1,nj2,nk1,nk2,n,s2,alfa_D,beta_D,alfa_B,beta_B,C_face)
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
