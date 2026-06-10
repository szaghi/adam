!< ADAM, PRISM coil source definition, FNL backend.

#include "fundal.H"

module adam_prism_fnl_coil_object
!< ADAM, PRISM coil source definition, FNL backend.

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
public :: prism_fnl_coil_object

type :: prism_fnl_coil_object
   real(R8P),    pointer :: A_gpu(:)=>null()               !< Current amplitude (A)
   real(R8P),    pointer :: f_gpu(:)=>null()               !< Current frequency, if AC (Hz)
   real(R8P),    pointer :: phase_gpu(:)=>null()           !< Current initial phase, if AC
   real(R8P),    pointer :: J_vec_gpu(:,:,:,:,:,:)=>null() !< Matrice contenente versori corrente spire (se assente = 0)
   contains
      ! public methods
      procedure, pass(self) :: copy_cpu_gpu !< Copy data from CPU to GPU.
      procedure, pass(self) :: copy_gpu_cpu !< Copy data from GPU to CPU.
      procedure, pass(self) :: initialize   !< Initialize class from global singletons.
endtype prism_fnl_coil_object

contains
   ! public methods
   subroutine copy_cpu_gpu(self, coil, grid, buf4D, buf6D, verbose)
   !< Copy data from CPU to GPU.
   class(prism_fnl_coil_object), intent(inout)           :: self         !< The field.
   class(prism_coil_object),     intent(in)              :: coil         !< Coils on host.
   type(grid_object),            intent(in)              :: grid         !< Grid (sibling realm component, threaded in).
   integer(I4P),                 intent(inout), optional :: buf4D(1:,                                &
                                                                  1-grid%ngc:,1-grid%ngc:,1-grid%ngc:&
                                                                  )      !< Buffer (host memory, device shape), rank 4D.
   real(R8P),                    intent(inout), optional :: buf6D(1:,                                 &
                                                                  1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,&
                                                                  1:,1:) !< Buffer (host memory, device shape), rank 6D.
   logical,                      intent(in),    optional :: verbose      !< Flag to activate verbose mode.
   logical                                               :: verbose_     !< Flag to activate verbose mode, local var.
   integer(I4P)                                          :: db4(2,4)     !< Device data bounds, rank 4.
   integer(I4P)                                          :: hb4(2,4)     !< Host   data bounds, rank 4.
   integer(I4P)                                          :: db6(2,6)     !< Device data bounds, rank 6.
   integer(I4P)                                          :: hb6(2,6)     !< Host   data bounds, rank 6.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (verbose_) call mpih_fnl%print_message('prism_fnl_coil_object%copy_cpu_gpu start')
   call dev_memcpy_to_device(src=coil%coil_amplitude    ,dst=self%A_gpu    )
   call dev_memcpy_to_device(src=coil%f                 ,dst=self%f_gpu    )
   call dev_memcpy_to_device(src=coil%phase             ,dst=self%phase_gpu)
   if (present(buf6D)) then
      db6(1,:) = lbound(self%j_vec_gpu ) ; db6(2,:) = ubound(self%j_vec_gpu )
      hb6(1,:) = lbound(coil%j_vec) ; hb6(2,:) = ubound(coil%j_vec)
      call dev_memcpy_to_device(bb=db6,ij=[1,5],tb=hb6,dst=self%j_vec_gpu,src=coil%j_vec,buf=buf6D)
   else
      call dev_assign_to_device(src=coil%j_vec,dst=self%j_vec_gpu,ij=[1,5])
   endif
   if (verbose_) call mpih_fnl%print_message('prism_fnl_coil_object%copy_cpu_gpu finish')
   endsubroutine copy_cpu_gpu

   subroutine copy_gpu_cpu(self, coil, grid, buf4D, buf6D, verbose)
   !< Copy data from GPU to CPU.
   class(prism_fnl_coil_object), intent(inout)           :: self         !< The field.
   class(prism_coil_object),     intent(inout)           :: coil         !< Coils on host.
   type(grid_object),            intent(in)              :: grid         !< Grid (sibling realm component, threaded in).
   integer(I4P),                 intent(inout), optional :: buf4D(1:,                                &
                                                                  1-grid%ngc:,1-grid%ngc:,1-grid%ngc:&
                                                                  )      !< Buffer (host memory, device shape), rank 4D.
   real(R8P),                    intent(inout), optional :: buf6D(1:,                                 &
                                                                  1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,&
                                                                  1:,1:) !< Buffer (host memory, device shape), rank 6D.
   logical,                      intent(in),    optional :: verbose      !< Flag to activate verbose mode.
   logical                                               :: verbose_     !< Flag to activate verbose mode, local var.
   integer(I4P)                                          :: db4(2,4)     !< Device data bounds, rank 4.
   integer(I4P)                                          :: hb4(2,4)     !< Host   data bounds, rank 4.
   integer(I4P)                                          :: db6(2,6)     !< Device data bounds, rank 6.
   integer(I4P)                                          :: hb6(2,6)     !< Host   data bounds, rank 6.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (verbose_) call mpih_fnl%print_message('prism_fnl_coil_object%copy_gpu_cpu start')
   call dev_memcpy_from_device(src=self%A_gpu    ,dst=coil%coil_amplitude)
   call dev_memcpy_from_device(src=self%f_gpu    ,dst=coil%f             )
   call dev_memcpy_from_device(src=self%phase_gpu,dst=coil%phase         )
   if (present(buf6D)) then
      db6(1,:) = lbound(self%j_vec_gpu) ; db6(2,:) = ubound(self%j_vec_gpu )
      hb6(1,:) = lbound(coil%j_vec    ) ; hb6(2,:) = ubound(coil%j_vec     )
      call dev_memcpy_from_device(bb=db6,ij=[1,5],tb=hb6,dst=coil%j_vec,src=self%j_vec_gpu,buf=buf6D)
   else
      call dev_assign_from_device(src=self%j_vec_gpu,dst=coil%j_vec,ij=[1,5])
   endif
   if (verbose_) call mpih_fnl%print_message('prism_fnl_coil_object%copy_gpu_cpu finish')
   endsubroutine copy_gpu_cpu

   subroutine initialize(self, coil, field, grid)
   !< Initialize class from program-scope `field` (adam_field_global) and `grid` (adam_grid_global) singletons.
   !< Requires `mpih_fnl` (adam_fnl_mpih_global), `field` and `grid` singletons to be ready.
   class(prism_fnl_coil_object), intent(inout)      :: self !< Coils.
   class(prism_coil_object),     intent(in), target :: coil !< Coils on host.
   type(field_object),           intent(in)       :: field !< Field (sibling realm component, threaded in).
   type(grid_object),            intent(in)       :: grid !< Grid (sibling realm component, threaded in).
   integer(I4P)                                     :: nv   !< Counter.
   integer(I4P)                                     :: ierr !< Error status.

   associate(ni=>grid%ni, nj=>grid%nj, nk=>grid%nk, ngc=>grid%ngc, nb=>field%nb, nc=>coil%total_coils_number)
   print '(A)', mpih_fnl%myrankstr//'prism_fnl_coil_object%initialize start'
   nv = size(coil%j_vec,dim=1)
   call dev_alloc(fptr_dev=self%A_gpu        ,ubounds=[nc                           ],lbounds=[0                      ],ierr=ierr)
   call dev_alloc(fptr_dev=self%f_gpu        ,ubounds=[nc                           ],lbounds=[0                      ],ierr=ierr)
   call dev_alloc(fptr_dev=self%phase_gpu    ,ubounds=[nc                           ],lbounds=[0                      ],ierr=ierr)
   call dev_alloc(fptr_dev=self%j_vec_gpu    ,ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,nv,nc],lbounds=[1,1-ngc,1-ngc,1-ngc,1,1],ierr=ierr)
   call self%copy_cpu_gpu(coil=coil, grid=grid)
   print '(A)', mpih_fnl%myrankstr//'prism_fnl_coil_object%initialize finish'
   endassociate
   endsubroutine initialize
endmodule adam_prism_fnl_coil_object
