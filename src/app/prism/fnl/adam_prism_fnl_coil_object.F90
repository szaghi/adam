!< ADAM, PRISM coil source definition, FNL backend.

#include "fundal.H"

module adam_prism_fnl_coil_object
!< ADAM, PRISM coil source definition, FNL backend.

! ADAM modules
use :: adam_field_object, only : field_object
use :: adam_global_grid, only: grid
use :: adam_fnl_mpih_object
! PRSIM modules
use :: adam_prism_coil_object
use :: adam_prism_parameters
! third party modules
use :: fundal
use :: penf

implicit none
private
public :: prism_fnl_coil_object

type :: prism_fnl_coil_object
   !< ADAM, PRISM coil source definition, FNL (GPU) backend.
   type(mpih_fnl_object)            :: mpih         !< MPI handler.
   type(prism_coil_object), pointer :: coil=>null() !< Coil common handler.
   ! device data
   real(R8P),    pointer :: A_gpu(:)=>null()               !< Current amplitude (A)
   real(R8P),    pointer :: f_gpu(:)=>null()               !< Current frequency, if AC (Hz)
   real(R8P),    pointer :: phase_gpu(:)=>null()           !< Current initial phase, if AC
   real(R8P),    pointer :: J_vec_gpu(:,:,:,:,:,:)=>null() !< Matrice contenente versori corrente spire (se assente = 0)
   integer(I4P), pointer :: coil_flag_gpu(:,:,:,:)=>null() !< Matrice contenente informazioni su quale spira pass.
   ! grid/field data replica for easy handling
   integer(I4P), pointer :: blocks_number=>null() !< Actual blocks number.
   integer(I4P), pointer :: nb=>null()            !< Total blocks number for MPI.
   integer(I4P), pointer :: ngc=>null()           !< Number of ghost cells.
   integer(I4P), pointer :: ni=>null()            !< Number of cells in i direction.
   integer(I4P), pointer :: nj=>null()            !< Number of cells in j direction.
   integer(I4P), pointer :: nk=>null()            !< Number of cells in k direction.
   contains
      ! public methods
      procedure, pass(self) :: copy_cpu_gpu !< Copy data from CPU to GPU.
      procedure, pass(self) :: copy_gpu_cpu !< Copy data from GPU to CPU.
      procedure, pass(self) :: initialize   !< Initialize class.
endtype prism_fnl_coil_object

contains
   ! public methods
   subroutine copy_cpu_gpu(self, verbose)
   !< Copy data from CPU to GPU.
   class(prism_fnl_coil_object), intent(inout)        :: self     !< The field.
   logical,                      intent(in), optional :: verbose  !< Flag to activate verbose mode.
   logical                                            :: verbose_ !< Flag to activate verbose mode, local var.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (verbose_) call self%mpih%print_message('prism_fnl_coil_object%copy_cpu_gpu start')
   call dev_assign_to_device(src=self%coil%A        ,dst=self%A_gpu                 )
   call dev_assign_to_device(src=self%coil%f        ,dst=self%f_gpu                 )
   call dev_assign_to_device(src=self%coil%phase    ,dst=self%phase_gpu             )
   call dev_assign_to_device(src=self%coil%j_vec    ,dst=self%j_vec_gpu    ,ij=[1,5])
   call dev_assign_to_device(src=self%coil%coil_flag,dst=self%coil_flag_gpu,ij=[1,4])
   if (verbose_) call self%mpih%print_message('prism_fnl_coil_object%copy_cpu_gpu finish')
   endsubroutine copy_cpu_gpu

   subroutine copy_gpu_cpu(self, verbose)
   !< Copy data from GPU to CPU.
   class(prism_fnl_coil_object), intent(inout)        :: self     !< The field.
   logical,                      intent(in), optional :: verbose  !< Flag to activate verbose mode.
   logical                                            :: verbose_ !< Flag to activate verbose mode, local var.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (verbose_) call self%mpih%print_message('prism_fnl_coil_object%copy_gpu_cpu start')
   call dev_assign_from_device(src=self%A_gpu        ,dst=self%coil%A                 )
   call dev_assign_from_device(src=self%f_gpu        ,dst=self%coil%f                 )
   call dev_assign_from_device(src=self%phase_gpu    ,dst=self%coil%phase             )
   call dev_assign_from_device(src=self%j_vec_gpu    ,dst=self%coil%j_vec    ,ij=[1,5])
   call dev_assign_from_device(src=self%coil_flag_gpu,dst=self%coil%coil_flag,ij=[1,4])
   if (verbose_) call self%mpih%print_message('prism_fnl_coil_object%copy_gpu_cpu finish')
   endsubroutine copy_gpu_cpu

   subroutine initialize(self, field, coil)
   !< Initialize class.
   class(prism_fnl_coil_object), intent(inout)      :: self  !< Coils.
   type(field_object),           intent(in), target :: field !< Field.
   class(prism_coil_object),     intent(in), target :: coil  !< Coils on host.
   integer(I4P)                                     :: nv    !< Counter.
   integer(I4P)                                     :: ierr  !< Error status.

   call self%mpih%initialize(do_mpi_init=.false.)
   associate(ni=>grid%ni,nj=>grid%nj,nk=>grid%nk,ngc=>grid%ngc,nb=>field%nb,nc=>coil%total_coils_number)
   print '(A)', self%mpih%myrankstr//'prism_fnl_coil_object%initialize start'
   self%coil          => coil
   self%blocks_number => field%blocks_number
   self%nb            => field%nb
   self%ngc           => field%ngc
   self%ni            => field%ni
   self%nj            => field%nj
   self%nk            => field%nk
   nv = size(self%coil%j_vec,dim=1)
   call dev_alloc(fptr_dev=self%j_vec_gpu    ,ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,nv,nc],lbounds=[1,1-ngc,1-ngc,1-ngc,1,1],ierr=ierr)
   call dev_alloc(fptr_dev=self%coil_flag_gpu,ubounds=[nb,ni+ngc,nj+ngc,nk+ngc      ],lbounds=[1,1-ngc,1-ngc,1-ngc    ],ierr=ierr)
   call self%copy_cpu_gpu
   print '(A)', self%mpih%myrankstr//'prism_fnl_coil_object%initialize finish'
   endassociate
   endsubroutine initialize
endmodule adam_prism_fnl_coil_object
