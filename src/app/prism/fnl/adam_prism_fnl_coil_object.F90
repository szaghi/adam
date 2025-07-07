!< ADAM, PRISM coil source definition, FNL (GPU) backend.
module adam_prism_fnl_coil_object
    !< ADAM, PRISM coil source definition, FNL (GPU) backend.

#include "fundal.H"

use adam_prism_coil_object
use adam_mpih_fnl_object
use fundal
use penf

implicit none
private
public :: prism_fnl_coil_object

type :: prism_fnl_coil_object
   !< ADAM, PRISM coil source definition, FNL (GPU) backend.
   type(mpih_fnl_object)            :: mpih                           !< MPI handler.
   type(prism_coil_object), pointer :: coil=>null()                   !< Coil common handler.
   real(R8P),               pointer :: A_gpu(:)=>null()               !< Current amplitude (A)
   real(R8P),               pointer :: f_gpu(:)=>null()               !< Current frequency, if AC (Hz)
   real(R8P),               pointer :: phase_gpu(:)=>null()           !< Current initial phase, if AC
   real(R8P),               pointer :: d_gpu(:)=>null()               !< Coil wire diameter
   real(R8P),               pointer :: J_vec_gpu(:,:,:,:,:)=>null()   !< Matrice contenente versori corrente spire (se assente = 0)
   integer(I4P),            pointer :: coil_flag_gpu(:,:,:,:)=>null() !< Matrice contenente informazioni su quale spira pass.
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
   class(prism_fnl_coil_object), intent(inout)        :: self                 !< The field.
   logical,                      intent(in), optional :: verbose              !< Flag to activate verbose mode.
   logical                                            :: verbose_             !< Flag to activate verbose mode, local var.
   real(R8P),    allocatable                          :: j_vec_t(:,:,:,:,:)   !< Transposed j_vec.
   integer(I4P), allocatable                          :: coil_flag_t(:,:,:,:) !< Transposed coil flag.
   integer(I4P)                                       :: nb,ngc,ni,nj,nk      !< Grid dimensions.
   integer(I4P)                                       :: i, j, k, b           !< Counter.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (verbose_) call self%mpih%print_message('prism_fnl_coil_object%copy_cpu_gpu start')
   call dev_assign_to_device(src=self%coil%A,         dst=self%A_gpu    )
   call dev_assign_to_device(src=self%coil%f,         dst=self%f_gpu    )
   call dev_assign_to_device(src=self%coil%phase,     dst=self%phase_gpu)
   call dev_assign_to_device(src=self%coil%d,         dst=self%d_gpu    )
   ! call dev_assign_to_device(src=self%coil%J_vec,     dst=self%J_vec_gpu,     transposed=.true.)
   ! call dev_assign_to_device(src=self%coil%coil_flag, dst=self%coil_flag_gpu, transposed=.true.)
   ngc = -lbound(self%coil%coil_flag, dim=1) + 1
   ni  =  ubound(self%coil%coil_flag, dim=1) - ngc
   nj  =  ubound(self%coil%coil_flag, dim=2) - ngc
   nk  =  ubound(self%coil%coil_flag, dim=3) - ngc
   nb  =  size(  self%coil%coil_flag, dim=4)
   allocate(j_vec_t(    1:nb,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:3))
   allocate(coil_flag_t(1:nb,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc    ))
   do b=1, nb
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               j_vec_t(    b,i,j,k,:) = self%coil%j_vec(    :,i,j,k,b)
               coil_flag_t(b,i,j,k  ) = self%coil%coil_flag(  i,j,k,b)
            enddo
         enddo
      enddo
   enddo
   call dev_memcpy_to_device(dst=self%j_vec_gpu,     src=j_vec_t    )
   call dev_memcpy_to_device(dst=self%coil_flag_gpu, src=coil_flag_t)
   if (verbose_) call self%mpih%print_message('prism_fnl_coil_object%copy_cpu_gpu finish')
   endsubroutine copy_cpu_gpu

   subroutine copy_gpu_cpu(self, verbose)
   !< Copy data from GPU to CPU.
   class(prism_fnl_coil_object), intent(inout)        :: self                 !< The field.
   logical,                      intent(in), optional :: verbose              !< Flag to activate verbose mode.
   logical                                            :: verbose_             !< Flag to activate verbose mode, local var.
   real(R8P),    allocatable                          :: j_vec_t(:,:,:,:,:)   !< Transposed j_vec.
   integer(I4P), allocatable                          :: coil_flag_t(:,:,:,:) !< Transposed coil flag.
   integer(I4P)                                       :: nb,ngc,ni,nj,nk,nv   !< Grid dimensions.
   integer(I4P)                                       :: i, j, k, b           !< Counter.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (verbose_) call self%mpih%print_message('prism_fnl_coil_object%copy_gpu_cpu start')

   nv  =  size(  self%coil%j_vec, dim=1)
   ngc = -lbound(self%coil%j_vec, dim=2) + 1
   ni  =  ubound(self%coil%j_vec, dim=2) - ngc
   nj  =  ubound(self%coil%j_vec, dim=3) - ngc
   nk  =  ubound(self%coil%j_vec, dim=4) - ngc
   nb  =  size(  self%coil%j_vec, dim=5)
   allocate(j_vec_t(    1:nb,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nv))
   allocate(coil_flag_t(1:nb,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc     ))
   call dev_memcpy_from_device(dst=j_vec_t,     src=self%j_vec_gpu    )
   call dev_memcpy_from_device(dst=coil_flag_t, src=self%coil_flag_gpu)
   do b=1, nb
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               self%coil%j_vec(    :,i,j,k,b) = j_vec_t(    b,i,j,k,:)
               self%coil%coil_flag(  i,j,k,b) = coil_flag_t(b,i,j,k  )
            enddo
         enddo
      enddo
   enddo
   if (verbose_) call self%mpih%print_message('prism_fnl_coil_object%copy_gpu_cpu finish')
   endsubroutine copy_gpu_cpu

   subroutine initialize(self, coil, blocks_number, nb, ngc, ni, nj, nk)
   !< Initialize class.
   class(prism_fnl_coil_object), intent(inout)      :: self          !< Coils.
   class(prism_coil_object),     intent(in), target :: coil          !< Coils on host.
   integer(I4P),                 intent(in), target :: blocks_number !< Actual blocks number.
   integer(I4P),                 intent(in), target :: nb            !< Maximum blocks number.
   integer(I4P),                 intent(in), target :: ngc           !< Number of ghost cells.
   integer(I4P),                 intent(in), target :: ni            !< Number of cells in i direction.
   integer(I4P),                 intent(in), target :: nj            !< Number of cells in j direction.
   integer(I4P),                 intent(in), target :: nk            !< Number of cells in k direction.
   integer(I4P)                                     :: ierr          !< Error status.

   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'prism_fnl_coil_object%initialize start'
   self%coil          => coil
   self%blocks_number => blocks_number
   self%nb            => nb
   self%ngc           => ngc
   self%ni            => ni
   self%nj            => nj
   self%nk            => nk
   call dev_alloc(fptr_dev=self%j_vec_gpu,     ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,3], lbounds=[1,1-ngc,1-ngc,1-ngc,1], ierr=ierr)
   call dev_alloc(fptr_dev=self%coil_flag_gpu, ubounds=[nb,ni+ngc,nj+ngc,nk+ngc  ], lbounds=[1,1-ngc,1-ngc,1-ngc  ], ierr=ierr)
   call self%copy_cpu_gpu
   print '(A)', self%mpih%myrankstr//'prism_fnl_coil_object%initialize finish'
   endsubroutine initialize
endmodule adam_prism_fnl_coil_object
