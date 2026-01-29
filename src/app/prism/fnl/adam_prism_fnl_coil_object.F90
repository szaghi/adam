!< ADAM, PRISM coil source definition, FNL backend.

#include "fundal.H"

module adam_prism_fnl_coil_object
!< ADAM, PRISM coil source definition, FNL backend.

! PRSIM modules
use adam_prism_coil_object
use adam_prism_parameters
! ADAM modules
use adam_fnl_mpih_object
! third party modules
use fundal
use penf

implicit none
private
public :: prism_fnl_coil_object
public :: compute_coils_current_dev

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
   integer(I4P)                                       :: nb,ngc,ni,nj,nk,nv   !< Grid dimensions.
   integer(I4P)                                       :: i, j, k, b           !< Counter.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (verbose_) call self%mpih%print_message('prism_fnl_coil_object%copy_cpu_gpu start')
   call dev_assign_to_device(src=self%coil%A,         dst=self%A_gpu    )
   call dev_assign_to_device(src=self%coil%f,         dst=self%f_gpu    )
   call dev_assign_to_device(src=self%coil%phase,     dst=self%phase_gpu)
   call dev_assign_to_device(src=self%coil%d,         dst=self%d_gpu    )
   ! call dev_assign_to_device(src=self%coil%J_vec,     dst=self%J_vec_gpu,     transposed=.true.)
   ! call dev_assign_to_device(src=self%coil%coil_flag, dst=self%coil_flag_gpu, transposed=.true.)
   nv  =  size(  self%coil%j_vec,     dim=1)
   ngc = -lbound(self%coil%coil_flag, dim=1) + 1
   ni  =  ubound(self%coil%coil_flag, dim=1) - ngc
   nj  =  ubound(self%coil%coil_flag, dim=2) - ngc
   nk  =  ubound(self%coil%coil_flag, dim=3) - ngc
   nb  =  size(  self%coil%coil_flag, dim=4)
   allocate(j_vec_t(    1:nb,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nv))
   allocate(coil_flag_t(1:nb,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc     ))
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
   integer(I4P)                                     :: nv            !< J_vec variables number.
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
   nv = size(  self%coil%j_vec, dim=1)
   call dev_alloc(fptr_dev=self%j_vec_gpu,     ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,nv], lbounds=[1,1-ngc,1-ngc,1-ngc,1], ierr=ierr)
   call dev_alloc(fptr_dev=self%coil_flag_gpu, ubounds=[nb,ni+ngc,nj+ngc,nk+ngc   ], lbounds=[1,1-ngc,1-ngc,1-ngc  ], ierr=ierr)
   call self%copy_cpu_gpu
   print '(A)', self%mpih%myrankstr//'prism_fnl_coil_object%initialize finish'
   endsubroutine initialize

   ! non TBP
   subroutine compute_coils_current_dev(ni, nj, nk, ngc, blocks_number,                               &
                                        time_s, td, A_gpu, f_gpu, phase_gpu, coil_flag_gpu, j_vec_gpu,&
                                        var_Jx, var_Jy, var_Jz, q_gpu)
   !< Compute current coils sources, device kernel.
   integer(I4P), intent(in)    :: ni                                     !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                                     !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                                     !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                                    !< Ghost grid number.
   integer(I4P), intent(in)    :: blocks_number                          !< Number of blocks.
   real(R8P),    intent(in)    :: time_s                                 !< Local time.
   real(R8P),    intent(in)    :: td                                     !< Delay coil start.
   real(R8P),    intent(in)    :: A_gpu(0:)                              !< Current amplitude (A)
   real(R8P),    intent(in)    :: f_gpu(0:)                              !< Current frequency, if AC (Hz)
   real(R8P),    intent(in)    :: phase_gpu(0:)                          !< Current initial phase, if AC
   real(R8P),    intent(in)    :: J_vec_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Matrice contenente versori corrente spire.
   integer(I4P), intent(in)    :: coil_flag_gpu(1:,1-ngc:,1-ngc:,1-ngc:) !< Matrice contenente informazioni su quale spira pass.
   integer(I4P), intent(in)    :: var_Jx, var_Jy, var_Jz                 !< Indices of current density components in q vector.
   real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Conservative variables on GPU.
   real(R8P)                   :: current_density                        !< Current density.
   real(R8P)                   :: g                                      !< Starting polynomial transitory of coils.
   integer(I4P)                :: w_, w_c_                               !< Step function coeff to avoid if in parallel regions.
   real(R8P)                   :: g_, f_                                 !< Current coefficients.
   integer(I4P)                :: coil_id                                !< Uniq coild ID.
   integer(I4P)                :: i,j,k,b                                !< Counter.

   g = 10._R8P*(time_s/td)**3 - 15._R8P*(time_s/td)**4 + 6._R8P*(time_s/td)**5
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,coil_flag_gpu,A_gpu,f_gpu,phase_gpu,j_vec_gpu)
   do b=1, blocks_number
   do k=1, nk
   do j=1, nj
   do i=1, ni
      coil_id = coil_flag_gpu(b,i,j,k)

      w_   = nint(sign(1._R8P,td-time_s) + 1._R8P)/2     ! = 1 if td>time, = 0                              if td<time
      w_c_ = 1_I4P - w_                                  ! = 0 if td>time, = 1                              if td<time
      g_   = w_ * g + w_c_                               ! = g if td>time, = 1                              if td<time
      f_   = w_c_ * 2._R8P*PI*f_gpu(coil_id)*(time_s-td) ! = 0 if td>time, = 2._R8P*PI*f(coil_id)*(time-td) if td<time
      current_density = g_ * A_gpu(coil_id) * cos(f_ + phase_gpu(coil_id)*PI/180.0_R8P)*j_vec_gpu(b,i,j,k,4)

      q_gpu(b,i,j,k,VAR_JX) = current_density * j_vec_gpu(b,i,j,k,1)
      q_gpu(b,i,j,k,VAR_JY) = current_density * j_vec_gpu(b,i,j,k,2)
      q_gpu(b,i,j,k,VAR_JZ) = current_density * j_vec_gpu(b,i,j,k,3)
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_coils_current_dev
endmodule adam_prism_fnl_coil_object
