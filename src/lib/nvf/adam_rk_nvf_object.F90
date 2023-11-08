!< ADAM, RK class NVF (NVF backend of [[rk_object]]).
module adam_rk_nvf_object
!< ADAM, RK class NVF (NVF backend of [[rk_object]]).

use adam_rk_object
use adam_memory_nvf_library
use adam_mpih_nvf_object
use adam_rk_nvf_kernels
use penf

implicit none
save
private
public :: rk_nvf_object

type :: rk_nvf_object
   !< RK NVF class definition.
   ! ADAM library objects
   type(rk_object), pointer :: rk=>null() !< RK common handler.
   ! ADAM NVF library objects
   type(mpih_nvf_object) :: mpih !< MPI handler.
   ! device data
   real(R8P), allocatable, device :: alph_gpu(:,:)         !< RK alpha coefficients.
   real(R8P), allocatable, device :: beta_gpu(:)           !< RK beta coefficients.
   real(R8P), allocatable, device :: gamm_gpu(:)           !< RK gamma coefficients.
   real(R8P), allocatable, device :: q_rk_gpu(:,:,:,:,:,:) !< Field cell centered variables,RK stages.
   contains
      ! public methods
      procedure, pass(self) :: assign_stage      !< Assign q to RK stage.
      procedure, pass(self) :: compute_stage     !< Compute RK stage.
      procedure, pass(self) :: compute_stage_ls  !< Compute RK stage, low storage scheme.
      procedure, pass(self) :: initialize        !< Initialize class.
      procedure, pass(self) :: initialize_stages !< Initialize RK stages.
      procedure, pass(self) :: update_q          !< Update RK q.
endtype rk_nvf_object
contains
   ! public methods
   subroutine assign_stage(self, s, q_gpu, phi_gpu)
   !< Assign q to RK stage.
   class(rk_nvf_object), intent(inout)                :: self        !< RK object.
   integer(I4P),         intent(in)                   :: s           !< Current stage number.
   real(R8P),            intent(in), device           :: q_gpu(1:     ,    &
                                                               1-self%rk%ngc:,&
                                                               1-self%rk%ngc:,&
                                                               1-self%rk%ngc:,&
                                                               1:)   !< Conservative variables.
   real(R8P),            intent(in), device, optional :: phi_gpu(1:,             &
                                                                 1-self%rk%ngc:, &
                                                                 1-self%rk%ngc:, &
                                                                 1-self%rk%ngc:, &
                                                                 1:) !< IB distance.

   associate(ni=>self%rk%ni, nj=>self%rk%nj, nk=>self%rk%nk, ngc=>self%rk%ngc, nv=>self%rk%nv, blocks_number=>self%rk%blocks_number)
      call rk_assign_stage_cuf(ni=ni,nj=nj,nk=nk,ngc=ngc,nv=nv,blocks_number=blocks_number,s=s,phi_gpu=phi_gpu,q_gpu=q_gpu, &
                               q_rk_gpu=self%q_rk_gpu)
      call self%mpih%check_cuda_error(error_code=CUDA_SILENT_ERR, msg='CUDA error in rk assign stage')
   endassociate
   endsubroutine assign_stage

   subroutine compute_stage(self, s, dt, phi_gpu)
   !< Compute RK stage.
   class(rk_nvf_object), intent(inout)                :: self        !< RK object.
   integer(I4P),         intent(in)                   :: s           !< Current stage number.
   real(R8P),            intent(in)                   :: dt          !< Current time step.
   real(R8P),            intent(in), device, optional :: phi_gpu(1:,             &
                                                                 1-self%rk%ngc:, &
                                                                 1-self%rk%ngc:, &
                                                                 1-self%rk%ngc:, &
                                                                 1:) !< IB distance.

   associate(ni=>self%rk%ni, nj=>self%rk%nj, nk=>self%rk%nk, ngc=>self%rk%ngc, nv=>self%rk%nv, blocks_number=>self%rk%blocks_number)
      call rk_compute_stage_cuf(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, s=s, dt=dt, alph=self%alph_gpu, &
                                phi_gpu=phi_gpu, q_rk_gpu=self%q_rk_gpu)
      call self%mpih%check_cuda_error(error_code=CUDA_SILENT_ERR, msg='CUDA error in rk compute stage')
   endassociate
   endsubroutine compute_stage

   subroutine compute_stage_ls(self, s, dt, phi_gpu, dq_gpu, q_gpu)
   !< Compute RK stage, low storage scheme.
   !< The first (only) stage is assumed to be the previous time step q solution.
   class(rk_nvf_object), intent(in)                      :: self        !< RK object.
   integer(I4P),         intent(in)                      :: s           !< Current RK stage.
   real(R8P),            intent(in)                      :: dt          !< Current time step.
   real(R8P),            intent(in),    device, optional :: phi_gpu(1:,             &
                                                                    1-self%rk%ngc:, &
                                                                    1-self%rk%ngc:, &
                                                                    1-self%rk%ngc:, &
                                                                    1:) !< IB distance.
   real(R8P),            intent(in),    device           :: dq_gpu(1:,       &
                                                                  1-self%rk%ngc:,&
                                                                  1-self%rk%ngc:,&
                                                                  1-self%rk%ngc:,&
                                                                  1:)   !< Conservative variables residuals.
   real(R8P),            intent(inout), device           :: q_gpu(1:,         &
                                                                  1-self%rk%ngc:,&
                                                                  1-self%rk%ngc:,&
                                                                  1-self%rk%ngc:,&
                                                                  1:)   !< Conservative variables stage.

   associate(ni=>self%rk%ni, nj=>self%rk%nj, nk=>self%rk%nk, ngc=>self%rk%ngc, nv=>self%rk%nv, blocks_number=>self%rk%blocks_number)
      call rk_compute_stage_ls_cuf(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, dt=dt, &
                                   ark=self%rk%ark(s), brk=self%rk%brk(s), crk=self%rk%crk(s),              &
                                   phi_gpu=phi_gpu, q_n_gpu=self%q_rk_gpu(:,:,:,:,:,1), dq_gpu=dq_gpu, q_rk_gpu=q_gpu)
      call self%mpih%check_cuda_error(error_code=CUDA_SILENT_ERR, msg='CUDA error in compute rk stage low storage scheme')
   endassociate
   endsubroutine compute_stage_ls

   subroutine initialize(self, rk, nb, ngc, ni, nj, nk, nv)
   !< Initialize class.
   class(rk_nvf_object), intent(inout)      :: self !< RK NVF object.
   type(rk_object),      intent(in), target :: rk   !< RK object.
   integer(I4P),         intent(in)         :: nb   !< Total blocks number for MPI.
   integer(I4P),         intent(in)         :: ngc  !< Number of ghost cells.
   integer(I4P),         intent(in)         :: ni   !< Number of cells in i direction.
   integer(I4P),         intent(in)         :: nj   !< Number of cells in j direction.
   integer(I4P),         intent(in)         :: nk   !< Number of cells in k direction.
   integer(I4P),         intent(in)         :: nv   !< Number of conservative variables.
   character(:), allocatable                :: msg_ !< Allocating message base.
   character(:), allocatable                :: msg  !< Allocating message.

   call self%mpih%initialize(do_mpi_init=.false.)
   call self%mpih%print_message('rk_nvf_object%initialize start')
   self%rk => rk
   msg_ = self%mpih%myrankstr//'rk_nvf_object%initialize'
   call assign_allocatable_gpu(lhs=self%alph_gpu, rhs=rk%alph, msg=msg_//' rk_alph ')
   call assign_allocatable_gpu(lhs=self%beta_gpu, rhs=rk%beta, msg=msg_//' rk_beta ')
   call assign_allocatable_gpu(lhs=self%gamm_gpu, rhs=rk%gamm, msg=msg_//' rk_gamm ')
   msg = msg_//' q_rk_gpu '
   select case(self%rk%scheme)
   case(RK_1, RK_2, RK_3) ! low storage, only stage 1 is necessary
   call alloc_var_gpu(var=self%q_rk_gpu,ulb=reshape([1,nb,1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,nv,1,1     ],[2,6]),msg=msg)
   case(RK_SSP_22, RK_SSP_33, RK_SSP_54)
   call alloc_var_gpu(var=self%q_rk_gpu,ulb=reshape([1,nb,1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,nv,1,rk%nrk],[2,6]),msg=msg)
   endselect
   self%q_rk_gpu = 0._R8P
   call self%mpih%print_message('rk_nvf_object%initialize finish')
   endsubroutine initialize

   subroutine initialize_stages(self, q_gpu)
   !< Initialize RK stages.
   class(rk_nvf_object), intent(inout)      :: self      !< RK object.
   real(R8P),            intent(in), device :: q_gpu(1:,             &
                                                     1-self%rk%ngc:, &
                                                     1-self%rk%ngc:, &
                                                     1-self%rk%ngc:, &
                                                     1:) !< Conservative variables.

   associate(ni=>self%rk%ni, nj=>self%rk%nj, nk=>self%rk%nk, ngc=>self%rk%ngc, nv=>self%rk%nv, blocks_number=>self%rk%blocks_number)
      call rk_initialize_stages_cuf(ni=ni,nj=nj,nk=nk,ngc=ngc,nv=nv,blocks_number=blocks_number,q_gpu=q_gpu,q_rk_gpu=self%q_rk_gpu)
      call self%mpih%check_cuda_error(error_code=CUDA_SILENT_ERR, msg='CUDA error in rk initialize stages')
   endassociate
   endsubroutine initialize_stages

   subroutine update_q(self, dt, phi_gpu, q_gpu)
   !< Update RK q.
   class(rk_nvf_object), intent(in)                      :: self      !< RK object.
   real(R8P),            intent(in)                      :: dt        !< Current time step.
   real(R8P),            intent(in),    device, optional :: phi_gpu(1:,             &
                                                                    1-self%rk%ngc:, &
                                                                    1-self%rk%ngc:, &
                                                                    1-self%rk%ngc:, &
                                                                    1:) !< IB distance.
   real(R8P),            intent(inout), device           :: q_gpu(1:     ,    &
                                                               1-self%rk%ngc:,&
                                                               1-self%rk%ngc:,&
                                                               1-self%rk%ngc:,&
                                                               1:)      !< Conservative variables.

   associate(ni=>self%rk%ni, nj=>self%rk%nj, nk=>self%rk%nk, ngc=>self%rk%ngc, nv=>self%rk%nv, blocks_number=>self%rk%blocks_number)
      call rk_update_q_cuf(ni=ni,nj=nj,nk=nk,ngc=ngc,nv=nv,blocks_number=blocks_number,nrk=self%rk%nrk,dt=dt,beta=self%beta_gpu, &
                           phi_gpu=phi_gpu,q_rk_gpu=self%q_rk_gpu,q_gpu=q_gpu)
      call self%mpih%check_cuda_error(error_code=CUDA_SILENT_ERR, msg='CUDA error in RK update q')
   endassociate
   endsubroutine update_q
endmodule adam_rk_nvf_object
