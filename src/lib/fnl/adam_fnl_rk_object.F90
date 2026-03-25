!< ADAM, RK class FNL (FNL backend of [[rk_object]]).
module adam_fnl_rk_object
!< ADAM, RK class FNL (FNL backend of [[rk_object]]).

! ADAM classes, libraries, parameters
use :: adam_rk_object
! ADAM singleton objects
use :: adam_global_field, only : field
use :: adam_global_grid,  only : grid
! ADAM FNL classes, libraries, parameters
use :: adam_fnl_rk_kernels
! ADAM FNL singleton objects
use :: adam_global_mpih_fnl, only : mpih_fnl
! third party modules
use :: fundal
use :: penf

implicit none
save
private
public :: rk_fnl_object

type :: rk_fnl_object
   !< RK FNL class definition.
   ! ADAM library objects
   type(rk_object), pointer :: rk=>null() !< RK common handler.
   ! device data
   real(R8P), pointer :: alph_gpu(:,:)=>null()         !< RK alpha coefficients.
   real(R8P), pointer :: beta_gpu(:)=>null()           !< RK beta coefficients.
   real(R8P), pointer :: gamm_gpu(:)=>null()           !< RK gamma coefficients.
   real(R8P), pointer :: q_rk_gpu(:,:,:,:,:,:)=>null() !< Field cell centered variables, RK stages.
   contains
      ! public methods
      procedure, pass(self) :: assign_stage      !< Assign q to RK stage.
      procedure, pass(self) :: compute_stage     !< Compute RK stage.
      procedure, pass(self) :: compute_stage_ls  !< Compute RK stage, low storage scheme.
      procedure, pass(self) :: initialize        !< Initialize class.
      procedure, pass(self) :: initialize_stages !< Initialize RK stages.
      procedure, pass(self) :: update_q          !< Update RK q.
endtype rk_fnl_object
contains
   ! public methods
   subroutine assign_stage(self, s, q_gpu, phi_gpu)
   !< Assign q to RK stage.
   class(rk_fnl_object), intent(inout)        :: self        !< RK object.
   integer(I4P),         intent(in)           :: s           !< Current stage number.
   real(R8P),            intent(in)           :: q_gpu(1:     ,        &
                                                       1-self%rk%ngc:, &
                                                       1-self%rk%ngc:, &
                                                       1-self%rk%ngc:, &
                                                       1:)   !< Conservative variables.
   real(R8P),            intent(in), optional :: phi_gpu(1:,             &
                                                         1-self%rk%ngc:, &
                                                         1-self%rk%ngc:, &
                                                         1-self%rk%ngc:, &
                                                         1:) !< IB distance.

   associate(ni=>grid%ni, nj=>grid%nj, nk=>grid%nk, ngc=>grid%ngc, nv=>field%nv, blocks_number=>field%blocks_number)
      call rk_assign_stage_dev(ni=ni,nj=nj,nk=nk,ngc=ngc,nv=nv,blocks_number=blocks_number,s=s,phi_gpu=phi_gpu,q_gpu=q_gpu, &
                               q_rk_gpu=self%q_rk_gpu)
   endassociate
   endsubroutine assign_stage

   subroutine compute_stage(self, s, dt, phi_gpu)
   !< Compute RK stage.
   class(rk_fnl_object), intent(inout)        :: self        !< RK object.
   integer(I4P),         intent(in)           :: s           !< Current stage number.
   real(R8P),            intent(in)           :: dt          !< Current time step.
   real(R8P),            intent(in), optional :: phi_gpu(1:,             &
                                                         1-self%rk%ngc:, &
                                                         1-self%rk%ngc:, &
                                                         1-self%rk%ngc:, &
                                                         1:) !< IB distance.

   associate(ni=>grid%ni, nj=>grid%nj, nk=>grid%nk, ngc=>grid%ngc, nv=>field%nv, blocks_number=>field%blocks_number)
      call rk_compute_stage_dev(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, s=s, dt=dt, alph=self%alph_gpu,&
                                phi_gpu=phi_gpu, q_rk_gpu=self%q_rk_gpu)
   endassociate
   endsubroutine compute_stage

   subroutine compute_stage_ls(self, s, dt, phi_gpu, dq_gpu, q_gpu)
   !< Compute RK stage, low storage scheme.
   !< The first (only) stage is assumed to be the previous time step q solution.
   class(rk_fnl_object), intent(in)           :: self        !< RK object.
   integer(I4P),         intent(in)           :: s           !< Current RK stage.
   real(R8P),            intent(in)           :: dt          !< Current time step.
   real(R8P),            intent(in), optional :: phi_gpu(1:,             &
                                                         1-self%rk%ngc:, &
                                                         1-self%rk%ngc:, &
                                                         1-self%rk%ngc:, &
                                                         1:) !< IB distance.
   real(R8P),            intent(in)           :: dq_gpu(1:,            &
                                                       1-self%rk%ngc:, &
                                                       1-self%rk%ngc:, &
                                                       1-self%rk%ngc:, &
                                                       1:)   !< Conservative variables residuals.
   real(R8P),            intent(inout)        :: q_gpu(1:,             &
                                                       1-self%rk%ngc:, &
                                                       1-self%rk%ngc:, &
                                                       1-self%rk%ngc:, &
                                                       1:)   !< Conservative variables stage.

   associate(ni=>grid%ni, nj=>grid%nj, nk=>grid%nk, ngc=>grid%ngc, nv=>field%nv, blocks_number=>field%blocks_number)
      call rk_compute_stage_ls_dev(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, dt=dt, &
                                   ark=self%rk%ark(s), brk=self%rk%brk(s), crk=self%rk%crk(s),              &
                                   phi_gpu=phi_gpu, q_n_gpu=self%q_rk_gpu(:,:,:,:,:,1), dq_gpu=dq_gpu, q_rk_gpu=q_gpu)
   endassociate
   endsubroutine compute_stage_ls

   subroutine initialize(self, rk, nb, ngc, ni, nj, nk, nv)
   !< Initialize class.
   !< Requires `mpih_fnl` (adam_global_mpih_fnl) to be initialized before calling.
   class(rk_fnl_object), intent(inout)      :: self !< RK FNL object.
   type(rk_object),      intent(in), target :: rk   !< RK object.
   integer(I4P),         intent(in)         :: nb   !< Total blocks number for MPI.
   integer(I4P),         intent(in)         :: ngc  !< Number of ghost cells.
   integer(I4P),         intent(in)         :: ni   !< Number of cells in i direction.
   integer(I4P),         intent(in)         :: nj   !< Number of cells in j direction.
   integer(I4P),         intent(in)         :: nk   !< Number of cells in k direction.
   integer(I4P),         intent(in)         :: nv   !< Number of conservative variables.
   integer(I4P)                             :: ierr !< Error status.
   integer(I4P)                             :: nrk  !< RK stages.

   call mpih_fnl%print_message('rk_fnl_object%initialize start')
   self%rk => rk
   call dev_assign_to_device(src=rk%alph, dst=self%alph_gpu)
   call dev_assign_to_device(src=rk%beta, dst=self%beta_gpu)
   call dev_assign_to_device(src=rk%gamm, dst=self%gamm_gpu)
   select case(self%rk%scheme)
   case(RK_1, RK_2, RK_3) ! low storage, only stage 1 is necessary
      nrk = 1
   case(RK_SSP_22, RK_SSP_33, RK_SSP_54)
      nrk = rk%nrk
   endselect
   call dev_alloc(fptr_dev=self%q_rk_gpu,                   &
                  ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,nv,nrk], &
                  lbounds=[1,1-ngc,1-ngc,1-ngc,1,1],        &
                  init_value=0._R8P,                        &
                  ierr=ierr)
   call mpih_fnl%print_message('rk_fnl_object%initialize finish')
   endsubroutine initialize

   subroutine initialize_stages(self, q_gpu)
   !< Initialize RK stages.
   class(rk_fnl_object), intent(inout) :: self      !< RK object.
   real(R8P),            intent(in)    :: q_gpu(1:,             &
                                                1-self%rk%ngc:, &
                                                1-self%rk%ngc:, &
                                                1-self%rk%ngc:, &
                                                1:) !< Conservative variables.

   associate(ni=>grid%ni, nj=>grid%nj, nk=>grid%nk, ngc=>grid%ngc, nv=>field%nv, blocks_number=>field%blocks_number)
      call rk_initialize_stages_dev(ni=ni,nj=nj,nk=nk,ngc=ngc,nv=nv,blocks_number=blocks_number,q_gpu=q_gpu,q_rk_gpu=self%q_rk_gpu)
   endassociate
   endsubroutine initialize_stages

   subroutine update_q(self, dt, phi_gpu, q_gpu)
   !< Update RK q.
   class(rk_fnl_object), intent(in)           :: self      !< RK object.
   real(R8P),            intent(in)           :: dt        !< Current time step.
   real(R8P),            intent(in), optional :: phi_gpu(1:,             &
                                                         1-self%rk%ngc:, &
                                                         1-self%rk%ngc:, &
                                                         1-self%rk%ngc:, &
                                                         1:) !< IB distance.
   real(R8P),            intent(inout)        :: q_gpu(1:,          &
                                                    1-self%rk%ngc:, &
                                                    1-self%rk%ngc:, &
                                                    1-self%rk%ngc:, &
                                                    1:)      !< Conservative variables.

   associate(ni=>grid%ni, nj=>grid%nj, nk=>grid%nk, ngc=>grid%ngc, nv=>field%nv, blocks_number=>field%blocks_number)
      call rk_update_q_dev(ni=ni,nj=nj,nk=nk,ngc=ngc,nv=nv,blocks_number=blocks_number,nrk=self%rk%nrk,dt=dt,beta=self%beta_gpu, &
                           phi_gpu=phi_gpu,q_rk_gpu=self%q_rk_gpu,q_gpu=q_gpu)
   endassociate
   endsubroutine update_q
endmodule adam_fnl_rk_object
