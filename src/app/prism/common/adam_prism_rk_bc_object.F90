!< PRISM, RK for BC integration.
module adam_prism_rk_bc_object
!< ADAM, RK-BC class definition.

! ADAM classes, libraries, parameters
use :: adam_rk_object
! ADAM singleton objects
use :: adam_global_field, only : field
use :: adam_global_grid,  only : grid
use :: adam_global_mpih,  only : mpih
! PRISM modules
use :: adam_prism_physics_object
! third party modules
use :: finer
use :: penf

implicit none
private
public :: prism_rk_bc_object

type :: prism_rk_bc_object
   !< RK class definition.
   character(:), pointer     :: scheme    !< RK scheme.
   integer(I4P)              :: nrk=3_I4P !< Runge-Kutta stages number.
   ! classic, Butcher schemes
   real(R8P), allocatable    :: ark(:)    !< Runge-Kutta low storage alpha coefficients.
   real(R8P), allocatable    :: brk(:)    !< Runge-Kutta low storage beta coefficients.
   real(R8P), allocatable    :: crk(:)    !< Runge-Kutta low storage beta coefficients.
   real(R8P), allocatable    :: alph(:,:) !< Runge-Kutta SSP alpha coefficients.
   real(R8P), allocatable    :: beta(:)   !< Runge-Kutta SSP beta coefficients.
   real(R8P), allocatable    :: gamm(:)   !< Runge-Kutta SSP gamma coefficients.
   ! symplectic (splitting) schemes
   real(R8P), allocatable :: ssa(:) !< Runge-Kutta sympletic-splitting part A coefficients.
   real(R8P), allocatable :: ssb(:) !< Runge-Kutta sympletic-splitting part B coefficients.
   ! RK data
   real(R8P), allocatable    :: q_bc_rk(:,:,:,:,:,:)
   real(R8P), allocatable    :: dq_bc_rk(:,:,:,:,:)
   ! grid/field data replica for easy handling
   integer(I4P),       pointer :: ngc=>null()  !< Number of ghost cells.
   integer(I4P),       pointer :: ni=>null()   !< Number of cells in i direction.
   integer(I4P),       pointer :: nj=>null()   !< Number of cells in j direction.
   integer(I4P),       pointer :: nk=>null()   !< Number of cells in k direction.
   integer(I4P),       pointer :: nv_c=>null() !< Number of variables.
   contains
      ! public methods
      procedure, pass(self) :: assign_stage      !< Assign q to RK stage.
      procedure, pass(self) :: compute_stage     !< Compute RK stage.
      !procedure, pass(self) :: compute_stage_ls  !< Compute RK stage, low storage scheme.
      procedure, pass(self) :: description       !< Return pretty-printed object description.
      procedure, pass(self) :: initialize        !< Initialize class.
      procedure, pass(self) :: initialize_stages !< Initialize RK stages.
      !procedure, pass(self) :: update_q          !< Update RK q.
endtype prism_rk_bc_object
contains
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(prism_rk_bc_object), intent(in)  :: self             !< RK object.
   character(len=:), allocatable :: desc             !< Description.
   character(len=1), parameter   :: NL=new_line('a') !< New line character.
   integer(I4P)                  :: s                !< Counter.

   desc =       mpih%myrankstr//'Runge-Kutta Boundary conditions scheme main data'//NL
   if (allocated(self%ark)) &
   desc = desc//mpih%myrankstr//'  ark:                             '//trim(str(self%ark                ))//NL
   if (allocated(self%brk)) &
   desc = desc//mpih%myrankstr//'  brk:                             '//trim(str(self%brk                ))//NL
   if (allocated(self%crk)) &
   desc = desc//mpih%myrankstr//'  crk:                             '//trim(str(self%crk                ))//NL
   if (allocated(self%alph)) then
   do s=1, self%nrk
   desc = desc//mpih%myrankstr//'  alph('//trim(str(s,.true.))//'): '//trim(str(self%alph(:,s)          ))//NL
   enddo
   endif
   if (allocated(self%beta)) &
   desc = desc//mpih%myrankstr//'  beta:                            '//trim(str(self%beta               ))//NL
   if (allocated(self%gamm)) &
   desc = desc//mpih%myrankstr//'  gamm:                            '//trim(str(self%gamm               ))//NL
   desc = desc//mpih%myrankstr//'  nrk:                             '//trim(str(self%nrk                ))
   endfunction description

   subroutine initialize(self, file_parameters, rk, physics)
   !< Initialize class.
   class(prism_rk_bc_object),   intent(inout)        :: self            !< RK object.
   type(file_ini),              intent(in), optional :: file_parameters !< Simulation parameters ini file handler.
   type(rk_object),             intent(in), target   :: rk              !< RK scheme
   type(prism_physics_object),  intent(in), target   :: physics         !< Physics object
   real(R8P)                                         :: w0, w1          !< Sympletic RK coefficients.

   call mpih%print_message('rk_object%initialize start')
   call associate_adam_data(rk=rk, physics=physics)
   select case(self%scheme)
   case(RK_1) ! 1 stage, 1st order, Euler
      self%nrk = 1
      allocate(self%ark(self%nrk), self%brk(self%nrk), self%crk(self%nrk))
      self%ark(1) = 1._R8P ; self%brk(1) = 0._R8P ; self%crk(1) = 1._R8P
   case(RK_2) ! 2 stages, low storage, 2nd order TVD
      self%nrk = 2
      allocate(self%ark(self%nrk), self%brk(self%nrk), self%crk(self%nrk))
      self%ark(1) = 1._R8P  ; self%brk(1) = 0._R8P  ; self%crk(1) = 1._R8P
      self%ark(2) = 0.5_R8P ; self%brk(2) = 0.5_R8P ; self%crk(2) = 0.5_R8P
   case(RK_3) ! 3 stages, low storage, 3rd order TVD
      self%nrk = 3
      allocate(self%ark(self%nrk), self%brk(self%nrk), self%crk(self%nrk))
      self%ark(1) = 1._R8P        ; self%brk(1) = 0._R8P        ; self%crk(1) = 1._R8P
      self%ark(2) = 0.75_R8P      ; self%brk(2) = 0.25_R8P      ; self%crk(2) = 0.25_R8P
      self%ark(3) = 1._R8P/3._R8P ; self%brk(3) = 2._R8P/3._R8P ; self%crk(3) = 2._R8P/3._R8P
   case(RK_SSP_22) ! 2 stages, 2nd order SSP
      self%nrk = 2
      allocate(self%alph(self%nrk,self%nrk), self%beta(self%nrk), self%gamm(self%nrk))
      self%alph = 0._R8P
      self%beta = 0._R8P
      self%gamm = 0._R8P

      self%beta(1) = 0.5_R8P
      self%beta(2) = 0.5_R8P

      self%alph(2,1) = 1._R8P

      self%gamm(2) = 1._R8P
   case(RK_SSP_33) ! 3 stages, 3rd order SSP
      self%nrk = 2
      allocate(self%alph(self%nrk,self%nrk), self%beta(self%nrk), self%gamm(self%nrk))
      self%alph = 0._R8P
      self%beta = 0._R8P
      self%gamm = 0._R8P

      self%beta(1) = 1._R8P/6._R8P
      self%beta(2) = 1._R8P/6._R8P
      self%beta(3) = 2._R8P/3._R8P

      self%alph(2,1) = 1._R8P
      self%alph(3,1) = 0.25_R8P ; self%alph(3,2) = 0.25_R8P

      self%gamm(2) = 1._R8P
      self%gamm(3) = 0.5_R8P
   case(RK_SSP_54) ! 5 stages, 4th order SSP
      self%nrk = 5
      allocate(self%alph(self%nrk,self%nrk), self%beta(self%nrk), self%gamm(self%nrk))
      self%alph = 0._R8P
      self%beta = 0._R8P
      self%gamm = 0._R8P

      self%beta(1) = 0.14681187618661_R8P
      self%beta(2) = 0.24848290924556_R8P
      self%beta(3) = 0.10425883036650_R8P
      self%beta(4) = 0.27443890091960_R8P
      self%beta(5) = 0.22600748319395_R8P

      self%alph(2,1) = 0.39175222700392_R8P
      self%alph(3,1) = 0.21766909633821_R8P ; self%alph(3,2) = 0.36841059262959_R8P
      self%alph(4,1) = 0.08269208670950_R8P ; self%alph(4,2) = 0.13995850206999_R8P ; self%alph(4,3) = 0.25189177424738_R8P
      self%alph(5,1) = 0.06796628370320_R8P ; self%alph(5,2) = 0.11503469844438_R8P ; self%alph(5,3) = 0.20703489864929_R8P
      self%alph(5,4) = 0.54497475021237_R8P

      self%gamm(2) = 0.39175222700392_R8P
      self%gamm(3) = 0.58607968896780_R8P
      self%gamm(4) = 0.47454236302687_R8P
      self%gamm(5) = 0.93501063100924_R8P
   case(RK_YOSHIDA)
      self%nrk = 4
      w0 = -1.702414383919315_R8P
      w1 =  1.351207191959658_R8P
      allocate(self%ssa(self%nrk), self%ssb(self%nrk-1))
      self%ssa = [w1/2.0_R8P,(w0+w1)/2.0_R8P,(w0+w1)/2.0_R8P,w1/2.0_R8P]
      self%ssb = [w1,w0,w1]
   case default
      !@TODO write error trap
   endselect

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv_c=>self%nv_c, nb=>field%nb, nrk=>self%nrk)
   select case(self%scheme)
   case(RK_1, RK_2, RK_3) ! low storage, only stage 1 is necessary

   case(RK_SSP_22, RK_SSP_33, RK_SSP_54)
      call allocate_variable(var=self%q_bc_rk,          &
                             ulb=reshape([1,nv_c,       &
                                          1-ngc,ni+ngc, &
                                          1-ngc,nj+ngc, &
                                          1-ngc,nk+ngc, &
                                          1,nb,         &
                                          1,nrk+1],[2,6]),  &
                             msg=mpih%myrankstr//'rk_object%initialize allocate q_bc_rk')
      call allocate_variable(var=self%dq_bc_rk,          &
                             ulb=reshape([1,nv_c,        &
                                          1-ngc,ni+ngc,  &
                                          1-ngc,nj+ngc,  &
                                          1-ngc,nk+ngc,  &
                                          1,nb           &
                                          ],[2,5]),      &
                             msg=mpih%myrankstr//'rk_object%initialize allocate dq_bc_rk')
   endselect
   endassociate
   print '(A)', self%description()
   call mpih%print_message('rk_object%initialize finish')
   contains
      subroutine associate_adam_data(rk, physics)
      !< Associate grid/physics/rk data pointers for easy handling.
      type(rk_object),            intent(in), target :: rk      !< The RK scheme.
      type(prism_physics_object), intent(in), target :: physics !< The physics.

      self%ni    => grid%ni
      self%nj    => grid%nj
      self%nk    => grid%nk
      self%ngc   => grid%ngc
      self%nv_c  => physics%nv_c
      self%scheme => rk%scheme
      endsubroutine associate_adam_data
   endsubroutine initialize

   subroutine initialize_stages(self, q)
   !< Initialize RK_bc stages.
   class(prism_rk_bc_object), intent(inout) :: self             !< RK object.
   real(R8P),                 intent(in)    :: q(1:,      &
                                             1-self%ngc:, &
                                             1-self%ngc:, &
                                             1-self%ngc:, &
                                             1:)                !< Conservative variables.
   integer(I4P)                             :: i, j, k, b, v, s !< Counter.

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>field%nv, nv_c=>self%nv_c, &
            blocks_number=>field%blocks_number)
   !$omp parallel do collapse(6) default(firstprivate) shared(q,self)
   do s=lbound(self%q_bc_rk,dim=6),ubound(self%q_bc_rk,dim=6)-1
      do b=1, blocks_number
         do k=1-ngc, nk+ngc
            do j=1-ngc, nj+ngc
               do i=1-ngc, ni+ngc
                  do v=1, nv_c
                     self%q_bc_rk(v,i,j,k,b,s) = q(v,i,j,k,b)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   !$omp parallel do collapse(6) default(firstprivate) shared(self)
   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               do v=1, nv_c
                  self%dq_bc_rk(v,i,j,k,b) = 0._R8P
               enddo
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endassociate
   endsubroutine initialize_stages

   subroutine compute_stage(self, s, dt, phi)
   !< Compute RK stage.
   class(prism_rk_bc_object), intent(inout) :: self               !< RK object.
   integer(I4P),      intent(in)            :: s                  !< Current stage number.
   real(R8P),         intent(in)            :: dt                 !< Current time step.
   real(R8P),         intent(in), optional  :: phi(1:,          &
                                                   1-self%ngc:, &
                                                   1-self%ngc:, &
                                                   1-self%ngc:, &
                                                   1:)            !< IB distance.
   integer(I4P)                             :: all_solids         !< Last phi index, all solids summary.
   integer(I4P)                             :: i, j, k, b, v, ss  !< Counter.

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>field%nv, nv_c=>self%nv_c, &
             blocks_number=>field%blocks_number)
   if (present(phi)) then
      all_solids = ubound(phi, dim=1)
      !$omp parallel do collapse(6) default(firstprivate) shared(phi,self)
      do ss=1, s-1
         do b=1, blocks_number
            do k=1-ngc, nk+ngc
               do j=1-ngc, nj+ngc
                  do i=1-ngc, ni+ngc
                     do v=1, nv_c
                        if (phi(all_solids,i,j,k,b) < 0._R8P) then
                           self%q_bc_rk(v,i,j,k,b,s) = self%q_bc_rk(v,i,j,k,b,s) + &
                                                         dt * self%alph(s,ss) * self%q_bc_rk(v,i,j,k,b,ss)
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
   else
      !$omp parallel do collapse(6) default(firstprivate) shared(self)
      do ss=1, s-1
         do b=1, blocks_number
            do k=1-ngc, nk+ngc
               do j=1-ngc, nj+ngc
                  do i=1-ngc, ni+ngc
                     do v=1, nv_c
                        self%q_bc_rk(v,i,j,k,b,s) = self%q_bc_rk(v,i,j,k,b,s) + &
                                                      dt * self%alph(s,ss) * self%q_bc_rk(v,i,j,k,b,ss)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
   endif
   endassociate
   endsubroutine compute_stage

   subroutine assign_stage(self, s, phi)
   !< Assign q to RK stage.
   class(prism_rk_bc_object), intent(inout)        :: self          !< RK object.
   integer(I4P),     intent(in)           :: s             !< Current stage number.
   real(R8P),        intent(in), optional :: phi(1:,          &
                                                 1-self%ngc:, &
                                                 1-self%ngc:, &
                                                 1-self%ngc:, &
                                                 1:)       !< IB distance.
   integer(I4P)                           :: all_solids    !< Last phi index, all solids summary.
   integer(I4P)                           :: i, j, k, b, v !< Counter.

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv_c=>self%nv_c, blocks_number=>field%blocks_number)
   if (present(phi)) then
      all_solids = ubound(phi, dim=1)
      !$omp parallel do collapse(5) default(firstprivate) shared(phi,q,self)
      do b=1, blocks_number
         do k=1-ngc, nk+ngc
            do j=1-ngc, nj+ngc
               do i=1-ngc, ni+ngc
                  do v=1, nv_c
                     if (phi(all_solids,i,j,k,b) < 0._R8P) then
                        self%q_bc_rk(v,i,j,k,b,s) = self%dq_bc_rk(v,i,j,k,b)
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
   else
      !$omp parallel do collapse(5) default(firstprivate) shared(q,self)
      do b=1, blocks_number
         do k=1-ngc, nk+ngc
            do j=1-ngc, nj+ngc
               do i=1-ngc, ni+ngc
                  do v=1, nv_c
                     self%q_bc_rk(v,i,j,k,b,s) = self%dq_bc_rk(v,i,j,k,b)
                  enddo
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
   endif
   endassociate
   endsubroutine assign_stage

endmodule adam_prism_rk_bc_object
