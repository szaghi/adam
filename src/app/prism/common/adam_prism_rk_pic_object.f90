!< PRISM, RK for BC integration.
module adam_prism_rk_pic_object
!< ADAM, RK-BC class definition.

!use adam_field_object
!use adam_grid_object
use adam_mpih_object
use adam_rk_object
use adam_prism_pic_object
use finer
use penf

implicit none
save
private
public :: prism_rk_pic_object

type :: prism_rk_pic_object
   !< RK class definition.
   type(mpih_object)         :: mpih      !< MPI handler.
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
   real(R8P), allocatable    :: q_pic_rk(:,:,:) !< RK stages for pic variables [particle, variable, stage].
   ! adam data replica for easy handling
   !type(field_object), pointer :: field=>null()           !< The field.
   !type(grid_object),  pointer :: grid=>null()            !< The grid.
   !integer(I4P),       pointer :: ngc=>null()             !< Number of ghost cells.
   !integer(I4P),       pointer :: ni=>null()              !< Number of cells in i direction.
   !integer(I4P),       pointer :: nj=>null()              !< Number of cells in j direction.
   !integer(I4P),       pointer :: nk=>null()              !< Number of cells in k direction.
   !integer(I4P),       pointer :: nb=>null()              !< Total blocks number for MPI.
   !integer(I4P),       pointer :: blocks_number=>null()   !< Actual blocks number.
   !integer(I4P),       pointer :: nv=>null()              !< Number of variables.
   integer(I4P),       pointer :: particle_number=>null() !< Number of particles.
   contains
      ! public methods
      !procedure, pass(self) :: assign_stage      !< Assign q to RK stage.
      procedure, pass(self) :: compute_stage     !< Compute RK stage.
      !procedure, pass(self) :: compute_stage_ls  !< Compute RK stage, low storage scheme.
      procedure, pass(self) :: description       !< Return pretty-printed object description.
      procedure, pass(self) :: initialize        !< Initialize class.
      procedure, pass(self) :: initialize_stages !< Initialize RK stages.
      !procedure, pass(self) :: update_q          !< Update RK q.
endtype prism_rk_pic_object
contains
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(prism_rk_pic_object), intent(in)  :: self             !< RK object.
   character(len=:), allocatable :: desc             !< Description.
   character(len=1), parameter   :: NL=new_line('a') !< New line character.
   integer(I4P)                  :: s                !< Counter.

   desc =       self%mpih%myrankstr//'Pic Runge-Kutta scheme main data'//NL
   if (allocated(self%ark)) &
   desc = desc//self%mpih%myrankstr//'  ark:                             '//trim(str(self%ark                ))//NL
   if (allocated(self%brk)) &
   desc = desc//self%mpih%myrankstr//'  brk:                             '//trim(str(self%brk                ))//NL
   if (allocated(self%crk)) &
   desc = desc//self%mpih%myrankstr//'  crk:                             '//trim(str(self%crk                ))//NL
   if (allocated(self%alph)) then
   do s=1, self%nrk
   desc = desc//self%mpih%myrankstr//'  alph('//trim(str(s,.true.))//'): '//trim(str(self%alph(:,s)          ))//NL
   enddo
   endif
   if (allocated(self%beta)) &
   desc = desc//self%mpih%myrankstr//'  beta:                            '//trim(str(self%beta               ))//NL
   if (allocated(self%gamm)) &
   desc = desc//self%mpih%myrankstr//'  gamm:                            '//trim(str(self%gamm               ))//NL
   desc = desc//self%mpih%myrankstr//'  nrk:                             '//trim(str(self%nrk                ))
   endfunction description

   subroutine initialize(self, file_parameters, rk, pic)
   !< Initialize class.
   class(prism_rk_pic_object),  intent(inout)        :: self            !< RK object.
   type(file_ini),              intent(in), optional :: file_parameters !< Simulation parameters ini file handler.
   !type(grid_object),           intent(in), target   :: grid            !< The grid.
   !type(field_object),          intent(in), target   :: field           !< The field.
   type(rk_object),             intent(in), target   :: rk              !< RK scheme
   type(prism_pic_object),  	  intent(in), target   :: pic         		!< Physics object
   real(R8P)                                         :: w0, w1          !< Sympletic RK coefficients.

   call self%mpih%initialize(do_mpi_init=.false.)
   call self%mpih%print_message('rk_object%initialize start')
   call associate_adam_data(rk=rk, pic=pic)
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

   !associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, particle_number=>self%particle_number, &
	!			 nb=>self%nb, nrk=>self%nrk)
	associate(particle_number=>self%particle_number, nrk=>self%nrk)
   select case(self%scheme)
   case(RK_1, RK_2, RK_3) ! low storage, only stage 1 is necessary

   case(RK_SSP_22, RK_SSP_33, RK_SSP_54)
      call allocate_variable(var=self%q_pic_rk,         		 &
                             ulb=reshape([1,particle_number, &
                                          1,8, 					 &
                                          1,nrk+1],[2,3]),   &
                             msg=self%mpih%myrankstr//'rk_pic_object%initialize allocate q_pic_rk')
   endselect
   endassociate
   print '(A)', self%description()
   call self%mpih%print_message('rk_object%initialize finish')
   contains
      subroutine associate_adam_data(rk, pic)
      !< Associate objects data to equation for easy handling.
      !type(grid_object),        intent(in), target :: grid   !< The grid.
      !type(field_object),       intent(in), target :: field  !< The field.
      type(rk_object),          intent(in), target :: rk     !< The RK scheme.
      type(prism_pic_object),   intent(in), target :: pic    !< Particle in Cell object.

      !self%grid          	=> grid
      !self%field         	=> field
      !self%blocks_number 	=> field%blocks_number
      !self%ni            	=> field%grid%ni
      !self%nj            	=> field%grid%nj
      !self%nk            	=> field%grid%nk
      !self%ngc           	=> field%grid%ngc
      !self%nb            	=> field%nb
      !self%nv            	=> field%nv
      self%particle_number => pic%particle_number
      self%scheme        	=> rk%scheme
      endsubroutine associate_adam_data
   endsubroutine initialize


   subroutine initialize_stages(self, q_pic)
   !< Initialize RK stages.
   class(prism_rk_pic_object), intent(inout) :: self         !< RK object.
   real(R8P),        intent(in)    				:: q_pic(1:,1:) !< Conservative variables.
   integer(I4P)                    				:: p, v, s 		 !< Counter.

   associate(particle_number => self%particle_number)
   !$omp parallel do collapse(6) default(firstprivate) shared(q_pic,self)
   do s=lbound(self%q_pic_rk,dim=3),ubound(self%q_pic_rk,dim=3)
      do p=1, particle_number
         do v=1, 8
            self%q_pic_rk(p,v,s) = q_pic(p,v)
         enddo
      enddo
   enddo
   !$omp end parallel do
   endassociate
   endsubroutine initialize_stages

	subroutine compute_stage(self, s, dt)!, phi)
   !< Compute RK stage.
   class(prism_rk_pic_object), intent(inout)    :: self     !< RK object.
   integer(I4P),     intent(in)           		:: s        !< Current stage number.
   real(R8P),        intent(in)           		:: dt       !< Current time step.
   !real(R8P),        intent(in), optional :: phi(1:,          &
   !                                              1-self%ngc:, &
   !                                              1-self%ngc:, &
   !                                              1-self%ngc:, &
   !                                              1:) !< IB distance.
   !integer(I4P)                           :: all_solids        !< Last phi index, all solids summary.
   integer(I4P)                           :: p, v, ss !< Counter.

   !associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>self%nv, blocks_number=>self%blocks_number)
	associate(particle_number => self%particle_number)
   !if (present(phi)) then
      !all_solids = ubound(phi, dim=1)
      !!$omp parallel do collapse(6) default(firstprivate) shared(phi,self)
      !do ss=1, s-1
      !   do b=1, blocks_number
      !      do k=1, nk
      !         do j=1, nj
      !            do i=1, ni
      !               do v=1, nv
      !                  if (phi(all_solids,i,j,k,b) < 0._R8P) then
      !                     self%q_rk(v,i,j,k,b,s) = self%q_rk(v,i,j,k,b,s) + dt * self%alph(s,ss) * self%q_rk(v,i,j,k,b,ss)
      !                  endif
      !               enddo
      !            enddo
      !         enddo
      !      enddo
      !   enddo
      !enddo
      !!$omp end parallel do
   !else
   !$omp parallel do collapse(3) default(firstprivate) shared(self)
   do ss=1, s-1
      do p=1, particle_number
         do v=1, 6
            self%q_pic_rk(p,v,s) = self%q_pic_rk(p,v,s) + dt * self%alph(s, ss) * self%q_pic_rk(p,v,ss)
         enddo
      enddo
   enddo
      !$omp end parallel do
   !endif
   endassociate
   endsubroutine compute_stage

	!subroutine assign_stage(self, s, q)!, phi)
   !!< Assign q to RK stage.
   !class(prism_rk_pic_object), intent(inout)        :: self          !< RK object.
   !integer(I4P),     intent(in)           :: s             !< Current stage number.
   !real(R8P),        intent(in)           :: q(1:     ,     &
   !                                            1-self%ngc:, &
   !                                            1-self%ngc:, &
   !                                            1-self%ngc:, &
   !                                            1:)         !< Conservative variables.
   !real(R8P),        intent(in), optional :: phi(1:,          &
   !                                              1-self%ngc:, &
   !                                              1-self%ngc:, &
   !                                              1-self%ngc:, &
   !                                              1:)       !< IB distance.
   !integer(I4P)                           :: all_solids    !< Last phi index, all solids summary.
   !integer(I4P)                           :: i, j, k, b, v !< Counter.
!
   !associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>self%nv, blocks_number=>self%blocks_number)
   !if (present(phi)) then
   !   all_solids = ubound(phi, dim=1)
   !   !$omp parallel do collapse(5) default(firstprivate) shared(phi,q,self)
   !   do b=1, blocks_number
   !      do k=1, nk
   !         do j=1, nj
   !            do i=1, ni
   !               do v=1, nv
   !                  if (phi(all_solids,i,j,k,b) < 0._R8P) then
   !                     self%q_rk(v,i,j,k,b,s) = q(v,i,j,k,b)
   !                  endif
   !               enddo
   !            enddo
   !         enddo
   !      enddo
   !   enddo
   !   !$omp end parallel do
   !else
   !   !$omp parallel do collapse(5) default(firstprivate) shared(q,self)
   !   do b=1, blocks_number
   !      do k=1, nk
   !         do j=1, nj
   !            do i=1, ni
   !               do v=1, nv
   !                  self%q_rk(v,i,j,k,b,s) = q(v,i,j,k,b)
   !               enddo
   !            enddo
   !         enddo
   !      enddo
   !   enddo
   !   !$omp end parallel do
   !endif
   !endassociate
   !endsubroutine assign_stage


endmodule adam_prism_rk_pic_object