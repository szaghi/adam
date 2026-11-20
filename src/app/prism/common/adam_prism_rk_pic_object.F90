!< PRISM, the particles' Runge-Kutta integrator (PIC), on the field stepper's SELAH tableau.
module adam_prism_rk_pic_object
!< PRISM, the particles' Runge-Kutta integrator (PIC), on the field stepper's SELAH tableau.

! ADAM classes, libraries, parameters
use :: adam_stepper_object, only : stepper_object, RK_SSP_22, RK_SSP_33, RK_SSP_54
! ADAM singleton objects
use :: adam_mpih_global,  only : mpih
! PRISM modules
use :: adam_prism_parameters
use :: adam_prism_pic_object
use :: adam_prism_physics_object, only : prism_physics_object, is_adim_model
! third party modules
use :: finer
use :: penf

implicit none
save
private
public :: prism_rk_pic_object

type :: prism_rk_pic_object
   !< RK class definition.
   character(:), allocatable :: scheme    !< RK scheme, the field stepper's.
   integer(I4P)              :: nrk=0_I4P !< Runge-Kutta stages number (0: no SSP scheme, the app refuses it).
   real(R8P)                 :: inv_eps_scale = 1._R8P/EPS0 !< Convert D to E in the particle force.
   ! the SSP Butcher tableau, the field stepper's (SELAH's)
   real(R8P), allocatable    :: alph(:,:) !< Butcher matrix A.
   real(R8P), allocatable    :: beta(:)   !< Butcher weights b.
   real(R8P), allocatable    :: gamm(:)   !< Butcher abscissas c.
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
      procedure, pass(self) :: assign_stage      !< Assign q to RK stage.
      procedure, pass(self) :: compute_stage     !< Compute RK stage.
      !procedure, pass(self) :: compute_stage_ls  !< Compute RK stage, low storage scheme.
      procedure, pass(self) :: description       !< Return pretty-printed object description.
      procedure, pass(self) :: initialize        !< Initialize class.
      procedure, pass(self) :: initialize_stages !< Initialize RK stages.
      procedure, pass(self) :: update_q_pic      !< Update RK q.
endtype prism_rk_pic_object
contains
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(prism_rk_pic_object), intent(in)  :: self   !< RK object.
   character(len=:), allocatable :: desc             !< Description.
   character(len=1), parameter   :: NL=new_line('a') !< New line character.
   integer(I4P)                  :: s                !< Counter.

   desc =       mpih%myrankstr//'Pic Runge-Kutta scheme main data'//NL
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

   subroutine initialize(self, ts, pic, physics)
   !< Initialize class: the particles advance in step with the field's stepper, with its SELAH Butcher tableau (A, b, c:
   !< bitwise the former `rk_object` tableau, SELAH's selah_test_adam_parity; #132 phase 6). A scheme other than SSP-22,
   !< -33, -54 leaves `nrk = 0` and nothing allocated: the backends refuse it for PIC.
   class(prism_rk_pic_object),  intent(inout)      :: self    !< RK object.
   type(stepper_object),        intent(in)         :: ts      !< The field's stepper (initialized).
   type(prism_pic_object),      intent(in), target :: pic     !< Particle in Cell object.
   type(prism_physics_object),  intent(in)         :: physics !< Physics.

   call mpih%print_message('rk_pic_object%initialize start')
   self%particle_number => pic%particle_number
   self%scheme = trim(ts%scheme)
   self%inv_eps_scale = 1._R8P/EPS0
   if (is_adim_model(physics%physical_model)) self%inv_eps_scale = 1._R8P
   select case(self%scheme)
   case(RK_SSP_22, RK_SSP_33, RK_SSP_54)
      self%nrk  = ts%stages()
      self%alph = ts%st%bt%a
      self%beta = ts%st%bt%b
      self%gamm = ts%st%bt%c
      call allocate_variable(var=self%q_pic_rk,                         &
                             ulb=reshape([1,8,                          &
                                          1,pic%particle_capacity,      &
                                          1,self%nrk+1],[2,3]),         &
                             msg=mpih%myrankstr//'rk_pic_object%initialize allocate q_pic_rk')
   case default
      self%nrk = 0_I4P
   endselect
   print '(A)', self%description()
   call mpih%print_message('rk_pic_object%initialize finish')
   endsubroutine initialize

   subroutine initialize_stages(self, q_pic)
   !< Initialize RK stages.
   class(prism_rk_pic_object), intent(inout) :: self         !< RK object.
   real(R8P),        intent(in)    				:: q_pic(1:,1:) !< Conservative variables.
   integer(I4P)                    				:: p, v, s 		 !< Counter.

   associate(particle_number => self%particle_number)
   !$omp parallel do collapse(3) default(firstprivate) shared(q_pic,self)
   do s=lbound(self%q_pic_rk,dim=3),ubound(self%q_pic_rk,dim=3)
      do p=1, particle_number
         do v=1, 8
            self%q_pic_rk(v,p,s) = q_pic(v,p)
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
   integer(I4P)                                 :: p, v, ss !< Counter.

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
   ! The ss reduction into q_pic_rk(:,p,s) must stay the inner sequential loop;
   ! only the particle/variable loops may be parallelised. Folding `ss` into the
   ! collapse races every thread on the same q_pic_rk(v,p,s). Serially
   ! bit-identical (per-point ss-order unchanged).
   !$omp parallel do collapse(2) default(firstprivate) shared(self)
   do p=1, particle_number
      do v=1, 6
         do ss=1, s-1
            self%q_pic_rk(v,p,s) = self%q_pic_rk(v,p,s) + dt * self%alph(s, ss) * self%q_pic_rk(v,p,ss)
         enddo
      enddo
   enddo
      !$omp end parallel do
   !endif
   endassociate
   endsubroutine compute_stage

	subroutine assign_stage(self, s, pic_fields)!, phi)
   !< Assign q to RK stage.
   class(prism_rk_pic_object), intent(inout) :: self                   !< RK object.
   integer(I4P),     intent(in)              :: s                      !< Current stage number.
   real(R8P),        intent(in)              :: pic_fields(1:,1:)      !< Pic fields.
   real(R8P)                                 :: v_p(3), B_p(3), E_p(3) !< Velocity and fields at particle position.
   real(R8P)                                 :: q_p, m_p               !< Particle charge and mass.
   real(R8P)                                 :: F_l(3), F_p(3)         !< Lorentz force and total force.
   !real(R8P),        intent(in), optional :: phi(1:,          &
   !                                              1-self%ngc:, &
   !                                              1-self%ngc:, &
   !                                              1-self%ngc:, &
   !                                              1:)       !< IB distance.
   !integer(I4P)                           :: all_solids    !< Last phi index, all solids summary.
   integer(I4P)                              :: p, v !< Counter.
!
   !associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>self%nv, blocks_number=>self%blocks_number)
   associate(particle_number => self%particle_number, q_pic_rk=>self%q_pic_rk)
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
      !$omp parallel do collapse(1) default(firstprivate) shared(pic_fields,self)
      do p=1, particle_number
         v_p = [q_pic_rk(4,p,s), q_pic_rk(5,p,s), q_pic_rk(6,p,s)]
         q_p = q_pic_rk(7,p,s)
         m_p = q_pic_rk(8,p,s)
         E_p = pic_fields(1:3,p)*self%inv_eps_scale
         B_p = pic_fields(4:6,p)
         F_l = crossproduct(a=v_p, b=B_p)
         F_p(1) = q_p*(E_p(1) + F_l(1))
         F_p(2) = q_p*(E_p(2) + F_l(2))
         F_p(3) = q_p*(E_p(3) + F_l(3))

         self%q_pic_rk(1,p,s) = v_p(1)
         self%q_pic_rk(2,p,s) = v_p(2)
         self%q_pic_rk(3,p,s) = v_p(3)
         self%q_pic_rk(4,p,s) = F_p(1)/m_p
         self%q_pic_rk(5,p,s) = F_p(2)/m_p
         self%q_pic_rk(6,p,s) = F_p(3)/m_p
         self%q_pic_rk(7,p,s) = 0._R8P
         self%q_pic_rk(8,p,s) = 0._R8P
      enddo
      !$omp end parallel do
   !endif
   endassociate
   endsubroutine assign_stage

   subroutine update_q_pic(self, dt, q_pic)!phi, q)
   !< Update RK q.
   class(prism_rk_pic_object), intent(in)           :: self             !< RK object.
   real(R8P),                  intent(in)           :: dt               !< Current time step.
   !real(R8P),        intent(in), optional :: phi(1:,          &
   !                                              1-self%ngc:, &
   !                                              1-self%ngc:, &
   !                                              1-self%ngc:, &
   !                                              1:)          !< IB distance.
   real(R8P),                  intent(inout)    	 :: q_pic(1:,1:)     !< Conservative variables.
   !integer(I4P)                           :: all_solids       !< Last phi index, all solids summary.
   integer(I4P)                                     :: p, v, s          !< Counter.

   !associate(nrk=>self%nrk, ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>self%nv, blocks_number=>self%blocks_number)
   associate(particle_number=>self%particle_number, nrk=>self%nrk)
   !if (present(phi)) then
   !   all_solids = ubound(phi, dim=1)
   !   !$omp parallel do collapse(6) default(firstprivate) shared(phi,q,self)
   !   do s=1, nrk
   !      do b=1, blocks_number
   !         do k=1, nk
   !            do j=1, nj
   !               do i=1, ni
   !                  do v=1, nv
   !                     if (phi(all_solids,i,j,k,b) < 0._R8P) then
   !                        q(v,i,j,k,b) = q(v,i,j,k,b) + dt * self%beta(s) * self%q_rk(v,i,j,k,b,s)
   !                     endif
   !                  enddo
   !               enddo
   !            enddo
   !         enddo
   !      enddo
   !   enddo
   !   !$omp end parallel do
   !else
   ! The stage sum over `s` into q_pic(:,p) must stay the inner sequential loop;
   ! only the particle/variable loops may be parallelised. Folding `s` into the
   ! collapse races every thread on the same q_pic(v,p). Serially bit-identical.
   !$omp parallel do collapse(2) default(firstprivate) shared(q_pic,self)
   do p=1, particle_number
      do v=1, 6
         do s=1, nrk
            q_pic(v,p) = q_pic(v,p) + dt * self%beta(s) * self%q_pic_rk(v,p,s)
         enddo
      enddo
   enddo
   !$omp end parallel do
   !endif
   endassociate
   endsubroutine update_q_pic

   function crossproduct(a, b) result(cross)
   real(R8P), intent(in) :: a(3)     !< Left hand side.
   real(R8P), intent(in) :: b(3)     !< Left hand side.
   real(R8P)             :: cross(3) !< Cross product.

   cross(1) = (a(2) * b(3)) - (a(3) * b(2))
   cross(2) = (a(3) * b(1)) - (a(1) * b(3))
   cross(3) = (a(1) * b(2)) - (a(2) * b(1))
   endfunction crossproduct

endmodule adam_prism_rk_pic_object
