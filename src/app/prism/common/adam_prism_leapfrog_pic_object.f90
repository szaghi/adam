!< ADAM, leapfrog class definition.
module adam_leapfrog_pic_object
!< ADAM, leapfrog class definition.

!< Considering the following ODE system:
!<
!< $$ U_t = R(t,U) $$
!<
!< where \(U_t = \frac{dU}{dt}\), *U* is the vector of *state* variables being a function of the time-like independent variable
!< *t*, *R* is the (vectorial) residual function, the leapfrog class scheme implemented (see [3]) is:
!<
!< $$ U^{n+2} = U^{n} + 2\Delta t \cdot R(t^{n+1}, U^{n+1}) $$
!<
!< Optionally, the Robert-Asselin-Williams (RAW) filter (see [3]) is applied to the computed integration steps:
!< $$ \Delta = \frac{\nu}{2}(U^{n} - 2 U^{n+1} + U^{n+2}) $$
!< $$ U^{n+1} = U^{n+1} + \Delta * \alpha $$
!< $$ U^{n+2} = U^{n+2} + \Delta * (\alpha-1) $$
!< Note that for \(\alpha=1\) the filter reverts back to the standard Robert-Asselin scheme.
!< The filter coefficients should be taken as \(\nu \in (0,1]\) and \(\alpha \in (0.5,1]\). The default values are
!<
!<  + \(\nu=0.01\)
!<  + \(\alpha=0.53\)
!<
!< @note The value of \(\Delta t\) must be provided, it not being computed by the integrator.
!<
!< The schemes are explicit. The filter coefficients \(\nu,\,\alpha \) define the actual scheme.
!<
!<#### Bibliography
!<
!< [1] *The integration of a low order spectral form of the primitive meteorological equations*, Robert, A. J., J. Meteor. Soc.
!< Japan,vol. 44, pages 237--245, 1966.
!<
!< [2] *Frequency filter for time integrations*, Asselin, R., Monthly Weather Review, vol. 100, pages 487--490, 1972.
!<
!< [3] *The RAW filter: An improvement to the Robert–Asselin filter in semi-implicit integrations*, Williams, P.D., Monthly
!< Weather Review, vol. 139(6), pages 1996--2007, June 2011.

use adam_field_object
use adam_grid_object
use adam_prism_pic_object
use adam_mpih_object
use finer
use penf

implicit none
save
private
public :: leapfrog_pic_object

character(len=8), parameter :: INI_SECTION_NAME="leapfrog" !< INI (config) file section name containing time configs.

type :: leapfrog_pic_object
   !< Leapforg class definition.
   type(mpih_object)         :: mpih                !< MPI handler.
   real(R8P)                 :: nu=0.01_R8P         !< Robert-Asselin filter coefficient.
   real(R8P)                 :: alpha=0.53_R8P      !< Robert-Asselin-Williams filter coefficient.
   logical                   :: is_filtered=.false. !< Flag to check if the integration if RAW filtered.
   real(R8P), allocatable    :: q_pic_old(:,:,:)    !< Pic variables, old time steps.
   ! Adam data replica for easy handling
   type(pic_object),   pointer :: pic=>null()              !< The PIC object.
   type(field_object), pointer :: field=>null()            !< The field.
   type(grid_object),  pointer :: grid=>null()             !< The grid.
   integer(I4P),       pointer :: particles_number=>null() !< Number of particles.
   integer(I4P),       pointer :: ngc=>null()              !< Number of ghost cells.
   integer(I4P),       pointer :: ni=>null()               !< Number of cells in i direction.
   integer(I4P),       pointer :: nj=>null()               !< Number of cells in j direction.
   integer(I4P),       pointer :: nk=>null()               !< Number of cells in k direction.
   integer(I4P),       pointer :: nb=>null()               !< Total blocks number for MPI.
   integer(I4P),       pointer :: blocks_number=>null()    !< Actual blocks number.
   integer(I4P),       pointer :: ns=>null()               !< Number of fluids specie.
   integer(I4P),       pointer :: nv=>null()               !< Number of conservative variables.
   contains
      ! public methods
      procedure, pass(self) :: assign_step    !< Assign q to old steps.
      procedure, pass(self) :: description    !< Return pretty-printed object description.
      procedure, pass(self) :: initialize     !< Initialize class.
      procedure, pass(self) :: integrate      !< Integrate.
      procedure, pass(self) :: load_from_file !< Load config from file.
endtype leapfrog_pic_object
contains
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(leapfrog_pic_object), intent(in)  :: self             !< Leapfrog object.
   character(len=:),       allocatable 	 :: desc             !< Description.
   character(len=1),       parameter   	 :: NL=new_line('a') !< New line character.
   integer(I4P)                        	 :: s                !< Counter.

   desc =       self%mpih%myrankstr//'Leapfrog pic scheme main data'//NL
   desc = desc//self%mpih%myrankstr//'  is RAW filtered: '//trim(str(self%is_filtered))//NL
   desc = desc//self%mpih%myrankstr//'  nu:              '//trim(str(self%nu         ))//NL
   desc = desc//self%mpih%myrankstr//'  alpha:           '//trim(str(self%alpha      ))
   endfunction description

   subroutine initialize(self, file_parameters, scheme, grid, field)
   !< Initialize class.
   class(leapfrog_pic_object),   intent(inout)        :: self            !< Leapfrog object.
   type(file_ini),           		intent(in), optional :: file_parameters !< Simulation parameters ini file handler.
   character(*),             		intent(in), optional :: scheme          !< Runge-Kutta scheme.
   type(grid_object),        	 	intent(in), target   :: grid            !< The grid.
   type(field_object),       	 	intent(in), target   :: field           !< The field.
	type(pic_object),         		intent(in), target   :: pic             !< The PIC object.

   call self%mpih%initialize(do_mpi_init=.false.)
   call self%mpih%print_message('leapfrog_pic_object%initialize start')
   call associate_adam_data(pic=pic) !grid=grid, field=field)
   if (present(file_parameters)) then
      call self%load_from_file(file_parameters=file_parameters)
   endif
   associate(particle_number=>self%particles_number)!, ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>self%nv, nb=>self%nb)
   call allocate_variable(var=self%q_pic_old,        					 &
                          ulb=reshape([1,particle_number,			 &
                                       1,8,								 &
													1,2], [2,3]),					 &
                          					msg=self%mpih%myrankstr//'leapfrog_pic_object%initialize allocate q_pic_old')
   endassociate
   print '(A)', self%description()
   call self%mpih%print_message('leapfrog_pic_object%initialize finish')
   contains
      subroutine associate_adam_data(grid, field, pic)
      !< Associate objects data to equation for easy handling.
      type(grid_object),          intent(in), target :: grid    !< The grid.
      type(field_object),         intent(in), target :: field   !< The field.
		type(pic_object),           intent(in), target :: pic     !< The PIC object.
      self%field         	 => field
      self%blocks_number 	 => field%blocks_number
      self%ni            	 => field%grid%ni
      self%nj            	 => field%grid%nj
      self%nk            	 => field%grid%nk
      self%ngc           	 => field%grid%ngc
      self%nb            	 => field%nb
      self%nv            	 => field%nv
		self%pic           	 => pic
		self%particles_number => pic%particles_number
      endsubroutine associate_adam_data
   endsubroutine initialize

	subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(leapfrog_pic_object), intent(inout)        :: self            !< Leapfrog object.
   type(file_ini),         	 intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,                	 intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                     	                   :: go_on_fail_     !< Go on if load fails.
   integer(I4P)                	                   :: error           !< Error status.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='is_filtered', val=self%is_filtered, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(is_filtered)')

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='nu', val=self%nu, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(nu)')

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='alpha', val=self%alpha, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(alpha)')
   endsubroutine load_from_file

   ! public methods
   subroutine assign_step(self, s, q_pic, phi)
   !< Assign q to leapfrog old step.
   class(leapfrog_pic_object), intent(inout)        :: self          !< Leapfrog object.
   integer(I4P),           	 intent(in)           :: s             !< Current step number.
   real(R8P),              	 intent(in)           :: q_pic(1:,1:)  !< Pic variables.
   real(R8P),              	 intent(in), optional :: phi(1:,          &
                           	                             1-self%ngc:, &
                           	                             1-self%ngc:, &
                           	                             1-self%ngc:, &
                           	                             1:)       !< IB distance.
   integer(I4P)            	                      :: all_solids    !< Last phi index, all solids summary.
   integer(I4P)            	                      :: i, j, k, b, v !< Counter.

   !associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>self%nv, blocks_number=>self%blocks_number)
   !if (present(phi)) then
   !   all_solids = ubound(phi, dim=1)
   !   !$omp parallel do collapse(5) default(firstprivate) shared(phi,q,self)
   !   do b=1, blocks_number
   !   do k=1, nk
   !   do j=1, nj
   !   do i=1, ni
   !   do v=1, nv
   !      if (phi(all_solids,i,j,k,b) < 0._R8P) self%q_old(v,i,j,k,b,s) = q(v,i,j,k,b)
   !   enddo
   !   enddo
   !   enddo
   !   enddo
   !   enddo
   !else
   !   !$omp parallel do collapse(5) default(firstprivate) shared(q,self)
   !   do b=1, blocks_number
   !   do k=1, nk
   !   do j=1, nj
   !   do i=1, ni
   !   do v=1, nv
   !      self%q_old(v,i,j,k,b,s) = q(v,i,j,k,b)
   !   enddo
   !   enddo
   !   enddo
   !   enddo
   !   enddo
   !endif
   !endassociate
	self%q_pic_old(:, :, s) = q_pic(:, :)
   endsubroutine assign_step

   subroutine integrate(self, dt, q_pic, dq_pic)
   !< Integrate.
   class(leapfrog_pic_object), intent(inout) :: self           !< Leapfrog object.
   real(R8P),              	 intent(in)    :: dt             !< Time step.
   real(R8P),              	 intent(inout) :: q_pic(1:, 1:)  !< Pic variables.
   real(R8P),              	 intent(in)    :: dq_pic(1:, 1:) !< Pic residuals.
   real(R8P)               	               :: filter         !< Filter field displacement.
   integer(I4P)            	               :: p, v			  !< Counter.

   associate(p=>self%particle_number, q_pic_old=>self%q_old)
   do p=1, particle_number
   	do v=1, 8
   	   q_pic_old(p,v,2) = q_pic_old(p,v,1) + 2._R8P * dt * dq_pic(p,v)
   	   q_pic_old(p,v,1) = q_pic(p,v)
   	   q_pic(p,v) = q_pic_old(p,v,2)
   	enddo
   enddo
   if (self%is_filtered) then
      do p=1, particle_number
      	do v=1, 8
      	   filter = (q_old(p,v,1) - (q_old(p,v,2) * 2._R8P) + q(p,v)) * self%nu * 0.5_R8P
      	   q_old(p,v,2) = q_old(p,v,2) + (filter * self%alpha)
      	   q(p,v) = q(p,v) + (filter * (self%alpha - 1._R8P))
      	enddo
      enddo
   endif
   endassociate
   endsubroutine integrate

endmodule adam_leapfrog_pic_object