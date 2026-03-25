!< ADAM, leapfrog class definition.
module adam_prism_leapfrog_pic_object
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

! ADAM singleton objects
use :: adam_mpih_global,      only : mpih
use :: adam_grid_global,      only : grid
use :: adam_prism_pic_object, only : prism_pic_object
! third party modules
use :: finer
use :: penf

implicit none
private
public :: prism_leapfrog_pic_object

character(len=8), parameter :: INI_SECTION_NAME="leapfrog" !< INI (config) file section name containing time configs.

type :: prism_leapfrog_pic_object
   !< Leapfrog class definition.
   real(R8P)                 :: nu=0.01_R8P         !< Robert-Asselin filter coefficient.
   real(R8P)                 :: alpha=0.53_R8P      !< Robert-Asselin-Williams filter coefficient.
   logical                   :: is_filtered=.false. !< Flag to check if the integration if RAW filtered.
   real(R8P), allocatable    :: q_pic_old(:,:,:)    !< Pic variables, old time steps.
   ! Prism data replica for easy handling
   type(prism_pic_object),   pointer :: pic=>null()              !< The PIC object.
   integer(I4P),             pointer :: particle_number=>null()  !< Number of particles.
   contains
      ! public methods
      procedure, pass(self) :: assign_step    !< Assign q to old steps.
      procedure, pass(self) :: description    !< Return pretty-printed object description.
      procedure, pass(self) :: initialize     !< Initialize class.
      procedure, pass(self) :: integrate      !< Integrate.
      procedure, pass(self) :: load_from_file !< Load config from file.
endtype prism_leapfrog_pic_object

contains
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(prism_leapfrog_pic_object), intent(in)  :: self             !< Leapfrog object.
   character(len=:),       allocatable 	       :: desc             !< Description.
   character(len=1),       parameter   	       :: NL=new_line('a') !< New line character.
   integer(I4P)                        	       :: s                !< Counter.

   desc =       mpih%myrankstr//'Leapfrog pic scheme main data'//NL
   desc = desc//mpih%myrankstr//'  is RAW filtered: '//trim(str(self%is_filtered))//NL
   desc = desc//mpih%myrankstr//'  nu:              '//trim(str(self%nu         ))//NL
   desc = desc//mpih%myrankstr//'  alpha:           '//trim(str(self%alpha      ))
   endfunction description

   subroutine initialize(self, file_parameters, scheme, pic)
   !< Initialize class.
   class(prism_leapfrog_pic_object),   intent(inout)        :: self            !< Leapfrog object.
   type(file_ini),           		      intent(in), optional :: file_parameters !< Simulation parameters ini file handler.
   character(*),             		      intent(in), optional :: scheme          !< Runge-Kutta scheme.
	type(prism_pic_object),   		      intent(in), target   :: pic             !< The PIC object.

   call mpih%print_message('leapfrog_pic_object%initialize start')
   call associate_adam_data(pic=pic)
   if (present(file_parameters)) then
      call self%load_from_file(file_parameters=file_parameters)
   endif
   associate(particle_number=>self%particle_number)!, ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>self%nv, nb=>self%nb)
   call allocate_variable(var=self%q_pic_old,        					 &
                          ulb=reshape([1,8,								 &
                                       1,particle_number,			 &
													1,2], [2,3]),					 &
                          					msg=mpih%myrankstr//'leapfrog_pic_object%initialize allocate q_pic_old')
   endassociate
   print '(A)', self%description()
   call mpih%print_message('leapfrog_pic_object%initialize finish')
   contains
      subroutine associate_adam_data(pic)
      !< Associate pic data pointers for easy handling.
		type(prism_pic_object), intent(in), target :: pic !< The PIC object.
		self%pic             => pic
		self%particle_number => pic%particle_number
      endsubroutine associate_adam_data
   endsubroutine initialize

	subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(prism_leapfrog_pic_object), intent(inout)        :: self            !< Leapfrog object.
   type(file_ini),         	 intent(in)                 :: file_parameters !< Simulation parameters ini file handler.
   logical,                	 intent(in), optional       :: go_on_fail      !< Go on if load fails.
   logical                     	                         :: go_on_fail_     !< Go on if load fails.
   integer(I4P)                	                         :: error           !< Error status.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='is_filtered', val=self%is_filtered, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(is_filtered)')

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='nu', val=self%nu, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(nu)')

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='alpha', val=self%alpha, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(alpha)')
   endsubroutine load_from_file

   ! public methods
   subroutine assign_step(self, s, q_pic, phi)
   !< Assign q to leapfrog old step.
   class(prism_leapfrog_pic_object), intent(inout)        :: self           !< Leapfrog object.
   integer(I4P),           	       intent(in)           :: s              !< Current step number.
   real(R8P),              	       intent(in)           :: q_pic(1:,1:)   !< Pic variables.
   real(R8P),              	       intent(in), optional :: phi(1:,      &
                           	                               1-grid%ngc:, &
                           	                               1-grid%ngc:, &
                           	                               1-grid%ngc:, &
                           	                               1:)            !< IB distance.
   integer(I4P)            	                            :: all_solids     !< Last phi index, all solids summary.
   integer(I4P)            	                            :: i, j, k, b, v  !< Counter.

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

   subroutine integrate(self, dt, q_pic, pic_fields)
   !< Integrate.
   class(prism_leapfrog_pic_object), intent(inout) :: self               !< Leapfrog object.
   real(R8P),              	 intent(in)          :: dt                 !< Time step.
   real(R8P),              	 intent(inout)       :: q_pic(1:, 1:)      !< Pic variables.
   real(R8P),              	 intent(in)          :: pic_fields(1:, 1:) !< Fields value at particle locations.
   real(R8P)               	                     :: filter             !< Filter field displacement.
   integer(I4P)            	                     :: p, v			       !< Counter.
   real(R8P)                                       :: v_star(3)          !< Auxiliary velocity v* Buneman-Boris
   real(R8P)                                       :: v_star_star(3)     !< Auxiliary velocity v** Buneman-Boris
   real(R8P)                                       :: t(3), s(3), w(3)   !< Auxiliary vector t, s, w Buneman-Boris

   associate(particle_number=>self%particle_number, q_pic_old=>self%q_pic_old)

   ! In ingresso ho: q_pic (velocità e posizione)            al tempo n
   !                 pic_fields (campi elettromagnetici)     al tempo n
   !                 q_pic_old(:,:,1) (velocità e posizione) al tempo n-1
   !                 q_pic_old(:,:,2) (velocità e posizione) al tempo n (tranne nella prima iterazione
   !                                                                     in cui è 0 da inizializzazione)

   do p=1, particle_number

      !Integrazione vettore velocità con schema di Buneman-Boris
      v_star(1) = q_pic_old(4,p,1) + dt * pic_fields(1,p) * q_pic(7,p) / q_pic(8,p)
      v_star(2) = q_pic_old(5,p,1) + dt * pic_fields(2,p) * q_pic(7,p) / q_pic(8,p)
      v_star(3) = q_pic_old(6,p,1) + dt * pic_fields(3,p) * q_pic(7,p) / q_pic(8,p)

      t(1) = dt * q_pic(7,p) / q_pic(8,p) * pic_fields(4,p)
      t(2) = dt * q_pic(7,p) / q_pic(8,p) * pic_fields(5,p)
      t(3) = dt * q_pic(7,p) / q_pic(8,p) * pic_fields(6,p)

      w = v_star + crossproduct(v_star, t)

      s = 2/(1._R8P + sq_norm(t)) * t

      v_star_star = v_star + crossproduct(w, s)

      !Aggiornamento vettore di appoggio e vettore velocità
      q_pic_old(4:6,p,1) = q_pic(4:6,p) !Salvo la velocità al tempo n, che userò nell'integrazione al tempo successivo e
                                        !nell'integrazione delle posizioni delle particelle

      q_pic_old(4:6,p,2) = v_star_star + dt * pic_fields(1:3,p) * q_pic(7,p) / q_pic(8,p) !Integro la velocità
                                                                                                    !al tempo n+1

      q_pic(4:6,p)       = q_pic_old(4:6,p,2) !Velocità al tempo n+1

      !Aggiornamento posizioni con schema leapfrog
      do v=1, 3
         q_pic_old(v,p,2) = q_pic_old(v,p,1) + 2._R8P * dt * q_pic_old(v+3_I4P,p,1) !Integro la posizione al tempo n+1
         q_pic_old(v,p,1) = q_pic(v,p) !Salvo la posizione al tempo n, che userò nell'integrazione al tempo successivo
         q_pic(v,p)       = q_pic_old(v,p,2) !Posizione al tempo n+1
      enddo
   enddo
   ! In uscita ho:   q_pic (velocità e posizione)            al tempo n+1
   !                 pic_fields (campi elettromagnetici)     al tempo n
   !                 q_pic_old(:,:,1) (velocità e posizione) al tempo n
   !                 q_pic_old(:,:,2) (velocità e posizione) al tempo n+1

   ! Lascio così in analogia a leapfrog_object, ma da verificare la sua implementazione
   if (self%is_filtered) then
      do p=1, particle_number
      	do v=1, 8
      	   filter = (q_pic_old(v,p,1) - (q_pic_old(v,p,2) * 2._R8P) + q_pic(v,p)) * self%nu * 0.5_R8P
      	   q_pic_old(v,p,2) = q_pic_old(v,p,2) + (filter * self%alpha)
      	   q_pic(v,p) = q_pic(v,p) + (filter * (self%alpha - 1._R8P))
      	enddo
      enddo
   endif
   endassociate
   endsubroutine integrate

   function dotproduct(a, b) result(dot)
   !< Compute the scalar (dot) product.
   real(R8P), intent(in) :: a(3) !< Left hand side.
   real(R8P), intent(in) :: b(3) !< Left hand side.
   real(R8P)             :: dot  !< Dot product.

   dot = (a(1) * b(1)) + (a(2) * b(2)) + (a(3) * b(3))
   endfunction dotproduct

   function crossproduct(a, b) result(cross)
   real(R8P), intent(in) :: a(3)     !< Left hand side.
   real(R8P), intent(in) :: b(3)     !< Left hand side.
   real(R8P)             :: cross(3) !< Cross product.

   cross(1) = (a(2) * b(3)) - (a(3) * b(2))
   cross(2) = (a(3) * b(1)) - (a(1) * b(3))
   cross(3) = (a(1) * b(2)) - (a(2) * b(1))
   endfunction crossproduct

   function sq_norm(a) result(sq)
   !< Return the square of the norm of vector.
   real(R8P), intent(in)  :: a(3)     !< Input vector
   real(R8P)              :: sq       !< Square norm of input

   sq = (a(1) * a(1)) + (a(2) * a(2)) + (a(3) * a(3))
   endfunction sq_norm

endmodule adam_prism_leapfrog_pic_object
