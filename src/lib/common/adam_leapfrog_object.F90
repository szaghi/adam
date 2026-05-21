!< ADAM, leapfrog class definition.
module adam_leapfrog_object
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
use adam_field_object, only : field_object
use adam_mpih_global,  only : mpih
use adam_grid_global,  only : grid
! third party modules
use finer
use penf

implicit none
save
private
public :: leapfrog_object

character(len=8), parameter :: INI_SECTION_NAME="leapfrog" !< INI (config) file section name containing time configs.

type :: leapfrog_object
   !< Leapforg class definition.
   real(R8P)                 :: nu=0.01_R8P        !< Robert-Asselin filter coefficient.
   real(R8P)                 :: alpha=0.53_R8P     !< Robert-Asselin-Williams filter coefficient.
   logical                   :: is_filtered=.false.!< Flag to check if the integration if RAW filtered.
   real(R8P), allocatable    :: q_old(:,:,:,:,:,:) !< Field cell centered variables, old time steps.
   ! grid data replica for easy handling
   integer(I4P),       pointer :: ngc=>null()           !< Number of ghost cells.
   integer(I4P),       pointer :: ni=>null()            !< Number of cells in i direction.
   integer(I4P),       pointer :: nj=>null()            !< Number of cells in j direction.
   integer(I4P),       pointer :: nk=>null()            !< Number of cells in k direction.
   contains
      ! public methods
      procedure, pass(self) :: assign_step    !< Assign q to old steps.
      procedure, pass(self) :: description    !< Return pretty-printed object description.
      procedure, pass(self) :: initialize     !< Initialize class.
      procedure, pass(self) :: integrate      !< Integrate.
      procedure, pass(self) :: load_from_file !< Load config from file.
endtype leapfrog_object
contains
   ! public methods
   subroutine assign_step(self, field, s, q, phi)
   !< Assign q to leapfrog old step.
   class(leapfrog_object), intent(inout)        :: self          !< Leapfrog object.
   type(field_object), intent(in) :: field !< Field (sibling realm component, threaded in).
   integer(I4P),           intent(in)           :: s             !< Current step number.
   real(R8P),              intent(in)           :: q(1:     ,     &
                                                     1-self%ngc:, &
                                                     1-self%ngc:, &
                                                     1-self%ngc:, &
                                                     1:)         !< Conservative variables.
   real(R8P),              intent(in), optional :: phi(1:,          &
                                                       1-self%ngc:, &
                                                       1-self%ngc:, &
                                                       1-self%ngc:, &
                                                       1:)       !< IB distance.
   integer(I4P)                                 :: all_solids    !< Last phi index, all solids summary.
   integer(I4P)                                 :: i, j, k, b, v !< Counter.

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>field%nv, blocks_number=>field%blocks_number)
   if (present(phi)) then
      all_solids = ubound(phi, dim=1)
      !$omp parallel do collapse(5) default(firstprivate) shared(phi,q,self)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
      do v=1, nv
         if (phi(all_solids,i,j,k,b) < 0._R8P) self%q_old(v,i,j,k,b,s) = q(v,i,j,k,b)
      enddo
      enddo
      enddo
      enddo
      enddo
   else
      !$omp parallel do collapse(5) default(firstprivate) shared(q,self)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
      do v=1, nv
         self%q_old(v,i,j,k,b,s) = q(v,i,j,k,b)
      enddo
      enddo
      enddo
      enddo
      enddo
   endif
   endassociate
   endsubroutine assign_step

   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(leapfrog_object), intent(in)  :: self             !< Leapfrog object.
   character(len=:),       allocatable :: desc             !< Description.
   character(len=1),       parameter   :: NL=new_line('a') !< New line character.
   integer(I4P)                        :: s                !< Counter.

   desc =       mpih%myrankstr//'Leapfrog scheme main data'//NL
   desc = desc//mpih%myrankstr//'  is RAW filtered: '//trim(str(self%is_filtered))//NL
   desc = desc//mpih%myrankstr//'  nu:              '//trim(str(self%nu         ))//NL
   desc = desc//mpih%myrankstr//'  alpha:           '//trim(str(self%alpha      ))
   endfunction description

   subroutine initialize(self, field, file_parameters, scheme)
   !< Initialize class.
   class(leapfrog_object),   intent(inout)        :: self            !< Leapfrog object.
   type(field_object), intent(in) :: field !< Field (sibling realm component, threaded in).
   type(file_ini),           intent(in), optional :: file_parameters !< Simulation parameters ini file handler.
   character(*),             intent(in), optional :: scheme          !< Runge-Kutta scheme.

   call mpih%print_message('leapfrog_object%initialize start')
   call associate_adam_data
   if (present(file_parameters)) then
      call self%load_from_file(file_parameters=file_parameters)
   endif
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>field%nv, nb=>field%nb)
   call allocate_variable(var=self%q_old,           &
                          ulb=reshape([1,nv,        &
                                       1-ngc,ni+ngc,&
                                       1-ngc,nj+ngc,&
                                       1-ngc,nk+ngc,&
                                       1,nb,        &
                                       1,2],[2,6]), &
                          msg=mpih%myrankstr//'leapfrog_object%initialize allocate q_old')
   endassociate
   print '(A)', self%description()
   call mpih%print_message('leapfrog_object%initialize finish')
   contains
      subroutine associate_adam_data
      !< Associate grid data pointers for easy handling.

      self%ni  => grid%ni
      self%nj  => grid%nj
      self%nk  => grid%nk
      self%ngc => grid%ngc
      endsubroutine associate_adam_data
   endsubroutine initialize

   subroutine integrate(self, field, dt, q, dq)
   !< Integrate.
   class(leapfrog_object), intent(inout) :: self          !< Leapfrog object.
   type(field_object), intent(in) :: field !< Field (sibling realm component, threaded in).
   real(R8P),              intent(in)    :: dt            !< Time step.
   real(R8P),              intent(inout) :: q(1:     ,     &
                                              1-self%ngc:, &
                                              1-self%ngc:, &
                                              1-self%ngc:, &
                                              1:)         !< Conservative variables.
   real(R8P),              intent(in)    :: dq(1:,         &
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1:)        !< Residuals.
   real(R8P)                             :: filter        !< Filter field displacement.
   integer(I4P)                          :: i, j, k, b, v !< Counter.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,nv=>field%nv,blocks_number=>field%blocks_number,q_old=>self%q_old)
   do b=1, blocks_number
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do v=1, nv
      q_old(v,i,j,k,b,2) = q_old(v,i,j,k,b,1) + 2._R8P * dt * dq(v,i,j,k,b)
      q_old(v,i,j,k,b,1) = q(v,i,j,k,b)
      q(v,i,j,k,b) = q_old(v,i,j,k,b,2)
   enddo
   enddo
   enddo
   enddo
   enddo
   if (self%is_filtered) then
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
      do v=1, nv
         filter = (q_old(v,i,j,k,b,1) - (q_old(v,i,j,k,b,2) * 2._R8P) + q(v,i,j,k,b)) * self%nu * 0.5_R8P
         q_old(v,i,j,k,b,2) = q_old(v,i,j,k,b,2) + (filter * self%alpha)
         q(v,i,j,k,b) = q(v,i,j,k,b) + (filter * (self%alpha - 1._R8P))
      enddo
      enddo
      enddo
      enddo
      enddo
   endif
   endassociate
   endsubroutine integrate

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(leapfrog_object), intent(inout)        :: self            !< Leapfrog object.
   type(file_ini),         intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,                intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                      :: go_on_fail_     !< Go on if load fails.
   integer(I4P)                                 :: error           !< Error status.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='is_filtered', val=self%is_filtered, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(is_filtered)')

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='nu', val=self%nu, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(nu)')

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='alpha', val=self%alpha, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(alpha)')
   endsubroutine load_from_file
endmodule adam_leapfrog_object
