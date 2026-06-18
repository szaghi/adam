!< ADAM, PRISM Particle-in-Cell class definition, CPU backend.
module adam_prism_pic_object
!< ADAM, PRISM Particle-in-Cell class definition, CPU backend.

! ADAM classes, libraries, parameters
use :: adam_field_object, only : field_object
use :: adam_grid_object,  only : grid_object
! ADAM singleton objects
use :: adam_mpih_global,  only : mpih
! PRISM modules
use :: adam_prism_parameters
! third party modules
use :: finer, only : file_ini
use :: penf,  only : I4P, R8P, str

implicit none
private
public :: INI_SECTION_NAME
public :: prism_pic_object
!public :: particle_weighting_interface
!public :: current_weighting_interface
!public :: field_weighting_interface
!public :: particle_weighting
!public :: current_weighting
!public :: field_weighting
public :: PLASMA_TYPE_PROBLEM
public :: SINGLE_PARTICLE_TYPE_PROBLEM
!public :: CIC_WEIGHTING_MODEL
!public :: NGP_WEIGHTING_MODEL
!public :: TSC_WEIGHTING_MODEL
!public :: ZEROD_FIELDS_WEIGHTING_MODEL
!public :: ONED_FIELDS_WEIGHTING_MODEL
public :: NUM_SCHEME_TIME_PIC_LEAPFROG
public :: NUM_SCHEME_TIME_PIC_RUNGE_KUTTA
public :: STANDARD_INITIALIZATION
public :: COHERENT_INITIALIZATION
public :: NEUMANN_BC
public :: DIRICHLET_BC
public :: ANALYTIC_BC
!public :: CIC_charge_weighting
!public :: NGP_charge_weighting
!public :: TSC_charge_weighting
!public :: CIC_current_weighting
!public :: NGP_current_weighting
!public :: TSC_current_weighting
!public :: zeroD_field_weighting
!public :: oneD_field_weighting

character(len=3 ), parameter :: INI_SECTION_NAME                 = 'PIC'             !< INI file section name for PIC configuration.
character(len=6 ), parameter :: PLASMA_TYPE_PROBLEM              = 'plasma'          !< Analyzing physical problem involving the presence of plasma
character(len=8 ), parameter :: STANDARD_INITIALIZATION          = 'standard'        !< Field initialization through elliptic solver wuth standard laplacian scheme (7-points)
character(len=8 ), parameter :: COHERENT_INITIALIZATION          = 'coherent'        !< Field initialization through elliptic solver wuth laplacian scheme coherent with centerd difference scheme (implemented only for 6th order)
character(len=7 ), parameter :: NEUMANN_BC                       = 'neumann'         !< Neumann boundary condition
character(len=9 ), parameter :: DIRICHLET_BC                     = 'dirichlet'       !< Dirichlet boundary condition
character(len=5 ), parameter :: ANALYTIC_BC                      = 'analytic'        !< Analytic boundary condition
character(len=15), parameter :: SINGLE_PARTICLE_TYPE_PROBLEM     = 'single_particle' !< Analyzing physical problem involving the presence of a single particle
character(len=3 ), parameter :: NGP_WEIGHTING_MODEL              = 'NGP'             !< NGP weighting model.
character(len=3 ), parameter :: CIC_WEIGHTING_MODEL              = 'CIC'             !< CIC weighting model.
character(len=3 ), parameter :: TSC_WEIGHTING_MODEL              = 'TSC'             !< TSC weighting model.
character(len=8 ), parameter :: GAUSSIAN_WEIGHTING_MODEL         = 'Gaussian'        !< Gaussian weighting model.
character(len=2 ), parameter :: ZEROD_FIELDS_WEIGHTING_MODEL     = '0D'              !< 0D field weighting.
character(len=2 ), parameter :: ONED_FIELDS_WEIGHTING_MODEL      = '1D'              !< 1D field weighting.
character(len=8 ), parameter :: NUM_SCHEME_TIME_PIC_LEAPFROG     = 'LEAPFROG'        !< Leapfrog numerical scheme for time operator.
character(len=11), parameter :: NUM_SCHEME_TIME_PIC_RUNGE_KUTTA  = 'RUNGE_KUTTA'     !< Runge-Kutta numerical scheme for time operator.

! PIC variables layout in q_pic array:
!q_pic(1) = x
!q_pic(2) = y
!q_pic(3) = z
!q_pic(4) = vx
!q_pic(5) = vy
!q_pic(6) = vz
!q_pic(7) = charge
!q_pic(8) = mass

type :: prism_pic_object
   real(R8P)                 :: plasma_density              !< Plasma density
   real(R8P)                 :: neutral_fraction = 0.0_R8P  !< Neutral fraction
   real(R8P)                 :: sigma = 0.0_R8P             !< Standard deviation for a gaussian weighting
   real(R8P)                 :: cutoff_sigma = 0._R8P       !< Gaussian cutoff.
   integer(I4P)              :: particle_number  = 0_I4P    !< Total number of particles.
	integer(I4P)				  :: n_ions = 0_I4P              !< Total ions number
	integer(I4P)				  :: n_electrons = 0_I4P         !< Total electrons number
	integer(I4P)				  :: n_neutrals = 0_I4P          !< Total neutrals number
   integer(I4P), allocatable :: neighbour_list(:,:)         !< Particle grid positions array.
   character(len=99)         :: problem_type                !< Type of problem analyzed
   character(len=99)         :: initialization              !< field initialization solver
   character(len=99)         :: bc_solver                   !< boundary conditions for the elliptic solver for initial conditions
   logical                   :: elliptic_correction=.false. !< elliptic correction for the initial fields
   character(len=99)         :: bc_correction               !< boundary conditions for the elliptic correction
   character(len=99)         :: particle_weighting_model    !< Particle weighting model.
   character(len=99)         :: current_weighting_model     !< Current weighting model.
   character(len=99)         :: field_weighting_model       !< Field weighting model.
   character(len=99)         :: scheme_time                 !< Numerical scheme for time operator [runge_kutta, leapfrog,...].
   !< Pointer (abstract) TBP.
   procedure(particle_weighting_interface), pass(self), pointer :: particle_weighting =>null() !< Particle weighting.
   procedure(current_weighting_interface),  pass(self), pointer :: current_weighting  =>null() !< Current weighting.
   procedure(field_weighting_interface),    pass(self), pointer :: field_weighting    =>null() !< field weighting.
contains
   procedure, pass(self) :: description                   !< Return pretty-printed object description.
   procedure, pass(self) :: initialize                    !< Initialize IC.
   procedure, pass(self) :: load_from_file                !< Load config from file.
   procedure, pass(self) :: particle_cartesian_grid_index !< Compute the grid index corresponding to a particle position.
   procedure, pass(self) :: CIC_charge_weighting          !< Cloud-in-Cell weighting of particle quantities to the grid.
   procedure, pass(self) :: NGP_charge_weighting          !< Nearest Grid Point weighting of particle quantities to the grid.
   procedure, pass(self) :: TSC_charge_weighting          !< Triangular Shaped Cloud weighting of particle quantities to the grid.
   procedure, pass(self) :: Gaussian_charge_weighting     !< Triangular Shaped Cloud weighting of particle quantities to the grid.
   procedure, pass(self) :: CIC_current_weighting         !< Cloud-in-Cell weighting of particle quantities to the grid.
   procedure, pass(self) :: NGP_current_weighting         !< Nearest Grid Point weighting of particle quantities to the grid.
   procedure, pass(self) :: TSC_current_weighting         !< Triangular Shaped Cloud weighting of particle quantities to the grid.
   procedure, pass(self) :: zeroD_field_weighting
   procedure, pass(self) :: oneD_field_weighting
endtype prism_pic_object

interface
   subroutine particle_weighting_interface(self, field, grid, q, q_pic, nv)
   import :: prism_pic_object, grid_object, field_object, I4P, R8P
   class(prism_pic_object), intent(inout) :: self                 !< External fields.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),       intent(in)    :: grid                 !< The grid.
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:, &
                                               1-grid%ngc:,    &
                                               1-grid%ngc:,1:)    !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Number of variables.
   endsubroutine particle_weighting_interface

   subroutine current_weighting_interface(self, field, grid, q, q_pic, nv)
   import :: prism_pic_object, grid_object, field_object, I4P, R8P
   class(prism_pic_object), intent(inout) :: self                 !< External fields.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),       intent(in)    :: grid                 !< The grid.
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:, &
                                               1-grid%ngc:,    &
                                               1-grid%ngc:,1:)    !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Number of variables.
   endsubroutine current_weighting_interface

   subroutine field_weighting_interface(self, field, grid, pic_fields, q, q_pic, nv)
   import :: prism_pic_object, grid_object, field_object, I4P, R8P
   class(prism_pic_object), intent(inout) :: self                 !< External fields.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),       intent(in)    :: grid                 !< The grid.
   real(R8P),               intent(inout) :: pic_fields(1:,1:)    !< Fields value at particle locations.
   real(R8P),               intent(in)    :: q(1:,1-grid%ngc:, &
                                               1-grid%ngc:,    &
                                               1-grid%ngc:,1:)    !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Number of variables.
   endsubroutine field_weighting_interface
endinterface

contains
   function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(prism_pic_object), intent(in) :: self             !< External fields.
   character(len=:), allocatable       :: desc             !< Description.
   character(len=1), parameter         :: NL=new_line('a') !< New line character.

   desc =            mpih%myrankstr//'PIC object description:'
   desc = desc//NL//mpih%myrankstr//'    Problem type: '//trim(self%problem_type)
   desc = desc//NL//mpih%myrankstr//'    Initialization type: '//trim(self%initialization)
   desc = desc//NL//mpih%myrankstr//'    Bc elliptic solver: '//trim(self%bc_solver)
   desc = desc//NL//mpih%myrankstr//'    Elliptic correction: '//trim(str(self%elliptic_correction))
   if (self%elliptic_correction) then
      desc = desc//NL//mpih%myrankstr//'    Bc elliptic correction: '//trim(self%bc_correction)
   endif
   if (self%problem_type == PLASMA_TYPE_PROBLEM) then
      desc = desc//NL//mpih%myrankstr//'    Input plasma density [m^(-3)]: '//trim(str(self%plasma_density))
      desc = desc//NL//mpih%myrankstr//'    Neutral fraction: '//trim(str(self%neutral_fraction))
      !desc = desc//NL//mpih%myrankstr//'    Total number of particles: '//trim(str(self%particle_number))
      !desc = desc//NL//mpih%myrankstr//'    of which ions: '//trim(str(self%n_ions))
      !desc = desc//NL//mpih%myrankstr//'    of which electrons: '//trim(str(self%n_electrons))
      !desc = desc//NL//mpih%myrankstr//'    of which neutrals: '//trim(str(self%n_neutrals))
   endif
   desc = desc//NL//mpih%myrankstr//'    Particle weighting model: '//trim(self%particle_weighting_model)
   desc = desc//NL//mpih%myrankstr//'    Current weighting model: '//trim(self%current_weighting_model)
   if (self%particle_weighting_model == GAUSSIAN_WEIGHTING_MODEL .or. &
         self%current_weighting_model == GAUSSIAN_WEIGHTING_MODEL) then
      desc = desc//NL//mpih%myrankstr//'    Sigma: '//trim(str(self%sigma))
      desc = desc//NL//mpih%myrankstr//'    Cutoff_sigma: '//trim(str(self%cutoff_sigma))
   endif
   desc = desc//NL//mpih%myrankstr//'    Field weighting model: '//trim(self%field_weighting_model)
   desc = desc//NL//mpih%myrankstr//'    Numerical scheme for time operator: '//trim(self%scheme_time)
   endfunction description

   subroutine initialize(self, field, grid, file_parameters)
   !< Initialize PIC.
   class(prism_pic_object), intent(inout) :: self            !< Pic object.
   type(field_object),      intent(in)    :: field           !< Field (sibling realm component, threaded in).
   type(grid_object),       intent(in)    :: grid            !< Grid (sibling realm component, threaded in).
   type(file_ini),          intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   real(R8P)                              :: domain_volume   !< Total volume of the computational domain where plasma is
                                                             !< present at t0
	character(len=:),        allocatable   :: desc
   character(len=1),        parameter     :: NL=new_line('a')

   print '(A)', mpih%myrankstr//'prism_pic_object%initialize start'

   call self%load_from_file(file_parameters=file_parameters)
   associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj, nk=>grid%nk, &
               emin=>grid%domain_emin, emax=>grid%domain_emax)

   if (self%problem_type == PLASMA_TYPE_PROBLEM) then
      domain_volume = (emax(1)-emin(1))*(emax(2)-emin(2))*(emax(3)-emin(3))
      self%particle_number = nint(self%plasma_density*domain_volume)
      self%n_neutrals = nint(self%neutral_fraction*real(self%particle_number,R8P))
	   self%n_ions = nint(real(self%particle_number-self%n_neutrals, R8P)/2.0_R8P)
	   self%n_electrons = self%n_ions
	   self%n_neutrals = self%particle_number - self%n_ions - self%n_electrons
   elseif (self%problem_type == SINGLE_PARTICLE_TYPE_PROBLEM) then
      self%particle_number = 1_I4P
   endif
   endassociate

   print '(A)', self%description()

   allocate(self%neighbour_list(4, self%particle_number))

   select case(self%particle_weighting_model)
   case(CIC_WEIGHTING_MODEL)
      self%particle_weighting => CIC_charge_weighting
   case(NGP_WEIGHTING_MODEL)
      self%particle_weighting => NGP_charge_weighting
   case(TSC_WEIGHTING_MODEL)
      self%particle_weighting => TSC_charge_weighting
   case(GAUSSIAN_WEIGHTING_MODEL)
      self%particle_weighting => Gaussian_charge_weighting
   case default
      call mpih%error_stop(msg=': invalid particle weighting model in prism_cpu_object%initialize')
   endselect

   select case(self%current_weighting_model)
   case(CIC_WEIGHTING_MODEL)
      self%current_weighting => CIC_current_weighting
   case(NGP_WEIGHTING_MODEL)
      self%current_weighting => NGP_current_weighting
   case(TSC_WEIGHTING_MODEL)
      self%current_weighting => TSC_current_weighting
   case default
      call mpih%error_stop(msg=': invalid current weighting model in prism_cpu_object%initialize')
   endselect

   select case(self%field_weighting_model)
   case(ZEROD_FIELDS_WEIGHTING_MODEL)
      self%field_weighting => zeroD_field_weighting
   case(ONED_FIELDS_WEIGHTING_MODEL)
      self%field_weighting => oneD_field_weighting
   !case(TSC_WEIGHTING_MODEL)
   !   current_weighting => TSC_current_weighting
   case default
      call mpih%error_stop(msg=': invalid field weighting model in prism_cpu_object%initialize')
   endselect

   print '(A)', mpih%myrankstr//'prism_pic_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load PIC configuration from file.
	class(prism_pic_object), intent(inout)   			 :: self             !< PIC object.
	type(file_ini),          intent(in)		  			 :: file_parameters  !< File handler.
   logical,                 intent(in), optional    :: go_on_fail      	!< Go on if load fails.
   logical                                          :: go_on_fail_     	!< Go on if load fails.
   integer(I4P)                                     :: error           	!< Error status.
   character(99)                                    :: buff       		!< Option character buffer.

	go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

	call file_parameters%get(section_name=INI_SECTION_NAME, option_name='particle_weighting_model', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(particle_weighting_model) from file')
   select case(trim(adjustl(buff)))
   case('CIC', 'cic', 'Cic')
      self%particle_weighting_model = CIC_WEIGHTING_MODEL
	case('NGP', 'ngp', 'Ngp')
		self%particle_weighting_model = NGP_WEIGHTING_MODEL
	case('TSC', 'tsc', 'Tsc')
		self%particle_weighting_model = TSC_WEIGHTING_MODEL
   case('GAUSSIAN', 'Gaussian', 'gaussian', 'GAUSS', 'gauss')
      self%particle_weighting_model = GAUSSIAN_WEIGHTING_MODEL
	case default
		call mpih%error_stop(msg=': invalid particle weighting model ['//trim(adjustl(buff))//'] in  &
      ['//INI_SECTION_NAME//'].(particle_weighting_model)')
	endselect

	call file_parameters%get(section_name=INI_SECTION_NAME, option_name='current_weighting_model', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(current_weighting_model) from file')
   select case(trim(adjustl(buff)))
   case('CIC', 'cic', 'Cic')
      self%current_weighting_model = CIC_WEIGHTING_MODEL
	case('NGP', 'ngp', 'Ngp')
		self%current_weighting_model = NGP_WEIGHTING_MODEL
	case('TSC', 'tsc', 'Tsc')
		self%current_weighting_model = TSC_WEIGHTING_MODEL
	case default
		call mpih%error_stop(msg=': invalid current weighting model ['//trim(adjustl(buff))//'] in  &
      ['//INI_SECTION_NAME//'].(current_weighting_model)')
	endselect

   if (self%particle_weighting_model == GAUSSIAN_WEIGHTING_MODEL .or. &
         self%current_weighting_model == GAUSSIAN_WEIGHTING_MODEL) then
         call file_parameters%get(section_name=INI_SECTION_NAME, option_name='sigma', &
         val=self%sigma, error=error)
         if (.not.go_on_fail_.and.error>0) &
         call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(sigma)')

         call file_parameters%get(section_name=INI_SECTION_NAME, option_name='cutoff_sigma', &
         val=self%cutoff_sigma, error=error)
         if (.not.go_on_fail_.and.error>0) &
         call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(cutoff_sigma)')
   endif

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='field_weighting_model', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(field_weighting_model) from file')
   select case(trim(adjustl(buff)))
   case('0D', '0d', '0_d', '0_D')
      self%field_weighting_model = ZEROD_FIELDS_WEIGHTING_MODEL
   case('1D', '1d', '1_d', '1_D')
      self%field_weighting_model = ONED_FIELDS_WEIGHTING_MODEL
   case default
      call mpih%error_stop(msg=': invalid field weighting model ['//trim(adjustl(buff))//'] in  &
      ['//INI_SECTION_NAME//'].(field_weighting_model)')
   endselect

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='scheme_time', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(scheme_time)')
   select case(trim(adjustl(buff)))
   case('LEAPFROG', 'leapfrog', 'Leapfrog')
      self%scheme_time = NUM_SCHEME_TIME_PIC_LEAPFROG
   case('RUNGE_KUTTA', 'runge_kutta', 'Runge_Kutta')
      self%scheme_time = NUM_SCHEME_TIME_PIC_RUNGE_KUTTA
   case default
      call mpih%print_message(msg='warning: numerical scheme "'//trim(adjustl(buff))//'" unknown. Revert back to RK scheme')
      self%scheme_time = NUM_SCHEME_TIME_PIC_RUNGE_KUTTA
   endselect

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='problem_type', &
   val=self%problem_type, error=error)
   if (.not.go_on_fail_.and.error>0) &
   call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(problem_type)')

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='initialization', &
   val=self%initialization, error=error)
   if (.not.go_on_fail_.and.error>0) &
   call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(initialization)')

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='bc_solver', &
   val=self%bc_solver, error=error)
   if (.not.go_on_fail_.and.error>0) &
   call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(bc_solver)')

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='elliptic_correction', &
   val=self%elliptic_correction, error=error)
   if (.not.go_on_fail_.and.error>0) &
   call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(elliptic_correction)')

   if (self%elliptic_correction) then
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='bc_correction', &
      val=self%bc_correction, error=error)
      if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(bc_correction)')
   endif

   if(self%problem_type == PLASMA_TYPE_PROBLEM) then
	   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='plasma_density', &
      val=self%plasma_density, error=error)
      if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(plasma_density)')

      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='neutral_fraction', &
      val=self%neutral_fraction, error=error)
      if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(neutral_fraction)')
   endif
   endsubroutine load_from_file

   subroutine particle_cartesian_grid_index(self, field, grid, q_pic)
   !< Compute the grid index corresponding to a particle position. Good for cartesian grids only.
   class(prism_pic_object), intent(inout) :: self               !< External fields.
   type(field_object),      intent(in)    :: field              !< The field.
   type(grid_object),       intent(in)    :: grid               !< Grid (sibling realm component, threaded in)
   real(R8P),               intent(in)    :: q_pic(1:,1:)       !< PIC variables.
   real(R8P)                              :: n, b               !< Counters
   integer(I4P)                           :: i_p, j_p, k_p, b_p !< Particle grid indices

   associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj,                                   &
            nk=>grid%nk, ngc=>grid%ngc, dx => field%dxyz(1,:), dy => field%dxyz(2,:), dz => field%dxyz(3,:), &
            np => self%particle_number, domain_emin => grid%domain_emin, domain_emax => grid%domain_emax,    &
            emin => field%emin, emax => field%emax, neighbour_list => self%neighbour_list)

   !Di sicuro va considerata una parte relativa alle particelle che escono dal dominio
   !Rivedi con Stefano, molto dipende se quei min max contano pure le gc. In tal caso a emin devi sommare ngc*dx o dx o dz (non dovrebbero contare)
   do n = 1, np
      do b = 1, blocks_number
         i_p = ceiling((q_pic(1,n) - emin(1,b)) / dx(b))
         j_p = ceiling((q_pic(2,n) - emin(2,b)) / dy(b))
         k_p = ceiling((q_pic(3,n) - emin(3,b)) / dz(b))
         b_p = b 
         if (i_p >= 1_I4P .and. i_p <= ni .and. &
             j_p >= 1_I4P .and. j_p <= nj .and. &
             k_p >= 1_I4P .and. k_p <= nk) then
            neighbour_list(1,n) = b_p
            neighbour_list(2,n) = i_p
            neighbour_list(3,n) = j_p
            neighbour_list(4,n) = k_p
            exit
         endif
      enddo
   enddo
   endassociate
   endsubroutine particle_cartesian_grid_index

   subroutine NGP_charge_weighting(self, field, grid, q, q_pic, nv)
   !!< Nearest Grid Point weighting of particle quantities to the grid.
   class(prism_pic_object), intent(inout) :: self                 !< External fields.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),       intent(in)    :: grid                 !< Grid (sibling realm component, threaded in).
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:,&
                                                  1-grid%ngc:,&
                                                  1-grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Number of variables.
   integer(I4P)                           :: n, i, j, k ,b        !< Particle counter
   integer(I4P)                           :: i_p, j_p, k_p, b_p   !< Particle grid indices
   real(R8P)                              :: wx, wy, wz           !< Weighting factors
   real(R8P)                              :: dx, dy, dz           !< Cell dimensions
   real(R8P)                              :: cell_coord(3)        !< Cell coordinates

   !Per iniziare, azzero tutte le cariche altrimenti vado a sommare le cariche del tempo precedente
   q(nv,:,:,:,:) = 0.0_R8P

   do n = 1, self%particle_number
      ! Get particle grid indices
      b_p = self%neighbour_list(1,n)
      i_p = self%neighbour_list(2,n)
      j_p = self%neighbour_list(3,n)
      k_p = self%neighbour_list(4,n)

      dx = field%dxyz(1,b_p)
      dy = field%dxyz(2,b_p)
      dz = field%dxyz(3,b_p)

      q(nv, i_p, j_p, k_p, b_p) = q(nv, i_p, j_p, k_p, b_p) + q_pic(7,n)/(dx*dy*dz)
   enddo
   endsubroutine NGP_charge_weighting

   subroutine CIC_charge_weighting(self, field, grid, q, q_PIC, nv)
   !< Cloud-in-Cell weighting of particle quantities to the grid.
   class(prism_pic_object), intent(inout) :: self                 !< External fields.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),       intent(in)    :: grid                 !< Grid (sibling realm component, threaded in).
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:, &
                                                  1-grid%ngc:, &
                                                  1-grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_PIC(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Number of variables.
   integer(I4P)                           :: n, i, j, k ,b        !< Particle counter
   integer(I4P)                           :: i_p, j_p, k_p, b_p   !< Particle grid indices
   real(R8P)                              :: dx, dy, dz           !< Grid spacing
   real(R8P)                              :: wx, wy, wz           !< Weighting factors
   real(R8P)                              :: cell_coord(3)        !< Cell coordinates

   associate(x_cell=>field%x_cell, y_cell=>field%y_cell, z_cell=>field%z_cell)

   q(nv,:,:,:,:) = 0.0_R8P

   do n = 1, self%particle_number
      ! Get particle grid indices
      b_p = self%neighbour_list(1,n)
      i_p = self%neighbour_list(2,n)
      j_p = self%neighbour_list(3,n)
      k_p = self%neighbour_list(4,n)

      dx = field%dxyz(1,b_p)
      dy = field%dxyz(2,b_p)
      dz = field%dxyz(3,b_p)

      do i = i_p-1, i_p+1
         do j = j_p-1, j_p+1
            do k = k_p-1, k_p+1
               cell_coord = [x_cell(i,b_p), y_cell(j,b_p), z_cell(k,b_p)]
               if (abs((q_pic(1,n) - cell_coord(1))/dx) <= 1.0_R8P) then
                  Wx = 1.0_R8P - abs((q_pic(1,n) - cell_coord(1))/dx)
               else
                  Wx = 0.0_R8P
               end if
               if (abs((q_pic(2,n) - cell_coord(2))/dy) <= 1.0_R8P) then
                  Wy = 1.0_R8P - abs((q_pic(2,n) - cell_coord(2))/dy)
               else
                  Wy = 0.0_R8P
               end if
               if (abs((q_pic(3,n) - cell_coord(3))/dz) <= 1.0_R8P) then
                  Wz = 1.0_R8P - abs((q_pic(3,n) - cell_coord(3))/dz)
               else
                  Wz = 0.0_R8P
               end if
               q(nv, i, j, k, b_p) = q(nv, i, j, k, b_p) + q_pic(7,n)/(dx*dy*dz) * Wx * Wy * Wz

               !Ok, ma va normalizzata e la carica nel vettore di stato va necessariamente azzerata a monte di ogni assegnazione
               !se scritta in questo modo
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine CIC_charge_weighting

   subroutine TSC_charge_weighting(self, field, grid, q, q_pic, nv)
   !!< Triangular Shaped Cloud weighting of particle quantities to the grid.
   class(prism_pic_object), intent(inout) :: self                 !< External fields.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),       intent(in)    :: grid                 !< Grid (sibling realm component, threaded in).
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:, &
                                                  1-grid%ngc:, &
                                                  1-grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Number of variables.
   real(R8P)                              :: n, i, j, k ,b        !< Particle counter
   real(R8P)                              :: i_p, j_p, k_p, b_p   !< Particle grid indices
   real(R8P)                              :: dx, dy, dz           !< Grid spacing
   real(R8P)                              :: wx, wy, wz           !< Weighting factors
   real(R8P)                              :: cell_coord(3)        !< Cell coordinates

   associate(x_cell=>field%x_cell, y_cell=>field%y_cell, z_cell=>field%z_cell)

   !Per iniziare, azzero tutte le cariche altrimenti vado a sommare le cariche del tempo precedente
   q(nv,:,:,:,:) = 0.0_R8P

   do n = 1, self%particle_number
      ! Get particle grid indices
      b_p = self%neighbour_list(1,n)
      i_p = self%neighbour_list(2,n)
      j_p = self%neighbour_list(3,n)
      k_p = self%neighbour_list(4,n)

      ! Qua va capito come gestire la questione dei blocchi multipli
      dx = field%dxyz(1,b_p)
      dy = field%dxyz(2,b_p)
      dz = field%dxyz(3,b_p)

      !Qua ci va sicuramente un if per le celle di confine, altrimenti darà errore quando arrivo alla frontiera

      do i = i_p-1, i_p+1
         do j = j_p-1, j_p+1
            do k = k_p-1, k_p+1
               cell_coord = [x_cell(i,b_p), y_cell(j,b_p), z_cell(k,b_p)]
               if (abs((q_pic(1,n) - cell_coord(1))/dx) <= 0.5_R8P) then
                  Wx = 0.75_R8P - ((q_pic(1,n) - cell_coord(1))/dx)**2
               elseif (abs((q_pic(1,n) - cell_coord(1))/dx) <= 1.5_R8P .and. abs((q_pic(1,n) - cell_coord(1))/dx) > 0.5_R8P) then
                  Wx = 0.5_R8P * (1.5_R8P - abs((q_pic(1,n) - cell_coord(1))/dx))**2
               else
                  Wx = 0.0_R8P
               end if
               if (abs((q_pic(2,n) - cell_coord(2))/dy) <= 0.5_R8P) then
                  Wy = 0.75_R8P - ((q_pic(2,n) - cell_coord(2))/dy)**2
               elseif (abs((q_pic(2,n) - cell_coord(2))/dy) <= 1.5_R8P .and. abs((q_pic(2,n) - cell_coord(2))/dy) > 0.5_R8P) then
                  Wy = 0.5_R8P * (1.5_R8P - abs((q_pic(2,n) - cell_coord(2))/dy))**2
               else
                  Wy = 0.0_R8P
               end if
               q(nv, i, j, k, b_p) = q(nv, i, j, k, b_p) + q_pic(7,n)/(dx*dy*dz) * Wx * Wy * Wz

               !Ok, ma va normalizzata e la carica nel vettore di stato va necessariamente azzerata a monte di ogni assegnazione
               !se scritta in questo modo
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine TSC_charge_weighting

   subroutine Gaussian_charge_weighting(self, field, grid, q, q_PIC, nv)
   !< Gaussian weighting of particle charge density to the grid.
   !<
   !< The Gaussian standard deviations are set equal to the local grid spacings:
   !<
   !<    sigma_x = dx
   !<    sigma_y = dy
   !<    sigma_z = dz
   !<
   !< The Gaussian distribution is truncated at +/- 4 sigma along each direction.
   !< The discrete weights are normalized particle by particle in order to preserve
   !< the total deposited charge.
   class(prism_pic_object), intent(inout) :: self                    !< PIC object.
   type(field_object),      intent(inout) :: field                   !< The field.
   type(grid_object),       intent(in)    :: grid                    !< Grid object.
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:, &
                                                   1-grid%ngc:, &
                                                   1-grid%ngc:,1:)   !< Field variables.
   real(R8P),               intent(in)    :: q_PIC(1:,1:)            !< PIC variables.
   integer(I4P),            intent(in)    :: nv                      !< Charge-density variable index.
   integer(I4P)                           :: n                       !< Particle counter.
   integer(I4P)                           :: i,j,k                   !< Grid counters.
   integer(I4P)                           :: i_p,j_p,k_p,b_p         !< Particle grid indices.
   integer(I4P)                           :: i_min,i_max             !< Weighting bounds.
   integer(I4P)                           :: j_min,j_max             !< Weighting bounds.
   integer(I4P)                           :: k_min,k_max             !< Weighting bounds.
   integer(I4P)                           :: ni_sigma                !< Support radius along x.
   integer(I4P)                           :: nj_sigma                !< Support radius along y.
   integer(I4P)                           :: nk_sigma                !< Support radius along z.
   real(R8P)                              :: dx,dy,dz                !< Grid spacings.
   real(R8P)                              :: sigma_x                 !< Gaussian standard deviation along x.
   real(R8P)                              :: sigma_y                 !< Gaussian standard deviation along y.
   real(R8P)                              :: sigma_z                 !< Gaussian standard deviation along z.
   real(R8P)                              :: rx,ry,rz                !< Normalized particle-cell distances.
   real(R8P)                              :: wx,wy,wz                !< One-dimensional Gaussian weights.
   real(R8P)                              :: weight                  !< Three-dimensional Gaussian weight.
   real(R8P)                              :: weight_sum              !< Discrete normalization factor.
   real(R8P)                              :: inverse_cell_volume     !< Inverse cell volume.

   associate(x_cell=>field%x_cell, y_cell=>field%y_cell, z_cell=>field%z_cell, sigma=>self%sigma, cutoff_sigma=>self%cutoff_sigma)

   ! Reset the charge density before depositing the current particle distribution.
   q(nv,:,:,:,:) = 0._R8P

   do n=1, self%particle_number
      ! Get particle block and grid indices.
      b_p = self%neighbour_list(1,n)
      i_p = self%neighbour_list(2,n)
      j_p = self%neighbour_list(3,n)
      k_p = self%neighbour_list(4,n)

      dx = field%dxyz(1,b_p)
      dy = field%dxyz(2,b_p)
      dz = field%dxyz(3,b_p)

      sigma_x = sigma
      sigma_y = sigma
      sigma_z = sigma

      inverse_cell_volume = 1._R8P / (dx*dy*dz)

      ! Number of cells required to cover the Gaussian cutoff.
      ni_sigma = ceiling(cutoff_sigma*sigma_x/dx, kind=I4P)
      nj_sigma = ceiling(cutoff_sigma*sigma_y/dy, kind=I4P)
      nk_sigma = ceiling(cutoff_sigma*sigma_z/dz, kind=I4P)

      ! Restrict the support to the locally available grid, including ghost cells.
      i_min = max(i_p-ni_sigma, lbound(q,dim=2))
      i_max = min(i_p+ni_sigma, ubound(q,dim=2))

      j_min = max(j_p-nj_sigma, lbound(q,dim=3))
      j_max = min(j_p+nj_sigma, ubound(q,dim=3))

      k_min = max(k_p-nk_sigma, lbound(q,dim=4))
      k_max = min(k_p+nk_sigma, ubound(q,dim=4))

      ! Compute the discrete normalization factor over the effective support.
      weight_sum = 0._R8P

      do k=k_min, k_max
         rz = (q_PIC(3,n)-z_cell(k,b_p))/sigma_z
         wz = exp(-0.5_R8P*rz*rz)
         do j=j_min, j_max
            ry = (q_PIC(2,n)-y_cell(j,b_p))/sigma_y
            wy = exp(-0.5_R8P*ry*ry)
            do i=i_min, i_max
               rx = (q_PIC(1,n)-x_cell(i,b_p))/sigma_x
               wx = exp(-0.5_R8P*rx*rx)
               weight_sum = weight_sum + wx*wy*wz
            enddo
         enddo
      enddo

      ! Deposit the particle charge density.
      if (weight_sum > tiny(1._R8P)) then
         do k=k_min, k_max
            rz = (q_PIC(3,n)-z_cell(k,b_p))/sigma_z
            wz = exp(-0.5_R8P*rz*rz)
            do j=j_min, j_max
               ry = (q_PIC(2,n)-y_cell(j,b_p))/sigma_y
               wy = exp(-0.5_R8P*ry*ry)
               do i=i_min, i_max
                  rx = (q_PIC(1,n)-x_cell(i,b_p))/sigma_x
                  wx = exp(-0.5_R8P*rx*rx)
                  weight = wx*wy*wz/weight_sum
                  q(nv,i,j,k,b_p) = q(nv,i,j,k,b_p) + q_PIC(7,n)*inverse_cell_volume*weight
               enddo
            enddo
         enddo
      endif
   enddo
   endassociate
   endsubroutine Gaussian_charge_weighting

   subroutine NGP_current_weighting(self, field, grid, q, q_pic, nv)
   !!< Nearest Grid Point weighting of particle quantities to the grid.
   class(prism_pic_object), intent(inout) :: self                 !< External fields.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),       intent(in)    :: grid                 !< Grid (sibling realm component, threaded in).
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:,&
                                                  1-grid%ngc:,&
                                                  1-grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Number of variables.
   real(R8P)                              :: n, i, j, k ,b        !< Particle counter
   real(R8P)                              :: i_p, j_p, k_p, b_p   !< Particle grid indices
   real(R8P)                              :: dx, dy, dz           !< Grid spacing
   real(R8P)                              :: wx, wy, wz           !< Weighting factors
   real(R8P)                              :: cell_coord(3)        !< Cell coordinates

   !Per iniziare, azzero tutte le correnti altrimenti vado a sommare le cariche del tempo precedente
   q((nv-3):(nv-1),:,:,:,:) = 0.0_R8P

   do n = 1, self%particle_number
      ! Get particle grid indices
      b_p = self%neighbour_list(1,n)
      i_p = self%neighbour_list(2,n)
      j_p = self%neighbour_list(3,n)
      k_p = self%neighbour_list(4,n)

      ! Qua va capito come gestire la questione dei blocchi multipli
      dx = field%dxyz(1,b_p)
      dy = field%dxyz(2,b_p)
      dz = field%dxyz(3,b_p)

      !Qua ci va sicuramente un if per le celle di confine, altrimenti darà errore quando arrivo alla frontiera
      q(nv-3, i_p, j_p, k_p, b_p) = q(nv-3, i_p, j_p, k_p, b_p) + q_pic(7,n)/(dx*dy*dz)*q_pic(4,n)
      q(nv-2, i_p, j_p, k_p, b_p) = q(nv-2, i_p, j_p, k_p, b_p) + q_pic(7,n)/(dx*dy*dz)*q_pic(5,n)
      q(nv-1, i_p, j_p, k_p, b_p) = q(nv-1, i_p, j_p, k_p, b_p) + q_pic(7,n)/(dx*dy*dz)*q_pic(6,n)
      !Ok, ma va normalizzata e la carica nel vettore di stato va necessariamente azzerata a monte di ogni assegnazione
      !se scritta in questo modo
   enddo
   endsubroutine NGP_current_weighting

   subroutine CIC_current_weighting(self, field, grid, q, q_pic, nv)
   !< Cloud-in-Cell weighting of particle quantities to the grid.
   class(prism_pic_object), intent(inout) :: self                 !< External fields.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),       intent(in)    :: grid                 !< Grid (sibling realm component, threaded in).
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:, &
                                                  1-grid%ngc:, &
                                                  1-grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Number of variables.
   real(R8P)                              :: n, i, j, k ,b        !< Particle counter
   real(R8P)                              :: i_p, j_p, k_p, b_p   !< Particle grid indices
   real(R8P)                              :: dx, dy, dz           !< Grid spacing
   real(R8P)                              :: wx, wy, wz           !< Weighting factors
   real(R8P)                              :: cell_coord(3)        !< Cell coordinates

   associate(x_cell=>field%x_cell, y_cell=>field%y_cell, z_cell=>field%z_cell)

   !Per iniziare, azzero tutte le correnti altrimenti vado a sommare le cariche del tempo precedente
   q((nv-3):(nv-1),:,:,:,:) = 0.0_R8P

   do n = 1, self%particle_number
      ! Get particle grid indices
      b_p = self%neighbour_list(1,n)
      i_p = self%neighbour_list(2,n)
      j_p = self%neighbour_list(3,n)
      k_p = self%neighbour_list(4,n)

      ! Qua va capito come gestire la questione dei blocchi multipli
      dx = field%dxyz(1,b_p)
      dy = field%dxyz(2,b_p)
      dz = field%dxyz(3,b_p)

      !Qua ci va sicuramente un if per le celle di confine, altrimenti darà errore quando arrivo alla frontiera

      do i = i_p-1, i_p+1
         do j = j_p-1, j_p+1
            do k = k_p-1, k_p+1
               cell_coord = [x_cell(i,b_p), y_cell(j,b_p), z_cell(k,b_p)]
               if (abs((q_pic(1,n) - cell_coord(1))/dx) <= 1.0_R8P) then
                  Wx = 1.0_R8P - abs((q_pic(1,n) - cell_coord(1))/dx)
               else
                  Wx = 0.0_R8P
               end if
               if (abs((q_pic(2,n) - cell_coord(2))/dy) <= 1.0_R8P) then
                  Wy = 1.0_R8P - abs((q_pic(2,n) - cell_coord(2))/dy)
               else
                  Wy = 0.0_R8P
               end if
               if (abs((q_pic(3,n) - cell_coord(3))/dz) <= 1.0_R8P) then
                  Wz = 1.0_R8P - abs((q_pic(3,n) - cell_coord(3))/dz)
               else
                  Wz = 0.0_R8P
               end if
               q(nv-3, i_p, j_p, k_p, b_p) = q(nv-3, i_p, j_p, k_p, b_p) + q_pic(7,n)/(dx*dy*dz)*q_pic(4,n)* Wx * Wy * Wz
               q(nv-2, i_p, j_p, k_p, b_p) = q(nv-2, i_p, j_p, k_p, b_p) + q_pic(7,n)/(dx*dy*dz)*q_pic(5,n)* Wx * Wy * Wz
               q(nv-1, i_p, j_p, k_p, b_p) = q(nv-1, i_p, j_p, k_p, b_p) + q_pic(7,n)/(dx*dy*dz)*q_pic(6,n)* Wx * Wy * Wz

               !Ok, ma va normalizzata e la carica nel vettore di stato va necessariamente azzerata a monte di ogni assegnazione
               !se scritta in questo modo
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine CIC_current_weighting

   subroutine TSC_current_weighting(self, field, grid, q, q_pic, nv)
   !!< Triangular Shaped Cloud weighting of particle quantities to the grid.
   class(prism_pic_object), intent(inout) :: self                 !< External fields.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),       intent(in)    :: grid                 !< Grid (sibling realm component, threaded in).
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:, &
                                                  1-grid%ngc:, &
                                                  1-grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Number of variables.
   real(R8P)                              :: n, i, j, k ,b        !< Particle counter
   real(R8P)                              :: i_p, j_p, k_p, b_p   !< Particle grid indices
   real(R8P)                              :: dx, dy, dz           !< Grid spacing
   real(R8P)                              :: wx, wy, wz           !< Weighting factors
   real(R8P)                              :: cell_coord(3)        !< Cell coordinates

   associate(x_cell=>field%x_cell, y_cell=>field%y_cell, z_cell=>field%z_cell)

   !Per iniziare, azzero tutte le correnti altrimenti vado a sommare le cariche del tempo precedente
   q((nv-3):(nv-1),:,:,:,:) = 0.0_R8P

   do n = 1, self%particle_number
      ! Get particle grid indices
      b_p = self%neighbour_list(1,n)
      i_p = self%neighbour_list(2,n)
      j_p = self%neighbour_list(3,n)
      k_p = self%neighbour_list(4,n)

      ! Qua va capito come gestire la questione dei blocchi multipli
      dx = field%dxyz(1,b_p)
      dy = field%dxyz(2,b_p)
      dz = field%dxyz(3,b_p)

      !Qua ci va sicuramente un if per le celle di confine, altrimenti darà errore quando arrivo alla frontiera

      do i = i_p-1, i_p+1
         do j = j_p-1, j_p+1
            do k = k_p-1, k_p+1
               cell_coord = [x_cell(i,b_p), y_cell(j,b_p), z_cell(k,b_p)]
               if (abs((q_pic(1,n) - cell_coord(1))/dx) <= 0.5_R8P) then
                  Wx = 0.75_R8P - ((q_pic(1,n) - cell_coord(1))/dx)**2
               elseif (abs((q_PIC(1,n) - cell_coord(1))/dx) <= 1.5_R8P .and. abs((q_pic(1,n) - cell_coord(1))/dx) > 0.5_R8P) then
                  Wx = 0.5_R8P * (1.5_R8P - abs((q_pic(1,n) - cell_coord(1))/dx))**2
               else
                  Wx = 0.0_R8P
               end if
               if (abs((q_pic(2,n) - cell_coord(2))/dy) <= 0.5_R8P) then
                  Wy = 0.75_R8P - ((q_pic(2,n) - cell_coord(2))/dy)**2
               elseif (abs((q_pic(2,n) - cell_coord(2))/dy) <= 1.5_R8P .and. abs((q_pic(2,n) - cell_coord(2))/dy) > 0.5_R8P) then
                  Wy = 0.5_R8P * (1.5_R8P - abs((q_pic(2,n) - cell_coord(2))/dy))**2
               else
                  Wy = 0.0_R8P
               end if
               q(nv-3, i_p, j_p, k_p, b_p) = q(nv-3, i_p, j_p, k_p, b_p) + q_pic(7,n)/(dx*dy*dz)*q_pic(4,n)* Wx * Wy * Wz
               q(nv-2, i_p, j_p, k_p, b_p) = q(nv-2, i_p, j_p, k_p, b_p) + q_pic(7,n)/(dx*dy*dz)*q_pic(5,n)* Wx * Wy * Wz
               q(nv-1, i_p, j_p, k_p, b_p) = q(nv-1, i_p, j_p, k_p, b_p) + q_pic(7,n)/(dx*dy*dz)*q_pic(6,n)* Wx * Wy * Wz

               !Ok, ma va normalizzata e la carica nel vettore di stato va necessariamente azzerata a monte di ogni assegnazione
               !se scritta in questo modo
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine TSC_current_weighting

   subroutine zeroD_field_weighting(self, field, grid, pic_fields, q, q_pic, nv)
   !< 0D spatial interpolation for the fields
   class(prism_pic_object), intent(inout) :: self                                             !< External fields.
   type(field_object),      intent(inout) :: field                                            !< The field.
   type(grid_object),       intent(in)    :: grid                                             !< Grid (sibling realm component, threaded in).
   real(R8P),               intent(inout) :: pic_fields(1:,1:)                                !< Fields value at particle locations
   real(R8P),               intent(in)    :: q(1:,1-grid%ngc:, &
                                                  1-grid%ngc:, &
                                                  1-grid%ngc:,1:)                             !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)                                     !< PIC variables.
   integer(I4P),            intent(in)    :: nv                                               !< Number of variables.
   real(R8P)                              :: n, i, j, k ,b                                    !< Particle counter
   real(R8P)                              :: i_p, j_p, k_p, block_p                           !< Particle grid indices
   real(R8P)                              :: v_p(3), D_p(3), B_p(3), F_p(3), F_l(3), m_p, q_p !< Particle scalars

   do n = 1, self%particle_number
      ! Get particle grid indices
      block_p = self%neighbour_list(1,n)
      i_p     = self%neighbour_list(2,n)
      j_p     = self%neighbour_list(3,n)
      k_p     = self%neighbour_list(4,n)

      D_p(1) = q(1, i_p, j_p, k_p, block_p)
      D_p(2) = q(2, i_p, j_p, k_p, block_p)
      D_p(3) = q(3, i_p, j_p, k_p, block_p)
      B_p(1) = q(4, i_p, j_p, k_p, block_p)
      B_p(2) = q(5, i_p, j_p, k_p, block_p)
      B_p(3) = q(6, i_p, j_p, k_p, block_p)

      !v_p = [q_pic(n,4), q_pic(n,5), q_pic(n,6)]
      !q_p = q_pic(n,7)
      !m_p = q_pic(n,8)
!
      !F_l = crossproduct(a=v_p, b=B_p)
      !F_p(1) = q_p*(D_p(1) + F_l(1))
      !F_p(2) = q_p*(D_p(2) + F_l(2))
      !F_p(3) = q_p*(D_p(3) + F_l(3))

      pic_fields(1,n) = D_p(1)
      pic_fields(2,n) = D_p(2)
      pic_fields(3,n) = D_p(3)
      pic_fields(4,n) = B_p(1)
      pic_fields(5,n) = B_p(2)
      pic_fields(6,n) = B_p(3)

   enddo
   endsubroutine zeroD_field_weighting

   subroutine oneD_field_weighting(self, field, grid, pic_fields, q, q_pic, nv)
   !< 1D spatial interpolation for the fields
   class(prism_pic_object), intent(inout) :: self                                             !< External fields.
   type(field_object),      intent(inout) :: field                                            !< The field.
   type(grid_object),       intent(in)    :: grid                                             !< Grid (sibling realm component, threaded in).
   real(R8P),               intent(inout) :: pic_fields(1:,1:)                                !< Fields value at particle locations
   real(R8P),               intent(in)    :: q(1:,1-grid%ngc:,&
                                                  1-grid%ngc:,&
                                                  1-grid%ngc:,1:)                             !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)                                     !< PIC variables.
   integer(I4P),            intent(in)    :: nv                                               !< Number of variables.
   real(R8P)                              :: n, i, j, k ,b                                    !< Particle counter
   real(R8P)                              :: i_p, j_p, k_p, block_p                           !< Particle grid indices
   real(R8P)                              :: x_p, y_p, z_p                                    !< Particle position scalar
   real(R8P)                              :: v_p(3), D_p(3), B_p(3), F_p(3), F_l(3), m_p, q_p !< Particle scalars
   real(R8P)                              :: dx, dy, dz                                       !< Grid spacing
   real(R8P)                              :: dDx_dx, dDx_dy, dDx_dz                           !<
   real(R8P)                              :: dDy_dx, dDy_dy, dDy_dz                           !<
   real(R8P)                              :: dDz_dx, dDz_dy, dDz_dz                           !<
   real(R8P)                              :: dBx_dx, dBx_dy, dBx_dz                           !<
   real(R8P)                              :: dBy_dx, dBy_dy, dBy_dz                           !<
   real(R8P)                              :: dBz_dx, dBz_dy, dBz_dz                           !<

   associate(x_cell=>field%x_cell, y_cell=>field%y_cell, z_cell=>field%z_cell)
   do n = 1, self%particle_number

      ! Get particle grid indices
      block_p = self%neighbour_list(1,n)
      i_p     = self%neighbour_list(2,n)
      j_p     = self%neighbour_list(3,n)
      k_p     = self%neighbour_list(4,n)

      ! Qua va capito come gestire la questione dei blocchi multipli

      !Interpolazione lineare dei campi nella posizione delle particelle
      !x
      if (x_cell(i_p,block_p) >= q_pic(1,n)) then !La particella è nella metà sinistra della cella
         dDx_dx = lininterp(x2=x_cell(i_p-1,block_p), y2=q(1,i_p-1,j_p,k_p,block_p), &
                            x1=x_cell(i_p,block_p), y1=q(1,i_p,j_p,k_p,block_p), xp=q_pic(1,n))
         dDy_dx = lininterp(x2=x_cell(i_p-1,block_p), y2=q(2,i_p-1,j_p,k_p,block_p), &
                            x1=x_cell(i_p,block_p), y1=q(2,i_p,j_p,k_p,block_p), xp=q_pic(1,n))
         dDz_dx = lininterp(x2=x_cell(i_p-1,block_p), y2=q(3,i_p-1,j_p,k_p,block_p), &
                            x1=x_cell(i_p,block_p), y1=q(3,i_p,j_p,k_p,block_p), xp=q_pic(1,n))
         dBx_dx = lininterp(x2=x_cell(i_p-1,block_p), y2=q(4,i_p-1,j_p,k_p,block_p), &
                            x1=x_cell(i_p,block_p), y1=q(4,i_p,j_p,k_p,block_p), xp=q_pic(1,n))
         dBy_dx = lininterp(x2=x_cell(i_p-1,block_p), y2=q(5,i_p-1,j_p,k_p,block_p), &
                            x1=x_cell(i_p,block_p), y1=q(5,i_p,j_p,k_p,block_p), xp=q_pic(1,n))
         dBz_dx = lininterp(x2=x_cell(i_p-1,block_p), y2=q(6,i_p-1,j_p,k_p,block_p), &
                            x1=x_cell(i_p,block_p), y1=q(6,i_p,j_p,k_p,block_p), xp=q_pic(1,n))
      else
         dDx_dx = lininterp(x2=x_cell(i_p+1,block_p), y2=q(1,i_p+1,j_p,k_p,block_p), &
                            x1=x_cell(i_p,block_p), y1=q(1,i_p,j_p,k_p,block_p), xp=q_pic(1,n))
         dDy_dx = lininterp(x2=x_cell(i_p+1,block_p), y2=q(2,i_p+1,j_p,k_p,block_p), &
                            x1=x_cell(i_p,block_p), y1=q(2,i_p,j_p,k_p,block_p), xp=q_pic(1,n))
         dDz_dx = lininterp(x2=x_cell(i_p+1,block_p), y2=q(3,i_p+1,j_p,k_p,block_p), &
                            x1=x_cell(i_p,block_p), y1=q(3,i_p,j_p,k_p,block_p), xp=q_pic(1,n))
         dBx_dx = lininterp(x2=x_cell(i_p+1,block_p), y2=q(4,i_p+1,j_p,k_p,block_p), &
                            x1=x_cell(i_p,block_p), y1=q(4,i_p,j_p,k_p,block_p), xp=q_pic(1,n))
         dBy_dx = lininterp(x2=x_cell(i_p+1,block_p), y2=q(5,i_p+1,j_p,k_p,block_p), &
                            x1=x_cell(i_p,block_p), y1=q(5,i_p,j_p,k_p,block_p), xp=q_pic(1,n))
         dBz_dx = lininterp(x2=x_cell(i_p+1,block_p), y2=q(6,i_p+1,j_p,k_p,block_p), &
                            x1=x_cell(i_p,block_p), y1=q(6,i_p,j_p,k_p,block_p), xp=q_pic(1,n))
      endif
      !y
      if (y_cell(j_p,block_p) >= q_pic(2,n)) then !La particella è nella metà sinistra della cella
         dDx_dy = lininterp(x2=y_cell(j_p-1,block_p), y2=q(1,i_p,j_p-1,k_p,block_p), &
                            x1=y_cell(j_p,block_p), y1=q(1,i_p,j_p,k_p,block_p), xp=q_pic(2,n))
         dDy_dy = lininterp(x2=y_cell(j_p-1,block_p), y2=q(2,i_p,j_p-1,k_p,block_p), &
                            x1=y_cell(j_p,block_p), y1=q(2,i_p,j_p,k_p,block_p), xp=q_pic(2,n))
         dDz_dy = lininterp(x2=y_cell(j_p-1,block_p), y2=q(3,i_p,j_p-1,k_p,block_p), &
                            x1=y_cell(j_p,block_p), y1=q(3,i_p,j_p,k_p,block_p), xp=q_pic(2,n))
         dBx_dy = lininterp(x2=y_cell(j_p-1,block_p), y2=q(4,i_p,j_p-1,k_p,block_p), &
                            x1=y_cell(j_p,block_p), y1=q(4,i_p,j_p,k_p,block_p), xp=q_pic(2,n))
         dBy_dy = lininterp(x2=y_cell(j_p-1,block_p), y2=q(5,i_p,j_p-1,k_p,block_p), &
                            x1=y_cell(j_p,block_p), y1=q(5,i_p,j_p,k_p,block_p), xp=q_pic(2,n))
         dBz_dy = lininterp(x2=y_cell(j_p-1,block_p), y2=q(6,i_p,j_p-1,k_p,block_p), &
                            x1=y_cell(j_p,block_p), y1=q(6,i_p,j_p,k_p,block_p), xp=q_pic(2,n))
      else
         dDx_dy = lininterp(x2=y_cell(j_p+1,block_p), y2=q(1,i_p,j_p+1,k_p,block_p), &
                            x1=y_cell(j_p,block_p), y1=q(1,i_p,j_p,k_p,block_p), xp=q_pic(2,n))
         dDy_dy = lininterp(x2=y_cell(j_p+1,block_p), y2=q(2,i_p,j_p+1,k_p,block_p), &
                            x1=y_cell(j_p,block_p), y1=q(2,i_p,j_p,k_p,block_p), xp=q_pic(2,n))
         dDz_dy = lininterp(x2=y_cell(j_p+1,block_p), y2=q(3,i_p,j_p+1,k_p,block_p), &
                            x1=y_cell(j_p,block_p), y1=q(3,i_p,j_p,k_p,block_p), xp=q_pic(2,n))
         dBx_dy = lininterp(x2=y_cell(j_p+1,block_p), y2=q(4,i_p,j_p+1,k_p,block_p), &
                            x1=y_cell(j_p,block_p), y1=q(4,i_p,j_p,k_p,block_p), xp=q_pic(2,n))
         dBy_dy = lininterp(x2=y_cell(j_p+1,block_p), y2=q(5,i_p,j_p+1,k_p,block_p), &
                            x1=y_cell(j_p,block_p), y1=q(5,i_p,j_p,k_p,block_p), xp=q_pic(2,n))
         dBz_dy = lininterp(x2=y_cell(j_p+1,block_p), y2=q(6,i_p,j_p+1,k_p,block_p), &
                            x1=y_cell(j_p,block_p), y1=q(6,i_p,j_p,k_p,block_p), xp=q_pic(2,n))
      endif
      !z
      if (z_cell(k_p,block_p) >= q_pic(3,n)) then !La particella è nella metà sinistra della cella
         dDx_dz = lininterp(x2=z_cell(k_p-1,block_p), y2=q(1,i_p,j_p,k_p-1,block_p), &
                            x1=z_cell(k_p,block_p), y1=q(1,i_p,j_p,k_p,block_p), xp=q_pic(3,n))
         dDy_dz = lininterp(x2=z_cell(k_p-1,block_p), y2=q(2,i_p,j_p,k_p-1,block_p), &
                            x1=z_cell(k_p,block_p), y1=q(2,i_p,j_p,k_p,block_p), xp=q_pic(3,n))
         dDz_dz = lininterp(x2=z_cell(k_p-1,block_p), y2=q(3,i_p,j_p,k_p-1,block_p), &
                            x1=z_cell(k_p,block_p), y1=q(3,i_p,j_p,k_p,block_p), xp=q_pic(3,n))
         dBx_dz = lininterp(x2=z_cell(k_p-1,block_p), y2=q(4,i_p,j_p,k_p-1,block_p), &
                            x1=z_cell(k_p,block_p), y1=q(4,i_p,j_p,k_p,block_p), xp=q_pic(3,n))
         dBy_dz = lininterp(x2=z_cell(k_p-1,block_p), y2=q(5,i_p,j_p,k_p-1,block_p), &
                            x1=z_cell(k_p,block_p), y1=q(5,i_p,j_p,k_p,block_p), xp=q_pic(3,n))
         dBz_dz = lininterp(x2=z_cell(k_p-1,block_p), y2=q(6,i_p,j_p,k_p-1,block_p), &
                            x1=z_cell(k_p,block_p), y1=q(6,i_p,j_p,k_p,block_p), xp=q_pic(3,n))
      else
         dDx_dz = lininterp(x2=z_cell(k_p+1,block_p), y2=q(1,i_p,j_p,k_p+1,block_p), &
                            x1=z_cell(k_p,block_p), y1=q(1,i_p,j_p,k_p,block_p), xp=q_pic(3,n))
         dDy_dz = lininterp(x2=z_cell(k_p+1,block_p), y2=q(2,i_p,j_p,k_p+1,block_p), &
                            x1=z_cell(k_p,block_p), y1=q(2,i_p,j_p,k_p,block_p), xp=q_pic(3,n))
         dDz_dz = lininterp(x2=z_cell(k_p+1,block_p), y2=q(3,i_p,j_p,k_p+1,block_p), &
                            x1=z_cell(k_p,block_p), y1=q(3,i_p,j_p,k_p,block_p), xp=q_pic(3,n))
         dBx_dz = lininterp(x2=z_cell(k_p+1,block_p), y2=q(4,i_p,j_p,k_p+1,block_p), &
                            x1=z_cell(k_p,block_p), y1=q(4,i_p,j_p,k_p,block_p), xp=q_pic(3,n))
         dBy_dz = lininterp(x2=z_cell(k_p+1,block_p), y2=q(5,i_p,j_p,k_p+1,block_p), &
                            x1=z_cell(k_p,block_p), y1=q(5,i_p,j_p,k_p,block_p), xp=q_pic(3,n))
         dBz_dz = lininterp(x2=z_cell(k_p+1,block_p), y2=q(6,i_p,j_p,k_p+1,block_p), &
                            x1=z_cell(k_p,block_p), y1=q(6,i_p,j_p,k_p,block_p), xp=q_pic(3,n))
      endif

      D_p(1) = q(1,i_p,j_p,k_p,block_p) + dDx_dx + dDx_dy + dDx_dz
      D_p(2) = q(2,i_p,j_p,k_p,block_p) + dDy_dx + dDy_dy + dDy_dz
      D_p(3) = q(3,i_p,j_p,k_p,block_p) + dDz_dx + dDz_dy + dDz_dz
      B_p(1) = q(4,i_p,j_p,k_p,block_p) + dBx_dx + dBx_dy + dBx_dz
      B_p(2) = q(5,i_p,j_p,k_p,block_p) + dBy_dx + dBy_dy + dBy_dz
      B_p(3) = q(6,i_p,j_p,k_p,block_p) + dBz_dx + dBz_dy + dBz_dz

      !v_p = [q_pic(n,4), q_pic(n,5), q_pic(n,6)]
      !q_p = q_PIC(n,7)
      !m_p = q_PIC(n,8)
      !F_l = crossproduct(a=v_p, b=B_p)
      !F_p(1) = q_p*(D_p(1) + F_l(1))
      !F_p(2) = q_p*(D_p(2) + F_l(2))
      !F_p(3) = q_p*(D_p(3) + F_l(3))

      pic_fields(1,n) = D_p(1)
      pic_fields(2,n) = D_p(2)
      pic_fields(3,n) = D_p(3)
      pic_fields(4,n) = B_p(1)
      pic_fields(5,n) = B_p(2)
      pic_fields(6,n) = B_p(3)

   enddo
   endassociate
   endsubroutine oneD_field_weighting

   function crossproduct(a, b) result(cross)
   real(R8P), intent(in) :: a(3)     !< Left hand side.
   real(R8P), intent(in) :: b(3)     !< Left hand side.
   real(R8P)             :: cross(3) !< Cross product.

   cross(1) = (a(2) * b(3)) - (a(3) * b(2))
   cross(2) = (a(3) * b(1)) - (a(1) * b(3))
   cross(3) = (a(1) * b(2)) - (a(2) * b(1))
   endfunction crossproduct

   function lininterp(x1, y1, x2, y2, xp) result(delta)
   real(R8P), intent(in) :: x1, x2
   real(R8P), intent(in) :: y1, y2
   real(R8P), intent(in) :: xp
   real(R8P)             :: delta
   real(R8P)             :: m

   m = (y2-y1)/(x2-x1)
   delta = m*(xp-x1)
   endfunction lininterp
endmodule adam_prism_pic_object
