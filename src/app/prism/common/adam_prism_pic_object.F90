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
use :: adam_prism_physics_object, only : prism_physics_object, is_adim_model
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
public :: UNIFORM_DOMAIN
public :: UNIFORM_CILINDER
public :: UNIFORM_CELL
!public :: CIC_WEIGHTING_MODEL
!public :: NGP_WEIGHTING_MODEL
!public :: TSC_WEIGHTING_MODEL
!public :: ZEROD_FIELDS_WEIGHTING_MODEL
!public :: ONED_FIELDS_WEIGHTING_MODEL
public :: NUM_SCHEME_TIME_PIC_LEAPFROG
public :: NUM_SCHEME_TIME_PIC_RUNGE_KUTTA
public :: COHERENT_INITIALIZATION
public :: effective_gaussian_weight
!public :: CIC_charge_weighting
!public :: NGP_charge_weighting
!public :: TSC_charge_weighting
!public :: CIC_current_weighting
!public :: NGP_current_weighting
!public :: TSC_current_weighting
!public :: zeroD_field_weighting
!public :: oneD_field_weighting

character(len=3 ), parameter :: INI_SECTION_NAME                 = 'PIC' !< INI file section name for PIC configuration.
character(len=6 ), parameter :: PLASMA_TYPE_PROBLEM              = 'plasma' !< Analyzing physical problem involving the presence of
   !< plasma
character(len=8 ), parameter :: COHERENT_INITIALIZATION          = 'coherent' !< Field initialization through elliptic solver with
   !< the coherent laplacian scheme.
character(len=15), parameter :: SINGLE_PARTICLE_TYPE_PROBLEM     = 'single_particle' !< Analyzing physical problem involving the
   !< presence of a single particle
character(len=3 ), parameter :: NGP_WEIGHTING_MODEL              = 'NGP'              !< NGP weighting model.
character(len=3 ), parameter :: CIC_WEIGHTING_MODEL              = 'CIC'              !< CIC weighting model.
character(len=3 ), parameter :: TSC_WEIGHTING_MODEL              = 'TSC'              !< TSC weighting model.
character(len=5 ), parameter :: CUBIC_WEIGHTING_MODEL            = 'cubic'            !< cubic order weighting model.
character(len=7 ), parameter :: QUARTIC_WEIGHTING_MODEL          = 'quartic'          !< quartic order weighting model.
character(len=7 ), parameter :: QUINTIC_WEIGHTING_MODEL          = 'quintic'          !< quintic order weighting model.
character(len=6 ), parameter :: SEXTIC_WEIGHTING_MODEL = 'sextic' !< Sextic B-spline weighting.
character(len=8 ), parameter :: GAUSSIAN_WEIGHTING_MODEL         = 'Gaussian'         !< Gaussian weighting model.
character(len=2 ), parameter :: ZEROD_FIELDS_WEIGHTING_MODEL     = '0D'               !< 0D field weighting.
character(len=2 ), parameter :: ONED_FIELDS_WEIGHTING_MODEL      = '1D'               !< 1D field weighting.
character(len=2 ), parameter :: TWOD_FIELDS_WEIGHTING_MODEL      = '2D'               !< 2D field weighting.
character(len=2 ), parameter :: THREED_FIELDS_WEIGHTING_MODEL    = '3D'               !< 3D field weighting.
character(len=2 ), parameter :: FOURD_FIELDS_WEIGHTING_MODEL     = '4D'               !< 4D field weighting.
character(len=2 ), parameter :: FIVED_FIELDS_WEIGHTING_MODEL     = '5D'               !< 5D field weighting.
character(len=2 ), parameter :: SIXD_FIELDS_WEIGHTING_MODEL = '6D' !< Sextic field weighting.
character(len=8 ), parameter :: NUM_SCHEME_TIME_PIC_LEAPFROG     = 'LEAPFROG' !< Leapfrog numerical scheme for time operator.
character(len=11), parameter :: NUM_SCHEME_TIME_PIC_RUNGE_KUTTA  = 'RUNGE_KUTTA' !< Runge-Kutta numerical scheme for time operator.
character(len=14), parameter :: UNIFORM_DOMAIN                   = 'Uniform_domain'   !<
character(len=16), parameter :: UNIFORM_CILINDER                 = 'Uniform_cilinder' !<
character(len=12), parameter :: UNIFORM_CELL                     = 'Uniform_cell'     !<
!character(len=32), parameter :: UNIFORM_BOX_SPACE_DISTRIBUTION   = 'Uniform_boxes_space_distribution'
!character(len=31), parameter :: UNIFORM_CELL_SPACE_DISTRIBUTION  = 'Uniform_cell_space_distribution'

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
   real(R8P)		           :: cilinder_length = 0.0_R8P   !<
   real(R8P)		           :: cilinder_radius = 0.0_R8P   !<
   real(R8P)                 :: cilinder_center(3)          !<
   character(len=1)          :: cilinder_axis               !<
   real(R8P)                 :: sigma = 0.0_R8P             !< Standard deviation for a gaussian weighting
   real(R8P)                 :: cutoff_sigma = 0._R8P       !< Gaussian cutoff.
   integer(I4P)              :: gaussian_support_cells = -1_I4P !< Gaussian support radius; negative selects automatic.
   integer(I4P)              :: particle_number  = 0_I4P    !< Total number of particles.
	integer(I4P)				  :: n_ions = 0_I4P              !< Total ions number
	integer(I4P)				  :: n_electrons = 0_I4P         !< Total electrons number
	integer(I4P)				  :: n_neutrals = 0_I4P          !< Total neutrals number
   integer(I4P), allocatable :: neighbour_list(:,:)         !< Particle grid positions array.
   character(len=99)         :: problem_type = ''           !< Type of problem analyzed
   character(len=99)         :: plasma_domain = ''          !< Domain of plasma at t0
   character(len=99)         :: initialization = COHERENT_INITIALIZATION !< Field initialization solver.
   logical                   :: elliptic_correction=.false. !< elliptic correction for the initial fields
   logical                   :: filter_deposition = .false. !< Binomial filter for deposition and gather weights.
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
   procedure, pass(self) :: cubic_charge_weighting        !< 3rd order weighting of particle quantities to the grid.
   procedure, pass(self) :: quartic_charge_weighting      !< 4th order weighting of particle quantities to the grid.
   procedure, pass(self) :: sextic_charge_weighting
   procedure, pass(self) :: sextic_current_weighting
   procedure, pass(self) :: quintic_charge_weighting      !< 5th order weighting of particle quantities to the grid.
   procedure, pass(self) :: Gaussian_charge_weighting     !< Gaussian Shaped Cloud weighting of particle quantities to the grid.
   procedure, pass(self) :: CIC_current_weighting         !< Cloud-in-Cell weighting of particle quantities to the grid.
   procedure, pass(self) :: NGP_current_weighting         !< Nearest Grid Point weighting of particle quantities to the grid.
   procedure, pass(self) :: TSC_current_weighting         !< Triangular Shaped Cloud weighting of particle quantities to the grid.
   procedure, pass(self) :: Gaussian_current_weighting    !< Gaussian Shaped Cloud weighting of particle quantities to the grid.
   procedure, pass(self) :: zeroD_field_weighting
   procedure, pass(self) :: oneD_field_weighting
   procedure, pass(self) :: twoD_field_weighting
   procedure, pass(self) :: threeD_field_weighting
   procedure, pass(self) :: fourD_field_weighting
   procedure, pass(self) :: sixD_field_weighting
   procedure, pass(self) :: fiveD_field_weighting
   procedure, pass(self) :: Gaussian_field_weighting
   procedure, pass(self) :: bspline_field_weighting
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
   implicit none
   class(prism_pic_object), intent(in) :: self             !< External fields.
   character(len=:), allocatable       :: desc             !< Description.
   character(len=1), parameter         :: NL=new_line('a') !< New line character.

   desc =            mpih%myrankstr//'PIC object description:'
   desc = desc//NL//mpih%myrankstr//'    Problem type: '//trim(self%problem_type)
   if (self%problem_type == PLASMA_TYPE_PROBLEM) then
      desc = desc//NL//mpih%myrankstr//'    Plasma initial domain: '//trim(self%plasma_domain)
      desc = desc//NL//mpih%myrankstr//'    Input plasma density [m^(-3)]: '//trim(str(self%plasma_density))
      desc = desc//NL//mpih%myrankstr//'    Neutral fraction: '//trim(str(self%neutral_fraction))
      if (self%plasma_domain == UNIFORM_CILINDER) then
			desc = desc//NL//mpih%myrankstr//'    	Cilinder radius: '//trim(str(self%cilinder_radius))
         desc = desc//NL//mpih%myrankstr//'    	Cilinder length: '//trim(str(self%cilinder_length))
         desc = desc//NL//mpih%myrankstr//'    	Cilinder x centre: '//trim(str(self%cilinder_center(1)))
         desc = desc//NL//mpih%myrankstr//'    	Cilinder y centre: '//trim(str(self%cilinder_center(2)))
         desc = desc//NL//mpih%myrankstr//'    	Cilinder z centre: '//trim(str(self%cilinder_center(3)))
         desc = desc//NL//mpih%myrankstr//'    	Cilinder axis: '//trim(self%cilinder_axis)
		endif
      !desc = desc//NL//mpih%myrankstr//'    Total number of particles: '//trim(str(self%particle_number))
      !desc = desc//NL//mpih%myrankstr//'    of which ions: '//trim(str(self%n_ions))
      !desc = desc//NL//mpih%myrankstr//'    of which electrons: '//trim(str(self%n_electrons))
      !desc = desc//NL//mpih%myrankstr//'    of which neutrals: '//trim(str(self%n_neutrals))
   endif
   desc = desc//NL//mpih%myrankstr//'    Initialization type: '//trim(self%initialization)
   desc = desc//NL//mpih%myrankstr//'    Elliptic correction: '//trim(str(self%elliptic_correction))
   desc = desc//NL//mpih%myrankstr//'    Filter deposition and gather: '//trim(str(self%filter_deposition))
   desc = desc//NL//mpih%myrankstr//'    Particle weighting model: '//trim(self%particle_weighting_model)
   desc = desc//NL//mpih%myrankstr//'    Current weighting model: '//trim(self%current_weighting_model)
   if (self%particle_weighting_model == GAUSSIAN_WEIGHTING_MODEL .or. &
         self%current_weighting_model == GAUSSIAN_WEIGHTING_MODEL .or. &
         self%field_weighting_model == GAUSSIAN_WEIGHTING_MODEL) then
      desc = desc//NL//mpih%myrankstr//'    Sigma: '//trim(str(self%sigma))
      desc = desc//NL//mpih%myrankstr//'    Cutoff_sigma: '//trim(str(self%cutoff_sigma))
      desc = desc//NL//mpih%myrankstr//'    Gaussian support cells: '//trim(str(self%gaussian_support_cells))
   endif
   desc = desc//NL//mpih%myrankstr//'    Field weighting model: '//trim(self%field_weighting_model)
   desc = desc//NL//mpih%myrankstr//'    Numerical scheme for time operator: '//trim(self%scheme_time)
   endfunction description

   subroutine initialize(self, field, grid, file_parameters, physics)
   !< Initialize PIC.
   implicit none
   class(prism_pic_object), intent(inout) :: self            !< Pic object.
   type(field_object),      intent(in)    :: field           !< Field (sibling realm component, threaded in).
   type(grid_object),       intent(in)    :: grid            !< Grid (sibling realm component, threaded in).
   type(file_ini),          intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   type(prism_physics_object), intent(in) :: physics
   real(R8P)                              :: domain_volume   !< Total volume of the computational domain where plasma is
                                                             !< present at t0
   integer(I4P)                           :: b
	character(len=:),        allocatable   :: desc
   character(len=1),        parameter     :: NL=new_line('a')

   print '(A)', mpih%myrankstr//'prism_pic_object%initialize start'

   call self%load_from_file(file_parameters=file_parameters)
   if (is_adim_model(physics%physical_model)) then
      if (self%problem_type == PLASMA_TYPE_PROBLEM) then
         self%plasma_density = self%plasma_density*physics%L0**3
         if (self%plasma_domain == UNIFORM_CILINDER) then
            self%cilinder_radius = self%cilinder_radius/physics%L0
            self%cilinder_length = self%cilinder_length/physics%L0
            self%cilinder_center = self%cilinder_center/physics%L0
         endif
      endif
      if (self%sigma > 0._R8P) self%sigma = self%sigma/physics%L0
   endif
   associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj, nk=>grid%nk, &
               emin=>grid%domain_emin, emax=>grid%domain_emax)

   if (self%problem_type == PLASMA_TYPE_PROBLEM) then
      if (self%plasma_domain==UNIFORM_DOMAIN) then
         domain_volume = (emax(1)-emin(1))*(emax(2)-emin(2))*(emax(3)-emin(3))
      elseif (self%plasma_domain==UNIFORM_CILINDER) then
         domain_volume = PI*self%cilinder_radius**2*self%cilinder_length
      else
         call mpih%error_stop(msg=': invalid plasma_domain for particle number computation')
      endif
      self%particle_number = nint(self%plasma_density*domain_volume)
      self%n_neutrals = nint(self%neutral_fraction*real(self%particle_number,R8P))
	   self%n_ions = nint(real(self%particle_number-self%n_neutrals, R8P)/2.0_R8P)
	   self%n_electrons = self%n_ions
	   self%n_neutrals = self%particle_number - self%n_ions - self%n_electrons
   elseif (self%problem_type == SINGLE_PARTICLE_TYPE_PROBLEM) then
      self%particle_number = 1_I4P
   endif
   if ((self%particle_weighting_model == GAUSSIAN_WEIGHTING_MODEL .or. &
        self%current_weighting_model == GAUSSIAN_WEIGHTING_MODEL .or. &
        self%field_weighting_model == GAUSSIAN_WEIGHTING_MODEL) .and. self%gaussian_support_cells < 0_I4P) then
      self%gaussian_support_cells = 0_I4P
      do b=1, blocks_number
         self%gaussian_support_cells = max(self%gaussian_support_cells, &
            ceiling(self%cutoff_sigma*self%sigma/minval(field%dxyz(:,b)), kind=I4P) + 1_I4P)
      enddo
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
   case(CUBIC_WEIGHTING_MODEL)
      self%particle_weighting => cubic_charge_weighting
   case(QUARTIC_WEIGHTING_MODEL)
      self%particle_weighting => quartic_charge_weighting
   case(QUINTIC_WEIGHTING_MODEL)
      self%particle_weighting => quintic_charge_weighting
   case(SEXTIC_WEIGHTING_MODEL)
      self%particle_weighting => sextic_charge_weighting
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
   case(CUBIC_WEIGHTING_MODEL)
      self%current_weighting => cubic_current_weighting
   case(QUARTIC_WEIGHTING_MODEL)
      self%current_weighting => quartic_current_weighting
   case(QUINTIC_WEIGHTING_MODEL)
      self%current_weighting => quintic_current_weighting
   case(SEXTIC_WEIGHTING_MODEL)
      self%current_weighting => sextic_current_weighting
   case(GAUSSIAN_WEIGHTING_MODEL)
      self%current_weighting => Gaussian_current_weighting
   case default
      call mpih%error_stop(msg=': invalid current weighting model in prism_cpu_object%initialize')
   endselect

   select case(self%field_weighting_model)
   case(ZEROD_FIELDS_WEIGHTING_MODEL)
      self%field_weighting => zeroD_field_weighting
   case(ONED_FIELDS_WEIGHTING_MODEL)
      self%field_weighting => oneD_field_weighting
   case(TWOD_FIELDS_WEIGHTING_MODEL)
      self%field_weighting => twoD_field_weighting
   case(THREED_FIELDS_WEIGHTING_MODEL)
      self%field_weighting => threeD_field_weighting
   case(FOURD_FIELDS_WEIGHTING_MODEL)
      self%field_weighting => fourD_field_weighting
   case(FIVED_FIELDS_WEIGHTING_MODEL)
      self%field_weighting => fiveD_field_weighting
   case(SIXD_FIELDS_WEIGHTING_MODEL)
      self%field_weighting => sixD_field_weighting
   case(GAUSSIAN_WEIGHTING_MODEL)
      self%field_weighting => Gaussian_field_weighting
   case default
      call mpih%error_stop(msg=': invalid field weighting model in prism_cpu_object%initialize')
   endselect

   print '(A)', mpih%myrankstr//'prism_pic_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load PIC configuration from file.
   implicit none
	class(prism_pic_object), intent(inout)   			 :: self             !< PIC object.
	type(file_ini),          intent(in)		  			 :: file_parameters  !< File handler.
   logical,                 intent(in), optional    :: go_on_fail      	!< Go on if load fails.
   logical                                          :: go_on_fail_     	!< Go on if load fails.
   integer(I4P)                                     :: error           	!< Error status.
   character(99)                                    :: buff       		!< Option character buffer.

	go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   buff = ''
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='problem_type', &
   val=buff, error=error)
   if (.not.go_on_fail_.and.error>0) &
   call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(problem_type)')
   self%problem_type = trim(adjustl(buff))

   if(self%problem_type == PLASMA_TYPE_PROBLEM) then
      buff = ''
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='plasma_domain', &
      val=buff, error=error)
      if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(plasma_domain)')
      buff = trim(adjustl(buff))

      select case(trim(adjustl(buff)))
      case('uniform', 'Uniform', 'UNIFORM', 'uniform_domain', 'all', 'full', 'Full', 'FULL')
         self%plasma_domain = UNIFORM_DOMAIN
      case('Cilinder', 'cilinder', 'CILINDER', 'uniform_cilinder', 'UNIFORM_CILINDER')
         self%plasma_domain = UNIFORM_CILINDER
      case default
         call mpih%error_stop(msg=': invalid initial domain for the plasma ['//trim(adjustl(buff))//'] in  &
         ['//INI_SECTION_NAME//'].(plasma_domain)')
      endselect

      if (self%plasma_domain == UNIFORM_CILINDER) then
	      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='cilinder_radius', &
         val=self%cilinder_radius, error=error)
         if (.not.go_on_fail_.and.error>0) &
         call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(cilinder_radius)')

	      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='cilinder_length', &
         val=self%cilinder_length, error=error)
         if (.not.go_on_fail_.and.error>0) &
         call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(cilinder_length)')

         call file_parameters%get(section_name=INI_SECTION_NAME, option_name='cilinder_x_center', &
         val=self%cilinder_center(1), error=error)
         if (.not.go_on_fail_.and.error>0) &
         call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(cilinder_x_center)')

         call file_parameters%get(section_name=INI_SECTION_NAME, option_name='cilinder_y_center', &
         val=self%cilinder_center(2), error=error)
         if (.not.go_on_fail_.and.error>0) &
         call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(cilinder_y_center)')

         call file_parameters%get(section_name=INI_SECTION_NAME, option_name='cilinder_z_center', &
         val=self%cilinder_center(3), error=error)
         if (.not.go_on_fail_.and.error>0) &
         call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(cilinder_z_center)')

	      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='cilinder_axis', &
         val=self%cilinder_axis, error=error)
         if (.not.go_on_fail_.and.error>0) &
         call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(cilinder_axis)')
      endif

	   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='plasma_density', &
      val=self%plasma_density, error=error)
      if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(plasma_density)')

      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='neutral_fraction', &
      val=self%neutral_fraction, error=error)
      if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(neutral_fraction)')
   endif

   buff = COHERENT_INITIALIZATION
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='initialization', &
   val=buff, error=error)
   if (error > 0) then
      self%initialization = COHERENT_INITIALIZATION
   else
      buff = trim(adjustl(buff))
      select case(trim(adjustl(buff)))
      case('', 'COHERENT', 'coherent', 'Coherent')
         self%initialization = COHERENT_INITIALIZATION
      case default
         call mpih%error_stop(msg=': invalid PIC field initialization ['//trim(buff)//']; only "'// &
                              COHERENT_INITIALIZATION//'" is supported')
      endselect
   endif

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='elliptic_correction', &
   val=self%elliptic_correction, error=error)
   if (.not.go_on_fail_.and.error>0) &
   call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(elliptic_correction)')

   self%filter_deposition = .false.
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='filter_deposition', &
                            val=self%filter_deposition, error=error)
   if (error > 0) self%filter_deposition = .false.

   buff = ''
	call file_parameters%get(section_name=INI_SECTION_NAME, option_name='particle_weighting_model', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(particle_weighting_model) from file')
   buff = trim(adjustl(buff))
   select case(trim(adjustl(buff)))
   case('CIC', 'cic', 'Cic')
      self%particle_weighting_model = CIC_WEIGHTING_MODEL
	case('NGP', 'ngp', 'Ngp')
		self%particle_weighting_model = NGP_WEIGHTING_MODEL
	case('TSC', 'tsc', 'Tsc')
		self%particle_weighting_model = TSC_WEIGHTING_MODEL
	case('cubic', '3rd', '3RD')
		self%particle_weighting_model = cubic_WEIGHTING_MODEL
	case('quartic', '4th', '4TH')
		self%particle_weighting_model = quartic_WEIGHTING_MODEL
	case('quintic', '5th', '5TH')
		self%particle_weighting_model = quintic_WEIGHTING_MODEL
   case('sextic', 'Sextic', 'SEXTIC', 'sixth', '6th', '6TH')
      self%particle_weighting_model = SEXTIC_WEIGHTING_MODEL
   case('GAUSSIAN', 'Gaussian', 'gaussian', 'GAUSS', 'gauss')
      self%particle_weighting_model = GAUSSIAN_WEIGHTING_MODEL
	case default
		call mpih%error_stop(msg=': invalid particle weighting model ['//trim(buff)//'] in  &
      ['//INI_SECTION_NAME//'].(particle_weighting_model)')
	endselect

   buff = ''
	call file_parameters%get(section_name=INI_SECTION_NAME, option_name='current_weighting_model', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(current_weighting_model) from file')
   buff = trim(adjustl(buff))
   select case(trim(adjustl(buff)))
   case('CIC', 'cic', 'Cic')
      self%current_weighting_model = CIC_WEIGHTING_MODEL
	case('NGP', 'ngp', 'Ngp')
		self%current_weighting_model = NGP_WEIGHTING_MODEL
	case('TSC', 'tsc', 'Tsc')
		self%current_weighting_model = TSC_WEIGHTING_MODEL
	case('cubic', '3rd', '3RD')
		self%current_weighting_model = cubic_WEIGHTING_MODEL
	case('quartic', '4th', '4TH')
		self%current_weighting_model = quartic_WEIGHTING_MODEL
	case('quintic', '5th', '5TH')
		self%current_weighting_model = quintic_WEIGHTING_MODEL
   case('sextic', 'Sextic', 'SEXTIC', 'sixth', '6th', '6TH')
      self%current_weighting_model = SEXTIC_WEIGHTING_MODEL
   case('GAUSSIAN', 'Gaussian', 'gaussian', 'GAUSS', 'gauss')
      self%current_weighting_model = GAUSSIAN_WEIGHTING_MODEL
	case default
		call mpih%error_stop(msg=': invalid current weighting model ['//trim(buff)//'] in  &
      ['//INI_SECTION_NAME//'].(current_weighting_model)')
	endselect


   buff = ''
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='field_weighting_model', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(field_weighting_model) from file')
   buff = trim(adjustl(buff))
   select case(trim(adjustl(buff)))
   case('0D', '0d', '0_d', '0_D')
      self%field_weighting_model = ZEROD_FIELDS_WEIGHTING_MODEL
   case('1D', '1d', '1_d', '1_D')
      self%field_weighting_model = ONED_FIELDS_WEIGHTING_MODEL
   case('2D', '2d', '2_d', '2_D')
      self%field_weighting_model = TWOD_FIELDS_WEIGHTING_MODEL
   case('3D', '3d', '3_d', '3_D')
      self%field_weighting_model = THREED_FIELDS_WEIGHTING_MODEL
   case('4D', '4d', '4_d', '4_D')
      self%field_weighting_model = FOURD_FIELDS_WEIGHTING_MODEL
   case('5D', '5d', '5_d', '5_D')
      self%field_weighting_model = FIVED_FIELDS_WEIGHTING_MODEL
   case('6D', '6d', '6_d', '6_D')
      self%field_weighting_model = SIXD_FIELDS_WEIGHTING_MODEL
   case('GAUSSIAN', 'Gaussian', 'gaussian', 'GAUSS', 'gauss')
      self%field_weighting_model = GAUSSIAN_WEIGHTING_MODEL
   case default
      call mpih%error_stop(msg=': invalid field weighting model ['//trim(buff)//'] in  &
      ['//INI_SECTION_NAME//'].(field_weighting_model)')
   endselect

   if (self%particle_weighting_model == GAUSSIAN_WEIGHTING_MODEL .or. &
         self%current_weighting_model == GAUSSIAN_WEIGHTING_MODEL .or. &
         self%field_weighting_model == GAUSSIAN_WEIGHTING_MODEL) then
         call file_parameters%get(section_name=INI_SECTION_NAME, option_name='sigma', &
         val=self%sigma, error=error)
         if (.not.go_on_fail_.and.error>0) &
         call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(sigma)')

         call file_parameters%get(section_name=INI_SECTION_NAME, option_name='cutoff_sigma', &
         val=self%cutoff_sigma, error=error)
         if (.not.go_on_fail_.and.error>0) &
         call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(cutoff_sigma)')
         self%gaussian_support_cells = -1_I4P
         call file_parameters%get(section_name=INI_SECTION_NAME, option_name='gaussian_support_cells', &
                                  val=self%gaussian_support_cells, error=error)
         if (error > 0) self%gaussian_support_cells = -1_I4P
         if (self%sigma <= 0._R8P .or. self%cutoff_sigma <= 0._R8P) &
            call mpih%error_stop(msg=': PIC Gaussian weighting requires positive sigma and cutoff_sigma')
   endif

   buff = ''
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='scheme_time', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(scheme_time)')
   buff = trim(adjustl(buff))
   select case(trim(adjustl(buff)))
   case('LEAPFROG', 'leapfrog', 'Leapfrog')
      self%scheme_time = NUM_SCHEME_TIME_PIC_LEAPFROG
   case('RUNGE_KUTTA', 'runge_kutta', 'Runge_Kutta')
      self%scheme_time = NUM_SCHEME_TIME_PIC_RUNGE_KUTTA
   case default
      call mpih%print_message(msg='warning: numerical scheme "'//trim(buff)//'" unknown. Revert back to RK scheme')
      self%scheme_time = NUM_SCHEME_TIME_PIC_RUNGE_KUTTA
   endselect
   endsubroutine load_from_file

   subroutine particle_cartesian_grid_index(self, field, grid, q_pic)
   !< Compute the grid index corresponding to a particle position. Good for cartesian grids only.
   class(prism_pic_object), intent(inout) :: self               !< External fields.
   type(field_object),      intent(in)    :: field              !< The field.
   type(grid_object),       intent(in)    :: grid               !< Grid (sibling realm component, threaded in)
   real(R8P),               intent(in)    :: q_pic(1:,1:)       !< PIC variables.
   integer(I4P)                           :: n, b               !< Counters
   integer(I4P)                           :: i_p, j_p, k_p, b_p !< Particle grid indices

   associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj,                                   &
            nk=>grid%nk, ngc=>grid%ngc, dx => field%dxyz(1,:), dy => field%dxyz(2,:), dz => field%dxyz(3,:), &
            np => self%particle_number, domain_emin => grid%domain_emin, domain_emax => grid%domain_emax,    &
            emin => field%emin, emax => field%emax, neighbour_list => self%neighbour_list)

   !Di sicuro va considerata una parte relativa alle particelle che escono dal dominio
   ! Rivedi con Stefano, molto dipende se quei min max contano pure le gc. In tal caso a emin devi sommare ngc*dx o dx o dz (non
   ! dovrebbero contare)
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
   implicit none
   class(prism_pic_object), intent(inout) :: self                 !< External fields.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),       intent(in)    :: grid                 !< Grid (sibling realm component, threaded in).
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:,&
                                                  1-grid%ngc:,&
                                                  1-grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Number of variables.

   call bspline_charge_weighting(self=self, field=field, grid=grid, q=q, q_pic=q_pic, nv=nv, order=0_I4P)
   endsubroutine NGP_charge_weighting

   subroutine CIC_charge_weighting(self, field, grid, q, q_PIC, nv)
   !< Cloud-in-Cell weighting of particle quantities to the grid.
   implicit none
   class(prism_pic_object), intent(inout) :: self                 !< External fields.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),       intent(in)    :: grid                 !< Grid (sibling realm component, threaded in).
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:, &
                                                  1-grid%ngc:, &
                                                  1-grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_PIC(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Number of variables.

   call bspline_charge_weighting(self=self, field=field, grid=grid, q=q, q_pic=q_pic, nv=nv, order=1_I4P)
   endsubroutine CIC_charge_weighting

   subroutine TSC_charge_weighting(self, field, grid, q, q_pic, nv)
   !!< Triangular Shaped Cloud weighting of particle quantities to the grid.
   implicit none
   class(prism_pic_object), intent(inout) :: self                 !< External fields.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),       intent(in)    :: grid                 !< Grid (sibling realm component, threaded in).
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:, &
                                                  1-grid%ngc:, &
                                                  1-grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Number of variables.

   call bspline_charge_weighting(self=self, field=field, grid=grid, q=q, q_pic=q_pic, nv=nv, order=2_I4P)
   endsubroutine TSC_charge_weighting

   subroutine cubic_charge_weighting(self, field, grid, q, q_pic, nv)
   !!< Cubic B-spline weighting of particle charge to the grid.
   implicit none
   class(prism_pic_object), intent(inout) :: self                 !< PIC object.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),       intent(in)    :: grid                 !< Grid (sibling realm component, threaded in).
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:, &
                                                  1-grid%ngc:, &
                                                  1-grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Charge variable index.

   call bspline_charge_weighting(self=self, field=field, grid=grid, q=q, q_pic=q_pic, nv=nv, order=3_I4P)
   endsubroutine cubic_charge_weighting

   subroutine quartic_charge_weighting(self, field, grid, q, q_pic, nv)
   !!< Quartic B-spline weighting of particle charge to the grid.
   implicit none
   class(prism_pic_object), intent(inout) :: self                 !< PIC object.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),       intent(in)    :: grid                 !< Grid (sibling realm component, threaded in).
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:, &
                                                  1-grid%ngc:, &
                                                  1-grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Charge variable index.

   call bspline_charge_weighting(self=self, field=field, grid=grid, q=q, q_pic=q_pic, nv=nv, order=4_I4P)
   endsubroutine quartic_charge_weighting

   subroutine quintic_charge_weighting(self, field, grid, q, q_pic, nv)
   !!< Quintic B-spline weighting of particle charge to the grid.
   implicit none
   class(prism_pic_object), intent(inout) :: self                 !< PIC object.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),       intent(in)    :: grid                 !< Grid (sibling realm component, threaded in).
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:, &
                                                  1-grid%ngc:, &
                                                  1-grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Charge variable index.

   call bspline_charge_weighting(self=self, field=field, grid=grid, q=q, q_pic=q_pic, nv=nv, order=5_I4P)
   endsubroutine quintic_charge_weighting

   subroutine sextic_charge_weighting(self, field, grid, q, q_pic, nv)
   !!< Sextic B-spline weighting of particle charge to the grid.
   implicit none
   class(prism_pic_object), intent(inout) :: self                 !< PIC object.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),       intent(in)    :: grid                 !< Grid (sibling realm component, threaded in).
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:, &
                                                  1-grid%ngc:, &
                                                  1-grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Charge variable index.

   call bspline_charge_weighting( &
        self=self, field=field, grid=grid, q=q, q_pic=q_pic, nv=nv, order=6_I4P)
   endsubroutine sextic_charge_weighting

   subroutine bspline_charge_weighting(self, field, grid, q, q_pic, nv, order)
   !< Deposit particle charge density with B-spline support and optional binomial filtering.
   implicit none
   class(prism_pic_object), intent(inout) :: self
   type(field_object),      intent(inout) :: field
   type(grid_object),       intent(in)    :: grid
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:, &
                                                  1-grid%ngc:, &
                                                  1-grid%ngc:,1:)
   real(R8P),               intent(in)    :: q_pic(1:,1:)
   integer(I4P),            intent(in)    :: nv
   integer(I4P),            intent(in)    :: order
   integer(I4P)                           :: n, i, j, k
   integer(I4P)                           :: i_p, j_p, k_p, b_p
   integer(I4P)                           :: i_min, i_max
   integer(I4P)                           :: j_min, j_max
   integer(I4P)                           :: k_min, k_max
   real(R8P)                              :: dx, dy, dz
   real(R8P)                              :: wx, wy, wz, weight
   real(R8P)                              :: prefactor

   associate(x_cell=>field%x_cell, y_cell=>field%y_cell, z_cell=>field%z_cell)

   q(nv,:,:,:,:) = 0.0_R8P

   do n = 1, self%particle_number
      b_p = self%neighbour_list(1,n)
      i_p = self%neighbour_list(2,n)
      j_p = self%neighbour_list(3,n)
      k_p = self%neighbour_list(4,n)

      dx = field%dxyz(1,b_p)
      dy = field%dxyz(2,b_p)
      dz = field%dxyz(3,b_p)
      prefactor = q_pic(7,n) / (dx * dy * dz)

      call set_bspline_stencil(order=order, x_p=q_pic(1,n), x_c=x_cell(i_p,b_p), &
                               i_p=i_p, i_min=i_min, i_max=i_max)
      call set_bspline_stencil(order=order, x_p=q_pic(2,n), x_c=y_cell(j_p,b_p), &
                               i_p=j_p, i_min=j_min, i_max=j_max)
      call set_bspline_stencil(order=order, x_p=q_pic(3,n), x_c=z_cell(k_p,b_p), &
                               i_p=k_p, i_min=k_min, i_max=k_max)

      if (self%filter_deposition) then
         i_min = i_min - 1_I4P ; i_max = i_max + 1_I4P
         j_min = j_min - 1_I4P ; j_max = j_max + 1_I4P
         k_min = k_min - 1_I4P ; k_max = k_max + 1_I4P
      endif

      i_min = max(i_min, lbound(q,dim=2)) ; i_max = min(i_max, ubound(q,dim=2))
      j_min = max(j_min, lbound(q,dim=3)) ; j_max = min(j_max, ubound(q,dim=3))
      k_min = max(k_min, lbound(q,dim=4)) ; k_max = min(k_max, ubound(q,dim=4))

      do k = k_min, k_max
         wz = effective_bspline_weight(order=order, r=(q_pic(3,n) - z_cell(k,b_p)) / dz, &
                                       filter=self%filter_deposition)
         do j = j_min, j_max
            wy = effective_bspline_weight(order=order, r=(q_pic(2,n) - y_cell(j,b_p)) / dy, &
                                          filter=self%filter_deposition)
            do i = i_min, i_max
               wx = effective_bspline_weight(order=order, r=(q_pic(1,n) - x_cell(i,b_p)) / dx, &
                                             filter=self%filter_deposition)
               weight = wx * wy * wz
               q(nv,i,j,k,b_p) = q(nv,i,j,k,b_p) + prefactor * weight
            enddo
         enddo
      enddo
   enddo

   endassociate
   endsubroutine bspline_charge_weighting

   subroutine Gaussian_charge_weighting(self, field, grid, q, q_PIC, nv)
   !< Gaussian weighting of particle charge density to the grid.
   implicit none
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
   real(R8P)                              :: cutoff_limit            !< Cutoff including roundoff tolerance.
   real(R8P)                              :: sx, sy, sz              !< Grid spacing in Gaussian-width units.

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
      cutoff_limit = cutoff_sigma + 64._R8P*epsilon(cutoff_sigma)*max(1._R8P, abs(cutoff_sigma))
      sx = dx / sigma_x
      sy = dy / sigma_y
      sz = dz / sigma_z

      inverse_cell_volume = 1._R8P / (dx*dy*dz)

      ! Fixed compact support radius for Gaussian loops.
      ni_sigma = self%gaussian_support_cells
      nj_sigma = self%gaussian_support_cells
      nk_sigma = self%gaussian_support_cells
      if (self%filter_deposition) then
         ni_sigma = ni_sigma + 1_I4P
         nj_sigma = nj_sigma + 1_I4P
         nk_sigma = nk_sigma + 1_I4P
      endif

      ! Restrict the support to the locally available grid, including ghost cells. !Parte da rivedere, perchè alla frontiera hai una
      ! distribuzione asimmetrica, che taglia e riscala di conseguenza
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
         wz = effective_gaussian_weight(r=rz, shift=sz, cutoff_limit=cutoff_limit, filter=self%filter_deposition)
         if (wz <= 0._R8P) cycle
         do j=j_min, j_max
            ry = (q_PIC(2,n)-y_cell(j,b_p))/sigma_y
            wy = effective_gaussian_weight(r=ry, shift=sy, cutoff_limit=cutoff_limit, filter=self%filter_deposition)
            if (wy <= 0._R8P) cycle
            do i=i_min, i_max
               rx = (q_PIC(1,n)-x_cell(i,b_p))/sigma_x
               wx = effective_gaussian_weight(r=rx, shift=sx, cutoff_limit=cutoff_limit, filter=self%filter_deposition)
               if (wx <= 0._R8P) cycle
               weight_sum = weight_sum + wx*wy*wz
            enddo
         enddo
      enddo

      ! Deposit the particle charge density.
      if (weight_sum > tiny(1._R8P)) then
         do k=k_min, k_max
            rz = (q_PIC(3,n)-z_cell(k,b_p))/sigma_z
            wz = effective_gaussian_weight(r=rz, shift=sz, cutoff_limit=cutoff_limit, filter=self%filter_deposition)
            if (wz <= 0._R8P) cycle
            do j=j_min, j_max
               ry = (q_PIC(2,n)-y_cell(j,b_p))/sigma_y
               wy = effective_gaussian_weight(r=ry, shift=sy, cutoff_limit=cutoff_limit, filter=self%filter_deposition)
               if (wy <= 0._R8P) cycle
               do i=i_min, i_max
                  rx = (q_PIC(1,n)-x_cell(i,b_p))/sigma_x
                  wx = effective_gaussian_weight(r=rx, shift=sx, cutoff_limit=cutoff_limit, filter=self%filter_deposition)
                  if (wx <= 0._R8P) cycle
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
   implicit none
   class(prism_pic_object), intent(inout) :: self                 !< External fields.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),       intent(in)    :: grid                 !< Grid (sibling realm component, threaded in).
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:,&
                                                  1-grid%ngc:,&
                                                  1-grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Number of variables.

   call bspline_current_weighting(self=self, field=field, grid=grid, q=q, q_pic=q_pic, nv=nv, order=0_I4P)
   endsubroutine NGP_current_weighting

   subroutine CIC_current_weighting(self, field, grid, q, q_pic, nv)
   !< Cloud-in-Cell weighting of particle quantities to the grid.
   implicit none
   class(prism_pic_object), intent(inout) :: self                 !< External fields.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),       intent(in)    :: grid                 !< Grid (sibling realm component, threaded in).
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:, &
                                                  1-grid%ngc:, &
                                                  1-grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Number of variables.

   call bspline_current_weighting(self=self, field=field, grid=grid, q=q, q_pic=q_pic, nv=nv, order=1_I4P)
   endsubroutine CIC_current_weighting

   subroutine TSC_current_weighting(self, field, grid, q, q_pic, nv)
   !< Triangular-Shaped-Cloud weighting of particle current density to the grid.
   implicit none
   class(prism_pic_object), intent(inout) :: self                 !< External fields.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),       intent(in)    :: grid                 !< Grid (sibling realm component, threaded in).
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:, &
                                                  1-grid%ngc:, &
                                                  1-grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Number of variables.

   call bspline_current_weighting(self=self, field=field, grid=grid, q=q, q_pic=q_pic, nv=nv, order=2_I4P)
   endsubroutine TSC_current_weighting

   subroutine cubic_current_weighting(self, field, grid, q, q_pic, nv)
   !< Cubic B-spline weighting of particle current density to the grid.
   implicit none
   class(prism_pic_object), intent(inout) :: self                 !< External fields.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),       intent(in)    :: grid                 !< Grid (sibling realm component, threaded in).
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:, &
                                                  1-grid%ngc:, &
                                                  1-grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Number of variables.

   call bspline_current_weighting(self=self, field=field, grid=grid, q=q, q_pic=q_pic, nv=nv, order=3_I4P)
   endsubroutine cubic_current_weighting

   subroutine quartic_current_weighting(self, field, grid, q, q_pic, nv)
   !< Quartic B-spline weighting of particle current density to the grid.
   implicit none
   class(prism_pic_object), intent(inout) :: self                 !< External fields.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),       intent(in)    :: grid                 !< Grid (sibling realm component, threaded in).
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:, &
                                                  1-grid%ngc:, &
                                                  1-grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Number of variables.

   call bspline_current_weighting(self=self, field=field, grid=grid, q=q, q_pic=q_pic, nv=nv, order=4_I4P)
   endsubroutine quartic_current_weighting

   subroutine quintic_current_weighting(self, field, grid, q, q_pic, nv)
   !< Quintic B-spline weighting of particle current density to the grid.
   implicit none
   class(prism_pic_object), intent(inout) :: self                 !< External fields.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),       intent(in)    :: grid                 !< Grid (sibling realm component, threaded in).
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:, &
                                                  1-grid%ngc:, &
                                                  1-grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Number of variables.

   call bspline_current_weighting(self=self, field=field, grid=grid, q=q, q_pic=q_pic, nv=nv, order=5_I4P)
   endsubroutine quintic_current_weighting

   subroutine sextic_current_weighting(self, field, grid, q, q_pic, nv)
   !< Sextic B-spline weighting of particle current density to the grid.
   implicit none
   class(prism_pic_object), intent(inout) :: self                 !< External fields.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),       intent(in)    :: grid                 !< Grid (sibling realm component, threaded in).
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:, &
                                                  1-grid%ngc:, &
                                                  1-grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Number of variables.

   call bspline_current_weighting( &
        self=self, field=field, grid=grid, q=q, q_pic=q_pic, nv=nv, order=6_I4P)
   endsubroutine sextic_current_weighting

   subroutine bspline_current_weighting(self, field, grid, q, q_pic, nv, order)
   !< Deposit particle current density with B-spline support and optional binomial filtering.
   implicit none
   class(prism_pic_object), intent(inout) :: self
   type(field_object),      intent(inout) :: field
   type(grid_object),       intent(in)    :: grid
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:, &
                                                  1-grid%ngc:, &
                                                  1-grid%ngc:,1:)
   real(R8P),               intent(in)    :: q_pic(1:,1:)
   integer(I4P),            intent(in)    :: nv
   integer(I4P),            intent(in)    :: order
   integer(I4P)                           :: n, i, j, k
   integer(I4P)                           :: i_p, j_p, k_p, b_p
   integer(I4P)                           :: i_min, i_max
   integer(I4P)                           :: j_min, j_max
   integer(I4P)                           :: k_min, k_max
   real(R8P)                              :: dx, dy, dz
   real(R8P)                              :: wx, wy, wz, weight
   real(R8P)                              :: prefactor
   real(R8P)                              :: jx, jy, jz

   associate(x_cell=>field%x_cell, y_cell=>field%y_cell, z_cell=>field%z_cell)

   q(nv-3:nv-1,:,:,:,:) = 0.0_R8P

   do n = 1, self%particle_number
      b_p = self%neighbour_list(1,n)
      i_p = self%neighbour_list(2,n)
      j_p = self%neighbour_list(3,n)
      k_p = self%neighbour_list(4,n)

      dx = field%dxyz(1,b_p)
      dy = field%dxyz(2,b_p)
      dz = field%dxyz(3,b_p)
      prefactor = q_pic(7,n) / (dx * dy * dz)
      jx = prefactor * q_pic(4,n)
      jy = prefactor * q_pic(5,n)
      jz = prefactor * q_pic(6,n)

      call set_bspline_stencil(order=order, x_p=q_pic(1,n), x_c=x_cell(i_p,b_p), &
                               i_p=i_p, i_min=i_min, i_max=i_max)
      call set_bspline_stencil(order=order, x_p=q_pic(2,n), x_c=y_cell(j_p,b_p), &
                               i_p=j_p, i_min=j_min, i_max=j_max)
      call set_bspline_stencil(order=order, x_p=q_pic(3,n), x_c=z_cell(k_p,b_p), &
                               i_p=k_p, i_min=k_min, i_max=k_max)

      if (self%filter_deposition) then
         i_min = i_min - 1_I4P ; i_max = i_max + 1_I4P
         j_min = j_min - 1_I4P ; j_max = j_max + 1_I4P
         k_min = k_min - 1_I4P ; k_max = k_max + 1_I4P
      endif

      i_min = max(i_min, lbound(q,dim=2)) ; i_max = min(i_max, ubound(q,dim=2))
      j_min = max(j_min, lbound(q,dim=3)) ; j_max = min(j_max, ubound(q,dim=3))
      k_min = max(k_min, lbound(q,dim=4)) ; k_max = min(k_max, ubound(q,dim=4))

      do k = k_min, k_max
         wz = effective_bspline_weight(order=order, r=(q_pic(3,n) - z_cell(k,b_p)) / dz, &
                                       filter=self%filter_deposition)
         do j = j_min, j_max
            wy = effective_bspline_weight(order=order, r=(q_pic(2,n) - y_cell(j,b_p)) / dy, &
                                          filter=self%filter_deposition)
            do i = i_min, i_max
               wx = effective_bspline_weight(order=order, r=(q_pic(1,n) - x_cell(i,b_p)) / dx, &
                                             filter=self%filter_deposition)
               weight = wx * wy * wz
               q(nv-3,i,j,k,b_p) = q(nv-3,i,j,k,b_p) + jx * weight
               q(nv-2,i,j,k,b_p) = q(nv-2,i,j,k,b_p) + jy * weight
               q(nv-1,i,j,k,b_p) = q(nv-1,i,j,k,b_p) + jz * weight
            enddo
         enddo
      enddo
   enddo

   endassociate
   endsubroutine bspline_current_weighting

   subroutine Gaussian_current_weighting(self, field, grid, q, q_PIC, nv)
   !< Gaussian weighting of particle current density to the grid.
   implicit none
   class(prism_pic_object), intent(inout) :: self                    !< PIC object.
   type(field_object),      intent(inout) :: field                   !< The field.
   type(grid_object),       intent(in)    :: grid                    !< Grid object.
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:,  &
                                                   1-grid%ngc:, &
                                                   1-grid%ngc:,1:)   !< Field variables.
   real(R8P),               intent(in)    :: q_PIC(1:,1:)            !< PIC variables.
   integer(I4P),            intent(in)    :: nv                      !< Index following Jx, Jy, Jz.
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
   real(R8P)                              :: weight                  !< Normalized three-dimensional weight.
   real(R8P)                              :: weight_sum              !< Discrete normalization factor.
   real(R8P)                              :: inverse_cell_volume     !< Inverse cell volume.
   real(R8P)                              :: current_prefactor(3)    !< Qp*vp/volume.
   real(R8P)                              :: cutoff_limit            !< Cutoff including roundoff tolerance.
   real(R8P)                              :: sx, sy, sz              !< Grid spacing in Gaussian-width units.

   associate(x_cell=>field%x_cell, y_cell=>field%y_cell, z_cell=>field%z_cell, sigma=>self%sigma, cutoff_sigma=>self%cutoff_sigma)

   ! Reset the current density before depositing the current particle distribution.
   q(nv-3:nv-1,:,:,:,:) = 0._R8P

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
      cutoff_limit = cutoff_sigma + 64._R8P*epsilon(cutoff_sigma)*max(1._R8P, abs(cutoff_sigma))
      sx = dx/sigma_x
      sy = dy/sigma_y
      sz = dz/sigma_z
      inverse_cell_volume = 1._R8P/(dx*dy*dz)
      current_prefactor(1) = q_PIC(7,n)*q_PIC(4,n)*inverse_cell_volume
      current_prefactor(2) = q_PIC(7,n)*q_PIC(5,n)*inverse_cell_volume
      current_prefactor(3) = q_PIC(7,n)*q_PIC(6,n)*inverse_cell_volume

      ! Fixed compact support radius for Gaussian loops.
      ni_sigma = self%gaussian_support_cells
      nj_sigma = self%gaussian_support_cells
      nk_sigma = self%gaussian_support_cells
      if (self%filter_deposition) then
         ni_sigma = ni_sigma + 1_I4P
         nj_sigma = nj_sigma + 1_I4P
         nk_sigma = nk_sigma + 1_I4P
      endif

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
         wz = effective_gaussian_weight(r=rz, shift=sz, cutoff_limit=cutoff_limit, filter=self%filter_deposition)
         if (wz <= 0._R8P) cycle
         do j=j_min, j_max
            ry = (q_PIC(2,n)-y_cell(j,b_p))/sigma_y
            wy = effective_gaussian_weight(r=ry, shift=sy, cutoff_limit=cutoff_limit, filter=self%filter_deposition)
            if (wy <= 0._R8P) cycle
            do i=i_min, i_max
               rx = (q_PIC(1,n)-x_cell(i,b_p))/sigma_x
               wx = effective_gaussian_weight(r=rx, shift=sx, cutoff_limit=cutoff_limit, filter=self%filter_deposition)
               if (wx <= 0._R8P) cycle
               weight_sum = weight_sum + wx*wy*wz
            enddo
         enddo
      enddo

      ! Deposit the particle current density.
      if (weight_sum > tiny(1._R8P)) then
         do k=k_min, k_max
            rz = (q_PIC(3,n)-z_cell(k,b_p))/sigma_z
            wz = effective_gaussian_weight(r=rz, shift=sz, cutoff_limit=cutoff_limit, filter=self%filter_deposition)
            if (wz <= 0._R8P) cycle
            do j=j_min, j_max
               ry = (q_PIC(2,n)-y_cell(j,b_p))/sigma_y
               wy = effective_gaussian_weight(r=ry, shift=sy, cutoff_limit=cutoff_limit, filter=self%filter_deposition)
               if (wy <= 0._R8P) cycle
               do i=i_min, i_max
                  rx = (q_PIC(1,n)-x_cell(i,b_p))/sigma_x
                  wx = effective_gaussian_weight(r=rx, shift=sx, cutoff_limit=cutoff_limit, filter=self%filter_deposition)
                  if (wx <= 0._R8P) cycle
                  weight = wx*wy*wz/weight_sum
                  q(nv-3,i,j,k,b_p) = q(nv-3,i,j,k,b_p) + &
                                       current_prefactor(1)*weight
                  q(nv-2,i,j,k,b_p) = q(nv-2,i,j,k,b_p) + &
                                       current_prefactor(2)*weight
                  q(nv-1,i,j,k,b_p) = q(nv-1,i,j,k,b_p) + &
                                       current_prefactor(3)*weight
               enddo
            enddo
         enddo
      endif
   enddo
   endassociate
   end subroutine Gaussian_current_weighting

      subroutine zeroD_field_weighting(self, field, grid, pic_fields, q, q_pic, nv)
   !< Zeroth-order spatial interpolation of cell-centered fields to particle locations.
   class(prism_pic_object), intent(inout) :: self
   type(field_object),      intent(inout) :: field
   type(grid_object),       intent(in)    :: grid
   real(R8P),               intent(inout) :: pic_fields(1:,1:)
   real(R8P),               intent(in)    :: q(1:,1-grid%ngc:, &
                                                  1-grid%ngc:, &
                                                  1-grid%ngc:,1:)
   real(R8P),               intent(in)    :: q_pic(1:,1:)
   integer(I4P),            intent(in)    :: nv

   call bspline_field_weighting(self       = self,       &
                                field      = field,      &
                                grid       = grid,       &
                                pic_fields = pic_fields, &
                                q          = q,          &
                                q_pic      = q_pic,      &
                                nv         = nv,         &
                                order      = 0_I4P)
   endsubroutine zeroD_field_weighting

   subroutine oneD_field_weighting(self, field, grid, pic_fields, q, q_pic, nv)
   !< First-order spatial interpolation of cell-centered fields to particle locations.
   class(prism_pic_object), intent(inout) :: self
   type(field_object),      intent(inout) :: field
   type(grid_object),       intent(in)    :: grid
   real(R8P),               intent(inout) :: pic_fields(1:,1:)
   real(R8P),               intent(in)    :: q(1:,1-grid%ngc:, &
                                                  1-grid%ngc:, &
                                                  1-grid%ngc:,1:)
   real(R8P),               intent(in)    :: q_pic(1:,1:)
   integer(I4P),            intent(in)    :: nv

   call bspline_field_weighting(self       = self,       &
                                field      = field,      &
                                grid       = grid,       &
                                pic_fields = pic_fields, &
                                q          = q,          &
                                q_pic      = q_pic,      &
                                nv         = nv,         &
                                order      = 1_I4P)
   endsubroutine oneD_field_weighting

   subroutine twoD_field_weighting(self, field, grid, pic_fields, q, q_pic, nv)
   !< Second-order spatial interpolation of cell-centered fields to particle locations.
   class(prism_pic_object), intent(inout) :: self
   type(field_object),      intent(inout) :: field
   type(grid_object),       intent(in)    :: grid
   real(R8P),               intent(inout) :: pic_fields(1:,1:)
   real(R8P),               intent(in)    :: q(1:,1-grid%ngc:, &
                                                  1-grid%ngc:, &
                                                  1-grid%ngc:,1:)
   real(R8P),               intent(in)    :: q_pic(1:,1:)
   integer(I4P),            intent(in)    :: nv

   call bspline_field_weighting(self       = self,       &
                                field      = field,      &
                                grid       = grid,       &
                                pic_fields = pic_fields, &
                                q          = q,          &
                                q_pic      = q_pic,      &
                                nv         = nv,         &
                                order      = 2_I4P)
   endsubroutine twoD_field_weighting

   subroutine threeD_field_weighting(self, field, grid, pic_fields, q, q_pic, nv)
   !< Third-order spatial interpolation of cell-centered fields to particle locations.
   class(prism_pic_object), intent(inout) :: self
   type(field_object),      intent(inout) :: field
   type(grid_object),       intent(in)    :: grid
   real(R8P),               intent(inout) :: pic_fields(1:,1:)
   real(R8P),               intent(in)    :: q(1:,1-grid%ngc:, &
                                                  1-grid%ngc:, &
                                                  1-grid%ngc:,1:)
   real(R8P),               intent(in)    :: q_pic(1:,1:)
   integer(I4P),            intent(in)    :: nv

   call bspline_field_weighting(self       = self,       &
                                field      = field,      &
                                grid       = grid,       &
                                pic_fields = pic_fields, &
                                q          = q,          &
                                q_pic      = q_pic,      &
                                nv         = nv,         &
                                order      = 3_I4P)
   endsubroutine threeD_field_weighting

   subroutine fourD_field_weighting(self, field, grid, pic_fields, q, q_pic, nv)
   !< Fourth-order spatial interpolation of cell-centered fields to particle locations.
   class(prism_pic_object), intent(inout) :: self
   type(field_object),      intent(inout) :: field
   type(grid_object),       intent(in)    :: grid
   real(R8P),               intent(inout) :: pic_fields(1:,1:)
   real(R8P),               intent(in)    :: q(1:,1-grid%ngc:, &
                                                  1-grid%ngc:, &
                                                  1-grid%ngc:,1:)
   real(R8P),               intent(in)    :: q_pic(1:,1:)
   integer(I4P),            intent(in)    :: nv

   call bspline_field_weighting(self       = self,       &
                                field      = field,      &
                                grid       = grid,       &
                                pic_fields = pic_fields, &
                                q          = q,          &
                                q_pic      = q_pic,      &
                                nv         = nv,         &
                                order      = 4_I4P)
   endsubroutine fourD_field_weighting

   subroutine fiveD_field_weighting(self, field, grid, pic_fields, q, q_pic, nv)
   !< Fifth-order spatial interpolation of cell-centered fields to particle locations.
   class(prism_pic_object), intent(inout) :: self
   type(field_object),      intent(inout) :: field
   type(grid_object),       intent(in)    :: grid
   real(R8P),               intent(inout) :: pic_fields(1:,1:)
   real(R8P),               intent(in)    :: q(1:,1-grid%ngc:, &
                                                  1-grid%ngc:, &
                                                  1-grid%ngc:,1:)
   real(R8P),               intent(in)    :: q_pic(1:,1:)
   integer(I4P),            intent(in)    :: nv

   call bspline_field_weighting(self       = self,       &
                                field      = field,      &
                                grid       = grid,       &
                                pic_fields = pic_fields, &
                                q          = q,          &
                                q_pic      = q_pic,      &
                                nv         = nv,         &
                                order      = 5_I4P)
   endsubroutine fiveD_field_weighting

   subroutine sixD_field_weighting(self, field, grid, pic_fields, q, q_pic, nv)
   !< Sixth-order spatial interpolation of cell-centered fields to particle locations.
   implicit none
   class(prism_pic_object), intent(inout) :: self
   type(field_object),      intent(inout) :: field
   type(grid_object),       intent(in)    :: grid
   real(R8P),               intent(inout) :: pic_fields(1:,1:)
   real(R8P),               intent(in)    :: q(1:,1-grid%ngc:, &
                                                  1-grid%ngc:, &
                                                  1-grid%ngc:,1:)
   real(R8P),               intent(in)    :: q_pic(1:,1:)
   integer(I4P),            intent(in)    :: nv

   call bspline_field_weighting(self       = self,       &
                                field      = field,      &
                                grid       = grid,       &
                                pic_fields = pic_fields, &
                                q          = q,          &
                                q_pic      = q_pic,      &
                                nv         = nv,         &
                                order      = 6_I4P)
   endsubroutine sixD_field_weighting

   subroutine Gaussian_field_weighting(self, field, grid, pic_fields, q, q_PIC, nv)
   !< Gaussian interpolation of cell-centered fields to particle locations.
   implicit none
   class(prism_pic_object), intent(inout) :: self                    !< PIC object.
   type(field_object),      intent(inout) :: field                   !< The field.
   type(grid_object),       intent(in)    :: grid                    !< Grid object.
   real(R8P),               intent(inout) :: pic_fields(1:,1:)       !< Fields at particle locations.
   real(R8P),               intent(in)    :: q(1:,1-grid%ngc:, &
                                                  1-grid%ngc:, &
                                                  1-grid%ngc:,1:)   !< Field variables.
   real(R8P),               intent(in)    :: q_PIC(1:,1:)            !< PIC variables.
   integer(I4P),            intent(in)    :: nv                      !< Number of variables.
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
   real(R8P)                              :: weight                  !< Normalized three-dimensional weight.
   real(R8P)                              :: weight_sum              !< Discrete normalization factor.
   real(R8P)                              :: cutoff_limit            !< Cutoff including roundoff tolerance.
   real(R8P)                              :: sx, sy, sz              !< Grid spacing in Gaussian-width units.

   associate(x_cell      => field%x_cell,       &
             y_cell      => field%y_cell,       &
             z_cell      => field%z_cell,       &
             sigma       => self%sigma,         &
             cutoff_sigma=> self%cutoff_sigma)

   do n=1, self%particle_number
      ! Get particle block and grid indices.
      b_p = self%neighbour_list(1,n)
      i_p = self%neighbour_list(2,n)
      j_p = self%neighbour_list(3,n)
      k_p = self%neighbour_list(4,n)

      ! Get local grid spacings.
      dx = field%dxyz(1,b_p)
      dy = field%dxyz(2,b_p)
      dz = field%dxyz(3,b_p)

      ! Isotropic Gaussian width.
      sigma_x = sigma
      sigma_y = sigma
      sigma_z = sigma
      cutoff_limit = cutoff_sigma + 64._R8P*epsilon(cutoff_sigma)*max(1._R8P, abs(cutoff_sigma))
      sx = dx / sigma_x
      sy = dy / sigma_y
      sz = dz / sigma_z

      ! Fixed compact support radius for Gaussian loops.
      ni_sigma = self%gaussian_support_cells
      nj_sigma = self%gaussian_support_cells
      nk_sigma = self%gaussian_support_cells
      if (self%filter_deposition) then
         ni_sigma = ni_sigma + 1_I4P
         nj_sigma = nj_sigma + 1_I4P
         nk_sigma = nk_sigma + 1_I4P
      endif

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
         wz = effective_gaussian_weight(r=rz, shift=sz, cutoff_limit=cutoff_limit, filter=self%filter_deposition)
         if (wz <= 0._R8P) cycle

         do j=j_min, j_max
            ry = (q_PIC(2,n)-y_cell(j,b_p))/sigma_y
            wy = effective_gaussian_weight(r=ry, shift=sy, cutoff_limit=cutoff_limit, filter=self%filter_deposition)
            if (wy <= 0._R8P) cycle

            do i=i_min, i_max
               rx = (q_PIC(1,n)-x_cell(i,b_p))/sigma_x
               wx = effective_gaussian_weight(r=rx, shift=sx, cutoff_limit=cutoff_limit, filter=self%filter_deposition)
               if (wx <= 0._R8P) cycle

               weight_sum = weight_sum + wx*wy*wz
            enddo
         enddo
      enddo

      ! Gaussian field gather.
      pic_fields(1:6,n) = 0._R8P

      if (weight_sum > tiny(1._R8P)) then
         do k=k_min, k_max
            rz = (q_PIC(3,n)-z_cell(k,b_p))/sigma_z
            wz = effective_gaussian_weight(r=rz, shift=sz, cutoff_limit=cutoff_limit, filter=self%filter_deposition)
            if (wz <= 0._R8P) cycle

            do j=j_min, j_max
               ry = (q_PIC(2,n)-y_cell(j,b_p))/sigma_y
               wy = effective_gaussian_weight(r=ry, shift=sy, cutoff_limit=cutoff_limit, filter=self%filter_deposition)
               if (wy <= 0._R8P) cycle

               do i=i_min, i_max
                  rx = (q_PIC(1,n)-x_cell(i,b_p))/sigma_x
                  wx = effective_gaussian_weight(r=rx, shift=sx, cutoff_limit=cutoff_limit, filter=self%filter_deposition)
                  if (wx <= 0._R8P) cycle

                  weight = wx*wy*wz/weight_sum

                  pic_fields(1:6,n) = pic_fields(1:6,n) + weight * &
                                      q(1:6,i,j,k,b_p)
               enddo
            enddo
         enddo
      endif
   enddo

   endassociate
   endsubroutine Gaussian_field_weighting

   subroutine bspline_field_weighting(self, field, grid, pic_fields, q, q_pic, nv, order)
   !< B-spline interpolation of cell-centered fields to particle locations.
   implicit none
   class(prism_pic_object), intent(inout) :: self
   type(field_object),      intent(inout) :: field
   type(grid_object),       intent(in)    :: grid
   real(R8P),               intent(inout) :: pic_fields(1:,1:)
   real(R8P),               intent(in)    :: q(1:,1-grid%ngc:, &
                                                  1-grid%ngc:, &
                                                  1-grid%ngc:,1:)
   real(R8P),               intent(in)    :: q_pic(1:,1:)
   integer(I4P),            intent(in)    :: nv
   integer(I4P),            intent(in)    :: order
   integer(I4P)                           :: n, i, j, k
   integer(I4P)                           :: i_p, j_p, k_p, block_p
   integer(I4P)                           :: i_min, i_max
   integer(I4P)                           :: j_min, j_max
   integer(I4P)                           :: k_min, k_max
   real(R8P)                              :: dx, dy, dz
   real(R8P)                              :: rx, ry, rz
   real(R8P)                              :: wx, wy, wz
   real(R8P)                              :: weight

   associate(x_cell=>field%x_cell, y_cell=>field%y_cell, z_cell=>field%z_cell)
   do n = 1, self%particle_number
      ! Get particle grid indices.
      block_p = self%neighbour_list(1,n)
      i_p     = self%neighbour_list(2,n)
      j_p     = self%neighbour_list(3,n)
      k_p     = self%neighbour_list(4,n)
      ! Get grid spacing.
      dx = field%dxyz(1,block_p)
      dy = field%dxyz(2,block_p)
      dz = field%dxyz(3,block_p)
      ! Select B-spline stencil along x.
      call set_bspline_stencil(order = order,               &
                               x_p   = q_pic(1,n),          &
                               x_c   = x_cell(i_p,block_p), &
                               i_p   = i_p,                 &
                               i_min = i_min,               &
                               i_max = i_max)
      ! Select B-spline stencil along y.
      call set_bspline_stencil(order = order,               &
                               x_p   = q_pic(2,n),          &
                               x_c   = y_cell(j_p,block_p), &
                               i_p   = j_p,                 &
                               i_min = j_min,               &
                               i_max = j_max)
      ! Select B-spline stencil along z.
      call set_bspline_stencil(order = order,               &
                               x_p   = q_pic(3,n),          &
                               x_c   = z_cell(k_p,block_p), &
                               i_p   = k_p,                 &
                               i_min = k_min,               &
                               i_max = k_max)
      if (self%filter_deposition) then
         i_min = i_min - 1_I4P ; i_max = i_max + 1_I4P
         j_min = j_min - 1_I4P ; j_max = j_max + 1_I4P
         k_min = k_min - 1_I4P ; k_max = k_max + 1_I4P
      endif
      i_min = max(i_min, lbound(q,dim=2)) ; i_max = min(i_max, ubound(q,dim=2))
      j_min = max(j_min, lbound(q,dim=3)) ; j_max = min(j_max, ubound(q,dim=3))
      k_min = max(k_min, lbound(q,dim=4)) ; k_max = min(k_max, ubound(q,dim=4))
      ! B-spline field gather.
      pic_fields(1:6,n) = 0.0_R8P
      do k = k_min, k_max
         rz = (q_pic(3,n) - z_cell(k,block_p)) / dz
         wz = effective_bspline_weight(order=order, r=rz, filter=self%filter_deposition)
         do j = j_min, j_max
            ry = (q_pic(2,n) - y_cell(j,block_p)) / dy
            wy = effective_bspline_weight(order=order, r=ry, filter=self%filter_deposition)
            do i = i_min, i_max
               rx = (q_pic(1,n) - x_cell(i,block_p)) / dx
               wx = effective_bspline_weight(order=order, r=rx, filter=self%filter_deposition)
               weight = wx * wy * wz
               pic_fields(1:6,n) = pic_fields(1:6,n) + weight * &
                                   q(1:6,i,j,k,block_p)
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine bspline_field_weighting

   pure subroutine set_bspline_stencil(order, x_p, x_c, i_p, i_min, i_max)
   !< Compute the one-dimensional B-spline stencil.
   implicit none
   integer(I4P), intent(in)  :: order
   real(R8P),    intent(in)  :: x_p, x_c
   integer(I4P), intent(in)  :: i_p
   integer(I4P), intent(out) :: i_min, i_max

   select case(order)
   case(0_I4P)
      i_min = i_p
      i_max = i_p
   case(1_I4P)
      if (x_p <= x_c) then
         i_min = i_p - 1_I4P
         i_max = i_p
      else
         i_min = i_p
         i_max = i_p + 1_I4P
      endif
   case(2_I4P)
      i_min = i_p - 1_I4P
      i_max = i_p + 1_I4P
   case(3_I4P)
      if (x_p <= x_c) then
         i_min = i_p - 2_I4P
         i_max = i_p + 1_I4P
      else
         i_min = i_p - 1_I4P
         i_max = i_p + 2_I4P
      endif
   case(4_I4P)
      i_min = i_p - 2_I4P
      i_max = i_p + 2_I4P
   case(5_I4P)
      if (x_p <= x_c) then
         i_min = i_p - 3_I4P
         i_max = i_p + 2_I4P
      else
         i_min = i_p - 2_I4P
         i_max = i_p + 3_I4P
      endif
   case(6_I4P)
      i_min = i_p - 3_I4P
      i_max = i_p + 3_I4P
   case default
      i_min = i_p
      i_max = i_p
   endselect
   endsubroutine set_bspline_stencil

   pure function effective_bspline_weight(order, r, filter) result(weight)
   !< Return the optionally binomial-filtered one-dimensional B-spline weight.
   implicit none
   integer(I4P), intent(in) :: order
   real(R8P),    intent(in) :: r
   logical,      intent(in) :: filter
   real(R8P)                :: weight
   real(R8P), parameter     :: tol = 64.0_R8P * epsilon(1.0_R8P)

   if (filter) then
      if (order == 0_I4P) then
         if (r > -0.5_R8P + tol .and. r <= 0.5_R8P + tol) then
            weight = 0.50_R8P
         elseif ((r > 0.5_R8P - tol .and. r <= 1.5_R8P + tol) .or. &
                 (r > -1.5_R8P - tol .and. r <= -0.5_R8P + tol)) then
            weight = 0.25_R8P
         else
            weight = 0.0_R8P
         endif
      else
         weight = 0.25_R8P*bspline_weight(order=order, r=r + 1.0_R8P) &
                + 0.50_R8P*bspline_weight(order=order, r=r)           &
                + 0.25_R8P*bspline_weight(order=order, r=r - 1.0_R8P)
      endif
   else
      weight = bspline_weight(order=order, r=r)
   endif
   endfunction effective_bspline_weight

   pure function gaussian_weight(r, cutoff_limit) result(weight)
   !< Return an unnormalized one-dimensional Gaussian weight with compact cutoff.
   implicit none
   real(R8P), intent(in) :: r
   real(R8P), intent(in) :: cutoff_limit
   real(R8P)             :: weight

   if (abs(r) <= cutoff_limit) then
      weight = exp(-0.5_R8P*r*r)
   else
      weight = 0.0_R8P
   endif
   endfunction gaussian_weight

   pure function effective_gaussian_weight(r, shift, cutoff_limit, filter) result(weight)
   !< Return the optionally binomial-filtered one-dimensional Gaussian weight.
   implicit none
   real(R8P), intent(in) :: r
   real(R8P), intent(in) :: shift
   real(R8P), intent(in) :: cutoff_limit
   logical,   intent(in) :: filter
   real(R8P)             :: weight

   if (filter) then
      weight = 0.25_R8P*gaussian_weight(r=r + shift, cutoff_limit=cutoff_limit) &
             + 0.50_R8P*gaussian_weight(r=r,         cutoff_limit=cutoff_limit) &
             + 0.25_R8P*gaussian_weight(r=r - shift, cutoff_limit=cutoff_limit)
   else
      weight = gaussian_weight(r=r, cutoff_limit=cutoff_limit)
   endif
   endfunction effective_gaussian_weight

   pure function bspline_weight(order, r) result(weight)
   !< Return the centered cardinal B-spline weight.
   implicit none
   integer(I4P), intent(in) :: order
   real(R8P),    intent(in) :: r
   real(R8P)                :: weight
   real(R8P)                :: a

   a = abs(r)

   select case(order)
   case(0_I4P)
      ! Nearest Grid Point, with the same tie convention as the ceiling-based
      ! particle cell lookup: the owning cell is (-1/2, +1/2].
      if (r > -0.5_R8P .and. r <= 0.5_R8P) then
         weight = 1.0_R8P
      else
         weight = 0.0_R8P
      endif
   case(1_I4P)
      ! Linear B-spline: support |r| <= 1.
      if (a <= 1.0_R8P) then
         weight = 1.0_R8P - a
      else
         weight = 0.0_R8P
      endif
   case(2_I4P)
      ! Quadratic B-spline: support |r| <= 3/2.
      if (a <= 0.5_R8P) then
         weight = 0.75_R8P - a*a
      elseif (a <= 1.5_R8P) then
         weight = 0.5_R8P * (1.5_R8P - a)**2
      else
         weight = 0.0_R8P
      endif
   case(3_I4P)
      ! Cubic B-spline: support |r| <= 2.
      if (a <= 1.0_R8P) then
         weight = (4.0_R8P - 6.0_R8P*a*a + 3.0_R8P*a**3) / 6.0_R8P
      elseif (a <= 2.0_R8P) then
         weight = (2.0_R8P - a)**3 / 6.0_R8P
      else
         weight = 0.0_R8P
      endif
   case(4_I4P)
      ! Quartic B-spline: support |r| <= 5/2.
      if (a <= 0.5_R8P) then
         weight = ((2.5_R8P - a)**4                          &
                  - 5.0_R8P  * (1.5_R8P - a)**4              &
                  + 10.0_R8P * (0.5_R8P - a)**4) / 24.0_R8P
      elseif (a <= 1.5_R8P) then
         weight = ((2.5_R8P - a)**4                          &
                  - 5.0_R8P * (1.5_R8P - a)**4) / 24.0_R8P
      elseif (a <= 2.5_R8P) then
         weight = (2.5_R8P - a)**4 / 24.0_R8P
      else
         weight = 0.0_R8P
      endif
   case(5_I4P)
      ! Quintic B-spline: support |r| <= 3.
      if (a <= 1.0_R8P) then
         weight = ((3.0_R8P - a)**5                          &
                  - 6.0_R8P  * (2.0_R8P - a)**5              &
                  + 15.0_R8P * (1.0_R8P - a)**5) / 120.0_R8P
      elseif (a <= 2.0_R8P) then
         weight = ((3.0_R8P - a)**5                          &
                  - 6.0_R8P * (2.0_R8P - a)**5) / 120.0_R8P
      elseif (a <= 3.0_R8P) then
         weight = (3.0_R8P - a)**5 / 120.0_R8P
      else
         weight = 0.0_R8P
      endif
   case(6_I4P)
      ! Sextic B-spline: support |r| <= 7/2.
      if (a <= 0.5_R8P) then
         weight = ((3.5_R8P - a)**6                         &
                  - 7.0_R8P  * (2.5_R8P - a)**6             &
                  + 21.0_R8P * (1.5_R8P - a)**6             &
                  - 35.0_R8P * (0.5_R8P - a)**6) / 720.0_R8P
      elseif (a <= 1.5_R8P) then
         weight = ((3.5_R8P - a)**6                         &
                  - 7.0_R8P  * (2.5_R8P - a)**6             &
                  + 21.0_R8P * (1.5_R8P - a)**6) / 720.0_R8P
      elseif (a <= 2.5_R8P) then
         weight = ((3.5_R8P - a)**6                         &
                  - 7.0_R8P * (2.5_R8P - a)**6) / 720.0_R8P
      elseif (a <= 3.5_R8P) then
         weight = (3.5_R8P - a)**6 / 720.0_R8P
      else
         weight = 0.0_R8P
      endif
   case default
      weight = 0.0_R8P
   endselect
   endfunction bspline_weight
endmodule adam_prism_pic_object
