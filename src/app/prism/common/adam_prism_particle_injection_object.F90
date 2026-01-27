!< ADAM, PRISM Particle-in-Cell class definition, CPU backend.
module adam_prism_particle_injection_object
!< ADAM, PRISM Particle-in-Cell class definition, CPU backend.
! ADAM modules
use :: adam_mpih_object, only : mpih_object
use adam_field_object, only : field_object
! PRISM modules
use :: adam_prism_parameters
use :: adam_prism_pic_object, only: prism_pic_object
! third party modules
use :: finer, only : file_ini
use :: penf,  only : I4P, R8P, str

implicit none
private
public :: INI_SECTION_NAME
public :: prism_particle_injection_object

character(len=18), parameter :: INI_SECTION_NAME = 'particle_injection'

character(len=33), parameter :: UNIFORM_DOMAIN_SPACE_DISTRIBUTION            = 'Uniform_domain_space_distribution'
character(len=32), parameter :: UNIFORM_BOX_SPACE_DISTRIBUTION               = 'Uniform_boxes_space_distribution'
character(len=31), parameter :: UNIFORM_CELL_SPACE_DISTRIBUTION              = 'Uniform_cell_space_distribution'
character(len=29), parameter :: SPACE_RANDOM_NUMBER_GENERATOR                = 'Space_random_number_generator'
character(len=30), parameter :: SPACE_LAYERED_NUMBER_GENERATOR               = 'Space_layered_number_generator'
character(len=18), parameter :: UNIFORM_MAXWELLIAN_VELOCITY_DISTRIBUTION     = 'Uniform_Maxwellian'
character(len=22), parameter :: NON_UNIFORM_MAXWELLIAN_VELOCITY_DISTRIBUTION = 'Non_uniform_Maxwellian'
character(len=10), parameter :: PURE_DRIFT_VELOCITY_DISTRIBUTION             = 'Pure_drift'
character(len=32), parameter :: VELOCITY_RANDOM_NUMBER_GENERATOR             = 'Velocity_random_number_generator'
character(len=33), parameter :: VELOCITY_LAYERED_NUMBER_GENERATOR            = 'Velocity_layered_number_generator'


type :: prism_particle_injection_object
   type(mpih_object)        :: mpih										 !< MPI handler.
   character(len=99)        :: space_distribution					 !< Particle space distribution type.
	character(len=99)        :: space_random_number_generator	 !< Type of random number generator for space distribution
	real(R8P)					 :: box_number = 0.0_R8P				 !< Number of boxes in which ensure charge neutrality
	logical						 :: space_pairing = .false.			 !< Enable space pairing of particles
	character(len=99)        :: velocity_distribution				 !< Particle velocity distribution type.
	real(R8P)					 :: T=0.0_R8P								 !< Plasma temperature (uniform)
	real(R8P)					 :: T_x=0.0_R8P							 !< Plasma temperature along x (non-uniform)
	real(R8P)					 :: T_y=0.0_R8P							 !< Plasma temperature along y (non-uniform)
	real(R8P)					 :: T_z=0.0_R8P							 !< Plasma temperature along z (non-uniform)
	character(len=99)        :: velocity_random_number_generator !< Type of random number generator for space distribution
	logical			 			 :: velocity_pairing = .false.		 !< Enable space pairing of particles
	real(R8P)					 :: v_drift_x=0.0_R8P					 !< Plasma drift velocity along x
	real(R8P)					 :: v_drift_y=0.0_R8P					 !< Plasma drift_velocity along y
	real(R8P)					 :: v_drift_z=0.0_R8P					 !< Plasma drift velocity along z
	logical        			 :: v_av_correction = .false.			 !< Flag to correct the average v.

   !< Pointer (abstract) TBP.
   !procedure(particle_space_injection_interface),	  	 	 pass(self), pointer :: particle_space_injection 	 => null() !< Particle space injection.
	!procedure(particle_velocity_injection_interface), 	 	 pass(self), pointer :: particle_velocity_injection => null() !< Particle velocity injection.
	!procedure(space_random_number_generator_interface), 	 pass(self), pointer :: space_rand_num_generator 	 => null() !< Space random number generator interface
	!procedure(velocity_random_number_generator_interface), pass(self), pointer :: velocity_rand_num_generator => null() !< Space random number generator interface
contains
   procedure, pass(self) :: description                   !< Return pretty-printed object description.
   procedure, pass(self) :: initialize                    !< Initialize IC.
   procedure, pass(self) :: load_from_file                !< Load config from file.
endtype prism_particle_injection_object

contains
   function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(prism_particle_injection_object), intent(in) :: self             !< External fields.
   character(len=:), allocatable                   	:: desc             !< Description.
   character(len=1), parameter                     	:: NL=new_line('a') !< New line character.

   desc =       self%mpih%myrankstr//'Particle injection object description:'
   desc = desc//NL//self%mpih%myrankstr//'    	Space initial distribution: '//self%space_distribution
	if (self%space_distribution == UNIFORM_BOX_SPACE_DISTRIBUTION) then
		desc = desc//NL//self%mpih%myrankstr//'    	Number of boxes: '//trim(str(self%box_number))
	endif
	desc = desc//NL//self%mpih%myrankstr//'    	Space random number generator: '//self%space_random_number_generator
	desc = desc//NL//self%mpih%myrankstr//'    	Space pairing: '//trim(str(self%space_pairing))
	desc = desc//NL//self%mpih%myrankstr//'    	Velocity initial distribution: '//self%velocity_distribution
	if (self%velocity_distribution == UNIFORM_MAXWELLIAN_VELOCITY_DISTRIBUTION) then
		desc = desc//NL//self%mpih%myrankstr//'    	Plasma Temperature: '//trim(str(self%T))
	elseif (self%velocity_distribution == NON_UNIFORM_MAXWELLIAN_VELOCITY_DISTRIBUTION) then
		desc = desc//NL//self%mpih%myrankstr//'    	Plasma x Temperature: '//trim(str(self%T_x))
		desc = desc//NL//self%mpih%myrankstr//'    	Plasma y Temperature: '//trim(str(self%T_y))
		desc = desc//NL//self%mpih%myrankstr//'    	Plasma z Temperature: '//trim(str(self%T_z))
	endif
	if (self%velocity_distribution == UNIFORM_MAXWELLIAN_VELOCITY_DISTRIBUTION .or. &
		self%velocity_distribution == NON_UNIFORM_MAXWELLIAN_VELOCITY_DISTRIBUTION) then
		desc = desc//NL//self%mpih%myrankstr//'    	Velocity random number generator: '//trim(self%velocity_random_number_generator)
		desc = desc//NL//self%mpih%myrankstr//'    	Velocity pairing: '//trim(str(self%velocity_pairing))
	endif
	desc = desc//NL//self%mpih%myrankstr//'    	Velocity averaging: '//trim(str(self%v_av_correction))
   endfunction description

   subroutine initialize(self, file_parameters)
   !< Initialize particle_injection.
   class(prism_particle_injection_object), intent(inout) :: self            !< External fields.
   type(file_ini),          intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'prism_particle_injection_object%initialize start'

   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%description()

	!select case(self%space_distribution)
   !case(UNIFORM_CELL_SPACE_DISTRIBUTION)
   !   self%particle_space_injection => uniform_cell_space_injection
   !case(UNIFORM_BOX_SPACE_DISTRIBUTION)
   !   self%particle_space_injection => uniform_box_space_injection
   !case(UNIFORM_DOMAIN_SPACE_DISTRIBUTION)
   !   self%particle_space_injection => uniform_domain_space_injection
   !case default
   !   call self%mpih%error_stop(msg=': invalid particle space injection model in prism_particle_injection_object%initialize')
   !endselect
	!select case(self%space_random_number_generator)
   !case(SPACE_RANDOM_NUMBER_GENERATOR)
   !   self%space_rand_num_generator => random_number_generator
   !case(SPACE_LAYERED_NUMBER_GENERATOR)
   !   self%space_rand_num_generator => layered_number_generator
   !case default
   !   call self%mpih%error_stop(msg=': invalid particle space random number generator in prism_particle_injection_object%initialize')
   !endselect
	!select case(self%velocity_distribution)
   !case(UNIFORM_MAXWELLIAN_VELOCITY_DISTRIBUTION)
   !   self%particle_velocity_injection => uniform_maxwellian_velocity_injection
   !case(NON_UNIFORM_MAXWELLIAN_VELOCITY_DISTRIBUTION)
   !   self%particle_velocity_injection => non_uniform_maxwellian_velocity_injection
   !!case(PURE_DRIFT_VELOCITY_DISTRIBUTION)
   !   !NON te lo scordare, ci vorrà un if nell'inizializzazione se non hai distribuzioni complesse (userai la funzione add drift velocity che scriverai)
   !case default
   !   call self%mpih%error_stop(msg=': invalid particle velocity injection model in prism_particle_injection_object%initialize')
   !endselect
	!select case(self%velocity_random_number_generator)
   !case(VELOCITY_RANDOM_NUMBER_GENERATOR)
   !   self%velocity_rand_num_generator => random_number_generator
   !case(VELOCITY_LAYERED_NUMBER_GENERATOR)
   !   self%velocity_rand_num_generator => layered_number_generator
   !case default
   !   call self%mpih%error_stop(msg=': invalid particle space random number generator in prism_particle_injection_object%initialize')
   !endselect
   print '(A)', self%mpih%myrankstr//'prism_particle_injection_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load PIC configuration from file.
	class(prism_particle_injection_object), intent(inout)   			 :: self             !< PIC object.
	type(file_ini),          intent(in)		  			 :: file_parameters  !< File handler.
   logical,                 intent(in), optional    :: go_on_fail      	!< Go on if load fails.
   logical                                          :: go_on_fail_     	!< Go on if load fails.
   integer(I4P)                                     :: error           	!< Error status.
   character(99)                                    :: buff       		!< Option character buffer.

	go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

	call file_parameters%get(section_name=INI_SECTION_NAME, option_name='space_distribution', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(space_distribution) from file')
   select case(trim(adjustl(buff)))
   case('Domain_uniform', 'domain_uniform', 'domain_Uniform')
      self%space_distribution = UNIFORM_DOMAIN_SPACE_DISTRIBUTION
   case('Box_uniform', 'box_uniform', 'box_Uniform')
      self%space_distribution = UNIFORM_BOX_SPACE_DISTRIBUTION
	case('Cell_uniform', 'cell_uniform', 'cell_Uniform')
      self%space_distribution = UNIFORM_CELL_SPACE_DISTRIBUTION
	case default
		call self%mpih%error_stop(msg=': invalid particle space distribution ['//trim(adjustl(buff))//'] in  & 
      ['//INI_SECTION_NAME//'].(space_distribution)')
	endselect

	call file_parameters%get(section_name=INI_SECTION_NAME, &
									 option_name='space_random_number_generator', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call self%mpih%error_stop(msg=': failed to load [' & 
												//INI_SECTION_NAME//'].(space_random_number_generator) from file')
   select case(trim(adjustl(buff)))
   case('Random', 'random', 'RANDOM')
      self%space_random_number_generator = SPACE_RANDOM_NUMBER_GENERATOR
   case('Layered', 'layered', 'LAYERED')
      self%space_random_number_generator = SPACE_LAYERED_NUMBER_GENERATOR
	case default
		call self%mpih%error_stop(msg=': invalid space random number generator ['//trim(adjustl(buff))//'] in  & 
      ['//INI_SECTION_NAME//'].(space_random_number_generator)')
	endselect

	if (self%space_distribution == UNIFORM_BOX_SPACE_DISTRIBUTION) then
		call file_parameters%get(section_name=INI_SECTION_NAME, option_name='box_number', &
   	val=self%box_number, error=error)
   	if (.not.go_on_fail_.and.error>0) & 
   	call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(box_number)')
	endif

	call file_parameters%get(section_name=INI_SECTION_NAME, option_name='space_pairing', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(space_pairing)')
   select case(trim(adjustl(buff)))
   case('NO', 'no', 'No', 'nO')
      self%space_pairing = .false.
	case('YES', 'yes', 'Yes')
		self%space_pairing = .true.
	case default
		call self%mpih%error_stop(msg=': invalid space pairing flag ['//trim(adjustl(buff))//'] in  & 
      ['//INI_SECTION_NAME//'].(space_pairing)')
	endselect

	call file_parameters%get(section_name=INI_SECTION_NAME, option_name='velocity_distribution', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(velocity_distribution) from file')
   select case(trim(adjustl(buff)))
   case('Uniform_Maxwellian', 'uniform_maxwellian', 'uniform_Maxwellian', & 
			'Maxwellian', 'maxwellian', 'Uniform', 'uniform')
      self%velocity_distribution = UNIFORM_MAXWELLIAN_VELOCITY_DISTRIBUTION
   case('Non_Uniform_Maxwellian', 'non_uniform_maxwellian', 'non_uniform_Maxwellian')
      self%velocity_distribution = NON_UNIFORM_MAXWELLIAN_VELOCITY_DISTRIBUTION
	case('Pure_drift', 'pure_drift', 'Drift', 'drift')
		self%velocity_distribution = PURE_DRIFT_VELOCITY_DISTRIBUTION
	case default
		call self%mpih%error_stop(msg=': invalid particle velocity distribution ['//trim(adjustl(buff))//'] in  & 
      ['//INI_SECTION_NAME//'].(velocity_distribution)')
	endselect

	if (self%velocity_distribution == UNIFORM_MAXWELLIAN_VELOCITY_DISTRIBUTION) then
		call file_parameters%get(section_name=INI_SECTION_NAME, option_name='Temperature', &
   	val=self%T, error=error)
   	if (.not.go_on_fail_.and.error>0) & 
   	call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(Temperature)')

	elseif (self%velocity_distribution == NON_UNIFORM_MAXWELLIAN_VELOCITY_DISTRIBUTION) then
		call file_parameters%get(section_name=INI_SECTION_NAME, option_name='Temperature_x', &
   	val=self%T_x, error=error)
   	if (.not.go_on_fail_.and.error>0) & 
   	call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(Temperature_x)')

		call file_parameters%get(section_name=INI_SECTION_NAME, option_name='Temperature_y', &
   	val=self%T_y, error=error)
   	if (.not.go_on_fail_.and.error>0) & 
   	call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(Temperature_y)')

		call file_parameters%get(section_name=INI_SECTION_NAME, option_name='Temperature_z', &
   	val=self%T_z, error=error)
   	if (.not.go_on_fail_.and.error>0) & 
   	call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(Temperature_z)')	
	endif

	if (self%velocity_distribution == UNIFORM_MAXWELLIAN_VELOCITY_DISTRIBUTION .or. &
		self%velocity_distribution == NON_UNIFORM_MAXWELLIAN_VELOCITY_DISTRIBUTION) then

		call file_parameters%get(section_name=INI_SECTION_NAME, &
										 option_name='velocity_random_number_generator', val=buff,error=error)
   	if (.not.go_on_fail_.and.error>0) &
   	   call self%mpih%error_stop(msg=': failed to load [' & 
													//INI_SECTION_NAME//'].(velocity_random_number_generator) from file')
   	select case(trim(adjustl(buff)))
   	case('Random', 'random', 'RANDOM')
   	   self%velocity_random_number_generator = VELOCITY_RANDOM_NUMBER_GENERATOR
   	case('Layered', 'layered', 'LAYERED')
   	   self%velocity_random_number_generator = VELOCITY_LAYERED_NUMBER_GENERATOR
		case default
			call self%mpih%error_stop(msg=': invalid velocity random number generator ['//trim(adjustl(buff))//'] in  & 
   	   ['//INI_SECTION_NAME//'].(velocity_random_number_generator)')
		endselect

		call file_parameters%get(section_name=INI_SECTION_NAME, option_name='velocity_pairing', val=buff,error=error)
   	if (.not.go_on_fail_.and.error>0) &
   	   call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(velocity_pairing)')
   	select case(trim(adjustl(buff)))
   	case('NO', 'no', 'No', 'nO')
   	   self%velocity_pairing = .false.
		case('YES', 'yes', 'Yes')
			self%velocity_pairing = .true.
		case default
			call self%mpih%error_stop(msg=': invalid velocity pairing flag ['//trim(adjustl(buff))//'] in  & 
   	   ['//INI_SECTION_NAME//'].(velocity_pairing)')
		endselect
	endif

	call file_parameters%get(section_name=INI_SECTION_NAME, option_name='v_drift_x', &
   		val=self%v_drift_x, error=error)
   if (.not.go_on_fail_.and.error>0) & 
   	call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(v_drift_x)')

	call file_parameters%get(section_name=INI_SECTION_NAME, option_name='v_drift_y', &
   		val=self%v_drift_y, error=error)
   if (.not.go_on_fail_.and.error>0) & 
   	call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(v_drift_y)')

	call file_parameters%get(section_name=INI_SECTION_NAME, option_name='v_drift_z', &
   		val=self%v_drift_z, error=error)
   if (.not.go_on_fail_.and.error>0) & 
   	call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(v_drift_z)')

	call file_parameters%get(section_name=INI_SECTION_NAME, option_name='v_av_correction', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(v_av_correction)')
   select case(trim(adjustl(buff)))
   case('NO', 'no', 'No', 'nO')
      self%v_av_correction = .false.
	case('YES', 'yes', 'Yes')
		self%v_av_correction = .true.
	case default
		call self%mpih%error_stop(msg=': invalid velocity average correction flag ['//trim(adjustl(buff))//'] in  & 
      ['//INI_SECTION_NAME//'].(v_av_correction)')
	endselect
   endsubroutine load_from_file

   endmodule adam_prism_particle_injection_object