!< ADAM, PRISM Particle-in-Cell class definition, CPU backend.
module adam_prism_particle_injection_object
!< ADAM, PRISM Particle-in-Cell class definition, CPU backend.

! ADAM classes, libraries, parameters
use :: adam_field_object, only : field_object
! ADAM singleton objects
use :: adam_mpih_global, only : mpih
use :: adam_grid_object, only : grid_object
! PRISM modules
use :: adam_prism_parameters
use :: adam_prism_pic_object, only : prism_pic_object, PLASMA_TYPE_PROBLEM, SINGLE_PARTICLE_TYPE_PROBLEM
! third party modules
use :: finer, only : file_ini
use :: penf,  only : I4P, R8P, str

implicit none
private
public :: INI_SECTION_NAME
public :: prism_particle_injection_object
public :: write_initial_injection_tab

character(len=18), parameter :: INI_SECTION_NAME = 'particle_injection'

character(len=33), parameter :: UNIFORM_DOMAIN_SPACE_DISTRIBUTION            = 'Uniform_domain_space_distribution'
character(len=32), parameter :: UNIFORM_BOX_SPACE_DISTRIBUTION               = 'Uniform_boxes_space_distribution'
character(len=31), parameter :: UNIFORM_CELL_SPACE_DISTRIBUTION              = 'Uniform_cell_space_distribution'
character(len=29), parameter :: SPACE_RANDOM_NUMBER_GENERATOR                = 'Space_random_number_generator'
character(len=30), parameter :: SPACE_LAYERED_NUMBER_GENERATOR               = 'Space_layered_number_generator'
character(len=18), parameter :: UNIFORM_MAXWELLIAN_VELOCITY_DISTRIBUTION     = 'Uniform_Maxwellian'
character(len=22), parameter :: NON_UNIFORM_MAXWELLIAN_VELOCITY_DISTRIBUTION = 'Non_uniform_Maxwellian'
character(len=32), parameter :: VELOCITY_RANDOM_NUMBER_GENERATOR             = 'Velocity_random_number_generator'
character(len=33), parameter :: VELOCITY_LAYERED_NUMBER_GENERATOR            = 'Velocity_layered_number_generator'

procedure(space_random_number_generator_interface),    pointer :: space_rand_num_generator 	  => null()
                                                                                                        !< Space random number
                                                                                                        !< generator interface
procedure(velocity_random_number_generator_interface), pointer :: velocity_rand_num_generator  => null()
                                                                                                         !< Space random number
                                                                                                         !< generator interface

type :: prism_particle_injection_object
   character(len=99)    :: space_distribution					!< Particle space distribution type.
	character(len=99)    :: space_random_number_generator		!< Type of random number generator for space distribution
	real(R8P)				:: box_number = 0.0_R8P					!< Number of boxes in which ensure charge neutrality
	logical					:: space_pairing = .false.				!< Enable space pairing of particles
	character(len=99)    :: velocity_distribution				!< Particle velocity distribution type.
	real(R8P)				:: T_i=0.0_R8P								!< Ionic plasma temperature (uniform)
	real(R8P)				:: T_i_x=0.0_R8P							!< Ionic plasma temperature along x (non-uniform)
	real(R8P)				:: T_i_y=0.0_R8P							!< Ionic plasma temperature along y (non-uniform)
	real(R8P)				:: T_i_z=0.0_R8P							!< Ionic plasma temperature along z (non-uniform)
	real(R8P)				:: T_e=0.0_R8P								!< Electronic plasma temperature (uniform)
	real(R8P)				:: T_e_x=0.0_R8P							!< Electronic plasma temperature along x (non-uniform)
	real(R8P)				:: T_e_y=0.0_R8P							!< Electronic plasma temperature along y (non-uniform)
	real(R8P)				:: T_e_z=0.0_R8P							!< Electronic plasma temperature along z (non-uniform)
	real(R8P)				:: T_n=0.0_R8P								!< Neutrals plasma temperature (uniform)
	real(R8P)				:: T_n_x=0.0_R8P							!< Neutrals plasma temperature along x (non-uniform)
	real(R8P)				:: T_n_y=0.0_R8P							!< Neutrals plasma temperature along y (non-uniform)
	real(R8P)				:: T_n_z=0.0_R8P							!< Neutrals plasma temperature along z (non-uniform)
	real(R8P)				:: x_position=0.0_R8P					!< x coordinate of the initial position of the particle in single particle problem
	real(R8P)				:: y_position=0.0_R8P					!< y coordinate of the initial position of the particle in single particle problem
	real(R8P)				:: z_position=0.0_R8P					!< z coordinate of the initial position of the particle in single particle problem
	real(R8P)				:: charge=0.0_R8P							!< charge of the particle in single particle problem
	real(R8P)				:: mass=0.0_R8P							!< mass of the particle in single particle problem
	character(len=99)    :: velocity_random_number_generator !< Type of random number generator for space distribution
	logical			 		:: velocity_pairing = .false.			!< Enable space pairing of particles
	real(R8P)				:: v_drift_x=0.0_R8P						!< Plasma drift velocity along x
	real(R8P)				:: v_drift_y=0.0_R8P						!< Plasma drift_velocity along y
	real(R8P)				:: v_drift_z=0.0_R8P						!< Plasma drift velocity along z
	logical        		:: v_av_correction = .false.			!< Flag to correct the average v.

   !< Pointer (abstract) TBP.
   procedure(particle_space_injection_interface),	  pass(self), pointer :: particle_space_injection 	  => null()
                                                                                                                 !< Particle space
                                                                                                                 !< injection.
	procedure(particle_velocity_injection_interface), pass(self), pointer :: particle_velocity_injection => null()
                                                                                                                !< Particle velocity
                                                                                                                !< injection.
contains
   procedure, pass(self) :: description    !< Return pretty-printed object description.
   procedure, pass(self) :: initialize     !< Initialize IC.
   procedure, pass(self) :: load_from_file !< Load config from file.
	procedure, pass(self) :: set_particle_initial_injection
	procedure, pass(self) :: uniform_domain_space_injection
	procedure, pass(self) :: uniform_cell_space_injection
	procedure, pass(self) :: uniform_maxwellian_velocity_injection
	procedure, pass(self) :: non_uniform_maxwellian_velocity_injection
	procedure, pass(self) :: single_particle_injection
endtype prism_particle_injection_object

interface
   subroutine particle_space_injection_interface(self, field, grid, pic, q_pic)
   import :: prism_particle_injection_object, field_object, grid_object, prism_pic_object, I4P, R8P
	class(prism_particle_injection_object), intent(inout) :: self
	type(field_object),                  	 intent(in) 	:: field
	type(grid_object),                   	 intent(in)	:: grid
	type(prism_pic_object),					 	 intent(in)		:: pic
	real(R8P),                           	 intent(inout) :: q_pic(1:,1:)
                                                                                                                              !< Num
                                                                                                                              !< ber
                                                                                                                              !< of
                                                                                                                              !< var
                                                                                                                              !< iab
                                                                                                                              !< les
                                                                                                                              !< .
   endsubroutine particle_space_injection_interface

	subroutine space_random_number_generator_interface(N, shuffled_list, i_numb, r_n)
	import :: I4P, R8P
	integer(I4P), intent(inout) :: shuffled_list(1:,1:)
	integer(I4P), intent(in) 	 :: i_numb
	integer(I4P), intent(in) 	 :: N
	real(R8P), intent(inout) 	 :: r_n(1:)
	endsubroutine space_random_number_generator_interface

	subroutine particle_velocity_injection_interface(self, field, grid, pic, q_pic)
   import :: prism_particle_injection_object, field_object, grid_object, prism_pic_object, I4P, R8P
	class(prism_particle_injection_object), intent(inout) :: self
	type(field_object),                  	 intent(in) 	:: field
	type(grid_object),                   	 intent(in)	:: grid
	type(prism_pic_object),					 	 intent(in)		:: pic
	real(R8P),                           	 intent(inout) :: q_pic(1:,1:)
                                                                                                                              !< Num
                                                                                                                              !< ber
                                                                                                                              !< of
                                                                                                                              !< var
                                                                                                                              !< iab
                                                                                                                              !< les
                                                                                                                              !< .
   endsubroutine particle_velocity_injection_interface

	subroutine velocity_random_number_generator_interface(N, shuffled_list, i_numb, r_n)
	import :: I4P, R8P
	integer(I4P), intent(inout) :: shuffled_list(1:,1:)
	integer(I4P), intent(in) 	 :: i_numb
	integer(I4P), intent(in) 	 :: N
	real(R8P), intent(inout) 	 :: r_n(1:)
	endsubroutine velocity_random_number_generator_interface
endinterface

contains
   function description(self, pic) result(desc)
   !< Return a pretty-formatted object description.
   class(prism_particle_injection_object), intent(in) :: self             !< External fields.
	type(prism_pic_object), 					 intent(in) :: pic				  !< PIC object
   character(len=:), allocatable                   	:: desc             !< Description.
   character(len=1), parameter                     	:: NL=new_line('a') !< New line character.

   desc =       mpih%myrankstr//'Particle injection object description:'
	if (pic%problem_type == PLASMA_TYPE_PROBLEM) then
   	desc = desc//NL//mpih%myrankstr//'    	Space initial distribution: '//self%space_distribution
		if (self%space_distribution == UNIFORM_BOX_SPACE_DISTRIBUTION) then
			desc = desc//NL//mpih%myrankstr//'    	Number of boxes: '//trim(str(self%box_number))
		endif
		desc = desc//NL//mpih%myrankstr//'    	Space random number generator: '//self%space_random_number_generator
		desc = desc//NL//mpih%myrankstr//'    	Space pairing: '//trim(str(self%space_pairing))
		desc = desc//NL//mpih%myrankstr//'    	Velocity initial distribution: '//self%velocity_distribution
		if (self%velocity_distribution == UNIFORM_MAXWELLIAN_VELOCITY_DISTRIBUTION) then
			desc = desc//NL//mpih%myrankstr//'    	Plasma  Ionic Temperature: '//trim(str(self%T_i))
			desc = desc//NL//mpih%myrankstr//'    	Plasma  Electronic Temperature: '//trim(str(self%T_e))
			desc = desc//NL//mpih%myrankstr//'    	Plasma  Neutrals Temperature: '//trim(str(self%T_n))
		elseif (self%velocity_distribution == NON_UNIFORM_MAXWELLIAN_VELOCITY_DISTRIBUTION) then
			desc = desc//NL//mpih%myrankstr//'    	Plasma x Ionic Temperature: '//trim(str(self%T_i_x))
			desc = desc//NL//mpih%myrankstr//'    	Plasma y Ionic Temperature: '//trim(str(self%T_i_y))
			desc = desc//NL//mpih%myrankstr//'    	Plasma z Ionic Temperature: '//trim(str(self%T_i_z))
			desc = desc//NL//mpih%myrankstr//'    	Plasma x Electronic Temperature: '//trim(str(self%T_e_x))
			desc = desc//NL//mpih%myrankstr//'    	Plasma y Electronic Temperature: '//trim(str(self%T_e_y))
			desc = desc//NL//mpih%myrankstr//'    	Plasma z Electronic Temperature: '//trim(str(self%T_e_z))
			desc = desc//NL//mpih%myrankstr//'    	Plasma x Neutrals Temperature: '//trim(str(self%T_n_x))
			desc = desc//NL//mpih%myrankstr//'    	Plasma y Neutrals Temperature: '//trim(str(self%T_n_y))
			desc = desc//NL//mpih%myrankstr//'    	Plasma z Neutrals Temperature: '//trim(str(self%T_n_z))
		endif
		if (self%velocity_distribution == UNIFORM_MAXWELLIAN_VELOCITY_DISTRIBUTION .or. &
			self%velocity_distribution == NON_UNIFORM_MAXWELLIAN_VELOCITY_DISTRIBUTION) then
			desc = desc//NL//mpih%myrankstr//'    	Velocity random number generator: '//trim(self%velocity_random_number_generator)
			desc = desc//NL//mpih%myrankstr//'    	Velocity pairing: '//trim(str(self%velocity_pairing))
		endif
		desc = desc//NL//mpih%myrankstr//'    	Velocity averaging: '//trim(str(self%v_av_correction))
	elseif (pic%problem_type == SINGLE_PARTICLE_TYPE_PROBLEM) then
		desc = desc//NL//mpih%myrankstr//'    	Particle x coordinate initial position: '//trim(str(self%x_position))
		desc = desc//NL//mpih%myrankstr//'    	Particle y coordinate initial position: '//trim(str(self%y_position))
		desc = desc//NL//mpih%myrankstr//'    	Particle z coordinate initial position: '//trim(str(self%z_position))
		desc = desc//NL//mpih%myrankstr//'    	Particle charge: '//trim(str(self%charge))
		desc = desc//NL//mpih%myrankstr//'    	Particle mass: '//trim(str(self%mass))
	endif
	desc = desc//NL//mpih%myrankstr//'    	Drift velocity along x: '//trim(str(self%v_drift_x))
	desc = desc//NL//mpih%myrankstr//'    	Drift velocity along y: '//trim(str(self%v_drift_y))
	desc = desc//NL//mpih%myrankstr//'    	Drift velocity along z: '//trim(str(self%v_drift_z))

   endfunction description

   subroutine initialize(self, file_parameters, pic)
   !< Initialize particle_injection.
   class(prism_particle_injection_object), intent(inout) :: self            !< External fields.
   type(file_ini),          					 intent(in)    :: file_parameters !< Simulation parameters ini file handler.
	type(prism_pic_object), 					 intent(in)		:: pic				 !< Pic object

   print '(A)', mpih%myrankstr//'prism_particle_injection_object%initialize start'

   call self%load_from_file(file_parameters=file_parameters, pic=pic)
   print '(A)', self%description(pic = pic)

	if (pic%problem_type == PLASMA_TYPE_PROBLEM) then
		select case(self%space_distribution)
   	case(UNIFORM_CELL_SPACE_DISTRIBUTION)
   	   self%particle_space_injection => uniform_cell_space_injection
   	!case(UNIFORM_BOX_SPACE_DISTRIBUTION)
   	!   self%particle_space_injection => uniform_box_space_injection
   	case(UNIFORM_DOMAIN_SPACE_DISTRIBUTION)
   	   self%particle_space_injection => uniform_domain_space_injection
   	case default
   	   call mpih%error_stop &
			(msg=': invalid particle space injection model in prism_particle_injection_object%initialize')
   	endselect
		select case(self%space_random_number_generator)
   	case(SPACE_RANDOM_NUMBER_GENERATOR)
   	   space_rand_num_generator => random_number_generator
   	case(SPACE_LAYERED_NUMBER_GENERATOR)
   	   space_rand_num_generator => layered_number_generator
   	case default
   	   call mpih%error_stop &
			(msg=': invalid particle space random number generator in prism_particle_injection_object%initialize')
   	endselect
		select case(self%velocity_distribution)
   	case(UNIFORM_MAXWELLIAN_VELOCITY_DISTRIBUTION)
   	   self%particle_velocity_injection => uniform_maxwellian_velocity_injection
   	case(NON_UNIFORM_MAXWELLIAN_VELOCITY_DISTRIBUTION)
   	   self%particle_velocity_injection => non_uniform_maxwellian_velocity_injection
   	case default
   	   call mpih%error_stop(msg=': invalid particle velocity injection model in prism_particle_injection_object%initialize')
   	endselect
		select case(self%velocity_random_number_generator)
   	case(VELOCITY_RANDOM_NUMBER_GENERATOR)
   	   velocity_rand_num_generator => random_number_generator
   	case(VELOCITY_LAYERED_NUMBER_GENERATOR)
   	   velocity_rand_num_generator => layered_number_generator
   	case default
   	   call mpih%error_stop &
			(msg=': invalid particle space random number generator in prism_particle_injection_object%initialize')
   	endselect
   	print '(A)', mpih%myrankstr//'prism_particle_injection_object%initialize finish'
	endif
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, pic, go_on_fail)
   !< Load PIC configuration from file.
	class(prism_particle_injection_object), intent(inout)   		 :: self             !< Particle injection object.
	type(file_ini),          					 intent(in)		  		 :: file_parameters  !< File handler.
	type(prism_pic_object), 					 intent(in)				 :: pic					!< PIC object
   logical,                 					 intent(in), optional :: go_on_fail      	!< Go on if load fails.
   logical                                          				 :: go_on_fail_     	!< Go on if load fails.
   integer(I4P)                                     				 :: error           	!< Error status.
   character(99)                                    				 :: buff       		!< Option character buffer.

	go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail
	if (pic%problem_type == PLASMA_TYPE_PROBLEM) then
		call file_parameters%get(section_name=INI_SECTION_NAME, option_name='space_distribution', val=buff,error=error)
   	if (.not.go_on_fail_.and.error>0) &
   	   call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(space_distribution) from file')
   	select case(trim(adjustl(buff)))
   	case('Domain_uniform', 'domain_uniform', 'domain_Uniform')
   	   self%space_distribution = UNIFORM_DOMAIN_SPACE_DISTRIBUTION
   	case('Box_uniform', 'box_uniform', 'box_Uniform')
   	   self%space_distribution = UNIFORM_BOX_SPACE_DISTRIBUTION
		case('Cell_uniform', 'cell_uniform', 'cell_Uniform')
   	   self%space_distribution = UNIFORM_CELL_SPACE_DISTRIBUTION
		case default
			call mpih%error_stop(msg=': invalid particle space distribution ['//trim(adjustl(buff))//'] in  &
   	   ['//INI_SECTION_NAME//'].(space_distribution)')
		endselect

		call file_parameters%get(section_name=INI_SECTION_NAME, &
										 option_name='space_random_number_generator', val=buff,error=error)
   	if (.not.go_on_fail_.and.error>0) &
   	   call mpih%error_stop(msg=': failed to load [' &
													//INI_SECTION_NAME//'].(space_random_number_generator) from file')
   	select case(trim(adjustl(buff)))
   	case('Random', 'random', 'RANDOM')
   	   self%space_random_number_generator = SPACE_RANDOM_NUMBER_GENERATOR
   	case('Layered', 'layered', 'LAYERED')
   	   self%space_random_number_generator = SPACE_LAYERED_NUMBER_GENERATOR
		case default
			call mpih%error_stop(msg=': invalid space random number generator ['//trim(adjustl(buff))//'] in  &
   	   ['//INI_SECTION_NAME//'].(space_random_number_generator)')
		endselect

		if (self%space_distribution == UNIFORM_BOX_SPACE_DISTRIBUTION) then
			call file_parameters%get(section_name=INI_SECTION_NAME, option_name='box_number', &
   		val=self%box_number, error=error)
   		if (.not.go_on_fail_.and.error>0) &
   		call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(box_number)')
		endif

		call file_parameters%get(section_name=INI_SECTION_NAME, option_name='space_pairing', val=buff,error=error)
   	if (.not.go_on_fail_.and.error>0) &
   	   call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(space_pairing)')
   	select case(trim(adjustl(buff)))
   	case('NO', 'no', 'No', 'nO')
   	   self%space_pairing = .false.
		case('YES', 'yes', 'Yes')
			self%space_pairing = .true.
		case default
			call mpih%error_stop(msg=': invalid space pairing flag ['//trim(adjustl(buff))//'] in  &
   	   ['//INI_SECTION_NAME//'].(space_pairing)')
		endselect

		call file_parameters%get(section_name=INI_SECTION_NAME, option_name='velocity_distribution', val=buff,error=error)
   	if (.not.go_on_fail_.and.error>0) &
   	   call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(velocity_distribution) from file')
   	select case(trim(adjustl(buff)))
   	case('Uniform_Maxwellian', 'uniform_maxwellian', 'uniform_Maxwellian', &
				'Maxwellian', 'maxwellian', 'Uniform', 'uniform')
   	   self%velocity_distribution = UNIFORM_MAXWELLIAN_VELOCITY_DISTRIBUTION
   	case('Non_Uniform_Maxwellian', 'non_uniform_maxwellian', 'non_uniform_Maxwellian')
   	   self%velocity_distribution = NON_UNIFORM_MAXWELLIAN_VELOCITY_DISTRIBUTION
		case default
			call mpih%error_stop(msg=': invalid particle velocity distribution ['//trim(adjustl(buff))//'] in  &
   	   ['//INI_SECTION_NAME//'].(velocity_distribution)')
		endselect

		if (self%velocity_distribution == UNIFORM_MAXWELLIAN_VELOCITY_DISTRIBUTION) then
			call file_parameters%get(section_name=INI_SECTION_NAME, option_name='Ionic_temperature', &
   		val=self%T_i, error=error)
   		if (.not.go_on_fail_.and.error>0) &
   		call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(Ionic_temperature)')

			call file_parameters%get(section_name=INI_SECTION_NAME, option_name='Electronic_temperature', &
   		val=self%T_e, error=error)
   		if (.not.go_on_fail_.and.error>0) &
   		call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(Electronic_temperature)')

			call file_parameters%get(section_name=INI_SECTION_NAME, option_name='Neutrals_temperature', &
   		val=self%T_n, error=error)
   		if (.not.go_on_fail_.and.error>0) &
   		call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(Neutrals_temperature)')

		elseif (self%velocity_distribution == NON_UNIFORM_MAXWELLIAN_VELOCITY_DISTRIBUTION) then
			call file_parameters%get(section_name=INI_SECTION_NAME, option_name='Ionic_temperature_x', &
   		val=self%T_i_x, error=error)
   		if (.not.go_on_fail_.and.error>0) &
   		call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(Ionic_temperature_x)')

			call file_parameters%get(section_name=INI_SECTION_NAME, option_name='Ionic_temperature_y', &
   		val=self%T_i_y, error=error)
   		if (.not.go_on_fail_.and.error>0) &
   		call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(Ionic_temperature_y)')

			call file_parameters%get(section_name=INI_SECTION_NAME, option_name='Ionic_temperature_z', &
   		val=self%T_i_z, error=error)
   		if (.not.go_on_fail_.and.error>0) &
   		call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(Ionic_temperature_z)')

			call file_parameters%get(section_name=INI_SECTION_NAME, option_name='Electronic_temperature_x', &
   		val=self%T_e_x, error=error)
   		if (.not.go_on_fail_.and.error>0) &
   		call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(Electronic_temperature_x)')

			call file_parameters%get(section_name=INI_SECTION_NAME, option_name='Electronic_temperature_y', &
   		val=self%T_e_y, error=error)
   		if (.not.go_on_fail_.and.error>0) &
   		call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(Electronic_temperature_y)')

			call file_parameters%get(section_name=INI_SECTION_NAME, option_name='Electronic_temperature_z', &
   		val=self%T_e_z, error=error)
   		if (.not.go_on_fail_.and.error>0) &
   		call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(Electronic_temperature_z)')

			call file_parameters%get(section_name=INI_SECTION_NAME, option_name='Neutrals_temperature_x', &
   		val=self%T_n_x, error=error)
   		if (.not.go_on_fail_.and.error>0) &
   		call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(Neutrals_temperature_x)')

			call file_parameters%get(section_name=INI_SECTION_NAME, option_name='Neutrals_temperature_y', &
   		val=self%T_n_y, error=error)
   		if (.not.go_on_fail_.and.error>0) &
   		call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(Neutrals_temperature_y)')

			call file_parameters%get(section_name=INI_SECTION_NAME, option_name='Neutrals_temperature_z', &
   		val=self%T_n_z, error=error)
   		if (.not.go_on_fail_.and.error>0) &
   		call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(Neutrals_temperature_z)')
		endif

		if (self%velocity_distribution == UNIFORM_MAXWELLIAN_VELOCITY_DISTRIBUTION .or. &
			self%velocity_distribution == NON_UNIFORM_MAXWELLIAN_VELOCITY_DISTRIBUTION) then

			call file_parameters%get(section_name=INI_SECTION_NAME, &
											 option_name='velocity_random_number_generator', val=buff,error=error)
   		if (.not.go_on_fail_.and.error>0) &
   		   call mpih%error_stop(msg=': failed to load [' &
														//INI_SECTION_NAME//'].(velocity_random_number_generator) from file')
   		select case(trim(adjustl(buff)))
   		case('Random', 'random', 'RANDOM')
   		   self%velocity_random_number_generator = VELOCITY_RANDOM_NUMBER_GENERATOR
   		case('Layered', 'layered', 'LAYERED')
   		   self%velocity_random_number_generator = VELOCITY_LAYERED_NUMBER_GENERATOR
			case default
				call mpih%error_stop(msg=': invalid velocity random number generator ['//trim(adjustl(buff))//'] in  &
   		   ['//INI_SECTION_NAME//'].(velocity_random_number_generator)')
			endselect

			call file_parameters%get(section_name=INI_SECTION_NAME, option_name='velocity_pairing', val=buff,error=error)
   		if (.not.go_on_fail_.and.error>0) &
   		   call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(velocity_pairing)')
   		select case(trim(adjustl(buff)))
   		case('NO', 'no', 'No', 'nO')
   		   self%velocity_pairing = .false.
			case('YES', 'yes', 'Yes')
				self%velocity_pairing = .true.
			case default
				call mpih%error_stop(msg=': invalid velocity pairing flag ['//trim(adjustl(buff))//'] in  &
   		   ['//INI_SECTION_NAME//'].(velocity_pairing)')
			endselect
		endif

		call file_parameters%get(section_name=INI_SECTION_NAME, option_name='v_av_correction', val=buff,error=error)
   	if (.not.go_on_fail_.and.error>0) &
   	   call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(v_av_correction)')
   	select case(trim(adjustl(buff)))
   	case('NO', 'no', 'No', 'nO')
   	   self%v_av_correction = .false.
		case('YES', 'yes', 'Yes')
			self%v_av_correction = .true.
		case default
			call mpih%error_stop(msg=': invalid velocity average correction flag ['//trim(adjustl(buff))//'] in  &
   	   ['//INI_SECTION_NAME//'].(v_av_correction)')
		endselect
	endif

	if (pic%problem_type == SINGLE_PARTICLE_TYPE_PROBLEM) then
		call file_parameters%get(section_name=INI_SECTION_NAME, option_name='x_position', &
   		val=self%x_position, error=error)
   	if (.not.go_on_fail_.and.error>0) &
   	call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(x_position)')

		call file_parameters%get(section_name=INI_SECTION_NAME, option_name='y_position', &
   		val=self%y_position, error=error)
   	if (.not.go_on_fail_.and.error>0) &
   	call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(y_position)')

		call file_parameters%get(section_name=INI_SECTION_NAME, option_name='z_position', &
   		val=self%z_position, error=error)
   	if (.not.go_on_fail_.and.error>0) &
   	call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(z_position)')

		call file_parameters%get(section_name=INI_SECTION_NAME, option_name='charge', &
   		val=self%charge, error=error)
   	if (.not.go_on_fail_.and.error>0) &
   	call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(charge)')

		call file_parameters%get(section_name=INI_SECTION_NAME, option_name='mass', &
   		val=self%mass, error=error)
   	if (.not.go_on_fail_.and.error>0) &
   	call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(mass)')
	endif

	call file_parameters%get(section_name=INI_SECTION_NAME, option_name='v_drift_x', &
   		val=self%v_drift_x, error=error)
   if (.not.go_on_fail_.and.error>0) &
   	call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(v_drift_x)')

	call file_parameters%get(section_name=INI_SECTION_NAME, option_name='v_drift_y', &
   		val=self%v_drift_y, error=error)
   if (.not.go_on_fail_.and.error>0) &
   	call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(v_drift_y)')

	call file_parameters%get(section_name=INI_SECTION_NAME, option_name='v_drift_z', &
   		val=self%v_drift_z, error=error)
   if (.not.go_on_fail_.and.error>0) &
   	call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(v_drift_z)')

   endsubroutine load_from_file

	subroutine set_particle_initial_injection(self, field, grid, pic, q_pic)
	class(prism_particle_injection_object), intent(inout) :: self
	type(field_object),                  	 intent(in) 	:: field
	type(grid_object),                   	 intent(in)	:: grid
	type(prism_pic_object),					 	 intent(inout)	:: pic
	real(R8P),                           	 intent(inout) :: q_pic(1:,1:)

	if (pic%problem_type == PLASMA_TYPE_PROBLEM) then
		!Setta posizione spaziale delle particelle e relativa tipologia
		call self%particle_space_injection(field=field, grid=grid, pic=pic, q_pic=q_pic)
		!Setta velocità iniziale delle particelle
		call self%particle_velocity_injection(field=field, grid=grid, pic=pic, q_pic=q_pic)
		!Definisci neighbour list
		call pic%particle_cartesian_grid_index(field=field, grid=grid, q_pic=q_pic)

		!Aggiungo qui successivamente interpolazione iniziale dei campi e spalmatura particelle (cariche e correnti) su griglia

	elseif(pic%problem_type == SINGLE_PARTICLE_TYPE_PROBLEM) then
	!Setta posizione velocità e caratteristiche della particella
	call self%single_particle_injection(q_pic=q_pic)
	!Definisci neighbour list
	call pic%particle_cartesian_grid_index(field=field, grid=grid, q_pic=q_pic)

	!Aggiungo qui successivamente interpolazione iniziale dei campi e spalmatura particelle (cariche e correnti) su griglia

	endif
	endsubroutine

	subroutine single_particle_injection(self, q_pic)
	class(prism_particle_injection_object), intent(inout) :: self
	real(R8P),                           	 intent(inout) :: q_pic(1:,1:)

	associate(x_p=>self%x_position, y_p=>self%y_position, z_p=>self%z_position, vx_p=>self%v_drift_x, &
			 vy_p=>self%v_drift_y,  vz_p=>self%v_drift_z, charge=>self%charge, mass=>self%mass)

	q_pic(1,1) = x_p
	q_pic(2,1) = y_p
	q_pic(3,1) = z_p
	q_pic(4,1) = vx_p
	q_pic(5,1) = vy_p
	q_pic(6,1) = vz_p
	q_pic(7,1) = charge
	q_pic(8,1) = mass
	endassociate
	endsubroutine single_particle_injection

   subroutine uniform_domain_space_injection(self, field, grid, pic, q_pic)
	class(prism_particle_injection_object), intent(inout) :: self
	type(field_object),                  	 intent(in) 	:: field
	type(grid_object),                       intent(in) 	:: grid !< Grid (sibling realm component, threaded in).
	type(prism_pic_object),					 	 intent(in)		:: pic
	real(R8P),                           	 intent(inout) :: q_pic(1:,1:)
	real(R8P)															:: r_n(3)
	real(R8P)															:: domain_volume
	real(R8P)															:: x_p, y_p, z_p
	integer(I4P)														:: i
	integer(I4P), allocatable										:: shuffled_list_ions(:,:)
	integer(I4P), allocatable										:: shuffled_list_electrons(:,:)
	integer(I4P), allocatable										:: shuffled_list_neutrals(:,:)
	character(len=:), allocatable		                   	:: desc
   character(len=1), parameter  		                   	:: NL=new_line('a')

	associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj, nk=>grid%nk,              &
      		ngc=>grid%ngc, dx=>field%dxyz(1,:), dy=>field%dxyz(2,:), dz=>field%dxyz(3,:),          &
      		np=>pic%particle_number, e_min=>grid%domain_emin, e_max=>grid%domain_emax,  		 &
				neutral_fraction=>pic%neutral_fraction, n_ions=>pic%n_ions, n_electrons=>pic%n_electrons, 	 &
				n_neutrals=>pic%n_neutrals)

	allocate(shuffled_list_ions	  (1:3,1:n_ions))
	allocate(shuffled_list_electrons(1:3,1:n_electrons))
	allocate(shuffled_list_neutrals (1:3,1:n_neutrals))
	shuffled_list_ions(:,:) 	  = 0_I4P
	shuffled_list_electrons(:,:) = 0_I4P
	shuffled_list_neutrals(:,:)  = 0_I4P

	!if(.not.space_pairing) then
		do i = 1, n_ions
			call space_rand_num_generator(N=n_ions, shuffled_list=shuffled_list_ions, i_numb=i, r_n=r_n)
			x_p = e_min(1)+r_n(1)*(e_max(1)-e_min(1))
			y_p = e_min(2)+r_n(2)*(e_max(2)-e_min(2))
			z_p = e_min(3)+r_n(3)*(e_max(3)-e_min(3))
			q_pic(1,i) = x_p
			q_pic(2,i) = y_p
			q_pic(3,i) = z_p
			q_pic(7,i) = -E_CHARGE
			q_pic(8,i) = E_MASS
			!Z0 = sqrt(-2.0_R8P*log(r_n(1)))*cos(2*PI*r_n(2))
			!Z1 = sqrt(-2.0_R8P*log(r_n(1)))*cos(2*PI*r_n(2))
		enddo
		do i = 1, n_electrons
			call space_rand_num_generator(N=n_electrons, shuffled_list=shuffled_list_electrons, i_numb=i, r_n=r_n)
			x_p = e_min(1)+r_n(1)*(e_max(1)-e_min(1))
			y_p = e_min(2)+r_n(2)*(e_max(2)-e_min(2))
			z_p = e_min(3)+r_n(3)*(e_max(3)-e_min(3))
			q_pic(1,i+n_ions) = x_p
			q_pic(2,i+n_ions) = y_p
			q_pic(3,i+n_ions) = z_p
			q_pic(7,i+n_ions) = E_CHARGE
			q_pic(8,i+n_ions) = E_MASS
		enddo
		do i = 1, n_neutrals
			call space_rand_num_generator(N=n_neutrals, shuffled_list=shuffled_list_neutrals, i_numb=i, r_n=r_n)
			x_p = e_min(1)+r_n(1)*(e_max(1)-e_min(1))
			y_p = e_min(2)+r_n(2)*(e_max(2)-e_min(2))
			z_p = e_min(3)+r_n(3)*(e_max(3)-e_min(3))
			q_pic(1,i+n_ions+n_electrons) = x_p
			q_pic(2,i+n_ions+n_electrons) = y_p
			q_pic(3,i+n_ions+n_electrons) = z_p
			q_pic(7,i+n_ions+n_electrons) = 0.0_R8P
			q_pic(8,i+n_ions+n_electrons) = E_MASS
		enddo
	!else

	!endif

	domain_volume = (e_max(1)-e_min(1))*(e_max(2)-e_min(2))*(e_max(3)-e_min(3))

	desc =       mpih%myrankstr//'Injected particles:'
   desc = desc//NL//mpih%myrankstr//'    	Electrons number: '//trim(str(n_electrons))
	desc = desc//NL//mpih%myrankstr//'    	Ions number: '//trim(str(n_ions))
	desc = desc//NL//mpih%myrankstr//'    	Neutrals number: '//trim(str(n_neutrals))
	desc = desc//NL//mpih%myrankstr//'    	Effective plasma density: '//trim(str(real(np)/domain_volume))
	print '(A)', desc
	endassociate
	endsubroutine uniform_domain_space_injection

	!subroutine uniform_box_space_injection(self, field, pic, q_pic)

	!endsubroutine uniform_box_space_injection

	subroutine uniform_cell_space_injection(self, field, grid, pic, q_pic)
	class(prism_particle_injection_object), intent(inout) :: self
	type(field_object),                  	 intent(in) 	:: field
	type(grid_object),                       intent(in) 	:: grid !< Grid (sibling realm component, threaded in).
	type(prism_pic_object),					 	 intent(in)		:: pic
	real(R8P),                           	 intent(inout) :: q_pic(1:,1:)
	real(R8P)															:: r_n(3)
	integer(I4P)														:: n_i_4c
	integer(I4P)														:: n_e_4c
	integer(I4P)														:: n_n_4c
	integer(I4P)														:: n_cells
	real(R8P)															:: x_min, y_min, z_min
	real(R8P)															:: x_max, y_max, z_max
	real(R8P)															:: x_p, y_p, z_p
	real(R8P)															:: deltax, deltay, deltaz
	real(R8P)															:: domain_volume
	integer(I4P)														:: i, j, k, b, p, n_p
	integer(I4P), allocatable										:: shuffled_list_ions(:,:)
	integer(I4P), allocatable										:: shuffled_list_electrons(:,:)
	integer(I4P), allocatable										:: shuffled_list_neutrals(:,:)
	character(len=:), allocatable		                   	:: desc
   character(len=1), parameter  		                   	:: NL=new_line('a')

	associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj, nk=>grid%nk,              &
      		ngc=>grid%ngc, dx=>field%dxyz(1,:), dy=>field%dxyz(2,:), dz=>field%dxyz(3,:),          &
      		np=>pic%particle_number, e_min=>grid%domain_emin, e_max=>grid%domain_emax,  		 &
				neutral_fraction=>pic%neutral_fraction, n_ions=>pic%n_ions, n_electrons=>pic%n_electrons, 	 &
				n_neutrals=>pic%n_neutrals)

	n_cells = ni*nj*nk
	n_i_4c = n_ions/n_cells
	n_e_4c = n_i_4c
	n_n_4c = n_neutrals/n_cells
	allocate(shuffled_list_ions	  (1:3,1:n_i_4c))
	allocate(shuffled_list_electrons(1:3,1:n_e_4c))
	allocate(shuffled_list_neutrals (1:3,1:n_n_4c))
	shuffled_list_ions(:,:) 	  = 0_I4P
	shuffled_list_electrons(:,:) = 0_I4P
	shuffled_list_neutrals(:,:)  = 0_I4P

	!if (.not.space_pairing) then
		p = 0_I4P
		do b=1,blocks_number
			deltax = dx(b)
			deltay = dy(b)
			deltaz = dz(b)
			do k=1,nk
				do j=1,nj
					do i=1,ni
						! Pensato per il monoblocco, eventualmente va esteso sul multiblocco
						x_min = e_min(1) + deltax * real(i-1,R8P)
						x_max = e_min(1) + deltax * real(i	,R8P)
						y_min = e_min(2) + deltay * real(j-1,R8P)
						y_max = e_min(2) + deltay * real(j	,R8P)
						z_min = e_min(3) + deltaz * real(k-1,R8P)
						z_max = e_min(3) + deltaz * real(k	,R8P)
						do n_p = 1, n_i_4c
							p = p + 1_I4P
							call space_rand_num_generator(N=n_i_4c, &
									shuffled_list=shuffled_list_ions, i_numb=n_p, r_n=r_n)
							x_p = x_min+r_n(1)*(x_max-x_min)
							y_p = y_min+r_n(2)*(y_max-y_min)
							z_p = z_min+r_n(3)*(z_max-z_min)
							q_pic(1,p) = x_p
							q_pic(2,p) = y_p
							q_pic(3,p) = z_p
							q_pic(7,p) = -E_CHARGE
							q_pic(8,p) = E_MASS
							!Z0 = sqrt(-2.0_R8P*log(r_n(1)))*cos(2*PI*r_n(2))
							!Z1 = sqrt(-2.0_R8P*log(r_n(1)))*cos(2*PI*r_n(2))
						enddo
					enddo
				enddo
			enddo
		enddo
		do b=1,blocks_number
			deltax = dx(b)
			deltay = dy(b)
			deltaz = dz(b)
			do k=1,nk
				do j=1,nj
					do i=1,ni
						! Pensato per il monoblocco, eventualmente va esteso sul multiblocco
						x_min = e_min(1) + deltax * real(i-1,R8P)
						x_max = e_min(1) + deltax * real(i	,R8P)
						y_min = e_min(2) + deltay * real(j-1,R8P)
						y_max = e_min(2) + deltay * real(j	,R8P)
						z_min = e_min(3) + deltaz * real(k-1,R8P)
						z_max = e_min(3) + deltaz * real(k	,R8P)
						do n_p = 1, n_e_4c
							p = p + 1_I4P
							call space_rand_num_generator(N=n_e_4c, &
									shuffled_list=shuffled_list_electrons, i_numb=n_p, r_n=r_n)
							x_p = x_min+r_n(1)*(x_max-x_min)
							y_p = y_min+r_n(2)*(y_max-y_min)
							z_p = z_min+r_n(3)*(z_max-z_min)
							q_pic(1,p) = x_p
							q_pic(2,p) = y_p
							q_pic(3,p) = z_p
							q_pic(7,p) = E_CHARGE
							q_pic(8,p) = E_MASS
						enddo
					enddo
				enddo
			enddo
		enddo
		do b=1,blocks_number
			deltax = dx(b)
			deltay = dy(b)
			deltaz = dz(b)
			do k=1,nk
				do j=1,nj
					do i=1,ni
						! Pensato per il monoblocco, eventualmente va esteso sul multiblocco
						x_min = e_min(1) + deltax * real(i-1,R8P)
						x_max = e_min(1) + deltax * real(i	,R8P)
						y_min = e_min(2) + deltay * real(j-1,R8P)
						y_max = e_min(2) + deltay * real(j	,R8P)
						z_min = e_min(3) + deltaz * real(k-1,R8P)
						z_max = e_min(3) + deltaz * real(k	,R8P)
						do n_p = 1, n_n_4c
							p = p + 1_I4P
							call space_rand_num_generator(N=n_n_4c, &
									shuffled_list=shuffled_list_neutrals, i_numb=n_p, r_n=r_n)
							x_p = x_min+r_n(1)*(x_max-x_min)
							y_p = y_min+r_n(2)*(y_max-y_min)
							z_p = z_min+r_n(3)*(z_max-z_min)
							q_pic(1,p) = x_p
							q_pic(2,p) = y_p
							q_pic(3,p) = z_p
							q_pic(7,p) = 0.0_R8P
							q_pic(8,p) = E_MASS
						enddo
					enddo
				enddo
			enddo
		enddo
	!else

	!endif

	domain_volume = (e_max(1)-e_min(1))*(e_max(2)-e_min(2))*(e_max(3)-e_min(3))

	desc =       mpih%myrankstr//'Injected particles:'
   desc = desc//NL//mpih%myrankstr//'    	Electrons number: '//trim(str(n_e_4c*n_cells))
	desc = desc//NL//mpih%myrankstr//'    	Ions number: '//trim(str(n_i_4c*n_cells))
	desc = desc//NL//mpih%myrankstr//'    	Neutrals number: '//trim(str(n_n_4c*n_cells))
	desc = desc//NL//mpih%myrankstr//'		Effective particle total number: ' &
													//trim(str((n_e_4c*n_cells+n_i_4c*n_cells+n_n_4c*n_cells)))
	desc = desc//NL//mpih%myrankstr//'		Effective plasma density: ' &
													//trim(str(real(n_e_4c*n_cells+n_i_4c*n_cells+n_n_4c*n_cells)/domain_volume))
	print '(A)', desc
	endassociate
	endsubroutine uniform_cell_space_injection

	subroutine uniform_maxwellian_velocity_injection(self, field, grid, pic, q_pic)
	class(prism_particle_injection_object), intent(inout) :: self
	type(field_object),                  	 intent(in) 	:: field
	type(grid_object),                       intent(in) 	:: grid !< Grid (sibling realm component, threaded in).
	type(prism_pic_object),					 	 intent(in)		:: pic
	real(R8P),                           	 intent(inout) :: q_pic(1:,1:)
	real(R8P)															:: r_n(4)
	real(R8P)															:: vx_p, vy_p, vz_p
	real(R8P)															:: Zx, Zy, Zz
	real(R8P)															:: v_t
	integer(I4P)														:: i
	integer(I4P), allocatable										:: shuffled_list_ions(:,:)
	integer(I4P), allocatable										:: shuffled_list_electrons(:,:)
	integer(I4P), allocatable										:: shuffled_list_neutrals(:,:)

	associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj, nk=>grid%nk,              &
      		ngc=>grid%ngc, dx=>field%dxyz(1,:), dy=>field%dxyz(2,:), dz=>field%dxyz(3,:),          &
      		np=>pic%particle_number, e_min=>grid%domain_emin, e_max=>grid%domain_emax,  		 &
				neutral_fraction=>pic%neutral_fraction, v_av_correction=>self%v_av_correction, 				 	 &
				v_drift_x=>self%v_drift_x, v_drift_y=>self%v_drift_y, v_drift_z=>self%v_drift_z, 				 &
				T_i=>self%T_i, T_e=>self%T_e, T_n=>self%T_n, 																 &
				n_ions=>pic%n_ions, n_electrons=>pic%n_electrons, n_neutrals=>pic%n_neutrals)

	allocate(shuffled_list_ions	  (1:4,1:n_ions))
	allocate(shuffled_list_electrons(1:4,1:n_electrons))
	allocate(shuffled_list_neutrals (1:4,1:n_neutrals))
	shuffled_list_ions(:,:) 	  = 0_I4P
	shuffled_list_electrons(:,:) = 0_I4P
	shuffled_list_neutrals(:,:)  = 0_I4P

	!if(.not.space_pairing) then
		do i = 1, n_ions
			call space_rand_num_generator(N=n_ions, shuffled_list=shuffled_list_ions, i_numb=i, r_n=r_n)
			v_t = sqrt(K_B*T_i/q_pic(8,i))
			Zx = sqrt(-2.0_R8P*log(r_n(1)))*cos(2*PI*r_n(2))
			Zy = sqrt(-2.0_R8P*log(r_n(1)))*sin(2*PI*r_n(2))
			Zz = sqrt(-2.0_R8P*log(r_n(3)))*cos(2*PI*r_n(4))
			vx_p = v_t*Zx
			vy_p = v_t*Zy
			vz_p = v_t*Zz
			q_pic(4,i) = vx_p
			q_pic(5,i) = vy_p
			q_pic(6,i) = vz_p
		enddo
		call add_drift_velocity(q_pic=q_pic(4,1:n_ions), v_drift=v_drift_x)
		call add_drift_velocity(q_pic=q_pic(5,1:n_ions), v_drift=v_drift_y)
		call add_drift_velocity(q_pic=q_pic(6,1:n_ions), v_drift=v_drift_z)
		if (v_av_correction) then
			call apply_vel_av_correction(q_pic=q_pic(4,1:n_ions), v_drift = v_drift_x)
			call apply_vel_av_correction(q_pic=q_pic(5,1:n_ions), v_drift = v_drift_y)
			call apply_vel_av_correction(q_pic=q_pic(6,1:n_ions), v_drift = v_drift_z)
		endif
		do i = 1, n_electrons
			call space_rand_num_generator(N=n_electrons, shuffled_list=shuffled_list_electrons, i_numb=i, r_n=r_n)
			v_t = sqrt(K_B*T_e/q_pic(8,i))
			Zx = sqrt(-2.0_R8P*log(r_n(1)))*cos(2*PI*r_n(2))
			Zy = sqrt(-2.0_R8P*log(r_n(1)))*sin(2*PI*r_n(2))
			Zz = sqrt(-2.0_R8P*log(r_n(3)))*cos(2*PI*r_n(4))
			vx_p = v_t*Zx
			vy_p = v_t*Zy
			vz_p = v_t*Zz
			q_pic(4,i+n_ions) = vx_p
			q_pic(5,i+n_ions) = vy_p
			q_pic(6,i+n_ions) = vz_p
		enddo
		call add_drift_velocity(q_pic=q_pic(4,n_ions+1:n_ions+n_electrons), v_drift=v_drift_x)
		call add_drift_velocity(q_pic=q_pic(5,n_ions+1:n_ions+n_electrons), v_drift=v_drift_y)
		call add_drift_velocity(q_pic=q_pic(6,n_ions+1:n_ions+n_electrons), v_drift=v_drift_z)
		if (v_av_correction) then
			call apply_vel_av_correction(q_pic=q_pic(4,n_ions+1:n_ions+n_electrons), v_drift = v_drift_x)
			call apply_vel_av_correction(q_pic=q_pic(5,n_ions+1:n_ions+n_electrons), v_drift = v_drift_y)
			call apply_vel_av_correction(q_pic=q_pic(6,n_ions+1:n_ions+n_electrons), v_drift = v_drift_z)
		endif
		do i = 1, n_neutrals
			call space_rand_num_generator(N=n_neutrals, shuffled_list=shuffled_list_neutrals, i_numb=i, r_n=r_n)
			v_t = sqrt(K_B*T_n/q_pic(8,i))
			Zx = sqrt(-2.0_R8P*log(r_n(1)))*cos(2*PI*r_n(2))
			Zy = sqrt(-2.0_R8P*log(r_n(1)))*sin(2*PI*r_n(2))
			Zz = sqrt(-2.0_R8P*log(r_n(3)))*cos(2*PI*r_n(4))
			vx_p = v_t*Zx
			vy_p = v_t*Zy
			vz_p = v_t*Zz
			q_pic(4,i+n_ions+n_electrons) = vx_p
			q_pic(5,i+n_ions+n_electrons) = vy_p
			q_pic(6,i+n_ions+n_electrons) = vz_p
		enddo
		call add_drift_velocity(q_pic=q_pic(4,n_ions+n_electrons+1:np), v_drift=v_drift_x)
		call add_drift_velocity(q_pic=q_pic(5,n_ions+n_electrons+1:np), v_drift=v_drift_y)
		call add_drift_velocity(q_pic=q_pic(6,n_ions+n_electrons+1:np), v_drift=v_drift_z)
		if (v_av_correction) then
			call apply_vel_av_correction(q_pic=q_pic(4,n_ions+n_electrons+1:np), v_drift = v_drift_x)
			call apply_vel_av_correction(q_pic=q_pic(5,n_ions+n_electrons+1:np), v_drift = v_drift_y)
			call apply_vel_av_correction(q_pic=q_pic(6,n_ions+n_electrons+1:np), v_drift = v_drift_z)
		endif
	!else

	!endif

	endassociate
	endsubroutine uniform_maxwellian_velocity_injection

	subroutine non_uniform_maxwellian_velocity_injection(self, field, grid, pic, q_pic)
	class(prism_particle_injection_object), intent(inout) :: self
	type(field_object),                  	 intent(in) 	:: field
	type(grid_object),                       intent(in) 	:: grid !< Grid (sibling realm component, threaded in).
	type(prism_pic_object),					 	 intent(in)		:: pic
	real(R8P),                           	 intent(inout) :: q_pic(1:,1:)
	real(R8P)															:: r_n(4)
	real(R8P)															:: vx_p, vy_p, vz_p
	real(R8P)															:: Zx, Zy, Zz
	real(R8P)															:: v_tx, v_ty, v_tz
	integer(I4P)														:: i
	integer(I4P), allocatable										:: shuffled_list_ions(:,:)
	integer(I4P), allocatable										:: shuffled_list_electrons(:,:)
	integer(I4P), allocatable										:: shuffled_list_neutrals(:,:)

	associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj, nk=>grid%nk,              &
      		ngc=>grid%ngc, dx=>field%dxyz(1,:), dy=>field%dxyz(2,:), dz=>field%dxyz(3,:),          &
      		np=>pic%particle_number, e_min=>grid%domain_emin, e_max=>grid%domain_emax,  		 &
				neutral_fraction=>pic%neutral_fraction, 																		 &
				T_i_x=>self%T_i_x, T_e_x=>self%T_e_x, T_n_x=>self%T_n_x, 												 &
				T_i_y=>self%T_i_y, T_e_y=>self%T_e_y, T_n_y=>self%T_n_y, 												 &
				T_i_z=>self%T_i_z, T_e_z=>self%T_e_z, T_n_z=>self%T_n_z, 												 &
				v_av_correction=>self%v_av_correction, v_drift_x=>self%v_drift_x, v_drift_y=>self%v_drift_y,  &
				v_drift_z=>self%v_drift_z, n_ions=>pic%n_ions, n_electrons=>pic%n_electrons, 						 &
				n_neutrals=>pic%n_neutrals)

	allocate(shuffled_list_ions	  (1:4,1:n_ions))
	allocate(shuffled_list_electrons(1:4,1:n_electrons))
	allocate(shuffled_list_neutrals (1:4,1:n_neutrals))
	shuffled_list_ions(:,:) 	  = 0_I4P
	shuffled_list_electrons(:,:) = 0_I4P
	shuffled_list_neutrals(:,:)  = 0_I4P

	!if(.not.velocity_pairing) then
		do i = 1, n_ions
			call velocity_rand_num_generator(N=n_ions, shuffled_list=shuffled_list_ions, i_numb=i, r_n=r_n)
			v_tx = sqrt(K_B*T_i_x/q_pic(8,i))
			v_ty = sqrt(K_B*T_i_y/q_pic(8,i))
			v_tz = sqrt(K_B*T_i_z/q_pic(8,i))
			Zx = sqrt(-2.0_R8P*log(r_n(1)))*cos(2*PI*r_n(2))
			Zy = sqrt(-2.0_R8P*log(r_n(1)))*sin(2*PI*r_n(2))
			Zz = sqrt(-2.0_R8P*log(r_n(3)))*cos(2*PI*r_n(4))
			vx_p = v_tx*Zx
			vy_p = v_ty*Zy
			vz_p = v_tz*Zz
			q_pic(4,i) = vx_p
			q_pic(5,i) = vy_p
			q_pic(6,i) = vz_p
		enddo
		call add_drift_velocity(q_pic=q_pic(4,1:n_ions), v_drift=v_drift_x)
		call add_drift_velocity(q_pic=q_pic(5,1:n_ions), v_drift=v_drift_y)
		call add_drift_velocity(q_pic=q_pic(6,1:n_ions), v_drift=v_drift_z)
		if (v_av_correction) then
			call apply_vel_av_correction(q_pic=q_pic(4,1:n_ions), v_drift = v_drift_x)
			call apply_vel_av_correction(q_pic=q_pic(5,1:n_ions), v_drift = v_drift_y)
			call apply_vel_av_correction(q_pic=q_pic(6,1:n_ions), v_drift = v_drift_z)
		endif
		do i = 1, n_electrons
			call velocity_rand_num_generator(N=n_electrons, shuffled_list=shuffled_list_electrons, i_numb=i, r_n=r_n)
			v_tx = sqrt(K_B*T_e_x/q_pic(8,i))
			v_ty = sqrt(K_B*T_e_y/q_pic(8,i))
			v_tz = sqrt(K_B*T_e_z/q_pic(8,i))
			Zx = sqrt(-2.0_R8P*log(r_n(1)))*cos(2*PI*r_n(2))
			Zy = sqrt(-2.0_R8P*log(r_n(1)))*sin(2*PI*r_n(2))
			Zz = sqrt(-2.0_R8P*log(r_n(3)))*cos(2*PI*r_n(4))
			vx_p = v_tx*Zx
			vy_p = v_ty*Zy
			vz_p = v_tz*Zz
			q_pic(4,i+n_ions) = vx_p
			q_pic(5,i+n_ions) = vy_p
			q_pic(6,i+n_ions) = vz_p
		enddo
		call add_drift_velocity(q_pic=q_pic(4,n_ions+1:n_ions+n_electrons), v_drift=v_drift_x)
		call add_drift_velocity(q_pic=q_pic(5,n_ions+1:n_ions+n_electrons), v_drift=v_drift_y)
		call add_drift_velocity(q_pic=q_pic(6,n_ions+1:n_ions+n_electrons), v_drift=v_drift_z)
		if (v_av_correction) then
			call apply_vel_av_correction(q_pic=q_pic(4,n_ions+1:n_ions+n_electrons), v_drift = v_drift_x)
			call apply_vel_av_correction(q_pic=q_pic(5,n_ions+1:n_ions+n_electrons), v_drift = v_drift_y)
			call apply_vel_av_correction(q_pic=q_pic(6,n_ions+1:n_ions+n_electrons), v_drift = v_drift_z)
		endif
		do i = 1, n_neutrals
			call velocity_rand_num_generator(N=n_neutrals, shuffled_list=shuffled_list_neutrals, i_numb=i, r_n=r_n)
			v_tx = sqrt(K_B*T_n_x/q_pic(8,i))
			v_ty = sqrt(K_B*T_n_y/q_pic(8,i))
			v_tz = sqrt(K_B*T_n_z/q_pic(8,i))
			Zx = sqrt(-2.0_R8P*log(r_n(1)))*cos(2*PI*r_n(2))
			Zy = sqrt(-2.0_R8P*log(r_n(1)))*sin(2*PI*r_n(2))
			Zz = sqrt(-2.0_R8P*log(r_n(3)))*cos(2*PI*r_n(4))
			vx_p = v_tx*Zx
			vy_p = v_ty*Zy
			vz_p = v_tz*Zz
			q_pic(4,i+n_ions+n_electrons) = vx_p
			q_pic(5,i+n_ions+n_electrons) = vy_p
			q_pic(6,i+n_ions+n_electrons) = vz_p
		enddo
		call add_drift_velocity(q_pic=q_pic(4,n_ions+n_electrons+1:np), v_drift=v_drift_x)
		call add_drift_velocity(q_pic=q_pic(5,n_ions+n_electrons+1:np), v_drift=v_drift_y)
		call add_drift_velocity(q_pic=q_pic(6,n_ions+n_electrons+1:np), v_drift=v_drift_z)
		if (v_av_correction) then
			call apply_vel_av_correction(q_pic=q_pic(4,n_ions+n_electrons+1:np), v_drift = v_drift_x)
			call apply_vel_av_correction(q_pic=q_pic(5,n_ions+n_electrons+1:np), v_drift = v_drift_y)
			call apply_vel_av_correction(q_pic=q_pic(6,n_ions+n_electrons+1:np), v_drift = v_drift_z)
		endif
	!else

	!endif

	endassociate
	endsubroutine non_uniform_maxwellian_velocity_injection

	subroutine write_initial_injection_tab(filename, q_pic, np)
	character(len=1), parameter  :: TAB = achar(9)
	character(len=*), intent(in) :: filename
	integer(I4P),     intent(in) :: np
	real(R8P),        intent(in) :: q_pic(1:, 1:)
	integer(I4P) 					  :: iu, i, j, ios, l

	l = size(q_pic, dim=1)
	open(newunit=iu, file=trim(filename), status='replace', action='write', &
	     form='formatted', iostat=ios)
	if (ios /= 0) then
	  write(*,'(a,i0)') 'write_dat_tab: errore open(), iostat=', ios
	  error stop
	end if
	do i = 1, np
	  write(iu,'(ES24.16,7(a,ES24.16))') q_pic(1,i), (TAB, q_pic(j,i), j=2,l)
	end do
	close(iu)
	endsubroutine write_initial_injection_tab

	subroutine random_number_generator(N, shuffled_list, i_numb, r_n)
	integer(I4P), intent(inout) :: shuffled_list(1:,1:)
	integer(I4P), intent(in) 	 :: i_numb !Non utilizzato qui
	integer(I4P), intent(in) 	 :: N !Numero di elementi Non utilizzato qui
	real(R8P), intent(inout) 	 :: r_n(1:) !Random numbers Non utilizzato qui

	call random_number(r_n)
	endsubroutine random_number_generator

	subroutine layered_number_generator(N, shuffled_list, i_numb, r_n)
	integer(I4P), intent(inout) :: shuffled_list(1:,1:)
	integer(I4P), intent(in) 	 :: i_numb
	integer(I4P), intent(in) 	 :: N  		 					!Numero di elementi
	real(R8P), intent(inout) 	 :: r_n(1:) 					!Random numbers
	integer(I4P) 					 :: index_list(N)
	integer(I4P) 					 :: w, h
	integer(I4P)					 :: n_rn

	n_rn = size(r_n)
	if (i_numb == 1) then
		do w = 1, N
			index_list(w) = w
		enddo
		do h = 1, n_rn
			if (h == 1) then
				shuffled_list(1,:) = index_list
			else
				shuffled_list(h,:) = fisher_yates_shuffle(index_list=index_list, nn=N)
			endif
			r_n(h) = (real(shuffled_list(h,i_numb),R8P)-0.5_R8P)/real(N,R8P)
		enddo
	else
		do h = 1, n_rn
			r_n(h) = (real(shuffled_list(h,i_numb),R8P)-0.5_R8P)/real(N,R8P)
		enddo
	endif
	endsubroutine layered_number_generator

	function fisher_yates_shuffle(index_list, nn) result(shuffled_list)
	integer(I4P), intent(in)  :: nn
	integer(I4P), intent(in)  :: index_list(nn)
	integer(I4P)				  :: shuffled_list(nn)
	integer(I4P) 				  :: ii, jj, tmp
	real(R8P) 	 				  :: u

	shuffled_list = index_list
	do ii = nn, 2, -1
      call random_number(u)     ! u in [0,1)
      jj = 1 + int(u*real(ii))    ! j in [1,i]
      tmp    			   = shuffled_list(ii)
      shuffled_list(ii) = shuffled_list(jj)
      shuffled_list(jj) = tmp
   enddo
	endfunction fisher_yates_shuffle

	subroutine add_drift_velocity(q_pic, v_drift)
	real(R8P), intent(inout) :: q_pic(1:)
	real(R8P), intent(in)	 :: v_drift

	q_pic(:) = q_pic(:) + v_drift
	endsubroutine add_drift_velocity

	subroutine apply_vel_av_correction(q_pic, v_drift)
	real(R8P), intent(inout) :: q_pic(1:)
	real(R8P), intent(in)	 :: v_drift
	real(R8P)					 :: v_av, correction

	v_av = sum(q_pic) /real(size(q_pic), kind=R8P)
	correction = (v_drift - v_av)
	q_pic(:) = q_pic(:) + correction
	endsubroutine apply_vel_av_correction

   endmodule adam_prism_particle_injection_object
