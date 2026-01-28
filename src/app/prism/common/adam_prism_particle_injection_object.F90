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

procedure(space_random_number_generator_interface), pointer :: space_rand_num_generator => null() !< Space random number generator interface

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
   procedure(particle_space_injection_interface), pass(self), pointer :: particle_space_injection => null() !< Particle space injection.
	!procedure(particle_velocity_injection_interface), 	 	 pass(self), pointer :: particle_velocity_injection => null() !< Particle velocity injection.
	!procedure(velocity_random_number_generator_interface), pass(self)?, pointer :: velocity_rand_num_generator => null() !< Space random number generator interface
contains
   procedure, pass(self) :: description                   !< Return pretty-printed object description.
   procedure, pass(self) :: initialize                    !< Initialize IC.
   procedure, pass(self) :: load_from_file                !< Load config from file.
endtype prism_particle_injection_object

interface
   subroutine particle_space_injection_interface(self, field, pic, q_pic)
   import :: prism_particle_injection_object, field_object, prism_pic_object, I4P, R8P
	class(prism_particle_injection_object), intent(inout) :: self 
	type(field_object),                  	 intent(in) 	:: field 
	type(prism_pic_object),					 	 intent(in)		:: pic
	real(R8P),                           	 intent(inout) :: q_pic(1:,1:)                                                         !< Number of variables.
   endsubroutine particle_space_injection_interface

	subroutine space_random_number_generator_interface(N, shuffled_list, i_numb, r_n)
	import :: I4P, R8P
	integer(I4P), intent(inout) :: shuffled_list(1:,1:)
	integer(I4P), intent(in) 	 :: i_numb
	integer(I4P), intent(in) 	 :: N 							!Numero di elementi
	real(R8P), intent(inout) 	 :: r_n(1:) 					!Random numbers
	endsubroutine space_random_number_generator_interface
endinterface

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

	select case(self%space_distribution)
   case(UNIFORM_CELL_SPACE_DISTRIBUTION)
      self%particle_space_injection => uniform_cell_space_injection
   !case(UNIFORM_BOX_SPACE_DISTRIBUTION)
   !   self%particle_space_injection => uniform_box_space_injection
   case(UNIFORM_DOMAIN_SPACE_DISTRIBUTION)
      self%particle_space_injection => uniform_domain_space_injection
   case default
      call self%mpih%error_stop &
		(msg=': invalid particle space injection model in prism_particle_injection_object%initialize')
   endselect
	select case(self%space_random_number_generator)
   case(SPACE_RANDOM_NUMBER_GENERATOR)
      space_rand_num_generator => random_number_generator
   case(SPACE_LAYERED_NUMBER_GENERATOR)
      space_rand_num_generator => layered_number_generator
   case default
      call self%mpih%error_stop & 
		(msg=': invalid particle space random number generator in prism_particle_injection_object%initialize')
   endselect
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

   subroutine uniform_domain_space_injection(self, field, pic, q_pic)
	class(prism_particle_injection_object), intent(inout) :: self 
	type(field_object),                  	 intent(in) 	:: field 
	type(prism_pic_object),					 	 intent(in)		:: pic
	real(R8P),                           	 intent(inout) :: q_pic(1:,1:)
	real(R8P)															:: r_n(3)
	integer(I4P)														:: n_ions
	integer(I4P)														:: n_electrons
	integer(I4P)														:: n_neutrals
	real(R8P)															:: x_p, y_p, z_p
	integer(I4P)														:: i
	integer(I4P), allocatable										:: shuffled_list_ions(:,:)
	integer(I4P), allocatable										:: shuffled_list_electrons(:,:)
	integer(I4P), allocatable										:: shuffled_list_neutrals(:,:)
	character(len=:), allocatable		                   	:: desc             
   character(len=1), parameter  		                   	:: NL=new_line('a') 

	associate(blocks_number=>field%blocks_number, ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, &
      		ngc=>field%grid%ngc, dx=>field%dxyz(1,:), dy=>field%dxyz(2,:), dz=>field%dxyz(3,:), 	 		 &
      		np=>pic%particle_number, e_min=>field%grid%domain_emin, e_max=>field%grid%domain_emax,  		 &
				neutral_fraction=>pic%neutral_fraction)

	n_neutrals = nint(neutral_fraction*real(np,R8P))
	n_ions = nint(real(np-n_neutrals, R8P)/2.0_R8P)
	n_electrons = n_ions
	n_neutrals = np - n_ions - n_electrons
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
	
	desc =       self%mpih%myrankstr//'Injected particles:'
   desc = desc//NL//self%mpih%myrankstr//'    	Electrons number: '//trim(str(n_electrons))
	desc = desc//NL//self%mpih%myrankstr//'    	Ions number: '//trim(str(n_ions))
	desc = desc//NL//self%mpih%myrankstr//'    	Neutrals number: '//trim(str(n_neutrals))
	endassociate
	endsubroutine uniform_domain_space_injection

	!subroutine uniform_box_space_injection(self, field, pic, q_pic)

	!endsubroutine uniform_box_space_injection

	subroutine uniform_cell_space_injection(self, field, pic, q_pic)
	class(prism_particle_injection_object), intent(inout) :: self 
	type(field_object),                  	 intent(in) 	:: field 
	type(prism_pic_object),					 	 intent(in)		:: pic
	real(R8P),                           	 intent(inout) :: q_pic(1:,1:)
	real(R8P)															:: r_n(3)
	integer(I4P)														:: n_ions, n_i_4c
	integer(I4P)														:: n_electrons, n_e_4c
	integer(I4P)														:: n_neutrals, n_n_4c
	integer(I4P)														:: n_cells
	real(R8P)															:: x_min, y_min, z_min
	real(R8P)															:: x_max, y_max, z_max
	real(R8P)															:: x_p, y_p, z_p
	real(R8P)															:: deltax, deltay, deltaz
	integer(I4P)														:: i, j, k, b, p, n_p
	integer(I4P), allocatable										:: shuffled_list_ions(:,:)
	integer(I4P), allocatable										:: shuffled_list_electrons(:,:)
	integer(I4P), allocatable										:: shuffled_list_neutrals(:,:)
	character(len=:), allocatable		                   	:: desc             
   character(len=1), parameter  		                   	:: NL=new_line('a') 

	associate(blocks_number=>field%blocks_number, ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, &
      		ngc=>field%grid%ngc, dx=>field%dxyz(1,:), dy=>field%dxyz(2,:), dz=>field%dxyz(3,:), 	 		 &
      		np=>pic%particle_number, e_min=>field%grid%domain_emin, e_max=>field%grid%domain_emax,  		 &
				neutral_fraction=>pic%neutral_fraction)

	n_neutrals = nint(neutral_fraction*real(np,R8P))
	n_ions = nint(real(np-n_neutrals, R8P)/2.0_R8P)
	n_electrons = n_ions
	n_neutrals = np - n_ions - n_electrons
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
	desc =       self%mpih%myrankstr//'Injected particles:'
   desc = desc//NL//self%mpih%myrankstr//'    	Electrons number: '//trim(str(n_e_4c*n_cells))
	desc = desc//NL//self%mpih%myrankstr//'    	Ions number: '//trim(str(n_i_4c*n_cells))
	desc = desc//NL//self%mpih%myrankstr//'    	Neutrals number: '//trim(str(n_n_4c*n_cells))
	endassociate
	endsubroutine uniform_cell_space_injection

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
	integer(I4P), intent(in) 	 :: N 							!Numero di elementi
	real(R8P), intent(inout) 	 :: r_n(1:) 					!Random numbers
	integer(I4P) 					 :: index_list(N)
	integer(I4P) 					 :: w

	if (i_numb == 1) then
		do w = 1, N 
			index_list(w) = w
		enddo
		shuffled_list(1,:) = index_list
		shuffled_list(2,:) = fisher_yates_shuffle(index_list=index_list, nn=N)
		shuffled_list(3,:) = fisher_yates_shuffle(index_list=index_list, nn=N)
		r_n(1) = (real(shuffled_list(1,i_numb),R8P)-0.5_R8P)/real(N,R8P)
		r_n(2) = (real(shuffled_list(2,i_numb),R8P)-0.5_R8P)/real(N,R8P)
		r_n(3) = (real(shuffled_list(3,i_numb),R8P)-0.5_R8P)/real(N,R8P)
	else
		r_n(1) = (real(shuffled_list(1,i_numb),R8P)-0.5_R8P)/real(N,R8P)
		r_n(2) = (real(shuffled_list(2,i_numb),R8P)-0.5_R8P)/real(N,R8P)
		r_n(3) = (real(shuffled_list(3,i_numb),R8P)-0.5_R8P)/real(N,R8P)
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
   endmodule adam_prism_particle_injection_object