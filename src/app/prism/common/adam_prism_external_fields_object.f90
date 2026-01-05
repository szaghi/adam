!< ADAM, PRISM external fields definition, CPU backend.
module adam_prism_external_fields_object
!< ADAM, PRISM external fields definition, CPU backend.
! ADAM modules
use adam_mpih_object,  only : mpih_object
use adam_field_object, only : field_object
! PRISM modules
use adam_prism_parameters
! third party modules
use finer
use penf

implicit none 
private
public :: INI_SECTION_NAME
public :: RMF
public :: MAGNETIC_NOZZLE
public :: RMF_AND_MAGNETIC_NOZZLE
public :: prism_external_fields_object
public :: add_external_fields_interface
public :: add_external_fields_rmf
!public :: add_external_fields_magnetic_nozzle
!public :: add_external_fields_rmf_and_magnetic_nozzle
public :: add_external_fields_none

character(len=15), parameter :: INI_SECTION_NAME        = 'external_fields'         !< INI (config) file section name containing external fields configs.
character(len=3),  parameter :: RMF                     = 'RMF'                     !< Rotating Magnetic Field.
character(len=15), parameter :: MAGNETIC_NOZZLE         = 'Magnetic_nozzle'         !< Magnetic Nozzle.
character(len=23), parameter :: RMF_AND_MAGNETIC_NOZZLE = 'RMF_and_magnetic_nozzle' !< Rotating Magnetic Field and Magnetic Nozzle.

type :: prism_external_fields_object
    !< PRISM external fields object.
   type(mpih_object) :: mpih                     !< MPI handler.
   character(len=99) :: external_field_applied   !< Field type.
   real(R8P)         :: RMF_frequency            !< Rotating magnetic field frequency.
   real(R8P)         :: RMF_B_amplitude          !< Rotating magnetic field amplitude.
	character(len=99) :: RMF_rotation_axis 		  !< Rotating magnetic field rotation axis (X, Y, Z).
	integer(I4P)      :: alpha                    !< RMF rotation axis coordinate 1
	integer(I4P)      :: beta                     !< RMF rotation axis coordinate 2
	integer(I4P)      :: gamma                    !< RMF rotation axis coordinate 3

contains
   procedure, pass(self) :: description                           !< Return pretty-printed object description.
   procedure, pass(self) :: initialize                            !< Initialize IC.
   procedure, pass(self) :: load_from_file                        !< Load config from file.
   procedure, pass(self) :: add_external_fields_rmf               !< Add rotating magnetic field to the field.
   !procedure, pass(self) :: add_external_fields_magnetic_nozzle  !< Add magnetic nozzle to the field.
   !procedure, pass(self) :: add_external_fields_rmf_and_magnetic_nozzle !< Add rotating magnetic field and magnetic nozzle to the field.
   procedure, pass(self) :: add_external_fields_none              !< No external field applied.

endtype prism_external_fields_object

interface
   subroutine add_external_fields_interface(self, field, time, dq)
   import :: prism_external_fields_object, field_object, I4P, R8P
   class(prism_external_fields_object), intent(inout) :: self                                                              !< External fields.
   type(field_object),                  intent(inout) :: field                                                             !< The field.
   real(R8P),                           intent(in)    :: time                                                              !< Current simulation time.
   real(R8P),                           intent(inout) :: dq(1:, 1-field%grid%ngc:,1-field%grid%ngc:,1-field%grid%ngc:,1:)  !< Primitive variables.
   endsubroutine add_external_fields_interface
endinterface

contains 
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(prism_external_fields_object), intent(in) :: self             !< External fields.
   character(len=:), allocatable                   :: desc             !< Description.
   character(len=1), parameter                     :: NL=new_line('a') !< New line character.
   desc =       self%mpih%myrankstr//'Applied external fields:'
   select case(self%external_field_applied)
   case(RMF)
   desc = desc//NL//self%mpih%myrankstr//'    Rotating magnetic field applied '
   desc = desc//NL//self%mpih%myrankstr//'    RMF frequency: '//trim(str(self%RMF_frequency))
   desc = desc//NL//self%mpih%myrankstr//'    RMF B amplitude: '//trim(str(self%RMF_B_amplitude))
   desc = desc//NL//self%mpih%myrankstr//'    RMF rotation axis: '//trim(self%RMF_rotation_axis)
   case(MAGNETIC_NOZZLE)
   desc = desc//NL//self%mpih%myrankstr//'    Magnetic nozzle applied '
   case(RMF_AND_MAGNETIC_NOZZLE)
   desc = desc//NL//self%mpih%myrankstr//'    Rotating magnetic field and magnetic nozzle applied '
   desc = desc//NL//self%mpih%myrankstr//'    RMF frequency: '//trim(str(self%RMF_frequency))
   desc = desc//NL//self%mpih%myrankstr//'    RMF B amplitude: '//trim(str(self%RMF_B_amplitude))
	desc = desc//NL//self%mpih%myrankstr//'    RMF rotation axis: '//trim(self%RMF_rotation_axis)
   case default
   desc = desc//NL//self%mpih%myrankstr//'    No external field applied'
   endselect
   endfunction description

   subroutine initialize(self, file_parameters) 
   !< Initialize external fields.
   class(prism_external_fields_object), intent(inout) :: self            !< External fields.
   type(file_ini),                      intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'prism_external_fields_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%description()
   print '(A)', self%mpih%myrankstr//'prism_external_fields_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(prism_external_fields_object), intent(inout)  :: self            !< External fields.
   type(file_ini),           intent(in)                :: file_parameters !< Simulation parameters ini file handler.
   logical,                  intent(in), optional      :: go_on_fail      !< Go on if load fails.
   logical                                             :: go_on_fail_     !< Go on if load fails.
   integer(I4P)                                        :: error           !< Error status.
   character(99)                                       :: buff_char       !< Option character buffer.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='external_fields_applied', &
                            val=buff_char, error=error)
	if (.not.go_on_fail_.and.error>0) &
	call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(external field applied)')
	self%external_field_applied = trim(buff_char)
	self%external_field_applied = trim(self%external_field_applied)

   selectcase(self%external_field_applied)

   case('RMF')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='RMF_frequency', &
   val=self%RMF_frequency, error=error)
   if (.not.go_on_fail_.and.error>0) & 
   call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(RMF_frequency)')

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='RMF_B_amplitude', &
   val=self%RMF_B_amplitude, error=error)
   if (.not.go_on_fail_.and.error>0) & 
   call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(RMF_B_amplitude)')

	call file_parameters%get(section_name=INI_SECTION_NAME, option_name='RMF_rotation_axis', &
                        val=buff_char, error=error)
   if (.not.go_on_fail_.and.error>0) &
   call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(RMF_rotation_axis)')
   self%RMF_rotation_axis = trim(buff_char)
   self%RMF_rotation_axis = trim(self%RMF_rotation_axis)
	select case(self%RMF_rotation_axis)
	case('X', 'x')
		self%alpha = 2_I4P
		self%beta  = 3_I4P
		self%gamma = 1_I4P
	case('Y', 'y')
		self%alpha = 3_I4P
		self%beta  = 1_I4P
		self%gamma = 2_I4P
	case('Z', 'z')
		self%alpha = 1_I4P
		self%beta  = 2_I4P
		self%gamma = 3_I4P
	endselect


   case('MAGNETIC_NOZZLE')


   case('RMF_AND_MAGNETIC_NOZZLE')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='RMF_frequency', &
   val=self%RMF_frequency, error=error)
   if (.not.go_on_fail_.and.error>0) & 
   call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(RMF_frequency)')

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='RMF_B_amplitude', &
   val=self%RMF_B_amplitude, error=error)
   if (.not.go_on_fail_.and.error>0) & 
   call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(RMF_B_amplitude)')

	call file_parameters%get(section_name=INI_SECTION_NAME, option_name='RMF_rotation_axis', &
                        val=buff_char, error=error)
   if (.not.go_on_fail_.and.error>0) &
   call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(RMF_rotation_axis)')
   self%RMF_rotation_axis = trim(buff_char)
   self%RMF_rotation_axis = trim(self%RMF_rotation_axis)
	select case(self%RMF_rotation_axis)
	case('X', 'x')
		self%alpha = 2_I4P
		self%beta  = 3_I4P
		self%gamma = 1_I4P
	case('Y', 'y')
		self%alpha = 3_I4P
		self%beta  = 1_I4P
		self%gamma = 2_I4P
	case('Z', 'z')
		self%alpha = 1_I4P
		self%beta  = 2_I4P
		self%gamma = 3_I4P
	endselect

   case default
   call self%mpih%print_message(msg='no external field applied')
   !self%external_field_applied = 'None'
   !self%external_field_applied = trim(self%external_field_applied)
   endselect

   endsubroutine load_from_file

   subroutine add_external_fields_rmf(self, field, time, dq)
   !< Add rotating magnetic field to the field.
   class(prism_external_fields_object), intent(inout) :: self                                                              !< External fields.
   type(field_object),                  intent(inout) :: field                                                             !< The field.
   real(R8P),                           intent(in)    :: time                                                              !< Current simulation time.
   real(R8P),                           intent(inout) :: dq(1:, 1-field%grid%ngc:,1-field%grid%ngc:,1-field%grid%ngc:,1:)  !< Primitive variables.
   real(R8P)                                          :: x_cell(1-field%grid%ngc:field%grid%ni+field%grid%ngc), &
                                                          y_cell(1-field%grid%ngc:field%grid%nj+field%grid%ngc), &
                                                          z_cell(1-field%grid%ngc:field%grid%nk+field%grid%ngc)             !< Vettori posizione centro celle del blocco b
	real(R8P) 										   :: dB_r, dB_theta 														!< Radial and azimuthal components of the rotating magnetic field
	real(R8P)										   :: theta, alfa, thetaabs									!< Angles in cylindrical coordinates
   integer(I4P)                                       :: b,i,j,k															!< Counters
	real(R8P)										   :: cell_coord(3)														!< Cell coordinates vector

   associate(blocks_number=>field%blocks_number, ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, ngc=>field%grid%ngc, &
		alpha=>self%alpha, beta=>self%beta, gamma=>self%gamma)	

   do b = 1, blocks_number
		call field%grid%cell_xyz(coordinates = field%coordinates(:,b), x_cell = x_cell, y_cell = y_cell, z_cell = z_cell)
         do i = 1, ni
            do j = 1, nj
               do k = 1, nk
					   cell_coord = [x_cell(i), y_cell(j), z_cell(k)]
                  theta = atan(cell_coord(beta)/cell_coord(alpha))
					   thetaabs = abs(atan(cell_coord(beta)/cell_coord(alpha)))
					   alfa = pi/2-thetaabs
					   if (cell_coord(alpha) > 0.0_R8P .and. cell_coord(beta) > 0.0_R8P) then ! 1 quadrante
                     theta = theta
                     dB_r = -(2*PI*self%RMF_frequency)*self%RMF_B_amplitude*sin(2*PI*self%RMF_frequency*time-theta)
	                  dB_theta = (2*PI*self%RMF_frequency)*self%RMF_B_amplitude*cos(2*PI*self%RMF_frequency*time-theta)
					   	dq(alpha+3_I4P,i,j,k,b) = dq(alpha+3_I4P,i,j,k,b) + dB_r*cos(thetaabs) - dB_theta*cos(alfa)
					   	dq(beta+3_I4P,i,j,k,b)  = dq(beta+3_I4P,i,j,k,b)  + dB_r*sin(thetaabs) + dB_theta*sin(alfa)
					   	dq(gamma,i,j,k,b) = -sqrt(cell_coord(alpha)**2 + cell_coord(beta)**2)*(2*PI*self%RMF_frequency)**2* &
					   		self%RMF_B_amplitude*cos(2*PI*self%RMF_frequency*time-theta)*EPS0
					   else if (cell_coord(alpha) < 0.0_R8P .and. cell_coord(beta) > 0.0_R8P) then ! 2 quadrante
                     theta = theta+PI
                     dB_r = -(2*PI*self%RMF_frequency)*self%RMF_B_amplitude*sin(2*PI*self%RMF_frequency*time-theta)
	                  dB_theta = (2*PI*self%RMF_frequency)*self%RMF_B_amplitude*cos(2*PI*self%RMF_frequency*time-theta)
					   	dq(alpha+3_I4P,i,j,k,b) = dq(alpha+3_I4P,i,j,k,b) - dB_r*cos(thetaabs) - dB_theta*cos(alfa)
					   	dq(beta+3_I4P,i,j,k,b)  = dq(beta+3_I4P,i,j,k,b)  + dB_r*sin(thetaabs) - dB_theta*sin(alfa)
					   	dq(gamma,i,j,k,b) = -sqrt(cell_coord(alpha)**2 + cell_coord(beta)**2)*(2*PI*self%RMF_frequency)**2* &
					   		self%RMF_B_amplitude*cos(2*PI*self%RMF_frequency*time-theta)*EPS0
					   else if (cell_coord(alpha) < 0.0_R8P .and. cell_coord(beta) < 0.0_R8P) then ! 3 quadrante
                     theta = theta+PI
                     dB_r = -(2*PI*self%RMF_frequency)*self%RMF_B_amplitude*sin(2*PI*self%RMF_frequency*time-theta)
	                  dB_theta = (2*PI*self%RMF_frequency)*self%RMF_B_amplitude*cos(2*PI*self%RMF_frequency*time-theta)
					   	dq(alpha+3_I4P,i,j,k,b) = dq(alpha+3_I4P,i,j,k,b) - dB_r*cos(thetaabs) + dB_theta*cos(alfa)
					   	dq(beta+3_I4P,i,j,k,b)  = dq(beta+3_I4P,i,j,k,b)  - dB_r*sin(thetaabs) - dB_theta*sin(alfa)
					   	dq(gamma,i,j,k,b) = -sqrt(cell_coord(alpha)**2 + cell_coord(beta)**2)*(2*PI*self%RMF_frequency)**2* &
					   		self%RMF_B_amplitude*cos(2*PI*self%RMF_frequency*time-theta)*EPS0
					   else if (cell_coord(alpha) > 0.0_R8P .and. cell_coord(beta) < 0.0_R8P) then ! 4 quadrante
                     theta = theta+2*PI
                     dB_r = -(2*PI*self%RMF_frequency)*self%RMF_B_amplitude*sin(2*PI*self%RMF_frequency*time-theta)
	                  dB_theta = (2*PI*self%RMF_frequency)*self%RMF_B_amplitude*cos(2*PI*self%RMF_frequency*time-theta)
					   	dq(alpha+3_I4P,i,j,k,b) = dq(alpha+3_I4P,i,j,k,b) + dB_r*cos(thetaabs) + dB_theta*cos(alfa)
					   	dq(beta+3_I4P,i,j,k,b)  = dq(beta+3_I4P,i,j,k,b)  - dB_r*sin(thetaabs) + dB_theta*sin(alfa)
					   	dq(gamma,i,j,k,b) = -sqrt(cell_coord(alpha)**2 + cell_coord(beta)**2)*(2*PI*self%RMF_frequency)**2* &
					   		self%RMF_B_amplitude*cos(2*PI*self%RMF_frequency*time-theta)*EPS0
					   end if
               enddo
            enddo
         enddo
   enddo
	endassociate
   endsubroutine add_external_fields_rmf

   subroutine add_external_fields_none(self, field, time, dq)
   !< Add rotating magnetic field to the field.
   class(prism_external_fields_object), intent(inout) :: self                                                              !< External fields.
   type(field_object),                  intent(inout) :: field                                                             !< The field.
   real(R8P),                           intent(in)    :: time                                                              !< Current simulation time.
   real(R8P),                           intent(inout) :: dq(1:, 1-field%grid%ngc:,1-field%grid%ngc:,1-field%grid%ngc:,1:)   !< Primitive variables.
   real(R8P)                                          :: x_cell(1-field%grid%ngc:field%grid%ni+field%grid%ngc), &
                                                         y_cell(1-field%grid%ngc:field%grid%nj+field%grid%ngc), &
                                                         z_cell(1-field%grid%ngc:field%grid%nk+field%grid%ngc)             !< Vettori posizione centro celle del blocco b
   integer(I4P)                                       :: b,i,j,k

   endsubroutine add_external_fields_none

endmodule adam_prism_external_fields_object