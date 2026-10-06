!< ADAM, PRISM external fields definition, CPU backend.
module adam_prism_external_fields_object
!< ADAM, PRISM external fields definition, CPU backend.

! ADAM singleton objects
use :: adam_mpih_global,  only : mpih
use :: adam_grid_object,  only : grid_object
use :: adam_field_object, only : field_object
use :: adam_prism_bc_object, only : BC_PEC
! PRISM modules
use :: adam_prism_parameters
! third party modules
use :: finer
use :: penf

implicit none
private
public :: INI_SECTION_NAME
public :: EF_TYPE_RMF
public :: EF_TYPE_UNIFORM_FIELD
!public :: EF_TYPE_MAGNETIC_NOZZLE
public :: EF_TYPE_NONE
!public :: EF_TYPE_RMF_AND_MAGNETIC_NOZZLE
public :: prism_external_fields_object
!public :: add_external_fields_interface
!public :: sub_external_fields_interface
!public :: add_external_fields_rmf
!public :: add_external_fields_magnetic_nozzle
!public :: add_external_fields_rmf_and_magnetic_nozzle
!public :: sub_external_fields_rmf
!public :: sub_external_fields_magnetic_nozzle
!public :: sub_external_fields_rmf_and_magnetic_nozzle

character(len=15), parameter :: INI_SECTION_NAME               ='external_fields'        !< INI (config) file section name.
character(len=15), parameter :: EF_TYPE_MAGNETIC_NOZZLE        ='magnetic_nozzle'        !< Magnetic Nozzle.
character(len=4),  parameter :: EF_TYPE_NONE                   ='none'                   !< Disable external field.
character(len=3),  parameter :: EF_TYPE_RMF                    ='RMF'                    !< Rotating Magnetic Field.
character(len=13), parameter :: EF_TYPE_UNIFORM_FIELD          ='Uniform_field'          !< Uniform D and B field.
character(len=23), parameter :: EF_TYPE_RMF_AND_MAGNETIC_NOZZLE='RMF_and_magnetic_nozzle'!< Rotating Magnetic Field/Magnetic Nozzle.

type :: prism_external_fields_object
   !< PRISM external fields object.
   character(len=99) :: ef_type           !< Field type.
   real(R8P)         :: RMF_frequency = 0._R8P !< Rotating magnetic field frequency.
   real(R8P)         :: RMF_B_amplitude = 0._R8P !< Rotating magnetic field amplitude.
	character(len=99) :: RMF_rotation_axis !< Rotating magnetic field rotation axis (X, Y, Z).
   real(R8P)         :: Uniform_D_amplitude = 0._R8P !< Uniform electric displacement field amplitude.
   real(R8P)         :: Uniform_B_amplitude = 0._R8P !< Uniform magnetic field amplitude.
   character(len=99) :: Uniform_direction   = 'x'    !< Uniform field direction (X, Y, Z).
	integer(I4P)      :: alpha             !< RMF rotation axis coordinate 1
	integer(I4P)      :: beta              !< RMF rotation axis coordinate 2
	integer(I4P)      :: gamm              !< RMF rotation axis coordinate 3
   integer(I4P)      :: uniform_axis = 1_I4P !< Uniform field direction index.
   integer(I4P), allocatable :: pec_faces(:,:) !< Physical PEC faces by local block, for FNL kernels.
   ! pointer methods
   procedure(add_external_fields_interface), pass(self), pointer :: add_external_fields=>null() !< Add external fields.
   procedure(sub_external_fields_interface), pass(self), pointer :: sub_external_fields=>null() !< Subtract external fields.
   contains
      ! public methods
      procedure, pass(self) :: description                           !< Return pretty-printed object description.
      procedure, pass(self) :: initialize                            !< Initialize IC.
      procedure, pass(self) :: load_from_file                        !< Load config from file.
      procedure, pass(self) :: add_external_fields_rmf               !< Add rotating magnetic field to the field.
      procedure, pass(self) :: add_external_fields_uniform           !< Add uniform external field to the field.
      !procedure, pass(self) :: add_external_fields_magnetic_nozzle  !< Add magnetic nozzle to the field.
      !procedure, pass(self) :: add_external_fields_rmf_and_magnetic_nozzle !< Add rotating magnetic field and magnetic nozzle to
      !the field.
      procedure, pass(self) :: sub_external_fields_rmf               !< Add rotating magnetic field to the field.
      procedure, pass(self) :: sub_external_fields_uniform           !< Subtract uniform external field from the field.
      !procedure, pass(self) :: sub_external_fields_magnetic_nozzle  !< Add magnetic nozzle to the field.
      !procedure, pass(self) :: sub_external_fields_rmf_and_magnetic_nozzle !< Add rotating magnetic field and magnetic nozzle to
      !the field.
endtype prism_external_fields_object

interface
   subroutine add_external_fields_interface(self, field, grid, time, dt, gamm, q)
   import :: prism_external_fields_object, grid_object, field_object, I4P, R8P
   class(prism_external_fields_object), intent(inout)        :: self                 !< External fields.
   type(field_object),                  intent(inout)        :: field                !< The field.
   type(grid_object),                  intent(in)              :: grid
                                                                                      !< Grid (sibling realm component, threaded
                                                                                      !< in).
   real(R8P),                           intent(in)           :: time                 !< Current simulation time.
   real(R8P),                           intent(in), optional :: dt                   !< Time step.
   real(R8P),                           intent(in), optional :: gamm                 !< Gamma values of RK SSP.
   real(R8P),                           intent(inout)        :: q(1:,1-grid%ngc:,&
                                                                     1-grid%ngc:,&
                                                                     1-grid%ngc:,1:) !< Primitive variables.
   endsubroutine add_external_fields_interface

   subroutine sub_external_fields_interface(self, field, grid, time, dt, gamm, q)
   import :: prism_external_fields_object, grid_object, field_object, I4P, R8P
   class(prism_external_fields_object), intent(inout)        :: self                 !< External fields.
   type(field_object),                  intent(inout)        :: field                !< The field.
   type(grid_object),                  intent(in)              :: grid
                                                                                      !< Grid (sibling realm component, threaded
                                                                                      !< in).
   real(R8P),                           intent(in)           :: time                 !< Current simulation time.
   real(R8P),                           intent(in), optional :: dt                   !< Time step.
   real(R8P),                           intent(in), optional :: gamm                 !< Gamma values of RK SSP.
   real(R8P),                           intent(inout)        :: q(1:,1-grid%ngc:,&
                                                                     1-grid%ngc:,&
                                                                     1-grid%ngc:,1:) !< Primitive variables.
   endsubroutine sub_external_fields_interface
endinterface

contains
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(prism_external_fields_object), intent(in) :: self             !< External fields.
   character(len=:), allocatable                   :: desc             !< Description.
   character(len=1), parameter                     :: NL=new_line('a') !< New line character.
   desc =       mpih%myrankstr//'Applied external fields:'
   select case(self%ef_type)
   case(EF_TYPE_RMF)
   desc = desc//NL//mpih%myrankstr//'    Rotating magnetic field applied '
   desc = desc//NL//mpih%myrankstr//'    RMF frequency: '//trim(str(self%RMF_frequency))
   desc = desc//NL//mpih%myrankstr//'    RMF B amplitude: '//trim(str(self%RMF_B_amplitude))
   desc = desc//NL//mpih%myrankstr//'    RMF rotation axis: '//trim(self%RMF_rotation_axis)
   case(EF_TYPE_UNIFORM_FIELD)
   desc = desc//NL//mpih%myrankstr//'    Uniform external field applied '
   desc = desc//NL//mpih%myrankstr//'    Uniform D amplitude: '//trim(str(self%Uniform_D_amplitude))
   desc = desc//NL//mpih%myrankstr//'    Uniform B amplitude: '//trim(str(self%Uniform_B_amplitude))
   desc = desc//NL//mpih%myrankstr//'    Uniform direction: '//trim(self%Uniform_direction)
   case(EF_TYPE_MAGNETIC_NOZZLE)
   desc = desc//NL//mpih%myrankstr//'    Magnetic nozzle applied '
   case(EF_TYPE_RMF_AND_MAGNETIC_NOZZLE)
   desc = desc//NL//mpih%myrankstr//'    Rotating magnetic field and magnetic nozzle applied '
   desc = desc//NL//mpih%myrankstr//'    RMF frequency: '//trim(str(self%RMF_frequency))
   desc = desc//NL//mpih%myrankstr//'    RMF B amplitude: '//trim(str(self%RMF_B_amplitude))
	desc = desc//NL//mpih%myrankstr//'    RMF rotation axis: '//trim(self%RMF_rotation_axis)
   case default
   desc = desc//NL//mpih%myrankstr//'    No external field applied'
   endselect
   endfunction description

   subroutine initialize(self, file_parameters, field, grid)
   !< Initialize external fields.
   class(prism_external_fields_object), intent(inout) :: self            !< External fields.
   type(file_ini),                      intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   type(field_object),                  intent(in)    :: field
   type(grid_object),                   intent(in)    :: grid
   integer(I4P)                                      :: b
   logical                                           :: faces(6)
   print '(A)', mpih%myrankstr//'prism_external_fields_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   if (allocated(self%pec_faces)) deallocate(self%pec_faces)
   allocate(self%pec_faces(6,field%blocks_number))
   do b=1, field%blocks_number
      call physical_pec_faces(field=field, grid=grid, b=b, pec_faces=faces)
      self%pec_faces(:,b) = merge(1_I4P, 0_I4P, faces)
   enddo

   select case(self%ef_type)
   case(EF_TYPE_RMF)
      self%add_external_fields => add_external_fields_rmf
      self%sub_external_fields => sub_external_fields_rmf
   case(EF_TYPE_UNIFORM_FIELD)
      self%add_external_fields => add_external_fields_uniform
      self%sub_external_fields => sub_external_fields_uniform
   !case(EF_TYPE_MAGNETIC_NOZZLE)
   !   self%add_external_fields => add_external_fields_magnetic_nozzle
   !case(EF_TYPE_RMF_AND_MAGNETIC_NOZZLE)
   !   self%add_external_fields => add_external_fields_rmf_and_magnetic_nozzle
   endselect
   print '(A)', self%description()
   print '(A)', mpih%myrankstr//'prism_external_fields_object%initialize finish'
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
	call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(external field applied)')

   select case(trim(adjustl(buff_char)))
   case('RMF', 'rmf', 'Rmf')
      self%ef_type = EF_TYPE_RMF
   case('Uniform_field', 'uniform_field', 'UNIFORM_FIELD', 'UniformField', 'uniformfield')
      self%ef_type = EF_TYPE_UNIFORM_FIELD
   case('Magnetic_nozzle', 'magnetic_nozzle', 'MAGNETIC_NOZZLE', 'MagneticNozzle', 'magneticnozzle')
      self%ef_type = EF_TYPE_MAGNETIC_NOZZLE
   case('RMF_and_magnetic_nozzle', 'rmf_and_magnetic_nozzle', 'RMF_AND_MAGNETIC_NOZZLE', &
        'Rmf_and_magnetic_nozzle', 'rmfAndMagneticNozzle', 'RMFAndMagneticNozzle')
      self%ef_type = EF_TYPE_RMF_AND_MAGNETIC_NOZZLE
   case default
      self%ef_type = EF_TYPE_NONE
   endselect

   selectcase(self%ef_type)
   case(EF_TYPE_RMF)
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='RMF_frequency', &
      val=self%RMF_frequency, error=error)
      if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(RMF_frequency)')

      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='RMF_B_amplitude', &
      val=self%RMF_B_amplitude, error=error)
      if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(RMF_B_amplitude)')

      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='RMF_rotation_axis', &
                           val=buff_char, error=error)
      if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(RMF_rotation_axis)')
      self%RMF_rotation_axis = trim(buff_char)
      self%RMF_rotation_axis = trim(self%RMF_rotation_axis)
      select case(self%RMF_rotation_axis)
      case('X', 'x')
         self%alpha = 2_I4P
         self%beta  = 3_I4P
         self%gamm  = 1_I4P
      case('Y', 'y')
         self%alpha = 3_I4P
         self%beta  = 1_I4P
         self%gamm  = 2_I4P
      case('Z', 'z')
         self%alpha = 1_I4P
         self%beta  = 2_I4P
         self%gamm  = 3_I4P
      endselect
   case(EF_TYPE_UNIFORM_FIELD)
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='Uniform_field_D_amplitude', &
                               val=self%Uniform_D_amplitude, error=error)
      if (error>0) call file_parameters%get(section_name=INI_SECTION_NAME, option_name='Uniform_D_amplitude', &
                                            val=self%Uniform_D_amplitude, error=error)
      if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(Uniform_field_D_amplitude)')

      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='Uniform_field_B_amplitude', &
                               val=self%Uniform_B_amplitude, error=error)
      if (error>0) call file_parameters%get(section_name=INI_SECTION_NAME, option_name='Uniform_B_amplitude', &
                                            val=self%Uniform_B_amplitude, error=error)
      if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(Uniform_field_B_amplitude)')

      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='Uniform_field_direction', &
                               val=buff_char, error=error)
      if (error>0) call file_parameters%get(section_name=INI_SECTION_NAME, option_name='Uniform_direction', &
                                            val=buff_char, error=error)
      if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(Uniform_field_direction)')
      self%Uniform_direction = trim(buff_char)
      select case(trim(adjustl(self%Uniform_direction)))
      case('X', 'x')
         self%uniform_axis = 1_I4P
      case('Y', 'y')
         self%uniform_axis = 2_I4P
      case('Z', 'z')
         self%uniform_axis = 3_I4P
      case default
         call mpih%error_stop(msg=': invalid ['//INI_SECTION_NAME//'].(Uniform_direction), expected x, y, or z')
      endselect
   case(EF_TYPE_MAGNETIC_NOZZLE)

   case(EF_TYPE_RMF_AND_MAGNETIC_NOZZLE)
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='RMF_frequency', &
      val=self%RMF_frequency, error=error)
      if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(RMF_frequency)')

      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='RMF_B_amplitude', &
      val=self%RMF_B_amplitude, error=error)
      if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(RMF_B_amplitude)')

      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='RMF_rotation_axis', &
                           val=buff_char, error=error)
      if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(RMF_rotation_axis)')
      self%RMF_rotation_axis = trim(buff_char)
      self%RMF_rotation_axis = trim(self%RMF_rotation_axis)
      select case(self%RMF_rotation_axis)
      case('X', 'x')
         self%alpha = 2_I4P
         self%beta  = 3_I4P
         self%gamm  = 1_I4P
      case('Y', 'y')
         self%alpha = 3_I4P
         self%beta  = 1_I4P
         self%gamm  = 2_I4P
      case('Z', 'z')
         self%alpha = 1_I4P
         self%beta  = 2_I4P
         self%gamm  = 3_I4P
      endselect
   case default
      call mpih%print_message(msg='no external field applied')
   endselect
   endsubroutine load_from_file

   !subroutine add_external_fields_rmf(self, field, time, dt, gamm, dq)
   !!< Add rotating magnetic field to the field.
   !class(prism_external_fields_object), intent(inout)           :: self                                                            
   ! !< External fields.
   !type(field_object),                  intent(inout)           :: field                                                           
   ! !< The field.
   !real(R8P),                           intent(in)              :: time                                                            
   ! !< Current simulation time.
   !real(R8P),                           intent(in), optional    :: dt                                                              
   ! !< Time step.
   !real(R8P),                           intent(in), optional    :: gamm                                                            
   ! !< Gamma values of RK SSP
   !real(R8P),                           intent(inout)           :: dq(1:, 1-field%grid%ngc:,1-field%grid%ngc:,1-field%grid%ngc:,1:)
   ! !< Primitive variables.
   !real(R8P)                                                    :: x_cell(1-field%grid%ngc:field%grid%ni+field%grid%ngc), &
   !                                                                y_cell(1-field%grid%ngc:field%grid%nj+field%grid%ngc), &
   !                                                                z_cell(1-field%grid%ngc:field%grid%nk+field%grid%ngc)           
   !  !< Vettori posizione centro celle del blocco b
 !real(R8P) 										                      :: dB_r, dB_theta 														          !< Radial and azimuthal components of the
 !rotating magnetic field
   !real(R8P)										                      :: time1															                !< Time at the next sub-step
	!real(R8P)										                      :: theta                  								                   !< Angle in cylindrical coordinates
   !integer(I4P)                                                 :: b,i,j,k															                !< Counters
 !real(R8P)										                      :: cell_coord(3)												                   !< Cell coordinates vector and scalar
 !variables
   !real(R8P)                                                    :: x, y, r, omega, phase, c, s
   !associate(blocks_number=>field%blocks_number, ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, ngc=>field%grid%ngc, &
	!	alpha=>self%alpha, beta=>self%beta, gamma=>self%gamma)
   !if (present(gamm)) then
   !   time1 = time + dt*gamm
   !else
   !   time1 = time
   !end if
   !do b = 1, blocks_number
	!	call field%grid%cell_xyz(coordinates = field%coordinates(:,b), x_cell = x_cell, y_cell = y_cell, z_cell = z_cell)
   !      do i = 1, ni
   !         do j = 1, nj
   !            do k = 1, nk
	!				   cell_coord = [x_cell(i), y_cell(j), z_cell(k)]
   !               x = cell_coord(alpha)
   !               y = cell_coord(beta)
   !               r = sqrt(x*x + y*y)
   !               theta = atan2(y, x)
   !               omega = 2.0_R8P*PI*self%RMF_frequency
   !               phase = omega*time1 - theta
   !               dB_r     = -omega*self%RMF_B_amplitude*sin(phase)
   !               dB_theta =  omega*self%RMF_B_amplitude*cos(phase)
   !               c = cos(theta)
   !               s = sin(theta)
   !               dq(alpha+3_I4P,i,j,k,b) = dB_r*c - dB_theta*s
   !               dq(beta +3_I4P,i,j,k,b) = dB_r*s + dB_theta*c
   !               dq(gamma,i,j,k,b) = -r*omega*omega*self%RMF_B_amplitude*sin(phase)*EPS0
   !            enddo
   !         enddo
   !      enddo
   !enddo
	!endassociate
   !endsubroutine add_external_fields_rmf

   subroutine add_external_fields_rmf(self, field, grid, time, dt, gamm, q)
   !< Add rotating magnetic field to the field.
   class(prism_external_fields_object), intent(inout)           :: self                 !< External fields.
   type(field_object),                  intent(inout)           :: field                !< The field.
   type(grid_object),                  intent(in)              :: grid
                                                                                      !< Grid (sibling realm component, threaded
                                                                                      !< in).
   real(R8P),                           intent(in)              :: time                 !< Current simulation time.
   real(R8P),                           intent(in), optional    :: dt                   !< Time step.
   real(R8P),                           intent(in), optional    :: gamm                 !< Gamma values of RK SSP
   real(R8P),                           intent(inout)           :: q(1:,1-grid%ngc:,&
                                                                        1-grid%ngc:,&
                                                                        1-grid%ngc:,1:) !< Primitive variables.
	real(R8P) 										                      :: B_r, B_theta 			 !< Radial and azimuthal components of the rotating magnetic field
   real(R8P)										                      :: time1					 !< Time at the next sub-step
	real(R8P)										                      :: theta                !< Angle in cylindrical coordinates
   integer(I4P)                                                 :: b,i,j,k					 !< Counters
   integer(I4P)                                                 :: i_sample, j_sample, k_sample
   logical                                                      :: pec_faces(6)
   real(R8P)                                                   :: sign_D(3), sign_B(3)
	real(R8P)										                      :: cell_coord(3)			 !< Cell coordinates vector and scalar variables
   real(R8P)                                                    :: x, y, r, omega, phase, c, s

   associate(blocks_number=>field%blocks_number, ni=>field%ni, nj=>grid%nj, nk=>grid%nk, ngc=>grid%ngc, &
		alpha=>self%alpha, beta=>self%beta, ef_gamma=>self%gamm, x_cell=>field%x_cell, y_cell=>field%y_cell, z_cell=>field%z_cell)
   if (present(gamm)) then
      time1 = time + dt*gamm
   else
      time1 = time + dt
   end if
   do b = 1, blocks_number
         call physical_pec_faces(field=field, grid=grid, b=b, pec_faces=pec_faces)
         do i = 1 - ngc, ni + ngc
            do j = 1 - ngc, nj + ngc
               do k = 1 - ngc, nk + ngc
                  call reflect_pec_ghost(i=i, j=j, k=k, ni=ni, nj=nj, nk=nk, pec_faces=pec_faces, &
                                         i_sample=i_sample, j_sample=j_sample, k_sample=k_sample, &
                                         sign_D=sign_D, sign_B=sign_B)
					   cell_coord = [x_cell(i_sample,b), y_cell(j_sample,b), z_cell(k_sample,b)]
                  x = cell_coord(alpha)
                  y = cell_coord(beta)
                  r = sqrt(x*x + y*y)
                  theta = atan2(y, x)
                  omega = 2.0_R8P*PI*self%RMF_frequency
                  phase = omega*time1 - theta
                  B_r     = self%RMF_B_amplitude*cos(phase)
                  B_theta = self%RMF_B_amplitude*sin(phase)
                  c = cos(theta)
                  s = sin(theta)
                  q(alpha+3_I4P,i,j,k,b) = q(alpha+3_I4P,i,j,k,b) + sign_B(alpha)*(B_r*c - B_theta*s)
                  q(beta +3_I4P,i,j,k,b) = q(beta +3_I4P,i,j,k,b) + sign_B(beta )*(B_r*s + B_theta*c)
                  q(ef_gamma,i,j,k,b) = q(ef_gamma,i,j,k,b) + &
                                         sign_D(ef_gamma)*r*omega*self%RMF_B_amplitude*cos(phase)*EPS0
               enddo
            enddo
         enddo
   enddo
	endassociate
   endsubroutine add_external_fields_rmf


   subroutine add_external_fields_uniform(self, field, grid, time, dt, gamm, q)
   !< Add uniform external electric displacement and magnetic field to the field.
   class(prism_external_fields_object), intent(inout)           :: self                 !< External fields.
   type(field_object),                  intent(inout)           :: field                !< The field.
   type(grid_object),                   intent(in)              :: grid                 !< Grid.
   real(R8P),                           intent(in)              :: time                 !< Current simulation time.
   real(R8P),                           intent(in), optional    :: dt                   !< Time step.
   real(R8P),                           intent(in), optional    :: gamm                 !< Gamma values of RK SSP.
   real(R8P),                           intent(inout)           :: q(1:,1-grid%ngc:,&
                                                                        1-grid%ngc:,&
                                                                        1-grid%ngc:,1:) !< Primitive variables.
   integer(I4P)                                                 :: b,i,j,k              !< Counters.
   integer(I4P)                                                 :: i_sample, j_sample, k_sample
   logical                                                      :: pec_faces(6)
   real(R8P)                                                   :: sign_D(3), sign_B(3)

   associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj, nk=>grid%nk, ngc=>grid%ngc, axis=>self%uniform_axis)
   if (present(dt)) continue
   if (present(gamm)) continue
   associate(time_unused=>time)
   endassociate
   do b = 1, blocks_number
      call physical_pec_faces(field=field, grid=grid, b=b, pec_faces=pec_faces)
      do i = 1 - ngc, ni + ngc
         do j = 1 - ngc, nj + ngc
            do k = 1 - ngc, nk + ngc
               call reflect_pec_ghost(i=i, j=j, k=k, ni=ni, nj=nj, nk=nk, pec_faces=pec_faces, &
                                      i_sample=i_sample, j_sample=j_sample, k_sample=k_sample, &
                                      sign_D=sign_D, sign_B=sign_B)
               q(axis,       i,j,k,b) = q(axis,       i,j,k,b) + sign_D(axis)*self%Uniform_D_amplitude
               q(axis+3_I4P, i,j,k,b) = q(axis+3_I4P, i,j,k,b) + sign_B(axis)*self%Uniform_B_amplitude
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine add_external_fields_uniform

   subroutine sub_external_fields_rmf(self, field, grid, time, dt, gamm, q)
   !< Add rotating magnetic field to the field.
   class(prism_external_fields_object), intent(inout)           :: self                 !< External fields.
   type(field_object),                  intent(inout)           :: field                !< The field.
   type(grid_object),                  intent(in)              :: grid
                                                                                      !< Grid (sibling realm component, threaded
                                                                                      !< in).
   real(R8P),                           intent(in)              :: time                 !< Current simulation time.
   real(R8P),                           intent(in), optional    :: dt                   !< Time step.
   real(R8P),                           intent(in), optional    :: gamm                 !< Gamma values of RK SSP
   real(R8P),                           intent(inout)           :: q(1:,1-grid%ngc:,&
                                                                        1-grid%ngc:,&
                                                                        1-grid%ngc:,1:) !< Primitive variables.
	real(R8P) 										                      :: B_r, B_theta 			 !< Radial and azimuthal components of the rotating magnetic field
   real(R8P)										                      :: time1					 !< Time at the next sub-step
	real(R8P)										                      :: theta                !< Angle in cylindrical coordinates
   integer(I4P)                                                 :: b,i,j,k					 !< Counters
   integer(I4P)                                                 :: i_sample, j_sample, k_sample
   logical                                                      :: pec_faces(6)
   real(R8P)                                                   :: sign_D(3), sign_B(3)
	real(R8P)										                      :: cell_coord(3)			 !< Cell coordinates vector and scalar variables
   real(R8P)                                                    :: x, y, r, omega, phase, c, s

   associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj, nk=>grid%nk, ngc=>grid%ngc, &
		alpha=>self%alpha, beta=>self%beta, ef_gamma=>self%gamm, x_cell=>field%x_cell, y_cell=>field%y_cell, z_cell=>field%z_cell)
   if (present(gamm)) then
      time1 = time + dt*gamm
   else
      time1 = time
   end if
   do b = 1, blocks_number
         call physical_pec_faces(field=field, grid=grid, b=b, pec_faces=pec_faces)
         do i = 1 - ngc, ni + ngc
            do j = 1 - ngc, nj + ngc
               do k = 1 - ngc, nk + ngc
                  call reflect_pec_ghost(i=i, j=j, k=k, ni=ni, nj=nj, nk=nk, pec_faces=pec_faces, &
                                         i_sample=i_sample, j_sample=j_sample, k_sample=k_sample, &
                                         sign_D=sign_D, sign_B=sign_B)
					   cell_coord = [x_cell(i_sample,b), y_cell(j_sample,b), z_cell(k_sample,b)]
                  x = cell_coord(alpha)
                  y = cell_coord(beta)
                  r = sqrt(x*x + y*y)
                  theta = atan2(y, x)
                  omega = 2.0_R8P*PI*self%RMF_frequency
                  phase = omega*time1 - theta
                  B_r     = self%RMF_B_amplitude*cos(phase)
                  B_theta = self%RMF_B_amplitude*sin(phase)
                  c = cos(theta)
                  s = sin(theta)
                  q(alpha+3_I4P,i,j,k,b) = q(alpha+3_I4P,i,j,k,b) - sign_B(alpha)*(B_r*c - B_theta*s)
                  q(beta +3_I4P,i,j,k,b) = q(beta +3_I4P,i,j,k,b) - sign_B(beta )*(B_r*s + B_theta*c)
                  q(ef_gamma,i,j,k,b) = q(ef_gamma,i,j,k,b) - &
                                         sign_D(ef_gamma)*r*omega*self%RMF_B_amplitude*cos(phase)*EPS0
               enddo
            enddo
         enddo
   enddo
	endassociate
   endsubroutine sub_external_fields_rmf

   subroutine sub_external_fields_uniform(self, field, grid, time, dt, gamm, q)
   !< Subtract uniform external electric displacement and magnetic field from the field.
   class(prism_external_fields_object), intent(inout)           :: self                 !< External fields.
   type(field_object),                  intent(inout)           :: field                !< The field.
   type(grid_object),                   intent(in)              :: grid                 !< Grid.
   real(R8P),                           intent(in)              :: time                 !< Current simulation time.
   real(R8P),                           intent(in), optional    :: dt                   !< Time step.
   real(R8P),                           intent(in), optional    :: gamm                 !< Gamma values of RK SSP.
   real(R8P),                           intent(inout)           :: q(1:,1-grid%ngc:,&
                                                                        1-grid%ngc:,&
                                                                        1-grid%ngc:,1:) !< Primitive variables.
   integer(I4P)                                                 :: b,i,j,k              !< Counters.
   integer(I4P)                                                 :: i_sample, j_sample, k_sample
   logical                                                      :: pec_faces(6)
   real(R8P)                                                   :: sign_D(3), sign_B(3)

   associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj, nk=>grid%nk, ngc=>grid%ngc, axis=>self%uniform_axis)
   if (present(dt)) continue
   if (present(gamm)) continue
   associate(time_unused=>time)
   endassociate
   do b = 1, blocks_number
      call physical_pec_faces(field=field, grid=grid, b=b, pec_faces=pec_faces)
      do i = 1 - ngc, ni + ngc
         do j = 1 - ngc, nj + ngc
            do k = 1 - ngc, nk + ngc
               call reflect_pec_ghost(i=i, j=j, k=k, ni=ni, nj=nj, nk=nk, pec_faces=pec_faces, &
                                      i_sample=i_sample, j_sample=j_sample, k_sample=k_sample, &
                                      sign_D=sign_D, sign_B=sign_B)
               q(axis,       i,j,k,b) = q(axis,       i,j,k,b) - sign_D(axis)*self%Uniform_D_amplitude
               q(axis+3_I4P, i,j,k,b) = q(axis+3_I4P, i,j,k,b) - sign_B(axis)*self%Uniform_B_amplitude
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine sub_external_fields_uniform

   pure subroutine physical_pec_faces(field, grid, b, pec_faces)
   !< Identify physical PEC faces of one block; internal block ghosts keep the unreflected field.
   type(field_object), intent(in)  :: field
   type(grid_object),  intent(in)  :: grid
   integer(I4P),      intent(in)  :: b
   logical,           intent(out) :: pec_faces(6)
   real(R8P)                      :: boundary_tol
   integer(I4P)                   :: axis

   boundary_tol = 16._R8P*epsilon(1._R8P)*max(1._R8P, maxval(abs(grid%domain_emin)), &
                                                    maxval(abs(grid%domain_emax)))
   do axis = 1, 3
      pec_faces(2*axis-1) = grid%bc_type(2*axis-1) == BC_PEC .and. &
                            abs(field%emin(axis,b)-grid%domain_emin(axis)) <= boundary_tol
      pec_faces(2*axis  ) = grid%bc_type(2*axis  ) == BC_PEC .and. &
                            abs(field%emax(axis,b)-grid%domain_emax(axis)) <= boundary_tol
   enddo
   endsubroutine physical_pec_faces

   pure subroutine reflect_pec_ghost(i, j, k, ni, nj, nk, pec_faces, i_sample, j_sample, k_sample, sign_D, sign_B)
   !< Apply the same PEC parity as prism_cpu_object%set_boundary_conditions to an external field.
   integer(I4P), intent(in)  :: i, j, k, ni, nj, nk
   logical,      intent(in)  :: pec_faces(6)
   integer(I4P), intent(out) :: i_sample, j_sample, k_sample
   real(R8P),    intent(out) :: sign_D(3), sign_B(3)
   integer(I4P)             :: cell_index(3), cell_count(3), sample_index(3), axis
   logical                  :: reflected

   cell_index = [i, j, k]
   cell_count = [ni, nj, nk]
   sample_index = cell_index
   sign_D = 1._R8P
   sign_B = 1._R8P
   do axis = 1, 3
      reflected = .false.
      if (cell_index(axis) < 1 .and. pec_faces(2*axis-1)) then
         sample_index(axis) = 1 - cell_index(axis)
         reflected = .true.
      elseif (cell_index(axis) > cell_count(axis) .and. pec_faces(2*axis)) then
         sample_index(axis) = 2*cell_count(axis) + 1 - cell_index(axis)
         reflected = .true.
      endif
      if (reflected) then
         sign_D = -sign_D
         sign_D(axis) = -sign_D(axis)
         sign_B(axis) = -sign_B(axis)
      endif
   enddo
   i_sample = sample_index(1)
   j_sample = sample_index(2)
   k_sample = sample_index(3)
   endsubroutine reflect_pec_ghost

endmodule adam_prism_external_fields_object
