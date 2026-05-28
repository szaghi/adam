!< ADAM, PRISM coil source definition, CPU backend.
module adam_prism_coil_object
!< ADAM, PRISM coil source definition, CPU backend.

! ADAM singleton objects
use :: adam_field_object, only : field_object
use :: adam_grid_object,  only : grid_object
use :: adam_mpih_global,  only : mpih
! PRISM modules
use :: adam_prism_physics_object, only : prism_physics_object, ADIM_EM_PHYSICAL_MODEL
use :: adam_prism_parameters
! third party modules
use :: finer
use :: penf
use :: stringifor

implicit none
private
public :: INI_SECTION_NAME
public :: COIL_TYPE_RECTANGULAR
public :: COIL_TYPE_CIRCULAR
public :: COIL_TYPE_SOLENOID
public :: CURRENT_TYPE_AC
public :: CURRENT_TYPE_DC
public :: prism_coil_object
public :: NORMAL_P_X
public :: NORMAL_M_X
public :: NORMAL_P_Y
public :: NORMAL_M_Y
public :: NORMAL_P_Z
public :: NORMAL_M_Z

character(len=11), parameter :: INI_SECTION_NAME="coils_input"        !< INI (config) file section name containing coils configs.
character(len=11), parameter :: COIL_TYPE_RECTANGULAR="rectangular"   !< Rectangular shape coil.
character(len=8),  parameter :: COIL_TYPE_CIRCULAR="circular"         !< Circular shape coil.
character(len=8),  parameter :: COIL_TYPE_SOLENOID="solenoid"         !< Solenoid shape coil.
character(len=10), parameter :: CURRENT_TYPE_DC="DC_current"          !< DC current.
character(len=10), parameter :: CURRENT_TYPE_AC="AC_current"          !< AC current
character(len=2),  parameter :: NORMAL_P_X="+x"                       !< Normal versor along positive x axis.
character(len=2),  parameter :: NORMAL_M_X="-x"                       !< Normal versor along negative x axis.
character(len=2),  parameter :: NORMAL_P_Y="+y"                       !< Normal versor along positive y axis.
character(len=2),  parameter :: NORMAL_M_Y="-y"                       !< Normal versor along negative y axis.
character(len=2),  parameter :: NORMAL_P_Z="+z"                       !< Normal versor along positive z axis.
character(len=2),  parameter :: NORMAL_M_Z="-z"                       !< Normal versor along negative z axis.

type :: prism_coil_object
   !< ADAM, PRISM coil source definition, CPU backend.
   character(len=99), allocatable :: coil_type(:)                          !< Coil type.
   character(len=99), allocatable :: current_type(:)                       !< Current type.
   character(len=2 ), allocatable :: normal(:)                             !< Versore normale alla spira, che identifica anche verso
   real(R8P),         allocatable :: A(:)                                  !< Current amplitude (A)
   real(R8P),         allocatable :: coil_amplitude(:)
                                                                           !< Current amplitude (A) for each coil, corrected for
                                                                           !< Gaussian current distribution.
   real(R8P),         allocatable :: f(:)                                  !< Current frequency, if AC (Hz)
   real(R8P),         allocatable :: phase(:)                              !< Current initial phase, if AC
   real(R8P),         allocatable :: x_center(:), y_center(:), z_center(:) !< Coil center
   real(R8P),         allocatable :: lx(:), ly(:)                          !< Rectangle's sizes (if rectangular coil)
   real(R8P),         allocatable :: r_coil(:)                             !< Circle's radius (if circular coil)
   real(R8P),         allocatable :: l_solenoid(:)                         !< Solenoid length (if solenoidal coil)
   real(R8P),         allocatable :: windings(:)                           !< Windings number (if solenoidal coil)
   real(R8P),         allocatable :: sigma(:)                              !< Gaussian current distribution sigma
   real(R8P),         allocatable :: J_vec(:,:,:,:,:,:)                    !< Matrice contenente versori corrente spire (se assente
   real(R8P)                      :: td                                    !< Delay di accensione della spira
   integer(I4P)                   :: circular_coils_number=0_I4P           !< Number of circular coils
   integer(I4P)                   :: rectangular_coils_number=0_I4P        !< Number of rectangular coils
   integer(I4P)                   :: total_coils_number=0_I4P              !< Number of coils
   type(string), allocatable      :: j_vec_name(:,:)                       !< J vec names.
   ! grid data replica for easy handling
   integer(I4P), pointer :: ngc=>null() !< Number of ghost cells.
   integer(I4P), pointer :: ni=>null()  !< Number of cells in i direction.
   integer(I4P), pointer :: nj=>null()  !< Number of cells in j direction.
   integer(I4P), pointer :: nk=>null()  !< Number of cells in k direction.
   contains
      ! public methods
      procedure, pass(self) :: allocate_coil                             !< Allocate coil data.
      procedure, pass(self) :: description                               !< Return pretty-printed object description.
      procedure, pass(self) :: initialize                                !< Initialize IC.
      procedure, pass(self) :: load_from_file                            !< Load config from file.
      procedure, pass(self) :: adimensionalize_coils_parameters          !< Adimensionalize coils parameters.
      !procedure, pass(self) :: set_coils                                 !< Set coil_object on PRISM fields.
      !procedure, pass(self) :: set_circular_coil                         !< Set circular coils on PRISM fields.
      !procedure, pass(self) :: set_rectangular_coil_quad_section_odd     !< Set rectangular coils on PRISM fields with quadratic
      !section and odd number of cell
      !procedure, pass(self) :: set_rectangular_coil_quad_section_odd_v2  !< Set rectangular coils on PRISM fields with quadratic
      !section and odd number of cell
      !procedure, pass(self) :: set_rectangular_coil_quad_section_even    !< Set rectangular coils on PRISM fields with quadratic
      !section and even number of cell
      !procedure, pass(self) :: set_rectangular_coil_quad_section_even_v2 !< Set rectangular coils on PRISM fields with quadratic
      !section and even number of cell
      !procedure, pass(self) :: set_rectangular_coil_circular_section     !< Set rectangular coils on PRISM fields with circular
      !section

endtype prism_coil_object

contains
   ! public methods
   subroutine allocate_coil(self, field, grid)
   !< Allocate coil data.
   class(prism_coil_object), intent(inout) :: self !< Coils.
   type(field_object), intent(in) :: field !< Field (sibling realm component, threaded in).
   type(grid_object),  intent(in), target :: grid !< Grid (sibling realm component, threaded in).
   integer(I4P)                            :: c    !< Counter.

   associate(ngc=>self%ngc,ni=>self%ni,nj=>self%nj,nk=>self%nk,nb=>field%nb,total_coils_number=>self%total_coils_number)

   allocate(self%r_coil              (0:total_coils_number)) ; self%r_coil = 0.0_R8P
   allocate(self%ly                  (0:total_coils_number)) ; self%ly = 0.0_R8P
   allocate(self%lx                  (0:total_coils_number)) ; self%lx = 0.0_R8P
   allocate(self%l_solenoid          (0:total_coils_number)) ; self%l_solenoid = 0.0_R8P
   allocate(self%windings            (0:total_coils_number)) ; self%windings = 0.0_R8P
   allocate(self%sigma               (0:total_coils_number)) ; self%sigma = 0.0_R8P
   allocate(self%x_center            (0:total_coils_number)) ; self%x_center = 0.0_R8P
   allocate(self%y_center            (0:total_coils_number)) ; self%y_center = 0.0_R8P
   allocate(self%z_center            (0:total_coils_number)) ; self%z_center = 0.0_R8P
   allocate(self%coil_type           (0:total_coils_number)) ; self%coil_type = ' '
   allocate(self%current_type        (0:total_coils_number)) ; self%current_type = ' '
   allocate(self%normal              (0:total_coils_number)) ; self%normal = ' '
   allocate(self%A                   (0:total_coils_number)) ; self%A = 0.0_R8P
   allocate(self%coil_amplitude      (0:total_coils_number)) ; self%coil_amplitude = 0.0_R8P
   allocate(self%f                   (0:total_coils_number)) ; self%f = 0.0_R8P
   allocate(self%phase               (0:total_coils_number)) ; self%phase = 0.0_R8P

   allocate(self%J_vec(3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb,total_coils_number)) ; self%J_vec = 0._R8P
   allocate(self%j_vec_name(3,total_coils_number))

   do c=1, self%total_coils_number
      self%j_vec_name(1,c) = 'coil_'//trim(strz(c,2))//'_j_vec_1'
      self%j_vec_name(2,c) = 'coil_'//trim(strz(c,2))//'_j_vec_2'
      self%j_vec_name(3,c) = 'coil_'//trim(strz(c,2))//'_j_vec_3'
   enddo

   endassociate
   endsubroutine allocate_coil

   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(prism_coil_object), intent(in) :: self             !< IC.
   character(len=:), allocatable        :: desc             !< Description.
   character(len=1), parameter          :: NL=new_line('a') !< New line character.
   integer(I4P)                         :: r                !< Counter.

   desc =       mpih%myrankstr//'Coils main data'
   if (self%total_coils_number > 0_I4P) then
      do r=1, self%total_coils_number
         desc = desc//NL//mpih%myrankstr//'  Coil('//trim(str(r,.true.))//')'
         select case(self%coil_type(r))
         case(COIL_TYPE_CIRCULAR)
         desc = desc//NL//mpih%myrankstr//'    Coil type: '//trim(self%coil_type(r))
         desc = desc//NL//mpih%myrankstr//'    Current type: '//trim(self%current_type(r))
         desc = desc//NL//mpih%myrankstr//'    Radius: '//trim(str(self%r_coil(r)))
         desc = desc//NL//mpih%myrankstr//'    Normal: '//trim(self%normal(r))
         desc = desc//NL//mpih%myrankstr//'    X_center: '//trim(str(self%x_center(r)))
         desc = desc//NL//mpih%myrankstr//'    Y_center: '//trim(str(self%y_center(r)))
         desc = desc//NL//mpih%myrankstr//'    Z_center: '//trim(str(self%z_center(r)))
         desc = desc//NL//mpih%myrankstr//'    Amplitude: '//trim(str(self%A(r)))
         desc = desc//NL//mpih%myrankstr//'    Frequency: '//trim(str(self%f(r)))
         desc = desc//NL//mpih%myrankstr//'    Phase: '//trim(str(self%phase(r)))
         desc = desc//NL//mpih%myrankstr//'    Sigma: '//trim(str(self%sigma(r)))
         case(COIL_TYPE_RECTANGULAR)
         desc = desc//NL//mpih%myrankstr//'    Coil type: '//trim(self%coil_type(r))
         desc = desc//NL//mpih%myrankstr//'    Current type: '//trim(self%current_type(r))
         desc = desc//NL//mpih%myrankstr//'    L1: '//trim(str(self%lx(r)))
         desc = desc//NL//mpih%myrankstr//'    L2: '//trim(str(self%ly(r)))
         desc = desc//NL//mpih%myrankstr//'    Normal: '//trim(self%normal(r))
         desc = desc//NL//mpih%myrankstr//'    X_center: '//trim(str(self%x_center(r)))
         desc = desc//NL//mpih%myrankstr//'    Y_center: '//trim(str(self%y_center(r)))
         desc = desc//NL//mpih%myrankstr//'    Z_center: '//trim(str(self%z_center(r)))
         desc = desc//NL//mpih%myrankstr//'    Amplitude: '//trim(str(self%A(r)))
         desc = desc//NL//mpih%myrankstr//'    Frequency: '//trim(str(self%f(r)))
         desc = desc//NL//mpih%myrankstr//'    Phase: '//trim(str(self%phase(r)))
         desc = desc//NL//mpih%myrankstr//'    Sigma: '//trim(str(self%sigma(r)))
         case(COIL_TYPE_SOLENOID)
         desc = desc//NL//mpih%myrankstr//'    Coil type: '//trim(self%coil_type(r))
         desc = desc//NL//mpih%myrankstr//'    Current type: '//trim(self%current_type(r))
         desc = desc//NL//mpih%myrankstr//'    Radius: '//trim(str(self%r_coil(r)))
         desc = desc//NL//mpih%myrankstr//'    Length: '//trim(str(self%l_solenoid(r)))
         desc = desc//NL//mpih%myrankstr//'    Windings: '//trim(str(self%windings(r)))
         desc = desc//NL//mpih%myrankstr//'    Normal: '//trim(self%normal(r))
         desc = desc//NL//mpih%myrankstr//'    X_center: '//trim(str(self%x_center(r)))
         desc = desc//NL//mpih%myrankstr//'    Y_center: '//trim(str(self%y_center(r)))
         desc = desc//NL//mpih%myrankstr//'    Z_center: '//trim(str(self%z_center(r)))
         desc = desc//NL//mpih%myrankstr//'    Amplitude: '//trim(str(self%A(r)))
         desc = desc//NL//mpih%myrankstr//'    Frequency: '//trim(str(self%f(r)))
         desc = desc//NL//mpih%myrankstr//'    Phase: '//trim(str(self%phase(r)))
         desc = desc//NL//mpih%myrankstr//'    Sigma: '//trim(str(self%sigma(r)))
         endselect
      enddo
   else
      desc = desc//mpih%myrankstr//'  No coils defined.'
   endif
   endfunction description

   subroutine initialize(self, field, grid, physics, file_parameters)
                                                             !Cfr ic%initialize, ma commentata parte descrizione perchè da
                                                             !implementare
   !< Initialize the equation.
   class(prism_coil_object),   intent(inout)      :: self            !< Coils.
   type(field_object),         intent(in)         :: field !< Field (sibling realm component, threaded in).
   type(grid_object),          intent(in), target :: grid !< Grid (sibling realm component, threaded in).
   type(prism_physics_object), intent(in)         :: physics         !< PRISM physics.
   type(file_ini),             intent(in)         :: file_parameters !< Simulation parameters ini file handler.

   print '(A)', mpih%myrankstr//'prism_coil_object%initialize start'
   call self%load_from_file(field=field, physics=physics, grid=grid, file_parameters=file_parameters)
   print '(A)', self%description()
   print '(A)', mpih%myrankstr//'prism_coil_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, field, grid, physics, file_parameters, go_on_fail)
   !< Load config from file.
   class(prism_coil_object),   intent(inout)        :: self            !< coils.
   type(field_object),         intent(in)           :: field           !< Field (sibling realm component, threaded in).
   type(grid_object),          intent(in), target   :: grid            !< Grid (sibling realm component, threaded in).
   type(prism_physics_object), intent(in)           :: physics         !< PRISM physics.
   type(file_ini),             intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,                    intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                          :: go_on_fail_     !< Go on if load fails.
   character(:), allocatable                        :: sname           !< Section name.
   integer(I4P)                                     :: i               !< Counter.
   integer(I4P)                                     :: error           !< Error status.
   character(99)                                    :: buff_char       !< Option character buffer.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='circular_coils_number', &
                            val=self%circular_coils_number, error=error)
   if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(circular_coils_number)')

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='rectangular_coils_number', &
                            val=self%rectangular_coils_number, error=error)
   if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(rectangular_coils_number)')

   self%total_coils_number = self%circular_coils_number + self%rectangular_coils_number

   if (self%total_coils_number==0_I4P) return

   call associate_adam_data

   call self%allocate_coil(field=field, grid=grid)

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='time_delay', val=self%td, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(time_delay)')

   do i=1, self%total_coils_number
      sname = INI_SECTION_NAME//'_coil_'//trim(str(i,.true.))

      call file_parameters%get(section_name=sname, option_name='coil_type', val=buff_char, error=error)
      if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(coil_type)')
      self%coil_type(i) = trim(buff_char)
      self%coil_type(i) = trim(self%coil_type(i))

      call file_parameters%get(section_name=sname, option_name='current_type', val=buff_char, error=error)
      if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(current_type)')
      self%current_type(i) = trim(buff_char)
      self%current_type(i) = trim(self%current_type(i))

      call file_parameters%get(section_name=sname, option_name='sigma', val=self%sigma(i), error=error)
      if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(sigma)')

      select case(self%coil_type(i))
      case(COIL_TYPE_CIRCULAR)
         call file_parameters%get(section_name=sname, option_name='r_coil', val=self%r_coil(i), error=error)
         if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(r_coil)')

         call file_parameters%get(section_name=sname, option_name='x_center', val=self%x_center(i), error=error)
         if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(x_center)')

         call file_parameters%get(section_name=sname, option_name='y_center', val=self%y_center(i), error=error)
         if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(y_center)')

         call file_parameters%get(section_name=sname, option_name='z_center', val=self%z_center(i), error=error)
         if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(z_center)')

         call file_parameters%get(section_name=sname, option_name='normal', val=buff_char, error=error)
         if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(normal)')
         self%normal(i) = trim(buff_char)
         self%normal(i) = trim(self%normal(i))

      case(COIL_TYPE_RECTANGULAR)
         call file_parameters%get(section_name=sname, option_name='lx', val=self%lx(i), error=error)
         if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(lx)')

         call file_parameters%get(section_name=sname, option_name='ly', val=self%ly(i), error=error)
         if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(ly)')

         call file_parameters%get(section_name=sname, option_name='x_center', val=self%x_center(i), error=error)
         if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(x_center)')

         call file_parameters%get(section_name=sname, option_name='y_center', val=self%y_center(i), error=error)
         if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(y_center)')

         call file_parameters%get(section_name=sname, option_name='z_center', val=self%z_center(i), error=error)
         if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(z_center)')

         call file_parameters%get(section_name=sname, option_name='normal', val=buff_char, error=error)
         if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(normal)')
         self%normal(i) = trim(buff_char)
         self%normal(i) = trim(self%normal(i))

      case(COIL_TYPE_SOLENOID)
         call file_parameters%get(section_name=sname, option_name='r_coil', val=self%r_coil(i), error=error)
         if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(r_coil)')

         call file_parameters%get(section_name=sname, option_name='l_solenoid', val=self%l_solenoid(i), error=error)
         if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(l_solenoid)')

         call file_parameters%get(section_name=sname, option_name='windings', val=self%windings(i), error=error)
         if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(windings)')

         call file_parameters%get(section_name=sname, option_name='x_center', val=self%x_center(i), error=error)
         if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(x_center)')

         call file_parameters%get(section_name=sname, option_name='y_center', val=self%y_center(i), error=error)
         if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(y_center)')

         call file_parameters%get(section_name=sname, option_name='z_center', val=self%z_center(i), error=error)
         if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(z_center)')

         call file_parameters%get(section_name=sname, option_name='normal', val=buff_char, error=error)
         if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(normal)')
         self%normal(i) = trim(buff_char)
         self%normal(i) = trim(self%normal(i))
      endselect

      select case(self%current_type(i))
      case(CURRENT_TYPE_DC)

         call file_parameters%get(section_name=sname, option_name='Amplitude', val=self%A(i), error=error)
         if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(Amplitude)')

      case(CURRENT_TYPE_AC)

         call file_parameters%get(section_name=sname, option_name='Amplitude', val=self%A(i), error=error)
         if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(Amplitude)')

         call file_parameters%get(section_name=sname, option_name='Frequency', val=self%f(i), error=error)
         if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(Frequency)')

         call file_parameters%get(section_name=sname, option_name='Phase', val=self%phase(i), error=error)
         if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(Phase)')

      endselect
   enddo

   if (physics%physical_model == ADIM_EM_PHYSICAL_MODEL) &
      call self%adimensionalize_coils_parameters(physics=physics)
   
   contains
      subroutine associate_adam_data
      !< Associate grid data pointers for easy handling.

      self%ni  => grid%ni
      self%nj  => grid%nj
      self%nk  => grid%nk
      self%ngc => grid%ngc
      endsubroutine associate_adam_data
   endsubroutine load_from_file

   subroutine adimensionalize_coils_parameters(self, physics)
   !< Adimensionalize coils parameters.
   class(prism_coil_object),   intent(inout) :: self    !< Coils
   type(prism_physics_object), intent(in)    :: physics !< PRISM physics.

   !Adimensionalization spatial parameters
   self%x_center(:)   = self%x_center(:)/physics%L0
   self%y_center(:)   = self%y_center(:)/physics%L0
   self%z_center(:)   = self%z_center(:)/physics%L0
   self%lx(:)         = self%lx(:)/physics%L0
   self%ly(:)         = self%ly(:)/physics%L0
   self%r_coil(:)     = self%r_coil(:)/physics%L0
   self%l_solenoid(:) = self%l_solenoid(:)/physics%L0
   self%sigma(:)      = self%sigma(:)/physics%L0
   !Adimensionalization of current parameters
   self%A(:) = self%A(:)/physics%J0
   !Adimensionalization of frequency parameters
   self%f(:) = self%f(:)*physics%T0
   !Adimensionalization of time parameters
   self%td = self%td/physics%T0
   endsubroutine adimensionalize_coils_parameters
endmodule adam_prism_coil_object
