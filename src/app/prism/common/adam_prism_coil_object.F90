!< ADAM, PRISM coil source definition, CPU backend.
module adam_prism_coil_object
!< ADAM, PRISM coil source definition, CPU backend.

! ADAM singleton objects
use :: adam_global_field, only : field
use :: adam_global_grid,  only : grid
use :: adam_global_mpih,  only : mpih
! PRISM modules
use :: adam_prism_physics_object, only : prism_physics_object
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
   integer(I4P),      allocatable :: coil_flag(:,:,:,:)                    !< Matrice contenente informazioni su quale spira pass pe
   integer(I4P)                   :: circular_coils_number=0_I4P           !< Number of circular coils
   integer(I4P)                   :: rectangular_coils_number=0_I4P        !< Number of rectangular coils
   integer(I4P)                   :: total_coils_number=0_I4P              !< Number of coils
   type(string), allocatable      :: j_vec_name(:,:)                       !< J vec names.
   type(string)                   :: coil_flag_name                        !< Coil flag name.
   ! grid data replica for easy handling
   integer(I4P), pointer :: ngc=>null()           !< Number of ghost cells.
   integer(I4P), pointer :: ni=>null()            !< Number of cells in i direction.
   integer(I4P), pointer :: nj=>null()            !< Number of cells in j direction.
   integer(I4P), pointer :: nk=>null()            !< Number of cells in k direction.
   contains
      ! public methods
      procedure, pass(self) :: allocate_coil                             !< Allocate coil data.
      procedure, pass(self) :: description                               !< Return pretty-printed object description.
      procedure, pass(self) :: initialize                                !< Initialize IC.
      procedure, pass(self) :: load_from_file                            !< Load config from file.
      !procedure, pass(self) :: set_coils                                 !< Set coil_object on PRISM fields.
      !procedure, pass(self) :: set_circular_coil                         !< Set circular coils on PRISM fields.
      !procedure, pass(self) :: set_rectangular_coil_quad_section_odd     !< Set rectangular coils on PRISM fields with quadratic section and odd number of cell
      !procedure, pass(self) :: set_rectangular_coil_quad_section_odd_v2  !< Set rectangular coils on PRISM fields with quadratic section and odd number of cell
      !procedure, pass(self) :: set_rectangular_coil_quad_section_even    !< Set rectangular coils on PRISM fields with quadratic section and even number of cell
      !procedure, pass(self) :: set_rectangular_coil_quad_section_even_v2 !< Set rectangular coils on PRISM fields with quadratic section and even number of cell
      !procedure, pass(self) :: set_rectangular_coil_circular_section     !< Set rectangular coils on PRISM fields with circular section

endtype prism_coil_object

contains
   ! public methods
   subroutine allocate_coil(self)
   !< Allocate coil data.
   class(prism_coil_object), intent(inout) :: self !< Coils.
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
   allocate(self%f                   (0:total_coils_number)) ; self%f = 0.0_R8P
   allocate(self%phase               (0:total_coils_number)) ; self%phase = 0.0_R8P

   allocate(self%coil_flag(1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb)) ; self%coil_flag = 0_I4P
   allocate(self%J_vec(4,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb,total_coils_number)) ; self%J_vec = 0._R8P
   self%coil_flag_name = 'coil_flag'
   allocate(self%j_vec_name(4,total_coils_number))
   do c=1, self%total_coils_number
      self%j_vec_name(1,c) = 'coil_'//trim(strz(c,2))//'_j_vec_1'
      self%j_vec_name(2,c) = 'coil_'//trim(strz(c,2))//'_j_vec_2'
      self%j_vec_name(3,c) = 'coil_'//trim(strz(c,2))//'_j_vec_3'
      self%j_vec_name(4,c) = 'coil_'//trim(strz(c,2))//'_f_Gauss'
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

   subroutine initialize(self, file_parameters) !Cfr ic%initialize, ma commentata parte descrizione perchè da implementare
   !< Initialize the equation.
   class(prism_coil_object), intent(inout) :: self            !< Coils.
   type(file_ini),           intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   print '(A)', mpih%myrankstr//'prism_coil_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%description()
   print '(A)', mpih%myrankstr//'prism_coil_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(prism_coil_object), intent(inout)        :: self            !< coils.
   type(file_ini),           intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,                  intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                        :: go_on_fail_     !< Go on if load fails.
   character(:), allocatable                      :: sname           !< Section name.
   integer(I4P)                                   :: i               !< Counter.
   integer(I4P)                                   :: error           !< Error status.
   character(99)                                  :: buff_char       !< Option character buffer.

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

   call self%allocate_coil

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
   contains
      subroutine associate_adam_data
      !< Associate grid data pointers for easy handling.

      self%ni  => grid%ni
      self%nj  => grid%nj
      self%nk  => grid%nk
      self%ngc => grid%ngc
      endsubroutine associate_adam_data
   endsubroutine load_from_file

endmodule adam_prism_coil_object

   !subroutine set_coils(self, physics, field)
   !!< Set initial conditions on PRISM fields.
   !class(prism_coil_object),     intent(inout) :: self    !< Coils
   !type(field_object),           intent(inout) :: field   !< Field object.
   !type(prism_physics_object),   intent(in)    :: physics !< Fluids physics.
   !integer(I4P)                                :: i       !< Counter.
   !if (self%total_coils_number >= 1_I4P) then
   !   do i=1, self%total_coils_number
   !      select case(self%coil_type(i))
   !      case(COIL_TYPE_CIRCULAR) !Caso spire circolari
   !         call self%set_circular_coil (physics = physics, field = field, n = i)
   !      case(COIL_TYPE_RECTANGULAR) !Caso spire rettangolari
   !         selectcase(self%coil_sec(i))
   !         case(CIRC_COIL_SECTION) !Caso spire rettangolari infinite
   !            call self%set_rectangular_coil_circular_section(physics = physics, field = field, n = i)
   !         case(QUAD_COIL_SECTION) !Caso spire rettangolari finite
   !            if (mod(self%d(i), 2.0_R8P) == 0.0_R8P) then
   !               call self%set_rectangular_coil_quad_section_even_v2(physics = physics, field = field, n = i)
   !            else
   !               call self%set_rectangular_coil_quad_section_odd_v2(physics = physics, field = field, n = i)
   !            endif
   !         endselect
   !      endselect
   !   enddo
   !endif
   !endsubroutine set_coils
   !
   !subroutine set_circular_coil(self, physics, field, n) !agli input aggiungo n del contatore per sapere a quale
   !                                                      !spira faccio riferimento
   !   !< Set coils on PRISM fields. La subroutine restituirà il vettore q contenuto in fields
   !   !< completo anche dei valori normalizzati delle correnti che passano per le celle (elementi 7,8,9)
   !   !< da calcolare poi tramite la funzione che assegna il valore della corrente compute_coils_current
   !   class(prism_coil_object),     intent(inout) :: self                                                                !< Coils
   !   type(field_object),           intent(inout) :: field                                                               !< Field object.
   !   type(prism_physics_object),   intent(in)    :: physics                                                             !< Fluids physiscs.
   !   integer(I4P),                 intent(in)    :: n                                                                   !< Coil number.
   !   !real(R8P),                    allocatable   :: flag(:,:,:,:)                                                      !< Flag per identificare se la spira passa per la cella
   !   real(R8P),                    allocatable   :: Gaussian(:,:,:,:)                                                   !< Matrice gaussiana per distribuzione corrente
   !   real(R8P)                                   :: dmax, dist                                                                !< Vincolo distanza massima dalla spira.
   !   real(R8P)                                   :: c_c(3)                                                              !< Vettore posizione centro spira
   !   real(R8P)                                   :: cell_coord(3)                                                       !< Vettore posizione centro cella
   !   integer(I4P)                                :: b,i,j,k                                                             !< Counter.
   !   !associo per dati su posizioni delle celle e contatori
   !   associate(blocks_number=>field%blocks_number, ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, ngc=>field%grid%ngc, &
   !            x_c => self%x_center(n), y_c => self%y_center(n), z_c => self%z_center(n), &
   !            dx => field%dxyz(1,:), dy => field%dxyz(2,:), dz => field%dxyz(3,:), r_coil => self%r_coil(n), &
   !            normal => self%normal(:,n), d => self%d(n), nb=>field%nb,current_distribution => self%current_distribution(n), &
   !            x_cell => field%x_cell, y_cell => field%y_cell, z_cell => field%z_cell, sigma => self%sigma(n))
!
   !   c_c = [ x_c, y_c, z_c ] !Vettore posizione centro spira
!
   !   !allocate(flag(ni,nj,nk,blocks_number))
   !   do b=1, blocks_number
   !      ! chiamo funzione che restituisce le coordinate delle varie celle che compongono il blocco
   !      !calcolo distanza massima dall'asse del filo della spira: somma di raggio del filo e metà della dimensione
   !      !massima della cella associata ai vettori dx dy e dz contenuti in field
   !      dmax = d/2
   !      do k=1, nk
   !         do j=1, nj
   !            do i=1, ni
   !               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]
   !               if ((dotproduct(a=(cell_coord-c_c),b=normal))**2 + (sqrt(sq_norm(cell_coord-c_c) - &
   !                  (dotproduct(a=(cell_coord-c_c),b=normal))**2) - r_coil)**2 <= (d/2)**2 .and. &
   !                  self%coil_flag(i,j,k,b) == 0_I4P) then
!
   !                  self%J_vec(1:3,i,j,k,b) = crossproduct(a=normal,b=(cell_coord-c_c))
!
   !                  !normalizzo per ottenere, alla fine il versore della corrente nella cella
   !                  !q(7:9,i,j,k,b) = q(7:9,i,j,k,b)/sqrt(sq_norm(q(7:9,i,j,k,b)))
   !                  self%J_vec(1:3,i,j,k,b) = self%J_vec(1:3,i,j,k,b)/sqrt(sq_norm(self%J_vec(1:3,i,j,k,b)))
!
   !                  !metto flag su quale spira passa per la cella
   !                  self%coil_flag(i,j,k,b) = n
!
   !                  selectcase (current_distribution)
   !                  case (GAUSS_CURRENT_DISTRIBUTION)
   !                     dist = sqrt((dotproduct(a=(cell_coord-c_c),b=normal))**2 + (sqrt(sq_norm(cell_coord-c_c) - &
   !                            (dotproduct(a=(cell_coord-c_c),b=normal))**2) - r_coil)**2)
   !                     self%J_vec(4,i,j,k,b) = gaussian_2D_ind(sigma = d/6, r = dist)
   !                  case (CONST_CURRENT_DISTRIBUTION)
   !                     self%J_vec(4,i,j,k,b) = 4.0/(PI*d**2)
   !                  endselect
   !               endif
   !            enddo
   !         enddo
   !      enddo
   !   enddo
   !   endassociate
   !endsubroutine set_circular_coil
!
   !subroutine set_rectangular_coil_quad_section_odd(self, physics, field, n)
   !class(prism_coil_object),     intent(inout) :: self                                                            !< Coils
   !type(field_object),           intent(inout) :: field                                                           !< Field object.
   !type(prism_physics_object),   intent(in)    :: physics                                                         !< Fluids physiscs.
   !integer(I4P),                 intent(in)    :: n                                                               !< Coil number.
   !integer(I4P),                 allocatable   :: flag(:,:,:,:)                                                   !< Flag per identificare se la spira passa per la cella
   !real(R8P),                    allocatable   :: Gaussian(:,:,:,:)                                               !< Matrice gaussiana per distribuzione corrente
   !real(R8P)                                   :: dmax                                                            !< Vincolo distanza massima dalla spira.
   !real(R8P)                                   :: c_c(3)                                                          !< Vettore posizione centro spira
   !real(R8P)                                   :: cell_coord(3)                                                   !< Vettore posizione centro cella
   !real(R8P)                                   :: vx(3), vy(3), vz(3)                                             !< Versori assi cartesiani
   !real(R8P)                                   :: V1(3), V2(3), V3(3), V4(3), V(4,3)                              !< Vertici rettangolo e relativa matrice
   !real(R8P)                                   :: V1_1(3), V2_1(3), V3_1(3), V4_1(3)
   !real(R8P)                                   :: V1_2(3), V2_2(3), V3_2(3), V4_2(3)
   !real(R8P)                                   :: v_l1(3), v_l2(3), vec(4,3)                                      !< Versori lati rettangolo (vale regola mano dx) e relativa matrice
   !real(R8P)                                   :: n1(3), d1, n2(3), d2, n3(3), d3, n4(3), d4                      !< Parametri piani diagonali perpendicolari a spira, per evitare sovrapposizioni
   !real(R8P)                                   :: kappa(3), K_rot(3,3), Id(3,3), theta                            !< Vettore, angolo e matrice di appoggio per formula Rodrigues + matrice identità
   !real(R8P)                                   :: Kquad(3,3), R(3,3)                                              !< Matrice K^2 e matrice rotazione tra vz e normale alla spira
   !real(R8P)                                   :: dist, prj_v(3),eps, d_real                                              !< Distanza punto retta e proiezione del punto sulla retta                                                                          !< Variabile utilizzata per definire direzione corrente
   !integer(I4P)                                :: b,i,j,k,w                                                       !< Counter.
   !integer(I4P)                                :: i1,j1,k1,w1                                                     !< Counter1.
   !integer(I4P)                                :: d_int
   !!associo per dati su posizioni delle celle e contatori
   !associate(blocks_number=>field%blocks_number, ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, ngc=>field%grid%ngc, &
   !         x_c => self%x_center(n), y_c => self%y_center(n), z_c => self%z_center(n),                                         &
   !         dx => field%dxyz(1,:), dy => field%dxyz(2,:), dz => field%dxyz(3,:), lx => self%lx(n), ly => self%ly(n),           &
   !         normal => self%normal(:,n), d => self%d(n), nb =>field%nb, current_distribution => self%current_distribution(n),   &
   !         x_cell => field%x_cell, y_cell => field%y_cell, z_cell => field%z_cell, sigma => self%sigma(n))
   !c_c = [ x_c, y_c, z_c ] !Vettore posizione centro spira
   !!vertici del rettangolo, lo costruisco come se avesse normale asse z e fosse centrato nell'origine;
   !!vertici in senso antiorario, partendo da in basso a sinistra; si ipotizza
   !!rotazione rispetto al centro, quindi, ad esempio, se la vuoi con normale lungo x
   !!con lato lungo lungo y AB sarà il lato corto (pensa a rotazione 3D). Con
   !!normale lungo y il problema non si pone
   !vx = [1._R8P,0._R8P,0._R8P]
   !vy = [0._R8P,1._R8P,0._R8P]
   !vz = [0._R8P,0._R8P,1._R8P]
   !V1 = [-lx/2, -ly/2, 0._R8P]
   !V2 = [+lx/2, -ly/2, 0._R8P]
   !V3 = [+lx/2, +ly/2, 0._R8P]
   !V4 = [-lx/2, +ly/2, 0._R8P]
   !! 1: interno
   !V1_1 = [-(lx/2 - (d-1.0)/2*dx(1)), -(ly/2 - (d-1.0)/2*dx(1)), 0._R8P]
   !V2_1 = [(lx/2  - (d-1.0)/2*dx(1)), -(ly/2 - (d-1.0)/2*dx(1)), 0._R8P]
   !V3_1 = [(lx/2  - (d-1.0)/2*dx(1)),  (ly/2 - (d-1.0)/2*dx(1)), 0._R8P]
   !V4_1 = [-(lx/2 - (d-1.0)/2*dx(1)), +(ly/2 - (d-1.0)/2*dx(1)), 0._R8P]
   !! 2: esterno
   !V1_2 = [-(lx/2 + (d-1.0)/2*dx(1)), -(ly/2 + (d-1.0)/2*dx(1)), 0._R8P]
   !V2_2 = [(lx/2  + (d-1.0)/2*dx(1)), -(ly/2 + (d-1.0)/2*dx(1)), 0._R8P]
   !V3_2 = [(lx/2  + (d-1.0)/2*dx(1)),  (ly/2 + (d-1.0)/2*dx(1)), 0._R8P]
   !V4_2 = [-(lx/2 + (d-1.0)/2*dx(1)), +(ly/2 + (d-1.0)/2*dx(1)), 0._R8P]
!
   !!            V4 ________________ V3     A
   !!              |                |       A
   !!              |                |       |
   !!              |                |       |
   !!              |                |       y
   !!              |________________|
   !!            V1                  V2     x --->>>
!
   !!calcolo rotazione tra i due vettori normali %OSS una rotazione di 180° dà problemi, qua ci metti un go to.
   !!isnan per il caso rotazione nulla rispetto
   !!a n // z
   !kappa = crossproduct(a=vz,b=normal)
!
   !if (any(normal /= vz )) then
   !    kappa = kappa/sqrt(sq_norm(kappa))
   !endif
!
   !theta = acos(dotproduct(a=vz,b=normal))
!
   !K_rot(1,1) = 0._R8P
   !K_rot(1,2) = -kappa(3)
   !K_rot(1,3) = kappa(2)
   !K_rot(2,1) = kappa(3)
   !K_rot(2,2) = 0._R8P
   !K_rot(2,3) = -kappa(1)
   !K_rot(3,1) = -kappa(2)
   !K_rot(3,2) = kappa(1)
   !K_rot(3,3) = 0._R8P
   !Kquad = matmul(K_rot,K_rot) !matrice K^2
   !Id(:,1) = vx
   !Id(:,2) = vy
   !Id(:,3) = vz
!
   !R = Id+sin(theta)*K_rot+(1-cos(theta))*Kquad !costruisco matrice di rotazione, ruoto i vettori posizione e poi traslo
   !                                            ! a seconda della posizione del centro richiesta
!
   !V1 = matmul(R,V1)+c_c
   !V2 = matmul(R,V2)+c_c
   !V3 = matmul(R,V3)+c_c
   !V4 = matmul(R,V4)+c_c
   !V(1,:) = V1
   !V(2,:) = V2
   !V(3,:) = V3
   !V(4,:) = V4 !genero matrice [V1; V2; V3; V4], quindi le colonne sono le coordinate x y z e le righe i vari vertici nell'ordine descritto prima
!
   !V1_1 = matmul(R,V1_1)+c_c
   !V2_1 = matmul(R,V2_1)+c_c
   !V3_1 = matmul(R,V3_1)+c_c
   !V4_1 = matmul(R,V4_1)+c_c
   !V1_2 = matmul(R,V1_2)+c_c
   !V2_2 = matmul(R,V2_2)+c_c
   !V3_2 = matmul(R,V3_2)+c_c
   !V4_2 = matmul(R,V4_2)+c_c
   !v_l1 = matmul(R,vx)
   !v_l1 = v_l1/sqrt(sq_norm(v_l1))
   !v_l2 = matmul(R,vy);
   !v_l2 = v_l2/sqrt(sq_norm(v_l2))
!
   !!matrice dei versori dei lati, generata ruotando i versori del rettangolo tramite la matrice di rotazione precedentemente calcolata
   !vec(1,:) = v_l1
   !vec(2,:) = v_l2
   !vec(3,:) = -v_l1
   !vec(4,:) = -v_l2
!
   !!calcolo piani perpendicolari alla spira su cui giacciono le due diagonali
   !!del rettangolo. Per convenzione, le normali puntano "verso i vertici" D e C
   !!piano 1, diagonale V1-V3 e normale verso V4
   !n1 = crossproduct(a=(V1_2-V1_1),b=normal)/sqrt(sq_norm(V1_2-V1_1)) !normale al piano
   !d1 = -dotproduct(a=n1,b=V1_1) !parametro d dell'equazione ax + by + cz + d1 = 0, con a b c coseni direttori della normale
   !n2 = crossproduct(a=(V2_1-V2_2),b=normal)/sqrt(sq_norm(V2_1-V2_2)) !normale al piano
   !d2 = -dotproduct(a=n2,b=V2_1) !parametro d dell'equazione ax + by + cz + d1 = 0, con a b c coseni direttori della normale
   !n3 = crossproduct(a=(V3_1-V3_2),b=normal)/sqrt(sq_norm(V3_1-V3_2)) !normale al piano
   !d3 = -dotproduct(a=n3,b=V3_1) !parametro d dell'equazione ax + by + cz + d1 = 0, con a b c coseni direttori della normale
   !n4 = crossproduct(a=(V4_2-V4_1),b=normal)/sqrt(sq_norm(V4_2-V4_1)) !normale al piano
   !d4 = -dotproduct(a=n4,b=V4_1) !parametro d dell'equazione ax + by + cz + d1 = 0, con a b c coseni direttori della normale
!
   !allocate(flag(1:ni,1:nj,1:nk,1:nb))
   !flag(:,:,:,:) = 0_I4P !inizializzo matrice flag a zero, per indicare che nessun lato passa per le celle
   !allocate(Gaussian(1:ni,1:nj,1:nk,1:nb))
   !Gaussian(:,:,:,:) = 0_I4P !inizializzo matrice gaussiana a zero
   !d_int = int(d/2)
   !d_real = real(d_int,R8P)
   !eps = 1*10e-10
   !!modificata per avere termine sorgente come Filippo, se infittiamo o aumentiamo sezione spira torna la precedente
   !do w = 1, 4, 2 !per ogni lato del rettangolo
   !   do b=1, blocks_number
   !      dmax = dx(b)*d_real+eps
   !      ! chiamo funzione che restituisce le coordinate delle varie celle che compongono il blocco
   !      do k=1+d_int, nk-d_int
   !         do j=1+d_int, nj-d_int
   !            do i=1+d_int, ni-d_int
   !               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]
   !               dist = sqrt(sq_norm(crossproduct(a=(cell_coord-V(w,:)),b=vec(w,:)))) !Distanza punto retta |(P-A) x v| / |v| con A punto sulla retta e v versore della retta
   !               prj_v = V(w,:)+dotproduct(a=(cell_coord-V(w,:)),b=vec(w,:))*vec(w,:); !Formula proiezione di un punto su una retta con A punto della retta A+[(P-A)*v]*v
   !               if (dist <= dx(b)/2 .and. &
   !                  prj_v(1) <= maxval(V(:,1)) + dmax .and. &
   !                  prj_v(2) <= maxval(V(:,2)) + dmax .and. prj_v(3) <= maxval(V(:,3)) + dmax .and. &
   !                  minval(V(:,1))-dmax <= prj_v(1) .and. minval(V(:,2))-dmax <= prj_v(2) .and. &
   !                  minval(V(:,3))-dmax <= prj_v(3) ) then
   !                  if (d_int /= 0_I4P) then
   !                     do i1 = -d_int, d_int
   !                        do j1 = -d_int, d_int
   !                           do k1 = -d_int, d_int
   !                              if (flag(i+i1,j+j1,k+k1,b) == 0_I4P) then
   !                                 flag(i+i1,j+j1,k+k1,b) = w
   !                              endif
   !                           enddo
   !                        enddo
   !                     enddo
   !                  else
   !                     if (flag(i,j,k,b) == 0_I4P) then
   !                        flag(i,j,k,b) = w
   !                     endif
   !                  endif
!
   !               endif
   !            enddo
   !         enddo
   !      enddo
   !   enddo
   !   do b=1, blocks_number
   !      do k = 1, nk
   !         do j = 1, nj
   !            do i = 1, ni
   !               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]
   !               dist = sqrt(sq_norm(crossproduct(a=(cell_coord-V(w,:)),b=vec(w,:))))
   !                !secondo if: per ogni lato verifico di essere dal "lato giusto" dei piani definiti dalle diagonali, al fine di
   !                !non avere sovrapposizioni in prossimità dei vertici
   !               if (w == 1) then
   !                  if ((dotproduct(a=n1,b=cell_coord)+d1 > eps .or. dotproduct(a=n2,b=cell_coord)+d2 > eps) .and. &
   !                       flag(i,j,k,b) == w) then !aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
   !                     flag(i,j,k,b) = 0_I4P
   !                     !print *, w
   !                  endif
   !               elseif (w == 3) then
   !                  if ((dotproduct(a=n3,b=cell_coord)+d3 < -eps .or. dotproduct(a=n4,b=cell_coord)+d4 < -eps) .and. &
   !                       flag(i,j,k,b) == w) then !aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
   !                     flag(i,j,k,b) = 0_I4P
   !                     !print *, w
   !                  endif
   !               endif
   !               if (flag(i,j,k,b) == w) then
   !                  selectcase (current_distribution)
   !                  case (GAUSS_CURRENT_DISTRIBUTION)
   !                     Gaussian(i,j,k,b) = gaussian_2D_ind(sigma = sigma*dx(b), r = dist)
   !                  case (CONST_CURRENT_DISTRIBUTION)
   !                     Gaussian(i,j,k,b) = 1/(d*dx(b))**2
   !                  endselect
   !               endif
   !            enddo
   !         enddo
   !      enddo
   !   enddo
   !enddo
!
   !do w = 2, 4, 2
   !   do b=1, blocks_number
   !      dmax = dx(b)*d_real+eps
   !      ! chiamo funzione che restituisce le coordinate delle varie celle che compongono il blocco
   !      do k=1+d_int, nk-d_int
   !         do j=1+d_int, nj-d_int
   !            do i=1+d_int, ni-d_int
   !               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]
   !               dist = sqrt(sq_norm(crossproduct(a=(cell_coord-V(w,:)),b=vec(w,:)))) !Distanza punto retta |(P-A) x v| / |v| con A punto sulla retta e v versore della retta
   !               prj_v = V(w,:)+dotproduct(a=(cell_coord-V(w,:)),b=vec(w,:))*vec(w,:); !Formula proiezione di un punto su una retta con A punto della retta A+[(P-A)*v]*v
   !               if (dist <= dx(b)/2 .and. &
   !               prj_v(1) <= maxval(V(:,1)) + dmax .and. &
   !               prj_v(2) <= maxval(V(:,2)) + dmax .and. prj_v(3) <= maxval(V(:,3)) + dmax .and. &
   !               minval(V(:,1))-dmax <= prj_v(1) .and. minval(V(:,2))-dmax <= prj_v(2) .and. &
   !               minval(V(:,3))-dmax <= prj_v(3) ) then
   !                  if (d_int /= 0_I4P) then
   !                     do i1 = -d_int, d_int
   !                        do j1 = -d_int, d_int
   !                           do k1 = -d_int, d_int
   !                              if (flag(i+i1,j+j1,k+k1,b) == 0_I4P) then
   !                                 flag(i+i1,j+j1,k+k1,b) = w
   !                              endif
   !                           enddo
   !                        enddo
   !                     enddo
   !                  else
   !                     if (flag(i,j,k,b) == 0_I4P) then
   !                        flag(i,j,k,b) = w
   !                     endif
   !                  endif
   !               endif
   !            enddo
   !         enddo
   !      enddo
   !   enddo
   !   do b=1, blocks_number
   !      do k = 1, nk
   !         do j = 1, nj
   !            do i = 1, ni
   !               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]
   !               dist = sqrt(sq_norm(crossproduct(a=(cell_coord-V(w,:)),b=vec(w,:))))
   !               if (w == 2) then
   !                  if ((dotproduct(a=n2,b=cell_coord)+d2 <= -eps .or. dotproduct(a=n3,b=cell_coord)+d3 >= eps) .and. &
   !                       flag(i,j,k,b) == w) then!aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
   !                     flag(i,j,k,b) = 0_I4P
   !                  endif
   !               elseif (w == 4) then
   !                  if ((dotproduct(a=n1,b=cell_coord)+d1 <= -eps .or. dotproduct(a=n4,b=cell_coord)+d4 >= eps) .and. &
   !                       flag(i,j,k,b) == w) then !aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
   !                     flag(i,j,k,b) = 0_I4P
   !                  endif
   !               endif
   !               if (flag(i,j,k,b) == w) then
   !                  selectcase (current_distribution)
   !                  case (GAUSS_CURRENT_DISTRIBUTION)
   !                     Gaussian(i,j,k,b) = gaussian_2D_ind(sigma = sigma*dx(b), r = dist)
   !                  case (CONST_CURRENT_DISTRIBUTION)
   !                     Gaussian(i,j,k,b) = 1/(d*dx(b))**2
   !                  endselect
   !               endif
   !            enddo
   !         enddo
   !      enddo
   !   enddo
   !enddo
   !!Ho un flag pari a 1 2 3 4 nelle celle per cui passa uno dei dati della spira. La direzione della corrente è
   !!Ceorente con quella dei versori dei lati precedentemente descritti
   !do b=1, blocks_number
   !   do k=1, nk
   !      do j=1, nj
   !         do i=1, ni
   !            if (flag(i,j,k,b) /= 0 .and. self%coil_flag(i,j,k,b) == 0_I4P) then
!
   !               self%J_vec(1:3,i,j,k,b) = vec(flag(i,j,k,b),:)
   !               self%J_vec(4,i,j,k,b) = Gaussian(i,j,k,b)
!
   !               if (abs(self%J_vec(1,i,j,k,b)) < 1.0e-10_R8P) then
   !                  self%J_vec(1,i,j,k,b) = 0._R8P
   !               endif
   !               if (abs(self%J_vec(2,i,j,k,b)) < 1.0e-10_R8P) then
   !                  self%J_vec(2,i,j,k,b) = 0._R8P
   !               endif
   !               if (abs(self%J_vec(3,i,j,k,b)) < 1.0e-10_R8P) then
   !                  self%J_vec(3,i,j,k,b) = 0._R8P
   !               endif
!
   !               self%coil_flag(i,j,k,b) = n
!
   !            endif
   !         enddo
   !      enddo
   !   enddo
   !enddo
   !endassociate
   !endsubroutine set_rectangular_coil_quad_section_odd
!
   !subroutine set_rectangular_coil_quad_section_even(self, physics, field, n)
   !class(prism_coil_object),     intent(inout) :: self                                                            !< Coils
   !type(field_object),           intent(inout) :: field                                                           !< Field object.
   !type(prism_physics_object),   intent(in)    :: physics                                                         !< Fluids physiscs.
   !integer(I4P),                 intent(in)    :: n                                                               !< Coil number.
   !integer(I4P),                 allocatable   :: flag(:,:,:,:)                                                   !< Flag per identificare se la spira passa per la cella
   !real(R8P),                    allocatable   :: Gaussian(:,:,:,:)                                               !< Matrice gaussiana per distribuzione corrente
   !real(R8P)                                   :: dmax                                                            !< Vincolo distanza massima dalla spira.
   !real(R8P)                                   :: c_c(3)                                                          !< Vettore posizione centro spira
   !real(R8P)                                   :: cell_coord(3)                                                   !< Vettore posizione centro cella
   !real(R8P)                                   :: vx(3), vy(3), vz(3)                                             !< Versori assi cartesiani
   !real(R8P)                                   :: V1(3), V2(3), V3(3), V4(3), V(4,3)                              !< Vertici rettangolo e relativa matrice
   !real(R8P)                                   :: V1_1(3), V2_1(3), V3_1(3), V4_1(3)
   !real(R8P)                                   :: V1_2(3), V2_2(3), V3_2(3), V4_2(3)
   !real(R8P)                                   :: v_l1(3), v_l2(3), vec(4,3)                                      !< Versori lati rettangolo (vale regola mano dx) e relativa matrice
   !real(R8P)                                   :: n1(3), d1, n2(3), d2, n3(3), d3, n4(3), d4                      !< Parametri piani diagonali perpendicolari a spira, per evitare sovrapposizioni
   !real(R8P)                                   :: kappa(3), K_rot(3,3), Id(3,3), theta                            !< Vettore, angolo e matrice di appoggio per formula Rodrigues + matrice identità
   !real(R8P)                                   :: Kquad(3,3), R(3,3)                                              !< Matrice K^2 e matrice rotazione tra vz e normale alla spira
   !real(R8P)                                   :: dist, prj_v(3),eps, d_real                                              !< Distanza punto retta e proiezione del punto sulla retta                                                                          !< Variabile utilizzata per definire direzione corrente
   !integer(I4P)                                :: b,i,j,k,w                                                       !< Counter.
   !integer(I4P)                                :: i1,j1,k1,w1                                                     !< Counter1.
   !integer(I4P)                                :: d_int
   !!associo per dati su posizioni delle celle e contatori
   !associate(blocks_number=>field%blocks_number, ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, ngc=>field%grid%ngc, &
   !         x_c => self%x_center(n), y_c => self%y_center(n), z_c => self%z_center(n),                                         &
   !         dx => field%dxyz(1,:), dy => field%dxyz(2,:), dz => field%dxyz(3,:), lx => self%lx(n), ly => self%ly(n),           &
   !         normal => self%normal(:,n), d => self%d(n), nb =>field%nb, current_distribution => self%current_distribution(n),   &
   !         x_cell => field%x_cell, y_cell => field%y_cell, z_cell => field%z_cell, sigma => self%sigma(n))
   !c_c = [ x_c, y_c, z_c ] !Vettore posizione centro spira
   !!vertici del rettangolo, lo costruisco come se avesse normale asse z e fosse centrato nell'origine;
   !!vertici in senso antiorario, partendo da in basso a sinistra; si ipotizza
   !!rotazione rispetto al centro, quindi, ad esempio, se la vuoi con normale lungo x
   !!con lato lungo lungo y AB sarà il lato corto (pensa a rotazione 3D). Con
   !!normale lungo y il problema non si pone
   !vx = [1._R8P,0._R8P,0._R8P]
   !vy = [0._R8P,1._R8P,0._R8P]
   !vz = [0._R8P,0._R8P,1._R8P]
   !V1 = [-lx/2, -ly/2, 0._R8P]
   !V2 = [+lx/2, -ly/2, 0._R8P]
   !V3 = [+lx/2, +ly/2, 0._R8P]
   !V4 = [-lx/2, +ly/2, 0._R8P]
   !! 1: interno
   !V1_1 = [-(lx/2 - (d)/2*dx(1)), -(ly/2 - (d)/2*dx(1)), 0._R8P]
   !V2_1 = [(lx/2  - (d)/2*dx(1)), -(ly/2 - (d)/2*dx(1)), 0._R8P]
   !V3_1 = [(lx/2  - (d)/2*dx(1)),  (ly/2 - (d)/2*dx(1)), 0._R8P]
   !V4_1 = [-(lx/2 - (d)/2*dx(1)), +(ly/2 - (d)/2*dx(1)), 0._R8P]
   !! 2: esterno
   !V1_2 = [-(lx/2 + (d)/2*dx(1)), -(ly/2 + (d)/2*dx(1)), 0._R8P]
   !V2_2 = [(lx/2  + (d)/2*dx(1)), -(ly/2 + (d)/2*dx(1)), 0._R8P]
   !V3_2 = [(lx/2  + (d)/2*dx(1)),  (ly/2 + (d)/2*dx(1)), 0._R8P]
   !V4_2 = [-(lx/2 + (d)/2*dx(1)), +(ly/2 + (d)/2*dx(1)), 0._R8P]
!
   !!            V4 ________________ V3     A
   !!              |                |       A
   !!              |                |       |
   !!              |                |       |
   !!              |                |       y
   !!              |________________|
   !!            V1                  V2     x --->>>
!
   !!calcolo rotazione tra i due vettori normali %OSS una rotazione di 180° dà problemi, qua ci metti un go to.
   !!isnan per il caso rotazione nulla rispetto
   !!a n // z
   !kappa = crossproduct(a=vz,b=normal)
!
   !if (any(normal /= vz )) then
   !    kappa = kappa/sqrt(sq_norm(kappa))
   !endif
!
   !theta = acos(dotproduct(a=vz,b=normal))
!
   !K_rot(1,1) = 0._R8P
   !K_rot(1,2) = -kappa(3)
   !K_rot(1,3) = kappa(2)
   !K_rot(2,1) = kappa(3)
   !K_rot(2,2) = 0._R8P
   !K_rot(2,3) = -kappa(1)
   !K_rot(3,1) = -kappa(2)
   !K_rot(3,2) = kappa(1)
   !K_rot(3,3) = 0._R8P
   !Kquad = matmul(K_rot,K_rot) !matrice K^2
   !Id(:,1) = vx
   !Id(:,2) = vy
   !Id(:,3) = vz
!
   !R = Id+sin(theta)*K_rot+(1-cos(theta))*Kquad !costruisco matrice di rotazione, ruoto i vettori posizione e poi traslo
   !                                            ! a seconda della posizione del centro richiesta
!
   !V1 = matmul(R,V1)+c_c
   !V2 = matmul(R,V2)+c_c
   !V3 = matmul(R,V3)+c_c
   !V4 = matmul(R,V4)+c_c
   !V(1,:) = V1
   !V(2,:) = V2
   !V(3,:) = V3
   !V(4,:) = V4 !genero matrice [V1; V2; V3; V4], quindi le colonne sono le coordinate x y z e le righe i vari vertici nell'ordine descritto prima
!
   !V1_1 = matmul(R,V1_1)+c_c
   !V2_1 = matmul(R,V2_1)+c_c
   !V3_1 = matmul(R,V3_1)+c_c
   !V4_1 = matmul(R,V4_1)+c_c
   !V1_2 = matmul(R,V1_2)+c_c
   !V2_2 = matmul(R,V2_2)+c_c
   !V3_2 = matmul(R,V3_2)+c_c
   !V4_2 = matmul(R,V4_2)+c_c
   !v_l1 = matmul(R,vx)
   !v_l1 = v_l1/sqrt(sq_norm(v_l1))
   !v_l2 = matmul(R,vy);
   !v_l2 = v_l2/sqrt(sq_norm(v_l2))
!
   !!matrice dei versori dei lati, generata ruotando i versori del rettangolo tramite la matrice di rotazione precedentemente calcolata
   !vec(1,:) = v_l1
   !vec(2,:) = v_l2
   !vec(3,:) = -v_l1
   !vec(4,:) = -v_l2
!
   !!calcolo piani perpendicolari alla spira su cui giacciono le due diagonali
   !!del rettangolo. Per convenzione, le normali puntano "verso i vertici" D e C
   !!piano 1, diagonale V1-V3 e normale verso V4
   !n1 = crossproduct(a=(V1_2-V1_1),b=normal)/sqrt(sq_norm(V1_2-V1_1)) !normale al piano
   !d1 = -dotproduct(a=n1,b=V1_1) !parametro d dell'equazione ax + by + cz + d1 = 0, con a b c coseni direttori della normale
   !n2 = crossproduct(a=(V2_1-V2_2),b=normal)/sqrt(sq_norm(V2_1-V2_2)) !normale al piano
   !d2 = -dotproduct(a=n2,b=V2_1) !parametro d dell'equazione ax + by + cz + d1 = 0, con a b c coseni direttori della normale
   !n3 = crossproduct(a=(V3_1-V3_2),b=normal)/sqrt(sq_norm(V3_1-V3_2)) !normale al piano
   !d3 = -dotproduct(a=n3,b=V3_1) !parametro d dell'equazione ax + by + cz + d1 = 0, con a b c coseni direttori della normale
   !n4 = crossproduct(a=(V4_2-V4_1),b=normal)/sqrt(sq_norm(V4_2-V4_1)) !normale al piano
   !d4 = -dotproduct(a=n4,b=V4_1) !parametro d dell'equazione ax + by + cz + d1 = 0, con a b c coseni direttori della normale
!
   !allocate(flag(1:ni,1:nj,1:nk,1:nb))
   !flag(:,:,:,:) = 0_I4P !inizializzo matrice flag a zero, per indicare che nessun lato passa per le celle
   !allocate(Gaussian(1:ni,1:nj,1:nk,1:nb))
   !Gaussian(:,:,:,:) = 0_I4P !inizializzo matrice gaussiana a zero
   !d_int = int(d/2) - 1_I4P
   !d_real = real(d_int,R8P)
   !eps = 1*10e-10
   !!modificata per avere termine sorgente come Filippo, se infittiamo o aumentiamo sezione spira torna la precedente
   !do w = 1, 4, 2 !per ogni lato del rettangolo
   !   do b=1, blocks_number
   !      dmax = dx(b)+eps
   !      ! chiamo funzione che restituisce le coordinate delle varie celle che compongono il blocco
   !      do k=1+d_int, nk-d_int
   !         do j=1+d_int, nj-d_int
   !            do i=1+d_int, ni-d_int
   !               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]
   !               !print *, cell_coord, 'Coordinata della cella'
   !               dist = sqrt(sq_norm(crossproduct(a=(cell_coord-V(w,:)),b=vec(w,:)))) !Distanza punto retta |(P-A) x v| / |v| con A punto sulla retta e v versore della retta
   !               prj_v = V(w,:)+dotproduct(a=(cell_coord-V(w,:)),b=vec(w,:))*vec(w,:); !Formula proiezione di un punto su una retta con A punto della retta A+[(P-A)*v]*v
   !               if (dist <= dmax .and. &
   !                  prj_v(1) <= maxval(V(:,1)) + dmax .and. &
   !                  prj_v(2) <= maxval(V(:,2)) + dmax .and. prj_v(3) <= maxval(V(:,3)) + dmax .and. &
   !                  minval(V(:,1))-dmax <= prj_v(1) .and. minval(V(:,2))-dmax <= prj_v(2) .and. &
   !                  minval(V(:,3))-dmax <= prj_v(3) ) then
   !                  if (d_int /= 0_I4P) then
   !                     do i1 = -d_int, d_int
   !                        do j1 = -d_int, d_int
   !                           do k1 = -d_int, d_int
   !                              if (flag(i+i1,j+j1,k+k1,b) == 0_I4P) then
   !                                 flag(i+i1,j+j1,k+k1,b) = w
   !                                 !print *, w
   !                              endif
   !                           enddo
   !                        enddo
   !                     enddo
   !                  else
   !                     if (flag(i,j,k,b) == 0_I4P) then
   !                        flag(i,j,k,b) = w
   !                        !print *, w
   !                     endif
   !                  endif
   !               endif
   !            enddo
   !         enddo
   !      enddo
   !   enddo
   !   do b=1, blocks_number
   !      do k = 1, nk
   !         do j = 1, nj
   !            do i = 1, ni
   !               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]
   !               dist = sqrt(sq_norm(crossproduct(a=(cell_coord-V(w,:)),b=vec(w,:))))
   !               !secondo if: per ogni lato verifico di essere dal "lato giusto" dei piani definiti dalle diagonali, al fine di
   !               !non avere sovrapposizioni in prossimità dei vertici
   !               if (w == 1) then
   !                  if ((dotproduct(a=n1,b=cell_coord)+d1 > eps .or. dotproduct(a=n2,b=cell_coord)+d2 > eps) .and. &
   !                       flag(i,j,k,b) == w) then !aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
   !                     flag(i,j,k,b) = 0_I4P
   !                     !print *, w
   !                  endif
   !               elseif (w == 3) then
   !                  if ((dotproduct(a=n3,b=cell_coord)+d3 < -eps .or. dotproduct(a=n4,b=cell_coord)+d4 < -eps) .and. &
   !                       flag(i,j,k,b) == w) then !aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
   !                     flag(i,j,k,b) = 0_I4P
   !                     !print *, w
   !                  endif
   !               endif
   !               if (flag(i,j,k,b) == w) then
   !                  selectcase (current_distribution)
   !                  case (GAUSS_CURRENT_DISTRIBUTION)
   !                     Gaussian(i,j,k,b) = gaussian_2D_ind(sigma = sigma*dx(b), r = dist)
   !                  case (CONST_CURRENT_DISTRIBUTION)
   !                     Gaussian(i,j,k,b) = 1/(d*dx(b))**2
   !                  endselect
   !               endif
   !            enddo
   !         enddo
   !      enddo
   !   enddo
   !enddo
!
   !do w = 2, 4, 2
   !   do b=1, blocks_number
   !      dmax = dx(b)+eps
   !      ! chiamo funzione che restituisce le coordinate delle varie celle che compongono il blocco
   !      do k=1+d_int, nk-d_int
   !         do j=1+d_int, nj-d_int
   !            do i=1+d_int, ni-d_int
   !               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]
   !               dist = sqrt(sq_norm(crossproduct(a=(cell_coord-V(w,:)),b=vec(w,:)))) !Distanza punto retta |(P-A) x v| / |v| con A punto sulla retta e v versore della retta
   !               prj_v = V(w,:)+dotproduct(a=(cell_coord-V(w,:)),b=vec(w,:))*vec(w,:); !Formula proiezione di un punto su una retta con A punto della retta A+[(P-A)*v]*v
   !               if (dist <= dmax .and. &
   !               prj_v(1) <= maxval(V(:,1)) + dmax .and. &
   !               prj_v(2) <= maxval(V(:,2)) + dmax .and. prj_v(3) <= maxval(V(:,3)) + dmax .and. &
   !               minval(V(:,1))-dmax <= prj_v(1) .and. minval(V(:,2))-dmax <= prj_v(2) .and. &
   !               minval(V(:,3))-dmax <= prj_v(3) ) then
   !                  if (d_int /= 0_I4P) then
   !                     do i1 = -d_int, d_int
   !                        do j1 = -d_int, d_int
   !                           do k1 = -d_int, d_int
   !                              if (flag(i+i1,j+j1,k+k1,b) == 0_I4P) then
   !                                 flag(i+i1,j+j1,k+k1,b) = w
   !                              endif
   !                           enddo
   !                        enddo
   !                     enddo
   !                  else
   !                     if (flag(i,j,k,b) == 0_I4P) then
   !                        flag(i,j,k,b) = w
   !                     endif
   !                  endif
   !               endif
   !            enddo
   !         enddo
   !      enddo
   !   enddo
   !   do b=1, blocks_number
   !      do k = 1, nk
   !         do j = 1, nj
   !            do i = 1, ni
   !               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]
   !               dist = sqrt(sq_norm(crossproduct(a=(cell_coord-V(w,:)),b=vec(w,:))))
   !               if (w == 2) then
   !                  if ((dotproduct(a=n2,b=cell_coord)+d2 <= -eps .or. dotproduct(a=n3,b=cell_coord)+d3 >= eps) .and. &
   !                       flag(i,j,k,b) == w) then!aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
   !                     flag(i,j,k,b) = 0_I4P
   !                  endif
   !               elseif (w == 4) then
   !                  if ((dotproduct(a=n1,b=cell_coord)+d1 <= -eps .or. dotproduct(a=n4,b=cell_coord)+d4 >= eps) .and. &
   !                       flag(i,j,k,b) == w) then !aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
   !                     flag(i,j,k,b) = 0_I4P
   !                  endif
   !               endif
   !               if (flag(i,j,k,b) == w) then
   !                  selectcase (current_distribution)
   !                  case (GAUSS_CURRENT_DISTRIBUTION)
   !                     Gaussian(i,j,k,b) = gaussian_2D_ind(sigma = sigma*dx(b), r = dist)
   !                  case (CONST_CURRENT_DISTRIBUTION)
   !                     Gaussian(i,j,k,b) = 1/(d*dx(b))**2
   !                  endselect
   !               endif
   !            enddo
   !         enddo
   !      enddo
   !   enddo
   !enddo
!
   !!Ho un flag pari a 1 2 3 4 nelle celle per cui passa uno dei dati della spira. La direzione della corrente è
   !!Ceorente con quella dei versori dei lati precedentemente descritti
!
   !do b=1, blocks_number
   !   do k=1, nk
   !      do j=1, nj
   !         do i=1, ni
   !            if (flag(i,j,k,b) /= 0 .and. self%coil_flag(i,j,k,b) == 0_I4P) then
!
   !               self%J_vec(1:3,i,j,k,b) = vec(flag(i,j,k,b),:)
   !               self%J_vec(4,i,j,k,b) = Gaussian(i,j,k,b)
!
   !               if (abs(self%J_vec(1,i,j,k,b)) < 1.0e-10_R8P) then
   !                  self%J_vec(1,i,j,k,b) = 0._R8P
   !               endif
   !               if (abs(self%J_vec(2,i,j,k,b)) < 1.0e-10_R8P) then
   !                  self%J_vec(2,i,j,k,b) = 0._R8P
   !               endif
   !               if (abs(self%J_vec(3,i,j,k,b)) < 1.0e-10_R8P) then
   !                  self%J_vec(3,i,j,k,b) = 0._R8P
   !               endif
!
   !               self%coil_flag(i,j,k,b) = n
!
   !            endif
   !         enddo
   !      enddo
   !   enddo
   !enddo
   !endassociate
   !endsubroutine set_rectangular_coil_quad_section_even
!
   !subroutine set_rectangular_coil_quad_section_even_v2(self, physics, field, n)
   !class(prism_coil_object),     intent(inout) :: self                                                            !< Coils
   !type(field_object),           intent(inout) :: field                                                           !< Field object.
   !type(prism_physics_object),   intent(in)    :: physics                                                         !< Fluids physiscs.
   !integer(I4P),                 intent(in)    :: n                                                               !< Coil number.
   !integer(I4P),                 allocatable   :: flag(:,:,:,:)                                                   !< Flag per identificare se la spira passa per la cella
   !real(R8P),                    allocatable   :: Gaussian(:,:,:,:)                                               !< Matrice gaussiana per distribuzione corrente
   !real(R8P)                                   :: dmax                                                            !< Vincolo distanza massima dalla spira.
   !real(R8P)                                   :: c_c(3)                                                          !< Vettore posizione centro spira
   !real(R8P)                                   :: cell_coord(3)                                                   !< Vettore posizione centro cella
   !real(R8P)                                   :: vx(3), vy(3), vz(3)                                             !< Versori assi cartesiani
   !real(R8P)                                   :: V1(3), V2(3), V3(3), V4(3), V(4,3)                              !< Vertici rettangolo e relativa matrice
   !real(R8P)                                   :: V1_1(3), V2_1(3), V3_1(3), V4_1(3)
   !real(R8P)                                   :: V1_2(3), V2_2(3), V3_2(3), V4_2(3)
   !real(R8P)                                   :: v_l1(3), v_l2(3), vec(4,3)                                      !< Versori lati rettangolo (vale regola mano dx) e relativa matrice
   !real(R8P)                                   :: n1(3), d1, n2(3), d2, n3(3), d3, n4(3), d4                      !< Parametri piani diagonali perpendicolari a spira, per evitare sovrapposizioni
   !real(R8P)                                   :: kappa(3), K_rot(3,3), Id(3,3), theta                            !< Vettore, angolo e matrice di appoggio per formula Rodrigues + matrice identità
   !real(R8P)                                   :: Kquad(3,3), R(3,3)                                              !< Matrice K^2 e matrice rotazione tra vz e normale alla spira
   !real(R8P)                                   :: dist, prj_v(3),eps, d_real                                              !< Distanza punto retta e proiezione del punto sulla retta                                                                          !< Variabile utilizzata per definire direzione corrente
   !integer(I4P)                                :: b,i,j,k,w                                                       !< Counter.
   !integer(I4P)                                :: i1,j1,k1,w1                                                     !< Counter1.
   !integer(I4P)                                :: d_int
   !integer(I4P)                                :: ind_vertex(4,2)
   !!associo per dati su posizioni delle celle e contatori
   !associate(blocks_number=>field%blocks_number, ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, ngc=>field%grid%ngc, &
   !         x_c => self%x_center(n), y_c => self%y_center(n), z_c => self%z_center(n),                                         &
   !         dx => field%dxyz(1,:), dy => field%dxyz(2,:), dz => field%dxyz(3,:), lx => self%lx(n), ly => self%ly(n),           &
   !         normal => self%normal(:,n), d => self%d(n), nb =>field%nb, current_distribution => self%current_distribution(n),   &
   !         x_cell => field%x_cell, y_cell => field%y_cell, z_cell => field%z_cell, sigma => self%sigma(n))
   !c_c = [ x_c, y_c, z_c ] !Vettore posizione centro spira
   !!vertici del rettangolo, lo costruisco come se avesse normale asse z e fosse centrato nell'origine;
   !!vertici in senso antiorario, partendo da in basso a sinistra; si ipotizza
   !!rotazione rispetto al centro, quindi, ad esempio, se la vuoi con normale lungo x
   !!con lato lungo lungo y AB sarà il lato corto (pensa a rotazione 3D). Con
   !!normale lungo y il problema non si pone
   !vx = [1._R8P,0._R8P,0._R8P]
   !vy = [0._R8P,1._R8P,0._R8P]
   !vz = [0._R8P,0._R8P,1._R8P]
   !V1 = [-lx/2, -ly/2, 0._R8P]
   !V2 = [+lx/2, -ly/2, 0._R8P]
   !V3 = [+lx/2, +ly/2, 0._R8P]
   !V4 = [-lx/2, +ly/2, 0._R8P]
!
   !ind_vertex(:,1) = [1,2,3,4]
   !ind_vertex(:,2) = [2,3,4,1]
!
   !! 1: interno
   !V1_1 = [-(lx/2 - (d)/2*dx(1)), -(ly/2 - (d)/2*dx(1)), 0._R8P]
   !V2_1 = [(lx/2  - (d)/2*dx(1)), -(ly/2 - (d)/2*dx(1)), 0._R8P]
   !V3_1 = [(lx/2  - (d)/2*dx(1)),  (ly/2 - (d)/2*dx(1)), 0._R8P]
   !V4_1 = [-(lx/2 - (d)/2*dx(1)), +(ly/2 - (d)/2*dx(1)), 0._R8P]
   !! 2: esterno
   !V1_2 = [-(lx/2 + (d)/2*dx(1)), -(ly/2 + (d)/2*dx(1)), 0._R8P]
   !V2_2 = [(lx/2  + (d)/2*dx(1)), -(ly/2 + (d)/2*dx(1)), 0._R8P]
   !V3_2 = [(lx/2  + (d)/2*dx(1)),  (ly/2 + (d)/2*dx(1)), 0._R8P]
   !V4_2 = [-(lx/2 + (d)/2*dx(1)), +(ly/2 + (d)/2*dx(1)), 0._R8P]
!
   !!            V4 ________________ V3     A
   !!              |                |       A
   !!              |                |       |
   !!              |                |       |
   !!              |                |       y
   !!              |________________|
   !!            V1                  V2     x --->>>
!
   !!calcolo rotazione tra i due vettori normali %OSS una rotazione di 180° dà problemi, qua ci metti un go to.
   !!isnan per il caso rotazione nulla rispetto
   !!a n // z
   !kappa = crossproduct(a=vz,b=normal)
!
   !if (any(normal /= vz )) then
   !    kappa = kappa/sqrt(sq_norm(kappa))
   !endif
!
   !theta = acos(dotproduct(a=vz,b=normal))
!
   !K_rot(1,1) = 0._R8P
   !K_rot(1,2) = -kappa(3)
   !K_rot(1,3) = kappa(2)
   !K_rot(2,1) = kappa(3)
   !K_rot(2,2) = 0._R8P
   !K_rot(2,3) = -kappa(1)
   !K_rot(3,1) = -kappa(2)
   !K_rot(3,2) = kappa(1)
   !K_rot(3,3) = 0._R8P
   !Kquad = matmul(K_rot,K_rot) !matrice K^2
   !Id(:,1) = vx
   !Id(:,2) = vy
   !Id(:,3) = vz
!
   !R = Id+sin(theta)*K_rot+(1-cos(theta))*Kquad !costruisco matrice di rotazione, ruoto i vettori posizione e poi traslo
   !                                            ! a seconda della posizione del centro richiesta
!
   !V1 = matmul(R,V1)+c_c
   !V2 = matmul(R,V2)+c_c
   !V3 = matmul(R,V3)+c_c
   !V4 = matmul(R,V4)+c_c
   !V(1,:) = V1
   !V(2,:) = V2
   !V(3,:) = V3
   !V(4,:) = V4 !genero matrice [V1; V2; V3; V4], quindi le colonne sono le coordinate x y z e le righe i vari vertici nell'ordine descritto prima
!
   !V1_1 = matmul(R,V1_1)+c_c
   !V2_1 = matmul(R,V2_1)+c_c
   !V3_1 = matmul(R,V3_1)+c_c
   !V4_1 = matmul(R,V4_1)+c_c
   !V1_2 = matmul(R,V1_2)+c_c
   !V2_2 = matmul(R,V2_2)+c_c
   !V3_2 = matmul(R,V3_2)+c_c
   !V4_2 = matmul(R,V4_2)+c_c
   !v_l1 = matmul(R,vx)
   !v_l1 = v_l1/sqrt(sq_norm(v_l1))
   !v_l2 = matmul(R,vy);
   !v_l2 = v_l2/sqrt(sq_norm(v_l2))
!
   !!matrice dei versori dei lati, generata ruotando i versori del rettangolo tramite la matrice di rotazione precedentemente calcolata
   !vec(1,:) = v_l1
   !vec(2,:) = v_l2
   !vec(3,:) = -v_l1
   !vec(4,:) = -v_l2
!
   !do i = 1, 4
   !   do j = 1, 3
   !      if (abs(vec(i,j)) < 1.0e-10_R8P) then
   !         vec(i,j) = 0._R8P
   !      endif
   !   enddo
   !enddo
!
   !!calcolo piani perpendicolari alla spira su cui giacciono le due diagonali
   !!del rettangolo. Per convenzione, le normali puntano "verso i vertici" D e C
   !!piano 1, diagonale V1-V3 e normale verso V4
   !n1 = crossproduct(a=(V1_2-V1_1),b=normal)/sqrt(sq_norm(V1_2-V1_1)) !normale al piano
   !d1 = -dotproduct(a=n1,b=V1_1) !parametro d dell'equazione ax + by + cz + d1 = 0, con a b c coseni direttori della normale
   !n2 = crossproduct(a=(V2_1-V2_2),b=normal)/sqrt(sq_norm(V2_1-V2_2)) !normale al piano
   !d2 = -dotproduct(a=n2,b=V2_1) !parametro d dell'equazione ax + by + cz + d1 = 0, con a b c coseni direttori della normale
   !n3 = crossproduct(a=(V3_1-V3_2),b=normal)/sqrt(sq_norm(V3_1-V3_2)) !normale al piano
   !d3 = -dotproduct(a=n3,b=V3_1) !parametro d dell'equazione ax + by + cz + d1 = 0, con a b c coseni direttori della normale
   !n4 = crossproduct(a=(V4_2-V4_1),b=normal)/sqrt(sq_norm(V4_2-V4_1)) !normale al piano
   !d4 = -dotproduct(a=n4,b=V4_1) !parametro d dell'equazione ax + by + cz + d1 = 0, con a b c coseni direttori della normale
!
   !allocate(flag(1:ni,1:nj,1:nk,1:nb))
   !flag(:,:,:,:) = 0_I4P !inizializzo matrice flag a zero, per indicare che nessun lato passa per le celle
   !allocate(Gaussian(1:ni,1:nj,1:nk,1:nb))
   !Gaussian(:,:,:,:) = 0_I4P !inizializzo matrice gaussiana a zero
   !d_int = int(d/2) - 1_I4P
   !d_real = real(d_int,R8P)
   !eps = 1*10e-10
   !!modificata per avere termine sorgente come Filippo, se infittiamo o aumentiamo sezione spira torna la precedente
   !do w = 1, 4!, 2 !per ogni lato del rettangolo
   !   do b=1, blocks_number
   !      dmax = dx(b)+eps
   !      ! chiamo funzione che restituisce le coordinate delle varie celle che compongono il blocco
   !      do k=1+d_int, nk-d_int
   !         do j=1+d_int, nj-d_int
   !            do i=1+d_int, ni-d_int
   !               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]
   !               !print *, cell_coord, 'Coordinata della cella'
   !               dist = sqrt(sq_norm(crossproduct(a=(cell_coord-V(w,:)),b=vec(w,:)))) !Distanza punto retta |(P-A) x v| / |v| con A punto sulla retta e v versore della retta
   !               prj_v = V(w,:)+dotproduct(a=(cell_coord-V(w,:)),b=vec(w,:))*vec(w,:) !Formula proiezione di un punto su una retta con A punto della retta A+[(P-A)*v]*v
   !               if (dist <= dmax .and. &
   !                  prj_v(1) <= maxval(V(:,1)) + dmax .and. &
   !                  prj_v(2) <= maxval(V(:,2)) + dmax .and. prj_v(3) <= maxval(V(:,3)) + dmax .and. &
   !                  minval(V(:,1))-dmax <= prj_v(1) .and. minval(V(:,2))-dmax <= prj_v(2) .and. &
   !                  minval(V(:,3))-dmax <= prj_v(3) ) then
   !                  if (d_int /= 0_I4P) then
   !                     do i1 = -d_int, d_int
   !                        do j1 = -d_int, d_int
   !                           do k1 = -d_int, d_int
   !                              !if (flag(i+i1,j+j1,k+k1,b) == 0_I4P) then
   !                                 flag(i+i1,j+j1,k+k1,b) = w
   !                                 !print *, w
   !                              !endif
   !                           enddo
   !                        enddo
   !                     enddo
   !                  else
   !                     !if (flag(i,j,k,b) == 0_I4P) then
   !                        flag(i,j,k,b) = w
   !                        !print *, w
   !                     !endif
   !                  endif
   !               endif
   !            enddo
   !         enddo
   !      enddo
   !   enddo
   !      !Ho un flag pari a 1 2 3 4 nelle celle per cui passa uno dei dati della spira. La direzione della corrente è
   !      !Coerente con quella dei versori dei lati precedentemente descritti
   !   do b=1, blocks_number
   !      do k = 1, nk
   !         do j = 1, nj
   !            do i = 1, ni
   !               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]
   !               prj_v = V(w,:)+dotproduct(a=(cell_coord-V(w,:)),b=vec(w,:))*vec(w,:) !Formula proiezione di un punto su una retta con A punto della retta A+[(P-A)*v]*v
!
   !               if (prj_v(1) <= maxval(V(:,1)) .and. prj_v(2) <= maxval(V(:,2)) .and. &
   !                   prj_v(3) <= maxval(V(:,3)) .and. minval(V(:,1)) <= prj_v(1) .and. &
   !                   minval(V(:,2)) <= prj_v(2) .and. minval(V(:,3)) <= prj_v(3)) then
!
   !                  dist = sqrt(sq_norm(crossproduct(a=(cell_coord-V(w,:)),b=vec(w,:))))
   !               else
   !                  dist = minval([sqrt(sq_norm(cell_coord-V(ind_vertex(w,1),:))), sqrt(sq_norm(cell_coord-V(ind_vertex(w,2),:)))])
   !               endif
   !               !secondo if: per ogni lato verifico di essere dal "lato giusto" dei piani definiti dalle diagonali, al fine di
   !               !non avere sovrapposizioni in prossimità dei vertici
   !               !if (w == 1) then
   !               !   if ((dotproduct(a=n1,b=cell_coord)+d1 > eps .or. dotproduct(a=n2,b=cell_coord)+d2 > eps) .and. &
   !               !        flag(i,j,k,b) == w) then !aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
   !               !      flag(i,j,k,b) = 0_I4P
   !               !      !print *, w
   !               !   endif
   !               !elseif (w == 3) then
   !               !   if ((dotproduct(a=n3,b=cell_coord)+d3 < -eps .or. dotproduct(a=n4,b=cell_coord)+d4 < -eps) .and. &
   !               !        flag(i,j,k,b) == w) then !aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
   !               !      flag(i,j,k,b) = 0_I4P
   !               !      !print *, w
   !               !   endif
   !               !endif
   !               if (flag(i,j,k,b) == w) then
   !                  selectcase (current_distribution)
   !                  case (GAUSS_CURRENT_DISTRIBUTION)
   !                     Gaussian(i,j,k,b) = gaussian_2D_ind(sigma = sigma*dx(b), r = dist)
   !                  case (CONST_CURRENT_DISTRIBUTION)
   !                     Gaussian(i,j,k,b) = 1/(d*dx(b))**2
   !                  endselect
   !                  self%J_vec(1:3,i,j,k,b) = self%J_vec(1:3,i,j,k,b)+vec(flag(i,j,k,b),:)*Gaussian(i,j,k,b)
   !                  self%J_vec(4,i,j,k,b) = 1._R8P
   !                  self%coil_flag(i,j,k,b) = n
   !               endif
   !            enddo
   !         enddo
   !      enddo
   !   enddo
   !   !do b=1, blocks_number
   !   !   do k=1, nk
   !   !      do j=1, nj
   !   !         do i=1, ni
   !   !            if (flag(i,j,k,b) == w) then! .and. self%coil_flag(i,j,k,b) == 0_I4P) then
   !   !               self%J_vec(1:3,i,j,k,b) = self%J_vec(1:3,i,j,k,b)+vec(flag(i,j,k,b),:)*Gaussian(i,j,k,b) !Se una cella è attraversata da più lati, sommo i contributi alla corrente
   !   !               self%J_vec(4,i,j,k,b) = 1._R8P
   !   !               !if (abs(self%J_vec(1,i,j,k,b)) < 1.0e-10_R8P) then
   !   !               !   self%J_vec(1,i,j,k,b) = 0._R8P
   !   !               !endif
   !   !               !if (abs(self%J_vec(2,i,j,k,b)) < 1.0e-10_R8P) then
   !   !               !   self%J_vec(2,i,j,k,b) = 0._R8P
   !   !               !endif
   !   !               !if (abs(self%J_vec(3,i,j,k,b)) < 1.0e-10_R8P) then
   !   !               !   self%J_vec(3,i,j,k,b) = 0._R8P
   !   !               !endif
   !   !               self%coil_flag(i,j,k,b) = n
   !   !            endif
   !   !         enddo
   !   !      enddo
   !   !   enddo
   !   !enddo
   !enddo
!
!
   !!do w = 2, 4, 2
   !!   do b=1, blocks_number
   !!      dmax = dx(b)+eps
   !!      ! chiamo funzione che restituisce le coordinate delle varie celle che compongono il blocco
   !!      do k=1+d_int, nk-d_int
   !!         do j=1+d_int, nj-d_int
   !!            do i=1+d_int, ni-d_int
   !!               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]
   !!               dist = sqrt(sq_norm(crossproduct(a=(cell_coord-V(w,:)),b=vec(w,:)))) !Distanza punto retta |(P-A) x v| / |v| con A punto sulla retta e v versore della retta
   !!               prj_v = V(w,:)+dotproduct(a=(cell_coord-V(w,:)),b=vec(w,:))*vec(w,:); !Formula proiezione di un punto su una retta con A punto della retta A+[(P-A)*v]*v
   !!               if (dist <= dmax .and. &
   !!               prj_v(1) <= maxval(V(:,1)) + dmax .and. &
   !!               prj_v(2) <= maxval(V(:,2)) + dmax .and. prj_v(3) <= maxval(V(:,3)) + dmax .and. &
   !!               minval(V(:,1))-dmax <= prj_v(1) .and. minval(V(:,2))-dmax <= prj_v(2) .and. &
   !!               minval(V(:,3))-dmax <= prj_v(3) ) then
   !!                  if (d_int /= 0_I4P) then
   !!                     do i1 = -d_int, d_int
   !!                        do j1 = -d_int, d_int
   !!                           do k1 = -d_int, d_int
   !!                              if (flag(i+i1,j+j1,k+k1,b) == 0_I4P) then
   !!                                 flag(i+i1,j+j1,k+k1,b) = w
   !!                              endif
   !!                           enddo
   !!                        enddo
   !!                     enddo
   !!                  else
   !!                     if (flag(i,j,k,b) == 0_I4P) then
   !!                        flag(i,j,k,b) = w
   !!                     endif
   !!                  endif
   !!               endif
   !!            enddo
   !!         enddo
   !!      enddo
   !!   enddo
   !!   do b=1, blocks_number
   !!      do k = 1, nk
   !!         do j = 1, nj
   !!            do i = 1, ni
   !!               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]
   !!               dist = sqrt(sq_norm(crossproduct(a=(cell_coord-V(w,:)),b=vec(w,:))))
   !!               if (w == 2) then
   !!                  if ((dotproduct(a=n2,b=cell_coord)+d2 <= -eps .or. dotproduct(a=n3,b=cell_coord)+d3 >= eps) .and. &
   !!                       flag(i,j,k,b) == w) then!aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
   !!                     flag(i,j,k,b) = 0_I4P
   !!                  endif
   !!               elseif (w == 4) then
   !!                  if ((dotproduct(a=n1,b=cell_coord)+d1 <= -eps .or. dotproduct(a=n4,b=cell_coord)+d4 >= eps) .and. &
   !!                       flag(i,j,k,b) == w) then !aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
   !!                     flag(i,j,k,b) = 0_I4P
   !!                  endif
   !!               endif
   !!               if (flag(i,j,k,b) == w) then
   !!                  selectcase (current_distribution)
   !!                  case (GAUSS_CURRENT_DISTRIBUTION)
   !!                     Gaussian(i,j,k,b) = gaussian_2D_ind(sigma = sigma*dx(b), r = dist)
   !!                  case (CONST_CURRENT_DISTRIBUTION)
   !!                     Gaussian(i,j,k,b) = 1/(d*dx(b))**2
   !!                  endselect
   !!               endif
   !!            enddo
   !!         enddo
   !!      enddo
   !!   enddo
   !!enddo
   !endassociate
   !endsubroutine set_rectangular_coil_quad_section_even_v2
!
   !subroutine set_rectangular_coil_quad_section_odd_v2(self, physics, field, n)
   !class(prism_coil_object),     intent(inout) :: self                                                            !< Coils
   !type(field_object),           intent(inout) :: field                                                           !< Field object.
   !type(prism_physics_object),   intent(in)    :: physics                                                         !< Fluids physiscs.
   !integer(I4P),                 intent(in)    :: n                                                               !< Coil number.
   !integer(I4P),                 allocatable   :: flag(:,:,:,:)                                                   !< Flag per identificare se la spira passa per la cella
   !real(R8P),                    allocatable   :: Gaussian(:,:,:,:)                                               !< Matrice gaussiana per distribuzione corrente
   !real(R8P)                                   :: dmax                                                            !< Vincolo distanza massima dalla spira.
   !real(R8P)                                   :: c_c(3)                                                          !< Vettore posizione centro spira
   !real(R8P)                                   :: cell_coord(3)                                                   !< Vettore posizione centro cella
   !real(R8P)                                   :: vx(3), vy(3), vz(3)                                             !< Versori assi cartesiani
   !real(R8P)                                   :: V1(3), V2(3), V3(3), V4(3), V(4,3)                              !< Vertici rettangolo e relativa matrice
   !real(R8P)                                   :: V1_1(3), V2_1(3), V3_1(3), V4_1(3)
   !real(R8P)                                   :: V1_2(3), V2_2(3), V3_2(3), V4_2(3)
   !real(R8P)                                   :: v_l1(3), v_l2(3), vec(4,3)                                      !< Versori lati rettangolo (vale regola mano dx) e relativa matrice
   !real(R8P)                                   :: n1(3), d1, n2(3), d2, n3(3), d3, n4(3), d4                      !< Parametri piani diagonali perpendicolari a spira, per evitare sovrapposizioni
   !real(R8P)                                   :: kappa(3), K_rot(3,3), Id(3,3), theta                            !< Vettore, angolo e matrice di appoggio per formula Rodrigues + matrice identità
   !real(R8P)                                   :: Kquad(3,3), R(3,3)                                              !< Matrice K^2 e matrice rotazione tra vz e normale alla spira
   !real(R8P)                                   :: dist, prj_v(3),eps, d_real                                              !< Distanza punto retta e proiezione del punto sulla retta                                                                          !< Variabile utilizzata per definire direzione corrente
   !integer(I4P)                                :: b,i,j,k,w                                                       !< Counter.
   !integer(I4P)                                :: i1,j1,k1,w1                                                     !< Counter1.
   !integer(I4P)                                :: d_int
   !integer(I4P)                                :: ind_vertex(4,2)
   !!associo per dati su posizioni delle celle e contatori
   !associate(blocks_number=>field%blocks_number, ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, ngc=>field%grid%ngc, &
   !         x_c => self%x_center(n), y_c => self%y_center(n), z_c => self%z_center(n),                                         &
   !         dx => field%dxyz(1,:), dy => field%dxyz(2,:), dz => field%dxyz(3,:), lx => self%lx(n), ly => self%ly(n),           &
   !         normal => self%normal(:,n), d => self%d(n), nb =>field%nb, current_distribution => self%current_distribution(n),   &
   !         x_cell => field%x_cell, y_cell => field%y_cell, z_cell => field%z_cell, sigma => self%sigma(n))
   !c_c = [ x_c, y_c, z_c ] !Vettore posizione centro spira
   !!vertici del rettangolo, lo costruisco come se avesse normale asse z e fosse centrato nell'origine;
   !!vertici in senso antiorario, partendo da in basso a sinistra; si ipotizza
   !!rotazione rispetto al centro, quindi, ad esempio, se la vuoi con normale lungo x
   !!con lato lungo lungo y AB sarà il lato corto (pensa a rotazione 3D). Con
   !!normale lungo y il problema non si pone
   !vx = [1._R8P,0._R8P,0._R8P]
   !vy = [0._R8P,1._R8P,0._R8P]
   !vz = [0._R8P,0._R8P,1._R8P]
   !V1 = [-lx/2, -ly/2, 0._R8P]
   !V2 = [+lx/2, -ly/2, 0._R8P]
   !V3 = [+lx/2, +ly/2, 0._R8P]
   !V4 = [-lx/2, +ly/2, 0._R8P]
!
   !ind_vertex(:,1) = [1,2,3,4]
   !ind_vertex(:,2) = [2,3,4,1]
!
   !! 1: interno
   !V1_1 = [-(lx/2 - (d)/2*dx(1)), -(ly/2 - (d)/2*dx(1)), 0._R8P]
   !V2_1 = [(lx/2  - (d)/2*dx(1)), -(ly/2 - (d)/2*dx(1)), 0._R8P]
   !V3_1 = [(lx/2  - (d)/2*dx(1)),  (ly/2 - (d)/2*dx(1)), 0._R8P]
   !V4_1 = [-(lx/2 - (d)/2*dx(1)), +(ly/2 - (d)/2*dx(1)), 0._R8P]
   !! 2: esterno
   !V1_2 = [-(lx/2 + (d)/2*dx(1)), -(ly/2 + (d)/2*dx(1)), 0._R8P]
   !V2_2 = [(lx/2  + (d)/2*dx(1)), -(ly/2 + (d)/2*dx(1)), 0._R8P]
   !V3_2 = [(lx/2  + (d)/2*dx(1)),  (ly/2 + (d)/2*dx(1)), 0._R8P]
   !V4_2 = [-(lx/2 + (d)/2*dx(1)), +(ly/2 + (d)/2*dx(1)), 0._R8P]
!
   !!            V4 ________________ V3     A
   !!              |                |       A
   !!              |                |       |
   !!              |                |       |
   !!              |                |       y
   !!              |________________|
   !!            V1                  V2     x --->>>
!
   !!calcolo rotazione tra i due vettori normali %OSS una rotazione di 180° dà problemi, qua ci metti un go to.
   !!isnan per il caso rotazione nulla rispetto
   !!a n // z
   !kappa = crossproduct(a=vz,b=normal)
!
   !if (any(normal /= vz )) then
   !    kappa = kappa/sqrt(sq_norm(kappa))
   !endif
!
   !theta = acos(dotproduct(a=vz,b=normal))
!
   !K_rot(1,1) = 0._R8P
   !K_rot(1,2) = -kappa(3)
   !K_rot(1,3) = kappa(2)
   !K_rot(2,1) = kappa(3)
   !K_rot(2,2) = 0._R8P
   !K_rot(2,3) = -kappa(1)
   !K_rot(3,1) = -kappa(2)
   !K_rot(3,2) = kappa(1)
   !K_rot(3,3) = 0._R8P
   !Kquad = matmul(K_rot,K_rot) !matrice K^2
   !Id(:,1) = vx
   !Id(:,2) = vy
   !Id(:,3) = vz
!
   !R = Id+sin(theta)*K_rot+(1-cos(theta))*Kquad !costruisco matrice di rotazione, ruoto i vettori posizione e poi traslo
   !                                            ! a seconda della posizione del centro richiesta
!
   !V1 = matmul(R,V1)+c_c
   !V2 = matmul(R,V2)+c_c
   !V3 = matmul(R,V3)+c_c
   !V4 = matmul(R,V4)+c_c
   !V(1,:) = V1
   !V(2,:) = V2
   !V(3,:) = V3
   !V(4,:) = V4 !genero matrice [V1; V2; V3; V4], quindi le colonne sono le coordinate x y z e le righe i vari vertici nell'ordine descritto prima
!
   !V1_1 = matmul(R,V1_1)+c_c
   !V2_1 = matmul(R,V2_1)+c_c
   !V3_1 = matmul(R,V3_1)+c_c
   !V4_1 = matmul(R,V4_1)+c_c
   !V1_2 = matmul(R,V1_2)+c_c
   !V2_2 = matmul(R,V2_2)+c_c
   !V3_2 = matmul(R,V3_2)+c_c
   !V4_2 = matmul(R,V4_2)+c_c
   !v_l1 = matmul(R,vx)
   !v_l1 = v_l1/sqrt(sq_norm(v_l1))
   !v_l2 = matmul(R,vy);
   !v_l2 = v_l2/sqrt(sq_norm(v_l2))
!
   !!matrice dei versori dei lati, generata ruotando i versori del rettangolo tramite la matrice di rotazione precedentemente calcolata
   !vec(1,:) = v_l1
   !vec(2,:) = v_l2
   !vec(3,:) = -v_l1
   !vec(4,:) = -v_l2
!
   !do i = 1, 4
   !   do j = 1, 3
   !      if (abs(vec(i,j)) < 1.0e-10_R8P) then
   !         vec(i,j) = 0._R8P
   !      endif
   !   enddo
   !enddo
!
   !!calcolo piani perpendicolari alla spira su cui giacciono le due diagonali
   !!del rettangolo. Per convenzione, le normali puntano "verso i vertici" D e C
   !!piano 1, diagonale V1-V3 e normale verso V4
   !n1 = crossproduct(a=(V1_2-V1_1),b=normal)/sqrt(sq_norm(V1_2-V1_1)) !normale al piano
   !d1 = -dotproduct(a=n1,b=V1_1) !parametro d dell'equazione ax + by + cz + d1 = 0, con a b c coseni direttori della normale
   !n2 = crossproduct(a=(V2_1-V2_2),b=normal)/sqrt(sq_norm(V2_1-V2_2)) !normale al piano
   !d2 = -dotproduct(a=n2,b=V2_1) !parametro d dell'equazione ax + by + cz + d1 = 0, con a b c coseni direttori della normale
   !n3 = crossproduct(a=(V3_1-V3_2),b=normal)/sqrt(sq_norm(V3_1-V3_2)) !normale al piano
   !d3 = -dotproduct(a=n3,b=V3_1) !parametro d dell'equazione ax + by + cz + d1 = 0, con a b c coseni direttori della normale
   !n4 = crossproduct(a=(V4_2-V4_1),b=normal)/sqrt(sq_norm(V4_2-V4_1)) !normale al piano
   !d4 = -dotproduct(a=n4,b=V4_1) !parametro d dell'equazione ax + by + cz + d1 = 0, con a b c coseni direttori della normale
!
   !allocate(flag(1:ni,1:nj,1:nk,1:nb))
   !flag(:,:,:,:) = 0_I4P !inizializzo matrice flag a zero, per indicare che nessun lato passa per le celle
   !allocate(Gaussian(1:ni,1:nj,1:nk,1:nb))
   !Gaussian(:,:,:,:) = 0_I4P !inizializzo matrice gaussiana a zero
   !d_int = int(d/2)
   !d_real = real(d_int,R8P)
   !eps = 1*10e-10
   !!modificata per avere termine sorgente come Filippo, se infittiamo o aumentiamo sezione spira torna la precedente
   !do w = 1, 4!, 2 !per ogni lato del rettangolo
   !   do b=1, blocks_number
   !      dmax = dx(b)*d_real+eps
   !      ! chiamo funzione che restituisce le coordinate delle varie celle che compongono il blocco
   !      do k=1+d_int, nk-d_int
   !         do j=1+d_int, nj-d_int
   !            do i=1+d_int, ni-d_int
   !               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]
   !               !print *, cell_coord, 'Coordinata della cella'
   !               dist = sqrt(sq_norm(crossproduct(a=(cell_coord-V(w,:)),b=vec(w,:)))) !Distanza punto retta |(P-A) x v| / |v| con A punto sulla retta e v versore della retta
   !               prj_v = V(w,:)+dotproduct(a=(cell_coord-V(w,:)),b=vec(w,:))*vec(w,:) !Formula proiezione di un punto su una retta con A punto della retta A+[(P-A)*v]*v
   !               if (dist <= dmax .and. &
   !                  prj_v(1) <= maxval(V(:,1)) + dmax .and. &
   !                  prj_v(2) <= maxval(V(:,2)) + dmax .and. prj_v(3) <= maxval(V(:,3)) + dmax .and. &
   !                  minval(V(:,1))-dmax <= prj_v(1) .and. minval(V(:,2))-dmax <= prj_v(2) .and. &
   !                  minval(V(:,3))-dmax <= prj_v(3) ) then
   !                  if (d_int /= 0_I4P) then
   !                     do i1 = -d_int, d_int
   !                        do j1 = -d_int, d_int
   !                           do k1 = -d_int, d_int
   !                              !if (flag(i+i1,j+j1,k+k1,b) == 0_I4P) then
   !                                 flag(i+i1,j+j1,k+k1,b) = w
   !                                 !print *, w
   !                              !endif
   !                           enddo
   !                        enddo
   !                     enddo
   !                  else
   !                     !if (flag(i,j,k,b) == 0_I4P) then
   !                        flag(i,j,k,b) = w
   !                        !print *, w
   !                     !endif
   !                  endif
   !               endif
   !            enddo
   !         enddo
   !      enddo
   !   enddo
   !      !Ho un flag pari a 1 2 3 4 nelle celle per cui passa uno dei dati della spira. La direzione della corrente è
   !      !Coerente con quella dei versori dei lati precedentemente descritti
   !   do b=1, blocks_number
   !      do k = 1, nk
   !         do j = 1, nj
   !            do i = 1, ni
   !               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]
   !               prj_v = V(w,:)+dotproduct(a=(cell_coord-V(w,:)),b=vec(w,:))*vec(w,:) !Formula proiezione di un punto su una retta con A punto della retta A+[(P-A)*v]*v
!
   !               if (prj_v(1) <= maxval(V(:,1))+dx(b) .and. prj_v(2) <= maxval(V(:,2))+dx(b) .and. &
   !                   prj_v(3) <= maxval(V(:,3))+dx(b) .and. minval(V(:,1))-dx(b) <= prj_v(1) .and. &
   !                   minval(V(:,2))-dx(b) <= prj_v(2) .and. minval(V(:,3))-dx(b) <= prj_v(3)) then
!
   !                 dist = sqrt(sq_norm(crossproduct(a=(cell_coord-V(w,:)),b=vec(w,:))))
   !               else
   !                  dist = minval([sqrt(sq_norm(cell_coord-V(ind_vertex(w,1),:))), sqrt(sq_norm(cell_coord-V(ind_vertex(w,2),:)))])
   !               endif
   !               !secondo if: per ogni lato verifico di essere dal "lato giusto" dei piani definiti dalle diagonali, al fine di
   !               !non avere sovrapposizioni in prossimità dei vertici
   !               !if (w == 1) then
   !               !   if ((dotproduct(a=n1,b=cell_coord)+d1 > eps .or. dotproduct(a=n2,b=cell_coord)+d2 > eps) .and. &
   !               !        flag(i,j,k,b) == w) then !aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
   !               !      flag(i,j,k,b) = 0_I4P
   !               !      !print *, w
   !               !   endif
   !               !elseif (w == 3) then
   !               !   if ((dotproduct(a=n3,b=cell_coord)+d3 < -eps .or. dotproduct(a=n4,b=cell_coord)+d4 < -eps) .and. &
   !               !        flag(i,j,k,b) == w) then !aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
   !               !      flag(i,j,k,b) = 0_I4P
   !               !      !print *, w
   !               !   endif
   !               !endif
   !               if (flag(i,j,k,b) == w) then
   !                  selectcase (current_distribution)
   !                  case (GAUSS_CURRENT_DISTRIBUTION)
   !                     Gaussian(i,j,k,b) = gaussian_2D_ind(sigma = sigma*dx(b), r = dist)
   !                  case (CONST_CURRENT_DISTRIBUTION)
   !                     Gaussian(i,j,k,b) = 1/(d*dx(b))**2
   !                  endselect
   !                  self%J_vec(1:3,i,j,k,b) = self%J_vec(1:3,i,j,k,b)+vec(flag(i,j,k,b),:)*Gaussian(i,j,k,b)
   !                  self%J_vec(4,i,j,k,b) = 1._R8P
   !                  self%coil_flag(i,j,k,b) = n
   !               endif
   !            enddo
   !         enddo
   !      enddo
   !   enddo
   !   !do b=1, blocks_number
   !   !   do k=1, nk
   !   !      do j=1, nj
   !   !         do i=1, ni
   !   !            if (flag(i,j,k,b) == w) then! .and. self%coil_flag(i,j,k,b) == 0_I4P) then
   !   !               self%J_vec(1:3,i,j,k,b) = self%J_vec(1:3,i,j,k,b)+vec(flag(i,j,k,b),:)*Gaussian(i,j,k,b) !Se una cella è attraversata da più lati, sommo i contributi alla corrente
   !   !               self%J_vec(4,i,j,k,b) = 1._R8P
   !   !               !if (abs(self%J_vec(1,i,j,k,b)) < 1.0e-10_R8P) then
   !   !               !   self%J_vec(1,i,j,k,b) = 0._R8P
   !   !               !endif
   !   !               !if (abs(self%J_vec(2,i,j,k,b)) < 1.0e-10_R8P) then
   !   !               !   self%J_vec(2,i,j,k,b) = 0._R8P
   !   !               !endif
   !   !               !if (abs(self%J_vec(3,i,j,k,b)) < 1.0e-10_R8P) then
   !   !               !   self%J_vec(3,i,j,k,b) = 0._R8P
   !   !               !endif
   !   !               self%coil_flag(i,j,k,b) = n
   !   !            endif
   !   !         enddo
   !   !      enddo
   !   !   enddo
   !   !enddo
   !enddo
!
!
   !!do w = 2, 4, 2
   !!   do b=1, blocks_number
   !!      dmax = dx(b)+eps
   !!      ! chiamo funzione che restituisce le coordinate delle varie celle che compongono il blocco
   !!      do k=1+d_int, nk-d_int
   !!         do j=1+d_int, nj-d_int
   !!            do i=1+d_int, ni-d_int
   !!               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]
   !!               dist = sqrt(sq_norm(crossproduct(a=(cell_coord-V(w,:)),b=vec(w,:)))) !Distanza punto retta |(P-A) x v| / |v| con A punto sulla retta e v versore della retta
   !!               prj_v = V(w,:)+dotproduct(a=(cell_coord-V(w,:)),b=vec(w,:))*vec(w,:); !Formula proiezione di un punto su una retta con A punto della retta A+[(P-A)*v]*v
   !!               if (dist <= dmax .and. &
   !!               prj_v(1) <= maxval(V(:,1)) + dmax .and. &
   !!               prj_v(2) <= maxval(V(:,2)) + dmax .and. prj_v(3) <= maxval(V(:,3)) + dmax .and. &
   !!               minval(V(:,1))-dmax <= prj_v(1) .and. minval(V(:,2))-dmax <= prj_v(2) .and. &
   !!               minval(V(:,3))-dmax <= prj_v(3) ) then
   !!                  if (d_int /= 0_I4P) then
   !!                     do i1 = -d_int, d_int
   !!                        do j1 = -d_int, d_int
   !!                           do k1 = -d_int, d_int
   !!                              if (flag(i+i1,j+j1,k+k1,b) == 0_I4P) then
   !!                                 flag(i+i1,j+j1,k+k1,b) = w
   !!                              endif
   !!                           enddo
   !!                        enddo
   !!                     enddo
   !!                  else
   !!                     if (flag(i,j,k,b) == 0_I4P) then
   !!                        flag(i,j,k,b) = w
   !!                     endif
   !!                  endif
   !!               endif
   !!            enddo
   !!         enddo
   !!      enddo
   !!   enddo
   !!   do b=1, blocks_number
   !!      do k = 1, nk
   !!         do j = 1, nj
   !!            do i = 1, ni
   !!               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]
   !!               dist = sqrt(sq_norm(crossproduct(a=(cell_coord-V(w,:)),b=vec(w,:))))
   !!               if (w == 2) then
   !!                  if ((dotproduct(a=n2,b=cell_coord)+d2 <= -eps .or. dotproduct(a=n3,b=cell_coord)+d3 >= eps) .and. &
   !!                       flag(i,j,k,b) == w) then!aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
   !!                     flag(i,j,k,b) = 0_I4P
   !!                  endif
   !!               elseif (w == 4) then
   !!                  if ((dotproduct(a=n1,b=cell_coord)+d1 <= -eps .or. dotproduct(a=n4,b=cell_coord)+d4 >= eps) .and. &
   !!                       flag(i,j,k,b) == w) then !aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
   !!                     flag(i,j,k,b) = 0_I4P
   !!                  endif
   !!               endif
   !!               if (flag(i,j,k,b) == w) then
   !!                  selectcase (current_distribution)
   !!                  case (GAUSS_CURRENT_DISTRIBUTION)
   !!                     Gaussian(i,j,k,b) = gaussian_2D_ind(sigma = sigma*dx(b), r = dist)
   !!                  case (CONST_CURRENT_DISTRIBUTION)
   !!                     Gaussian(i,j,k,b) = 1/(d*dx(b))**2
   !!                  endselect
   !!               endif
   !!            enddo
   !!         enddo
   !!      enddo
   !!   enddo
   !!enddo
   !endassociate
   !endsubroutine set_rectangular_coil_quad_section_odd_v2
!
   !subroutine set_rectangular_coil_circular_section(self, physics, field, n)
   !class(prism_coil_object),     intent(inout) :: self                                                            !< Coils
   !type(field_object),           intent(inout) :: field                                                           !< Field object.
   !type(prism_physics_object),   intent(in)    :: physics                                                         !< Fluids physiscs.
   !integer(I4P),                 intent(in)    :: n                                                               !< Coil number.
   !integer(I4P),                 allocatable   :: flag(:,:,:,:)                                                   !< Flag per identificare se la spira passa per la cella
   !real(R8P),                    allocatable   :: Gaussian(:,:,:,:)                                               !< Matrice gaussiana per distribuzione corrente
   !real(R8P)                                   :: dmax                                                            !< Vincolo distanza massima dalla spira.
   !real(R8P)                                   :: c_c(3)                                                          !< Vettore posizione centro spira
   !real(R8P)                                   :: cell_coord(3)                                                   !< Vettore posizione centro cella
   !real(R8P)                                   :: vx(3), vy(3), vz(3)                                             !< Versori assi cartesiani
   !real(R8P)                                   :: V1(3), V2(3), V3(3), V4(3), V(4,3)                              !< Vertici rettangolo e relativa matrice
   !real(R8P)                                   :: V1_1(3), V2_1(3), V3_1(3), V4_1(3)
   !real(R8P)                                   :: V1_2(3), V2_2(3), V3_2(3), V4_2(3)
   !real(R8P)                                   :: v_l1(3), v_l2(3), vec(4,3)                                      !< Versori lati rettangolo (vale regola mano dx) e relativa matrice
   !real(R8P)                                   :: n1(3), d1, n2(3), d2, n3(3), d3, n4(3), d4                      !< Parametri piani diagonali perpendicolari a spira, per evitare sovrapposizioni
   !real(R8P)                                   :: kappa(3), K_rot(3,3), Id(3,3), theta                            !< Vettore, angolo e matrice di appoggio per formula Rodrigues + matrice identità
   !real(R8P)                                   :: Kquad(3,3), R(3,3)                                              !< Matrice K^2 e matrice rotazione tra vz e normale alla spira
   !real(R8P)                                   :: dist, prj_v(3),eps, d_real , dist_sec                                             !< Distanza punto retta e proiezione del punto sulla retta                                                                          !< Variabile utilizzata per definire direzione corrente
   !integer(I4P)                                :: b,i,j,k,w                                                       !< Counter.
   !integer(I4P)                                :: i1,j1,k1,w1                                                     !< Counter1.
   !integer(I4P)                                :: d_int
   !!associo per dati su posizioni delle celle e contatori
   !associate(blocks_number=>field%blocks_number, ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, ngc=>field%grid%ngc, &
   !         x_c => self%x_center(n), y_c => self%y_center(n), z_c => self%z_center(n),                                         &
   !         dx => field%dxyz(1,:), dy => field%dxyz(2,:), dz => field%dxyz(3,:), lx => self%lx(n), ly => self%ly(n),           &
   !         normal => self%normal(:,n), d => self%d(n), nb =>field%nb, current_distribution => self%current_distribution(n),   &
   !         x_cell => field%x_cell, y_cell => field%y_cell, z_cell => field%z_cell, sigma => self%sigma(n))
   !c_c = [ x_c, y_c, z_c ] !Vettore posizione centro spira
   !!vertici del rettangolo, lo costruisco come se avesse normale asse z e fosse centrato nell'origine;
   !!vertici in senso antiorario, partendo da in basso a sinistra; si ipotizza
   !!rotazione rispetto al centro, quindi, ad esempio, se la vuoi con normale lungo x
   !!con lato lungo lungo y AB sarà il lato corto (pensa a rotazione 3D). Con
   !!normale lungo y il problema non si pone
   !vx = [1._R8P,0._R8P,0._R8P]
   !vy = [0._R8P,1._R8P,0._R8P]
   !vz = [0._R8P,0._R8P,1._R8P]
   !V1 = [-lx/2, -ly/2, 0._R8P]
   !V2 = [+lx/2, -ly/2, 0._R8P]
   !V3 = [+lx/2, +ly/2, 0._R8P]
   !V4 = [-lx/2, +ly/2, 0._R8P]
!
   !! 1: interno
   !V1_1 = [-(lx/2 - d/2), -(ly/2 - d/2), 0._R8P]
   !V2_1 = [(lx/2  - d/2), -(ly/2 - d/2), 0._R8P]
   !V3_1 = [(lx/2  - d/2),  (ly/2 - d/2), 0._R8P]
   !V4_1 = [-(lx/2 - d/2), +(ly/2 - d/2), 0._R8P]
   !! 2: esterno
   !V1_2 = [-(lx/2 + d/2), -(ly/2 + d/2), 0._R8P]
   !V2_2 = [(lx/2  + d/2), -(ly/2 + d/2), 0._R8P]
   !V3_2 = [(lx/2  + d/2),  (ly/2 + d/2), 0._R8P]
   !V4_2 = [-(lx/2 + d/2), +(ly/2 + d/2), 0._R8P]
!
   !!            V4 ________________ V3     A
   !!              |                |       A
   !!              |                |       |
   !!              |                |       |
   !!              |                |       y
   !!              |________________|
   !!            V1                  V2     x --->>>
!
   !!calcolo rotazione tra i due vettori normali %OSS una rotazione di 180° dà problemi, qua ci metti un go to.
   !!isnan per il caso rotazione nulla rispetto
   !!a n // z
   !kappa = crossproduct(a=vz,b=normal)
!
   !if (any(normal /= vz )) then
   !    kappa = kappa/sqrt(sq_norm(kappa))
   !endif
!
   !theta = acos(dotproduct(a=vz,b=normal))
!
   !K_rot(1,1) = 0._R8P
   !K_rot(1,2) = -kappa(3)
   !K_rot(1,3) = kappa(2)
   !K_rot(2,1) = kappa(3)
   !K_rot(2,2) = 0._R8P
   !K_rot(2,3) = -kappa(1)
   !K_rot(3,1) = -kappa(2)
   !K_rot(3,2) = kappa(1)
   !K_rot(3,3) = 0._R8P
   !Kquad = matmul(K_rot,K_rot) !matrice K^2
   !Id(:,1) = vx
   !Id(:,2) = vy
   !Id(:,3) = vz
!
   !R = Id+sin(theta)*K_rot+(1-cos(theta))*Kquad !costruisco matrice di rotazione, ruoto i vettori posizione e poi traslo
   !                                            ! a seconda della posizione del centro richiesta
!
   !V1 = matmul(R,V1)+c_c
   !V2 = matmul(R,V2)+c_c
   !V3 = matmul(R,V3)+c_c
   !V4 = matmul(R,V4)+c_c
   !V(1,:) = V1
   !V(2,:) = V2
   !V(3,:) = V3
   !V(4,:) = V4 !genero matrice [V1; V2; V3; V4], quindi le colonne sono le coordinate x y z e le righe i vari vertici nell'ordine descritto prima
!
   !V1_1 = matmul(R,V1_1)+c_c
   !V2_1 = matmul(R,V2_1)+c_c
   !V3_1 = matmul(R,V3_1)+c_c
   !V4_1 = matmul(R,V4_1)+c_c
   !V1_2 = matmul(R,V1_2)+c_c
   !V2_2 = matmul(R,V2_2)+c_c
   !V3_2 = matmul(R,V3_2)+c_c
   !V4_2 = matmul(R,V4_2)+c_c
   !v_l1 = matmul(R,vx)
   !v_l1 = v_l1/sqrt(sq_norm(v_l1))
   !v_l2 = matmul(R,vy);
   !v_l2 = v_l2/sqrt(sq_norm(v_l2))
!
   !!matrice dei versori dei lati, generata ruotando i versori del rettangolo tramite la matrice di rotazione precedentemente calcolata
   !vec(1,:) = v_l1
   !vec(2,:) = v_l2
   !vec(3,:) = -v_l1
   !vec(4,:) = -v_l2
!
   !!calcolo piani perpendicolari alla spira su cui giacciono le due diagonali
   !!del rettangolo. Per convenzione, le normali puntano "verso i vertici" D e C
   !!piano 1, diagonale V1-V3 e normale verso V4
   !n1 = crossproduct(a=(V1_2-V1_1),b=normal)/sqrt(sq_norm(V1_2-V1_1)) !normale al piano
   !d1 = -dotproduct(a=n1,b=V1_1) !parametro d dell'equazione ax + by + cz + d1 = 0, con a b c coseni direttori della normale
   !n2 = crossproduct(a=(V2_1-V2_2),b=normal)/sqrt(sq_norm(V2_1-V2_2)) !normale al piano
   !d2 = -dotproduct(a=n2,b=V2_1) !parametro d dell'equazione ax + by + cz + d1 = 0, con a b c coseni direttori della normale
   !n3 = crossproduct(a=(V3_1-V3_2),b=normal)/sqrt(sq_norm(V3_1-V3_2)) !normale al piano
   !d3 = -dotproduct(a=n3,b=V3_1) !parametro d dell'equazione ax + by + cz + d1 = 0, con a b c coseni direttori della normale
   !n4 = crossproduct(a=(V4_2-V4_1),b=normal)/sqrt(sq_norm(V4_2-V4_1)) !normale al piano
   !d4 = -dotproduct(a=n4,b=V4_1) !parametro d dell'equazione ax + by + cz + d1 = 0, con a b c coseni direttori della normale
!
   !allocate(flag(1:ni,1:nj,1:nk,1:nb))
   !flag(:,:,:,:) = 0_I4P !inizializzo matrice flag a zero, per indicare che nessun lato passa per le celle
   !allocate(Gaussian(1:ni,1:nj,1:nk,1:nb))
   !Gaussian(:,:,:,:) = 0_I4P !inizializzo matrice gaussiana a zero
   !d_int = 10_I4P
   !eps = 1*10e-10_R8P
   !do w = 1, 4, 2
   !   do b=1, blocks_number
   !      dmax = d/2+eps
   !      ! chiamo funzione che restituisce le coordinate delle varie celle che compongono il blocco
   !      do k=1+d_int, nk-d_int
   !         do j=1+d_int, nj-d_int
   !            do i=1+d_int, ni-d_int
   !               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]
   !               dist = sqrt(sq_norm(crossproduct(a=(cell_coord-V(w,:)),b=vec(w,:)))) !Distanza punto retta |(P-A) x v| / |v| con A punto sulla retta e v versore della retta
   !               prj_v = V(w,:)+dotproduct(a=(cell_coord-V(w,:)),b=vec(w,:))*vec(w,:); !Formula proiezione di un punto su una retta con A punto della retta A+[(P-A)*v]*v
   !               !Con questo if inziale traccio tutta la "linea mediana" della spira
   !               if (dist <= dx(b)/2 .and. &
   !                  prj_v(1) <= maxval(V(:,1)) + dmax .and. &
   !                  prj_v(2) <= maxval(V(:,2)) + dmax .and. prj_v(3) <= maxval(V(:,3)) + dmax .and. &
   !                  minval(V(:,1))-dmax <= prj_v(1) .and. minval(V(:,2))-dmax <= prj_v(2) .and. &
   !                  minval(V(:,3))-dmax <= prj_v(3) ) then
   !                  if (d_int /= 0_I4P) then !Qua ci va algoritmo per indicare la sezione.
   !                     do i1 = -d_int, d_int
   !                        do j1 = -d_int, d_int
   !                           do k1 = -d_int, d_int
   !                              dist_sec = sqrt((x_cell(i,b)-x_cell(i+i1,b))**2.0_R8P + (y_cell(j,b)-y_cell(j+j1,b))**2.0_R8P + &
   !                                                 (z_cell(k,b)-z_cell(k+k1,b))**2.0_R8P)
   !                              if (flag(i+i1,j+j1,k+k1,b) == 0_I4P .and. dist_sec <= d/2+eps) then
   !                                 flag(i+i1,j+j1,k+k1,b) = w
   !                              endif
   !                           enddo
   !                        enddo
   !                     enddo
   !                  else
   !                     if (flag(i,j,k,b) == 0_I4P) then
   !                        flag(i,j,k,b) = w
   !                     endif
   !                  endif
   !               endif
   !            enddo
   !         enddo
   !      enddo
   !   enddo
   !   do b=1, blocks_number
   !      do k = 1, nk
   !         do j = 1, nj
   !            do i = 1, ni
   !               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]
   !               dist = sqrt(sq_norm(crossproduct(a=(cell_coord-V(w,:)),b=vec(w,:))))
   !               if (w == 1) then
   !                  if ((dotproduct(a=n1,b=cell_coord)+d1 > eps .or. dotproduct(a=n2,b=cell_coord)+d2 > eps) .and. &
   !                       flag(i,j,k,b) == w) then !aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
   !                     flag(i,j,k,b) = 0_I4P
   !                     !print *, w
   !                  endif
   !               elseif (w == 3) then
   !                  if ((dotproduct(a=n3,b=cell_coord)+d3 < -eps .or. dotproduct(a=n4,b=cell_coord)+d4 < -eps) .and. &
   !                       flag(i,j,k,b) == w) then !aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
   !                     flag(i,j,k,b) = 0_I4P
   !                     !print *, w
   !                  endif
   !               endif
   !               if (flag(i,j,k,b) == w) then
   !                  selectcase (current_distribution)
   !                  case (GAUSS_CURRENT_DISTRIBUTION)
   !                     Gaussian(i,j,k,b) = gaussian_2D_ind(sigma = d/6, r = dist)
   !                  case (CONST_CURRENT_DISTRIBUTION)
   !                     Gaussian(i,j,k,b) = 1/(PI*d**2/4)
   !                  endselect
   !               endif
   !            enddo
   !         enddo
   !      enddo
   !   enddo
   !enddo
   !do w = 2, 4, 2
   !   do b=1, blocks_number
   !      !dmax = 0.0_R8P!*d_real+eps
   !      ! chiamo funzione che restituisce le coordinate delle varie celle che compongono il blocco
   !      do k=1+d_int, nk-d_int
   !         do j=1+d_int, nj-d_int
   !            do i=1+d_int, ni-d_int
   !               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]
   !               dist = sqrt(sq_norm(crossproduct(a=(cell_coord-V(w,:)),b=vec(w,:)))) !Distanza punto retta |(P-A) x v| / |v| con A punto sulla retta e v versore della retta
   !               prj_v = V(w,:)+dotproduct(a=(cell_coord-V(w,:)),b=vec(w,:))*vec(w,:); !Formula proiezione di un punto su una retta con A punto della retta A+[(P-A)*v]*v
   !               if (dist <= dx(b)/2 .and. &
   !                 prj_v(1) <= maxval(V(:,1)) + dmax .and. &
   !                 prj_v(2) <= maxval(V(:,2)) + dmax .and. prj_v(3) <= maxval(V(:,3)) + dmax .and. &
   !                 minval(V(:,1))-dmax <= prj_v(1) .and. minval(V(:,2))-dmax <= prj_v(2) .and. &
   !                 minval(V(:,3))-dmax <= prj_v(3) ) then
   !                  if (d_int /= 0_I4P) then !Qua ci va algoritmo per indicare la sezione.
   !                     !Per alleggerire il calcolo itero sul cubo che racchiude la cella considerando
   !                     !solo una parte degli indici. Magari calcola quanti indici in funzione del rapporto
   !                     !d/dx. La cella "centrale" è indicata dagli indici i j k b
   !                     do i1 = -d_int, d_int
   !                        do j1 = -d_int, d_int
   !                           do k1 = -d_int, d_int
   !                              dist_sec = sqrt((x_cell(i,b)-x_cell(i+i1,b))**2.0_R8P + (y_cell(j,b)-y_cell(j+j1,b))**2.0_R8P + &
   !                                                 (z_cell(k,b)-z_cell(k+k1,b))**2.0_R8P)
   !                              if (flag(i+i1,j+j1,k+k1,b) == 0_I4P .and. dist_sec <= d/2+eps) then
   !                                 flag(i+i1,j+j1,k+k1,b) = w
   !                              endif
   !                           enddo
   !                        enddo
   !                     enddo
   !                  else
   !                     if (flag(i,j,k,b) == 0_I4P) then
   !                        flag(i,j,k,b) = w
   !                     endif
   !                  endif
   !               endif
   !            enddo
   !         enddo
   !      enddo
   !   enddo
   !   do b=1, blocks_number
   !      do k = 1, nk
   !         do j = 1, nj
   !            do i = 1, ni
   !               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]
   !               dist = sqrt(sq_norm(crossproduct(a=(cell_coord-V(w,:)),b=vec(w,:))))
   !               if (w == 2) then
   !                  if ((dotproduct(a=n2,b=cell_coord)+d2 <= -eps .or. dotproduct(a=n3,b=cell_coord)+d3 >= eps) .and. &
   !                       flag(i,j,k,b) == w) then!aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
   !                     flag(i,j,k,b) = 0_I4P
   !                     !print *, w
   !                  endif
   !               elseif (w == 4) then
   !                  if ((dotproduct(a=n1,b=cell_coord)+d1 <= -eps .or. dotproduct(a=n4,b=cell_coord)+d4 >= eps) .and. &
   !                       flag(i,j,k,b) == w) then !aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
   !                     flag(i,j,k,b) = 0_I4P
   !                     !print *, w
   !                  endif
   !               endif
   !               if (flag(i,j,k,b) == w) then
   !                  selectcase (current_distribution)
   !                  case (GAUSS_CURRENT_DISTRIBUTION)
   !                     Gaussian(i,j,k,b) = gaussian_2D_ind(sigma = d/6, r = dist)
   !                  case (CONST_CURRENT_DISTRIBUTION)
   !                     Gaussian(i,j,k,b) = 1/(PI*d**2/4)
   !                  endselect
   !               endif
   !            enddo
   !         enddo
   !      enddo
   !   enddo
   !enddo
   !!Ho un flag pari a 1 2 3 4 nelle celle per cui passa uno dei dati della spira. La direzione della corrente è
   !!Ceorente con quella dei versori dei lati precedentemente descritti
   !do b=1, blocks_number
   !   do k=1, nk
   !      do j=1, nj
   !         do i=1, ni
   !            if (flag(i,j,k,b) /= 0 .and. self%coil_flag(i,j,k,b) == 0_I4P) then
!
   !               self%J_vec(1:3,i,j,k,b) = vec(flag(i,j,k,b),:)
   !               self%J_vec(4,i,j,k,b) = Gaussian(i,j,k,b)
!
   !               if (abs(self%J_vec(1,i,j,k,b)) < 1.0e-10_R8P) then
   !                  self%J_vec(1,i,j,k,b) = 0._R8P
   !               endif
   !               if (abs(self%J_vec(2,i,j,k,b)) < 1.0e-10_R8P) then
   !                  self%J_vec(2,i,j,k,b) = 0._R8P
   !               endif
   !               if (abs(self%J_vec(3,i,j,k,b)) < 1.0e-10_R8P) then
   !                  self%J_vec(3,i,j,k,b) = 0._R8P
   !               endif
!
   !               self%coil_flag(i,j,k,b) = n
!
   !            endif
   !         enddo
   !      enddo
   !   enddo
   !enddo
   !endassociate
   !endsubroutine set_rectangular_coil_circular_section
!
   !function dotproduct(a, b) result(dot)
   !!< Compute the scalar (dot) product.
   !real(R8P), intent(in) :: a(3) !< Left hand side.
   !real(R8P), intent(in) :: b(3) !< Left hand side.
   !real(R8P)             :: dot  !< Dot product.
!
   !dot = (a(1) * b(1)) + (a(2) * b(2)) + (a(3) * b(3))
   !endfunction dotproduct
!
   !function crossproduct(a, b) result(cross)
   !real(R8P), intent(in) :: a(3)     !< Left hand side.
   !real(R8P), intent(in) :: b(3)     !< Left hand side.
   !real(R8P)             :: cross(3) !< Cross product.
!
   !cross(1) = (a(2) * b(3)) - (a(3) * b(2))
   !cross(2) = (a(3) * b(1)) - (a(1) * b(3))
   !cross(3) = (a(1) * b(2)) - (a(2) * b(1))
   !endfunction crossproduct
!
   !function sq_norm(a) result(sq)
   !!< Return the square of the norm of vector.
   !real(R8P), intent(in)  :: a(3)     !< Input vector
   !real(R8P)              :: sq       !< Square norm of input
!
   !sq = (a(1) * a(1)) + (a(2) * a(2)) + (a(3) * a(3))
   !endfunction sq_norm
!
   !function gaussian_2D_ind(sigma, r) result(f)
   !! Compute the 2D Gaussian function for indipendent variables x,y linked with the distance from the cell centre to the
   !! coil wire in the plain perpendicular to the coil direction. Mu is 0 (obviously)
   !real(R8P), intent(in)  :: sigma, r  !< Standard deviation and distance from the coil wire
   !real(R8P)              :: f         !< Density of probability function value
!
   !f = exp(-0.5 * (r / sigma)**2) / (2.0 * PI * sigma**2)
   !endfunction gaussian_2D_ind
