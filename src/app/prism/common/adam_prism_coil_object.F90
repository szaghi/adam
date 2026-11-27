!< ADAM, PRISM coil source definition, CPU backend.
module adam_prism_coil_object
!< ADAM, PRISM coil source definition, CPU backend.

! ADAM singleton objects
use :: adam_field_object, only : field_object
use :: adam_grid_object,  only : grid_object
use :: adam_ini_library, only : ini_check, ini_count_sections, ini_get
use :: adam_mpih_global,  only : mpih
! PRISM modules
use :: adam_prism_physics_object, only : prism_physics_object, is_adim_model
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
public :: COIL_TYPE_HELICON
public :: CURRENT_TYPE_AC
public :: CURRENT_TYPE_DC
public :: prism_coil_object
public :: NORMAL_P_X
public :: NORMAL_M_X
public :: NORMAL_P_Y
public :: NORMAL_M_Y
public :: NORMAL_P_Z
public :: NORMAL_M_Z
public :: COIL_SLOT_COMMITTED
public :: COIL_SLOT_STAGE

character(len=11), parameter :: INI_SECTION_NAME="coils_input"        !< INI (config) file section name containing coils configs.
character(len=11), parameter :: COIL_TYPE_RECTANGULAR="rectangular"   !< Rectangular shape coil.
character(len=8),  parameter :: COIL_TYPE_CIRCULAR="circular"         !< Circular shape coil.
character(len=8),  parameter :: COIL_TYPE_SOLENOID="solenoid"         !< Solenoid shape coil.
character(len=8),  parameter :: COIL_TYPE_HELICON="helicon"           !< Helicon shape coil.
character(len=10), parameter :: CURRENT_TYPE_DC="DC_current"          !< DC current.
character(len=10), parameter :: CURRENT_TYPE_AC="AC_current"          !< AC current
character(len=2),  parameter :: NORMAL_P_X="+x"                       !< Normal versor along positive x axis.
character(len=2),  parameter :: NORMAL_M_X="-x"                       !< Normal versor along negative x axis.
character(len=2),  parameter :: NORMAL_P_Y="+y"                       !< Normal versor along positive y axis.
character(len=2),  parameter :: NORMAL_M_Y="-y"                       !< Normal versor along negative y axis.
character(len=2),  parameter :: NORMAL_P_Z="+z"                       !< Normal versor along positive z axis.
character(len=2),  parameter :: NORMAL_M_Z="-z"                       !< Normal versor along negative z axis.
integer(I4P),      parameter :: COIL_SLOT_COMMITTED=1_I4P              !< Plasma-current slot of the committed state q.
integer(I4P),      parameter :: COIL_SLOT_STAGE=2_I4P                  !< Plasma-current slot of the stage state.

type :: prism_coil_object
   !< ADAM, PRISM coil source definition, CPU backend.
   character(len=99), allocatable :: coil_type(:)                                !< Coil type.
   character(len=99), allocatable :: current_type(:)                             !< Current type.
   character(len=2 ), allocatable :: normal(:) !< Versore normale alla spira, che identifica anche verso
   real(R8P),         allocatable :: A(:)                                        !< Current amplitude (A)
   real(R8P),         allocatable :: coil_amplitude(:)      
                                                                                 !< Current amplitude (A) for each coil, corrected
                                                                                 !< for Gaussian current distribution.
   real(R8P),         allocatable :: f(:)                                        !< Current frequency, if AC (Hz)
   real(R8P),         allocatable :: phase(:)                                    !< Current initial phase, if AC
   real(R8P),         allocatable :: x_center(:), y_center(:), z_center(:)       !< Coil center
   real(R8P),         allocatable :: lx(:), ly(:)                                !< Rectangle's sizes (if rectangular coil)
   real(R8P),         allocatable :: r_coil(:)                                   !< Circle's radius (if circular coil)
   real(R8P),         allocatable :: l_solenoid(:)                               !< Solenoid length (if solenoidal coil)
   real(R8P),         allocatable :: windings(:)                                 !< Windings number (if solenoidal coil)
   integer(I4P),      allocatable :: N_points(:) !< Number of points for coil representation (for helicon coil setting)
   real(R8P),         allocatable :: x_points(:,:), y_points(:,:), z_points(:,:) !< Coil points coordinates (for helicon coil
      !< setting)
   real(R8P),         allocatable :: rectangular_current_flux(:,:)
                                                                                 !< Per-side unit-current flux for rectangular
                                                                                 !< coils.
   real(R8P),         allocatable :: helicon_current_flux(:,:)                   !< Per-segment unit-current flux for helicon coils.
   real(R8P),         allocatable :: sigma(:)                                    !< Gaussian current distribution sigma
   real(R8P),         allocatable :: J_vec(:,:,:,:,:,:) !< Matrice contenente versori corrente spire (se assente
   real(R8P)                      :: td                                          !< Delay di accensione della spira
   integer(I4P)                   :: circular_coils_number=0_I4P                 !< Number of circular coils
   integer(I4P)                   :: rectangular_coils_number=0_I4P              !< Number of rectangular coils
   integer(I4P)                   :: solenoid_coils_number=0_I4P                 !< Number of solenoid coils
   integer(I4P)                   :: helicon_coils_number=0_I4P                  !< Number of helicon coils
   integer(I4P)                   :: total_coils_number=0_I4P                    !< Number of coils
   type(string), allocatable      :: j_vec_name(:,:)                             !< J vec names.
   ! Plasma current deposited on the coil support. The coil stamp clears the J rows of the support cells and writes the
   ! coil current there; where the PIC deposit also put current, the plasma part is kept here and summed by the stamp
   ! (detect_plasma and capture_plasma after the deposit, restore_plasma inside the stamp), so the stamp stays idempotent.
   integer(I4P), allocatable      :: support(:,:)                                !< Cells (b,i,j,k) with J_vec /= 0: interior
                                                                                 !< first, then ghost.
   integer(I4P)                   :: support_interior=0_I4P                      !< Number of interior support cells.
   logical                        :: support_valid=.false.                       !< Support built from the current J_vec.
   real(R8P),    allocatable      :: plasma_j(:,:,:)                             !< Plasma current on the support (3,n,slot).
   logical                        :: plasma_j_active(2)=.false.                  !< Slot holds a plasma current.
   integer(I4P)                   :: overlap_first_it=-1_I4P                     !< First iteration with plasma on coil cells.
   integer(I4P)                   :: overlap_last_it=-1_I4P                      !< Last iteration with plasma on coil cells.
   character(len=99)              :: overlap_first_context=''                    !< Call site of the first overlap.
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
      procedure, pass(self) :: load_coils_number                         !< Load coils number from file.
      procedure, pass(self) :: load_from_file                            !< Load config from file.
      procedure, pass(self) :: adimensionalize_coils_parameters          !< Adimensionalize coils parameters.
      procedure, pass(self) :: build_support                             !< Build the coil support and its plasma slots.
      procedure, pass(self) :: invalidate_support                        !< Mark the support stale (J_vec recomputed).
      procedure, pass(self) :: record_overlap                            !< Record plasma deposited on coil cells.
      procedure, pass(self) :: capture_plasma                            !< Keep the plasma current of the coil support.
      procedure, pass(self) :: detect_plasma                             !< Plasma on the interior coil support.
      procedure, pass(self) :: restore_plasma                            !< Write the kept plasma current back.
      procedure, pass(self) :: support_cells                             !< Interior cells where J_vec /= 0.
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
   type(field_object), intent(in)          :: field !< Field (sibling realm component, threaded in).
   type(grid_object),  intent(in), target  :: grid !< Grid (sibling realm component, threaded in).
   integer(I4P)                            :: c    !< Counter.

   associate(ngc=>self%ngc,ni=>self%ni,nj=>self%nj,nk=>self%nk,nb=>field%nb,total_coils_number=>self%total_coils_number)

   allocate(self%r_coil              (0:total_coils_number))    ; self%r_coil = 0.0_R8P
   allocate(self%ly                  (0:total_coils_number))    ; self%ly = 0.0_R8P
   allocate(self%lx                  (0:total_coils_number))    ; self%lx = 0.0_R8P
   allocate(self%l_solenoid          (0:total_coils_number))    ; self%l_solenoid = 0.0_R8P
   allocate(self%windings            (0:total_coils_number))    ; self%windings = 0.0_R8P
   allocate(self%N_points            (0:total_coils_number))    ; self%N_points = 0_I4P
   allocate(self%x_points            (0:total_coils_number,50)) ; self%x_points = 0.0_R8P
   allocate(self%y_points            (0:total_coils_number,50)) ; self%y_points = 0.0_R8P
   allocate(self%z_points            (0:total_coils_number,50)) ; self%z_points = 0.0_R8P
   allocate(self%rectangular_current_flux(4, 0:total_coils_number)) ; self%rectangular_current_flux = 0.0_R8P
   allocate(self%helicon_current_flux    (50,0:total_coils_number)) ; self%helicon_current_flux     = 0.0_R8P
   allocate(self%sigma               (0:total_coils_number))    ; self%sigma = 0.0_R8P
   allocate(self%x_center            (0:total_coils_number))    ; self%x_center = 0.0_R8P
   allocate(self%y_center            (0:total_coils_number))    ; self%y_center = 0.0_R8P
   allocate(self%z_center            (0:total_coils_number))    ; self%z_center = 0.0_R8P
   allocate(self%coil_type           (0:total_coils_number))    ; self%coil_type = ' '
   allocate(self%current_type        (0:total_coils_number))    ; self%current_type = ' '
   allocate(self%normal              (0:total_coils_number))    ; self%normal = ' '
   allocate(self%A                   (0:total_coils_number))    ; self%A = 0.0_R8P
   allocate(self%coil_amplitude      (0:total_coils_number))    ; self%coil_amplitude = 0.0_R8P
   allocate(self%f                   (0:total_coils_number))    ; self%f = 0.0_R8P
   allocate(self%phase               (0:total_coils_number))    ; self%phase = 0.0_R8P

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
   integer(I4P)                         :: r, p             !< Counter.

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
         case(COIL_TYPE_HELICON)
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
         desc = desc//NL//mpih%myrankstr//'    N_points: '//trim(str(self%N_points(r)))
         do p=1, self%N_points(r)
            desc = desc//NL//mpih%myrankstr//'      Point '//trim(str(p,.true.))//': ('// &
               trim(str(self%x_points(r,p)))//', '//trim(str(self%y_points(r,p)))//', '// &
               trim(str(self%z_points(r,p)))//')'
         enddo
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

   subroutine load_coils_number(self, file_parameters)
   !< Load the number of coils of each kind from file, and their total.
   !<
   !< Separate from `load_from_file` because it needs no grid or field: the solver reads it before the blocks budget
   !< is computed, since the coils current density is a block-sized field per coil.
   class(prism_coil_object), intent(inout)        :: self            !< Coils.
   type(file_ini),           intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   integer(I4P)                                   :: counts(4)       !< Coil sections by type.


   ! each count defaults to the number of [coils_input_coil_N] sections of that coil_type (issue #80, F5)
   call count_coil_sections
   call ini_get(file_parameters, INI_SECTION_NAME, 'circular_coils_number', self%circular_coils_number, &
                default=counts(1))
   call ini_get(file_parameters, INI_SECTION_NAME, 'rectangular_coils_number', self%rectangular_coils_number, &
                default=counts(2))
   call ini_get(file_parameters, INI_SECTION_NAME, 'solenoid_coils_number', self%solenoid_coils_number, &
                default=counts(3))
   call ini_get(file_parameters, INI_SECTION_NAME, 'helicon_coils_number', self%helicon_coils_number, &
                default=counts(4))
   self%total_coils_number = self%circular_coils_number + self%rectangular_coils_number &
                              + self%solenoid_coils_number + self%helicon_coils_number
   call ini_check()
   contains
      subroutine count_coil_sections
      !< Count the coil sections by coil_type.
      character(99) :: buff_char !< Option character buffer.
      integer(I4P)  :: c         !< Counter.

      counts = 0_I4P
      do c=1, ini_count_sections(file_parameters, INI_SECTION_NAME//'_coil_')
         buff_char = ''
         call ini_get(file_parameters, INI_SECTION_NAME//'_coil_'//trim(str(c, .true.)), 'coil_type', buff_char)
         select case(trim(buff_char))
         case(COIL_TYPE_CIRCULAR)    ; counts(1) = counts(1) + 1_I4P
         case(COIL_TYPE_RECTANGULAR) ; counts(2) = counts(2) + 1_I4P
         case(COIL_TYPE_SOLENOID)    ; counts(3) = counts(3) + 1_I4P
         case(COIL_TYPE_HELICON)     ; counts(4) = counts(4) + 1_I4P
         endselect
      enddo
      endsubroutine count_coil_sections
   endsubroutine load_coils_number

   subroutine load_from_file(self, field, grid, physics, file_parameters)
   !< Load config from file.
   class(prism_coil_object),   intent(inout)        :: self            !< coils.
   type(field_object),         intent(in)           :: field           !< Field (sibling realm component, threaded in).
   type(grid_object),          intent(in), target   :: grid            !< Grid (sibling realm component, threaded in).
   type(prism_physics_object), intent(in)           :: physics         !< PRISM physics.
   type(file_ini),             intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   character(:), allocatable                        :: sname           !< Section name.
   integer(I4P)                                     :: i,j             !< Counter.
   integer(I4P)                                     :: error           !< Error status.
   character(99)                                    :: buff_char       !< Option character buffer.


   call self%load_coils_number(file_parameters=file_parameters)
   if (self%total_coils_number==0_I4P) return

   call associate_adam_data

   call self%allocate_coil(field=field, grid=grid)

   call ini_get(file_parameters, INI_SECTION_NAME, 'time_delay', self%td)

   call ini_check()
   do i=1, self%total_coils_number
      sname = INI_SECTION_NAME//'_coil_'//trim(str(i,.true.))

      call ini_get(file_parameters, sname, 'coil_type', buff_char)
      self%coil_type(i) = trim(buff_char)
      self%coil_type(i) = trim(self%coil_type(i))

      call ini_get(file_parameters, sname, 'current_type', buff_char)
      self%current_type(i) = trim(buff_char)
      self%current_type(i) = trim(self%current_type(i))

      call ini_get(file_parameters, sname, 'sigma', self%sigma(i))

      call ini_check()
      select case(self%coil_type(i))
      case(COIL_TYPE_CIRCULAR)
         call ini_get(file_parameters, sname, 'r_coil', self%r_coil(i))

         call ini_get(file_parameters, sname, 'x_center', self%x_center(i))

         call ini_get(file_parameters, sname, 'y_center', self%y_center(i))

         call ini_get(file_parameters, sname, 'z_center', self%z_center(i))

         call ini_get(file_parameters, sname, 'normal', buff_char)
         self%normal(i) = trim(buff_char)
         self%normal(i) = trim(self%normal(i))

         call ini_check()
      case(COIL_TYPE_RECTANGULAR)
         call ini_get(file_parameters, sname, 'lx', self%lx(i))

         call ini_get(file_parameters, sname, 'ly', self%ly(i))

         call ini_get(file_parameters, sname, 'x_center', self%x_center(i))

         call ini_get(file_parameters, sname, 'y_center', self%y_center(i))

         call ini_get(file_parameters, sname, 'z_center', self%z_center(i))

         call ini_get(file_parameters, sname, 'normal', buff_char)
         self%normal(i) = trim(buff_char)
         self%normal(i) = trim(self%normal(i))

         call ini_check()
      case(COIL_TYPE_SOLENOID)
         call ini_get(file_parameters, sname, 'r_coil', self%r_coil(i))

         call ini_get(file_parameters, sname, 'l_solenoid', self%l_solenoid(i))

         call ini_get(file_parameters, sname, 'windings', self%windings(i))

         call ini_get(file_parameters, sname, 'x_center', self%x_center(i))

         call ini_get(file_parameters, sname, 'y_center', self%y_center(i))

         call ini_get(file_parameters, sname, 'z_center', self%z_center(i))

         call ini_get(file_parameters, sname, 'normal', buff_char)
         self%normal(i) = trim(buff_char)
         self%normal(i) = trim(self%normal(i))

         call ini_check()
      case(COIL_TYPE_HELICON)
         call ini_get(file_parameters, sname, 'x_center', self%x_center(i))

         call ini_get(file_parameters, sname, 'y_center', self%y_center(i))

         call ini_get(file_parameters, sname, 'z_center', self%z_center(i))

         call ini_get(file_parameters, sname, 'r_coil', self%r_coil(i))

         call ini_get(file_parameters, sname, 'normal', buff_char)
         self%normal(i) = trim(buff_char)
         self%normal(i) = trim(self%normal(i))

         call ini_get(file_parameters, sname, 'N_points', self%N_points(i))
         call ini_check()
         if (self%N_points(i) > 50.0_R8P) then
            call mpih%error_stop(msg=': too many points for ['//sname//'].(N_points). Max allowed is 50.')
         end if

         do j=1, self%N_points(i)
            call ini_get(file_parameters, sname, 'x_point_'//trim(str(j,.true.)), self%x_points(i,j))
            call ini_get(file_parameters, sname, 'y_point_'//trim(str(j,.true.)), self%y_points(i,j))
            call ini_get(file_parameters, sname, 'z_point_'//trim(str(j,.true.)), self%z_points(i,j))
         enddo
         call ini_check()
      endselect

      select case(self%current_type(i))
      case(CURRENT_TYPE_DC)

         call ini_get(file_parameters, sname, 'Amplitude', self%A(i))

         call ini_check()
      case(CURRENT_TYPE_AC)

         call ini_get(file_parameters, sname, 'Amplitude', self%A(i))

         call ini_get(file_parameters, sname, 'Frequency', self%f(i))

         call ini_get(file_parameters, sname, 'Phase', self%phase(i))

         call ini_check()
      endselect
   enddo

   if (is_adim_model(physics%physical_model)) &
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
   self%x_points(:,:) = self%x_points(:,:)/physics%L0
   self%y_points(:,:) = self%y_points(:,:)/physics%L0
   self%z_points(:,:) = self%z_points(:,:)/physics%L0
   !Adimensionalization of current parameters
   self%A(:) = self%A(:)/physics%I0
   !Adimensionalization of frequency parameters
   self%f(:) = self%f(:)*physics%T0
   !Adimensionalization of time parameters
   self%td = self%td/physics%T0
   endsubroutine adimensionalize_coils_parameters

   subroutine build_support(self, blocks_number)
   !< Build the coil support from the current J_vec and (re)allocate its plasma slots, empty.
   class(prism_coil_object), intent(inout) :: self          !< Coils.
   integer(I4P),             intent(in)    :: blocks_number !< Number of blocks in use.

   call self%support_cells(blocks_number=blocks_number, cells=self%support, interior=self%support_interior)
   if (allocated(self%plasma_j)) deallocate(self%plasma_j)
   allocate(self%plasma_j(1:3, 1:size(self%support, dim=2), 1:2))
   self%plasma_j = 0._R8P
   self%plasma_j_active = .false.
   self%support_valid = .true.
   endsubroutine build_support

   subroutine invalidate_support(self)
   !< Mark the support stale: J_vec was recomputed (initialization, AMR), the support and the kept plasma current with it.
   class(prism_coil_object), intent(inout) :: self !< Coils.

   self%support_valid = .false.
   self%plasma_j_active = .false.
   endsubroutine invalidate_support

   subroutine record_overlap(self, it, context)
   !< Record that a PIC deposit put plasma on coil cells: a warning at the first time (repeated at the end of the run by
   !< prism_common_object%report_plasma_on_coils), the iteration range kept.
   class(prism_coil_object), intent(inout) :: self    !< Coils.
   integer(I4P),             intent(in)    :: it      !< Current iteration.
   character(*),             intent(in)    :: context !< Call site.

   if (self%overlap_first_it < 0_I4P) then
      self%overlap_first_it = it
      self%overlap_first_context = context
      if (mpih%myrank == 0_I4P) &
         call mpih%print_message('warning: '//trim(context)//' deposited plasma on coil cells at iteration '// &
                                 trim(str(it, .true.))//': plasma and coil currents are summed there (repeated at '// &
                                 'the end of the run)')
   endif
   self%overlap_last_it = it
   endsubroutine record_overlap

   subroutine detect_plasma(self, q, var_jx, var_jy, var_jz, check_current, check_charge, blocks_number, &
                            overlap_current, overlap_charge)
   !< Right after a PIC deposit into q: plasma current (and charge) on the interior cells of the coil support (the deposit
   !< owns them after its fold; the ghost ones may hold an earlier state). Local to this rank.
   class(prism_coil_object), intent(inout) :: self            !< Coils.
   real(R8P),                intent(in)    :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Deposited field.
   integer(I4P),             intent(in)    :: var_jx          !< Index of Jx.
   integer(I4P),             intent(in)    :: var_jy          !< Index of Jy.
   integer(I4P),             intent(in)    :: var_jz          !< Index of Jz.
   logical,                  intent(in)    :: check_current   !< The J rows hold a fresh deposit.
   logical,                  intent(in)    :: check_charge    !< The charge row holds a fresh deposit.
   integer(I4P),             intent(in)    :: blocks_number   !< Number of blocks in use.
   logical,                  intent(out)   :: overlap_current !< Plasma current on the support.
   logical,                  intent(out)   :: overlap_charge  !< Plasma charge on the support.
   integer(I4P)                            :: nv_q            !< Width of q (charge is the last row).
   integer(I4P)                            :: c, b, i, j, k   !< Counters.

   overlap_current = .false.
   overlap_charge  = .false.
   if (self%total_coils_number <= 0_I4P) return
   if (.not.self%support_valid) call self%build_support(blocks_number=blocks_number)
   nv_q = int(size(q, dim=1), I4P)
   do c=1, self%support_interior
      b = self%support(1,c) ; i = self%support(2,c) ; j = self%support(3,c) ; k = self%support(4,c)
      if (check_current) overlap_current = overlap_current .or. (q(var_jx,i,j,k,b) /= 0._R8P) .or. &
                                           (q(var_jy,i,j,k,b) /= 0._R8P) .or. (q(var_jz,i,j,k,b) /= 0._R8P)
      if (check_charge) overlap_charge = overlap_charge .or. (q(nv_q,i,j,k,b) /= 0._R8P)
   enddo
   endsubroutine detect_plasma

   subroutine capture_plasma(self, q, slot, var_jx, var_jy, var_jz, active)
   !< Keep the J rows of the support cells (interior and ghost, the ghosts refreshed by the caller) in `slot`, or empty
   !< the slot (no plasma current on the support: the stamp then works exactly as before).
   class(prism_coil_object), intent(inout) :: self   !< Coils.
   real(R8P),                intent(in)    :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Deposited field.
   integer(I4P),             intent(in)    :: slot   !< COIL_SLOT_COMMITTED or COIL_SLOT_STAGE.
   integer(I4P),             intent(in)    :: var_jx !< Index of Jx.
   integer(I4P),             intent(in)    :: var_jy !< Index of Jy.
   integer(I4P),             intent(in)    :: var_jz !< Index of Jz.
   logical,                  intent(in)    :: active !< Plasma current on the support (of some rank).
   integer(I4P)                            :: c, b, i, j, k !< Counters.

   if (.not.self%support_valid) return
   self%plasma_j_active(slot) = active
   if (.not.active) return
   do c=1, size(self%support, dim=2)
      b = self%support(1,c) ; i = self%support(2,c) ; j = self%support(3,c) ; k = self%support(4,c)
      self%plasma_j(1,c,slot) = q(var_jx,i,j,k,b)
      self%plasma_j(2,c,slot) = q(var_jy,i,j,k,b)
      self%plasma_j(3,c,slot) = q(var_jz,i,j,k,b)
   enddo
   endsubroutine capture_plasma

   subroutine restore_plasma(self, q, slot, var_jx, var_jy, var_jz)
   !< Inside the coil stamp, after the support was cleared: write back the plasma current kept in `slot` (nothing when the
   !< slot is empty), the coil current is then added on top.
   class(prism_coil_object), intent(in)    :: self   !< Coils.
   real(R8P),                intent(inout) :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field.
   integer(I4P),             intent(in)    :: slot   !< COIL_SLOT_COMMITTED or COIL_SLOT_STAGE.
   integer(I4P),             intent(in)    :: var_jx !< Index of Jx.
   integer(I4P),             intent(in)    :: var_jy !< Index of Jy.
   integer(I4P),             intent(in)    :: var_jz !< Index of Jz.
   integer(I4P)                            :: c, b, i, j, k !< Counters.

   if (.not.self%support_valid) return
   if (.not.self%plasma_j_active(slot)) return
   do c=1, size(self%support, dim=2)
      b = self%support(1,c) ; i = self%support(2,c) ; j = self%support(3,c) ; k = self%support(4,c)
      q(var_jx,i,j,k,b) = self%plasma_j(1,c,slot)
      q(var_jy,i,j,k,b) = self%plasma_j(2,c,slot)
      q(var_jz,i,j,k,b) = self%plasma_j(3,c,slot)
   enddo
   endsubroutine restore_plasma

   subroutine support_cells(self, blocks_number, cells, interior)
   !< Cells (b,i,j,k) where some coil has J_vec /= 0, the interior ones first and then the ghost ones, in a fixed order
   !< (the FNL twin pushes the same list).
   class(prism_coil_object),  intent(in)  :: self          !< Coils.
   integer(I4P),              intent(in)  :: blocks_number !< Number of blocks in use.
   integer(I4P), allocatable, intent(out) :: cells(:,:)    !< Support cells (4,n).
   integer(I4P),              intent(out) :: interior      !< Number of interior support cells (the first ones).
   integer(I4P)                           :: n, b, i, j, k !< Counters.
   integer(I4P)                           :: sweep         !< Count, then fill.
   integer(I4P)                           :: part          !< 1 interior cells, 2 ghost cells.
   logical                                :: is_interior   !< The cell is an interior one.

   interior = 0_I4P
   if (self%total_coils_number <= 0_I4P .or. (.not.allocated(self%J_vec))) then
      allocate(cells(1:4,1:0))
      return
   endif
   n = 0_I4P
   do sweep=1, 2
      if (sweep == 2) allocate(cells(1:4,1:n))
      n = 0_I4P
      do part=1, 2
         do b=1, blocks_number
            do k=1-self%ngc, self%nk+self%ngc
               do j=1-self%ngc, self%nj+self%ngc
                  do i=1-self%ngc, self%ni+self%ngc
                     is_interior = i >= 1 .and. i <= self%ni .and. j >= 1 .and. j <= self%nj .and. &
                                   k >= 1 .and. k <= self%nk
                     if (is_interior .neqv. (part == 1)) cycle
                     if (any(self%J_vec(:,i,j,k,b,:) /= 0._R8P)) then
                        n = n + 1_I4P
                        if (sweep == 2) cells(:,n) = [b, i, j, k]
                     endif
                  enddo
               enddo
            enddo
         enddo
         if (part == 1) interior = n
      enddo
   enddo
   endsubroutine support_cells
endmodule adam_prism_coil_object
