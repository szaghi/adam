!< ADAM, FLUME numerics class definition.
module adam_flume_numerics_object
!< ADAM, FLUME numerics class definition.
!<
!< The time scheme is not duplicated here: it is the library `[runge_kutta].(scheme)`.
!<
!< `reflux` switches the Berger-Colella correction at AMR coarse-fine faces: `.true.` for every production run;
!< `.false.` is a diagnostic (the negative control of the conservation test, issue #35 V3), since a run on a grid
!< with coarse-fine faces is then not conservative.

! ADAM singleton objects
use :: adam_mpih_global,      only : mpih
! FLUME modules
use :: adam_flume_parameters, only : RECON_CHARACTERISTIC, RECON_CONSERVATIVE, SCHEME_SPACE_WENO, strip_control
! third party modules
use :: finer,                 only : file_ini
use :: penf,                  only : I4P

implicit none
private
public :: flume_numerics_object

character(len=8), parameter :: INI_SECTION_NAME="numerics" !< INI (config) file section name containing numerics configs.

type :: flume_numerics_object
   !< FLUME numerics class definition.
   character(:), allocatable :: scheme_space             !< Spatial scheme.
   character(:), allocatable :: reconstruction_variables !< Variables reconstructed at cell interfaces.
   logical                   :: reflux=.true.            !< Berger-Colella reflux at AMR coarse-fine faces.
   contains
      ! public methods
      procedure, pass(self) :: description    !< Return pretty-printed object description.
      procedure, pass(self) :: initialize     !< Initialize numerics.
      procedure, pass(self) :: load_from_file !< Load config from file.
endtype flume_numerics_object

contains
   ! public methods
   function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(flume_numerics_object), intent(in) :: self             !< Numerics.
   character(len=:), allocatable            :: desc             !< Description.
   character(len=1), parameter              :: NL=new_line('a') !< New line character.

   desc =       mpih%myrankstr//'Numerics main data'//NL
   desc = desc//mpih%myrankstr//'  scheme_space:             '//self%scheme_space//NL
   desc = desc//mpih%myrankstr//'  reconstruction_variables: '//self%reconstruction_variables//NL
   desc = desc//mpih%myrankstr//'  reflux:                   '//trim(merge('.true. ', '.false.', self%reflux))
   endfunction description

   subroutine initialize(self, file_parameters)
   !< Initialize numerics.
   class(flume_numerics_object), intent(inout) :: self            !< Numerics.
   type(file_ini),               intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   print '(A)', mpih%myrankstr//'flume_numerics_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%description()
   print '(A)', mpih%myrankstr//'flume_numerics_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters)
   !< Load config from file; every key is required and an unknown value is fatal.
   class(flume_numerics_object), intent(inout) :: self            !< Numerics.
   type(file_ini),               intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   character(999)                              :: buff            !< Option value buffer.
   integer(I4P)                                :: error           !< Error status.

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='scheme_space', val=buff, error=error)
   if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(scheme_space)')
   self%scheme_space = trim(adjustl(strip_control(buff)))
   select case(self%scheme_space)
   case(SCHEME_SPACE_WENO)
   case default
      call mpih%error_stop(msg=': unknown ['//INI_SECTION_NAME//'].(scheme_space) "'//self%scheme_space// &
                               '"; expected one of '//SCHEME_SPACE_WENO)
   endselect

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='reconstruction_variables', val=buff, error=error)
   if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(reconstruction_variables)')
   self%reconstruction_variables = trim(adjustl(strip_control(buff)))
   select case(self%reconstruction_variables)
   case(RECON_CHARACTERISTIC, RECON_CONSERVATIVE)
   case default
      call mpih%error_stop(msg=': unknown ['//INI_SECTION_NAME//'].(reconstruction_variables) "'// &
                               self%reconstruction_variables//'"; expected one of '//RECON_CHARACTERISTIC//', '// &
                               RECON_CONSERVATIVE)
   endselect

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='reflux', val=self%reflux, error=error)
   if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(reflux)')
   endsubroutine load_from_file
endmodule adam_flume_numerics_object
