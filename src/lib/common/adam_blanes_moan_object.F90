!< ADAM, Blanes-Moan solver class definition.
module adam_blanes_moan_object
!< ADAM, Blanes-Moan solver class definition.

use adam_global_field, only: field
use adam_global_mpih, only: mpih
use adam_global_grid, only: grid
use finer
use penf

implicit none
save
private
public :: blanesmoan_object
public :: BM_SCHEME4
public :: BM_SCHEME6

character(len=12), parameter :: BM_SCHEME4="blanes-moan4" !< blanes-Moan 4th order scheme name.
character(len=12), parameter :: BM_SCHEME6="blanes-moan6" !< blanes-Moan 6th order scheme name.

character(len=11), parameter :: INI_SECTION_NAME="blanes_moan" !< INI (config) file section name containing configs.

type :: blanesmoan_object
   !< Leapforg class definition.
   character(:), allocatable :: scheme    !< Scheme name [blanes-moan4, blanes-moan6].
   integer(I4P)              :: nc=4_I4P  !< Number of coefficients.
   real(R8P),    allocatable :: a(:),b(:) !< Splitted residuals coefficients.
   ! grid data replica for easy handling
   integer(I4P),       pointer :: ngc=>null()           !< Number of ghost cells.
   integer(I4P),       pointer :: ni=>null()            !< Number of cells in i direction.
   integer(I4P),       pointer :: nj=>null()            !< Number of cells in j direction.
   integer(I4P),       pointer :: nk=>null()            !< Number of cells in k direction.
   contains
      ! public methods
      procedure, pass(self) :: description    !< Return pretty-printed object description.
      procedure, pass(self) :: initialize     !< Initialize class.
      procedure, pass(self) :: load_from_file !< Load config from file.
endtype blanesmoan_object
contains
   ! public methods
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(blanesmoan_object), intent(in) :: self             !< Blanes-Moan object.
   character(len=:),       allocatable  :: desc             !< Description.
   character(len=1),       parameter    :: NL=new_line('a') !< New line character.

   desc =       mpih%myrankstr//'Blanes-Moan scheme main data'                  //NL
   desc = desc//mpih%myrankstr//'  number of coefficients: '//trim(str(self%nc))//NL
   desc = desc//mpih%myrankstr//'  a coefficients:         '//trim(str(self%a ))//NL
   desc = desc//mpih%myrankstr//'  b coefficients:         '//trim(str(self%b ))
   endfunction description

   subroutine initialize(self, file_parameters, scheme)
   !< Initialize class.
   class(blanesmoan_object), intent(inout)        :: self            !< Blanes-Moan object.
   type(file_ini),           intent(in), optional :: file_parameters !< Simulation parameters ini file handler.
   character(*),             intent(in), optional :: scheme          !< Scheme name.
   real(R8P)                                      :: theta           !< Temporary coefficient.

   call mpih%print_message('blanesmoan_object%initialize start')
   call associate_adam_data
   if (present(file_parameters)) then
      call self%load_from_file(file_parameters=file_parameters)
   elseif (present(scheme)) then
      self%scheme = trim(adjustl(scheme))
   else
      call mpih%error_stop(msg=': failed to initialize blanesmoan object, file parameters or scheme name must be passed')
   endif
   ! associate(a=>self%a,b=>self%b)
   select case(self%scheme)
   case(BM_SCHEME4)
      self%nc = 4
      call allocate_variable(var=self%a, ulb=[1,self%nc], msg=mpih%myrankstr//'blanesmoan_object%initialize allocate a')
      call allocate_variable(var=self%b, ulb=[1,self%nc], msg=mpih%myrankstr//'blanesmoan_object%initialize allocate b')
      self%a(1) =  0.0792036964311957_R8P
      self%a(2) =  0.3531729060497740_R8P
      self%a(3) = -0.0420650803577195_R8P
      self%a(4) =  1.0_R8P - 2.0_R8P*(self%a(1) + self%a(2) + self%a(3))

      self%b(1) =  0.209515106613362_R8P
      self%b(2) = -0.143851773179818_R8P
      self%b(3) =  0.5_R8P - (self%b(1) + self%b(2))
      self%b(4) =  0.0_R8P
   case(BM_SCHEME6)
      self%nc = 6
      call allocate_variable(var=self%a, ulb=[1,self%nc], msg=mpih%myrankstr//'blanesmoan_object%initialize allocate a')
      call allocate_variable(var=self%b, ulb=[1,self%nc], msg=mpih%myrankstr//'blanesmoan_object%initialize allocate b')
      self%a(1) =  0.0502627644003922_R8P
      self%a(2) =  0.4135143004283440_R8P
      self%a(3) =  0.0450798897943977_R8P
      self%a(4) = -0.1880548538195690_R8P
      self%a(5) =  0.5419606784507800_R8P
      self%a(6) =  1.0_R8P - 2.0_R8P*(self%a(1) + self%a(2) + self%a(3) + self%a(4) + self%a(5))

      self%b(1) =  0.148816447901042_R8P
      self%b(2) = -0.132385865767784_R8P
      self%b(3) =  0.067307604692185_R8P
      self%b(4) =  0.432666402578175_R8P
      self%b(5) =  0.5_R8P - (self%b(1) + self%b(2) + self%b(3) + self%b(4))
      self%b(6) =  0.0_R8P
   endselect
   ! endassociate
   print '(A)', self%description()
   call mpih%print_message('blanesmoan_object%initialize finish')
   contains
      subroutine associate_adam_data
      !< Associate grid data pointers for easy handling.

      self%ni  => grid%ni
      self%nj  => grid%nj
      self%nk  => grid%nk
      self%ngc => grid%ngc
      endsubroutine associate_adam_data
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(blanesmoan_object), intent(inout)        :: self            !< Blanes-Moan object.
   type(file_ini),           intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,                  intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                        :: go_on_fail_     !< Go on if load fails.
   integer(I4P)                                   :: error           !< Error status.
   character(99)                                  :: buff            !< Buffer variable.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='scheme', val=buff, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(scheme)')
   self%scheme = trim(adjustl(buff))
   endsubroutine load_from_file
endmodule adam_blanes_moan_object
