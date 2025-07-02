module adam_prism_physics_object

use adam_mpih_object, only : mpih_object
!use adam_nasto_eos_object, only : nasto_eos_object
use finer, only : file_ini
use penf, only : I4P, R8P, str

implicit none
private
public :: prism_physics_object

character(len=7), parameter :: INI_SECTION_NAME='physics' !< INI file section name containing fluid physics.

integer(I4P), parameter :: VAR_DX = 1_I4P
integer(I4P), parameter :: VAR_DY = 2_I4P
integer(I4P), parameter :: VAR_DZ = 3_I4P
integer(I4P), parameter :: VAR_BX = 4_I4P
integer(I4P), parameter :: VAR_BY = 5_I4P
integer(I4P), parameter :: VAR_BZ = 6_I4P
integer(I4P), parameter :: VAR_JX = 7_I4P
integer(I4P), parameter :: VAR_JY = 8_I4P
integer(I4P), parameter :: VAR_JZ = 9_I4P

type :: prism_physics_object
   !< PRISM physics class definition.
   type(mpih_object)                   :: mpih                          !< MPI handler.
   integer(I4P)                        :: nv=0_I4P                      !< Number of variables (see below).
   logical                             :: D_divergence_cleaner=.false.  !< Enable electric field divergence cleaning.
   logical                             :: B_divergence_cleaner=.false.  !< Enable magnetic field divergence cleaning.
   real(R8P)                           :: chi                           !< Coefficiente velocità
                                                                        !< trasporto errore divergenza campi
   real(R8P)                           :: eta                           !< Coefficiente correzione parabolica
                                                                        !< trasporto errore divergenza campi

   !integer(I4P)                        :: ns=1_I4P     !< Number of species.
   !integer(I4P)                        :: nv_aux=9_I4P !< Number of auxiliary variables (rns+r+u+v+w+p+g=ns+6).
   !integer(I4P)                        :: np=5_I4P     !< Number of 1D primitive variables (rns+r+un+p+g=ns+4).
   !type(nasto_eos_object), allocatable :: eos(:)       !< Equations of state of each specie [1:ns].

   contains
      ! public methods
      !procedure, pass(self) :: conservative2primitive !< Return primitive variables from conservative ones.
      procedure, pass(self) :: description            !< Return pretty-printed object description.
      procedure, pass(self) :: initialize             !< Initialize physics.
      procedure, pass(self) :: load_from_file         !< Load config from file.
      !procedure, pass(self) :: primitive2conservative !< Return conservative variables from primitive ones.
endtype prism_physics_object

contains

pure function description(self) result(desc)
!< Return a pretty-formatted object description.
class(prism_physics_object), intent(in) :: self             !< Physics.
character(len=:), allocatable           :: desc             !< Description.
character(len=1), parameter             :: NL=new_line('a') !< New line character.
!integer(I4P)                            :: s                !< Counter.

desc =       self%mpih%myrankstr//'Physics main data:'                                     //NL
desc = desc//self%mpih%myrankstr//'  prism_physics_object nv:     '//trim(str(self%nv    ))//NL
desc = desc//self%mpih%myrankstr//'  D divergence correction:     '//trim(str(self%D_divergence_cleaner    ))//NL
desc = desc//self%mpih%myrankstr//'  B divergence correction:     '//trim(str(self%B_divergence_cleaner    ))//NL
!desc = desc//self%mpih%myrankstr//'  prism_physics_object VAR_DX: '//trim(str(VAR_DX     ))//NL
!desc = desc//self%mpih%myrankstr//'  prism_physics_object VAR_DY: '//trim(str(VAR_DY     ))//NL
!desc = desc//self%mpih%myrankstr//'  prism_physics_object VAR_DZ: '//trim(str(VAR_DZ     ))//NL
!desc = desc//self%mpih%myrankstr//'  prism_physics_object VAR_BX: '//trim(str(VAR_BX     ))//NL
!desc = desc//self%mpih%myrankstr//'  prism_physics_object VAR_BY: '//trim(str(VAR_BY     ))//NL
!desc = desc//self%mpih%myrankstr//'  prism_physics_object VAR_BZ: '//trim(str(VAR_BZ     ))//NL
!desc = desc//self%mpih%myrankstr//'  prism_physics_object VAR_JX: '//trim(str(VAR_JX     ))//NL
!desc = desc//self%mpih%myrankstr//'  prism_physics_object VAR_JY: '//trim(str(VAR_JY     ))//NL
!desc = desc//self%mpih%myrankstr//'  prism_physics_object VAR_JZ: '//trim(str(VAR_JZ     ))//NL

endfunction description

subroutine initialize(self, file_parameters)
!< Initialize the equation.
class(prism_physics_object), intent(inout) :: self            !< Physics.
type(file_ini),              intent(in)    :: file_parameters !< Simulation parameters ini file handler.

call self%mpih%initialize(do_mpi_init=.false.)
print '(A)', self%mpih%myrankstr//'prism_physics_object%initialize start'
call self%load_from_file(file_parameters=file_parameters)
!self%nv     = self%nv
! initialize named index of q_aux array
print '(A)', self%description()
print '(A)', self%mpih%myrankstr//'prism_physics_object%initialize finish'
endsubroutine initialize

subroutine load_from_file(self, file_parameters, go_on_fail)
!< Load config from file.
class(prism_physics_object), intent(inout)        :: self            !< Physics.
type(file_ini),              intent(in)           :: file_parameters !< Simulation parameters ini file handler.
logical,                     intent(in), optional :: go_on_fail      !< Go on if load fails.
logical                                           :: go_on_fail_     !< Go on if load fails.
integer(I4P)                                      :: s               !< Counter.
integer(I4P)                                      :: error           !< Error status.

go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

call file_parameters%get(section_name=INI_SECTION_NAME, option_name='nv', val=self%nv, error=error)
if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(nv)')

call file_parameters%get(section_name=INI_SECTION_NAME, option_name='D_divergence_cleaner', &
                         val=self%D_divergence_cleaner,error=error)
if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load [' &
                                          //INI_SECTION_NAME//'].(D_divergence_cleaner)')

call file_parameters%get(section_name=INI_SECTION_NAME, option_name='B_divergence_cleaner', &
                         val=self%B_divergence_cleaner,error=error)
if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load [' &
                                          //INI_SECTION_NAME//'].(B_divergence_cleaner)')

call file_parameters%get(section_name=INI_SECTION_NAME, option_name='chi', val=self%chi, error=error)
if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(chi)')

call file_parameters%get(section_name=INI_SECTION_NAME, option_name='eta', val=self%eta, error=error)
if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(eta)')

if(self%D_divergence_cleaner .and. .not.self%B_divergence_cleaner) then
   self%nv = 10_I4P
   !call self%mpih%error_stop(msg=': D_divergence_cleaner is true but nv /= 10')
elseif(self%D_divergence_cleaner .and. self%B_divergence_cleaner .and. self%nv /= 11_I4P) then
   self%nv = 11_I4P
   !call self%mpih%error_stop(msg=': D and B_divergence_cleaner is true but nv /= 11')
endif

if(.not.self%D_divergence_cleaner .and. .not.self%B_divergence_cleaner .and. self%nv /= 9_I4P) then
   call self%mpih%error_stop(msg=': D and B_divergence_cleaner are false but nv /= 9')
elseif(.not.self%D_divergence_cleaner .and. self%B_divergence_cleaner) then
   call self%mpih%error_stop(msg=': D_divergence_cleaner is false but B_divergence_cleaner is true')
elseif(self%D_divergence_cleaner .and. self%chi <= 1._R8P) then
   call self%mpih%error_stop(msg=': If cleaner is true chi cannot be lower than 1.0')
endif

print *, self%nv
print *, self%chi
print *, self%eta
endsubroutine load_from_file

endmodule adam_prism_physics_object
