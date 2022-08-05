!< ADAM, Physics class definition, CPU backend.
module adam_physics_cpu_object
!< ADAM, Physics class definition, CPU backend.

use adam_eos_ic_cpu_object, only : eos_ic_cpu_object
use adam_mpih_object,       only : mpih_object
use finer
use penf

implicit none
private
public :: physics_cpu_object
public :: IR
public :: IU
public :: IV
public :: IW
public :: IG
public :: IP

character(len=7), parameter :: INI_SECTION_NAME="physics" !< INI (config) file section name containing physics configs.

! Named indexes of q_aux variables, to be initialized when ns is known.
integer(I4P) :: IR = 2_I4P
integer(I4P) :: IU = 3_I4P
integer(I4P) :: IV = 4_I4P
integer(I4P) :: IW = 5_I4P
integer(I4P) :: IG = 6_I4P
integer(I4P) :: IP = 7_I4P

type :: physics_cpu_object
   !< Physics class definition, CPU backend.
   type(mpih_object)                    :: mpih         !< MPI handler.
   integer(I4P)                         :: ns=1_I4P     !< Number of species.
   integer(I4P)                         :: nv=5_I4P     !< Number of variables.
   integer(I4P)                         :: nv_aux=7_I4P !< Number of auxiliary variables.
   type(eos_ic_cpu_object), allocatable :: eos(:)       !< Equations of state (cp, cv...) of each specie [1:ns].
   contains
      ! public methods
      procedure, pass(self) :: description    !< Return pretty-printed object description.
      procedure, pass(self) :: initialize     !< Initialize physics.
      procedure, pass(self) :: load_from_file !< Load config from file.
endtype physics_cpu_object

contains
   ! public methods
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(physics_cpu_object), intent(in) :: self             !< Physics.
   character(len=:), allocatable         :: desc             !< Description.
   character(len=1), parameter           :: NL=new_line('a') !< New line character.
   integer(I4P)                          :: s                !< Counter.

   desc =       self%mpih%myrankstr//'physics_cpu_object ns:     '//trim(str(self%ns    ))//NL
   desc = desc//self%mpih%myrankstr//'physics_cpu_object nv:     '//trim(str(self%nv    ))//NL
   desc = desc//self%mpih%myrankstr//'physics_cpu_object nv_aux: '//trim(str(self%nv_aux))//NL
   desc = desc//self%mpih%myrankstr//'physics_cpu_object IR:     '//trim(str(IR         ))//NL
   desc = desc//self%mpih%myrankstr//'physics_cpu_object IU:     '//trim(str(IU         ))//NL
   desc = desc//self%mpih%myrankstr//'physics_cpu_object IV:     '//trim(str(IV         ))//NL
   desc = desc//self%mpih%myrankstr//'physics_cpu_object IW:     '//trim(str(IW         ))//NL
   desc = desc//self%mpih%myrankstr//'physics_cpu_object IG:     '//trim(str(IG         ))//NL
   desc = desc//self%mpih%myrankstr//'physics_cpu_object IP:     '//trim(str(IP         ))
   do s=1, self%ns
      desc = desc//NL//self%mpih%myrankstr//'physics_cpu_object cp('//trim(str(s,.true.))//'):  '//trim(str(self%eos(s)%cp))
      desc = desc//NL//self%mpih%myrankstr//'physics_cpu_object cv('//trim(str(s,.true.))//'):  '//trim(str(self%eos(s)%cv))
      desc = desc//NL//self%mpih%myrankstr//'physics_cpu_object  g('//trim(str(s,.true.))//'):  '//trim(str(self%eos(s)%g ))
      desc = desc//NL//self%mpih%myrankstr//'physics_cpu_object  R('//trim(str(s,.true.))//'):  '//trim(str(self%eos(s)%R ))
   enddo
   endfunction description

   subroutine initialize(self, file_parameters)
   !< Initialize the equation.
   class(physics_cpu_object), intent(inout) :: self            !< Physics.
   type(file_ini),            intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'physics_cpu_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   self%nv     = self%ns - 1 + 5
   self%nv_aux = self%ns     + 6
   ! initialize named index of q_aux array
   IR = self%ns + 1
   IU = self%ns + 2
   IV = self%ns + 3
   IW = self%ns + 4
   IG = self%ns + 5
   IP = self%ns + 6
   print '(A)', self%description()
   print '(A)', self%mpih%myrankstr//'physics_cpu_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(physics_cpu_object), intent(inout)        :: self            !< Physics.
   type(file_ini),            intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,                   intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                         :: go_on_fail_     !< Go on if load fails.
   character(:), allocatable                       :: sname           !< Section name.
   real(R8P)                                       :: cp, cv          !< Constant specific heats.
   integer(I4P)                                    :: s               !< Counter.
   integer(I4P)                                    :: error           !< Error status.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='ns', val=self%ns, error=error)
   ! if (.not.go_on_fail_.and.error>0) error stop self%mpih%myrankstr//'error: failed to load ['//INI_SECTION_NAME//'].(ns)'
   allocate(self%eos(1:self%ns))
   do s=1, self%ns
      sname = INI_SECTION_NAME//'_specie_'//trim(str(s,.true.))
      call file_parameters%get(section_name=sname, option_name='cp', val=cp, error=error)
      ! if (.not.go_on_fail_.and.error>0) error stop self%mpih%myrankstr//'error: failed to load ['//sname//'].(cp)'
      call file_parameters%get(section_name=sname, option_name='cv', val=cv, error=error)
      ! if (.not.go_on_fail_.and.error>0) error stop self%mpih%myrankstr//'error: failed to load ['//sname//'].(cv)'
      call self%eos(s)%initialize(cp=cp, cv=cv)
   enddo
   endsubroutine load_from_file
endmodule adam_physics_cpu_object
