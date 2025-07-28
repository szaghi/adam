!< ADAM, PATCH Boundary Conditions class definition, common CPU backend.
module adam_patch_bc_object
!< ADAM, PATCH Boundary Conditions class definition, CPU backend.

! ADAM modules
use adam_mpih_object
! third party modules
use finer
use penf

implicit none
private
public :: patch_bc_object
public :: BC_DIRICHLET

character(len=8), parameter :: INI_SECTION_NAMES(6)=["bc_x_min", "bc_x_max", &
                                                     "bc_y_min", "bc_y_max", &
                                                     "bc_z_min", "bc_z_max"] !< INI (config) file section name containing BC configs.

integer(I4P), parameter :: BC_DIRICHLET     = 1_I4P              !< Dirichlet, phi = c (0).
character(9), parameter :: BC_DIRICHLET_STR = 'dirichlet'        !< Dirichlet, string input.
character(9), parameter :: BC_TYPE_STR(1)   = [BC_DIRICHLET_STR] !< BC types list, string cast.

type :: patch_bc_object
   !< Boundary Conditions class definition, CPU backend.
   type(mpih_object)      :: mpih       !< MPI handler.
   integer(I4P)           :: bc_type(6) !< Boundary condition type.
   real(R8P), allocatable :: q(:,:)     !< Poisson field variable at BC.
   contains
      ! public methods
      procedure, pass(self) :: description    !< Return pretty-printed object description.
      procedure, pass(self) :: initialize     !< Initialize BC.
      procedure, pass(self) :: load_from_file !< Load config from file.
endtype patch_bc_object

contains
   ! public methods
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(patch_bc_object), intent(in) :: self             !< BC.
   character(len=:), allocatable      :: desc             !< Description.
   character(len=1), parameter        :: NL=new_line('a') !< New line character.
   integer(I4P)                       :: f, v             !< Counter.

   desc =       self%mpih%myrankstr//'BC main data'
   do f=1, 6
      desc = desc//NL//self%mpih%myrankstr//'  face "'//trim(str(f,.true.))//'" BC type "'//trim(BC_TYPE_STR(self%bc_type(f)))//'"'
      do v=1, size(self%q, dim=1)
         desc = desc//NL//self%mpih%myrankstr//'    q('//trim(str(v,.true.))//'): '//trim(str(self%q(v,f)))
      enddo
   enddo
   endfunction description

   subroutine initialize(self, file_parameters)
   !< Initialize the equation.
   class(patch_bc_object), intent(inout) :: self            !< BC.
   type(file_ini),         intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   call self%mpih%initialize
   print '(A)', self%mpih%myrankstr//'patch_bc_object%initialize start'
   allocate(self%q(6,6)) ; self%q = 0._R8P
   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%description()
   print '(A)', self%mpih%myrankstr//'patch_bc_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(patch_bc_object), intent(inout)        :: self            !< BC.
   type(file_ini),         intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,                intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                      :: go_on_fail_     !< Go on if load fails.
   character(8)                                 :: sname           !< Section name.
   character(99)                                :: buff_c          !< Character buffer.
   integer(I4P)                                 :: f               !< Counter.
   integer(I4P)                                 :: error           !< Error status.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail
   do f=1, 6
      sname = INI_SECTION_NAMES(f)
      call file_parameters%get(section_name=sname, option_name='type', val=buff_c, error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(type)')
      select case(trim(adjustl(buff_c)))
      case(BC_DIRICHLET_STR)
         self%bc_type(f) = BC_DIRICHLET
         call file_parameters%get(section_name=sname, option_name='phi', val=self%q(1,f), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(phi)')
      endselect
   enddo
   endsubroutine load_from_file
endmodule adam_patch_bc_object
