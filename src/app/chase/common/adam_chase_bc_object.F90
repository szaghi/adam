!< ADAM, CHASE Boundary Conditions class definition, common CPU backend.
module adam_chase_bc_object
!< ADAM, CHASE Boundary Conditions class definition, CPU backend.

! ADAM modules
use adam_mpih_object
! third party modules
use finer
use penf

implicit none
private
public :: chase_bc_object
public :: BC_EXTRAPOLATION
public :: BC_INFLOW
public :: BC_WALL_INVISCID

character(len=8), parameter :: INI_SECTION_NAMES(6)=["bc_x_min", "bc_x_max", &
                                                     "bc_y_min", "bc_y_max", &
                                                     "bc_z_min", "bc_z_max"] !< INI (config) file section name containing BC configs.

integer(I4P),  parameter :: BC_EXTRAPOLATION     = 1_I4P                  !< Extrapolation.
integer(I4P),  parameter :: BC_INFLOW            = 2_I4P                  !< Supersonic inflow.
integer(I4P),  parameter :: BC_WALL_INVISCID     = 3_I4P                  !< Inviscid wall.
character(13), parameter :: BC_EXTRAPOLATION_STR = 'extrapolation'        !< Extrapolation, string input.
character(13), parameter :: BC_INFLOW_STR        = 'inflow       '        !< Supersonic inflow, string input.
character(13), parameter :: BC_WALL_INVISCID_STR = 'wall-inviscid'        !< Inviscid wall, string input.
character(13), parameter :: BC_TYPE_STR(3)       = [BC_EXTRAPOLATION_STR,&
                                                    BC_INFLOW_STR,       &
                                                    BC_WALL_INVISCID_STR] !< BC types list, string cast.

type :: chase_bc_object
   !< Boundary Conditions class definition, CPU backend.
   type(mpih_object)      :: mpih       !< MPI handler.
   integer(I4P)           :: bc_type(6) !< Boundary condition type.
   real(R8P), allocatable :: q(:,:)     !< Primitive variables (Dx,Dy,Dz,Bx,By,Bz,Jx,Jy,Jz) at BC.
   contains
      ! public methods
      procedure, pass(self) :: description    !< Return pretty-printed object description.
      procedure, pass(self) :: initialize     !< Initialize BC.
      procedure, pass(self) :: load_from_file !< Load config from file.
endtype chase_bc_object

contains
   ! public methods
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(chase_bc_object), intent(in) :: self             !< BC.
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
   class(chase_bc_object), intent(inout) :: self            !< BC.
   type(file_ini),         intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   call self%mpih%initialize
   print '(A)', self%mpih%myrankstr//'chase_bc_object%initialize start'
   allocate(self%q(6,6)) ; self%q = 0._R8P
   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%description()
   print '(A)', self%mpih%myrankstr//'chase_bc_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(chase_bc_object), intent(inout)        :: self            !< BC.
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
      case(BC_EXTRAPOLATION_STR)
         self%bc_type(f) = BC_EXTRAPOLATION
      case(BC_INFLOW_STR)
         self%bc_type(f) = BC_INFLOW
         call file_parameters%get(section_name=sname, option_name='r', val=self%q(1,f), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(r)')
         call file_parameters%get(section_name=sname, option_name='u', val=self%q(2,f), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(u)')
         call file_parameters%get(section_name=sname, option_name='v', val=self%q(3,f), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(v)')
         call file_parameters%get(section_name=sname, option_name='w', val=self%q(4,f), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(w)')
         call file_parameters%get(section_name=sname, option_name='p', val=self%q(5,f), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(p)')
         call file_parameters%get(section_name=sname, option_name='s', val=self%q(6,f), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(s)')
      case(BC_WALL_INVISCID_STR)
         self%bc_type(f) = BC_WALL_INVISCID
      endselect
   enddo
   endsubroutine load_from_file
endmodule adam_chase_bc_object
