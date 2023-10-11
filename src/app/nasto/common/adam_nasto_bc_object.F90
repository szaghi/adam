!< ADAM, NASTO Boundary Conditions class definition, CPU backend.
module adam_nasto_bc_object
!< ADAM, NASTO Boundary Conditions class definition, CPU backend.

use adam_mpih_object
use finer
use penf

implicit none
private
public :: nasto_bc_object
public :: BC_EXTRAPOLATION
public :: BC_INFLOW
public :: BC_WALL_INVISCID

character(len=8), parameter :: INI_SECTION_NAMES(6)=["bc_x_min", "bc_x_max", &
                                                     "bc_y_min", "bc_y_max", &
                                                     "bc_z_min", "bc_z_max"] !< INI (config) file section name containing BC configs.

integer(I4P), parameter :: BC_EXTRAPOLATION   = 1_I4P !< Extrapolation.
integer(I4P), parameter :: BC_INFLOW          = 2_I4P !< Supersonic inflow.
integer(I4P), parameter :: BC_WALL_INVISCID   = 3_I4P !< Inviscid wall.

type :: nasto_bc_object
   !< Boundary Conditions class definition, CPU backend.
   type(mpih_object)           :: mpih         !< MPI handler.
   integer(I4P)                :: bc_type(6)   !< Boundary condition type.
   real(R8P)                   :: q(6,6)       !< Primitive variables (r,u,v,w,p,s with s being the specie index) at BC.
   contains
      ! public methods
      procedure, pass(self) :: initialize     !< Initialize BC.
      procedure, pass(self) :: load_from_file !< Load config from file.
endtype nasto_bc_object

contains
   ! public methods
   subroutine initialize(self, file_parameters)
   !< Initialize the equation.
   class(nasto_bc_object), intent(inout) :: self            !< BC.
   type(file_ini),         intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'nasto_bc_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%mpih%myrankstr//'nasto_bc_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(nasto_bc_object), intent(inout)        :: self            !< BC.
   type(file_ini),         intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,                intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                      :: go_on_fail_     !< Go on if load fails.
   character(8)                                 :: sname           !< Section name.
   character(99)                                :: buff_c          !< Character buffer.
   integer(I4P)                                 :: b               !< Counter.
   integer(I4P)                                 :: error           !< Error status.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail
   do b=1, 6
      sname = INI_SECTION_NAMES(b)
      call file_parameters%get(section_name=sname, option_name='type', val=buff_c, error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(type)')
      select case(trim(adjustl(buff_c)))
      case('extrapolation')
         self%bc_type(b) = BC_EXTRAPOLATION
      case('inflow')
         self%bc_type(b) = BC_INFLOW
         call file_parameters%get(section_name=sname, option_name='r', val=self%q(1,b), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(r)')
         call file_parameters%get(section_name=sname, option_name='u', val=self%q(2,b), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(u)')
         call file_parameters%get(section_name=sname, option_name='v', val=self%q(3,b), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(v)')
         call file_parameters%get(section_name=sname, option_name='w', val=self%q(4,b), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(w)')
         call file_parameters%get(section_name=sname, option_name='p', val=self%q(5,b), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(p)')
         call file_parameters%get(section_name=sname, option_name='s', val=self%q(6,b), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(s)')
      case('wall-inviscid')
         self%bc_type(b) = BC_WALL_INVISCID
      endselect
   enddo
   endsubroutine load_from_file
endmodule adam_nasto_bc_object
