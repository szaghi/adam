!< ADAM, Boundary Conditions class definition, CPU backend.
module adam_bc_cpu_object
!< ADAM, Boundary Conditions class definition, CPU backend.

use adam_grid_object, only : grid_object
use adam_mpih_object, only : mpih_object
use finer
use penf

implicit none
private
public :: bc_cpu_object
public :: BC_EXTRAPOLATION
public :: BC_INFLOW

character(len=8), parameter :: INI_SECTION_NAMES(6)=["bc_x_min", "bc_x_max", &
                                                     "bc_y_min", "bc_y_max", &
                                                     "bc_z_min", "bc_z_max"] !< INI (config) file section name containing BC configs.

integer(I4P), parameter :: BC_EXTRAPOLATION   = 1_I4P !< Extrapolation.
integer(I4P), parameter :: BC_INFLOW          = 2_I4P !< Supersonic inflow.

type :: bc_cpu_object
   !< Boundary Conditions class definition, CPU backend.
   type(mpih_object)           :: mpih         !< MPI handler.
   type(grid_object),  pointer :: grid=>null() !< The grid.
   real(R8P)                   :: q(6,6)       !< Primitive variables (r,u,v,w,p,s with s being the speccie index) at BC.
   contains
      ! public methods
      procedure, pass(self) :: initialize     !< Initialize BC.
      procedure, pass(self) :: load_from_file !< Load config from file.
endtype bc_cpu_object

contains
   ! public methods
   subroutine initialize(self, file_parameters, grid)
   !< Initialize the equation.
   class(bc_cpu_object), intent(inout)      :: self            !< BC.
   type(file_ini),       intent(in)         :: file_parameters !< Simulation parameters ini file handler.
   type(grid_object),    intent(in), target :: grid            !< The grid.

   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'bc_cpu_object%initialize start'
   self%grid => grid
   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%mpih%myrankstr//'bc_cpu_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(bc_cpu_object), intent(inout)        :: self            !< BC.
   type(file_ini),       intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,              intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                    :: go_on_fail_     !< Go on if load fails.
   character(8)                               :: sname           !< Section name.
   integer(I4P)                               :: bc_type         !< BC type.
   integer(I4P)                               :: b               !< Counter.
   integer(I4P)                               :: error           !< Error status.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail
   do b=1,6
      sname = INI_SECTION_NAMES(b)
      call file_parameters%get(section_name=sname, option_name='type', val=bc_type, error=error)
      ! if (.not.go_on_fail_.and.error>0) error stop self%mpih%myrankstr//'error: failed to load ['//sname//'].(type)'
      select case(bc_type)
      case(BC_INFLOW)
         call file_parameters%get(section_name=sname, option_name='r', val=self%q(1,b), error=error)
         ! if (.not.go_on_fail_.and.error>0) error stop self%mpih%myrankstr//'error: failed to load ['//sname//'].(r)'
         call file_parameters%get(section_name=sname, option_name='u', val=self%q(2,b), error=error)
         ! if (.not.go_on_fail_.and.error>0) error stop self%mpih%myrankstr//'error: failed to load ['//sname//'].(u)'
         call file_parameters%get(section_name=sname, option_name='v', val=self%q(3,b), error=error)
         ! if (.not.go_on_fail_.and.error>0) error stop self%mpih%myrankstr//'error: failed to load ['//sname//'].(v)'
         call file_parameters%get(section_name=sname, option_name='w', val=self%q(4,b), error=error)
         ! if (.not.go_on_fail_.and.error>0) error stop self%mpih%myrankstr//'error: failed to load ['//sname//'].(w)'
         call file_parameters%get(section_name=sname, option_name='p', val=self%q(5,b), error=error)
         ! if (.not.go_on_fail_.and.error>0) error stop self%mpih%myrankstr//'error: failed to load ['//sname//'].(p)'
         call file_parameters%get(section_name=sname, option_name='s', val=self%q(6,b), error=error)
         ! if (.not.go_on_fail_.and.error>0) error stop self%mpih%myrankstr//'error: failed to load ['//sname//'].(s)'
      endselect
   enddo
   endsubroutine load_from_file
endmodule adam_bc_cpu_object
