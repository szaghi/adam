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

character(len=8), parameter :: INI_SECTION_NAMES(6)=["bc_x_min", "bc_x_max", &
                                                     "bc_y_min", "bc_y_max", &
                                                     "bc_z_min", "bc_z_max"] !< INI (config) file section name containing BC configs.

type :: bc_cpu_object
   !< Boundary Conditions class definition, CPU backend.
   type(mpih_object)           :: mpih         !< MPI handler.
   type(grid_object),  pointer :: grid=>null() !< The grid.
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
   character(99)                              :: char_buff       !< Character buffer.
   integer(I4P)                               :: error           !< Error status.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   do i_bc=1,6
      sname = INI_SECTION_NAMES(i_bc)
      call file_parameters%get(section_name=sname, option_name='type', val=char_buff, error=error)
      call self%file_input%get(section_name=sname, option_name='type', val=bc_type_item)
      n_vars = BC_VARS_NUMBER(bc_type_item)
      do i_var=1,n_vars
          call self%file_input%get(section_name=sname, option_name="var"//trim(str(i_var,.true.)), val=buf_R8)
          self%bc_vars(i_var, i_bc) = buf_R8
      enddo
   enddo

   if (.not.go_on_fail_.and.error>0) error stop self%mpih%myrankstr//'error: failed to load ['//INI_SECTION_NAME//'].(runge_kutta)'

   endsubroutine load_from_file
endmodule adam_bc_cpu_object
