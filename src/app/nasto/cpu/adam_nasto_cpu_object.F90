!< ADAM, Navier-Stokes equations system class definition, CPU backend.
module adam_nasto_cpu_object
!< ADAM, Navier-Stokes equations system class definition, CPU backend.

! use adam_base_cpu_object
use adam_nasto_common_object
use penf
use MPI

implicit none
private
public :: nasto_cpu_object

type, extends(nasto_common_object) :: nasto_cpu_object
   !< Navier-Stokes equations system class definition, CPU backend.
   ! ADAM library objects
   ! type(base_cpu_object) :: base_cpu !< The base CPU handler.
   contains
      ! auxiliary methods
      procedure, pass(self) :: initialize !< Initialize the equation.
      ! numerical methods
      procedure, pass(self) :: simulate   !< Perform the simulation.
endtype nasto_cpu_object
contains
   subroutine initialize(self, filename)
   !< Initialize the equation.
   class(nasto_cpu_object), intent(inout) :: self     !< The equation.
   character(*),            intent(in)    :: filename !< Input file name.
   endsubroutine initialize

   subroutine simulate(self, filename)
   !< Perform the simulation.
   class(nasto_cpu_object), intent(inout) :: self     !< The equation.
   character(*),            intent(in)    :: filename !< Input file name.
   endsubroutine simulate
endmodule adam_nasto_cpu_object
