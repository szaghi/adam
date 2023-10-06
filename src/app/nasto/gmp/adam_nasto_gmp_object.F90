!< ADAM, Navier-Stokes equations system class definition, GPU-OpenMP (GMP) backend.
module adam_nasto_gmp_object
!< ADAM, Navier-Stokes equations system class definition, GPU-OpenMP (GMP) backend.

! use adam_base_gmp_object
use adam_nasto_common_object
use penf
use MPI

implicit none
private
public :: nasto_gmp_object

type, extends(nasto_common_object) :: nasto_gmp_object
   !< Navier-Stokes equations system class definition, GPU-OpenMP (GMP) backend.
   ! ADAM library objects
   ! type(base_gmp_object) :: base_gmp !< The base GMP handler.
   contains
      ! auxiliary methods
      procedure, pass(self) :: initialize !< Initialize the equation.
      ! numerical methods
      procedure, pass(self) :: simulate   !< Perform the simulation.
endtype nasto_gmp_object
contains
   subroutine initialize(self, filename)
   !< Initialize the equation.
   class(nasto_gmp_object), intent(inout) :: self     !< The equation.
   character(*),            intent(in)    :: filename !< Input file name.
   endsubroutine initialize

   subroutine simulate(self, filename)
   !< Perform the simulation.
   class(nasto_gmp_object), intent(inout) :: self     !< The equation.
   character(*),            intent(in)    :: filename !< Input file name.
   endsubroutine simulate
endmodule adam_nasto_gmp_object
