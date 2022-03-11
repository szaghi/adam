!< ADAM, MPI class definition.
module adam_base_mpi_object
!< ADAM, MPI class definition: minimal object to handle MPI communications.

use PENF
use MPI
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

implicit none
private
public :: base_mpi_object

type :: base_mpi_object
   !< Base MPI class definition.
   !<
   !< Provide methods for MPI communications handling.
   integer(I4P)              :: myrank=0_I4P       !< MPI rank process.
   character(:), allocatable :: myrankstr          !< MPI rank process stringified.
   integer(I4P)              :: procs_number=1_I4P !< Number of MPI processes.
   integer(I4P)              :: error=0_I4P        !< Error traping flag.
   contains
      ! public methods
      procedure, pass(self) :: initialize !< Initialize MPI data.
endtype base_mpi_object

contains
   subroutine initialize(self)
   !< Initialize MPI data.
   class(base_mpi_object), intent(inout) :: self !< The MPI object.

   call MPI_COMM_SIZE(MPI_COMM_WORLD, self%procs_number, self%error)
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, self%error)
   self%myrankstr = '[myrank-'//trim(strz(self%myrank,6))//']'
   endsubroutine initialize
endmodule adam_base_mpi_object
