!< ADAM, MPI handling library.
module adam_mpi_lib
!< ADAM, MPI handling library.

use PENF
use MPI

implicit none
save
private
public :: myrank
public :: myrankstr
public :: procs_number
public :: error
public :: mpi_initialize

integer(I4P)              :: myrank=0_I4P       !< MPI rank process.
character(:), allocatable :: myrankstr          !< MPI rank process stringified.
integer(I4P)              :: procs_number=1_I4P !< Number of MPI processes.
integer(I4P)              :: error=0_I4P        !< Error traping flag.

contains
   subroutine mpi_initialize(do_mpi_init)
   !< Initialize MPI data.
   logical, intent(in), optional :: do_mpi_init !< Flag to activate MPI init call.

   if (present(do_mpi_init)) then
      if (do_mpi_init) call MPI_INIT(error)
   endif
   call MPI_COMM_SIZE(MPI_COMM_WORLD, procs_number, error)
   call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, error)
   myrankstr = '[myrank-'//trim(strz(myrank,6))//']'
   endsubroutine mpi_initialize
endmodule adam_mpi_lib
