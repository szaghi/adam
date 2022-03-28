!< ADAM, MPI handler object.
module adam_mpih_object
!< ADAM, MPI handler object.

use PENF
use MPI

implicit none
save
private
public :: mpih_object

type :: mpih_object
   integer(I4P)              :: myrank=0_I4P       !< MPI rank process.
   character(:), allocatable :: myrankstr          !< MPI rank process stringified.
   integer(I4P)              :: procs_number=1_I4P !< Number of MPI processes.
   integer(I4P)              :: error=0_I4P        !< Error traping flag.
   contains
      ! public methods
      procedure, pass(self) :: initialize !< Initialize MPI handler data.
endtype mpih_object

contains
   subroutine initialize(self, do_mpi_init)
   !< Initialize MPI handler data.
   class(mpih_object) , intent(inout)        :: self        !< MPI handler.
   logical,             intent(in), optional :: do_mpi_init !< Flag to activate MPI init call.

   if (present(do_mpi_init)) then
      if (do_mpi_init) call MPI_INIT(self%error)
   endif
   call MPI_COMM_SIZE(MPI_COMM_WORLD, self%procs_number, self%error)
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, self%error)
   self%myrankstr = '[myrank-'//trim(strz(self%myrank,6))//']'
   endsubroutine initialize
endmodule adam_mpih_object
