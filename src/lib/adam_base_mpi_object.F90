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
      procedure, pass(self) :: compute_blocks_number !< Compute maximum blocks number available.
      procedure, pass(self) :: initialize            !< Initialize MPI data.
endtype base_mpi_object

contains
   subroutine compute_blocks_number(self, memory_avail, block_weight, fields_number, nb, nodes_number, verbose)
   !< Compute maximum blocks number available.
   class(base_mpi_object), intent(inout)        :: self          !< The MPI object.
   real(R8P),              intent(in)           :: memory_avail  !< Memory available for single MPI process.
   integer(I4P),           intent(in)           :: block_weight  !< Single block weight.
   integer(I4P),           intent(in)           :: fields_number !< Fields number.
   integer(I4P),           intent(out)          :: nb            !< Maximum blocks number for single MPI process.
   integer(I8P),           intent(out)          :: nodes_number  !< Maximum blocks number for all MPI processes (nodes).
   logical,                intent(in), optional :: verbose       !< Flag to activate verbose output.
   integer(I4P)                                 :: size_of_real  !< Size of real.
   real(R8P)                                    :: save_factor   !< Factor to avoid GPU completely full.

   size_of_real = storage_size(1._R8P)/8._R8P
   save_factor = 0.6_R8P
   nb = nint(save_factor * memory_avail*1e9 / (fields_number * block_weight * size_of_real))
   nodes_number  = nb * self%procs_number
   if (present(verbose)) then
      if (verbose) then
         print '(A)', self%myrankstr//'blocks number for single MPI [nb]: '//trim(str(nb))
         print '(A)', self%myrankstr//'blocks number for all MPI [nodes_number]: '//trim(str(nodes_number))
      endif
   endif
   endsubroutine compute_blocks_number

   subroutine initialize(self, do_mpi_init)
   !< Initialize MPI data.
   class(base_mpi_object), intent(inout)        :: self        !< The MPI object.
   logical,                intent(in), optional :: do_mpi_init !< Flag to activate MPI init call.

   if (present(do_mpi_init)) then
      if (do_mpi_init) call MPI_INIT(self%error)
   endif
   call MPI_COMM_SIZE(MPI_COMM_WORLD, self%procs_number, self%error)
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, self%error)
   self%myrankstr = '[myrank-'//trim(strz(self%myrank,6))//']'
   endsubroutine initialize
endmodule adam_base_mpi_object
