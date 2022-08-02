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
   real(R8P)                 :: timing(1:2)        !< Tic toc timing.
   integer(I4P)              :: tictoc=1_I4P       !< Next is tic or toc?
   contains
      ! public methods
      procedure, pass(self) :: barrier       !< Handy MPI barrier wrapper.
      procedure, pass(self) :: initialize    !< Initialize MPI handler data.
      procedure, pass(self) :: tictoc_timing !< Return the last tic toc timing.
endtype mpih_object

contains
   subroutine barrier(self, tictoc, timing, single)
   !< Handy MPI barrier wrapper.
   class(mpih_object) , intent(inout)         :: self    !< MPI handler.
   logical,             intent(in),  optional :: tictoc  !< Activate tic toc timing between 2 barrier calls.
   real(R8P),           intent(out), optional :: timing  !< Current timing.
   logical,             intent(in),  optional :: single  !< Single tictoc for one-shot timing.
   logical                                    :: tictoc_ !< Activate tic toc timing between 2 barrier calls, local var.
   logical                                    :: single_ !< Single tictoc for one-shot timing, local var.

   tictoc_ = .false. ; if (present(tictoc)) tictoc_ = tictoc
   single_ = .false. ; if (present(single)) single_ = single
   call MPI_BARRIER(MPI_COMM_WORLD, self%error)
   if (tictoc_) then
      self%timing(self%tictoc) = MPI_WTIME()
      if (present(timing)) timing = self%timing(self%tictoc)
      if (.not.single_) then
         if (self%tictoc==1) then
            self%tictoc = 2
         else
            self%tictoc = 1
         endif
      endif
   endif
   endsubroutine barrier

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

   function tictoc_timing(self) result(timing)
   !< Return the last tic toc timing.
   class(mpih_object) , intent(in) :: self   !< MPI handler.
   real(R8P)                       :: timing !< Last tic toc timing.

   timing = self%timing(2) - self%timing(1)
   endfunction tictoc_timing
endmodule adam_mpih_object
