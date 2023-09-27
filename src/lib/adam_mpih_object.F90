!< ADAM, MPI handler object.
module adam_mpih_object
!< ADAM, MPI handler object.

use adam_memory_cpu_lib
use penf
use mpi
use, intrinsic :: iso_c_binding, only : C_LONG
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

implicit none
save
private
public :: mpih_object

type :: mpih_object
   integer(I4P)              :: myrank=0_I4P        !< MPI rank process.
   character(:), allocatable :: myrankstr           !< MPI rank process stringified.
   integer(I4P)              :: procs_number=1_I4P  !< Number of MPI processes.
   real(R8P)                 :: memory_avail=0._R8P !< CPU memory available (GB) for each process.
   integer(I4P)              :: error=0_I4P         !< Error traping flag.
   real(R8P)                 :: timing(1:2)         !< Tic toc timing.
   integer(I4P)              :: tictoc=1_I4P        !< Next is tic or toc?
   contains
      ! public methods
      procedure, pass(self) :: abort         !< Handy MPI abort wrapper.
      procedure, pass(self) :: barrier       !< Handy MPI barrier wrapper.
      procedure, pass(self) :: error_stop    !< Stop run with error output.
      procedure, pass(self) :: finalize      !< Handy MPI finalize wrapper.
      procedure, pass(self) :: initialize    !< Initialize MPI handler data.
      procedure, pass(self) :: tictoc_timing !< Return the last tic toc timing.
endtype mpih_object

contains
   subroutine abort(self, error_code, msg)
   !< Handy MPI abort wrapper.
   class(mpih_object) , intent(inout)        :: self        !< MPI handler.
   integer(I4P),        intent(in), optional :: error_code  !< Abort error code.
   character(*),        intent(in), optional :: msg         !< Error message.
   character(:), allocatable                 :: msg_        !< Error message, local variable.
   integer(I4P)                              :: error_code_ !< Abort error code, local variable.

   msg_        = ''   ; if (present(msg))        msg_        = msg
   error_code_ = -101 ; if (present(error_code)) error_code_ = error_code
   if (msg_ /='') write(stderr, '(A)') self%myrankstr//'abort '//msg_
   call MPI_ABORT(MPI_COMM_WORLD, error_code_, self%error)
   stop
   endsubroutine abort

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

   subroutine error_stop(self, msg)
   !< Stop run with error output.
   class(mpih_object), intent(inout)        :: self !< MPI handler.
   character(*),       intent(in), optional :: msg  !< Error message.
   character(:), allocatable                :: msg_ !< Error message, local variable.

   msg_ = '' ; if (present(msg)) msg_ = msg
   write(stderr, '(A)') self%myrankstr//'error stop '//msg_
   call self%finalize
   stop
   endsubroutine error_stop

   subroutine finalize(self)
   !< Handy MPI finalize wrapper.
   class(mpih_object) , intent(inout) :: self !< MPI handler.

   call MPI_FINALIZE(self%error)
   endsubroutine finalize

   subroutine initialize(self, do_mpi_init)
   !< Initialize MPI handler data.
   class(mpih_object) , intent(inout)        :: self                !< MPI handler.
   logical,             intent(in), optional :: do_mpi_init         !< Flag to activate MPI init call.
   integer(C_LONG)                           :: mem_free, mem_total !< CPU memory.

   if (present(do_mpi_init)) then
      if (do_mpi_init) call MPI_INIT(self%error)
   endif
   call MPI_COMM_SIZE(MPI_COMM_WORLD, self%procs_number, self%error)
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, self%error)
   self%myrankstr = '[myrank-'//trim(strz(self%myrank,6))//']'
   call cpuMemGetInfo(mem_total, mem_free)
   self%memory_avail = real(mem_total, R8P)/1e9/self%procs_number
   endsubroutine initialize

   function tictoc_timing(self) result(timing)
   !< Return the last tic toc timing.
   class(mpih_object) , intent(in) :: self   !< MPI handler.
   real(R8P)                       :: timing !< Last tic toc timing.

   timing = self%timing(2) - self%timing(1)
   endfunction tictoc_timing
endmodule adam_mpih_object
