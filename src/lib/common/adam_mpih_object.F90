!< ADAM, MPI handler class definition.
module adam_mpih_object
!< ADAM, MPI handler class definition.

! third party modules
use :: penf
! sdk modules
use :: mpi
! intrinsic modules
use, intrinsic :: iso_c_binding,   only : C_LONG
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

implicit none
private
public :: mpih_object

type :: mpih_object
   !< MPI handler class.
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
      procedure, pass(self) :: description   !< Return pretty-printed object description.
      procedure, pass(self) :: error_stop    !< Stop run with error output.
      procedure, pass(self) :: finalize      !< Handy MPI finalize wrapper.
      procedure, pass(self) :: initialize    !< Initialize MPI handler data.
      procedure, pass(self) :: print_message !< Print a message on stdout with rank prefix.
      procedure, pass(self) :: tictoc_timing !< Return the last tic toc timing.
      procedure, pass(self) :: tic           !< Start a tic toc timing.
      procedure, pass(self) :: toc           !< Stop  a tic toc timing.
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

   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(mpih_object) , intent(in) :: self             !< MPI handler.
   character(len=:), allocatable   :: desc             !< Description.
   character(len=1), parameter     :: NL=new_line('a') !< New line character.

   desc =       self%myrankstr//'MPIH main data'//NL
   desc = desc//self%myrankstr//'  myrank:            '//trim(str(self%myrank      ))//NL
   desc = desc//self%myrankstr//'  procs_number:      '//trim(str(self%procs_number))//NL
   desc = desc//self%myrankstr//'  memory_avail [GB]: '//trim(str(self%memory_avail))
   endfunction description

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

   subroutine initialize(self, do_mpi_init, do_device_init, myrankstr_char_length, verbose)
   !< Initialize MPI handler data.
   class(mpih_object) , intent(inout)        :: self                   !< MPI handler.
   logical,             intent(in), optional :: do_mpi_init            !< Flag to activate MPI init call.
   logical,             intent(in), optional :: do_device_init         !< Flag to activate device init call (used by backends).
   integer(I4P),        intent(in), optional :: myrankstr_char_length  !< MPI ID string length.
   logical,             intent(in), optional :: verbose                !< Trigger verbose output.
   logical                                   :: verbose_               !< Trigger verbose output, local variable.
   integer(I4P)                              :: myrankstr_char_length_ !< MPI ID string length, local variable.
   integer(C_LONG)                           :: mem_free, mem_total    !< CPU memory.
   logical                                   :: is_initialized         !< Flag to check if MPI has been inizialied.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   myrankstr_char_length_ = 5 ; if (present(myrankstr_char_length)) myrankstr_char_length_ = myrankstr_char_length

   call MPI_INITIALIZED(is_initialized, self%error)
   if (.not.is_initialized) call MPI_INIT(self%error)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, self%procs_number, self%error)
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, self%error)
   self%myrankstr = '[mpi-'//trim(strz(self%myrank,myrankstr_char_length_))//']'
   if (verbose_) call self%print_message('mpih_object%initialize start')
   call get_memory_info(mem_free=mem_free, mem_total=mem_total)
   self%memory_avail = real(mem_total, R8P)/1e6/self%procs_number
   if (verbose_) then
      print '(A)', self%description()
      call self%print_message('mpih_object%initialize finish')
   endif
   endsubroutine initialize

   subroutine print_message(self, msg)
   !< Print a message on stdout with rank prefix.
   class(mpih_object) , intent(in) :: self !< MPI handler.
   character(*),        intent(in) :: msg  !< Message to print.

   print '(A)', self%myrankstr//trim(adjustl(msg))
   endsubroutine print_message

   function tictoc_timing(self) result(timing)
   !< Return the last tic toc timing.
   class(mpih_object) , intent(in) :: self   !< MPI handler.
   real(R8P)                       :: timing !< Last tic toc timing.

   timing = self%timing(2) - self%timing(1)
   endfunction tictoc_timing

   subroutine tic(self)
   !< Start a tic toc timing.
   class(mpih_object) , intent(inout) :: self !< MPI handler.

   self%timing(1) = MPI_WTIME()
   endsubroutine tic

   function toc(self) result(timing)
   !< Stop a tic toc timing.
   class(mpih_object) , intent(inout) :: self !< MPI handler.
   real(R8P)                          :: timing !< Tic toc timing.

   self%timing(2) = MPI_WTIME()
   timing = self%tictoc_timing()
   endfunction toc
endmodule adam_mpih_object
