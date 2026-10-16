!< ADAM, MPI handler class definition.
module adam_mpih_object
!< ADAM, MPI handler class definition.

! third party modules
use :: penf
! sdk modules
use :: mpi
! intrinsic modules
use, intrinsic :: iso_c_binding,   only : C_INT, C_LONG
use, intrinsic :: iso_fortran_env, only : output_unit, stderr=>error_unit

implicit none
private

interface
   function c_usleep(useconds) bind(C, name='usleep') result(ierr)
   !< POSIX usleep: suspend the calling process for `useconds` microseconds.
   import :: C_INT
   integer(C_INT), value :: useconds !< Microseconds.
   integer(C_INT)        :: ierr     !< Status.
   endfunction c_usleep
endinterface
public :: mpih_object

type :: mpih_object
   !< MPI handler class.
   integer(I4P)              :: myrank=0_I4P        !< MPI rank process.
   character(:), allocatable :: myrankstr           !< MPI rank process stringified.
   integer(I4P)              :: procs_number=1_I4P  !< Number of MPI processes (MPI_COMM_WORLD).
   integer(I4P)              :: node_comm=0_I4P     !< Node-local (shared-memory) communicator.
   integer(I4P)              :: node_procs_number=1_I4P !< Number of MPI processes sharing THIS node's RAM.
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

   function cgroup_memory_limit() result(limit)
   !< Return the memory limit (bytes) of this process' cgroup, the tightest along its path to the root; -1 if none is
   !< set or there is no cgroup memory controller.
   !<
   !< The controller is located through /proc/self/cgroup: a v1 line `<id>:<controllers>:<path>` whose controller list
   !< holds `memory` (limit file `/sys/fs/cgroup/memory<path>/memory.limit_in_bytes`), else the v2 line `0::<path>`
   !< (`/sys/fs/cgroup<path>/memory.max`). The ancestors are read too: Slurm sets the --mem limit on the job and step
   !< cgroups while the task cgroup itself is unlimited (measured on medooza, cgroup v1).
   integer(I8P)              :: limit                              !< Memory limit (bytes), -1 if none.
   integer(I8P), parameter   :: UNLIMITED=4611686018427387904_I8P  !< 2**62: v1 reports "no limit" as about 2**63.
   character(999)            :: line                               !< Line buffer.
   character(:), allocatable :: path                               !< Cgroup path, walked up to the root.
   character(:), allocatable :: root                               !< Controller mount point.
   character(:), allocatable :: leaf                               !< Limit file name.
   integer(I8P)              :: value                              !< Limit read at one level.
   integer(I4P)              :: file_unit                          !< File unit.
   integer(I4P)              :: ios                                !< I/O status.
   integer(I4P)              :: c1, c2                             !< Positions of the two colons of a line.
   logical                   :: is_v1                              !< A v1 memory controller line was found.

   limit = -1_I8P
   path = ''
   is_v1 = .false.
   open(newunit=file_unit, file='/proc/self/cgroup', status='old', action='read', iostat=ios)
   if (ios /= 0) return
   do
      read(file_unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      c1 = index(line, ':')
      if (c1 == 0) cycle
      c2 = index(line(c1+1:), ':') + c1
      if (c2 == c1) cycle
      if (index(','//line(c1+1:c2-1)//',', ',memory,') > 0) then
         path = trim(line(c2+1:))
         is_v1 = .true.
         exit
      elseif (line(1:c2) == '0::') then
         path = trim(line(c2+1:))
      endif
   enddo
   close(file_unit)
   if (path == '') return
   if (is_v1) then
      root = '/sys/fs/cgroup/memory'
      leaf = '/memory.limit_in_bytes'
   else
      root = '/sys/fs/cgroup'
      leaf = '/memory.max'
   endif
   do
      if (path == '/') path = '' ! the root level
      open(newunit=file_unit, file=root//path//leaf, status='old', action='read', iostat=ios)
      if (ios == 0) then
         read(file_unit, *, iostat=ios) value ! v2 "max" (no limit) does not read as an integer
         if (ios == 0 .and. value > 0_I8P .and. value < UNLIMITED) then
            if (limit < 0_I8P .or. value < limit) limit = value
         endif
         close(file_unit)
      endif
      if (path == '') exit
      path = path(1:index(path, '/', back=.true.)-1)
   enddo
   endfunction cgroup_memory_limit

   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(mpih_object) , intent(in) :: self             !< MPI handler.
   character(len=:), allocatable   :: desc             !< Description.
   character(len=1), parameter     :: NL=new_line('a') !< New line character.

   desc =       self%myrankstr//'MPIH main data'//NL
   desc = desc//self%myrankstr//'  myrank:            '//trim(str(self%myrank      ))//NL
   desc = desc//self%myrankstr//'  procs_number:      '//trim(str(self%procs_number))//NL
   desc = desc//self%myrankstr//'  node procs_number: '//trim(str(self%node_procs_number))//NL
   desc = desc//self%myrankstr//'  memory_avail [GB]: '//trim(str(self%memory_avail))
   endfunction description

   subroutine error_stop(self, msg)
   !< Stop the run with error output, aborting every rank.
   !<
   !< The error may be rank-local (a check only some ranks evaluate): `MPI_ABORT` on `MPI_COMM_WORLD` terminates all
   !< the ranks. `MPI_FINALIZE` would block the failing rank until the others finalize, while they block forever in
   !< their next collective (issue #43). The message is flushed first and the rank sleeps 0.5 s before the abort, so the
   !< launcher forwards it (Open MPI kills the ranks at once, and output not yet forwarded is lost: measured). Without
   !< MPI (not yet initialized, or already finalized) the process stops alone.
   class(mpih_object), intent(inout)        :: self           !< MPI handler.
   character(*),       intent(in), optional :: msg            !< Error message.
   character(:), allocatable                :: msg_           !< Error message, local variable.
   logical                                  :: is_initialized !< MPI initialized.
   logical                                  :: is_finalized   !< MPI finalized.
   integer(I4P)                             :: ierr           !< MPI error status.
   integer(C_INT)                           :: cerr           !< usleep status.

   msg_ = '' ; if (present(msg)) msg_ = msg
   write(stderr, '(A)') self%myrankstr//'error stop '//msg_
   flush(stderr)
   flush(output_unit)
   call MPI_INITIALIZED(is_initialized, ierr)
   call MPI_FINALIZED(is_finalized, ierr)
   if (is_initialized .and. .not.is_finalized) then
      cerr = c_usleep(500000_C_INT)
      call MPI_ABORT(MPI_COMM_WORLD, 1_I4P, ierr)
   endif
   ! issue #25: plain `stop` exits with code 0 -- every error_stop looked SUCCESSFUL to
   ! shells/harnesses/CI. A failed run must report a nonzero exit code.
   stop 1
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
   integer(I8P)                              :: mem_limit              !< Cgroup memory limit (bytes), -1 if none.
   logical                                   :: is_initialized         !< Flag to check if MPI has been inizialied.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   myrankstr_char_length_ = 5 ; if (present(myrankstr_char_length)) myrankstr_char_length_ = myrankstr_char_length

   call MPI_INITIALIZED(is_initialized, self%error)
   if (.not.is_initialized) call MPI_INIT(self%error)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, self%procs_number, self%error)
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, self%error)
   self%myrankstr = '[mpi-'//trim(strz(self%myrank,myrankstr_char_length_))//']'
   if (verbose_) call self%print_message('mpih_object%initialize start')
   ! Node-local rank count. get_memory_info reads /proc/meminfo (PENF
   ! penf_allocatable_memory.F90:3799), i.e. THIS node's RAM -- a quantity shared
   ! only by the ranks resident on this node. Dividing it by the WORLD size
   ! under-estimates every rank's share by the number of nodes: at 2 nodes x 4
   ! ranks each rank computed mem/8 when its true share is mem/4, so adding nodes
   ! bought no capacity at all. Split by MPI_COMM_TYPE_SHARED and divide by that
   ! communicator's size instead -- the same primitive the device backends already
   ! use to bind one GPU per node-local rank (adam_mpih_nvf_object.F90:58).
   call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, self%node_comm, self%error)
   call MPI_COMM_SIZE(self%node_comm, self%node_procs_number, self%error)
   if (self%node_procs_number < 1_I4P) self%node_procs_number = 1_I4P
   call get_memory_info(mem_free=mem_free, mem_total=mem_total)
   ! The node's RAM is not the job's: under a batch scheduler (Slurm --mem) or a container the process lives in a memory
   ! cgroup whose limit /proc/meminfo does not show, and a pool sized from MemTotal is killed by the OOM handler (medooza,
   ! 2026-10-16: 1 TB nodes, 2 ranks -> 425 GB per rank in a 64 GB job). Cap by the cgroup limit; none -> unchanged.
   mem_limit = cgroup_memory_limit()
   if (mem_limit > 0_I8P) mem_total = min(mem_total, int(mem_limit/1024_I8P, kind(mem_total))) ! meminfo is in kB
   self%memory_avail = real(mem_total, R8P)/1e6/self%node_procs_number
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
