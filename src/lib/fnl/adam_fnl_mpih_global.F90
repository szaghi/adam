!< ADAM, global FNL MPI handler singleton — single program-scope mpih_fnl_object instance.
module adam_fnl_mpih_global
!< ADAM, global FNL MPI handler singleton — single program-scope mpih_fnl_object instance.
!<
!< Provides the GPU-aware MPI handler (FUNDAL mpih_object) as a program-scope singleton,
!< mirroring the CPU-side adam_mpih_global pattern.
!<
!< Unlike the six per-realm FNL objects (field/ib/rk/weno/coil/fwlayer), the MPI handler
!< is genuinely program-scope: it wraps the process rank/size, the bound GPU device, and
!< the local communicator — all one-per-process, not one-per-realm. FUNDAL's mpih_object
!< even associates its device members (`mydev`, `local_comm`, ...) to fundal_env global
!< variables, so there is exactly one device context regardless of realm count. Promoting
!< it to a per-realm value component duplicates the handler and re-runs device init per
!< realm, corrupting the CUDA primary context on the multi-realm path. True singleton here.
!<
!< Init-once discipline: FUNDAL's `initialize` is `intent(out)` (it wipes self) and runs the
!< device init + communicator split under `do_device_init`. It must therefore be called
!< exactly once per process. `is_initialized` guards that: the first realm initializes the
!< singleton; subsequent realms skip. See `prism_fnl_object%initialize_prism`.
!<```fortran
!< if (.not. mpih_fnl_is_initialized) then
!<    call mpih_fnl%initialize(do_mpi_init=..., do_device_init=.true., verbose=.true.)
!<    mpih_fnl_is_initialized = .true.
!<    mpih_fnl_device_sharers = count_device_sharers()
!< endif
!<```
!<
!< `mpih_fnl_device_sharers` is the number of node-local ranks bound to this process's device, itself included. The
!< apps divide the device capacity budget by it (and by the realm count), so that ranks sharing one GPU do not each
!< size their block pool for the whole device (issue #96, 6.2).

! ADAM FNL classes
use :: adam_fnl_mpih_object, only : mpih_fnl_object
! external libraries
use :: fundal, only : dev_is_host_fallback
use :: mpi
use :: penf,   only : I4P

implicit none
private
public :: count_device_sharers
public :: mpih_fnl
public :: mpih_fnl_device_sharers
public :: mpih_fnl_is_initialized

type(mpih_fnl_object), target :: mpih_fnl                        !< Program-scope GPU MPI handler singleton.
integer(I4P)                  :: mpih_fnl_device_sharers=1_I4P   !< Node-local ranks bound to this process's device.
logical                       :: mpih_fnl_is_initialized=.false. !< Init-once guard (FUNDAL initialize is intent(out)).

contains
   function count_device_sharers() result(sharers)
   !< Count the node-local ranks bound to the same physical device as this process, itself included.
   !<
   !< Collective on `mpih_fnl%local_comm`; call it once, after `mpih_fnl%initialize(do_device_init=.true.)`.
   !<
   !< A device is identified by the device list the process sees (`CUDA_VISIBLE_DEVICES`, `ROCR_VISIBLE_DEVICES`,
   !< `HIP_VISIBLE_DEVICES`) plus the number FUNDAL picked in it, `mydev = mod(local_rank, devs_number)`. Ranks
   !< launched with the same list share device `mydev` when they picked the same number (`mpirun -np 2` on one GPU);
   !< ranks given different lists do not (per-task binding, `srun --gpu-bind`). A launcher that shows every rank a
   !< distinct GPU under the same list (e.g. per-task cgroups, all ranks seeing `CUDA_VISIBLE_DEVICES=0`) makes the
   !< count err high: a smaller budget, never an oversubscribed device.
   integer(I4P)                    :: sharers     !< Ranks bound to this process's device.
   integer(I4P),       parameter   :: KEY_LEN=512 !< Device key length.
   character(20),      parameter   :: VISIBLE(3)=[character(20) :: 'CUDA_VISIBLE_DEVICES', 'ROCR_VISIBLE_DEVICES', &
                                                                   'HIP_VISIBLE_DEVICES']
                                                  !< Device-list environment variables.
   character(KEY_LEN)              :: key         !< This process's device key.
   character(KEY_LEN), allocatable :: keys(:)     !< Device keys of the node-local ranks.
   character(KEY_LEN)              :: list        !< One device list.
   character(12)                   :: dev         !< Device number, stringified.
   integer(I4P)                    :: ranks       !< Node-local ranks.
   integer(I4P)                    :: i           !< Counter.
   integer(I4P)                    :: ierr        !< MPI error flag.

   sharers = 1_I4P
   if (dev_is_host_fallback() .or. mpih_fnl%devs_number < 1_I4P) return
   key = ''
   do i = 1, size(VISIBLE, dim=1)
      call get_environment_variable(trim(VISIBLE(i)), list)
      key = trim(key)//trim(VISIBLE(i))//'='//trim(list)//';'
   enddo
   write(dev, '(I0)') mpih_fnl%mydev
   key = trim(key)//'device='//trim(dev)
   call MPI_COMM_SIZE(mpih_fnl%local_comm, ranks, ierr)
   allocate(keys(ranks))
   call MPI_ALLGATHER(key, KEY_LEN, MPI_CHARACTER, keys, KEY_LEN, MPI_CHARACTER, mpih_fnl%local_comm, ierr)
   sharers = count(keys == key)
   endfunction count_device_sharers
endmodule adam_fnl_mpih_global
