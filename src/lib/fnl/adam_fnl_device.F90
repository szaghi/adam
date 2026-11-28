!< ADAM FNL, the device binding of the process: HERALD's (`mpih%device`), then FUNDAL's `dev_init`.
module adam_fnl_device
!< ADAM FNL, the device binding of the process: HERALD's (`mpih%device`), then FUNDAL's `dev_init`.
!<
!< The FNL backend's single handler is `mpih` (`adam_mpih_global`, HERALD's `mpih_object`, issue #158): there is no FNL
!< handler and no module state here. `bind_device` binds the rank to its device, once per process, from every FNL entry
!< point (the backends' `initialize_process_context`, the FNL test drivers):
!<
!< - HERALD binds the device, `mod(node rank, devices)`, and counts the ranks of the node that share it
!<   (`mpih%device%sharers`, issue #96). HERALD's `initialize` is once-only: a call asking for the device after `mpih`
!<   was initialized without it would be ignored without a message, so when `mpih` is already up (a test driver, the
!<   CPU layer) the device is bound through `mpih%device%initialize`, collective on `mpih%node_comm`;
!< - FUNDAL's `dev_init(local_rank=node rank)` is still called: it fills `fundal_env`, which FUNDAL's data layer reads
!<   (`mydev`, the default device of every `dev_alloc`/`dev_free`; `devtype`; `myhos`), and binds the same device by the
!<   same rule. The two must agree on every rank, or the run stops. Which of the two binders survives is a later issue.
!<
!< The FNL block budget keeps FUNDAL's device memory (`dev_memory_total`), bit for bit; HERALD's is logged next to it.

! ADAM classes
use :: adam_mpih_global, only : mpih
! external libraries
use :: fundal, only : dev_init, dev_memory_total, devs_number, mydev
use :: penf,   only : I4P, str

implicit none
private
public :: bind_device

contains
   subroutine bind_device(verbose)
   !< Bind this rank to its device, once per process: HERALD's binding (with `mpih`'s initialization if nobody did it),
   !< then FUNDAL's `dev_init` on the same node rank; stop if the two picked different devices. Collective on the node.
   logical, intent(in), optional :: verbose !< Print the handler's description (default .false.).
   logical                       :: verbose_ !< Print the description, local variable.
   integer(I4P)                  :: error    !< Error of HERALD's binding.
   character(:), allocatable     :: msg      !< Error or warning of HERALD's binding.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (mpih%device%devices_number > 0_I4P) return ! already bound (a second entry point of the same process)
   if (.not.mpih%is_initialized) then
      call mpih%initialize(do_device_init=.true., verbose=verbose_)
   else
      call mpih%device%initialize(node_comm=mpih%node_comm, node_myrank=mpih%node_myrank, error=error, msg=msg)
      if (error /= 0_I4P) call mpih%error_stop(msg='bind_device: '//msg)
      if (len(msg) > 0) call mpih%warning(msg=msg)
      if (verbose_) print '(A)', mpih%device%description(prefix=mpih%prefix())
   endif
   call dev_init(local_rank=mpih%node_myrank)
   if ((.not.mpih%device%is_host) .and. mydev /= mpih%device%mydev) &
      call mpih%error_stop(msg='bind_device: FUNDAL bound device '//trim(str(mydev, no_sign=.true.))//', HERALD device '// &
                               trim(str(mpih%device%mydev, no_sign=.true.))//': the two binders disagree')
   call mpih%print_message('device: FUNDAL mydev '//trim(str(mydev, no_sign=.true.))//                                &
                           ' of '//trim(str(devs_number, no_sign=.true.))//                                            &
                           ', dev_memory_total '//trim(str(dev_memory_total, no_sign=.true.))//' B; HERALD mydev '//     &
                           trim(str(mpih%device%mydev, no_sign=.true.))//                                              &
                           ', memory_total '//trim(str(mpih%device%memory_total, no_sign=.true.))//' B'//              &
                           '; shared by '//trim(str(mpih%device%sharers, no_sign=.true.))//' rank(s)')
   endsubroutine bind_device
endmodule adam_fnl_device
