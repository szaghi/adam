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
!< it to a per-realm value component (the C.3 over-application, issue #13 D.4a) duplicated
!< the handler and re-ran device init per realm, corrupting the CUDA primary context on the
!< multi-realm path. Reverted to a true singleton here.
!<
!< Init-once discipline: FUNDAL's `initialize` is `intent(out)` (it wipes self) and runs the
!< device init + communicator split under `do_device_init`. It must therefore be called
!< exactly once per process. `is_initialized` guards that: the first realm initializes the
!< singleton; subsequent realms skip. See `prism_fnl_object%initialize_prism`.
!<```fortran
!< if (.not. mpih_fnl_is_initialized) then
!<    call mpih_fnl%initialize(do_mpi_init=..., do_device_init=.true., verbose=.true.)
!<    mpih_fnl_is_initialized = .true.
!< endif
!<```

! ADAM FNL classes
use :: adam_fnl_mpih_object, only : mpih_fnl_object

implicit none
private
public :: mpih_fnl
public :: mpih_fnl_is_initialized

type(mpih_fnl_object), target :: mpih_fnl                      !< Program-scope GPU MPI handler singleton.
logical                       :: mpih_fnl_is_initialized=.false. !< Init-once guard (FUNDAL initialize is intent(out)).
endmodule adam_fnl_mpih_global
