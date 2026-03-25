!< ADAM, global FNL MPI handler singleton — single program-scope mpih_fnl_object instance.
module adam_global_mpih_fnl
!< ADAM, global FNL MPI handler singleton — single program-scope mpih_fnl_object instance.
!<
!< Provides the GPU-aware MPI handler (FUNDAL mpih_object) as a program-scope singleton,
!< mirroring the CPU-side adam_global_mpih pattern.
!<
!< Requires explicit initialization before any FNL object is constructed.
!< GPU applications call it with device init enabled:
!<```fortran
!< call mpih_fnl%initialize(do_mpi_init=.true., do_device_init=.true., verbose=.true.)
!<```
!< Non-GPU callers (e.g. unit tests) may pass `do_mpi_init=.false.` if MPI is already up.

! ADAM FNL classes
use :: adam_fnl_mpih_object, only : mpih_fnl_object

implicit none
private
public :: mpih_fnl

type(mpih_fnl_object), target :: mpih_fnl !< Program-scope GPU MPI handler singleton.
endmodule adam_global_mpih_fnl
