!< ADAM, FNL program-scope singletons — convenience re-export of all adam_fnl_*_global modules.
module adam_fnl_globals
!< ADAM, FNL program-scope singletons — convenience re-export of all adam_fnl_*_global modules.
!<
!< Provides a single USE point for the FNL (OpenACC GPU) program-scope singleton:
!<```fortran
!< use :: adam_fnl_globals, only : mpih_fnl
!<```
!< Only mpih_fnl remains a program-scope shim — the base lib/fnl objects read it
!< directly (no realm self in scope). The former field/weno/ib/rk shims were
!< dropped in the consolidation pass; those objects are now per-realm value
!< components reached through self%.

! ADAM FNL global singletons
use :: adam_fnl_mpih_global, only : mpih_fnl

implicit none
private
public :: mpih_fnl
endmodule adam_fnl_globals
