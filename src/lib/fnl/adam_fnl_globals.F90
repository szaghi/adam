!< ADAM, FNL program-scope singletons — convenience re-export of all adam_fnl_*_global modules.
module adam_fnl_globals
!< ADAM, FNL program-scope singletons — convenience re-export of all adam_fnl_*_global modules.
!<
!< Provides a single USE point for all FNL (OpenACC GPU) program-scope singleton objects:
!<```fortran
!< use :: adam_fnl_globals, only : mpih_fnl, field_fnl, weno_fnl, ib_fnl, rk_fnl
!<```
!< instead of listing every individual adam_fnl_*_global module.

! ADAM FNL global singletons
use :: adam_fnl_mpih_global,  only : mpih_fnl
use :: adam_fnl_field_global, only : field_fnl
use :: adam_fnl_weno_global,  only : weno_fnl
use :: adam_fnl_ib_global,    only : ib_fnl
use :: adam_fnl_rk_global,    only : rk_fnl

implicit none
private
public :: mpih_fnl
public :: field_fnl
public :: weno_fnl
public :: ib_fnl
public :: rk_fnl
endmodule adam_fnl_globals
