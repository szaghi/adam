!< ADAM, common program-scope singletons — convenience re-export of all adam_*_global modules.
module adam_globals
!< ADAM, common program-scope singletons — convenience re-export of all adam_*_global modules.
!<
!< Provides a single USE point for all CPU-side program-scope singleton objects:
!<```fortran
!< use :: adam_globals, only : mpih
!<```
!< instead of listing every individual adam_*_global module.

! ADAM global singletons
use :: adam_mpih_global,  only : mpih

implicit none
private
public :: mpih
endmodule adam_globals
