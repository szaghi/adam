!< ADAM, common program-scope singletons — convenience re-export of all adam_*_global modules.
module adam_globals
!< ADAM, common program-scope singletons — convenience re-export of all adam_*_global modules.
!<
!< Provides a single USE point for all CPU-side program-scope singleton objects:
!<```fortran
!< use :: adam_globals, only : mpih, grid, field, maps, weno, ib, rk
!<```
!< instead of listing every individual adam_*_global module.

! ADAM global singletons
use :: adam_mpih_global,  only : mpih
use :: adam_grid_global,  only : grid
use :: adam_field_global, only : field
use :: adam_maps_global,  only : maps
use :: adam_weno_global,  only : weno
use :: adam_ib_global,    only : ib
use :: adam_rk_global,    only : rk

implicit none
private
public :: mpih
public :: grid
public :: field
public :: maps
public :: weno
public :: ib
public :: rk
endmodule adam_globals
