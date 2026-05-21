!< ADAM, common program-scope singletons — convenience re-export of all adam_*_global modules.
module adam_globals
!< ADAM, common program-scope singletons — convenience re-export of all adam_*_global modules.
!<
!< Provides a single USE point for all CPU-side program-scope singleton objects:
!<```fortran
!< use :: adam_globals, only : mpih, grid, field, maps, tree, weno, ib, rk
!<```
!< instead of listing every individual adam_*_global module.

! ADAM global singletons
use :: adam_mpih_global,  only : mpih
use :: adam_adam_global,  only : adam
use :: adam_grid_global,  only : grid
use :: adam_field_global, only : field
use :: adam_maps_global,  only : maps
use :: adam_tree_global,  only : tree

implicit none
private
public :: mpih
public :: adam
public :: grid
public :: field
public :: maps
public :: tree
endmodule adam_globals
