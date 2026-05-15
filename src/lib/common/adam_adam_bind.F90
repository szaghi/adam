!< ADAM, binder for legacy grid/field/tree/maps shims (Step 1 of forest-of-trees migration, issue #10).
module adam_adam_bind
!< Binder module that aliases the four legacy shim singletons
!< (`grid`, `field`, `tree`, `maps`) into the corresponding value
!< components of the program-scope `adam` (the `adam_adam_global`
!< singleton). Must be called once after `adam%initialize`.
!<
!< This logic lives outside the four shim modules to break a circular
!< dependency: `adam_field_object`/`adam_tree_object`/`adam_maps_object`
!< already import their sibling shims, and `adam_adam_object` imports
!< those object modules — putting the binder inside the shims would close
!< the cycle. Keeping it here makes the dependency strictly
!<```
!<   adam_*_global  →  adam_*_object
!<   adam_adam_object  →  adam_*_object
!<   adam_adam_global  →  adam_adam_object
!<   adam_adam_bind    →  adam_adam_global, adam_*_global, adam_adam_object
!<```
!< with no cycles.

! ADAM classes, libraries, parameters
use :: adam_adam_object,  only : adam_object
! ADAM singleton objects
use :: adam_adam_global,  only : adam
use :: adam_field_global, only : field
use :: adam_grid_global,  only : grid
use :: adam_maps_global,  only : maps
use :: adam_tree_global,  only : tree

implicit none
private
public :: bind_globals_to_adam

contains
   subroutine bind_globals_to_adam
   !< Alias the four legacy shim singletons into `adam`'s value components.
   !< Idempotent: safe to call more than once.

   grid  => adam%grid
   field => adam%field
   tree  => adam%tree
   maps  => adam%maps
   endsubroutine bind_globals_to_adam
endmodule adam_adam_bind
