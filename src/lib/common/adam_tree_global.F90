!< ADAM, global tree singleton — pointer shim into adam%tree (Step 1 of forest-of-trees migration, issue #10).
module adam_tree_global
!< ADAM, global tree singleton.
!<
!< Compatibility shim: storage lives in `adam%tree` (the `adam_adam_global`
!< singleton). See `adam_grid_global` for the rationale and the dependency
!< note. Bound by `adam_adam_bind`.

! ADAM classes, libraries, parameters
use :: adam_tree_object, only : tree_object

implicit none
private
public :: tree

type(tree_object), pointer :: tree => null() !< Program-scope tree shim, bound by adam_adam_bind.
endmodule adam_tree_global
