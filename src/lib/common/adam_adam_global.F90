!< ADAM, global adam singleton — single program-scope adam_object instance.
module adam_adam_global
!< ADAM, global adam singleton — single program-scope adam_object instance.
!<
!< Step 1 of the forest-of-trees migration (issue #10): adam_object absorbs
!< grid/tree/field/maps as value components. Exposing the one adam_object as
!< a module-scope singleton gives the grid/field/tree/maps shim modules a
!< stable anchor to alias into (`grid => adam%grid`, etc.) without touching
!< the 28 existing consumers of those shims.

! ADAM classes, libraries, parameters
use :: adam_adam_object, only : adam_object

implicit none
private
public :: adam

type(adam_object), target :: adam !< Program-scope adam singleton.
endmodule adam_adam_global
