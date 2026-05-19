!< ADAM, global adam shim — pointer shim into realm(N)%adam (C.3 closure of issue #10).
module adam_adam_global
!< ADAM, global adam shim — pointer shim into realm(N)%adam (C.3 closure of issue #10).
!<
!< History: was the **owner** of `adam_singleton` after Step 1 of the
!< forest-of-trees migration (the grid/tree/field/maps shims aliased into
!< its value components). C.3-closure (2026-05-19, follow-up to D.3 / issue
!< #13) promotes `adam_object` to a value component on `realm_object`,
!< making this module a **pointer shim** like the seven other
!< `adam_*_global` modules — `adam` now aliases the last realm initialized
!< by the forest (`realm(N)`), via the `bind_globals_to_adam` binder.
!<
!< For N=1 (current rmf usage) the shim aliases `realm(1)%adam` and the 28
!< consumer files that read `adam%...` through this module see exactly the
!< same data they did before. For N>1 (the rmf-2realm use case under #13),
!< the shim aliases `realm(N)%adam` after `forest%initialize` finishes —
!< matching the documented C.3 behaviour of the other shims, and matching
!< the Phase C driver comment ("the shims alias realm(N)").

! ADAM classes, libraries, parameters
use :: adam_adam_object, only : adam_object

implicit none
private
public :: adam

type(adam_object), pointer :: adam => null() !< Program-scope adam shim, bound by adam_adam_bind to the last initialized realm's adam value component.
endmodule adam_adam_global
