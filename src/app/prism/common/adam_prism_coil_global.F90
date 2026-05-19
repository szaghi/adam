!< ADAM, global PRISM coil shim — pointer into prism_common_object's value component (PRISM C.3 closure, issue #13).
module adam_prism_coil_global
!< ADAM, global PRISM coil shim — pointer into the current realm's coil value component.
!<
!< Was a singleton owner until the PRISM C.3 closure (issue #13 D.4b follow-up).
!< Now a pointer shim aliasing `realm(N)%coil` after startup and rebound per realm
!< by `prism_common_object%bind_my_globals_forest` on the multi-realm path.

! PRISM modules
use :: adam_prism_coil_object, only : prism_coil_object

implicit none
private
public :: coil

type(prism_coil_object), pointer :: coil => null() !< Program-scope PRISM coil shim, bound by adam_prism_common_bind.
endmodule adam_prism_coil_global
