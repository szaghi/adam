!< ADAM, global PRISM fWLayer shim — pointer into prism_common_object's value component (PRISM C.3 closure, issue #13).
module adam_prism_fWLayer_global
!< ADAM, global PRISM fWLayer shim — pointer into the current realm's fWLayer value component.
!<
!< Was a singleton owner until the PRISM C.3 closure (issue #13 D.4b follow-up).
!< Now a pointer shim aliasing `realm(N)%fWLayer` after startup and rebound per realm
!< by `prism_common_object%bind_my_globals_forest` on the multi-realm path.

! PRISM modules
use :: adam_prism_fWLayer_object, only : prism_fWLayer_object

implicit none
private
public :: fWLayer

type(prism_fWLayer_object), pointer :: fWLayer => null() !< Program-scope PRISM fWLayer shim, bound by adam_prism_common_bind.
endmodule adam_prism_fWLayer_global
