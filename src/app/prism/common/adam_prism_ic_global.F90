!< ADAM, global PRISM initial conditions shim — pointer into prism_common_object's value component (PRISM C.3 closure, issue #13).
module adam_prism_ic_global
!< ADAM, global PRISM initial conditions shim — pointer into the current realm's ic value component.
!<
!< Was a singleton owner until the PRISM C.3 closure (issue #13 D.4b follow-up).
!< Now a pointer shim aliasing `realm(N)%ic` after startup and rebound per realm
!< by `prism_common_object%bind_my_globals_forest` on the multi-realm path.

! PRISM modules
use :: adam_prism_ic_object, only : prism_ic_object

implicit none
private
public :: ic

type(prism_ic_object), pointer :: ic => null() !< Program-scope PRISM ic shim, bound by adam_prism_common_bind.
endmodule adam_prism_ic_global
