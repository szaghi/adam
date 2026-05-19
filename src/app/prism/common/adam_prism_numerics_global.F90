!< ADAM, global PRISM numerics shim — pointer into prism_common_object's value component (PRISM C.3 closure, issue #13).
module adam_prism_numerics_global
!< ADAM, global PRISM numerics shim — pointer into the current realm's numerics value component.
!<
!< Was a singleton owner until the PRISM C.3 closure (issue #13 D.4b follow-up).
!< Now a pointer shim aliasing `realm(N)%numerics` after startup and rebound per realm
!< by `prism_common_object%bind_my_globals_forest` on the multi-realm path.

! PRISM modules
use :: adam_prism_numerics_object, only : prism_numerics_object

implicit none
private
public :: numerics

type(prism_numerics_object), pointer :: numerics => null() !< Program-scope PRISM numerics shim, bound by adam_prism_common_bind.
endmodule adam_prism_numerics_global
