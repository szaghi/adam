!< ADAM, global PRISM physics shim — pointer into prism_common_object's value component (PRISM C.3 closure, issue #13).
module adam_prism_physics_global
!< ADAM, global PRISM physics shim — pointer into the current realm's physics value component.
!<
!< Was a singleton owner until the PRISM C.3 closure (issue #13 D.4b follow-up).
!< Now a pointer shim aliasing `realm(N)%physics` after startup and rebound per realm
!< by `prism_common_object%bind_my_globals_forest` on the multi-realm path.

! PRISM modules
use :: adam_prism_physics_object, only : prism_physics_object

implicit none
private
public :: physics

type(prism_physics_object), pointer :: physics => null() !< Program-scope PRISM physics shim, bound by adam_prism_common_bind.
endmodule adam_prism_physics_global
