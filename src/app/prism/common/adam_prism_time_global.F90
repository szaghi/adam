!< ADAM, global PRISM time shim — pointer into prism_common_object's value component (PRISM C.3 closure, issue #13).
module adam_prism_time_global
!< ADAM, global PRISM time shim — pointer into the current realm's time value component.
!<
!< Was a singleton owner until the PRISM C.3 closure (issue #13 D.4b follow-up).
!< Now a pointer shim aliasing `realm(N)%time` after startup and rebound per realm
!< by `prism_common_object%bind_my_globals_forest` on the multi-realm path.

! PRISM modules
use :: adam_prism_time_object, only : prism_time_object

implicit none
private
public :: time

type(prism_time_object), pointer :: time => null() !< Program-scope PRISM time shim, bound by adam_prism_common_bind.
endmodule adam_prism_time_global
