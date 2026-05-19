!< ADAM, global PRISM RK-BC shim — pointer into prism_common_object's value component (PRISM C.3 closure, issue #13).
module adam_prism_rk_bc_global
!< ADAM, global PRISM RK-BC shim — pointer into the current realm's rk_bc value component.
!<
!< Was a singleton owner until the PRISM C.3 closure (issue #13 D.4b follow-up).
!< Now a pointer shim aliasing `realm(N)%rk_bc` after startup and rebound per realm
!< by `prism_common_object%bind_my_globals_forest` on the multi-realm path.

! PRISM modules
use :: adam_prism_rk_bc_object, only : prism_rk_bc_object

implicit none
private
public :: rk_bc

type(prism_rk_bc_object), pointer :: rk_bc => null() !< Program-scope PRISM rk_bc shim, bound by adam_prism_common_bind.
endmodule adam_prism_rk_bc_global
