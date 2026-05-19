!< ADAM, global PRISM RK-PIC shim — pointer into prism_common_object's value component (PRISM C.3 closure, issue #13).
module adam_prism_rk_pic_global
!< ADAM, global PRISM RK-PIC shim — pointer into the current realm's rk_pic value component.
!<
!< Was a singleton owner until the PRISM C.3 closure (issue #13 D.4b follow-up).
!< Now a pointer shim aliasing `realm(N)%rk_pic` after startup and rebound per realm
!< by `prism_common_object%bind_my_globals_forest` on the multi-realm path.

! PRISM modules
use :: adam_prism_rk_pic_object, only : prism_rk_pic_object

implicit none
private
public :: rk_pic

type(prism_rk_pic_object), pointer :: rk_pic => null() !< Program-scope PRISM rk_pic shim, bound by adam_prism_common_bind.
endmodule adam_prism_rk_pic_global
