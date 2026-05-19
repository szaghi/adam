!< ADAM, global PRISM PIC shim — pointer into prism_common_object's value component (PRISM C.3 closure, issue #13).
module adam_prism_pic_global
!< ADAM, global PRISM PIC shim — pointer into the current realm's pic value component.
!<
!< Was a singleton owner until the PRISM C.3 closure (issue #13 D.4b follow-up).
!< Now a pointer shim aliasing `realm(N)%pic` after startup and rebound per realm
!< by `prism_common_object%bind_my_globals_forest` on the multi-realm path.

! PRISM modules
use :: adam_prism_pic_object, only : prism_pic_object

implicit none
private
public :: pic

type(prism_pic_object), pointer :: pic => null() !< Program-scope PRISM pic shim, bound by adam_prism_common_bind.
endmodule adam_prism_pic_global
