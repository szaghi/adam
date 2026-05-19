!< ADAM, global PRISM external fields shim — pointer into prism_common_object's value component (PRISM C.3 closure, issue #13).
module adam_prism_external_fields_global
!< ADAM, global PRISM external fields shim — pointer into the current realm's external_fields value component.
!<
!< Was a singleton owner until the PRISM C.3 closure (issue #13 D.4b follow-up).
!< Now a pointer shim aliasing `realm(N)%external_fields` after startup and rebound per realm
!< by `prism_common_object%bind_my_globals_forest` on the multi-realm path.

! PRISM modules
use :: adam_prism_external_fields_object, only : prism_external_fields_object

implicit none
private
public :: external_fields

type(prism_external_fields_object), pointer :: external_fields => null() !< Program-scope PRISM external_fields shim, bound by adam_prism_common_bind.
endmodule adam_prism_external_fields_global
