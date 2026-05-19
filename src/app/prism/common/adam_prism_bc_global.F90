!< ADAM, global PRISM boundary conditions shim — pointer into prism_common_object's value component (PRISM C.3 closure, issue #13).
module adam_prism_bc_global
!< ADAM, global PRISM boundary conditions shim — pointer into the current realm's bc value component.
!<
!< Was a singleton owner (`type(prism_bc_object), target :: bc`) until the PRISM C.3
!< closure (issue #13 D.4b follow-up). For N=1 (single-realm rmf) the shim aliases
!< `realm(1)%bc` and the dozens of consumer sites that read `bc%...` see the same
!< data they did before. For N>1 (rmf-2realm) the shim is rebound per realm by
!< `prism_common_object%bind_my_globals_forest` before each realm's per-step TBPs.

! PRISM modules
use :: adam_prism_bc_object, only : prism_bc_object

implicit none
private
public :: bc

type(prism_bc_object), pointer :: bc => null() !< Program-scope PRISM bc shim, bound by adam_prism_common_bind.
endmodule adam_prism_bc_global
