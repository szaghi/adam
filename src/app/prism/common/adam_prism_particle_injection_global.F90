!< ADAM, global PRISM particle injection shim — pointer into prism_common_object's value component (PRISM C.3 closure, issue #13).
module adam_prism_particle_injection_global
!< ADAM, global PRISM particle injection shim — pointer into the current realm's particle_injection value component.
!<
!< Was a singleton owner until the PRISM C.3 closure (issue #13 D.4b follow-up).
!< Now a pointer shim aliasing `realm(N)%particle_injection` after startup and rebound per realm
!< by `prism_common_object%bind_my_globals_forest` on the multi-realm path.

! PRISM modules
use :: adam_prism_particle_injection_object, only : prism_particle_injection_object

implicit none
private
public :: particle_injection

type(prism_particle_injection_object), pointer :: particle_injection => null() !< Program-scope PRISM particle_injection shim, bound by adam_prism_common_bind.
endmodule adam_prism_particle_injection_global
