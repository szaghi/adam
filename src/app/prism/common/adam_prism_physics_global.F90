!< ADAM, global PRISM physics singleton — single program-scope prism_physics_object instance.
module adam_prism_physics_global
!< ADAM, global PRISM physics singleton — single program-scope prism_physics_object instance.

! PRISM modules
use :: adam_prism_physics_object, only : prism_physics_object

implicit none
private
public :: physics

type(prism_physics_object), target :: physics !< Program-scope PRISM physics singleton.
endmodule adam_prism_physics_global
