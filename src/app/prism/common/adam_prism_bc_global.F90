!< ADAM, global PRISM boundary conditions singleton — single program-scope prism_bc_object instance.
module adam_prism_bc_global
!< ADAM, global PRISM boundary conditions singleton — single program-scope prism_bc_object instance.

! PRISM modules
use :: adam_prism_bc_object, only : prism_bc_object

implicit none
private
public :: bc

type(prism_bc_object), target :: bc !< Program-scope PRISM boundary conditions singleton.
endmodule adam_prism_bc_global
