!< ADAM, global PRISM RK-BC singleton — single program-scope prism_rk_bc_object instance.
module adam_prism_rk_bc_global
!< ADAM, global PRISM RK-BC singleton — single program-scope prism_rk_bc_object instance.

! PRISM modules
use :: adam_prism_rk_bc_object, only : prism_rk_bc_object

implicit none
private
public :: rk_bc

type(prism_rk_bc_object), target :: rk_bc !< Program-scope PRISM RK-BC singleton.
endmodule adam_prism_rk_bc_global
