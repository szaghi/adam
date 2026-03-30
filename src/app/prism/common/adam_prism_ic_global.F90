!< ADAM, global PRISM initial conditions singleton — single program-scope prism_ic_object instance.
module adam_prism_ic_global
!< ADAM, global PRISM initial conditions singleton — single program-scope prism_ic_object instance.

! PRISM modules
use :: adam_prism_ic_object, only : prism_ic_object

implicit none
private
public :: ic

type(prism_ic_object), target :: ic !< Program-scope PRISM initial conditions singleton.
endmodule adam_prism_ic_global
