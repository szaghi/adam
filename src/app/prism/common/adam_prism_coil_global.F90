!< ADAM, global PRISM coil singleton — single program-scope prism_coil_object instance.
module adam_prism_coil_global
!< ADAM, global PRISM coil singleton — single program-scope prism_coil_object instance.

! PRISM modules
use :: adam_prism_coil_object, only : prism_coil_object

implicit none
private
public :: coil

type(prism_coil_object), target :: coil !< Program-scope PRISM coil singleton.
endmodule adam_prism_coil_global
