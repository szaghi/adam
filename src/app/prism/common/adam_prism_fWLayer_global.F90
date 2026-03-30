!< ADAM, global PRISM fWLayer singleton — single program-scope prism_fWLayer_object instance.
module adam_prism_fWLayer_global
!< ADAM, global PRISM fWLayer singleton — single program-scope prism_fWLayer_object instance.

! PRISM modules
use :: adam_prism_fWLayer_object, only : prism_fWLayer_object

implicit none
private
public :: fWLayer

type(prism_fWLayer_object), target :: fWLayer !< Program-scope PRISM fWLayer singleton.
endmodule adam_prism_fWLayer_global
