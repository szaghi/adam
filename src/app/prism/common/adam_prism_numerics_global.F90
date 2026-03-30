!< ADAM, global PRISM numerics singleton — single program-scope prism_numerics_object instance.
module adam_prism_numerics_global
!< ADAM, global PRISM numerics singleton — single program-scope prism_numerics_object instance.

! PRISM modules
use :: adam_prism_numerics_object, only : prism_numerics_object

implicit none
private
public :: numerics

type(prism_numerics_object), target :: numerics !< Program-scope PRISM numerics singleton.
endmodule adam_prism_numerics_global
