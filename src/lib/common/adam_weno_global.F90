!< ADAM, global WENO singleton — single program-scope weno_object instance.
module adam_weno_global
!< ADAM, global WENO singleton — single program-scope weno_object instance.

! ADAM classes, libraries, parameters
use :: adam_weno_object, only : weno_object

implicit none
private
public :: weno

type(weno_object), target :: weno !< Program-scope WENO singleton.
endmodule adam_weno_global
