!< ADAM, global PRISM time singleton — single program-scope prism_time_object instance.
module adam_prism_time_global
!< ADAM, global PRISM time singleton — single program-scope prism_time_object instance.

! PRISM modules
use :: adam_prism_time_object, only : prism_time_object

implicit none
private
public :: time

type(prism_time_object), target :: time !< Program-scope PRISM time singleton.
endmodule adam_prism_time_global
