!< ADAM, global RK singleton — single program-scope rk_object instance.
module adam_rk_global
!< ADAM, global RK singleton — single program-scope rk_object instance.

! ADAM classes, libraries, parameters
use :: adam_rk_object, only : rk_object

implicit none
private
public :: rk

type(rk_object), target :: rk !< Program-scope RK singleton.
endmodule adam_rk_global
