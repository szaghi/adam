!< ADAM, global IB singleton — single program-scope ib_object instance.
module adam_ib_global
!< ADAM, global IB singleton — single program-scope ib_object instance.

! ADAM classes, libraries, parameters
use :: adam_ib_object, only : ib_object

implicit none
private
public :: ib

type(ib_object), target :: ib !< Program-scope IB singleton.
endmodule adam_ib_global
