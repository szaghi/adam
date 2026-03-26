!< ADAM, global FNL IB singleton — single program-scope ib_fnl_object instance.
module adam_fnl_ib_global
!< ADAM, global FNL IB singleton — single program-scope ib_fnl_object instance.
!<
!< Requires explicit initialization after the CPU-side ib global singleton and
!< the field_fnl global are populated:
!<```fortran
!< ib = self%ib   ! populate CPU singleton
!< call ib_fnl%initialize()
!<```

! ADAM FNL classes
use :: adam_fnl_ib_object, only : ib_fnl_object

implicit none
private
public :: ib_fnl

type(ib_fnl_object), target :: ib_fnl !< Program-scope IB FNL singleton.
endmodule adam_fnl_ib_global
