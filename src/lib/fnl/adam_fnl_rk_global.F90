!< ADAM, global FNL RK singleton — single program-scope rk_fnl_object instance.
module adam_fnl_rk_global
!< ADAM, global FNL RK singleton — single program-scope rk_fnl_object instance.
!<
!< Requires explicit initialization after the CPU-side rk global singleton and
!< the field/grid globals are populated:
!<```fortran
!< rk = self%rk   ! populate CPU singleton
!< call rk_fnl%initialize()
!<```

! ADAM FNL classes
use :: adam_fnl_rk_object, only : rk_fnl_object

implicit none
private
public :: rk_fnl

type(rk_fnl_object), target :: rk_fnl !< Program-scope RK FNL singleton.
endmodule adam_fnl_rk_global
