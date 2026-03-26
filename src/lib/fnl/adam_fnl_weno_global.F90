!< ADAM, global FNL WENO singleton — single program-scope weno_fnl_object instance.
module adam_fnl_weno_global
!< ADAM, global FNL WENO singleton — single program-scope weno_fnl_object instance.
!<
!< Requires explicit initialization after the CPU-side weno global singleton is populated:
!<```fortran
!< weno = self%weno   ! populate CPU singleton
!< call weno_fnl%initialize()
!<```

! ADAM FNL classes
use :: adam_fnl_weno_object, only : weno_fnl_object

implicit none
private
public :: weno_fnl

type(weno_fnl_object), target :: weno_fnl !< Program-scope WENO FNL singleton.
endmodule adam_fnl_weno_global
