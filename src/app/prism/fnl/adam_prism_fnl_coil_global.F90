!< ADAM, global PRISM FNL coil singleton — single program-scope prism_fnl_coil_object instance.
module adam_prism_fnl_coil_global
!< ADAM, global PRISM FNL coil singleton — single program-scope prism_fnl_coil_object instance.
!<
!< Requires explicit initialization after the CPU-side coil object is populated:
!<```fortran
!< call coil_fnl%initialize(coil=self%coil)
!<```

! PRISM FNL classes
use :: adam_prism_fnl_coil_object, only : prism_fnl_coil_object

implicit none
private
public :: coil_fnl

type(prism_fnl_coil_object), target :: coil_fnl !< Program-scope coil FNL singleton.
endmodule adam_prism_fnl_coil_global
