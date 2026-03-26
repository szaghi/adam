!< ADAM, global PRISM FNL fWLayer singleton — single program-scope prism_fnl_fwlayer_object instance.
module adam_prism_fnl_fwlayer_global
!< ADAM, global PRISM FNL fWLayer singleton — single program-scope prism_fnl_fwlayer_object instance.
!<
!< Requires explicit initialization after the CPU-side fwlayer object is populated:
!<```fortran
!< call fwlayer_fnl%initialize(fwlayer=self%fwlayer)
!<```

! PRISM FNL classes
use :: adam_prism_fnl_fwlayer_object, only : prism_fnl_fwlayer_object

implicit none
private
public :: fwlayer_fnl

type(prism_fnl_fwlayer_object), target :: fwlayer_fnl !< Program-scope fWLayer FNL singleton.
endmodule adam_prism_fnl_fwlayer_global
