!< ADAM PRISM FNL library, entry for all PRISM FNL (and common) classes and libraries.
module adam_prism_fnl_library
!< ADAM PRISM FNL library, entry for all PRISM FNL (and common) classes and libraries.

! ADAM FNL classes, libraries, parameters
use :: adam_fnl_library
! PRISM common classes, libraries, parameters
use :: adam_prism_common_library
! PRISM FNL classes, libraries, parameters
use :: adam_prism_fnl_coil_object
use :: adam_prism_fnl_external_fields_kernels
use :: adam_prism_fnl_fwlayer_object
! PRISM FNL singleton objects
use :: adam_prism_fnl_globals

implicit none
public

endmodule adam_prism_fnl_library
