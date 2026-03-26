!< ADAM, PRISM FNL program-scope singletons — convenience re-export of all adam_prism_fnl_*_global modules.
module adam_prism_fnl_globals
!< ADAM, PRISM FNL program-scope singletons — convenience re-export of all adam_prism_fnl_*_global modules.
!<
!< Provides a single USE point for all PRISM FNL program-scope singleton objects:
!<```fortran
!< use :: adam_prism_fnl_globals, only : coil_fnl, fwlayer_fnl
!<```
!< instead of listing every individual adam_prism_fnl_*_global module.

! PRISM FNL global singletons
use :: adam_prism_fnl_coil_global,    only : coil_fnl
use :: adam_prism_fnl_fwlayer_global, only : fwlayer_fnl

implicit none
private
public :: coil_fnl
public :: fwlayer_fnl
endmodule adam_prism_fnl_globals
