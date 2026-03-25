!< ADAM FNL library, entry for all FNL (and common) classes and libraries.
module adam_fnl_library
!< ADAM FNL library, entry for all FNL (and common) classes and libraries.

! ADAM modules
use :: adam_common_library
use :: adam_global_mpih_fnl
use :: adam_fnl_fdv_operators_library
use :: adam_fnl_field_kernels
use :: adam_fnl_field_object
use :: adam_fnl_ib_kernels
use :: adam_fnl_ib_object
use :: adam_fnl_maps_object
use :: adam_fnl_mpih_object
use :: adam_fnl_rk_kernels
use :: adam_fnl_rk_object
use :: adam_fnl_weno_kernels
use :: adam_fnl_weno_object

implicit none
public

endmodule adam_fnl_library
