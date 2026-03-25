!< ADAM common library, entry for all common classes and libraries.
module adam_common_library
!< ADAM common library, entry for all common classes and libraries.

! ADAM classes, libraries, parameters
use :: adam_adam_object
use :: adam_amr_object
use :: adam_blanes_moan_object
use :: adam_cfm_object
use :: adam_eos_ic_object
use :: adam_equation_object
use :: adam_fdv_operators_library
use :: adam_field_object
use :: adam_flail_object
use :: adam_grid_object
use :: adam_ib_object
use :: adam_io_object
use :: adam_leapfrog_object
use :: adam_maps_object
use :: adam_mpih_object
use :: adam_parameters
use :: adam_refinement_plan_object
use :: adam_rk_object
use :: adam_riemann_euler_library
use :: adam_slices_object
use :: adam_tree_node_object
use :: adam_tree_bucket_object
use :: adam_tree_object
use :: adam_weno_object
! ADAM singleton objects
use :: adam_mpih_global
use :: adam_grid_global
use :: adam_field_global
use :: adam_maps_global

implicit none
public

endmodule adam_common_library
