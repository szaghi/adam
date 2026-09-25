!< ADAM, FLUME common library: re-export of every backend-independent FLUME module.
module adam_flume_common_library
!< ADAM, FLUME common library: re-export of every backend-independent FLUME module.

! FLUME modules
use :: adam_flume_bc_object
use :: adam_flume_common_object
use :: adam_flume_diagnostics_object
use :: adam_flume_euler_library
use :: adam_flume_ic_object
use :: adam_flume_mhd_library
use :: adam_flume_mhd_object
use :: adam_flume_numerics_object
use :: adam_flume_parameters
use :: adam_flume_physics_object
use :: adam_flume_time_object

implicit none
public
endmodule adam_flume_common_library
