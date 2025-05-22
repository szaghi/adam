!< ADAM PRISM common library, entry for all PRISM common classes and libraries.
module adam_prism_common_library
!< ADAM PRISM common library, entry for all common classes and libraries.
! use adam_nasto_common_library
! use adam_nasto_fnl_cns_kernels
! use adam_nasto_fnl_kernels


use adam_prism_bc_object
use adam_prism_coil_object
use adam_prism_common_object
use adam_prism_ic_object
use adam_prism_io_object
use adam_prism_parameters
use adam_prism_physics_object
use adam_prism_time_object

! da aggiungere queste due (con la seconda da verificare trovandosi fuori da common in realtà)!!!


use adam_riemann_maxwell_library

implicit none
public

endmodule adam_prism_common_library
