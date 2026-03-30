!< ADAM, PRISM program-scope singletons — convenience re-export of all adam_prism_*_global modules.
module adam_prism_globals
!< ADAM, PRISM program-scope singletons — convenience re-export of all adam_prism_*_global modules.
!<
!< Provides a single USE point for all PRISM program-scope singleton objects:
!<```fortran
!< use :: adam_prism_globals, only : bc, coil, external_fields, fWLayer, ic, leapfrog_pic, &
!<                                   numerics, particle_injection, physics, pic, rk_bc,     &
!<                                   rk_pic, time
!<```
!< instead of listing every individual adam_prism_*_global module.

! PRISM global singletons
use :: adam_prism_bc_global,                 only : bc
use :: adam_prism_coil_global,               only : coil
use :: adam_prism_external_fields_global,    only : external_fields
use :: adam_prism_fWLayer_global,            only : fWLayer
use :: adam_prism_ic_global,                 only : ic
use :: adam_prism_leapfrog_pic_global,       only : leapfrog_pic
use :: adam_prism_numerics_global,           only : numerics
use :: adam_prism_particle_injection_global, only : particle_injection
use :: adam_prism_physics_global,            only : physics
use :: adam_prism_pic_global,                only : pic
use :: adam_prism_rk_bc_global,              only : rk_bc
use :: adam_prism_rk_pic_global,             only : rk_pic
use :: adam_prism_time_global,               only : time

implicit none
private
public :: bc
public :: coil
public :: external_fields
public :: fWLayer
public :: ic
public :: leapfrog_pic
public :: numerics
public :: particle_injection
public :: physics
public :: pic
public :: rk_bc
public :: rk_pic
public :: time
endmodule adam_prism_globals
