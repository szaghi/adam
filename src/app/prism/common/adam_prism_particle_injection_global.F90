!< ADAM, global PRISM particle injection singleton — single program-scope prism_particle_injection_object instance.
module adam_prism_particle_injection_global
!< ADAM, global PRISM particle injection singleton — single program-scope prism_particle_injection_object instance.

! PRISM modules
use :: adam_prism_particle_injection_object, only : prism_particle_injection_object

implicit none
private
public :: particle_injection

type(prism_particle_injection_object), target :: particle_injection !< Program-scope PRISM particle injection singleton.
endmodule adam_prism_particle_injection_global
