!< ADAM, global maps singleton — single program-scope maps_object instance.
module adam_global_maps
!< ADAM, global maps singleton — single program-scope maps_object instance.

! ADAM classes, libraries, parameters
use :: adam_maps_object, only : maps_object

implicit none
private
public :: maps

type(maps_object), target :: maps !< Program-scope maps singleton.
endmodule adam_global_maps
