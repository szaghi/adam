!< ADAM, global field singleton — single program-scope field_object instance.
module adam_global_field
!< ADAM, global field singleton — single program-scope field_object instance.
use :: adam_field_object, only: field_object
implicit none
private
public :: field
type(field_object), target :: field  !< Program-scope field singleton.
endmodule adam_global_field
