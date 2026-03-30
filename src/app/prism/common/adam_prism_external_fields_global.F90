!< ADAM, global PRISM external fields singleton — single program-scope prism_external_fields_object instance.
module adam_prism_external_fields_global
!< ADAM, global PRISM external fields singleton — single program-scope prism_external_fields_object instance.

! PRISM modules
use :: adam_prism_external_fields_object, only : prism_external_fields_object

implicit none
private
public :: external_fields

type(prism_external_fields_object), target :: external_fields !< Program-scope PRISM external fields singleton.
endmodule adam_prism_external_fields_global
