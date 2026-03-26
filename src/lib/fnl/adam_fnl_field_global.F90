!< ADAM, global FNL field singleton — single program-scope field_fnl_object instance.
module adam_fnl_field_global
!< ADAM, global FNL field singleton — single program-scope field_fnl_object instance.
!<
!< Provides the GPU field handler (FNL backend) as a program-scope singleton,
!< mirroring the CPU-side adam_field_global pattern.
!<
!< Requires explicit initialization before any FNL kernel is invoked:
!<```fortran
!< call field_fnl%initialize(verbose=.true.)
!<```

! ADAM FNL classes
use :: adam_fnl_field_object, only : field_fnl_object

implicit none
private
public :: field_fnl

type(field_fnl_object), target :: field_fnl !< Program-scope field FNL singleton.
endmodule adam_fnl_field_global
