!< ADAM, global field singleton — pointer shim into adam%field (Step 1 of forest-of-trees migration, issue #10).
module adam_field_global
!< ADAM, global field singleton.
!<
!< Compatibility shim: storage lives in `adam%field` (the `adam_adam_global`
!< singleton). See `adam_grid_global` for the rationale and the dependency
!< note. Bound by `adam_adam_bind`.

! ADAM classes, libraries, parameters
use :: adam_field_object, only : field_object

implicit none
private
public :: field

type(field_object), pointer :: field => null() !< Program-scope field shim, bound by adam_adam_bind.
endmodule adam_field_global
