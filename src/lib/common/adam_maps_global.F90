!< ADAM, global maps singleton — pointer shim into adam%maps (Step 1 of forest-of-trees migration, issue #10).
module adam_maps_global
!< ADAM, global maps singleton.
!<
!< Compatibility shim: storage lives in `adam%maps` (the `adam_adam_global`
!< singleton). See `adam_grid_global` for the rationale and the dependency
!< note. Bound by `adam_adam_bind`.

! ADAM classes, libraries, parameters
use :: adam_maps_object, only : maps_object

implicit none
private
public :: maps

type(maps_object), pointer :: maps => null() !< Program-scope maps shim, bound by adam_adam_bind.
endmodule adam_maps_global
