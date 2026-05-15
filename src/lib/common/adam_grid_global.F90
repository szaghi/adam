!< ADAM, global grid singleton — pointer shim into adam%grid (Step 1 of forest-of-trees migration, issue #10).
module adam_grid_global
!< ADAM, global grid singleton.
!<
!< Historically this module owned a `type(grid_object) :: grid` instance.
!< Step 1 of the forest-of-trees migration moves grid/tree/field/maps inside
!< `adam_object` as value components, so the true storage lives at
!< `adam%grid` (the `adam_adam_global` singleton). This module is now a
!< compatibility shim: `grid` is a pointer aliased into `adam%grid` by
!< `adam_adam_bind`, which must be called once after `adam%initialize`.
!<
!< The 28 consumer files that `use adam_grid_global, only: grid` see no
!< API change — `grid%ngc`, `grid%ni`, ... continue to work because the
!< pointer transparently dereferences. Specification expressions like
!< `q(1:, 1-grid%ngc:, ...)` remain valid as long as the binding has run
!< before the procedure is called.
!<
!< The shim deliberately does NOT `use adam_adam_object` — that would create
!< a circular dependency, because `adam_field_object`/`adam_tree_object`/
!< `adam_maps_object` import `grid` here, and `adam_adam_object` imports
!< those. The binder lives in a separate module `adam_adam_bind`.

! ADAM classes, libraries, parameters
use :: adam_grid_object, only : grid_object

implicit none
private
public :: grid

type(grid_object), pointer :: grid => null() !< Program-scope grid shim, bound by adam_adam_bind.
endmodule adam_grid_global
