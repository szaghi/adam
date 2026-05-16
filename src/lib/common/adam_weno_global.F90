!< ADAM, global WENO singleton — pointer shim into equation%weno (Step 2 of forest-of-trees migration, issue #10).
module adam_weno_global
!< ADAM, global WENO singleton.
!<
!< Historically this module owned a `type(weno_object) :: weno` instance.
!< Step 2 of the forest-of-trees migration moves weno/ib/rk inside the
!< `realm_object` as value components, so the true storage lives at
!< `equation%weno` (the running solver instance — `prism` for the PRISM
!< app, etc.). This module is now a compatibility shim: `weno` is a
!< pointer aliased into `<equation>%weno` by inline pointer assignment
!< in realm_object%initialize.
!<
!< The shim deliberately does NOT `use adam_realm_object` — that would
!< create a circular dependency, because the equation module's transitive
!< object dependencies (ib_object, rk_object, etc.) already import the
!< sibling shims. Unlike the Step 1 grid/field/tree/maps shims (which use
!< the standalone `adam_adam_bind` module), the equation-side binding
!< cannot live in a separate module — a binder taking `class(realm_object)`
!< would force exactly the cycle this comment guards against. The binding
!< is therefore inlined at the only call site that needs it.

! ADAM classes, libraries, parameters
use :: adam_weno_object, only : weno_object

implicit none
private
public :: weno

type(weno_object), pointer :: weno => null() !< Program-scope WENO shim, bound inside realm_object%initialize.
endmodule adam_weno_global
