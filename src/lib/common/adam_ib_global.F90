!< ADAM, global IB singleton — pointer shim into equation%ib (Step 2 of forest-of-trees migration, issue #10).
module adam_ib_global
!< ADAM, global IB singleton.
!<
!< Compatibility shim: storage lives in `equation%ib` (the running solver
!< instance). See `adam_weno_global` for the rationale and the dependency
!< note. Bound by inline pointer assignment in equation_object%initialize.

! ADAM classes, libraries, parameters
use :: adam_ib_object, only : ib_object

implicit none
private
public :: ib

type(ib_object), pointer :: ib => null() !< Program-scope IB shim, bound inside equation_object%initialize.
endmodule adam_ib_global
