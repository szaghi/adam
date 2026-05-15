!< ADAM, global RK singleton — pointer shim into equation%rk (Step 2 of forest-of-trees migration, issue #10).
module adam_rk_global
!< ADAM, global RK singleton.
!<
!< Compatibility shim: storage lives in `equation%rk` (the running solver
!< instance). See `adam_weno_global` for the rationale and the dependency
!< note. Bound by inline pointer assignment in equation_object%initialize.

! ADAM classes, libraries, parameters
use :: adam_rk_object, only : rk_object

implicit none
private
public :: rk

type(rk_object), pointer :: rk => null() !< Program-scope RK shim, bound inside equation_object%initialize.
endmodule adam_rk_global
