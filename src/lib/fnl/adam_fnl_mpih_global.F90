!< ADAM, global FNL MPI handler shim — pointer into prism_fnl_object's value component (PRISM FNL C.3 closure, issue #13).
module adam_fnl_mpih_global
!< ADAM, global FNL MPI handler shim — pointer into the current realm's mpih_fnl value component.
!<
!< Was a singleton owner (`type(mpih_fnl_object), target :: mpih_fnl`) until the
!< PRISM FNL C.3 closure (issue #13 D.4a follow-up — unblocks the per-substage
!< error_stop guard in prism_fnl_object). For N=1 the shim aliases
!< `realm(1)%mpih_fnl`; reads through `mpih_fnl%...` see the same data the
!< singleton used to own. For N>1 (rmf-2realm) the shim is rebound per realm
!< by `prism_fnl_object%bind_my_globals_forest` before each per-substage TBP.
!<
!< This shim is shared with NASTO-FNL by code path but NASTO is out of scope
!< for issue #13 — see #10/#13 PRISM-only discipline.

! ADAM FNL classes
use :: adam_fnl_mpih_object, only : mpih_fnl_object

implicit none
private
public :: mpih_fnl

type(mpih_fnl_object), pointer :: mpih_fnl => null() !< Program-scope mpih_fnl shim, bound by prism_fnl_object.
endmodule adam_fnl_mpih_global
