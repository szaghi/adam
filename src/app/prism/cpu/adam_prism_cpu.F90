!< ADAM, Maxwell application solver, CPU backend.
program adam_prism_cpu
!< ADAM, Maxwell application solver, CPU backend.
!<
!< Phase C of issue #10: the driver is now organized as a **forest of
!< realms**. For v1 the forest holds a single realm (N=1) and behaves
!< identically to the legacy `prism%simulate(filename)` entry point;
!< Phase D will land N>1 once a concrete multi-realm use case is
!< specified.
!<
!< The legacy seven `adam_*_global` shims (`grid`, `field`, ...) continue
!< to alias correctly for N=1: the last call to `realm(is)%initialize_forest`
!< inside `forest%initialize` lands on `realm(N)`, which for N=1 is
!< `realm(1)`. Phase D will need an explicit rebind-to-realm(1) after the
!< per-realm init loop; with N=1 the implicit bind suffices.

! ADAM common library — forest orchestrator
use :: adam_forest_object,    only : forest_object
! PRISM modules
use :: adam_prism_cpu_object, only : prism_cpu_object

implicit none

type(prism_cpu_object) :: realm(1)        !< Single-realm array (v1; Phase D extends to N>1).
type(forest_object)    :: forest          !< Orchestrator that drives the realm array.
integer                :: na              !< Number of command line arguments.
character(999)         :: input_file_name !< Input file name.

na = command_argument_count()
if (na==0) then
   call forest%simulate(realm=realm, filename='input.ini')
else
   call get_command_argument(1, input_file_name)
   call forest%simulate(realm=realm, filename=trim(adjustl(input_file_name)))
endif
endprogram adam_prism_cpu
