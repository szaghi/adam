!< RIPPLE, scalar linear advection on the ADAM library, FNL (OpenACC) backend.
program adam_ripple_fnl
!< RIPPLE, scalar linear advection on the ADAM library, FNL (OpenACC) backend.
!<
!< Usage: `adam_ripple_fnl [input.ini]`. The CPU driver with the FNL realm: the forest does not know the backend.

! ADAM classes, libraries, parameters
use :: adam_forest_object,     only : forest_object
use :: adam_realm_object,      only : realm_pointer
! RIPPLE modules
use :: adam_ripple_fnl_object, only : ripple_fnl_object

implicit none

type(ripple_fnl_object), target :: realm(1)          !< The realm.
type(realm_pointer)             :: realm_pointers(1) !< A pointer to it: what the forest takes (issue #96).
type(forest_object)             :: forest            !< Orchestrator.
character(999)                  :: input_file_name   !< Input file name.

input_file_name = 'input.ini'
if (command_argument_count() > 0) call get_command_argument(1, input_file_name)
realm_pointers(1)%ptr => realm(1)
call forest%simulate(realm=realm_pointers, filename=trim(adjustl(input_file_name)))
endprogram adam_ripple_fnl
