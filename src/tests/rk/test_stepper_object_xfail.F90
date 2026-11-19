!< Expected-failure test of the SELAH adapter: an unknown `[runge_kutta]` scheme stops the run (ADAM #132 box 1.2).
program test_stepper_object_xfail
!< Expected-failure test of the SELAH adapter: an unknown `[runge_kutta]` scheme stops the run (ADAM #132 box 1.2).
!<
!< `stepper_object%initialize` with the scheme of the first argument (default `runge-kutta-4`, a name ADAM never had)
!< must stop with exit code 1 and the message "unknown Runge-Kutta scheme" on stderr, before allocating anything. Run:
!< `exe/test_stepper_object_xfail; echo $?` prints the message and 1. A known scheme as argument (e.g. runge-kutta-ssp-33)
!< initializes and exits 0, the control of the test.
use adam_field_object,   only : field_object
use adam_grid_object,    only : grid_object
use adam_mpih_global,    only : mpih
use adam_stepper_object, only : stepper_object
use penf,                only : I4P

implicit none
type(field_object)   :: field  !< Field.
type(grid_object)    :: grid   !< Grid.
type(stepper_object) :: ts     !< Stepper.
character(64)        :: scheme !< Scheme name.

mpih%myrankstr = ''
scheme = 'runge-kutta-4'
if (command_argument_count() > 0) call get_command_argument(1, scheme)
field%nv = 3 ; field%nb = 2 ; field%blocks_number = 1
grid%ni = 4 ; grid%nj = 4 ; grid%nk = 4 ; grid%ngc = 2
call ts%initialize(field=field, grid=grid, scheme=trim(scheme))
print '(A)', 'test_stepper_object_xfail: "'//trim(scheme)//'" initialized'
call ts%destroy
endprogram test_stepper_object_xfail
