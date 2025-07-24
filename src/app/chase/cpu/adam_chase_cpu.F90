!< ADAM, Maxwell application solver, CPU backend.
program adam_chase_cpu
!< ADAM, Maxwell application solver, CPU backend.

use adam_chase_cpu_object, only : chase_cpu_object

implicit none

type(chase_cpu_object) :: chase           !< Navier-Stokes equations system class.
integer                :: na              !< Number of command line arguments.
character(999)         :: input_file_name !< Input file name.

na = command_argument_count()
if (na==0) then
   call chase%simulate(filename='input.ini')
else
   call get_command_argument(1, input_file_name)
   call chase%simulate(filename=trim(adjustl(input_file_name)))
endif
endprogram adam_chase_cpu
