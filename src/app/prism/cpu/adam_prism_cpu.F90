!< ADAM, Maxwell application solver, CPU backend.
program adam_prism_cpu
!< ADAM, Maxwell application solver, CPU backend.

use adam_prism_cpu_object, only : prism_cpu_object

implicit none

type(prism_cpu_object) :: prism           !< Navier-Stokes equations system class.
integer                :: na              !< Number of command line arguments.
character(999)         :: input_file_name !< Input file name.

na = command_argument_count()
if (na==0) then
   call prism%simulate(filename='input.ini')
else
   call get_command_argument(1, input_file_name)
   call prism%simulate(filename=trim(adjustl(input_file_name)))
endif
endprogram adam_prism_cpu