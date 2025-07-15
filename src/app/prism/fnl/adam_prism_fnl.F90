!< ADAM, Maxwell application solver, GPU (FNL) backend.
program adam_prism_fnl
!< ADAM, Maxwell application solver, GPU (FNL) backend.

! PRSIM modules
use adam_prism_fnl_object, only : prism_fnl_object

implicit none

type(prism_fnl_object) :: prism           !< Navier-Stokes equations system class.
integer                :: na              !< Number of command line arguments.
character(999)         :: input_file_name !< Input file name.

na = command_argument_count()
if (na==0) then
   call prism%simulate(filename='input.ini')
else
   call get_command_argument(1, input_file_name)
   call prism%simulate(filename=trim(adjustl(input_file_name)))
endif
endprogram adam_prism_fnl
