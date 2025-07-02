!< ADAM, Navier-Stokes application solver, GPU (FNL) backend.
program adam_nasto_fnl
!< ADAM, Navier-Stokes application solver, GPU (FNL) backend.

use adam_nasto_fnl_object, only : nasto_fnl_object

implicit none

type(nasto_fnl_object) :: nasto           !< Navier-Stokes equations system class.
integer                :: na              !< Number of command line arguments.
character(999)         :: input_file_name !< Input file name.

na = command_argument_count()
if (na==0) then
   call nasto%simulate(filename='input.ini')
else
   call get_command_argument(1, input_file_name)
   call nasto%simulate(filename=trim(adjustl(input_file_name)))
endif
endprogram adam_nasto_fnl
