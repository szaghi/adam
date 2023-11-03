!< ADAM, Navier-Stokes application solver, GPU (NVF) backend.
program adam_mhd_nvf
!< ADAM, Navier-Stokes application solver, GPU (NVF) backend.

use adam_mhd_nvf_object, only : mhd_nvf_object

implicit none

type(mhd_nvf_object) :: mhd           !< Navier-Stokes equations system class.
integer                :: na              !< Number of command line arguments.
character(999)         :: input_file_name !< Input file name.

na = command_argument_count()
if (na==0) then
   call mhd%simulate(filename='input.ini')
else
   call get_command_argument(1, input_file_name)
   call mhd%simulate(filename=trim(adjustl(input_file_name)))
endif
endprogram adam_mhd_nvf
