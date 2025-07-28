!< ADAM, PATCH Poisson solver with Adpative mesh Refinement for HPC computing, CPU backend.
program adam_patch_cpu
!< ADAM, PATCH Poisson solver with Adpative mesh Refinement for HPC computing, CPU backend.

use adam_patch_cpu_object, only : patch_cpu_object

implicit none

type(patch_cpu_object) :: patch           !< Navier-Stokes equations system class.
integer                :: na              !< Number of command line arguments.
character(999)         :: input_file_name !< Input file name.

na = command_argument_count()
if (na==0) then
   call patch%simulate(filename='input.ini')
else
   call get_command_argument(1, input_file_name)
   call patch%simulate(filename=trim(adjustl(input_file_name)))
endif
endprogram adam_patch_cpu
