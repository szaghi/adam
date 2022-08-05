!< ADAM, Euler solver, CPU backend.
program adam_euler_cpu
!< ADAM, Euler solver, CPU backend.

use adam_equation_euler_cpu_object, only : equation_euler_cpu_object

implicit none
type(equation_euler_cpu_object) :: euler           !< Euler equations solver.
character(99)                   :: file_parameters !< Input file of simulation parameters.

call get_command_argument(1,file_parameters)
call euler%solve(filename=trim(adjustl(file_parameters)))
endprogram adam_euler_cpu
