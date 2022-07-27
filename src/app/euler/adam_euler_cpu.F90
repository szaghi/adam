!< ADAM, Euler solver, CPU backend.
program adam_euler_cpu
!< ADAM, Euler solver, CPU backend.

use adam_equation_euler_cpu_object, only : equation_euler_cpu_object

implicit none
type(equation_euler_cpu_object) :: euler !< Euler equations solver.

! call euler%solve(filename='adam_euler_sphere_shock.ini')

endprogram adam_euler_cpu
