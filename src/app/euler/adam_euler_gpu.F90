!< ADAM, Euler solver, GPU backend.
program adam_euler_gpu
!< ADAM, Euler solver, GPU backend.

use adam_equation_euler_gpu_object

implicit none
type(equation_euler_gpu_object) :: euler !< Euler solver.

call euler%run(filename='adam_euler_sphere_shock.ini')

endprogram adam_euler_gpu
