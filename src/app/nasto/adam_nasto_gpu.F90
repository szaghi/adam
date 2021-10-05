!< ADAM, ADAM for Navier-Stokes equation, GPU backend.
program adam_nasto_gpu
!< ADAM, ADAM for Navier-Stokes equation, GPU backend.

use adam_equation_nasto_gpu_object
use MPI

implicit none
type(equation_nasto_gpu_object) :: nasto !< Navier-Stokes equations system.

! call nasto%run(filename='adam_nasto_sphere_shock.ini')
call MPI_INIT(nasto%error)
call nasto%run(filename='adam_nasto_sphere_shock-work-small.ini')

endprogram adam_nasto_gpu
