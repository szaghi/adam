!< ADAM, ADAM for Euler equation, GPU backend.
program adam_euler_gpu
!< ADAM, ADAM for Euler equation, GPU backend.

use adam_adam_object, only : adam_object
use adam_equation_euler_gpu_object, only : equation_euler_gpu_object, BC_EXTRAPOLATION, BC_INFLOW
use adam_parameters
use FINER, only : file_ini
use PENF
use CUDAFOR
use MPI

implicit none

type(adam_object)               :: adam            !< ADAM.
type(equation_euler_gpu_object) :: euler           !< Euler equations system.
type(file_ini)                  :: file_parameters !< Simulation parameters file handler.
integer(I4P)                    :: t               !< Counter.
real(R8P)                       :: time            !< Time.
real(R8P)                       :: time_max        !< Maximum time of integration.
integer(I4P)                    :: t_max           !< Maximum time iteration.
integer(I4P)                    :: n_save          !< Frequency of saving output.
character(999)                  :: output_basename !< Output base name.
real(R8P)                       :: timing(1:2)     !< Tic toc timing.

call initialize(filename='src/app/euler/adam_euler.ini')

print '(A)', 'initial status'

call euler%adam%print_status

! call euler%refine_uniform(refinement_levels=euler%adam%tree%iu_ref_levels, do_blocks_reorder=.false.)

! call euler%adam%prune(ijkl_prune=euler%adam%tree%ijkl_prune, do_blocks_reorder=.false.)

call euler%set_initial_conditions(file_parameters=file_parameters)

call euler%save_hdf5(output_basename=output_basename, t=0, time=0._R8P)

print '(A)', 'refined/pruned status'

call euler%adam%print_status

time = 0._R8P
t = 0
call MPI_BARRIER(MPI_COMM_WORLD, adam%error) ; timing(1) = MPI_Wtime()
integration: do
   t = t + 1
   ! call euler%amr_update(iterations=1)
   call euler%compute_dt
   if (time + euler%dt > time_max) euler%dt = time_max - time
   if (mod(t,1)==0.and.adam%myrank==0) call euler%print_progress(t=t, time=time, time_max=time_max)
   call euler%integrate(t=time)
   if (mod(t,n_save)==0) call euler%save_hdf5(output_basename=output_basename, t=t, time=time)
   time = time + euler%dt
   if (time>=time_max.or.t>=t_max) exit integration
enddo integration
call MPI_BARRIER(MPI_COMM_WORLD, adam%error) ; timing(2) = MPI_Wtime()
print '(A, F18.10)', 'timing: ', (timing(2) - timing(1))/t

call euler%save_hdf5(output_basename=output_basename, t=t, time=time)

call adam%finalize
contains
   subroutine initialize(filename)
   !< Parse parameters file getting simulation input data.
   character(*), intent(in)  :: filename             !< Parameters file name.
   type(file_ini)            :: file_adam_euler      !< Adam Euler input file handler.
   character(999)            :: file_parameters_name !< Name of file parameters.

   call file_adam_euler%initialize(filename=trim(filename))
   call file_adam_euler%load
   call file_adam_euler%get(section_name='adam_euler', option_name='file_parameters_name', val=file_parameters_name)

   call file_parameters%initialize(filename=trim(file_parameters_name))
   call file_parameters%load

   call adam%initialize(file_parameters=file_parameters)
   call adam%print_status
   call euler%initialize(adam=adam, file_parameters=file_parameters)

   call file_parameters%get(section_name='simulation', option_name='time_max'       , val=time_max       )
   call file_parameters%get(section_name='simulation', option_name='t_max'          , val=t_max          )
   call file_parameters%get(section_name='simulation', option_name='output_basename', val=output_basename)
   call file_parameters%get(section_name='simulation', option_name='n_save'         , val=n_save         )
   endsubroutine initialize
endprogram adam_euler_gpu
