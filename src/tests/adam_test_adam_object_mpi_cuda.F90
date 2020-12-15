!< ADAM, test ADAM class in MPI env.
program adam_test_adam_object_mpi
!< ADAM, test ADAM class in MPI env.

use adam_adam_object
use adam_parameters
use adam_field_gpu_object
use PENF

implicit none

type(adam_object)      :: adam            !< ADAM.
type(field_gpu_object) :: field_gpu       !< GPU field.
integer(I4P)           :: l, t            !< Counter.
logical                :: is_grid_changed !< Flag to check grid changes.
real(R8P)              :: time            !< Time.
integer(I8P)           :: timing(0:2)     !< Tic toc timing.
integer(I4P)           :: n_iter          !< Number of iterations

print '(A)', 'initialize ADAM'
call adam%initialize(max_level=6, emin=[0._R8P,0._R8P,0._R8P], emax=[1._R8P,1._R8P,1._R8P], nb=40000)
call field_gpu%initialize(field_cpu=adam%field)

do l=1, 10
   print '(A)', 'refine ADAM at level '//trim(str(l))
   call adam%tree%mark_all_nodes(mark=TO_BE_REFINED)
   call adam%amr_update(do_blocks_reorder=.false.)
enddo

print*,'n_blocks: ',adam%tree%nodes_number

print '(A)', 'set initial conditions'
call adam%field%set_initial_conditions
call adam%save_hdf5(basename='sphere-'//trim(strz(0,9)))

call field_gpu%copy_cpu_gpu
time = 0._R8P
n_iter = 25
call system_clock(timing(1))
do t=1, n_iter
   if (mod(t,1)==0) print '(A)', 'track iteration '//trim(str(t, .true.))
   call field_gpu%rk_integrate(t=time, Dt=0.1_R8P)
   call field_gpu%copy_gpu_cpu
   !RIMETTEREcall adam%save_hdf5(basename='sphere-'//trim(strz(t,9)))
   time = time + 0.1_R8P
enddo
call system_clock(timing(2), timing(0))
print '(A, F8.3)', 'timing: ', real(timing(2) - timing(1))/ timing(0) / n_iter

print '(A)', 'finalize ADAM'
call adam%finalize

endprogram adam_test_adam_object_mpi
