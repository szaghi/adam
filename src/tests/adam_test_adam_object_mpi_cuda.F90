!< ADAM, test ADAM class in MPI env.
program adam_test_adam_object_mpi
!< ADAM, test ADAM class in MPI env.

use adam_adam_object
use adam_parameters
!GPUuse adam_field_gpu_object
use adam_tree_node_object
use PENF

implicit none

type(adam_object)      :: adam            !< ADAM.
!GPUtype(field_gpu_object) :: field_gpu       !< GPU field.
integer(I4P)           :: l, t            !< Counter.
logical                :: is_grid_changed !< Flag to check grid changes.
real(R8P)              :: time            !< Time.
integer(I8P)           :: timing(0:2)     !< Tic toc timing.
integer(I4P)           :: n_iter          !< Number of iterations
type(tree_node_object), pointer :: node

print '(A)', 'initialize ADAM'
call adam%initialize(max_level=7, emin=[0._R8P,0._R8P,0._R8P], emax=[1._R8P,1._R8P,1._R8P], nb=40000, nodes_number=16*40000_I8P)
!GPUcall field_gpu%initialize(field_cpu=adam%field)

do l=1, 6
   print '(A)', 'refine ADAM at level '//trim(str(l))
   print *, 'blocks_number: ',adam%field%blocks_number
   call adam%tree%mark_all_nodes(mark=TO_BE_REFINED)
   call adam%amr_update(do_blocks_reorder=.false., do_mpi_redistribute=.true. )
enddo

!node => adam%tree%node(code=15_I8P)
!node%refinement_needed = TO_BE_REFINED
!call adam%amr_update(do_blocks_reorder=.false.)
print*,'n_blocks: ',adam%tree%nodes_number

print '(A)', 'set initial conditions'
call adam%field%set_initial_conditions
!call adam%save_hdf5(basename='sphere-'//trim(strz(0,9)), with_ghost=.false., with_cell_morton=.true.)
!call adam%save_vtk(basename='sphere-'//trim(strz(0,9)), with_ghost=.true.)

!call field_gpu%copy_cpu_gpu
!call update_ghost_gpu_u(local_map_ghost_gpu=field_gpu%local_map_ghost_gpu, u_s=field_gpu%u_gpu)
!call update_ghost_mpi_gpu_u(comm_map_recv_ptr_ghost=field_gpu%comm_map_recv_ptr_ghost, &
!                            comm_map_send_ptr_ghost=field_gpu%comm_map_send_ptr_ghost, &
!                            comm_map_recv_ghost_gpu=field_gpu%comm_map_recv_ghost_gpu, &
!                            comm_map_send_ghost_gpu=field_gpu%comm_map_send_ghost_gpu, &
!                            recv_buffer_ghost_gpu=field_gpu%recv_buffer_ghost_gpu,     &
!                            send_buffer_ghost_gpu=field_gpu%send_buffer_ghost_gpu,     &
!                            u_s=field_gpu%u_gpu, procs_number=adam%procs_number)
!GPUcall field_gpu%copy_gpu_cpu
!call adam%save_hdf5(basename='sphere-'//trim(strz(1,9)), with_ghost=.false., with_cell_morton=.true.)
!call adam%save_vtk(basename='sphere-'//trim(strz(1,9)), with_ghost=.true.)
!call adam%finalize

time = 0._R8P
n_iter = 10
call system_clock(timing(1))
do t=1, n_iter
   if (mod(t,1)==0) print '(A)', 'track iteration '//trim(str(t, .true.))
   !GPUcall field_gpu%rk_integrate(t=time, Dt=0.1_R8P)
   call adam%field%rk_integrate(t=time, Dt=0.1_R8P)
   !call field_gpu%copy_gpu_cpu
   !call adam%save_hdf5(basename='sphere-'//trim(strz(t,9)), with_ghost=.false., with_cell_morton=.true.)
   !call adam%save_vtk(basename='sphere-'//trim(strz(t,9)), with_ghost=.true.)
   time = time + 0.1_R8P
enddo
call system_clock(timing(2), timing(0))
print '(A, F8.3)', 'timing: ', real(timing(2) - timing(1))/ timing(0) / n_iter
!GPUcall field_gpu%copy_gpu_cpu

print '(A)', 'finalize ADAM'
call adam%finalize

endprogram adam_test_adam_object_mpi
