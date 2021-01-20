!< ADAM, test ADAM class in MPI env.
program adam_test_adam_object_mpi_cuda
!< ADAM, test ADAM class in MPI env.

use adam_adam_object
use adam_equation_convect1D_cpu_object
use adam_equation_laplace_cpu_object
use adam_parameters
use PENF

implicit none

type(adam_object)                   :: adam            !< ADAM.
type(equation_convect1D_cpu_object) :: convect         !< 1D convection equation.
type(equation_laplace_cpu_object)   :: laplace         !< Laplace equation.
type(tree_node_object), pointer     :: node            !< Tree node pointer.
integer(I4P)                        :: l, t            !< Counter.
logical                             :: is_grid_changed !< Flag to check grid changes.
real(R8P)                           :: time            !< Time.
integer(I8P)                        :: timing(0:2)     !< Tic toc timing.
integer(I4P)                        :: n_iter          !< Number of iterations.
integer(I4P)                        :: n_save          !< Frequency of saving output.
!GPUtype(field_gpu_object) :: field_gpu       !< GPU field.

print '(A)', 'initialize ADAM'
call adam%initialize(max_level=7,                                              &
                     emin=[0._R8P,0._R8P,0._R8P], emax=[1._R8P,1._R8P,1._R8P], &
                     ni=8_I4P, nj=8_I4P, nk=8_I4P, gc=[2_I4P,2_I4P,2_I4P],     &
                     ! bc_type=[BC_INFLOW,       BC_EXTRAPOLATION,               &
                     ! bc_type=[BC_EXTRAPOLATION,BC_EXTRAPOLATION,               &
                     bc_type=[BC_PERIODIC,BC_PERIODIC,                         &
                              BC_EXTRAPOLATION,BC_EXTRAPOLATION,               &
                              BC_EXTRAPOLATION,BC_EXTRAPOLATION],              &
                     nb=40000, nodes_number=16*40000_I8P)

call laplace%initialize(field=adam%field)
!GPUcall field_gpu%initialize(field_cpu=adam%field)

do l=1, 2
   print '(A)', 'refine ADAM at level '//trim(str(l))
   print *, 'blocks_number: ',adam%field%blocks_number
   call adam%tree%mark_all_nodes(mark=TO_BE_REFINED)
   call laplace%update_ghost(q=adam%field%q)
   call adam%amr_update(do_blocks_reorder=.false., do_mpi_redistribute=.true.)
enddo

! node => adam%tree%node(code=0_I8P)
! node%refinement_needed = TO_BE_REFINED
! call adam%amr_update(do_blocks_reorder=.false.)
! print*,'n_blocks: ',adam%tree%nodes_number

print*,' BC faces number: ', size(adam%tree%local_map_bc_face, dim=1)
print*,' BC edges number: ', size(adam%tree%local_map_bc_edge, dim=1)
print*,' BC corners number: ', size(adam%tree%local_map_bc_corner, dim=1)

print '(A)', 'set initial conditions'
call laplace%set_initial_conditions

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

time = 0._R8P
n_iter = 100
n_save = 1
call system_clock(timing(1))
do t=1, n_iter
   if (mod(t,1)==0) print '(A)', 'track iteration '//trim(str(t, .true.))
   !GPUcall field_gpu%rk_integrate(t=time, Dt=0.1_R8P)

   call laplace%mark_by_grad_q
   call laplace%update_ghost(q=adam%field%q)
   call adam%amr_update(is_marked_by_field=.true., do_blocks_reorder=.false.)
   print*, 'blocks number: ', adam%tree%nodes_number

   call laplace%integrate(t=time, Dt=0.006_R8P)
   !call field_gpu%copy_gpu_cpu
   if (mod(t,n_save)==0) call adam%save_hdf5(basename='sphere-'//trim(strz(t,9)), with_ghost=.false., with_cell_morton=.true.)
   time = time + 0.2_R8P
enddo
call system_clock(timing(2), timing(0))
print '(A, F8.3)', 'timing: ', real(timing(2) - timing(1))/ timing(0) / n_iter
!GPUcall field_gpu%copy_gpu_cpu

print '(A)', 'finalize ADAM'
call adam%finalize

endprogram adam_test_adam_object_mpi_cuda
