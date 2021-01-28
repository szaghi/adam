!< ADAM, ADAM for convective 1D equation, GPU backend.
program adam_convect1D_gpu
!< ADAM, ADAM for convective 1D equation, GPU backend.

use adam_adam_object
use adam_equation_convect1D_gpu_object
use adam_parameters
use PENF
use MPI

implicit none

type(adam_object)                   :: adam            !< ADAM.
type(equation_convect1D_gpu_object) :: convect_gpu     !< 1D convection equation, GPU backend.
integer(I4P)                        :: l, t            !< Counter.
logical                             :: is_grid_changed !< Flag to check grid changes.
real(R8P)                           :: time            !< Time.
integer(I4P)                        :: n_iter=100      !< Number of iterations.
integer(I4P)                        :: n_save=150       !< Frequency of saving output.
integer(I8P)                        :: timing(0:2)     !< Tic toc timing.

print '(A)', '1D convection equation integration'
call adam%initialize(max_level=9,                                              &
                     emin=[0._R8P,0._R8P,0._R8P], emax=[1._R8P,1._R8P,1._R8P], &
                     ni=16_I4P, nj=16_I4P, nk=16_I4P, gc=[2_I4P,2_I4P,2_I4P],     &
                     bc_type=[BC_PERIODIC,BC_PERIODIC,                         &
                              BC_EXTRAPOLATION,BC_EXTRAPOLATION,               &
                              BC_EXTRAPOLATION,BC_EXTRAPOLATION],              &
                     nb=40000, nodes_number=16*15000_I8P)

call convect_gpu%initialize(field=adam%field)

print '(A)', 'create 2 levels of refinement'
do l=1, 2
   print '(A)', 'refine ADAM at level '//trim(str(l))
   print *, 'blocks_number: ',adam%field%blocks_number
   call adam%tree%mark_all_nodes(mark=TO_BE_REFINED)
   call adam%amr_update(do_blocks_reorder=.false., do_mpi_redistribute=.true.)
enddo
print '(A)', 'set initial conditions'
call convect_gpu%set_initial_conditions
call adam%save_hdf5(basename='convect-'//trim(strz(0,9)), with_ghost=.false., with_cell_morton=.true.)
call convect_gpu%copy_cpu_gpu

do l=1,5 ! n_blocks = 8464
!do l=1,6 ! n_blocks = 53824
   call convect_gpu%mark_by_grad_q
   call convect_gpu%update_ghost_gpu(q_gpu=convect_gpu%q_gpu)
   call convect_gpu%copy_gpu_cpu
   call adam%amr_update(is_marked_by_field=.true., do_blocks_reorder=.false.)
   call convect_gpu%copy_cpu_gpu
   if (mod(t,1)==0.and.adam%myrank==0) print '(A)', 'blocks number: '//trim(str(adam%tree%nodes_number, .true.))
enddo

time = 0._R8P
call MPI_BARRIER(MPI_COMM_WORLD, adam%error) ; call system_clock(timing(1))
do t=1, n_iter
   if (mod(t,1)==0.and.adam%myrank==0) print '(A)', 'track iteration '//trim(str(t, .true.))
   call convect_gpu%integrate(t=time, Dt=0.00006_R8P)
   if (mod(t,n_save)==0) then
      call convect_gpu%copy_gpu_cpu
      call adam%save_hdf5(basename='convect-'//trim(strz(t,9)), with_ghost=.false., with_cell_morton=.true.)
   endif
   time = time + 0.2_R8P
enddo
call MPI_BARRIER(MPI_COMM_WORLD, adam%error) ; call system_clock(timing(2), timing(0))
print '(A, F8.3)', 'timing: ', real(timing(2) - timing(1))/ timing(0) / n_iter

call convect_gpu%copy_gpu_cpu
call adam%save_hdf5(basename='convect-'//trim(strz(t,9)), with_ghost=.false., with_cell_morton=.true.)
print*,'maxval :',maxval(adam%field%q(:,:,:,:,1:adam%field%blocks_number))
print*,'minval :',minval(adam%field%q(:,:,:,:,1:adam%field%blocks_number))

call adam%finalize
endprogram adam_convect1D_gpu
