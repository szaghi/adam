!< ADAM, ADAM for convective 1D equation, GPU backend.
program adam_convect1D_gpu
!< ADAM, ADAM for convective 1D equation, GPU backend.

use adam_adam_object
use adam_equation_convect1D_gpu_object
use adam_parameters
use PENF

implicit none

type(adam_object)                   :: adam            !< ADAM.
type(equation_convect1D_gpu_object) :: convect         !< 1D convection equation, GPU backend.
integer(I4P)                        :: l, t            !< Counter.
logical                             :: is_grid_changed !< Flag to check grid changes.
real(R8P)                           :: time            !< Time.
integer(I4P)                        :: n_iter=100      !< Number of iterations.
integer(I4P)                        :: n_save=1        !< Frequency of saving output.

print '(A)', '1D convection equation integration'
call adam%initialize(max_level=7,                                              &
                     emin=[0._R8P,0._R8P,0._R8P], emax=[1._R8P,1._R8P,1._R8P], &
                     ni=8_I4P, nj=8_I4P, nk=8_I4P, ngc=2_I4P,                  &
                     bc_type=[BC_PERIODIC,BC_PERIODIC,                         &
                              BC_EXTRAPOLATION,BC_EXTRAPOLATION,               &
                              BC_EXTRAPOLATION,BC_EXTRAPOLATION],              &
                     nb=5000, nodes_number=16*5000_I8P)

call convect%initialize(field=adam%field)

print '(A)', 'create 2 levels of refinement'
do l=1, 2
   call adam%tree%mark_all_nodes(mark=TO_BE_REFINED)
   call adam%amr_update(do_blocks_reorder=.false., do_mpi_redistribute=.true.)
   print '(A)', 'refine ADAM at level '//trim(str(l))
   print *, 'blocks number: ',adam%tree%nodes_number
enddo
print '(A)', 'set initial conditions'
call convect%set_initial_conditions
print '(A)', 'save initial conditions'
call convect%copy_cpu_gpu
call convect%update_ghost_gpu(q_gpu=convect%q_gpu)
call convect%copy_gpu_cpu
call adam%save_hdf5(basename='convect-'//trim(strz(0,9)), q=convect%field%q, q_name=['T'], &
                    with_cell_morton=.true.)
print '(A)', 'copy CPU to GPU'
time = 0._R8P
do t=1, n_iter
   if (mod(t,1)==0.and.adam%myrank==0) print '(A)', 'track iteration '//trim(str(t, .true.))
   call convect%mark_by_grad_q(grad_tol=9.2_R8P, delta_fine=0.004_R8P, delta_coarse=0.08_R8P)
   call convect%update_ghost_gpu(q_gpu=convect%q_gpu)
   call convect%copy_gpu_cpu
   call adam%amr_update(is_marked_by_field=.true., do_blocks_reorder=.false.)
   call convect%copy_cpu_gpu
   if (mod(t,1)==0.and.adam%myrank==0) print '(A)', 'blocks number: '//trim(str(adam%tree%nodes_number, .true.))
   call convect%integrate(t=time, Dt=0.006_R8P)
   if (mod(t,n_save)==0) then
      call convect%copy_gpu_cpu
      call adam%save_hdf5(basename='convect-'//trim(strz(t,9)), q=convect%field%q, q_name=['T'], with_cell_morton=.true.)
   endif
   time = time + 0.2_R8P
enddo

call adam%finalize
endprogram adam_convect1D_gpu
