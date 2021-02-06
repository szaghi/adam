!< ADAM, ADAM for Laplace equation.
program adam_laplace
!< ADAM, ADAM for Laplace equation.

use adam_adam_object
use adam_equation_laplace_cpu_object
use adam_parameters
use PENF

implicit none

type(adam_object)                 :: adam            !< ADAM.
type(equation_laplace_cpu_object) :: laplace         !< Laplace equation.
integer(I4P)                      :: l, t            !< Counter.
logical                           :: is_grid_changed !< Flag to check grid changes.
real(R8P)                         :: time            !< Time.
integer(I4P)                      :: n_iter=100      !< Number of iterations.
integer(I4P)                      :: n_save=1        !< Frequency of saving output.

print '(A)', 'Laplace equation integration'
call adam%initialize(max_level=7,                                              &
                     emin=[0._R8P,0._R8P,0._R8P], emax=[1._R8P,1._R8P,1._R8P], &
                     ni=8_I4P, nj=8_I4P, nk=8_I4P, ngc=2_I4P,                  &
                     bc_type=[BC_INFLOW,BC_EXTRAPOLATION,                      &
                              BC_EXTRAPOLATION,BC_EXTRAPOLATION,               &
                              BC_EXTRAPOLATION,BC_EXTRAPOLATION],              &
                     nb=5000, nodes_number=16*5000_I8P)
call laplace%initialize(field=adam%field)
print '(A)', 'create 2 levels of refinement'
do l=1, 2
   print '(A)', 'refine ADAM at level '//trim(str(l))
   print *, 'blocks_number: ',adam%field%blocks_number
   call adam%tree%mark_all_nodes(mark=TO_BE_REFINED)
   call laplace%update_ghost(q=adam%field%q)
   call adam%amr_update(do_blocks_reorder=.false., do_mpi_redistribute=.true.)
enddo
print '(A)', 'set initial conditions'
call laplace%set_initial_conditions
call adam%save_hdf5(basename='laplace-'//trim(strz(0,9)), q=laplace%field%q, q_name=['T'], with_cell_morton=.true.)
time = 0._R8P
do t=1, n_iter
   if (mod(t,1)==0.and.adam%myrank==0) print '(A)', 'track iteration '//trim(str(t, .true.))
   call laplace%mark_by_grad_q
   call laplace%update_ghost(q=adam%field%q)
   call adam%amr_update(is_marked_by_field=.true., do_blocks_reorder=.false.)
   if (mod(t,1)==0.and.adam%myrank==0) print '(A)', 'blocks number: '//trim(str(adam%tree%nodes_number, .true.))
   call laplace%integrate(t=time, Dt=0.000001_R8P)
   if (mod(t,n_save)==0) &
      call adam%save_hdf5(basename='laplace-'//trim(strz(t,9)), q=laplace%field%q, q_name=['T'], with_cell_morton=.true.)
   time = time + 0.2_R8P
enddo
call adam%finalize
endprogram adam_laplace
