!< ADAM, test ADAM for Euler equation.
program adam_test_euler
!< ADAM, test ADAM for Euler equation.

use adam_adam_object
use adam_equation_euler_cpu_object
use adam_parameters
use PENF

implicit none

type(adam_object)               :: adam            !< ADAM.
type(equation_euler_cpu_object) :: euler           !< Euler equations system.
integer(I4P)                    :: l, t            !< Counter.
logical                         :: is_grid_changed !< Flag to check grid changes.
real(R8P)                       :: time            !< Time.
integer(I4P)                    :: n_iter=300      !< Number of iterations.
integer(I4P)                    :: n_save=10       !< Frequency of saving output.

print '(A)', 'Laplace equation integration'
call adam%initialize(max_level=7,                                              &
                     emin=[0._R8P,0._R8P,0._R8P], emax=[1._R8P,1._R8P,1._R8P], &
                     ni=8_I4P, nj=8_I4P, nk=8_I4P, gc=[2_I4P,2_I4P,2_I4P],     &
                     bc_type=[BC_EXTRAPOLATION,BC_EXTRAPOLATION,               &
                              BC_EXTRAPOLATION,BC_EXTRAPOLATION,               &
                              BC_EXTRAPOLATION,BC_EXTRAPOLATION],              &
                     nb=15000, nv=5, nodes_number=16*15000_I8P)

call euler%initialize(field=adam%field, ns=1)
print '(A)', 'create 2 levels of refinement'
do l=1, 2
   print '(A)', 'refine ADAM at level '//trim(str(l))
   print *, 'blocks_number: ',adam%field%blocks_number
   call adam%tree%mark_all_nodes(mark=TO_BE_REFINED)
   call euler%update_ghost(q=adam%field%q)
   call adam%amr_update(do_blocks_reorder=.false., do_mpi_redistribute=.true.)
enddo
print '(A)', 'set initial conditions'
call euler%set_initial_conditions
do t=1, 10
   if (mod(t,1)==0.and.adam%myrank==0) print '(A)', 'track iteration '//trim(str(t, .true.))
   call euler%mark_by_grad_rho(grad_tol=10.0_R8P)
   call euler%update_ghost(q=adam%field%q)
   call adam%amr_update(is_marked_by_field=.true., do_blocks_reorder=.false., is_grid_changed=is_grid_changed)
   if (.not.is_grid_changed) exit
enddo
call adam%save_hdf5(basename='euler-'//trim(strz(0,9)), with_ghost=.false., with_cell_morton=.true.)
time = 0._R8P
do t=1, n_iter
   if (mod(t,1000)==0) then
      call euler%mark_by_grad_rho(grad_tol=10.0_R8P)
      call euler%update_ghost(q=adam%field%q)
      call adam%amr_update(is_marked_by_field=.true., do_blocks_reorder=.false.)
   endif
   call euler%compute_dt
   if (mod(t,1)==0.and.adam%myrank==0) then
      print '(A)',       'blocks number: '//trim(str(adam%tree%nodes_number, .true.))
      print '(A,F12.7)', 'time step:     ', euler%dt
      print '(A,F12.7)', 'time:          ', time
      print '(A)',       't:             '//trim(str(t,.true.))
   endif
   call euler%integrate(t=time)
   if (mod(t,n_save)==0) call adam%save_hdf5(basename='euler-'//trim(strz(t,9)), with_ghost=.false., with_cell_morton=.true.)
   time = time + euler%dt
enddo

call adam%finalize
endprogram adam_test_euler
