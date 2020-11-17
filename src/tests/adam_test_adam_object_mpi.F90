!< ADAM, test ADAM class in MPI env.
program adam_test_adam_object_mpi
!< ADAM, test ADAM class in MPI env.

use adam_adam_object
use adam_parameters
use adam_tree_node_object
use PENF

implicit none

type(adam_object) :: adam            !< ADAM.
type(tree_node_object), pointer :: node            !< ADAM.
integer(I4P)      :: l, t, st        !< Counter.
logical           :: is_grid_changed !< Flag to check grid changes.
real(R8P)         :: time            !< Time.

print '(A)', 'initialize ADAM'
call adam%initialize(max_level=6, emin=[0._R8P,0._R8P,0._R8P], emax=[1._R8P,1._R8P,1._R8P], nb=500000)

do l=1, 2
   print '(A)', 'refine ADAM at level '//trim(str(l))
   call adam%tree%mark_all_nodes(mark=TO_BE_REFINED)
   call adam%amr_update(do_mpi_redistribute=.false.) ! no need to pass is_marked_by_tree=.true.
                                                     ! because the marks are the same in all processes
enddo
print '(A)', 'redistribute ADAM nodes/blocks between processes, load balancing'
call adam%mpi_redistribute(print_mpi_stats=.false.)
print '(A)', 'set initial conditions'
call adam%field%set_initial_conditions
call adam%save_hdf5(basename='sphere-initial')

time = 0._R8P
do t=1, 25
   ! if (mod(t,1)==0) call adam%save_hdf5(basename='sphere-'//trim(strz(t,9)))
   if (mod(t,1)==0) print '(A)', 'track iteration '//trim(str(t, .true.))
   ! call adam%field%rk_integrate(t=time, Dt=0.1_R8P)
   ! sub_iteration_loop : do st=1, 10
      call adam%save_hdf5(basename='sphere-'//trim(strz(t,9)))
      ! call adam%field%mark_by_u_value(u_value=0.5_R8P)
      call adam%field%mark_sphere(center=[0.2_R8P+t*0.05_R8P,0.5_R8P,0.5_R8P], radius=0.1_R8P)
      call adam%amr_update(is_marked_by_field=.true.)
      ! if (.not.is_grid_changed) exit sub_iteration_loop
   ! enddo sub_iteration_loop
   time = time + 0.1_R8P
enddo
call adam%save_hdf5(basename='sphere-'//trim(strz(t,9)))

print '(A)', 'finalize ADAM'
call adam%finalize

endprogram adam_test_adam_object_mpi
