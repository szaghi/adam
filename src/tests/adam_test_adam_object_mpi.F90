!< ADAM, test ADAM class in MPI env.
program adam_test_adam_object_mpi
!< ADAM, test ADAM class in MPI env.

use adam_adam_object
use adam_parameters
use adam_tree_node_object
use PENF

implicit none

type(adam_object)               :: adam            !< ADAM.
type(tree_node_object), pointer :: node            !< Node tree.
integer(I4P)                    :: l, t, st, b     !< Counter.
logical                         :: is_grid_changed !< Flag to check grid changes.
integer(I8P)                    :: timing(0:2)     !< Tic toc timing.
real(R8P)                       :: residual        !< Global residual.

print '(A)', 'initialize ADAM'
call adam%initialize(max_level=8, emin=[0._R8P,0._R8P,0._R8P], emax=[1._R8P,1._R8P,1._R8P], nb=50000)

do l=1, 2
   print '(A)', 'refine ADAM at level '//trim(str(l))
   call adam%tree%mark_all_nodes(mark=TO_BE_REFINED)
   call adam%amr_update
enddo
print '(A)', 'set initial conditions'

! print*, 'Reorder as tree saw'
! do while(adam%tree%loop(node=node))
!    print*, 'code,block,proc: ', node%code, node%block_index, node%myrank
! enddo
! print*, 'Reorder as field saw'
! do b=1, adam%field%blocks_number
!    print*, 'code,block,proc: ', adam%field%code(b), b, adam%field%myrank
! enddo
! call adam%finalize

call adam%field%set_initial_conditions
call adam%save_hdf5(basename='sphere-'//trim(strz(0,9)))

! do t=1, 3
!    print '(A)', 'track iteration '//trim(str(t, .true.))
!    ! call adam%field%rk_integrate(t=t*0.1_R8P, Dt=0.1_R8P)
!    sub_iteration_loop : do st=1, 1
!       ! call adam%field%mark_sphere(center=[0.2_R8P+t*0.05_R8P,0.5_R8P,0.5_R8P], radius=0.1_R8P)
!       call adam%field%mark_sphere(center=[0.2_R8P,0.5_R8P,0.5_R8P], radius=0.1_R8P)
!       call adam%amr_update(do_blocks_reorder=.false., is_marked_by_field=.true., is_grid_changed=is_grid_changed)
!       if (.not.is_grid_changed) exit sub_iteration_loop
!    enddo sub_iteration_loop
!    ! call adam%save_hdf5(basename='sphere-'//trim(strz(t,9)))
! enddo

call system_clock(timing(1))
do t=1, 4
   call adam%save_hdf5(basename='sphere-'//trim(strz(t,9)))
   call adam%field%rk_integrate(t=t*0.1_R8P, Dt=0.1_R8P, residual=residual)
   print '(A)', 'integrate iteration: '//trim(str(t, .true.))//' global residual: '//trim(str(residual))
enddo
call system_clock(timing(2), timing(0))
print '(A, F8.3)', 'timing: ', real(timing(2) - timing(1))/ timing(0)/10
call adam%save_hdf5(basename='sphere-'//trim(strz(t,9)))

print '(A)', 'finalize ADAM'
call adam%finalize

endprogram adam_test_adam_object_mpi
