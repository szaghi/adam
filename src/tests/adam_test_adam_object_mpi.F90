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
integer(I4P)                    :: l, t, st        !< Counter.
logical                         :: is_grid_changed !< Flag to check grid changes.

print '(A)', 'initialize ADAM'
call adam%initialize(max_level=8, emin=[0._R8P,0._R8P,0._R8P], emax=[1._R8P,1._R8P,1._R8P], nb=500000)

do l=1, 1
   print '(A)', 'refine ADAM at level '//trim(str(l))
   call adam%tree%mark_all_nodes(mark=TO_BE_REFINED)
   call adam%amr_update
enddo
print '(A)', 'set initial conditions'

!!!call adam%field%set_initial_conditions
!!!call adam%field%rk_integrate(t=t*0.1_R8P, Dt=0.1_R8P)
!!!! adam%field%u_s(:,:,:,1:adam%field%blocks_number,1) = adam%field%u(:,:,:,1:adam%field%blocks_number)
!!!! call adam%field%update_ghost(s=1)
!!!! adam%field%u(:,:,:,1:adam%field%blocks_number) = adam%field%u_s(:,:,:,1:adam%field%blocks_number,1)
!!!! call adam%save_vtk(basename='sphere-'//trim(strz(0,9)))
!!!call adam%save_hdf5(basename='sphere-'//trim(strz(0,9)))
!!!call adam%finalize

node => adam%tree%node(code=1_I8P)
node%refinement_needed = TO_BE_REFINED
call adam%amr_update
call adam%field%set_initial_conditions
call adam%field%rk_integrate(t=t*0.1_R8P, Dt=0._R8P)
adam%field%u_s(:,:,:,1:adam%field%blocks_number,1) = adam%field%u(:,:,:,1:adam%field%blocks_number)
call adam%field%update_ghost(s=1)
adam%field%u(:,:,:,1:adam%field%blocks_number) = adam%field%u_s(:,:,:,1:adam%field%blocks_number,1)
call adam%save_hdf5(basename='sphere-'//trim(strz(1,9)))
call adam%save_vtk(basename='sphere-'//trim(strz(1,9)))

call adam%finalize

do t=1, 25
   print '(A)', 'track iteration '//trim(str(t, .true.))
   call adam%field%rk_integrate(t=t*0.1_R8P, Dt=0.1_R8P)
   sub_iteration_loop : do st=1, 1
      ! call adam%field%mark_sphere(center=[0.2_R8P+t*0.05_R8P,0.5_R8P,0.5_R8P], radius=0.1_R8P)
      call adam%field%mark_sphere(center=[0.2_R8P,0.5_R8P,0.5_R8P], radius=0.1_R8P)
      call adam%amr_update(is_marked_by_field=.true., is_grid_changed=is_grid_changed)
      if (.not.is_grid_changed) exit sub_iteration_loop
   enddo sub_iteration_loop
   call adam%save_hdf5(basename='sphere-'//trim(strz(t,9)))
enddo

print '(A)', 'finalize ADAM'
call adam%finalize

endprogram adam_test_adam_object_mpi
