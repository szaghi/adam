!< ADAM, test ADAM class tracking a STL surface.
program adam_test_track_surface_stl
!< ADAM, test ADAM class tracking a STL surface.

use adam_adam_object
use adam_parameters
use adam_tree_node_object
use PENF

implicit none

type(adam_object) :: adam            !< ADAM.
integer(I4P)      :: l, t, st        !< Counter.
logical           :: is_grid_changed !< Flag to check grid changes.
type(tree_node_object), pointer :: node            !< Node tree.

print '(A)', 'initialize ADAM'
call adam%initialize(max_level=6, emin=[0._R8P,0._R8P,0._R8P], emax=[1._R8P,1._R8P,1._R8P], nb=50000)

do l=1, 2
   print '(A)', 'refine ADAM at level '//trim(str(l))
   call adam%tree%mark_all_nodes(mark=TO_BE_REFINED)
   call adam%amr_update(do_blocks_reorder=.false.)
enddo

call adam%tree%load_surface_stl(file_name='src/tests/space_shuttle.stl')

do t=1, 1
   print '(A)', 'refine close to STL iteration '//trim(str(t, .true.))
   sub_iteration_loop : do st=1, 2
      call adam%tree%mark_surface_stl(surface_stl=adam%tree%surface_stl)
      call adam%amr_update(is_grid_changed=is_grid_changed)
      call adam%tree%compute_surface_stl_distance(surface_stl=adam%tree%surface_stl)
      if (.not.is_grid_changed) exit sub_iteration_loop
   enddo sub_iteration_loop
   call adam%tree%compute_surface_stl_distance(surface_stl=adam%tree%surface_stl, from_cell=.true., cell_distance=adam%field%ls)
   call adam%save_hdf5(basename='space_shuttle-'//trim(strz(t,9)), with_ls=.true.)
enddo

! do while(adam%tree%loop(node=node))
!    print*, "node: ", node%code, " distance: ", node%surface_stl_distance
! enddo

print '(A)', 'finalize ADAM'
call adam%finalize

endprogram adam_test_track_surface_stl
