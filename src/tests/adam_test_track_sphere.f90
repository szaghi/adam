!< ADAM, test track a sphere.
program adam_test_track_sphere
!< ADAM, test track a sphere.

use adam_objects
use PENF, only : R8P, I8P, I4P, str, strz

implicit none

type(tree_object)               :: tree                 !< The tree.
type(tree_node_object), pointer :: node                 !< Pointer to node.
type(field_object)              :: field                !< Field.
integer(I8P), allocatable       :: block_to_refine(:)   !< List of field blocks to be refined.
integer(I8P), allocatable       :: block_refined(:,:)   !< List of field refined blocks with Morton code.
integer(I8P), allocatable       :: block_to_derefine(:) !< List of field blocks to be derefined.
integer(I8P), allocatable       :: block_derefined(:,:) !< List of field derefined blocks with Morton code.
integer(I4P)                    :: t, st                !< Counter.

print '(A)', 'sphere tracking, initialize grid'
call tree%initialize
call field%initialize(nb=190000, emin=[0._R8P,0._R8P,0._R8P], emax=[1._R8P,1._R8P,1._R8P])
do t=1,2
   call tree%mark_all_nodes(mark=NODE_TO_BE_REFINED)
   call tree%adapt(block_to_refine=block_to_refine, block_refined=block_refined, &
                   block_to_derefine=block_to_derefine, block_derefined=block_derefined)
   call field%adapt(ratio=tree%ratio, block_to_refine=block_to_refine, block_refined=block_refined, &
                    block_to_derefine=block_to_derefine, block_derefined=block_derefined)
enddo

print '(A)', 'catch first position'
do t=1,8
   print*, ''
   print '(A)', 'track iteration '//trim(str(t, .true.))
   call mark_sphere_nodes(tree=tree, field=field, center=[0.2_R8P,0.5_R8P,0.5_R8P], radius=0.1_R8P)
   call tree%adapt(block_to_refine=block_to_refine, block_refined=block_refined, &
                   block_to_derefine=block_to_derefine, block_derefined=block_derefined)
   print '(A)', 'nodes refined n. '//trim(str(size(block_refined(1,:), dim=1),.true.))
   print '(A)', 'nodes derefined n. '//trim(str(size(block_derefined(1,:), dim=1),.true.))
   call field%adapt(ratio=tree%ratio, block_to_refine=block_to_refine, block_refined=block_refined, &
                    block_to_derefine=block_to_derefine, block_derefined=block_derefined)
enddo

print*, ''
print '(A)', 'move sphere'
do t=1,10
   print*, ''
   print '(A)', 'track iteration '//trim(str(t, .true.))//' position x='//trim(str(0.2_R8P + t*0.05_R8P))
   sub_iteration_loop : do st=1, 10
      print '(A)', '  track su-iteration '//trim(str(st, .true.))
      call mark_sphere_nodes(tree=tree, field=field, center=[0.2_R8P + t*0.05_R8P,0.5_R8P,0.5_R8P], radius=0.1_R8P)
      call tree%adapt(block_to_refine=block_to_refine, block_refined=block_refined, &
                      block_to_derefine=block_to_derefine, block_derefined=block_derefined)
      print '(A)', '  nodes refined n. '//trim(str(size(block_refined(1,:), dim=1),.true.))
      print '(A)', '  nodes derefined n. '//trim(str(size(block_derefined(1,:), dim=1),.true.))
      if (size(block_refined(1,:), dim=1)==0_I4P.and.size(block_derefined(1,:), dim=1)==0_I4P) exit sub_iteration_loop
      call field%adapt(ratio=tree%ratio, block_to_refine=block_to_refine, block_refined=block_refined, &
                       block_to_derefine=block_to_derefine, block_derefined=block_derefined)
   enddo sub_iteration_loop
   call field_save_vtk(tree=tree, field=field, basename='sphere-t-'//trim(strz(t,3)), directory='sphere/')
enddo
endprogram adam_test_track_sphere
