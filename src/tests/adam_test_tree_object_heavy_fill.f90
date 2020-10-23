!< ADAM, test tree class.
program adam_test_tree_object_heavy_fill
!< ADAM, test tree class.

use adam_objects
use PENF, only : I8P, I4P, str

implicit none

type(tree_object)               :: tree                 !< Tree.
type(tree_node_object), pointer :: node                 !< Pointer to node.
integer(I8P)                    :: code                 !< Tree node code.
integer(I4P)                    :: level                !< Refinement level counter.
integer(I8P), allocatable       :: block_to_refine(:)   !< List of field blocks to be refined.
integer(I8P), allocatable       :: block_refined(:,:)   !< List of field refined blocks with Morton code.
integer(I8P), allocatable       :: block_to_derefine(:) !< List of field blocks to be derefined.
integer(I8P), allocatable       :: block_derefined(:,:) !< List of field derefined blocks with Morton code.
integer(I8P)                    :: timing(0:2)          !< Tic toc timing.

print '(A)', 'initialize tree'
call tree%initialize(nodes_number=10_I8P**7)

print '(A)', 'fill tree with some levels'
do level=1, 7
   call system_clock(timing(1))
   call tree%mark_all_nodes(mark=NODE_TO_BE_REFINED)
   call tree%adapt(block_to_refine=block_to_refine, block_refined=block_refined, &
                   block_to_derefine=block_to_derefine, block_derefined=block_derefined)
   call system_clock(timing(2), timing(0))
   print '(A, F8.3)', 'level '//trim(str(level,.true.))//' timing: ', real(timing(2) - timing(1))/ timing(0)
enddo

print '(A)', 'get some nodes'
call system_clock(timing(1))
do code=1_I8P, 10_I8P**7, 400_I8P ! 25000 gets
   node => tree%node(code=code)
enddo
call system_clock(timing(2), timing(0))
print '(A, F8.3)', 'timing: ', real(timing(2) - timing(1))/ timing(0)

print '(A)', 'add new level'
call system_clock(timing(1))
call tree%mark_all_nodes(mark=NODE_TO_BE_REFINED)
call tree%adapt(block_to_refine=block_to_refine, block_refined=block_refined, &
                block_to_derefine=block_to_derefine, block_derefined=block_derefined)
call system_clock(timing(2), timing(0))
print '(A, F8.3)', 'timing: ', real(timing(2) - timing(1))/ timing(0)

print '(A)', 'remove some nodes'
call system_clock(timing(1))
do code=31_I8P, 10_I8P**7, 400_I8P ! 25000 adds
   call tree%remove_node(code=code)
enddo
call system_clock(timing(2), timing(0))
print '(A, F8.3)', 'timing: ', real(timing(2) - timing(1))/ timing(0)
endprogram adam_test_tree_object_heavy_fill
