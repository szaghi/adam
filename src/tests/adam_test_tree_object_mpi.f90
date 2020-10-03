!< ADAM, test tree class.
program adam_test_tree_object_mpi
!< ADAM, test tree class.

use adam_objects
use PENF, only : I8P, I4P, str

implicit none

type(tree_object)               :: tree        !< Tree.
type(tree_node_object), pointer :: tree_node   !< Pointer to node.
integer(I8P)                    :: code        !< Tree node code.
! integer(I8P)                    :: offset      !< Tree node code offset.
! integer(I8P)                    :: content     !< Tree node content.
! integer(I8P)                    :: max_content !< Maximum content value.
integer(I4P)                    :: l!, i, j, k  !< Counter.

print '(A)', 'initialize tree'
call tree%initialize(ratio=8_I4P, max_level=2_I4P)
print*, ''
print '(A)', 'test codes comparison'
print '(A,L1)', '2  < 3  (T): ', tree%lower(2_I8P, 3_I8P)
print '(A,L1)', '2  > 3  (F): ', tree%greater(2_I8P, 3_I8P)
print '(A,L1)', '3  < 2  (F): ', tree%lower(3_I8P, 2_I8P)
print '(A,L1)', '3  > 2  (T): ', tree%greater(3_I8P, 2_I8P)
print '(A,L1)', '2  < 22 (F): ', tree%lower(2_I8P, 22_I8P)
print '(A,L1)', '2  > 22 (T): ', tree%greater(2_I8P, 22_I8P)
print '(A,L1)', '22 < 2  (T): ', tree%lower(22_I8P, 2_I8P)
print '(A,L1)', '22 > 2  (F): ', tree%greater(22_I8P, 2_I8P)
print*, ''

print '(A)', 'test uniform refinement'
print '(A)', 'add ancestor node to the tree, Morton code -1'
call tree%add_node(code=-1_I8P, content=-1_I8P)
do l=1, 3
   ! third level should not be done because max refinement level has been set to 2
   ! the level 2 should be printed twice
   print '(A)', 'create children of level '//trim(str(l))
   call tree%refine(force_all=.true.)
   print '(A)', 'loop in tree'
   do while(tree%loop(code=code))
      call tree%print_code_topology(code=code)
   enddo
   print*, ''
enddo

print '(A)', 'test non uniform refinement'
call tree%initialize(ratio=8_I4P, max_level=2_I4P)
print '(A)', 'add ancestor node to the tree, Morton code -1'
call tree%add_node(code=-1_I8P, content=-1_I8P)
print '(A)', 'create children of level 1'
call tree%refine(force_all=.true.)
print '(A)', 'refine nodes 2, 3, 7'
tree_node => tree%node(code=2_I8P) ; tree_node%refinement_needed = TO_BE_REFINED
tree_node => tree%node(code=3_I8P) ; tree_node%refinement_needed = TO_BE_REFINED
tree_node => tree%node(code=7_I8P) ; tree_node%refinement_needed = TO_BE_REFINED
call tree%refine()
print '(A)', 'loop in tree'
do while(tree%loop(code=code))
   call tree%print_code_topology(code=code)
enddo
print*, ''

endprogram adam_test_tree_object_mpi
