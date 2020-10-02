!< ADAM, test tree class.
program adam_test_tree_object_heavy_fill
!< ADAM, test tree class.

use adam_objects
use PENF, only : I8P, I4P, str

implicit none

type(tree_object) :: tree        !< Tree.
integer(I8P)      :: code        !< Tree node code.
integer(I8P)      :: content     !< Tree node content.
integer(I8P)      :: timing(0:2) !< Tic toc timing.

print '(A)', 'initialize tree'
call tree%initialize(nodes_number=10_I4P**7)

print '(A)', 'fill tree'
call system_clock(timing(1))
do code=0_I8P, 10_I8P**7
   call tree%add_node(code=code, content=content)
enddo
call system_clock(timing(2), timing(0))
print '(A, F8.3)', 'timing: ', real(timing(2) - timing(1))/ timing(0)

print '(A)', 'get some nodes'
call system_clock(timing(1))
do code=1_I8P, 10_I8P**7, 400_I8P ! 25000 gets
   content = tree%node_content(code=code)
enddo
call system_clock(timing(2), timing(0))
print '(A, F8.3)', 'timing: ', real(timing(2) - timing(1))/ timing(0)

print '(A)', 'add some nodes'
call system_clock(timing(1))
do code=7_I8P, 10_I8P**7, 400_I8P ! 25000 adds
   call tree%add_node(code=code, content=code-1_I8P)
enddo
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
