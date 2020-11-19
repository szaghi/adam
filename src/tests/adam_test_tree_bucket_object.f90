!< ADAM, test tree bucket class.
program adam_test_tree_bucket_object
!< ADAM, test tree bucket class.

use adam_tree_bucket_object
use adam_tree_node_object
use PENF, only : I8P, str

implicit none

type(tree_bucket_object)        :: tree_bucket       !< Tree bucket.
type(tree_node_object), pointer :: tree_node=>null() !< Tree bucket node pointer.
integer(I8P)                    :: code              !< Tree node code.
integer(I8P)                    :: max_code          !< Maximum code value.

print '(A)', 'initialize tree bucket'

call tree_bucket%add_node(code=1_I8P)
call tree_bucket%add_node(code=2_I8P)
call tree_bucket%add_node(code=3_I8P)
call tree_bucket%add_node(code=4_I8P)
call tree_bucket%add_node(code=5_I8P)
call tree_bucket%add_node(code=6_I8P)

print '(A)', 'loop into tree bucket'
do while(tree_bucket%loop(code=code))
   print '(A)', ' node: "'//trim(str(code))
enddo

print '(A,L1)', 'tree bucket has code "3"?', tree_bucket%has_code(code=3_I8P)
print '(A)', 'getting node 3rd pointer'
tree_node => tree_bucket%node(code=3_I8P)
print '(A)', ' node: "3" ='//trim(str(tree_node%refinement_needed))

print '(A)', 'removing node 3rd'
call tree_bucket%remove_node(code=3_I8P)
print '(A,L1)', 'tree bucket has "2" code?', tree_bucket%has_code(code=2_I8P)
print '(A,L1)', 'tree bucket has "3" code?', tree_bucket%has_code(code=3_I8P)
print '(A,L1)', 'tree bucket has "4" code?', tree_bucket%has_code(code=4_I8P)

max_code = 0
call tree_bucket%traverse(iterator=iterator_max)
print '(A)', 'maximum code value = '//trim(str(max_code))

contains
   subroutine iterator_max(node, done)
   !< Iterator that computes the max of codes.
   type(tree_node_object), pointer, intent(in)  :: node !< Actual node pointer in the tree bucket.
   logical,                         intent(out) :: done !< Flag to set to true to stop traversing.

   max_code = max(max_code, node%code)
   done = .false.
   endsubroutine iterator_max
endprogram adam_test_tree_bucket_object
