!< ADAM, test tree node class.
program adam_test_tree_node
!< ADAM, test tree node class.

use adam_tree_node_object
use PENF, only : I8P
implicit none

type(tree_node_object)          :: tree_node                 !< Tree node.
type(tree_node_object), pointer :: tree_node_pointer=>null() !< Tree node pointer.

print '(A)', 'initialize tree node'
call tree_node%initialize(code=33_I8P)

print '(A)', 'destroy tree node'
call tree_node%destroy

print '(A)', 'initialize tree node pointer'
allocate(tree_node_pointer)
call tree_node_pointer%initialize(code=3_I8P)

print '(A)', 'destroy tree node pointer'
call destroy_tree_node(node=tree_node_pointer)
endprogram adam_test_tree_node
