!< ADAM, test dictionary node class.
program adam_test_dictionary_node
!< ADAM, test dictionary node class.

use adam_dictionary_node_object
use PENF, only : I8P
implicit none

type(dictionary_node_object) :: dictionary_node !< Dictionary node.

print '(A)', 'initialize dictionary node'
call dictionary_node%initialize(key=key_str(l=4,  tijk=[3,5,6], bijk=[1023,34054,35667]), content=-2_I8P)

print '(A)', 'destroy dictionary node'
call dictionary_node%destroy
endprogram adam_test_dictionary_node
