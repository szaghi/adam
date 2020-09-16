!< ADAM, test dictionary class.
program adam_text_dictionary
!< ADAM, test dictionary class.

use adam_objects
use PENF, only : FI4P, I4P

implicit none

type(dictionary_object)               :: dictionary              !< Dictionary.
type(dictionary_node_object), pointer :: dictionary_node=>null() !< Dictionary node pointer.
character(len=KEY_LEN)                :: key                     !< Dictionary node key.
integer(I4P)                          :: content                 !< Dictionary node content.
integer(I4P)                          :: max_content             !< Maximum content value.

print '(A)', 'initialize dictionary'

call dictionary%add_node(key=repeat('a', KEY_LEN), content=1)
call dictionary%add_node(key=repeat('b', KEY_LEN), content=2)
call dictionary%add_node(key=repeat('c', KEY_LEN), content=3)
call dictionary%add_node(key=repeat('d', KEY_LEN), content=4)
call dictionary%add_node(key=repeat('e', KEY_LEN), content=5)
call dictionary%add_node(key=repeat('f', KEY_LEN), content=6)

print '(A)', 'loop into dictionary'
do while(dictionary%loop(key=key, content=content))
   print '(A,'//FI4P//')', ' node: "'//key//'" =', content
enddo

print '(A,L1)', 'dictionary has "b" key = ', dictionary%has_key(key=repeat('b', KEY_LEN))
print '(A)', 'getting b-node pointer'
dictionary_node => dictionary%node(key=repeat('b', KEY_LEN))
print '(A,'//FI4P//')', ' node: "'//repeat('b', KEY_LEN)//'" =', dictionary_node%content

print '(A)', 'removing b-node'
call dictionary%remove_node(key=repeat('b', KEY_LEN))
print '(A,L1)', 'dictionary has "b" key = ', dictionary%has_key(key=repeat('b', KEY_LEN))
print '(A,L1)', 'dictionary has "c" key = ', dictionary%has_key(key=repeat('c', KEY_LEN))
print '(A,L1)', 'dictionary has "h" key = ', dictionary%has_key(key=repeat('h', KEY_LEN))

max_content = 0
call dictionary%traverse(iterator=iterator_max)
print '(A,'//FI4P//')', 'maximum content value = ', max_content

contains
   subroutine iterator_max(key, content, done)
   !< Iterator that computes the max of contents.
   character(len=*),  intent(in)  :: key     !< The node key.
   integer(I4P),      intent(in)  :: content !< The generic content.
   logical,           intent(out) :: done    !< Flag to set to true to stop traversing.

   max_content = max(max_content, content)
   done = .false.
   endsubroutine iterator_max
endprogram adam_text_dictionary
