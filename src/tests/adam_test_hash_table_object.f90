!< ADAM, test hash table class.
program adam_test_hash_table_object
!< ADAM, test hash table class.

use adam_objects
use PENF, only : FI8P, I8P, FI4P, I4P

implicit none

type(hash_table_object) :: hash_table  !< Hash table.
character(len=KEY_LEN)  :: key         !< Dictionary node key.
integer(I8P)            :: content     !< Dictionary node content.
integer(I8P)            :: max_content !< Maximum content value.
integer(I4P)            :: b           !< Counter.

print '(A)', 'initialize hash table'
call hash_table%initialize

print '(A,'//FI4P//')', 'Add key "'//repeat('a', KEY_LEN)//'" hashed in bucket: ', hash_table%hash(key=repeat('a', KEY_LEN))
call hash_table%add_node(key=repeat('a', KEY_LEN), content=1_I8P)

print '(A,'//FI4P//')', 'Add key "'//repeat('b', KEY_LEN)//'" hashed in bucket: ', hash_table%hash(key=repeat('b', KEY_LEN))
call hash_table%add_node(key=repeat('b', KEY_LEN), content=2_I8P)

print '(A,'//FI4P//')', 'Add key "'//repeat('c', KEY_LEN)//'" hashed in bucket: ', hash_table%hash(key=repeat('c', KEY_LEN))
call hash_table%add_node(key=repeat('c', KEY_LEN), content=3_I8P)

print '(A,'//FI4P//')', 'Add key "'//repeat('d', KEY_LEN)//'" hashed in bucket: ', hash_table%hash(key=repeat('d', KEY_LEN))
call hash_table%add_node(key=repeat('d', KEY_LEN), content=4_I8P)

print '(A,'//FI4P//')', 'Add key "'//repeat('e', KEY_LEN)//'" hashed in bucket: ', hash_table%hash(key=repeat('e', KEY_LEN))
call hash_table%add_node(key=repeat('e', KEY_LEN), content=5_I8P)

print '(A,'//FI4P//')', 'Add key "'//repeat('f', KEY_LEN)//'" hashed in bucket: ', hash_table%hash(key=repeat('f', KEY_LEN))
call hash_table%add_node(key=repeat('f', KEY_LEN), content=6_I8P)

print '(A,'//FI4P//')', 'Hash table nodes number:', hash_table%nodes_number

print '(A)', 'loop into hash table'
do b=1, hash_table%buckets_number
   do while(hash_table%bucket(b)%loop(key=key, content=content))
      print '(A,'//FI8P//')', ' node: "'//key//'" =', content
   enddo
enddo

max_content = 0
call hash_table%traverse(iterator=iterator_max)
print '(A,'//FI8P//')', 'maximum content value = ', max_content

print '(A,'//FI4P//')', 'Remove key "'//repeat('d', KEY_LEN)//'" hashed in bucket: ', hash_table%hash(key=repeat('d', KEY_LEN))
call hash_table%remove_node(key=repeat('d', KEY_LEN))

print '(A,'//FI4P//')', 'Hash table nodes number:', hash_table%nodes_number

print '(A)', 'loop into hash table'
do b=1, hash_table%buckets_number
   do while(hash_table%bucket(b)%loop(key=key, content=content))
      print '(A,'//FI8P//')', ' node: "'//key//'" =', content
   enddo
enddo

print '(A)', 'destroy hash table'
call hash_table%destroy
print '(A,'//FI4P//')', 'Hash table nodes number:', hash_table%nodes_number
print '(A)', 're-initialize hash table'
call hash_table%initialize

print '(A,'//FI4P//')', 'Add key "'//repeat('z', KEY_LEN)//'" hashed in bucket: ', hash_table%hash(key=repeat('z', KEY_LEN))
call hash_table%add_node(key=repeat('z', KEY_LEN), content=1_I8P)

print '(A,'//FI4P//')', 'Hash table nodes number:', hash_table%nodes_number

print '(A)', 'loop into hash table'
do b=1, hash_table%buckets_number
   do while(hash_table%bucket(b)%loop(key=key, content=content))
      print '(A,'//FI8P//')', ' node: "'//key//'" =', content
   enddo
enddo

contains
   subroutine iterator_max(node, done)
   !< Iterator that computes the max of contents.
   type(dictionary_node_object), pointer, intent(in)  :: node !< Actual node pointer in the dictionary.
   logical,                               intent(out) :: done !< Flag to set to true to stop traversing.

   max_content = max(max_content, node%content)
   done = .false.
   endsubroutine iterator_max
endprogram adam_test_hash_table_object
