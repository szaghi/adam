!< ADAM, test hash table class.
program adam_test_hash_table_object
!< ADAM, test hash table class.

use adam_objects
use PENF, only : FI8P, I8P, FI4P, I4P, str

implicit none

type(hash_table_object) :: hash_table                !< Hash table.
character(len=KEY_LEN)  :: key                       !< Hash table node key.
integer(I8P)            :: content                   !< Hash table node content.
character(len=KEY_LEN)  :: keys(6)                   !< Hash table node keys.
integer(I8P)            :: max_content               !< Maximum content value.
integer(I4P)            :: b, k, l, tijk(3), bijk(3) !< Counter.

print '(A)', 'initialize hash table'
call hash_table%initialize

print '(A)', 'prepare keys'
keys(1) = key_str(l=4,  tijk=[3,5,6     ], bijk=[1023,34054,35667])
keys(2) = key_str(l=3,  tijk=[25,5,9    ], bijk=[13,3454,30567   ])
keys(3) = key_str(l=1,  tijk=[43,8,6    ], bijk=[123,354,35667   ])
keys(4) = key_str(l=4,  tijk=[23,5,123  ], bijk=[1023,3454,35667 ])
keys(5) = key_str(l=7,  tijk=[73,14355,6], bijk=[30,344,35667    ])
keys(6) = key_str(l=12, tijk=[143,5,1126], bijk=[120,34054,30667 ])
do k=1, size(keys, dim=1)
   print '(A,'//FI4P//')', 'key "'//keys(k)//'" hashed in bucket: ', hash_table%hash(key=keys(k))
enddo
print '(A)', 'convert back keys'
do k=1, size(keys, dim=1)
   call key_int(key=keys(k), l=l, tijk=tijk, bijk=bijk)
   print '(A)', 'level: '//trim(str(l))//' tijk:'//trim(str(tijk))//' bijk:'//trim(str(bijk))
enddo

do k=1, size(keys, dim=1)
   call hash_table%add_node(key=keys(k), content=int(k, I8P))
enddo

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

print '(A,'//FI4P//')', 'Remove key "'//keys(3)//'" hashed in bucket: ', hash_table%hash(key=keys(3))
call hash_table%remove_node(key=keys(3))

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

print '(A,'//FI4P//')', 'Add key "'//keys(4)//'" hashed in bucket: ', hash_table%hash(key=keys(4))
call hash_table%add_node(key=keys(4), content=1_I8P)

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
