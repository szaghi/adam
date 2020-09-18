!< ADAM, test hash table class.
program adam_test_hash_table_object
!< ADAM, test hash table class.

use adam_objects
use PENF, only : FI8P, I8P, FI4P, I4P, str, R8P

implicit none

type(hash_table_object) :: hash_table  !< Hash table.
character(len=KEY_LEN)  :: key         !< Hash table node key.
integer(I8P)            :: content     !< Hash table node content.
integer(I4P)            :: l           !< Counter.
integer(I4P)            :: bi, bj, bk  !< Counter.
integer(I8P)            :: timing(0:2) !< Tic toc timing.

print '(A)', 'initialize hash table'
call hash_table%initialize(buckets_number=10352717)

print '(A)', 'fill hash table'
call system_clock(timing(1))
do l=1, 10
   print*, l
   do bk=1, 100
      do bj=1, 100
         do bi=1, 100
            key = key_str(l=l, tijk=[0,0,0], bijk=[bi,bj,bk])
            call hash_table%add_node(key=key, content=int(bi+bj+bk+l, I8P))
         enddo
      enddo
   enddo
enddo
call system_clock(timing(2), timing(0))
print '(A, F8.3)', 'timing: ', real(timing(2) - timing(1))/ timing(0)

print '(A)', 'get some keys'
call system_clock(timing(1))
do bk=1, 10
   do bj=1, 100, 2
      do bi=1, 20, 4
         key = key_str(l=5, tijk=[0,0,0], bijk=[bi,bj,bk])
         content = hash_table%node_content(key=key)
      enddo
   enddo
enddo
call system_clock(timing(2), timing(0))
print '(A, F8.3)', 'timing: ', real(timing(2) - timing(1))/ timing(0)

print '(A)', 'add some keys'
call system_clock(timing(1))
do bk=1, 10
   do bj=2, 100, 2
      do bi=3, 20, 4
         key = key_str(l=5, tijk=[0,0,0], bijk=[bi,bj,bk])
         call hash_table%add_node(key=key, content=int(bi+bj+bk+2, I8P))
      enddo
   enddo
enddo
call system_clock(timing(2), timing(0))
print '(A, F8.3)', 'timing: ', real(timing(2) - timing(1))/ timing(0)

print '(A)', 'remove some keys'
call system_clock(timing(1))
do bk=1, 10
   do bj=2, 100, 2
      do bi=3, 20, 4
         key = key_str(l=5, tijk=[0,0,0], bijk=[bi,bj,bk])
         call hash_table%remove_node(key=key)
      enddo
   enddo
enddo
call system_clock(timing(2), timing(0))
print '(A, F8.3)', 'timing: ', real(timing(2) - timing(1))/ timing(0)
endprogram adam_test_hash_table_object
