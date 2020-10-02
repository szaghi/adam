program adam_test_app

! TODO shift dell'ancestor e find finest_code

use adam_objects
use PENF, only : I8P, I4P, str

implicit none

type(tree_object)               :: tree        !< The tree.
type(tree_node_object), pointer :: tree_node   !< Pointer to node.
integer(I8P)                    :: code        !< Counter.

print '(A)', 'initialize tree'
call tree%initialize

print '(A)', 'add ancestor node -1'
call tree%add_node(code=-1_I8P, content=-1_I8P)
print '(A)', 'loop in tree'
do while(tree%loop(node=tree_node))
   print '(A)', 'code: '//trim(str(tree_node%code))//' content:'//trim(str(tree_node%content))
enddo
print*, ''

print '(A)', 'refine level 0'
call tree%refine(force_all=.true.)
print '(A)', 'loop in tree'
do while(tree%loop(node=tree_node))
   print '(A)', 'code: '//trim(str(tree_node%code))//' content:'//trim(str(tree_node%content))
enddo
print*, ''

print '(A)', 'refine level 1'
call tree%refine(force_all=.true.)
print '(A)', 'loop in tree'
do while(tree%loop(node=tree_node))
   print '(A)', 'code: '//trim(str(tree_node%code))//' content:'//trim(str(tree_node%content))
enddo
print*, ''
!type code_object
!   integer(I8P) :: code
!   integer(I8P) :: finest_code
!endtype code_object

!call octree%initialize
!key = key_str(l=0_I4P, tijk=[0_I4P,0_I4P,0_I4P], bijk=[0_I4P,0_I4P,0_I4P])
!call octree%add_node(key=key, content=int(0, I8P))

!do l=1, 2
!   do b=1, hash_table%buckets_number
!      do while(hash_table%bucket(b)%loop(key=key, content=content, code=code))
!         do i=0,7
!            call octree%add_node(code=child(code=code, i=i, ratio=8_I4P), content=content)
!         enddo
!         call octree%remove(key=key)
!      enddo
!   enddo
!enddo
!allocate(codes(       1:octree%nodes_number))
!allocate(finest_codes(1:octree%nodes_number))
!c = 0
!do b=1, hash_table%buckets_number
!   do while(hash_table%bucket(b)%loop(key=key, content=content, code=code))
!      c = c + 1
!      finest_codes(c) = finest_code
!             codes(c) = code
!   enddo
!enddo

!contains
!   subroutine iterator_max(node, done)
!   !< Iterator that computes the max of contents.
!   type(dictionary_node_object), pointer, intent(in)  :: node !< Actual node pointer in the dictionary.
!   logical,                               intent(out) :: done !< Flag to set to true to stop traversing.

!   max_content = max(max_content, node%content)
!   done = .false.
!   endsubroutine iterator_max
endprogram adam_test_app
