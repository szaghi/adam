!< ADAM, test tree class.
program adam_test_tree_object
!< ADAM, test tree class.

use adam_objects
use PENF, only : I8P, I4P, str
use vtk_fortran, only : vtk_file

implicit none

type(tree_object)               :: tree               !< Tree.
type(tree_node_object), pointer :: tree_node          !< Pointer to node.
integer(I8P), allocatable       :: block_to_refine(:) !< List of field blocks to be refined.
integer(I8P), allocatable       :: block_refined(:)   !< List of field refined blocks.
integer(I8P)                    :: code               !< Tree node code.
integer(I8P)                    :: offset             !< Tree node code offset.
integer(I8P)                    :: content            !< Tree node content.
integer(I8P)                    :: max_code           !< Maximum code value.
integer(I4P)                    :: l, i, j, k         !< Counter.

print '(A)', 'initialize tree'
call tree%initialize
print '(A)', 'testing buckets number calculation for 10**7 nodes: '//trim(str(tree%prime_buckets_number(nodes_number=10**7)))
print '(A)', 'hash codes'
do code=0_I8P, 83_I8P
   print '(A)', 'code "'//trim(str(code))//'" hashed in bucket: '//trim(str(tree%hash(code=code)))
enddo

print*, ''
print '(A)', 'add nodes to the tree'
call tree%refine(force_all=.true., block_to_refine=block_to_refine, block_refined=block_refined)
call tree%refine(force_all=.true., block_to_refine=block_to_refine, block_refined=block_refined)
print '(A)', 'tree nodes number:', trim(str(tree%nodes_number))
print '(A)', 'loop into tree'
do while(tree%loop(node=tree_node))
   print '(A)', ' code: '//trim(str(tree_node%code))
enddo
print '(A)', 'tree global nodes codes min/max: '//trim(str(minval(tree%code(1,:))))//'/'//trim(str(maxval(tree%code(2,:))))

print*, ''
max_code = 0
call tree%traverse(iterator=iterator_max)
print '(A)', 'maximum code value = '//trim(str(max_code))

print '(A)', 'Remove code "7" hashed in bucket: '//trim(str(tree%hash(code=7_I8P)))
call tree%remove_node(code=7_I8P)
print '(A)', 'tree nodes number:'//trim(str(tree%nodes_number))
print '(A,L1)', 'tree has code "7"?', tree%has_code(code=7_I8P)
print '(A)', 'loop into tree'
do while(tree%loop(code=code))
   print '(A)', ' code: '//trim(str(code))
enddo

print*, ''
print '(A)', 'destroy tree'
call tree%destroy
print '(A)', 'tree nodes number:'//trim(str(tree%nodes_number))
print '(A)', 're-initialize tree'
call tree%initialize
print '(A)', 'tree nodes number:'//trim(str(tree%nodes_number))
print '(A)', 'loop into tree'
do while(tree%loop(code=code))
   print '(A)', ' code: '//trim(str(code))
enddo

print*, ''
print '(A)', 'destroy tree'
call tree%destroy
print '(A)', 're-initialize tree'
call tree%initialize(nodes_number=100)
call tree%refine(force_all=.true., block_to_refine=block_to_refine, block_refined=block_refined)
call tree%refine(force_all=.true., block_to_refine=block_to_refine, block_refined=block_refined)
print '(A)', 'resize tree'
call tree%resize(nodes_number=500)
print '(A)', 'loop into tree'
do while(tree%loop(code=code))
   print '(A)', ' code: '//trim(str(code))
enddo

print*, ''
print '(A)', 'test Morton ordering'
call tree%destroy
call tree%initialize(ratio=8_I4P)
print*, ''
print '(A)', 'first/last at levels:'
do l=1, 17
   print '(A)', 'l: '//trim(str(l))//' first: '//trim(str(tree%first_at_level(level=l)))//&
                                     ' last: '//trim(str(tree%last_at_level(level=l)))
enddo
print*, ''
print '(A)', 'from Morton code to coordinates and viceversa:'
do code=0, tree%last_at_level(level=2), tree%ratio
   do offset=0, tree%ratio-1
      call tree%morton_to_coordinates(code=code+offset, i=i, j=j, k=k, l=l)
      print '(A)', 'code:'//trim(str(code+offset))//' i,j,k,l: '//trim(str([i,j,k,l]))//' c-check: '//&
                   trim(str(tree%coordinates_to_morton(i=i, j=j, k=k, l=l)))
   enddo
   print*, ''
enddo
print '(A)', 'loop in tree'
do while(tree%loop(node=tree_node))
   call tree%print_code_topology(code=tree_node%code, whole=.true.)
enddo
print*, ''
do l=1, 2
   print '(A)', 'create children of level '//trim(str(l))
   call tree%refine(force_all=.true., block_to_refine=block_to_refine, block_refined=block_refined)
   print '(A)', 'loop in tree'
   do while(tree%loop(node=tree_node))
      call tree%print_code_topology(code=tree_node%code, whole=.true.)
   enddo
   print*, ''
enddo
print*, ''
print '(A)', 'first common parent'
print '(A)', 'codes: 23, 7  first common parent: '//trim(str(tree%first_common_parent(code1=23_I8P, code2=7_I8P )))
print '(A)', 'codes: 7, 23  first common parent: '//trim(str(tree%first_common_parent(code1=7_I8P,  code2=23_I8P)))
print '(A)', 'codes: 29, 23 first common parent: '//trim(str(tree%first_common_parent(code1=29_I8P, code2=23_I8P)))
print '(A)', 'codes: 13, 14 first common parent: '//trim(str(tree%first_common_parent(code1=13_I8P, code2=14_I8P)))
print '(A)', 'codes: 48, 54 first common parent: '//trim(str(tree%first_common_parent(code1=48_I8P, code2=54_I8P)))
contains
   subroutine iterator_max(node, done)
   !< Iterator that computes the max of contents.
   type(tree_node_object), pointer, intent(in)  :: node !< Actual node pointer in the dictionary.
   logical,                         intent(out) :: done !< Flag to set to true to stop traversing.

   max_code = max(max_code, node%code)
   done = .false.
   endsubroutine iterator_max
endprogram adam_test_tree_object
