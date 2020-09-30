!< ADAM, test dictionary node class.
program adam_test_dictionary_node
!< ADAM, test dictionary node class.

use adam_tree_topology
use MORTIF, only : demorton2D
use PENF, only : I4P, I8P, str
implicit none

integer(I4P) :: coord(3) !< Counter.
integer(I4P) :: l        !< Counter.
integer(I8P) :: c        !< Counter.

print '(A)', 'test Morton code functions'
print '(A)', 'first/last at levels:'
do l=1, 17
   print*, 'l:', l, ' first=', first_at_level(level=l, ratio=4_I4P), ' last=', last_at_level(level=l, ratio=4_I4P)
enddo
print '(A)', 'from Morton code to coordinates and viceversa:'
do c=0, 83, 4
   call morton_to_coordinates(code=c  , ratio=4_I4P, i=coord(1), j=coord(2), l=coord(3))
   print '(A,I3,A,3(I3,1X),A,I3)', 'c=', c  , ' i,j,l=', coord, ' c-check=', &
         coordinates_to_morton(i=coord(1), j=coord(2), l=coord(3), ratio=4_I4P)
   call morton_to_coordinates(code=c+1, ratio=4_I4P, i=coord(1), j=coord(2), l=coord(3))
   print '(A,I3,A,3(I3,1X),A,I3)', 'c=', c+1, ' i,j,l=', coord, ' c-check=', &
         coordinates_to_morton(i=coord(1), j=coord(2), l=coord(3), ratio=4_I4P)
   call morton_to_coordinates(code=c+2, ratio=4_I4P, i=coord(1), j=coord(2), l=coord(3))
   print '(A,I3,A,3(I3,1X),A,I3)', 'c=', c+2, ' i,j,l=', coord, ' c-check=', &
         coordinates_to_morton(i=coord(1), j=coord(2), l=coord(3), ratio=4_I4P)
   call morton_to_coordinates(code=c+3, ratio=4_I4P, i=coord(1), j=coord(2), l=coord(3))
   print '(A,I3,A,3(I3,1X),A,I3)', 'c=', c+3, ' i,j,l=', coord, ' c-check=', &
         coordinates_to_morton(i=coord(1), j=coord(2), l=coord(3), ratio=4_I4P)
   print*, ''
enddo
print '(A)', 'parent computation given Morton code'
do c=0, 83, 4
   print '(A,I3,A,I3)', 'c=', c  , ' parent=', parent(code=c  , ratio=4_I4P)
   print '(A,I3,A,I3)', 'c=', c+1, ' parent=', parent(code=c+1, ratio=4_I4P)
   print '(A,I3,A,I3)', 'c=', c+2, ' parent=', parent(code=c+2, ratio=4_I4P)
   print '(A,I3,A,I3)', 'c=', c+3, ' parent=', parent(code=c+3, ratio=4_I4P)
   print*, ''
enddo
print '(A)', 'child computation given Morton code'
do c=0, 83, 4
   print '(A,I3,A,I3)', 'c=', c  , ' child=', child(code=c  , i=0, ratio=4_I4P)
   print '(A,I3,A,I3)', 'c=', c+1, ' child=', child(code=c+1, i=0, ratio=4_I4P)
   print '(A,I3,A,I3)', 'c=', c+2, ' child=', child(code=c+2, i=0, ratio=4_I4P)
   print '(A,I3,A,I3)', 'c=', c+3, ' child=', child(code=c+3, i=0, ratio=4_I4P)
   print*, ''
enddo
print '(A)', 'local child and siblings computation given Morton code'
do c=0, 83, 4
   print '(A,I3,A,I1,A,3I4)', 'c=',  c  , ' child_local=', child_local(code=c  , ratio=4_I4P), ' sib=', &
      siblings(code=c  , ratio=4_I4P)
   print '(A,I3,A,I1,A,3I4)', 'c=',  c+1, ' child_local=', child_local(code=c+1, ratio=4_I4P), ' sib=', &
      siblings(code=c+1, ratio=4_I4P)
   print '(A,I3,A,I1,A,3I4)', 'c=',  c+2, ' child_local=', child_local(code=c+2, ratio=4_I4P), ' sib=', &
      siblings(code=c+2, ratio=4_I4P)
   print '(A,I3,A,I1,A,3I4)', 'c=',  c+3, ' child_local=', child_local(code=c+3, ratio=4_I4P), ' sib=', &
      siblings(code=c+3, ratio=4_I4P)
   print*, ''
enddo
print '(A)', 'path to root computation given Morton code'
do c=0, 83, 4
   print '(A,I3,A,3I4)', 'c=', c  , ' path=', path(code=c  , ratio=4_I4P)
   print '(A,I3,A,3I4)', 'c=', c+1, ' path=', path(code=c+1, ratio=4_I4P)
   print '(A,I3,A,3I4)', 'c=', c+2, ' path=', path(code=c+2, ratio=4_I4P)
   print '(A,I3,A,3I4)', 'c=', c+3, ' path=', path(code=c+3, ratio=4_I4P)
   print*, ''
enddo
print '(A)', 'level computation given Morton code'
do c=0, 83, 4
   print '(A,I3,A,I3)', 'c=', c  , ' level=', level(code=c  , ratio=4_I4P)
   print '(A,I3,A,I3)', 'c=', c+1, ' level=', level(code=c+1, ratio=4_I4P)
   print '(A,I3,A,I3)', 'c=', c+2, ' level=', level(code=c+2, ratio=4_I4P)
   print '(A,I3,A,I3)', 'c=', c+3, ' level=', level(code=c+3, ratio=4_I4P)
   print*, ''
enddo
endprogram adam_test_dictionary_node
