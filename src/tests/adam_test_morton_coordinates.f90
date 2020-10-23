!< ADAM, test Morton coordinates class.
program adam_test_morton_coordinates
!< ADAM, test Morton coordinates class.
!< Validate against "test-morton-coordinates.py".

use PENF, only : I8P, I4P, str
use MORTIF, only : demorton3D

implicit none

integer(I8P)      :: code       !< Tree node code.
integer(I4P)      :: l, i, j, k !< Counter.

do code=0, 20000
   ! call tree%morton_to_coordinates(code=code, i=i, j=j, k=k, l=l)
   call demorton3D(code=code, i=i, j=j, k=k)
   print '(A)', trim(str(code,.true.))//' ('//trim(str(i,.true.))//', '//&
                                              trim(str(j,.true.))//', '//&
                                              trim(str(k,.true.))//')'
enddo

endprogram adam_test_morton_coordinates
