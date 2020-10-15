!< ADAM, test field class.
program adam_test_field_object
!< ADAM, test field class.

use adam_objects
use PENF, only : R8P, I8P, I4P, str

implicit none

type(tree_object)               :: tree               !< The tree.
type(tree_node_object), pointer :: node               !< Pointer to node.
type(field_object)              :: field              !< Field.
integer(I8P)                    :: code               !< Counter.
integer(I8P), allocatable       :: block_to_refine(:) !< List of field blocks to be refined.
integer(I8P), allocatable       :: block_refined(:)   !< List of field refined blocks.
real(R8P)                       :: emin(3), emax(3)

print '(A)', 'initialize'
call tree%initialize

emin = 0._R8P
emax = 1._R8P
call field%initialize(nb=100, emin=emin, emax=emax)

print '(A)', 'loop in tree'
do while(tree%loop(node=node))
   call tree%print_code_topology(code=node%code, block_index=.true., coordinates=.true.)
enddo
print*, ''

print '(A)', 'refine level 0'
call tree%refine(force_all=.true., block_to_refine=block_to_refine, block_refined=block_refined)
call field%refine(ratio=tree%ratio, block_to_refine=block_to_refine, block_refined=block_refined)
print '(A)', 'loop in tree'
do while(tree%loop(node=node))
   call tree%print_code_topology(code=node%code, block_index=.true., coordinates=.true.)
enddo
print*, ''

print '(A)', 'blocks index to refine'
print '(A)', trim(str(block_to_refine))
print '(A)', 'blocks index refined'
print '(A)', trim(str(block_refined))
print*, ''
print '(A)', 'refine level 1'
node => tree%node(code=3_I8P)
node%refinement_needed = TO_BE_REFINED
call tree%sanitize
call tree%refine(block_to_refine=block_to_refine, block_refined=block_refined)
call field%refine(ratio=tree%ratio, block_to_refine=block_to_refine, block_refined=block_refined)
print '(A)', 'loop in tree'
do while(tree%loop(node=node))
   call tree%print_code_topology(code=node%code, block_index=.true., coordinates=.true.)
enddo
print*, ''
print '(A)', 'blocks index to refine'
print '(A)', trim(str(block_to_refine))
print '(A)', 'blocks refined'
print '(A)', trim(str(block_refined))
print*, ''
print '(A)', 'refine index level 2'
node => tree%node(code=37_I8P)
node%refinement_needed = TO_BE_REFINED
call tree%sanitize
call tree%refine(block_to_refine=block_to_refine, block_refined=block_refined)
call field%refine(ratio=tree%ratio, block_to_refine=block_to_refine, block_refined=block_refined)
print '(A)', 'loop in tree'
do while(tree%loop(node=node))
   call tree%print_code_topology(code=node%code, block_index=.true., coordinates=.true.)
enddo
print*, ''
print '(A)', 'blocks index to refine'
print '(A)', trim(str(block_to_refine))
print '(A)', 'blocks index refined'
print '(A)', trim(str(block_refined))
call field%save_vtk(basename='adam')
endprogram adam_test_field_object
