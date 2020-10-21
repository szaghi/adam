!< ADAM, test field class.
program adam_test_field_object
!< ADAM, test field class.

use adam_objects
use PENF, only : R8P, I8P, I4P, str

implicit none

type(tree_object)               :: tree                 !< The tree.
type(tree_node_object), pointer :: node                 !< Pointer to node.
type(field_object)              :: field                !< Field.
integer(I8P), allocatable       :: block_to_refine(:)   !< List of field blocks to be refined.
integer(I8P), allocatable       :: block_refined(:,:)   !< List of field refined blocks with Morton code.
integer(I8P), allocatable       :: block_to_derefine(:) !< List of field blocks to be derefined.
integer(I8P), allocatable       :: block_derefined(:,:) !< List of field derefined blocks with Morton code.
real(R8P)                       :: emin(3), emax(3)     !< Domain extents.
integer(I4P)                    :: c                    !< Counter.

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
call tree%mark_all_nodes_to_be_refined
call tree%adapt(block_to_refine=block_to_refine, block_refined=block_refined, &
                block_to_derefine=block_to_derefine, block_derefined=block_derefined)
print*, ''
print '(A)', 'blocks index to refine'
print '(A)', trim(str(block_to_refine))
print '(A)', 'blocks index refined'
print '(A)', trim(str(block_refined(2,:)))
call field%adapt(ratio=tree%ratio, block_to_refine=block_to_refine, block_refined=block_refined, &
                 block_to_derefine=block_to_derefine, block_derefined=block_derefined)
print '(A)', 'loop in tree'
do while(tree%loop(node=node))
   call tree%print_code_topology(code=node%code, block_index=.true., coordinates=.true.)
enddo
print*, ''

print '(A)', 'refine level 1'
node => tree%node(code=3_I8P)
node%refinement_needed = NODE_TO_BE_REFINED
call tree%adapt(block_to_refine=block_to_refine, block_refined=block_refined, &
                block_to_derefine=block_to_derefine, block_derefined=block_derefined)
print*, ''
print '(A)', 'blocks index to refine'
print '(A)', trim(str(block_to_refine))
print '(A)', 'blocks refined'
print '(A)', trim(str(block_refined(2,:)))
call field%adapt(ratio=tree%ratio, block_to_refine=block_to_refine, block_refined=block_refined, &
                 block_to_derefine=block_to_derefine, block_derefined=block_derefined)
print '(A)', 'loop in tree'
do while(tree%loop(node=node))
   call tree%print_code_topology(code=node%code, block_index=.true., coordinates=.true.)
enddo
print*, ''
print '(A)', 'refine index level 2'
node => tree%node(code=37_I8P)
node%refinement_needed = NODE_TO_BE_REFINED
call tree%adapt(block_to_refine=block_to_refine, block_refined=block_refined, &
                block_to_derefine=block_to_derefine, block_derefined=block_derefined)
print*, ''

print '(A)', 'blocks index to refine'
print '(A)', trim(str(block_to_refine))
print '(A)', 'blocks index refined'
print '(A)', trim(str(block_refined(2,:)))
call field%adapt(ratio=tree%ratio, block_to_refine=block_to_refine, block_refined=block_refined, &
                 block_to_derefine=block_to_derefine, block_derefined=block_derefined)
print '(A)', 'loop in tree'
do while(tree%loop(node=node))
   call tree%print_code_topology(code=node%code, block_index=.true., coordinates=.true.)
enddo
print*, ''

print '(A)', 'refine index level 3'
node => tree%node(code=307_I8P)
node%refinement_needed = NODE_TO_BE_REFINED
call tree%adapt(block_to_refine=block_to_refine, block_refined=block_refined, &
                block_to_derefine=block_to_derefine, block_derefined=block_derefined)
print*, ''
print '(A)', 'blocks index to refine'
print '(A)', trim(str(block_to_refine))
print '(A)', 'blocks index refined'
print '(A)', trim(str(block_refined(2,:)))
call field%adapt(ratio=tree%ratio, block_to_refine=block_to_refine, block_refined=block_refined, &
                 block_to_derefine=block_to_derefine, block_derefined=block_derefined)
print '(A)', 'loop in tree'
do while(tree%loop(node=node))
   call tree%print_code_topology(code=node%code, block_index=.true., coordinates=.true.)
enddo
print*, ''
call field_save_vtk(tree=tree, field=field, basename='adam-finest')

print*, ''
print '(A)', 'derefine octants 2464-2471'
do c=0, tree%ratio-1
   node => tree%node(code=2464_I8P+c)
   node%refinement_needed = NODE_TO_BE_DEREFINED
enddo
call tree%adapt(block_to_refine=block_to_refine, block_refined=block_refined, &
                block_to_derefine=block_to_derefine, block_derefined=block_derefined)
print*, ''
print '(A)', 'blocks index to derefine'
print '(A)', trim(str(block_to_derefine))
print '(A)', 'blocks index derefined'
print '(A)', trim(str(block_derefined(2,:)))
call field%adapt(ratio=tree%ratio, block_to_refine=block_to_refine, block_refined=block_refined, &
                 block_to_derefine=block_to_derefine, block_derefined=block_derefined)
print '(A)', 'loop in tree'
do while(tree%loop(node=node))
   call tree%print_code_topology(code=node%code, block_index=.true., coordinates=.true.)
enddo
print*, ''

call field_save_vtk(tree=tree, field=field, basename='adam-derefined')
endprogram adam_test_field_object
