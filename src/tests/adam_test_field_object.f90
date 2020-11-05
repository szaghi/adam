!< ADAM, test field class.
program adam_test_field_object
!< ADAM, test field class.

use adam_objects
use PENF

implicit none

type(tree_object)               :: tree             !< The tree.
type(tree_node_object), pointer :: node             !< Pointer to node.
type(field_object)              :: field            !< Field.
real(R8P)                       :: emin(3), emax(3) !< Domain extents.
integer(I4P)                    :: c                !< Counter.

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
call tree%mark_all_nodes(mark=NODE_TO_BE_REFINED)
call tree%adapt
print*, ''
print '(A)', 'blocks index to refine'
print '(A)', trim(str(tree%block_to_refine(1,:), .true.))
print '(A)', 'blocks index refined'
print '(A)', trim(str(tree%block_refined(2,:), .true.))
call field%adapt(ratio=tree%ratio,                                                       &
                 block_to_refine=tree%block_to_refine, block_refined=tree%block_refined, &
                 block_to_derefine=tree%block_to_derefine, block_derefined=tree%block_derefined)
print '(A)', 'loop in tree'
do while(tree%loop(node=node))
   call tree%print_code_topology(code=node%code, block_index=.true., coordinates=.true.)
enddo
print*, ''

print '(A)', 'refine level 1'
node => tree%node(code=3_I8P)
node%refinement_needed = NODE_TO_BE_REFINED
call tree%adapt
print*, ''
print '(A)', 'blocks index to refine'
print '(A)', trim(str(tree%block_to_refine(1,:),.true.))
print '(A)', 'blocks refined'
print '(A)', trim(str(tree%block_refined(2,:),.true.))
call field%adapt(ratio=tree%ratio,                                                       &
                 block_to_refine=tree%block_to_refine, block_refined=tree%block_refined, &
                 block_to_derefine=tree%block_to_derefine, block_derefined=tree%block_derefined)
print '(A)', 'loop in tree'
do while(tree%loop(node=node))
   call tree%print_code_topology(code=node%code, block_index=.true., coordinates=.true.)
enddo
print*, ''
print '(A)', 'refine index level 2'
node => tree%node(code=37_I8P)
node%refinement_needed = NODE_TO_BE_REFINED
call tree%adapt
print*, ''

print '(A)', 'blocks index to refine'
print '(A)', trim(str(tree%block_to_refine(1,:),.true.))
print '(A)', 'blocks index refined'
print '(A)', trim(str(tree%block_refined(2,:), .true.))
call field%adapt(ratio=tree%ratio,                                                       &
                 block_to_refine=tree%block_to_refine, block_refined=tree%block_refined, &
                 block_to_derefine=tree%block_to_derefine, block_derefined=tree%block_derefined)
print '(A)', 'loop in tree'
do while(tree%loop(node=node))
   call tree%print_code_topology(code=node%code, block_index=.true., coordinates=.true.)
enddo
print*, ''

print '(A)', 'refine index level 3'
node => tree%node(code=307_I8P)
node%refinement_needed = NODE_TO_BE_REFINED
call tree%adapt
print*, ''
print '(A)', 'blocks index to refine'
print '(A)', trim(str(tree%block_to_refine(1,:),.true.))
print '(A)', 'blocks index refined'
print '(A)', trim(str(tree%block_refined(2,:), .true.))
call field%adapt(ratio=tree%ratio,                                                       &
                 block_to_refine=tree%block_to_refine, block_refined=tree%block_refined, &
                 block_to_derefine=tree%block_to_derefine, block_derefined=tree%block_derefined)
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
call tree%adapt
print*, ''
print '(A)', 'blocks index to derefine'
print '(A)', trim(str(tree%block_to_derefine(:), .true.))
print '(A)', 'blocks index derefined'
print '(A)', trim(str(tree%block_derefined(2,:), .true.))
call field%adapt(ratio=tree%ratio,                                                       &
                 block_to_refine=tree%block_to_refine, block_refined=tree%block_refined, &
                 block_to_derefine=tree%block_to_derefine, block_derefined=tree%block_derefined)
print '(A)', 'loop in tree'
do while(tree%loop(node=node))
   call tree%print_code_topology(code=node%code, block_index=.true., coordinates=.true.)
enddo
print*, ''

call field_save_vtk(tree=tree, field=field, basename='adam-derefined')
endprogram adam_test_field_object
