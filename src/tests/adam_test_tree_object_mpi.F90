!< ADAM, test tree class.
program adam_test_tree_object_mpi
!< ADAM, test tree class.

use adam_objects
use PENF
#ifdef _MPI_
use MPI
#endif

implicit none

type(tree_object)               :: tree     !< Tree.
type(tree_node_object), pointer :: node     !< Pointer to node.
type(field_object)              :: field    !< Field.
integer(I8P)                    :: code     !< Tree node code.
integer(I8P), allocatable       :: codes(:) !< Tree node codes list.
integer(I4P)                    :: myrank   !< Rank of current process.
integer(I4P)                    :: l, t, st !< Counter.
#ifdef _MPI_
integer(I4P)                    :: error    !< Error traping flag.

call MPI_INIT(error)
call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, error)
#endif

print '(A)', 'initialize tree'
call tree%initialize
call field%initialize(nb=190000, emin=[0._R8P,0._R8P,0._R8P], emax=[1._R8P,1._R8P,1._R8P])

do l=1, 2
   print '(A)', 'create children of level '//trim(str(l))
   call tree%mark_all_nodes(mark=NODE_TO_BE_REFINED)
   call tree%adapt
   call field%adapt(ratio=tree%ratio,                                                       &
                    block_to_refine=tree%block_to_refine, block_refined=tree%block_refined, &
                    block_to_derefine=tree%block_to_derefine, block_derefined=tree%block_derefined)
   call tree%mpi_redistribute
   call tree%print_mpi_stats
   call field%redistribute(comm_map_send=tree%comm_map_send,        &
                           comm_map_recv=tree%comm_map_recv,        &
                           comm_map_send_ptr=tree%comm_map_send_ptr,&
                           comm_map_recv_ptr=tree%comm_map_recv_ptr,&
                           local_map=tree%local_map,                &
                           coordinates=tree%block_coordinates)
enddo

print '(A)', 'sphere tracking'
do t=1,10
   print '(A)', 'track iteration '//trim(str(t, .true.))//' position x='//trim(str(0.2_R8P + t*0.05_R8P))
   sub_iteration_loop : do st=1, 10
      call mark_sphere_nodes(tree=tree, field=field, center=[0.2_R8P + t*0.05_R8P,0.5_R8P,0.5_R8P], radius=0.1_R8P)
      call tree%adapt
      call field%adapt(ratio=tree%ratio,                         &
                       block_to_refine=tree%block_to_refine,     &
                       block_refined=tree%block_refined,         &
                       block_to_derefine=tree%block_to_derefine, &
                       block_derefined=tree%block_derefined)
      call tree%mpi_redistribute
      call tree%print_mpi_stats
      call field%redistribute(comm_map_send=tree%comm_map_send,         &
                              comm_map_recv=tree%comm_map_recv,         &
                              comm_map_send_ptr=tree%comm_map_send_ptr, &
                              comm_map_recv_ptr=tree%comm_map_recv_ptr, &
                              local_map=tree%local_map,                 &
                              coordinates=tree%block_coordinates)
   enddo sub_iteration_loop
   call field_save_vtk(tree=tree, field=field, basename='sphere-t-'//trim(strz(t,3)), directory='sphere/')
enddo

#ifdef _MPI_
call MPI_FINALIZE(error)
#endif

endprogram adam_test_tree_object_mpi
