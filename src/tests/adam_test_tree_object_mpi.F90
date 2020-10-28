!< ADAM, test tree class.
program adam_test_tree_object_mpi
!< ADAM, test tree class.

use adam_objects
#ifdef _MPI_
use MPI
#endif
use PENF, only : R8P, I8P, I4P, str, strz

implicit none

type(tree_object)               :: tree                 !< Tree.
type(tree_node_object), pointer :: node                 !< Pointer to node.
type(field_object)              :: field                !< Field.
integer(I8P)                    :: code                 !< Tree node code.
integer(I8P), allocatable       :: codes(:)             !< Tree node codes list.
integer(I8P), allocatable       :: block_to_refine(:)   !< List of field blocks to be refined.
integer(I8P), allocatable       :: block_refined(:,:)   !< List of field refined blocks with Morton code.
integer(I8P), allocatable       :: block_to_derefine(:) !< List of field blocks to be derefined.
integer(I8P), allocatable       :: block_derefined(:,:) !< List of field derefined blocks with Morton code.
integer(I4P)                    :: error                !< Error traping flag.
integer(I4P)                    :: myrank               !< Rank of current process.
integer(I8P)                    :: nodes_for_mpi        !< Number of nodes for each MPI process.
integer(I4P)                    :: l, t                 !< Counter.
integer(I4P), allocatable       :: coordinates(:,:)     !< Coordinates (ijkl,nb) of redistributed nodes.

#ifdef _MPI_
call MPI_INIT(error)
call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, error)
#endif

   if (myrank==0) &
   print '(A)', 'initialize tree'
   call tree%initialize
   call field%initialize(nb=20000, emin=[0._R8P,0._R8P,0._R8P], emax=[1._R8P,1._R8P,1._R8P])

   do l=1, 1
      if (myrank==0) &
      print '(A)', 'create children of level '//trim(str(l))
      call tree%mark_all_nodes(mark=NODE_TO_BE_REFINED)
      call tree%adapt(block_to_refine=block_to_refine, block_refined=block_refined, &
                      block_to_derefine=block_to_derefine, block_derefined=block_derefined)
      call field%adapt(ratio=tree%ratio, block_to_refine=block_to_refine, block_refined=block_refined, &
                       block_to_derefine=block_to_derefine, block_derefined=block_derefined)
      call tree%redistribute(coordinates=coordinates)
      call field%redistribute(comm_map_send=tree%comm_map_send,        &
                              comm_map_recv=tree%comm_map_recv,        &
                              comm_map_send_ptr=tree%comm_map_send_ptr,&
                              comm_map_recv_ptr=tree%comm_map_recv_ptr,&
                              local_map=tree%local_map,                &
                              coordinates=coordinates)
   enddo
      ! do code=32_I8P, 39_I8P
      !    node => tree%node(code=code)
      !    node%refinement_needed = NODE_TO_BE_DEREFINED
      ! enddo
      ! call tree%adapt(block_to_refine=block_to_refine, block_refined=block_refined, &
      !                 block_to_derefine=block_to_derefine, block_derefined=block_derefined)
      ! call field%adapt(ratio=tree%ratio, block_to_refine=block_to_refine, block_refined=block_refined, &
      !                  block_to_derefine=block_to_derefine, block_derefined=block_derefined)

      ! call tree%redistribute(coordinates=coordinates)
      ! call field%redistribute(comm_map_send=tree%comm_map_send,        &
      !                         comm_map_recv=tree%comm_map_recv,        &
      !                         comm_map_send_ptr=tree%comm_map_send_ptr,&
      !                         comm_map_recv_ptr=tree%comm_map_recv_ptr,&
      !                         local_map=tree%local_map,                &
      !                         coordinates=coordinates)
   ! ! if (myrank==0) then
      ! print*, 'cazzo send'
      ! print*, trim(str(tree%comm_map_n_send,.true.))
      ! if (allocated(tree%comm_map_send)) &
      ! print*, trim(str(tree%comm_map_send  ,.true.))
      ! print*, 'cazzo recv'
      ! print*, trim(str(tree%comm_map_n_recv,.true.))
      ! if (allocated(tree%comm_map_recv)) &
      ! print*, trim(str(tree%comm_map_recv  ,.true.))
      ! print*, 'cazzo u'
      ! print*, trim(str(int(field%u(1,1,1,:)),.true.))
      ! print*, 'cazzo blocks_number'
      ! print*, trim(str(field%blocks_number,.true.))
      ! print*, 'cazzo local'
      ! if (allocated(tree%local_map)) &
      ! print*, trim(str(tree%local_map(:,1),.true.))
   ! ! endif
! call MPI_FINALIZE(error)
! stop

   if (myrank==0) &
   print '(A)', 'sphere tracking'
   do t=1,2
      if (myrank==0) &
      print '(A)', '  track iteration '//trim(str(t, .true.))//repeat('-',30)
      call mark_sphere_nodes(tree=tree, field=field, center=[0.5_R8P,0.5_R8P,0.5_R8P], radius=0.2_R8P)
      call tree%adapt(block_to_refine=block_to_refine, block_refined=block_refined, &
                      block_to_derefine=block_to_derefine, block_derefined=block_derefined)
      call field%adapt(ratio=tree%ratio, block_to_refine=block_to_refine, block_refined=block_refined, &
                       block_to_derefine=block_to_derefine, block_derefined=block_derefined)
      call tree%redistribute(coordinates=coordinates)
   call MPI_BARRIER(MPI_COMM_WORLD, error)
   print*, 'cazzooiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii'
      print*, 'cazzo n send '//trim(str(tree%comm_map_n_send,.true.))
      print*, 'cazzo n recv '//trim(str(tree%comm_map_n_recv,.true.))
      if (allocated(tree%comm_map_send)) &
      print*, 'cazzo map send '//trim(str(tree%comm_map_send,.true.))
      if (allocated(tree%comm_map_recv)) &
      print*, 'cazzo map recv '//trim(str(tree%comm_map_recv,.true.))
      if (allocated(tree%comm_map_send_ptr)) &
      print*, 'cazzo map send_ptr '//trim(str(tree%comm_map_send_ptr,.true.))
      if (allocated(tree%comm_map_recv_ptr)) &
      print*, 'cazzo map recv_ptr '//trim(str(tree%comm_map_recv_ptr,.true.))
      call field%redistribute(comm_map_send=tree%comm_map_send, &
                              comm_map_recv=tree%comm_map_recv, &
                              comm_map_send_ptr=tree%comm_map_send_ptr,&
                              comm_map_recv_ptr=tree%comm_map_recv_ptr,&
                              local_map=tree%local_map,                &
                              coordinates=coordinates)

   enddo
      call field_save_vtk(tree=tree, field=field, basename='sphere-t-'//trim(strz(t,3)), directory='sphere/')

#ifdef _MPI_
call MPI_FINALIZE(error)
#endif

endprogram adam_test_tree_object_mpi
