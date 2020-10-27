!< ADAM, test tree class.
program adam_test_tree_object_mpi
!< ADAM, test tree class.

use adam_objects
#ifdef _MPI_
use MPI
#endif
use PENF, only : R8P, I8P, I4P, str

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
integer(I4P)                    :: mpi_number           !< Number of MPI processes.
integer(I8P)                    :: nodes_for_mpi        !< Number of nodes for each MPI process.
integer(I4P)                    :: l, n, p, t           !< Counter.

#ifdef _MPI_
call MPI_INIT(error)
call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, error)
call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_number, error)
tree%procs_number = mpi_number
#endif

! if (myrank==0) then
   ! print '(A)', 'initialize tree'
   call tree%initialize
   call field%initialize(nb=190000, emin=[0._R8P,0._R8P,0._R8P], emax=[1._R8P,1._R8P,1._R8P])
   ! print*, ''
   ! print '(A)', 'test codes comparison'
   ! print '(A,L1)', '2  < 3  (T): ', tree%lower(2_I8P, 3_I8P)
   ! print '(A,L1)', '2  > 3  (F): ', tree%greater(2_I8P, 3_I8P)
   ! print '(A,L1)', '3  < 2  (F): ', tree%lower(3_I8P, 2_I8P)
   ! print '(A,L1)', '3  > 2  (T): ', tree%greater(3_I8P, 2_I8P)
   ! print '(A,L1)', '2  < 22 (F): ', tree%lower(2_I8P, 22_I8P)
   ! print '(A,L1)', '2  > 22 (T): ', tree%greater(2_I8P, 22_I8P)
   ! print '(A,L1)', '22 < 2  (T): ', tree%lower(22_I8P, 2_I8P)
   ! print '(A,L1)', '22 > 2  (F): ', tree%greater(22_I8P, 2_I8P)
   ! print*, ''

   do l=1, 2
      ! print '(A)', 'create children of level '//trim(str(l))
      call tree%mark_all_nodes(mark=NODE_TO_BE_REFINED)
      call tree%adapt(block_to_refine=block_to_refine, block_refined=block_refined, &
                      block_to_derefine=block_to_derefine, block_derefined=block_derefined)
      call field%adapt(ratio=tree%ratio, block_to_refine=block_to_refine, block_refined=block_refined, &
                       block_to_derefine=block_to_derefine, block_derefined=block_derefined)
      call tree%redistribute
      ! print '(A)', '  codes: '//trim(str(tree%codes(),.true.))
      ! if (myrank==0) then
      !    do while(tree%loop(node=node))
      !       print '(A)', '  code: '//trim(str(node%code,.true.))// &
      !                    '  myrank: '//trim(str(node%myrank,.true.))// &
      !                    '  myrank_new: '//trim(str(node%myrank_new,.true.))
      !    enddo
      ! endif

      ! node => tree%node(code=6_I8P)
      ! node%refinement_needed = NODE_TO_BE_REFINED
      ! call tree%adapt(block_to_refine=block_to_refine, block_refined=block_refined, &
      !                 block_to_derefine=block_to_derefine, block_derefined=block_derefined)
      ! call field%adapt(ratio=tree%ratio, block_to_refine=block_to_refine, block_refined=block_refined, &
      !                  block_to_derefine=block_to_derefine, block_derefined=block_derefined)
      ! call tree%redistribute
      ! if (myrank==0) then
      !    print '(A)', ''
      !    print '(A)', ''
      !    print '(A)', '  codes: '//trim(str(tree%codes(),.true.))
      !    do while(tree%loop(node=node))
      !       print '(A)', '  code: '//trim(str(node%code,.true.))// &
      !                    '  myrank: '//trim(str(node%myrank,.true.))// &
      !                    '  myrank_new: '//trim(str(node%myrank_new,.true.))
      !    enddo
      ! endif

   enddo

   ! print '(A)', 'sphere tracking'
   do t=1,3
      print*, ''
      print '(A)', '  track iteration '//trim(str(t, .true.))
      call mark_sphere_nodes(tree=tree, field=field, center=[0.2_R8P,0.5_R8P,0.5_R8P], radius=0.1_R8P)
      call tree%adapt(block_to_refine=block_to_refine, block_refined=block_refined, &
                      block_to_derefine=block_to_derefine, block_derefined=block_derefined)
      ! print '(A)', '  nodes refined n. '//trim(str(size(block_refined(1,:), dim=1),.true.))
      ! print '(A)', '  nodes derefined n. '//trim(str(size(block_derefined(1,:), dim=1),.true.))
      call field%adapt(ratio=tree%ratio, block_to_refine=block_to_refine, block_refined=block_refined, &
                       block_to_derefine=block_to_derefine, block_derefined=block_derefined)
      call tree%redistribute
      ! print '(A)', '  codes: '//trim(str(tree%codes(),.true.))
      if (t==3.and.myrank==0) then
         do while(tree%loop(node=node))
            print '(A)', '  code: '//trim(str(node%code,.true.))// &
                         '  myrank: '//trim(str(node%myrank,.true.))// &
                         '  myrank_new: '//trim(str(node%myrank_new,.true.))
         enddo
      endif
   enddo

! endif

! #ifdef _MPI_
! call MPI_BCAST(nodes_for_mpi, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, error)
! #endif
! if (myrank/=0) then
!    allocate(codes(nodes_for_mpi))
!    print '(A)', 'myrank: '//trim(str(myrank))//' nodes for each process: '//trim(str(nodes_for_mpi))
! endif
! #ifdef _MPI_
! call MPI_SCATTER(codes, int(nodes_for_mpi,I4P), MPI_INTEGER8, codes, int(nodes_for_mpi,I4P), MPI_INTEGER8, 0, MPI_COMM_WORLD, error)
! #endif
! if (myrank/=0) then
!    print '(A)', 'myrank: '//trim(str(myrank))//' my codes: '//trim(str(codes))
! endif

#ifdef _MPI_
call MPI_FINALIZE(error)
#endif

endprogram adam_test_tree_object_mpi
