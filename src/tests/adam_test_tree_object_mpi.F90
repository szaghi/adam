!< ADAM, test tree class.
program adam_test_tree_object_mpi
!< ADAM, test tree class.

use adam_objects
#ifdef _MPI_
use MPI
#endif
use PENF, only : I8P, I4P, str

implicit none

type(tree_object)               :: tree               !< Tree.
type(tree_node_object), pointer :: tree_node          !< Pointer to node.
integer(I8P)                    :: code               !< Tree node code.
integer(I8P), allocatable       :: codes(:)           !< Tree node codes list.
integer(I8P), allocatable       :: block_to_refine(:) !< List of field blocks to be refined.
integer(I8P), allocatable       :: block_refined(:)   !< List of field refined blocks.
integer(I4P)                    :: error              !< Error traping flag.
integer(I4P)                    :: myrank             !< Rank of current process.
integer(I4P)                    :: mpi_number         !< Number of MPI processes.
integer(I8P)                    :: nodes_for_mpi      !< Number of nodes for each MPI process.
integer(I4P)                    :: l, n, p            !< Counter.

#ifdef _MPI_
call MPI_INIT(error)
call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, error)
call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_number, error)
#endif

if (myrank==0) then
   print '(A)', 'initialize tree'
   call tree%initialize(ratio=8_I4P, max_level=2_I4P)
   print*, ''
   print '(A)', 'test codes comparison'
   print '(A,L1)', '2  < 3  (T): ', tree%lower(2_I8P, 3_I8P)
   print '(A,L1)', '2  > 3  (F): ', tree%greater(2_I8P, 3_I8P)
   print '(A,L1)', '3  < 2  (F): ', tree%lower(3_I8P, 2_I8P)
   print '(A,L1)', '3  > 2  (T): ', tree%greater(3_I8P, 2_I8P)
   print '(A,L1)', '2  < 22 (F): ', tree%lower(2_I8P, 22_I8P)
   print '(A,L1)', '2  > 22 (T): ', tree%greater(2_I8P, 22_I8P)
   print '(A,L1)', '22 < 2  (T): ', tree%lower(22_I8P, 2_I8P)
   print '(A,L1)', '22 > 2  (F): ', tree%greater(22_I8P, 2_I8P)
   print*, ''

   print '(A)', 'test uniform refinement'
   do l=1, 3
      ! third level should not be done because max refinement level has been set to 2
      ! the level 2 should be printed twice
      print '(A)', 'create children of level '//trim(str(l))
      call tree%refine(force_all=.true., block_to_refine=block_to_refine, block_refined=block_refined)
      print '(A)', 'loop in tree'
      do while(tree%loop(code=code))
         call tree%print_code_topology(code=code)
      enddo
      print*, ''
   enddo

   print '(A)', 'test non uniform refinement'
   call tree%initialize(ratio=8_I4P, max_level=2_I4P)
   print '(A)', 'create children of level 1'
   call tree%refine(force_all=.true., block_to_refine=block_to_refine, block_refined=block_refined)
   print '(A)', 'refine nodes 2, 3, 7'
   tree_node => tree%node(code=2_I8P) ; tree_node%refinement_needed = TO_BE_REFINED
   tree_node => tree%node(code=3_I8P) ; tree_node%refinement_needed = TO_BE_REFINED
   tree_node => tree%node(code=7_I8P) ; tree_node%refinement_needed = TO_BE_REFINED
   call tree%refine(block_to_refine=block_to_refine, block_refined=block_refined)
   print '(A)', 'loop in tree'
   do while(tree%loop(code=code))
      call tree%print_code_topology(code=code)
   enddo
   print*, ''
   codes = tree%codes()
   print '(A)', 'list of sorted codes actually stored'
   do code=1, size(codes, dim=1)
      print*, codes(code)
   enddo
   print*, ''

   print '(A)', 'test MPI load balancing'
   print '(A)', 'MPI processes number: '//trim(str(mpi_number))
   print '(A)', 'distribute nodes to other processes'
   nodes_for_mpi = tree%nodes_number / mpi_number
   p = 0
   do n=1, tree%nodes_number, nodes_for_mpi
      if (n+nodes_for_mpi-1<=tree%nodes_number) &
      print '(A)', 'rank: '//trim(str(p))//' codes: '//trim(str(codes(n:n+nodes_for_mpi-1)))
      p = p + 1
   enddo
   print*, ''
endif

#ifdef _MPI_
call MPI_BCAST(nodes_for_mpi, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, error)
#endif
if (myrank/=0) then
   allocate(codes(nodes_for_mpi))
   print '(A)', 'myrank: '//trim(str(myrank))//' nodes for each process: '//trim(str(nodes_for_mpi))
endif
#ifdef _MPI_
call MPI_SCATTER(codes, int(nodes_for_mpi,I4P), MPI_INTEGER8, codes, int(nodes_for_mpi,I4P), MPI_INTEGER8, 0, MPI_COMM_WORLD, error)
#endif
if (myrank/=0) then
   print '(A)', 'myrank: '//trim(str(myrank))//' my codes: '//trim(str(codes))
endif

#ifdef _MPI_
call MPI_FINALIZE(error)
#endif

endprogram adam_test_tree_object_mpi
