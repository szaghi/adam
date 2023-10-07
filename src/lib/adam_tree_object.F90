!< ADAM, tree class definition.
module adam_tree_object
!< ADAM, tree class definition.
!< The tree data structure is organized as a hash table and can be arranged as octree or quadtree simply adopting the proper
!< refinement ratio, namely 8 or 4 ratio respectively.

!< The tree exploits the Morton order for the linearization of the tree. A prototype of quadtree is represented below.
!<
!<            j^
!<    L=0      |
!<      L=1    |
!<        L=2  |
!<          L=3|
!<             |------------------------------------------------
!<          7  |  10   11  |  14   15  |  58   59  |  62   63  |
!<        3----|    10     |    11     |    14     |    15     |
!<          6  |  8    9   |  12   13  |  56   57  |  60   61  |
!<      1------|---------- 2 ----------|---------- 3 ----------|
!<          5  |  34   45  |  38   39  |  50   51  |  54   55  |
!<        2----|     8     |     9     |    12     |    13     |
!<          4  |  32   33  |  36   37  |  48   49  |  52   53  |
!<    0--------|-----------|--------- -1 ----------|-----------|
!<          3  |  10   11  |  14   15  |  26   27  |  30   31  |
!<        1----|     2     |     3     |     6     |     7     |
!<          2  |  8    9   |  12   13  |  24   25  |  28   29  |
!<      0------|---------- 0 ----------|---------- 1 ----------|
!<          1  |  2    3   |  6    7   |  18   19  |  22   23  |
!<        0----|     0     |     1     |     4     |     5     |
!<          0  |  0    1   |  4    5   |  16   17  |  20   21  |
!<             O------------------------------------------------------>
!<                0 |  1   |  2 |  3   |  4 |  5   |  6 |  7     L=3  i
!<                  0      |    1      |    2      |    3        L=2
!<                         0           |           1             L=1
!<                                     0                         L=0
!<
!< In the above representation each refinement level is represented *alone*, namely without taking into account the existence of
!< previous levels. However, the real hierarchy can be taken into account by simply adding the offset derived from the previous
!< levels, thus the numbering becomes as below.
!<
!<            j^
!<    L=0      |
!<      L=1    |
!<        L=2  |
!<          L=3|
!<             |------------------------------------------------
!<          7  |  62   63  |  66   67  |  78   79  |  82   63  |
!<        3----|    14     |    15     |    18     |    19     |
!<          6  |  60   61  |  64   65  |  76   77  |  80   81  |
!<      1------|---------- 2 ----------|---------- 3 ----------|
!<          5  |  54   55  |  58   59  |  70   71  |  74   75  |
!<        2----|    12     |    13     |    16     |    17     |
!<          4  |  52   53  |  56   57  |  68   69  |  72   73  |
!<    0--------|-----------|--------- -1 ----------|-----------|
!<          3  |  30   31  |  34   35  |  46   47  |  50   51  |
!<        1----|     6     |     7     |    10     |    11     |
!<          2  |  28   29  |  32   33  |  44   45  |  48   49  |
!<      0------|---------- 0 ----------|---------- 1 ----------|
!<          1  |  22   23  |  26   27  |  38   39  |  42   43  |
!<        0----|     4     |     5     |     8     |     9     |
!<          0  |  20   21  |  24   25  |  36   37  |  40   41  |
!<             O------------------------------------------------------>
!<                0 |  1   |  2 |  3   |  4 |  5   |  6 |  7     L=3  i
!<                  0      |    1      |    2      |    3        L=2
!<                         0           |           1             L=1
!<                                     0                         L=0
!<
!< This last numbering is the complete Morton order where the Morton code ideintifing a node entails all the spatial information,
!< the refinement level L and the spatial coordinatates IJ.
!< Sometimes it is convenient to use the representation where level L is not encoded into the Morton order, but the conversion is as
!< as adding the offset of previous level, namely using the funcion `first_at_level`.
!< Note that the ancestor node has Morton code -1, it is the ancestor of all nodes, it is **Adam**.
!< For the octree case (3D case) the ordering is equivalent, simply the local children fall in [0,7] instead of [0,3], namely the
!< local numbering is always [0, ratio-1], as represented below.
!<
!<  ^
!< /|\Z
!<  |
!<  |                            G----------*----------H
!<  |                           /|         /|         /|
!<  |                          / |        / |        / |
!<  |                         /  |       /  |       /  |
!<  |                        /   |      /   |      /   |
!<  |                       /    *-----/----*-----/----*
!<  |            6 <-------/--+ /|    /    /|    /   +---------> 7
!<  |                     *----------*----------*    / |
!<  |                    /|   /  |  /|   /  |  /|   /  |
!<  |                   / |  /   | / |  /   | / |  /   |
!<  |                  /  | /    E/--|-/----*/--|-/----F
!<  |         2 <-----/---|/--+ //   |/    //   |/  +---------> 3
!<  |                /    *-----/----*-----/----*    /
!<  |    4 <--------/--+ /|   //    /|   //   +------------> 5
!<  |              C----------*----------D    / |  /
!<  |              |   /  | / |   /  | / |   /  | /
!<  |              |  /   |/  |  /   |/  |  /   |/
!<  |              | /    *---|-/----*---|-/----*
!<  |      0 <-----|/--+ /    |/    /    |/  +------------> 1
!<  |              *----------*----------*    /
!<  |              |   /      |   /      |   /
!<  |              |  /       |  /       |  /
!<  |   _ Y        | /        | /        | /
!<  |   /|         |/         |/         |/
!<  |  /           A----------*----------B
!<  | /
!<  |/                                                    X
!<  o------------------------------------------------------------------->
!<
!< Edges numeration:
!<  AC = fec-7
!<  BD = fec-8
!<  EG = fec-9
!<  FH = fec-10
!<  AE = fec-11
!<  BF = fec-12
!<  CG = fec-13
!<  DH = fec-14
!<  AB = fec-15
!<  EF = fec-16
!<  CD = fec-17
!<  GH = fec-18
!< Corners numeration:
!<  A  = fec-19
!<  B  = fec-20
!<  E  = fec-21
!<  F  = fec-22
!<  C  = fec-23
!<  D  = fec-24
!<  G  = fec-25
!<  H  = fec-26

use adam_grid_object, only : grid_object
use adam_mpih_object, only : mpih_object
use adam_parameters
use adam_tree_node_object, only : tree_node_object
use adam_tree_bucket_object, only : tree_bucket_object, iterator_interface
use FINER, only : file_ini
use FOSSIL
use MORTIF
use PENF
use VECFOR
use MPI
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

implicit none
private
public :: tree_object

! tree node parameters
integer(I4P), parameter :: NODE_LESS_REFINED = -1_I4P      !< More refined node type.
integer(I4P), parameter :: NODE_STANDARD     =  0_I4P      !< Standard node type.
integer(I4P), parameter :: NODE_MORE_REFINED =  1_I4P      !< More refined node type.
integer(I4P), parameter :: NODE_BOUNDARY_CONDITION = 2_I4P !< Boundary condition node type.
! tree parameters
integer(I8P), parameter :: TREE_BUCKETS_NUMBER_DEF = 9973_I8P    !< Default number of buckets of hash table.
real(R8P),    parameter :: TREE_MAX_LOAD = 0.9_R8P               !< Maximum load of hash table buckets.
integer(I4P), parameter :: TREE_MAX_SANITIZE_ITERATIONS = 10_I4P !< Default number of tree sanitize iterations.

type :: tree_object
   !< Tree class definition.
   ! tree data
   type(mpih_object)                     :: mpih                        !< MPI handler.
   type(grid_object), pointer            :: grid=>null()                !< Grid data.
   type(tree_bucket_object), allocatable :: bucket(:)                   !< Tree buckets.
   integer(I8P)                          :: buckets_number=0_I8P        !< Number of buckets used.
   integer(I4P)                          :: nodes_number=0_I4P          !< Number of nodes actually stored, namely the tree length.
   real(R8P)                             :: max_load=TREE_MAX_LOAD      !< Maximum load of tree buckets.
   integer(I4P)                          :: ratio=8_I4P                 !< Refinement ratio.
   integer(I4P)                          :: max_level=12_I4P            !< Maximum refinement level.
   logical                               :: is_initialized_=.false.     !< Initialization status.
   integer(I4P)                          :: ijkl_prune(4)=[-1,-1,-1,-1] !< IJKL prune indexes for simple-initial-forest.
   integer(I4P)                          :: iu_ref_levels=-1_I4P        !< Initial uniform refinement levels.
   ! AMR data
   integer(I8P)              :: n_my_derefine=0_I8P      !< Number of my nodes to be derefined.
   integer(I8P)              :: n_my_refine=0_I8P        !< Number of my nodes to be refined.
   integer(I8P)              :: last_block_index=0_I8P   !< Last block index in the field array.
   integer(I8P), allocatable :: node_to_refine(:)        !< List of nodes to be refined.
   integer(I8P), allocatable :: node_to_derefine(:)      !< List of nodes to be derefined.
   integer(I8P), allocatable :: block_to_refine(:,:)     !< List of field blocks to be refined.
   integer(I8P), allocatable :: block_refined(:,:)       !< List of field refined blocks with Morton code.
   integer(I8P), allocatable :: block_to_derefine(:)     !< List of field blocks to be derefined.
   integer(I8P), allocatable :: block_derefined(:,:)     !< List of field derefined blocks with Morton code.
   integer(I4P), allocatable :: block_coordinates(:,:)   !< Block coordinates of redistributed blocks [4,blocks_number].
   integer(I8P), allocatable :: block_code(:)            !< Block Morton code of redistributed blocks [blocks_number].
   integer(I8P), allocatable :: local_map(:,:)           !< Local map, list block index changes of my nodes.
   integer(I8P), allocatable :: local_map_ghost(:,:)     !< Local map for ghost cells updating [fec_number, 4].
   integer(I8P), allocatable :: local_map_bc_face(:,:)   !< Local map for face BC ghost cells.
   integer(I8P), allocatable :: local_map_bc_edge(:,:)   !< Local map for edge BC ghost cells.
   integer(I8P), allocatable :: local_map_bc_corner(:,:) !< Local map for corner BC ghost cells.
   ! MPI data of nodes
   integer(I4P)              :: my_nodes_number=0_I4P     !< Number of my nodes, keep_nodes_number + recv_nodes_number.
   integer(I4P)              :: send_nodes_number=0_I4P   !< Number of nodes to be sent.
   integer(I4P)              :: recv_nodes_number=0_I4P   !< Number of nodes to be received.
   integer(I4P)              :: keep_nodes_number=0_I4P   !< Number of nodes to be keept.
   integer(I4P)              :: inner_blocks_number=0_I4P !< Number of inner blocks where I need fecs.
   integer(I4P), allocatable :: inner_outer_block_map(:)  !< Inner/outer blocks map.
   integer(I4P), allocatable :: comm_map_n_send(:)        !< Communication map, number of blocks to send [procs_number].
   integer(I4P), allocatable :: comm_map_n_recv(:)        !< Communication map, number of blocks to recv [procs_number].
   integer(I4P), allocatable :: comm_map_send_ptr(:)      !< Communication map, pointers in list to send [procs_number+1].
   integer(I4P), allocatable :: comm_map_recv_ptr(:)      !< Communication map, pointers in list to recv [procs_number+1].
   integer(I8P), allocatable :: comm_map_send(:)          !< Communication map, blocks to send [sum(comm_map_n_send)].
   integer(I8P), allocatable :: comm_map_recv(:)          !< Communication map, blocks to receive [sum(comm_map_n_recv)].
   ! MPI data of ghost cells
   integer(I4P), allocatable :: comm_map_n_send_ghost(:)   !< Communication map, number of ghost celss to send [procs_number].
   integer(I4P), allocatable :: comm_map_n_recv_ghost(:)   !< Communication map, number of ghost celss to recv [procs_number].
   integer(I4P), allocatable :: comm_map_send_ptr_ghost(:) !< Communication map, pointers in list to send [procs_number+1].
   integer(I4P), allocatable :: comm_map_recv_ptr_ghost(:) !< Communication map, pointers in list to recv [procs_number+1].
   integer(I8P), allocatable :: comm_map_send_ghost(:,:)   !< Communication map, `fec` information [fec_number, 5].
   integer(I8P), allocatable :: comm_map_recv_ghost(:,:)   !< Communication map, `fec` information [fec_number, 5].
   ! STL surfaces data
   type(surface_stl_object)  :: surface_stl !< STL surface.

   integer(I8P), allocatable :: temp_array_i8(:)
   integer(I8P), allocatable :: codes_analyzed(:)    !< List of codes analyzed.
   contains
      ! public methods
      procedure, pass(self) :: adapt                        !< Adapt tree accordingly to refine/derefine necessity.
      procedure, pass(self) :: blocks_reorder               !< Reorder blocks indexes in field.
      procedure, pass(self) :: codes                        !< Return the list of (sorted) codes actually stored in the tree.
      procedure, pass(self) :: compute_surface_stl_distance !< Compute signed distance of nodes from a STL surface.
      procedure, pass(self) :: get_closest_block            !< Get the closest block to a given point.
      procedure, pass(self) :: get_closest_cells            !< Get the closest cells to a given point.
      procedure, pass(self) :: hash                         !< Hash the key.
      procedure, pass(self) :: has_code                     !< Check if the code is present in the tree.
      procedure, pass(self) :: initialize                   !< Initialize the tree.
      procedure, pass(self) :: load_nodes                   !< Load nodes data, used for restart.
      procedure, pass(self) :: load_from_ini_file           !< Load object data from INI file.
      procedure, pass(self) :: load_surface_stl             !< Load surface from STL file.
      procedure, pass(self) :: loop                         !< Sentinel while-loop on nodes returning the code.
      procedure, pass(self) :: make_local_maps_bc           !< Make local maps of boundary conditions.
      procedure, pass(self) :: make_neighborhood            !< Make neighborhood all whole tree and store it in nodes.
      procedure, pass(self) :: mark_all_nodes               !< Mark all nodes to be refined, derefined, ecc.
      procedure, pass(self) :: mark_sphere                  !< Mark nodes to be refined/derefined by sphere distance.
      procedure, pass(self) :: mark_surface_stl             !< Mark all nodes inside a surface defined by STL triangulation.
      procedure, pass(self) :: max_cell_delta               !< Return the maximum cell delta given a comparison distance.
      procedure, pass(self) :: node                         !< Return a pointer to a node.
      procedure, pass(self) :: prime_buckets_number         !< Return the buckets number as nearest prime number given nodes number.
      procedure, pass(self) :: print_status                 !< Print status of main data.
      procedure, pass(self) :: prune                        !< Prune nodes.
      procedure, pass(self) :: resize                       !< Resize the tree.
      procedure, pass(self) :: save_nodes                   !< Save nodes data, used for restart.
      procedure, pass(self) :: save_local_map               !< Save local map.
      procedure, pass(self) :: traverse                     !< Traverse tree calling the iterator procedure.
      ! MPI methods
      procedure, pass(self) :: import_refinements_needed  !< Import refinements needed status changed externally.
      procedure, pass(self) :: make_comm_local_maps       !< Make communication/local maps.
      procedure, pass(self) :: make_comm_local_maps_ghost !< Make communication/local maps of ghost cells.
      procedure, pass(self) :: mpi_gather_nodes_data      !< Gather nodes data between MPI processes.
      procedure, pass(self) :: mpi_print_stats            !< Print MPI stats.
      procedure, pass(self) :: mpi_redistribute           !< Redistribute nodes to MPI processes, load balancing.
      ! Morton ordering methods
      generic               :: coordinates_to_morton => &
                               coordinates3D_to_morton, &
                               coordinates2D_to_morton !< Return the Morton code given space-level coordinates.
      generic               :: morton_to_coordinates => &
                               morton_to_coordinates3D, &
                               morton_to_coordinates2D !< Return the space-level coordinates given Morton code.
      procedure, pass(self) :: all_siblings            !< Return all siblings Morton code given Morton code.
      procedure, pass(self) :: child                   !< Return the i-th child given Morton code.
      procedure, pass(self) :: child_local             !< Return the child index in the local numbering.
      procedure, pass(self) :: children                !< Return the children list given Morton code.
      procedure, pass(self) :: finest_at_level         !< Return the finest node code at given level.
      procedure, pass(self) :: first_at_level          !< Return the first node code at given level.
      procedure, pass(self) :: first_common_parent     !< Return the first common parent given two codes.
      procedure, pass(self) :: get_neighbor            !< Return the neighbor in a given face of given Morton code.
      procedure, pass(self) :: get_neighbor_all        !< Return the neighbor in a given face/edge/corner of given Morton code.
      procedure, pass(self) :: greater                 !< Return true if code is greater than other.
      procedure, pass(self) :: is_greater_by_level     !< Return true if code is greater than other, comparing level.
      procedure, pass(self) :: is_greater_by_pos       !< Return true if code is greater than other, comparing position
      procedure, pass(self) :: is_lower_by_level       !< Return true if code is lower than other, comparing level.
      procedure, pass(self) :: is_lower_by_pos         !< Return true if code is lower than other, comparing position.
      procedure, pass(self) :: last_at_level           !< Return the last node code at given level.
      procedure, pass(self) :: level                   !< Return the refinement level given the code.
      procedure, pass(self) :: lower                   !< Return true if code is lower than other.
      procedure, pass(self) :: parent                  !< Return the parent given Morton code.
      procedure, pass(self) :: parent_at_level         !< Return the parent given Morton code at given level.
      procedure, pass(self) :: path                    !< Return the path codes, the list of codes from given node to root.
      procedure, pass(self) :: print_code_topology     !< Print all code topology data.
      procedure, pass(self) :: siblings                !< Return the siblings Morton code given Morton code.
      ! private methods
      procedure, pass(self), private :: add_node                !< Add a node pointer to the tree.
      procedure, pass(self), private :: coordinates3D_to_morton !< Return the Morton code given ijkl coordinates.
      procedure, pass(self), private :: coordinates2D_to_morton !< Return the Morton code given ijl coordinates.
      procedure, pass(self), private :: derefine                !< Derefine nodes.
      procedure, pass(self), private :: empty                   !< Empty tree, i.e. remove all nodes, without destroy the tree.
      procedure, pass(self), private :: morton_to_coordinates3D !< Return the ijkl coordinates given Morton code.
      procedure, pass(self), private :: morton_to_coordinates2D !< Return the ijkl coordinates given Morton code.
      procedure, pass(self), private :: refine                  !< Refine nodes.
      procedure, pass(self), private :: remove_node             !< Remove a node from the tree, given the key.
      procedure, pass(self), private :: sanitize                !< Sanitize the tree.
endtype tree_object

abstract interface
   subroutine is_greater_interface(self, lhs, rhs, res)
   !< Return true if code is greater than other.
   import :: tree_object, I8P
   class(tree_object), intent(in)  :: self !< The tree.
   integer(I8P),       intent(in)  :: lhs  !< Left hand side of code comparison.
   integer(I8P),       intent(in)  :: rhs  !< Right hand side of code comparison.
   logical,            intent(out) :: res  !< Comparison result.
   endsubroutine is_greater_interface

   subroutine is_lower_interface(self, lhs, rhs, res)
   !< Return true if code is lower than other.
   import :: tree_object, I8P
   class(tree_object), intent(in)  :: self !< The tree.
   integer(I8P),       intent(in)  :: lhs  !< Left hand side of code comparison.
   integer(I8P),       intent(in)  :: rhs  !< Right hand side of code comparison.
   logical,            intent(out) :: res  !< Comparison result.
   endsubroutine is_lower_interface
endinterface

contains
   ! public methods
   subroutine adapt(self)
   !< Adapt tree accordingly to refine/derefine necessity.
   class(tree_object), intent(inout) :: self        !< The tree.
   real(R8P)                         :: timing(1:2) !< Tic toc timing.

   timing(1) = MPI_WTIME()
   call self%sanitize
   timing(2) = MPI_WTIME()
   ! print '(A, F18.10)', self%mpih%myrankstr//'step timing TREE%ADAPT%SANITIZE: ', timing(2) - timing(1)

   timing(1) = MPI_WTIME()
   call self%refine
   timing(2) = MPI_WTIME()
   ! print '(A, F18.10)', self%mpih%myrankstr//'step timing TREE%ADAPT%REFINE: ', timing(2) - timing(1)

   timing(1) = MPI_WTIME()
   call self%derefine
   timing(2) = MPI_WTIME()
   ! print '(A, F18.10)', self%mpih%myrankstr//'step timing TREE%ADAPT%DEREFINE: ', timing(2) - timing(1)

   timing(1) = MPI_WTIME()
   call self%make_neighborhood
   timing(2) = MPI_WTIME()
   ! print '(A, F18.10)', self%mpih%myrankstr//'step timing TREE%ADAPT%MAKE_NEIGHBORHOOD: ', timing(2) - timing(1)
   endsubroutine adapt

   function codes(self, only_mine, sort_by_level)
   !< Return the list of (sorted) codes actually stored in the tree.
   class(tree_object), intent(in)           :: self           !< The tree.
   logical,            intent(in), optional :: only_mine      !< If true return only the nodes of myrank process.
   logical,            intent(in), optional :: sort_by_level  !< If true sort codes be level instead of position.
   integer(I8P), allocatable                :: codes(:)       !< List of codes.
   logical                                  :: only_mine_     !< If true return only the nodes of myrank process.
   logical                                  :: sort_by_level_ !< If true sort codes be level instead of position.
   type(tree_node_object), pointer          :: node_ptr       !< Pointer to current node.
   integer(I8P)                             :: c              !< Counter.
   integer(I8P), allocatable                :: work(:)        !< Working memory for sorting codes list.

   only_mine_ = .false. ; if (present(only_mine)) only_mine_ = only_mine
   sort_by_level_ = .false. ; if (present(sort_by_level)) sort_by_level_ = sort_by_level
   allocate(codes(self%nodes_number))
   c = 0
   do while(self%loop(node_ptr=node_ptr))
      if (only_mine_.and.self%mpih%myrank/=node_ptr%myrank) cycle
      c = c + 1
      codes(c) = node_ptr%code
   enddo
   if (c < self%nodes_number) then
      work = codes(1:c)
      !RIMETTERE SENZA
      if(allocated(codes)) deallocate(codes)
      !RIMETTERE SENZA
      call move_alloc(from=work, to=codes)
   endif
   allocate(work((size(codes,dim=1)+1)/2))
   ! sort codes
   if (sort_by_level_) then
      call mergesort(array=codes, is_lower=is_lower_by_level, is_greater=is_greater_by_level)
   else
      call mergesort(array=codes, is_lower=is_lower_by_pos,   is_greater=is_greater_by_pos  )
   endif
   contains
      recursive subroutine mergesort(array, is_lower, is_greater)
      !< Sort input array by means of mergesort algorithm.
      integer(I8P), intent(inout)     :: array(:)   !< Array to be sorted.
      procedure(is_lower_interface)   :: is_lower   !< Lower comparison function.
      procedure(is_greater_interface) :: is_greater !< Greater comparison function.
      integer(I8P)                    :: half       !< Half size counter.
      logical                         :: is_great   !< Result of is greater.

      half = (size(array) + 1) / 2
      if (size(array) < 2) then
         continue
      else if (size(array) == 2) then
         ! if (self%greater(array(1), array(2))) call swap_element(array(1), array(2))
         call is_greater(self=self, lhs=array(1), rhs=array(2), res=is_great)
         if (is_great) call swap_element(array(1), array(2))
      else
         call mergesort(array=array( : half),    is_lower=is_lower, is_greater=is_greater)
         call mergesort(array=array(half + 1 :), is_lower=is_lower, is_greater=is_greater)
         ! if (self%greater(array(half), array(half + 1))) then
         call is_greater(self=self, lhs=array(half), rhs=array(half + 1), res=is_great)
         if (is_great) then
            work(1 : half) = array(1 : half)
            call merge_array(A=work(1 : half), B=array(half + 1:), C=array, is_lower=is_lower)
         endif
      endif
      endsubroutine mergesort

      subroutine merge_array(A, B, C, is_lower)
      !< Merge arrays A/B in C.
      integer(I8P), target, intent(in)    :: A(:), B(:) !< Input arrays.
      integer(I8P), target, intent(inout) :: C(:)       !< Output array.
      procedure(is_lower_interface)       :: is_lower   !< Lower comparison function.
      logical                             :: is_low     !< Result of is greater.
      integer(I8P)                        :: i, j, k    !< Counter.

      if (size(A) + size(B) > size(C)) stop

      i = 1 ; j = 1
      do k = 1, size(C)
         if (i <= size(A) .and. j <= size(B)) then
            ! if (self%lower(A(i), B(j))) then
            call is_lower(self=self, lhs=A(i), rhs=B(j), res=is_low)
            if (is_low) then
               C(k) = A(i)
               i = i + 1
            else
               C(k) = B(j)
               j = j + 1
            endif
         else if (i <= size(A)) then
            C(k) = A(i)
            i = i + 1
         else if (j <= size(B)) then
            C(k) = B(j)
            j = j + 1
         endif
      enddo
      endsubroutine merge_array

      subroutine swap_element(x, y)
      !< Swap array element.
      integer(I8P), intent(inout) :: x, y !< Array element to be swaped.
      integer(I8P)                :: tmp  !< Temporary memory.

      tmp = x ; x = y ; y = tmp
      endsubroutine swap_element
   endfunction codes

   subroutine compute_surface_stl_distance(self, surface_stl, from_cell, cell_distance)
   !< Compute signed distance of blocks/cells from a STL surface.
   !<
   !< Distance from blocks are computed for all blocks, there is not need to exchange data between processes.
   !< On the contrary, distance from cells is computed only for my blocks, each process must pass the distance
   !< to each field.
   class(tree_object),       intent(inout)         :: self                   !< The tree.
   type(surface_stl_object), intent(in)            :: surface_stl            !< STL surface.
   logical,                  intent(in),  optional :: from_cell              !< Distance from cells instead of blocks.
   real(R8P), allocatable,   intent(out), optional :: cell_distance(:,:,:,:) !< Distance from cells.
   logical                                         :: from_cell_             !< Distance from cells instead of blocks, local var.
   type(tree_node_object), pointer                 :: node_ptr               !< Pointer to current node.
   real(R8P)                                       :: block_center(3)        !< block center coordinates.
   real(R8P)                                       :: distance(0:8)          !< Distances between block and sphere.
   integer(I4P)                                    :: i,j,k,l                !< Counter.
   real(R8P)                                       :: emin(3), emax(3)       !< Node extents.
   type(vector_R8P)                                :: point(0:8)             !< Vector point coordinates.
   real(R8P), allocatable                          :: x_cell(:)              !< X cell coordinates.
   real(R8P), allocatable                          :: y_cell(:)              !< Y cell coordinates.
   real(R8P), allocatable                          :: z_cell(:)              !< Z cell coordinates.

   from_cell_ = .false. ; if (present(from_cell)) from_cell_ = from_cell
   associate(ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk, ngc=>self%grid%ngc)
   if (from_cell_.and.present(cell_distance)) then
      ! compute distance from cell, distance from blocks must be already computed
      allocate(cell_distance(1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:self%my_nodes_number))
      ! cell_distance = huge(0._R8P)
      cell_distance = 10._R8P
      allocate(x_cell(1-ngc:ni+ngc))
      allocate(y_cell(1-ngc:nj+ngc))
      allocate(z_cell(1-ngc:nk+ngc))
      do while(self%loop(node_ptr=node_ptr))
         if (node_ptr%myrank==self%mpih%myrank) then
            if (node_ptr%surface_stl_distance<epsilon(0._R8P)) then
               ! compute cell distance only in blocks where is STL surface
               call self%morton_to_coordinates(code=node_ptr%code, i=i, j=j, k=k, l=l)
               call self%grid%compute_metrics(coordinates=[i,j,k,l], x_cell=x_cell, y_cell=y_cell, z_cell=z_cell)
               do k=1-ngc, nk+ngc
                  do j=1-ngc, nj+ngc
                     do i=1-ngc, ni+ngc
                        point(0) = x_cell(i) * ex_R8P + y_cell(j) * ey_R8P + z_cell(k) * ez_R8P
                        cell_distance(i,j,k,node_ptr%block_index) = surface_stl%distance(point=point(0),   &
                                                                                         is_signed=.true., &
                                                                                         sign_algorithm='ray_intersections')
                     enddo
                  enddo
               enddo
            endif
         endif
      enddo
   else
      do while(self%loop(node_ptr=node_ptr))
         if (node_ptr%myrank==self%mpih%myrank) then
            if (node_ptr%surface_stl_distance<huge(0._R8P)/2._R8P) cycle ! distance already computed for this node
            call self%morton_to_coordinates(code=node_ptr%code, i=i, j=j, k=k, l=l)
            call self%grid%compute_metrics(coordinates=[i,j,k,l], emin=emin, emax=emax)
            block_center = (emax + emin) / 2._R8P
            point(0) = block_center(1) * ex_R8P +  block_center(2) * ey_R8P + block_center(3) * ez_R8P
            point(1) =         emin(1) * ex_R8P +          emin(2) * ey_R8P +         emin(3) * ez_R8P
            point(2) =         emax(1) * ex_R8P +          emin(2) * ey_R8P +         emin(3) * ez_R8P
            point(3) =         emin(1) * ex_R8P +          emax(2) * ey_R8P +         emin(3) * ez_R8P
            point(4) =         emax(1) * ex_R8P +          emax(2) * ey_R8P +         emin(3) * ez_R8P
            point(5) =         emin(1) * ex_R8P +          emin(2) * ey_R8P +         emax(3) * ez_R8P
            point(6) =         emax(1) * ex_R8P +          emin(2) * ey_R8P +         emax(3) * ez_R8P
            point(7) =         emin(1) * ex_R8P +          emax(2) * ey_R8P +         emax(3) * ez_R8P
            point(8) =         emax(1) * ex_R8P +          emax(2) * ey_R8P +         emax(3) * ez_R8P
            distance(0) = surface_stl%distance(point=point(0), is_signed=.true., sign_algorithm='ray_intersections')
            distance(1) = surface_stl%distance(point=point(1), is_signed=.true., sign_algorithm='ray_intersections')
            distance(2) = surface_stl%distance(point=point(2), is_signed=.true., sign_algorithm='ray_intersections')
            distance(3) = surface_stl%distance(point=point(3), is_signed=.true., sign_algorithm='ray_intersections')
            distance(4) = surface_stl%distance(point=point(4), is_signed=.true., sign_algorithm='ray_intersections')
            distance(5) = surface_stl%distance(point=point(5), is_signed=.true., sign_algorithm='ray_intersections')
            distance(6) = surface_stl%distance(point=point(6), is_signed=.true., sign_algorithm='ray_intersections')
            distance(7) = surface_stl%distance(point=point(7), is_signed=.true., sign_algorithm='ray_intersections')
            distance(8) = surface_stl%distance(point=point(8), is_signed=.true., sign_algorithm='ray_intersections')
            if (maxval(distance(0:8),dim=1)*minval(distance(0:8),dim=1) < 0._R8P) distance(0) = 0._R8P
            node_ptr%surface_stl_distance = distance(0)
         endif
      enddo
   endif
   endassociate
   endsubroutine compute_surface_stl_distance

   function get_closest_block(self, point) result(code)
   !< Get the closest block to a given point.
   class(tree_object), intent(in) :: self        !< The tree.
   real(R8P),          intent(in) :: point(3)    !< Point xyz coordinates.
   integer(I8P)                   :: code        !< The Morton code.
   integer(I4P)                   :: ijkl(4)     !< Indexes.
   integer(I8P), allocatable      :: path(:)     !< Path codes from node to root.
   integer(I4P)                   :: c           !< Counter.
   logical                        :: block_found !< Flag to check if block has been found.

   ! find the block closest to the point in the finest refinement level
   ijkl(1:3) = self%grid%get_closest_block(point=point, level=self%max_level)
   code = self%coordinates_to_morton(i=ijkl(1), j=ijkl(2), k=ijkl(3), l=self%max_level)
   if (.not.self%has_code(code=code)) then
      ! finest-level block does not exist, return the finest living into its parents-path
      path = self%path(code=code)
      block_found = .false.
      do c=lbound(path, dim=1), ubound(path, dim=1)
         code = path(c)
         if (self%has_code(code=code)) then
            block_found = .true.
            exit
          endif
      enddo
      if (.not.block_found) then
         print '(A)', self%mpih%myrankstr//'ERROR: tree%get_closest block failed, path: '//str(path)//' point: '//str(point)
      endif
   endif
   endfunction get_closest_block

   subroutine get_closest_cells(self, point, code, ijk, v, xyz)
   !< Get the closest cells to a given point.
   class(tree_object), intent(in)            :: self       !< The tree.
   real(R8P),          intent(in)            :: point(3)   !< Point xyz coordinates.
   integer(I8P),       intent(in)            :: code       !< The Morton code of the closest block.
   integer(I4P),       intent(out)           :: ijk(3,8)   !< Closest cells indexes.
   integer(I4P),       intent(out), optional :: v          !< Closest vertex index.
   real(R8P),          intent(out), optional :: xyz(3,8)   !< Closest cells center-coordinates.
   integer(I4P)                              :: v_         !< Closest vertex index, local var.
   integer(I4P)                              :: ijkl(4)    !< Block indexes.
   real(R8P)                                 :: emin(3)    !< Minimum block extent.
   real(R8P)                                 :: dxyz(3)    !< Space steps.
   integer(I4P)                              :: vijk(3)    !< Vertex indices.
   integer(I4P)                              :: c,cc,i,j,k !< Counter.
   ! parameters
   ! integer(I4P), parameter :: vi(0:1,0:1,0:1) = reshape([0,1,2,3,4,5,6,7],[2,2,2]) !< Vertices indexes.
   ! integer(I4P), parameter :: ci(0:7)         = [-1,+1,-1,+1,-1,+1,-1,+1]          !< Cell-vertices x increments.
   ! integer(I4P), parameter :: cj(0:7)         = [-1,-1,+1,+1,-1,-1,+1,+1]          !< Cell-vertices y increments.
   ! integer(I4P), parameter :: ck(0:7)         = [-1,-1,-1,-1,+1,+1,+1,+1]          !< Cell-vertices z increments.

   ! get metrics of the closest block
   call self%morton_to_coordinates(code=code, i=ijkl(1), j=ijkl(2), k=ijkl(3), l=ijkl(4))
   call self%grid%compute_metrics(coordinates=ijkl, emin=emin, dx=dxyz(1), dy=dxyz(2), dz=dxyz(3))
   ijk(1,1) = floor((point(1) - 0.5_R8P*dxyz(1) - emin(1)) / dxyz(1), I4P) + 1
   ijk(2,1) = floor((point(2) - 0.5_R8P*dxyz(2) - emin(2)) / dxyz(2), I4P) + 1
   ijk(3,1) = floor((point(3) - 0.5_R8P*dxyz(3) - emin(3)) / dxyz(3), I4P) + 1
   cc = 0
   do k=0,1
     do j=0,1
       do i=0,1
          cc = cc + 1
          ijk(1,cc) = ijk(1,1) + i
          ijk(2,cc) = ijk(2,1) + j
          ijk(3,cc) = ijk(3,1) + k
       enddo
     enddo
   enddo
   ! find closest cell, put it into first element of ijk array
   ! ijk(1,1) = ceiling((point(1) - emin(1)) / dxyz(1), I4P)
   ! ijk(2,1) = ceiling((point(2) - emin(2)) / dxyz(2), I4P)
   ! ijk(3,1) = ceiling((point(3) - emin(3)) / dxyz(3), I4P)
   ! ! find closest vertex into the closest cell
   ! vijk(1) = nint((point(1) - (emin(1) + (ijk(1,1)-1) * dxyz(1))) / dxyz(1), I4P)
   ! vijk(2) = nint((point(2) - (emin(2) + (ijk(2,1)-1) * dxyz(2))) / dxyz(2), I4P)
   ! vijk(3) = nint((point(3) - (emin(3) + (ijk(3,1)-1) * dxyz(3))) / dxyz(3), I4P)
   ! v_= vi(vijk(1), vijk(2), vijk(3))
   ! ! find other 7 cells surrounding the vertex
   !                                ijk(1,2) = ijk(1,1) + ci(v_) ;  ijk(1,3) = ijk(1,1)          ;  ijk(1,4) = ijk(1,1) + ci(v_)
   !                                ijk(2,2) = ijk(2,1)          ;  ijk(2,3) = ijk(2,1) + cj(v_) ;  ijk(2,4) = ijk(2,1) + cj(v_)
   !                                ijk(3,2) = ijk(3,1)          ;  ijk(3,3) = ijk(3,1)          ;  ijk(3,4) = ijk(3,1)

   ! ijk(1,5) = ijk(1,1)          ; ijk(1,6) = ijk(1,1) + ci(v_) ;  ijk(1,7) = ijk(1,1)          ;  ijk(1,8) = ijk(1,1) + ci(v_)
   ! ijk(2,5) = ijk(2,1)          ; ijk(2,6) = ijk(2,1)          ;  ijk(2,7) = ijk(2,1) + cj(v_) ;  ijk(2,8) = ijk(2,1) + cj(v_)
   ! ijk(3,5) = ijk(3,1) + ck(v_) ; ijk(3,6) = ijk(3,1) + ck(v_) ;  ijk(3,7) = ijk(3,1) + ck(v_) ;  ijk(3,8) = ijk(3,1) + ck(v_)
   ! if (present(v)) v = v_
   if (present(xyz)) then
      do cc=1,8
         do c=1,3
            xyz(c,cc) = emin(c) + (ijk(c,cc)-0.5_R8P) * dxyz(c)
         enddo
      enddo
   endif
   endsubroutine get_closest_cells

   function has_code(self, code)
   !< Check if the key is present in the tree.
   class(tree_object), intent(in) :: self     !< The tree.
   integer(I8P),       intent(in) :: code     !< The Morton code.
   logical                        :: has_code !< Check result.

   has_code = .false.
   if (self%is_initialized_) has_code = self%bucket(self%hash(code=code))%has_code(code=code)
   endfunction has_code

   elemental function hash(self, code) result(bucket)
   !< Hash the key.
   class(tree_object), intent(in) :: self   !< The tree.
   integer(I8P),       intent(in) :: code   !< The Morton code.
   integer(I4P)                   :: bucket !< Bucket index corresponding to the key.

   bucket = modulo(code, int(self%buckets_number, I8P)) + 1
   endfunction hash

   subroutine initialize(self, grid, file_parameters, max_load, nodes_number, buckets_number, ratio, max_level, add_adam, &
                         iu_ref_levels, i_prune, j_prune, k_prune, l_prune)
   !< Initialize the tree.
   class(tree_object), intent(inout)           :: self            !< The tree.
   type(grid_object),  intent(in), target      :: grid            !< Grid data.
   type(file_ini),     intent(inout), optional :: file_parameters !< INI file handler.
   real(R8P),          intent(in),    optional :: max_load        !< Maximum load of tree buckets.
   integer(I8P),       intent(in),    optional :: nodes_number    !< Nodes number to be stored in the tree.
   integer(I8P),       intent(in),    optional :: buckets_number  !< Number of buckets for initialize the tree.
   integer(I4P),       intent(in),    optional :: ratio           !< Refinement ratio.
   integer(I4P),       intent(in),    optional :: max_level       !< Maximum refinement level.
   integer(I4P),       intent(in),    optional :: iu_ref_levels   !< Uniform initial refinement.
   integer(I4P),       intent(in),    optional :: i_prune         !< Pruning along x.
   integer(I4P),       intent(in),    optional :: j_prune         !< Pruning along y.
   integer(I4P),       intent(in),    optional :: k_prune         !< Pruning along z.
   integer(I4P),       intent(in),    optional :: l_prune         !< Pruning level.
   logical,            intent(in),    optional :: add_adam        !< Add ADAM node, the ancestor of all nodes.
   logical                                     :: add_adam_       !< Add ADAM node, the ancestor of all nodes, local var.

   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'tree%initialize start'
   self%grid => grid
   if (present(file_parameters)) call self%load_from_ini_file(file_parameters)

   ! parameters explicitely passed ovveride ones file-passed
   add_adam_ = .true. ; if (present(add_adam)) add_adam_ = add_adam
   ! tree data
   if (present(max_load)) self%max_load = max_load
   if (present(nodes_number)) then
      self%buckets_number = self%prime_buckets_number(nodes_number=nodes_number)
   else
      self%buckets_number = TREE_BUCKETS_NUMBER_DEF ; if (present(buckets_number)) self%buckets_number = buckets_number
   endif
   allocate(self%bucket(1:self%buckets_number))
   ! allocate(self%bucket(1:4133))
   if (present(ratio)) then
      if (ratio==8_I8P.or.ratio==4_I8P) then
         self%ratio = ratio
      else
         write(stderr, '(A)') self%mpih%myrankstr//'ADAM-ERROR: tree ratio must be 8 o 4'
         call MPI_FINALIZE(self%mpih%error)
         stop
      endif
   endif
   if (present(max_level)) self%max_level = max_level
   self%is_initialized_ = .true.
   if (add_adam_) call self%add_node(code=-1_I8P) ! add ADAM node, the ancestor of all nodes

   if(present(iu_ref_levels)) self%iu_ref_levels = iu_ref_levels
   if(present(i_prune)) self%ijkl_prune(1)       = i_prune
   if(present(j_prune)) self%ijkl_prune(2)       = j_prune
   if(present(k_prune)) self%ijkl_prune(3)       = k_prune
   if(present(l_prune)) self%ijkl_prune(4)       = l_prune

   ! MPI data
   allocate(self%comm_map_n_send(0:self%mpih%procs_number-1))
   allocate(self%comm_map_n_recv(0:self%mpih%procs_number-1))
   allocate(self%comm_map_send_ptr(0:self%mpih%procs_number))
   allocate(self%comm_map_recv_ptr(0:self%mpih%procs_number))
   self%comm_map_n_send = 0_I4P
   self%comm_map_n_recv = 0_I4P
   self%comm_map_send_ptr = 0_I4P
   self%comm_map_recv_ptr = 0_I4P
   call self%mpi_redistribute
   ! MPI ghost data
   allocate(self%comm_map_n_send_ghost(0:self%mpih%procs_number-1))
   allocate(self%comm_map_n_recv_ghost(0:self%mpih%procs_number-1))
   allocate(self%comm_map_send_ptr_ghost(0:self%mpih%procs_number))
   allocate(self%comm_map_recv_ptr_ghost(0:self%mpih%procs_number))

   call self%make_neighborhood
   print '(A)', self%mpih%myrankstr//'tree%initialize finish'
   endsubroutine initialize

   subroutine load_nodes(self, file_name)
   !< Load nodes data, used for restart.
   !<
   !< Note: the tree is made empty before lading nodes data.
   class(tree_object), intent(inout) :: self             !< The tree.
   character(*),       intent(in)    :: file_name        !< Output file name.
   integer(I8P)                      :: codes_number     !< Number of codes actually stored.
   integer(I4P)                      :: file_unit        !< Output file unit.
   logical                           :: file_exist       !< Flag to check file's existance.
   integer(I8P)                      :: node_code        !< Morton code of current node.
   integer(I4P)                      :: node_rank        !< MPI rank of current node.
   integer(I8P)                      :: node_block_index !< Block index of current node.
   integer(I8P)                      :: c                !< Counter.

   inquire(file=trim(adjustl(file_name)), exist=file_exist)
   if (file_exist) then
      print '(A)', self%mpih%myrankstr//'load tree nodes from file '//trim(adjustl(file_name))
      open(newunit=file_unit, file=trim(adjustl(file_name)), form='UNFORMATTED', access='STREAM')
      read(unit=file_unit) codes_number
      print '(A)', self%mpih%myrankstr//'tree nodes number '//trim(str(codes_number))
      if (codes_number > 0) then
         call self%empty
         do c=1, codes_number
            read(unit=file_unit) node_code, node_rank, node_block_index
            call self%add_node(code=node_code, rank=node_rank, block_index=node_block_index)
         enddo
      endif
      close(file_unit)
      call self%make_neighborhood
      print '(A)', self%mpih%myrankstr//'load tree nodes from file '//trim(adjustl(file_name))//' completed'
   else
      write(stderr, '(A)') self%mpih%myrankstr//'ERROR: file "'//trim(adjustl(file_name))//'" does not exist!'
   endif
   endsubroutine load_nodes

   subroutine load_from_ini_file(self, file_parameters)
   !< Load object data from INI file.
   class(tree_object), intent(inout) :: self            !< The tree.
   type(file_ini),     intent(inout) :: file_parameters !< INI file handler.
   integer(I4P)                      :: buff_I4P        !< I4P buffer.

   call file_parameters%get(section_name='amr', option_name='max_level'    , val=buff_I4P) ; self%max_level     = buff_I4P
   call file_parameters%get(section_name='amr', option_name='iu_ref_levels', val=buff_I4P) ; self%iu_ref_levels = buff_I4P
   call file_parameters%get(section_name='amr', option_name='i_prune'      , val=buff_I4P) ; self%ijkl_prune(1) = buff_I4P
   call file_parameters%get(section_name='amr', option_name='j_prune'      , val=buff_I4P) ; self%ijkl_prune(2) = buff_I4P
   call file_parameters%get(section_name='amr', option_name='k_prune'      , val=buff_I4P) ; self%ijkl_prune(3) = buff_I4P
   call file_parameters%get(section_name='amr', option_name='l_prune'      , val=buff_I4P) ; self%ijkl_prune(4) = buff_I4P
   endsubroutine load_from_ini_file

   subroutine load_surface_stl(self, file_name)
   !< Load surface from STL file and compute signed distance of nodes from it.
   class(tree_object), intent(inout) :: self      !< The tree.
   character(*),       intent(in)    :: file_name !< STL file name.
   type(file_stl_object)             :: file_stl  !< STL file handler.

   print '(A)', self%mpih%myrankstr//'ADAM: load STL file: '//trim(adjustl(file_name))
   call file_stl%load_from_file(facet=self%surface_stl%facet, file_name=trim(adjustl(file_name)), guess_format=.true.)
   call self%surface_stl%analize(aabb_refinement_levels=3)
   call self%surface_stl%sanitize
   print '(A)', self%mpih%myrankstr//'ADAM: compute distance from STL surface'
   call self%compute_surface_stl_distance(surface_stl=self%surface_stl)
   endsubroutine load_surface_stl

   function loop(self, code, node_ptr) result(again)
   !< Sentinel while-loop on nodes returning the code (for tree looping).
   class(tree_object),     intent(in)                     :: self      !< The tree bucket.
   integer(I8P),           intent(out), optional          :: code      !< The Morton code.
   type(tree_node_object), intent(out), optional, pointer :: node_ptr  !< Pointer to current node.
   logical                                                :: again     !< Sentinel flag to contine the loop.
   integer(I4P), save                                     :: b=1_I4P   !< Bucket counter.
   type(tree_node_object), pointer, save                  :: p=>null() !< Pointer to current node.

   again = .false.
   if (self%nodes_number>0) then
      do
         if (b>self%buckets_number) then
            b = 1
            p => null()
            again = .false.
            return
         else
            if (self%bucket(b)%nodes_number>0) then
               if (.not.associated(p)) then
                  p => self%bucket(b)%head
                  if (present(code)) code = p%code
                  if (present(node_ptr)) node_ptr => p
                  again = .true.
                  return
               elseif (associated(p%next)) then
                  p => p%next
                  if (present(code)) code = p%code
                  if (present(node_ptr)) node_ptr => p
                  again = .true.
                  return
               else
                  b = b + 1
                  p => null()
               endif
            else
               b = b + 1
               p => null()
            endif
         endif
      enddo
   endif
   endfunction loop

   subroutine make_local_maps_bc(self)
   !< Make local maps of boundary conditions.
   class(tree_object), intent(inout) :: self                  !< The tree.
   type(tree_node_object), pointer   :: node_ptr              !< Pointer to current node.
   integer(I8P), allocatable         :: neighbor(:)           !< List of code neighbors.
   type(tree_node_object), pointer   :: neigh                 !< Pointer to node neighbor.
   integer(I4P)                      :: neighbor_type         !< Neighbors type.
   integer(I4P)                      :: neighbor_bc_fec       !< Neighbors fec for BC.
   integer(I4P)                      :: my_fec_number         !< Number of faces/edges/corners ghost cells locally exchanged.
   integer(I4P)                      :: fec, mf, rf, sf, n, p !< Counter.
   integer(I8P)                      :: ijk_min_max_delta(1:9) !< IJK min/max/delta.
   integer(I4P)                      :: fec_bc_faces_number
   integer(I4P)                      :: fec_bc_edges_number
   integer(I4P)                      :: fec_bc_corners_number
   integer(I4P)                      :: fec_bc_type

   if (allocated(self%local_map_bc_face))   deallocate(self%local_map_bc_face)
   if (allocated(self%local_map_bc_edge))   deallocate(self%local_map_bc_edge)
   if (allocated(self%local_map_bc_corner)) deallocate(self%local_map_bc_corner)
   fec_bc_faces_number = 0_I4P
   fec_bc_edges_number = 0_I4P
   fec_bc_corners_number = 0_I4P
   do while(self%loop(node_ptr=node_ptr))
      if (node_ptr%myrank==self%mpih%myrank) then
         do fec=1, 26
            ! call self%get_neighbor_all(code=node_ptr%code,          &
            !                            face=fec,                    &
            !                            neighbor=neighbor,           &
            !                            neighbor_type=neighbor_type, &
            !                            neighbor_bc_fec=neighbor_bc_fec)
            neighbor        = node_ptr%neighbor(fec)%codes
            neighbor_type   = node_ptr%neighbor(fec)%ntype
            neighbor_bc_fec = node_ptr%neighbor(fec)%bc_fec
            if (neighbor_type == NODE_BOUNDARY_CONDITION) then
               if (fec<=6)             fec_bc_faces_number   = fec_bc_faces_number   + 1_I4P
               if (fec>=7.and.fec<=18) fec_bc_edges_number   = fec_bc_edges_number   + 1_I4P
               if (fec>=19)            fec_bc_corners_number = fec_bc_corners_number + 1_I4P
            endif
         enddo
      endif
   enddo
   if (fec_bc_faces_number > 0.or.fec_bc_edges_number > 0.or.fec_bc_corners_number > 0) then
      if (fec_bc_faces_number > 0)   allocate(self%local_map_bc_face(fec_bc_faces_number, 12))
      if (fec_bc_edges_number > 0)   allocate(self%local_map_bc_edge(fec_bc_edges_number, 12))
      if (fec_bc_corners_number > 0) allocate(self%local_map_bc_corner(fec_bc_corners_number, 12))
      fec_bc_faces_number = 0_I4P
      fec_bc_edges_number = 0_I4P
      fec_bc_corners_number = 0_I4P
      do while(self%loop(node_ptr=node_ptr))
         if (node_ptr%myrank==self%mpih%myrank) then
            do fec=1, 26
               ! call self%get_neighbor_all(code=node_ptr%code,          &
               !                            face=fec,                    &
               !                            neighbor=neighbor,           &
               !                            neighbor_type=neighbor_type, &
               !                            neighbor_bc_fec=neighbor_bc_fec)
               neighbor        = node_ptr%neighbor(fec)%codes
               neighbor_type   = node_ptr%neighbor(fec)%ntype
               neighbor_bc_fec = node_ptr%neighbor(fec)%bc_fec
               if (neighbor_type == NODE_BOUNDARY_CONDITION) then
                  select case(neighbor_bc_fec)
                  case(1,7,9,11,13,19,21,23,25)
                     ! BC i min
                     fec_bc_type = self%grid%bc_type(1)
                  case(2,8,10,12,14,20,22,24,26)
                     ! BC i max
                     fec_bc_type = self%grid%bc_type(2)
                  case(3,15,17)
                     ! BC j min
                     fec_bc_type = self%grid%bc_type(3)
                  case(4,16,18)
                     ! BC j max
                     fec_bc_type = self%grid%bc_type(4)
                  case(5)
                     ! BC k min
                     fec_bc_type = self%grid%bc_type(5)
                  case(6)
                     ! BC k max
                     fec_bc_type = self%grid%bc_type(6)
                  endselect
                  call compute_ijk_min_max_delta(fec=fec, neighbor_bc_fec=neighbor_bc_fec, ijk_min_max_delta=ijk_min_max_delta)
                  if (fec<=6) then
                     fec_bc_faces_number = fec_bc_faces_number + 1_I4P
                     self%local_map_bc_face(fec_bc_faces_number,1)    = node_ptr%block_index
                     self%local_map_bc_face(fec_bc_faces_number,2)    = neighbor_bc_fec
                     self%local_map_bc_face(fec_bc_faces_number,3:11) = ijk_min_max_delta
                     self%local_map_bc_face(fec_bc_faces_number,12)   = fec_bc_type
                  endif
                  if (fec>=7.and.fec<=18) then
                     fec_bc_edges_number = fec_bc_edges_number + 1_I4P
                     self%local_map_bc_edge(fec_bc_edges_number,1)    = node_ptr%block_index
                     self%local_map_bc_edge(fec_bc_edges_number,2)    = neighbor_bc_fec
                     self%local_map_bc_edge(fec_bc_edges_number,3:11) = ijk_min_max_delta
                     self%local_map_bc_edge(fec_bc_edges_number,12)   = fec_bc_type
                  endif
                  if (fec>=19) then
                     fec_bc_corners_number = fec_bc_corners_number + 1_I4P
                     self%local_map_bc_corner(fec_bc_corners_number,1)    = node_ptr%block_index
                     self%local_map_bc_corner(fec_bc_corners_number,2)    = neighbor_bc_fec
                     self%local_map_bc_corner(fec_bc_corners_number,3:11) = ijk_min_max_delta
                     self%local_map_bc_corner(fec_bc_corners_number,12)   = fec_bc_type
                  endif
               endif
            enddo
         endif
      enddo
   endif
   contains
      subroutine compute_ijk_min_max_delta(fec, neighbor_bc_fec, ijk_min_max_delta)
      !< Compute IJK min/max/delta.
      integer(I4P), intent(in)  :: fec                    !< Current fec number.
      integer(I4P)              :: neighbor_bc_fec        !< Neighbors fec for BC.
      integer(I8P), intent(out) :: ijk_min_max_delta(1:9) !< IJK min/max/delta.
      integer(I4P)              :: ijkmin(3)              !< Lower limit of ijk indexes.
      integer(I4P)              :: ijkmax(3)              !< Upper limit of ijk indexes.
      integer(I4P)              :: ijkdelta(3)            !< Delta offset for ghost-inner cells mapping same refinement.
      integer(I4P)              :: ijkdelta_global(3)     !< Delta offset for ghost-inner cells mapping same refinement.
      integer(I4P)              :: nijk(3)                !< Ni, Nj , Nk stored in array.
      integer(I4P)              :: i                      !< Counter.

      associate(ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk, ngc=>self%grid%ngc)
         nijk = [ni, nj, nk]

         ijkdelta = FEC_TO_DELTA(1:3, fec)
         ijkdelta_global = FEC_TO_DELTA(1:3, neighbor_bc_fec)

         do i=1, 3
            if     (ijkdelta(i)==1) then
               ijkmin(i) = nijk(i) + 1
               ijkmax(i) = nijk(i) + ngc
            elseif (ijkdelta(i)==-1) then
               ijkmin(i) = 0
               ijkmax(i) = 1 - ngc
            elseif (ijkdelta(i)==0) then
               ijkmin(i) = 1
               ijkmax(i) = nijk(i)
            endif
         enddo
         ijk_min_max_delta = [ijkmin, ijkmax, ijkdelta_global]
      endassociate
      endsubroutine
   endsubroutine make_local_maps_bc

   subroutine make_neighborhood(self)
   !< Make neighborhood all whole tree and store it in nodes.
   class(tree_object), intent(inout) :: self     !< The tree.
   type(tree_node_object), pointer   :: node_ptr !< Pointer to current node.
   integer(I4P)                      :: fec      !< Counter.

   do while(self%loop(node_ptr=node_ptr))
      if (node_ptr%i_am_new) then
         do fec=1, 26
            call self%get_neighbor_all(code             = node_ptr%code,                  &
                                       face             = fec,                            &
                                       neighbor         = node_ptr%neighbor(fec)%codes,   &
                                       neighbor_type    = node_ptr%neighbor(fec)%ntype,   &
                                       neighbor_portion = node_ptr%neighbor(fec)%portion, &
                                       neighbor_bc_fec  = node_ptr%neighbor(fec)%bc_fec)
         enddo
         node_ptr%i_am_new = .false.
      else
         do fec=1, 26
            if (allocated(node_ptr%neighbor(fec)%codes)) then
               if (.not.self%has_code(code=node_ptr%neighbor(fec)%codes(1))) then
                  call self%get_neighbor_all(code             = node_ptr%code,                  &
                                             face             = fec,                            &
                                             neighbor         = node_ptr%neighbor(fec)%codes,   &
                                             neighbor_type    = node_ptr%neighbor(fec)%ntype,   &
                                             neighbor_portion = node_ptr%neighbor(fec)%portion, &
                                             neighbor_bc_fec  = node_ptr%neighbor(fec)%bc_fec)
               endif
            endif
         enddo
      endif
   enddo
   endsubroutine make_neighborhood

   function max_cell_delta(self, distance) result(delta)
   !< Return the maximum cell delta given a comparison distance.
   class(tree_object), intent(in) :: self     !< The field.
   real(R8P),          intent(in) :: distance !< Comparison distance.
   real(R8P)                      :: delta    !< Maximum cell delta admissible.

   if (abs(distance) < epsilon(0._R8P)) then
      ! delta = 0.001_R8P
      delta = 0.005_R8P
   else
      delta = huge(0._R8P)
   endif
   endfunction max_cell_delta

   subroutine mark_all_nodes(self, mark)
   !< Mark all nodes to be refined, derefined, ecc.
   class(tree_object), intent(inout) :: self     !< The tree.
   integer(I4P),       intent(in)    :: mark     !< Mark to be imposed [TO_BE_REFINED,...].
   type(tree_node_object), pointer   :: node_ptr !< Pointer to current node.

   do while(self%loop(node_ptr=node_ptr))
      node_ptr%refinement_needed = mark
   enddo
   endsubroutine mark_all_nodes

   subroutine mark_sphere(self, center, radius, threshold)
   !< Mark all nodes inside a sphere to be refined.
   class(tree_object), intent(inout)        :: self             !< The tree.
   real(R8P),          intent(in)           :: center(3)        !< Sphere center coordinates [x,y,z].
   real(R8P),          intent(in)           :: radius           !< Sphere radius.
   real(R8P),          intent(in), optional :: threshold        !< Threshold for sphere proximity.
   real(R8P)                                :: threshold_       !< Threshold for sphere proximity, local var.
   type(tree_node_object), pointer          :: node_ptr         !< Pointer to current node.
   real(R8P)                                :: block_center(3)  !< block center coordinates.
   real(R8P)                                :: block_diagonal   !< block diagonal.
   real(R8P)                                :: distance(0:8)    !< Distances between block and sphere.
   real(R8P)                                :: max_delta        !< Max cell delta.
   integer(I4P)                             :: i,j,k,l          !< Counter.
   real(R8P)                                :: emin(3), emax(3) !< Node extents.

   threshold_ = 2.2_R8P ; if (present(threshold)) threshold_ = threshold
   do while(self%loop(node_ptr=node_ptr))
      call self%morton_to_coordinates(code=node_ptr%code, i=i, j=j, k=k, l=l)
      call self%grid%compute_metrics(coordinates=[i,j,k,l], emin=emin, emax=emax)
      block_center = (emax + emin) / 2._R8P
      block_diagonal = sqrt((emax(1) - emin(1))**2 + &
                            (emax(2) - emin(2))**2 + &
                            (emax(3) - emin(3))**2)

      associate (ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk)
      distance(0) = sphere_distance(point=block_center)
      distance(1) = sphere_distance(point=[emin(1), emin(2), emin(3)])
      distance(2) = sphere_distance(point=[emax(1), emin(2), emin(3)])
      distance(3) = sphere_distance(point=[emin(1), emax(2), emin(3)])
      distance(4) = sphere_distance(point=[emax(1), emax(2), emin(3)])
      distance(5) = sphere_distance(point=[emin(1), emin(2), emax(3)])
      distance(6) = sphere_distance(point=[emax(1), emin(2), emax(3)])
      distance(7) = sphere_distance(point=[emin(1), emax(2), emax(3)])
      distance(8) = sphere_distance(point=[emax(1), emax(2), emax(3)])
      if (maxval(distance(0:8),dim=1)*minval(distance(0:8),dim=1) < 0._R8P) then
         distance(0) = 0._R8P
      endif

      max_delta = self%max_cell_delta(distance=distance(0))

      if (block_diagonal/min(ni,nj,nk) > max_delta) then
         node_ptr%refinement_needed = TO_BE_REFINED
      elseif (block_diagonal/min(ni,nj,nk) * threshold_ < max_delta) then
         node_ptr%refinement_needed = TO_BE_DEREFINED
      endif
      endassociate
   enddo
   contains
      pure function sphere_distance(point)
      !< Return the distance from a point to the sphere surface, with sign.
      real(R8P), intent(in) :: point(3)        !< Point coordinates.
      real(R8P)             :: sphere_distance !< Distance from sphere surface.

      sphere_distance = sqrt((center(1) - point(1))**2 + &
                             (center(2) - point(2))**2 + &
                             (center(3) - point(3))**2) - radius
      endfunction sphere_distance
   endsubroutine mark_sphere

   subroutine mark_surface_stl(self, surface_stl, threshold)
   !< Mark all nodes inside a surface defined by STL triangulation.
   class(tree_object),       intent(inout)        :: self             !< The tree.
   type(surface_stl_object), intent(in)           :: surface_stl      !< STL surface.
   real(R8P),                intent(in), optional :: threshold        !< Threshold for sphere proximity.
   real(R8P)                                      :: threshold_       !< Threshold for sphere proximity, local var.
   type(tree_node_object), pointer                :: node_ptr         !< Pointer to current node.
   real(R8P)                                      :: block_diagonal   !< block diagonal.
   real(R8P)                                      :: max_delta        !< Max cell delta.
   integer(I4P)                                   :: i,j,k,l          !< Counter.
   real(R8P)                                      :: emin(3), emax(3) !< Node extents.

   threshold_ = 2.2_R8P ; if (present(threshold)) threshold_ = threshold
   do while(self%loop(node_ptr=node_ptr))
      if (node_ptr%myrank==self%mpih%myrank) then
         call self%morton_to_coordinates(code=node_ptr%code, i=i, j=j, k=k, l=l)
         call self%grid%compute_metrics(coordinates=[i,j,k,l], emin=emin, emax=emax)
         block_diagonal = sqrt((emax(1) - emin(1))**2 + &
                               (emax(2) - emin(2))**2 + &
                               (emax(3) - emin(3))**2)
         associate (ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk)
            max_delta = self%max_cell_delta(distance=node_ptr%surface_stl_distance)
            if (block_diagonal/min(ni,nj,nk) > max_delta) then
               node_ptr%refinement_needed = TO_BE_REFINED
            elseif (block_diagonal/min(ni,nj,nk) * threshold_ < max_delta) then
               node_ptr%refinement_needed = TO_BE_DEREFINED
            endif
         endassociate
      endif
   enddo
   endsubroutine mark_surface_stl

   function node(self, code) result(p)
   !< Return a pointer to a node in the tree.
   class(tree_object), intent(in)  :: self !< The tree.
   integer(I8P),       intent(in)  :: code !< The Morton code.
   type(tree_node_object), pointer :: p    !< Pointer to node queried.

   p => null()
   if (self%is_initialized_) p => self%bucket(self%hash(code=code))%node(code=code)
   endfunction node

   elemental function prime_buckets_number(self, nodes_number) result(buckets_number)
   !< Return the buckets number as the nearest prime number given nodes number.
   !<
   !< @note The balanced buckets number is computing considering the tree load defined in `self` and using the
   !< Sieve of Eratoshenes for findining the nearest prime number.
   class(tree_object), intent(in) :: self           !< The tree.
   integer(I8P),       intent(in) :: nodes_number   !< Nodes number to be stored in the tree.
   integer(I8P)                   :: buckets_number !< Well balanced, prime buckets number.
   logical, allocatable           :: is_prime(:)    !< List of prime numbers up to buckets number.
   integer(I8P)                   :: b              !< Counter.

   buckets_number = int((1._R8P / self%max_load) * nodes_number)
   allocate(is_prime(buckets_number))
   is_prime = .true.
   is_prime(1) = .false.
   do b=2, int(sqrt(real(buckets_number, R8P)))
      if (is_prime(b)) is_prime(b*b:buckets_number:b) = .false.
   enddo
   b = buckets_number
   do while(.not.is_prime(b))
      b = b - 1
   enddo
   buckets_number = b
   endfunction prime_buckets_number

   subroutine print_status(self)
   !< Print status of main data.
   class(tree_object), intent(in) :: self !< The tree.

   print '(A)', self%mpih%myrankstr//'tree status of main data'
   print '(A)', self%mpih%myrankstr//'  buckets number: '//trim(str(self%buckets_number))
   print '(A)', self%mpih%myrankstr//'  nodes number:   '//trim(str(self%nodes_number  ))
   print '(A)', self%mpih%myrankstr//'  max load:       '//trim(str(self%max_load      ))
   print '(A)', self%mpih%myrankstr//'  ratio:          '//trim(str(self%ratio         ))
   print '(A)', self%mpih%myrankstr//'  max level:      '//trim(str(self%max_level     ))
   print '(A)', self%mpih%myrankstr//'  ijkl_prune:     '//trim(str(self%ijkl_prune    ))
   print '(A)', self%mpih%myrankstr//'  iu_ref_levels:  '//trim(str(self%iu_ref_levels ))
   print '(A)', self%mpih%myrankstr//''
   endsubroutine print_status

   subroutine prune(self, ijkl_prune)
   !< Prune nodes.
   class(tree_object), intent(inout) :: self          !< The tree.
   integer(I4P),       intent(inout) :: ijkl_prune(4) !< Maximum coordinates after which the prune operates.
   type(tree_node_object), pointer   :: node_ptr      !< Pointer to current node.
   integer(I4P)                      :: ijkl(4)       !< Coordinates counter.
   integer(I4P)                      :: i             !< Counter.

   if (ijkl_prune(4)>0) then
      do i=1, 3
         ijkl_prune(i) = min(ijkl_prune(i), 2**ijkl_prune(4) - 1)
      enddo
      self%ijkl_prune = ijkl_prune
      print '(A)', self%mpih%myrankstr//'prune tree with IJKL max: '//trim(str(self%ijkl_prune))
      do while(self%loop(node_ptr=node_ptr))
         call self%morton_to_coordinates(code=node_ptr%code, i=ijkl(1), j=ijkl(2), k=ijkl(3), l=ijkl(4))
         if (ijkl(4)/=self%ijkl_prune(4)) then
            print '(A)', self%mpih%myrankstr//'ERROR: cannot prune nodes at different prune-level, node: '//&
                         trim(str(node_ptr%code))
         endif
         if (ijkl(1)>self%ijkl_prune(1).or.ijkl(2)>self%ijkl_prune(2).or.ijkl(3)>self%ijkl_prune(3)) then
            call self%remove_node(code=node_ptr%code)
         endif
      enddo
   endif
   endsubroutine prune

   subroutine remove_node(self, code)
   !< Remove a node from the tree, given the code.
   class(tree_object), intent(inout) :: self !< The tree.
   integer(I8P),       intent(in)    :: code !< The Morton code.
   integer(I4P)                      :: b    !< Bucket index, namely hashed key.

   if (self%is_initialized_) then
      if (self%has_code(code=code)) then
         b = self%hash(code=code)
         call self%bucket(b)%remove_node(code=code)
         self%nodes_number = self%nodes_number - 1
      endif
   endif
   endsubroutine remove_node

   subroutine resize(self, nodes_number, max_load)
   !< Resize the tree.
   class(tree_object), intent(inout)        :: self         !< The tree.
   integer(I8P),       intent(in)           :: nodes_number !< Nodes number to be stored in the tree.
   real(R8P),          intent(in), optional :: max_load     !< Maximum load of tree buckets.
   type(tree_object)                        :: swap         !< Temporary (swap) tree.
   type(tree_node_object), pointer          :: node_ptr     !< Pointer to node.

   if (self%is_initialized_) then
      if (present(max_load)) self%max_load = max_load
      if (self%nodes_number > int((1._R8P/self%max_load)*nodes_number, I4P)) return ! new size too small, cannot previous nodes
      call swap%initialize(grid=self%grid, max_load=self%max_load, nodes_number=nodes_number, ratio=self%ratio, add_adam=.false.)
      do while(self%loop(node_ptr=node_ptr)) ! re-hash all codes
         call swap%add_node(code=node_ptr%code,                           &
                            refinement_needed=node_ptr%refinement_needed, &
                            rank=node_ptr%myrank,                         &
                            block_index=node_ptr%block_index)
      enddo
      !RIMETTERE SENZA
      if(allocated(self%bucket)) deallocate(self%bucket)
      !RIMETTERE SENZA
      call move_alloc(from=swap%bucket, to=self%bucket)
      self%buckets_number = swap%buckets_number
      self%nodes_number   = swap%nodes_number
   else
      print '(A)', self%mpih%myrankstr//'ERROR: tree is not initialized, cannot be resized'
      stop
   endif
   endsubroutine resize

   subroutine save_local_map(self, file_name)
   !< Save local map.
   class(tree_object), intent(in)  :: self      !< The tree.
   character(*),       intent(in)  :: file_name !< Output file name.
   integer(I4P)                    :: file_unit !< Output file unit.
   integer(I4P)                    :: i1,i2     !< Counter.

   if (self%is_initialized_) then
      if (allocated(self%local_map_ghost)) then
         open(newunit=file_unit, file=trim(adjustl(file_name)), form='FORMATTED')
         do i1=lbound(self%local_map_ghost, dim=1), ubound(self%local_map_ghost, dim=1)
            ! do i2=lbound(self%local_map_ghost, dim=2), ubound(self%local_map_ghost, dim=2)
               write(unit=file_unit, fmt='(13(I6,1X))') self%local_map_ghost(i1, :)
            ! enddo
         enddo
         close(file_unit)
      endif
   endif
   endsubroutine save_local_map

   subroutine save_nodes(self, file_name)
   !< Save nodes data, used for restart.
   class(tree_object), intent(in)  :: self            !< The tree.
   character(*),       intent(in)  :: file_name       !< Output file name.
   integer(I8P), allocatable       :: actual_codes(:) !< List of actually stored codes.
   integer(I8P)                    :: codes_number    !< Number of codes actually stored.
   integer(I4P)                    :: file_unit       !< Output file unit.
   type(tree_node_object), pointer :: node_ptr        !< Pointer to current node.
   integer(I8P)                    :: c               !< Counter.

   if (self%is_initialized_) then
      actual_codes = self%codes()
      codes_number = size(actual_codes, dim=1)
      if (codes_number > 0) then
         print '(A)', self%mpih%myrankstr//'save tree nodes in file  '//trim(adjustl(file_name))
         open(newunit=file_unit, file=trim(adjustl(file_name)), form='UNFORMATTED', access='STREAM')
         write(unit=file_unit) codes_number
         do c=1, codes_number
            node_ptr => self%node(code=actual_codes(c))
            write(unit=file_unit) node_ptr%code, node_ptr%myrank, node_ptr%block_index
         enddo
         close(file_unit)
      endif
   endif
   endsubroutine save_nodes

   subroutine traverse(self, iterator)
   !< Traverse tree calling the iterator procedure.
   class(tree_object), intent(in) :: self     !< The tree.
   procedure(iterator_interface)  :: iterator !< The (key) iterator procedure to call for each node.
   integer(I8P)                   :: b        !< Counter.

   if (self%is_initialized_) then
      do b=1_I8P, self%buckets_number
         call self%bucket(b)%traverse(iterator)
      enddo
   endif
   endsubroutine traverse

   ! MPI methods
   subroutine blocks_reorder(self)
   !< Reorder blocks indexes in field.
   class(tree_object), intent(inout) :: self                !< The tree.
   type(tree_node_object), pointer   :: node_ptr            !< Pointer to current node.
   integer(I4P)                      :: outer_blocks_number !< Number of outer blocks where I need fecs.
   integer(I4P)                      :: inner_blocks_number !< Number of inner blocks where I need fecs.
   logical                           :: is_inner_block      !< Flag to check if a neighbor block is inner or not.
   integer(I8P), allocatable         :: neighbor(:)         !< List of code neighbors.
   type(tree_node_object), pointer   :: neigh               !< Pointer to node neighbor.
   integer(I4P)                      :: neighbor_type       !< Neighbors type.
   integer(I4P), allocatable         :: inner_block_map(:)  !< Inner blocks map.
   integer(I4P), allocatable         :: outer_block_map(:)  !< Outer blocks map.
   integer(I4P)                      :: fec, n              !< Counter.

   outer_blocks_number = 0
   inner_blocks_number = 0
   allocate(inner_block_map(self%my_nodes_number))
   allocate(outer_block_map(self%my_nodes_number))
   if (allocated(self%inner_outer_block_map)) deallocate(self%inner_outer_block_map)
   allocate(self%inner_outer_block_map(self%my_nodes_number))
   do while(self%loop(node_ptr=node_ptr))
      if (self%mpih%myrank == node_ptr%myrank) then
         is_inner_block = .true.
         do fec=1, 26
            ! call self%get_neighbor_all(code=node_ptr%code, face=fec, neighbor=neighbor, neighbor_type=neighbor_type)
            neighbor      = node_ptr%neighbor(fec)%codes
            neighbor_type = node_ptr%neighbor(fec)%ntype
            if (neighbor_type /= NODE_BOUNDARY_CONDITION) then
               do n=1, size(neighbor, dim=1)
                  neigh => self%node(code=neighbor(n))
                  if (self%mpih%myrank /= neigh%myrank) is_inner_block = .false.
               enddo
            endif
         enddo
         if (is_inner_block) then
            inner_blocks_number = inner_blocks_number + 1
            inner_block_map(inner_blocks_number) = node_ptr%block_index
            node_ptr%block_index_new = inner_blocks_number
         else
            outer_blocks_number = outer_blocks_number + 1
            outer_block_map(outer_blocks_number) = node_ptr%block_index
            node_ptr%block_index_new = -outer_blocks_number
         endif
      endif
   enddo
   self%inner_blocks_number = inner_blocks_number
   self%inner_outer_block_map(1:inner_blocks_number ) = inner_block_map(1:inner_blocks_number)
   self%inner_outer_block_map(inner_blocks_number+1:) = outer_block_map(1:outer_blocks_number)
   do while(self%loop(node_ptr=node_ptr))
      if (self%mpih%myrank == node_ptr%myrank) then
         if (node_ptr%block_index_new < 0) then
            node_ptr%block_index = -node_ptr%block_index_new + inner_blocks_number
         else
            node_ptr%block_index =  node_ptr%block_index_new
         endif
      endif
   enddo
   call self%mpi_gather_nodes_data(node_member='block_index')
   endsubroutine blocks_reorder

   subroutine import_refinements_needed(self, refinements_needed_all, disp_count)
   !< Import refinements needed status changed externally.
   class(tree_object),        intent(inout) :: self                      !< The tree.
   integer(I4P), allocatable, intent(in)    :: refinements_needed_all(:) !< Refinements needed of all blocks.
   integer(I4P), allocatable, intent(in)    :: disp_count(:)             !< Displacement of blocks that are received from process.
   type(tree_node_object), pointer          :: node_ptr                  !< Pointer to current node.
   integer(I8P)                             :: b                         !< Counter.
   integer(I4P)                             :: myrank                    !< Counter.

   do while(self%loop(node_ptr=node_ptr))
      myrank = node_ptr%myrank
      b = node_ptr%block_index
      node_ptr%refinement_needed = refinements_needed_all(disp_count(myrank)+b)
   enddo
   endsubroutine import_refinements_needed

   subroutine make_comm_local_maps(self)
   !< Make communication/local maps.
   !<```comm_map_send     = [ 17, 511, 92, 3, 54, 56, 11, 12...] (block index).
   !<                          |   |       |  ||       |
   !<   comm_map_send_ptr = [  0,  1,  3,  4,  4, 6, (8)]        (pointer to comm_map_send)
   !<   comm_map_recv     = [ 23, 4, 51, 69, 145, 2, 72, 16, 6]  (block index).
   !<                         |       |  ||          |       |
   !<   comm_map_recv_prt = [ 0,  2,  3,  3,  6, 8, (9)]         (pointer to comm_map_recv)```
   class(tree_object), intent(inout) :: self                 !< The tree.
   type(tree_node_object), pointer   :: node_ptr             !< Pointer to current node.
   integer(I8P), allocatable         :: codes_sorted(:)      !< List of (sorted) codes.
   integer(I4P), allocatable         :: comm_map_send_ctr(:) !< Communication map, counters in list to send [procs_number+1].
   integer(I4P), allocatable         :: comm_map_recv_ctr(:) !< Communication map, counters in list to recv [procs_number+1].
   integer(I8P)                      :: my_nodes_number      !< Number of my nodes.
   integer(I4P)                      :: n_send               !< Number of nodes that I have to send.
   integer(I4P)                      :: n_recv               !< Number of nodes that I have to receive.
   integer(I4P)                      :: n_keep               !< Number of nodes that I have to keep.
   integer(I8P)                      :: c, l                 !< Counter.
   integer(I4P)                      :: p                    !< Counter.

   ! initialize communication maps
   self%comm_map_n_send = 0_I4P
   self%comm_map_n_recv = 0_I4P
   if (allocated(self%comm_map_send)) deallocate(self%comm_map_send)
   if (allocated(self%comm_map_recv)) deallocate(self%comm_map_recv)
   if (allocated(self%local_map    )) deallocate(self%local_map    )

   ! compute the number of blocks to send/receive
   my_nodes_number = 0_I8P
   do while(self%loop(node_ptr=node_ptr))
      if (node_ptr%myrank==self%mpih%myrank) my_nodes_number = my_nodes_number + 1_I8P
      if     (is_node_to_send(n=node_ptr)) then
         ! I have this node that must be sent to node_ptr%myrank_new
         self%comm_map_n_send(node_ptr%myrank_new) = self%comm_map_n_send(node_ptr%myrank_new) + 1
      elseif (is_node_to_receive(n=node_ptr)) then
         ! node_ptr%rank has this node that must be sent to me
         self%comm_map_n_recv(node_ptr%myrank) = self%comm_map_n_recv(node_ptr%myrank) + 1
      endif
   enddo
   n_send = sum(self%comm_map_n_send, dim=1)
   n_recv = sum(self%comm_map_n_recv, dim=1)
   n_keep = my_nodes_number - n_send
   self%my_nodes_number = n_keep + n_recv
   if (n_keep > 0_I4P) allocate(self%local_map(n_keep,2))

   ! allocate communications maps
   if (n_send>0_I4P) then
      allocate(self%comm_map_send(n_send))
      self%comm_map_send = 0_I8P
   endif
   if (n_recv>0_I4P) then
      allocate(self%comm_map_recv(n_recv))
      self%comm_map_recv = 0_I8P
   endif

   ! compute communication maps pointers/counters
   self%comm_map_send_ptr = 0_I4P
   self%comm_map_recv_ptr = 0_I4P
   do p=1, self%mpih%procs_number
      self%comm_map_send_ptr(p) = self%comm_map_send_ptr(p-1) + self%comm_map_n_send(p-1)
      self%comm_map_recv_ptr(p) = self%comm_map_recv_ptr(p-1) + self%comm_map_n_recv(p-1)
   enddo
   comm_map_send_ctr = self%comm_map_send_ptr
   comm_map_recv_ctr = self%comm_map_recv_ptr

   ! populate communication maps
   codes_sorted = self%codes() ! sorted list of Morton codes
   ! send map
   if (n_send>0_I4P) then
      do c=1, size(codes_sorted, dim=1) ! sort blocks in Morton order
         node_ptr => self%node(code=codes_sorted(c))
         if (is_node_to_send(n=node_ptr)) then
            self%comm_map_send(comm_map_send_ctr(node_ptr%myrank_new)+1) = node_ptr%block_index
            comm_map_send_ctr(node_ptr%myrank_new) = comm_map_send_ctr(node_ptr%myrank_new) + 1
         endif
      enddo
   endif
   ! receive map
   if (n_recv>0_I4P) then
      do c=1, size(codes_sorted, dim=1) ! sort blocks in Morton order
         node_ptr => self%node(code=codes_sorted(c))
         if (is_node_to_receive(n=node_ptr)) then
            self%comm_map_recv(comm_map_recv_ctr(node_ptr%myrank)+1) = node_ptr%block_index_new
            comm_map_recv_ctr(node_ptr%myrank) = comm_map_recv_ctr(node_ptr%myrank) + 1
         endif
      enddo
   endif
   ! local maps
   if (n_keep > 0_I4P) then
      l = 0
      do c=1, size(codes_sorted, dim=1) ! sort blocks in Morton order
         node_ptr => self%node(code=codes_sorted(c))
         if (is_node_to_keep(n=node_ptr)) then
            l = l + 1
            self%local_map(l,1) = node_ptr%block_index_new
            self%local_map(l,2) = node_ptr%block_index
         endif
      enddo
   endif
   contains
      function is_node_to_keep(n)
      !< Check if node `n` must be kept.
      type(tree_node_object), intent(in), pointer :: n               !< Pointer to current node.
      logical                                     :: is_node_to_keep !< Check result.

      is_node_to_keep = ((self%mpih%myrank == n%myrank).and.(n%myrank == n%myrank_new))
      endfunction is_node_to_keep

      function is_node_to_send(n)
      !< Check if node `n` must be sent.
      type(tree_node_object), intent(in), pointer :: n               !< Pointer to current node.
      logical                                     :: is_node_to_send !< Check result.

      is_node_to_send = ((self%mpih%myrank == n%myrank).and.(n%myrank /= n%myrank_new))
      endfunction is_node_to_send

      function is_node_to_receive(n)
      !< Check if node `n` must be received.
      type(tree_node_object), intent(in), pointer :: n                  !< Pointer to current node.
      logical                                     :: is_node_to_receive !< Check result.

      is_node_to_receive = ((self%mpih%myrank == n%myrank_new).and.(n%myrank /= n%myrank_new))
      endfunction is_node_to_receive
   endsubroutine make_comm_local_maps

   subroutine make_comm_local_maps_ghost(self)
   !< Make communication/local maps of ghost cells.
   class(tree_object), intent(inout) :: self                       !< The tree.
   type(tree_node_object), pointer   :: node_ptr                   !< Pointer to current node.
   integer(I8P), allocatable         :: neighbor(:)                !< List of code neighbors.
   type(tree_node_object), pointer   :: neigh                      !< Pointer to node neighbor.
   integer(I4P)                      :: neighbor_type              !< Neighbors type.
   integer(I4P)                      :: neighbor_portion           !< Neighbors portion.
   integer(I4P)                      :: my_fec_number              !< Number of faces/edges/corners ghost cells locally exchanged.
   integer(I4P)                      :: recv_fec_number            !< Number of faces/edges/corners ghost cells externally received.
   integer(I4P)                      :: send_fec_number            !< Number of faces/edges/corners ghost cells externally sent.
   integer(I4P)                      :: fec, mf, rf, sf, n, p      !< Counter.
   integer(I4P)                      :: weight_reduction           !< Neighbor weight reduction for send/recv at different level.
   integer(I4P), allocatable         :: comm_map_send_ctr_ghost(:) !< Communication map, counters to send.
   integer(I4P), allocatable         :: comm_map_recv_ctr_ghost(:) !< Communication map, counters to recv.

   if (allocated(self%local_map_ghost    )) deallocate(self%local_map_ghost    )
   if (allocated(self%comm_map_send_ghost)) deallocate(self%comm_map_send_ghost)
   if (allocated(self%comm_map_recv_ghost)) deallocate(self%comm_map_recv_ghost)
   my_fec_number = 0
   recv_fec_number = 0
   send_fec_number = 0
   self%comm_map_n_recv_ghost = 0
   self%comm_map_n_send_ghost = 0
   do while(self%loop(node_ptr=node_ptr))
      do fec=1, 26
         weight_reduction = 2 ** count(FEC_TO_DELTA(:, fec)==0_I4P, dim=1)
         ! call self%get_neighbor_all(code=node_ptr%code, face=fec, neighbor=neighbor, neighbor_type=neighbor_type)
         neighbor      = node_ptr%neighbor(fec)%codes
         neighbor_type = node_ptr%neighbor(fec)%ntype
         if (neighbor_type /= NODE_BOUNDARY_CONDITION) then
            do n=1, size(neighbor, dim=1)
               neigh => self%node(code=neighbor(n))
               if     ((self%mpih%myrank == neigh%myrank).and.(self%mpih%myrank == node_ptr%myrank)) then
                  my_fec_number = my_fec_number + 1
               elseif ((self%mpih%myrank /= neigh%myrank).and.(self%mpih%myrank == node_ptr%myrank)) then
                  ! when receiving from same or less refined than me the size of the message is full, when receiving from
                  ! more refined the message is an averaged portion (reduced size)
                  recv_fec_number = recv_fec_number + 1
                  if (neighbor_type==NODE_MORE_REFINED) then
                     self%comm_map_n_recv_ghost(neigh%myrank) = self%comm_map_n_recv_ghost(neigh%myrank) + &
                                                                self%grid%weight_neighbor(fec)/ weight_reduction
                  else
                     self%comm_map_n_recv_ghost(neigh%myrank) = self%comm_map_n_recv_ghost(neigh%myrank) + &
                                                                self%grid%weight_neighbor(fec)
                  endif
               elseif ((self%mpih%myrank == neigh%myrank).and.(self%mpih%myrank /= node_ptr%myrank)) then
                  ! when sending to same or more refined than me the size of the message is full, when sending to less
                  ! refined the message is an averaged portion (reduced size)
                  send_fec_number = send_fec_number + 1
                  if (neighbor_type==NODE_MORE_REFINED) then
                     self%comm_map_n_send_ghost(node_ptr%myrank) = self%comm_map_n_send_ghost(node_ptr%myrank) + &
                                                                   self%grid%weight_neighbor(fec)/ weight_reduction
                  else
                     self%comm_map_n_send_ghost(node_ptr%myrank) = self%comm_map_n_send_ghost(node_ptr%myrank) + &
                                                                   self%grid%weight_neighbor(fec)
                  endif
               endif
            enddo
         endif
      enddo
   enddo
   ! create pointer for each process in the send/recv buffers
   self%comm_map_send_ptr_ghost = 0_I4P
   self%comm_map_recv_ptr_ghost = 0_I4P
   do p=1, self%mpih%procs_number
      self%comm_map_send_ptr_ghost(p) = self%comm_map_send_ptr_ghost(p-1) + self%comm_map_n_send_ghost(p-1)
      self%comm_map_recv_ptr_ghost(p) = self%comm_map_recv_ptr_ghost(p-1) + self%comm_map_n_recv_ghost(p-1)
   enddo
   ! create counter
   comm_map_send_ctr_ghost = self%comm_map_send_ptr_ghost
   comm_map_recv_ctr_ghost = self%comm_map_recv_ptr_ghost
   ! populate maps
   if (my_fec_number>0  ) allocate(self%local_map_ghost(    1:my_fec_number,  1:13))
   if (send_fec_number>0) allocate(self%comm_map_send_ghost(1:send_fec_number,1:15))
   if (recv_fec_number>0) allocate(self%comm_map_recv_ghost(1:recv_fec_number,1:15))
   sf = 0
   rf = 0
   mf = 0
   do while(self%loop(node_ptr=node_ptr))
      do fec=1, 26
         weight_reduction = 2 ** count(FEC_TO_DELTA(:, fec)==0_I4P, dim=1)
         ! call self%get_neighbor_all(code=node_ptr%code, face=fec, neighbor=neighbor, &
                                    ! neighbor_type=neighbor_type, neighbor_portion=neighbor_portion)
         neighbor         = node_ptr%neighbor(fec)%codes
         neighbor_type    = node_ptr%neighbor(fec)%ntype
         neighbor_portion = node_ptr%neighbor(fec)%portion
         if (neighbor_type /= NODE_BOUNDARY_CONDITION) then
            do n=1, size(neighbor, dim=1)
               neigh => self%node(code=neighbor(n))
               if     ((self%mpih%myrank == neigh%myrank).and.(self%mpih%myrank == node_ptr%myrank)) then
                  mf = mf + 1
                  self%local_map_ghost(mf, 1) = node_ptr%block_index
                  self%local_map_ghost(mf, 2) = neigh%block_index
                  self%local_map_ghost(mf, 3) = fec
                  if     (neighbor_type==NODE_STANDARD) then
                     self%local_map_ghost(mf, 4) = 0
                  elseif (neighbor_type==NODE_MORE_REFINED) then
                     self%local_map_ghost(mf, 4) = n
                  elseif (neighbor_type==NODE_LESS_REFINED) then
                     self%local_map_ghost(mf, 4) = -neighbor_portion
                  endif
                  call compute_ijk_min_max_delta(fec=fec, portion=self%local_map_ghost(mf, 4), &
                                                 ijk_min_max_delta=self%local_map_ghost(mf, 5:))
               elseif ((self%mpih%myrank /= neigh%myrank).and.(self%mpih%myrank == node_ptr%myrank)) then
                  rf = rf + 1
                  self%comm_map_recv_ghost(rf, 15) = comm_map_recv_ctr_ghost(neigh%myrank)
                  self%comm_map_recv_ghost(rf, 1) = node_ptr%block_index
                  self%comm_map_recv_ghost(rf, 2) = neigh%block_index
                  self%comm_map_recv_ghost(rf, 3) = neigh%myrank
                  self%comm_map_recv_ghost(rf, 4) = fec
                  if     (neighbor_type==NODE_STANDARD) then
                     self%comm_map_recv_ghost(rf, 5) = 0
                     comm_map_recv_ctr_ghost(neigh%myrank) = comm_map_recv_ctr_ghost(neigh%myrank) + &
                                                             self%grid%weight_neighbor(fec)
                  elseif (neighbor_type==NODE_MORE_REFINED) then
                     self%comm_map_recv_ghost(rf, 5) = n
                     comm_map_recv_ctr_ghost(neigh%myrank) = comm_map_recv_ctr_ghost(neigh%myrank) + &
                                                             self%grid%weight_neighbor(fec) / weight_reduction
                  elseif (neighbor_type==NODE_LESS_REFINED) then
                     self%comm_map_recv_ghost(rf, 5) = -neighbor_portion
                     comm_map_recv_ctr_ghost(neigh%myrank) = comm_map_recv_ctr_ghost(neigh%myrank) + &
                                                             self%grid%weight_neighbor(fec)
                  endif
                  call compute_ijk_min_max_delta(fec=fec, portion=self%comm_map_recv_ghost(rf, 5), &
                                                 ijk_min_max_delta=self%comm_map_recv_ghost(rf, 6:))
               elseif ((self%mpih%myrank == neigh%myrank).and.(self%mpih%myrank /= node_ptr%myrank)) then
                  sf = sf + 1
                  self%comm_map_send_ghost(sf, 15) = comm_map_send_ctr_ghost(node_ptr%myrank)
                  self%comm_map_send_ghost(sf, 1) = node_ptr%block_index
                  self%comm_map_send_ghost(sf, 2) = neigh%block_index
                  self%comm_map_send_ghost(sf, 3) = node_ptr%myrank
                  self%comm_map_send_ghost(sf, 4) = fec
                  if     (neighbor_type==NODE_STANDARD) then
                     self%comm_map_send_ghost(sf, 5) = 0
                     comm_map_send_ctr_ghost(node_ptr%myrank) = comm_map_send_ctr_ghost(node_ptr%myrank) + &
                                                                self%grid%weight_neighbor(fec)
                  elseif (neighbor_type==NODE_MORE_REFINED) then
                     self%comm_map_send_ghost(sf, 5) = n
                     comm_map_send_ctr_ghost(node_ptr%myrank) = comm_map_send_ctr_ghost(node_ptr%myrank) + &
                                                                self%grid%weight_neighbor(fec) / weight_reduction
                  elseif (neighbor_type==NODE_LESS_REFINED) then
                     self%comm_map_send_ghost(sf, 5) = -neighbor_portion
                     comm_map_send_ctr_ghost(node_ptr%myrank) = comm_map_send_ctr_ghost(node_ptr%myrank) + &
                                                                self%grid%weight_neighbor(fec)
                  endif
                  call compute_ijk_min_max_delta(fec=fec, portion=self%comm_map_send_ghost(sf, 5), &
                                                 ijk_min_max_delta=self%comm_map_send_ghost(sf, 6:))
               endif
            enddo
         endif
      enddo
   enddo
   contains
      subroutine compute_ijk_min_max_delta(fec, portion, ijk_min_max_delta)
      !< Compute IJK min/max/delta.
      integer(I4P), intent(in)  :: fec                    !< Current fec number.
      integer(I8P), intent(in)  :: portion                !< Current portion.
      integer(I8P), intent(out) :: ijk_min_max_delta(1:9) !< IJK min/max/delta.
      integer(I4P)              :: abs_portion            !< ABS of current portion.
      integer(I4P)              :: portion_cur            !< Current portion of fec updated (0=>whole fec).
      integer(I4P)              :: ijkmin(3)              !< Lower limit of ijk indexes.
      integer(I4P)              :: ijkmax(3)              !< Upper limit of ijk indexes.
      integer(I4P)              :: ijkdelta(3)            !< Delta offset for ghost-inner cells mapping same refinement.
      integer(I4P)              :: delta(3)               !< Neighbor delta of current fec.
      integer(I4P)              :: nijk(3)                !< Ni, Nj , Nk stored in array.
      integer(I4P)              :: portion_array(3)       !< Portion array binary useful for less refined case.
      integer(I4P)              :: i                      !< Counter.

      associate(ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk, ngc=>self%grid%ngc)
      nijk = [ni, nj, nk]

      abs_portion = abs(portion)
      delta = FEC_TO_DELTA(1:3, fec)

      if     (portion==0) then
         do i=1, 3
            if     (delta(i)==1) then
               ijkmin(i) = nijk(i) + 1
               ijkmax(i) = nijk(i) + ngc
               ijkdelta(i) = -nijk(i)
            elseif (delta(i)==-1) then
               ijkmin(i) = 1 - ngc
               ijkmax(i) = 0
               ijkdelta(i) = nijk(i)
            elseif (delta(i)==0) then
               ijkmin(i) = 1
               ijkmax(i) = nijk(i)
               ijkdelta(i) = 0
            endif
         enddo
      elseif (portion>0) then
         do i=1, 3
            if     (delta(i)==1) then
               ijkmin(i) = nijk(i) + 1
               ijkmax(i) = nijk(i) + ngc
               ijkdelta(i) = -1 - 2 * nijk(i)
            elseif (delta(i)==-1) then
               ijkmin(i) = 1 - ngc
               ijkmax(i) = 0
               ijkdelta(i) = -1 + nijk(i)
            elseif (delta(i)==0) then
               portion_cur = mod(abs_portion-1, 2)
               abs_portion = 2 * portion_cur + (abs_portion - 1) / 2 + 1
               ijkmin(i) = portion_cur * nijk(i) / 2 + 1
               ijkmax(i) = ijkmin(i) + nijk(i) / 2 - 1
               ijkdelta(i) = -2 * ijkmin(i) + 1
            endif
         enddo
      else
         portion_array(1) = mod(abs_portion-1,2)
         portion_array(2) = mod((abs_portion-1)/2,2)
         portion_array(3) = mod((abs_portion-1)/4,2)
         do i=1, 3
            if     (delta(i)==1) then
               ijkmin(i) = nijk(i) / 2 * portion_array(i) + 1
               ijkmax(i) = ijkmin(i) + ngc / 2 - 1
               ijkdelta(i) = - nijk(i) * portion_array(i) + nijk(i) - 1
            elseif (delta(i)==-1) then
               ijkmin(i) = nijk(i) / 2 + nijk(i) / 2 * portion_array(i) + 1 - ngc / 2
               ijkmax(i) = ijkmin(i) + ngc / 2 - 1
               ijkdelta(i) = - nijk(i) -  nijk(i) * portion_array(i) - 1
            elseif (delta(i)==0) then
               ijkmin(i) = nijk(i) / 2 * portion_array(i) + 1
               ijkmax(i) = ijkmin(i) + nijk(i) / 2 - 1
               ijkdelta(i) = - nijk(i) * portion_array(i) - 1
            endif
         enddo
      endif
      ijk_min_max_delta = [ijkmin, ijkmax, ijkdelta]
      endassociate
      endsubroutine
   endsubroutine make_comm_local_maps_ghost

   subroutine mpi_gather_nodes_data(self, node_member)
   !< Gather nodes data status between MPI processes.
   class(tree_object), intent(inout) :: self           !< The tree.
   character(*),       intent(in)    :: node_member    !< Node member to be shared.
   type(tree_node_object), pointer   :: node_ptr       !< Pointer to current node.
   integer(I8P), allocatable         :: send_buffer(:) !< Send buffer nodes data.
   integer(I8P), allocatable         :: recv_buffer(:) !< Recv buffer nodes data.
   integer(I4P), allocatable         :: recv_count(:)  !< Number of nodes that are received from each process.
   integer(I4P), allocatable         :: disp_count(:)  !< Displacement of nodes that are received from each process.
   integer(I8P)                      :: n, p           !< Counter.

   allocate(send_buffer(self%my_nodes_number * 2)) ! [Morton code, refinement_needed]
   allocate(recv_buffer(self%nodes_number    * 2)) ! [Morton code, refinement_needed]
   allocate(recv_count(0:self%mpih%procs_number - 1))
   allocate(disp_count(0:self%mpih%procs_number - 1))
   ! populating receive counters and send buffer
   recv_count = 0_I4P
   n = 0_I8P
   do while(self%loop(node_ptr=node_ptr))
      recv_count(node_ptr%myrank) = recv_count(node_ptr%myrank) + 2
      if (self%mpih%myrank == node_ptr%myrank) then
         n = n + 1
         send_buffer(n) = node_ptr%code
         n = n + 1
         select case(trim(node_member))
         case('refinement_needed')
            send_buffer(n) = node_ptr%refinement_needed
         case('block_index')
            send_buffer(n) = node_ptr%block_index
         endselect
      endif
   enddo
   ! computing displacement counts
   disp_count = 0_I4P
   do p=1, self%mpih%procs_number - 1
      disp_count(p) = disp_count(p-1) + recv_count(p-1)
   enddo

   call MPI_ALLGATHERV(send_buffer, self%my_nodes_number * 2, MPI_INTEGER8, &
                       recv_buffer, recv_count, disp_count, MPI_INTEGER8, MPI_COMM_WORLD, self%mpih%error)

   ! update nodes data
   do n=1, self%nodes_number*2, 2
      node_ptr => self%node(code=recv_buffer(n))
      select case(trim(node_member))
      case('refinement_needed')
         node_ptr%refinement_needed = int(recv_buffer(n+1), I4P)
      case('block_index')
         node_ptr%block_index = recv_buffer(n+1)
      endselect
   enddo
   endsubroutine mpi_gather_nodes_data

   subroutine mpi_print_stats(self)
   !< Print MPI stats.
   class(tree_object), intent(in) :: self !< The tree.
   integer(I4P)                   :: p    !< Counter.

   do p=0, self%mpih%procs_number-1
      print '(A)', self%mpih%myrankstr//' send to: '//trim(str(p,.true.))//' blocks n.: '//&
                   trim(str(self%comm_map_n_send(p),.true.))
   enddo
   if (allocated(self%comm_map_send)) &
      print '(A)', self%mpih%myrankstr//' blocks sent: '//trim(str(self%comm_map_send,.true.))
   do p=0, self%mpih%procs_number-1
      print '(A)', self%mpih%myrankstr//' recv from: '//trim(str(p,.true.))//' blocks n.: '//&
                   trim(str(self%comm_map_n_recv(p),.true.))
   enddo
   if (allocated(self%comm_map_recv)) &
      print '(A)', self%mpih%myrankstr//' blocks recv:  '//trim(str(self%comm_map_recv,.true.))
   if (allocated(self%local_map)) &
      print '(A)', self%mpih%myrankstr//' blocks kept n.: '//trim(str(size(self%local_map(:,1),dim=1),.true.))
   endsubroutine mpi_print_stats

   subroutine mpi_redistribute(self)
   !< Redistribute nodes to processes, load balancing.
   !<
   !< The nodes are distributed among all process exploiting the Morton ordering spatiality. Simply, the sorted list of codes are
   !< splitted in chunk of `nodes_number/procs_number` balancing the workload. However, the algorithm checks if the splits fall
   !< among siblings that cannot be splitted: if this scenario happens the siblings are placed in the same process for preserving
   !< the spatiality. The algorithm is sophisticated enough to place the siblings alternatively to the *left* or *right* process
   !< accordingly to where the split falls, namely to the *right* if the split falls in the first half of siblings or to the *left*
   !< if it falls in the second half.
   class(tree_object), intent(inout) :: self             !< The tree.
   integer(I8P), allocatable         :: codes_sorted(:)  !< List of (sorted) codes.
   type(tree_node_object), pointer   :: node_ptr         !< Pointer to current node.
   integer(I4P)                      :: p                !< Processes counter.
   integer(I4P)                      :: cl               !< Local child counter.
   integer(I8P)                      :: c                !< Codes counter.
   integer(I4P)                      :: i, j, k, l       !< Coordinates.
   integer(I8P)                      :: block_index_new  !< New block index counter.
   integer(I8P)                      :: my_codes_number  !< Number of codes for each process for a balanced workload.
   integer(I4P)                      :: child_local_code !< Local numbering.
   integer(I8P)                      :: n_keep           !< Number of keept nodes.
   integer(I8P)                      :: n_recv           !< Number of nodes that I have to receive.
   real(R8P)                         :: timing(2)        !< Tic toc timing.

   codes_sorted = self%codes() ! sorted list of codes
   my_codes_number = nint(real(size(codes_sorted, dim=1),R8P) / self%mpih%procs_number)
   ! initialize process rank and my codes number
   p = 0_I4P
   block_index_new = 0_I8P
   ! loop over all codes
   c = 1
   timing(1) = MPI_WTIME()
   do while(c<=size(codes_sorted, dim=1))
      if (block_index_new > my_codes_number.and.p < self%mpih%procs_number-1) then ! I would like to split...
         if (can_split()) then
            ! I am lucky, the split does not separate siblings
            p = p + 1_I4P
            block_index_new = 1_I8P
            node_ptr => self%node(code=codes_sorted(c))
            node_ptr%myrank_new = p
            node_ptr%block_index_new = block_index_new
         else
            ! I am not lucky, the split would separate siblings
            child_local_code = self%child_local(code=codes_sorted(c))
            if (child_local_code > self%ratio/2 -1) then
               ! the split falls in the second half of siblings list, place all nodes in the current process
               do cl=child_local_code, self%ratio-1
                  block_index_new = block_index_new + 1
                  node_ptr => self%node(code=codes_sorted(c+cl-child_local_code))
                  node_ptr%myrank_new = p
                  node_ptr%block_index_new = block_index_new
               enddo
            else
               ! the split falls in the first half of siblings list, place all nodes in the next process
               p = p + 1_I4P
               block_index_new = 0_I8P
               do cl=0, self%ratio-1
                  block_index_new = block_index_new + 1
                  node_ptr => self%node(code=codes_sorted(c+cl-child_local_code))
                  node_ptr%myrank_new = p
                  node_ptr%block_index_new = block_index_new
               enddo
            endif
            ! update codes counter skipping all current siblings
            c = c + self%ratio - 1 - child_local_code
         endif
      else ! no split, keeping to place nodes in the current process
         block_index_new = block_index_new + 1
         node_ptr => self%node(code=codes_sorted(c))
         node_ptr%myrank_new = p
         node_ptr%block_index_new = block_index_new
      endif
      c = c + 1
   enddo
   timing(2) = MPI_WTIME()
   ! print '(A, F18.10)', self%mpih%myrankstr//'step timing TREE%MPI_REDISTRIBUTE%newrank_loop: ', timing(2) - timing(1)

   ! create communication/local maps
   call self%make_comm_local_maps
   ! update tree status and compute coordinates of redistributed nodes
   n_keep = 0_I8P ; if (allocated(self%local_map    )) n_keep = size(self%local_map,     dim=1)
   n_recv = 0_I8P ; if (allocated(self%comm_map_recv)) n_recv = size(self%comm_map_recv, dim=1)
   if (allocated(self%block_coordinates)) deallocate(self%block_coordinates) ; allocate(self%block_coordinates(4, n_keep + n_recv))
   if (allocated(self%block_code)) deallocate(self%block_code) ; allocate(self%block_code(n_keep + n_recv))
   timing(1) = MPI_WTIME()
   do while(self%loop(node_ptr=node_ptr))
      node_ptr%myrank = node_ptr%myrank_new
      node_ptr%block_index = node_ptr%block_index_new
      if (node_ptr%myrank == self%mpih%myrank) then
         call self%morton_to_coordinates(code=node_ptr%code, i=i, j=j, k=k, l=l)
         self%block_coordinates(1, node_ptr%block_index) = i
         self%block_coordinates(2, node_ptr%block_index) = j
         self%block_coordinates(3, node_ptr%block_index) = k
         self%block_coordinates(4, node_ptr%block_index) = l
         self%block_code(node_ptr%block_index) = node_ptr%code
      endif
   enddo
   timing(2) = MPI_WTIME()
   ! print '(A, F18.10)', self%mpih%myrankstr//'step timing TREE%MPI_REDISTRIBUTE%create_coordinate: ', timing(2) - timing(1)

   self%last_block_index = n_keep + n_recv
   contains
      function can_split()
      !< Return true if the split can be done.
      !<
      !< The split is not allowed if all siblings exist and the previous code in the ordered list is one of my siblings.
      logical                   :: can_split   !< Result of test.
      integer(I8P), allocatable :: siblings(:) !< List of siblings
      integer(I4P)              :: s           !< Counter.

      can_split = .true.
      if (c==1_I8P) return
      siblings = self%siblings(code=codes_sorted(c))
      if (all([(self%has_code(code=siblings(s)), s=1,self%ratio-1)])) then ! if all siblings exist
         can_split = .not.(findloc(siblings, codes_sorted(c-1),dim=1) > 0) ! if my predecessor is a sibling the split is not allowed
      endif
      endfunction can_split
   endsubroutine mpi_redistribute

   ! Morton ordering methods
   pure function all_siblings(self, code) result(siblings)
   !< Return all siblings Morton code given Morton code (included into the list).
   class(tree_object), intent(in) :: self                   !< The tree.
   integer(I8P),       intent(in) :: code                   !< Morton code.
   integer(I8P)                   :: siblings(1:self%ratio) !< Siblings Morton codes [1:ratio].
   integer(I4P)                   :: local                  !< Local child code [0,ratio-1].
   integer(I4P)                   :: start                  !< Start code in the sibblings.
   integer(I4P)                   :: l, s                   !< Counter.

   if (code==-1) then
      siblings = -1_I8P ! anceestor of all has not siblings
   else
      local = self%child_local(code=code)
      start = code - local + 1
      s = 0
      do l=0, self%ratio - 1
         s = s + 1
         siblings(s) = start + l - 1
      enddo
   endif
   endfunction all_siblings

   elemental function child(self, code, i)
   !< Return the i-th child given Morton code.
   class(tree_object), intent(in) :: self  !< The tree.
   integer(I8P),       intent(in) :: code  !< Morton code.
   integer(I4P),       intent(in) :: i     !< Child index [0, ratio-1].
   integer(I8P)                   :: child !< Child Morton code.

   child = self%ratio * code + self%ratio + i
   endfunction child

   elemental function child_local(self, code) result(child)
   !< Return the child index in the local numbering.
   class(tree_object), intent(in) :: self  !< The tree.
   integer(I8P),       intent(in) :: code  !< Morton code.
   integer(I8P)                   :: child !< Child Morton code.

   child = 0
   if (code>0) child = int(code + self%ratio - ((code + self%ratio)/self%ratio)*self%ratio, I4P)
   endfunction child_local

   pure function children(self, code)
   !< Return the children given Morton code.
   class(tree_object), intent(in) :: self        !< The tree.
   integer(I8P),       intent(in) :: code        !< Morton code.
   integer(I8P), allocatable      :: children(:) !< Children Morton code.

   children = [self%child(code=code, i=0), &
               self%child(code=code, i=1), &
               self%child(code=code, i=2), &
               self%child(code=code, i=3), &
               self%child(code=code, i=4), &
               self%child(code=code, i=5), &
               self%child(code=code, i=6), &
               self%child(code=code, i=7)]
   endfunction children

   elemental function finest_at_level(self, code, level) result(finest)
   !< Return the inest node code at given level, namely the last child at a given level.
   class(tree_object), intent(in) :: self       !< The tree.
   integer(I8P),       intent(in) :: code       !< Morton code.
   integer(I4P),       intent(in) :: level      !< Refinement level.
   integer(I8P)                   :: finest     !< Finest Morton code.
   integer(I4P)                   :: code_level !< Level of the given code.
   integer(I4P)                   :: l          !< Counter.

   finest = code
   code_level = self%level(code=code)
   if (code_level<level) then
      finest = code
      do l=code_level+1, level
         finest = self%child(code=finest, i=self%ratio-1)
      enddo
   endif
   endfunction finest_at_level

   elemental function first_at_level(self, level) result(code)
   !< Return the first node code at given level.
   class(tree_object), intent(in) :: self  !< The tree.
   integer(I4P),       intent(in) :: level !< Refinement level.
   integer(I8P)                   :: code  !< Morton code.
   integer(I4P)                   :: l     !< Counter.

   code = -1
   if (level>0) then
      code = 0
      do l=2, level
         code = self%child(code=code, i=0_I4P)
      enddo
   endif
   endfunction first_at_level

   elemental function first_common_parent(self, code1, code2) result(fc_parent)
   !< Return the first common parent given two codes.
   class(tree_object), intent(in) :: self         !< The tree.
   integer(I8P),       intent(in) :: code1, code2 !< Morton codes.
   integer(I8P)                   :: fc_parent    !< First common parent.
   integer(I4P)                   :: level(2)     !< Levels of two codes.
   integer(I8P)                   :: parent(2)    !< Parents of two codes.
   logical                        :: is_found     !< Flag for searching exit.

   parent(1) = self%parent(code=code1)
   parent(2) = self%parent(code=code2)
   level(1) = self%level(code=parent(1))
   level(2) = self%level(code=parent(2))
   if     (level(1)>level(2)) then
      do while(self%level(code=parent(1))==level(2))
         parent(1) = self%parent(code=parent(1))
      enddo
   elseif (level(1)<level(2)) then
      do while(self%level(code=parent(2))==level(1))
         parent(2) = self%parent(code=parent(2))
      enddo
   endif
   is_found = parent(1) == parent(2)
   do while(.not.is_found)
      parent(1) = self%parent(code=parent(1))
      parent(2) = self%parent(code=parent(2))
      is_found = parent(1) == parent(2)
   enddo
   fc_parent = parent(1)
   endfunction first_common_parent

   subroutine get_neighbor(self, code, face, neighbor, neighbor_type)
   !< Return the neighbor in a given face of given Morton code.
   !<
   !< We define *direct neighbor* the neighbor of given code in the given face at the same level of the given code
   !< either if it exists or not.
   class(tree_object), intent(in)               :: self                   !< The tree.
   integer(I8P),       intent(in)               :: code                   !< Morton code.
   integer(I4P),       intent(in)               :: face                   !< Face queried.
   integer(I8P),       intent(out), allocatable :: neighbor(:)            !< Neighbors codes list, [1] or [ratio/2].
   integer(I4P),       intent(out)              :: neighbor_type          !< Type of neighbor.
   integer(I8P)                                 :: direct_neighbor        !< Morton code of direct neighbor.
   integer(I8P)                                 :: direct_neighbor_parent !< Morton code of direct neighbor parent.
   integer(I4P)                                 :: i_dn(4)                !< I coordinate of direct neighbor, or its 4 children.
   integer(I4P)                                 :: j_dn(4)                !< J coordinate of direct neighbor, or its 4 children.
   integer(I4P)                                 :: k_dn(4)                !< K coordinate of direct neighbor, or its 4 children.
   integer(I4P)                                 :: l_dn                   !< L coordinate of direct neighbor.
   integer(I4P)                                 :: i_dn_offset            !< I coordinate offset of direct neighbor l+1.
   integer(I4P)                                 :: j_dn_offset            !< J coordinate offset of direct neighbor l+1.
   integer(I4P)                                 :: k_dn_offset            !< K coordinate offset of direct neighbor l+1.

   ! compute coordinates of code
   call self%morton_to_coordinates(code=code, i=i_dn(1), j=j_dn(1), k=k_dn(1), l=l_dn)

   ! compute coordinates of direct neighbor and check if it falls outside the ancestor, in case
   ! it is a boundary condition node
   select case(face)
   case(1_I4P)
      i_dn(1) = i_dn(1) - 1
      if (i_dn(1) < 0) then
         neighbor_type = NODE_BOUNDARY_CONDITION
         return
      endif
   case(2_I4P)
      i_dn(1) = i_dn(1) + 1
      if (i_dn(1) > 2**l_dn - 1) then
         neighbor_type = NODE_BOUNDARY_CONDITION
         return
      endif
   case(3_I4P)
      j_dn(1) = j_dn(1) - 1
      if (j_dn(1) < 0) then
         neighbor_type = NODE_BOUNDARY_CONDITION
         return
      endif
   case(4_I4P)
      j_dn(1) = j_dn(1) + 1
      if (j_dn(1) > 2**l_dn - 1) then
         neighbor_type = NODE_BOUNDARY_CONDITION
         return
      endif
   case(5_I4P)
      k_dn(1) = k_dn(1) - 1
      if (k_dn(1) < 0) then
         neighbor_type = NODE_BOUNDARY_CONDITION
         return
      endif
   case(6_I4P)
      k_dn(1) = k_dn(1) + 1
      if (k_dn(1) > 2**l_dn - 1) then
         neighbor_type = NODE_BOUNDARY_CONDITION
         return
      endif
   endselect

   ! compute direct neighbor code
   select case(self%ratio)
   case(4_I4P)
      direct_neighbor = self%coordinates_to_morton(i=i_dn(1), j=j_dn(1), l=l_dn)
   case(8_I4P)
      direct_neighbor = self%coordinates_to_morton(i=i_dn(1), j=j_dn(1), k=k_dn(1), l=l_dn)
   endselect

   ! direct neighbor is not a sibling, check if it exists
   if (self%has_code(code=direct_neighbor)) then
      neighbor = [direct_neighbor]
      neighbor_type = NODE_STANDARD
      return
   endif

   ! direct neighbor does not exist, check if its parent exists
   direct_neighbor_parent = self%parent(code=direct_neighbor)
   if (self%has_code(code=direct_neighbor_parent)) then
      neighbor = [direct_neighbor_parent]
      neighbor_type = NODE_LESS_REFINED
      return
   endif

   ! direct neighbor parent does not exist, its children must exist or 2-1 rule has been broken
   ! using ijk coordinates at level l one can find the ijk coordinates of neighbor at level l+1
   l_dn = l_dn + 1
   i_dn_offset = (i_dn(1) - 1) * 2 + 1
   j_dn_offset = (j_dn(1) - 1) * 2 + 1
   k_dn_offset = (k_dn(1) - 1) * 2 + 1
   select case(face)
   case(1_I4P)
      i_dn(1) = i_dn_offset + 2
      j_dn(1) = j_dn_offset + 1
      k_dn(1) = k_dn_offset + 1

      i_dn(2) = i_dn_offset + 2
      j_dn(2) = j_dn_offset + 2
      k_dn(2) = k_dn_offset + 1

      i_dn(3) = i_dn_offset + 2
      j_dn(3) = j_dn_offset + 1
      k_dn(3) = k_dn_offset + 2

      i_dn(4) = i_dn_offset + 2
      j_dn(4) = j_dn_offset + 2
      k_dn(4) = k_dn_offset + 2
   case(2_I4P)
      i_dn(1) = i_dn_offset + 1
      j_dn(1) = j_dn_offset + 1
      k_dn(1) = k_dn_offset + 1

      i_dn(2) = i_dn_offset + 1
      j_dn(2) = j_dn_offset + 2
      k_dn(2) = k_dn_offset + 1

      i_dn(3) = i_dn_offset + 1
      j_dn(3) = j_dn_offset + 1
      k_dn(3) = k_dn_offset + 2

      i_dn(4) = i_dn_offset + 1
      j_dn(4) = j_dn_offset + 2
      k_dn(4) = k_dn_offset + 2
   case(3_I4P)
      i_dn(1) = i_dn_offset + 1
      j_dn(1) = j_dn_offset + 2
      k_dn(1) = k_dn_offset + 1

      i_dn(2) = i_dn_offset + 2
      j_dn(2) = j_dn_offset + 2
      k_dn(2) = k_dn_offset + 1

      i_dn(3) = i_dn_offset + 1
      j_dn(3) = j_dn_offset + 2
      k_dn(3) = k_dn_offset + 2

      i_dn(4) = i_dn_offset + 2
      j_dn(4) = j_dn_offset + 2
      k_dn(4) = k_dn_offset + 2
   case(4_I4P)
      i_dn(1) = i_dn_offset + 1
      j_dn(1) = j_dn_offset + 1
      k_dn(1) = k_dn_offset + 1

      i_dn(2) = i_dn_offset + 2
      j_dn(2) = j_dn_offset + 1
      k_dn(2) = k_dn_offset + 1

      i_dn(3) = i_dn_offset + 1
      j_dn(3) = j_dn_offset + 1
      k_dn(3) = k_dn_offset + 2

      i_dn(4) = i_dn_offset + 2
      j_dn(4) = j_dn_offset + 1
      k_dn(4) = k_dn_offset + 2
   case(5_I4P)
      i_dn(1) = i_dn_offset + 1
      j_dn(1) = j_dn_offset + 1
      k_dn(1) = k_dn_offset + 2

      i_dn(2) = i_dn_offset + 2
      j_dn(2) = j_dn_offset + 1
      k_dn(2) = k_dn_offset + 2

      i_dn(3) = i_dn_offset + 1
      j_dn(3) = j_dn_offset + 2
      k_dn(3) = k_dn_offset + 2

      i_dn(4) = i_dn_offset + 2
      j_dn(4) = j_dn_offset + 2
      k_dn(4) = k_dn_offset + 2
   case(6_I4P)
      i_dn(1) = i_dn_offset + 1
      j_dn(1) = j_dn_offset + 1
      k_dn(1) = k_dn_offset + 1

      i_dn(2) = i_dn_offset + 2
      j_dn(2) = j_dn_offset + 1
      k_dn(2) = k_dn_offset + 1

      i_dn(3) = i_dn_offset + 1
      j_dn(3) = j_dn_offset + 2
      k_dn(3) = k_dn_offset + 1

      i_dn(4) = i_dn_offset + 2
      j_dn(4) = j_dn_offset + 2
      k_dn(4) = k_dn_offset + 1
   endselect
   select case(self%ratio)
   case(4_I4P)
      neighbor = [self%coordinates_to_morton(i=i_dn(1), j=j_dn(1), l=l_dn), &
                  self%coordinates_to_morton(i=i_dn(2), j=j_dn(2), l=l_dn)]
      neighbor_type = NODE_MORE_REFINED
   case(8_I4P)
      neighbor = [self%coordinates_to_morton(i=i_dn(1), j=j_dn(1), k=k_dn(1), l=l_dn), &
                  self%coordinates_to_morton(i=i_dn(2), j=j_dn(2), k=k_dn(2), l=l_dn), &
                  self%coordinates_to_morton(i=i_dn(3), j=j_dn(3), k=k_dn(3), l=l_dn), &
                  self%coordinates_to_morton(i=i_dn(4), j=j_dn(4), k=k_dn(4), l=l_dn)]
      neighbor_type = NODE_MORE_REFINED
   endselect
   endsubroutine get_neighbor

   subroutine get_neighbor_all(self, code, face, neighbor, neighbor_type, neighbor_portion, neighbor_bc_fec)
   !< Return the neighbor in a given face/edge/corner of given Morton code.
   !<
   !< The direction `fec` is organized as: faces=[1,6], edges=[7,18], corners=[19,26].
   !<
   !< We define *direct neighbor* the neighbor of given code in the given face at the same level of the given code
   !< either if it exists or not.
   class(tree_object), intent(in)               :: self                        !< The tree.
   integer(I8P),       intent(in)               :: code                        !< Morton code.
   integer(I4P),       intent(in)               :: face                        !< Face queried.
   integer(I8P),       intent(out), allocatable :: neighbor(:)                 !< Neighbors codes list, [1] or [ratio/2].
   integer(I4P),       intent(out)              :: neighbor_type               !< Type of neighbor.
   integer(I4P),       intent(out), optional    :: neighbor_portion            !< Neighbors portion.
   integer(I4P),       intent(out), optional    :: neighbor_bc_fec             !< Neighbors fec for BC.
   integer(I8P)                                 :: direct_neighbor             !< Morton code of direct neighbor.
   integer(I8P)                                 :: direct_neighbor_parent      !< Morton code of direct neighbor parent.
   integer(I8P)                                 :: direct_neighbor_first_child !< Morton code of direct neighbor first child.
   integer(I4P)                                 :: i, j, k, l                  !< Counter.
   integer(I4P)                                 :: delta(3)                    !< Counter.
   integer(I4P)                                 :: ijk(3)                      !< Counter.
   integer(I4P)                                 :: ijkmin(3), ijkmax(3)        !< Counter.
   integer(I4P)                                 :: cl, cl_neighbor             !< Counter.
   integer(I4P)                                 :: cl_array(3)                 !< Counter.
   integer(I4P)                                 :: ijk_bc(3)                   !< Counter.
   integer(I4P)                                 :: ijk_size(3)                 !< Counter.

   if (present(neighbor_portion)) neighbor_portion = 1

   ! compute coordinates of code
   call self%morton_to_coordinates(code=code, i=i, j=j, k=k, l=l)

   ! compute coordinates of direct neighbor and check if it falls outside the ancestor, in case
   ! it is a boundary condition node
   delta = FEC_TO_DELTA(1:3, face)
   ijk = [i, j, k] + delta

   if (all(self%ijkl_prune>=0)) then
      ijk_size(1:3) = 2**(l - self%ijkl_prune(4)) * (self%ijkl_prune(1:3) + 1)
   else
      ijk_size = 2**l
   endif

   ! initialize ijk of direct neighbor node as a standard node
   neighbor_type = NODE_STANDARD
   ijk_bc = 0 ! standard node has ijk_bc = 0; BC node has ijk_bc = +-1
   do i=1, 3
      if (ijk(i)<0.or.ijk(i)>ijk_size(i) - 1) then
         ! ijk of neighbor node is outside ADAM, it is a BC except for periodic BC that is indeed a standard node
         ! check for periodic BC
         if (self%grid%is_ijk_periodic(i)) then
            ! direction ijk(i) is periodic, preserve the standard node initialization because direct neighbor is a standard node
            ijk(i) = modulo(ijk(i)+ijk_size(i), ijk_size(i)) ! ijk of neighbor node set to the direct neighbor in opposite direction
         else
            ! direction ijk(i) is not periodic, this is an actual BC
            ijk_bc(i) = sign(1, ijk(i))
            neighbor_type = NODE_BOUNDARY_CONDITION
         endif
      endif
   enddo
   if (neighbor_type == NODE_BOUNDARY_CONDITION) then
      if (present(neighbor_bc_fec)) neighbor_bc_fec = DELTA_TO_FEC(ijk_bc(1),ijk_bc(2),ijk_bc(3))
      return
   endif

   ! compute direct neighbor code
   direct_neighbor = self%coordinates_to_morton(i=ijk(1), j=ijk(2), k=ijk(3), l=l)

   ! direct neighbor is not a sibling, check if it exists
   if (self%has_code(code=direct_neighbor)) then
      neighbor = [direct_neighbor]
      neighbor_type = NODE_STANDARD
      return
   endif

   ! direct neighbor does not exist, check if its parent exists
   direct_neighbor_parent = self%parent(code=direct_neighbor)
   if (self%has_code(code=direct_neighbor_parent)) then
      neighbor = [direct_neighbor_parent]
      neighbor_type = NODE_LESS_REFINED
      if (present(neighbor_portion)) then
          neighbor_portion = self%child_local(code=direct_neighbor) + 1
      endif
      return
   endif

   ! direct neighbor parent does not exist, check its children
   ! using ijk coordinates at level l one can find the ijk coordinates of neighbor at level l+1
   delta = -delta ! thats magic
   allocate(neighbor(0))
   ! build the ratio/2 indexes of the children of direct_neighbor that are face-to-face with node
   do i=1, 3
      if     (delta(i)==1) then
         ijkmin(i) = 1
         ijkmax(i) = 1
      elseif (delta(i)==-1) then
         ijkmin(i) = 0
         ijkmax(i) = 0
      elseif (delta(i)==0) then
         ijkmin(i) = 0
         ijkmax(i) = 1
      endif
   enddo
   ! in ijkmin/ijkmax there are the ratio/2 deltas in binary notation that select the ratio/2 children of
   ! my interest
   do k=ijkmin(3), ijkmax(3)
      do j=ijkmin(2), ijkmax(2)
         do i=ijkmin(1), ijkmax(1)
            ! the first child (global) Morton code of direct neighbor has been already computed thus the other
            ! Morton codes are simply computed adding the local numeration to the first child code
            neighbor = [neighbor, self%child(code=direct_neighbor, i=i + 2*j + 4*k)]
            !RIMETTEREneighbor = [neighbor, [self%child(code=direct_neighbor, i=i + 2*j + 4*k)]]
         enddo
      enddo
   enddo
   neighbor_type = NODE_MORE_REFINED
   return
   endsubroutine get_neighbor_all

   elemental function greater(self, lhs, rhs) result(res)
   !< Return true if code is greater than other.
   class(tree_object), intent(in) :: self !< The tree.
   integer(I8P),       intent(in) :: lhs  !< Left hand side of code comparison.
   integer(I8P),       intent(in) :: rhs  !< Right hand side of code comparison.
   logical                        :: res  !< Comparison result.

   res = self%finest_at_level(code=lhs, level=self%max_level) > self%finest_at_level(code=rhs, level=self%max_level)
   endfunction greater

   subroutine is_greater_by_level(self, lhs, rhs, res)
   !< Return true if code is greater than other, comparing level.
   class(tree_object), intent(in)  :: self !< The tree.
   integer(I8P),       intent(in)  :: lhs  !< Left hand side of code comparison.
   integer(I8P),       intent(in)  :: rhs  !< Right hand side of code comparison.
   logical,            intent(out) :: res  !< Comparison result.

   res = self%level(code=lhs) > self%level(code=rhs)
   endsubroutine is_greater_by_level

   subroutine is_greater_by_pos(self, lhs, rhs, res)
   !< Return true if code is greater than other, comparing position.
   class(tree_object), intent(in)  :: self !< The tree.
   integer(I8P),       intent(in)  :: lhs  !< Left hand side of code comparison.
   integer(I8P),       intent(in)  :: rhs  !< Right hand side of code comparison.
   logical,            intent(out) :: res  !< Comparison result.

   res = self%finest_at_level(code=lhs, level=self%max_level) > self%finest_at_level(code=rhs, level=self%max_level)
   endsubroutine is_greater_by_pos

   subroutine is_lower_by_level(self, lhs, rhs, res)
   !< Return true if code is lower than other, comparing level.
   class(tree_object), intent(in)  :: self !< The tree.
   integer(I8P),       intent(in)  :: lhs  !< Left hand side of code comparison.
   integer(I8P),       intent(in)  :: rhs  !< Right hand side of code comparison.
   logical,            intent(out) :: res  !< Comparison result.

   res = self%level(code=lhs) < self%level(code=rhs)
   endsubroutine is_lower_by_level

   subroutine is_lower_by_pos(self, lhs, rhs, res)
   !< Return true if code is lower than other, comparing position.
   class(tree_object), intent(in)  :: self !< The tree.
   integer(I8P),       intent(in)  :: lhs  !< Left hand side of code comparison.
   integer(I8P),       intent(in)  :: rhs  !< Right hand side of code comparison.
   logical,            intent(out) :: res  !< Comparison result.

   res = self%finest_at_level(code=lhs, level=self%max_level) < self%finest_at_level(code=rhs, level=self%max_level)
   endsubroutine is_lower_by_pos

   elemental function last_at_level(self, level) result(code)
   !< Return the last node code at given level.
   class(tree_object), intent(in) :: self  !< The tree.
   integer(I4P),       intent(in) :: level !< Refinement level.
   integer(I8P)                   :: code  !< Morton code.

   code = self%first_at_level(level=level) + self%ratio**level - 1
   endfunction last_at_level

   elemental function level(self, code)
   !< Return the refinement level given the code.
   class(tree_object), intent(in) :: self  !< The tree.
   integer(I8P),       intent(in) :: code  !< Morton code.
   integer(I4P)                   :: level !< Refinement level.
   integer(I8P)                   :: c     !< Counter.

   if (code<=-1) then
      level = 0 ! ancestor of all has level 0
   else
      level = 1
      c = code
      do while(c>=self%ratio)
         c = (c - self%ratio) / self%ratio
         if (c>=0) level = level + 1
      enddo
   endif
   endfunction level

   elemental function lower(self, lhs, rhs) result(res)
   !< Return true if code is lower than other.
   class(tree_object), intent(in) :: self !< The tree.
   integer(I8P),       intent(in) :: lhs  !< Left hand side of code comparison.
   integer(I8P),       intent(in) :: rhs  !< Right hand side of code comparison.
   logical                        :: res  !< Comparison result.

   res = self%finest_at_level(code=lhs, level=self%max_level) < self%finest_at_level(code=rhs, level=self%max_level)
   endfunction lower

   elemental function parent(self, code)
   !< Return the parent given Morton code.
   class(tree_object), intent(in) :: self   !< The tree.
   integer(I8P),       intent(in) :: code   !< Morton code.
   integer(I8P)                   :: parent !< Parent Morton code.

   parent = -2_I8P ! ancestor of all has not parent
   if (code>self%ratio-1) parent = (code - self%ratio) / self%ratio
   endfunction parent

   elemental function parent_at_level(self, code, level) result(parent)
   !< Return the parent given Morton code at a given level.
   class(tree_object), intent(in) :: self     !< The tree.
   integer(I8P),       intent(in) :: code     !< Morton code.
   integer(I4P),       intent(in) :: level    !< Refinement level.
   integer(I8P)                   :: parent   !< Parent Morton code.
   integer(I8P), allocatable      :: path_(:) !< Parent Morton code.

   parent = -2_I8P ! ancestor of all nodes
   if (level>=1) then
      path_ = self%path(code=code)
      parent = path_(size(path_) - level + 1)
   endif
   endfunction parent_at_level

   pure function path(self, code)
   !< Return the path codes, the list of codes from given node to root.
   class(tree_object), intent(in) :: self    !< The tree.
   integer(I8P),       intent(in) :: code    !< Morton code.
   integer(I8P), allocatable      :: path(:) !< Path codes from node to root.
   integer(I8P), allocatable      :: temp(:) !< Temporary path list.
   integer(I8P)                   :: c       !< Counter.

   path = [code]
   c = code
   do while(self%level(code=c)>1)
      allocate(temp(1:size(path)+1))
      temp(1:size(path)) = path
      temp(size(path)+1) = self%parent(code=c)
      !RIMETTERE SENZA
      if(allocated(path)) deallocate(path)
      !RIMETTERE SENZA
      call move_alloc(from=temp,to=path)
      c = self%parent(code=c)
   enddo
   !RIMETTERE SENZA
   !EVITEREIif(allocated(temp)) deallocate(temp)
   !RIMETTERE SENZA
   endfunction path

   subroutine print_code_topology(self, code,  &
                                  coordinates, &
                                  level,       &
                                  parent,      &
                                  parents,     &
                                  path,        &
                                  child,       &
                                  child_local_code, &
                                  finest,      &
                                  siblings,    &
                                  neighbor,    &
                                  block_index, &
                                  whole)
   !< Print all code topology data.
   class(tree_object), intent(in)           :: self             !< The tree.
   integer(I8P),       intent(in)           :: code             !< Morton code.
   logical,            intent(in), optional :: coordinates      !< Coordinates.
   logical,            intent(in), optional :: level            !< Level of node.
   logical,            intent(in), optional :: parent           !< Parent of code.
   logical,            intent(in), optional :: parents          !< Parents list.
   logical,            intent(in), optional :: path             !< Path from node to parent of first level.
   logical,            intent(in), optional :: child            !< (First) Child of code.
   logical,            intent(in), optional :: child_local_code !< Local child-numbering of code.
   logical,            intent(in), optional :: finest           !< Finest Morton code.
   logical,            intent(in), optional :: siblings         !< Siblings of code.
   logical,            intent(in), optional :: neighbor         !< Neighbor of code.
   logical,            intent(in), optional :: block_index      !< Block index in the field array.
   logical,            intent(in), optional :: whole            !< Whole topology data.
   logical                                  :: coordinates_     !< Coordinates, local var.
   logical                                  :: level_           !< Level of node, local var.
   logical                                  :: parent_          !< Parent of code, local var.
   logical                                  :: child_           !< (First) Child of code, local var.
   logical                                  :: child_local_     !< Local child-numbering of code, local var.
   logical                                  :: finest_          !< Finest Morton code, local var.
   logical                                  :: siblings_        !< Siblings of code, local var.
   logical                                  :: path_            !< Path from node to parent of first level, local var.
   logical                                  :: parents_         !< Parents list, local var.
   logical                                  :: neighbor_        !< Neighbor of code, local var.
   logical                                  :: block_index_     !< Block index in the field array.
   logical                                  :: whole_           !< Whole topology data, local var.
   character(:), allocatable                :: topology         !< Topology string.
   character(:), allocatable                :: parents_str      !< Parents string list.
   character(:), allocatable                :: neighbors_str    !< Neighbors string.
   integer(I8P), allocatable                :: neighbors(:)     !< Neighbor of code.
   integer(I4P)                             :: neighbor_type    !< Type of neighbor.
   type(tree_node_object), pointer          :: node_ptr         !< Node pointer.
   integer(I4P)                             :: f                !< Counter.
   integer(I4P)                             :: i, j, k, l       !< Counter.

   coordinates_ = .false. ; if (present(coordinates     )) coordinates_ = coordinates
   level_       = .false. ; if (present(level           )) level_       = level
   parent_      = .false. ; if (present(parent          )) parent_      = parent
   child_       = .false. ; if (present(child           )) child_       = child
   child_local_ = .false. ; if (present(child_local_code)) child_local_ = child_local_code
   finest_      = .false. ; if (present(finest          )) finest_      = finest
   siblings_    = .false. ; if (present(siblings        )) siblings_    = siblings
   path_        = .false. ; if (present(path            )) path_        = path
   parents_     = .false. ; if (present(parents         )) parents_     = parents
   neighbor_    = .false. ; if (present(neighbor        )) neighbor_    = neighbor
   block_index_ = .false. ; if (present(block_index     )) block_index_ = block_index
   whole_       = .false. ; if (present(whole           )) whole_       = whole

   node_ptr => self%node(code=code)

   parents_str = ''
   do l=self%level(code=code)-1, 0, -1
      parents_str = parents_str//' pal_'//trim(str(l,.true.))//' '//trim(str(self%parent_at_level(code=code, level=l),.true.))
   enddo

   neighbors_str = ''
   do f=1, 6
      if (self%ratio==4_I4P.and.f>4) exit
      ! call self%get_neighbor_all(code=code, neighbor=neighbors, face=f, neighbor_type=neighbor_type)
      neighbors        = node_ptr%neighbor(f)%codes
      neighbor_type    = node_ptr%neighbor(f)%ntype
      if (allocated(neighbors)) then
         neighbors_str = neighbors_str//' f_'//trim(str(f,.true.))//' '//trim(str(neighbors,.true.))//&
                         ' type-'//trim(str(neighbor_type,.true.))
      else
         neighbors_str = neighbors_str//' f-'//trim(str(f,.true.))//' type_'//trim(str(neighbor_type,.true.))
      endif
   enddo

   call self%morton_to_coordinates(code=node_ptr%code, i=i, j=j, k=k, l=l)
   topology = ' code: '//trim(str(code))
   if (coordinates_.or.whole_) topology = topology//' coordinates: '//trim(str([i,j,k,l],.true.))
   if (level_.or.whole_      ) topology = topology//' level: '//trim(str(self%level(code=code),.true.))
   if (parent_.or.whole_     ) topology = topology//' parent: '//trim(str(self%parent(code=code),.true.))
   if (parents_.or.whole_    ) topology = topology//' parents: '//parents_str
   if (path_.or.whole_       ) topology = topology//' path: '//trim(str(self%path(code=code),.true.))
   if (child_.or.whole_      ) topology = topology//' child: '//trim(str(self%child(code=code, i=0),.true.))
   if (child_local_.or.whole_) topology = topology//' child_local: '//trim(str(self%child_local(code=code),.true.))
   if (finest_.or.whole_     ) topology = topology//' finest: '//trim(str(self%finest_at_level(code=code, &
                                                                                               level=self%max_level),.true.))
   if (siblings_.or.whole_   ) topology = topology//' siblings: '//trim(str(self%siblings(code=code),.true.))
   if (neighbor_.or.whole_   ) topology = topology//' neighbor: '//neighbors_str
   if (block_index_.or.whole_) topology = topology//' block_index: '//trim(str(node_ptr%block_index,.true.))

   print '(A)', self%mpih%myrankstr//topology
   endsubroutine print_code_topology

   pure function siblings(self, code)
   !< Return the siblings Morton code given Morton code.
   class(tree_object), intent(in) :: self                     !< The tree.
   integer(I8P),       intent(in) :: code                     !< Morton code.
   integer(I8P)                   :: siblings(1:self%ratio-1) !< Siblings Morton codes [1:ratio-1].
   integer(I4P)                   :: local                    !< Local child code [0,ratio-1].
   integer(I4P)                   :: start                    !< Start code in the sibblings.
   integer(I4P)                   :: l, s                     !< Counter.

   if (code==-1) then
      siblings = -1_I8P ! anceestor of all has not siblings
   else
      local = self%child_local(code=code)
      start = code - local + 1
      s = 0
      do l=0, self%ratio - 1
         if (l/=local) then
            s = s + 1
            siblings(s) = start + l - 1
         endif
      enddo
   endif
   endfunction siblings

   ! private methods
   subroutine add_node(self, code, refinement_needed, rank, block_index, update_last_block_index)
   !< Add a node pointer to the tree.
   !<
   !< @note If a node with the same key is already in the tree, it is removed and the new one will replace it.
   class(tree_object), intent(inout)        :: self                     !< The tree.
   integer(I8P),       intent(in)           :: code                     !< The Morton code.
   integer(I4P),       intent(in), optional :: refinement_needed        !< Flag for refinement/derefinement algorithm.
   integer(I4P),       intent(in), optional :: rank                     !< MPI rank process.
   integer(I8P),       intent(in), optional :: block_index              !< Block index in the field array.
   logical,            intent(in), optional :: update_last_block_index  !< Update or not last block index.
   logical                                  :: update_last_block_index_ !< Update or not last block index, local var.
   integer(I4P)                             :: b                        !< Bucket index, namely hashed key.
   integer(I4P)                             :: rank_                    !< MPI rank process, local variable.

   if (.not.self%is_initialized_) then
      print '(A)', self%mpih%myrankstr//'ERROR: cannot add a node a non initialized tree'
   endif

   rank_ = 0 ; if (present(rank)) rank_ = rank

   ! if the code is not already in the tree update the nodes number otherwise not
   if (.not.self%has_code(code=code)) self%nodes_number = self%nodes_number + 1
   b = self%hash(code=code)
   call self%bucket(b)%add_node(code=code, refinement_needed=refinement_needed, &
                                myrank=rank_, block_index=block_index)

   if (self%mpih%myrank == rank_) then
      update_last_block_index_ = .true. ; if (present(update_last_block_index)) update_last_block_index_ = update_last_block_index
      if (update_last_block_index_) self%last_block_index = self%last_block_index + 1
   endif
   endsubroutine add_node

   function coordinates2D_to_morton(self, i, j, l) result(code)
   !< Return the Morton code given ijl coordinates.
   class(tree_object), intent(in) :: self  !< The tree.
   integer(I4P),       intent(in) :: i     !< I coordinate.
   integer(I4P),       intent(in) :: j     !< J coordinate.
   integer(I4P),       intent(in) :: l     !< L coordinate.
   integer(I8P)                   :: code  !< Morton code.

   code = self%first_at_level(level=l) + morton2D(i=i, j=j)
   endfunction coordinates2D_to_morton

   function coordinates3D_to_morton(self, i, j, k, l) result(code)
   !< Return the Morton code given ijkl coordinates.
   class(tree_object), intent(in) :: self  !< The tree.
   integer(I4P),       intent(in) :: i     !< I coordinate.
   integer(I4P),       intent(in) :: j     !< J coordinate.
   integer(I4P),       intent(in) :: k     !< K coordinate.
   integer(I4P),       intent(in) :: l     !< L coordinate.
   integer(I8P)                   :: code  !< Morton code.

   code = self%first_at_level(level=l) + morton3D(i=i, j=j, k=k)
   endfunction coordinates3D_to_morton

   subroutine derefine(self)
   !< Derefine nodes.
   class(tree_object), intent(inout) :: self             !< The tree.
   type(tree_node_object), pointer   :: first_child      !< Pointer to first child node.
   type(tree_node_object), pointer   :: node_ptr         !< Pointer to node.
   integer(I8P)                      :: derefined_number !< Number of derefined blocks.
   integer(I8P)                      :: mn               !< Counter.
   integer(I8P)                      :: n                !< Counter.
   integer(I4P)                      :: i                !< Counter.

   if (allocated(self%block_to_derefine)) deallocate(self%block_to_derefine)
   if (allocated(self%block_derefined)) deallocate(self%block_derefined)
   derefined_number = size(self%node_to_derefine, dim=1)
   if (derefined_number>0) then
      allocate(self%block_to_derefine(self%n_my_derefine))
      allocate(self%block_derefined(2, self%n_my_derefine/self%ratio))
      mn = -7
      do n=1, derefined_number, self%ratio
         first_child => self%node(code=self%node_to_derefine(n))
         if (self%mpih%myrank == first_child%myrank) then
            mn = mn + 8
            self%block_derefined(1,(mn-1)/self%ratio+1) = self%parent(code=first_child%code)
            self%block_derefined(2,(mn-1)/self%ratio+1) = first_child%block_index
         endif
         call self%add_node(code=self%parent(code=first_child%code),              &
                                             rank=first_child%myrank,             &
                                             block_index=first_child%block_index, &
                                             update_last_block_index=.false.)
         do i=0, self%ratio - 1
            node_ptr => self%node(code=self%node_to_derefine(n+i))
            if (self%mpih%myrank == node_ptr%myrank) then
               self%block_to_derefine(mn+i) = node_ptr%block_index
            endif
            call self%remove_node(code=self%node_to_derefine(n+i))
         enddo
      enddo
   endif
   endsubroutine derefine

   subroutine empty(self)
   !< Empty tree, i.e. remove all nodes, without destroy the tree.
   class(tree_object), intent(inout) :: self !< The tree.
   integer(I4P)                      :: b    !< Counter.

   do b=lbound(self%bucket, dim=1), ubound(self%bucket, dim=1)
      call self%bucket(b)%destroy
   enddo
   self%nodes_number = 0_I4P
   endsubroutine empty

   subroutine morton_to_coordinates2D(self, code, i, j, l)
   !< Return the ijkl coordinates given Morton code.
   class(tree_object), intent(in)  :: self     !< The tree.
   integer(I8P),       intent(in)  :: code     !< Morton code.
   integer(I4P),       intent(out) :: i        !< I coordinate.
   integer(I4P),       intent(out) :: j        !< J coordinate.
   integer(I4P),       intent(out) :: l        !< L coordinate.
   integer(I4P)                    :: ij(2)    !< IJ local coordinates.
   integer(I8P), allocatable       :: path_(:) !< Path from code to root.
   integer(I4P)                    :: p        !< Counter.
   integer(I8P)                    :: c        !< Counter.

   i = 0_I4P
   j = 0_I4P
   l = self%level(code=code)
   path_ = self%path(code=code)
   do p=1, size(path_, dim=1)
      c = path_(p)
      call demorton2D(code=self%child_local(code=c), i=ij(1), j=ij(2))
      i = i + ij(1) * 2**(p-1)
      j = j + ij(2) * 2**(p-1)
   enddo
   endsubroutine morton_to_coordinates2D

   subroutine morton_to_coordinates3D(self, code, i, j, k, l)
   !< Return the ijkl coordinates given Morton code.
   class(tree_object), intent(in)  :: self     !< The tree.
   integer(I8P),       intent(in)  :: code     !< Morton code.
   integer(I4P),       intent(out) :: i        !< I coordinate.
   integer(I4P),       intent(out) :: j        !< J coordinate.
   integer(I4P),       intent(out) :: k        !< K coordinate.
   integer(I4P),       intent(out) :: l        !< L coordinate.
   integer(I4P)                    :: ijk(3)   !< IJK local coordinates.
   integer(I8P), allocatable       :: path_(:) !< Path from code to root.
   integer(I4P)                    :: p        !< Counter.
   integer(I8P)                    :: c        !< Counter.

   i = 0_I4P
   j = 0_I4P
   k = 0_I4P
   l = self%level(code=code)
   path_ = self%path(code=code)
   do p=1, size(path_, dim=1)
      c = path_(p)
      call demorton3D(code=self%child_local(code=c), i=ijk(1), j=ijk(2), k=ijk(3))
      i = i + ijk(1) * 2**(p-1)
      j = j + ijk(2) * 2**(p-1)
      k = k + ijk(3) * 2**(p-1)
   enddo
   endsubroutine morton_to_coordinates3D

   subroutine refine(self)
   !< Refine nodes.
   class(tree_object), intent(inout) :: self           !< The tree.
   type(tree_node_object), pointer   :: parent         !< Pointer to parent node.
   integer(I8P)                      :: refined_number !< Number of nodes to be refined.
   integer(I8P)                      :: mn             !< Counter.
   integer(I8P)                      :: n              !< Counter.
   integer(I4P)                      :: i              !< Counter.

   if (allocated(self%block_to_refine)) deallocate(self%block_to_refine)
   if (allocated(self%block_refined)) deallocate(self%block_refined)
   refined_number = size(self%node_to_refine, dim=1)
   if (refined_number>0) then
      allocate(self%block_to_refine(2, self%n_my_refine))
      allocate(self%block_refined(2, self%ratio*self%n_my_refine))

      self%block_refined = -100

      mn = 0
      do n=1, refined_number
         parent => self%node(code=self%node_to_refine(n))
         if (parent%myrank == self%mpih%myrank) then
            mn = mn + 1
            self%block_to_refine(1,mn) = parent%block_index
            self%block_to_refine(2,mn) = parent%myrank
         endif
         call self%add_node(code=self%child(code=parent%code, i=0), rank=parent%myrank, &
                            block_index=parent%block_index, update_last_block_index=.false.)
         if (parent%myrank == self%mpih%myrank) then
            self%block_refined(1, (mn-1)*self%ratio+1) = self%child(code=parent%code, i=0)
            self%block_refined(2, (mn-1)*self%ratio+1) = parent%block_index
         endif
         do i=1, self%ratio-1
            if (parent%myrank == self%mpih%myrank) then
               self%block_refined(1, (mn-1)*self%ratio+1+i) = self%child(code=parent%code, i=i)
               self%block_refined(2, (mn-1)*self%ratio+1+i) = self%last_block_index + 1
            endif
            call self%add_node(code=self%child(code=parent%code, i=i), rank=parent%myrank, &
                               block_index=self%last_block_index+1)
         enddo
         call self%remove_node(code=parent%code)
      enddo
   endif
   endsubroutine refine

   subroutine sanitize(self, iterations_number)
   !< Sanitize the tree.
   class(tree_object),        intent(inout)        :: self                 !< The tree.
   integer(I4P),              intent(in), optional :: iterations_number    !< Sanitazie iterations number.
   integer(I4P)                                    :: iterations_number_   !< Sanitazie iterations number.
   type(tree_node_object), pointer                 :: node_ptr             !< Pointer to node.
   type(tree_node_object), pointer                 :: sibling              !< Pointer to node sibling.
   integer(I8P)                                    :: code                 !< Code.
   integer(I8P), allocatable                       :: siblings(:)          !< List of code siblings, excluded the quering code.
   integer(I8P), allocatable                       :: all_siblings(:)      !< List of code siblings, included the quering code.
   integer(I8P), allocatable                       :: neighbor(:)          !< List of code neighbors.
   type(tree_node_object), pointer                 :: neigh                !< Pointer to node neighbor.
   integer(I4P)                                    :: neighbor_type        !< Neighbors type.
   logical                                         :: is_sanitize_complete !< Flag for finishing sanitize.
   logical                                         :: can_be_derefined     !< Flag for checking derefinement possibility.
   integer(I8P), allocatable                       :: codes_analyzed(:)    !< List of codes analyzed.
   integer(I4P)                                    :: new_level            !< New level counter.
   integer(I4P)                                    :: new_level_n          !< Neighbor new level counter.
   integer(I4P)                                    :: s, sib, f, n         !< Counter.
   real(R8P)                                       :: timing(4)            !< Tic toc timing.

   iterations_number_ = TREE_MAX_SANITIZE_ITERATIONS ; if (present(iterations_number)) iterations_number_ = iterations_number

   timing(1) = MPI_WTIME()
   min_max_check_loop : do while(self%loop(node_ptr=node_ptr))
      new_level = self%level(code=node_ptr%code) + node_ptr%refinement_needed
      if ((new_level > self%max_level).or.(new_level < 0)) then
         node_ptr%refinement_needed = TO_NOT_TOUCH
      endif
   enddo min_max_check_loop
   timing(2) = MPI_WTIME()
   ! print '(A, F18.10)', self%mpih%myrankstr//'step timing TREE%SANITIZE%min_max_check_loop: ', timing(2) - timing(1)

   !MEMORYLEAK LEVAREallocate(self%temp_array_i8(10000))

   timing(1) = MPI_WTIME()
   sanitize_loop : do s=1, iterations_number_
      print '(A)', self%mpih%myrankstr//'sanitize iteration '//trim(str(s,.true.))
      is_sanitize_complete = .true.

      ! check for the sanity of derefinement
      self%n_my_derefine = 0
      if (allocated(self%node_to_derefine)) deallocate(self%node_to_derefine) ; allocate(self%node_to_derefine(0))
      if (allocated(codes_analyzed)) deallocate(codes_analyzed) ; allocate(codes_analyzed(0))
      !MEMORYLEAKLEVAREif (allocated(self%codes_analyzed)) deallocate(self%codes_analyzed) ; allocate(self%codes_analyzed(0))
      timing(3) = MPI_WTIME()
      derefine_loop : do while(self%loop(node_ptr=node_ptr))
         ! check if I want to be derefined and I have not been analyzed yet
         if (node_ptr%refinement_needed == TO_BE_DEREFINED) then
            if (findloc(codes_analyzed, node_ptr%code, dim=1)==0) then ! avoid to re-analyze already confirmed siblingsi to derefine
            !MEMORYLEAKLEVAREif (findloc(self%codes_analyzed, node_ptr%code, dim=1)==0) then ! avoid to re-analyze already confirmed siblingsi to derefine
               ! check sibling for derefinement
               can_be_derefined = .true.
               code = node_ptr%code
               siblings = self%siblings(code=code)
               sibs_check_loop : do sib=1, self%ratio -1
                  if (.not.self%has_code(code=siblings(sib))) then
                     can_be_derefined = .false.
                     exit sibs_check_loop
                  endif
                  sibling => self%node(code=siblings(sib))
                  if (sibling%refinement_needed /= TO_BE_DEREFINED) then
                     can_be_derefined = .false.
                     exit sibs_check_loop
                  endif
               enddo sibs_check_loop
               if (can_be_derefined) then
                  all_siblings = self%all_siblings(code=code)

                  self%node_to_derefine = [self%node_to_derefine, all_siblings]
                  codes_analyzed = [codes_analyzed, all_siblings]

                  if (self%mpih%myrank==node_ptr%myrank) self%n_my_derefine = self%n_my_derefine + 8
               else
                  is_sanitize_complete = .false.
                  node_ptr%refinement_needed = TO_NOT_TOUCH
                  do sib=1, self%ratio -1
                     if (self%has_code(code=siblings(sib))) then
                        sibling => self%node(code=siblings(sib))
                        if (sibling%refinement_needed == TO_BE_DEREFINED) then
                           ! due some of your siblings you cannot be derefined, you need to be altered
                           sibling%refinement_needed = TO_NOT_TOUCH
                        endif
                     endif
                  enddo
               endif
            endif
         endif
      enddo derefine_loop
      timing(4) = MPI_WTIME()
      ! print '(A, F18.10)', self%mpih%myrankstr//'step timing TREE%SANITIZE%sanitize_loop%derefine_loop: ', timing(4) - timing(3)

      ! check for the sanity of refinement (2:1 rule)
      timing(3) = MPI_WTIME()
      refine_loop : do while(self%loop(node_ptr=node_ptr))
         new_level = self%level(code=node_ptr%code) + node_ptr%refinement_needed
         face_loop : do f=1, 26
            ! call self%get_neighbor_all(code=node_ptr%code, face=f, neighbor=neighbor, neighbor_type=neighbor_type)
            neighbor         = node_ptr%neighbor(f)%codes
            neighbor_type    = node_ptr%neighbor(f)%ntype
            if (neighbor_type /= NODE_BOUNDARY_CONDITION) then
               neighbor_loop : do n=1, size(neighbor, dim=1)
                  ! check level
                  neigh => self%node(code=neighbor(n))
                  new_level_n = self%level(code=neighbor(n)) + neigh%refinement_needed
                  if (new_level_n > new_level + 1) then
                     ! a neighbour want to be refined 2 levels more than me, I have to refine more too
                     is_sanitize_complete = .false.
                     if     (new_level_n - new_level == 3) then ! node want to derefine, but it must be refined
                        node_ptr%refinement_needed = 1
                     elseif (new_level_n - new_level == 2) then
                        node_ptr%refinement_needed = node_ptr%refinement_needed + 1
                     else
                        print '(A)', self%mpih%myrankstr//'SOMETHING WENT TERRIBLY WRONG. EXIT!'
                        print '(A)', self%mpih%myrankstr//'REFINEMENT NEEDED '//trim(str(node_ptr%refinement_needed,.true.))
                        print '(A)', self%mpih%myrankstr//'SANITIZE ITERATIONS '//trim(str(s,.true.))
                        stop
                     endif
                     new_level = self%level(code=node_ptr%code) + node_ptr%refinement_needed
                  endif
               enddo neighbor_loop
            endif
         enddo face_loop

        if (node_ptr%refinement_needed > 1) then
           print '(A)', self%mpih%myrankstr//'CANNOT REFINE TWICE IN A ROW. SOMETHING WENT TERRIBLY WRONG. EXIT!'
           print '(A)', self%mpih%myrankstr//'SANITIZE ITERATIONS '//trim(str(s,.true.))
           stop
        endif

        new_level = self%level(code=node_ptr%code) + node_ptr%refinement_needed
        if (new_level > self%max_level) then
           print '(A)', self%mpih%myrankstr//'CANNOT REFINE MORE. SOMETHING WENT TERRIBLY WRONG. EXIT!'
           stop
        endif
      enddo refine_loop
      timing(4) = MPI_WTIME()
      ! print '(A, F18.10)', self%mpih%myrankstr//'step timing TREE%SANITIZE%sanitize_loop%refine_loop: ', timing(4) - timing(3)

      if (is_sanitize_complete) exit sanitize_loop
   enddo sanitize_loop
   timing(2) = MPI_WTIME()
   ! print '(A, F18.10)', self%mpih%myrankstr//'step timing TREE%SANITIZE%sanitize_loop: ', timing(2) - timing(1)

   if (.not.is_sanitize_complete) then
      print '(A)', self%mpih%myrankstr//'SANITZE CANNOT BE COMPLETED. SOMETHING WENT TERRIBLY WRONG. EXIT!'
      stop
   endif

   ! update to_refine list
   self%n_my_refine = 0
   if (allocated(self%node_to_refine)) deallocate(self%node_to_refine) ; allocate(self%node_to_refine(0))
   timing(1) = MPI_WTIME()
   do while(self%loop(node_ptr=node_ptr))
      if (node_ptr%refinement_needed==TO_BE_REFINED) then
         self%node_to_refine = [self%node_to_refine, [node_ptr%code]]
         if (self%mpih%myrank==node_ptr%myrank) self%n_my_refine = self%n_my_refine + 1
      endif
   enddo
   timing(2) = MPI_WTIME()
   ! print '(A, F18.10)', self%mpih%myrankstr//'step timing TREE%SANITIZE%to_refine_loop: ', timing(2) - timing(1)
   endsubroutine sanitize
endmodule adam_tree_object
