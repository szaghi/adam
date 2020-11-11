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
!<  |                            *----------*----------*
!<  |                           /|         /|         /|
!<  |                          / |        / |        / |
!<  |                         /  |       /  |       /  |
!<  |                        /   |      /   |      /   |
!<  |                       /    *-----/----*-----/----*
!<  |            6 <-------/--+ /|    /    /|    /   +---------> 7
!<  |                     *----------*----------*    / |
!<  |                    /|   /  |  /|   /  |  /|   /  |
!<  |                   / |  /   | / |  /   | / |  /   |
!<  |                  /  | /    */--|-/----*/--|-/----*
!<  |         2 <-----/---|/--+ //   |/    //   |/  +---------> 3
!<  |                /    *-----/----*-----/----*    /
!<  |    4 <--------/--+ /|   //    /|   //   +------------> 5
!<  |              *----------*----------*    / |  /
!<  |              |   /  | / |   /  | / |   /  | /
!<  |              |  /   |/  |  /   |/  |  /   |/
!<  |              | /    *---|-/----*---|-/----*
!<  |      0 <-----|/--+ /    |/    /    |/  +------------> 1
!<  |              *----------*----------*    /
!<  |              |   /      |   /      |   /
!<  |              |  /       |  /       |  /
!<  |   _ Y        | /        | /        | /
!<  |   /|         |/         |/         |/
!<  |  /           *----------*----------*
!<  | /
!<  |/                                                    X
!<  o------------------------------------------------------------------->

use adam_grid_object
use adam_parameters
use adam_tree_node_object
use adam_tree_bucket_object
use MORTIF
use PENF
#ifdef _MPI_
use MPI
#endif
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

implicit none
private
public :: tree_object

! tree node parameters
integer(I4P), parameter :: NODE_STANDARD = 0_I4P           !< Standard node type.
integer(I4P), parameter :: NODE_MORE_REFINED = 1_I4P       !< More refined node type.
integer(I4P), parameter :: NODE_BOUNDARY_CONDITION = 2_I4P !< Boundary condition node type.
! tree parameters
integer(I8P), parameter :: TREE_BUCKETS_NUMBER_DEF = 9973_I8P    !< Default number of buckets of hash table.
real(R8P),    parameter :: TREE_MAX_LOAD = 0.9_R8P               !< Maximum load of hash table buckets.
integer(I4P), parameter :: TREE_MAX_SANITIZE_ITERATIONS = 10_I4P !< Default number of tree sanitize iterations.

type :: tree_object
   !< Tree class definition.
   ! tree data
   type(grid_object), pointer            :: grid=>null()            !< Grid data.
   type(tree_bucket_object), allocatable :: bucket(:)               !< Tree buckets.
   integer(I8P)                          :: buckets_number=0_I8P    !< Number of buckets used.
   integer(I4P)                          :: nodes_number=0_I4P      !< Number of nodes actually stored, namely the tree length.
   real(R8P)                             :: max_load=TREE_MAX_LOAD  !< Maximum load of tree buckets.
   integer(I4P)                          :: ratio=8_I4P             !< Refinement ratio.
   integer(I4P)                          :: max_level=12_I4P        !< Maximum refinement level.
   logical                               :: is_initialized_=.false. !< Initialization status.
   ! AMR data
   integer(I8P)              :: last_block_index=0_I8P !< Last block index in the field array.
   integer(I8P), allocatable :: node_to_refine(:)      !< List of nodes to be refined.
   integer(I8P), allocatable :: node_to_derefine(:)    !< List of nodes to be derefined.
   integer(I8P), allocatable :: block_to_refine(:,:)   !< List of field blocks to be refined.
   integer(I8P), allocatable :: block_refined(:,:)     !< List of field refined blocks with Morton code.
   integer(I8P), allocatable :: block_to_derefine(:)   !< List of field blocks to be derefined.
   integer(I8P), allocatable :: block_derefined(:,:)   !< List of field derefined blocks with Morton code.
   integer(I4P), allocatable :: block_coordinates(:,:) !< Block coordinates of redistributed blocks [4,blocks_number].
   ! MPI data
   integer(I4P)              :: procs_number=1_I4P      !< MPI Number of processes.
   integer(I4P)              :: myrank=0_I4P            !< MPI rank process.
   integer(I4P)              :: my_nodes_number=0_I4P   !< Number of my nodes, keep_nodes_number + recv_nodes_number.
   integer(I4P)              :: send_nodes_number=0_I4P !< Number of nodes to be sent.
   integer(I4P)              :: recv_nodes_number=0_I4P !< Number of nodes to be received.
   integer(I4P)              :: keep_nodes_number=0_I4P !< Number of nodes to be keept.
   integer(I4P), allocatable :: comm_map_n_send(:)      !< Communication map, number of blocks to send [procs_number].
   integer(I4P), allocatable :: comm_map_n_recv(:)      !< Communication map, number of blocks to recv [procs_number].
   integer(I4P), allocatable :: comm_map_send_ptr(:)    !< Communication map, pointers in list to send [procs_number+1].
   integer(I4P), allocatable :: comm_map_recv_ptr(:)    !< Communication map, pointers in list to recv [procs_number+1].
   integer(I8P), allocatable :: comm_map_send(:)        !< Communication map, blocks to send [sum(comm_map_n_send)].
   integer(I8P), allocatable :: comm_map_recv(:)        !< Communication map, blocks to receive [sum(comm_map_n_recv)].
   integer(I8P), allocatable :: local_map(:,:)          !< Local map, list block index changes of my nodes.
   contains
      ! public methods
      procedure, pass(self) :: adapt                !< Adapt tree accordingly to refine/derefine necessity.
      procedure, pass(self) :: codes                !< Return the list of (sorted) codes actually stored in the tree.
      procedure, pass(self) :: destroy              !< Destroy the tree.
      procedure, pass(self) :: loop                 !< Sentinel while-loop on nodes returning the code.
      procedure, pass(self) :: hash                 !< Hash the key.
      procedure, pass(self) :: has_code             !< Check if the code is present in the tree.
      procedure, pass(self) :: initialize           !< Initialize the tree.
      procedure, pass(self) :: mark_all_nodes       !< Mark all nodes to be refined, derefined....
      procedure, pass(self) :: node                 !< Return a pointer to a node.
      procedure, pass(self) :: prime_buckets_number !< Return the buckets number as nearest prime number given nodes number.
      procedure, pass(self) :: resize               !< Resize the tree.
      procedure, pass(self) :: sanitize             !< Sanitize the tree.
      procedure, pass(self) :: traverse             !< Traverse tree calling the iterator procedure.
      ! MPI methods
      procedure, pass(self) :: import_refinements_needed     !< Import refinements needed status changed externally.
      procedure, pass(self) :: make_comm_local_maps          !< Make communication/local maps.
      procedure, pass(self) :: mpi_gather_refinements_needed !< Gather refinements needed status between MPI processes.
      procedure, pass(self) :: mpi_print_stats               !< Print MPI stats.
      procedure, pass(self) :: mpi_redistribute              !< Redistribute nodes to MPI processes, load balancing.
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
      procedure, pass(self) :: last_at_level           !< Return the last node code at given level.
      procedure, pass(self) :: level                   !< Return the refinement level given the code.
      procedure, pass(self) :: lower                   !< Return true if code is lower than other.
      procedure, pass(self) :: get_neighbor            !< Return the neighbor in a given face of given Morton code.
      procedure, pass(self) :: greater                 !< Return true if code is greater than other.
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
      procedure, pass(self), private :: morton_to_coordinates3D !< Return the ijkl coordinates given Morton code.
      procedure, pass(self), private :: morton_to_coordinates2D !< Return the ijkl coordinates given Morton code.
      procedure, pass(self), private :: refine                  !< Refine nodes.
      procedure, pass(self), private :: remove_node             !< Remove a node from the tree, given the key.
      ! operators
      generic :: assignment(=) => tree_assign_tree      !< Overload `=`.
      procedure, pass(lhs), private :: tree_assign_tree !< Operator `=`.
endtype tree_object

contains
   ! public methods
   subroutine adapt(self)
   !< Adapt tree accordingly to refine/derefine necessity.
   class(tree_object), intent(inout) :: self !< The tree.

   call self%sanitize
   call self%refine
   call self%derefine
   endsubroutine adapt

   function codes(self, only_mine)
   !< Return the list of (sorted) codes actually stored in the tree.
   class(tree_object), intent(in)           :: self       !< The tree.
   logical,            intent(in), optional :: only_mine  !< If true return only the nodes of myrank process.
   integer(I8P), allocatable                :: codes(:)   !< List of codes.
   logical                                  :: only_mine_ !< If true return only the nodes of myrank process.
   type(tree_node_object), pointer          :: node       !< Pointer to current node.
   integer(I8P)                             :: c          !< Counter.
   integer(I8P), allocatable                :: work(:)    !< Working memory for sorting codes list.

   only_mine_ = .false. ; if (present(only_mine)) only_mine_ = only_mine
   allocate(codes(self%nodes_number))
   c = 0
   do while(self%loop(node=node))
      if (only_mine_.and.self%myrank/=node%myrank) cycle
      c = c + 1
      codes(c) = node%code
   enddo
   if (c < self%nodes_number) then
      work = codes(1:c)
      call move_alloc(from=work, to=codes)
   endif
   allocate(work((size(codes,dim=1)+1)/2))
   call mergesort(array=codes)
   contains
      recursive subroutine mergesort(array)
      !< Sort input array by means of mergesort algorithm.
      integer(I8P), intent(inout) :: array(:) !< Array to be sorted.
      integer(I8P)                :: half     !< Half size counter.

      half = (size(array) + 1) / 2
      if (size(array) < 2) then
         continue
      else if (size(array) == 2) then
         if (self%greater(array(1), array(2))) call swap_element(array(1), array(2))
      else
         call mergesort(array( : half))
         call mergesort(array(half + 1 :))
         if (self%greater(array(half), array(half + 1))) then
            work(1 : half) = array(1 : half)
            call merge_array(work(1 : half), array(half + 1:), array)
         endif
      endif
      endsubroutine mergesort

      subroutine merge_array(A, B, C)
      !< Merge arrays A/B in C.
      integer(I8P), target, intent(in)    :: A(:), B(:) !< Input arrays.
      integer(I8P), target, intent(inout) :: C(:)       !< Output array.
      integer(I8P)                        :: i, j, k    !< Counter.

      if (size(A) + size(B) > size(C)) stop

      i = 1 ; j = 1
      do k = 1, size(C)
         if (i <= size(A) .and. j <= size(B)) then
            if (self%lower(A(i), B(j))) then
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

   subroutine destroy(self)
   !< Destroy the tree.
   class(tree_object), intent(inout) :: self  !< The tree.
   type(tree_object)                 :: fresh !< Fresh tree.

   self = fresh
   endsubroutine destroy

   function loop(self, code, node) result(again)
   !< Sentinel while-loop on nodes returning the code (for tree looping).
   class(tree_object),     intent(in)                     :: self      !< The tree bucket.
   integer(I8P),           intent(out), optional          :: code      !< The Morton code.
   type(tree_node_object), intent(out), optional, pointer :: node      !< Pointer to current node.
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
                  if (present(node)) node => p
                  again = .true.
                  return
               elseif (associated(p%next)) then
                  p => p%next
                  if (present(code)) code = p%code
                  if (present(node)) node => p
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

   subroutine initialize(self, grid, max_load, nodes_number, buckets_number, ratio, max_level, add_adam)
   !< Initialize the tree.
   class(tree_object), intent(inout)        :: self           !< The tree.
   type(grid_object),  intent(in), target   :: grid           !< Grid data.
   real(R8P),          intent(in), optional :: max_load       !< Maximum load of tree buckets.
   integer(I8P),       intent(in), optional :: nodes_number   !< Nodes number to be stored in the tree.
   integer(I8P),       intent(in), optional :: buckets_number !< Number of buckets for initialize the tree.
   integer(I4P),       intent(in), optional :: ratio          !< Refinement ratio.
   integer(I4P),       intent(in), optional :: max_level      !< Maximum refinement level.
   logical,            intent(in), optional :: add_adam       !< Add ADAM node, the ancestor of all nodes.
   logical                                  :: add_adam_      !< Add ADAM node, the ancestor of all nodes, local var.
#ifdef _MPI_
   integer(I4P)                             :: error          !< Error traping flag.
#endif

   call self%destroy
   self%grid => grid
   add_adam_ = .true. ; if (present(add_adam)) add_adam_ = add_adam
   ! tree data
   if (present(max_load)) self%max_load = max_load
   if (present(nodes_number)) then
      self%buckets_number = self%prime_buckets_number(nodes_number=nodes_number)
   else
      self%buckets_number = TREE_BUCKETS_NUMBER_DEF ; if (present(buckets_number)) self%buckets_number = buckets_number
   endif
   allocate(self%bucket(1:self%buckets_number))
   if (present(ratio)) then
      if (ratio==8_I8P.or.ratio==4_I8P) then
         self%ratio = ratio
      else
         write(stderr, '(A)') 'ADAM-ERROR: tree ratio must be 8 o 4'
#ifdef _MPI_
   call MPI_FINALIZE(error)
#endif
        stop
      endif
   endif
   if (present(max_level)) self%max_level = max_level
   self%is_initialized_ = .true.
   if (add_adam_) call self%add_node(code=-1_I8P) ! add ADAM node, the ancestor of all nodes
   ! MPI data
#ifdef _MPI_
   call MPI_COMM_SIZE(MPI_COMM_WORLD, self%procs_number, error)
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, error)
#endif
   allocate(self%comm_map_n_send(0:self%procs_number-1))
   allocate(self%comm_map_n_recv(0:self%procs_number-1))
   allocate(self%comm_map_send_ptr(0:self%procs_number))
   allocate(self%comm_map_recv_ptr(0:self%procs_number))
   self%comm_map_n_send = 0_I4P
   self%comm_map_n_recv = 0_I4P
   self%comm_map_send_ptr = 0_I4P
   self%comm_map_recv_ptr = 0_I4P
   call self%mpi_redistribute
   endsubroutine initialize

   subroutine mark_all_nodes(self, mark)
   !< Mark all nodes to be refined.
   class(tree_object), intent(inout) :: self !< The tree.
   integer(I4P),       intent(in)    :: mark !< Mark to be imposed [TO_BE_REFINED,...]
   type(tree_node_object), pointer   :: node !< Pointer to current node.

   do while(self%loop(node=node))
      node%refinement_needed = mark
   enddo
   endsubroutine mark_all_nodes

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
   type(tree_node_object), pointer          :: node         !< Pointer to node.

   if (self%is_initialized_) then
      if (present(max_load)) self%max_load = max_load
      if (self%nodes_number > int((1._R8P/self%max_load)*nodes_number, I4P)) return ! new size too small, cannot previous nodes
      call swap%initialize(grid=self%grid, max_load=self%max_load, nodes_number=nodes_number, ratio=self%ratio, add_adam=.false.)
      do while(self%loop(node=node)) ! re-hash all codes
         call swap%add_node(code=node%code,                           &
                            refinement_needed=node%refinement_needed, &
                            myrank=node%myrank,                       &
                            block_index=node%block_index)
      enddo
      call move_alloc(from=swap%bucket, to=self%bucket)
      self%buckets_number = swap%buckets_number
      self%nodes_number   = swap%nodes_number
   else
      print '(A)', 'ERROR: tree is not initialized, cannot be resized'
      stop
   endif
   endsubroutine resize

   subroutine sanitize(self, iterations_number)
   !< Sanitize the tree.
   class(tree_object),        intent(inout)        :: self                 !< The tree.
   integer(I4P),              intent(in), optional :: iterations_number    !< Sanitazie iterations number.
   integer(I4P)                                    :: iterations_number_   !< Sanitazie iterations number.
   type(tree_node_object), pointer                 :: node                 !< Pointer to node.
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

   iterations_number_ = TREE_MAX_SANITIZE_ITERATIONS ; if (present(iterations_number)) iterations_number_ = iterations_number

   min_max_check_loop : do while(self%loop(node=node))
      new_level = self%level(code=node%code) + node%refinement_needed
      if ((new_level > self%max_level).or.(new_level < 0)) then
         node%refinement_needed = TO_NOT_TOUCH
      endif
   enddo min_max_check_loop

   sanitize_loop : do s=1, iterations_number_
      is_sanitize_complete = .true.

      ! check for the sanity of derefinement
      if (allocated(self%node_to_derefine)) deallocate(self%node_to_derefine) ; allocate(self%node_to_derefine(0))
      if (allocated(codes_analyzed))   deallocate(codes_analyzed)   ; allocate(codes_analyzed(0))
      derefine_loop : do while(self%loop(node=node))
         ! check if I want to be derefined and I have not been analyzed yet
         if (node%refinement_needed == TO_BE_DEREFINED) then
            if (findloc(codes_analyzed, node%code, dim=1)==0) then ! avoid to re-analyze already confirmed siblingsi to derefine
               ! check sibling for derefinement
               can_be_derefined = .true.
               code = node%code
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
               else
                  is_sanitize_complete = .false.
                  node%refinement_needed = TO_NOT_TOUCH
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

     ! check for the sanity of refinement (2:1 rule)
     refine_loop : do while(self%loop(node=node))
        new_level = self%level(code=node%code) + node%refinement_needed
        face_loop : do f=1, 6
           call self%get_neighbor(code=node%code, face=f, neighbor=neighbor, neighbor_type=neighbor_type)
           if (neighbor_type /= NODE_BOUNDARY_CONDITION) then
              neighbor_loop : do n=1, size(neighbor, dim=1)
                 ! check level
                 neigh => self%node(code=neighbor(n))
                 new_level_n = self%level(code=neighbor(n)) + neigh%refinement_needed
                 if (new_level_n > new_level + 1) then
                    ! a neighbour want to be refined 2 levels more than me, I have to refine more too
                    is_sanitize_complete = .false.
                    if     (new_level_n - new_level == 3) then ! node want to derefine, but it must be refined
                       node%refinement_needed = 1
                    elseif (new_level_n - new_level == 2) then
                       node%refinement_needed = node%refinement_needed + 1
                    else
                       print '(A)',  'SOMETHING WENT TERRIBLY WRONG. EXIT!'
                       print '(A)',  'REFINEMENT NEEDED '//trim(str(node%refinement_needed,.true.))
                       print '(A)',  'SANITIZE ITERATIONS '//trim(str(s,.true.))
                       stop
                    endif
                    new_level = self%level(code=node%code) + node%refinement_needed
                 endif
              enddo neighbor_loop
           endif
        enddo face_loop

        if (node%refinement_needed > 1) then
           print '(A)',  'CANNOT REFINE TWICE IN A ROW. SOMETHING WENT TERRIBLY WRONG. EXIT!'
           print '(A)',  'SANITIZE ITERATIONS '//trim(str(s,.true.))
           stop
        endif

        new_level = self%level(code=node%code) + node%refinement_needed
        if (new_level > self%max_level) then
           print '(A)',  'CANNOT REFINE MORE. SOMETHING WENT TERRIBLY WRONG. EXIT!'
           stop
        endif
     enddo refine_loop
     if (is_sanitize_complete) exit sanitize_loop
   enddo sanitize_loop

   if (.not.is_sanitize_complete) then
      print '(A)',  'SANITZE CANNOT BE COMPLETED. SOMETHING WENT TERRIBLY WRONG. EXIT!'
      stop
   endif

   ! update to_refine list
   if (allocated(self%node_to_refine)) deallocate(self%node_to_refine) ; allocate(self%node_to_refine(0))
   do while(self%loop(node=node))
      if (node%refinement_needed==TO_BE_REFINED) self%node_to_refine = [self%node_to_refine, [node%code]]
   enddo
   endsubroutine sanitize

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
   subroutine import_refinements_needed(self, refinements_needed_all, disp_count)
   !< Import refinements needed status changed externally.
   class(tree_object),        intent(inout) :: self                      !< The tree.
   integer(I4P), allocatable, intent(in)    :: refinements_needed_all(:) !< Refinements needed of all blocks.
   integer(I4P), allocatable, intent(in)    :: disp_count(:)             !< Displacement of blocks that are received from process.
   type(tree_node_object), pointer          :: node                      !< Pointer to current node.
   integer(I8P)                             :: b                         !< Counter.
   integer(I4P)                             :: myrank                    !< Counter.

   do while(self%loop(node=node))
      myrank = node%myrank
      b = node%block_index
      node%refinement_needed = refinements_needed_all(disp_count(myrank)+b)
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
   type(tree_node_object), pointer   :: node                 !< Pointer to current node.
   integer(I8P), allocatable         :: codes(:)             !< List of (sorted) codes.
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
   do while(self%loop(node=node))
      if (node%myrank==self%myrank) my_nodes_number = my_nodes_number + 1_I8P
      if     (is_node_to_send(n=node)) then
         ! I have this node that must be sent to node%myrank_new
         self%comm_map_n_send(node%myrank_new) = self%comm_map_n_send(node%myrank_new) + 1
      elseif (is_node_to_receive(n=node)) then
         ! node%rank has this node that must be sent to me
         self%comm_map_n_recv(node%myrank) = self%comm_map_n_recv(node%myrank) + 1
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
   do p=1, self%procs_number
      self%comm_map_send_ptr(p) = self%comm_map_send_ptr(p-1) + self%comm_map_n_send(p-1)
      self%comm_map_recv_ptr(p) = self%comm_map_recv_ptr(p-1) + self%comm_map_n_recv(p-1)
   enddo
   comm_map_send_ctr = self%comm_map_send_ptr
   comm_map_recv_ctr = self%comm_map_recv_ptr

   ! populate communication maps
   codes = self%codes() ! sorted list of Morton codes
   ! send map
   if (n_send>0_I4P) then
      do c=1, size(codes, dim=1) ! sort blocks in Morton order
         node => self%node(code=codes(c))
         if (is_node_to_send(n=node)) then
            self%comm_map_send(comm_map_send_ctr(node%myrank_new)+1) = node%block_index
            comm_map_send_ctr(node%myrank_new) = comm_map_send_ctr(node%myrank_new) + 1
         endif
      enddo
   endif
   ! receive map
   if (n_recv>0_I4P) then
      do c=1, size(codes, dim=1) ! sort blocks in Morton order
         node => self%node(code=codes(c))
         if (is_node_to_receive(n=node)) then
            self%comm_map_recv(comm_map_recv_ctr(node%myrank)+1) = node%block_index_new
            comm_map_recv_ctr(node%myrank) = comm_map_recv_ctr(node%myrank) + 1
         endif
      enddo
   endif
   ! local maps
   if (n_keep > 0_I4P) then
      l = 0
      do c=1, size(codes, dim=1) ! sort blocks in Morton order
         node => self%node(code=codes(c))
         if (is_node_to_keep(n=node)) then
            l = l + 1
            self%local_map(l,1) = node%block_index_new
            self%local_map(l,2) = node%block_index
         endif
      enddo
   endif
   contains
      function is_node_to_keep(n)
      !< Check if node `n` must be kept.
      type(tree_node_object), intent(in), pointer :: n               !< Pointer to current node.
      logical                                     :: is_node_to_keep !< Check result.

      is_node_to_keep = ((self%myrank == n%myrank).and.(n%myrank == n%myrank_new))
      endfunction is_node_to_keep

      function is_node_to_send(n)
      !< Check if node `n` must be sent.
      type(tree_node_object), intent(in), pointer :: n               !< Pointer to current node.
      logical                                     :: is_node_to_send !< Check result.

      is_node_to_send = ((self%myrank == n%myrank).and.(n%myrank /= n%myrank_new))
      endfunction is_node_to_send

      function is_node_to_receive(n)
      !< Check if node `n` must be received.
      type(tree_node_object), intent(in), pointer :: n                  !< Pointer to current node.
      logical                                     :: is_node_to_receive !< Check result.

      is_node_to_receive = ((self%myrank == n%myrank_new).and.(n%myrank /= n%myrank_new))
      endfunction is_node_to_receive
   endsubroutine make_comm_local_maps

   subroutine mpi_gather_refinements_needed(self)
   !< Gather refinements needed status between MPI processes.
   class(tree_object), intent(inout) :: self           !< The tree.
   type(tree_node_object), pointer   :: node           !< Pointer to current node.
   integer(I8P), allocatable         :: send_buffer(:) !< Send buffer nodes data.
   integer(I8P), allocatable         :: recv_buffer(:) !< Recv buffer nodes data.
   integer(I4P), allocatable         :: recv_count(:)  !< Number of nodes that are received from each process.
   integer(I4P), allocatable         :: disp_count(:)  !< Displacement of nodes that are received from each process.
   integer(I8P)                      :: n, p           !< Counter.
#ifdef _MPI_
   integer(I4P)                      :: error          !< Error traping flag.
#endif

   allocate(send_buffer(self%my_nodes_number * 2)) ! [Morton code, refinement_needed]
   allocate(recv_buffer(self%nodes_number    * 2)) ! [Morton code, refinement_needed]
   allocate(recv_count(0:self%procs_number - 1))
   allocate(disp_count(0:self%procs_number - 1))
   ! populating receive counters and send buffer
   recv_count = 0_I4P
   n = 0_I8P
   do while(self%loop(node=node))
      recv_count(node%myrank) = recv_count(node%myrank) + 2
      if (self%myrank == node%myrank) then
         n = n + 1 ; send_buffer(n) = node%code
         n = n + 1 ; send_buffer(n) = node%refinement_needed - 100_I8P ! shift to negative for resolving conflicts with Morton code
      endif
   enddo
   ! computing displacement counts
   disp_count = 0_I4P
   do p=1, self%procs_number - 1
      disp_count(p) = disp_count(p-1) + recv_count(p-1)
   enddo

#ifdef _MPI_
   call MPI_ALLGATHERV(send_buffer, self%my_nodes_number * 2, MPI_INTEGER8, &
                       recv_buffer, recv_count, disp_count, MPI_INTEGER8, MPI_COMM_WORLD, error)
#endif

   ! update nodes data
   do while(self%loop(node=node))
      if (self%myrank == node%myrank) cycle ! my nodes are already updated
      n = findloc(recv_buffer, node%code, dim=1) + 1
      node%refinement_needed = int(recv_buffer(n) + 100_I8P, I4P)
   enddo
   endsubroutine mpi_gather_refinements_needed

   subroutine mpi_print_stats(self)
   !< Print MPI stats.
   class(tree_object), intent(in) :: self !< The tree.
   integer(I4P)                   :: p    !< Counter.

   do p=0, self%procs_number-1
      print '(A)', ' myrank: '//trim(str(self%myrank,.true.))//&
                   ' send to: '//trim(str(p,.true.))//' blocks n.: '//trim(str(self%comm_map_n_send(p),.true.))
   enddo
   if (allocated(self%comm_map_send)) &
      print '(A)', ' myrank: '//trim(str(self%myrank,.true.))//&
                   ' blocks sent: '//trim(str(self%comm_map_send,.true.))
   do p=0, self%procs_number-1
      print '(A)', ' myrank: '//trim(str(self%myrank,.true.))//&
                   ' recv from: '//trim(str(p,.true.))//' blocks n.: '//trim(str(self%comm_map_n_recv(p),.true.))
   enddo
   if (allocated(self%comm_map_recv)) &
      print '(A)', ' myrank: '//trim(str(self%myrank,.true.))//&
                   ' blocks recv:  '//trim(str(self%comm_map_recv,.true.))
   if (allocated(self%local_map)) &
      print '(A)', ' myrank: '//trim(str(self%myrank,.true.))//&
                   ' blocks kept n.: '//trim(str(size(self%local_map(:,1),dim=1),.true.))
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
   integer(I8P), allocatable         :: codes(:)         !< List of (sorted) codes.
   type(tree_node_object), pointer   :: node             !< Pointer to current node.
   integer(I4P)                      :: p                !< Processes counter.
   integer(I4P)                      :: cl               !< Local child counter.
   integer(I8P)                      :: c                !< Codes counter.
   integer(I4P)                      :: i, j, k, l       !< Coordinates.
   integer(I8P)                      :: block_index_new  !< New block index counter.
   integer(I8P)                      :: my_codes_number  !< Number of codes for each process for a balanced workload.
   integer(I4P)                      :: child_local      !< Local numbering.
   integer(I8P)                      :: n_keep           !< Number of keept nodes.
   integer(I8P)                      :: n_recv           !< Number of nodes that I have to receive.

   codes = self%codes() ! sorted list of codes
   my_codes_number = nint(real(size(codes, dim=1),R8P) / self%procs_number)
   ! initialize process rank and my codes number
   p = 0_I4P
   block_index_new = 0_I8P
   ! loop over all codes
   c = 1
   do while(c<=size(codes, dim=1))
      if (block_index_new > my_codes_number.and.p < self%procs_number-1) then ! I would like to split...
         if (can_split()) then
            ! I am lucky, the split does not separate siblings
            p = p + 1_I4P
            block_index_new = 1_I8P
            node => self%node(code=codes(c))
            node%myrank_new = p
            node%block_index_new = block_index_new
         else
            ! I am not lucky, the split would separate siblings
            child_local = self%child_local(code=codes(c))
            if (child_local > self%ratio/2 -1) then
               ! the split falls in the second half of siblings list, place all nodes in the current process
               do cl=child_local, self%ratio-1
                  block_index_new = block_index_new + 1
                  node => self%node(code=codes(c+cl-child_local))
                  node%myrank_new = p
                  node%block_index_new = block_index_new
               enddo
            else
               ! the split falls in the first half of siblings list, place all nodes in the next process
               p = p + 1_I4P
               block_index_new = 0_I8P
               do cl=0, self%ratio-1
                  block_index_new = block_index_new + 1
                  node => self%node(code=codes(c+cl-child_local))
                  node%myrank_new = p
                  node%block_index_new = block_index_new
               enddo
            endif
            ! update codes counter skipping all current siblings
            c = c + self%ratio - 1 - child_local
         endif
      else ! no split, keeping to place nodes in the current process
         block_index_new = block_index_new + 1
         node => self%node(code=codes(c))
         node%myrank_new = p
         node%block_index_new = block_index_new
      endif
      c = c + 1
   enddo
   ! create communication/local maps
   call self%make_comm_local_maps
   ! update tree status and compute coordinates of redistributed nodes
   n_keep = 0_I8P ; if (allocated(self%local_map    )) n_keep = size(self%local_map,     dim=1)
   n_recv = 0_I8P ; if (allocated(self%comm_map_recv)) n_recv = size(self%comm_map_recv, dim=1)
   if (allocated(self%block_coordinates)) deallocate(self%block_coordinates) ; allocate(self%block_coordinates(n_keep + n_recv, 4))
   do while(self%loop(node=node))
      node%myrank = node%myrank_new
      node%block_index = node%block_index_new
      if (node%myrank == self%myrank) then
         select case(self%ratio)
         case(4_I4P)
            call self%morton_to_coordinates(code=node%code, i=i, j=j, l=l)
         case(8_I4P)
            call self%morton_to_coordinates(code=node%code, i=i, j=j, k=k, l=l)
         endselect
         self%block_coordinates(node%block_index, 1) = i
         self%block_coordinates(node%block_index, 2) = j
         self%block_coordinates(node%block_index, 3) = k
         self%block_coordinates(node%block_index, 4) = l
      endif
   enddo
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
      siblings = self%siblings(code=codes(c))
      if (all([(self%has_code(code=siblings(s)), s=1,self%ratio-1)])) then ! if all siblings exist
         can_split = .not.(findloc(siblings, codes(c-1),dim=1) > 0) ! if my predecessor is a sibling the split is not allowed
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

   select case(self%ratio)
   case(4_I4P)
      children = [self%child(code=code, i=0), &
                  self%child(code=code, i=1), &
                  self%child(code=code, i=2), &
                  self%child(code=code, i=3)]
   case(8_I4P)
      children = [self%child(code=code, i=0), &
                  self%child(code=code, i=1), &
                  self%child(code=code, i=2), &
                  self%child(code=code, i=3), &
                  self%child(code=code, i=4), &
                  self%child(code=code, i=5), &
                  self%child(code=code, i=6), &
                  self%child(code=code, i=7)]
   endselect
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

   code = 0
   if (level>1) then
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

   if (code==-1) then
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
   select case(self%ratio)
   case(4_I4P)
      call self%morton_to_coordinates(code=code, i=i_dn(1), j=j_dn(1), l=l_dn)
   case(8_I4P)
      call self%morton_to_coordinates(code=code, i=i_dn(1), j=j_dn(1), k=k_dn(1), l=l_dn)
   endselect

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
      neighbor_type = NODE_STANDARD
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

   elemental function greater(self, lhs, rhs) result(res)
   !< Return true if code is greater than other.
   class(tree_object), intent(in) :: self !< The tree.
   integer(I8P),       intent(in) :: lhs  !< Left hand side of code comparison.
   integer(I8P),       intent(in) :: rhs  !< Right hand side of code comparison.
   logical                        :: res  !< Comparison result.

   res = self%finest_at_level(code=lhs, level=self%max_level) > self%finest_at_level(code=rhs, level=self%max_level)
   endfunction greater

   elemental function parent(self, code)
   !< Return the parent given Morton code.
   class(tree_object), intent(in) :: self   !< The tree.
   integer(I8P),       intent(in) :: code   !< Morton code.
   integer(I8P)                   :: parent !< Parent Morton code.

   parent = -1 ! ancestor of all has not parent
   if (code>self%ratio-1) parent = (code - self%ratio) / self%ratio
   endfunction parent

   elemental function parent_at_level(self, code, level) result(parent)
   !< Return the parent given Morton code at a given level.
   class(tree_object), intent(in) :: self     !< The tree.
   integer(I8P),       intent(in) :: code     !< Morton code.
   integer(I4P),       intent(in) :: level    !< Refinement level.
   integer(I8P)                   :: parent   !< Parent Morton code.
   integer(I8P), allocatable      :: path_(:) !< Parent Morton code.

   parent = -1_I8P ! ancestor of all nodes
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
      call move_alloc(from=temp,to=path)
      c = self%parent(code=c)
   enddo
   endfunction path

   subroutine print_code_topology(self, code,  &
                                  coordinates, &
                                  level,       &
                                  parent,      &
                                  parents,     &
                                  path,        &
                                  child,       &
                                  child_local, &
                                  finest,      &
                                  siblings,    &
                                  neighbor,    &
                                  block_index, &
                                  whole)
   !< Print all code topology data.
   class(tree_object), intent(in)           :: self          !< The tree.
   integer(I8P),       intent(in)           :: code          !< Morton code.
   logical,            intent(in), optional :: coordinates   !< Coordinates.
   logical,            intent(in), optional :: level         !< Level of node.
   logical,            intent(in), optional :: parent        !< Parent of code.
   logical,            intent(in), optional :: parents       !< Parents list.
   logical,            intent(in), optional :: path          !< Path from node to parent of first level.
   logical,            intent(in), optional :: child         !< (First) Child of code.
   logical,            intent(in), optional :: child_local   !< Local child-numbering of code.
   logical,            intent(in), optional :: finest        !< Finest Morton code.
   logical,            intent(in), optional :: siblings      !< Siblings of code.
   logical,            intent(in), optional :: neighbor      !< Neighbor of code.
   logical,            intent(in), optional :: block_index   !< Block index in the field array.
   logical,            intent(in), optional :: whole         !< Whole topology data.
   logical                                  :: coordinates_  !< Coordinates, local var.
   logical                                  :: level_        !< Level of node, local var.
   logical                                  :: parent_       !< Parent of code, local var.
   logical                                  :: child_        !< (First) Child of code, local var.
   logical                                  :: child_local_  !< Local child-numbering of code, local var.
   logical                                  :: finest_       !< Finest Morton code, local var.
   logical                                  :: siblings_     !< Siblings of code, local var.
   logical                                  :: path_         !< Path from node to parent of first level, local var.
   logical                                  :: parents_      !< Parents list, local var.
   logical                                  :: neighbor_     !< Neighbor of code, local var.
   logical                                  :: block_index_  !< Block index in the field array.
   logical                                  :: whole_        !< Whole topology data, local var.
   character(:), allocatable                :: topology      !< Topology string.
   character(:), allocatable                :: parents_str   !< Parents string list.
   character(:), allocatable                :: neighbors_str !< Neighbors string.
   integer(I8P), allocatable                :: neighbors(:)  !< Neighbor of code.
   integer(I4P)                             :: neighbor_type !< Type of neighbor.
   type(tree_node_object), pointer          :: node          !< Node pointer.
   integer(I4P)                             :: f             !< Counter.
   integer(I4P)                             :: i, j, k, l    !< Counter.

   coordinates_ = .false. ; if (present(coordinates))  coordinates_ = coordinates
   level_       = .false. ; if (present(level      ))  level_       = level
   parent_      = .false. ; if (present(parent     ))  parent_      = parent
   child_       = .false. ; if (present(child      ))  child_       = child
   child_local_ = .false. ; if (present(child_local))  child_local_ = child_local
   finest_      = .false. ; if (present(finest     ))  finest_      = finest
   siblings_    = .false. ; if (present(siblings   ))  siblings_    = siblings
   path_        = .false. ; if (present(path       ))  path_        = path
   parents_     = .false. ; if (present(parents    ))  parents_     = parents
   neighbor_    = .false. ; if (present(neighbor   ))  neighbor_    = neighbor
   block_index_ = .false. ; if (present(block_index))  block_index_ = block_index
   whole_       = .false. ; if (present(whole      ))  whole_       = whole

   node => self%node(code=code)

   parents_str = ''
   do l=self%level(code=code)-1, 0, -1
      parents_str = parents_str//' pal_'//trim(str(l,.true.))//' '//trim(str(self%parent_at_level(code=code, level=l),.true.))
   enddo

   neighbors_str = ''
   do f=1, 6
      if (self%ratio==4_I4P.and.f>4) exit
      call self%get_neighbor(code=code, neighbor=neighbors, face=f, neighbor_type=neighbor_type)
      if (allocated(neighbors)) then
         neighbors_str = neighbors_str//' f_'//trim(str(f,.true.))//' '//trim(str(neighbors,.true.))//&
                         ' type-'//trim(str(neighbor_type,.true.))
      else
         neighbors_str = neighbors_str//' f-'//trim(str(f,.true.))//' type_'//trim(str(neighbor_type,.true.))
      endif
   enddo

   select case(self%ratio)
   case(4_I4P)
      call self%morton_to_coordinates(code=node%code, i=i, j=j, l=l)
   case(8_I4P)
      call self%morton_to_coordinates(code=node%code, i=i, j=j, k=k, l=l)
   endselect
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
   if (block_index_.or.whole_) topology = topology//' block_index: '//trim(str(node%block_index,.true.))

   print '(A)', topology
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
   subroutine add_node(self, code, refinement_needed, myrank, block_index, update_last_block_index)
   !< Add a node pointer to the tree.
   !<
   !< @note If a node with the same key is already in the tree, it is removed and the new one will replace it.
   class(tree_object), intent(inout)        :: self                     !< The tree.
   integer(I8P),       intent(in)           :: code                     !< The Morton code.
   integer(I4P),       intent(in), optional :: refinement_needed        !< Flag for refinement/derefinement algorithm.
   integer(I4P),       intent(in), optional :: myrank                   !< MPI rank process.
   integer(I8P),       intent(in), optional :: block_index              !< Block index in the field array.
   logical,            intent(in), optional :: update_last_block_index  !< Update or not last block index.
   logical                                  :: update_last_block_index_ !< Update or not last block index, local var.
   integer(I4P)                             :: b                        !< Bucket index, namely hashed key.

   if (.not.self%is_initialized_) then
      print '(A)', 'ERROR: cannot add a node a non initialized tree'
   endif
   ! if the code is not already in the tree update the nodes number otherwise not
   if (.not.self%has_code(code=code)) self%nodes_number = self%nodes_number + 1
   b = self%hash(code=code)
   call self%bucket(b)%add_node(code=code, refinement_needed=refinement_needed, &
                                myrank=myrank, block_index=block_index)
   update_last_block_index_ = .true. ; if (present(update_last_block_index)) update_last_block_index_ = update_last_block_index
   if (update_last_block_index_) self%last_block_index = self%last_block_index + 1
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
   type(tree_node_object), pointer   :: node             !< Pointer to node.
   integer(I8P)                      :: derefined_number !< Number of derefined blocks.
   integer(I8P)                      :: n                !< Counter.
   integer(I4P)                      :: i                !< Counter.

   if (allocated(self%block_to_derefine)) deallocate(self%block_to_derefine)
   if (allocated(self%block_derefined)) deallocate(self%block_derefined)
   derefined_number = size(self%node_to_derefine, dim=1)
   if (derefined_number>0) then
      allocate(self%block_to_derefine(derefined_number))
      allocate(self%block_derefined(2, derefined_number/self%ratio))
      if (allocated(self%node_to_derefine)) then
         do n=1, size(self%node_to_derefine, dim=1), self%ratio
            first_child => self%node(code=self%node_to_derefine(n))
            self%block_derefined(1,(n-1)/self%ratio+1) = self%parent(code=first_child%code)
            self%block_derefined(2,(n-1)/self%ratio+1) = first_child%block_index
            call self%add_node(code=self%parent(code=first_child%code),              &
                                                myrank=first_child%myrank,           &
                                                block_index=first_child%block_index, &
                                                update_last_block_index=.false.)
            do i=0, self%ratio - 1
               node => self%node(code=self%node_to_derefine(n+i))
               self%block_to_derefine(n+i) = node%block_index
               call self%remove_node(code=self%node_to_derefine(n+i))
            enddo
         enddo
      endif
   endif
   endsubroutine derefine

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
   integer(I8P)                      :: n              !< Counter.
   integer(I4P)                      :: i              !< Counter.

   if (allocated(self%block_to_refine)) deallocate(self%block_to_refine)
   if (allocated(self%block_refined)) deallocate(self%block_refined)
   refined_number = size(self%node_to_refine, dim=1)
   if (refined_number>0) then
      allocate(self%block_to_refine(2, refined_number))
      allocate(self%block_refined(2, self%ratio*refined_number))
      do n=1, refined_number
         parent => self%node(code=self%node_to_refine(n))
         self%block_to_refine(1,n) = parent%block_index
         self%block_to_refine(2,n) = parent%myrank
         call self%add_node(code=self%child(code=parent%code, i=0), myrank=parent%myrank, &
                            block_index=parent%block_index, update_last_block_index=.false.)
         self%block_refined(1, (n-1)*self%ratio+1) = self%child(code=parent%code, i=0)
         self%block_refined(2, (n-1)*self%ratio+1) = parent%block_index
         do i=1, self%ratio-1
            self%block_refined(1, (n-1)*self%ratio+1+i) = self%child(code=parent%code, i=i)
            self%block_refined(2, (n-1)*self%ratio+1+i) = self%last_block_index + 1
            call self%add_node(code=self%child(code=parent%code, i=i), myrank=parent%myrank, &
                               block_index=self%last_block_index+1)
         enddo
         call self%remove_node(code=parent%code)
      enddo
   endif
   endsubroutine refine

   ! operators
   ! =
   subroutine tree_assign_tree(lhs, rhs)
   !< Operator `=`.
   class(tree_object), intent(inout) :: lhs !< Left hand side.
   type(tree_object),  intent(in)    :: rhs !< Right hand side.
   integer(I4P)                      :: b   !< Counter.

   ! tree data
   lhs%grid => rhs%grid
   if (allocated(rhs%bucket)) then
      lhs%bucket = rhs%bucket
   else
      if (allocated(lhs%bucket)) then
         do b=lbound(lhs%bucket, dim=1), ubound(lhs%bucket, dim=1)
            call lhs%bucket(b)%destroy
         enddo
         deallocate(lhs%bucket)
      endif
   endif
   lhs%buckets_number = rhs%buckets_number
   lhs%nodes_number = rhs%nodes_number
   lhs%max_load = rhs%max_load
   lhs%ratio = rhs%ratio
   lhs%max_level = rhs%max_level
   lhs%is_initialized_ = rhs%is_initialized_
   ! AMR data
   lhs%last_block_index = rhs%last_block_index
   if (allocated(rhs%node_to_refine)) then
      lhs%node_to_refine = rhs%node_to_refine
   else
      if (allocated(lhs%node_to_refine)) deallocate(lhs%node_to_refine)
   endif
   if (allocated(rhs%node_to_derefine)) then
      lhs%node_to_derefine = rhs%node_to_derefine
   else
      if (allocated(lhs%node_to_derefine)) deallocate(lhs%node_to_derefine)
   endif
   if (allocated(rhs%block_to_refine)) then
      lhs%block_to_refine = rhs%block_to_refine
   else
      if (allocated(lhs%block_to_refine)) deallocate(lhs%block_to_refine)
   endif
   if (allocated(rhs%block_refined)) then
      lhs%block_refined = rhs%block_refined
   else
      if (allocated(lhs%block_refined)) deallocate(lhs%block_refined)
   endif
   if (allocated(rhs%block_to_derefine)) then
      lhs%block_to_derefine = rhs%block_to_derefine
   else
      if (allocated(lhs%block_to_derefine)) deallocate(lhs%block_to_derefine)
   endif
   if (allocated(rhs%block_derefined)) then
      lhs%block_derefined = rhs%block_derefined
   else
      if (allocated(lhs%block_derefined)) deallocate(lhs%block_derefined)
   endif
   if (allocated(rhs%block_coordinates)) then
      lhs%block_coordinates = rhs%block_coordinates
   else
      if (allocated(lhs%block_coordinates)) deallocate(lhs%block_coordinates)
   endif
   ! MPI data
   lhs%procs_number = rhs%procs_number
   lhs%myrank = rhs%myrank
   lhs%my_nodes_number = rhs%my_nodes_number
   lhs%send_nodes_number = rhs%send_nodes_number
   lhs%recv_nodes_number = rhs%recv_nodes_number
   lhs%keep_nodes_number = rhs%keep_nodes_number
   if (allocated(rhs%comm_map_n_send)) then
      lhs%comm_map_n_send = rhs%comm_map_n_send
   else
      if (allocated(lhs%comm_map_n_send)) deallocate(lhs%comm_map_n_send)
   endif
   if (allocated(rhs%comm_map_n_recv)) then
      lhs%comm_map_n_recv = rhs%comm_map_n_recv
   else
      if (allocated(lhs%comm_map_n_recv)) deallocate(lhs%comm_map_n_recv)
   endif
   if (allocated(rhs%comm_map_send_ptr)) then
      lhs%comm_map_send_ptr = rhs%comm_map_send_ptr
   else
      if (allocated(lhs%comm_map_send_ptr)) deallocate(lhs%comm_map_send_ptr)
   endif
   if (allocated(rhs%comm_map_recv_ptr)) then
      lhs%comm_map_recv_ptr = rhs%comm_map_recv_ptr
   else
      if (allocated(lhs%comm_map_recv_ptr)) deallocate(lhs%comm_map_recv_ptr)
   endif
   if (allocated(rhs%comm_map_send)) then
      lhs%comm_map_send = rhs%comm_map_send
   else
      if (allocated(lhs%comm_map_send)) deallocate(lhs%comm_map_send)
   endif
   if (allocated(rhs%comm_map_recv)) then
      lhs%comm_map_recv = rhs%comm_map_recv
   else
      if (allocated(lhs%comm_map_recv)) deallocate(lhs%comm_map_recv)
   endif
   if (allocated(rhs%local_map)) then
      lhs%local_map = rhs%local_map
   else
      if (allocated(lhs%local_map)) deallocate(lhs%local_map)
   endif
   endsubroutine tree_assign_tree
endmodule adam_tree_object
