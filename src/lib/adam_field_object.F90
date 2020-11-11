!< ADAM, field class definition.
module adam_field_object
!< ADAM, field class definition.

!< A structured block is composed of hexahedron finite volumes with quadrilateral faces using the
!< following internal numeration for nodes and faces:
!<```
!< /|\Z
!<  |                            F(4)         _ F(6)
!<  |                            /|\          /!
!<  |                        7    |          /    8
!<  |                         *------------------*
!<  |                        /|   |        /    /|
!<  |                       / |   |       /    / |
!<  |                      /  |   |      /    /  |
!<  |                     /   |   |     /    /   |
!<  |                    /    |   |    +    /    |
!<  |                   /     |   |        /     |
!<  |                  /      |   +       /      |
!<  |                 /      3|          /       |4
!<  |                /        * --------/--------*
!<  |      F(1)<----/----+   /         /        /
!<  |              *------------------*    +-------->F(2)
!<  |             5|       /          |6      /
!<  |              |      /           |      /
!<  |              |     /        +   |     /
!<  |              |    /         |   |    /
!<  |              |   /      +   |   |   /
!<  |              |  /      /    |   |  /
!<  |              | /      /     |   | /
!<  |              |/      /      |   |/
!<  |              *------------------*
!<  |             1      /        |    2
!<  |                   /        \|/
!<  |   _ Y           |/_       F(3)
!<  |   /|         F(5)
!<  |  /
!<  | /
!<  |/                                                    X
!<  O----------------------------------------------------->
!<```
!< Each hexadron cells is faces-connected to its neighboring, thus the cells build a structured block with implicit
!< connectivity, e.g. in 2D space a FriVolous block could be as the following:
!<```
!<                 _ J
!<                 /|                          _____
!<               5+ ...*----*----*----*----*...     |
!<               /    /    /    /    /    /         |
!<              /    /    /    /    /    /          |
!<            4+ ...*----*----*----*----*...        |
!<            /    /    /    /    /    /            |
!<           /    /    /    /    /    /             |
!<         3+ ...*----*----*----*----*...           |  Structured block of 4x4 Finite Volumes
!<         /    /    / FV /    /    /               |
!<        /    /    /    /    /    /                |
!<      2+ ...*----*----*----*----*...              |
!<      /    /    /    /    /    /                  |
!<     /    /    /    /    /    /                   |
!<   1+ ...*----*----*----*----*...                 |
!<   /     .    .    .    .    .                    |
!<  /      .    .    .    .    .               _____
!< O-------+----+----+----+----+-------------> I
!<         1    2    3    4    5
!<```

use adam_grid_object
use adam_parameters
use PENF
#ifdef _MPI_
use MPI
#endif

implicit none
private
public :: field_object

type :: field_object
   !< Field class definition.
   type(grid_object), pointer :: grid=>null()            !< Grid data.
   integer(I4P)               :: nv=1_I4P                !< Number of field variables.
   integer(I4P)               :: block_weight=0_I4P      !< Block weight, `cells_number * variables_number`.
   integer(I4P)               :: nb=0_I4P                !< Number of all blocks that can be stored.
   integer(I4P)               :: blocks_number=0_I4P     !< Number of blocks actually stored.
   integer(I8P), allocatable  :: code(:)                 !< Morton codes [nb].
   integer(I4P), allocatable  :: coordinates(:,:)        !< Coordinates IJKL for each block [nb,4].
   real(R8P),    allocatable  :: emin(:,:)               !< Coordinates of minimum abscissa of each block [3,nb].
   real(R8P),    allocatable  :: emax(:,:)               !< Coordinates of maximum abscissa of each block [3,nb].
   real(R8P),    allocatable  :: u(:,:,:,:)              !< Field cell centered variables [ni+gc12,nj+gc34,nk+gc56,nv,nb].
   real(R8P),    allocatable  :: u_new(:,:,:,:)          !< Field cell centered variables, buffer memory.
   logical                    :: is_initialized_=.false. !< Initialization status.
   ! MPI data
   integer(I4P)               :: myrank=0_I4P              !< MPI rank process.
   integer(I4P)               :: procs_number=1_I4P        !< Number of processes.
   integer(I4P), allocatable  :: blocks_numbers(:)         !< Number of blocks actually stored in all processes.
   integer(I4P), allocatable  :: refinements_needed(:)     !< Refinements needed of my blocks.
   integer(I4P), allocatable  :: refinements_needed_all(:) !< Refinements needed of all blocks.
   integer(I4P), allocatable  :: disp_count(:)             !< Displacement of blocks that are received from process.
   contains
      ! public methods
      procedure, pass(self) :: adapt                         !< Adapt field accordingly to refine/derefine necessity.
      procedure, pass(self) :: compute_emin_emax             !< Compute emin/emax of each block.
      procedure, pass(self) :: compute_xyz                   !< Compute grids coordinates from grids extents emin/emax.
      procedure, pass(self) :: destroy                       !< Destroy the field.
      procedure, pass(self) :: initialize                    !< Initialize the field.
      procedure, pass(self) :: mark_sphere                   !< Mark blocks to be refined/derefined by sphere distance.
      procedure, pass(self) :: max_cell_delta                !< Return the maximum cell delta given a comparison distance.
      procedure, pass(self) :: mpi_gather_refinements_needed !< Gather blocks refinement needed status between MPI processes.
      procedure, pass(self) :: mpi_redistribute              !< Redistribute blocks to processes.
      ! private methods
      procedure, pass(self), private :: derefine !< Derefine blocks.
      procedure, pass(self), private :: refine   !< Refine blocks.
      ! operators
      generic :: assignment(=) => field_assign_field      !< Overload `=`.
      procedure, pass(lhs), private :: field_assign_field !< Operator `=`.
endtype field_object

contains
   ! public methods
   subroutine adapt(self, ratio, block_to_refine, block_refined, block_to_derefine, block_derefined)
   !< Adapt field accordingly to refine/derefine necessity.
   class(field_object),       intent(inout) :: self                 !< The field.
   integer(I4P),              intent(in)    :: ratio                !< Refinement ratio.
   integer(I8P), allocatable, intent(in)    :: block_to_refine(:,:) !< List of field blocks to be refined.
   integer(I8P), allocatable, intent(in)    :: block_refined(:,:)   !< List of field refined blocks with Morton code.
   integer(I8P), allocatable, intent(in)    :: block_to_derefine(:) !< List of field blocks to be derefined.
   integer(I8P), allocatable, intent(in)    :: block_derefined(:,:) !< List of field derefined blocks with Morton code.

   call self%refine(  ratio=ratio, block_to_refine=block_to_refine,     block_refined=block_refined    )
   call self%derefine(ratio=ratio, block_to_derefine=block_to_derefine, block_derefined=block_derefined)
   endsubroutine adapt

   subroutine compute_emin_emax(self)
   !< Compute emin/emax of each block.
   class(field_object), intent(inout) :: self          !< The field.
   real(R8P)                          :: dx, dy, dz    !< Domain delta space.
   real(R8P)                          :: dxl, dyl, dzl !< Local delta space.
   integer(I4P)                       :: i, j, k, l, b !< Counter.

   dx = self%grid%domain_emax(1) - self%grid%domain_emin(1)
   dy = self%grid%domain_emax(2) - self%grid%domain_emin(2)
   dz = self%grid%domain_emax(3) - self%grid%domain_emin(3)
   do b=1, self%blocks_number
      i = self%coordinates(b,1)
      j = self%coordinates(b,2)
      k = self%coordinates(b,3)
      l = self%coordinates(b,4)
      dxl = dx / 2**l
      dyl = dy / 2**l
      dzl = dz / 2**l
      self%emin(1,b) = i * dxl ; self%emax(1,b) = self%emin(1,b) + dxl
      self%emin(2,b) = j * dyl ; self%emax(2,b) = self%emin(2,b) + dyl
      self%emin(3,b) = k * dzl ; self%emax(3,b) = self%emin(3,b) + dzl
   enddo
   endsubroutine compute_emin_emax

   pure function compute_xyz(self, b, axis) result(xyz)
   !< Compute grids coordinates from grids extents emin/emax of b-th block.
   class(field_object), intent(in) :: self   !< The field.
   integer(I4P),        intent(in) :: b      !< Block index.
   character(1),        intent(in) :: axis   !< Axis direction queried ['x','y','z'].
   real(R8P), allocatable          :: xyz(:) !< Grid coordinates.
   real(R8P)                       :: dxyz   !< Space delta.
   integer(I4P)                    :: i      !< Counter.

   associate(emin=>self%emin, emax=>self%emax, ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk, &
             gc1=>self%grid%gc1, gc2=>self%grid%gc2, gc3=>self%grid%gc3, gc4=>self%grid%gc4, gc5=>self%grid%gc5, gc6=>self%grid%gc6)
      select case(axis)
      case('x')
         allocate(xyz(0:ni))
         dxyz = (emax(1,b) - emin(1,b)) / ni
         do i=0, ni
            xyz(i) = emin(1,b) + i * dxyz
         enddo
      case('y')
         allocate(xyz(0:nj))
         dxyz = (emax(2,b) - emin(2,b)) / nj
         do i=0, nj
            xyz(i) = emin(2,b) + i * dxyz
         enddo
      case('z')
         allocate(xyz(0:nk))
         dxyz = (emax(3,b) - emin(3,b)) / nk
         do i=0, nk
            xyz(i) = emin(3,b) + i * dxyz
         enddo
      endselect
   endassociate
   endfunction compute_xyz

   subroutine destroy(self)
   !< Destroy field.
   class(field_object), intent(inout) :: self  !< The field.
   type(field_object)                 :: fresh !< Fresh field.

   self = fresh
   endsubroutine destroy

   subroutine initialize(self, grid, nv, nb)
   !< Initialize field.
   class(field_object), intent(inout)        :: self    !< The field.
   type(grid_object),   intent(in), target   :: grid    !< Grid data.
   integer(I4P),        intent(in), optional :: nv      !< Number of field variables.
   integer(I4P),        intent(in), optional :: nb      !< Number of all blocks that can be stored.
#ifdef _MPI_
   integer(I4P)                              :: error   !< Error traping flag.
#endif

   call self%destroy
   self%grid => grid
   if (present(nv)) self%nv  = nv
   self%block_weight = (self%grid%gc1+self%grid%ni+self%grid%gc2)* &
                       (self%grid%gc3+self%grid%nj+self%grid%gc4)* &
                       (self%grid%gc5+self%grid%nk+self%grid%gc6)*self%nv
   if (present(nb)) self%nb  = nb
   if (nb>0) then

      allocate(self%code(nb))
      self%code    = -2_I8P
      self%code(1) = -1_I8P ! first block is assumed to be ADAM

      allocate(self%coordinates(nb,4))

      allocate(self%emin(3,nb))
      allocate(self%emax(3,nb))
      self%emin(:,1) = self%grid%domain_emin
      self%emax(:,1) = self%grid%domain_emax

      allocate(self%u(1-self%grid%gc1:self%grid%ni+self%grid%gc2, &
                      1-self%grid%gc3:self%grid%nj+self%grid%gc4, &
                      1-self%grid%gc5:self%grid%nk+self%grid%gc6, 1:self%nb))
      allocate(self%u_new(1-self%grid%gc1:self%grid%ni+self%grid%gc2, &
                          1-self%grid%gc3:self%grid%nj+self%grid%gc4, &
                          1-self%grid%gc5:self%grid%nk+self%grid%gc6, 1:self%nb))
      self%u = 0._R8P
      self%u_new = 0._R8P
   endif
   self%is_initialized_ = .true.
   ! MPI data
#ifdef _MPI_
   call MPI_COMM_SIZE(MPI_COMM_WORLD, self%procs_number, error)
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, error)
   allocate(self%blocks_numbers(0:self%procs_number-1))
#endif
   endsubroutine initialize

   subroutine mark_sphere(self, center, radius, threshold)
   !< Mark blocks to be refined/derefined by sphere distance.
   class(field_object),       intent(inout)        :: self            !< The field.
   real(R8P),                 intent(in)           :: center(3)       !< Sphere center coordinates [x,y,z].
   real(R8P),                 intent(in)           :: radius          !< Sphere radius.
   real(R8P),                 intent(in), optional :: threshold       !< Threshold for sphere proximity.
   real(R8P)                                       :: threshold_      !< Threshold for sphere proximity, local var.
   real(R8P)                                       :: block_center(3) !< block center coordinates.
   real(R8P)                                       :: block_diagonal  !< block diagonal.
   real(R8P)                                       :: distance(0:8)   !< Distances between block and sphere.
   real(R8P)                                       :: max_cell_delta  !< Max cell delta.
   integer(I4P)                                    :: b               !< Counter.

   threshold_ = 2.2_R8P ; if (present(threshold)) threshold_ = threshold
   if (allocated(self%refinements_needed)) deallocate(self%refinements_needed)
   allocate(self%refinements_needed(self%blocks_number))
   do b=1, self%blocks_number
      block_center = (self%emax(:,b) + self%emin(:,b)) / 2._R8P
      block_diagonal = sqrt((self%emax(1,b) - self%emin(1,b))**2 + &
                            (self%emax(2,b) - self%emin(2,b))**2 + &
                            (self%emax(3,b) - self%emin(3,b))**2)

      associate (emin=>self%emin(:,b), emax=>self%emax(:,b), ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk)
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

      max_cell_delta = self%max_cell_delta(distance=distance(0))

      if (block_diagonal/min(ni,nj,nk) > max_cell_delta) then
         self%refinements_needed(b) = TO_BE_REFINED
      elseif (block_diagonal/min(ni,nj,nk) * threshold_ < max_cell_delta) then
         self%refinements_needed(b) = TO_BE_DEREFINED
      else
         self%refinements_needed(b) = TO_NOT_TOUCH
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

   function max_cell_delta(self, distance) result(delta)
   !< Return the maximum cell delta given a comparison distance.
   class(field_object), intent(in) :: self     !< The field.
   real(R8P),           intent(in) :: distance !< Comparison distance.
   real(R8P)                       :: delta    !< Maximum cell delta admissible.

   if (abs(distance) < epsilon(0._R8P)) then
      ! delta = 0.001_R8P
      delta = 0.005_R8P
   else
      delta = huge(0._R8P)
   endif
   endfunction max_cell_delta

   subroutine mpi_gather_refinements_needed(self)
   !< Gather blocks refinement needed status between MPI processes.
   class(field_object),       intent(inout) :: self          !< The field.
   integer(I4P), allocatable                :: recv_count(:) !< Number of blocks that are received from process.
   integer(I8P)                             :: p             !< Counter.
#ifdef _MPI_
   integer(I4P)                             :: error         !< Error traping flag.
#endif

   ! computing received blocks
   allocate(recv_count(0:self%procs_number - 1))
   call MPI_ALLGATHER(self%blocks_number, 1_I4P, MPI_INTEGER, &
                      recv_count, 1_I4P, MPI_INTEGER, MPI_COMM_WORLD, error)

   ! computing displacement counts
   if (allocated(self%disp_count)) deallocate(self%disp_count)
   allocate(self%disp_count(0:self%procs_number - 1))
   self%disp_count = 0_I4P
   do p=1, self%procs_number - 1
      self%disp_count(p) = self%disp_count(p-1) + recv_count(p-1)
   enddo

   if (allocated(self%refinements_needed_all)) deallocate(self%refinements_needed_all)
   allocate(self%refinements_needed_all(sum(recv_count, dim=1)))
#ifdef _MPI_
   call MPI_ALLGATHERV(self%refinements_needed, self%blocks_number, MPI_INTEGER, &
                       self%refinements_needed_all, recv_count, self%disp_count, MPI_INTEGER, MPI_COMM_WORLD, error)
#endif
   endsubroutine mpi_gather_refinements_needed

   subroutine mpi_redistribute(self, comm_map_send, comm_map_recv, comm_map_send_ptr, comm_map_recv_ptr, local_map, coordinates)
   !< Redistribute blocks to processes.
   class(field_object),       intent(inout) :: self                   !< The field.
   integer(I8P), allocatable, intent(in)    :: comm_map_send(:)       !< Comm map, blocks to send [sum(comm_map_n_send)].
   integer(I8P), allocatable, intent(in)    :: comm_map_recv(:)       !< Comm map, blocks to receive [sum(comm_map_n_recv)].
   integer(I4P), allocatable, intent(in)    :: comm_map_send_ptr(:)   !< Comm map, pointers in list to send [procs_number+1].
   integer(I4P), allocatable, intent(in)    :: comm_map_recv_ptr(:)   !< Comm map, pointers in list to recv [procs_number+1].
   integer(I4P), allocatable, intent(in)    :: coordinates(:,:)       !< Coordinates (ijkl,nb) of redistributed nodes.
   integer(I8P), allocatable                :: local_map(:,:)         !< Local map, list block index changes of my blocks.
   real(R8P),    allocatable                :: send_buffer(:)         !< Send buffer of field cell centered variables.
   real(R8P),    allocatable                :: recv_buffer(:)         !< Recv buffer of field cell centered variables.
   integer(I8P)                             :: send_size, send_offset !< Total size of send buffer.
   integer(I8P)                             :: recv_size, recv_offset !< Total size of recv buffer.
   integer(I4P)                             :: n_keep                 !< Number of keept blocks.
   integer(I4P)                             :: b, bi, p               !< Counter.
   integer(I4P)                             :: ptr_start, ptr_end     !< Counter.
   integer(I4P)                             :: n_recv, n_send         !< Counter.
#ifdef _MPI_
   integer(I4P), allocatable                :: req_recv(:)            !< MPI request receive flags.
   integer(I4P)                             :: error                  !< Error traping flag.

   allocate(req_recv(0:self%procs_number-1))
   req_recv = MPI_REQUEST_NULL
#endif

   send_size = 0_I8P ; if (allocated(comm_map_send)) send_size = size(comm_map_send, dim=1) * self%block_weight
   recv_size = 0_I8P ; if (allocated(comm_map_recv)) recv_size = size(comm_map_recv, dim=1) * self%block_weight
   n_keep    = 0_I8P ; if (allocated(local_map    )) n_keep    = size(local_map    , dim=1)
   if (send_size > 0_I8P) allocate(send_buffer(send_size))
   if (recv_size > 0_I8P) allocate(recv_buffer(recv_size))

   if (send_size > 0_I8P) then
      send_offset = 1
      do b=1, size(comm_map_send, dim=1)
         bi = comm_map_send(b)
         send_buffer(send_offset:send_offset + self%block_weight - 1) = reshape(self%u(:,:,:,bi),[self%block_weight])
         send_offset = send_offset + self%block_weight
      enddo
   endif

   do p=0, self%procs_number - 1_I4P
      ptr_start = comm_map_recv_ptr(p)   * self%block_weight + 1
      ptr_end   = comm_map_recv_ptr(p+1) * self%block_weight
      n_recv    = ptr_end - ptr_start + 1
      if (n_recv > 0) then
#ifdef _MPI_
         call MPI_IRECV(recv_buffer(ptr_start), n_recv, MPI_REAL8, p, 100, &
                        MPI_COMM_WORLD, req_recv(p), error)
#endif
      endif
   enddo

   do p=0, self%procs_number - 1_I4P
      ptr_start = comm_map_send_ptr(p)   * self%block_weight + 1
      ptr_end   = comm_map_send_ptr(p+1) * self%block_weight
      n_send    = ptr_end - ptr_start + 1
      if (n_send > 0) then
#ifdef _MPI_
         call MPI_SEND(send_buffer(ptr_start), n_send, MPI_REAL8, p, 100, &
                       MPI_COMM_WORLD, error)
#endif
      endif
   enddo

#ifdef _MPI_
   call MPI_WAITALL(self%procs_number, req_recv, MPI_STATUSES_IGNORE, error)
#endif

   if (recv_size > 0_I8P) then
      recv_offset = 1
      do b=1, size(comm_map_recv, dim=1)
          bi = comm_map_recv(b)
          self%u_new(:,:,:,bi) = reshape(recv_buffer(recv_offset:recv_offset + self%block_weight -1),&
                                         [self%grid%gc1+self%grid%ni+self%grid%gc2, &
                                          self%grid%gc3+self%grid%nj+self%grid%gc4, &
                                          self%grid%gc5+self%grid%nk+self%grid%gc6])
          recv_offset = recv_offset + self%block_weight
      enddo
   endif

   do b=1, n_keep
      self%u_new(:,:,:,local_map(b,1)) = self%u(:,:,:,local_map(b,2))
   enddo
   self%u = self%u_new
   self%blocks_number = n_keep  + recv_size / self%block_weight
   self%coordinates(1:self%blocks_number,:) = coordinates
   call self%compute_emin_emax
   endsubroutine mpi_redistribute

   ! private methods
   pure subroutine derefine(self, ratio, block_to_derefine, block_derefined)
   !< Derefine blocks.
   class(field_object),       intent(inout) :: self                 !< The field.
   integer(I4P),              intent(in)    :: ratio                !< Refinement ratio.
   integer(I8P), allocatable, intent(in)    :: block_to_derefine(:) !< List of blocks to be derefined.
   integer(I8P), allocatable, intent(in)    :: block_derefined(:,:) !< List of derefined blocks with Morton code.
   real(R8P)                                :: dx, dy, dz           !< Space deltas.
   integer(I4P)                             :: b, ib                !< Counter.
   integer(I4P)                             :: ic1, ic2, ic3, ic4   !< Counter.
   integer(I4P)                             :: ic5, ic6, ic7, ic8   !< Counter.

   if (allocated(block_derefined)) then
      do b=1, size(block_derefined, dim=2)
         ib = block_derefined(2,b)

         ic1 = block_to_derefine((b-1)*ratio+1)
         ic2 = block_to_derefine((b-1)*ratio+2)
         ic3 = block_to_derefine((b-1)*ratio+3)
         ic4 = block_to_derefine((b-1)*ratio+4)
         ic5 = block_to_derefine((b-1)*ratio+5)
         ic6 = block_to_derefine((b-1)*ratio+6)
         ic7 = block_to_derefine((b-1)*ratio+7)
         ic8 = block_to_derefine((b-1)*ratio+8)

         self%u(:,:,:,ib) = block_derefined(1,b) ; self%code(ib) = block_derefined(1,b)

      enddo

      do b=1, size(block_derefined, dim=2)
         ib = block_derefined(2,b)

         ic1 = block_to_derefine((b-1)*ratio+1)

         dx = self%emax(1,ic1) - self%emin(1,ic1)
         dy = self%emax(2,ic1) - self%emin(2,ic1)
         dz = self%emax(3,ic1) - self%emin(3,ic1)

         self%emin(1,ib) = self%emin(1,ic1)
         self%emin(2,ib) = self%emin(2,ic1)
         self%emin(3,ib) = self%emin(3,ic1)

         self%emax(1,ib) = self%emin(1,ic1) + 2 * dx
         self%emax(2,ib) = self%emin(2,ic1) + 2 * dy
         self%emax(3,ib) = self%emin(3,ic1) + 2 * dz
      enddo
   endif
   endsubroutine derefine

   pure subroutine refine(self, ratio, block_to_refine, block_refined)
   !< Refine blocks.
   class(field_object),       intent(inout) :: self                 !< The field.
   integer(I4P),              intent(in)    :: ratio                !< Refinement ratio.
   integer(I8P), allocatable, intent(in)    :: block_to_refine(:,:) !< List of blocks to be refined.
   integer(I8P), allocatable, intent(in)    :: block_refined(:,:)   !< List of refined blocks with Morton code.
   real(R8P)                                :: dx, dy, dz           !< Space deltas.
   integer(I4P)                             :: b, i, j, k           !< Spatial counter.
   integer(I4P)                             :: ib, ic, ii           !< Counter.
   integer(I4P)                             :: ic1, ic2, ic3, ic4   !< Counter.
   integer(I4P)                             :: ic5, ic6, ic7, ic8   !< Counter.

   if (allocated(block_to_refine)) then
      do b=1, size(block_to_refine, dim=2)
         if (self%myrank /= block_to_refine(2,b)) cycle
         ib = block_to_refine(1,b)

         ic1 = block_refined(2,(b-1)*ratio+1)
         ic2 = block_refined(2,(b-1)*ratio+2)
         ic3 = block_refined(2,(b-1)*ratio+3)
         ic4 = block_refined(2,(b-1)*ratio+4)
         ic5 = block_refined(2,(b-1)*ratio+5)
         ic6 = block_refined(2,(b-1)*ratio+6)
         ic7 = block_refined(2,(b-1)*ratio+7)
         ic8 = block_refined(2,(b-1)*ratio+8)

         self%u(:,:,:,ic1) = block_refined(1,(b-1)*ratio+1) ; self%code(ic1) = block_refined(1,(b-1)*ratio+1)
         self%u(:,:,:,ic2) = block_refined(1,(b-1)*ratio+2) ; self%code(ic2) = block_refined(1,(b-1)*ratio+2)
         self%u(:,:,:,ic3) = block_refined(1,(b-1)*ratio+3) ; self%code(ic3) = block_refined(1,(b-1)*ratio+3)
         self%u(:,:,:,ic4) = block_refined(1,(b-1)*ratio+4) ; self%code(ic4) = block_refined(1,(b-1)*ratio+4)
         self%u(:,:,:,ic5) = block_refined(1,(b-1)*ratio+5) ; self%code(ic5) = block_refined(1,(b-1)*ratio+5)
         self%u(:,:,:,ic6) = block_refined(1,(b-1)*ratio+6) ; self%code(ic6) = block_refined(1,(b-1)*ratio+6)
         self%u(:,:,:,ic7) = block_refined(1,(b-1)*ratio+7) ; self%code(ic7) = block_refined(1,(b-1)*ratio+7)
         self%u(:,:,:,ic8) = block_refined(1,(b-1)*ratio+8) ; self%code(ic8) = block_refined(1,(b-1)*ratio+8)
      enddo

      do b=1, size(block_to_refine, dim=2)
         if (self%myrank /= block_to_refine(2,b)) cycle
         ib = block_to_refine(1,b)

         dx = self%emax(1,ib) - self%emin(1,ib)
         dy = self%emax(2,ib) - self%emin(2,ib)
         dz = self%emax(3,ib) - self%emin(3,ib)

         ii = 1
         do k=0, 1
            do j=0, 1
               do i=0, 1
                  ic = block_refined(2,(b-1)*ratio + ii)

                  self%emin(1,ic) = self%emin(1,ib) + i * dx/2
                  self%emin(2,ic) = self%emin(2,ib) + j * dy/2
                  self%emin(3,ic) = self%emin(3,ib) + k * dz/2

                  self%emax(1,ic) = self%emin(1,ic) + dx/2
                  self%emax(2,ic) = self%emin(2,ic) + dy/2
                  self%emax(3,ic) = self%emin(3,ic) + dz/2

                  ii = ii + 1
               enddo
            enddo
         enddo
      enddo
   endif
   endsubroutine refine

   ! operators
   ! =
   subroutine field_assign_field(lhs, rhs)
   !< Operator `=`.
   class(field_object), intent(inout) :: lhs !< Left hand side.
   type(field_object),  intent(in)    :: rhs !< Right hand side.

   lhs%grid => rhs%grid
   lhs%nv            = rhs%nv
   lhs%block_weight  = rhs%block_weight
   lhs%nb            = rhs%nb
   lhs%blocks_number = rhs%blocks_number
   if (allocated(rhs%code)) then
      lhs%code = rhs%code
   else
      if (allocated(lhs%code)) deallocate(lhs%code)
   endif
   if (allocated(rhs%coordinates)) then
      lhs%coordinates = rhs%coordinates
   else
      if (allocated(lhs%coordinates)) deallocate(lhs%coordinates)
   endif
   if (allocated(rhs%emin)) then
      lhs%emin = rhs%emin
   else
      if (allocated(lhs%emin)) deallocate(lhs%emin)
   endif
   if (allocated(rhs%emax)) then
      lhs%emax = rhs%emax
   else
      if (allocated(lhs%emax)) deallocate(lhs%emax)
   endif
   if (allocated(rhs%u)) then
      lhs%u = rhs%u
   else
      if (allocated(lhs%u)) deallocate(lhs%u)
   endif
   if (allocated(rhs%u_new)) then
      lhs%u_new = rhs%u_new
   else
      if (allocated(lhs%u_new)) deallocate(lhs%u_new)
   endif
   lhs%is_initialized_ = rhs%is_initialized_
   ! MPI data
   lhs%myrank        = rhs%myrank
   lhs%procs_number  = rhs%procs_number
   if (allocated(rhs%blocks_numbers)) then
      lhs%blocks_numbers = rhs%blocks_numbers
   else
      if (allocated(lhs%blocks_numbers)) deallocate(lhs%blocks_numbers)
   endif
   endsubroutine field_assign_field
endmodule adam_field_object
