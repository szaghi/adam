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

use PENF, only : I8P, I4P, R8P, str
#ifdef _MPI_
use MPI
#endif

implicit none
private
public :: field_object

type :: field_object
   integer(I4P)              :: ni=4_I4P                !< Number of cells in i direction.
   integer(I4P)              :: nj=4_I4P                !< Number of cells in j direction.
   integer(I4P)              :: nk=4_I4P                !< Number of cells in k direction.
   integer(I4P)              :: gc1=0_I4P               !< Number of ghost cells in i- direction for boundary conditions.
   integer(I4P)              :: gc2=0_I4P               !< Number of ghost cells in i+ direction for boundary conditions.
   integer(I4P)              :: gc3=0_I4P               !< Number of ghost cells in j- direction for boundary conditions.
   integer(I4P)              :: gc4=0_I4P               !< Number of ghost cells in j+ direction for boundary conditions.
   integer(I4P)              :: gc5=0_I4P               !< Number of ghost cells in k- direction for boundary conditions.
   integer(I4P)              :: gc6=0_I4P               !< Number of ghost cells in k+ direction for boundary conditions.
   integer(I4P)              :: nv=1_I4P                !< Number of field variables.
   integer(I4P)              :: block_weight=0_I4P      !< Block weight, `cells_number * variables_number`.
   integer(I4P)              :: nb=0_I4P                !< Number of all blocks that can be stored.
   integer(I4P)              :: blocks_number=0_I4P     !< Number of blocks actually stored.
   integer(I8P), allocatable :: code(:)                 !< Morton codes [nb].
   integer(I4P), allocatable :: coordinates(:,:)        !< Coordinates IJKL for each block [nb,4].
   real(R8P)                 :: domain_emin(3)          !< Coordinates of minimum abscissa of whole domain [3].
   real(R8P)                 :: domain_emax(3)          !< Coordinates of maximum abscissa of whole domain [3].
   real(R8P),    allocatable :: emin(:,:)               !< Coordinates of minimum abscissa of each block [3,nb].
   real(R8P),    allocatable :: emax(:,:)               !< Coordinates of maximum abscissa of each block [3,nb].
   real(R8P),    allocatable :: u(:,:,:,:)              !< Field cell centered variables [ni+gc12,nj+gc34,nk+gc56,nv,nb].
   real(R8P),    allocatable :: u_new(:,:,:,:)          !< Field cell centered variables, buffer memory.
   logical                   :: is_initialized_=.false. !< Initialization status.
   ! MPI data
   integer(I4P) :: myrank=0_I4P       !< MPI rank process.
   integer(I4P) :: procs_number=1_I4P !< Number of processes.
   contains
      ! public methods
      procedure, pass(self) :: adapt             !< Adapt field accordingly to refine/derefine necessity.
      procedure, pass(self) :: compute_emin_emax !< Compute emin/emax of each block.
      procedure, pass(self) :: compute_xyz       !< Compute grids coordinates from grids extents emin/emax.
      procedure, pass(self) :: destroy           !< Destroy the field.
      procedure, pass(self) :: initialize        !< Initialize the field.
      procedure, pass(self) :: max_cell_delta    !< Return the maximum cell delta given a comparison distance.
      procedure, pass(self) :: redistribute      !< Redistribute blocks to processes.
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
   integer(I8P), allocatable, intent(in)    :: block_to_refine(:)   !< List of field blocks to be refined.
   integer(I8P), allocatable, intent(in)    :: block_refined(:,:)   !< List of field refined blocks with Morton code.
   integer(I8P), allocatable, intent(in)    :: block_to_derefine(:) !< List of field blocks to be derefined.
   integer(I8P), allocatable, intent(in)    :: block_derefined(:,:) !< List of field derefined blocks with Morton code.

   call self%refine(  ratio=ratio, block_to_refine=block_to_refine,     block_refined=block_refined    )
   call self%derefine(ratio=ratio, block_to_derefine=block_to_derefine, block_derefined=block_derefined)
   endsubroutine adapt

   pure subroutine compute_emin_emax(self)
   !< Compute emin/emax of each block.
   class(field_object), intent(inout) :: self          !< The field.
   real(R8P)                          :: dx, dy, dz    !< Domain delta space.
   real(R8P)                          :: dxl, dyl, dzl !< Local delta space.
   integer(I4P)                       :: i, j, k, l, b !< Counter.

   dx = self%domain_emax(1) - self%domain_emin(1)
   dy = self%domain_emax(2) - self%domain_emin(2)
   dz = self%domain_emax(3) - self%domain_emin(3)
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

   associate(emin=>self%emin, emax=>self%emax, ni=>self%ni, nj=>self%nj, nk=>self%nk, &
             gc1=>self%gc1, gc2=>self%gc2, gc3=>self%gc3,  gc4=>self%gc4, gc5=>self%gc5, gc6=>self%gc6)
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

   elemental subroutine destroy(self)
   !< Destroy field.
   class(field_object), intent(inout) :: self  !< The field.
   type(field_object)                 :: fresh !< Fresh field.

   self = fresh
   endsubroutine destroy

   subroutine initialize(self, ni, nj, nk, gc, nv, nb, emin, emax)
   !< Initialize field.
   class(field_object), intent(inout)        :: self    !< The field.
   integer(I4P),        intent(in), optional :: ni      !< Number of cells in X direction.
   integer(I4P),        intent(in), optional :: nj      !< Number of cells in Y direction.
   integer(I4P),        intent(in), optional :: nk      !< Number of cells in Z direction.
   integer(I4P),        intent(in), optional :: gc(6)   !< Number of ghost cells in each direction.
   integer(I4P),        intent(in), optional :: nv      !< Number of field variables.
   integer(I4P),        intent(in), optional :: nb      !< Number of all blocks that can be stored.
   real(R8P),           intent(in), optional :: emin(3) !< Coordinates of minium abscissa.
   real(R8P),           intent(in), optional :: emax(3) !< Coordinates of maxium abscissa.
   integer(I4P)                              :: error   !< Error traping flag.

   call self%destroy
   if (present(ni)) self%ni  = ni
   if (present(nj)) self%nj  = nj
   if (present(nk)) self%nk  = nk
   if (present(gc)) self%gc1 = gc(1)
   if (present(gc)) self%gc2 = gc(2)
   if (present(gc)) self%gc3 = gc(3)
   if (present(gc)) self%gc4 = gc(4)
   if (present(gc)) self%gc5 = gc(5)
   if (present(gc)) self%gc6 = gc(6)
   if (present(nv)) self%nv  = nv
   if (present(nb)) self%nb  = nb
   if (nb>0) then
      self%block_weight = (self%gc1+self%ni+self%gc2)* &
                          (self%gc3+self%nj+self%gc4)* &
                          (self%gc5+self%nk+self%gc6)*self%nv

      allocate(self%code(nb))
      allocate(self%coordinates(nb,4))
      self%code    = -2_I8P
      self%code(1) = -1_I8P ! first block is assumed to be ADAM

      allocate(self%emin(3,nb))
      allocate(self%emax(3,nb))
      if (present(emin)) then
         self%domain_emin = emin
         self%emin(:,1) = emin
      else
         self%domain_emin = 0._R8P
         self%emin(:,1) = 0._R8P
      endif
      if (present(emax)) then
         self%domain_emax = emax
         self%emax(:,1) = emax
      else
         self%domain_emax = 1._R8P
         self%emax(:,1) = 1._R8P
      endif

      allocate(self%u(1-self%gc1:self%ni+self%gc2, &
                      1-self%gc3:self%nj+self%gc4, &
                      1-self%gc5:self%nk+self%gc6, 1:self%nb))
      allocate(self%u_new(1-self%gc1:self%ni+self%gc2, &
                          1-self%gc3:self%nj+self%gc4, &
                          1-self%gc5:self%nk+self%gc6, 1:self%nb))
      self%u = 0._R8P
   endif
   self%is_initialized_ = .true.
   ! MPI data
#ifdef _MPI_
   call MPI_COMM_SIZE(MPI_COMM_WORLD, self%procs_number, error)
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, error)
#endif
   endsubroutine initialize

   function max_cell_delta(self, distance) result(delta)
   !< Return the maximum cell delta given a comparison distance.
   class(field_object), intent(in) :: self     !< The field.
   real(R8P),           intent(in) :: distance !< Comparison distance.
   real(R8P)                       :: delta    !< Maximum cell delta admissible.

   if (abs(distance) < epsilon(0._R8P)) then
      delta = 0.005_R8P
   else
      delta = huge(0._R8P)
   endif
   endfunction max_cell_delta

   subroutine redistribute(self, comm_map_send, comm_map_recv, comm_map_send_ptr, comm_map_recv_ptr, local_map, coordinates)
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
   integer(I4P)                             :: error                  !< Error traping flag.
   integer(I4P), allocatable                :: req_recv(:)            !< MPI request receive flags.

   allocate(req_recv(0:self%procs_number-1))
   req_recv = MPI_REQUEST_NULL

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

   call MPI_BARRIER(MPI_COMM_WORLD, error)
   print*, 'cazzooooooooooooooooooooooooooooooooooooooooooooooooooooooooo'
   do p=0, self%procs_number - 1_I4P
      ptr_start = comm_map_recv_ptr(p)   * self%block_weight + 1
      ptr_end   = comm_map_recv_ptr(p+1) * self%block_weight
      n_recv    = ptr_end - ptr_start + 1
      print*, 'cazzo recv ', n_recv, ' da ', p
      if (n_recv > 0) then
         call MPI_IRECV(recv_buffer(ptr_start), n_recv, MPI_REAL8, p, 100, &
                        MPI_COMM_WORLD, req_recv(p), error)
      endif
   enddo

   do p=0, self%procs_number - 1_I4P
      ptr_start = comm_map_send_ptr(p)   * self%block_weight + 1
      ptr_end   = comm_map_send_ptr(p+1) * self%block_weight
      n_send    = ptr_end - ptr_start + 1
      print*, 'cazzo send ', n_send, ' a ', p
      if (n_send > 0) then
         call MPI_SEND(send_buffer(ptr_start), n_send, MPI_REAL8, p, 100, &
                       MPI_COMM_WORLD, error)
      endif
   enddo

   print*, ' cazzo prima wait'
   call MPI_WAITALL(self%procs_number, req_recv, MPI_STATUSES_IGNORE, error)
   print*, ' cazzo dopo wait'

   if (recv_size > 0_I8P) then
      recv_offset = 1
      do b=1, size(comm_map_recv, dim=1)
          bi = comm_map_recv(b)
          self%u_new(:,:,:,bi) = reshape(recv_buffer(recv_offset:recv_offset + self%block_weight -1),[self%gc1+self%ni+self%gc2, &
                                                                                                      self%gc3+self%nj+self%gc4, &
                                                                                                      self%gc5+self%nk+self%gc6])
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
   endsubroutine redistribute

   ! privatec methods
   pure subroutine derefine(self, ratio, block_to_derefine, block_derefined)
   !< Derefine blocks.
   class(field_object),       intent(inout) :: self                 !< The field.
   integer(I4P),              intent(in)    :: ratio                !< Refinement ratio.
   integer(I8P), allocatable, intent(in)    :: block_to_derefine(:) !< List of blocks to be derefined.
   integer(I8P), allocatable, intent(in)    :: block_derefined(:,:) !< List of derefined blocks with Morton code.
   real(R8P)                                :: dx, dy, dz           !< Space deltas.
   integer(I4P)                             :: b, i, j, k           !< Spatial counter.
   integer(I4P)                             :: ib, ic, ii           !< Counter.
   integer(I4P)                             :: ic1, ic2, ic3, ic4   !< Counter.
   integer(I4P)                             :: ic5, ic6, ic7, ic8   !< Counter.

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
   endsubroutine derefine

   pure subroutine refine(self, ratio, block_to_refine, block_refined)
   !< Refine blocks.
   class(field_object),       intent(inout) :: self               !< The field.
   integer(I4P),              intent(in)    :: ratio              !< Refinement ratio.
   integer(I8P), allocatable, intent(in)    :: block_to_refine(:) !< List of blocks to be refined.
   integer(I8P), allocatable, intent(in)    :: block_refined(:,:) !< List of refined blocks with Morton code.
   real(R8P)                                :: dx, dy, dz         !< Space deltas.
   integer(I4P)                             :: b, i, j, k         !< Spatial counter.
   integer(I4P)                             :: ib, ic, ii         !< Counter.
   integer(I4P)                             :: ic1, ic2, ic3, ic4 !< Counter.
   integer(I4P)                             :: ic5, ic6, ic7, ic8 !< Counter.

   self%blocks_number = maxval(block_refined(2,:))
   do b=1, size(block_to_refine, dim=1)
      ib = block_to_refine(b)

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

   do b=1, size(block_to_refine, dim=1)
      ib = block_to_refine(b)

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
   endsubroutine refine

   ! operators
   ! =
   pure subroutine field_assign_field(lhs, rhs)
   !< Operator `=`.
   class(field_object), intent(inout) :: lhs !< Left hand side.
   type(field_object),  intent(in)    :: rhs !< Right hand side.

   lhs%ni            = rhs%ni
   lhs%nj            = rhs%nj
   lhs%nk            = rhs%nk
   lhs%gc1           = rhs%gc1
   lhs%gc2           = rhs%gc2
   lhs%gc3           = rhs%gc3
   lhs%gc4           = rhs%gc4
   lhs%gc5           = rhs%gc5
   lhs%gc6           = rhs%gc6
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

   lhs%domain_emin = rhs%domain_emin
   lhs%domain_emax = rhs%domain_emax

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
   endsubroutine field_assign_field
endmodule adam_field_object
