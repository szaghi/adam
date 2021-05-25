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

use adam_grid_object, only : grid_object
use adam_parameters
use FINER, only : file_ini
use PENF
use MPI

implicit none
private
public :: field_object

type :: field_object
   !< Field class definition.
   ! mesh related data, unrelated to field equations
   type(grid_object), pointer :: grid=>null()             !< Grid data.
   integer(I4P)               :: nv=1_I4P                 !< Number of field variables.
   integer(I4P)               :: block_weight=0_I4P       !< Block weight, `cells_number * variables_number`.
   integer(I4P)               :: nb=0_I4P                 !< Number of all blocks that can be stored.
   integer(I4P)               :: blocks_number=0_I4P      !< Number of blocks actually stored.
   integer(I8P), allocatable  :: code(:)                  !< Morton codes [nb].
   integer(I4P), allocatable  :: coordinates(:,:)         !< Coordinates IJKL for each block [nb,4].
   real(R8P),    allocatable  :: emin(:,:)                !< Coordinates of minimum abscissa of each block [3,nb].
   real(R8P),    allocatable  :: emax(:,:)                !< Coordinates of maximum abscissa of each block [3,nb].
   real(R8P),    allocatable  :: dxyz(:,:)                !< Space steps of each block [3,nb].
   real(R8P),    allocatable  :: x_node(:,:)              !< X node coordinates.
   real(R8P),    allocatable  :: y_node(:,:)              !< Y node coordinates.
   real(R8P),    allocatable  :: z_node(:,:)              !< Z node coordinates.
   real(R8P),    allocatable  :: x_cell(:,:)              !< X cell coordinates.
   real(R8P),    allocatable  :: y_cell(:,:)              !< Y cell coordinates.
   real(R8P),    allocatable  :: z_cell(:,:)              !< Z cell coordinates.
   integer(I8P), allocatable  :: local_map_ghost(:,:)     !< Local map for ghost cells updating.
   integer(I8P), allocatable  :: local_map_bc_face(:,:)   !< Local map for face BC ghost cells.
   integer(I8P), allocatable  :: local_map_bc_edge(:,:)   !< Local map for edge BC ghost cells.
   integer(I8P), allocatable  :: local_map_bc_corner(:,:) !< Local map for corner BC ghost cells.
   ! MPI data, unrelated to field equations
   integer(I4P)              :: error=0_I4P                !< Error traping flag.
   integer(I4P)              :: myrank=0_I4P               !< MPI rank process.
   integer(I4P)              :: procs_number=1_I4P         !< Number of processes.
   integer(I4P), allocatable :: blocks_numbers(:)          !< Number of blocks actually stored in all processes.
   integer(I4P), allocatable :: refinements_needed(:)      !< Refinements needed of my blocks.
   integer(I4P), allocatable :: refinements_needed_all(:)  !< Refinements needed of all blocks.
   integer(I4P), allocatable :: disp_count(:)              !< Displacement of blocks that are received from process.
   integer(I4P)              :: inner_blocks_number=0_I4P  !< Number of inner blocks where I need fecs.
   integer(I4P), allocatable :: req_send_recv(:)           !< MPI request receive flags.
   integer(I4P), allocatable :: comm_map_n_send_ghost(:)   !< Communication map, number of cells to send [procs_number].
   integer(I4P), allocatable :: comm_map_n_recv_ghost(:)   !< Communication map, number of cells to recv [procs_number].
   integer(I4P), allocatable :: comm_map_send_ptr_ghost(:) !< Communication map, pointers in list to send [procs_number+1].
   integer(I4P), allocatable :: comm_map_recv_ptr_ghost(:) !< Communication map, pointers in list to recv [procs_number+1].
   integer(I4P), allocatable :: comm_map_send_ptr_ghost_s(:) !< Communication map, pointers in list to send, single var.
   integer(I4P), allocatable :: comm_map_recv_ptr_ghost_s(:) !< Communication map, pointers in list to recv, single var.
   integer(I8P), allocatable :: comm_map_send_ghost(:,:)   !< Communication map, `fec` information [fec_number, 15].
   integer(I8P), allocatable :: comm_map_recv_ghost(:,:)   !< Communication map, `fec` information [fec_number, 15].
   integer(I8P), allocatable :: comm_map_send_ghost_s(:,:)   !< Communication map, `fec` information [fec_number, 15], single var.
   integer(I8P), allocatable :: comm_map_recv_ghost_s(:,:)   !< Communication map, `fec` information [fec_number, 15], single var.
   ! MPI data, related to field equations
   real(R8P), allocatable :: send_buffer_ghost(:)   !< Send buffer of ghost cells.
   real(R8P), allocatable :: recv_buffer_ghost(:)   !< Receive buffer of ghost cells.
   real(R8P), allocatable :: send_buffer_ghost_s(:) !< Send buffer of ghost cells, single var.
   real(R8P), allocatable :: recv_buffer_ghost_s(:) !< Receive buffer of ghost cells, single var.
   ! field equations data
   real(R8P), allocatable :: q(     :,:,:,:,:) !< Field cell centered variables.
   real(R8P), allocatable :: q_work(:,:,:,:,:) !< Field cell centered variables, working buffer memory.
   contains
      ! public methods
      procedure, pass(self) :: adapt                         !< Adapt field accordingly to refine/derefine necessity.
      procedure, pass(self) :: blocks_reorder                !< Reorder blocks indexes in field.
      procedure, pass(self) :: compute_metrics               !< Compute metrics of each block.
      procedure, pass(self) :: destroy                       !< Destroy the field.
      procedure, pass(self) :: initialize                    !< Initialize the field.
      procedure, pass(self) :: load_from_ini_file            !< Load object data from INI file.
      procedure, pass(self) :: mark_all_blocks               !< Mark all blocks to be refined, derefined, ecc.
      procedure, pass(self) :: mark_sphere                   !< Mark blocks to be refined/derefined by sphere distance.
      procedure, pass(self) :: mpi_gather_refinements_needed !< Gather blocks refinement needed status between MPI processes.
      procedure, pass(self) :: mpi_redistribute              !< Redistribute blocks to processes.
      procedure, pass(self) :: prepare_comm_local_ghost      !< Prepare communication and local maps/buffers for ghosts update.
      procedure, pass(self) :: prepare_local_bc              !< Prepare local maps for boundary conditions.
      procedure, pass(self) :: print_status                  !< Print status of main data.
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

   subroutine blocks_reorder(self, inner_outer_block_map, inner_blocks_number)
   !< Reorder blocks indexes in field.
   class(field_object), intent(inout) :: self                     !< The field.
   integer(I4P),        intent(in)    :: inner_outer_block_map(:) !< Inner/outer blocks map.
   integer(I4P),        intent(in)    :: inner_blocks_number      !< Number of inner blocks where I need fecs.
   integer(I4P), allocatable          :: coordinates_new(:,:)     !< Temporary coordinates array.
   integer(I8P), allocatable          :: code_new(:)              !< Temporary Morton codes.
   integer(I4P)                       :: b                        !< Counter.

   allocate(coordinates_new(4,self%blocks_number))
   allocate(code_new(self%blocks_number))
   do b=1, self%blocks_number
      self%q_work(:,:,:,:,b) = self%q(:,:,:,:,inner_outer_block_map(b))
      coordinates_new(:,b) = self%coordinates(:,inner_outer_block_map(b))
      code_new(b) = self%code(inner_outer_block_map(b))
   enddo
   do b=1, self%blocks_number
      self%q(:,:,:,:,b) = self%q_work(:,:,:,:,b)
      self%coordinates(:,b) = coordinates_new(:,b)
      self%code(b) = code_new(b)
   enddo
   self%inner_blocks_number = inner_blocks_number
   call self%compute_metrics
   endsubroutine blocks_reorder

   subroutine compute_metrics(self)
   !< Compute metrics of each block.
   class(field_object), intent(inout) :: self !< The field.
   integer(I4P)                       :: b    !< Counter.

   do b=1, self%blocks_number
      call self%grid%compute_metrics(coordinates=self%coordinates(:,b),                                         &
                                     emin=self%emin(:,b), emax=self%emax(:,b),                                  &
                                     dx=self%dxyz(1,b), dy=self%dxyz(2,b), dz=self%dxyz(3,b),                   &
                                     x_node=self%x_node(:,b), y_node=self%y_node(:,b), z_node=self%z_node(:,b), &
                                     x_cell=self%x_cell(:,b), y_cell=self%y_cell(:,b), z_cell=self%z_cell(:,b))
   enddo
   endsubroutine compute_metrics

   subroutine destroy(self)
   !< Destroy field.
   class(field_object), intent(inout) :: self  !< The field.
   type(field_object)                 :: fresh !< Fresh field.

   self = fresh
   endsubroutine destroy

   subroutine initialize(self, grid, file_parameters, nv, nb)
   !< Initialize field.
   class(field_object), intent(inout)           :: self            !< The field.
   type(grid_object),   intent(in), target      :: grid            !< Grid data.
   type(file_ini),      intent(inout), optional :: file_parameters !< INI file handler.
   integer(I4P),        intent(in),    optional :: nv              !< Number of field variables.
   integer(I4P),        intent(in),    optional :: nb              !< Number of all blocks that can be stored.

   call self%destroy
   self%grid => grid
   if (present(file_parameters)) call self%load_from_ini_file(file_parameters)

   ! parameters explicitely passed ovveride ones file-passed
   if (present(nv)) self%nv  = nv
   self%block_weight = (self%grid%ngc+self%grid%ni+self%grid%ngc)* &
                       (self%grid%ngc+self%grid%nj+self%grid%ngc)* &
                       (self%grid%ngc+self%grid%nk+self%grid%ngc)*self%nv
   if (present(nb)) self%nb  = nb
   if (self%nb>0) then

      allocate(self%code(self%nb))
      self%code    = -2_I8P
      self%code(1) = -1_I8P ! first block is assumed to be ADAM

      allocate(self%coordinates(4, self%nb))

      allocate(self%emin(3,self%nb))
      allocate(self%emax(3,self%nb))
      allocate(self%dxyz(3,self%nb))
      allocate(self%x_cell(1-self%grid%ngc:self%grid%ni+self%grid%ngc,self%nb))
      allocate(self%y_cell(1-self%grid%ngc:self%grid%nj+self%grid%ngc,self%nb))
      allocate(self%z_cell(1-self%grid%ngc:self%grid%nk+self%grid%ngc,self%nb))
      allocate(self%x_node(0-self%grid%ngc:self%grid%ni+self%grid%ngc,self%nb))
      allocate(self%y_node(0-self%grid%ngc:self%grid%nj+self%grid%ngc,self%nb))
      allocate(self%z_node(0-self%grid%ngc:self%grid%nk+self%grid%ngc,self%nb))
      self%emin(:,1) = self%grid%domain_emin
      self%emax(:,1) = self%grid%domain_emax

      allocate(     self%q(1:self%nv,                                  &
                           1-self%grid%ngc:self%grid%ni+self%grid%ngc, &
                           1-self%grid%ngc:self%grid%nj+self%grid%ngc, &
                           1-self%grid%ngc:self%grid%nk+self%grid%ngc, 1:self%nb))
      allocate(self%q_work(1:self%nv,                                  &
                           1-self%grid%ngc:self%grid%ni+self%grid%ngc, &
                           1-self%grid%ngc:self%grid%nj+self%grid%ngc, &
                           1-self%grid%ngc:self%grid%nk+self%grid%ngc, 1:self%nb))
      self%q = 0._R8P
      self%q_work = 0._R8P
   endif

   ! MPI data
   call MPI_COMM_SIZE(MPI_COMM_WORLD, self%procs_number, self%error)
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, self%error)
   allocate(self%blocks_numbers(0:self%procs_number-1))
   allocate(self%req_send_recv(0:self%procs_number*2-1))
   endsubroutine initialize

   subroutine load_from_ini_file(self, file_parameters)
   !< Load object data from INI file.
   class(field_object), intent(inout) :: self            !< The field.
   type(file_ini),      intent(inout) :: file_parameters !< INI file handler.
   integer(I4P)                       :: buff_I4P        !< I4P buffer.

   call file_parameters%get(section_name='field', option_name='nv', val=buff_I4P) ; self%nv = buff_I4P
   call file_parameters%get(section_name='field', option_name='nb', val=buff_I4P) ; self%nb = buff_I4P
   endsubroutine load_from_ini_file

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
   integer(I8P)                                    :: b               !< Counter.

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

         max_cell_delta = max_cell_delta_dist(distance=distance(0))

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

      function max_cell_delta_dist(distance) result(delta)
      !< Return the maximum cell delta given a comparison distance.
      real(R8P), intent(in) :: distance !< Comparison distance.
      real(R8P)             :: delta    !< Maximum cell delta admissible.

      if (abs(distance) < epsilon(0._R8P)) then
         delta = 0.001_R8P
      else
         delta = huge(0._R8P)
      endif
      endfunction max_cell_delta_dist
   endsubroutine mark_sphere

   subroutine mark_all_blocks(self, mark)
   !< Mark all blocks to be refined, derefined, ecc.
   class(field_object), intent(inout) :: self !< The tree.
   integer(I4P),        intent(in)    :: mark !< Mark to be imposed [TO_BE_REFINED,...].
   integer(I8P)                       :: b    !< Counter.

   if (allocated(self%refinements_needed)) deallocate(self%refinements_needed)
   allocate(self%refinements_needed(self%blocks_number))
   do b=1, self%blocks_number
      self%refinements_needed(b) = mark
   enddo
   endsubroutine mark_all_blocks

   subroutine mpi_gather_refinements_needed(self)
   !< Gather blocks refinement needed status between MPI processes.
   class(field_object), intent(inout) :: self          !< The field.
   integer(I4P), allocatable          :: recv_count(:) !< Number of blocks that are received from process.
   integer(I8P)                       :: p             !< Counter.

   ! computing received blocks
   allocate(recv_count(0:self%procs_number - 1))
   call MPI_ALLGATHER(self%blocks_number, 1_I4P, MPI_INTEGER, &
                      recv_count, 1_I4P, MPI_INTEGER, MPI_COMM_WORLD, self%error)

   ! computing displacement counts
   if (allocated(self%disp_count)) deallocate(self%disp_count)
   allocate(self%disp_count(0:self%procs_number - 1))
   self%disp_count = 0_I4P
   do p=1, self%procs_number - 1
      self%disp_count(p) = self%disp_count(p-1) + recv_count(p-1)
   enddo

   if (allocated(self%refinements_needed_all)) deallocate(self%refinements_needed_all)
   allocate(self%refinements_needed_all(sum(recv_count, dim=1)))
   call MPI_ALLGATHERV(self%refinements_needed, self%blocks_number, MPI_INTEGER, &
                       self%refinements_needed_all, recv_count, self%disp_count, MPI_INTEGER, MPI_COMM_WORLD, self%error)
   endsubroutine mpi_gather_refinements_needed

   subroutine mpi_redistribute(self, comm_map_send, comm_map_recv, comm_map_send_ptr, comm_map_recv_ptr, &
                               local_map, coordinates, code)
   !< Redistribute blocks to processes.
   !< @TODO: Morton codes are not yet redistributed, must be fixed.
   class(field_object),       intent(inout) :: self                   !< The field.
   integer(I8P), allocatable, intent(in)    :: comm_map_send(:)       !< Comm map, blocks to send [sum(comm_map_n_send)].
   integer(I8P), allocatable, intent(in)    :: comm_map_recv(:)       !< Comm map, blocks to receive [sum(comm_map_n_recv)].
   integer(I4P), allocatable, intent(in)    :: comm_map_send_ptr(:)   !< Comm map, pointers in list to send [procs_number+1].
   integer(I4P), allocatable, intent(in)    :: comm_map_recv_ptr(:)   !< Comm map, pointers in list to recv [procs_number+1].
   integer(I8P), allocatable, intent(in)    :: local_map(:,:)         !< Local map, list block index changes of my blocks.
   integer(I4P), allocatable, intent(in)    :: coordinates(:,:)       !< Coordinates of redistributed nodes [nb, ijkl].
   integer(I8P), allocatable, intent(in)    :: code(:)                !< Morton code of redistributed nodes [nb].
   real(R8P),    allocatable                :: send_buffer(:)         !< Send buffer of field cell centered variables.
   real(R8P),    allocatable                :: recv_buffer(:)         !< Recv buffer of field cell centered variables.
   integer(I8P)                             :: send_size, send_offset !< Total size of send buffer.
   integer(I8P)                             :: recv_size, recv_offset !< Total size of recv buffer.
   integer(I4P)                             :: n_keep                 !< Number of keept blocks.
   integer(I4P)                             :: b, bi, p               !< Counter.
   integer(I4P)                             :: ptr_start, ptr_end     !< Counter.
   integer(I4P)                             :: n_recv, n_send         !< Counter.
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
         send_buffer(send_offset:send_offset + self%block_weight - 1) = reshape(self%q(:,:,:,:,bi),[self%block_weight])
         send_offset = send_offset + self%block_weight
      enddo
   endif

   do p=0, self%procs_number - 1_I4P
      ptr_start = comm_map_recv_ptr(p)   * self%block_weight + 1
      ptr_end   = comm_map_recv_ptr(p+1) * self%block_weight
      n_recv    = ptr_end - ptr_start + 1
      if (n_recv > 0) then
         call MPI_IRECV(recv_buffer(ptr_start), n_recv, MPI_REAL8, p, 100, MPI_COMM_WORLD, req_recv(p), self%error)
      endif
   enddo

   do p=0, self%procs_number - 1_I4P
      ptr_start = comm_map_send_ptr(p)   * self%block_weight + 1
      ptr_end   = comm_map_send_ptr(p+1) * self%block_weight
      n_send    = ptr_end - ptr_start + 1
      if (n_send > 0) then
         call MPI_SEND(send_buffer(ptr_start), n_send, MPI_REAL8, p, 100, MPI_COMM_WORLD, self%error)
      endif
   enddo

   call MPI_WAITALL(self%procs_number, req_recv, MPI_STATUSES_IGNORE, self%error)

   if (recv_size > 0_I8P) then
      recv_offset = 1
      do b=1, size(comm_map_recv, dim=1)
          bi = comm_map_recv(b)
          self%q_work(:,:,:,:,bi) = reshape(recv_buffer(recv_offset:recv_offset + self%block_weight -1),&
                                            [self%nv,                                                   &
                                             self%grid%ngc+self%grid%ni+self%grid%ngc,                  &
                                             self%grid%ngc+self%grid%nj+self%grid%ngc,                  &
                                             self%grid%ngc+self%grid%nk+self%grid%ngc])
          recv_offset = recv_offset + self%block_weight
      enddo
   endif

   do b=1, n_keep
      self%q_work(:,:,:,:,local_map(b,1)) = self%q(:,:,:,:,local_map(b,2))
   enddo
   self%blocks_number = n_keep  + recv_size / self%block_weight
   self%q(:,:,:,:,1:self%blocks_number) = self%q_work(:,:,:,:,1:self%blocks_number)
   self%coordinates(:, 1:self%blocks_number) = coordinates
   self%code(1:self%blocks_number) = code
   call self%compute_metrics
   endsubroutine mpi_redistribute

   subroutine prepare_comm_local_ghost(self,                    &
                                       local_map_ghost,         &
                                       comm_map_n_send_ghost  , &
                                       comm_map_n_recv_ghost  , &
                                       comm_map_send_ptr_ghost, &
                                       comm_map_recv_ptr_ghost, &
                                       comm_map_send_ghost    , &
                                       comm_map_recv_ghost)
   !< Prepare communication and local maps for ghosts update and send/receive buffer.
   class(field_object),       intent(inout) :: self                       !< The field.
   integer(I8P), allocatable, intent(in)    :: local_map_ghost(:,:)       !< Local map for ghost cells updating.
   integer(I4P), allocatable, intent(in)    :: comm_map_n_send_ghost(:)   !< Communication map, number of ghost celss to send.
   integer(I4P), allocatable, intent(in)    :: comm_map_n_recv_ghost(:)   !< Communication map, number of ghost celss to recv.
   integer(I4P), allocatable, intent(in)    :: comm_map_send_ptr_ghost(:) !< Communication map, pointers in list to send.
   integer(I4P), allocatable, intent(in)    :: comm_map_recv_ptr_ghost(:) !< Communication map, pointers in list to recv.
   integer(I8P), allocatable, intent(in)    :: comm_map_send_ghost(:,:)   !< Communication map, `fec` information.
   integer(I8P), allocatable, intent(in)    :: comm_map_recv_ghost(:,:)   !< Communication map, `fec` information.

   !OLD---------------------------------------------------------------------------------------------
   !OLD---------------------------------------------------------------------------------------------
   !OLDcall assign_allocatable(lhs=self%local_map_ghost        , rhs=local_map_ghost        )
   !OLDcall assign_allocatable(lhs=self%comm_map_n_send_ghost  , rhs=comm_map_n_send_ghost  )
   !OLDcall assign_allocatable(lhs=self%comm_map_n_recv_ghost  , rhs=comm_map_n_recv_ghost  )
   !OLDcall assign_allocatable(lhs=self%comm_map_send_ptr_ghost, rhs=comm_map_send_ptr_ghost)
   !OLDcall assign_allocatable(lhs=self%comm_map_recv_ptr_ghost, rhs=comm_map_recv_ptr_ghost)
   !OLDcall assign_allocatable(lhs=self%comm_map_send_ghost    , rhs=comm_map_send_ghost    )
   !OLDcall assign_allocatable(lhs=self%comm_map_recv_ghost    , rhs=comm_map_recv_ghost    )
   !OLDif (allocated(self%send_buffer_ghost)) deallocate(self%send_buffer_ghost)
   !OLDif (allocated(self%recv_buffer_ghost)) deallocate(self%recv_buffer_ghost)

   !OLDif (allocated(self%comm_map_send_ghost)) self%comm_map_send_ghost(:,15) = self%comm_map_send_ghost(:,15) * self%nv
   !OLDif (allocated(self%comm_map_recv_ghost)) self%comm_map_recv_ghost(:,15) = self%comm_map_recv_ghost(:,15) * self%nv
   !OLDif (allocated(self%comm_map_send_ptr_ghost)) self%comm_map_send_ptr_ghost = self%comm_map_send_ptr_ghost * self%nv
   !OLDif (allocated(self%comm_map_recv_ptr_ghost)) self%comm_map_recv_ptr_ghost = self%comm_map_recv_ptr_ghost * self%nv
   !OLDif (allocated(self%comm_map_n_send_ghost)) allocate(self%send_buffer_ghost(sum(self%comm_map_n_send_ghost, dim=1)*self%nv))
   !OLDif (allocated(self%comm_map_n_recv_ghost)) allocate(self%recv_buffer_ghost(sum(self%comm_map_n_recv_ghost, dim=1)*self%nv))
   !OLD---------------------------------------------------------------------------------------------
   !OLD---------------------------------------------------------------------------------------------

   call assign_allocatable(lhs=self%local_map_ghost, rhs=local_map_ghost)

   call assign_allocatable(lhs=self%comm_map_n_send_ghost, rhs=comm_map_n_send_ghost  )
   call assign_allocatable(lhs=self%comm_map_n_recv_ghost, rhs=comm_map_n_recv_ghost  )

   ! Nv maps
   call assign_allocatable(lhs=self%comm_map_send_ptr_ghost, rhs=comm_map_send_ptr_ghost)
   call assign_allocatable(lhs=self%comm_map_recv_ptr_ghost, rhs=comm_map_recv_ptr_ghost)
   call assign_allocatable(lhs=self%comm_map_send_ghost    , rhs=comm_map_send_ghost    )
   call assign_allocatable(lhs=self%comm_map_recv_ghost    , rhs=comm_map_recv_ghost    )
   if (allocated(self%send_buffer_ghost)) deallocate(self%send_buffer_ghost)
   !if (allocated(self%send_buffer_ghost))   then
   !    deallocate(self%send_buffer_ghost)
   !    !print*,'deallocate send_buffer_ghost'
   !endif
   if (allocated(self%recv_buffer_ghost)) deallocate(self%recv_buffer_ghost)
   !if (allocated(self%recv_buffer_ghost))   then
   !    deallocate(self%recv_buffer_ghost)
   !    !print*,'deallocate recv_buffer_ghost'
   !endif

   if (allocated(self%comm_map_send_ghost)) self%comm_map_send_ghost(:,15) = self%comm_map_send_ghost(:,15) * self%nv
   if (allocated(self%comm_map_recv_ghost)) self%comm_map_recv_ghost(:,15) = self%comm_map_recv_ghost(:,15) * self%nv
   if (allocated(self%comm_map_send_ptr_ghost)) self%comm_map_send_ptr_ghost = self%comm_map_send_ptr_ghost * self%nv
   if (allocated(self%comm_map_recv_ptr_ghost)) self%comm_map_recv_ptr_ghost = self%comm_map_recv_ptr_ghost * self%nv
   if (allocated(self%comm_map_n_send_ghost)) then
       allocate(self%send_buffer_ghost(sum(self%comm_map_n_send_ghost, dim=1)*self%nv))
       print*,'allocating send_buffer_ghost',size(self%send_buffer_ghost)
       print*,'allocating send_buffer_ghost/2',self%comm_map_n_send_ghost
   endif
   if (allocated(self%comm_map_n_recv_ghost)) then
       allocate(self%recv_buffer_ghost(sum(self%comm_map_n_recv_ghost, dim=1)*self%nv))
       print*,'allocating recv_buffer_ghost',size(self%recv_buffer_ghost)
   endif

   ! single variable maps
   !RIMETTEREcall assign_allocatable(lhs=self%comm_map_send_ptr_ghost_s, rhs=comm_map_send_ptr_ghost)
   !RIMETTEREcall assign_allocatable(lhs=self%comm_map_recv_ptr_ghost_s, rhs=comm_map_recv_ptr_ghost)
   !RIMETTEREcall assign_allocatable(lhs=self%comm_map_send_ghost_s    , rhs=comm_map_send_ghost    )
   !RIMETTEREcall assign_allocatable(lhs=self%comm_map_recv_ghost_s    , rhs=comm_map_recv_ghost    )
   !RIMETTEREif (allocated(self%send_buffer_ghost_s)) deallocate(self%send_buffer_ghost_s)
   !RIMETTEREif (allocated(self%recv_buffer_ghost_s)) deallocate(self%recv_buffer_ghost_s)

   !RIMETTEREif (allocated(self%comm_map_n_send_ghost)) allocate(self%send_buffer_ghost_s(sum(self%comm_map_n_send_ghost, dim=1)))
   !RIMETTEREif (allocated(self%comm_map_n_recv_ghost)) allocate(self%recv_buffer_ghost_s(sum(self%comm_map_n_recv_ghost, dim=1)))
   endsubroutine prepare_comm_local_ghost

   subroutine prepare_local_bc(self, local_map_bc_face, local_map_bc_edge, local_map_bc_corner)
   !< Prepare local maps of boundary conditions.
   class(field_object),       intent(inout) :: self                     !< The field.
   integer(I8P), allocatable, intent(in)    :: local_map_bc_face(:,:)   !< Local map for face BC ghost cells.
   integer(I8P), allocatable, intent(in)    :: local_map_bc_edge(:,:)   !< Local map for edge BC ghost cells.
   integer(I8P), allocatable, intent(in)    :: local_map_bc_corner(:,:) !< Local map for corner BC ghost cells.

   call assign_allocatable(lhs=self%local_map_bc_face  , rhs=local_map_bc_face  )
   call assign_allocatable(lhs=self%local_map_bc_edge  , rhs=local_map_bc_edge  )
   call assign_allocatable(lhs=self%local_map_bc_corner, rhs=local_map_bc_corner)
   endsubroutine prepare_local_bc

   subroutine print_status(self)
   !< Print status of main data.
   class(field_object), intent(in) :: self !< The field.

   print '(A)', 'field status of main data'
   print '(A)', '  field variables number (nv): '//trim(str(self%nv           ))
   print '(A)', '  all blocks number (nb):      '//trim(str(self%nb           ))
   print '(A)', '  blocks number:               '//trim(str(self%blocks_number))
   print '(A)', '  block weight:                '//trim(str(self%block_weight ))
   print '(A)', ''
   endsubroutine print_status

   ! private methods
   subroutine derefine(self, ratio, block_to_derefine, block_derefined)
   !< Derefine blocks.
   !<
   !< Note: blocks number is not updated: mpi redistribute does it. This is dangerous...
   class(field_object),       intent(inout) :: self                 !< The field.
   integer(I4P),              intent(in)    :: ratio                !< Refinement ratio.
   integer(I8P), allocatable, intent(in)    :: block_to_derefine(:) !< List of blocks to be derefined.
   integer(I8P), allocatable, intent(in)    :: block_derefined(:,:) !< List of derefined blocks with Morton code.
   real(R8P)                                :: dx, dy, dz           !< Space deltas.
   integer(I4P)                             :: b, ib                !< Counter.
   integer(I4P)                             :: ic1, ic2, ic3, ic4   !< Counter.
   integer(I4P)                             :: ic5, ic6, ic7, ic8   !< Counter.
   integer(I4P)                             :: iii, jjj, kkk        !< Counter.
   integer(I4P)                             :: i, j, k              !< Counter.

   associate(ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk, q=>self%q, q_work=>self%q_work)
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

         do k=1, nk/2
            do j=1, nj/2
               do i=1, ni/2
                  kkk = (k - 1) * 2 + 1
                  jjj = (j - 1) * 2 + 1
                  iii = (i - 1) * 2 + 1

                  q_work(:,i,     j,     k     ,ib) = (q(:,iii,jjj,  kkk  ,ic1) + q(:,iii+1,jjj,  kkk  ,ic1) + &
                                                       q(:,iii,jjj+1,kkk  ,ic1) + q(:,iii+1,jjj+1,kkk  ,ic1) + &
                                                       q(:,iii,jjj,  kkk+1,ic1) + q(:,iii+1,jjj,  kkk+1,ic1) + &
                                                       q(:,iii,jjj+1,kkk+1,ic1) + q(:,iii+1,jjj+1,kkk+1,ic1)) / 8._R8P

                  q_work(:,i+ni/2,j,     k     ,ib) = (q(:,iii,jjj,  kkk  ,ic2) + q(:,iii+1,jjj,  kkk  ,ic2) + &
                                                       q(:,iii,jjj+1,kkk  ,ic2) + q(:,iii+1,jjj+1,kkk  ,ic2) + &
                                                       q(:,iii,jjj,  kkk+1,ic2) + q(:,iii+1,jjj,  kkk+1,ic2) + &
                                                       q(:,iii,jjj+1,kkk+1,ic2) + q(:,iii+1,jjj+1,kkk+1,ic2)) / 8._R8P

                  q_work(:,i,     j+nj/2,k     ,ib) = (q(:,iii,jjj,  kkk  ,ic3) + q(:,iii+1,jjj,  kkk  ,ic3) + &
                                                       q(:,iii,jjj+1,kkk  ,ic3) + q(:,iii+1,jjj+1,kkk  ,ic3) + &
                                                       q(:,iii,jjj,  kkk+1,ic3) + q(:,iii+1,jjj,  kkk+1,ic3) + &
                                                       q(:,iii,jjj+1,kkk+1,ic3) + q(:,iii+1,jjj+1,kkk+1,ic3)) / 8._R8P

                  q_work(:,i+ni/2,j+nj/2,k     ,ib) = (q(:,iii,jjj,  kkk  ,ic4) + q(:,iii+1,jjj,  kkk  ,ic4) + &
                                                       q(:,iii,jjj+1,kkk  ,ic4) + q(:,iii+1,jjj+1,kkk  ,ic4) + &
                                                       q(:,iii,jjj,  kkk+1,ic4) + q(:,iii+1,jjj,  kkk+1,ic4) + &
                                                       q(:,iii,jjj+1,kkk+1,ic4) + q(:,iii+1,jjj+1,kkk+1,ic4)) / 8._R8P

                  q_work(:,i,     j,     k+nk/2,ib) = (q(:,iii,jjj,  kkk  ,ic5) + q(:,iii+1,jjj,  kkk  ,ic5) + &
                                                       q(:,iii,jjj+1,kkk  ,ic5) + q(:,iii+1,jjj+1,kkk  ,ic5) + &
                                                       q(:,iii,jjj,  kkk+1,ic5) + q(:,iii+1,jjj,  kkk+1,ic5) + &
                                                       q(:,iii,jjj+1,kkk+1,ic5) + q(:,iii+1,jjj+1,kkk+1,ic5)) / 8._R8P

                  q_work(:,i+ni/2,j,     k+nk/2,ib) = (q(:,iii,jjj,  kkk  ,ic6) + q(:,iii+1,jjj,  kkk  ,ic6) + &
                                                       q(:,iii,jjj+1,kkk  ,ic6) + q(:,iii+1,jjj+1,kkk  ,ic6) + &
                                                       q(:,iii,jjj,  kkk+1,ic6) + q(:,iii+1,jjj,  kkk+1,ic6) + &
                                                       q(:,iii,jjj+1,kkk+1,ic6) + q(:,iii+1,jjj+1,kkk+1,ic6)) / 8._R8P

                  q_work(:,i,     j+nj/2,k+nk/2,ib) = (q(:,iii,jjj,  kkk  ,ic7) + q(:,iii+1,jjj,  kkk  ,ic7) + &
                                                       q(:,iii,jjj+1,kkk  ,ic7) + q(:,iii+1,jjj+1,kkk  ,ic7) + &
                                                       q(:,iii,jjj,  kkk+1,ic7) + q(:,iii+1,jjj,  kkk+1,ic7) + &
                                                       q(:,iii,jjj+1,kkk+1,ic7) + q(:,iii+1,jjj+1,kkk+1,ic7)) / 8._R8P

                  q_work(:,i+ni/2,j+nj/2,k+nk/2,ib) = (q(:,iii,jjj,  kkk  ,ic8) + q(:,iii+1,jjj,  kkk  ,ic8) + &
                                                       q(:,iii,jjj+1,kkk  ,ic8) + q(:,iii+1,jjj+1,kkk  ,ic8) + &
                                                       q(:,iii,jjj,  kkk+1,ic8) + q(:,iii+1,jjj,  kkk+1,ic8) + &
                                                       q(:,iii,jjj+1,kkk+1,ic8) + q(:,iii+1,jjj+1,kkk+1,ic8)) / 8._R8P
               enddo
            enddo
         enddo

         q(:,1:ni,1:nj,1:nk,ib) = q_work(:,1:ni,1:nj,1:nk,ib)

         self%code(ib) = block_derefined(1,b)
      enddo
   endif
   endassociate
   endsubroutine derefine

   subroutine refine(self, ratio, block_to_refine, block_refined)
   !< Refine blocks.
   !<
   !< Note: blocks number is not updated: mpi redistribute does it. This is dangerous...
   class(field_object),       intent(inout) :: self                      !< The field.
   integer(I4P),              intent(in)    :: ratio                     !< Refinement ratio.
   integer(I8P), allocatable, intent(in)    :: block_to_refine(:,:)      !< List of blocks to be refined.
   integer(I8P), allocatable, intent(in)    :: block_refined(:,:)        !< List of refined blocks with Morton code.
   real(R8P)                                :: dx, dy, dz                !< Space deltas.
   integer(I4P)                             :: b, i, j, k                !< Spatial counter.
   integer(I4P)                             :: ib, ic, ii, ic_local      !< Counter.
   integer(I4P)                             :: i_fine, j_fine, k_fine    !< Counter.
   integer(I4P)                             :: i_delta, j_delta, k_delta !< Counter.
   integer(I4P)                             :: ic1, ic2, ic3, ic4        !< Counter.
   integer(I4P)                             :: ic5, ic6, ic7, ic8        !< Counter.

   associate(ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk, q=>self%q, q_work=>self%q_work)
   if (allocated(block_to_refine)) then
      do b=1, size(block_to_refine, dim=2)
         if (self%myrank /= block_to_refine(2,b)) cycle
         ib = block_to_refine(1,b)

         q_work(:,:,:,:,ib) = q(:,:,:,:,ib)

         do ic_local=1, 8
            ic = block_refined(2,(b-1)*ratio+ic_local)
            ic1 = mod(ic_local - 1, 2)
            ic2 = mod((ic_local - 1)/2, 2)
            ic3 = mod((ic_local - 1)/4, 2)
            do k=1+nk/2*ic3,nk/2+nk/2*ic3
               do j=1+nj/2*ic2,nj/2+nj/2*ic2
                  do i=1+ni/2*ic1,ni/2+ni/2*ic1
                     k_fine = mod(k - 1, nk/2) * 2 + 1
                     j_fine = mod(j - 1, nj/2) * 2 + 1
                     i_fine = mod(i - 1, ni/2) * 2 + 1
                     q(:,i_fine:i_fine+1,j_fine:j_fine+1,k_fine:k_fine+1,ic) = 0._R8P
                     do k_delta=0,1
                     do j_delta=0,1
                     do i_delta=0,1
                     q(:,i_fine,  j_fine,  k_fine,  ic) = q(:,i_fine,j_fine,k_fine,ic) +   &
                                                          (0.25_R8P + i_delta * 0.5_R8P) * &
                                                          (0.25_R8P + j_delta * 0.5_R8P) * &
                                                          (0.25_R8P + k_delta * 0.5_R8P) * &
                                                          q_work(:,i+i_delta-1, j+j_delta-1, k+k_delta-1,ib)
                     q(:,i_fine+1,j_fine,  k_fine,  ic) = q(:,i_fine+1,j_fine,k_fine,ic) + &
                                                          (0.75_R8P - i_delta * 0.5_R8P) * &
                                                          (0.25_R8P + j_delta * 0.5_R8P) * &
                                                          (0.25_R8P + k_delta * 0.5_R8P) * &
                                                          q_work(:,i+i_delta,   j+j_delta-1, k+k_delta-1,ib)
                     q(:,i_fine,  j_fine+1,k_fine,  ic) = q(:,i_fine,j_fine+1,k_fine,ic) + &
                                                          (0.25_R8P + i_delta * 0.5_R8P) * &
                                                          (0.75_R8P - j_delta * 0.5_R8P) * &
                                                          (0.25_R8P + k_delta * 0.5_R8P) * &
                                                          q_work(:,i+i_delta-1, j+j_delta  , k+k_delta-1,ib)
                     q(:,i_fine+1,j_fine+1,k_fine,  ic) = q(:,i_fine+1,j_fine+1,k_fine,ic) + &
                                                          (0.75_R8P - i_delta * 0.5_R8P) *   &
                                                          (0.75_R8P - j_delta * 0.5_R8P) *   &
                                                          (0.25_R8P + k_delta * 0.5_R8P) *   &
                                                          q_work(:,i+i_delta,   j+j_delta  , k+k_delta-1,ib)
                     q(:,i_fine,  j_fine,  k_fine+1,ic) = q(:,i_fine,j_fine,k_fine+1,ic) + &
                                                          (0.25_R8P + i_delta * 0.5_R8P) * &
                                                          (0.25_R8P + j_delta * 0.5_R8P) * &
                                                          (0.75_R8P - k_delta * 0.5_R8P) * &
                                                          q_work(:,i+i_delta-1, j+j_delta-1, k+k_delta  ,ib)
                     q(:,i_fine+1,j_fine,  k_fine+1,ic) = q(:,i_fine+1,j_fine,k_fine+1,ic) + &
                                                          (0.75_R8P - i_delta * 0.5_R8P) *   &
                                                          (0.25_R8P + j_delta * 0.5_R8P) *   &
                                                          (0.75_R8P - k_delta * 0.5_R8P) *   &
                                                          q_work(:,i+i_delta,   j+j_delta-1, k+k_delta  ,ib)
                     q(:,i_fine,  j_fine+1,k_fine+1,ic) = q(:,i_fine,j_fine+1,k_fine+1,ic) + &
                                                          (0.25_R8P + i_delta * 0.5_R8P) *   &
                                                          (0.75_R8P - j_delta * 0.5_R8P) *   &
                                                          (0.75_R8P - k_delta * 0.5_R8P) *   &
                                                          q_work(:,i+i_delta-1, j+j_delta  , k+k_delta  ,ib)
                     q(:,i_fine+1,j_fine+1,k_fine+1,ic) = q(:,i_fine+1,j_fine+1,k_fine+1,ic) + &
                                                          (0.75_R8P - i_delta * 0.5_R8P) *     &
                                                          (0.75_R8P - j_delta * 0.5_R8P) *     &
                                                          (0.75_R8P - k_delta * 0.5_R8P) *     &
                                                          q_work(:,i+i_delta,   j+j_delta  , k+k_delta  ,ib)
                     enddo
                     enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo

         ic1 = block_refined(2,(b-1)*ratio+1)
         ic2 = block_refined(2,(b-1)*ratio+2)
         ic3 = block_refined(2,(b-1)*ratio+3)
         ic4 = block_refined(2,(b-1)*ratio+4)
         ic5 = block_refined(2,(b-1)*ratio+5)
         ic6 = block_refined(2,(b-1)*ratio+6)
         ic7 = block_refined(2,(b-1)*ratio+7)
         ic8 = block_refined(2,(b-1)*ratio+8)
         self%code(ic1) = block_refined(1,(b-1)*ratio+1)
         self%code(ic2) = block_refined(1,(b-1)*ratio+2)
         self%code(ic3) = block_refined(1,(b-1)*ratio+3)
         self%code(ic4) = block_refined(1,(b-1)*ratio+4)
         self%code(ic5) = block_refined(1,(b-1)*ratio+5)
         self%code(ic6) = block_refined(1,(b-1)*ratio+6)
         self%code(ic7) = block_refined(1,(b-1)*ratio+7)
         self%code(ic8) = block_refined(1,(b-1)*ratio+8)
      enddo
   endif
   endassociate
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
   call assign_allocatable(lhs=lhs%code, rhs=rhs%code)
   call assign_allocatable(lhs=lhs%coordinates, rhs=rhs%coordinates)
   call assign_allocatable(lhs=lhs%emin, rhs=rhs%emin)
   call assign_allocatable(lhs=lhs%emax, rhs=rhs%emax)
   call assign_allocatable(lhs=lhs%dxyz, rhs=rhs%dxyz)
   call assign_allocatable(lhs=lhs%x_cell, rhs=rhs%x_cell)
   call assign_allocatable(lhs=lhs%y_cell, rhs=rhs%y_cell)
   call assign_allocatable(lhs=lhs%z_cell, rhs=rhs%z_cell)
   call assign_allocatable(lhs=lhs%x_node, rhs=rhs%x_node)
   call assign_allocatable(lhs=lhs%y_node, rhs=rhs%y_node)
   call assign_allocatable(lhs=lhs%z_node, rhs=rhs%z_node)
   call assign_allocatable(lhs=lhs%local_map_ghost, rhs=rhs%local_map_ghost)
   ! MPI data, unrelated to field equations
   lhs%error        = rhs%error
   lhs%myrank       = rhs%myrank
   lhs%procs_number = rhs%procs_number
   call assign_allocatable(lhs=lhs%blocks_numbers, rhs=rhs%blocks_numbers)
   call assign_allocatable(lhs=lhs%refinements_needed, rhs=rhs%refinements_needed)
   call assign_allocatable(lhs=lhs%refinements_needed_all, rhs=rhs%refinements_needed_all)
   call assign_allocatable(lhs=lhs%disp_count, rhs=rhs%disp_count)
   lhs%inner_blocks_number = rhs%inner_blocks_number
   call assign_allocatable(lhs=lhs%local_map_bc_face, rhs=rhs%local_map_bc_face)
   call assign_allocatable(lhs=lhs%local_map_bc_edge, rhs=rhs%local_map_bc_edge)
   call assign_allocatable(lhs=lhs%local_map_bc_corner, rhs=rhs%local_map_bc_corner)
   call assign_allocatable(lhs=lhs%req_send_recv, rhs=rhs%req_send_recv)
   call assign_allocatable(lhs=lhs%comm_map_n_send_ghost, rhs=rhs%comm_map_n_send_ghost)
   call assign_allocatable(lhs=lhs%comm_map_n_recv_ghost, rhs=rhs%comm_map_n_recv_ghost)
   call assign_allocatable(lhs=lhs%comm_map_send_ptr_ghost, rhs=rhs%comm_map_send_ptr_ghost)
   call assign_allocatable(lhs=lhs%comm_map_recv_ptr_ghost, rhs=rhs%comm_map_recv_ptr_ghost)
   call assign_allocatable(lhs=lhs%comm_map_send_ghost, rhs=rhs%comm_map_send_ghost)
   call assign_allocatable(lhs=lhs%comm_map_recv_ghost, rhs=rhs%comm_map_recv_ghost)
   ! MPI data, related to field equations
   call assign_allocatable(lhs=lhs%send_buffer_ghost, rhs=rhs%send_buffer_ghost)
   call assign_allocatable(lhs=lhs%recv_buffer_ghost, rhs=rhs%recv_buffer_ghost)
   ! field equations data
   call assign_allocatable(lhs=lhs%q, rhs=rhs%q)
   call assign_allocatable(lhs=lhs%q_work, rhs=rhs%q_work)
   endsubroutine field_assign_field
endmodule adam_field_object
