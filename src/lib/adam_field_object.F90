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
   ! mesh related data, unrelated to field equations
   type(grid_object), pointer :: grid=>null()         !< Grid data.
   integer(I4P)               :: nv=1_I4P             !< Number of field variables.
   integer(I4P)               :: block_weight=0_I4P   !< Block weight, `cells_number * variables_number`.
   integer(I4P)               :: nb=0_I4P             !< Number of all blocks that can be stored.
   integer(I4P)               :: blocks_number=0_I4P  !< Number of blocks actually stored.
   integer(I8P), allocatable  :: code(:)              !< Morton codes [nb].
   integer(I4P), allocatable  :: coordinates(:,:)     !< Coordinates IJKL for each block [nb,4].
   real(R8P),    allocatable  :: emin(:,:)            !< Coordinates of minimum abscissa of each block [3,nb].
   real(R8P),    allocatable  :: emax(:,:)            !< Coordinates of maximum abscissa of each block [3,nb].
   real(R8P),    allocatable  :: x_node(:,:)          !< X coordinates of [0-gc1:ni+gc2,nb].
   real(R8P),    allocatable  :: y_node(:,:)          !< Y coordinates of [0-gc3:nj+gc4,nb].
   real(R8P),    allocatable  :: z_node(:,:)          !< Z coordinates of [0-gc5:nk+gc6,nb].
   real(R8P),    allocatable  :: x_cell(:,:)          !< X coordinates of [1-gc1:ni+gc2,nb].
   real(R8P),    allocatable  :: y_cell(:,:)          !< Y coordinates of [1-gc3:nj+gc4,nb].
   real(R8P),    allocatable  :: z_cell(:,:)          !< Z coordinates of [1-gc5:nk+gc6,nb].
   integer(I8P), allocatable  :: local_map_ghost(:,:) !< Local map for ghost cells updating.
   ! MPI data, unrelated to field equations
   integer(I4P)              :: error                      !< Error traping flag.
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
   integer(I8P), allocatable :: comm_map_send_ghost(:,:)   !< Communication map, `fec` information [fec_number, 5].
   integer(I8P), allocatable :: comm_map_recv_ghost(:,:)   !< Communication map, `fec` information [fec_number, 5].
   ! MPI data, related to field equations
   real(R8P), allocatable    :: send_buffer_ghost(:) !< Send buffer of ghost cells.
   real(R8P), allocatable    :: recv_buffer_ghost(:) !< Receive buffer of ghost cells.
   ! RK data, related to field equations
   real(R8P) :: alph(3,3) = reshape([0._R8P, 1._R8P, 0.25_R8P, &
                                     0._R8P, 0._R8P, 0.25_R8P, &
                                     0._R8P, 0._R8P,0._R8P], [3,3]) !< RK alpha coefficients.
   real(R8P) :: beta(3) = [1._R8P/6._R8P, &
                           1._R8P/6._R8P, &
                           2._R8P/3._R8P]                           !< RK beta coefficients.
   real(R8P) :: gamm(3) = [0._R8P, &
                           1._R8P, &
                           0._R8P]                                  !< RK gamma coefficients.
   ! field equations data
   real(R8P), allocatable  :: u(     :,:,:,:  ) !< Field cell centered variables [ni+gc12,nj+gc34,nk+gc56,nv,nb].
   real(R8P), allocatable  :: u_work(:,:,:,:  ) !< Field cell centered variables, buffer memory.
   real(R8P), allocatable  :: u_s(   :,:,:,:,:) !< RK field stages.
   contains
      ! public methods
      procedure, pass(self) :: adapt                         !< Adapt field accordingly to refine/derefine necessity.
      procedure, pass(self) :: blocks_reorder                !< Reorder blocks indexes in field.
      procedure, pass(self) :: compute_metrics               !< Compute metrics of each block.
      procedure, pass(self) :: destroy                       !< Destroy the field.
      procedure, pass(self) :: initialize                    !< Initialize the field.
      procedure, pass(self) :: mark_by_u_value               !< Mark blocks to be refined/derefined by a `u` value.
      procedure, pass(self) :: mark_all_blocks               !< Mark all blocks to be refined, derefined, ecc.
      procedure, pass(self) :: mark_sphere                   !< Mark blocks to be refined/derefined by sphere distance.
      procedure, pass(self) :: max_cell_delta                !< Return the maximum cell delta given a comparison distance.
      procedure, pass(self) :: mpi_gather_refinements_needed !< Gather blocks refinement needed status between MPI processes.
      procedure, pass(self) :: mpi_redistribute              !< Redistribute blocks to processes.
      procedure, pass(self) :: prepare_comm_local_ghost      !< Prepare communication and local maps/buffers for ghosts update.
      procedure, pass(self) :: rk_integrate                  !< Runge Kutta integration of field.
      procedure, pass(self) :: compute_residuals             !< Compute residuals of field.
      procedure, pass(self) :: set_initial_conditions        !< Set initial conditions of field.
      procedure, pass(self) :: update_ghost                  !< Update ghost cells.
      procedure, pass(self) :: update_ghost_mpi              !< Update ghost cells within other processes.
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
      self%u_work(:,:,:,b) = self%u(:,:,:,inner_outer_block_map(b))
      coordinates_new(:,b) = self%coordinates(:,inner_outer_block_map(b))
      code_new(b) = self%code(inner_outer_block_map(b))
   enddo
   do b=1, self%blocks_number
      self%u(:,:,:,b) = self%u_work(:,:,:,b)
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

   subroutine initialize(self, grid, nv, nb)
   !< Initialize field.
   class(field_object), intent(inout)        :: self    !< The field.
   type(grid_object),   intent(in), target   :: grid    !< Grid data.
   integer(I4P),        intent(in), optional :: nv      !< Number of field variables.
   integer(I4P),        intent(in), optional :: nb      !< Number of all blocks that can be stored.

   call self%destroy
   self%grid => grid
   if (present(nv)) self%nv  = nv
   self%block_weight = (self%grid%gc1+self%grid%ni+self%grid%gc2)* &
                       (self%grid%gc3+self%grid%nj+self%grid%gc4)* &
                       (self%grid%gc5+self%grid%nk+self%grid%gc6)*self%nv
   if (present(nb)) self%nb  = nb
   if (self%nb>0) then

      allocate(self%code(self%nb))
      self%code    = -2_I8P
      self%code(1) = -1_I8P ! first block is assumed to be ADAM

      allocate(self%coordinates(4, self%nb))

      allocate(self%emin(3,self%nb))
      allocate(self%emax(3,self%nb))
      allocate(self%x_cell(1-self%grid%gc1:self%grid%ni+self%grid%gc2,self%nb))
      allocate(self%y_cell(1-self%grid%gc3:self%grid%nj+self%grid%gc4,self%nb))
      allocate(self%z_cell(1-self%grid%gc5:self%grid%nk+self%grid%gc6,self%nb))
      allocate(self%x_node(0-self%grid%gc1:self%grid%ni+self%grid%gc2,self%nb))
      allocate(self%y_node(0-self%grid%gc3:self%grid%nj+self%grid%gc4,self%nb))
      allocate(self%z_node(0-self%grid%gc5:self%grid%nk+self%grid%gc6,self%nb))
      self%emin(:,1) = self%grid%domain_emin
      self%emax(:,1) = self%grid%domain_emax

      allocate(self%u(1-self%grid%gc1:self%grid%ni+self%grid%gc2, &
                      1-self%grid%gc3:self%grid%nj+self%grid%gc4, &
                      1-self%grid%gc5:self%grid%nk+self%grid%gc6, 1:self%nb))
      allocate(self%u_work(1-self%grid%gc1:self%grid%ni+self%grid%gc2, &
                           1-self%grid%gc3:self%grid%nj+self%grid%gc4, &
                           1-self%grid%gc5:self%grid%nk+self%grid%gc6, 1:self%nb))
      allocate(self%u_s(1-self%grid%gc1:self%grid%ni+self%grid%gc2, &
                        1-self%grid%gc3:self%grid%nj+self%grid%gc4, &
                        1-self%grid%gc5:self%grid%nk+self%grid%gc6, 1:self%nb, 1:3))
      self%u = 0._R8P
      self%u_work = 0._R8P
   endif
   ! MPI data
#ifdef _MPI_
   call MPI_COMM_SIZE(MPI_COMM_WORLD, self%procs_number, self%error)
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, self%error)
   allocate(self%blocks_numbers(0:self%procs_number-1))
   allocate(self%req_send_recv(0:self%procs_number*2-1))
#endif
   endsubroutine initialize

   subroutine mark_by_u_value(self, u_value, threshold)
   !< Mark blocks to be refined/derefined by a `u` value.
   class(field_object), intent(inout)        :: self       !< The field.
   real(R8P),           intent(in)           :: u_value    !< `u` value to be tracked.
   real(R8P),           intent(in), optional :: threshold  !< Threshold for sphere proximity.
   real(R8P)                                 :: threshold_ !< Threshold for sphere proximity, local var.
   real(R8P)                                 :: u_mean     !< `u` mean value of block.
   integer(I4P)                              :: b          !< Counter.

   threshold_ = 1.5_R8P ; if (present(threshold)) threshold_ = threshold
   if (allocated(self%refinements_needed)) deallocate(self%refinements_needed)
   allocate(self%refinements_needed(self%blocks_number))
   do b=1, self%blocks_number

      associate (ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk)
         u_mean = sum(self%u(1:ni,1:nj,1:nk,b)) / (ni*nj*nk)
         ! u_mean = maxval(self%u(1:ni,1:nj,1:nk,b))
      endassociate

      if (abs(u_mean - u_value) < 0.05) then
      ! if ((u_mean - u_value) > 0.47) then
         self%refinements_needed(b) = TO_BE_REFINED
      elseif (abs(u_mean - u_value) * threshold_ > 0.05) then
      ! elseif ((u_mean - u_value) * threshold_ < 0.47) then
         self%refinements_needed(b) = TO_BE_DEREFINED
      else
         self%refinements_needed(b) = TO_NOT_TOUCH
      endif
   enddo
   endsubroutine mark_by_u_value

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

   function max_cell_delta(self, distance) result(delta)
   !< Return the maximum cell delta given a comparison distance.
   class(field_object), intent(in) :: self     !< The field.
   real(R8P),           intent(in) :: distance !< Comparison distance.
   real(R8P)                       :: delta    !< Maximum cell delta admissible.

   if (abs(distance) < epsilon(0._R8P)) then
      ! delta = 0.001_R8P
      delta = 0.001_R8P
   else
      delta = huge(0._R8P)
   endif
   endfunction max_cell_delta

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
#ifdef _MPI_
   call MPI_ALLGATHERV(self%refinements_needed, self%blocks_number, MPI_INTEGER, &
                       self%refinements_needed_all, recv_count, self%disp_count, MPI_INTEGER, MPI_COMM_WORLD, self%error)
#endif
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
#ifdef _MPI_
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
         call MPI_IRECV(recv_buffer(ptr_start), n_recv, MPI_REAL8, p, 100, MPI_COMM_WORLD, req_recv(p), self%error)
#endif
      endif
   enddo

   do p=0, self%procs_number - 1_I4P
      ptr_start = comm_map_send_ptr(p)   * self%block_weight + 1
      ptr_end   = comm_map_send_ptr(p+1) * self%block_weight
      n_send    = ptr_end - ptr_start + 1
      if (n_send > 0) then
#ifdef _MPI_
         call MPI_SEND(send_buffer(ptr_start), n_send, MPI_REAL8, p, 100, MPI_COMM_WORLD, self%error)
#endif
      endif
   enddo

#ifdef _MPI_
   call MPI_WAITALL(self%procs_number, req_recv, MPI_STATUSES_IGNORE, self%error)
#endif

   if (recv_size > 0_I8P) then
      recv_offset = 1
      do b=1, size(comm_map_recv, dim=1)
          bi = comm_map_recv(b)
          self%u_work(:,:,:,bi) = reshape(recv_buffer(recv_offset:recv_offset + self%block_weight -1),&
                                          [self%grid%gc1+self%grid%ni+self%grid%gc2,                  &
                                           self%grid%gc3+self%grid%nj+self%grid%gc4,                  &
                                           self%grid%gc5+self%grid%nk+self%grid%gc6])
          recv_offset = recv_offset + self%block_weight
      enddo
   endif

   do b=1, n_keep
      self%u_work(:,:,:,local_map(b,1)) = self%u(:,:,:,local_map(b,2))
   enddo
   self%blocks_number = n_keep  + recv_size / self%block_weight
   self%u(:,:,:,1:self%blocks_number) = self%u_work(:,:,:,1:self%blocks_number)
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
   integer(I4P)                             :: n_recv, n_send             !< Counter.

   if (allocated(self%local_map_ghost        )) deallocate(self%local_map_ghost        )
   if (allocated(self%comm_map_n_send_ghost  )) deallocate(self%comm_map_n_send_ghost  )
   if (allocated(self%comm_map_n_recv_ghost  )) deallocate(self%comm_map_n_recv_ghost  )
   if (allocated(self%comm_map_send_ptr_ghost)) deallocate(self%comm_map_send_ptr_ghost)
   if (allocated(self%comm_map_recv_ptr_ghost)) deallocate(self%comm_map_recv_ptr_ghost)
   if (allocated(self%comm_map_send_ghost    )) deallocate(self%comm_map_send_ghost    )
   if (allocated(self%comm_map_recv_ghost    )) deallocate(self%comm_map_recv_ghost    )
   if (allocated(self%send_buffer_ghost      )) deallocate(self%send_buffer_ghost      )
   if (allocated(self%recv_buffer_ghost      )) deallocate(self%recv_buffer_ghost      )

   if (allocated(local_map_ghost        )) self%local_map_ghost         = local_map_ghost
   if (allocated(comm_map_n_send_ghost  )) self%comm_map_n_send_ghost   = comm_map_n_send_ghost
   if (allocated(comm_map_n_recv_ghost  )) self%comm_map_n_recv_ghost   = comm_map_n_recv_ghost
   if (allocated(comm_map_send_ptr_ghost)) self%comm_map_send_ptr_ghost = comm_map_send_ptr_ghost
   if (allocated(comm_map_recv_ptr_ghost)) self%comm_map_recv_ptr_ghost = comm_map_recv_ptr_ghost
   if (allocated(comm_map_send_ghost    )) self%comm_map_send_ghost     = comm_map_send_ghost
   if (allocated(comm_map_recv_ghost    )) self%comm_map_recv_ghost     = comm_map_recv_ghost

   n_send = 0_I8P ; if (allocated(self%comm_map_n_send_ghost)) n_send = sum(self%comm_map_n_send_ghost, dim=1)
   n_recv = 0_I8P ; if (allocated(self%comm_map_n_recv_ghost)) n_recv = sum(self%comm_map_n_recv_ghost, dim=1)
   if (n_send > 0_I8P) allocate(self%send_buffer_ghost(n_send))
   if (n_recv > 0_I8P) allocate(self%recv_buffer_ghost(n_recv))
   endsubroutine prepare_comm_local_ghost

   subroutine rk_integrate(self, t, Dt, do_ghost_syncro, residual)
   !< Runge Kutta integration of field.
   class(field_object), intent(inout)         :: self             !< The field.
   real(R8P),           intent(in)            :: t                !< Time.
   real(R8P),           intent(in)            :: Dt               !< Time step.
   logical,             intent(in),  optional :: do_ghost_syncro  !< Flag to do syncrous ghost update.
   real(R8P),           intent(out), optional :: residual         !< Global residual.
   logical                                    :: do_ghost_syncro_ !< Flag to do syncrous ghost update, local var.
   integer(I4P)                               :: b, s, ss         !< Counter.

   do_ghost_syncro_ = .true. ; if (present(do_ghost_syncro)) do_ghost_syncro_ = do_ghost_syncro
   associate(alph=>self%alph, beta=>self%beta, gamm=>self%gamm,    &
             ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk, &
             blocks_number=>self%blocks_number,                    &
             inner_blocks_number=>self%inner_blocks_number,        &
             u=>self%u, u_s=>self%u_s)
   do s=1, 3
      u_s(1:ni,1:nj,1:nk,1:blocks_number,s) = u(1:ni,1:nj,1:nk,1:blocks_number)
      do ss=1, s - 1
         u_s(1:ni,1:nj,1:nk,1:blocks_number,s) = u_s(1:ni,1:nj,1:nk,1:blocks_number,s ) + &
                                                (u_s(1:ni,1:nj,1:nk,1:blocks_number,ss) * (Dt * alph(s, ss)))
      enddo
      if (do_ghost_syncro_) then
         call self%update_ghost(s=s)
         call self%update_ghost_mpi(s=s)
         call self%compute_residuals(s=s, t=t + gamm(s) * Dt, &
                                     block_start=1, block_end=blocks_number)
      else
         call self%update_ghost(s=s)
         call self%update_ghost_mpi(s=s, step=1)
         call self%update_ghost_mpi(s=s, step=2)
         call self%compute_residuals(s=s, t=t + gamm(s) * Dt, &
                                     block_start=1, block_end=inner_blocks_number)
         call self%update_ghost_mpi(s=s, step=3)
         call self%compute_residuals(s=s, t=t + gamm(s) * Dt, &
                                     block_start=inner_blocks_number+1, block_end=blocks_number)
      endif
      if (present(residual).and.s==3) then
         residual = 0._R8P
         do b=1, blocks_number
            residual = residual + sum(u_s(1:ni,1:nj,1:nk,b,s))/ni/nj/nk
         enddo
#ifdef _MPI_
         call MPI_ALLREDUCE(MPI_IN_PLACE, residual, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, self%error)
#endif
      endif
   enddo
   do s=1, 3
      u(1:ni,1:nj,1:nk,1:blocks_number) = u(  1:ni,1:nj,1:nk,1:blocks_number) + &
                                          u_s(1:ni,1:nj,1:nk,1:blocks_number,s) * Dt * beta(s)
   enddo
   endassociate
   endsubroutine rk_integrate

   subroutine compute_residuals(self, s, t, block_start, block_end)
   !< Compute residuals of field.
   class(field_object), intent(inout) :: self        !< The field.
   integer(I4P),        intent(in)    :: s           !< Current stage.
   real(R8P),           intent(in)    :: t           !< Time.
   integer(I4P),        intent(in)    :: block_start !< Index of block to start residuals computation.
   integer(I4P),        intent(in)    :: block_end   !< Index of block to end   residuals computation.
   integer(I4P)                       :: b, i, j, k  !< Counter.

   do b=block_start, block_end
      do k=1, self%grid%nk
         do j=1, self%grid%nj
            do i=1, self%grid%ni
               self%u_work(i,j,k,b) = self%u_s(i+1,j,  k,  b,s) + self%u_s(i-1,j,  k,  b,s) + &
                                      self%u_s(i,  j+1,k,  b,s) + self%u_s(i,  j-1,k,  b,s) + &
                                      self%u_s(i,  j,  k+1,b,s) + self%u_s(i,  j,  k-1,b,s) - &
                             6._R8P * self%u_s(i,  j,  k,  b,s)
            enddo
         enddo
      enddo
   enddo
   self%u_s(:,:,:,block_start:block_end,s) = self%u_work(:,:,:,block_start:block_end)
   endsubroutine compute_residuals

   subroutine set_initial_conditions(self)
   !< Set initial conditions of field.
   class(field_object), intent(inout) :: self    !< The field.
   integer(I4P)                       :: b       !< Counter.
   integer(I4P)                       :: i, j, k !< Counter.
   real(R8P)                          :: a       !< Gaussian amplitude.
   real(R8P)                          :: sigma_x !< Gaussian x variance.
   real(R8P)                          :: sigma_y !< Gaussian y variance.
   real(R8P)                          :: sigma_z !< Gaussian z variance.
   real(R8P)                          :: x_0     !< Gaussian x center.
   real(R8P)                          :: y_0     !< Gaussian y center.
   real(R8P)                          :: z_0     !< Gaussian z center.

   a = 1.0_R8P
   x_0 = (self%grid%domain_emax(1) - self%grid%domain_emin(1)) / 2.0_R8P
   y_0 = (self%grid%domain_emax(2) - self%grid%domain_emin(2)) / 2.0_R8P
   z_0 = (self%grid%domain_emax(3) - self%grid%domain_emin(3)) / 2.0_R8P
   sigma_x = 0.2_R8P
   sigma_y = 0.2_R8P
   sigma_z = 0.2_R8P
   associate(blocks_number=>self%blocks_number,                             &
             u=>self%u,                                                     &
             ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk,          &
             gc1=>self%grid%gc1, gc2=>self%grid%gc2, gc3=>self%grid%gc3,    &
             gc4=>self%grid%gc4, gc5=>self%grid%gc5, gc6=>self%grid%gc6,    &
             x_cell=>self%x_cell, y_cell=>self%y_cell, z_cell=>self%z_cell)
   do b=1, blocks_number
      do k=1-gc5, nk+gc6
         do j=1-gc3, nj+gc4
            do i=1-gc1, ni+gc2
               u(i,j,k,b) = a * exp(-((x_cell(i,b) - x_0)**2/(2 * sigma_x**2)+&
                                      (y_cell(j,b) - y_0)**2/(2 * sigma_y**2)+&
                                      (z_cell(k,b) - z_0)**2/(2 * sigma_z**2)))
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine set_initial_conditions

   subroutine update_ghost(self, s)
   !< Update ghost cells.
   class(field_object), intent(inout) :: self          !< The field.
   integer(I4P),        intent(in)    :: s             !< Stage.
   integer(I4P)                       :: i, j, k, mf   !< Counter.
   integer(I4P)                       :: iii, jjj, kkk !< Counter.
   integer(I4P)                       :: fec           !< Direction where ghost cells are updated, faces/edges/corners.
   integer(I4P)                       :: portion       !< Portion of fec updated (0=>whole fec).
   integer(I4P)                       :: b_recv        !< Index of receiving block.
   integer(I4P)                       :: b_send        !< Index of sending block.
   integer(I4P)                       :: ijkmin(3)     !< Lower limit of ijk indexes.
   integer(I4P)                       :: ijkmax(3)     !< Upper limit of ijk indexes.
   integer(I4P)                       :: ijkdelta(3)   !< Delta offset for ghost-inner cells mapping same refinement.

   if (.not.allocated(self%local_map_ghost)) return
   associate(u_s=>self%u_s)
   do mf=1, size(self%local_map_ghost, dim=1)
      b_recv   = self%local_map_ghost(mf, 1)
      b_send   = self%local_map_ghost(mf, 2)
      fec      = self%local_map_ghost(mf, 3)
      portion  = self%local_map_ghost(mf, 4)
      ijkmin   = self%local_map_ghost(mf, 5:7)
      ijkmax   = self%local_map_ghost(mf, 8:10)
      ijkdelta = self%local_map_ghost(mf, 11:13)
      if     (portion==0) then
         ! receiving from a block with the same refinement
         do k=ijkmin(3), ijkmax(3)
            do j=ijkmin(2), ijkmax(2)
               do i=ijkmin(1), ijkmax(1)
                  u_s(i,j,k,b_recv,s) = u_s(i+ijkdelta(1),j+ijkdelta(2),k+ijkdelta(3),b_send,s)
               enddo
            enddo
         enddo
      elseif (portion>0) then
         ! receiving from a block finer than me
         do k=ijkmin(3), ijkmax(3)
            do j=ijkmin(2), ijkmax(2)
               do i=ijkmin(1), ijkmax(1)
                  kkk = 2 * k + ijkdelta(3)
                  jjj = 2 * j + ijkdelta(2)
                  iii = 2 * i + ijkdelta(1)
                  u_s(i,j,k,b_recv,s) = (u_s(iii,jjj,  kkk,  b_send,s) + u_s(iii+1,jjj,  kkk,  b_send,s) + &
                                         u_s(iii,jjj+1,kkk,  b_send,s) + u_s(iii+1,jjj+1,kkk,  b_send,s) + &
                                         u_s(iii,jjj,  kkk+1,b_send,s) + u_s(iii+1,jjj,  kkk+1,b_send,s) + &
                                         u_s(iii,jjj+1,kkk+1,b_send,s) + u_s(iii+1,jjj+1,kkk+1,b_send,s)) / 8._R8P
               enddo
            enddo
         enddo
      else
         ! receiving from a block coarser than me
         do k=ijkmin(3), ijkmax(3)
            do j=ijkmin(2), ijkmax(2)
               do i=ijkmin(1), ijkmax(1)
                  kkk = 2 * k + ijkdelta(3)
                  jjj = 2 * j + ijkdelta(2)
                  iii = 2 * i + ijkdelta(1)
                  u_s(iii,  jjj,  kkk  ,b_recv,s) = u_s(i,j,k,b_send,s)
                  u_s(iii+1,jjj,  kkk  ,b_recv,s) = u_s(i,j,k,b_send,s)
                  u_s(iii,  jjj+1,kkk  ,b_recv,s) = u_s(i,j,k,b_send,s)
                  u_s(iii+1,jjj+1,kkk  ,b_recv,s) = u_s(i,j,k,b_send,s)
                  u_s(iii,  jjj,  kkk+1,b_recv,s) = u_s(i,j,k,b_send,s)
                  u_s(iii+1,jjj,  kkk+1,b_recv,s) = u_s(i,j,k,b_send,s)
                  u_s(iii,  jjj+1,kkk+1,b_recv,s) = u_s(i,j,k,b_send,s)
                  u_s(iii+1,jjj+1,kkk+1,b_recv,s) = u_s(i,j,k,b_send,s)
               enddo
            enddo
         enddo
      endif
   enddo
   endassociate
   endsubroutine update_ghost

   subroutine update_ghost_mpi(self, s, step)
   !< Update ghost cells within other processes.
   class(field_object), intent(inout)        :: self                       !< The field.
   integer(I4P),        intent(in)           :: s                          !< Stage.
   integer(I4P),        intent(in), optional :: step                       !< Step to be perfordmed.
   logical                                   :: steps(3)                   !< Steps to be performed.
   integer(I4P)                              :: i, j, k                    !< Counter.
   integer(I4P)                              :: iii, jjj, kkk              !< Counter.
   integer(I4P)                              :: fec, mf, rf, sf, n, p      !< Counter.
   integer(I4P), allocatable                 :: comm_map_send_ctr_ghost(:) !< Communication map, counters to send [procs_number+1].
   integer(I4P), allocatable                 :: comm_map_recv_ctr_ghost(:) !< Communication map, counters to recv [procs_number+1].
   integer(I4P)                              :: portion                    !< Portion of fec updated (0=>whole fec).
   integer(I4P)                              :: b_recv                     !< Index of receiving block.
   integer(I4P)                              :: b_send                     !< Index of sending block.
   integer(I4P)                              :: delta(3)                   !< Neighbor delta of current fec.
   integer(I4P)                              :: ijkmin(3)                  !< Lower limit of ijk indexes.
   integer(I4P)                              :: ijkmax(3)                  !< Upper limit of ijk indexes.
   integer(I4P)                              :: ijkdelta(3)                !< Delta offset of ghost-inner cells map same refinement.
   integer(I4P)                              :: ptr_start, ptr_end         !< Counter.
   integer(I4P)                              :: n_recv, n_send             !< Counter.
   integer(I4P)                              :: recv_rank                  !< Rank of receiving block.
   integer(I4P)                              :: send_rank                  !< Rank of sending block.

   associate(u_s=>self%u_s)

   steps = .true.
   if (present(step)) then
      steps = .false.
      steps(step) = .true.
   endif

   if (steps(1)) then
#ifdef _MPI_
      self%req_send_recv = MPI_REQUEST_NULL
#endif

      if ((.not.allocated(self%comm_map_recv_ghost)).and.&
          (.not.allocated(self%comm_map_send_ghost))) return

      comm_map_send_ctr_ghost = self%comm_map_send_ptr_ghost

      ! populate send buffer
      do sf=1, size(self%comm_map_send_ghost, dim=1)
         ! b_ghost   =     comm_map_send_ghost(sf, 1) ! block-index
         b_send    = self%comm_map_send_ghost(sf, 2) ! neighbor-block-index of block
         send_rank = self%comm_map_send_ghost(sf, 3)
         fec       = self%comm_map_send_ghost(sf, 4)
         portion   = self%comm_map_send_ghost(sf, 5)
         ijkmin    = self%comm_map_send_ghost(sf, 6:8)
         ijkmax    = self%comm_map_send_ghost(sf, 9:11)
         ijkdelta  = self%comm_map_send_ghost(sf, 12:14)
         if (portion==0_I4P) then
            ! sending to a block at my level
            do k=ijkmin(3), ijkmax(3)
               do j=ijkmin(2), ijkmax(2)
                  do i=ijkmin(1), ijkmax(1)
                     self%send_buffer_ghost(comm_map_send_ctr_ghost(send_rank)+1) = &
                        u_s(i+ijkdelta(1),j+ijkdelta(2),k+ijkdelta(3),b_send,s)
                     comm_map_send_ctr_ghost(send_rank) = comm_map_send_ctr_ghost(send_rank) + 1
                  enddo
               enddo
            enddo
         elseif (portion<0_I4P) then ! Beware! This is < 0 because the reference is the receiver
            ! sending to a block finer than me
            do k=ijkmin(3), ijkmax(3)
               do j=ijkmin(2), ijkmax(2)
                  do i=ijkmin(1), ijkmax(1)
                     do n=1,8
                        self%send_buffer_ghost(comm_map_send_ctr_ghost(send_rank)+1) = &
                           u_s(i,j,k,b_send,s)
                        comm_map_send_ctr_ghost(send_rank) = comm_map_send_ctr_ghost(send_rank) + 1
                     enddo
                  enddo
               enddo
            enddo
         else
            ! sending to a block coarser than me, loop is over the coarser grid
            do k=ijkmin(3), ijkmax(3)
               do j=ijkmin(2), ijkmax(2)
                  do i=ijkmin(1), ijkmax(1)
                     kkk = 2 * k + ijkdelta(3)
                     jjj = 2 * j + ijkdelta(2)
                     iii = 2 * i + ijkdelta(1)
                     self%send_buffer_ghost(comm_map_send_ctr_ghost(send_rank)+1) =        &
                        (u_s(iii,jjj,  kkk,  b_send,s) + u_s(iii+1,jjj,  kkk,  b_send,s) + &
                         u_s(iii,jjj+1,kkk,  b_send,s) + u_s(iii+1,jjj+1,kkk,  b_send,s) + &
                         u_s(iii,jjj,  kkk+1,b_send,s) + u_s(iii+1,jjj,  kkk+1,b_send,s) + &
                         u_s(iii,jjj+1,kkk+1,b_send,s) + u_s(iii+1,jjj+1,kkk+1,b_send,s)) / 8._R8P
                     comm_map_send_ctr_ghost(send_rank) = comm_map_send_ctr_ghost(send_rank) + 1
                  enddo
               enddo
            enddo
         endif
      enddo
   endif

   if (steps(2)) then
      ! receive
      do p=0, self%procs_number - 1_I4P
         ptr_start = self%comm_map_recv_ptr_ghost(p) + 1
         ptr_end   = self%comm_map_recv_ptr_ghost(p+1)
         n_recv    = ptr_end - ptr_start + 1
         if (n_recv > 0) then
#ifdef _MPI_
            call MPI_IRECV(self%recv_buffer_ghost(ptr_start), n_recv, MPI_REAL8, p, 100, MPI_COMM_WORLD, &
                           self%req_send_recv(p), self%error)
#endif
         endif
      enddo

      ! send
      do p=0, self%procs_number - 1_I4P
         ptr_start = self%comm_map_send_ptr_ghost(p) + 1
         ptr_end   = self%comm_map_send_ptr_ghost(p+1)
         n_send    = ptr_end - ptr_start + 1
         if (n_send > 0) then
#ifdef _MPI_
            call MPI_ISEND(self%send_buffer_ghost(ptr_start), n_send, MPI_REAL8, p, 100, MPI_COMM_WORLD, &
                           self%req_send_recv(p+self%procs_number), self%error)
#endif
         endif
      enddo
   endif

   if (steps(3)) then
      comm_map_recv_ctr_ghost = self%comm_map_recv_ptr_ghost
#ifdef _MPI_
      call MPI_WAITALL(self%procs_number * 2, self%req_send_recv, MPI_STATUSES_IGNORE, self%error)
#endif

      call MPI_BARRIER(MPI_COMM_WORLD, self%error)

      ! retrive from receive buffer
      do rf=1, size(self%comm_map_recv_ghost, dim=1)
         b_recv    = self%comm_map_recv_ghost(rf, 1) ! block-index
         ! b_recv    =     comm_map_recv_ghost(rf, 2) ! neighbor-block-index of block
         recv_rank = self%comm_map_recv_ghost(rf, 3)
         fec       = self%comm_map_recv_ghost(rf, 4)
         portion   = self%comm_map_recv_ghost(rf, 5)
         ijkmin    = self%comm_map_recv_ghost(rf, 6:8)
         ijkmax    = self%comm_map_recv_ghost(rf, 9:11)
         ijkdelta  = self%comm_map_recv_ghost(rf, 12:14)
         if (portion==0_I4P) then
            do k=ijkmin(3), ijkmax(3)
               do j=ijkmin(2), ijkmax(2)
                  do i=ijkmin(1), ijkmax(1)
                     u_s(i,j,k,b_recv,s) = self%recv_buffer_ghost(comm_map_recv_ctr_ghost(recv_rank)+1)
                     comm_map_recv_ctr_ghost(recv_rank) = comm_map_recv_ctr_ghost(recv_rank) + 1
                  enddo
               enddo
            enddo
         elseif (portion>0_I4P) then
            ! receiving from a block finer than me
            do k=ijkmin(3), ijkmax(3)
               do j=ijkmin(2), ijkmax(2)
                  do i=ijkmin(1), ijkmax(1)
                     u_s(i,j,k,b_recv,s) = self%recv_buffer_ghost(comm_map_recv_ctr_ghost(recv_rank)+1)
                     comm_map_recv_ctr_ghost(recv_rank) = comm_map_recv_ctr_ghost(recv_rank) + 1
                  enddo
               enddo
            enddo
         else
            ! receiving from a block coarser than me
            do k=ijkmin(3), ijkmax(3)
               do j=ijkmin(2), ijkmax(2)
                  do i=ijkmin(1), ijkmax(1)
                     kkk = 2 * k + ijkdelta(3)
                     jjj = 2 * j + ijkdelta(2)
                     iii = 2 * i + ijkdelta(1)
                     u_s(iii,  jjj,  kkk  ,b_recv,s) = self%recv_buffer_ghost(comm_map_recv_ctr_ghost(recv_rank)+1)
                     comm_map_recv_ctr_ghost(recv_rank) = comm_map_recv_ctr_ghost(recv_rank) + 1
                     u_s(iii+1,jjj,  kkk  ,b_recv,s) = self%recv_buffer_ghost(comm_map_recv_ctr_ghost(recv_rank)+1)
                     comm_map_recv_ctr_ghost(recv_rank) = comm_map_recv_ctr_ghost(recv_rank) + 1
                     u_s(iii,  jjj+1,kkk  ,b_recv,s) = self%recv_buffer_ghost(comm_map_recv_ctr_ghost(recv_rank)+1)
                     comm_map_recv_ctr_ghost(recv_rank) = comm_map_recv_ctr_ghost(recv_rank) + 1
                     u_s(iii+1,jjj+1,kkk  ,b_recv,s) = self%recv_buffer_ghost(comm_map_recv_ctr_ghost(recv_rank)+1)
                     comm_map_recv_ctr_ghost(recv_rank) = comm_map_recv_ctr_ghost(recv_rank) + 1
                     u_s(iii,  jjj,  kkk+1,b_recv,s) = self%recv_buffer_ghost(comm_map_recv_ctr_ghost(recv_rank)+1)
                     comm_map_recv_ctr_ghost(recv_rank) = comm_map_recv_ctr_ghost(recv_rank) + 1
                     u_s(iii+1,jjj,  kkk+1,b_recv,s) = self%recv_buffer_ghost(comm_map_recv_ctr_ghost(recv_rank)+1)
                     comm_map_recv_ctr_ghost(recv_rank) = comm_map_recv_ctr_ghost(recv_rank) + 1
                     u_s(iii,  jjj+1,kkk+1,b_recv,s) = self%recv_buffer_ghost(comm_map_recv_ctr_ghost(recv_rank)+1)
                     comm_map_recv_ctr_ghost(recv_rank) = comm_map_recv_ctr_ghost(recv_rank) + 1
                     u_s(iii+1,jjj+1,kkk+1,b_recv,s) = self%recv_buffer_ghost(comm_map_recv_ctr_ghost(recv_rank)+1)
                     comm_map_recv_ctr_ghost(recv_rank) = comm_map_recv_ctr_ghost(recv_rank) + 1
                  enddo
               enddo
            enddo
         endif
      enddo
   endif

   endassociate
   endsubroutine update_ghost_mpi

   ! private methods
   subroutine derefine(self, ratio, block_to_derefine, block_derefined)
   !< Derefine blocks.
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

   associate(ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk, u_work=>self%u_work)
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

                  u_work(i,     j,     k     ,ib) = (self%u(iii,jjj,  kkk,  ic1) + self%u(iii+1,jjj,  kkk,  ic1) + &
                                                     self%u(iii,jjj+1,kkk,  ic1) + self%u(iii+1,jjj+1,kkk,  ic1) + &
                                                     self%u(iii,jjj,  kkk+1,ic1) + self%u(iii+1,jjj,  kkk+1,ic1) + &
                                                     self%u(iii,jjj+1,kkk+1,ic1) + self%u(iii+1,jjj+1,kkk+1,ic1)) / 8._R8P

                  u_work(i+ni/2,j,     k     ,ib) = (self%u(iii,jjj,  kkk,  ic2) + self%u(iii+1,jjj,  kkk,  ic2) + &
                                                     self%u(iii,jjj+1,kkk,  ic2) + self%u(iii+1,jjj+1,kkk,  ic2) + &
                                                     self%u(iii,jjj,  kkk+1,ic2) + self%u(iii+1,jjj,  kkk+1,ic2) + &
                                                     self%u(iii,jjj+1,kkk+1,ic2) + self%u(iii+1,jjj+1,kkk+1,ic2)) / 8._R8P

                  u_work(i,     j+nj/2,k     ,ib) = (self%u(iii,jjj,  kkk,  ic3) + self%u(iii+1,jjj,  kkk,  ic3) + &
                                                     self%u(iii,jjj+1,kkk,  ic3) + self%u(iii+1,jjj+1,kkk,  ic3) + &
                                                     self%u(iii,jjj,  kkk+1,ic3) + self%u(iii+1,jjj,  kkk+1,ic3) + &
                                                     self%u(iii,jjj+1,kkk+1,ic3) + self%u(iii+1,jjj+1,kkk+1,ic3)) / 8._R8P

                  u_work(i+ni/2,j+nj/2,k     ,ib) = (self%u(iii,jjj,  kkk,  ic4) + self%u(iii+1,jjj,  kkk,  ic4) + &
                                                     self%u(iii,jjj+1,kkk,  ic4) + self%u(iii+1,jjj+1,kkk,  ic4) + &
                                                     self%u(iii,jjj,  kkk+1,ic4) + self%u(iii+1,jjj,  kkk+1,ic4) + &
                                                     self%u(iii,jjj+1,kkk+1,ic4) + self%u(iii+1,jjj+1,kkk+1,ic4)) / 8._R8P

                  u_work(i,     j,     k+nk/2,ib) = (self%u(iii,jjj,  kkk,  ic5) + self%u(iii+1,jjj,  kkk,  ic5) + &
                                                     self%u(iii,jjj+1,kkk,  ic5) + self%u(iii+1,jjj+1,kkk,  ic5) + &
                                                     self%u(iii,jjj,  kkk+1,ic5) + self%u(iii+1,jjj,  kkk+1,ic5) + &
                                                     self%u(iii,jjj+1,kkk+1,ic5) + self%u(iii+1,jjj+1,kkk+1,ic5)) / 8._R8P

                  u_work(i+ni/2,j,     k+nk/2,ib) = (self%u(iii,jjj,  kkk,  ic6) + self%u(iii+1,jjj,  kkk,  ic6) + &
                                                     self%u(iii,jjj+1,kkk,  ic6) + self%u(iii+1,jjj+1,kkk,  ic6) + &
                                                     self%u(iii,jjj,  kkk+1,ic6) + self%u(iii+1,jjj,  kkk+1,ic6) + &
                                                     self%u(iii,jjj+1,kkk+1,ic6) + self%u(iii+1,jjj+1,kkk+1,ic6)) / 8._R8P

                  u_work(i,     j+nj/2,k+nk/2,ib) = (self%u(iii,jjj,  kkk,  ic7) + self%u(iii+1,jjj,  kkk,  ic7) + &
                                                     self%u(iii,jjj+1,kkk,  ic7) + self%u(iii+1,jjj+1,kkk,  ic7) + &
                                                     self%u(iii,jjj,  kkk+1,ic7) + self%u(iii+1,jjj,  kkk+1,ic7) + &
                                                     self%u(iii,jjj+1,kkk+1,ic7) + self%u(iii+1,jjj+1,kkk+1,ic7)) / 8._R8P

                  u_work(i+ni/2,j+nj/2,k+nk/2,ib) = (self%u(iii,jjj,  kkk,  ic8) + self%u(iii+1,jjj,  kkk,  ic8) + &
                                                     self%u(iii,jjj+1,kkk,  ic8) + self%u(iii+1,jjj+1,kkk,  ic8) + &
                                                     self%u(iii,jjj,  kkk+1,ic8) + self%u(iii+1,jjj,  kkk+1,ic8) + &
                                                     self%u(iii,jjj+1,kkk+1,ic8) + self%u(iii+1,jjj+1,kkk+1,ic8)) / 8._R8P
               enddo
            enddo
         enddo

         self%u(1:ni,1:nj,1:nk,ib) = u_work(1:ni,1:nj,1:nk,ib)

         self%code(ib) = block_derefined(1,b)

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
   endassociate
   endsubroutine derefine

   subroutine refine(self, ratio, block_to_refine, block_refined)
   !< Refine blocks.
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

   associate(ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk,       &
             gc1=>self%grid%gc1, gc2=>self%grid%gc2, gc3=>self%grid%gc3, &
             gc4=>self%grid%gc4, gc5=>self%grid%gc5, gc6=>self%grid%gc6)
   if (allocated(block_to_refine)) then
      do b=1, size(block_to_refine, dim=2)
         if (self%myrank /= block_to_refine(2,b)) cycle
         ib = block_to_refine(1,b)

         self%u_work(:,:,:,ib) = self%u(:,:,:,ib)

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
                     self%u(i_fine:i_fine+1,j_fine:j_fine+1,k_fine:k_fine+1,ic) = 0._R8P
                     do k_delta=0,1
                     do j_delta=0,1
                     do i_delta=0,1
                     self%u(i_fine,  j_fine,  k_fine,  ic) = self%u(i_fine,j_fine,k_fine,ic) + &
                                                             (0.25_R8P + i_delta * 0.5_R8P) *  &
                                                             (0.25_R8P + j_delta * 0.5_R8P) *  &
                                                             (0.25_R8P + k_delta * 0.5_R8P) *  &
                                                             self%u_work(i+i_delta-1, j+j_delta-1, k+k_delta-1,ib)
                     self%u(i_fine+1,j_fine,  k_fine,  ic) = self%u(i_fine+1,j_fine,k_fine,ic) + &
                                                             (0.75_R8P - i_delta * 0.5_R8P) *    &
                                                             (0.25_R8P + j_delta * 0.5_R8P) *    &
                                                             (0.25_R8P + k_delta * 0.5_R8P) *    &
                                                             self%u_work(i+i_delta,   j+j_delta-1, k+k_delta-1,ib)
                     self%u(i_fine,  j_fine+1,k_fine,  ic) = self%u(i_fine,j_fine+1,k_fine,ic) + &
                                                             (0.25_R8P + i_delta * 0.5_R8P) *    &
                                                             (0.75_R8P - j_delta * 0.5_R8P) *    &
                                                             (0.25_R8P + k_delta * 0.5_R8P) *    &
                                                             self%u_work(i+i_delta-1, j+j_delta  , k+k_delta-1,ib)
                     self%u(i_fine+1,j_fine+1,k_fine,  ic) = self%u(i_fine+1,j_fine+1,k_fine,ic) + &
                                                             (0.75_R8P - i_delta * 0.5_R8P) *      &
                                                             (0.75_R8P - j_delta * 0.5_R8P) *      &
                                                             (0.25_R8P + k_delta * 0.5_R8P) *      &
                                                             self%u_work(i+i_delta,   j+j_delta  , k+k_delta-1,ib)
                     self%u(i_fine,  j_fine,  k_fine+1,ic) = self%u(i_fine,j_fine,k_fine+1,ic) + &
                                                             (0.25_R8P + i_delta * 0.5_R8P) *    &
                                                             (0.25_R8P + j_delta * 0.5_R8P) *    &
                                                             (0.75_R8P - k_delta * 0.5_R8P) *    &
                                                             self%u_work(i+i_delta-1, j+j_delta-1, k+k_delta  ,ib)
                     self%u(i_fine+1,j_fine,  k_fine+1,ic) = self%u(i_fine+1,j_fine,k_fine+1,ic) + &
                                                             (0.75_R8P - i_delta * 0.5_R8P) *      &
                                                             (0.25_R8P + j_delta * 0.5_R8P) *      &
                                                             (0.75_R8P - k_delta * 0.5_R8P) *      &
                                                             self%u_work(i+i_delta,   j+j_delta-1, k+k_delta  ,ib)
                     self%u(i_fine,  j_fine+1,k_fine+1,ic) = self%u(i_fine,j_fine+1,k_fine+1,ic) + &
                                                             (0.25_R8P + i_delta * 0.5_R8P) *      &
                                                             (0.75_R8P - j_delta * 0.5_R8P) *      &
                                                             (0.75_R8P - k_delta * 0.5_R8P) *      &
                                                             self%u_work(i+i_delta-1, j+j_delta  , k+k_delta  ,ib)
                     self%u(i_fine+1,j_fine+1,k_fine+1,ic) = self%u(i_fine+1,j_fine+1,k_fine+1,ic) + &
                                                             (0.75_R8P - i_delta * 0.5_R8P) *        &
                                                             (0.75_R8P - j_delta * 0.5_R8P) *        &
                                                             (0.75_R8P - k_delta * 0.5_R8P) *        &
                                                             self%u_work(i+i_delta,   j+j_delta  , k+k_delta  ,ib)

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

      ! do b=1, size(block_to_refine, dim=2)
      !    if (self%myrank /= block_to_refine(2,b)) cycle
      !    ib = block_to_refine(1,b)

      !    dx = self%emax(1,ib) - self%emin(1,ib)
      !    dy = self%emax(2,ib) - self%emin(2,ib)
      !    dz = self%emax(3,ib) - self%emin(3,ib)

      !    ii = 1
      !    do k=0, 1
      !       do j=0, 1
      !          do i=0, 1
      !             ic = block_refined(2,(b-1)*ratio + ii)

      !             self%emin(1,ic) = self%emin(1,ib) + i * dx/2
      !             self%emin(2,ic) = self%emin(2,ib) + j * dy/2
      !             self%emin(3,ic) = self%emin(3,ib) + k * dz/2

      !             self%emax(1,ic) = self%emin(1,ic) + dx/2
      !             self%emax(2,ic) = self%emin(2,ic) + dy/2
      !             self%emax(3,ic) = self%emin(3,ic) + dz/2

      !             ii = ii + 1
      !          enddo
      !       enddo
      !    enddo
      ! enddo
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
   if (allocated(rhs%x_cell)) then
      lhs%x_cell = rhs%x_cell
   else
      if (allocated(lhs%x_cell)) deallocate(lhs%x_cell)
   endif
   if (allocated(rhs%y_cell)) then
      lhs%y_cell = rhs%y_cell
   else
      if (allocated(lhs%y_cell)) deallocate(lhs%y_cell)
   endif
   if (allocated(rhs%z_cell)) then
      lhs%z_cell = rhs%z_cell
   else
      if (allocated(lhs%z_cell)) deallocate(lhs%z_cell)
   endif
   if (allocated(rhs%x_node)) then
      lhs%x_node = rhs%x_node
   else
      if (allocated(lhs%x_node)) deallocate(lhs%x_node)
   endif
   if (allocated(rhs%y_node)) then
      lhs%y_node = rhs%y_node
   else
      if (allocated(lhs%y_node)) deallocate(lhs%y_node)
   endif
   if (allocated(rhs%z_node)) then
      lhs%z_node = rhs%z_node
   else
      if (allocated(lhs%z_node)) deallocate(lhs%z_node)
   endif
   if (allocated(rhs%u)) then
      lhs%u = rhs%u
   else
      if (allocated(lhs%u)) deallocate(lhs%u)
   endif
   if (allocated(rhs%u_work)) then
      lhs%u_work = rhs%u_work
   else
      if (allocated(lhs%u_work)) deallocate(lhs%u_work)
   endif
   if (allocated(rhs%local_map_ghost)) then
      lhs%local_map_ghost = rhs%local_map_ghost
   else
      if (allocated(lhs%local_map_ghost)) deallocate(lhs%local_map_ghost)
   endif
   ! MPI data
   lhs%myrank        = rhs%myrank
   lhs%procs_number  = rhs%procs_number
   if (allocated(rhs%blocks_numbers)) then
      lhs%blocks_numbers = rhs%blocks_numbers
   else
      if (allocated(lhs%blocks_numbers)) deallocate(lhs%blocks_numbers)
   endif
   ! RK data
   lhs%alph = rhs%alph
   lhs%beta = rhs%beta
   lhs%gamm = rhs%gamm
   endsubroutine field_assign_field
endmodule adam_field_object
