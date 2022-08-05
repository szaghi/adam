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
!<  /      .    .    .    .    .               _____|
!< O-------+----+----+----+----+-------------> I
!<         1    2    3    4    5
!<```

use adam_grid_object, only : grid_object
use adam_mpih_object, only : mpih_object
use adam_memory_cpu_lib
use adam_parameters
use FINER, only : file_ini
use PENF
use MPI
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

implicit none
private
public :: field_object

type :: field_object
   !< Field class definition.
   ! mesh related data, unrelated to field equations
   type(mpih_object)          :: mpih                !< MPI handler.
   type(grid_object), pointer :: grid=>null()        !< Grid data.
   integer(I4P)               :: nv=1_I4P            !< Number of field variables.
   integer(I4P)               :: block_weight=0_I4P  !< Block weight, `cells_number * variables_number`.
   integer(I4P)               :: nb=0_I4P            !< Number of all blocks that can be stored.
   integer(I4P)               :: blocks_number=0_I4P !< Number of blocks actually stored.
   integer(I8P), allocatable  :: code(:)             !< Morton codes [nb].
   integer(I4P), allocatable  :: coordinates(:,:)    !< Coordinates IJKL for each block [nb,4].
   real(R8P),    allocatable  :: emin(:,:)           !< Coordinates of minimum abscissa of each block [3,nb].
   real(R8P),    allocatable  :: emax(:,:)           !< Coordinates of maximum abscissa of each block [3,nb].
   real(R8P),    allocatable  :: dxyz(:,:)           !< Space steps of each block [3,nb].
   real(R8P),    allocatable  :: x_node(:,:)         !< X node coordinates.
   real(R8P),    allocatable  :: y_node(:,:)         !< Y node coordinates.
   real(R8P),    allocatable  :: z_node(:,:)         !< Z node coordinates.
   real(R8P),    allocatable  :: x_cell(:,:)         !< X cell coordinates.
   real(R8P),    allocatable  :: y_cell(:,:)         !< Y cell coordinates.
   real(R8P),    allocatable  :: z_cell(:,:)         !< Z cell coordinates.
   ! MPI data, unrelated to field equations
   integer(I4P), allocatable :: blocks_numbers(:)         !< Number of blocks actually stored in all processes.
   integer(I4P), allocatable :: refinements_needed(:)     !< Refinements needed of my blocks.
   integer(I4P), allocatable :: refinements_needed_all(:) !< Refinements needed of all blocks.
   integer(I4P), allocatable :: disp_count(:)             !< Displacement of blocks that are received from process.
   integer(I4P)              :: inner_blocks_number=0_I4P !< Number of inner blocks where I need fecs.
   integer(I4P), allocatable :: req_send_recv(:)          !< MPI request receive flags.
   ! field equations data
   real(R8P), allocatable :: q(     :,:,:,:,:) !< Field cell centered variables.
   real(R8P), allocatable :: q_work(:,:,:,:,:) !< Field cell centered variables, working buffer memory.
   contains
      ! public methods
      procedure, pass(self) :: adapt                         !< Adapt field accordingly to refine/derefine necessity.
      procedure, pass(self) :: blocks_reorder                !< Reorder blocks indexes in field.
      procedure, pass(self) :: compute_metrics               !< Compute metrics of each block.
      procedure, pass(self) :: description                   !< Return pretty-printed object description.
      procedure, pass(self) :: do_caxis_intersect            !< Return true if a block is intersected by coordinate-axis.
      procedure, pass(self) :: do_cplane_intersect           !< Return true if a block is intersected by coordinate-plane.
      procedure, pass(self) :: do_ray_intersect              !< Return true if a block is intersected by ray.
      procedure, pass(self) :: initialize                    !< Initialize the field.
      procedure, pass(self) :: load_blocks                   !< Load blocks data, used for restarting.
      procedure, pass(self) :: load_from_ini_file            !< Load object data from INI file.
      procedure, pass(self) :: mark_all_blocks               !< Mark all blocks to be refined, derefined, ecc.
      procedure, pass(self) :: mark_sphere                   !< Mark blocks to be refined/derefined by sphere distance.
      procedure, pass(self) :: mpi_gather_refinements_needed !< Gather blocks refinement needed status between MPI processes.
      procedure, pass(self) :: mpi_redistribute              !< Redistribute blocks to processes.
      procedure, pass(self) :: save_blocks                   !< Save blocks data, used for restarting.
      ! private methods
      procedure, pass(self), private :: derefine !< Derefine blocks.
      procedure, pass(self), private :: refine   !< Refine blocks.
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
   if (allocated(coordinates_new)) deallocate(coordinates_new)
   if (allocated(code_new)) deallocate(code_new)
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

   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(field_object), intent(in) :: self             !< The field.
   character(len=:), allocatable   :: desc             !< Description.
   character(len=1), parameter     :: NL=new_line('a') !< New line character.

   desc =       self%mpih%myrankstr//'field main data'                                               //NL
   desc = desc//self%mpih%myrankstr//'  field variables number (nv): '//trim(str(self%nv           ))//NL
   desc = desc//self%mpih%myrankstr//'  all blocks number (nb):      '//trim(str(self%nb           ))//NL
   desc = desc//self%mpih%myrankstr//'  blocks number:               '//trim(str(self%blocks_number))//NL
   desc = desc//self%mpih%myrankstr//'  block weight:                '//trim(str(self%block_weight ))//NL
   desc = desc//self%mpih%myrankstr//'  q shape:                     '//trim(str(shape(self%q)     ))
   endfunction description

   function do_caxis_intersect(self, b, caxis_origin, caxis_direction, caxis_block_indexes) result(do_intersect)
   !< Return true if a block is intersected by coordinate-axis.
   class(field_object), intent(inout)         :: self                   !< The field.
   integer(I4P),        intent(in)            :: b                      !< Block index.
   real(R8P),           intent(in)            :: caxis_origin(3)        !< Coordinate-axis origin.
   real(R8P),           intent(in)            :: caxis_direction(3)     !< Coordinate-axis direction.
   integer(I4P),        intent(out), optional :: caxis_block_indexes(3) !< Block-local indexes of caxis intersection.
   logical                                    :: do_intersect           !< Test result.

   do_intersect = .false.
   if     (nint(caxis_direction(1))==1) then
      if ((caxis_origin(2) >= self%emin(2,b)).and.(caxis_origin(2) <= self%emax(2,b)).and. &
          (caxis_origin(3) >= self%emin(3,b)).and.(caxis_origin(3) <= self%emax(3,b))) then
         do_intersect = .true.
         if (present(caxis_block_indexes)) then
            caxis_block_indexes(2) = ceiling((caxis_origin(2) - self%emin(2,b)) / self%dxyz(2,b), I4P)
            caxis_block_indexes(3) = ceiling((caxis_origin(3) - self%emin(3,b)) / self%dxyz(3,b), I4P)
         endif
      endif
   elseif (nint(caxis_direction(2))==1) then
      if ((caxis_origin(1) >= self%emin(1,b)).and.(caxis_origin(1) <= self%emax(1,b)).and. &
          (caxis_origin(3) >= self%emin(3,b)).and.(caxis_origin(3) <= self%emax(3,b))) then
         do_intersect = .true.
         if (present(caxis_block_indexes)) then
            caxis_block_indexes(1) = ceiling((caxis_origin(1) - self%emin(1,b)) / self%dxyz(1,b), I4P)
            caxis_block_indexes(3) = ceiling((caxis_origin(3) - self%emin(3,b)) / self%dxyz(3,b), I4P)
         endif
      endif
   elseif (nint(caxis_direction(3))==1) then
      if ((caxis_origin(2) >= self%emin(2,b)).and.(caxis_origin(2) <= self%emax(2,b)).and. &
          (caxis_origin(1) >= self%emin(1,b)).and.(caxis_origin(1) <= self%emax(1,b))) then
         do_intersect = .true.
         if (present(caxis_block_indexes)) then
            caxis_block_indexes(2) = ceiling((caxis_origin(2) - self%emin(2,b)) / self%dxyz(2,b), I4P)
            caxis_block_indexes(1) = ceiling((caxis_origin(1) - self%emin(1,b)) / self%dxyz(1,b), I4P)
         endif
      endif
   endif
   endfunction do_caxis_intersect

   function do_cplane_intersect(self, b, cplane_origin, cplane_normal, cplane_block_indexes) result(do_intersect)
   !< Return true if a block is intersected by coordinate-plane.
   class(field_object), intent(inout)         :: self                    !< The field.
   integer(I4P),        intent(in)            :: b                       !< Block index.
   real(R8P),           intent(in)            :: cplane_origin(3)        !< Coordinate-plane origin.
   real(R8P),           intent(in)            :: cplane_normal(3)        !< Coordinate-plane normal.
   integer(I4P),        intent(out), optional :: cplane_block_indexes(3) !< Block-local indexes of cplane intersection.
   logical                                    :: do_intersect            !< Test result.

   do_intersect = self%grid%do_cplane_intersect(emin=self%emin(:,b),         &
                                                emax=self%emax(:,b),         &
                                                dxyz=self%dxyz(:,b),         &
                                                cplane_origin=cplane_origin, &
                                                cplane_normal=cplane_normal, &
                                                cplane_block_indexes=cplane_block_indexes)
   endfunction do_cplane_intersect

   function do_ray_intersect(self, b, ray_origin, ray_direction) result(do_intersect)
   !< Return true if a block is intersected by ray from ray_origin and oriented as ray_direction vector.
   class(field_object), intent(inout) :: self             !< The field.
   integer(I4P),        intent(in)    :: b                !< Block index.
   real(R8P),           intent(in)    :: ray_origin(3)    !< Ray origin.
   real(R8P),           intent(in)    :: ray_direction(3) !< Ray direction.
   logical                            :: do_intersect     !< Test result.
   logical                            :: must_return      !< Flag to check when to return from procedure.
   real(R8P)                          :: tmin, tmax       !< Minimum maximum ray intersections with box slabs.

   do_intersect = .false.
   must_return = .false.
   tmin = 0._R8P
   tmax = MaxR8P
   call check_slab(bmin=self%emin(1,b), bmax=self%emax(1,b), o=ray_origin(1), d=ray_direction(1), &
                   must_return=must_return, tmin=tmin, tmax=tmax)
   if (must_return) return
   call check_slab(bmin=self%emin(2,b), bmax=self%emax(2,b), o=ray_origin(2), d=ray_direction(2), &
                   must_return=must_return, tmin=tmin, tmax=tmax)
   if (must_return) return
   call check_slab(bmin=self%emin(3,b), bmax=self%emax(3,b), o=ray_origin(3), d=ray_direction(3), &
                   must_return=must_return, tmin=tmin, tmax=tmax)
   if (must_return) return
   ! ray intersects all 3 slabs
   do_intersect = .true.
   contains
      subroutine check_slab(bmin, bmax, o, d, must_return, tmin, tmax)
      !< Perform ray intersection check in a direction-split fashion over slabs.
      real(R8P), intent(in)    :: bmin        !< Box minimum bound in the current direction.
      real(R8P), intent(in)    :: bmax        !< Box maximum bound in the current direction.
      real(R8P), intent(in)    :: o           !< Ray origin in the current direction.
      real(R8P), intent(in)    :: d           !< Ray slope in the current direction.
      logical,   intent(inout) :: must_return !< Flag to check when to return from procedure.
      real(R8P), intent(inout) :: tmin, tmax  !< Minimum maximum ray intersections with box slabs.
      real(R8P)                :: ood, t1, t2 !< Intersection coefficients.
      real(R8P)                :: tmp         !< Temporary buffer.

      if ((d) < 1.E-16_R8P) then
         ! ray is parallel to slab, no hit if origin not within slab
         if ((o < bmin).or.(o > bmax)) then
            must_return = .true.
            return
         endif
      else
         ! compute intersection t value of ray with near and far plane of slab
         ood = 1._R8P / d
         t1 = (bmin - o) * ood
         t2 = (bmax - o) * ood
         ! make t1 be intersection with near plane, t2 with far plane
         if (t1 > t2) then
            tmp = t1
            t1 = t2
            t2 = tmp
         endif
         ! compute the intersection of slab intersection intervals
         if (t1 > tmin) tmin = t1
         if (t2 > tmax) tmax = t2
         ! exit with no collision as soon as slab intersection becomes empty
         if (tmin > tmax) then
            must_return = .true.
            return
         endif
      endif
      endsubroutine check_slab
   endfunction do_ray_intersect

   subroutine initialize(self, grid, file_parameters, nv, nb)
   !< Initialize field.
   class(field_object), intent(inout)           :: self            !< The field.
   type(grid_object),   intent(in), target      :: grid            !< Grid data.
   type(file_ini),      intent(inout), optional :: file_parameters !< INI file handler.
   integer(I4P),        intent(in),    optional :: nv              !< Number of field variables.
   integer(I4P),        intent(in),    optional :: nb              !< Number of all blocks that can be stored.

   call self%mpih%initialize
   print '(A)', self%mpih%myrankstr//'field%initialize start'
   self%grid => grid
   if (present(file_parameters)) call self%load_from_ini_file(file_parameters)
   ! parameters explicitely passed ovveride ones file-passed
   if (present(nv)) self%nv  = nv
   self%block_weight = (self%grid%ngc+self%grid%ni+self%grid%ngc)* &
                       (self%grid%ngc+self%grid%nj+self%grid%ngc)* &
                       (self%grid%ngc+self%grid%nk+self%grid%ngc)*self%nv
   if (present(nb)) self%nb  = nb
   if (self%nb>0) then
      call alloc_var_cpu(var=self%code,  &
                         ulb=[1,self%nb],&
                         msg=self%mpih%myrankstr//'field%initialize(code) ', verbose=.true.)
      self%code    = -2_I8P
      self%code(1) = -1_I8P ! first block is assumed to be ADAM
      call alloc_var_cpu(var=self%coordinates,                &
                         ulb=reshape([1,4, 1,self%nb],[2,2]), &
                         msg=self%mpih%myrankstr//'field%initialize(coordinates) ', verbose=.true.)
      call alloc_var_cpu(var=self%emin,                       &
                         ulb=reshape([1,3, 1,self%nb],[2,2]), &
                         msg=self%mpih%myrankstr//'field%initialize(emin) ', verbose=.true.)
      call alloc_var_cpu(var=self%emax,                       &
                         ulb=reshape([1,3, 1,self%nb],[2,2]), &
                         msg=self%mpih%myrankstr//'field%initialize(emax) ', verbose=.true.)
      self%emin(:,1) = self%grid%domain_emin
      self%emax(:,1) = self%grid%domain_emax
      call alloc_var_cpu(var=self%dxyz,                       &
                         ulb=reshape([1,3, 1,self%nb],[2,2]), &
                         msg=self%mpih%myrankstr//'field%initialize(dxyz) ', verbose=.true.)
      call alloc_var_cpu(var=self%x_cell,                                         &
                         ulb=reshape([1-self%grid%ngc,self%grid%ni+self%grid%ngc, &
                                      1,self%nb],[2,2]),                          &
                         msg=self%mpih%myrankstr//'field%initialize(x_cell) ', verbose=.true.)
      call alloc_var_cpu(var=self%y_cell,                                         &
                         ulb=reshape([1-self%grid%ngc,self%grid%nj+self%grid%ngc, &
                                      1,self%nb],[2,2]),                          &
                         msg=self%mpih%myrankstr//'field%initialize(y_cell) ', verbose=.true.)
      call alloc_var_cpu(var=self%z_cell,                                         &
                         ulb=reshape([1-self%grid%ngc,self%grid%nk+self%grid%ngc, &
                                      1,self%nb],[2,2]),                          &
                         msg=self%mpih%myrankstr//'field%initialize(z_cell) ', verbose=.true.)
      call alloc_var_cpu(var=self%x_node,                                         &
                         ulb=reshape([0-self%grid%ngc,self%grid%ni+self%grid%ngc, &
                                      1,self%nb],[2,2]),                          &
                         msg=self%mpih%myrankstr//'field%initialize(x_node) ', verbose=.true.)
      call alloc_var_cpu(var=self%y_node,                                         &
                         ulb=reshape([0-self%grid%ngc,self%grid%nj+self%grid%ngc, &
                                      1,self%nb],[2,2]),                          &
                         msg=self%mpih%myrankstr//'field%initialize(y_node) ', verbose=.true.)
      call alloc_var_cpu(var=self%z_node,                                         &
                         ulb=reshape([0-self%grid%ngc,self%grid%nk+self%grid%ngc, &
                                      1,self%nb],[2,2]),                          &
                         msg=self%mpih%myrankstr//'field%initialize(z_node) ', verbose=.true.)
      call alloc_var_cpu(var=self%q,                                              &
                         ulb=reshape([1,self%nv,                                  &
                                      1-self%grid%ngc,self%grid%ni+self%grid%ngc, &
                                      1-self%grid%ngc,self%grid%nj+self%grid%ngc, &
                                      1-self%grid%ngc,self%grid%nk+self%grid%ngc, &
                                      1,self%nb],[2,5]),                          &
                         msg=self%mpih%myrankstr//'field%initialize(q) ', verbose=.true.)
      call alloc_var_cpu(var=self%q_work,                                         &
                         ulb=reshape([1,self%nv,                                  &
                                      1-self%grid%ngc,self%grid%ni+self%grid%ngc, &
                                      1-self%grid%ngc,self%grid%nj+self%grid%ngc, &
                                      1-self%grid%ngc,self%grid%nk+self%grid%ngc, &
                                      1,self%nb],[2,5]),                          &
                         msg=self%mpih%myrankstr//'field%initialize(q_work) ', verbose=.true.)
      self%q = 0._R8P
      self%q_work = 0._R8P
   endif
   call alloc_var_cpu(var=self%blocks_numbers, &
                      ulb=[0,self%mpih%procs_number-1],  &
                      msg=self%mpih%myrankstr//'field%initialize(blocks_numbers) ', verbose=.true.)
   call alloc_var_cpu(var=self%req_send_recv,  &
                      ulb=[0,self%mpih%procs_number*2-1],&
                      msg=self%mpih%myrankstr//'field%initialize(req_send_recv) ', verbose=.true.)
   print '(A)', self%description()
   print '(A)', self%mpih%myrankstr//'field%initialize finish'
   endsubroutine initialize

   subroutine load_blocks(self, basename)
   !< Load blocks data, used for restarting.
   !<
   !< Note: blocks memory must be already initialized with enough memory (proper nv,ni,nj,nk,ngc and nb>=blocks number).
   class(field_object), intent(inout) :: self                !< The field.
   character(*),        intent(in)    :: basename            !< Output base name.
   integer(I4P)                       :: file_unit           !< Output file unit.
   logical                            :: file_exist          !< Flag to check file's existance.
   integer(I4P)                       :: blocks_number       !< Blocks number.
   integer(I4P)                       :: nv, ni, nj, nk, ngc !< Dimensions.
   integer(I4P)                       :: b                   !< Counter.

   inquire(file=trim(adjustl(basename))//'-proc'//trim(strz(self%mpih%myrank,6))//'.fbd', exist=file_exist)
   if (file_exist) then
      print '(A)', self%mpih%myrankstr//'load field blocks from file '//trim(adjustl(basename))//&
                   '-proc'//trim(strz(self%mpih%myrank,6))//'.fbd'
      open(newunit=file_unit,                                                             &
           file=trim(adjustl(basename))//'-proc'//trim(strz(self%mpih%myrank,6))//'.fbd', &
           form='UNFORMATTED',                                                            &
           access='STREAM')
      read(unit=file_unit) nv, ni, nj, nk, ngc
      if (nv==self%nv.and.ni==self%grid%ni.and.nj==self%grid%nj.and.nk==self%grid%nk.and.ngc==self%grid%ngc) then
         read(unit=file_unit) blocks_number
         print '(A)', self%mpih%myrankstr//'field blocks number '//trim(str(blocks_number))
         if (blocks_number <= self%nb) then
            self%blocks_number = blocks_number
            do b=1, self%blocks_number
               read(unit=file_unit) self%code(b)
               read(unit=file_unit) self%coordinates(1:4,b)
               read(unit=file_unit) self%q(1:self%nv,                                  &
                                           1-self%grid%ngc:self%grid%ni+self%grid%ngc, &
                                           1-self%grid%ngc:self%grid%nj+self%grid%ngc, &
                                           1-self%grid%ngc:self%grid%nk+self%grid%ngc,b)
            enddo
            call self%compute_metrics
         else
            write(stderr, '(A)') self%mpih%myrankstr//'ERROR: blocks number to read "'//trim(str(blocks_number))//&
                                 '" is greater than blocks allocated "'//trim(str(self%nb))//'"!'
         endif
      else
         write(stderr, '(A)') self%mpih%myrankstr//'ERROR: blocks dimensions to read "'//trim(str([nv, ni, nj, nk, ngc]))//&
                              '" are different than blocks allocated "'//trim(str([self%nv,                                &
                                                                                   self%grid%ni,                           &
                                                                                   self%grid%nj,                           &
                                                                                   self%grid%nk,                           &
                                                                                   self%grid%ngc]))//'"!'
      endif
      close(file_unit)
      print '(A)', self%mpih%myrankstr//'load field blocks from file '//&
                   trim(adjustl(basename))//'-proc'//trim(strz(self%mpih%myrank,6))//'.fbd completed'
   else
      write(stderr, '(A)') self%mpih%myrankstr//'WARNING: file "'//&
                           trim(adjustl(basename))//'-proc'//trim(strz(self%mpih%myrank,6))//'.fbd" does not exist!'
   endif
   endsubroutine load_blocks

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
   allocate(recv_count(0:self%mpih%procs_number - 1))
   call MPI_ALLGATHER(self%blocks_number, 1_I4P, MPI_INTEGER, recv_count, 1_I4P, MPI_INTEGER, MPI_COMM_WORLD, self%mpih%error)

   ! computing displacement counts
   if (allocated(self%disp_count)) deallocate(self%disp_count)
   allocate(self%disp_count(0:self%mpih%procs_number - 1))
   self%disp_count = 0_I4P
   do p=1, self%mpih%procs_number - 1
      self%disp_count(p) = self%disp_count(p-1) + recv_count(p-1)
   enddo

   if (allocated(self%refinements_needed_all)) deallocate(self%refinements_needed_all)
   allocate(self%refinements_needed_all(sum(recv_count, dim=1)))
   call MPI_ALLGATHERV(self%refinements_needed, self%blocks_number, MPI_INTEGER, &
                       self%refinements_needed_all, recv_count, self%disp_count, MPI_INTEGER, MPI_COMM_WORLD, self%mpih%error)
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

   allocate(req_recv(0:self%mpih%procs_number-1))
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

   do p=0, self%mpih%procs_number - 1_I4P
      ptr_start = comm_map_recv_ptr(p)   * self%block_weight + 1
      ptr_end   = comm_map_recv_ptr(p+1) * self%block_weight
      n_recv    = ptr_end - ptr_start + 1
      if (n_recv > 0) then
         call MPI_IRECV(recv_buffer(ptr_start), n_recv, MPI_REAL8, p, 100, MPI_COMM_WORLD, req_recv(p), self%mpih%error)
      endif
   enddo

   do p=0, self%mpih%procs_number - 1_I4P
      ptr_start = comm_map_send_ptr(p)   * self%block_weight + 1
      ptr_end   = comm_map_send_ptr(p+1) * self%block_weight
      n_send    = ptr_end - ptr_start + 1
      if (n_send > 0) then
         call MPI_SEND(send_buffer(ptr_start), n_send, MPI_REAL8, p, 100, MPI_COMM_WORLD, self%mpih%error)
      endif
   enddo

   call MPI_WAITALL(self%mpih%procs_number, req_recv, MPI_STATUSES_IGNORE, self%mpih%error)

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

   subroutine save_blocks(self, basename)
   !< Save blocks data, used for restarting.
   class(field_object), intent(in) :: self      !< The field.
   character(*),        intent(in) :: basename  !< Output base name.
   integer(I4P)                    :: file_unit !< Output file unit.
   integer(I4P)                    :: b         !< Counter.

   if (self%blocks_number > 0) then
      open(newunit=file_unit,                                                             &
           file=trim(adjustl(basename))//'-proc'//trim(strz(self%mpih%myrank,6))//'.fbd', &
           form='UNFORMATTED',                                                            &
           access='STREAM')
      write(unit=file_unit) self%nv,  self%grid%ni,  self%grid%nj,  self%grid%nk,  self%grid%ngc
      write(unit=file_unit) self%blocks_number
      do b=1, self%blocks_number
         write(unit=file_unit) self%code(b)
         write(unit=file_unit) self%coordinates(1:4,b)
         write(unit=file_unit) self%q(1:self%nv,                                  &
                                      1-self%grid%ngc:self%grid%ni+self%grid%ngc, &
                                      1-self%grid%ngc:self%grid%nj+self%grid%ngc, &
                                      1-self%grid%ngc:self%grid%nk+self%grid%ngc,b)
      enddo
      close(file_unit)
   endif
   endsubroutine save_blocks

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
         if (self%mpih%myrank /= block_to_refine(2,b)) cycle
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
endmodule adam_field_object
