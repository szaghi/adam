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

! ADAM classes, libraries, parameters
use :: adam_maps_object
use :: adam_refinement_plan_object, only : refinement_plan_object
use :: adam_parameters
! ADAM singleton objects
use :: adam_global_grid, only : grid
use :: adam_global_maps, only : maps
use :: adam_global_mpih, only : mpih
! third party modules
use :: finer, only : file_ini
use :: penf
! sdk modules
use :: mpi
! intrinsic modules
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

implicit none
private
public :: field_object

character(len=5), parameter :: INI_SECTION_NAME="field" !< INI (config) file section name containing configs.

type :: field_object
   !< Field class definition.
   ! ADAM objects
   type(maps_object), pointer :: maps=>null() !< The maps.
   ! field data dimensions
   integer(I4P) :: nv=0_I4P               !< Number of field variables.
   integer(I4P) :: nv_pic=0_I4P           !< Number of PIC (Particle In Cell) variables.
   integer(I4P) :: block_weight=0_I4P     !< Block weight, `cells_number * variables_number`.
   integer(I4P) :: block_weight_pic=0_I4P !< Block weight, `particles_number * pic_variables_number`.
   integer(I4P) :: nb=0_I4P               !< Number of all blocks that can be stored.
   integer(I4P) :: blocks_number=0_I4P    !< Number of blocks actually stored.
   integer(I4P) :: np=0_I4P               !< Number of all particles that can be stored for each block.
   ! mesh related data, unrelated to field equations
   integer(I8P), allocatable :: code(:)             !< Morton codes [nb].
   integer(I4P), allocatable :: coordinates(:,:)    !< Coordinates IJKL for each block [4,nb].
   integer(I4P), allocatable :: particles_number(:) !< Number of partcles actually stored in each block [nb].
   real(R8P),    allocatable :: emin(:,:)           !< Coordinates of minimum abscissa of each block [3,nb].
   real(R8P),    allocatable :: emax(:,:)           !< Coordinates of maximum abscissa of each block [3,nb].
   real(R8P),    allocatable :: dxyz(:,:)           !< Space steps of each block [3,nb].
   real(R8P),    allocatable :: x_node(:,:)         !< X node coordinates [0-ngc:ni+ngc,nb].
   real(R8P),    allocatable :: y_node(:,:)         !< Y node coordinates [0-ngc:nj+ngc,nb].
   real(R8P),    allocatable :: z_node(:,:)         !< Z node coordinates [0-ngc:nk+ngc,nb].
   real(R8P),    allocatable :: x_cell(:,:)         !< X cell coordinates [1-ngc:ni+ngc,nb].
   real(R8P),    allocatable :: y_cell(:,:)         !< Y cell coordinates [1-ngc:nj+ngc,nb].
   real(R8P),    allocatable :: z_cell(:,:)         !< Z cell coordinates [1-ngc:nk+ngc,nb].
   ! MPI data, unrelated to field equations
   integer(I4P), allocatable :: req_send_recv(:)          !< MPI request send/receive flags.
   integer(I4P), allocatable :: blocks_numbers(:)         !< Number of blocks actually stored in all processes.
   integer(I4P), allocatable :: refinements_needed(:)     !< Refinements needed of my blocks.
   integer(I4P), allocatable :: refinements_needed_all(:) !< Refinements needed of all blocks.
   integer(I4P), allocatable :: disp_count(:)             !< Displacement of blocks that are received from process.
   integer(I4P)              :: inner_blocks_number=0_I4P !< Number of inner blocks where I need fecs.
   ! field equations data
   real(R8P), allocatable :: q_work(:,:,:,:,:) !< Field cell centered variables, working buffer memory.
   real(R8P), allocatable :: residuals(:)      !< Field residuals, normalized.
   ! grid data replica for easy handling
   integer(I4P), pointer :: ngc=>null() !< Number of ghost cells.
   integer(I4P), pointer :: ni=>null()  !< Number of cells in i direction.
   integer(I4P), pointer :: nj=>null()  !< Number of cells in j direction.
   integer(I4P), pointer :: nk=>null()  !< Number of cells in k direction.
   contains
      ! public methods
      procedure, pass(self) :: adapt                         !< Adapt field accordingly to refine/derefine necessity.
      procedure, pass(self) :: blocks_reorder                !< Reorder blocks indexes in field.
      procedure, pass(self) :: compute_metrics               !< Compute metrics of each block.
      procedure, pass(self) :: compute_normL2_residuals      !< Compute L2 norm of residuals.
      procedure, pass(self) :: description                   !< Return pretty-printed object description.
      procedure, pass(self) :: do_caxis_intersect            !< Return true if a block is intersected by coordinate-axis.
      procedure, pass(self) :: do_cplane_intersect           !< Return true if a block is intersected by coordinate-plane.
      procedure, pass(self) :: do_ray_intersect              !< Return true if a block is intersected by ray.
      procedure, pass(self) :: initialize                    !< Initialize the field.
      procedure, pass(self) :: load_blocks                   !< Load blocks data, used for restarting.
      procedure, pass(self) :: load_from_ini_file            !< Load object data from INI file.
      procedure, pass(self) :: mark_all_blocks               !< Mark all blocks to be refined, derefined, ecc.
      procedure, pass(self) :: mark_sphere                   !< Mark blocks to be refined/derefined by sphere distance.
      procedure, pass(self) :: get_refinements_needed        !< Return gathered refinement flags and displacements.
      procedure, pass(self) :: mpi_gather_refinements_needed !< Gather blocks refinement needed status between MPI processes.
      procedure, pass(self) :: mpi_redistribute              !< Redistribute blocks to processes.
      procedure, pass(self) :: save_blocks                   !< Save blocks data, used for restarting.
      procedure, pass(self) :: update_ghost_local            !< Update ghosts locally.
      procedure, pass(self) :: update_ghost_mpi              !< Update ghosts MPI.
      ! private methods
      procedure, pass(self), private :: derefine1D         !< Derefine blocks, 1D.
      procedure, pass(self), private :: derefine2D         !< Derefine blocks, 2D.
      procedure, pass(self), private :: derefine3D         !< Derefine blocks, 3D.
      procedure, pass(self), private :: refine1D           !< Refine blocks, 1D.
      procedure, pass(self), private :: refine2D           !< Refine blocks, 2D.
      procedure, pass(self), private :: refine3D           !< Refine blocks, 3D.
      procedure, pass(self), private :: update_coordinates !< Update coordinates using the updated data in maps.
endtype field_object

contains
   ! public methods
   subroutine adapt(self, plan, q)
   !< Adapt field accordingly to refine/derefine necessity.
   class(field_object),           intent(inout) :: self !< The field.
   type(refinement_plan_object),  intent(in)    :: plan !< Refinement plan produced by tree%adapt.
   real(R8P),                     intent(inout) :: q(1:,          &
                                                     1-grid%ngc:, &
                                                     1-grid%ngc:, &
                                                     1-grid%ngc:, &
                                                     1:)          !< Field cell centered variables.

   select case(plan%ratio)
   case(2_I4P)
      call self%refine1D(  ratio=plan%ratio, block_to_refine=plan%block_to_refine,     block_refined=plan%block_refined,     q=q)
      call self%derefine1D(ratio=plan%ratio, block_to_derefine=plan%block_to_derefine, block_derefined=plan%block_derefined, q=q)
   case(4_I4P)
      call self%refine2D(  ratio=plan%ratio, block_to_refine=plan%block_to_refine,     block_refined=plan%block_refined,     q=q)
      call self%derefine2D(ratio=plan%ratio, block_to_derefine=plan%block_to_derefine, block_derefined=plan%block_derefined, q=q)
   case(8_I4P)
      call self%refine3D(  ratio=plan%ratio, block_to_refine=plan%block_to_refine,     block_refined=plan%block_refined,     q=q)
      call self%derefine3D(ratio=plan%ratio, block_to_derefine=plan%block_to_derefine, block_derefined=plan%block_derefined, q=q)
   endselect
   endsubroutine adapt

   subroutine blocks_reorder(self, q)
   !< Reorder blocks indexes in field.
   class(field_object), intent(inout) :: self                 !< The field.
   real(R8P),           intent(inout) :: q(1:,              &
                                           1-grid%ngc:,&
                                           1-grid%ngc:,&
                                           1-grid%ngc:,&
                                           1:)                !< Field cell centered variables.
   integer(I4P), allocatable          :: coordinates_new(:,:) !< Temporary coordinates array.
   integer(I8P), allocatable          :: code_new(:)          !< Temporary Morton codes.
   integer(I4P)                       :: b                    !< Counter.

   allocate(coordinates_new(4,self%blocks_number))
   allocate(code_new(self%blocks_number))
   do b=1, self%blocks_number
      self%q_work(:,:,:,:,b) = q(:,:,:,:,self%maps%inner_outer_block_map(b))
      coordinates_new(:,b) = self%coordinates(:,self%maps%inner_outer_block_map(b))
      code_new(b) = self%code(self%maps%inner_outer_block_map(b))
   enddo
   do b=1, self%blocks_number
      q(:,:,:,:,b) = self%q_work(:,:,:,:,b)
      self%coordinates(:,b) = coordinates_new(:,b)
      self%code(b) = code_new(b)
   enddo
   self%inner_blocks_number = self%maps%inner_blocks_number
   call self%compute_metrics
   if (allocated(coordinates_new)) deallocate(coordinates_new)
   if (allocated(code_new)) deallocate(code_new)
   endsubroutine blocks_reorder

   subroutine compute_metrics(self)
   !< Compute metrics of each block.
   class(field_object), intent(inout) :: self  !< The field.
   integer(I4P)                       :: b     !< Counter.

   do b=1, self%blocks_number
      call grid%compute_metrics(coordinates=self%coordinates(:,b),                                         &
                                     emin=self%emin(:,b), emax=self%emax(:,b),                                  &
                                     dx=self%dxyz(1,b), dy=self%dxyz(2,b), dz=self%dxyz(3,b),                   &
                                     x_node=self%x_node(:,b), y_node=self%y_node(:,b), z_node=self%z_node(:,b), &
                                     x_cell=self%x_cell(:,b), y_cell=self%y_cell(:,b), z_cell=self%z_cell(:,b))
   enddo
   endsubroutine compute_metrics

   subroutine compute_normL2_residuals(self, dq, norm)
   !< Compute L2 norm of residuals.
   class(field_object), intent(in)    :: self       !< The field.
   real(R8P),           intent(in)    :: dq(1:,               &
                                            1-grid%ngc:, &
                                            1-grid%ngc:, &
                                            1-grid%ngc:, &
                                            1:)     !< Residuals.
   real(R8P),           intent(inout) :: norm(1:)   !< Residuals norm.
   integer(I4P)                       :: b, i, j, k !< Counter.

   norm   = 0._R8P
   do b=1, self%blocks_number
      do k=1, grid%nk
         do j=1, grid%nj
            do i=1, grid%ni
               norm(:) = norm(:) + dq(:, i, j, k, b)**2
            enddo
         enddo
      enddo
   enddo
   endsubroutine compute_normL2_residuals

   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(field_object), intent(in) :: self             !< The field.
   character(len=:), allocatable   :: desc             !< Description.
   character(len=1), parameter     :: NL=new_line('a') !< New line character.
   integer(I4P)                    :: b                !< Counter.

   desc =       mpih%myrankstr//'field main data'                                                      //NL
   desc = desc//mpih%myrankstr//'  field variables number (nv):     '//trim(str(self%nv              ))//NL
   desc = desc//mpih%myrankstr//'  PIC   variables number (nv_pic): '//trim(str(self%nv_pic          ))//NL
   desc = desc//mpih%myrankstr//'  all blocks number (nb):          '//trim(str(self%nb              ))//NL
   desc = desc//mpih%myrankstr//'  blocks number:                   '//trim(str(self%blocks_number   ))//NL
   desc = desc//mpih%myrankstr//'  block weight:                    '//trim(str(self%block_weight    ))//NL
   desc = desc//mpih%myrankstr//'  block weigh pic                  '//trim(str(self%block_weight_pic))
   do b=1, self%blocks_number
   desc = desc//NL//mpih%myrankstr//'  coordinates b='//trim(strz(b,9)) //'             '//trim(str(self%coordinates(:,b)))
   desc = desc//NL//mpih%myrankstr//'  dxyz        b='//trim(strz(b,9)) //'             '//trim(str(self%dxyz       (:,b)))
   enddo
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

   do_intersect = grid%do_cplane_intersect(emin=self%emin(:,b),         &
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

   subroutine initialize(self, nb, file_parameters, nv, verbose)
   !< Initialize field.
   class(field_object),    intent(inout)        :: self            !< The field.
   type(file_ini),         intent(inout)        :: file_parameters !< INI file handler.
   integer(I4P),           intent(in)           :: nb              !< Number of all blocks that can be stored.
   integer(I4P),           intent(in), optional :: nv              !< Number of field variables.
   logical,                intent(in), optional :: verbose         !< Trigger verbose output.
   logical                                      :: verbose_        !< Trigger verbose output, local variable.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (allocated(self%req_send_recv)) deallocate(self%req_send_recv)
   allocate(self%req_send_recv(0:mpih%procs_number*2-1))
   if (verbose_) call mpih%print_message('field_object%initialize start')
   self%maps => maps
   self%nb = nb
   call self%load_from_ini_file(file_parameters, nv)
   self%block_weight = (grid%ngc+grid%ni+grid%ngc)* &
                       (grid%ngc+grid%nj+grid%ngc)* &
                       (grid%ngc+grid%nk+grid%ngc)*self%nv
   if (self%nb>0) then
      call allocate_variable(var=self%code, ulb=[1,self%nb],&
                             msg=mpih%myrankstr//'field_object%initialize(code) ', verbose=verbose_)
      self%code    = -2_I8P
      self%code(1) = -1_I8P ! first block is assumed to be ADAM
      call allocate_variable(var=self%coordinates, ulb=reshape([1,4, 1,self%nb],[2,2]), &
                             msg=mpih%myrankstr//'field_object%initialize(coordinates) ', verbose=verbose_)
      call allocate_variable(var=self%particles_number, ulb=[1,self%nb],&
                             msg=mpih%myrankstr//'field_object%initialize(particles_number) ', verbose=verbose_)
      call allocate_variable(var=self%emin, ulb=reshape([1,3, 1,self%nb],[2,2]), &
                             msg=mpih%myrankstr//'field_object%initialize(emin) ', verbose=verbose_)
      call allocate_variable(var=self%emax, ulb=reshape([1,3, 1,self%nb],[2,2]), &
                             msg=mpih%myrankstr//'field_object%initialize(emax) ', verbose=verbose_)
      self%emin(:,1) = grid%domain_emin
      self%emax(:,1) = grid%domain_emax
      call allocate_variable(var=self%dxyz, ulb=reshape([1,3, 1,self%nb],[2,2]), &
                             msg=mpih%myrankstr//'field_object%initialize(dxyz) ', verbose=verbose_)
      call allocate_variable(var=self%x_cell,                                         &
                             ulb=reshape([1-grid%ngc,grid%ni+grid%ngc, &
                                          1,self%nb],[2,2]),                          &
                             msg=mpih%myrankstr//'field_object%initialize(x_cell) ', verbose=verbose_)
      call allocate_variable(var=self%y_cell,                                         &
                             ulb=reshape([1-grid%ngc,grid%nj+grid%ngc, &
                                          1,self%nb],[2,2]),                          &
                             msg=mpih%myrankstr//'field_object%initialize(y_cell) ', verbose=verbose_)
      call allocate_variable(var=self%z_cell,                                         &
                             ulb=reshape([1-grid%ngc,grid%nk+grid%ngc, &
                                          1,self%nb],[2,2]),                          &
                             msg=mpih%myrankstr//'field_object%initialize(z_cell) ', verbose=verbose_)
      call allocate_variable(var=self%x_node,                                         &
                             ulb=reshape([0-grid%ngc,grid%ni+grid%ngc, &
                                          1,self%nb],[2,2]),                          &
                             msg=mpih%myrankstr//'field_object%initialize(x_node) ', verbose=verbose_)
      call allocate_variable(var=self%y_node,                                         &
                             ulb=reshape([0-grid%ngc,grid%nj+grid%ngc, &
                                          1,self%nb],[2,2]),                          &
                             msg=mpih%myrankstr//'field_object%initialize(y_node) ', verbose=verbose_)
      call allocate_variable(var=self%z_node,                                         &
                             ulb=reshape([0-grid%ngc,grid%nk+grid%ngc, &
                                          1,self%nb],[2,2]),                          &
                             msg=mpih%myrankstr//'field_object%initialize(z_node) ', verbose=verbose_)
      ! if (present(q)) then
      !    call allocate_variable(var=q,                                                   &
      !                           ulb=reshape([1,self%nv,                                  &
      !                                        1-grid%ngc,grid%ni+grid%ngc, &
      !                                        1-grid%ngc,grid%nj+grid%ngc, &
      !                                        1-grid%ngc,grid%nk+grid%ngc, &
      !                                        1,self%nb],[2,5]),                          &
      !                           msg=mpih%myrankstr//'field_object%initialize(q) ', verbose=verbose_)
      ! endif
      call allocate_variable(var=self%q_work,                                         &
                             ulb=reshape([1,self%nv,                                  &
                                          1-grid%ngc,grid%ni+grid%ngc, &
                                          1-grid%ngc,grid%nj+grid%ngc, &
                                          1-grid%ngc,grid%nk+grid%ngc, &
                                          1,self%nb],[2,5]),                          &
                             msg=mpih%myrankstr//'field_object%initialize(q_work) ', verbose=verbose_)
      call allocate_variable(var=self%residuals, &
                             ulb=[1,self%nv],    &
                             msg=mpih%myrankstr//'field_object%initialize(residuals) ', verbose=verbose_)
      self%q_work    = 0._R8P
      self%residuals = 0._R8P
      self%block_weight_pic = 0_I4P
   endif
   call allocate_variable(var=self%blocks_numbers, ulb=[0,mpih%procs_number-1], &
                          msg=mpih%myrankstr//'field_object%initialize(blocks_numbers) ', verbose=verbose_)
   self%ngc => grid%ngc
   self%ni  => grid%ni
   self%nj  => grid%nj
   self%nk  => grid%nk
   call self%update_coordinates
   call self%compute_metrics
   if (verbose_) print '(A)', self%description()
   if (verbose_) call mpih%print_message('field_object%initialize finish')
   endsubroutine initialize

   subroutine load_blocks(self, basename, q)
   !< Load blocks data, used for restarting.
   !<
   !< Note: blocks memory must be already initialized with enough memory (proper nv,ni,nj,nk,ngc and nb>=blocks number).
   class(field_object), intent(inout) :: self                            !< The field.
   character(*),        intent(in)    :: basename                        !< Output base name.
   real(R8P),           intent(inout) :: q(1:,              &
                                           1-grid%ngc:,&
                                           1-grid%ngc:,&
                                           1-grid%ngc:,&
                                           1:)                           !< Field cell centered variables.
   character(:), allocatable          :: filename                        !< Output file name.
   integer(I4P)                       :: file_unit                       !< Output file unit.
   logical                            :: file_exist                      !< Flag to check file's existance.
   integer(I4P)                       :: blocks_number                   !< Blocks number.
   integer(I4P)                       :: nv, ni, nj, nk, ngc, nv_pic, np !< Dimensions.
   integer(I4P)                       :: b                               !< Counter.

   filename = trim(adjustl(basename))//'-proc'//trim(strz(mpih%myrank,6))//'.fbd'
   inquire(file=filename, exist=file_exist)
   if (file_exist) then
      call mpih%print_message('field_object%load_blocks from file '//filename)
      open(newunit=file_unit, file=filename, form='UNFORMATTED', access='STREAM')
      read(unit=file_unit) nv, ni, nj, nk, ngc, nv_pic, np
      if (nv==self%nv.and.ni==grid%ni.and.nj==grid%nj.and.nk==grid%nk.and.ngc==grid%ngc) then
         read(unit=file_unit) blocks_number
         call mpih%print_message('field blocks number '//trim(str(blocks_number)))
         if (blocks_number <= self%nb) then
            self%blocks_number = blocks_number
            do b=1, self%blocks_number
               read(unit=file_unit) self%code(b)
               read(unit=file_unit) self%coordinates(1:4,b)
               read(unit=file_unit) q(1:self%nv,                                  &
                                      1-grid%ngc:grid%ni+grid%ngc, &
                                      1-grid%ngc:grid%nj+grid%ngc, &
                                      1-grid%ngc:grid%nk+grid%ngc,b)
            enddo
            call self%compute_metrics
         else
            call mpih%abort(error_code=-102, msg='ERROR: blocks number to read "'//trim(str(blocks_number))//&
                                                      '" is greater than blocks allocated "'//trim(str(self%nb))//'"!')
         endif
      else
         call mpih%abort(error_code=-103,                                                            &
                              msg='ERROR: blocks dimensions to read "'//trim(str([nv, ni, nj, nk, ngc]))//&
                                  '" are different than blocks allocated "'//trim(str([self%nv,           &
                                                                                       grid%ni,      &
                                                                                       grid%nj,      &
                                                                                       grid%nk,      &
                                                                                       grid%ngc]))//'"!')
      endif
      close(file_unit)
      call mpih%print_message('field_object%load_blocks from file '//filename//' completed')
   else
      call mpih%abort(error_code=-104, msg='WARNING: file "'//filename//'" does not exist!')
   endif
   endsubroutine load_blocks

   subroutine load_from_ini_file(self, file_parameters, nv, go_on_fail)
   !< Load object data from INI file.
   class(field_object), intent(inout)        :: self            !< The field.
   type(file_ini),      intent(inout)        :: file_parameters !< INI file handler.
   integer(I4P),        intent(in), optional :: nv              !< Number of field variables.
   logical,             intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                   :: go_on_fail_     !< Go on if load fails.
   integer(I4P)                              :: error           !< Error status.
   integer(I4P)                              :: buff_I4P        !< I4P buffer.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   if (present(nv)) then
      self%nv = nv
   else
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='nv', val=buff_I4P, error=error)
      if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(nv)')
      self%nv = buff_I4P
   endif
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

      associate (emin=>self%emin(:,b), emax=>self%emax(:,b), ni=>grid%ni, nj=>grid%nj, nk=>grid%nk)
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
   allocate(recv_count(0:mpih%procs_number - 1))
   call MPI_ALLGATHER(self%blocks_number, 1_I4P, MPI_INTEGER, recv_count, 1_I4P, MPI_INTEGER, MPI_COMM_WORLD, mpih%error)

   ! computing displacement counts
   if (allocated(self%disp_count)) deallocate(self%disp_count)
   allocate(self%disp_count(0:mpih%procs_number - 1))
   self%disp_count = 0_I4P
   do p=1, mpih%procs_number - 1
      self%disp_count(p) = self%disp_count(p-1) + recv_count(p-1)
   enddo

   if (allocated(self%refinements_needed_all)) deallocate(self%refinements_needed_all)
   allocate(self%refinements_needed_all(sum(recv_count, dim=1)))
   call MPI_ALLGATHERV(self%refinements_needed, self%blocks_number, MPI_INTEGER, &
                       self%refinements_needed_all, recv_count, self%disp_count, MPI_INTEGER, MPI_COMM_WORLD, mpih%error)
   endsubroutine mpi_gather_refinements_needed

   subroutine get_refinements_needed(self, flags, disp)
   !< Return gathered refinement flags and process displacements.
   !< Provides an accessor so callers do not reach into field internals.
   !< Call mpi_gather_refinements_needed before invoking this procedure.
   class(field_object),       intent(in)  :: self     !< The field.
   integer(I4P), allocatable, intent(out) :: flags(:) !< Refinements needed of all blocks.
   integer(I4P), allocatable, intent(out) :: disp(:)  !< Displacement of received blocks per process.

   flags = self%refinements_needed_all
   disp  = self%disp_count
   endsubroutine get_refinements_needed

   subroutine mpi_redistribute(self, q)
   !< Redistribute blocks to processes.
   !< @TODO: Morton codes are not yet redistributed, must be fixed.
   class(field_object), intent(inout) :: self                   !< The field.
   real(R8P),           intent(inout) :: q(1:,              &
                                           1-grid%ngc:,&
                                           1-grid%ngc:,&
                                           1-grid%ngc:,&
                                           1:)                  !< Field cell centered variables.
   real(R8P),    allocatable          :: send_buffer(:)         !< Send buffer of field cell centered variables.
   real(R8P),    allocatable          :: recv_buffer(:)         !< Recv buffer of field cell centered variables.
   integer(I8P)                       :: send_size, send_offset !< Total size of send buffer.
   integer(I8P)                       :: recv_size, recv_offset !< Total size of recv buffer.
   integer(I4P)                       :: n_keep                 !< Number of keept blocks.
   integer(I4P)                       :: b, bi, p               !< Counter.
   integer(I4P)                       :: ptr_start, ptr_end     !< Counter.
   integer(I4P)                       :: n_recv, n_send         !< Counter.
   integer(I4P), allocatable          :: req_recv(:)            !< MPI request receive flags.
   integer(I4P)                       :: bwt                    !< Block weight total.

   allocate(req_recv(0:mpih%procs_number-1))
   req_recv = MPI_REQUEST_NULL
   associate(bw=>self%block_weight, bw_pic=>self%block_weight_pic)
   bwt = bw + bw_pic
   send_size = 0_I8P ; if (allocated(self%maps%comm_map_send)) send_size = size(self%maps%comm_map_send, dim=1) * bwt
   recv_size = 0_I8P ; if (allocated(self%maps%comm_map_recv)) recv_size = size(self%maps%comm_map_recv, dim=1) * bwt
   n_keep    = 0_I8P ; if (allocated(self%maps%local_map    )) n_keep    = size(self%maps%local_map    , dim=1)
   if (send_size > 0_I8P) allocate(send_buffer(send_size))
   if (recv_size > 0_I8P) allocate(recv_buffer(recv_size))

   if (send_size > 0_I8P) then
      send_offset = 1
      do b=1, size(self%maps%comm_map_send, dim=1)
         bi = self%maps%comm_map_send(b)
         send_buffer(send_offset:send_offset+bw-1) = reshape(q(:,:,:,:,bi),[bw])
         send_offset = send_offset + bw
      enddo
   endif

   do p=0, mpih%procs_number - 1_I4P
      ptr_start = self%maps%comm_map_recv_ptr(p)   * bwt + 1
      ptr_end   = self%maps%comm_map_recv_ptr(p+1) * bwt
      n_recv    = ptr_end - ptr_start + 1
      if (n_recv > 0) then
         call MPI_IRECV(recv_buffer(ptr_start), n_recv, MPI_REAL8, p, 100, MPI_COMM_WORLD, req_recv(p), mpih%error)
      endif
   enddo

   do p=0, mpih%procs_number - 1_I4P
      ptr_start = self%maps%comm_map_send_ptr(p)   * bwt + 1
      ptr_end   = self%maps%comm_map_send_ptr(p+1) * bwt
      n_send    = ptr_end - ptr_start + 1
      if (n_send > 0) then
         call MPI_SEND(send_buffer(ptr_start), n_send, MPI_REAL8, p, 100, MPI_COMM_WORLD, mpih%error)
      endif
   enddo

   call MPI_WAITALL(mpih%procs_number, req_recv, MPI_STATUSES_IGNORE, mpih%error)

   if (recv_size > 0_I8P) then
      recv_offset = 1
      do b=1, size(self%maps%comm_map_recv, dim=1)
          bi = self%maps%comm_map_recv(b)
          self%q_work(:,:,:,:,bi) = reshape(recv_buffer(recv_offset:recv_offset + bw -1),&
                                            [self%nv,                                    &
                                             grid%ngc+grid%ni+grid%ngc,   &
                                             grid%ngc+grid%nj+grid%ngc,   &
                                             grid%ngc+grid%nk+grid%ngc])
          recv_offset = recv_offset + bw
      enddo
   endif
   do b=1, n_keep
      self%q_work(:,:,:,:,self%maps%local_map(b,1)) = q(:,:,:,:,self%maps%local_map(b,2))
   enddo
   self%blocks_number = n_keep  + recv_size / bwt
   q(:,:,:,:,1:self%blocks_number) = self%q_work(:,:,:,:,1:self%blocks_number)
   call self%update_coordinates
   call self%maps%get_block_layout(code=self%code(1:self%blocks_number))
   call self%compute_metrics
   endassociate
   endsubroutine mpi_redistribute

   subroutine save_blocks(self, basename, q)
   !< Save blocks data, used for restarting.
   class(field_object), intent(in) :: self      !< The field.
   character(*),        intent(in) :: basename  !< Output base name.
   real(R8P),           intent(in) :: q(1:,              &
                                        1-grid%ngc:,&
                                        1-grid%ngc:,&
                                        1-grid%ngc:,&
                                        1:)     !< Field cell centered variables.
   integer(I4P)                    :: file_unit !< Output file unit.
   integer(I4P)                    :: b         !< Counter.

   if (self%blocks_number > 0) then
      open(newunit=file_unit,                                                             &
           file=trim(adjustl(basename))//'-proc'//trim(strz(mpih%myrank,6))//'.fbd', &
           form='UNFORMATTED',                                                            &
           access='STREAM')
      write(unit=file_unit) self%nv, grid%ni, grid%nj, grid%nk, grid%ngc, self%nv_pic, self%np
      write(unit=file_unit) self%blocks_number
      do b=1, self%blocks_number
         write(unit=file_unit) self%code(b)
         write(unit=file_unit) self%coordinates(1:4,b)
         write(unit=file_unit) q(1:self%nv,                                  &
                                 1-grid%ngc:grid%ni+grid%ngc, &
                                 1-grid%ngc:grid%nj+grid%ngc, &
                                 1-grid%ngc:grid%nk+grid%ngc,b)
      enddo
      close(file_unit)
   endif
   endsubroutine save_blocks

   subroutine update_ghost_local(self, q)
   !< Update (local) ghost cells.
   class(field_object), intent(in)    :: self              !< The base backend.
   real(R8P),           intent(inout) :: q(1:,              &
                                           1-grid%ngc:,&
                                           1-grid%ngc:,&
                                           1-grid%ngc:,&
                                           1:)             !< Field component to be updated.
   integer(I4P)                       :: ic, jc, kc, mf, v !< Counter.
   integer(I4P)                       :: b_recv            !< Index of receiving block.
   integer(I4P)                       :: b_send            !< Index of sending block.
   integer(I4P)                       :: i_recv            !< I recv index.
   integer(I4P)                       :: j_recv            !< J recv index.
   integer(I4P)                       :: k_recv            !< K recv index.
   integer(I4P)                       :: i_send            !< I send index.
   integer(I4P)                       :: j_send            !< J send index.
   integer(I4P)                       :: k_send            !< K send index.
   integer(I4P)                       :: one_or_eight      !< Flag triggering 8 cells mean.

   if (.not.allocated(self%maps%local_map_ghost_cell)) return
   associate(local_map_ghost_cell=>self%maps%local_map_ghost_cell)
   do mf=1, size(local_map_ghost_cell, dim=1)
      do v=1, size(q, dim=1)
         b_send       = local_map_ghost_cell(mf,1)
         b_recv       = local_map_ghost_cell(mf,2)
         i_send       = local_map_ghost_cell(mf,3)
         j_send       = local_map_ghost_cell(mf,4)
         k_send       = local_map_ghost_cell(mf,5)
         i_recv       = local_map_ghost_cell(mf,6)
         j_recv       = local_map_ghost_cell(mf,7)
         k_recv       = local_map_ghost_cell(mf,8)
         one_or_eight = local_map_ghost_cell(mf,9)
         if (one_or_eight==1) then
            q(v,i_recv,j_recv,k_recv,b_recv) = q(v,i_send,j_send,k_send,b_send)
         else
            q(v,i_recv,j_recv,k_recv,b_recv) = 0._R8P
            do kc=0,1 ; do jc=0,1 ; do ic=0,1
               q(v,i_recv,j_recv,k_recv,b_recv) = q(v,i_recv,   j_recv,   k_recv,   b_recv) + &
                                                  q(v,i_send+ic,j_send+jc,k_send+kc,b_send)
            enddo ; enddo ; enddo
            q(v,i_recv,j_recv,k_recv,b_recv) = q(v,i_recv,j_recv,k_recv,b_recv) * 0.125_R8P
         endif
      enddo
   enddo
   endassociate
   endsubroutine update_ghost_local

   subroutine update_ghost_mpi(self, q, step)
   !< Update ghost cells within other processes.
   class(field_object), intent(inout)        :: self                               !< The base backend.
   real(R8P),           intent(inout)        :: q(1:,              &
                                                  1-grid%ngc:,&
                                                  1-grid%ngc:,&
                                                  1-grid%ngc:,&
                                                  1:)                              !< Field component to be updated.
   integer(I4P),        intent(in), optional :: step                               !< Step to be perfordmed in asyncronous comp.
   logical                                   :: do_step(3)                         !< Steps performed in async comp.
   integer(I4P)                              :: p                                  !< Counter.
   integer(I4P)                              :: ptr_start, ptr_end                 !< Counter.
   integer(I4P)                              :: n_recv, n_send                     !< Counter.
   integer(I4P)                              :: ic, jc, kc, sf                     !< Counter.
   integer(I4P)                              :: b_send,i_send,j_send,k_send,v_send !< Send indexes.
   integer(I4P)                              :: c_recv                             !< Counter.
   integer(I4P)                              :: one_or_eight                       !< Flag triggering 8 cells mean.
   integer(I4P)                              :: rf                                 !< Counter.
   integer(I4P)                              :: b_recv,i_recv,j_recv,k_recv,v_recv !< Receive indexes.
   integer(I4P)                              :: c_send                             !< Counter.

   associate(procs_number=>mpih%procs_number,                       &
             error=>mpih%error,                                     &
             req_send_recv=>self%req_send_recv,                     &
             comm_map_send_ptr_ghost=>self%maps%comm_map_send_ptr_ghost, &
             comm_map_recv_ptr_ghost=>self%maps%comm_map_recv_ptr_ghost, &
             recv_buffer_ghost=>self%maps%recv_buffer_ghost,             &
             send_buffer_ghost=>self%maps%send_buffer_ghost,             &
             ngc=>grid%ngc)
   do_step = .true.
   if (present(step)) then
      do_step = .false.
      do_step(step) = .true.
   endif

   if (do_step(1)) then
      req_send_recv = MPI_REQUEST_NULL
      if (allocated(self%maps%comm_map_send_ghost_cell)) then
         do sf=1, size(self%maps%comm_map_send_ghost_cell, dim=1)
            b_send       = self%maps%comm_map_send_ghost_cell(sf,1)
            i_send       = self%maps%comm_map_send_ghost_cell(sf,2)
            j_send       = self%maps%comm_map_send_ghost_cell(sf,3)
            k_send       = self%maps%comm_map_send_ghost_cell(sf,4)
            v_send       = self%maps%comm_map_send_ghost_cell(sf,5)
            c_recv       = self%maps%comm_map_send_ghost_cell(sf,6)
            one_or_eight = self%maps%comm_map_send_ghost_cell(sf,7)
            if (one_or_eight==1) then
               send_buffer_ghost(c_recv) = q(v_send,i_send,j_send,k_send,b_send)
            else
               send_buffer_ghost(c_recv) = 0._R8P
               do kc=0,1 ; do jc=0,1 ; do ic=0,1
                  send_buffer_ghost(c_recv) = send_buffer_ghost(c_recv) + q(v_send,i_send+ic,j_send+jc,k_send+kc,b_send)
               enddo ; enddo ; enddo
               send_buffer_ghost(c_recv) = send_buffer_ghost(c_recv) * 0.125_R8P
            endif
         enddo
      endif
   endif

   if (do_step(2)) then
      ! receive
      do p=0, procs_number - 1_I4P
         ptr_start = comm_map_recv_ptr_ghost(p) + 1
         ptr_end   = comm_map_recv_ptr_ghost(p+1)
         n_recv    = ptr_end - ptr_start + 1
         if (n_recv > 0) then
            call MPI_IRECV(recv_buffer_ghost(ptr_start), n_recv, MPI_REAL8, p, 100, MPI_COMM_WORLD, req_send_recv(p), error)
         endif
      enddo
      ! send
      do p=0, procs_number - 1_I4P
         ptr_start = comm_map_send_ptr_ghost(p) + 1
         ptr_end   = comm_map_send_ptr_ghost(p+1)
         n_send    = ptr_end - ptr_start + 1
         if (n_send > 0) then
            call MPI_ISEND(send_buffer_ghost(ptr_start), n_send, MPI_REAL8, p, 100, MPI_COMM_WORLD, &
                           req_send_recv(p+procs_number), error)
         endif
      enddo
   endif

   if (do_step(3)) then
      call MPI_WAITALL(procs_number * 2, req_send_recv, MPI_STATUSES_IGNORE, error)
      call MPI_Barrier(MPI_COMM_WORLD, error)
      !RIMETTERE SENZA
      if (allocated(self%maps%comm_map_recv_ghost_cell)) then
         do rf=1, size(self%maps%comm_map_recv_ghost_cell, dim=1)
            c_send = self%maps%comm_map_recv_ghost_cell(rf,1)
            b_recv = self%maps%comm_map_recv_ghost_cell(rf,2)
            i_recv = self%maps%comm_map_recv_ghost_cell(rf,3)
            j_recv = self%maps%comm_map_recv_ghost_cell(rf,4)
            k_recv = self%maps%comm_map_recv_ghost_cell(rf,5)
            v_recv = self%maps%comm_map_recv_ghost_cell(rf,6)
            q(v_recv,i_recv,j_recv,k_recv,b_recv) = recv_buffer_ghost(c_send)
         enddo
      endif
   endif
   call MPI_Barrier(MPI_COMM_WORLD, error)
   endassociate
   endsubroutine update_ghost_mpi

   ! private methods
   subroutine derefine1D(self, ratio, block_to_derefine, block_derefined, q)
   !< Derefine blocks.
   !<
   !< Note: blocks number is not updated: mpi redistribute does it. This is dangerous...
   class(field_object),       intent(inout) :: self                 !< The field.
   integer(I4P),              intent(in)    :: ratio                !< Refinement ratio.
   integer(I8P), allocatable, intent(in)    :: block_to_derefine(:) !< List of blocks to be derefined.
   integer(I8P), allocatable, intent(in)    :: block_derefined(:,:) !< List of derefined blocks with Morton code.
   real(R8P),                 intent(inout) :: q(1:,              &
                                                 1-grid%ngc:,&
                                                 1-grid%ngc:,&
                                                 1-grid%ngc:,&
                                                 1:)                !< Field cell centered variables.
   integer(I4P)                             :: b, ib                !< Counter.
   integer(I4P)                             :: ic(ratio)            !< Counter.
   integer(I4P)                             :: iii                  !< Counter.
   integer(I4P)                             :: i, j, k              !< Counter.

   associate(ni=>grid%ni, nj=>grid%nj, nk=>grid%nk, q_work=>self%q_work)
   if (allocated(block_derefined)) then
      do b=1, size(block_derefined, dim=2)
         ib = block_derefined(2,b)
         do i=1, ratio
            ic(i) = block_to_derefine((b-1)*ratio+i)
         enddo
         do k=1, nk
            do j=1, nj
               do i=1, ni/2
                  iii = (i-1)*2+1
                  q_work(:,i,     j,k,ib) = (q(:,iii,j,k,ic(1)) + q(:,iii+1,j,k,ic(1))) / 2._R8P
                  q_work(:,i+ni/2,j,k,ib) = (q(:,iii,j,k,ic(2)) + q(:,iii+1,j,k,ic(2))) / 2._R8P
               enddo
            enddo
         enddo
         q(:,1:ni,1:nj,1:nk,ib) = q_work(:,1:ni,1:nj,1:nk,ib)
         self%code(ib) = block_derefined(1,b)
      enddo
   endif
   endassociate
   endsubroutine derefine1D

   subroutine derefine2D(self, ratio, block_to_derefine, block_derefined, q)
   !< Derefine blocks.
   !<
   !< Note: blocks number is not updated: mpi redistribute does it. This is dangerous...
   class(field_object),       intent(inout) :: self                 !< The field.
   integer(I4P),              intent(in)    :: ratio                !< Refinement ratio.
   integer(I8P), allocatable, intent(in)    :: block_to_derefine(:) !< List of blocks to be derefined.
   integer(I8P), allocatable, intent(in)    :: block_derefined(:,:) !< List of derefined blocks with Morton code.
   real(R8P),                 intent(inout) :: q(1:,              &
                                                 1-grid%ngc:,&
                                                 1-grid%ngc:,&
                                                 1-grid%ngc:,&
                                                 1:)                !< Field cell centered variables.
   integer(I4P)                             :: b, ib                !< Counter.
   integer(I4P)                             :: ic(ratio)            !< Counter.
   integer(I4P)                             :: iii, jjj             !< Counter.
   integer(I4P)                             :: i, j, k              !< Counter.

   associate(ni=>grid%ni, nj=>grid%nj, nk=>grid%nk, q_work=>self%q_work)
   if (allocated(block_derefined)) then
      do b=1, size(block_derefined, dim=2)
         ib = block_derefined(2,b)
         do i=1, ratio
            ic(i) = block_to_derefine((b-1)*ratio+i)
         enddo
         do k=1, nk
            do j=1, nj/2
               do i=1, ni/2
                  jjj = (j-1)*2+1
                  iii = (i-1)*2+1
                  q_work(:,i,     j,     k,ib) = (q(:,iii,jjj,  k,ic(1)) + q(:,iii+1,jjj,  k,ic(1)) + &
                                                  q(:,iii,jjj+1,k,ic(1)) + q(:,iii+1,jjj+1,k,ic(1))) / 4._R8P
                  q_work(:,i+ni/2,j,     k,ib) = (q(:,iii,jjj,  k,ic(2)) + q(:,iii+1,jjj,  k,ic(2)) + &
                                                  q(:,iii,jjj+1,k,ic(2)) + q(:,iii+1,jjj+1,k,ic(2))) / 4._R8P
                  q_work(:,i,     j+nj/2,k,ib) = (q(:,iii,jjj,  k,ic(3)) + q(:,iii+1,jjj,  k,ic(3)) + &
                                                  q(:,iii,jjj+1,k,ic(3)) + q(:,iii+1,jjj+1,k,ic(3))) / 4._R8P
                  q_work(:,i+ni/2,j+nj/2,k,ib) = (q(:,iii,jjj,  k,ic(4)) + q(:,iii+1,jjj,  k,ic(4)) + &
                                                  q(:,iii,jjj+1,k,ic(4)) + q(:,iii+1,jjj+1,k,ic(4))) / 4._R8P
               enddo
            enddo
         enddo
         q(:,1:ni,1:nj,1:nk,ib) = q_work(:,1:ni,1:nj,1:nk,ib)
         self%code(ib) = block_derefined(1,b)
      enddo
   endif
   endassociate
   endsubroutine derefine2D

   subroutine derefine3D(self, ratio, block_to_derefine, block_derefined, q)
   !< Derefine blocks.
   !<
   !< Note: blocks number is not updated: mpi redistribute does it. This is dangerous...
   class(field_object),       intent(inout) :: self                 !< The field.
   integer(I4P),              intent(in)    :: ratio                !< Refinement ratio.
   integer(I8P), allocatable, intent(in)    :: block_to_derefine(:) !< List of blocks to be derefined.
   integer(I8P), allocatable, intent(in)    :: block_derefined(:,:) !< List of derefined blocks with Morton code.
   real(R8P),                 intent(inout) :: q(1:,              &
                                                 1-grid%ngc:,&
                                                 1-grid%ngc:,&
                                                 1-grid%ngc:,&
                                                 1:)                !< Field cell centered variables.
   integer(I4P)                             :: b, ib                !< Counter.
   integer(I4P)                             :: ic(ratio)            !< Counter.
   integer(I4P)                             :: iii, jjj, kkk        !< Counter.
   integer(I4P)                             :: i, j, k              !< Counter.

   associate(ni=>grid%ni, nj=>grid%nj, nk=>grid%nk, q_work=>self%q_work)
   if (allocated(block_derefined)) then
      do b=1, size(block_derefined, dim=2)
         ib = block_derefined(2,b)
         do i=1, ratio
            ic(i) = block_to_derefine((b-1)*ratio+i)
         enddo
         do k=1, nk/2
            do j=1, nj/2
               do i=1, ni/2
                  kkk = (k - 1) * 2 + 1
                  jjj = (j - 1) * 2 + 1
                  iii = (i - 1) * 2 + 1
                  q_work(:,i,     j,     k     ,ib) = (q(:,iii,jjj,  kkk  ,ic(1)) + q(:,iii+1,jjj,  kkk  ,ic(1)) + &
                                                       q(:,iii,jjj+1,kkk  ,ic(1)) + q(:,iii+1,jjj+1,kkk  ,ic(1)) + &
                                                       q(:,iii,jjj,  kkk+1,ic(1)) + q(:,iii+1,jjj,  kkk+1,ic(1)) + &
                                                       q(:,iii,jjj+1,kkk+1,ic(1)) + q(:,iii+1,jjj+1,kkk+1,ic(1))) / 8._R8P
                  q_work(:,i+ni/2,j,     k     ,ib) = (q(:,iii,jjj,  kkk  ,ic(2)) + q(:,iii+1,jjj,  kkk  ,ic(2)) + &
                                                       q(:,iii,jjj+1,kkk  ,ic(2)) + q(:,iii+1,jjj+1,kkk  ,ic(2)) + &
                                                       q(:,iii,jjj,  kkk+1,ic(2)) + q(:,iii+1,jjj,  kkk+1,ic(2)) + &
                                                       q(:,iii,jjj+1,kkk+1,ic(2)) + q(:,iii+1,jjj+1,kkk+1,ic(2))) / 8._R8P
                  q_work(:,i,     j+nj/2,k     ,ib) = (q(:,iii,jjj,  kkk  ,ic(3)) + q(:,iii+1,jjj,  kkk  ,ic(3)) + &
                                                       q(:,iii,jjj+1,kkk  ,ic(3)) + q(:,iii+1,jjj+1,kkk  ,ic(3)) + &
                                                       q(:,iii,jjj,  kkk+1,ic(3)) + q(:,iii+1,jjj,  kkk+1,ic(3)) + &
                                                       q(:,iii,jjj+1,kkk+1,ic(3)) + q(:,iii+1,jjj+1,kkk+1,ic(3))) / 8._R8P
                  q_work(:,i+ni/2,j+nj/2,k     ,ib) = (q(:,iii,jjj,  kkk  ,ic(4)) + q(:,iii+1,jjj,  kkk  ,ic(4)) + &
                                                       q(:,iii,jjj+1,kkk  ,ic(4)) + q(:,iii+1,jjj+1,kkk  ,ic(4)) + &
                                                       q(:,iii,jjj,  kkk+1,ic(4)) + q(:,iii+1,jjj,  kkk+1,ic(4)) + &
                                                       q(:,iii,jjj+1,kkk+1,ic(4)) + q(:,iii+1,jjj+1,kkk+1,ic(4))) / 8._R8P
                  q_work(:,i,     j,     k+nk/2,ib) = (q(:,iii,jjj,  kkk  ,ic(5)) + q(:,iii+1,jjj,  kkk  ,ic(5)) + &
                                                       q(:,iii,jjj+1,kkk  ,ic(5)) + q(:,iii+1,jjj+1,kkk  ,ic(5)) + &
                                                       q(:,iii,jjj,  kkk+1,ic(5)) + q(:,iii+1,jjj,  kkk+1,ic(5)) + &
                                                       q(:,iii,jjj+1,kkk+1,ic(5)) + q(:,iii+1,jjj+1,kkk+1,ic(5))) / 8._R8P
                  q_work(:,i+ni/2,j,     k+nk/2,ib) = (q(:,iii,jjj,  kkk  ,ic(6)) + q(:,iii+1,jjj,  kkk  ,ic(6)) + &
                                                       q(:,iii,jjj+1,kkk  ,ic(6)) + q(:,iii+1,jjj+1,kkk  ,ic(6)) + &
                                                       q(:,iii,jjj,  kkk+1,ic(6)) + q(:,iii+1,jjj,  kkk+1,ic(6)) + &
                                                       q(:,iii,jjj+1,kkk+1,ic(6)) + q(:,iii+1,jjj+1,kkk+1,ic(6))) / 8._R8P
                  q_work(:,i,     j+nj/2,k+nk/2,ib) = (q(:,iii,jjj,  kkk  ,ic(7)) + q(:,iii+1,jjj,  kkk  ,ic(7)) + &
                                                       q(:,iii,jjj+1,kkk  ,ic(7)) + q(:,iii+1,jjj+1,kkk  ,ic(7)) + &
                                                       q(:,iii,jjj,  kkk+1,ic(7)) + q(:,iii+1,jjj,  kkk+1,ic(7)) + &
                                                       q(:,iii,jjj+1,kkk+1,ic(7)) + q(:,iii+1,jjj+1,kkk+1,ic(7))) / 8._R8P
                  q_work(:,i+ni/2,j+nj/2,k+nk/2,ib) = (q(:,iii,jjj,  kkk  ,ic(8)) + q(:,iii+1,jjj,  kkk  ,ic(8)) + &
                                                       q(:,iii,jjj+1,kkk  ,ic(8)) + q(:,iii+1,jjj+1,kkk  ,ic(8)) + &
                                                       q(:,iii,jjj,  kkk+1,ic(8)) + q(:,iii+1,jjj,  kkk+1,ic(8)) + &
                                                       q(:,iii,jjj+1,kkk+1,ic(8)) + q(:,iii+1,jjj+1,kkk+1,ic(8))) / 8._R8P
               enddo
            enddo
         enddo
         q(:,1:ni,1:nj,1:nk,ib) = q_work(:,1:ni,1:nj,1:nk,ib)
         self%code(ib) = block_derefined(1,b)
      enddo
   endif
   endassociate
   endsubroutine derefine3D

   subroutine refine1D(self, ratio, block_to_refine, block_refined, q)
   !< Refine blocks.
   !<
   !< Note: blocks number is not updated: mpi redistribute does it. This is dangerous...
   class(field_object),       intent(inout) :: self                 !< The field.
   integer(I4P),              intent(in)    :: ratio                !< Refinement ratio.
   integer(I8P), allocatable, intent(in)    :: block_to_refine(:,:) !< List of blocks to be refined.
   integer(I8P), allocatable, intent(in)    :: block_refined(:,:)   !< List of refined blocks with Morton code.
   real(R8P),                 intent(inout) :: q(1:,              &
                                                 1-grid%ngc:,&
                                                 1-grid%ngc:,&
                                                 1-grid%ngc:,&
                                                 1:)                !< Field cell centered variables.
   integer(I4P)                             :: b, i, j, k           !< Spatial counter.
   integer(I4P)                             :: ib, ic, ic_local     !< Counter.
   integer(I4P)                             :: i_fine               !< Counter.
   integer(I4P)                             :: i_delta              !< Counter.
   integer(I4P)                             :: ic1, ic2             !< Counter.

   associate(ni=>grid%ni, nj=>grid%nj, nk=>grid%nk, q_work=>self%q_work)
   if (allocated(block_to_refine)) then
      do b=1, size(block_to_refine, dim=2)
         if (mpih%myrank /= block_to_refine(2,b)) cycle
         ib = block_to_refine(1,b)

         q_work(:,:,:,:,ib) = q(:,:,:,:,ib)

         do ic_local=1, ratio
            ic = block_refined(2,(b-1)*ratio+ic_local)
            ic1 = mod(ic_local - 1, 2)
            do k=1,nk
               do j=1,nj
                  do i=1+ni/2*ic1,ni/2+ni/2*ic1
                     i_fine = mod(i - 1, ni/2) * 2 + 1
                     q(:,i_fine:i_fine+1,j,k,ic) = 0._R8P
                     do i_delta=0,1
                     ! q(:,i_fine,  j,k,ic) = q(:,i_fine,  j,k,ic) + 0.5_R8P * q_work(:,i+i_delta-1,j,k,ib)
                     ! q(:,i_fine+1,j,k,ic) = q(:,i_fine+1,j,k,ic) + 0.5_R8P * q_work(:,i+i_delta,  j,k,ib)
                     q(:,i_fine,  j,k,ic) = q(:,i_fine  ,j,k,ic) + (0.25_R8P + i_delta * 0.5_R8P) * q_work(:,i+i_delta-1,j,k,ib)
                     q(:,i_fine+1,j,k,ic) = q(:,i_fine+1,j,k,ic) + (0.75_R8P - i_delta * 0.5_R8P) * q_work(:,i+i_delta  ,j,k,ib)
                     enddo
                  enddo
               enddo
            enddo
         enddo

         ic1 = block_refined(2,(b-1)*ratio+1)
         ic2 = block_refined(2,(b-1)*ratio+2)
         self%code(ic1) = block_refined(1,(b-1)*ratio+1)
         self%code(ic2) = block_refined(1,(b-1)*ratio+2)
      enddo
   endif
   endassociate
   endsubroutine refine1D

   subroutine refine2D(self, ratio, block_to_refine, block_refined, q)
   !< Refine blocks.
   !<
   !< Note: blocks number is not updated: mpi redistribute does it. This is dangerous...
   class(field_object),       intent(inout) :: self                 !< The field.
   integer(I4P),              intent(in)    :: ratio                !< Refinement ratio.
   integer(I8P), allocatable, intent(in)    :: block_to_refine(:,:) !< List of blocks to be refined.
   integer(I8P), allocatable, intent(in)    :: block_refined(:,:)   !< List of refined blocks with Morton code.
   real(R8P),                 intent(inout) :: q(1:,              &
                                                 1-grid%ngc:,&
                                                 1-grid%ngc:,&
                                                 1-grid%ngc:,&
                                                 1:)                !< Field cell centered variables.
   integer(I4P)                             :: b, i, j, k           !< Spatial counter.
   integer(I4P)                             :: ib, ic, ic_local     !< Counter.
   integer(I4P)                             :: i_fine, j_fine       !< Counter.
   integer(I4P)                             :: i_delta, j_delta     !< Counter.
   integer(I4P)                             :: ic1, ic2, ic3, ic4   !< Counter.

   associate(ni=>grid%ni, nj=>grid%nj, nk=>grid%nk, q_work=>self%q_work)
   if (allocated(block_to_refine)) then
      do b=1, size(block_to_refine, dim=2)
         if (mpih%myrank /= block_to_refine(2,b)) cycle
         ib = block_to_refine(1,b)

         q_work(:,:,:,:,ib) = q(:,:,:,:,ib)

         do ic_local=1, ratio
            ic = block_refined(2,(b-1)*ratio+ic_local)
            ic1 = mod( ic_local - 1   , 2)
            ic2 = mod((ic_local - 1)/2, 2)
            do k=1,nk
               do j=1+nj/2*ic2,nj/2+nj/2*ic2
                  do i=1+ni/2*ic1,ni/2+ni/2*ic1
                     j_fine = mod(j - 1, nj/2) * 2 + 1
                     i_fine = mod(i - 1, ni/2) * 2 + 1
                     q(:,i_fine:i_fine+1,j_fine:j_fine+1,k,ic) = 0._R8P
                     do j_delta=0,1
                     do i_delta=0,1
                     q(:,i_fine,  j_fine,  k,ic) = q(:,i_fine,j_fine,k,ic) +   &
                                                   (0.25_R8P + i_delta * 0.5_R8P) * &
                                                   (0.25_R8P + j_delta * 0.5_R8P) * &
                                                   q_work(:,i+i_delta-1, j+j_delta-1,k,ib)
                     q(:,i_fine+1,j_fine,  k,ic) = q(:,i_fine+1,j_fine,k,ic) + &
                                                   (0.75_R8P - i_delta * 0.5_R8P) * &
                                                   (0.25_R8P + j_delta * 0.5_R8P) * &
                                                   q_work(:,i+i_delta,   j+j_delta-1,k,ib)
                     q(:,i_fine,  j_fine+1,k,ic) = q(:,i_fine,j_fine+1,k,ic) + &
                                                   (0.25_R8P + i_delta * 0.5_R8P) * &
                                                   (0.75_R8P - j_delta * 0.5_R8P) * &
                                                   q_work(:,i+i_delta-1, j+j_delta  ,k,ib)
                     q(:,i_fine+1,j_fine+1,k,ic) = q(:,i_fine+1,j_fine+1,k,ic) + &
                                                   (0.75_R8P - i_delta * 0.5_R8P) * &
                                                   (0.75_R8P - j_delta * 0.5_R8P) * &
                                                   q_work(:,i+i_delta,   j+j_delta  ,k,ib)
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
         self%code(ic1) = block_refined(1,(b-1)*ratio+1)
         self%code(ic2) = block_refined(1,(b-1)*ratio+2)
         self%code(ic3) = block_refined(1,(b-1)*ratio+3)
         self%code(ic4) = block_refined(1,(b-1)*ratio+4)
      enddo
   endif
   endassociate
   endsubroutine refine2D

   subroutine refine3D(self, ratio, block_to_refine, block_refined, q)
   !< Refine blocks.
   !<
   !< Note: blocks number is not updated: mpi redistribute does it. This is dangerous...
   class(field_object),       intent(inout) :: self                      !< The field.
   integer(I4P),              intent(in)    :: ratio                     !< Refinement ratio.
   integer(I8P), allocatable, intent(in)    :: block_to_refine(:,:)      !< List of blocks to be refined.
   integer(I8P), allocatable, intent(in)    :: block_refined(:,:)        !< List of refined blocks with Morton code.
   real(R8P),                 intent(inout) :: q(1:,              &
                                                 1-grid%ngc:,&
                                                 1-grid%ngc:,&
                                                 1-grid%ngc:,&
                                                 1:)                     !< Field cell centered variables.
   integer(I4P)                             :: b, i, j, k                !< Spatial counter.
   integer(I4P)                             :: ib, ic, ic_local          !< Counter.
   integer(I4P)                             :: i_fine, j_fine, k_fine    !< Counter.
   integer(I4P)                             :: i_delta, j_delta, k_delta !< Counter.
   integer(I4P)                             :: ic1, ic2, ic3, ic4        !< Counter.
   integer(I4P)                             :: ic5, ic6, ic7, ic8        !< Counter.

   associate(ni=>grid%ni, nj=>grid%nj, nk=>grid%nk, q_work=>self%q_work)
   if (allocated(block_to_refine)) then
      do b=1, size(block_to_refine, dim=2)
         if (mpih%myrank /= block_to_refine(2,b)) cycle
         ib = block_to_refine(1,b)

         q_work(:,:,:,:,ib) = q(:,:,:,:,ib)

         do ic_local=1, ratio
            ic = block_refined(2,(b-1)*ratio+ic_local)
            ic1 = mod( ic_local - 1   , 2)
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
   endsubroutine refine3D

   subroutine update_coordinates(self)
   !< Update coordinates using the updated data in maps (that in turn is updated by tree).
   class(field_object), intent(inout) :: self !< The field.

   if (self%blocks_number>0) call self%maps%get_block_layout(coordinates=self%coordinates(:, 1:self%blocks_number))
   endsubroutine update_coordinates
endmodule adam_field_object
