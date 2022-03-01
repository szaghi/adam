!< ADAM, grid class definition.
module adam_grid_object
!< ADAM, grid class definition.

use adam_parameters
use FINER, only : file_ini
use PENF

implicit none
private
public :: grid_object

! grid parameters
integer(I4P), parameter :: MAX_REF_LEVELS = 100_I4P !< Maximum refinement levels.

type :: grid_object
   !< Grid class definition.
   real(R8P)                 :: domain_emin(3)=[0._R8P,0._R8P,0._R8P] !< Coordinates of minimum abscissa of whole domain.
   real(R8P)                 :: domain_emax(3)=[1._R8P,1._R8P,1._R8P] !< Coordinates of maximum abscissa of whole domain.
   integer(I4P)              :: ni=16_I4P                             !< Number of cells in i direction.
   integer(I4P)              :: nj=16_I4P                             !< Number of cells in j direction.
   integer(I4P)              :: nk=16_I4P                             !< Number of cells in k direction.
   integer(I4P)              :: ngc=2_I4P                             !< Number of ghost cells.
   integer(I4P)              :: weight_neighbor(26)=0_I4P             !< Weight of neighbors (cells number).
   integer(I4P)              :: bc_type(6)=0_I4P                      !< Type of boundary conditions in the 6 faces of grid.
   logical                   :: is_ijk_periodic(3)=.false.            !< Flag to indicate if the direction i, j or k is periodic.
   real(R8P),    allocatable :: block_dxyz(:,:)                       !< Blocks space steps for each level [3,MAX_REF_LEVELS].
   real(R8P),    allocatable :: cell_dxyz(:,:)                        !< Cells  space steps for each level [3,MAX_REF_LEVELS].
   integer(I4P), allocatable :: nb_max(:)                             !< Number of maximum blocks in each direction for each level.
   real(R8P),    allocatable :: lin_space_x(:,:)                      !< Lin. space x for each level [0-ngc:ni+ngc,MAX_REF_LEVELS].
   real(R8P),    allocatable :: lin_space_y(:,:)                      !< Lin. space y for each level [0-ngc:nj+ngc,MAX_REF_LEVELS].
   real(R8P),    allocatable :: lin_space_z(:,:)                      !< Lin. space z for each level [0-ngc:nk+ngc,MAX_REF_LEVELS].
   contains
      ! public methods
      procedure, pass(self) :: block_emin              !< Return block emin given its coordinates.
      procedure, pass(self) :: block_emax              !< Return block emax given its coordinates.
      procedure, pass(self) :: cell_xyz                !< Return cells xyz abscissa given block coordinates.
      procedure, pass(self) :: compute_metrics         !< Compute metrics of a block.
      procedure, pass(self) :: compute_weight_neighbor !< Compute weight of neighbors.
      procedure, pass(self) :: destroy                 !< Destroy the field.
      procedure, pass(self) :: do_cplane_intersect     !< Return true if a block is intersected by coordinate-plane.
      procedure, pass(self) :: get_closest_block       !< Get the closest block to a given point at a given level.
      procedure, pass(self) :: initialize              !< Initialize the field.
      procedure, pass(self) :: load_from_ini_file      !< Load object data from INI file.
      procedure, pass(self) :: node_xyz                !< Return nodes xyz abscissa given block coordinates.
      procedure, pass(self) :: print_status            !< Print status of main data.
      ! operators
      generic :: assignment(=) => grid_assign_grid      !< Overload `=`.
      procedure, pass(lhs), private :: grid_assign_grid !< Operator `=`.
endtype grid_object

contains
   ! public methods
   function block_emin(self, coordinates) result(emin)
   !< Return block emin given its coordinates.
   class(grid_object), intent(in) :: self           !< The grid.
   integer(I4P),       intent(in) :: coordinates(4) !< Block coordinates.
   real(R8P)                      :: emin(3)        !< Min abscissa of block.

   emin(:) = coordinates(:) * self%block_dxyz(:, coordinates(4))
   endfunction block_emin

   function block_emax(self, coordinates) result(emax)
   !< Return block emax given its coordinates.
   class(grid_object), intent(in) :: self           !< The grid.
   integer(I4P),       intent(in) :: coordinates(4) !< Block coordinates.
   real(R8P)                      :: emax(3)        !< Max abscissa of block.

   emax(:) = self%block_emin(coordinates) + self%block_dxyz(:, coordinates(4))
   endfunction block_emax

   subroutine node_xyz(self, coordinates, x_node, y_node, z_node)
   !< Return nodes xyz abscissa given block coordinates.
   class(grid_object), intent(in)            :: self                                !< The grid.
   integer(I4P),       intent(in)            :: coordinates(4)                      !< Block coordinates.
   real(R8P),          intent(out), optional :: x_node(0-self%ngc:self%ni+self%ngc) !< X coordinates.
   real(R8P),          intent(out), optional :: y_node(0-self%ngc:self%nj+self%ngc) !< Y coordinates.
   real(R8P),          intent(out), optional :: z_node(0-self%ngc:self%nk+self%ngc) !< Z coordinates.
   real(R8P)                                 :: emin(3)                             !< Min abscissa of block.
   integer(I4P)                              :: i, j, k                             !< Counter.

   emin = self%block_emin(coordinates)
   if (present(x_node)) x_node(:) = emin(1) + self%lin_space_x(:,coordinates(4))
   if (present(y_node)) y_node(:) = emin(2) + self%lin_space_y(:,coordinates(4))
   if (present(z_node)) z_node(:) = emin(3) + self%lin_space_z(:,coordinates(4))
   ! if (present(x_node)) then
   !    do i=0-self%ngc, self%ni+self%ngc
   !       x_node(i) = emin(1) + i * self%cell_dxyz(1,coordinates(4))
   !    enddo
   ! endif
   ! if (present(y_node)) then
   !    do j=0-self%ngc, self%nj+self%ngc
   !       y_node(j) = emin(2) + j * self%cell_dxyz(2,coordinates(4))
   !    enddo
   ! endif
   ! if (present(z_node)) then
   !    do k=0-self%ngc, self%nk+self%ngc
   !       z_node(k) = emin(3) + k * self%cell_dxyz(3,coordinates(4))
   !    enddo
   ! endif
   endsubroutine node_xyz

   subroutine cell_xyz(self, coordinates, x_cell, y_cell, z_cell)
   !< Return cells xyz abscissa given block coordinates.
   class(grid_object), intent(in)            :: self                                !< The grid.
   integer(I4P),       intent(in)            :: coordinates(4)                      !< Block coordinates.
   real(R8P),          intent(out), optional :: x_cell(1-self%ngc:self%ni+self%ngc) !< X coordinates.
   real(R8P),          intent(out), optional :: y_cell(1-self%ngc:self%nj+self%ngc) !< Y coordinates.
   real(R8P),          intent(out), optional :: z_cell(1-self%ngc:self%nk+self%ngc) !< Z coordinates.
   real(R8P)                                 :: emin(3)                             !< Min abscissa of block.
   integer(I4P)                              :: i, j, k                             !< Counter.

   emin = self%block_emin(coordinates)
   if (present(x_cell)) x_cell(:) = emin(1) + self%lin_space_x(1-self%ngc:self%ni+self%ngc,coordinates(4))
   if (present(y_cell)) y_cell(:) = emin(2) + self%lin_space_y(1-self%ngc:self%nj+self%ngc,coordinates(4))
   if (present(z_cell)) z_cell(:) = emin(3) + self%lin_space_z(1-self%ngc:self%nk+self%ngc,coordinates(4))
   ! if (present(x_cell)) then
   !    do i=1-self%ngc, self%ni+self%ngc
   !       x_cell(i) = emin(1) + (i-0.5_R8P) * self%cell_dxyz(1,coordinates(4))
   !    enddo
   ! endif
   ! if (present(y_cell)) then
   !    do j=1-self%ngc, self%nj+self%ngc
   !       y_cell(j) = emin(2) + (j-0.5_R8P) * self%cell_dxyz(2,coordinates(4))
   !    enddo
   ! endif
   ! if (present(z_cell)) then
   !    do k=1-self%ngc, self%nk+self%ngc
   !       z_cell(k) = emin(3) + (k-0.5_R8P) * self%cell_dxyz(3,coordinates(4))
   !    enddo
   ! endif
   endsubroutine cell_xyz

   subroutine compute_metrics(self, coordinates,      &
                              dx, dy, dz,             &
                              emin, emax,             &
                              x_node, y_node, z_node, &
                              x_cell, y_cell, z_cell)
   !< Compute metrics of a block.
   class(grid_object), intent(in)            :: self                                 !< The grid.
   integer(I4P),       intent(in)            :: coordinates(4)                       !< Block coordinates.
   real(R8P),          intent(out), optional :: dx, dy, dz                           !< Space steps.
   real(R8P),          intent(out), optional :: emin(3), emax(3)                     !< Min/max abscissa of block.
   real(R8P),          intent(out), optional :: x_node(0-self%ngc:self%ni+self%ngc)  !< X coordinates.
   real(R8P),          intent(out), optional :: y_node(0-self%ngc:self%nj+self%ngc)  !< Y coordinates.
   real(R8P),          intent(out), optional :: z_node(0-self%ngc:self%nk+self%ngc)  !< Z coordinates.
   real(R8P),          intent(out), optional :: x_cell(1-self%ngc:self%ni+self%ngc)  !< X coordinates.
   real(R8P),          intent(out), optional :: y_cell(1-self%ngc:self%nj+self%ngc)  !< Y coordinates.
   real(R8P),          intent(out), optional :: z_cell(1-self%ngc:self%nk+self%ngc)  !< Z coordinates.

   if (present(dx)) dx = self%cell_dxyz(1,coordinates(4))
   if (present(dy)) dy = self%cell_dxyz(2,coordinates(4))
   if (present(dz)) dz = self%cell_dxyz(3,coordinates(4))
   if (present(emin)) emin = self%block_emin(coordinates)
   if (present(emax)) emax = self%block_emax(coordinates)
   call self%node_xyz(coordinates=coordinates, x_node=x_node, y_node=y_node, z_node=z_node)
   call self%cell_xyz(coordinates=coordinates, x_cell=x_cell, y_cell=y_cell, z_cell=z_cell)
   endsubroutine compute_metrics

   subroutine compute_weight_neighbor(self)
   !< Compute weight of neighbors.
   class(grid_object), intent(inout) :: self    !< The grid.
   integer(I4P)                      :: nijk(3) !< Ni, nj, nk in array.
   integer(I4P)                      :: ngc     !< Ghost cells number.
   integer(I4P)                      :: fec     !< Counter.
   integer(I4P)                      :: i       !< Counter.

   self%weight_neighbor = 1_I4P
   nijk = [self%ni, self%nj, self%nk]
   ngc = self%ngc
   do fec=1, 26
      do i=1, 3
         self%weight_neighbor(fec) = self%weight_neighbor(fec) * (  abs(fec_to_delta(i,fec))  * ngc + &
                                                                 (1-abs(fec_to_delta(i,fec))) * nijk(i))
      enddo
   enddo
   endsubroutine compute_weight_neighbor

   elemental subroutine destroy(self)
   !< Destroy field.
   class(grid_object), intent(inout) :: self  !< The grid.
   type(grid_object)                 :: fresh !< Fresh grid.

   self = fresh
   endsubroutine destroy

   function do_cplane_intersect(self, emin, emax, dxyz, cplane_origin, cplane_normal, cplane_block_indexes) result(do_intersect)
   !< Return true if a block is intersected by coordinate-plane.
   class(grid_object), intent(inout)         :: self                    !< The grid.
   real(R8P),          intent(in)            :: emin(3), emax(3)        !< Block extents.
   real(R8P),          intent(in)            :: dxyz(3)                 !< Block space steps.
   real(R8P),          intent(in)            :: cplane_origin(3)        !< Coordinate-plane origin.
   real(R8P),          intent(in)            :: cplane_normal(3)        !< Coordinate-plane normal.
   integer(I4P),       intent(out), optional :: cplane_block_indexes(3) !< Block-local indexes of cplane intersection.
   logical                                   :: do_intersect            !< Test result.

   do_intersect = .false.
   if     (nint(cplane_normal(1))==1) then
      if ((cplane_origin(1) >= emin(1)).and.(cplane_origin(1) <= emax(1))) then
         do_intersect = .true.
         if (present(cplane_block_indexes)) then
            cplane_block_indexes(1) = ceiling((cplane_origin(1) - emin(1)) / dxyz(1), I4P)
         endif
      endif
   elseif (nint(cplane_normal(2))==1) then
      if ((cplane_origin(2) >= emin(2)).and.(cplane_origin(2) <= emax(2))) then
         do_intersect = .true.
         if (present(cplane_block_indexes)) then
            cplane_block_indexes(2) = ceiling((cplane_origin(2) - emin(2)) / dxyz(2), I4P)
         endif
      endif
   elseif (nint(cplane_normal(3))==1) then
      if ((cplane_origin(3) >= emin(3)).and.(cplane_origin(3) <= emax(3))) then
         do_intersect = .true.
         if (present(cplane_block_indexes)) then
            cplane_block_indexes(3) = ceiling((cplane_origin(3) - emin(3)) / dxyz(3), I4P)
         endif
      endif
   endif
   endfunction do_cplane_intersect

   function get_closest_block(self, point, level) result(ijk)
   !< Get the closest block to a given point at a given level.
   class(grid_object), intent(in) :: self     !< The grid.
   real(R8P),          intent(in) :: point(3) !< Point xyz coordinates.
   integer(I4P),       intent(in) :: level    !< Refinement level.
   integer(I4P)                   :: ijk(3)   !< Indexes of the closest (living or not) block.

   associate(nb_max=>self%nb_max(level), emin=>self%domain_emin, dxyz=>self%block_dxyz(:,level))
      ! ijk(:) = min(nb_max, max(1, ceiling((point(:) - emin(:)) / dxyz(:), I4P)))
      ijk(:) = int((point(:) - emin(:)) / dxyz(:), I4P)
      if (any(ijk)<0.or.any(ijk)>2**level-1) then
         print '(A)', 'ERROR: grid%get_closest block failed ijk: '//str(ijk)//&
               ' level:'//str(level)//' point:'//str(point)
      endif
   endassociate
   endfunction get_closest_block

   subroutine initialize(self, file_parameters, ni, nj, nk, ngc, emin, emax, bc_type)
   !< Initialize field.
   class(grid_object), intent(inout)           :: self            !< The grid.
   type(file_ini),     intent(inout), optional :: file_parameters !< INI file handler.
   integer(I4P),       intent(in),    optional :: ni              !< Number of cells in X direction.
   integer(I4P),       intent(in),    optional :: nj              !< Number of cells in Y direction.
   integer(I4P),       intent(in),    optional :: nk              !< Number of cells in Z direction.
   integer(I4P),       intent(in),    optional :: ngc             !< Number of ghost cells.
   real(R8P),          intent(in),    optional :: emin(3)         !< Coordinates of minium abscissa.
   real(R8P),          intent(in),    optional :: emax(3)         !< Coordinates of maxium abscissa.
   integer(I4P),       intent(in),    optional :: bc_type(6)      !< Type of boundary conditions in the 6 faces of grid.
   integer(I4P)                                :: i, j, k, l      !< Counter.
   integer(I4P)                                :: nijk(3)         !< Cells number.

   call self%destroy
   if (present(file_parameters)) call self%load_from_ini_file(file_parameters)

   ! parameters explicitely passed ovveride ones file-passed
   if (present(emin)) self%domain_emin = emin
   if (present(emax)) self%domain_emax = emax
   if (present(ni))  self%ni = ni
   if (present(nj))  self%nj = nj
   if (present(nk))  self%nk = nk
   if (present(ngc)) self%ngc = ngc
   if (present(bc_type)) self%bc_type = bc_type

   call self%compute_weight_neighbor
   if (any(self%bc_type(1:2)==BC_PERIODIC)) self%is_ijk_periodic(1) = .true.
   if (any(self%bc_type(3:4)==BC_PERIODIC)) self%is_ijk_periodic(2) = .true.
   if (any(self%bc_type(5:6)==BC_PERIODIC)) self%is_ijk_periodic(3) = .true.

   nijk = [self%ni, self%nj, self%nk]
   allocate(self%block_dxyz(3,                           MAX_REF_LEVELS))
   allocate(self%cell_dxyz( 3,                           MAX_REF_LEVELS))
   allocate(self%nb_max(                                 MAX_REF_LEVELS))
   allocate(self%lin_space_x(0-self%ngc:self%ni+self%ngc,MAX_REF_LEVELS))
   allocate(self%lin_space_y(0-self%ngc:self%nj+self%ngc,MAX_REF_LEVELS))
   allocate(self%lin_space_z(0-self%ngc:self%nk+self%ngc,MAX_REF_LEVELS))
   do l=1, MAX_REF_LEVELS
      self%nb_max(l) = 2**l
      self%block_dxyz(:,l) = (self%domain_emax(:) - self%domain_emin(:)) / self%nb_max(l)
      self%cell_dxyz(:,l) = self%block_dxyz(:,l) / nijk(:)
      do i=0-self%ngc, self%ni+self%ngc
         self%lin_space_x(i,l) = i * self%cell_dxyz(1,l)
      enddo
      do j=0-self%ngc, self%nj+self%ngc
         self%lin_space_y(j,l) = j * self%cell_dxyz(2,l)
      enddo
      do k=0-self%ngc, self%nk+self%ngc
         self%lin_space_z(k,l) = k * self%cell_dxyz(3,l)
      enddo
   enddo
   endsubroutine initialize

   subroutine load_from_ini_file(self, file_parameters)
   !< Load object data from INI file.
   class(grid_object), intent(inout) :: self            !< The grid.
   type(file_ini),     intent(inout) :: file_parameters !< INI file handler.
   integer(I4P)                      :: buff_I4P        !< I4P buffer.
   real(R8P)                         :: buff_R8P        !< R8P buffer.
   logical                           :: buff_LOG        !< LOG buffer.

   call file_parameters%get(section_name='grid', option_name='ni'           , val=buff_I4P) ; self%ni                 = buff_I4P
   call file_parameters%get(section_name='grid', option_name='nj'           , val=buff_I4P) ; self%nj                 = buff_I4P
   call file_parameters%get(section_name='grid', option_name='nk'           , val=buff_I4P) ; self%nk                 = buff_I4P
   call file_parameters%get(section_name='grid', option_name='ngc'          , val=buff_I4P) ; self%ngc                = buff_I4P
   call file_parameters%get(section_name='grid', option_name='emin_x'       , val=buff_R8P) ; self%domain_emin(1)     = buff_R8P
   call file_parameters%get(section_name='grid', option_name='emin_y'       , val=buff_R8P) ; self%domain_emin(2)     = buff_R8P
   call file_parameters%get(section_name='grid', option_name='emin_z'       , val=buff_R8P) ; self%domain_emin(3)     = buff_R8P
   call file_parameters%get(section_name='grid', option_name='emax_x'       , val=buff_R8P) ; self%domain_emax(1)     = buff_R8P
   call file_parameters%get(section_name='grid', option_name='emax_y'       , val=buff_R8P) ; self%domain_emax(2)     = buff_R8P
   call file_parameters%get(section_name='grid', option_name='emax_z'       , val=buff_R8P) ; self%domain_emax(3)     = buff_R8P
   call file_parameters%get(section_name='grid', option_name='bc_type_1'    , val=buff_I4P) ; self%bc_type(1)         = buff_I4P
   call file_parameters%get(section_name='grid', option_name='bc_type_2'    , val=buff_I4P) ; self%bc_type(2)         = buff_I4P
   call file_parameters%get(section_name='grid', option_name='bc_type_3'    , val=buff_I4P) ; self%bc_type(3)         = buff_I4P
   call file_parameters%get(section_name='grid', option_name='bc_type_4'    , val=buff_I4P) ; self%bc_type(4)         = buff_I4P
   call file_parameters%get(section_name='grid', option_name='bc_type_5'    , val=buff_I4P) ; self%bc_type(5)         = buff_I4P
   call file_parameters%get(section_name='grid', option_name='bc_type_6'    , val=buff_I4P) ; self%bc_type(6)         = buff_I4P
   call file_parameters%get(section_name='grid', option_name='is_i_periodic', val=buff_LOG) ; self%is_ijk_periodic(1) = buff_LOG
   call file_parameters%get(section_name='grid', option_name='is_j_periodic', val=buff_LOG) ; self%is_ijk_periodic(2) = buff_LOG
   call file_parameters%get(section_name='grid', option_name='is_k_periodic', val=buff_LOG) ; self%is_ijk_periodic(3) = buff_LOG
   endsubroutine load_from_ini_file

   subroutine print_status(self)
   !< Print status of main data.
   class(grid_object), intent(in) :: self !< The field.

   print '(A)',          'grid status of main data'
   print '(A)',          '  domain minimum extent: '//trim(str(self%domain_emin    ))
   print '(A)',          '  domain maximum extent: '//trim(str(self%domain_emax    ))
   print '(A)',          '  ni:                    '//trim(str(self%ni             ))
   print '(A)',          '  nj:                    '//trim(str(self%nj             ))
   print '(A)',          '  nk:                    '//trim(str(self%nk             ))
   print '(A)',          '  ngc:                   '//trim(str(self%ngc            ))
   print '(A)',          '  boundary conditions:   '//trim(str(self%bc_type        ))
   print '(A,3(L1,1X))', '  IJK periodic:          ',          self%is_ijk_periodic
   print '(A)',          ''
   endsubroutine print_status

   ! operators
   ! =
   pure subroutine grid_assign_grid(lhs, rhs)
   !< Operator `=`.
   class(grid_object), intent(inout) :: lhs !< Left hand side.
   type(grid_object),  intent(in)    :: rhs !< Right hand side.

   lhs%domain_emin     = rhs%domain_emin
   lhs%domain_emax     = rhs%domain_emax
   lhs%ni              = rhs%ni
   lhs%nj              = rhs%nj
   lhs%nk              = rhs%nk
   lhs%ngc             = rhs%ngc
   lhs%weight_neighbor = rhs%weight_neighbor
   lhs%bc_type         = rhs%bc_type
   lhs%is_ijk_periodic = rhs%is_ijk_periodic
   endsubroutine grid_assign_grid
endmodule adam_grid_object
