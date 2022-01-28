!< ADAM, grid class definition.
module adam_grid_object
!< ADAM, grid class definition.

use adam_parameters
use FINER, only : file_ini
use PENF

implicit none
private
public :: grid_object

type :: grid_object
   !< Grid class definition.
   real(R8P)    :: domain_emin(3)=[0._R8P,0._R8P,0._R8P] !< Coordinates of minimum abscissa of whole domain.
   real(R8P)    :: domain_emax(3)=[1._R8P,1._R8P,1._R8P] !< Coordinates of maximum abscissa of whole domain.
   integer(I4P) :: ni=16_I4P                             !< Number of cells in i direction.
   integer(I4P) :: nj=16_I4P                             !< Number of cells in j direction.
   integer(I4P) :: nk=16_I4P                             !< Number of cells in k direction.
   integer(I4P) :: ngc=2_I4P                             !< Number of ghost cells.
   integer(I4P) :: weight_neighbor(26)=0_I4P             !< Weight of neighbors (cells number).
   integer(I4P) :: bc_type(6)=0_I4P                      !< Type of boundary conditions in the 6 faces of grid.
   logical      :: is_ijk_periodic(3)=.false.            !< Flag to indicate if the direction i, j or k is periodic.
   contains
      ! public methods
      procedure, pass(self) :: compute_metrics         !< Compute metrics of a block.
      procedure, pass(self) :: compute_weight_neighbor !< Compute weight of neighbors.
      procedure, pass(self) :: destroy                 !< Destroy the field.
      procedure, pass(self) :: do_cplane_intersect     !< Return true if a block is intersected by coordinate-plane.
      procedure, pass(self) :: initialize              !< Initialize the field.
      procedure, pass(self) :: load_from_ini_file      !< Load object data from INI file.
      procedure, pass(self) :: print_status            !< Print status of main data.
      ! operators
      generic :: assignment(=) => grid_assign_grid      !< Overload `=`.
      procedure, pass(lhs), private :: grid_assign_grid !< Operator `=`.
endtype grid_object

contains
   ! public methods
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
   real(R8P)                                 :: emin_(3), emax_(3)                   !< Min/max abscissa of block, local var.
   real(R8P)                                 :: x_node_(0-self%ngc:self%ni+self%ngc) !< X coordinates, local var.
   real(R8P)                                 :: y_node_(0-self%ngc:self%nj+self%ngc) !< Y coordinates, local var.
   real(R8P)                                 :: z_node_(0-self%ngc:self%nk+self%ngc) !< Z coordinates, local var.
   real(R8P)                                 :: x_cell_(1-self%ngc:self%ni+self%ngc) !< X coordinates, local var.
   real(R8P)                                 :: y_cell_(1-self%ngc:self%nj+self%ngc) !< Y coordinates, local var.
   real(R8P)                                 :: z_cell_(1-self%ngc:self%nk+self%ngc) !< Z coordinates, local var.
   real(R8P)                                 :: dx_, dy_, dz_                        !< Space steps, local var.
   integer(I4P)                              :: i, j, k, l                           !< Counter.

   i = coordinates(1)
   j = coordinates(2)
   k = coordinates(3)
   l = coordinates(4)
   dx_ = (self%domain_emax(1) - self%domain_emin(1)) / 2**l
   dy_ = (self%domain_emax(2) - self%domain_emin(2)) / 2**l
   dz_ = (self%domain_emax(3) - self%domain_emin(3)) / 2**l
   emin_(1) = i * dx_ ; emax_(1) = emin_(1) + dx_
   emin_(2) = j * dy_ ; emax_(2) = emin_(2) + dy_
   emin_(3) = k * dz_ ; emax_(3) = emin_(3) + dz_
   dx_ = dx_ / self%ni
   dy_ = dy_ / self%nj
   dz_ = dz_ / self%nk
   do i=0-self%ngc, self%ni+self%ngc
      x_node_(i) = emin_(1) + i * dx_
   enddo
   do j=0-self%ngc, self%nj+self%ngc
      y_node_(j) = emin_(2) + j * dy_
   enddo
   do k=0-self%ngc, self%nk+self%ngc
      z_node_(k) = emin_(3) + k * dz_
   enddo
   do i=1-self%ngc, self%ni+self%ngc
      x_cell_(i) = x_node_(i-1) + dx_ * 0.5_R8P
   enddo
   do j=1-self%ngc, self%nj+self%ngc
      y_cell_(j) = y_node_(j-1) + dy_ * 0.5_R8P
   enddo
   do k=1-self%ngc, self%nk+self%ngc
      z_cell_(k) = z_node_(k-1) + dz_ * 0.5_R8P
   enddo
   if (present(dx)) dx = dx_
   if (present(dy)) dy = dy_
   if (present(dz)) dz = dz_
   if (present(emin)) emin = emin_
   if (present(emax)) emax = emax_
   if (present(x_node)) x_node = x_node_
   if (present(y_node)) y_node = y_node_
   if (present(z_node)) z_node = z_node_
   if (present(x_cell)) x_cell = x_cell_
   if (present(y_cell)) y_cell = y_cell_
   if (present(z_cell)) z_cell = z_cell_
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
