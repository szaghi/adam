!< ADAM, grid class definition.
module adam_grid_object
!< ADAM, grid class definition.

! ADAM classes, libraries, parameters
use :: adam_parameters
! ADAM singleton objects
use :: adam_global_mpih, only : mpih
! third party modules
use :: finer, only : file_ini
use :: penf
! sdk modules
use :: mpi

implicit none
private
public :: grid_object

! grid parameters
integer(I4P), parameter :: MAX_REF_LEVELS = 30_I4P !< Maximum refinement levels.

character(len=4), parameter :: INI_SECTION_NAME="grid" !< INI (config) file section name containing configs.

type :: grid_object
   !< Grid class definition.
   real(R8P)                 :: domain_emin(3)=[0._R8P,0._R8P,0._R8P] !< Coordinates of minimum abscissa of whole domain.
   real(R8P)                 :: domain_emax(3)=[1._R8P,1._R8P,1._R8P] !< Coordinates of maximum abscissa of whole domain.
   integer(I4P)              :: ni=16_I4P                             !< Number of cells in i direction.
   integer(I4P)              :: nj=16_I4P                             !< Number of cells in j direction.
   integer(I4P)              :: nk=16_I4P                             !< Number of cells in k direction.
   integer(I4P)              :: ngc=2_I4P                             !< Number of ghost cells.
   integer(I4P)              :: block_weight=0._R8P                   !< Weight of single block.
   integer(I4P)              :: weight_neighbor(26)=0_I4P             !< Weight of neighbors (cells number).
   integer(I4P)              :: bc_type(6)=0_I4P                      !< Type of boundary conditions in the 6 faces of grid.
   logical                   :: is_ijk_periodic(3)=.false.            !< Flag to indicate if the direction i, j or k is periodic.
   real(R8P),    allocatable :: block_dxyz(:,:)                       !< Blocks space steps for each level [3,0:MAX_REF_LEVELS].
   real(R8P),    allocatable :: cell_dxyz(:,:)                        !< Cells  space steps for each level [3,0:MAX_REF_LEVELS].
   integer(I4P), allocatable :: nb_max(:)                             !< Number of maximum blocks in each direction for each level.
   real(R8P),    allocatable :: lin_space_x(:,:)                      !< Lin. space x for each level [0-ngc:ni+ngc,MAX_REF_LEVELS].
   real(R8P),    allocatable :: lin_space_y(:,:)                      !< Lin. space y for each level [0-ngc:nj+ngc,MAX_REF_LEVELS].
   real(R8P),    allocatable :: lin_space_z(:,:)                      !< Lin. space z for each level [0-ngc:nk+ngc,MAX_REF_LEVELS].
   logical                   :: null_xyz(3)=[.false.,.false.,.false.] !< Handy tags for nullify equations in some direction.
   contains
      ! public methods
      procedure, pass(self) :: block_emin              !< Return block emin given its coordinates.
      procedure, pass(self) :: block_emax              !< Return block emax given its coordinates.
      procedure, pass(self) :: cell_xyz                !< Return cells xyz abscissa given block coordinates.
      procedure, pass(self) :: compute_metrics         !< Compute metrics of a block.
      procedure, pass(self) :: compute_weight_neighbor !< Compute weight of neighbors.
      procedure, pass(self) :: description             !< Return pretty-printed object description.
      procedure, nopass     :: do_cplane_intersect     !< Return true if a block is intersected by coordinate-plane.
      procedure, pass(self) :: fec_bc_type             !< Return BC type of given fec.
      procedure, pass(self) :: get_closest_block       !< Get the closest block to a given point at a given level.
      procedure, pass(self) :: initialize              !< Initialize the field.
      procedure, pass(self) :: load_from_ini_file      !< Load object data from INI file.
      procedure, pass(self) :: node_xyz                !< Return nodes xyz abscissa given block coordinates.
      procedure, pass(self) :: set_bc_type             !< Set grid boundary conditions accordingly with app'equation object.
endtype grid_object

contains
   ! public methods
   function block_emin(self, coordinates) result(emin)
   !< Return block emin given its coordinates.
   class(grid_object), intent(in) :: self           !< The grid.
   integer(I4P),       intent(in) :: coordinates(4) !< Block coordinates.
   real(R8P)                      :: emin(3)        !< Min abscissa of block.

   emin(1:3) = self%domain_emin(1:3) + coordinates(1:3) * self%block_dxyz(1:3, coordinates(4))
   endfunction block_emin

   function block_emax(self, coordinates) result(emax)
   !< Return block emax given its coordinates.
   class(grid_object), intent(in) :: self           !< The grid.
   integer(I4P),       intent(in) :: coordinates(4) !< Block coordinates.
   real(R8P)                      :: emax(3)        !< Max abscissa of block.

   emax(:) = self%block_emin(coordinates) + self%block_dxyz(:, coordinates(4))
   endfunction block_emax

   subroutine cell_xyz(self, coordinates, x_cell, y_cell, z_cell)
   !< Return cells xyz abscissa given block coordinates.
   class(grid_object), intent(in)            :: self                                !< The grid.
   integer(I4P),       intent(in)            :: coordinates(4)                      !< Block coordinates.
   real(R8P),          intent(out), optional :: x_cell(1-self%ngc:self%ni+self%ngc) !< X coordinates.
   real(R8P),          intent(out), optional :: y_cell(1-self%ngc:self%nj+self%ngc) !< Y coordinates.
   real(R8P),          intent(out), optional :: z_cell(1-self%ngc:self%nk+self%ngc) !< Z coordinates.
   integer(I4P)                              :: nijk(3)                             !< Cells number.
   real(R8P)                                 :: emin(3)                             !< Min abscissa of block.

   ! Cells number in each direction.
   nijk = [self%ni, self%nj, self%nk]

   emin = self%block_emin(coordinates)
   associate(l => coordinates(4))
   if (present(x_cell)) x_cell(:) = emin(1) - self%block_dxyz(1,l)/nijk(1)*0.5_R8P + self%lin_space_x(1-self%ngc:self%ni+self%ngc,l)
   if (present(y_cell)) y_cell(:) = emin(2) - self%block_dxyz(2,l)/nijk(2)*0.5_R8P + self%lin_space_y(1-self%ngc:self%nj+self%ngc,l)
   if (present(z_cell)) z_cell(:) = emin(3) - self%block_dxyz(3,l)/nijk(3)*0.5_R8P + self%lin_space_z(1-self%ngc:self%nk+self%ngc,l)
   endassociate
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
         self%weight_neighbor(fec) = self%weight_neighbor(fec) * (  abs(FEC_TO_DELTA(i,fec))  * ngc + &
                                                                 (1-abs(FEC_TO_DELTA(i,fec))) * nijk(i))
      enddo
   enddo
   endsubroutine compute_weight_neighbor

   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(grid_object), intent(in) :: self             !< The grid.
   character(len=:), allocatable  :: desc             !< Description.
   character(len=1), parameter    :: NL=new_line('a') !< New line character.
   integer(I4P)                   :: l                !< Counter.

   desc =       mpih%myrankstr//'grid main data'                                            //NL
   desc = desc//mpih%myrankstr//'  domain minimum extent: '//trim(str(self%domain_emin    ))//NL
   desc = desc//mpih%myrankstr//'  domain maximum extent: '//trim(str(self%domain_emax    ))//NL
   desc = desc//mpih%myrankstr//'  ni:                    '//trim(str(self%ni             ))//NL
   desc = desc//mpih%myrankstr//'  nj:                    '//trim(str(self%nj             ))//NL
   desc = desc//mpih%myrankstr//'  nk:                    '//trim(str(self%nk             ))//NL
   desc = desc//mpih%myrankstr//'  ngc:                   '//trim(str(self%ngc            ))//NL
   desc = desc//mpih%myrankstr//'  boundary conditions:   '//trim(str(self%bc_type        ))//NL
   desc = desc//mpih%myrankstr//'  block weight:          '//trim(str(self%block_weight   ))//NL
   desc = desc//mpih%myrankstr//'  IJK periodic:          '//str(self%is_ijk_periodic(1))//' '//&
                                                                  str(self%is_ijk_periodic(2))//' '//&
                                                                  str(self%is_ijk_periodic(3))//NL
   desc = desc//mpih%myrankstr//'  nullify xyz:           '//str(self%null_xyz(1))//' '//&
                                                                  str(self%null_xyz(2))//' '//&
                                                                  str(self%null_xyz(3))
   do l=0, ubound(self%cell_dxyz, dim=2)
   desc = desc//NL//mpih%myrankstr//'  dxyz l='//trim(strz(l,2)) //'             '//trim(str(self%cell_dxyz(:,l)))
   enddo
   endfunction description

   function do_cplane_intersect(emin, emax, dxyz, cplane_origin, cplane_normal, cplane_block_indexes) result(do_intersect)
   !< Return true if a block is intersected by coordinate-plane.
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

   function fec_bc_type(self, fec) result(bc_type)
   !< Return BC type of given fec.
   class(grid_object), intent(in) :: self    !< The grid.
   integer(I4P),       intent(in) :: fec     !< Current fec.
   integer(I4P)                   :: bc_type !< BC type.

   select case(fec)
   case(1,7,9,11,13,19,21,23,25)
      ! BC i min
      bc_type = self%bc_type(1)
   case(2,8,10,12,14,20,22,24,26)
      ! BC i max
      bc_type = self%bc_type(2)
   case(3,15,17)
      ! BC j min
      bc_type = self%bc_type(3)
   case(4,16,18)
      ! BC j max
      bc_type = self%bc_type(4)
   case(5)
      ! BC k min
      bc_type = self%bc_type(5)
   case(6)
      ! BC k max
      bc_type = self%bc_type(6)
   endselect
   endfunction fec_bc_type

   function get_closest_block(self, point, level) result(ijk)
   !< Get the closest block to a given point at a given level.
   class(grid_object), intent(inout) :: self     !< The grid.
   real(R8P),          intent(in)    :: point(3) !< Point xyz coordinates.
   integer(I4P),       intent(in)    :: level    !< Refinement level.
   integer(I4P)                      :: ijk(3)   !< Indexes of the closest (living or not) block.

   associate(nb_max=>self%nb_max(level), emin=>self%domain_emin, dxyz=>self%block_dxyz(:,level))
      ijk(:) = int((point(:) - emin(:)) / dxyz(:), I4P)
      if (any(ijk<0).or.any(ijk>2**level-1)) then
         call mpih%abort(error_code=-110, msg='ERROR: grid_object%get_closest block failed ijk: '//str(ijk)//&
                                                   ' level:'//str(level)//' point:'//str(point))
      endif
   endassociate
   endfunction get_closest_block

   subroutine initialize(self, file_parameters, verbose)
   !< Initialize field.
   !< Note: bc_type is not initialized, must be set separately by equation app.
   class(grid_object), intent(inout)        :: self            !< The grid.
   type(file_ini),     intent(inout)        :: file_parameters !< INI file handler.
   logical,            intent(in), optional :: verbose         !< Flag to activate verbose output.
   logical                                  :: verbose_        !< Trigger verbose output, local variable.
   integer(I4P)                             :: i, j, k, l      !< Counter.
   integer(I4P)                             :: nijk(3)         !< Cells number.
   integer(I4P)                             :: ratio           !< AMR ratio.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (verbose_) call mpih%print_message('grid_object%initialize start')
   call self%load_from_ini_file(file_parameters)
   call file_parameters%get(section_name='amr', option_name='ratio', val=ratio)

   self%block_weight = (self%ngc+self%ni+self%ngc) * (self%ngc+self%nj+self%ngc) * (self%ngc+self%nk+self%ngc)

   call self%compute_weight_neighbor

   nijk = [self%ni, self%nj, self%nk]
   allocate(self%block_dxyz(3,                           0:MAX_REF_LEVELS))
   allocate(self%cell_dxyz( 3,                           0:MAX_REF_LEVELS))
   allocate(self%nb_max(                                 0:MAX_REF_LEVELS))
   allocate(self%lin_space_x(0-self%ngc:self%ni+self%ngc,0:MAX_REF_LEVELS))
   allocate(self%lin_space_y(0-self%ngc:self%nj+self%ngc,0:MAX_REF_LEVELS))
   allocate(self%lin_space_z(0-self%ngc:self%nk+self%ngc,0:MAX_REF_LEVELS))
   do l=0, MAX_REF_LEVELS
      self%nb_max(l) = 2**l
      select case(ratio)
      case(2_I4P)
         self%block_dxyz(1,l) = (self%domain_emax(1) - self%domain_emin(1)) / self%nb_max(l)
         self%block_dxyz(2,l) = (self%domain_emax(2) - self%domain_emin(2))
         self%block_dxyz(3,l) = (self%domain_emax(3) - self%domain_emin(3))
      case(4_I4P)
         self%block_dxyz(1,l) = (self%domain_emax(1) - self%domain_emin(1)) / self%nb_max(l)
         self%block_dxyz(2,l) = (self%domain_emax(2) - self%domain_emin(2)) / self%nb_max(l)
         self%block_dxyz(3,l) = (self%domain_emax(3) - self%domain_emin(3))
      case(8_I4P)
         self%block_dxyz(:,l) = (self%domain_emax(:) - self%domain_emin(:)) / self%nb_max(l)
      endselect
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
   if (verbose_) print '(A)', self%description()
   if (verbose_) call mpih%print_message('grid_object%initialize finish')
   endsubroutine initialize

   subroutine load_from_ini_file(self, file_parameters, go_on_fail)
   !< Load object data from INI file.
   class(grid_object), intent(inout)        :: self            !< The grid.
   type(file_ini),     intent(inout)        :: file_parameters !< INI file handler.
   logical,            intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                  :: go_on_fail_     !< Go on if load fails.
   integer(I4P)                             :: buff_I4P        !< I4P buffer.
   real(R8P)                                :: buff_R8P        !< R8P buffer.
   logical                                  :: buff_log        !< Logical buffer.
   integer(I4P)                             :: error           !< Error status.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='ni', val=buff_I4P, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(ni)')
   self%ni=buff_I4P
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='nj', val=buff_I4P, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(nj)')
   self%nj=buff_I4P
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='nk', val=buff_I4P, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(nk)')
   self%nk=buff_I4P
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='ngc', val=buff_I4P, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(ngc)')
   self%ngc=buff_I4P
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='emin_x', val=buff_R8P, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(emin_x)')
   self%domain_emin(1)=buff_R8P
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='emin_y', val=buff_R8P, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(emin_y)')
   self%domain_emin(2)=buff_R8P
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='emin_z', val=buff_R8P, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(emin_z)')
   self%domain_emin(3)=buff_R8P
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='emax_x', val=buff_R8P, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(emax_x)')
   self%domain_emax(1)=buff_R8P
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='emax_y', val=buff_R8P, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(emax_y)')
   self%domain_emax(2)=buff_R8P
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='emax_z', val=buff_R8P, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(emax_z)')
   self%domain_emax(3)=buff_R8P
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='null_x', val=buff_log, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(null_x)')
   self%null_xyz(1)=buff_log
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='null_y', val=buff_log, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(null_y)')
   self%null_xyz(2)=buff_log
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='null_z', val=buff_log, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(null_z)')
   self%null_xyz(3)=buff_log
   endsubroutine load_from_ini_file

   subroutine node_xyz(self, coordinates, x_node, y_node, z_node)
   !< Return nodes xyz abscissa given block coordinates.
   class(grid_object), intent(in)            :: self                                !< The grid.
   integer(I4P),       intent(in)            :: coordinates(4)                      !< Block coordinates.
   real(R8P),          intent(out), optional :: x_node(0-self%ngc:self%ni+self%ngc) !< X coordinates.
   real(R8P),          intent(out), optional :: y_node(0-self%ngc:self%nj+self%ngc) !< Y coordinates.
   real(R8P),          intent(out), optional :: z_node(0-self%ngc:self%nk+self%ngc) !< Z coordinates.
   real(R8P)                                 :: emin(3)                             !< Min abscissa of block.

   emin = self%block_emin(coordinates)
   if (present(x_node)) x_node(:) = emin(1) + self%lin_space_x(:,coordinates(4))
   if (present(y_node)) y_node(:) = emin(2) + self%lin_space_y(:,coordinates(4))
   if (present(z_node)) z_node(:) = emin(3) + self%lin_space_z(:,coordinates(4))
   endsubroutine node_xyz

   subroutine set_bc_type(self, bc_type)
   !< Set grid boundary conditions accordingly with app'equation object.
   class(grid_object), intent(inout) :: self       !< The grid.
   integer(I4P),       intent(in)    :: bc_type(6) !< Type of boundary conditions in the 6 faces of grid.

   self%bc_type = bc_type
   if (any(self%bc_type(1:2)==BC_PERIODIC)) self%is_ijk_periodic(1) = .true.
   if (any(self%bc_type(3:4)==BC_PERIODIC)) self%is_ijk_periodic(2) = .true.
   if (any(self%bc_type(5:6)==BC_PERIODIC)) self%is_ijk_periodic(3) = .true.
   endsubroutine set_bc_type
endmodule adam_grid_object
