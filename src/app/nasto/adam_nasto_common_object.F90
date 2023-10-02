!< ADAM, Navier-Stokes equations system class definition, common data to all backends.
module adam_nasto_common_object
!< ADAM, Navier-Stokes equations system class definition, common data to all backends.

use adam_adam_object
use adam_amr_marker_object
use adam_field_object
use adam_grid_object
use adam_mpih_object
use adam_slice_object
use adam_nasto_ic_object
use adam_nasto_physics_object
use adam_nasto_parameters
use FiNeR
use PENF
use ISO_C_BINDING

implicit none
private
public :: nasto_common_object

type :: nasto_common_object
   !< Navier-Stokes equations system class definition, common data to all backends.

   ! ADAM library objects
   type(mpih_object)           :: mpih          !< MPI handler.
   type(adam_object)           :: adam          !< ADAM.
   type(field_object), pointer :: field=>null() !< The field.
   type(grid_object),  pointer :: grid=>null()  !< The grid.
   ! grid/field data replica for easy handling
   integer(I4P), pointer :: ngc=>null()           !< Number of ghost cells.
   integer(I4P), pointer :: ni=>null()            !< Number of cells in i direction.
   integer(I4P), pointer :: nj=>null()            !< Number of cells in j direction.
   integer(I4P), pointer :: nk=>null()            !< Number of cells in k direction.
   integer(I4P), pointer :: nb=>null()            !< Total blocks number for MPI.
   integer(I4P), pointer :: blocks_number=>null() !< Actual blocks number.
   integer(I4P), pointer :: nv=>null()            !< Number of variables.
   integer(I4P), pointer :: nv_aux=>null()        !< Number of auxiliary variables.
   ! IO data
   type(file_ini)                  :: file_input         !< Nasto input file handler.
   logical                         :: save_memory_status !< Flag to activate memory status saving.
   integer(I4P)                    :: slices_number=0    !< Number of slices to be save.
   type(slice_object), allocatable :: slice(:)           !< Slices data.
   ! numerics data
   real(R8P), allocatable    :: fd_coeff1(:)            !< First order derivatives coeffs.
   real(R8P), allocatable    :: fd_coeff2(:)            !< Second order derivatives coeffs.
   real(R8P), allocatable    :: fd_conv(:,:)            !< Second order derivatives coeffs.
   integer(I4P)              :: visc_scheme=2_I4P       !< Laplacian viscosity scheme.
   integer(I4P)              :: euler_scheme=2_I4P      !< Centered (1) or weno (2) euler scheme.
   integer(I4P)              :: central_order=4_I4P     !< Centered euler scheme order.
   integer(I4P)              :: lmax=2_I4P              !< Central convective half stencil.
   integer(I4P)              :: weno_n_ror=1_I4P        !< Number of ror weno (1 disables ror).
   integer(I4P), allocatable :: weno_schemes(:)         !< Weno schemes.
   real(R8P)                 :: ror_threshold           !< Threshold for ror check.
   integer(I4P)              :: ror_n_indexes=2_I4P     !< Number of variables checked by ror.
   integer(I4P), allocatable :: ror_indexes(:)          !< Ror variable indexes.
   integer(I4P)              :: enable_ror_stats=0_I4P  !< Ror stats (0=disable, 1=enable)
   integer(I4P), allocatable :: ror_stats(:,:,:,:,:)    !< ROR statistics.
   integer(I4P)              :: reduction_extent        !< Length of stencil to consider to reduce weno order close to solids.
   integer(I4P)              :: reduced_order           !< Weno reduced order close to solids.
   integer(I4P), allocatable :: order_modify(:,:,:,:,:) !< Modified order close to solids.
   integer(I4P)              :: iweno=2_I4P             !< WENO order.
   integer(I4P)              :: visc_order=4_I4P        !< Laplacian viscosity order.
   integer(I4P)              :: nrk=4_I4P               !< Runge-Kutta stages number.
   real(R8P), allocatable    :: ark(:)                  !< Runge-Kutta alpha coefficients.
   real(R8P), allocatable    :: brk(:)                  !< Runge-Kutta beta coefficients.
   real(R8P), allocatable    :: crk(:)                  !< Runge-Kutta beta coefficients.
   ! Immersed Boundary (IB) data
   character(999)           :: solid_name     !< Name of solid off file.
   integer(I4P)             :: solid_bc_type  !< Solid bc.
   integer(I4P)             :: n_solids=0     !< Number of solids (only 1 supported now).
   type(c_ptr), allocatable :: ptree(:)       !< CGAL trees for solids.
   real(R8P),   allocatable :: phi(:,:,:,:,:) !< Distance function.
   ! Adaptive Mesh Refinement (AMR) data
   integer(I4P)                         :: amr_iters      !< AMR number of iterations.
   integer(I4P)                         :: amr_frequency  !< AMR time step interval.
   integer(I4P)                         :: amr_n_markers  !< AMR number of markers.
   type(amr_marker_object), allocatable :: amr_markers(:) !< AMR array of marker objects.
   ! simulation iterations data
   integer(I4P)   :: it = 0           !< Time iteration counter.
   real(R8P)      :: time = 0.0_R8P   !< Time.
   logical        :: restart=.false.  !< Restart flag.
   character(999) :: restart_basename !< Restart file basename.
   integer(I4P)   :: restart_save     !< Iteration interval between subsequent restart saves.
   real(R8P)      :: time_max         !< Maximum time of run.
   integer(I4P)   :: t_max            !< Maximum number of iterations of run.
   real(R8P)      :: time_save        !< Time interval between subsequent saves.
   integer(I4P)   :: n_save           !< Iteration interval between subsequent saves.
   character(999) :: output_basename  !< Output file basename.
   real(R8P)      :: CFL              !< CFL time limit.
   real(R8P)      :: dt=0.0001_R8P    !< Maximum time step accordingly to CFL criterion.
   ! Boundary Conditions (BC) data
   integer(I4P)              :: bc_type(6)                     !< Boundary condition type.
   real(R8P)                 :: bc_vars(BC_VARS_NUMBER_MAX, 6) !< Variables' array for boundary conditions.
   integer(I4P), allocatable :: bcs_type(:)                    !< Immersed boundary condition type.
   real(R8P),    allocatable :: bcs_vars(:, :)                 !< Variables' array for immersed boundary conditions.
   ! Fields data: see nasto parameters definition for the arrangement of conservative and auxiliary variables
   real(R8P), allocatable :: q_aux(:,:,:,:,:) !< Auxiliary cell centered variables.
   ! fluids physics
   integer(I4P)                            :: ns=1_I4P          !< Number of fluid species.
   type(nasto_physics_object), allocatable :: fluids_physics(:) !< Fluids physiscs.
   ! initial conditions
   type(nasto_ic_object) :: ic !< Initial Conditions.
   ! integer(I4P) :: ic_type                     !< Initial condition type.
   ! real(R8P)    :: ic_vars(IC_VARS_NUMBER_MAX) !< Variables' array for initial conditions.
   ! real(R8P)    :: ic_vars(12                ) !< Variables' array for initial conditions.
   contains
      procedure, pass(self) :: allocate_common            !< Allocate common data.
      procedure, pass(self) :: initialize_common          !< Initialize the equation common data.
      procedure, pass(self) :: initialize_fd_coefficients !< Initialize Finite Difference coefficients.
      procedure, pass(self) :: runge_kutta_initialize     !< Initialize Runge-Kutta data.
endtype nasto_common_object

contains
   subroutine allocate_common(self)
   !< Allocate common data.
   class(nasto_common_object), intent(inout) :: self !< The equation.

   associate(nv=>self%nv, ns=>self%ns, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk, &
             nb=>self%nb, nrk=>self%nrk, nv_aux=>self%nv_aux, n_solids=>self%n_solids, iweno=>self%iweno)
   allocate(self%q_aux(1:nv_aux, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nb))
   self%q_aux = 0._R8P
   if (self%n_solids > 0) then
      allocate(self%phi(1:nb, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:n_solids))
      self%phi = -1._R8P
   endif
   allocate(self%order_modify(1:nb, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:3))
   self%order_modify = 0
   if (self%enable_ror_stats > 0) allocate(self%ror_stats(1:nb, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:3))
   endassociate
   endsubroutine allocate_common

   subroutine initialize_common(self, filename, nv, nb, nodes_number, do_mpi_init)
   !< Initialize the equation common data.
   class(nasto_common_object), intent(inout)        :: self         !< The equation.
   character(*),               intent(in)           :: filename     !< Input file name.
   integer(I8P),               intent(in)           :: nodes_number !< Allocated nodes on tree.
   integer(I4P),               intent(in)           :: nb           !< Number of allocated blocks.
   integer(I4P),               intent(in)           :: nv           !< Number of evolved variables.
   logical,                    intent(in), optional :: do_mpi_init  !< Flag to activate MPI init call.

   call self%mpih%initialize(do_mpi_init=do_mpi_init)
   print '(A)', self%mpih%myrankstr//'nasto_common_object%initialize start'
   call self%file_input%initialize(filename=trim(filename))
   call self%file_input%load
   call self%adam%grid%initialize(file_parameters=self%file_input, verbose=.true.)
   self%bc_type = self%adam%grid%bc_type
   call self%adam%initialize(file_parameters=self%file_input, &
                             do_tree_init=.true.,             &
                             do_field_init=.true.,            &
                             nv=nv, nb=nb, nodes_number=nodes_number)
   call associate_adam_data(grid=self%adam%grid, field=self%adam%field)
   call self%adam%refine_uniform(refinement_levels=self%adam%tree%iu_ref_levels, do_blocks_reorder=.false.)
   call self%adam%prune(ijkl_prune=self%adam%tree%ijkl_prune, do_blocks_reorder=.false.)
   call load_amr_from_ini_file
   call load_bc_from_ini_file
   call self%ic%initialize(file_parameters=self%file_input)
   call load_physics_from_ini_file
   call load_schemes_from_ini_file
   call load_slices_from_ini_file
   call load_solids_from_ini_file
   call load_timing_from_ini_file
   call self%initialize_fd_coefficients
   call self%runge_kutta_initialize
   call self%allocate_common
   print '(A)', self%mpih%myrankstr//'nasto_common_object%initialize finish'
   contains
      subroutine associate_adam_data(grid, field)
      !< Associate main ADAM data to equation for easy handling.
      type(grid_object),  intent(in), target :: grid  !< The grid.
      type(field_object), intent(in), target :: field !< The field.

      self%grid          => grid
      self%field         => field
      self%blocks_number => field%blocks_number
      self%ni            => field%grid%ni
      self%nj            => field%grid%nj
      self%nk            => field%grid%nk
      self%ngc           => field%grid%ngc
      self%nb            => field%nb
      self%nv            => field%nv

      allocate(self%nv_aux) ; self%nv_aux = 9
      endsubroutine associate_adam_data

      subroutine load_amr_from_ini_file
      !< Parse AMR setting from input file.
      integer(I4P)   :: buf_I4   !< I4 buffer.
      real(R8P)      :: buf_R8   !< R8 buffer.
      character(999) :: sname    !< Section name.
      integer(I4P)   :: i_marker !< Counter.
      integer(I4P)   :: mode     !< AMR mode.

      call self%file_input%get(section_name='amr', option_name='frequency', val=buf_I4) ; self%amr_frequency = buf_I4
      call self%file_input%get(section_name='amr', option_name='iters',     val=buf_I4) ; self%amr_iters     = buf_I4
      call self%file_input%get(section_name='amr', option_name='n_markers', val=buf_I4) ; self%amr_n_markers = buf_I4
      allocate(self%amr_markers(self%amr_n_markers))
      do i_marker=1,self%amr_n_markers
         sname = 'amr_marker_'//trim(str(i_marker,.true.))
         call self%file_input%get(section_name=sname, option_name='mode', val=mode)
         self%amr_markers(i_marker)%mode = mode
         call self%file_input%get(section_name=sname, option_name='delta_fine', val=buf_R8)
         self%amr_markers(i_marker)%delta_fine = buf_R8
         call self%file_input%get(section_name=sname, option_name='delta_coarse', val=buf_R8)
         self%amr_markers(i_marker)%delta_coarse = buf_R8
         if(mode == 1) then
            call self%file_input%get(section_name=sname, option_name='solid', val=buf_I4)
            self%amr_markers(i_marker)%solid = buf_I4
         elseif(mode == 2) then
            call self%file_input%get(section_name=sname, option_name='var', val=buf_I4)
            self%amr_markers(i_marker)%ivar = buf_I4
            call self%file_input%get(section_name=sname, option_name='tol', val=buf_R8)
            self%amr_markers(i_marker)%tol = buf_R8
         endif
      enddo
      endsubroutine load_amr_from_ini_file

      subroutine load_bc_from_ini_file
      !< Parse boundary conditions setting from input file.
      real(R8P)      :: buf_R8       !< R8 buffer.
      character(999) :: sname        !< Section name.
      character(999) :: snames_bc(6) !< Section names bc.
      integer(I4P)   :: bc_type_item !< Boundary condition type element.
      integer(I4P)   :: i_var, i_bc  !< Counter.
      integer(I4P)   :: n_vars       !< Number of vars.

      snames_bc(:) = ["bc_x_min", "bc_x_max", "bc_y_min", "bc_y_max", "bc_z_min", "bc_z_max"]
      do i_bc=1,6
         sname = snames_bc(i_bc)
         call self%file_input%get(section_name=sname, option_name='type', val=bc_type_item)
         n_vars = BC_VARS_NUMBER(bc_type_item)
         do i_var=1,n_vars
             call self%file_input%get(section_name=sname, option_name="var"//trim(str(i_var,.true.)), val=buf_R8)
             self%bc_vars(i_var, i_bc) = buf_R8
         enddo
      enddo
      endsubroutine load_bc_from_ini_file

      subroutine load_physics_from_ini_file
      !< Parse fluids physics setting from input file.
      integer(I4P) :: s !< Counter.

      call self%file_input%get(section_name='physics', option_name='ns', val=self%ns, error=self%mpih%error)
      if (self%mpih%error>0) call self%mpih%error_stop(msg=': failed to load [physics].(ns)')
      if (self%ns<1) call self%mpih%error_stop(msg=': error [physics].(ns) must be >=1')

      allocate(self%fluids_physics(1:self%ns))
      do s=1, self%ns
         call self%fluids_physics(s)%load_from_file(fini=self%file_input, s=s)
      enddo
      endsubroutine load_physics_from_ini_file

      subroutine load_schemes_from_ini_file
      !< Parse schemes setting from input file.
      integer(I4P)              :: buf_I4          !< I4 buffer.
      integer(I4P), allocatable :: buf_array_I4(:) !< I4 buffer array.
      real(R8P)                 :: buf_R8          !< R8 buffer.
      integer(I4P)              :: i               !< Counter.
      character(128)            :: oname           !< Option name buffer.

      call self%file_input%get(section_name='schemes', option_name='euler_scheme', val=buf_I4) ; self%euler_scheme = buf_I4
      call self%file_input%get(section_name='schemes', option_name='central_order' , val=buf_I4) ; self%central_order  = buf_I4
      self%lmax  = (self%central_order)/2
      call self%file_input%get(section_name='schemes', option_name='weno_n_ror' , val=buf_I4) ; self%weno_n_ror  = buf_I4
      allocate(buf_array_I4(1:self%weno_n_ror))
      self%weno_schemes = buf_array_I4
      do i=1,self%weno_n_ror
         oname = 'weno_schemes_'//trim(str(i,.true.))
         call self%file_input%get(section_name='schemes', option_name=trim(oname), val=buf_I4)
         self%weno_schemes(i) = buf_I4
      enddo
      deallocate(buf_array_I4)
      self%iweno = IWENO_FROM_SCHEME(self%weno_schemes(1))
      call self%file_input%get(section_name='schemes', option_name='ror_threshold' , val=buf_R8) ; self%ror_threshold = buf_R8
      call self%file_input%get(section_name='schemes', option_name='ror_n_indexes' , val=buf_I4) ; self%ror_n_indexes = buf_I4
      allocate(buf_array_I4(1:self%ror_n_indexes))
      self%ror_indexes = buf_array_I4
      do i=1,self%ror_n_indexes
         oname = 'ror_indexes_'//trim(str(i,.true.))
         call self%file_input%get(section_name='schemes', option_name=trim(oname), val=buf_I4)
         self%ror_indexes(i) = buf_I4
      enddo
      deallocate(buf_array_I4)
      call self%file_input%get(section_name='schemes', option_name='enable_ror_stats' , val=buf_I4) ; self%enable_ror_stats = buf_I4
      call self%file_input%get(section_name='schemes', option_name='reduction_extent' , val=buf_I4) ; self%reduction_extent = buf_I4
      call self%file_input%get(section_name='schemes', option_name='reduced_order'    , val=buf_I4) ; self%reduced_order    = buf_I4
      call self%file_input%get(section_name='schemes', option_name='visc_scheme' , val=buf_I4) ; self%visc_scheme = buf_I4
      call self%file_input%get(section_name='schemes', option_name='visc_order'  , val=buf_I4) ; self%visc_order  = buf_I4
      endsubroutine load_schemes_from_ini_file

      subroutine load_slices_from_ini_file
      !< Parse slices setting from input file.
      character(999) :: buf_CHAR   !< String buffer.
      integer(I4P)   :: buf_I4     !< I4 buffer.
      real(R8P)      :: buf_R8     !< R8 buffer.
      character(999) :: sname      !< Section name.
      real(R8P)      :: dxyz(3)    !< Space steps.
      integer(I4P)   :: i, j, k, s !< Counter.

      call self%file_input%get(section_name="time", option_name="slices_number", val=buf_I4) ; self%slices_number = buf_I4
      if (self%slices_number > 0) then
         allocate(self%slice(self%slices_number))
         do s=1, self%slices_number
            sname = 'slice_'//trim(str(s,.true.))
            call self%file_input%get(section_name=sname, option_name='slice_itype',  val=buf_CHAR)
            self%slice(s)%slice_itype  =buf_CHAR
            call self%file_input%get(section_name=sname, option_name='slice_save',   val=buf_I4)
            self%slice(s)%slice_save   =buf_I4
            call self%file_input%get(section_name=sname, option_name='slice_ni',     val=buf_I4)
            self%slice(s)%slice_nijk(1)=buf_I4
            call self%file_input%get(section_name=sname, option_name='slice_nj',     val=buf_I4)
            self%slice(s)%slice_nijk(2)=buf_I4
            call self%file_input%get(section_name=sname, option_name='slice_nk',     val=buf_I4)
            self%slice(s)%slice_nijk(3)=buf_I4
            call self%file_input%get(section_name=sname, option_name='slice_emin_x', val=buf_R8)
            self%slice(s)%slice_emin(1)=buf_R8
            call self%file_input%get(section_name=sname, option_name='slice_emin_y', val=buf_R8)
            self%slice(s)%slice_emin(2)=buf_R8
            call self%file_input%get(section_name=sname, option_name='slice_emin_z', val=buf_R8)
            self%slice(s)%slice_emin(3)=buf_R8
            call self%file_input%get(section_name=sname, option_name='slice_emax_x', val=buf_R8)
            self%slice(s)%slice_emax(1)=buf_R8
            call self%file_input%get(section_name=sname, option_name='slice_emax_y', val=buf_R8)
            self%slice(s)%slice_emax(2)=buf_R8
            call self%file_input%get(section_name=sname, option_name='slice_emax_z', val=buf_R8)
            self%slice(s)%slice_emax(3)=buf_R8
            allocate(self%slice(s)%slice_points(3,self%slice(s)%slice_nijk(1),&
                                                  self%slice(s)%slice_nijk(2),&
                                                  self%slice(s)%slice_nijk(3)))
            dxyz(1) = (self%slice(s)%slice_emax(1) - self%slice(s)%slice_emin(1)) / self%slice(s)%slice_nijk(1)
            dxyz(2) = (self%slice(s)%slice_emax(2) - self%slice(s)%slice_emin(2)) / self%slice(s)%slice_nijk(2)
            dxyz(3) = (self%slice(s)%slice_emax(3) - self%slice(s)%slice_emin(3)) / self%slice(s)%slice_nijk(3)
            do k=1, self%slice(s)%slice_nijk(3)
               do j=1, self%slice(s)%slice_nijk(2)
                  do i=1, self%slice(s)%slice_nijk(1)
                     self%slice(s)%slice_points(1,i,j,k) = self%slice(s)%slice_emin(1) + (i - 0.5_R8P) * dxyz(1)
                     self%slice(s)%slice_points(2,i,j,k) = self%slice(s)%slice_emin(2) + (j - 0.5_R8P) * dxyz(2)
                     self%slice(s)%slice_points(3,i,j,k) = self%slice(s)%slice_emin(3) + (k - 0.5_R8P) * dxyz(3)
                  enddo
               enddo
            enddo
         enddo
      endif
      endsubroutine load_slices_from_ini_file

      subroutine load_solids_from_ini_file
      !< Parse immersed boundary solids setting from input file.
      character(999) :: buf_CHAR !< String buffer.
      integer(I4P)   :: buf_I4   !< I4 buffer.
      character(999) :: sname    !< Section name.
      integer(I4P)   :: i_solid  !< Counter.

      call self%file_input%get(section_name='solids', option_name='n_solids', val=buf_I4)
      self%n_solids = buf_I4
      allocate(self%bcs_type(self%n_solids))
      allocate(self%bcs_vars(BCS_VARS_NUMBER_MAX, self%n_solids))
      allocate(self%ptree(self%n_solids))
      do i_solid=1,self%n_solids
         sname = 'solid_'//trim(str(i_solid,.true.))
         call self%file_input%get(section_name=sname, option_name='name', val=buf_CHAR)
         self%solid_name = buf_CHAR
         call self%file_input%get(section_name=sname, option_name='bcs_type', val=buf_I4)
         self%solid_bc_type = buf_I4
         self%bcs_type(i_solid) = self%solid_bc_type
         ! RIMETTERE CGAL
         ! call cgal_polyhedron_read(self%ptree(i_solid), self%solid_name)
      enddo
      endsubroutine load_solids_from_ini_file

      subroutine load_timing_from_ini_file
      !< Parse timing setting from input file.
      logical        :: buf_BOOL !< Logical buffer.
      character(999) :: buf_CHAR !< String buffer.
      integer(I4P)   :: buf_I4   !< I4 buffer.
      real(R8P)      :: buf_R8   !< R8 buffer.

      call self%file_input%get(section_name="time", option_name="restart",         val=buf_BOOL) ; self%restart          = buf_BOOL
      call self%file_input%get(section_name="time", option_name="restart_basename",val=buf_CHAR) ; self%restart_basename = buf_CHAR
      call self%file_input%get(section_name="time", option_name="restart_save",    val=buf_I4)   ; self%restart_save     = buf_I4
      call self%file_input%get(section_name="time", option_name="time_max",        val=buf_R8)   ; self%time_max         = buf_R8
      call self%file_input%get(section_name="time", option_name="t_max",           val=buf_I4)   ; self%t_max            = buf_I4
      call self%file_input%get(section_name="time", option_name="time_save",       val=buf_R8)   ; self%time_save        = buf_R8
      call self%file_input%get(section_name="time", option_name="n_save",          val=buf_I4)   ; self%n_save           = buf_I4
      call self%file_input%get(section_name="time", option_name="output_basename", val=buf_CHAR) ; self%output_basename  = buf_CHAR
      call self%file_input%get(section_name='time', option_name='CFL',             val=buf_R8)   ; self%CFL              = buf_R8

      call self%file_input%get(section_name='equation', option_name='save_memory_status', val=buf_BOOL)
      self%save_memory_status = buf_BOOL
      endsubroutine load_timing_from_ini_file
   endsubroutine initialize_common

   subroutine initialize_fd_coefficients(self)
   !< Initialize Finite Difference coefficients.
   class(nasto_common_object), intent(inout) :: self !< The equation.

   allocate(self%fd_conv(4,4), self%fd_coeff1(3), self%fd_coeff2(0:3))
   associate(fd_conv=>self%fd_conv, fd_coeff1=>self%fd_coeff1,fd_coeff2=>self%fd_coeff2)
   ! Coefficients for computation of convective terms
   fd_conv(1,1) = 1._R8P/2._R8P

   fd_conv(1,2) =  2._R8P/3._R8P
   fd_conv(2,2) = -1._R8P/12._R8P

   fd_conv(1,3) =  3._R8P/4._R8P
   fd_conv(2,3) = -3._R8P/20._R8P
   fd_conv(3,3) =  1._R8P/60._R8P

   fd_conv(1,4) =  4._R8P/5._R8P
   fd_conv(2,4) = -1._R8P/5._R8P
   fd_conv(3,4) =  4._R8P/105._R8P
   fd_conv(4,4) = -1._R8P/280._R8P

   ! Coefficients for computation of viscous terms
   select case (self%visc_order/2)
   case (1)
    fd_coeff1(1) = 0.5_R8P
   case (2)
    fd_coeff1(1) = 2._R8P/3._R8P
    fd_coeff1(2) = -1._R8P/12._R8P
   case (3)
    fd_coeff1(1) = 0.75_R8P
    fd_coeff1(2) = -0.15_R8P
    fd_coeff1(3) = 1._R8P/60._R8P
   end select

   select case (self%visc_order/2)
   case (1)
    fd_coeff2(0) = -2._R8P
    fd_coeff2(1) =  1._R8P
   case (2)
    fd_coeff2(0) = -2.5_R8P
    fd_coeff2(1) = 4._R8P/3._R8P
    fd_coeff2(2) = -1._R8P/12._R8P
   case (3)
    fd_coeff2(0) = -245._R8P/90._R8P
    fd_coeff2(1) = 1.5_R8P
    fd_coeff2(2) = -0.15_R8P
    fd_coeff2(3) = 1._R8P/90._R8P
   endselect
   endassociate
   endsubroutine initialize_fd_coefficients

   subroutine runge_kutta_initialize(self)
   !< Initialize Runge-Kutta data.
   class(nasto_common_object), intent(inout) :: self !< The equation.

   call self%file_input%get(section_name='time', option_name='nrk', val=self%nrk)
   allocate(self%ark(self%nrk), self%brk(self%nrk), self%crk(self%nrk))
   select case(self%nrk)
      case(1_I4P) ! Eulero
         self%ark(1) = 1d0  ; self%brk(1) = 0d0; self%crk(1) = 1d0
      case(2_I4P) ! secondo ordine TVD
         self%ark(1) = 1d0    ; self%brk(1) = 0d0  ; self%crk(1) = 1d0
         self%ark(2) = 0.5d0  ; self%brk(2) = 0.5d0; self%crk(2) = 0.5d0
      case(3_I4P) ! terzo ordine TVD
         self%ark(1) = 1d0     ; self%brk(1) = 0d0     ; self%crk(1) = 1d0
         self%ark(2) = 0.75d0  ; self%brk(2) = 0.25d0  ; self%crk(2) = 0.25d0
         self%ark(3) = 1d0/3d0 ; self%brk(3) = 2d0/3d0 ; self%crk(3) = 2d0/3d0
   endselect
   endsubroutine runge_kutta_initialize
endmodule adam_nasto_common_object
