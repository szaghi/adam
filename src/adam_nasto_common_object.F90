!< ADAM, Navier-Stokes equations system class definition, common data to all backends.
module adam_nasto_common_object
!< ADAM, Navier-Stokes equations system class definition, common data to all backends.

use adam_adam_object
use adam_amr_marker_object
use adam_field_object
use adam_grid_object
use adam_mpih_object
use adam_slice_object
! use adam_parameters
use FiNeR
use PENF
use ISO_C_BINDING

implicit none
private
public :: nasto_common_object

! the following 2 parameters are already defined into adam_parameters module, but the compiler is too much stupid (see the use
! module commented above) and do not accept to take them from the module...
integer(I4P), parameter :: IC_VARS_NUMBER_MAX = 12
integer(I4P), parameter :: BC_VARS_NUMBER_MAX = 5

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
   integer(I4P)              :: ns=1_I4P                !< Number of fluid species.
   integer(I4P)              :: visc_law=0_I4P          !< Diffusivity type (0=constant, 1=power, 2=Sutherland).
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
   ! Initial Conditions (IC) data
   integer(I4P) :: ic_type                     !< Initial condition type.
   real(R8P)    :: ic_vars(IC_VARS_NUMBER_MAX) !< Variables' array for initial conditions.
   ! real(R8P)    :: ic_vars(12                ) !< Variables' array for initial conditions.
   ! Boundary Conditions (BC) data
   integer(I4P)              :: bc_type(6)                     !< Boundary condition type.
   real(R8P)                 :: bc_vars(BC_VARS_NUMBER_MAX, 6) !< Variables' array for boundary conditions.
   integer(I4P), allocatable :: bcs_type(:)                    !< Immersed boundary condition type.
   real(R8P),    allocatable :: bcs_vars(:, :)                 !< Variables' array for immersed boundary conditions.
   ! physics data
   real(R8P) :: Lewis=1._R8P        !< Lewis number.
   real(R8P) :: Zeldovich=1060._R8P !< Zeldovich number.
   real(R8P) :: Damkohler=1800._R8P !< Damkohler number.
   real(R8P) :: gamma_fluid=1.4_R8P !< Gamma.
   real(R8P) :: R_star=287._R8P     !< Gas constant.
   real(R8P) :: cv_star=714._R8P    !< Constant volume specific heat.
   real(R8P) :: mu_star=0.001_R8P   !< Dynamic viscosity.
   real(R8P) :: cp_star=1000._R8P   !< Constant pressure specific heat.
   real(R8P) :: k_star=0.0013_R8P   !< Thermal diffusivity.
   real(R8P) :: dha_star=10000._R8P !< Entalpy formation.
   ! Fields data: see nasto parameters definition for the arrangement of conservative and auxiliary variables
   real(R8P), allocatable :: q_aux(:,:,:,:,:) !< Auxiliary cell centered variables.
   contains
      procedure, pass(self) :: initialize !< Initialize the equation.
endtype nasto_common_object

contains
   subroutine initialize(self, filename)
   !< Initialize the equation.
   class(nasto_common_object), intent(inout) :: self     !< The equation.
   character(*),               intent(in)    :: filename !< Input file name.

   call self%mpih%initialize(do_mpi_init=.true.)
   print '(A)', self%mpih%myrankstr//'nasto_common_object%initialize start'
   call self%file_input%initialize(filename=trim(filename))
   call self%file_input%load
   call self%adam%grid%initialize(file_parameters=self%file_input, verbose=.true.)
   self%bc_type = self%adam%grid%bc_type
   endsubroutine initialize
endmodule adam_nasto_common_object
