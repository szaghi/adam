!< ADAM, Euler equations system class definition, GPU backend.
module adam_equation_nasto_gpu_object
!< ADAM, Euler equations system class definition, GPU backend.

use adam_adam_object
use adam_base_gpu_object
use adam_field_object
use adam_grid_object
use adam_parameters
use adam_tree_node_object, only : tree_node_object
use adam_weno_library_gpu
use FiNeR
use PENF
use MPI
use CUDAFOR
use cgal_wrappers
use ISO_C_BINDING
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use memorysaver

implicit none
private
public :: equation_nasto_gpu_object

integer(I4P), parameter :: IC_UNIFORM         = 1_I4P
integer(I4P), parameter :: IC_LEFTRIGHT       = 2_I4P
integer(I4P), parameter :: IC_FLAME           = 3_I4P
integer(I4P), parameter :: IC_VARS_NUMBER(3)  = [6, 11, 12]
integer(I4P), parameter :: IC_VARS_NUMBER_MAX = 12 !maxval(IC_VARS_NUMBER)    !< Maximum number of variables needed for IC.

integer(I4P), parameter :: BC_EXTRAPOLATION   = 1_I4P
integer(I4P), parameter :: BC_INFLOW          = 2_I4P
integer(I4P), parameter :: BC_VARS_NUMBER(2)  = [0, 5]
integer(I4P), parameter :: BC_VARS_NUMBER_MAX = 5 !maxval(BC_VARS_NUMBER)    !< Maximum number of variables needed for BC.

integer(I4P), parameter :: BCS_VISCOUS  = 1_I4P
integer(I4P), parameter :: BCS_EULER = 2_I4P
integer(I4P), parameter :: BCS_VARS_NUMBER(2)  = [3, 0]
integer(I4P), parameter :: BCS_VARS_NUMBER_MAX = 3 !maxval(BCS_VARS_NUMBER)  !< Maximum number of variables needed for BCS.

! q_aux indexes
integer(I4P), parameter :: IRHO  = 1
integer(I4P), parameter :: IU    = 2
integer(I4P), parameter :: IV    = 3
integer(I4P), parameter :: IW    = 4
integer(I4P), parameter :: IYA   = 5
integer(I4P), parameter :: ITEM  = 6
integer(I4P), parameter :: IPRES = 7
integer(I4P), parameter :: IENTA = 8
integer(I4P), parameter :: ISOUN = 9

type :: amr_marker_obj
   !< AMR marker object.
   integer(I4P) :: mode         !< Marker mode.
   integer(I4P) :: solid        !< Solid number.
   real(R8P)    :: delta_fine   !< Fine cell space step.
   real(R8P)    :: delta_coarse !< Coarse cell space step.
   integer(I4P) :: ivar         !< ivar.
   real(R8P)    :: tol          !< Tolerance.
endtype amr_marker_obj

type :: slice_obj
   !< Slice object.
   character(99)          :: slice_itype           !< Slice interpolation type.
   integer(I4P)           :: slice_save            !< Iteration interval between subsequent data-slice saves.
   integer(I4P)           :: slice_nijk(3)         !< Slice number of points.
   real(R8P)              :: slice_emin(3)         !< Slice minimum extents.
   real(R8P)              :: slice_emax(3)         !< Slice maximum extents.
   real(R8P), allocatable :: slice_points(:,:,:,:) !< Slice points coordinates [3,ni,nj,nk].
endtype slice_obj

type :: equation_nasto_gpu_object
   !< Flame single-step class definition, GPU backend.
   !<
   !< Arrenhius form for reaction term.
   !<
   !< The conservative variables are arranged as follows:
   !<```
   !< q(1): rho
   !< q(2): rho * u
   !< q(3): rho * v
   !< q(4): rho * w
   !< q(5): rho * E
   !< q(6): rho * Ya, specific density of specie
   !<```
   !< q_aux(1): rho
   !< q_aux(2): u
   !< q_aux(3): v
   !< q_aux(4): w
   !< q_aux(5): ya
   !< q_aux(6): tem
   !< q_aux(7): pressure
   !< q_aux(8): entalpy
   !< q_aux(9): sound speed
   !<```
   type(file_ini)              :: file_input              !< Nasto input file handler.
   type(adam_object)           :: adam                    !< ADAM.
   type(field_object), pointer :: field=>null()           !< The field.
   type(grid_object),  pointer :: grid=>null()            !< The grid.
   integer(I4P)                :: ngc                     !< Number of ghost cells.
   integer(I4P)                :: ni                      !< Number of cells in i direction.
   integer(I4P)                :: nj                      !< Number of cells in j direction.
   integer(I4P)                :: nk                      !< Number of cells in k direction.
   integer(I4P)                :: nb                      !< Total blocks number for MPI.
   integer(I4P), pointer       :: blocks_number=>null()   !< Actual blocks number.
   integer(I4P)                :: nv                      !< Number of variables.
   integer(I4P)                :: nv_aux                  !< Number of auxiliary variables.
   type(base_gpu_object)       :: base_gpu                !< The base GPU handler.
   integer(I4P)                :: myrank=0_I4P            !< MPI rank process.
   integer(I4P)                :: procs_number=1_I4P      !< Number of MPI processes.
   integer(I4P)                :: error=0_I4P             !< Error traping flag.
   integer(I4P)                :: field_gpu_number=12_I4P !< Number of nv fields used in memory.
   ! equation data
   real(R8P), allocatable      :: fd_coeff1(:)       !< First order derivatives coeffs.
   real(R8P), allocatable      :: fd_coeff2(:)       !< Second order derivatives coeffs.
   real(R8P), allocatable      :: fd_conv(:,:)       !< Second order derivatives coeffs.
   integer(I4P)                :: visc_scheme=2_I4P  !< Laplacian viscosity scheme.
   integer(I4P)                :: euler_scheme=2_I4P !< Centered euler scheme scheme.
   integer(I4P)                :: visc_order=4_I4P   !< Laplacian viscosity order.
   integer(I4P)                :: euler_order=4_I4P  !< Centered euler scheme order.
   integer(I4P)                :: ns=1_I4P           !< Number of fluid species.
   integer(I4P)                :: iweno=2_I4P        !< WENO order.
   integer(I4P)                :: lmax=2_I4P         !< Central convective half stencil.
   integer(I4P)                :: visc_law=0_I4P     !< Diffusivity type (0=constant, 1=power, 2=Sutherland)
   ! Runge-Kutta data
   integer(I4P)                      :: nrk=4_I4P       !< Runge-Kutta stages number.
   real(R8P), allocatable            :: ark(:)          !< RK alpha coefficients.
   real(R8P), allocatable            :: brk(:)          !< RK beta coefficients.
   ! Immersed boundary
   character(999)                    :: solid_name      !< Name of solid off file.
   integer(I4P)                      :: solid_bc_type   !< Solid bc.
   integer(I4P)                      :: n_solids=0      !< Number of solids (only 1 supported now).
   type(c_ptr), allocatable          :: ptree(:)        !< CGAL trees for solids.
   real(R8P),    allocatable         :: phi(:,:,:,:,:)  !< Distance function.
   ! AMR
   integer(I4P)                      :: amr_iters       !< AMR number of iterations.
   integer(I4P)                      :: amr_frequency   !< AMR time step interval.
   integer(I4P)                      :: amr_n_markers   !< AMR number of markers.
   type(amr_marker_obj), allocatable :: amr_markers(:)  !< AMR array of marker objects.
   ! Time
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
   ! Slices
   integer(I4P)                 :: slices_number=0 !< Number of slices to be save.
   type(slice_obj), allocatable :: slice(:)        !< Slices data.
   ! Initial conditions
   integer(I4P)           :: ic_type                     !< Initial condition type.
   real(R8P)              :: ic_vars(IC_VARS_NUMBER_MAX) !< Variables' array for initial conditions.
   ! Boundary conditions
   integer(I4P)              :: bc_type(6)                     !< Boundary condition type.
   real(R8P)                 :: bc_vars(BC_VARS_NUMBER_MAX, 6) !< Variables' array for boundary conditions.
   integer(I4P), allocatable :: bcs_type(:)                    !< Immersed boundary condition type.
   real(R8P), allocatable    :: bcs_vars(:, :)                 !< Variables' array for immersed boundary conditions.
   ! Physics
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
   ! Fields
   real(R8P), allocatable         :: q_aux(:,:,:,:,:)        !< Auxiliary cell centered variables.
   real(R8P), allocatable, device :: dq_gpu(:,:,:,:,:)       !< Eikonal right hand side.
   real(R8P), allocatable, device :: fl_gpu(:,:,:,:,:)       !< Residuals.
   real(R8P), allocatable, device :: flx_gpu(:,:,:,:,:)      !< Fluxes along x.
   real(R8P), allocatable, device :: fly_gpu(:,:,:,:,:)      !< Fluxes along y.
   real(R8P), allocatable, device :: flz_gpu(:,:,:,:,:)      !< Fluxes along z.
   real(R8P), allocatable, device :: prhs_gpu(:,:,:,:,:)     !< Prhs for Runge-Kutta.
   real(R8P), allocatable, device :: dxyz_gpu(:,:)           !< Space steps.
   real(R8P), allocatable, device :: fd_coeff1_gpu(:)        !< First order derivatives coeffs.
   real(R8P), allocatable, device :: fd_coeff2_gpu(:)        !< Second order derivatives coeffs.
   real(R8P), allocatable, device :: fd_conv_gpu(:,:)        !< Second order derivatives coeffs.
   real(R8P), allocatable, device :: x_cell_gpu(:,:)         !< Cell positions in x.
   real(R8P), allocatable, device :: y_cell_gpu(:,:)         !< Cell positions in y.
   real(R8P), allocatable, device :: z_cell_gpu(:,:)         !< Cell positions in z.
   real(R8P), allocatable, device :: q_aux_gpu(:,:,:,:,:)    !< Auxiliary cell centered variables.
   real(R8P), allocatable, device :: q_gpu(:,:,:,:,:)        !< Field cell centered variables.
   real(R8P), allocatable, device :: q_invert_gpu(:,:,:,:,:) !< Field cell with boundary set on immersed bodies.
   real(R8P), allocatable, device :: gplus_x(:,:,:,:,:)      !< For weno-x
   real(R8P), allocatable, device :: gminus_x(:,:,:,:,:)     !< For weno-x
   real(R8P), allocatable, device :: gplus_y(:,:,:,:,:)      !< For weno-y
   real(R8P), allocatable, device :: gminus_y(:,:,:,:,:)     !< For weno-y
   real(R8P), allocatable, device :: gplus_z(:,:,:,:,:)      !< For weno-z
   real(R8P), allocatable, device :: gminus_z(:,:,:,:,:)     !< For weno-z
   real(R8P), allocatable, device :: phi_gpu(:,:,:,:,:)      !< Distance function on GPU.
   real(R8P), allocatable, device :: bc_vars_gpu(:, :)       !< Variables' array for boundary conditions on GPU.
   real(R8P), allocatable, device :: bcs_vars_gpu(:, :)      !< Variables' array for immersed boundary on GPU.
   contains
      ! public methods
      procedure, pass(self) :: amr_update              !< Do AMR update.
      procedure, pass(self) :: move_phi                !< Move phi and the actual ptree representation.
      procedure, pass(self) :: compute_aux             !< Compute auxiliary variables.
      procedure, pass(self) :: compute_dt              !< Compute time step.
      procedure, pass(self) :: copy_cpu_gpu            !< Copy data from CPU to GPU.
      procedure, pass(self) :: copy_gpu_cpu            !< Copy data from GPU to CPU.
      procedure, pass(self) :: destroy                 !< Destroy the equation.
      procedure, pass(self) :: fd_initialize           !< Initialize Finite Difference Coefficients.
      procedure, pass(self) :: initialize              !< Initialize the equation.
      procedure, pass(self) :: load_restart_files      !< Load restart files.
      procedure, pass(self) :: mark_by_grad_var        !< Mark blocks to be refined/derefined by a `grad(var)` value.
      procedure, pass(self) :: mark_by_geo             !< Mark blocks to be refined/derefined by a `grad(var)` value.
      procedure, pass(self) :: integrate               !< Runge Kutta integration of equation.
      procedure, pass(self) :: parse_input             !< Parse input file.
      procedure, pass(self) :: print_progress          !< Print simulation progress.
      procedure, pass(self) :: refine_uniform          !< Refine all blocks uniformly.
      procedure, pass(self) :: run                     !< Run all.
      procedure, pass(self) :: update_phi              !< Refine all blocks uniformly.
      procedure, pass(self) :: runge_kutta_initialize  !< Initialize Runge-Kutta data.
      procedure, pass(self) :: save_simulation_data    !< Save all simulation data.
      procedure, pass(self) :: save_restart_files      !< Save restart files.
      procedure, pass(self) :: save_hdf5               !< Save simulation data in HDF5 format.
      procedure, pass(self) :: save_slices             !< Save simulation data slices.
      procedure, pass(self) :: set_boundary_conditions !< Set boundary conditions of equation.
      procedure, pass(self) :: set_initial_conditions  !< Set initial conditions of equation.
      procedure, pass(self) :: update_ghost_gpu        !< Update ghost cells and set boundary conditions.
      procedure, pass(self) :: update_ghost_fluxes_gpu !< Update fluxes cells and set boundary conditions.
      procedure, pass(self) :: compute_residuals_gpu   !< Compute residuals.
      ! operators
      generic :: assignment(=) => eq_assign_eq      !< Overload `=`.
      procedure, pass(lhs), private :: eq_assign_eq !< Operator `=`.
endtype equation_nasto_gpu_object

contains
   ! public methods
   subroutine amr_update(self)
   !< Do AMR update.
   class(equation_nasto_gpu_object), intent(inout) :: self                 !< The equation.
   integer(I4P)                                    :: iterations_          !< Number of AMR iterations, local var.
   logical                                         :: is_grid_changed      !< Flag to check grid changes for each marker.
   logical                                         :: is_grid_changed_all  !< Flag to check grid changes for each iter.
   integer(I4P)                                    :: b, i, j, k, i_marker !< Counter.
   type(amr_marker_obj)                            :: amr_marker           !< Current amr marker.

   amr: do i=1, self%amr_iters
      is_grid_changed_all = .false.
      do i_marker=1, self%amr_n_markers
         amr_marker = self%amr_markers(i_marker)
         call self%update_ghost_gpu(q_gpu=self%q_gpu)
         if(amr_marker%mode == 1) then ! marker "geo"
            call self%mark_by_geo(delta_fine=amr_marker%delta_fine, delta_coarse=amr_marker%delta_coarse)
         elseif(amr_marker%mode == 2) then ! marker "grad"
            call self%mark_by_grad_var(grad_tol=amr_marker%tol, delta_fine=amr_marker%delta_fine, &
                                       delta_coarse=amr_marker%delta_coarse, ivar=amr_marker%ivar)
         endif
         call self%copy_gpu_cpu() ! needed for adam%amr_update
         call self%adam%amr_update(is_marked_by_field=.true., do_blocks_reorder=.false., is_grid_changed=is_grid_changed)
         if(self%n_solids > 0) call self%update_phi()
         call self%copy_cpu_gpu
         is_grid_changed_all = is_grid_changed_all.or.is_grid_changed
      enddo
      if (.not.is_grid_changed_all) then
          print '(A)','AMR Grid stabilized after : '//trim(str(i))//' AMR iterations'
          exit amr
      endif
   enddo amr
   endsubroutine amr_update

   subroutine compute_aux(self, q_gpu, q_aux_gpu)
   !< Compute auxiliary variables.
   class(equation_nasto_gpu_object), intent(in)          :: self          !< The equation.
   real(R8P),                        intent(in),  device :: q_gpu(1:,                    &
                                                                  1-self%ngc:,&
                                                                  1-self%ngc:,&
                                                                  1-self%ngc:,&
                                                                  1:)     !< Conservative variables.
   real(R8P),                        intent(out), device :: q_aux_gpu(1:,                    &
                                                                      1-self%ngc:,&
                                                                      1-self%ngc:,&
                                                                      1-self%ngc:,&
                                                                      1:) !< Auxiliary variables.

   associate(blocks_number=>self%field%blocks_number, ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, ns=>self%ns)
      call compute_aux_cuf(ni=ni, nj=nj, nk=nk, ngc=ngc, ns=ns, blocks_number=blocks_number, &
                           gamma_fluid=self%gamma_fluid, dha_star=self%dha_star, cv_star=self%cv_star, R_star=self%R_star, &
                           q_gpu=q_gpu, q_aux_gpu=q_aux_gpu)
   endassociate
   endsubroutine compute_aux

   subroutine compute_dt(self)
   !< Compute maximum time step accordingly to CFL stabilty criterion.
   class(equation_nasto_gpu_object), intent(inout) :: self !< The equation.
   real(R8P)                                       :: umax !< Maximum speed of waves propagation.
   integer(I4P)                                    :: b    !< Counter.

   associate(blocks_number=>self%field%blocks_number, dxyz=>self%field%dxyz, ni=>self%ni, nj=>self%nj, nk=>self%nk, &
             ngc=>self%ngc, ns=>self%ns, q=>self%field%q, dt=>self%dt, CFL=>self%CFL, mu_star=>self%mu_star)
      call self%compute_aux(q_gpu=self%q_gpu, q_aux_gpu=self%q_aux_gpu)
      dt = huge(1._R8P)
      do b=1, self%field%blocks_number
         call compute_umax_cuf(b, ni=ni, nj=nj, nk=nk, ngc=ngc, ns=ns,   &
                               dx=dxyz(1,b), dy=dxyz(2,b), dz=dxyz(3,b), &
                               q_aux_gpu=self%q_aux_gpu, umax=umax, mu_star=mu_star)
         !dt = min(dt, minval(dxyz(:,b)) / umax * CFL)
         dt = min(dt, CFL / umax )
      enddo
      call MPI_ALLREDUCE(MPI_IN_PLACE, dt, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, self%error)
   endassociate
   endsubroutine compute_dt

   subroutine copy_cpu_gpu(self)
   !< Copy data from CPU to GPU.
   class(equation_nasto_gpu_object), intent(inout) :: self        !< The base backend.
   real(R8P), allocatable                          :: dxyz_t(:,:) !< Space steps transposed.
   integer(I4P)                                    :: i, b        !< Counter.

   call self%base_gpu%copy_transpose_cpu_gpu(q_cpu=self%field%q, q_gpu=self%q_gpu)
   call self%base_gpu%copy_cpu_gpu
   self%x_cell_gpu = self%base_gpu%x_cell_gpu
   self%y_cell_gpu = self%base_gpu%y_cell_gpu
   self%z_cell_gpu = self%base_gpu%z_cell_gpu
   self%dxyz_gpu   = self%base_gpu%dxyz_gpu
   endsubroutine copy_cpu_gpu

   subroutine copy_gpu_cpu(self, compute_q_aux)
   !< Copy data from GPU to CPU.
   class(equation_nasto_gpu_object), intent(inout)        :: self          !< The base backend.
   logical,                          intent(in), optional :: compute_q_aux !< Flag to compute auxiliary variables.

   call self%base_gpu%copy_transpose_gpu_cpu(nv=self%nv, q_gpu=self%q_gpu, q_cpu=self%field%q)
   if (present(compute_q_aux)) then
      if (compute_q_aux) then
         call self%compute_aux(q_gpu=self%q_gpu, q_aux_gpu=self%q_aux_gpu)
         call self%base_gpu%copy_transpose_gpu_cpu(nv=self%nv_aux, q_gpu=self%q_aux_gpu, q_cpu=self%q_aux)
      endif
   endif
   endsubroutine copy_gpu_cpu

   subroutine destroy(self)
   !< Destroy the equation.
   class(equation_nasto_gpu_object), intent(inout) :: self  !< The equation.
   type(equation_nasto_gpu_object)                 :: fresh !< Fresh equation.

   self = fresh
   endsubroutine destroy

   subroutine fd_initialize(self)
   !< Initialize Finite-Difference coefficients.
   class(equation_nasto_gpu_object), intent(inout) :: self !< The equation.

   allocate(self%fd_conv(4,4), self%fd_coeff1(3), self%fd_coeff2(0:3))

   ! Coefficients for finite difference computations
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
   endsubroutine fd_initialize

   subroutine initialize(self, filename)
   !< Initialize the equation.
   class(equation_nasto_gpu_object), intent(inout) :: self              !< The equation.
   character(*),                     intent(in)    :: filename          !< Input file name.
   integer(I4P)                                    :: i, j, k, v, s     !< Counter.
   integer(I4P)                                    :: n_solids          !< Number of immersed bodies.
   integer(I4P)                                    :: n_vars            !< Number of ic/bc vars.
   integer(I4P)                                    :: iu_ref_levels     !< Uniform refinement initial.
   integer(I4P)                                    :: ni, nj, nk, ngc   !< Number of cells along directions.
   integer(I4P)                                    :: mode              !< AMR mode.
   integer(I4P)                                    :: i_prune           !< Pruning along x.
   integer(I4P)                                    :: j_prune           !< Pruning along y.
   integer(I4P)                                    :: k_prune           !< Pruning along z.
   integer(I4P)                                    :: l_prune           !< Pruning level.
   integer(I4P)                                    :: max_level         !< Max refinement level.
   integer(I8P)                                    :: nodes_number      !< Allocated nodes on tree.
   integer(I4P)                                    :: nb                !< Number of allocated blocks.
   integer(I4P)                                    :: nv                !< Number of evolved variables.
   character(999)                                  :: sname             !< Section name.
   character(999)                                  :: oname             !< Option name.
   character(999)                                  :: snames_bc(6)      !< Section names bc.
   integer(I4P)                                    :: bc_type(6)        !< Boundary condition type array.
   integer(I4P)                                    :: bc_type_item      !< Boundary condition type element.
   integer(I4P)                                    :: i_var, i_bc       !< Counter.
   integer(I4P)                                    :: i_solid, i_marker !< Counter.
   real(R8P)                                       :: gpu_memory        !< Available GPU memory.
   real(R8P)                                       :: emin(3), emax(3)  !< Domain dimension.
   real(R8P)                                       :: dxyz(3)           !< Space steps.
   logical                                         :: buf_BOOL          !< Logical buffer.
   integer(I4P)                                    :: buf_I4            !< I4 buffer.
   real(R8P)                                       :: buf_R8            !< R8 buffer.
   character(999)                                  :: buf_CHAR          !< String buffer.

   ! call self%destroy

   print '(A)', 'assign device and get device memory'
   call self%base_gpu%initialize(gpu_memory)

   print '(A)', 'initialize MPI'
   call MPI_COMM_SIZE(MPI_COMM_WORLD, self%procs_number, self%error)
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, self%error)

   print '(A)', 'parse main input'
   call self%parse_input(avail_memory=gpu_memory, filename=filename, nb=nb, nodes_number=nodes_number, nv=nv,      &
                         ni=ni, nj=nj, nk=nk, ngc=ngc, max_level=max_level, bc_type=bc_type, emin=emin, emax=emax, &
                         iu_ref_levels=iu_ref_levels, i_prune=i_prune, j_prune=j_prune, k_prune=k_prune, l_prune=l_prune)
   print '(A)', 'ni, nj, nk, ngc, nv:               '//trim(str([ni, nj, nk, ngc, nv]))
   print '(A)', 'BC types:                          '//trim(str(bc_type))
   print '(A)', 'emin, emax:                        '//trim(str([emin, emax]))
   print '(A)', 'nb, nodes_number:                  '//trim(str([nb, nodes_number]))
   print '(A)', 'max_level, iu_ref_levels:          '//trim(str([max_level, iu_ref_levels]))
   print '(A)', 'i_prune, j_prune, k_prune, l_prune:'//trim(str([i_prune, j_prune, k_prune, l_prune]))

   print '(A)', 'initialize ADAM (tree, field)'
   call self%adam%initialize(ni=ni, nj=nj, nk=nk, ngc=ngc, bc_type=bc_type,                                        &
                             emin=emin, emax=emax, nv=nv, nb=nb, nodes_number=nodes_number, max_level=max_level,   &
                             iu_ref_levels=iu_ref_levels, i_prune=i_prune, j_prune=j_prune, k_prune=k_prune, l_prune=l_prune)

   print '(A)', 'do uniform refine if any'
   call self%adam%refine_uniform(refinement_levels=self%adam%tree%iu_ref_levels, do_blocks_reorder=.false.)

   print '(A)', 'do ijk prune if any'
   call self%adam%prune(ijkl_prune=self%adam%tree%ijkl_prune, do_blocks_reorder=.false.)

   self%ni     =  self%adam%field%grid%ni
   self%nj     =  self%adam%field%grid%nj
   self%nk     =  self%adam%field%grid%nk
   self%ngc    =  self%adam%field%grid%ngc
   self%nb     =  self%adam%field%nb
   self%nv     =  self%adam%field%nv
   self%nv_aux =  9

   print '(A)', 'allocate GPU variables'
   call self%base_gpu%alloc(field=self%adam%field, nv_aux=self%nv_aux)
   self%field         => self%base_gpu%field
   self%grid          => self%base_gpu%field%grid
   self%blocks_number => self%base_gpu%field%blocks_number

   print '(A)', 'parse numerical schemes input setting'
   call self%file_input%get(section_name='schemes', option_name='euler_scheme', val=buf_I4) ; self%euler_scheme = buf_I4
   call self%file_input%get(section_name='schemes', option_name='euler_order' , val=buf_I4) ; self%euler_order  = buf_I4
   self%iweno = (self%euler_order+1)/2
   self%lmax  = (self%euler_order)/2
   call self%file_input%get(section_name='schemes', option_name='visc_scheme' , val=buf_I4) ; self%visc_scheme = buf_I4
   call self%file_input%get(section_name='schemes', option_name='visc_order'  , val=buf_I4) ; self%visc_order  = buf_I4
   call self%fd_initialize
   self%fd_coeff1_gpu = self%fd_coeff1
   self%fd_coeff2_gpu = self%fd_coeff2
   self%fd_conv_gpu   = self%fd_conv

   print '(A)', 'parse physical input setting'
   call self%file_input%get(section_name='physics', option_name='cp',        val=buf_R8) ; self%cp_star   = buf_R8
   call self%file_input%get(section_name='physics', option_name='cv',        val=buf_R8) ; self%cv_star   = buf_R8
   call self%file_input%get(section_name='physics', option_name='mu',        val=buf_R8) ; self%mu_star   = buf_R8
   call self%file_input%get(section_name='physics', option_name='k',         val=buf_R8) ; self%k_star    = buf_R8
   call self%file_input%get(section_name='physics', option_name='dha',       val=buf_R8) ; self%dha_star  = buf_R8
   call self%file_input%get(section_name='physics', option_name='Zeldovich', val=buf_R8) ; self%Zeldovich = buf_R8
   call self%file_input%get(section_name='physics', option_name='Damkohler', val=buf_R8) ; self%Damkohler = buf_R8
   call self%file_input%get(section_name='physics', option_name='Lewis',     val=buf_R8) ; self%Lewis     = buf_R8
   call self%file_input%get(section_name='physics', option_name='visc_law',  val=buf_I4) ; self%visc_law  = buf_I4
   self%gamma_fluid = self%cp_star/self%cv_star
   self%R_star      = self%cp_star-self%cv_star

   print '(A)', 'parse AMR input setting'
   call self%file_input%get(section_name='amr', option_name='frequency', val=buf_I4) ; self%amr_frequency = buf_I4
   call self%file_input%get(section_name='amr', option_name='iters',     val=buf_I4) ; self%amr_iters = buf_I4
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

   print '(A)', 'parse timing input setting'
   call self%file_input%get(section_name="time", option_name="restart",         val=buf_BOOL) ; self%restart          = buf_BOOL
   call self%file_input%get(section_name="time", option_name="restart_basename",val=buf_CHAR) ; self%restart_basename = buf_CHAR
   call self%file_input%get(section_name="time", option_name="restart_save",    val=buf_I4)   ; self%restart_save     = buf_I4
   call self%file_input%get(section_name="time", option_name="time_max",        val=buf_R8)   ; self%time_max         = buf_R8
   call self%file_input%get(section_name="time", option_name="t_max",           val=buf_I4)   ; self%t_max            = buf_I4
   call self%file_input%get(section_name="time", option_name="time_save",       val=buf_R8)   ; self%time_save        = buf_R8
   call self%file_input%get(section_name="time", option_name="n_save",          val=buf_I4)   ; self%n_save           = buf_I4
   call self%file_input%get(section_name="time", option_name="output_basename", val=buf_CHAR) ; self%output_basename  = buf_CHAR
   call self%file_input%get(section_name='time', option_name='CFL',             val=buf_R8)   ; self%CFL              = buf_R8
   call self%file_input%get(section_name="time", option_name="slices_number",   val=buf_I4)   ; self%slices_number    = buf_I4
   call self%runge_kutta_initialize

   if (self%slices_number > 0) then
      allocate(self%slice(self%slices_number))
      print '(A)', 'parse slices input setting'
      do s=1, self%slices_number
         sname = 'slice_'//trim(str(s,.true.))
         call self%file_input%get(section_name=sname, option_name='slice_itype', val=buf_CHAR) ; self%slice(s)%slice_itype=buf_CHAR
         call self%file_input%get(section_name=sname, option_name='slice_save', val=buf_I4)    ; self%slice(s)%slice_save =buf_I4
         call self%file_input%get(section_name=sname, option_name='slice_ni', val=buf_I4)      ; self%slice(s)%slice_nijk(1)=buf_I4
         call self%file_input%get(section_name=sname, option_name='slice_nj', val=buf_I4)      ; self%slice(s)%slice_nijk(2)=buf_I4
         call self%file_input%get(section_name=sname, option_name='slice_nk', val=buf_I4)      ; self%slice(s)%slice_nijk(3)=buf_I4
         call self%file_input%get(section_name=sname, option_name='slice_emin_x', val=buf_R8)  ; self%slice(s)%slice_emin(1)=buf_R8
         call self%file_input%get(section_name=sname, option_name='slice_emin_y', val=buf_R8)  ; self%slice(s)%slice_emin(2)=buf_R8
         call self%file_input%get(section_name=sname, option_name='slice_emin_z', val=buf_R8)  ; self%slice(s)%slice_emin(3)=buf_R8
         call self%file_input%get(section_name=sname, option_name='slice_emax_x', val=buf_R8)  ; self%slice(s)%slice_emax(1)=buf_R8
         call self%file_input%get(section_name=sname, option_name='slice_emax_y', val=buf_R8)  ; self%slice(s)%slice_emax(2)=buf_R8
         call self%file_input%get(section_name=sname, option_name='slice_emax_z', val=buf_R8)  ; self%slice(s)%slice_emax(3)=buf_R8
         allocate(self%slice(s)%slice_points(3,self%slice(s)%slice_nijk(1),self%slice(s)%slice_nijk(2),self%slice(s)%slice_nijk(3)))
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

   print '(A)', 'parse initial conditions input setting'
   call self%file_input%get(section_name="initial_conditions", option_name='ic_type', val=buf_I4) ; self%ic_type = buf_I4
   n_vars = IC_VARS_NUMBER(self%ic_type)
   do i_var=1,n_vars
      oname = "var"//trim(str(i_var,.true.))
      call self%file_input%get(section_name="initial_conditions", option_name=oname, val=buf_R8)
      self%ic_vars(i_var) = buf_R8
   enddo

   print '(A)', 'parse boundary conditions input setting'
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
   self%bc_vars_gpu = self%bc_vars

   print '(A)', 'parse Immersed Boundary input setting'
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
   self%bcs_vars_gpu = self%bcs_vars

   print '(A)', 'allocate large arrays'
   associate(nv=>self%nv, ns=>self%ns, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk, &
             nb=>self%nb, nrk=>self%nrk, nv_aux=>self%nv_aux, n_solids=>self%n_solids, iweno=>self%iweno)
   ! CPU data
   allocate(self%q_aux(1:nv_aux, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nb))
   ! GPU data
   allocate(self%q_gpu(1:nb,        1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nv))
   allocate(self%q_aux_gpu(1:nb,    1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nv_aux))
   allocate(self%fl_gpu(1:nb,       1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nv))
   allocate(self%flx_gpu(1:nb,      1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nv))
   allocate(self%fly_gpu(1:nb,      1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nv))
   allocate(self%flz_gpu(1:nb,      1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nv))
   allocate(self%prhs_gpu(1:nb,     1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nv))
   ! ! debug restart
   ! allocate(self%pbuffer(1:nb,      1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nv))
   ! ! debug restart
   allocate(self%dq_gpu(1:nb,       1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nv))
   allocate(self%q_invert_gpu(1:nb, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nv))
   self%prhs_gpu = 0._R8P
   self%fl_gpu   = 0._R8P
   self%flx_gpu  = 0._R8P
   self%fly_gpu  = 0._R8P
   self%flz_gpu  = 0._R8P
   if (self%n_solids > 0) then
      allocate(self%phi(1:nb,     1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:n_solids))
      allocate(self%phi_gpu(1:nb, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:n_solids))
      self%phi     = -1._R8P
      self%phi_gpu = self%phi
   endif
   allocate(self%gplus_x (nv, 2*iweno, nj, nk, nb))
   allocate(self%gminus_x(nv, 2*iweno, nj, nk, nb))
   allocate(self%gplus_y (nv, 2*iweno, ni, nk, nb))
   allocate(self%gminus_y(nv, 2*iweno, ni, nk, nb))
   allocate(self%gplus_z (nv, 2*iweno, ni, nj, nb))
   allocate(self%gminus_z(nv, 2*iweno, ni, nj, nb))

   allocate(self%dxyz_gpu(1:nb, 1:3)) !TODO
   endassociate

   endsubroutine initialize

   subroutine integrate(self, t, do_ghost_syncro, residual)
   !< Runge Kutta integration of field.
   class(equation_nasto_gpu_object), intent(inout)         :: self             !< The equation.
   real(R8P),                        intent(in)            :: t                !< Time.
   logical,                          intent(in),  optional :: do_ghost_syncro  !< Flag to do syncrous ghost update.
   real(R8P),                        intent(out), optional :: residual         !< Global residual.
   logical                                                 :: do_ghost_syncro_ !< Flag to do syncrous ghost update, local var.
   integer(I4P)                                            :: s                !< Counter.
   integer(I4P)                                            :: i_eikonal        !< Counter.
   integer(I4P), parameter                                 :: n_eikonal=2      !< Counter.
   real(R8P)                                               :: t_s
   real(R8P)                                               :: qnrk
   integer(I4P)                                            :: iermpi

   integer(I4P)                                     :: error                         !< Error trapping flag for CUDAFortran.

   do_ghost_syncro_ = .true. ; if (present(do_ghost_syncro)) do_ghost_syncro_ = do_ghost_syncro
   associate(dt=>self%dt, ni=>self%ni, nj=>self%nj, nk=>self%nk,                                         &
             ngc=>self%ngc, nv=>self%nv, nrk=>self%nrk, ns=>self%ns, blocks_number=>self%blocks_number,  &
             inner_blocks_number=>self%field%inner_blocks_number, n_solids=>self%n_solids, bcs_type=>self%bcs_type(1))

   do s=1, nrk
      call MPI_Barrier(MPI_COMM_WORLD, iermpi)
      t_s = t + dt*(self%ark(s)+self%brk(s))
      call compute_rk_prhs_gpu_cuf(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number,    &
                                   dt=dt, s=s, q_gpu=self%q_gpu, prhs_gpu=self%prhs_gpu,                &
                                   fl_gpu=self%fl_gpu, phi_gpu=self%phi_gpu, qnrk=dt*self%brk(s))

      ! ! debug restart
      ! if (self%itt == 51) then
      ! print*, ' cazzo 1', s
      ! self%pbuffer = self%prhs_gpu
      ! print '(A)', 'debug-restart prhs_gpu[1,ni-3:ni,nj  ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 1)))
      ! print '(A)', 'debug-restart prhs_gpu[2,ni-3:ni,nj  ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 2)))
      ! print '(A)', 'debug-restart prhs_gpu[3,ni-3:ni,nj  ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 3)))
      ! print '(A)', 'debug-restart prhs_gpu[4,ni-3:ni,nj  ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 4)))
      ! print '(A)', 'debug-restart prhs_gpu[5,ni-3:ni,nj  ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 5)))
      ! self%pbuffer = self%fl_gpu
      ! print '(A)', 'debug-restart fl_gpu[1,ni-3:ni,nj    ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 1)))
      ! print '(A)', 'debug-restart fl_gpu[2,ni-3:ni,nj    ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 2)))
      ! print '(A)', 'debug-restart fl_gpu[3,ni-3:ni,nj    ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 3)))
      ! print '(A)', 'debug-restart fl_gpu[4,ni-3:ni,nj    ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 4)))
      ! print '(A)', 'debug-restart fl_gpu[5,ni-3:ni,nj    ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 5)))
      ! endif
      ! ! debug restart

      call self%update_ghost_gpu(q_gpu=self%q_gpu)

      ! ! debug restart
      ! if (self%itt == 51) then
      ! print*, ' cazzo 2', s
      ! self%pbuffer = self%q_gpu
      ! print '(A)', 'debug-restart q_gpu[1,ni-3:ni,0     ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 0, 4, 1)))
      ! print '(A)', 'debug-restart q_gpu[1,ni-3:ni,nj    ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 1)))
      ! print '(A)', 'debug-restart q_gpu[2,ni-3:ni,0     ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 0, 4, 2)))
      ! print '(A)', 'debug-restart q_gpu[2,ni-3:ni,nj    ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 2)))
      ! print '(A)', 'debug-restart q_gpu[3,ni-3:ni,0     ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 0, 4, 3)))
      ! print '(A)', 'debug-restart q_gpu[3,ni-3:ni,nj    ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 3)))
      ! print '(A)', 'debug-restart q_gpu[4,ni-3:ni,0     ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 0, 4, 4)))
      ! print '(A)', 'debug-restart q_gpu[4,ni-3:ni,nj    ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 4)))
      ! print '(A)', 'debug-restart q_gpu[5,ni-3:ni,0     ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 0, 4, 5)))
      ! print '(A)', 'debug-restart q_gpu[5,ni-3:ni,nj    ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 5)))
      ! endif
      ! ! debug restart

      ! ------------------------------------------------------------------------------------------------
      ! ANDREA IB
      ! ------------------------------------------------------------------------------------------------
      if (n_solids > 0) then
         do i_eikonal=1,n_eikonal
            call MPI_Barrier(MPI_COMM_WORLD, iermpi)

            error = cudaGetLastError()
            if (error /= cudaSuccess) then
               print*,'BEFORE EIK POST FRA CUDA ERROR ',cudaGetErrorString(error)
               call MPI_Abort(MPI_COMM_WORLD, -15,error) ; STOP
            endif
            call evolve_eikonal_q_gpu_cuf(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, &
                                          phi_gpu=self%phi_gpu,                                             &
                                          dx_gpu=self%dxyz_gpu(:,1),                                        &
                                          dy_gpu=self%dxyz_gpu(:,2),                                        &
                                          dz_gpu=self%dxyz_gpu(:,3),                                        &
                                          dq_gpu=self%dq_gpu,                                               &
                                          q_gpu=self%q_gpu)

      ! ! debug restart
      ! if (self%itt == 51) then
      ! print*, ' cazzo 3', s
      ! self%pbuffer = self%q_gpu
      ! print '(A)', 'debug-restart q_gpu[1,ni-3:ni,0     ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 0, 4, 1)))
      ! print '(A)', 'debug-restart q_gpu[1,ni-3:ni,nj    ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 1)))
      ! print '(A)', 'debug-restart q_gpu[2,ni-3:ni,0     ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 0, 4, 2)))
      ! print '(A)', 'debug-restart q_gpu[2,ni-3:ni,nj    ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 2)))
      ! print '(A)', 'debug-restart q_gpu[3,ni-3:ni,0     ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 0, 4, 3)))
      ! print '(A)', 'debug-restart q_gpu[3,ni-3:ni,nj    ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 3)))
      ! print '(A)', 'debug-restart q_gpu[4,ni-3:ni,0     ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 0, 4, 4)))
      ! print '(A)', 'debug-restart q_gpu[4,ni-3:ni,nj    ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 4)))
      ! print '(A)', 'debug-restart q_gpu[5,ni-3:ni,0     ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 0, 4, 5)))
      ! print '(A)', 'debug-restart q_gpu[5,ni-3:ni,nj    ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 5)))
      ! endif
      ! ! debug restart

            error = cudaGetLastError()
            if (error /= cudaSuccess) then
               print*,'AFTER EIK POST FRA CUDA ERROR ',cudaGetErrorString(error)
               call MPI_Abort(MPI_COMM_WORLD, -15,error) ; STOP
            endif
            call self%update_ghost_gpu(q_gpu=self%q_gpu)
         enddo
         call invert_eikonal_field(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number,        &
                                   q_gpu=self%q_gpu(:,:,:,:,:), q_invert_gpu=self%q_invert_gpu(:,:,:,:,:),  &
                                   phi_gpu=self%phi_gpu, bcs_type=bcs_type)
      else
         ! added for restart debug...
         self%q_invert_gpu = self%q_gpu
      endif

      ! ! debug restart
      ! if (self%itt == 51) then
      ! print*, ' cazzo 4', s
      ! self%pbuffer = self%q_invert_gpu
      ! print '(A)', 'debug-restart q_invert_gpu[1,ni-3:ni,0     ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 0, 4, 1)))
      ! print '(A)', 'debug-restart q_invert_gpu[1,ni-3:ni,nj    ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 1)))
      ! print '(A)', 'debug-restart q_invert_gpu[2,ni-3:ni,0     ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 0, 4, 2)))
      ! print '(A)', 'debug-restart q_invert_gpu[2,ni-3:ni,nj    ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 2)))
      ! print '(A)', 'debug-restart q_invert_gpu[3,ni-3:ni,0     ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 0, 4, 3)))
      ! print '(A)', 'debug-restart q_invert_gpu[3,ni-3:ni,nj    ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 3)))
      ! print '(A)', 'debug-restart q_invert_gpu[4,ni-3:ni,0     ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 0, 4, 4)))
      ! print '(A)', 'debug-restart q_invert_gpu[4,ni-3:ni,nj    ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 4)))
      ! print '(A)', 'debug-restart q_invert_gpu[5,ni-3:ni,0     ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 0, 4, 5)))
      ! print '(A)', 'debug-restart q_invert_gpu[5,ni-3:ni,nj    ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 5)))
      ! endif
      ! ! debug restart

      call self%compute_aux(q_gpu=self%q_invert_gpu, q_aux_gpu=self%q_aux_gpu)

      ! ! debug restart
      ! if (self%itt == 51) then
      ! print*, ' cazzo 4-bis', s
      ! self%pbuffer = self%q_aux_gpu(:,:,:,:,1:5)
      ! print '(A)', 'debug-restart q_aux_gpu[1,ni-3:ni,0   ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 0, 4, 1)))
      ! print '(A)', 'debug-restart q_aux_gpu[1,ni-3:ni,nj  ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 1)))
      ! print '(A)', 'debug-restart q_aux_gpu[2,ni-3:ni,0   ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 0, 4, 2)))
      ! print '(A)', 'debug-restart q_aux_gpu[2,ni-3:ni,nj  ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 2)))
      ! print '(A)', 'debug-restart q_aux_gpu[3,ni-3:ni,0   ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 0, 4, 3)))
      ! print '(A)', 'debug-restart q_aux_gpu[3,ni-3:ni,nj  ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 3)))
      ! print '(A)', 'debug-restart q_aux_gpu[4,ni-3:ni,0   ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 0, 4, 4)))
      ! print '(A)', 'debug-restart q_aux_gpu[4,ni-3:ni,nj  ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 4)))
      ! print '(A)', 'debug-restart q_aux_gpu[5,ni-3:ni,0   ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 0, 4, 5)))
      ! print '(A)', 'debug-restart q_aux_gpu[5,ni-3:ni,nj  ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 5)))
      ! self%pbuffer = self%q_aux_gpu(:,:,:,:,4:8)
      ! print '(A)', 'debug-restart q_aux_gpu[7,ni-3:ni,0   ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 0, 4, 3)))
      ! print '(A)', 'debug-restart q_aux_gpu[7,ni-3:ni,nj  ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 3)))
      ! print '(A)', 'debug-restart q_aux_gpu[8,ni-3:ni,0   ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 0, 4, 4)))
      ! print '(A)', 'debug-restart q_aux_gpu[8,ni-3:ni,nj  ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 4)))
      ! print '(A)', 'debug-restart q_aux_gpu[5,ni-3:ni,0   ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 0, 4, 5)))
      ! print '(A)', 'debug-restart q_aux_gpu[5,ni-3:ni,nj  ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 5)))
      ! endif
      ! ! debug restart

      call MPI_Barrier(MPI_COMM_WORLD, iermpi)
      call self%compute_residuals_gpu(ni=ni, nj=nj, nk=nk, ngc=ngc, ns=ns, blocks_number=blocks_number,                      &
                                      dx_gpu = self%dxyz_gpu(:,1), dy_gpu = self%dxyz_gpu(:,2), dz_gpu = self%dxyz_gpu(:,3), &
                                      q_aux_gpu = self%q_aux_gpu,  phi_gpu = self%phi_gpu,      fl_gpu = self%fl_gpu,        &
                                      flx_gpu   = self%flx_gpu,    fly_gpu = self%fly_gpu,      flz_gpu = self%flz_gpu,      &
                                      fd_conv_gpu = self%fd_conv_gpu, fd_coeff1_gpu = self%fd_coeff1_gpu,                    &
                                      fd_coeff2_gpu = self%fd_coeff2_gpu,                                                    &
                                      gminus_x = self%gminus_x, gminus_y = self%gminus_y, gminus_z = self%gminus_z,          &
                                      gplus_x = self%gplus_x,   gplus_y  = self%gplus_y,  gplus_z  = self%gplus_z,           &
                                      euler_scheme = self%euler_scheme, visc_scheme = self%visc_scheme,                      &
                                      lmax = self%lmax, iweno = self%iweno, visc_order = self%visc_order,                    &
                                      visc_law = self%visc_law,                                                              &
                                      cp_star       = self%cp_star,  cv_star = self%cv_star, gamma_fluid = self%gamma_fluid, &
                                      R_star        = self%R_star,   mu_star = self%mu_star, k_star      = self%k_star,      &
                                      dha_star      = self%dha_star, Lewis   = self%Lewis,   Zeldovich   = self%Zeldovich,   &
                                      Damkohler     = self%Damkohler)
      !call self%set_bc_rhs(q_gpu=self%fl_gpu(:,:,:,:,:), q_aux_gpu=self%q_aux_gpu)

      ! ! debug restart
      ! if (self%itt == 51) then
      ! print*, ' cazzo 5', s
      ! self%pbuffer = self%fl_gpu
      ! print '(A)', 'debug-restart fl_gpu[rho*v,ni-3:ni,0     ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 0, 4, 3)))
      ! print '(A)', 'debug-restart fl_gpu[rho*v,ni-3:ni,nj    ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 3)))
      ! endif
      ! ! debug restart

      call compute_rk_linear_gpu_cuf(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number,      &
                                     dt=dt, q_gpu=self%q_gpu, prhs_gpu=self%prhs_gpu,                       &
                                     fl_gpu=self%fl_gpu, phi_gpu=self%phi_gpu, qnrk=dt*self%ark(s))

      ! ! debug restart
      ! if (self%itt == 51) then
      ! print*, ' cazzo 6', s
      ! self%pbuffer = self%q_gpu
      ! print '(A)', 'debug-restart q_gpu[rho*v,ni-3:ni,0   ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 0, 4, 3)))
      ! print '(A)', 'debug-restart q_gpu[rho*v,ni-3:ni,nj  ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 3)))
      ! endif
      ! ! debug restart

   enddo
   endassociate
   endsubroutine integrate

   subroutine load_restart_files(self, t, time)
   !< Save restart files.
   class(equation_nasto_gpu_object), intent(inout) :: self !< The equation.
   integer(I4P),                     intent(out)   :: t    !< Time iteration.
   real(R8P),                        intent(out)   :: time !< Time.

   call self%adam%load_restart_files(basename=self%restart_basename, t=t, time=time)
   call self%adam%make_comm_local_maps_ghost_bc
   call self%copy_cpu_gpu
   ! note: update ghost has been commented because it is not necessary
   ! call self%update_ghost_gpu(q_gpu=self%q_gpu)
   endsubroutine load_restart_files

   subroutine mark_by_geo(self, delta_fine, delta_coarse, threshold, do_init)
   !< Mark blocks to be refined/derefined by a `grad(rho)` value.
   class(equation_nasto_gpu_object), intent(inout)        :: self           !< The equation.
   real(R8P),                        intent(in)           :: delta_fine     !< Maximum cell delta in fine grids.
   real(R8P),                        intent(in)           :: delta_coarse   !< Minimum cell delta in coarse grids.
   real(R8P),                        intent(in), optional :: threshold      !< Threshold for sphere proximity.
   real(R8P)                                              :: threshold_     !< Threshold for sphere proximity, local var.
   real(R8P)                                              :: max_cell_delta !< Maximum cell delta.
   real(R8P)                                              :: distance       !< Value (max) of gradient of rho.
   integer(I4P)                                           :: b              !< Counter.
   logical, optional                , intent(in)          :: do_init
   logical                                                :: do_init_

   do_init_ = .true.    ; if (present(do_init)) do_init_ = do_init
   threshold_ = 2.2_R8P ; if (present(threshold)) threshold_ = threshold
   if(do_init_) then
       self%field%refinements_needed = [(TO_BE_DEREFINED,b=1,self%blocks_number)]

       !if (allocated(self%field%refinements_needed)) deallocate(self%field%refinements_needed)
       !allocate(self%field%refinements_needed(self%blocks_number))
       !self%field%refinements_needed(:) = TO_BE_DEREFINED
   endif
   associate (ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, &
              blocks_number=>self%blocks_number, ns=>self%ns, dxyz=>self%field%dxyz, phi=>self%phi)
      do b=1, blocks_number
         distance = 1._R8P
         if(maxval(phi(b,:,:,:,1))*minval(phi(b,:,:,:,1)) < 0._R8P) then
            distance = 0._R8P
         endif
         max_cell_delta = max_cell_delta_dist(distance=distance)

         if (maxval(dxyz(:,b)) > max_cell_delta) then
            self%field%refinements_needed(b) = TO_BE_REFINED
         elseif (maxval(dxyz(:,b)) * threshold_ < max_cell_delta) then
            self%field%refinements_needed(b) = max(self%field%refinements_needed(b), TO_BE_DEREFINED)
         else
            self%field%refinements_needed(b) = max(self%field%refinements_needed(b), TO_NOT_TOUCH)
         endif
      enddo
   endassociate
   contains
      function max_cell_delta_dist(distance) result(delta)
      !< Return the maximum cell delta given a comparison distance.
      real(R8P),          intent(in) :: distance !< Comparison distance.
      real(R8P)                      :: delta    !< Maximum cell delta admissible.

      if (abs(distance) < epsilon(0._R8P)) then
         delta = delta_fine
      else
         delta = delta_coarse
         !delta = huge(0._R8P)
      endif
      endfunction max_cell_delta_dist
   endsubroutine mark_by_geo

   subroutine mark_by_grad_var(self, grad_tol, delta_fine, delta_coarse, ivar, threshold, do_init)
   !< Mark blocks to be refined/derefined by a `grad(rho)` value.
   class(equation_nasto_gpu_object), intent(inout)        :: self           !< The equation.
   real(R8P),                        intent(in)           :: grad_tol       !< Gradiend tolerance value.
   real(R8P),                        intent(in)           :: delta_fine     !< Maximum cell delta in fine grids.
   real(R8P),                        intent(in)           :: delta_coarse   !< Minimum cell delta in coarse grids.
   integer(I4P),                     intent(in), optional :: ivar            !< Variable for marking.
   integer(I4P)                                           :: ivar_           !< Variable for marking (local var).
   real(R8P),                        intent(in), optional :: threshold      !< Threshold for sphere proximity.
   real(R8P)                                              :: threshold_     !< Threshold for sphere proximity, local var.
   real(R8P)                                              :: max_cell_delta !< Maximum cell delta.
   real(R8P)                                              :: grad_var       !< Value (max) of gradient of var.
   integer(I4P)                                           :: b              !< Counter.
   logical, optional                , intent(in)          :: do_init
   logical                                                :: do_init_

   ivar_     = IRHO     ; if (present(ivar)) ivar_ = ivar
   do_init_ = .true.    ; if (present(do_init)) do_init_ = do_init
   threshold_ = 2.2_R8P ; if (present(threshold)) threshold_ = threshold
   if(do_init_) self%field%refinements_needed = [(TO_BE_DEREFINED,b=1,self%blocks_number)]
   associate (ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, &
              blocks_number=>self%blocks_number, ns=>self%ns, dxyz=>self%field%dxyz)
      call self%compute_aux(q_gpu=self%q_gpu, q_aux_gpu=self%q_aux_gpu)
      do b=1, blocks_number
         grad_var = gradient_cuf(b=b, ni=ni, nj=nj, nk=nk, ngc=ngc, &
                                 dx=dxyz(1,b), dy=dxyz(2,b), dz=dxyz(3,b), q_gpu=self%q_aux_gpu, ivar=ivar_)

         max_cell_delta = max_cell_delta_grad(grad=grad_var)
         !max_cell_delta = 0.1

         if (maxval(dxyz(:,b)) > max_cell_delta) then
            self%field%refinements_needed(b) = TO_BE_REFINED
         elseif (maxval(dxyz(:,b)) * threshold_ < max_cell_delta) then
            self%field%refinements_needed(b) = max(self%field%refinements_needed(b), TO_BE_DEREFINED)
         else
            self%field%refinements_needed(b) = max(self%field%refinements_needed(b), TO_NOT_TOUCH)
         endif
      enddo
   endassociate
   contains
      function gradient_cuf(b, ni, nj, nk, ngc, dx, dy, dz, q_gpu, ivar) result(gradient)
      !< Gradient done by CUF threads.
      integer(I4P), intent(in)         :: b                                 !< Block index.
      integer(I4P), intent(in)         :: ni                                !< Grid cells number in I direction.
      integer(I4P), intent(in)         :: nj                                !< Grid cells number in J direction.
      integer(I4P), intent(in)         :: nk                                !< Grid cells number in K direction.
      integer(I4P), intent(in)         :: ngc                               !< Ghost cells number.
      real(R8P),    intent(in)         :: dx                                !< X space step.
      real(R8P),    intent(in)         :: dy                                !< Y space step.
      real(R8P),    intent(in)         :: dz                                !< Z space step.
      real(R8P),    intent(in), device :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field component to which apply gradient.
      integer(I4P), intent(in)         :: ivar                              !< Ghost cells number.
      real(R8P)                        :: gradient                          !< Maximum gradient of q.
      real(R8P)                        :: grad                              !< Current gradient of q.
      integer(I4P)                     :: i, j, k                           !< Counter.
      integer(I4P)                     :: iercuda                           !< Error trapping flag for CUDAFortran.
      real(R8P), parameter             :: tol=1.e-12                        !< Gradient denominator tolerance.

      gradient = 0._R8P
      !$cuf kernel do(3) <<<*,*>>>
      do k=1, nk
         do j=1, nj
            do i=1, ni
               grad = sqrt(((q_gpu(b,i+1,j,k,ivar) - q_gpu(b,i-1,j,k,ivar))/(2*dx))**2 + &
                           ((q_gpu(b,i,j+1,k,ivar) - q_gpu(b,i,j-1,k,ivar))/(2*dy))**2 + &
                           ((q_gpu(b,i,j,k+1,ivar) - q_gpu(b,i,j,k-1,ivar))/(2*dz))**2)
               grad = grad/(abs(q_gpu(b,i,j,k,ivar))+tol)
               gradient = max(gradient, grad)
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
      endfunction gradient_cuf

      function max_cell_delta_grad(grad) result(delta)
      !< Return the maximum cell delta given a gradient tollerance.
      real(R8P), intent(in) :: grad  !< Gradient value.
      real(R8P)             :: delta !< Maximum cell delta admissible.

      if (grad > grad_tol) then
         delta = delta_fine
      else
         delta = delta_coarse
      endif
      endfunction max_cell_delta_grad
   endsubroutine mark_by_grad_var

   subroutine move_phi(self, velocity)
   !< Move phi and the actual ptree representation.
   class(equation_nasto_gpu_object), intent(inout) :: self        !< The equation.
   real(R8P),                        intent(in)    :: velocity(3) !< Velocity of the movement.

   call move_phi_cuf(ni            = self%ni,            &
                     nj            = self%nj,            &
                     nk            = self%nk,            &
                     ngc           = self%ngc,           &
                     blocks_number = self%blocks_number, &
                     velocity      = velocity,           &
                     phi_gpu       = self%phi_gpu,       &
                     dphi_gpu      = self%dq_gpu)
   endsubroutine move_phi

   subroutine move_phi_cuf(ni, nj, nk, ngc, blocks_number, velocity, phi_gpu, dphi_gpu)
   !< Move phi and the actual ptree representation.
   integer(I4P), intent(in)            :: ni                                   !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj                                   !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk                                   !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: ngc                                  !< Ghost grid number.
   integer(I4P), intent(in)            :: blocks_number                        !< Number of blocks.
   real(R8P),    intent(in)            :: velocity(3)                          !< Velocity of the movement.
   real(R8P),    intent(inout), device ::  phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Distance function.
   real(R8P),    intent(inout), device :: dphi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Distance function gradient.
   real(R8P)                           :: n_phi_x, n_phi_y, n_phi_z, n_phi     !< Eikonal direction.
   integer(I4P)                        :: b, i, j, k, v                        !< Counter.
   integer(I4P)                        :: iercuda                              !< Error trapping flag for CUDAFortran.

   n_phi_x = velocity(1)
   n_phi_y = velocity(2)
   n_phi_z = velocity(3)
   n_phi = abs(n_phi_x) + abs(n_phi_y) + abs(n_phi_z) + 10e-12
   n_phi = 0.9_R8P / n_phi
   n_phi_x = n_phi_x * n_phi
   n_phi_y = n_phi_y * n_phi
   n_phi_z = n_phi_z * n_phi

   !$cuf kernel do(4) <<<*,*>>>
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      do v=1, 1
         dphi_gpu(b,i,j,k,v) = 0._R8P
      enddo
      if (n_phi_x > 0._R8P) then
         do v=1, 1
            dphi_gpu(b,i,j,k,v) = dphi_gpu(b,i,j,k,v) + abs(n_phi_x) * (phi_gpu(b,i,j,k,v) - phi_gpu(b,i-1,j,k,v))
         enddo
      else
         do v=1, 1
            dphi_gpu(b,i,j,k,v) = dphi_gpu(b,i,j,k,v) + abs(n_phi_x) * (phi_gpu(b,i,j,k,v) - phi_gpu(b,i+1,j,k,v))
         enddo
      endif
      if (n_phi_y > 0._R8P) then
         do v=1, 1
            dphi_gpu(b,i,j,k,v) = dphi_gpu(b,i,j,k,v) + abs(n_phi_y) * (phi_gpu(b,i,j,k,v) - phi_gpu(b,i,j-1,k,v))
         enddo
      else
         do v=1, 1
            dphi_gpu(b,i,j,k,v) = dphi_gpu(b,i,j,k,v) + abs(n_phi_y) * (phi_gpu(b,i,j,k,v) - phi_gpu(b,i,j+1,k,v))
         enddo
      endif
      if (n_phi_z > 0._R8P) then
         do v=1, 1
            dphi_gpu(b,i,j,k,v) = dphi_gpu(b,i,j,k,v) + abs(n_phi_z) * (phi_gpu(b,i,j,k,v) - phi_gpu(b,i,j,k-1,v))
         enddo
      else
         do v=1, 1
            dphi_gpu(b,i,j,k,v) = dphi_gpu(b,i,j,k,v) + abs(n_phi_z) * (phi_gpu(b,i,j,k,v) - phi_gpu(b,i,j,k+1,v))
         enddo
      endif
   enddo
   enddo
   enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()

   !$cuf kernel do(4) <<<*,*>>>
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      do v=1, 1
         phi_gpu(b,i,j,k,v) = phi_gpu(b,i,j,k,v) - dphi_gpu(b,i,j,k,v)
      enddo
   enddo
   enddo
   enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine move_phi_cuf

   subroutine parse_input(self, avail_memory, filename, nb, nodes_number, nv, ni, nj, nk, ngc, &
                          max_level, emin, emax, bc_type, iu_ref_levels, i_prune, j_prune, k_prune, l_prune)
   class(equation_nasto_gpu_object), intent(inout) :: self             !< The equation.
   real(R8P)                       , intent(in)    :: avail_memory     !< Available memory
   character(*)                    , intent(in)    :: filename         !< Input file name.
   integer(I4P)                    , intent(out)   :: nv               !< Number of variables.
   integer(I4P)                    , intent(out)   :: nb               !< Number of allocated blocks.
   integer(I8P)                    , intent(out)   :: nodes_number     !< Number of tree nodes.
   integer(I4P)                    , intent(out)   :: ni, nj, nk, ngc  !< Grid dimensions.
   integer(I4P)                    , intent(out)   :: max_level        !< Max tree level.
   integer(I4P)                    , intent(out)   :: bc_type(6)       !< Boundary conditions flags.
   real(R8P)                       , intent(out)   :: emin(3), emax(3) !< Domain sizes.
   integer(I4P)                    , intent(out)   :: iu_ref_levels    !< Domain sizes.
   integer(I4P)                    , intent(out)   :: i_prune          !< Pruning along x.
   integer(I4P)                    , intent(out)   :: j_prune          !< Pruning along y.
   integer(I4P)                    , intent(out)   :: k_prune          !< Pruning along z.
   integer(I4P)                    , intent(out)   :: l_prune          !< Pruning level.
   integer(I4P)                                    :: ns               !< Number of species and variables.
   real(R8P)                                       :: block_weight     !< Number of points per block.
   real(R8P), parameter                            :: save_factor=0.8  !< Factor to avoid GPU completely full.
   integer(I4P)                                    :: buf_I4           !< Integer buffer.
   real(R8P)                                       :: buf_R8           !< Real buffer.
   integer(I4P)                                    :: size_of_real     !< Size of real.

   call self%file_input%initialize(filename=trim(filename))
   call self%file_input%load
   call self%file_input%get(section_name='grid', option_name='ni', val=ni)
   call self%file_input%get(section_name='grid', option_name='nj', val=nj)
   call self%file_input%get(section_name='grid', option_name='nk', val=nk)
   call self%file_input%get(section_name='grid', option_name='ngc', val=ngc)
   call self%file_input%get(section_name='physics', option_name='ns', val=ns)

   nv           = ns + 4
   block_weight = (ngc+ni+ngc) * (ngc+nj+ngc) * (ngc+nk+ngc) * nv
   size_of_real = storage_size(emin(1))/8.
   nb           = nint(save_factor*avail_memory*1e9/(self%field_gpu_number*block_weight*size_of_real))
   nodes_number = nb*self%procs_number
   print '(A)', 'available places for blocks [nb]: '//trim(str(nb))

   call self%file_input%get(section_name='bc_x_min', option_name='type', val=buf_I4) ; bc_type(1) = buf_I4
   call self%file_input%get(section_name='bc_x_max', option_name='type', val=buf_I4) ; bc_type(2) = buf_I4
   call self%file_input%get(section_name='bc_y_min', option_name='type', val=buf_I4) ; bc_type(3) = buf_I4
   call self%file_input%get(section_name='bc_y_max', option_name='type', val=buf_I4) ; bc_type(4) = buf_I4
   call self%file_input%get(section_name='bc_z_min', option_name='type', val=buf_I4) ; bc_type(5) = buf_I4
   call self%file_input%get(section_name='bc_z_max', option_name='type', val=buf_I4) ; bc_type(6) = buf_I4

   call self%file_input%get(section_name='grid', option_name='emin_x', val=buf_R8) ; emin(1) = buf_R8
   call self%file_input%get(section_name='grid', option_name='emin_y', val=buf_R8) ; emin(2) = buf_R8
   call self%file_input%get(section_name='grid', option_name='emin_z', val=buf_R8) ; emin(3) = buf_R8
   call self%file_input%get(section_name='grid', option_name='emax_x', val=buf_R8) ; emax(1) = buf_R8
   call self%file_input%get(section_name='grid', option_name='emax_y', val=buf_R8) ; emax(2) = buf_R8
   call self%file_input%get(section_name='grid', option_name='emax_z', val=buf_R8) ; emax(3) = buf_R8
   call self%file_input%get(section_name='amr', option_name='max_level', val=max_level)
   call self%file_input%get(section_name='amr', option_name='iu_ref_levels', val=iu_ref_levels)
   call self%file_input%get(section_name='amr', option_name='i_prune', val=i_prune)
   call self%file_input%get(section_name='amr', option_name='j_prune', val=j_prune)
   call self%file_input%get(section_name='amr', option_name='k_prune', val=k_prune)
   call self%file_input%get(section_name='amr', option_name='l_prune', val=l_prune)

   !TODO
   ! Check all input parameters

   endsubroutine parse_input

   subroutine print_progress(self, t, time, t_max, time_max)
   !< Print simulation progress.
   class(equation_nasto_gpu_object), intent(in) :: self     !< The equation.
   integer(I4P),                     intent(in) :: t        !< Time iteration.
   real(R8P),                        intent(in) :: time     !< Time.
   integer(I4P),                     intent(in) :: t_max    !< Maximum time iteration.
   real(R8P),                        intent(in) :: time_max !< Maximum time of integration.

   print '(A)', ''
   print '(A)', 'Iteration[rank] = '//trim(str(t))//'['//trim(str(self%myrank))//']'
   if (self%myrank==0) then
         print '(A)', 't:             '//trim(str(t,.true.))
         print '(A)', 'blocks number: '//trim(str(self%adam%tree%nodes_number, .true.))
         print '(A)', 'time step:     '//trim(str(self%dt, .true.))
         print '(A)', 'time:          '//trim(str(time, .true.))
      if (t_max <= 0) then
         print '(A)', 'progress:      '//trim(str(int(time/time_max * 100), .true.))//'%'
      else
         print '(A)', 'progress:      '//trim(str(int((t*1._R8P)/t_max * 100), .true.))//'%'
      endif
   endif
   print '(A)', ''
   endsubroutine print_progress

   subroutine refine_uniform(self, refinement_levels)
   !< Refine all blocks uniformly.
   class(equation_nasto_gpu_object), intent(inout) :: self              !< The equation.
   integer(I4P),                     intent(in)    :: refinement_levels !< Number of refinement to be performed.
   integer(I4P)                                    :: l                 !< Counter.

   do l=1, refinement_levels
      call self%adam%tree%mark_all_nodes(mark=TO_BE_REFINED)
      call self%adam%amr_update(do_blocks_reorder=.false., do_mpi_redistribute=.true.)
   enddo
   endsubroutine

   subroutine run(self, filename)
   class(equation_nasto_gpu_object), intent(inout) :: self             !< The equation.
   character(*),                     intent(in)    :: filename         !< Input file name.
   real(R8P)                                       :: timing(1:2)      !< Tic toc timing.
   real(R8P)                                       :: timing_step(1:2) !< Tic toc timing.

   call MPI_INIT(self%error)

   call self%initialize(filename=filename)
   if (self%restart) then
      print '(A)', 'restart simulation from "'//trim(self%restart_basename)//'" files'
      call self%load_restart_files(t=self%it, time=self%time)
      print '(A)', 'restart [t, time]: '//trim(str(self%it))//', '//trim(str(self%time))
   else
      call self%set_initial_conditions()
      self%time = 0._R8P
      self%it = 0
   endif
   if (self%n_solids > 0) call self%update_phi()

   call self%save_simulation_data
   call MPI_BARRIER(MPI_COMM_WORLD, self%error) ; timing(1) = MPI_Wtime()

   integration: do
      call MPI_BARRIER(MPI_COMM_WORLD, self%error) ; timing_step(1) = MPI_Wtime()
      self%it = self%it + 1

      if(mod(self%it,self%amr_frequency) == 0) then
         call self%amr_update()
         call MPI_BARRIER(MPI_COMM_WORLD, self%error) ; timing_step(2) = MPI_Wtime()
         print '(A, F18.10)', 'step timing (AMR): ', timing_step(2) - timing_step(1)
      endif

      call self%compute_dt()
      if ((self%t_max <= 0).and.(self%time + self%dt > self%time_max)) self%dt = self%time_max - self%time

      call self%integrate(t=self%time)

      self%time = self%time + self%dt
      call self%print_progress(t=self%it, time=self%time, t_max=self%t_max, time_max=self%time_max)

      call self%save_simulation_data
      call MPI_BARRIER(MPI_COMM_WORLD, self%error) ; timing_step(2) = MPI_Wtime()
      print '(A, F18.10)', 'step timing (save data): ', timing_step(2) - timing_step(1)

      if (((self%t_max <= 0).and.(self%time >= self%time_max)).or.((self%it>=self%t_max).and.(self%t_max > 0))) exit integration

      call MPI_BARRIER(MPI_COMM_WORLD, self%error) ; timing_step(2) = MPI_Wtime()
      print '(A, F18.10)', 'step timing: ', timing_step(2) - timing_step(1)
   enddo integration

   call MPI_BARRIER(MPI_COMM_WORLD, self%error) ; timing(2) = MPI_Wtime()
   print '(A, F18.10)', 'averaged timing: ', (timing(2) - timing(1))/self%it

   call self%save_simulation_data
   call self%adam%finalize
   endsubroutine run

   subroutine runge_kutta_initialize(self)
   !< Initialize Runge-Kutta data.
   class(equation_nasto_gpu_object), intent(inout) :: self !< The equation.

   call self%file_input%get(section_name='time', option_name='nrk', val=self%nrk)
   allocate(self%ark(self%nrk), self%brk(self%nrk))
   select case(self%nrk)
   case(3_I4P)
      self%ark(:)=[8._R8P  /15._R8P, 5._R8P  /12._R8P, 3._R8P  /4._R8P]
      self%brk(:)=[0._R8P, -17._R8P/60._R8P , -5._R8P /12._R8P]
   case(4_I4P)
      self%ark(:) = [8._R8P/17._R8P,17._R8P /60._R8P,5._R8P /12._R8P,3._R8P/4._R8P]
      self%brk(:) = [0._R8P,-15._R8P/68._R8P,-17._R8P/60._R8P,-5._R8P/12._R8P]
   endselect
   endsubroutine runge_kutta_initialize

   subroutine save_simulation_data(self)
   !< Save all simulation data.
   class(equation_nasto_gpu_object), intent(inout) :: self                     !< The equation.
   logical                                         :: is_update_ghost_gpu_done !< Flag to minimize ghosts-update-calls for IO.
   logical                                         :: is_alo_slice_to_save     !< Flag to check if at least one slice is to save.
   integer(I4P)                                    :: s                        !< Slices counter.

   ! update ghost cells if necessary
   is_update_ghost_gpu_done = .false.
   is_alo_slice_to_save     = .false.
   if (self%slices_number>0) then
      do s=1, self%slices_number
         if (mod(self%it,self%slice(s)%slice_save)==0.or.self%it==self%t_max.or.&
            (((self%t_max <= 0).and.(self%time >= self%time_max)).or.((self%it>=self%t_max).and.(self%t_max > 0)))) then
            is_alo_slice_to_save = .true.
            if (.not.is_update_ghost_gpu_done) then
               is_update_ghost_gpu_done = .true.
               call self%update_ghost_gpu(q_gpu=self%q_gpu)
               call self%update_ghost_gpu(q_gpu=self%q_aux_gpu)
            endif
         endif
      enddo
   endif
   ! copy GPU data to CPU
   if (mod(self%it,self%n_save)==0.or.self%it==self%t_max.or.& ! HDF5 output
      (mod(self%it,self%restart_save)==0).or.                & ! restart output
      is_alo_slice_to_save.or.                               & ! slices output
      (((self%t_max <= 0).and.(self%time >= self%time_max)).or.((self%it>=self%t_max).and.(self%t_max > 0)))) then
      call self%copy_gpu_cpu(compute_q_aux=.true.)
   endif
   ! save data
   call self%save_hdf5
   call self%save_restart_files
   call self%save_slices
   endsubroutine save_simulation_data

   subroutine save_hdf5(self, output_basename)
   !< Save simulation data in HDF5 format.
   class(equation_nasto_gpu_object), intent(inout)        :: self             !< The equation.
   character(*),                     intent(in), optional :: output_basename  !< Output basename.
   character(:), allocatable                              :: output_basename_ !< Output basename, local var.

   if (mod(self%it,self%n_save)==0.or.self%it==self%t_max.or.&
      (((self%t_max <= 0).and.(self%time >= self%time_max)).or.((self%it>=self%t_max).and.(self%t_max > 0)))) then
      print '(A)', 'save HDF5 files t: '//trim(str(self%it,.true.))//', time: '//trim(str(self%time,.true.))
      output_basename_ = trim(self%output_basename)//'-'//trim(strz(self%it,9))
      if (present(output_basename)) output_basename_ = trim(output_basename)
      call self%adam%save_hdf5(basename=trim(output_basename_),                                  &
                               q=self%field%q,                                                   &
                               q_aux=self%q_aux,                                                 &
                               q_name=['rho','rhu','rhv','rhw','rhe'],                           &
                               q_aux_name=['rhob','u','v','w','ya','tem','pres','ental','csp'],  &
                               with_cell_morton=.true.)
   endif
   endsubroutine save_hdf5

   subroutine save_restart_files(self)
   !< Save restart files.
   class(equation_nasto_gpu_object), intent(inout) :: self !< The equation.

   if (mod(self%it,self%restart_save)==0) then
      print '(A)', 'save restart files t: '//trim(str(self%it,.true.))//', time: '//trim(str(self%time,.true.))
      call self%adam%save_restart_files(basename=self%restart_basename, t=self%it, time=self%time)
      call self%save_hdf5(output_basename=self%restart_basename)
   endif
   endsubroutine save_restart_files

   subroutine save_slices(self)
   !< Save simulation data slices.
   class(equation_nasto_gpu_object), intent(inout) :: self !< The equation.
   integer(I4P)                                    :: s    !< Slices counter.

   if (self%slices_number>0) then
      do s=1, self%slices_number
         if (self%it>0) then
            if (mod(self%it,self%slice(s)%slice_save)==0.or.self%it==self%t_max.or.&
               (((self%t_max <= 0).and.(self%time >= self%time_max)).or.((self%it>=self%t_max).and.(self%t_max > 0)))) then
               print '(A)', 'save slice n: '//trim(str(s,.true.))//&
                     ', t: '//trim(str(self%it,.true.))//', time: '//trim(str(self%time,.true.))
               call self%adam%save_slice(points=self%slice(s)%slice_points,                               &
                                         itype=trim(self%slice(s)%slice_itype),                           &
                                         basename=trim(self%output_basename)//                            &
                                                  '-slice_'//trim(strz(s,2))//'-'//trim(strz(self%it,9)), &
                                         q=self%field%q,                                                  &
                                         q_name=['rho','rhu','rhv','rhw','rhe'],                          &
                                         phi=self%phi(:,:,:,:,1))
            endif
         endif
      enddo
   endif
   endsubroutine save_slices

   subroutine set_boundary_conditions(self, q_gpu)
   !< Set boundary conditions of equation.
   class(equation_nasto_gpu_object), intent(in)            :: self                  !< The equation.
   real(R8P),                        intent(inout), device :: q_gpu(1:,         &
                                                                    1-self%ngc:,&
                                                                    1-self%ngc:,&
                                                                    1-self%ngc:,1:) !< Conservative variables.

   if (allocated(self%base_gpu%local_map_bc_crown_gpu)) call set_bc(nv=self%nv, ngc=self%ngc,                          &
                                                                    local_map_bc=self%base_gpu%local_map_bc_crown_gpu, &
                                                                    x_cell_gpu=self%x_cell_gpu,                        &
                                                                    y_cell_gpu=self%y_cell_gpu,                        &
                                                                    z_cell_gpu=self%z_cell_gpu,                        &
                                                                    fec_1_6_array_gpu=self%base_gpu%fec_1_6_array_gpu, &
                                                                    q_bc_vars_gpu=self%bc_vars_gpu,                  &
                                                                    gamma_fluid=self%gamma_fluid,                      &
                                                                    dha_star=self%dha_star,                            &
                                                                    cv_star=self%cv_star,                              &
                                                                    R_star=self%R_star)
   contains
      subroutine set_bc(nv, ngc, local_map_bc,                    &
                        x_cell_gpu, y_cell_gpu, z_cell_gpu,       &
                        fec_1_6_array_gpu, q_bc_vars_gpu,         &
                        gamma_fluid, dha_star,  cv_star, R_star)
      integer(I4P), intent(in)         :: nv                      !< Number of variables.
      integer(I4P), intent(in)         :: ngc                     !< Ghost cells number.
      integer(I8P), intent(in), device :: local_map_bc(:,:,:)     !< Local map for BC ghost cells.
      real(R8P),    intent(in), device :: x_cell_gpu(1:,1-ngc:)   !< Conservative variables.
      real(R8P),    intent(in), device :: y_cell_gpu(1:,1-ngc:)   !< Conservative variables.
      real(R8P),    intent(in), device :: z_cell_gpu(1:,1-ngc:)   !< Conservative variables.
      integer(I4P), intent(in), device :: fec_1_6_array_gpu(:)    !< Local map for BC ghost cells.
      real(R8P),    intent(in), device :: q_bc_vars_gpu(:,:)      !< Boundary variables.
      integer(I4P)                     :: b                       !< Counter.
      integer(I4P)                     :: c, i, j, k, v           !< Counter.
      integer(I4P)                     :: idelta                  !< IJK delta step for extrapolation.
      integer(I4P)                     :: jdelta                  !< IJK delta step for extrapolation.
      integer(I4P)                     :: kdelta                  !< IJK delta step for extrapolation.
      integer(I4P)                     :: bc_type                 !< Boundary condition type.
      integer(I4P)                     :: crown                   !< Crown counter.
      integer(I4P)                     :: iercuda                 !< Error trapping flag for CUDAFortran.
      integer(I4P)                     :: fec                     !< Boundary fec (1 to 26).
      integer(I4P)                     :: fec_1_6                 !< Boundary fec (1 to 6).
      real(R8P)                        :: gamma_fluid             !< Gamma.
      real(R8P)                        :: dha_star                !< Entalpy formation.
      real(R8P)                        :: cv_star                 !< Constant volume specific heat.
      real(R8P)                        :: R_star                  !< Gas constant.

      do crown=1, ngc
         !$cuf kernel do(1) <<<*,*>>>
         do c=1, size(local_map_bc, dim=1)
            b = local_map_bc(c, 1 ,crown)
            if (b>0) then
               i       = local_map_bc(c, 2 ,crown)
               j       = local_map_bc(c, 3 ,crown)
               k       = local_map_bc(c, 4 ,crown)
               idelta  = local_map_bc(c, 5 ,crown)
               jdelta  = local_map_bc(c, 6 ,crown)
               kdelta  = local_map_bc(c, 7 ,crown)
               bc_type = local_map_bc(c, 8 ,crown)
               fec     = local_map_bc(c, 9 ,crown)
               fec_1_6 = fec_1_6_array_gpu(fec)
               if (bc_type == BC_EXTRAPOLATION) then
                  do v=1, nv
                     q_gpu(b,i,j,k,v) = q_gpu(b,i-idelta,j-jdelta,k-kdelta,v)
                  enddo
               else if (bc_type == BC_INFLOW) then
                   q_gpu(b,i,j,k,1) = q_bc_vars_gpu(1, fec_1_6)
                   q_gpu(b,i,j,k,2) = q_bc_vars_gpu(1, fec_1_6)* q_bc_vars_gpu(2, fec_1_6)
                   q_gpu(b,i,j,k,3) = q_bc_vars_gpu(1, fec_1_6)* q_bc_vars_gpu(3, fec_1_6)
                   q_gpu(b,i,j,k,4) = q_bc_vars_gpu(1, fec_1_6)* q_bc_vars_gpu(4, fec_1_6)
                   q_gpu(b,i,j,k,5) = q_bc_vars_gpu(1, fec_1_6)*                         &
                       (cv_star*q_bc_vars_gpu(5, fec_1_6)/(q_bc_vars_gpu(1, fec_1_6)*R_star)+ &
                       0.5_R8P*(q_bc_vars_gpu(2, fec_1_6)**2+q_bc_vars_gpu(3, fec_1_6)**2+q_bc_vars_gpu(4, fec_1_6)**2))
               endif
            endif
         enddo
         !@cuf iercuda=cudaDeviceSynchronize()
      enddo
      endsubroutine set_bc
   endsubroutine set_boundary_conditions

   subroutine set_initial_conditions(self)
   !< Set initial conditions of field.
   class(equation_nasto_gpu_object), intent(inout) :: self       !< The equation.
   integer(I4P)                                    :: b, i, j, k !< Counter.
   real(R8P)                                       :: x_split    !< Scalar.
   real(R8P)                                       :: uu, vv, ww !< Scalar.
   real(R8P)                                       :: rn         !< Scalar.

   associate(blocks_number=>self%blocks_number, q=>self%field%q, ni=>self%ni, nj=>self%nj, nk=>self%nk,        &
             ngc=>self%ngc, x_cell=>self%field%x_cell, y_cell=>self%field%y_cell, z_cell=>self%field%z_cell,   &
             gamma_fluid=>self%gamma_fluid, R_star=>self%R_star, cv_star=>self%cv_star, cp_star=>self%cp_star, &
             dha_star=>self%dha_star, ic_vars => self%ic_vars)

      if(self%ic_type == IC_UNIFORM) then
         do b=1, blocks_number
            do k=1, nk
               do j=1, nj
                  do i=1, ni
                     q(1,i,j,k,b) = ic_vars(1)
                     call random_number(rn) ; rn = rn*ic_vars(6) ; uu = ic_vars(2)+rn
                     q(2,i,j,k,b) = ic_vars(1)*uu
                     call random_number(rn) ; rn = rn*ic_vars(6) ; vv = ic_vars(3)+rn
                     q(3,i,j,k,b) = ic_vars(1)*vv
                     call random_number(rn) ; rn = rn*ic_vars(6) ; ww = ic_vars(4)+rn
                     q(4,i,j,k,b) = ic_vars(1)*ww
                     q(5,i,j,k,b) = ic_vars(1)*(cv_star*ic_vars(5)/(ic_vars(1)*R_star)+&
                         0.5_R8P*(uu**2+vv**2+ww**2))
                  enddo
               enddo
            enddo
         enddo
      elseif(self%ic_type == IC_LEFTRIGHT) then
         x_split = ic_vars(11)
         do b=1, blocks_number
            do k=1, nk
               do j=1, nj
                  do i=1, ni
                     if(x_cell(i,b) <= x_split) then
                        q(1,i,j,k,b) = ic_vars(1)
                        q(2,i,j,k,b) = ic_vars(1)*ic_vars(2)
                        q(3,i,j,k,b) = ic_vars(1)*ic_vars(3)
                        q(4,i,j,k,b) = ic_vars(1)*ic_vars(4)
                        q(5,i,j,k,b) = ic_vars(1)*(cv_star*ic_vars(5)/(ic_vars(1)*R_star)+ &
                            0.5_R8P*(ic_vars(2)**2+ic_vars(3)**2+ic_vars(4)**2))
                     else
                        q(1,i,j,k,b) = ic_vars(6)
                        q(2,i,j,k,b) = ic_vars(6)*ic_vars(7)
                        q(3,i,j,k,b) = ic_vars(6)*ic_vars(8)
                        q(4,i,j,k,b) = ic_vars(6)*ic_vars(9)
                        q(5,i,j,k,b) = ic_vars(6)*(cv_star*ic_vars(10)/(ic_vars(6)*R_star)+&
                            0.5_R8P*(ic_vars(7)**2+ic_vars(8)**2+ic_vars(9)**2))
                     endif
                  enddo
               enddo
            enddo
         enddo
      elseif(self%ic_type == IC_FLAME) then
      endif
   endassociate
   call self%copy_cpu_gpu
   endsubroutine set_initial_conditions

   subroutine update_ghost_gpu(self, q_gpu, step)
   !< Update ghost cells.
   !< If not specified all steps are perfermod, syncronous computation
   class(equation_nasto_gpu_object), intent(inout)         :: self            !< The equation.
   real(R8P),                        intent(inout), device :: q_gpu(1:,         &
                                                                    1-self%ngc:,&
                                                                    1-self%ngc:,&
                                                                    1-self%ngc:,&
                                                                    1:)       !< Conservative variables.
   integer(I4P),                     intent(in), optional  :: step            !< Step to be perfordmed in asyncronous comp.
   logical                                                 :: do_local_update !< Flag for triggering local update.
   logical                                                 :: do_set_bc       !< Flag for triggering setting bc.

   ! perform local update if step is not speficied or if first step is selected
   do_local_update = .false.
   do_set_bc       = .false.
   if (.not.present(step)) then
      do_local_update = .true.
      do_set_bc       = .true.
   else
      if (step==1) do_local_update = .true.
      if (step==3) do_set_bc       = .true.
   endif

   if (do_local_update) call self%base_gpu%update_ghost_local_gpu(q_gpu=q_gpu)
                        call self%base_gpu%update_ghost_mpi_gpu(q_gpu=q_gpu, step=step)
   if (do_set_bc)       call self%set_boundary_conditions(q_gpu=q_gpu)
   endsubroutine update_ghost_gpu

   subroutine update_ghost_fluxes_gpu(self, flx_gpu, fly_gpu, flz_gpu, step)
   !< Update ghost cells.
   !< If not specified all steps are perfermod, syncronous computation
   class(equation_nasto_gpu_object), intent(inout)         :: self            !< The equation.
   real(R8P),                        intent(inout), device :: flx_gpu(1:,         &
                                                                      1-self%ngc:,&
                                                                      1-self%ngc:,&
                                                                      1-self%ngc:,&
                                                                      1:)       !< Conservative variables.
   real(R8P),                        intent(inout), device :: fly_gpu(1:,         &
                                                                      1-self%ngc:,&
                                                                      1-self%ngc:,&
                                                                      1-self%ngc:,&
                                                                      1:)       !< Conservative variables.
   real(R8P),                        intent(inout), device :: flz_gpu(1:,         &
                                                                      1-self%ngc:,&
                                                                      1-self%ngc:,&
                                                                      1-self%ngc:,&
                                                                      1:)       !< Conservative variables.
   integer(I4P),                     intent(in), optional  :: step            !< Step to be perfordmed in asyncronous comp.
   logical                                                 :: do_local_update !< Flag for triggering local update.

   ! perform local update if step is not speficied or if first step is selected
   do_local_update = .false.
   if (.not.present(step)) then
      do_local_update = .true.
   else
      if (step==1) do_local_update = .true.
   endif

   if (do_local_update) call self%base_gpu%update_ghost_fluxes_local_gpu(flx_gpu=flx_gpu, fly_gpu=fly_gpu, flz_gpu=flz_gpu)
   !TODO                     call self%base_gpu%update_ghost_fluxes_mpi_gpu(q_gpu=q_gpu, step=step)
   endsubroutine update_ghost_fluxes_gpu

   subroutine update_phi(self)
   !< Update x/y/z_cell_gpu
   class(equation_nasto_gpu_object), intent(inout) :: self                      !< The equation.
   integer(I4P)                                    :: b, i, j, k, ib            !< Counter.
   real(R8P)                                       :: query_x, query_y, query_z !< Query point coordinates.
   real(R8P)                                       :: near_x, near_y, near_z    !< Nearest point coordinates.
   real(R8P)                                       :: distance                  !< Distance from solid.
   logical                                         :: inside                    !< Inside/outside boolean.

   associate(blocks_number=>self%blocks_number, ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, &
             x_cell=>self%field%x_cell, y_cell=>self%field%y_cell, z_cell=>self%field%z_cell,         &
             ptree => self%ptree, phi=>self%phi, phi_gpu=>self%phi_gpu, n_solids=>self%n_solids)
   print '(A)', 'update IB distance'
   do ib=1,n_solids
      do b=1,blocks_number
         do i=1-ngc,ni+ngc
         do j=1-ngc,nj+ngc
         do k=1-ngc,nk+ngc
            query_x = x_cell(i,b)
            query_y = y_cell(j,b)
            query_z = z_cell(k,b)

            ! p = your point
            ! c = centre of the cube
            ! s = half size of the cube
            ! r = the point we are looking for
            !RIMETTEREv = p - c
            !RIMETTEREm = maxval(abs((v)))
            !RIMETTEREr = c + ( v / m * s )

            ! RIMETTERE CGAL
            ! if(query_y < 31.2 .or. query_y > 32.5 .or. query_x < 19.5 .or. query_x > 21.5) then
            !     distance = -100.
            ! else
            !     call polyhedron_closest(ptree(ib),query_x,query_y,query_z,near_x,near_y,near_z)
            !     distance = sqrt((near_x-query_x)**2+(near_y-query_y)**2+(near_z-query_z)**2)
            !     inside   = cgal_polyhedron_inside(ptree(ib),query_x,query_y,query_z)
            !     if(.not.inside) distance = - distance
            ! endif

            !if(inside) print*,'Point inside!!!!!!!!!!!!!!!!'
            ! RIMETTERE CGAL

            distance = - (sqrt((query_x-10.)**2+(query_y-10.)**2+(query_z-10.)**2)-1.)

            phi(b,i,j,k,ib) = distance
         enddo
         enddo
         enddo
      enddo
   enddo
   phi_gpu = phi
   endassociate
   endsubroutine update_phi

   ! operators
   ! =
   subroutine eq_assign_eq(lhs, rhs)
   !< Operator `=`.
   class(equation_nasto_gpu_object), intent(inout) :: lhs !< Left hand side.
   type(equation_nasto_gpu_object),  intent(in)    :: rhs !< Right hand side.

   ! TODO to be written
   endsubroutine eq_assign_eq

   subroutine flame_find_x_v_cuf(ni, nj, nk, ngc, nv, blocks_number, dt, &
       tem_stabil, x_cell_gpu, y_cell_gpu, z_cell_gpu, q_aux_gpu, x_min, x_min_old, velrel)

   !< Advance q_gpu by means of RK stages.
   integer(I4P), intent(in)            :: ni                                     !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj                                     !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk                                     !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: ngc                                    !< Ghost grid number.
   integer(I4P), intent(in)            :: nv                                     !< Number of conservative varibales.
   integer(I4P), intent(in)            :: blocks_number                          !< Number of blocks.
   real(R8P),    intent(in)            :: tem_stabil                             !< Time step.
   real(R8P),    intent(in)            :: Dt                                     !< Time step.
   real(R8P),    intent(in), device    :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Conservative variables.
   real(R8P),    intent(in), device    :: x_cell_gpu(1:,1-ngc:)                  !< Conservative variables.
   real(R8P),    intent(in), device    :: y_cell_gpu(1:,1-ngc:)                  !< Conservative variables.
   real(R8P),    intent(in), device    :: z_cell_gpu(1:,1-ngc:)                  !< Conservative variables.
   real(R8P),    intent(inout)         :: x_min                                  !< Conservative variables.
   real(R8P),    intent(inout)         :: x_min_old                              !< Conservative variables.
   real(R8P),    intent(inout)         :: velrel                                 !< Conservative variables.
   integer(I4P)                        :: i, j, k, b, s, v                       !< Counter.
   integer(I4P)                        :: iercuda                                !< Error trapping flag for CUDAFortran.
   integer(I4P)                        :: ierr                                   !< Error trapping flag for CUDAFortran.

   x_min = huge(1.0_R8P)

   !$cuf kernel do(4) <<<*,*>>>
   do k=1-ngc, nk+ngc
      do j=1-ngc, nj+ngc
         do i=1-ngc, ni+ngc
            do b=1, blocks_number
               if(q_aux_gpu(b,i,j,k,6) > tem_stabil) then
                   x_min = min(x_min, x_cell_gpu(b,i))
               endif
            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   call MPI_ALLREDUCE(MPI_IN_PLACE, x_min, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
   endsubroutine flame_find_x_v_cuf

   subroutine compute_aux_cuf(ni, nj, nk, ngc, ns, blocks_number, gamma_fluid, dha_star, R_star, cv_star, q_gpu,  q_aux_gpu)
   !< Compute auxiliary variables by means of CUF threads.
   integer(I4P), intent(in)          :: ni                                     !< Grid cells number in I direction.
   integer(I4P), intent(in)          :: nj                                     !< Grid cells number in J direction.
   integer(I4P), intent(in)          :: nk                                     !< Grid cells number in K direction.
   integer(I4P), intent(in)          :: ngc                                    !< Ghost cells number.
   integer(I4P), intent(in)          :: ns                                     !< Number of fluid species.
   integer(I4P), intent(in)          :: blocks_number                          !< Number of blocks.
   real(R8P),    intent(in)          :: gamma_fluid                            !< Gamma fluid.
   real(R8P),    intent(in)          :: dha_star                               !< Entalpy fluid.
   real(R8P),    intent(in),  device ::     q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Conservative variables.
   real(R8P),    intent(out), device :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Auxiliary variables.
   integer(I4P)                      :: b, i, j, k, s                          !< Counter.
   real(R8P)                         :: rho, uuu, vvv, www, rhe, rya, yya, tem !< Scalars.
   real(R8P), intent(in)             :: R_star, cv_star                        !< Scalars.
   integer(I4P)                      :: iercuda                                !< Error trapping flag for CUDAFortran.

   !$cuf kernel do(4) <<<*,*>>>
   do k=1-ngc, nk+ngc
      do j=1-ngc, nj+ngc
         do i=1-ngc, ni+ngc
            do b=1, blocks_number
               rho = q_gpu(b,i,j,k,1)
               uuu = q_gpu(b,i,j,k,2)/rho
               vvv = q_gpu(b,i,j,k,3)/rho
               www = q_gpu(b,i,j,k,4)/rho
               rhe = q_gpu(b,i,j,k,5)
               if(ns == 2) then
                   rya = q_gpu(b,i,j,k,ns+4)
               else
                   rya = 0._R8P
               endif
               yya = rya/rho
               tem = ((rhe-rya*dha_star)/rho-0.5*(uuu**2+vvv**2+www**2))/cv_star

               q_aux_gpu(b,i,j,k,1) = rho                            ! density
               q_aux_gpu(b,i,j,k,2) = uuu                            ! velocity x
               q_aux_gpu(b,i,j,k,3) = vvv                            ! velocity y
               q_aux_gpu(b,i,j,k,4) = www                            ! velocity z
               q_aux_gpu(b,i,j,k,5) = yya                            ! mass fraction
               q_aux_gpu(b,i,j,k,6) = tem                            ! temperature
               q_aux_gpu(b,i,j,k,7) = R_star*rho*tem                 ! pressure
               q_aux_gpu(b,i,j,k,8) = rhe/rho+R_star*tem             ! entalpy
               q_aux_gpu(b,i,j,k,9) = sqrt(gamma_fluid*R_star*tem)   ! sound speed
            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine compute_aux_cuf

   subroutine evolve_eikonal_q_gpu_cuf(ni, nj, nk, ngc, nv, phi_gpu, dx_gpu, dy_gpu, dz_gpu, blocks_number, dq_gpu, q_gpu)
   !< Initialize RK stage with q_gpu.
   integer(I4P), intent(in)            :: ni                                  !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj                                  !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk                                  !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: ngc                                 !< Ghost cells number.
   integer(I4P), intent(in)            :: nv                                  !< Number of conservative varibales.
   integer(I4P), intent(in)            :: blocks_number                       !< Number of blocks.
   real(R8P),    intent(in),    device :: phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Distance function.
   real(R8P),    intent(in),    device :: dx_gpu(1:)                          !< X space steps.
   real(R8P),    intent(in),    device :: dy_gpu(1:)                          !< Y space steps.
   real(R8P),    intent(in),    device :: dz_gpu(1:)                          !< Z space steps.
   real(R8P),    intent(inout), device :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !<
   real(R8P),    intent(inout), device :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Conservative variables.
   integer(I4P)                        :: i, j, k, b, v                       !< Counter.
   integer(I4P)                        :: iercuda                             !< Error trapping flag for CUDAFortran.
   type(dim3)                          :: grid, tBlock                        !< CUDA grid and block.

   integer :: error

      error = cudaGetLastError() ; if(error /= cudaSuccess) then
         print*,'BEFORE DQ PRE FRA CUDA ERROR ',cudaGetErrorString(error)
         print*,'THINGS: ',blocks_number, ni
         call MPI_Abort(MPI_COMM_WORLD, -15,error) ; STOP
      endif
   tBlock = dim3(32,8,1) ; grid = dim3(ceiling(real(blocks_number)/tBlock%x),ceiling(real(ni)/tBlock%y),1)
   if(blocks_number > 0) then
   call compute_eikonal_dq_gpu<<<grid, tBlock>>>(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, &
                                                 phi_gpu=phi_gpu, dx_gpu=dx_gpu, dy_gpu=dy_gpu, dz_gpu=dz_gpu,     &
                                                 dq_gpu=dq_gpu, q_gpu=q_gpu)
                                         endif
      error = cudaGetLastError() ; if(error /= cudaSuccess) then
         print*,'AFTER DQ POST FRA CUDA ERROR ',cudaGetErrorString(error)
         print*,'THINGS: ',blocks_number, ni
         call MPI_Abort(MPI_COMM_WORLD, -15,error) ; STOP
      endif

   !$cuf kernel do(4) <<<*,*>>>
   do k=1, nk
      do j=1, nj
         do i=1,ni
            do b=1, blocks_number
               do v=1, nv
                  if (phi_gpu(b,i,j,k,1) > 0._R8P) then
                     q_gpu(b,i,j,k,v) = q_gpu(b,i,j,k,v) - dq_gpu(b,i,j,k,v)
                  endif
               enddo
            enddo
         enddo
      enddo
   enddo
   endsubroutine evolve_eikonal_q_gpu_cuf

   attributes(global) subroutine compute_eikonal_dq_gpu(ni, nj, nk, ngc, nv, blocks_number, &
                                                        phi_gpu, dx_gpu, dy_gpu, dz_gpu,    &
                                                        q_gpu, dq_gpu)
   !< Initialize RK stage with q_gpu.
   integer(I4P), intent(in), value     :: ni                                  !< Grid cells number in I direction.
   integer(I4P), intent(in), value     :: nj                                  !< Grid cells number in J direction.
   integer(I4P), intent(in), value     :: nk                                  !< Grid cells number in K direction.
   integer(I4P), intent(in), value     :: ngc                                 !< Ghost cells number.
   integer(I4P), intent(in), value     :: nv                                  !< Number of conservative varibales.
   integer(I4P), intent(in), value     :: blocks_number                       !< Number of blocks.
   real(R8P),    intent(in),    device :: phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Distance function.
   real(R8P),    intent(in),    device :: dx_gpu(1:)                          !< X space steps.
   real(R8P),    intent(in),    device :: dy_gpu(1:)                          !< Y space steps.
   real(R8P),    intent(in),    device :: dz_gpu(1:)                          !< Z space steps.
   real(R8P),    intent(in),    device :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Conservative variables.
   real(R8P),    intent(inout), device :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !<
   integer(I4P)                        :: i, j, k, b, v                       !< Counter.
   real(R8P)                           :: n_phi_x, n_phi_y, n_phi_z, n_phi    !<

   b = blockDim%x * (blockIdx%x - 1) + threadIdx%x
   i = blockDim%y * (blockIdx%y - 1) + threadIdx%y
   if (b > blocks_number .or. i > ni) return

   do k=1, nk
      do j=1, nj
         if (phi_gpu(b,i,j,k,1) > 0._R8P) then
            n_phi_x = (phi_gpu(b,i+1,j,k,1) - phi_gpu(b,i-1,j,k,1) ) ! / (2 * dx_gpu(b))
            n_phi_y = (phi_gpu(b,i,j+1,k,1) - phi_gpu(b,i,j-1,k,1) ) ! / (2 * dy_gpu(b))
            n_phi_z = (phi_gpu(b,i,j,k+1,1) - phi_gpu(b,i,j,k-1,1) ) ! / (2 * dz_gpu(b))
            n_phi = abs(n_phi_x) + abs(n_phi_y) + abs(n_phi_z) + 10e-12
            n_phi = 0.9_R8P / n_phi
            n_phi_x = n_phi_x * n_phi
            n_phi_y = n_phi_y * n_phi
            n_phi_z = n_phi_z * n_phi
            do v=1, nv
               dq_gpu(b,i,j,k,v) = 0._R8P
            enddo
            if (n_phi_x > 0._R8P) then
               do v=1, nv
                  dq_gpu(b,i,j,k,v) = dq_gpu(b,i,j,k,v) + abs(n_phi_x) * (q_gpu(b,i,j,k,v) - q_gpu(b,i-1,j,k,v))
               enddo
            else
               do v=1, nv
                  dq_gpu(b,i,j,k,v) = dq_gpu(b,i,j,k,v) + abs(n_phi_x) * (q_gpu(b,i,j,k,v) - q_gpu(b,i+1,j,k,v))
               enddo
            endif
            if (n_phi_y > 0._R8P) then
               do v=1, nv
                  dq_gpu(b,i,j,k,v) = dq_gpu(b,i,j,k,v) + abs(n_phi_y) * (q_gpu(b,i,j,k,v) - q_gpu(b,i,j-1,k,v))
               enddo
            else
               do v=1, nv
                  dq_gpu(b,i,j,k,v) = dq_gpu(b,i,j,k,v) + abs(n_phi_y) * (q_gpu(b,i,j,k,v) - q_gpu(b,i,j+1,k,v))
               enddo
            endif
            if (n_phi_z > 0._R8P) then
               do v=1, nv
                  dq_gpu(b,i,j,k,v) = dq_gpu(b,i,j,k,v) + abs(n_phi_z) * (q_gpu(b,i,j,k,v) - q_gpu(b,i,j,k-1,v))
               enddo
            else
               do v=1, nv
                  dq_gpu(b,i,j,k,v) = dq_gpu(b,i,j,k,v) + abs(n_phi_z) * (q_gpu(b,i,j,k,v) - q_gpu(b,i,j,k+1,v))
               enddo
            endif
         endif
      enddo
   enddo
   endsubroutine compute_eikonal_dq_gpu

   subroutine compute_residuals_gpu(self, ni, nj, nk, ngc, ns, blocks_number,   &
                                    dx_gpu , dy_gpu , dz_gpu ,                  &
                                    q_aux_gpu ,  phi_gpu ,      fl_gpu ,        &
                                    flx_gpu   ,    fly_gpu ,      flz_gpu ,     &
                                    fd_conv_gpu , fd_coeff1_gpu ,               &
                                    fd_coeff2_gpu ,                             &
                                    gminus_x , gminus_y , gminus_z ,            &
                                    gplus_x ,   gplus_y  ,  gplus_z  ,          &
                                    euler_scheme , visc_scheme ,                &
                                    lmax , iweno , visc_order ,                 &
                                    visc_law ,                                  &
                                    cp_star  , cv_star , gamma_fluid,           &
                                    R_star   , mu_star , k_star     ,           &
                                    dha_star , Lewis   ,  Zeldovich ,           &
                                    Damkohler    )
   !< Compute residuals of equation.
   class(equation_nasto_gpu_object), intent(inout) :: self                      !< The equation.
   integer(I4P), intent(in)            :: ni                                    !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj                                    !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk                                    !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: ngc                                   !< Ghost cells number.
   integer(I4P), intent(in)            :: ns                                    !< Number of species.
   integer(I4P), intent(in)            :: blocks_number                         !< Number of blocks.
   real(R8P),    intent(in),    device :: dx_gpu(1:)                            !< X space steps.
   real(R8P),    intent(in),    device :: dy_gpu(1:)                            !< Y space steps.
   real(R8P),    intent(in),    device :: dz_gpu(1:)                            !< Z space steps.
   real(R8P),    intent(in),    device :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout), device :: phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Conservative variables.
   real(R8P),    intent(inout), device :: fl_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Positive fluxes.
   real(R8P),    intent(inout), device :: flx_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Positive fluxes.
   real(R8P),    intent(inout), device :: fly_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Positive fluxes.
   real(R8P),    intent(inout), device :: flz_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Positive fluxes.
   real(R8P),    intent(in),    device :: fd_conv_gpu(1:,1:)                    !< Convective coefficients.
   real(R8P),    intent(in),    device :: fd_coeff1_gpu(1:)                     !< Derivatives first coefficients.
   real(R8P),    intent(in),    device :: fd_coeff2_gpu(0:)                     !< Derivatives second coefficients.
   real(R8P),    intent(inout), device :: gplus_x(1:,1:,1:,1:,1:)               !< Auxiliary variables.
   real(R8P),    intent(inout), device :: gminus_x(1:,1:,1:,1:,1:)              !< Auxiliary variables.
   real(R8P),    intent(inout), device :: gplus_y(1:,1:,1:,1:,1:)               !< Auxiliary variables.
   real(R8P),    intent(inout), device :: gminus_y(1:,1:,1:,1:,1:)              !< Auxiliary variables.
   real(R8P),    intent(inout), device :: gplus_z(1:,1:,1:,1:,1:)               !< Auxiliary variables.
   real(R8P),    intent(inout), device :: gminus_z(1:,1:,1:,1:,1:)              !< Auxiliary variables.
   integer(I4P),  intent(in)           :: euler_scheme                          !< Euler scheme.
   integer(I4P),  intent(in)           :: visc_scheme                           !< Diffusive terms scheme.
   integer(I4P), intent(in)            :: lmax                                  !< Conservative stencil size.
   integer(I4P), intent(in)            :: iweno                                 !< Weno order.
   integer(I4P), intent(in)            :: visc_order                            !< Diffusive terms' order.
   integer(I4P), intent(in)            :: visc_law                              !< Viscosity temperature law.
   real(R8P),    intent(in)            :: cp_star                               !< Constant pressure specific heat.
   real(R8P),    intent(in)            :: cv_star                               !< Constant volume specific heat.
   real(R8P),    intent(in)            :: R_star                                !< Gas number.
   real(R8P),    intent(in)            :: gamma_fluid                           !< Gas gamma..
   real(R8P),    intent(in)            :: mu_star                               !< Viscosity.
   real(R8P),    intent(in)            :: k_star                                !< Thermal diffusion.
   real(R8P),    intent(in)            :: dha_star                              !< Formation entalpy.
   real(R8P),    intent(in)            :: Lewis                                 !< Lewis number
   real(R8P),    intent(in)            :: Zeldovich                             !< Zeldovich number.
   real(R8P),    intent(in)            :: Damkohler                             !< Damkohler number.
   real(R8P)                           :: ib_eps                                !< Tolerance immersed boundary delta ratio.
   integer(I4P)                        :: iercuda                               !< Error trapping flag for CUDAFortran.
   type(dim3)                          :: grid, tBlock                          !< CUDA grid and block.
   integer(I4P)                        :: error                               !< Error trapping flag for CUDAFortran.

      error = cudaGetLastError()
      if(error /= cudaSuccess) then
         print*,'FRA-1 CUDA ERROR ',cudaGetErrorString(error)
         call MPI_Abort(MPI_COMM_WORLD, -15,error)
         STOP
      endif

   if(blocks_number > 0) then
       if(euler_scheme == 1) then
          tBlock = dim3(32,8,1) ; grid = dim3(ceiling(real(blocks_number)/tBlock%x),ceiling(real(nj)/tBlock%y),1)
          call euler_x_central_kernel<<<grid, tBlock>>>(q_aux_gpu, flx_gpu, fd_conv_gpu, dx_gpu,     &
                                                        blocks_number, ni, nj, nk, ngc, ns+4, lmax)
          tBlock = dim3(32,8,1) ; grid = dim3(ceiling(real(blocks_number)/tBlock%x),ceiling(real(ni)/tBlock%y),1)
          call euler_y_central_kernel<<<grid, tBlock>>>(q_aux_gpu, fly_gpu, fd_conv_gpu, dy_gpu,     &
                                                        blocks_number, ni, nj, nk, ngc, ns+4, lmax)
          tBlock = dim3(32,8,1) ; grid = dim3(ceiling(real(blocks_number)/tBlock%x),ceiling(real(ni)/tBlock%y),1)
          call euler_z_central_kernel<<<grid, tBlock>>>(q_aux_gpu, flz_gpu, fd_conv_gpu, dz_gpu,     &
                                                        blocks_number, ni, nj, nk, ngc, ns+4, lmax)
       else
          tBlock = dim3(32,8,1) ; grid = dim3(ceiling(real(blocks_number)/tBlock%x),ceiling(real(nj)/tBlock%y),1)
          call euler_x_kernel<<<grid, tBlock>>>(q_aux_gpu, flx_gpu, gplus_x, gminus_x,                                       &
                                                blocks_number, ni, nj, nk, ngc, ns+4, iweno, dha_star, gamma_fluid, R_star, cv_star)
          tBlock = dim3(32,8,1) ; grid = dim3(ceiling(real(blocks_number)/tBlock%x),ceiling(real(ni)/tBlock%y),1)
          call euler_y_kernel<<<grid, tBlock>>>(q_aux_gpu, fly_gpu, gplus_y, gminus_y,                                        &
                                                blocks_number, ni, nj, nk, ngc, ns+4, iweno, dha_star, gamma_fluid, R_star, cv_star)
          tBlock = dim3(32,8,1) ; grid = dim3(ceiling(real(blocks_number)/tBlock%x),ceiling(real(ni)/tBlock%y),1)
          call euler_z_kernel<<<grid, tBlock>>>(q_aux_gpu, flz_gpu, gplus_z, gminus_z,                                        &
                                                blocks_number, ni, nj, nk, ngc, ns+4, iweno, dha_star, gamma_fluid, R_star, cv_star)
       endif
   endif

   ! ! debug restart
   ! if (self%itt == 51) then
   ! print*, ' cazzo residual euler 1'
   ! self%pbuffer = flx_gpu
   ! print '(A)', 'debug-restart flx_gpu[rho*v,ni-3:ni,1   ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 1, 4, 3)))
   ! print '(A)', 'debug-restart flx_gpu[rho*v,ni-3:ni,nj  ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 3)))
   ! self%pbuffer = fly_gpu
   ! print '(A)', 'debug-restart fly_gpu[rho*v,ni-3:ni,1   ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 1, 4, 3)))
   ! print '(A)', 'debug-restart fly_gpu[rho*v,ni-3:ni,nj  ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 3)))
   ! self%pbuffer = flz_gpu
   ! print '(A)', 'debug-restart flz_gpu[rho*v,ni-3:ni,1   ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 1, 4, 3)))
   ! print '(A)', 'debug-restart flz_gpu[rho*v,ni-3:ni,nj  ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 3)))
   ! endif
   ! ! debug restart

      error = cudaGetLastError()
      if(error /= cudaSuccess) then
         print*,'FRA-2 CUDA ERROR ',cudaGetErrorString(error)
         call MPI_Abort(MPI_COMM_WORLD, -15,error)
         STOP
      endif

   if(visc_scheme == 1) then
      !if(mu_star > 0.) call viscous_cuf(ni, nj, nk, ngc, blocks_number, ivis, visc_type, fd_coeff1_gpu, fd_coeff2_gpu, &
      !   gamma_fluid, Prandtl, q_coeff, Lewis, Zeldovich, Damkohler, dha, q_aux_gpu, dx_gpu, dy_gpu, dz_gpu, fl_gpu)
   else
      if(mu_star > 0.) call viscous_part(blocks_number, ni, nj, nk, ngc, ns+4, mu_star, k_star, &
                                         q_aux_gpu, flx_gpu, fly_gpu, flz_gpu,                           &
                                         dx_gpu, dy_gpu, dz_gpu)
   endif
      !@cuf iercuda=cudaDeviceSynchronize()

   ! ! debug restart
   ! if (self%itt == 51) then
   ! print*, ' cazzo residual visc 2'
   ! self%pbuffer = flx_gpu
   ! print '(A)', 'debug-restart flx_gpu[rho*v,ni-3:ni,1   ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 1, 4, 3)))
   ! print '(A)', 'debug-restart flx_gpu[rho*v,ni-3:ni,nj  ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 3)))
   ! self%pbuffer = fly_gpu
   ! print '(A)', 'debug-restart fly_gpu[rho*v,ni-3:ni,1   ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 1, 4, 3)))
   ! print '(A)', 'debug-restart fly_gpu[rho*v,ni-3:ni,nj  ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 3)))
   ! self%pbuffer = flz_gpu
   ! print '(A)', 'debug-restart flz_gpu[rho*v,ni-3:ni,1   ,nk/2]'//trim(str(self%pbuffer(296, 6:8, 1, 4, 3)))
   ! print '(A)', 'debug-restart flz_gpu[rho*v,ni-3:ni,nj  ,nk/2]'//trim(str(self%pbuffer(167, 6:8, 8, 4, 3)))
   ! endif
   ! ! debug restart

      error = cudaGetLastError()
      if(error /= cudaSuccess) then
         print*,'FRA-3 CUDA ERROR ',cudaGetErrorString(error)
         call MPI_Abort(MPI_COMM_WORLD, -15,error)
         STOP
      endif

   !RIMETTEREcall self%update_ghost_fluxes_gpu(flx_gpu, fly_gpu, flz_gpu)

   ib_eps = 1.e-12_R8P
   call compute_flux_diff(blocks_number, ni, nj, nk, ngc, ns+4, &
                          fl_gpu, flx_gpu, fly_gpu, flz_gpu, phi_gpu, &
                          dx_gpu, dy_gpu, dz_gpu, ib_eps)
      error = cudaGetLastError()
      if(error /= cudaSuccess) then
         print*,'FRA-4 CUDA ERROR ',cudaGetErrorString(error)
         call MPI_Abort(MPI_COMM_WORLD, -15,error)
         STOP
      endif

   endsubroutine compute_residuals_gpu

   attributes(global) subroutine euler_x_central_kernel(q_aux_gpu, flx_gpu, fd_conv_gpu, dx_gpu, &
                                                        blocks_number, ni, nj, nk, ngc, nv, lmax)

   real(R8P), intent(in), device     :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::   flx_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(in), device     :: fd_conv_gpu(1:, 1:)
   real(R8P), intent(in), device     :: dx_gpu(1:)
   integer, intent(in), value        :: blocks_number, ni, nj, nk, ngc, nv, lmax

   real(R8P)                         :: rhom, uui, vvi, wwi, ppi, enti, yai
   real(R8P)                         :: uuip, vvip, wwip, ppip, entip, yaip
   real(R8P)                         :: ft1, ft2, ft3, ft4, ft5, ft6, ft7
   real(R8P)                         :: uvs1, uvs2, uvs3, uvs4, uvs5, uvs6, uvs7, uv_part
   integer                           :: b, i, j, k, l, v, m

   b = blockDim%x * (blockIdx%x - 1) + threadIdx%x
   j = blockDim%y * (blockIdx%y - 1) + threadIdx%y
   if(b > blocks_number .or. j > nj) return

   do k=1,nk

      do i=0,ni ! loop on faces
         ft1 = 0._R8P ; ft2 = 0._R8P ; ft3 = 0._R8P ; ft4 = 0._R8P ; ft5 = 0._R8P ; ft6 = 0._R8P ; ft7 = 0._R8P
         do l=1,lmax
             uvs1 = 0._R8P ; uvs2 = 0._R8P ; uvs3 = 0._R8P ; uvs4 = 0._R8P ; uvs5 = 0._R8P ; uvs6 = 0._R8P ; uvs7 = 0._R8P
             do m=0,l-1
                 rhom    = q_aux_gpu(b,i-m,j,k,1) + q_aux_gpu(b,i-m+l,j,k,1)

                 uui     = q_aux_gpu(b,i-m,j,k,2)
                 vvi     = q_aux_gpu(b,i-m,j,k,3)
                 wwi     = q_aux_gpu(b,i-m,j,k,4)
                 ppi     = q_aux_gpu(b,i-m,j,k,7)
                 enti    = q_aux_gpu(b,i-m,j,k,8)
                 yai     = q_aux_gpu(b,i-m,j,k,5)

                 uuip    = q_aux_gpu(b,i-m+l,j,k,2)
                 vvip    = q_aux_gpu(b,i-m+l,j,k,3)
                 wwip    = q_aux_gpu(b,i-m+l,j,k,4)
                 ppip    = q_aux_gpu(b,i-m+l,j,k,7)
                 entip   = q_aux_gpu(b,i-m+l,j,k,8)
                 yaip    = q_aux_gpu(b,i-m+l,j,k,5)

                 uv_part = (uui+uuip) * rhom
                 uvs1    = uvs1 + uv_part * (2._R8P)
                 uvs2    = uvs2 + uv_part * (uui+uuip)
                 uvs3    = uvs3 + uv_part * (vvi+vvip)
                 uvs4    = uvs4 + uv_part * (wwi+wwip)
                 uvs5    = uvs5 + uv_part * (enti+entip)
                 uvs6    = uvs6 + (2._R8P)*(ppi+ppip)
                 uvs7    = uvs7 + uv_part * (yai+yaip)
             enddo
             ft1 = ft1 + fd_conv_gpu(l,lmax)*uvs1
             ft2 = ft2 + fd_conv_gpu(l,lmax)*uvs2
             ft3 = ft3 + fd_conv_gpu(l,lmax)*uvs3
             ft4 = ft4 + fd_conv_gpu(l,lmax)*uvs4
             ft5 = ft5 + fd_conv_gpu(l,lmax)*uvs5
             ft6 = ft6 + fd_conv_gpu(l,lmax)*uvs6
             ft7 = ft7 + fd_conv_gpu(l,lmax)*uvs7
         enddo
!
         flx_gpu(b,i,j,k,1) = 0.25_R8P*ft1
         flx_gpu(b,i,j,k,2) = 0.25_R8P*ft2 + 0.5_R8P*ft6
         flx_gpu(b,i,j,k,3) = 0.25_R8P*ft3
         flx_gpu(b,i,j,k,4) = 0.25_R8P*ft4
         flx_gpu(b,i,j,k,5) = 0.25_R8P*ft5
         if(nv == 6) flx_gpu(b,i,j,k,6) = 0.25_R8P*ft7
      enddo

   enddo
   endsubroutine euler_x_central_kernel

   attributes(global) subroutine euler_y_central_kernel(q_aux_gpu, fly_gpu, fd_conv_gpu, dy_gpu, &
                                                        blocks_number, ni, nj, nk, ngc, nv, lmax)

   real(R8P), intent(in), device     :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::   fly_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(in), device     :: fd_conv_gpu(1:, 1:)
   real(R8P), intent(in), device     :: dy_gpu(1:)
   integer, intent(in), value        :: blocks_number, ni, nj, nk, ngc, nv, lmax

   real(R8P)                         :: rhom, uui, vvi, wwi, ppi, enti, yai
   real(R8P)                         :: uuip, vvip, wwip, ppip, entip, yaip
   real(R8P)                         :: ft1, ft2, ft3, ft4, ft5, ft6, ft7
   real(R8P)                         :: uvs1, uvs2, uvs3, uvs4, uvs5, uvs6, uvs7, uv_part
   integer                           :: b, i, j, k, l, v, m

   b = blockDim%x * (blockIdx%x - 1) + threadIdx%x
   i = blockDim%y * (blockIdx%y - 1) + threadIdx%y
   if(b > blocks_number .or. i > ni) return

   do k=1,nk

      do j=0,nj ! loop on faces
         ft1 = 0._R8P ; ft2 = 0._R8P ; ft3 = 0._R8P ; ft4 = 0._R8P ; ft5 = 0._R8P ; ft6 = 0._R8P ; ft7 = 0._R8P
         do l=1,lmax
             uvs1 = 0._R8P ; uvs2 = 0._R8P ; uvs3 = 0._R8P ; uvs4 = 0._R8P ; uvs5 = 0._R8P ; uvs6 = 0._R8P ; uvs7 = 0._R8P
             do m=0,l-1
                 rhom    = q_aux_gpu(b,i,j-m,k,1) + q_aux_gpu(b,i,j-m+l,k,1)

                 uui     = q_aux_gpu(b,i,j-m,k,2)
                 vvi     = q_aux_gpu(b,i,j-m,k,3)
                 wwi     = q_aux_gpu(b,i,j-m,k,4)
                 ppi     = q_aux_gpu(b,i,j-m,k,7)
                 enti    = q_aux_gpu(b,i,j-m,k,8)
                 yai     = q_aux_gpu(b,i,j-m,k,5)

                 uuip    = q_aux_gpu(b,i,j-m+l,k,2)
                 vvip    = q_aux_gpu(b,i,j-m+l,k,3)
                 wwip    = q_aux_gpu(b,i,j-m+l,k,4)
                 ppip    = q_aux_gpu(b,i,j-m+l,k,7)
                 entip   = q_aux_gpu(b,i,j-m+l,k,8)
                 yaip    = q_aux_gpu(b,i,j-m+l,k,5)

                 uv_part = (vvi+vvip) * rhom
                 uvs1    = uvs1 + uv_part * (2._R8P)
                 uvs2    = uvs2 + uv_part * (uui+uuip)
                 uvs3    = uvs3 + uv_part * (vvi+vvip)
                 uvs4    = uvs4 + uv_part * (wwi+wwip)
                 uvs5    = uvs5 + uv_part * (enti+entip)
                 uvs6    = uvs6 + (2._R8P)*(ppi+ppip)
                 uvs7    = uvs7 + uv_part * (yai+yaip)
             enddo
             ft1 = ft1 + fd_conv_gpu(l,lmax)*uvs1
             ft2 = ft2 + fd_conv_gpu(l,lmax)*uvs2
             ft3 = ft3 + fd_conv_gpu(l,lmax)*uvs3
             ft4 = ft4 + fd_conv_gpu(l,lmax)*uvs4
             ft5 = ft5 + fd_conv_gpu(l,lmax)*uvs5
             ft6 = ft6 + fd_conv_gpu(l,lmax)*uvs6
             ft7 = ft7 + fd_conv_gpu(l,lmax)*uvs7
         enddo
!
         fly_gpu(b,i,j,k,1) = 0.25_R8P*ft1
         fly_gpu(b,i,j,k,2) = 0.25_R8P*ft2
         fly_gpu(b,i,j,k,3) = 0.25_R8P*ft3 + 0.5_R8P*ft6
         fly_gpu(b,i,j,k,4) = 0.25_R8P*ft4
         fly_gpu(b,i,j,k,5) = 0.25_R8P*ft5
         if(nv == 6) fly_gpu(b,i,j,k,6) = 0.25_R8P*ft7
      enddo

   enddo
   endsubroutine euler_y_central_kernel

   attributes(global) subroutine euler_z_central_kernel(q_aux_gpu, flz_gpu, fd_conv_gpu, dz_gpu, &
                                                        blocks_number, ni, nj, nk, ngc, nv, lmax)

   real(R8P), intent(in), device     :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::   flz_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(in), device     :: fd_conv_gpu(1:, 1:)
   real(R8P), intent(in), device     :: dz_gpu(1:)
   integer, intent(in), value        :: blocks_number, ni, nj, nk, ngc, nv, lmax

   real(R8P)                         :: rhom, uui, vvi, wwi, ppi, enti, yai
   real(R8P)                         :: uuip, vvip, wwip, ppip, entip, yaip
   real(R8P)                         :: ft1, ft2, ft3, ft4, ft5, ft6, ft7
   real(R8P)                         :: uvs1, uvs2, uvs3, uvs4, uvs5, uvs6, uvs7, uv_part
   integer                           :: b, i, j, k, l, v, m

   b = blockDim%x * (blockIdx%x - 1) + threadIdx%x
   i = blockDim%y * (blockIdx%y - 1) + threadIdx%y
   if(b > blocks_number .or. i > ni) return

   do j=1,nj

      do k=0,nk ! loop on faces
         ft1 = 0._R8P ; ft2 = 0._R8P ; ft3 = 0._R8P ; ft4 = 0._R8P ; ft5 = 0._R8P ; ft6 = 0._R8P ; ft7 = 0._R8P
         do l=1,lmax
             uvs1 = 0._R8P ; uvs2 = 0._R8P ; uvs3 = 0._R8P ; uvs4 = 0._R8P ; uvs5 = 0._R8P ; uvs6 = 0._R8P ; uvs7 = 0._R8P
             do m=0,l-1
                 rhom    = q_aux_gpu(b,i,j,k-m,1) + q_aux_gpu(b,i,j,k-m+l,1)

                 uui     = q_aux_gpu(b,i,j,k-m,2)
                 vvi     = q_aux_gpu(b,i,j,k-m,3)
                 wwi     = q_aux_gpu(b,i,j,k-m,4)
                 ppi     = q_aux_gpu(b,i,j,k-m,7)
                 enti    = q_aux_gpu(b,i,j,k-m,8)
                 yai     = q_aux_gpu(b,i,j,k-m,5)

                 uuip    = q_aux_gpu(b,i,j,k-m+l,2)
                 vvip    = q_aux_gpu(b,i,j,k-m+l,3)
                 wwip    = q_aux_gpu(b,i,j,k-m+l,4)
                 ppip    = q_aux_gpu(b,i,j,k-m+l,7)
                 entip   = q_aux_gpu(b,i,j,k-m+l,8)
                 yaip    = q_aux_gpu(b,i,j,k-m+l,5)

                 uv_part = (wwi+wwip) * rhom
                 uvs1    = uvs1 + uv_part * (2._R8P)
                 uvs2    = uvs2 + uv_part * (uui+uuip)
                 uvs3    = uvs3 + uv_part * (vvi+vvip)
                 uvs4    = uvs4 + uv_part * (wwi+wwip)
                 uvs5    = uvs5 + uv_part * (enti+entip)
                 uvs6    = uvs6 + (2._R8P)*(ppi+ppip)
                 uvs7    = uvs7 + uv_part * (yai+yaip)
             enddo
             ft1 = ft1 + fd_conv_gpu(l,lmax)*uvs1
             ft2 = ft2 + fd_conv_gpu(l,lmax)*uvs2
             ft3 = ft3 + fd_conv_gpu(l,lmax)*uvs3
             ft4 = ft4 + fd_conv_gpu(l,lmax)*uvs4
             ft5 = ft5 + fd_conv_gpu(l,lmax)*uvs5
             ft6 = ft6 + fd_conv_gpu(l,lmax)*uvs6
             ft7 = ft7 + fd_conv_gpu(l,lmax)*uvs7
         enddo
!
         flz_gpu(b,i,j,k,1) = 0.25_R8P*ft1
         flz_gpu(b,i,j,k,2) = 0.25_R8P*ft2
         flz_gpu(b,i,j,k,3) = 0.25_R8P*ft3
         flz_gpu(b,i,j,k,4) = 0.25_R8P*ft4 + 0.5_R8P*ft6
         flz_gpu(b,i,j,k,5) = 0.25_R8P*ft5
         if(nv == 6) flz_gpu(b,i,j,k,6) = 0.25_R8P*ft7
      enddo

   enddo
   endsubroutine euler_z_central_kernel

   attributes(global) subroutine euler_x_kernel(q_aux_gpu, flx_gpu, gplus, gminus,  &
                                                blocks_number, ni, nj, nk, ngc, nv, iweno, dha_star, gamma_fluid, R_star, cv_star)

   real(R8P), intent(in), device     :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::   flx_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::     gplus(1:, 1:, 1:, 1:, 1:)
   real(R8P), intent(inout), device  ::    gminus(1:, 1:, 1:, 1:, 1:)
   integer, intent(in), value        :: blocks_number, ni, nj, nk, ngc, nv, iweno
   real(R8P), intent(in), value      :: dha_star, gamma_fluid, R_star, cv_star
   integer                           :: b, i, j, k, l, ll, m, mm, v
   ! here 6 is used instead of nv to help the compiler to use registers instead of global memory
   real(R8P)                         :: er(5,5), el(5,5), ev(5), evmax(5), ghat(5), gl(5), gr(5), fi(5), vi(5)
   !real(R8P)                         :: er(nv,nv), el(nv,nv), ev(nv), evmax(nv), ghat(nv), gl(nv), gr(nv), fi(nv), vi(nv)
   real(R8P)                         :: uu, vv, ww, h, ya, qq, c, ci, b1, b2
   real(R8P)                         :: gc, wc

   b = blockDim%x * (blockIdx%x - 1) + threadIdx%x
   j = blockDim%y * (blockIdx%y - 1) + threadIdx%y
   if(b > blocks_number .or. j > nj) return

   do k=1,nk

      do i=0,ni ! loop on faces

         ! Compute Roe average
         call compute_roe_average(q_aux_gpu, dha_star, gamma_fluid, R_star, &
            ngc, b, i, j, k, i+1, j, k, uu, vv, ww, h, ya, qq, c, ci, b1, b2)

         ! Compute right and left eigenvectors matrices (at Roe state)
         er(1,1) = 1._R8P ;  er(1,2) = uu-c    ; er(1,3) = vv     ; er(1,4) = ww     ; er(1,5) = h-uu*c
         er(2,1) = 1._R8P ;  er(2,2) = uu      ; er(2,3) = vv     ; er(2,4) = ww     ; er(2,5) = qq
         er(3,1) = 1._R8P ;  er(3,2) = uu+c    ; er(3,3) = vv     ; er(3,4) = ww     ; er(3,5) = h+uu*c
         er(4,1) = 0._R8P ;  er(4,2) = 0._R8P  ; er(4,3) = 1._R8P ; er(4,4) = 0._R8P ; er(4,5) = 0._R8P
         er(5,1) = 0._R8P ;  er(5,2) = 0._R8P  ; er(5,3) = 0._R8P ; er(5,4) = 1._R8P ; er(5,5) = 0._R8P

         el(1,1) =  0.5_R8P*(b1+uu*ci) ; el(1,2) = 1._R8P-b1 ; el(1,3) =  0.5_R8P*(b1-uu*ci)
         el(2,1) = -0.5_R8P*(b2*uu+ci) ; el(2,2) = b2*uu     ; el(2,3) = -0.5_R8P*(b2*uu-ci)
         el(3,1) = -0.5_R8P*(b2*vv   ) ; el(3,2) = b2*vv     ; el(3,3) = -0.5_R8P*(b2*vv   )
         el(4,1) = -0.5_R8P*(b2*ww   ) ; el(4,2) = b2*ww     ; el(4,3) = -0.5_R8P*(b2*ww   )
         el(5,1) =  0.5_R8P*b2         ; el(5,2) = -b2       ; el(5,3) =  0.5_R8P*b2

         el(1,4) = -vv                 ; el(1,5) = -ww
         el(2,4) = 0._R8P              ; el(2,5) = 0._R8P
         el(3,4) = 1._R8P              ; el(3,5) = 0._R8P
         el(4,4) = 0._R8P              ; el(4,5) = 1._R8P
         el(5,4) = 0._R8P              ; el(5,5) = 0._R8P

         if(nv == 6) then
         ! Compute right and left eigenvectors matrices (at Roe state)
                                                                                                          er(1,6) = ya
                                                                                                          er(2,6) = 0._R8P
                                                                                                          er(3,6) = ya
                                                                                                          er(4,6) = -vv/dha_star
                                                                                                          er(5,6) = -ww/dha_star
         er(6,1) = 0._R8P ;  er(6,2) = 0._R8P  ; er(6,3) = 0._R8P ; er(6,4) = 0._R8P ; er(6,5) = 1._R8P ; er(6,6) = 1._R8P/dha_star

         el(6,1) = -0.5_R8P*b2*dha_star     ; el(6,2) = b2*dha_star    ; el(6,3) = -0.5_R8P*b2*dha_star

                                                               el(1,6) = -ya*dha_star*b1-vv**2-ww**2
                                                               el(2,6) = uu*ya*dha_star*b2
                                                               el(3,6) = vv*(1._R8P+ya*dha_star*b2)
                                                               el(4,6) = ww*(1._R8P+ya*dha_star*b2)
                                                               el(5,6) = -ya*dha_star*b2
         el(6,4) = 0._R8P              ; el(6,5) = 0._R8P    ; el(6,6) = dha_star*(1._R8P+ya*dha_star*b2)
         endif

         ! Find max eigenvalues on the stencil
         do m=1,nv  ! loop on characteristic fields
            evmax(m) = -1._R8P
         enddo
         do l=1,2*iweno ! LLF
            ll = i + l - iweno
            uu = q_aux_gpu(b,ll,j,k,2)
            c  = q_aux_gpu(b,ll,j,k,9)
            ev(1) = abs(uu-c) ; ev(2) = abs(uu) ; ev(3) = abs(uu+c) ; ev(4) = ev(2) ; ev(5) = ev(2)
            if(nv == 6) ev(6) = ev(2)
            do m=1,nv
                evmax(m) = max(ev(m),evmax(m))
            enddo
         enddo

         ! Decompose fluxes as + and -
         do l=1,2*iweno ! loop over the stencil centered at face i
            ll = i + l - iweno
            vi(1) = q_aux_gpu(b,ll,j,k,1)
            vi(2) = vi(1)*q_aux_gpu(b,ll,j,k,2)
            vi(3) = vi(1)*q_aux_gpu(b,ll,j,k,3)
            vi(4) = vi(1)*q_aux_gpu(b,ll,j,k,4)
            vi(5) = vi(1)*(cv_star*q_aux_gpu(b,ll,j,k,6)+&
                0.5_R8P*(q_aux_gpu(b,ll,j,k,2)**2+q_aux_gpu(b,ll,j,k,3)**2+q_aux_gpu(b,ll,j,k,4)**2)+&
                q_aux_gpu(b,ll,j,k,5)*dha_star)
            if(nv==6) vi(6) = vi(1)*q_aux_gpu(b,ll,j,k,5)
            fi(1) = vi(2)
            fi(2) = fi(1) * q_aux_gpu(b,ll,j,k,2) + q_aux_gpu(b,ll,j,k,7)
            fi(3) = fi(1) * q_aux_gpu(b,ll,j,k,3)
            fi(4) = fi(1) * q_aux_gpu(b,ll,j,k,4)
            fi(5) = fi(1) * vi(5) / vi(1) + q_aux_gpu(b,ll,j,k,7)*q_aux_gpu(b,ll,j,k,2)
            if(nv==6) fi(6) = fi(1) * q_aux_gpu(b,ll,j,k,5)
            do m=1,nv
               wc = 0._R8P
               gc = 0._R8P
               do mm=1,nv
                  wc = wc + el(mm,m) * vi(mm)
                  gc = gc + el(mm,m) * fi(mm)
               enddo
               gplus (m,l,j,k,b) = 0.5_R8P * (gc + evmax(m) * wc)
               gminus(m,l,j,k,b) = gc - gplus(m,l,j,k,b)
            enddo
         enddo

         ! Reconstruction of the + and - fluxes
         call weno_reconstruction(nv, gplus(1,1,j,k,b), gminus(1,1,j,k,b), gl, gr, iweno)

         ! Reassemble + and - characteristic fluxes
         do m=1,nv
            ghat(m) = gl(m) + gr(m)
         enddo

         ! Return to conservative fluxes
         do m=1,nv
            flx_gpu(b,i,j,k,m) = 0._R8P
            do mm=1,nv
               flx_gpu(b,i,j,k,m) = flx_gpu(b,i,j,k,m) + er(mm,m) * ghat(mm)
            enddo
         enddo

      enddo
   enddo
   endsubroutine euler_x_kernel

   subroutine compute_flux_diff(blocks_number, ni, nj, nk, ngc, nv, &
                                fl_gpu, flx_gpu, fly_gpu, flz_gpu, phi_gpu, &
                                dx_gpu, dy_gpu, dz_gpu, ib_eps)
       integer(I4P), intent(in)          :: blocks_number, ni, nj, nk, ngc, nv
       real(R8P), intent(inout), device  ::  fl_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
       real(R8P), intent(in)   , device  :: flx_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
       real(R8P), intent(in)   , device  :: fly_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
       real(R8P), intent(in)   , device  :: flz_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
       real(R8P), intent(in)   , device  :: phi_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
       real(R8P), intent(in)   , device  :: dx_gpu(1:), dy_gpu(1:), dz_gpu(1:)
       real(R8P), intent(in)             :: ib_eps
       real(R8P)                         :: delta_x, delta_y, delta_z, dx_locale, dy_locale, dz_locale
       integer(I4P)                      :: b, i, j, k, v, iercuda

      !$cuf kernel do(4) <<<*,*>>>
      do k=1,nk
      do j=1,nj
      do i=1,ni
      do b=1,blocks_number

         dx_locale = dx_gpu(b)
         ! Update net flux (procedura alternativa all'interpolazione proposta nel paper, utilizza dx_locale).
         if(phi_gpu(b,i,j,k,1)<0.) then
             if(phi_gpu(b,i+1,j,k,1)*phi_gpu(b,i-1,j,k,1)<0) then
                 if(phi_gpu(b,i+1,j,k,1)>0.) then
                     delta_x = -phi_gpu(b,i,j,k,1)/(phi_gpu(b,i+1,j,k,1)-phi_gpu(b,i,j,k,1)+ib_eps)*dx_gpu(b)
                     dx_locale = dx_gpu(b)/2 + delta_x
                 else !if(phi_gpu(b,i-1,j,k,1)>0) then
                     delta_x = -phi_gpu(b,i,j,k,1)/(phi_gpu(b,i-1,j,k,1)-phi_gpu(b,i,j,k,1)+ib_eps)*dx_gpu(b)
                     dx_locale = dx_gpu(b)/2 + delta_x
                 endif
             endif
         endif

         dy_locale = dy_gpu(b)
         if(phi_gpu(b,i,j,k,1)<0.) then
             if(phi_gpu(b,i,j+1,k,1)*phi_gpu(b,i,j-1,k,1)<0) then
                 if(phi_gpu(b,i,j+1,k,1)>0.) then
                     delta_y = -phi_gpu(b,i,j,k,1)/(phi_gpu(b,i,j+1,k,1)-phi_gpu(b,i,j,k,1)+ib_eps)*dy_gpu(b)
                     dy_locale = dy_gpu(b)/2 + delta_y
                 else !if(phi_gpu(b,i-1,j,k,1)>0) then
                     delta_y = -phi_gpu(b,i,j,k,1)/(phi_gpu(b,i,j-1,k,1)-phi_gpu(b,i,j,k,1)+ib_eps)*dy_gpu(b)
                     dy_locale = dy_gpu(b)/2 + delta_y
                 endif
             endif
         endif

         dz_locale = dz_gpu(b)
         if(phi_gpu(b,i,j,k,1)<0.) then
             if(phi_gpu(b,i,j,k+1,1)*phi_gpu(b,i,j,k-1,1)<0) then
                 if(phi_gpu(b,i,j,k+1,1)>0.) then
                     delta_z = -phi_gpu(b,i,j,k,1)/(phi_gpu(b,i,j,k+1,1)-phi_gpu(b,i,j,k,1)+ib_eps)*dz_gpu(b)
                     dz_locale = dz_gpu(b)/2 + delta_z
                 else !if(phi_gpu(b,i,j,k-1,1)>0) then
                     delta_z = -phi_gpu(b,i,j,k,1)/(phi_gpu(b,i,j,k-1,1)-phi_gpu(b,i,j,k,1)+ib_eps)*dz_gpu(b)
                     dz_locale = dz_gpu(b)/2 + delta_z
                 endif
             endif
         endif

         do v=1,nv
            fl_gpu(b,i,j,k,v) = - (flx_gpu(b,i,j,k,v)-flx_gpu(b,i-1,j,k,v))/dx_locale &
                                - (fly_gpu(b,i,j,k,v)-fly_gpu(b,i,j-1,k,v))/dy_locale &
                                - (flz_gpu(b,i,j,k,v)-flz_gpu(b,i,j,k-1,v))/dz_locale
         enddo

      enddo
      enddo
      enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine compute_flux_diff

   subroutine viscous_part(blocks_number, ni, nj, nk, ngc, nv, mu_star, k_star, &
                           q_aux_gpu, flx_gpu, fly_gpu, flz_gpu,                         &
                           dx_gpu, dy_gpu, dz_gpu)

   integer, intent(in)               :: blocks_number, ni, nj, nk, ngc, nv
   real(R8P), intent(in)             :: mu_star, k_star
   real(R8P), intent(in), device     :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::   flx_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::   fly_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::   flz_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(in), device     ::    dx_gpu(1:), dy_gpu(1:), dz_gpu(1:)
   integer                           :: b, i, j, k, v, iercuda
   real(R8P)                         :: du_dx, dv_dx, dw_dx, du_dy, dv_dy, dw_dy, du_dz, dv_dz, dw_dz
   real(R8P)                         :: dx_locale, dy_locale, dz_locale
   real(R8P)                         :: delta_x, delta_y, delta_z
   real(R8P)                         :: sigq, sigl
   real(R8P)                         :: tau_1_1, tau_2_1, tau_3_1, dT_dx
   real(R8P)                         :: tau_1_2, tau_2_2, tau_3_2, dT_dy
   real(R8P)                         :: tau_1_3, tau_2_3, tau_3_3, dT_dz
   real(R8P)                         :: vel_u, vel_v, vel_w
   real(R8P)                         :: mu, k_coeff
   real(R8P), parameter              :: ib_eps=1.e-12_R8P

   mu       = mu_star
   k_coeff  = k_star

   !$cuf kernel do(3) <<<*,*>>>
   do k=1,nk
      do j=1,nj
         do b=1,blocks_number
            do i=0,ni ! loop on faces
                du_dx = (q_aux_gpu(b,i+1,j,k,2)-q_aux_gpu(b,i,j,k,2))/dx_gpu(b)
                dv_dx = (q_aux_gpu(b,i+1,j,k,3)-q_aux_gpu(b,i,j,k,3))/dx_gpu(b)
                dw_dx = (q_aux_gpu(b,i+1,j,k,4)-q_aux_gpu(b,i,j,k,4))/dx_gpu(b)

                du_dy = (q_aux_gpu(b,i+1,j+1,k,2) - q_aux_gpu(b,i+1,j-1,k,2)+ &
                         q_aux_gpu(b,i,j+1,k,2)   - q_aux_gpu(b,i,j-1,k,2) )*0.25_R8P/dy_gpu(b)
                dv_dy = (q_aux_gpu(b,i+1,j+1,k,3) - q_aux_gpu(b,i+1,j-1,k,3)+ &
                         q_aux_gpu(b,i,j+1,k,3)   - q_aux_gpu(b,i,j-1,k,3) )*0.25_R8P/dy_gpu(b)

                du_dz = (q_aux_gpu(b,i+1,j,k+1,2) - q_aux_gpu(b,i+1,j,k-1,2)+ &
                         q_aux_gpu(b,i,j,k+1,2)   - q_aux_gpu(b,i,j,k-1,2) )*0.25_R8P/dz_gpu(b)
                dw_dz = (q_aux_gpu(b,i+1,j,k+1,4) - q_aux_gpu(b,i+1,j,k-1,4)+ &
                         q_aux_gpu(b,i,j,k+1,4)   - q_aux_gpu(b,i,j,k-1,4) )*0.25_R8P/dz_gpu(b)

                vel_u = 0.5*(q_aux_gpu(b,i,j,k,2) + q_aux_gpu(b,i+1,j,k,2))
                vel_v = 0.5*(q_aux_gpu(b,i,j,k,3) + q_aux_gpu(b,i+1,j,k,3))
                vel_w = 0.5*(q_aux_gpu(b,i,j,k,4) + q_aux_gpu(b,i+1,j,k,4))

                tau_1_1 = 2.0*mu*(du_dx-1./3.*(du_dx+dv_dy+dw_dz))
                tau_2_1 = mu*(dv_dx+du_dy)
                tau_3_1 = mu*(dw_dx+du_dz)

                dT_dx = (q_aux_gpu(b,i+1,j,k,6)-q_aux_gpu(b,i,j,k,6))/dx_gpu(b)

                sigq = k_coeff*dT_dx
                sigl = vel_u*tau_1_1+vel_v*tau_2_1+vel_w*tau_3_1

                flx_gpu(b,i,j,k,2) = flx_gpu(b,i,j,k,2) - tau_1_1
                flx_gpu(b,i,j,k,3) = flx_gpu(b,i,j,k,3) - tau_2_1
                flx_gpu(b,i,j,k,4) = flx_gpu(b,i,j,k,4) - tau_3_1
                flx_gpu(b,i,j,k,5) = flx_gpu(b,i,j,k,5) - sigq + sigl
            enddo

         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()

   !$cuf kernel do(3) <<<*,*>>>
   do k=1,nk
      do i=1,ni
         do b=1,blocks_number
            do j=0,nj ! loop on faces
                du_dy = (q_aux_gpu(b,i,j+1,k,2)-q_aux_gpu(b,i,j,k,2))/dy_gpu(b)
                dv_dy = (q_aux_gpu(b,i,j+1,k,3)-q_aux_gpu(b,i,j,k,3))/dy_gpu(b)
                dw_dy = (q_aux_gpu(b,i,j+1,k,4)-q_aux_gpu(b,i,j,k,4))/dy_gpu(b)

                du_dx = (q_aux_gpu(b,i+1,j+1,k,2) - q_aux_gpu(b,i-1,j+1,k,2)+ &
                         q_aux_gpu(b,i+1,j,k,2)   - q_aux_gpu(b,i-1,j,k,2) )*0.25_R8P/dx_gpu(b)
                dv_dx = (q_aux_gpu(b,i+1,j+1,k,3) - q_aux_gpu(b,i-1,j+1,k,3)+ &
                         q_aux_gpu(b,i+1,j,k,3)   - q_aux_gpu(b,i-1,j,k,3) )*0.25_R8P/dx_gpu(b)

                dv_dz = (q_aux_gpu(b,i+1,j,k+1,3) - q_aux_gpu(b,i-1,j,k+1,3)+ &
                         q_aux_gpu(b,i+1,j,k,3)   - q_aux_gpu(b,i-1,j,k,3) )*0.25_R8P/dz_gpu(b)
                dw_dz = (q_aux_gpu(b,i+1,j,k+1,4) - q_aux_gpu(b,i-1,j,k+1,4)+ &
                         q_aux_gpu(b,i+1,j,k,4)   - q_aux_gpu(b,i-1,j,k,4) )*0.25_R8P/dz_gpu(b)

                vel_u = 0.5*(q_aux_gpu(b,i,j,k,2) + q_aux_gpu(b,i,j+1,k,2))
                vel_v = 0.5*(q_aux_gpu(b,i,j,k,3) + q_aux_gpu(b,i,j+1,k,3))
                vel_w = 0.5*(q_aux_gpu(b,i,j,k,4) + q_aux_gpu(b,i,j+1,k,4))

                tau_1_2 = mu*(du_dy+dv_dx)
                tau_2_2 = 2.0*mu*(dv_dy-1./3.*(du_dx+dv_dy+dw_dz))
                tau_3_2 = mu*(dw_dy+dv_dz)

                dT_dy = (q_aux_gpu(b,i,j+1,k,6)-q_aux_gpu(b,i,j,k,6))/dy_gpu(b)

                sigq = k_coeff*dT_dy
                sigl = vel_u*tau_1_2+vel_v*tau_2_2+vel_w*tau_3_2

                fly_gpu(b,i,j,k,2) = fly_gpu(b,i,j,k,2) - tau_1_2
                fly_gpu(b,i,j,k,3) = fly_gpu(b,i,j,k,3) - tau_2_2
                fly_gpu(b,i,j,k,4) = fly_gpu(b,i,j,k,4) - tau_3_2
                fly_gpu(b,i,j,k,5) = fly_gpu(b,i,j,k,5) - sigq + sigl
            enddo

         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()

   !$cuf kernel do(3) <<<*,*>>>
   do j=1,nj
      do i=1,ni
         do b=1,blocks_number
            do k=0,nk ! loop on faces
                du_dz = (q_aux_gpu(b,i,j,k+1,2)-q_aux_gpu(b,i,j,k,2))/dz_gpu(b)
                dv_dz = (q_aux_gpu(b,i,j,k+1,3)-q_aux_gpu(b,i,j,k,3))/dz_gpu(b)
                dw_dz = (q_aux_gpu(b,i,j,k+1,4)-q_aux_gpu(b,i,j,k,4))/dz_gpu(b)

                du_dx = (q_aux_gpu(b,i+1,j,k+1,2) - q_aux_gpu(b,i-1,j,k+1,2)+ &
                         q_aux_gpu(b,i+1,j,k,2)   - q_aux_gpu(b,i-1,j,k,2) )*0.25_R8P/dx_gpu(b)
                dw_dx = (q_aux_gpu(b,i+1,j,k+1,4) - q_aux_gpu(b,i-1,j,k+1,4)+ &
                         q_aux_gpu(b,i+1,j,k,4)   - q_aux_gpu(b,i-1,j,k,4) )*0.25_R8P/dx_gpu(b)

                dv_dy = (q_aux_gpu(b,i,j+1,k+1,3) - q_aux_gpu(b,i,j-1,k+1,3)+ &
                         q_aux_gpu(b,i,j+1,k,3)   - q_aux_gpu(b,i,j-1,k,3) )*0.25_R8P/dy_gpu(b)
                dw_dy = (q_aux_gpu(b,i,j+1,k+1,4) - q_aux_gpu(b,i,j-1,k+1,4)+ &
                         q_aux_gpu(b,i,j+1,k,4)   - q_aux_gpu(b,i,j-1,k,4) )*0.25_R8P/dy_gpu(b)

                vel_u = 0.5*(q_aux_gpu(b,i,j,k,2) + q_aux_gpu(b,i,j,k+1,2))
                vel_v = 0.5*(q_aux_gpu(b,i,j,k,3) + q_aux_gpu(b,i,j,k+1,3))
                vel_w = 0.5*(q_aux_gpu(b,i,j,k,4) + q_aux_gpu(b,i,j,k+1,4))

                tau_1_3 = mu*(du_dz+dw_dx)
                tau_2_3 = mu*(dv_dz+dw_dy)
                tau_3_3 = 2.0*mu*(dw_dz-1./3.*(du_dx+dv_dy+dw_dz))

                dT_dz = (q_aux_gpu(b,i,j,k+1,6)-q_aux_gpu(b,i,j,k,6))/dz_gpu(b)

                sigq = k_coeff*dT_dz
                sigl = vel_u*tau_1_3+vel_v*tau_2_3+vel_w*tau_3_3

                flz_gpu(b,i,j,k,2) = flz_gpu(b,i,j,k,2) - tau_1_3
                flz_gpu(b,i,j,k,3) = flz_gpu(b,i,j,k,3) - tau_2_3
                flz_gpu(b,i,j,k,4) = flz_gpu(b,i,j,k,4) - tau_3_3
                flz_gpu(b,i,j,k,5) = flz_gpu(b,i,j,k,5) - sigq + sigl
            enddo

         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()

   endsubroutine viscous_part

   attributes(global) subroutine euler_y_kernel(q_aux_gpu, fly_gpu, gplus, gminus, &
                                                blocks_number, ni, nj, nk, ngc, nv, iweno, dha_star, gamma_fluid, R_star, cv_star)

   real(R8P), intent(in), device     :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::   fly_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::     gplus(1:, 1:, 1:, 1:, 1:)
   real(R8P), intent(inout), device  ::    gminus(1:, 1:, 1:, 1:, 1:)
   integer, intent(in), value        :: blocks_number, ni, nj, nk, ngc, nv, iweno
   real(R8P), intent(in), value      :: dha_star, gamma_fluid, R_star, cv_star
   integer                           :: b, i, j, k, l, ll, m, mm, v
   real(R8P)                         :: er(5,5), el(5,5), ev(5), evmax(5), ghat(5), gl(5), gr(5), fi(5), vi(5)
   !real(R8P)                         :: er(nv,nv), el(nv,nv), ev(nv), evmax(nv), ghat(nv), gl(nv), gr(nv), fi(nv), vi(nv)
   real(R8P)                         :: uu, vv, ww, h, ya, qq, c, ci, b1, b2
   real(R8P)                         :: gc, wc

   b = blockDim%x * (blockIdx%x - 1) + threadIdx%x
   i = blockDim%y * (blockIdx%y - 1) + threadIdx%y
   if(b > blocks_number .or. i > ni) return

   do k=1,nk

      do j=0,nj ! loop on faces

         ! Compute Roe average
         call compute_roe_average(q_aux_gpu, dha_star, gamma_fluid, R_star, &
            ngc, b, i, j, k, i, j+1, k, uu, vv, ww, h, ya, qq, c, ci, b1, b2)

         ! Compute right and left eigenvectors matrices (at Roe state)
         er(1,1) = 1._R8P ;  er(1,2) = uu      ; er(1,3) = vv-c   ; er(1,4) = ww     ; er(1,5) = h-vv*c
         er(2,1) = 1._R8P ;  er(2,2) = uu      ; er(2,3) = vv     ; er(2,4) = ww     ; er(2,5) = qq
         er(3,1) = 1._R8P ;  er(3,2) = uu      ; er(3,3) = vv+c   ; er(3,4) = ww     ; er(3,5) = h+vv*c
         er(4,1) = 0._R8P ;  er(4,2) = 1._R8P  ; er(4,3) = 0._R8P ; er(4,4) = 0._R8P ; er(4,5) = 0._R8P
         er(5,1) = 0._R8P ;  er(5,2) = 0._R8P  ; er(5,3) = 0._R8P ; er(5,4) = 1._R8P ; er(5,5) = 0._R8P

         el(1,1) =  0.5_R8P*(b1+vv*ci) ; el(1,2) = 1._R8P-b1 ; el(1,3) = 0.5_R8P*(b1-vv*ci)
         el(2,1) = -0.5_R8P*(b2*uu)    ; el(2,2) = b2*uu     ; el(2,3) = -0.5_R8P*(b2*uu)
         el(3,1) = -0.5_R8P*(b2*vv+ci) ; el(3,2) = b2*vv     ; el(3,3) = -0.5_R8P*(b2*vv-ci)
         el(4,1) = -0.5_R8P*(b2*ww)    ; el(4,2) = b2*ww     ; el(4,3) = -0.5_R8P*(b2*ww)
         el(5,1) =  0.5_R8P*b2         ; el(5,2) = -b2       ; el(5,3) = 0.5_R8P*b2

         el(1,4) = -ww                 ; el(1,5) = uu
         el(2,4) = 0._R8P              ; el(2,5) = -1._R8P
         el(3,4) = 0._R8P              ; el(3,5) = 0._R8P
         el(4,4) = 1._R8P              ; el(4,5) = 0._R8P
         el(5,4) = 0._R8P              ; el(5,5) = 0._R8P

         if(nv == 6) then
         ! Compute right and left eigenvectors matrices (at Roe state)
                                                                                                          er(1,6) = ya
                                                                                                          er(2,6) = 0._R8P
                                                                                                          er(3,6) = ya
                                                                                                          er(4,6) = -uu/dha_star
                                                                                                          er(5,6) = -ww/dha_star
         er(6,1) = 0._R8P ;  er(6,2) = 0._R8P  ; er(6,3) = 0._R8P ; er(6,4) = 0._R8P ; er(6,5) = 1._R8P ; er(6,6) = 1._R8P/dha_star

         el(6,1) = -0.5_R8P*dha_star*b2     ; el(6,2) = dha_star*b2    ; el(6,3) = -0.5_R8P*dha_star*b2

                                                               el(1,6) = -ya*dha_star*b1-uu**2-ww**2
                                                               el(2,6) = uu*(1._R8P+ya*dha_star*b2)
                                                               el(3,6) = ya*dha_star*vv*b2
                                                               el(4,6) = ww*(1._R8P+ya*dha_star*b2)
                                                               el(5,6) = -ya*dha_star*b2
         el(6,4) = 0._R8P              ; el(6,5) = 0._R8P    ; el(6,6) = dha_star*(1._R8P+ya*dha_star*b2)
         endif

         ! Find max eigenvalues on the stencil
         do m=1,nv  ! loop on characteristic fields
            evmax(m) = -1._R8P
         enddo
         do l=1,2*iweno ! LLF
            ll = j + l - iweno
            uu = q_aux_gpu(b,i,ll,k,3)
            c  = q_aux_gpu(b,i,ll,k,9)
            ev(1) = abs(uu-c) ; ev(2) = abs(uu) ; ev(3) = abs(uu+c) ; ev(4) = ev(2) ; ev(5) = ev(2) ;
            if(nv == 6) ev(6) = ev(2)
            do m=1,nv
                evmax(m) = max(ev(m),evmax(m))
            enddo
         enddo

         ! Decompose fluxes as + and -
         do l=1,2*iweno ! loop over the stencil centered at face i
            ll = j + l - iweno
            vi(1) = q_aux_gpu(b,i,ll,k,1)
            vi(2) = vi(1)*q_aux_gpu(b,i,ll,k,2)
            vi(3) = vi(1)*q_aux_gpu(b,i,ll,k,3)
            vi(4) = vi(1)*q_aux_gpu(b,i,ll,k,4)
            vi(5) = vi(1)*(cv_star*q_aux_gpu(b,i,ll,k,6)+&
                0.5_R8P*(q_aux_gpu(b,i,ll,k,2)**2+q_aux_gpu(b,i,ll,k,3)**2+q_aux_gpu(b,i,ll,k,4)**2)+&
                q_aux_gpu(b,i,ll,k,5)*dha_star)
            if(nv == 6) vi(6) = vi(1)*q_aux_gpu(b,i,ll,k,5)
            fi(1) = vi(3)
            fi(2) = fi(1) * q_aux_gpu(b,i,ll,k,2)
            fi(3) = fi(1) * q_aux_gpu(b,i,ll,k,3) + q_aux_gpu(b,i,ll,k,7)
            fi(4) = fi(1) * q_aux_gpu(b,i,ll,k,4)
            fi(5) = fi(1) * vi(5) / vi(1) + q_aux_gpu(b,i,ll,k,7)*q_aux_gpu(b,i,ll,k,3)
            if(nv == 6) fi(6) = fi(1) * q_aux_gpu(b,i,ll,k,5)
            do m=1,nv
               wc = 0._R8P
               gc = 0._R8P
               do mm=1,nv
                  wc = wc + el(mm,m) * vi(mm)
                  gc = gc + el(mm,m) * fi(mm)
               enddo
               gplus (m,l,i,k,b) = 0.5_R8P * (gc + evmax(m) * wc)
               gminus(m,l,i,k,b) = gc - gplus(m,l,i,k,b)
            enddo
         enddo

         ! Reconstruction of the + and - fluxes
         call weno_reconstruction(nv, gplus(1,1,i,k,b), gminus(1,1,i,k,b), gl, gr, iweno)

         ! Reassemble + and - characteristic fluxes
         do m=1,nv
            ghat(m) = gl(m) + gr(m)
         enddo

         ! Return to conservative fluxes
         do m=1,nv
            fly_gpu(b,i,j,k,m) = 0._R8P
            do mm=1,nv
               fly_gpu(b,i,j,k,m) = fly_gpu(b,i,j,k,m) + er(mm,m) * ghat(mm)
            enddo
         enddo

      enddo
   enddo
   endsubroutine euler_y_kernel

   attributes(global) subroutine euler_z_kernel(q_aux_gpu, flz_gpu, gplus, gminus, &
                                                blocks_number, ni, nj, nk, ngc, nv, iweno, dha_star, gamma_fluid, R_star, cv_star)

   real(R8P), intent(in), device     :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::   flz_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::     gplus(1:, 1:, 1:, 1:, 1:)
   real(R8P), intent(inout), device  ::    gminus(1:, 1:, 1:, 1:, 1:)
   integer, intent(in), value        :: blocks_number, ni, nj, nk, ngc, nv, iweno
   real(R8P), intent(in), value      :: dha_star, gamma_fluid, R_star, cv_star
   integer                           :: b, i, j, k, l, ll, m, mm, v
   real(R8P)                         :: er(5,5), el(5,5), ev(5), evmax(5), ghat(5), gl(5), gr(5), fi(5), vi(5)
   !real(R8P)                         :: er(nv,nv), el(nv,nv), ev(nv), evmax(nv), ghat(nv), gl(nv), gr(nv), fi(nv), vi(nv)
   real(R8P)                         :: uu, vv, ww, h, ya, qq, c, ci, b1, b2
   real(R8P)                         :: gc, wc

   b = blockDim%x * (blockIdx%x - 1) + threadIdx%x
   i = blockDim%y * (blockIdx%y - 1) + threadIdx%y
   if(b > blocks_number .or. i > ni) return

   do j=1,nj

      do k=0,nk ! loop on faces

         ! Compute Roe average
         call compute_roe_average(q_aux_gpu, dha_star, gamma_fluid, R_star, &
            ngc, b, i, j, k, i, j, k+1, uu, vv, ww, h, ya, qq, c, ci, b1, b2)

         ! Compute right and left eigenvectors matrices (at Roe state)
         er(1,1) = 1._R8P ;  er(1,2) = uu      ; er(1,3) = vv     ; er(1,4) = ww-c    ; er(1,5) = h-ww*c
         er(2,1) = 1._R8P ;  er(2,2) = uu      ; er(2,3) = vv     ; er(2,4) = ww      ; er(2,5) = qq
         er(3,1) = 1._R8P ;  er(3,2) = uu      ; er(3,3) = vv     ; er(3,4) = ww+c    ; er(3,5) = h+ww*c
         er(4,1) = 0._R8P ;  er(4,2) = 1._R8P  ; er(4,3) = 0._R8P ; er(4,4) = 0._R8P  ; er(4,5) = 0._R8P
         er(5,1) = 0._R8P ;  er(5,2) = 0._R8P  ; er(5,3) = 1._R8P ; er(5,4) = 0._R8P  ; er(5,5) = 0._R8P

         el(1,1) = 0.5_R8P*(b1+ww*ci)          ; el(1,2) = 1._R8P-b1 ; el(1,3) = 0.5_R8P*(b1-ww*ci)
         el(2,1) = -0.5_R8P*(b2*uu)            ; el(2,2) = b2*uu     ; el(2,3) = -0.5_R8P*(b2*uu)
         el(3,1) = -0.5_R8P*(b2*vv)            ; el(3,2) = b2*vv     ; el(3,3) = -0.5_R8P*(b2*vv)
         el(4,1) = -0.5_R8P*(b2*ww+ci)         ; el(4,2) = b2*ww     ; el(4,3) = -0.5_R8P*(b2*ww-ci)
         el(5,1) = 0.5_R8P*b2                  ; el(5,2) = -b2       ; el(5,3) = 0.5_R8P*b2

         el(1,4) = -uu                         ; el(1,5) = -vv
         el(2,4) = 1._R8P                      ; el(2,5) = 0._R8P
         el(3,4) = 0._R8P                      ; el(3,5) = 1._R8P
         el(4,4) = 0._R8P                      ; el(4,5) = 0._R8P
         el(5,4) = 0._R8P                      ; el(5,5) = 0._R8P

         if(nv == 6) then
         ! Compute right and left eigenvectors matrices (at Roe state)
                                                                                                           er(1,6) = ya
                                                                                                           er(2,6) = 0._R8P
                                                                                                           er(3,6) = ya
                                                                                                           er(4,6) = -uu/dha_star
                                                                                                           er(5,6) = -vv/dha_star
         er(6,1) = 0._R8P ;  er(6,2) = 0._R8P  ; er(6,3) = 0._R8P ; er(6,4) = 0._R8P  ; er(6,5) = 1._R8P ; er(6,6) = 1._R8P/dha_star

         el(6,1) = -0.5_R8P*dha_star*b2             ; el(6,2) = dha_star*b2    ; el(6,3) = -0.5_R8P*dha_star*b2

                                                                       el(1,6) = -ya*dha_star*b1-uu**2-vv**2
                                                                       el(2,6) = uu*(1+ya*dha_star*b2)
                                                                       el(3,6) = vv*(1+ya*dha_star*b2)
                                                                       el(4,6) = ya*dha_star*ww*b2
                                                                       el(5,6) = -ya*dha_star*b2
         el(6,4) = 0._R8P                      ; el(6,5) = 0._R8P    ; el(6,6) = dha_star*(1._R8P+ya*dha_star*b2)
         endif

         ! Find max eigenvalues on the stencil
         do m=1,nv  ! loop on characteristic fields
            evmax(m) = -1._R8P
         enddo
         do l=1,2*iweno ! LLF
            ll = k + l - iweno
            uu = q_aux_gpu(b,i,j,ll,4)
            c  = q_aux_gpu(b,i,j,ll,9)
            ev(1) = abs(uu-c) ; ev(2) = abs(uu) ; ev(3) = abs(uu+c) ; ev(4) = ev(2) ; ev(5) = ev(2)
            if(nv == 6) ev(6) = ev(2)
            do m=1,nv
                evmax(m) = max(ev(m),evmax(m))
            enddo
         enddo

         ! Decompose fluxes as + and -
         do l=1,2*iweno ! loop over the stencil centered at face i
            ll = k + l - iweno
            vi(1) = q_aux_gpu(b,i,j,ll,1)
            vi(2) = vi(1)*q_aux_gpu(b,i,j,ll,2)
            vi(3) = vi(1)*q_aux_gpu(b,i,j,ll,3)
            vi(4) = vi(1)*q_aux_gpu(b,i,j,ll,4)
            vi(5) = vi(1)*(cv_star*q_aux_gpu(b,i,j,ll,6)+&
                0.5_R8P*(q_aux_gpu(b,i,j,ll,2)**2+q_aux_gpu(b,i,j,ll,3)**2+q_aux_gpu(b,i,j,ll,4)**2)+&
                q_aux_gpu(b,i,j,ll,5)*dha_star)
            if(nv == 6) vi(6) = vi(1)*q_aux_gpu(b,i,j,ll,5)
            fi(1) = vi(4)
            fi(2) = fi(1) * q_aux_gpu(b,i,j,ll,2)
            fi(3) = fi(1) * q_aux_gpu(b,i,j,ll,3)
            fi(4) = fi(1) * q_aux_gpu(b,i,j,ll,4) + q_aux_gpu(b,i,j,ll,7)
            fi(5) = fi(1) * vi(5) / vi(1) + q_aux_gpu(b,i,j,ll,7)*q_aux_gpu(b,i,j,ll,4)
            if(nv == 6) fi(6) = fi(1) * q_aux_gpu(b,i,j,ll,5)
            do m=1,nv
               wc = 0._R8P
               gc = 0._R8P
               do mm=1,nv
                  wc = wc + el(mm,m) * vi(mm)
                  gc = gc + el(mm,m) * fi(mm)
               enddo
               gplus (m,l,i,j,b) = 0.5_R8P * (gc + evmax(m) * wc)
               gminus(m,l,i,j,b) = gc - gplus(m,l,i,j,b)
            enddo
         enddo

         ! Reconstruction of the + and - fluxes
         call weno_reconstruction(nv, gplus(1,1,i,j,b), gminus(1,1,i,j,b), gl, gr, iweno)

         ! Reassemble + and - characteristic fluxes
         do m=1,nv
            ghat(m) = gl(m) + gr(m)
         enddo

         ! Return to conservative fluxes
         do m=1,nv
            flz_gpu(b,i,j,k,m) = 0._R8P
            do mm=1,nv
               flz_gpu(b,i,j,k,m) = flz_gpu(b,i,j,k,m) + er(mm,m) * ghat(mm)
            enddo
         enddo

      enddo
   enddo
   endsubroutine euler_z_kernel

   attributes(device) subroutine compute_roe_average(q_aux_gpu, dha_star, gamma_fluid, R_star, &
      ngc, b, i, j, k, ip, jp, kp, uu, vv, ww, h, ya, qq, c, ci, b1, b2)

   implicit none
   real(R8P), intent(in), device  :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(in)          :: dha_star, gamma_fluid, R_star
   integer(I4P), intent(in)       :: ngc, b, i, j, k, ip, jp, kp
   real(R8P), intent(out)         :: uu, vv, ww, h, ya, qq, c, ci, b1, b2
   real(R8P)                      :: ri, up, vp, wp, hp, yap, r, rp1, cc
   ! Left state (node i)
   ri        =  1._R8P/q_aux_gpu(b,i,j,k,1)
   uu        =  q_aux_gpu(b,i,j,k,2)
   vv        =  q_aux_gpu(b,i,j,k,3)
   ww        =  q_aux_gpu(b,i,j,k,4)
   h         =  q_aux_gpu(b,i,j,k,8)
   ya        =  q_aux_gpu(b,i,j,k,5)
   ! Right state (node i+1)
   up        =  q_aux_gpu(b,ip,jp,kp,2)
   vp        =  q_aux_gpu(b,ip,jp,kp,3)
   wp        =  q_aux_gpu(b,ip,jp,kp,4)
   hp        =  q_aux_gpu(b,ip,jp,kp,8)
   yap       =  q_aux_gpu(b,ip,jp,kp,5)
   ! Average state
   r         =  sqrt(q_aux_gpu(b,ip,jp,kp,1)*ri)
   rp1       =  1._R8P/(r+1._R8P)
   uu        =  (r*up+uu)*rp1
   vv        =  (r*vp+vv)*rp1
   ww        =  (r*wp+ww)*rp1
   h         =  (r*hp+h)*rp1
   ya        =  (r*yap+ya)*rp1
   qq        =  0.5_R8P * (uu*uu+vv*vv+ww*ww)
   cc        =  (gamma_fluid-1._R8P) * (h - qq - ya*dha_star)
   !ERRATODIREIcc        =  gamma_fluid * (gamma_fluid-1._R8P) * (h - qq - ya*dha_star)
   c         =  sqrt(cc)
   ci        =  1._R8P/c
   b2        = (gamma_fluid-1)/cc  ! alias 1/(cp*theta)
   b1        = b2 * qq             ! alias q/(cp*theta)

   endsubroutine compute_roe_average

   attributes(device) subroutine weno_reconstruction(nvar,vp,vm,vminus,vplus,iweno)

   implicit none
   integer, intent(in)                             :: nvar, iweno
   !real(R8P), dimension(nvar,2*iweno), intent(in)  :: vm,vp
   real(R8P), dimension(1:nvar,1:*) :: vm,vp
   real(R8P), dimension(nvar), intent(out)         :: vminus,vplus

   real(R8P), dimension(-1:4) :: dwe               ! linear weights
   real(R8P), dimension(-1:4) :: alfp,alfm         ! alpha_l
   real(R8P), dimension(-1:4) :: alfp_map,alfm_map ! alpha_l
   real(R8P), dimension(-1:4) :: betap,betam       ! beta_l
   real(R8P), dimension(-1:4) :: omp,omm           ! WENO weights
   integer                    :: r,i,j,k,l,m
   real(R8P)                  :: c0,c1,c2,c3,c4,d0,d1,d2,d3,d4,summ,sump
   real(R8P)                  :: x,y,y2

   if (iweno==1) then ! Godunov

       i = iweno ! index of intermediate node to perform reconstruction

       vminus(1:nvar) = vp(1:nvar,i)
       vplus (1:nvar) = vm(1:nvar,i+1)

   elseif (iweno==2) then ! WENO-3

       i = iweno ! index of intermediate node to perform reconstruction

       dwe(1)   = 2._R8P/3._R8P
       dwe(0)   = 1._R8P/3._R8P

       do m=1,nvar

           betap(0)  = (vp(m,i  )-vp(m,i-1))**2
           betap(1)  = (vp(m,i+1)-vp(m,i  ))**2
           betam(0)  = (vm(m,i+2)-vm(m,i+1))**2
           betam(1)  = (vm(m,i+1)-vm(m,i  ))**2

           sump = 0._R8P
           summ = 0._R8P
           do l=0,1
               alfp(l) = dwe(l)/(0.000001_R8P+betap(l))**2
               alfm(l) = dwe(l)/(0.000001_R8P+betam(l))**2
               sump = sump + alfp(l)
               summ = summ + alfm(l)
           enddo
           do l=0,1
               omp(l) = alfp(l)/sump
               omm(l) = alfm(l)/summ
           enddo

           vminus(m) = omp(0) *(-vp(m,i-1)+3*vp(m,i  )) + omp(1) *( vp(m,i  )+ vp(m,i+1))
           vplus(m)  = omm(0) *(-vm(m,i+2)+3*vm(m,i+1)) + omm(1) *( vm(m,i  )+ vm(m,i+1))

       enddo

       do m=1,nvar
           vminus(m) = 0.5_R8P*vminus(m)
           vplus(m)  = 0.5_R8P*vplus(m)
       enddo

     elseif (iweno==3) then ! WENO-5
!
      i = iweno ! index of intermediate node to perform reconstruction
!
      dwe( 0) = 1._R8P/10._R8P
      dwe( 1) = 6._R8P/10._R8P
      dwe( 2) = 3._R8P/10._R8P
!     JS
      d0 = 13._R8P/12._R8P
      d1 = 1._R8P/4._R8P
!     Weights for polynomial reconstructions
      c0 = 1._R8P/3._R8P
      c1 = 5._R8P/6._R8P
      c2 =-1._R8P/6._R8P
      c3 =-7._R8P/6._R8P
      c4 =11._R8P/6._R8P
!
      do m=1,nvar
!
       betap(2) = d0*(     vp(m,i)-2._R8P*vp(m,i+1)+vp(m,i+2))**2+d1*(3._R8P*vp(m,i)-4._R8P*vp(m,i+1)+vp(m,i+2))**2
       betap(1) = d0*(     vp(m,i-1)-2._R8P*vp(m,i)+vp(m,i+1))**2+d1*(     vp(m,i-1)-vp(m,i+1) )**2
       betap(0) = d0*(     vp(m,i)-2._R8P*vp(m,i-1)+vp(m,i-2))**2+d1*(3._R8P*vp(m,i)-4._R8P*vp(m,i-1)+vp(m,i-2))**2
!
       betam(2) = d0*(     vm(m,i+1)-2._R8P*vm(m,i)+vm(m,i-1))**2+d1*(3._R8P*vm(m,i+1)-4._R8P*vm(m,i)+vm(m,i-1))**2
       betam(1) = d0*(     vm(m,i+2)-2._R8P*vm(m,i+1)+vm(m,i))**2+d1*(     vm(m,i+2)-vm(m,i) )**2
       betam(0) = d0*(     vm(m,i+1)-2._R8P*vm(m,i+2)+vm(m,i+3))**2+d1*(3._R8P*vm(m,i+1)-4._R8P*vm(m,i+2)+vm(m,i+3))**2
!
       sump = 0._R8P
       summ = 0._R8P
       do l=0,2
        alfp(l) = dwe(  l)/(0.000001_R8P+betap(l))**2
        alfm(l) = dwe(  l)/(0.000001_R8P+betam(l))**2
        sump = sump + alfp(l)
        summ = summ + alfm(l)
       enddo
       do l=0,2
        omp(l) = alfp(l)/sump
        omm(l) = alfm(l)/summ
       enddo
!
       vminus(m)   = omp(2)*(c0*vp(m,i  )+c1*vp(m,i+1)+c2*vp(m,i+2)) + &
         & omp(1)*(c2*vp(m,i-1)+c1*vp(m,i  )+c0*vp(m,i+1)) + omp(0)*(c0*vp(m,i-2)+c3*vp(m,i-1)+c4*vp(m,i  ))
       vplus(m)   = omm(2)*(c0*vm(m,i+1)+c1*vm(m,i  )+c2*vm(m,i-1)) +  &
         & omm(1)*(c2*vm(m,i+2)+c1*vm(m,i+1)+c0*vm(m,i  )) + omm(0)*(c0*vm(m,i+3)+c3*vm(m,i+2)+c4*vm(m,i+1))
!
      enddo ! end of m-loop
!
   elseif (iweno==4) then ! WENO-7
!
      i = iweno ! index of intermediate node to perform reconstruction
!
      dwe( 0) = 1._R8P/35._R8P
      dwe( 1) = 12._R8P/35._R8P
      dwe( 2) = 18._R8P/35._R8P
      dwe( 3) = 4._R8P/35._R8P
!
!     JS weights
      d1 = 1._R8P/36._R8P
      d2 = 13._R8P/12._R8P
      d3 = 781._R8P/720._R8P
!
      do m=1,nvar
!
       betap(3)= d1*(-11*vp(m,  i)+18*vp(m,i+1)- 9*vp(m,i+2)+ 2*vp(m,i+3))**2+&
       &  d2*(  2*vp(m,  i)- 5*vp(m,i+1)+ 4*vp(m,i+2)-   vp(m,i+3))**2+ &
       & d3*(   -vp(m,  i)+ 3*vp(m,i+1)- 3*vp(m,i+2)+   vp(m,i+3))**2
       betap(2)= d1*(- 2*vp(m,i-1)- 3*vp(m,i  )+ 6*vp(m,i+1)-   vp(m,i+2))**2+&
       &  d2*(    vp(m,i-1)- 2*vp(m,i  )+   vp(m,i+1)             )**2+&
       &  d3*(   -vp(m,i-1)+ 3*vp(m,i  )- 3*vp(m,i+1)+   vp(m,i+2))**2
       betap(1)= d1*(    vp(m,i-2)- 6*vp(m,i-1)+ 3*vp(m,i  )+ 2*vp(m,i+1))**2+&
       &  d2*( vp(m,i-1)- 2*vp(m,i  )+   vp(m,i+1))**2+ &
       &  d3*(   -vp(m,i-2)+ 3*vp(m,i-1)- 3*vp(m,i  )+   vp(m,i+1))**2
       betap(0)= d1*(- 2*vp(m,i-3)+ 9*vp(m,i-2)-18*vp(m,i-1)+11*vp(m,i  ))**2+&
       &  d2*(-   vp(m,i-3)+ 4*vp(m,i-2)- 5*vp(m,i-1)+ 2*vp(m,i  ))**2+&
       &  d3*(   -vp(m,i-3)+ 3*vp(m,i-2)- 3*vp(m,i-1)+   vp(m,i  ))**2
!
       betam(3)= d1*(-11*vm(m,i+1)+18*vm(m,i  )- 9*vm(m,i-1)+ 2*vm(m,i-2))**2+&
       &  d2*(  2*vm(m,i+1)- 5*vm(m,i  )+ 4*vm(m,i-1)-   vm(m,i-2))**2+&
       &  d3*(   -vm(m,i+1)+ 3*vm(m,i  )- 3*vm(m,i-1)+   vm(m,i-2))**2
       betam(2)= d1*(- 2*vm(m,i+2)- 3*vm(m,i+1)+ 6*vm(m,i  )-   vm(m,i-1))**2+&
       &  d2*(    vm(m,i+2)- 2*vm(m,i+1)+   vm(m,i  )             )**2+&
       &  d3*(   -vm(m,i+2)+ 3*vm(m,i+1)- 3*vm(m,i  )+   vm(m,i-1))**2
       betam(1)= d1*(    vm(m,i+3)- 6*vm(m,i+2)+ 3*vm(m,i+1)+ 2*vm(m,i  ))**2+&
       &  d2*(                 vm(m,i+2)- 2*vm(m,i+1)+   vm(m,i  ))**2+&
       &  d3*(   -vm(m,i+3)+ 3*vm(m,i+2)- 3*vm(m,i+1)+   vm(m,i  ))**2
       betam(0)= d1*(- 2*vm(m,i+4)+ 9*vm(m,i+3)-18*vm(m,i+2)+11*vm(m,i+1))**2+&
       &  d2*(-   vm(m,i+4)+ 4*vm(m,i+3)- 5*vm(m,i+2)+ 2*vm(m,i+1))**2+&
       &  d3*(   -vm(m,i+4)+ 3*vm(m,i+3)- 3*vm(m,i+2)+   vm(m,i+1))**2
!
       sump = 0._R8P
       summ = 0._R8P
       do l=0,3
        alfp(l) = dwe(  l)/(0.000001_R8P+betap(l))**2
        alfm(l) = dwe(  l)/(0.000001_R8P+betam(l))**2
        sump = sump + alfp(l)
        summ = summ + alfm(l)
       enddo
       do l=0,3
        omp(l) = alfp(l)/sump
        omm(l) = alfm(l)/summ
       enddo
!
       vminus(m)   = omp(3)*( 6*vp(m,i  )+26*vp(m,i+1)-10*vp(m,i+2)+ 2*vp(m,i+3))+&
        omp(2)*(-2*vp(m,i-1)+14*vp(m,i  )+14*vp(m,i+1)- 2*vp(m,i+2))+&
        omp(1)*( 2*vp(m,i-2)-10*vp(m,i-1)+26*vp(m,i  )+ 6*vp(m,i+1))+&
        omp(0)*(-6*vp(m,i-3)+26*vp(m,i-2)-46*vp(m,i-1)+50*vp(m,i  ))
       vplus(m)   =  omm(3)*( 6*vm(m,i+1)+26*vm(m,i  )-10*vm(m,i-1)+ 2*vm(m,i-2))+&
        omm(2)*(-2*vm(m,i+2)+14*vm(m,i+1)+14*vm(m,i  )- 2*vm(m,i-1))+&
        omm(1)*( 2*vm(m,i+3)-10*vm(m,i+2)+26*vm(m,i+1)+ 6*vm(m,i  ))+&
        omm(0)*(-6*vm(m,i+4)+26*vm(m,i+3)-46*vm(m,i+2)+50*vm(m,i+1))
!
      enddo ! end of m-loop
!
      vminus = vminus/24._R8P
      vplus  = vplus /24._R8P
!
   else
      write(*,*) 'Error! WENO scheme not implemented'
      stop
   endif

   endsubroutine weno_reconstruction

   subroutine compute_rk_linear_gpu_cuf(ni, nj, nk, ngc, nv, blocks_number, dt, q_gpu, prhs_gpu, fl_gpu, phi_gpu, qnrk)

   !< Initialize RK stage with q_gpu.
   integer(I4P), intent(in)            :: ni                                     !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj                                     !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk                                     !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: ngc                                    !< Ghost cells number.
   integer(I4P), intent(in)            :: nv                                     !< Number of conservative varibales.
   integer(I4P), intent(in)            :: blocks_number                          !< Number of blocks.
   real(R8P),    intent(in)            :: dt                                     !< Time step.
   real(R8P),    intent(in)            :: qnrk                                   !< Time step.
   real(R8P),    intent(inout), device ::   q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Conservative field.
   real(R8P),    intent(in),    device ::  fl_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Conservative field.
   real(R8P),    intent(in),    device :: phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Conservative field.
   real(R8P),    intent(inout), device :: prhs_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< RK stage.
   integer(I4P)                        :: i, j, k, b, v, ss                      !< Counter.
   integer(I4P)                        :: iercuda                                !< Error trapping flag for CUDAFortran.
   integer(I4P)                        :: error                                !< Error trapping flag for CUDAFortran.

   !$cuf kernel do(5) <<<*,*>>>
   do v=1, nv
      do k=1, nk
         do j=1, nj
            do i=1, ni
               do b=1, blocks_number
                  if(phi_gpu(b,i,j,k,1) < 0.) then
                     q_gpu(b,i,j,k,v) = prhs_gpu(b,i,j,k,v) + qnrk * fl_gpu(b,i,j,k,v)
                  endif
               enddo
            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine compute_rk_linear_gpu_cuf

   subroutine compute_rk_prhs_gpu_cuf(ni, nj, nk, ngc, nv, blocks_number, dt, s, q_gpu, prhs_gpu, fl_gpu, phi_gpu, qnrk)
   !< Initialize RK stage with q_gpu.
   integer(I4P), intent(in)            :: ni                                     !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj                                     !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk                                     !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: ngc                                    !< Ghost cells number.
   integer(I4P), intent(in)            :: nv                                     !< Number of conservative varibales.
   integer(I4P), intent(in)            :: blocks_number                          !< Number of blocks.
   real(R8P),    intent(in)            :: dt                                     !< Time step.
   real(R8P),    intent(in)            :: qnrk                                   !< Time step.
   integer(I4P), intent(in)            :: s                                      !< Stage to initialize.
   real(R8P),    intent(in),    device ::   q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Conservative field.
   real(R8P),    intent(in),    device ::  fl_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Conservative field.
   real(R8P),    intent(in),    device :: phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Conservative field.
   real(R8P),    intent(inout), device :: prhs_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< RK stage.
   integer(I4P)                        :: i, j, k, b, v, ss                      !< Counter.
   integer(I4P)                        :: iercuda                                !< Error trapping flag for CUDAFortran.

   !$cuf kernel do(5) <<<*,*>>>
   do v=1, nv
      do k=1, nk
         do j=1, nj
            do i=1, ni
               do b=1, blocks_number
                  if(phi_gpu(b,i,j,k,1) < 0.) then
                     prhs_gpu(b,i,j,k,v) = q_gpu(b,i,j,k,v) + qnrk * fl_gpu(b,i,j,k,v)
                  endif
               enddo
            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine compute_rk_prhs_gpu_cuf

   subroutine invert_eikonal_field(ni, nj, nk, ngc, nv, blocks_number, q_gpu, q_invert_gpu, phi_gpu, bcs_type)
   !< Initialize RK stage with q_gpu.
   integer(I4P), intent(in)            :: ni                                        !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj                                        !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk                                        !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: ngc                                       !< Ghost cells number.
   integer(I4P), intent(in)            :: nv                                        !< Number of conservative varibales.
   integer(I4P), intent(in)            :: blocks_number                             !< Number of blocks.
   integer(I4P), intent(in)            :: bcs_type                                  !< Immersed boundary type.
   real(R8P),    intent(in),    device ::  q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)        !< Conservative field.
   real(R8P),    intent(in),    device ::  phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Distance field.
   real(R8P),    intent(inout), device ::  q_invert_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Inverted internal field.
   integer(I4P)                        :: i, j, k, b, v, ss                         !< Counter.
   integer(I4P)                        :: iercuda                                   !< Error trapping flag for CUDAFortran.
   real(R8P)                           :: n_phi_x, n_phi_y, n_phi_z                 !< Distance function normals.
   real(R8P)                           :: n_phi_mod, un_mod                         !< Distance abs normal and normal velocity.

   if(bcs_type == BCS_VISCOUS) then
      !$cuf kernel do(4) <<<*,*>>>
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               do b=1, blocks_number
                  if(phi_gpu(b,i,j,k,1) < 0) then
                      do v=1,nv
                         q_invert_gpu(b,i,j,k,v) = q_gpu(b,i,j,k,v)
                      enddo
                  else
                      q_invert_gpu(b,i,j,k,1) =   q_gpu(b,i,j,k,1)
                      q_invert_gpu(b,i,j,k,2) = - q_gpu(b,i,j,k,2)
                      q_invert_gpu(b,i,j,k,3) = - q_gpu(b,i,j,k,3)
                      q_invert_gpu(b,i,j,k,4) = - q_gpu(b,i,j,k,4)
                      q_invert_gpu(b,i,j,k,5) =   q_gpu(b,i,j,k,5)
                      if(nv == 6) q_invert_gpu(b,i,j,k,nv) = q_gpu(b,i,j,k,nv)
                  endif
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   elseif(bcs_type == BCS_EULER) then
      !$cuf kernel do(4) <<<*,*>>>
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               do b=1, blocks_number
                  if(phi_gpu(b,i,j,k,1) < 0) then
                     do v=1,nv
                        q_invert_gpu(b,i,j,k,v) = q_gpu(b,i,j,k,v)
                     enddo
                  else
                     n_phi_x = phi_gpu(b,i+1,j,k,1)-phi_gpu(b,i-1,j,k,1)
                     n_phi_y = phi_gpu(b,i,j+1,k,1)-phi_gpu(b,i,j-1,k,1)
                     n_phi_z = phi_gpu(b,i,j,k+1,1)-phi_gpu(b,i,j,k-1,1)
                     n_phi_mod = sqrt(n_phi_x**2+n_phi_y**2+n_phi_z**2)
                     n_phi_x = n_phi_x/n_phi_mod
                     n_phi_y = n_phi_y/n_phi_mod
                     n_phi_z = n_phi_z/n_phi_mod
                     un_mod = q_gpu(b,i,j,k,2)*n_phi_x+q_gpu(b,i,j,k,3)*n_phi_y+q_gpu(b,i,j,k,4)*n_phi_z

                     q_invert_gpu(b,i,j,k,1) = q_gpu(b,i,j,k,1)
                     q_invert_gpu(b,i,j,k,2) = q_gpu(b,i,j,k,2) - 2*un_mod*n_phi_x
                     q_invert_gpu(b,i,j,k,3) = q_gpu(b,i,j,k,3) - 2*un_mod*n_phi_y
                     q_invert_gpu(b,i,j,k,4) = q_gpu(b,i,j,k,4) - 2*un_mod*n_phi_z
                     q_invert_gpu(b,i,j,k,5) = q_gpu(b,i,j,k,5)
                     if(nv == 6) q_invert_gpu(b,i,j,k,6) = q_gpu(b,i,j,k,6)
                  endif
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   endif
   endsubroutine invert_eikonal_field

   subroutine compute_umax_cuf(b, ni, nj, nk, ngc, ns, dx, dy, dz, q_aux_gpu, umax, mu_star)
   !< Compute maximum speed by means of CUF threads.
   integer(I4P), intent(in)         :: b                                     !< Block index.
   integer(I4P), intent(in)         :: ni                                    !< Grid cells number in I direction.
   integer(I4P), intent(in)         :: nj                                    !< Grid cells number in J direction.
   integer(I4P), intent(in)         :: nk                                    !< Grid cells number in K direction.
   integer(I4P), intent(in)         :: ngc                                   !< Ghost cells number.
   integer(I4P), intent(in)         :: ns                                    !< Number of species.
   real(R8P),    intent(in)         :: dx                                    !< X space step.
   real(R8P),    intent(in)         :: dy                                    !< Y space step.
   real(R8P),    intent(in)         :: dz                                    !< Z space step.
   real(R8P),    intent(in)         :: mu_star                               !< Z space step.
   real(R8P),    intent(in), device :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary varibales.
   real(R8P),    intent(out)        :: umax                                  !< Maximum speed.
   real(R8P)                        :: ss                                    !< Speed of sound.
   integer(I4P)                     :: i, j, k                               !< Counter.
   integer(I4P)                     :: iercuda                               !< Error trapping flag for CUDAFortran.
   real(R8P)                        :: dx_locale, dy_locale, dz_locale

   umax = 0._R8P
   !$cuf kernel do(3) <<<*,*>>>
   do k=1, nk
      do j=1, nj
         do i=1, ni
            ss = q_aux_gpu(b,i,j,k,9)
            !umax = max(umax, abs(q_aux_gpu(b,i,j,k,2)) + ss, &
            !                 abs(q_aux_gpu(b,i,j,k,3)) + ss, &
            !                 abs(q_aux_gpu(b,i,j,k,4)) + ss)
            dx_locale = dx/2.
            dy_locale = dy/2.
            dz_locale = dz/2.
            umax = max(umax, (abs(q_aux_gpu(b,i,j,k,2)) + ss)/dx_locale +      &
                              2.*mu_star/(q_aux_gpu(b,i,j,k,1))/dx_locale**2 + &
                             (abs(q_aux_gpu(b,i,j,k,3)) + ss)/dy_locale +      &
                              2.*mu_star/(q_aux_gpu(b,i,j,k,1))/dy_locale**2 + &
                             (abs(q_aux_gpu(b,i,j,k,4)) + ss)/dz_locale +      &
                              2.*mu_star/(q_aux_gpu(b,i,j,k,1))/dz_locale**2)
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine compute_umax_cuf

endmodule adam_equation_nasto_gpu_object
