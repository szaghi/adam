!< ADAM, Navier-Stokes equations system class definition, GPU (NVF) backend.
module adam_nasto_nvf_object
!< ADAM, Navier-Stokes equations system class definition, GPU (NVF) backend.

use adam_amr_object
use adam_base_gpu_object
use adam_ib_object
use adam_memory_cpu_lib
use adam_memory_gpu_lib
use adam_parameters
use adam_nasto_bc_object
use adam_nasto_common_object
use adam_nasto_schemes_object
use adam_nasto_nvf_kernels
use penf
use MPI
use CUDAFOR
use ISO_C_BINDING

implicit none
private
public :: nasto_nvf_object

type, extends(nasto_common_object) :: nasto_nvf_object
   !< Navier-Stokes equations system class definition, GPU (NVF) backend.
   ! ADAM library objects
   type(base_gpu_object) :: base_gpu !< The base GPU handler.
   ! GPU data
   integer(I4P), allocatable, device :: ror_schemes_gpu(:)         !< ROR WENO schemes (GPU).
   integer(I4P), allocatable, device :: ror_ivar_gpu(:)            !< ROR variables indexes (GPU).
   integer(I4P), allocatable, device :: ror_stats_gpu(:,:,:,:,:)   !< ROR statistics (GPU.)
   real(R8P),    allocatable, device :: dq_gpu(:,:,:,:,:)          !< Eikonal right hand side.
   real(R8P),    allocatable, device :: fl_gpu(:,:,:,:,:)          !< Residuals.
   real(R8P),    allocatable, device :: flx_gpu(:,:,:,:,:)         !< Fluxes along x.
   real(R8P),    allocatable, device :: fly_gpu(:,:,:,:,:)         !< Fluxes along y.
   real(R8P),    allocatable, device :: flz_gpu(:,:,:,:,:)         !< Fluxes along z.
   real(R8P),    allocatable, device :: prhs_gpu(:,:,:,:,:)        !< Prhs for Runge-Kutta.
   real(R8P),    allocatable, device :: fd_coeff1_gpu(:)           !< Diffusive fluxes integration coefficients, first order.
   real(R8P),    allocatable, device :: fd_coeff2_gpu(:)           !< Diffusive fluxes integration coefficients, second order.
   real(R8P),    allocatable, device :: fc_coeff_gpu(:,:)          !< Convective fluxes integration coefficients.
   real(R8P),    allocatable, device :: q_aux_gpu(:,:,:,:,:)       !< Auxiliary cell centered variables.
   real(R8P),    allocatable, device :: q_gpu(:,:,:,:,:)           !< Field cell centered variables.
   real(R8P),    allocatable, device :: q_old_gpu(:,:,:,:,:)       !< Field cell centered variables (old iteration).
   real(R8P),    allocatable, device :: q_invert_gpu(:,:,:,:,:)    !< Field cell with boundary set on immersed bodies.
   real(R8P),    allocatable, device :: gplus_x_gpu(:,:,:,:,:)     !< Positive fluxes for weno-x.
   real(R8P),    allocatable, device :: gminus_x_gpu(:,:,:,:,:)    !< Negative fluxes for weno-x.
   real(R8P),    allocatable, device :: gplus_y_gpu(:,:,:,:,:)     !< Positive fluxes for weno-y.
   real(R8P),    allocatable, device :: gminus_y_gpu(:,:,:,:,:)    !< Negative fluxes for weno-y.
   real(R8P),    allocatable, device :: gplus_z_gpu(:,:,:,:,:)     !< Positive fluxes for weno-z.
   real(R8P),    allocatable, device :: gminus_z_gpu(:,:,:,:,:)    !< Negative fluxes for weno-z.
   real(R8P),    allocatable, device :: phi_gpu(:,:,:,:,:)         !< Distance function on GPU.
   real(R8P),    allocatable, device :: q_bc_vars_gpu(:, :)        !< Variables array for boundary conditions on GPU.
   real(R8P),    allocatable, device :: q_bcs_vars_gpu(:, :)       !< Variables array for immersed boundary on GPU.
   integer(I4P), allocatable, device :: cell_scheme_gpu(:,:,:,:,:) !< Modified order close to solids (GPU variable).
   contains
      ! auxiliary methods
      procedure, pass(self) :: allocate_gpu            !< Allocate GPU data.
      procedure, pass(self) :: check_cuda_error        !< Check if CUDA error occurs and abort in case.
      procedure, nopass     :: compute_cuda_dimensions !< Compute CUDA grid dimensions for GPU parallel computations.
      procedure, pass(self) :: copy_cpu_gpu            !< Copy data from CPU to GPU.
      procedure, pass(self) :: copy_gpu_cpu            !< Copy data from GPU to CPU.
      procedure, pass(self) :: destroy                 !< Destroy the equation.
      procedure, pass(self) :: initialize              !< Initialize the equation.
      ! AMR methods
      procedure, pass(self) :: amr_update       !< Do AMR update.
      procedure, pass(self) :: mark_by_grad_var !< Mark blocks to be refined/derefined by a `grad(var)` value.
      procedure, pass(self) :: mark_by_geo      !< Mark blocks to be refined/derefined by a `grad(var)` value.
      procedure, pass(self) :: move_phi         !< Move phi and the actual ptree representation.
      procedure, pass(self) :: refine_uniform   !< Refine all blocks uniformly.
      procedure, pass(self) :: update_phi       !< Update IB distance.
      ! eikonal methods
      procedure, pass(self) :: evolve_eikonal_q_gpu !< Evolve eikonal equation over q.
      procedure, pass(self) :: invert_eikonal_q_gpu !< Invert eikonal equation over q.
      ! IO methods
      procedure, pass(self) :: load_restart_files   !< Load restart files.
      procedure, pass(self) :: save_simulation_data !< Save all simulation data.
      procedure, pass(self) :: save_restart_files   !< Save restart files.
      procedure, pass(self) :: save_hdf5            !< Save simulation data in HDF5 format.
      ! IC/BC
      procedure, pass(self) :: set_boundary_conditions !< Set boundary conditions of equation.
      procedure, pass(self) :: set_initial_conditions  !< Set initial conditions of equation.
      procedure, pass(self) :: update_ghost_gpu        !< Update ghost cells and set boundary conditions.
      procedure, pass(self) :: update_ghost_fluxes_gpu !< Update fluxes cells and set boundary conditions.
      ! numerical methods
      procedure, pass(self) :: compute_dt        !< Compute time step.
      procedure, pass(self) :: compute_q_aux_gpu !< Compute auxiliary variables.
      procedure, pass(self) :: compute_residuals !< Compute residuals.
      procedure, pass(self) :: compute_rk_q_gpu  !< Compute RK approximation over q.
      procedure, pass(self) :: integrate         !< Perform one step integration.
      procedure, pass(self) :: simulate          !< Perform the simulation.
endtype nasto_nvf_object

contains
   ! auxiliary methods
   subroutine allocate_gpu(self)
   !< Allocate common data.
   class(nasto_nvf_object), intent(inout) :: self !< The equation.

   print '(A)', self%mpih%myrankstr//'nasto_nvf_object%allocate_gpu start'
   ! allocate by CPU data copy
   self%ror_schemes_gpu = self%schemes%ror_schemes
   self%ror_ivar_gpu    = self%schemes%ror_ivar

   self%fc_coeff_gpu  = self%schemes%fc_coeff
   self%fd_coeff1_gpu = self%schemes%fd_coeff1
   self%fd_coeff2_gpu = self%schemes%fd_coeff2

   self%q_bc_vars_gpu  = self%bc%q
   self%q_bcs_vars_gpu = self%ib%q

   if (self%ib%solids_number > 0) self%phi_gpu = self%phi

   ! allocate standalone
   associate(nv=>self%nv, ns=>self%ns, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk, &
             nb=>self%nb, nv_aux=>self%nv_aux, iweno=>self%schemes%iweno)

   ! call alloc_var_gpu(var=self%q_gpu, msg=self%mpih%myrankstr//'equation_nasto_gpu%alloc(q_gpu) ', verbose=.false.,&
   !                    ulb=reshape([1,nb,1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,nv],[2,5]))

   !< @NOTE gplus e gminus hanno Nb e Nv invertiti rispetto a tutti gli altri array GPU, errore o voluto?
                                          allocate(self%gplus_x_gpu (    1:nv, 1:2*iweno,    1:nj,         1:nk,         1:nb    ))
                                          allocate(self%gminus_x_gpu(    1:nv, 1:2*iweno,    1:nj,         1:nk,         1:nb    ))
                                          allocate(self%gplus_y_gpu (    1:nv, 1:2*iweno,    1:ni,         1:nk,         1:nb    ))
                                          allocate(self%gminus_y_gpu(    1:nv, 1:2*iweno,    1:ni,         1:nk,         1:nb    ))
                                          allocate(self%gplus_z_gpu (    1:nv, 1:2*iweno,    1:ni,         1:nj,         1:nb    ))
                                          allocate(self%gminus_z_gpu(    1:nv, 1:2*iweno,    1:ni,         1:nj,         1:nb    ))
                                          allocate(self%q_gpu(           1:nb, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nv    ))
                                          allocate(self%q_aux_gpu(       1:nb, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nv_aux))
                                          allocate(self%q_old_gpu(       1:nb, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nv    ))
                                          allocate(self%q_invert_gpu(    1:nb, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nv    ))
                                          allocate(self%fl_gpu(          1:nb, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nv    ))
                                          allocate(self%flx_gpu(         1:nb, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nv    ))
                                          allocate(self%fly_gpu(         1:nb, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nv    ))
                                          allocate(self%flz_gpu(         1:nb, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nv    ))
                                          allocate(self%dq_gpu(          1:nb, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nv    ))
                                          allocate(self%prhs_gpu(        1:nb, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nv    ))
                                          allocate(self%cell_scheme_gpu( 1:nb, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:3     ))
   if (self%schemes%enable_ror_stats > 0) allocate(self%ror_stats_gpu(   1:nb, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:3     ))
   self%q_gpu           = 0._R8P
   self%q_aux_gpu       = 0._R8P
   self%q_old_gpu       = 0._R8P
   self%q_invert_gpu    = 0._R8P
   self%fl_gpu          = 0._R8P
   self%flx_gpu         = 0._R8P
   self%fly_gpu         = 0._R8P
   self%flz_gpu         = 0._R8P
   self%dq_gpu          = 0._R8P
   self%prhs_gpu        = 0._R8P
   self%cell_scheme_gpu = self%schemes%iweno
   endassociate
   print '(A)', self%mpih%myrankstr//'nasto_nvf_object%allocate_gpu finish'
   endsubroutine allocate_gpu

   subroutine check_cuda_error(self, error_code, msg)
   !< Check if CUDA error occurs and abort in case.
   class(nasto_nvf_object), intent(inout)        :: self       !< The equation.
   integer(I4P),            intent(in), optional :: error_code !< Abort error code.
   character(*),            intent(in), optional :: msg        !< Error message.
   character(:), allocatable                     :: msg_       !< Error message, local variable.

   self%mpih%error = cudaGetLastError()
   if (self%mpih%error /= cudaSuccess) then
      msg_ = cudaGetErrorString(self%mpih%error) ; if (present(msg)) msg_ = '"'//msg//'": '//msg_
      call self%mpih%abort(error_code=error_code, msg=msg_)
   endif
   endsubroutine

   subroutine compute_cuda_dimensions(grid_x, grid_y, grid, tBlock)
   !< Compute CUDA grid dimensions for GPU parallel computations.
   integer(I4P), intent(in)  :: grid_x, grid_y !< CUDA grid xy dimensions.
   type(dim3), intent(inout) :: grid, tBlock   !< CUDA grid and block.

   tBlock = dim3(32,8,1)
   grid = dim3(ceiling(real(grid_x)/tBlock%x), ceiling(real(grid_y)/tBlock%y), 1)
   endsubroutine compute_cuda_dimensions

   subroutine copy_cpu_gpu(self)
   !< Copy data from CPU to GPU.
   class(nasto_nvf_object), intent(inout) :: self !< The equation.

   call self%base_gpu%copy_transpose_cpu_gpu(nv=self%nv, q_cpu=self%field%q, q_gpu=self%q_gpu)
   call self%base_gpu%copy_cpu_gpu(verbose=.false.)
   endsubroutine copy_cpu_gpu

   subroutine copy_gpu_cpu(self, compute_q_aux)
   !< Copy data from GPU to CPU.
   class(nasto_nvf_object), intent(inout)        :: self          !< The equation.
   logical,                 intent(in), optional :: compute_q_aux !< Flag to compute auxiliary variables.

   call self%base_gpu%copy_transpose_gpu_cpu(nv=self%nv, q_gpu=self%q_gpu, q_cpu=self%field%q)
   if (present(compute_q_aux)) then
      if (compute_q_aux) then
         call self%compute_q_aux_gpu(q_gpu=self%q_gpu, q_aux_gpu=self%q_aux_gpu)
         call self%base_gpu%copy_transpose_gpu_cpu(nv=self%nv_aux, q_gpu=self%q_aux_gpu, q_cpu=self%q_aux)
      endif
   endif
   endsubroutine copy_gpu_cpu

   subroutine destroy(self)
   !< Destroy the equation.
   class(nasto_nvf_object), intent(inout) :: self  !< The equation.
   type(nasto_nvf_object)                 :: fresh !< Fresh equation.

   ! TODO to be implemented
   endsubroutine destroy

   subroutine initialize(self, filename)
   !< Initialize the equation.
   class(nasto_nvf_object), intent(inout) :: self         !< The equation.
   character(*),            intent(in)    :: filename     !< Input file name.

   call self%base_gpu%initialize_gpu(do_mpi_init=.true.)
   print '(A)', self%base_gpu%mpih%myrankstr//'nasto_nvf_object%initialize start'
   call self%initialize_common(filename=filename, memory_avail=self%base_gpu%memory_avail)
   call self%base_gpu%initialize(field=self%adam%field, nv_aux=self%nv_aux, verbose=.false.)
   call self%allocate_gpu
   print '(A)', self%mpih%myrankstr//'nasto_nvf_object%initialize finish'
   endsubroutine initialize

   ! AMR methods
   subroutine amr_update(self)
   !< Do AMR update.
   class(nasto_nvf_object), intent(inout) :: self                 !< The equation.
   integer(I4P)                           :: iterations_          !< Number of AMR iterations, local var.
   logical                                :: is_grid_changed      !< Flag to check grid changes for each marker.
   logical                                :: is_grid_changed_all  !< Flag to check grid changes for each iter.
   integer(I4P)                           :: b, i, j, k, i_marker !< Counter.
   type(amr_marker_object)                :: amr_marker           !< Current amr marker.

   amr: do i=1, self%amr%iters
      is_grid_changed_all = .false.
      do i_marker=1, self%amr%markers_number
         amr_marker = self%amr%markers(i_marker)
         call self%update_ghost_gpu(q_gpu=self%q_gpu)
         select case(amr_marker%mode)
         case(AMR_GEO)
            call self%mark_by_geo(delta_fine=amr_marker%delta_fine, delta_coarse=amr_marker%delta_coarse)
         case(AMR_GRAD)
            select case(amr_marker%field)
            case(1)
               call self%mark_by_grad_var(grad_tol=amr_marker%tol, delta_fine=amr_marker%delta_fine, &
                                          delta_coarse=amr_marker%delta_coarse, ivar=amr_marker%ivar)
            case(2)
               call self%mark_by_grad_var(grad_tol=amr_marker%tol, delta_fine=amr_marker%delta_fine, &
                                          delta_coarse=amr_marker%delta_coarse, ivar=amr_marker%ivar)
            endselect
         endselect
         call self%copy_gpu_cpu ! needed for adam%amr_update
         call self%adam%amr_update(is_marked_by_field=.true., do_blocks_reorder=.false., is_grid_changed=is_grid_changed)
         if (self%ib%solids_number > 0) call self%update_phi()
         call self%copy_cpu_gpu
         is_grid_changed_all = is_grid_changed_all.or.is_grid_changed
      enddo
      if (.not.is_grid_changed_all) then
          print '(A)', self%mpih%myrankstr//'AMR Grid stabilized after : '//trim(str(i))//' AMR iterations'
          exit amr
       elseif (i==self%amr%iters) then
          print '(A)', self%mpih%myrankstr//'AMR Grid is NOT stabilized after : '//trim(str(i))//' AMR iterations'
      endif
   enddo amr
   endsubroutine amr_update

   subroutine mark_by_geo(self, delta_fine, delta_coarse, threshold, do_init)
   !< Mark blocks to be refined/derefined by a `grad(rho)` value.
   class(nasto_nvf_object), intent(inout)        :: self           !< The equation.
   real(R8P),               intent(in)           :: delta_fine     !< Maximum cell delta in fine grids.
   real(R8P),               intent(in)           :: delta_coarse   !< Minimum cell delta in coarse grids.
   real(R8P),               intent(in), optional :: threshold      !< Threshold for sphere proximity.
   real(R8P)                                     :: threshold_     !< Threshold for sphere proximity, local var.
   real(R8P)                                     :: max_cell_delta !< Maximum cell delta.
   real(R8P)                                     :: distance       !< Value (max) of gradient of rho.
   integer(I4P)                                  :: b              !< Counter.
   logical, optional,       intent(in)           :: do_init
   logical                                       :: do_init_

   do_init_ = .true.    ; if (present(do_init)) do_init_ = do_init
   threshold_ = 2.2_R8P ; if (present(threshold)) threshold_ = threshold
   if (do_init_) self%field%refinements_needed = [(TO_BE_DEREFINED,b=1,self%blocks_number)]
   associate (ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, &
              blocks_number=>self%blocks_number, ns=>self%ns, dxyz=>self%field%dxyz, phi=>self%phi)
      do b=1, blocks_number
         distance = 1._R8P
         if (maxval(phi(b,:,:,:,1))*minval(phi(b,:,:,:,1)) < 0._R8P) then
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
      endif
      endfunction max_cell_delta_dist
   endsubroutine mark_by_geo

   subroutine mark_by_grad_var(self, grad_tol, delta_fine, delta_coarse, ivar, threshold, do_init)
   !< Mark blocks to be refined/derefined by a `grad(rho)` value.
   class(nasto_nvf_object), intent(inout)        :: self           !< The equation.
   real(R8P),               intent(in)           :: grad_tol       !< Gradiend tolerance value.
   real(R8P),               intent(in)           :: delta_fine     !< Maximum cell delta in fine grids.
   real(R8P),               intent(in)           :: delta_coarse   !< Minimum cell delta in coarse grids.
   integer(I4P),            intent(in), optional :: ivar           !< Variable for marking.
   real(R8P),               intent(in), optional :: threshold      !< Threshold for sphere proximity.
   logical,                 intent(in), optional :: do_init        !< Re-initialize refinements queries.
   integer(I4P)                                  :: ivar_          !< Variable for marking (local var).
   logical                                       :: do_init_       !< Re-initialize refinements queries, local var.
   real(R8P)                                     :: threshold_     !< Threshold for sphere proximity, local var.
   real(R8P)                                     :: max_cell_delta !< Maximum cell delta.
   real(R8P)                                     :: grad_var       !< Value (max) of gradient of var.
   integer(I4P)                                  :: b              !< Counter.

   ivar_     = 1_R4P    ; if (present(ivar)) ivar_ = ivar
   do_init_ = .true.    ; if (present(do_init)) do_init_ = do_init
   threshold_ = 2.2_R8P ; if (present(threshold)) threshold_ = threshold
   if(do_init_) self%field%refinements_needed = [(TO_BE_DEREFINED,b=1,self%blocks_number)]
   associate (ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, &
              blocks_number=>self%blocks_number, ns=>self%ns, dxyz=>self%field%dxyz)
      call self%update_ghost_gpu(q_gpu=self%q_gpu)
      call self%compute_q_aux_gpu(q_gpu=self%q_gpu, q_aux_gpu=self%q_aux_gpu)
      do b=1, blocks_number
         call compute_q_gradient_cuf(b=b, ni=ni, nj=nj, nk=nk, ngc=ngc, &
                                     dx=dxyz(1,b), dy=dxyz(2,b), dz=dxyz(3,b), q_gpu=self%q_aux_gpu, ivar=ivar_, gradient=grad_var)
         max_cell_delta = max_cell_delta_grad(grad=grad_var)
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
   class(nasto_nvf_object), intent(inout) :: self        !< The equation.
   real(R8P),               intent(in)    :: velocity(3) !< Velocity of the movement.

   call move_phi_cuf(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                     velocity=velocity, phi_gpu=self%phi_gpu, dphi_gpu=self%dq_gpu)
   endsubroutine move_phi

   subroutine refine_uniform(self, refinement_levels)
   !< Refine all blocks uniformly.
   class(nasto_nvf_object), intent(inout) :: self              !< The equation.
   integer(I4P),            intent(in)    :: refinement_levels !< Number of refinement to be performed.
   integer(I4P)                           :: l                 !< Counter.

   do l=1, refinement_levels
      call self%adam%tree%mark_all_nodes(mark=TO_BE_REFINED)
      call self%adam%amr_update(do_blocks_reorder=.false., do_mpi_redistribute=.true.)
   enddo
   endsubroutine

   subroutine update_phi(self)
   !< Update IB distance.
   class(nasto_nvf_object), intent(inout) :: self                      !< The equation.
   integer(I4P)                           :: b, i, j, k, ib, l         !< Counter.
   real(R8P)                              :: query_x, query_y, query_z !< Query point coordinates.
   real(R8P)                              :: near_x, near_y, near_z    !< Nearest point coordinates.
   real(R8P)                              :: distance                  !< Distance from solid.
   logical                                :: inside                    !< Inside/outside boolean.

   associate(blocks_number=>self%blocks_number, ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc,         &
             x_cell=>self%field%x_cell, y_cell=>self%field%y_cell, z_cell=>self%field%z_cell,                 &
             ptree => self%ptree, phi=>self%phi, phi_gpu=>self%phi_gpu, solids_number=>self%ib%solids_number, &
             ib_reduction_extent => self%schemes%ib_reduction_extent, ib_reduced_order => self%schemes%ib_reduced_order)
   print '(A)', self%mpih%myrankstr//'update IB distance start'
   do ib=1, solids_number
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

            distance = - (sqrt((query_x-10._R8P)**2+(query_y-10._R8P)**2+(query_z-10._R8P)**2)-1.0_R8P)

            phi(b,i,j,k,ib) = distance
         enddo
         enddo
         enddo
      enddo

      self%schemes%cell_scheme = self%schemes%iweno

      if (ib_reduction_extent > 0) then
          print*,'Reduced order close to solids. Extent/base/reduced-order: ',ib_reduction_extent, &
             self%schemes%ror_schemes(1), ib_reduced_order

          do b=1,blocks_number
             do j=1,nj
             do k=1,nk
             idir: do i=0,ni
             do l=1,ib_reduction_extent
                if(phi(b,i+l,j,k,ib)>0 .or. phi(b,i-l+1,j,k,ib)>0) then
                    self%schemes%cell_scheme(b,i,j,k,1) = ib_reduced_order
                    cycle idir
                endif
             enddo
             enddo idir
             enddo
             enddo
          enddo
          do b=1,blocks_number
             do i=1,ni
             do k=1,nk
             jdir: do j=0,nj
             do l=1,ib_reduction_extent
                if(phi(b,i,j+l,k,ib)>0 .or. phi(b,i,j-l+1,k,ib)>0) then
                    self%schemes%cell_scheme(b,i,j,k,2) = ib_reduced_order
                    cycle jdir
                endif
             enddo
             enddo jdir
             enddo
             enddo
          enddo
          do b=1,blocks_number
             do j=1,nj
             do i=1,ni
             kdir: do k=0,nk
             do l=1,ib_reduction_extent
                if(phi(b,i,j,k+l,ib)>0 .or. phi(b,i,j,k-l+1,ib)>0) then
                    self%schemes%cell_scheme(b,i,j,k,3) = ib_reduced_order
                    cycle kdir
                endif
             enddo
             enddo kdir
             enddo
             enddo
          enddo
      endif
   enddo
   ! TODO: next assignments should be reduced only to blocks_number.
   phi_gpu = phi
   self%cell_scheme_gpu = self%schemes%cell_scheme
   endassociate
   print '(A)', self%mpih%myrankstr//'update IB distance finish'
   endsubroutine update_phi

   ! eikonal methods
   subroutine evolve_eikonal_q_gpu(self)
   !< Evolve eikonal equation over q.
   class(nasto_nvf_object), intent(inout) :: self         !< The equation.
   type(dim3)                             :: grid, tBlock !< CUDA grid and block.

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>self%nv, blocks_number=>self%blocks_number,  &
             phi_gpu=>self%phi_gpu,                                                                                 &
             dx_gpu=>self%base_gpu%dxyz_gpu(:,1),                                                                   &
             dy_gpu=>self%base_gpu%dxyz_gpu(:,2),                                                                   &
             dz_gpu=>self%base_gpu%dxyz_gpu(:,3),                                                                   &
             dq_gpu=>self%dq_gpu, q_gpu=>self%q_gpu)
   if (blocks_number > 0) then
      call self%compute_cuda_dimensions(grid_x=blocks_number, grid_y=ni, grid=grid, tBlock=tBlock)
      call compute_eikonal_dq_gpu<<<grid, tBlock>>>(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, &
                                                    phi_gpu=phi_gpu, dx_gpu=dx_gpu, dy_gpu=dy_gpu, dz_gpu=dz_gpu,     &
                                                    dq_gpu=dq_gpu, q_gpu=q_gpu)
   endif

   call evolve_eikonal_q_gpu_cuf(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, &
                                 phi_gpu=self%phi_gpu,                                             &
                                 dx_gpu=self%base_gpu%dxyz_gpu(:,1),                               &
                                 dy_gpu=self%base_gpu%dxyz_gpu(:,2),                               &
                                 dz_gpu=self%base_gpu%dxyz_gpu(:,3),                               &
                                 dq_gpu=self%dq_gpu,                                               &
                                 q_gpu=self%q_gpu)
   endassociate
   endsubroutine evolve_eikonal_q_gpu

   subroutine invert_eikonal_q_gpu(self)
   !< Invert eikonal equation over q.
   class(nasto_nvf_object), intent(inout) :: self !< The equation.

   call invert_eikonal_q_gpu_cuf(BCS_VISCOUS=BCS_VISCOUS, BCS_EULER=BCS_EULER,                                                   &
                                 ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, nv=self%nv, blocks_number=self%blocks_number, &
                                 q_gpu=self%q_gpu(:,:,:,:,:), q_invert_gpu=self%q_invert_gpu(:,:,:,:,:),                         &
                                 phi_gpu=self%phi_gpu, bcs_type=self%ib%bc_type(1))
   endsubroutine invert_eikonal_q_gpu

   ! IO methods
   subroutine load_restart_files(self, t, time)
   !< Save restart files.
   class(nasto_nvf_object), intent(inout) :: self !< The equation.
   integer(I4P),            intent(out)   :: t    !< Time iteration.
   real(R8P),               intent(out)   :: time !< Time.

   call self%adam%load_restart_files(basename=self%io%restart_basename, t=t, time=time)
   call self%adam%make_comm_local_maps_ghost_bc
   call self%copy_cpu_gpu
   endsubroutine load_restart_files

   subroutine save_simulation_data(self)
   !< Save all simulation data.
   class(nasto_nvf_object), intent(inout) :: self !< The equation.

   if ((self%time%is_to_save(it_save=self%io%it_save)).or. &
       (mod(self%time%it,self%io%restart_save)==0).or.     &
       (self%slices%is_to_save(it=self%time%it,it_max=self%time%it_max,time=self%time%time,time_max=self%time%time_max))) then
      call self%update_ghost_gpu(q_gpu=self%q_gpu)
      call self%copy_gpu_cpu(compute_q_aux=.true.)

      if (self%time%is_to_save(it_save=self%io%it_save)) call self%save_hdf5
      if (mod(self%time%it,self%io%restart_save)==0) call self%save_restart_files
      if (self%slices%is_to_save(it=self%time%it,it_max=self%time%it_max,time=self%time%time,time_max=self%time%time_max))&
         call self%slices%save_mat(basename=self%io%output_basename, &
                                   it=self%time%it,                  &
                                   it_max=self%time%it_max,          &
                                   time=self%time%time,              &
                                   time_max=self%time%time_max,      &
                                   adam=self%adam,                   &
                                   q=self%field%q,                   &
                                   q_name=['rho','rhu','rhv','rhw','rhe'])
   endif
   endsubroutine save_simulation_data

   subroutine save_hdf5(self, output_basename)
   !< Save simulation data in HDF5 format.
   class(nasto_nvf_object), intent(inout)        :: self             !< The equation.
   character(*),            intent(in), optional :: output_basename  !< Output basename.
   character(:), allocatable                     :: output_basename_ !< Output basename, local var.

   call self%mpih%barrier(tictoc=.true.)
   print '(A)', self%mpih%myrankstr//'save HDF5 files t: '//trim(str(self%time%it,.true.))//', time: '//&
                trim(str(self%time%time,.true.))
   output_basename_ = trim(self%io%output_basename)//'-'//trim(strz(self%time%it,9))
   if (present(output_basename)) output_basename_ = trim(output_basename)
   call self%adam%save_hdf5(basename=trim(output_basename_),                                  &
                            q=self%field%q,                                                   &
                            q_aux=self%q_aux,                                                 &
                            q_name=['rho','rhu','rhv','rhw','rhe'],                           &
                            q_aux_name=['rhob','u','v','w','ya','tem','pres','ental','csp'],  &
                            with_cell_morton=.true.)
   call self%mpih%barrier(tictoc=.true.)
   ! print '(A, F18.10)', self%mpih%myrankstr//'step timing (save HDF5): ', self%mpih%tictoc_timing()
   endsubroutine save_hdf5

   subroutine save_restart_files(self)
   !< Save restart files.
   class(nasto_nvf_object), intent(inout) :: self !< The equation.

   call self%mpih%barrier(tictoc=.true.)
   print '(A)', self%mpih%myrankstr//'save restart files t: '//trim(str(self%time%it,.true.))//', time: '//&
                trim(str(self%time%time,.true.))
   call self%adam%save_restart_files(basename=self%io%restart_basename, t=self%time%it, time=self%time%time)
   call self%save_hdf5(output_basename=self%io%restart_basename)
   call self%mpih%barrier(tictoc=.true.)
   ! print '(A, F18.10)', self%mpih%myrankstr//'step timing (save restart): ', self%mpih%tictoc_timing()
   endsubroutine save_restart_files

   ! IC/BC
   subroutine set_boundary_conditions(self, q_gpu)
   !< Set boundary conditions of equation.
   class(nasto_nvf_object), intent(in)            :: self                  !< The equation.
   real(R8P),               intent(inout), device :: q_gpu(1:,         &
                                                           1-self%ngc:,&
                                                           1-self%ngc:,&
                                                           1-self%ngc:,1:) !< Conservative variables.

   if (allocated(self%base_gpu%local_map_bc_crown_gpu)) &
      call set_bc_q_gpu_cuf(BC_EXTRAPOLATION=BC_EXTRAPOLATION, BC_INFLOW=BC_INFLOW, &
                            nv=self%nv, ngc=self%ngc, cv=self%physics%eos(1)%cv, R=self%physics%eos(1)%R, &
                            local_map_bc_gpu=self%base_gpu%local_map_bc_crown_gpu,                        &
                            fec_1_6_array_gpu=self%base_gpu%fec_1_6_array_gpu,                            &
                            q_bc_vars_gpu=self%q_bc_vars_gpu,                                             &
                            q_gpu=q_gpu)
   endsubroutine set_boundary_conditions

   subroutine set_initial_conditions(self)
   !< Set initial conditions of field.
   class(nasto_nvf_object), intent(inout) :: self !< The equation.

   call self%ic%set_initial_conditions(physics=self%physics, field=self%field)
   call self%copy_cpu_gpu
   endsubroutine set_initial_conditions

   subroutine update_ghost_gpu(self, q_gpu, step)
   !< Update ghost cells.
   !< If not specified all steps are perfermod, syncronous computation
   class(nasto_nvf_object), intent(inout)         :: self            !< The equation.
   real(R8P),               intent(inout), device :: q_gpu(1:,         &
                                                           1-self%ngc:,&
                                                           1-self%ngc:,&
                                                           1-self%ngc:,&
                                                           1:)       !< Conservative variables.
   integer(I4P),            intent(in), optional  :: step            !< Step to be perfordmed in asyncronous comp.
   logical                                        :: do_local_update !< Flag for triggering local update.
   logical                                        :: do_set_bc       !< Flag for triggering setting bc.

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
   class(nasto_nvf_object), intent(inout)         :: self            !< The equation.
   real(R8P),               intent(inout), device :: flx_gpu(1:,         &
                                                             1-self%ngc:,&
                                                             1-self%ngc:,&
                                                             1-self%ngc:,&
                                                             1:)     !< Conservative variables.
   real(R8P),               intent(inout), device :: fly_gpu(1:,         &
                                                             1-self%ngc:,&
                                                             1-self%ngc:,&
                                                             1-self%ngc:,&
                                                             1:)     !< Conservative variables.
   real(R8P),               intent(inout), device :: flz_gpu(1:,         &
                                                             1-self%ngc:,&
                                                             1-self%ngc:,&
                                                             1-self%ngc:,&
                                                             1:)     !< Conservative variables.
   integer(I4P),            intent(in), optional  :: step            !< Step to be perfordmed in asyncronous comp.
   logical                                        :: do_local_update !< Flag for triggering local update.

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

   ! numerical methods
   subroutine compute_q_aux_gpu(self, q_gpu, q_aux_gpu)
   !< Compute auxiliary variables.
   class(nasto_nvf_object), intent(in)          :: self          !< The equation.
   real(R8P),               intent(in),  device :: q_gpu(1:,         &
                                                         1-self%ngc:,&
                                                         1-self%ngc:,&
                                                         1-self%ngc:,&
                                                         1:)     !< Conservative variables.
   real(R8P),               intent(out), device :: q_aux_gpu(1:,         &
                                                             1-self%ngc:,&
                                                             1-self%ngc:,&
                                                             1-self%ngc:,&
                                                             1:) !< Auxiliary variables.

   call compute_q_aux_cuf(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, ns=self%ns, blocks_number=self%blocks_number, &
                          R=self%physics%eos(1)%R, cv=self%physics%eos(1)%cv, g=self%physics%eos(1)%g,                    &
                          dha=self%physics%eos(1)%dha, q_gpu=q_gpu, q_aux_gpu=q_aux_gpu)
   endsubroutine compute_q_aux_gpu

   subroutine compute_dt(self)
   !< Compute maximum time step accordingly to CFL stabilty criterion.
   class(nasto_nvf_object), intent(inout) :: self !< The equation.
   real(R8P)                              :: umax !< Maximum speed of waves propagation.
   integer(I4P)                           :: b    !< Counter.

   call self%compute_q_aux_gpu(q_gpu=self%q_gpu, q_aux_gpu=self%q_aux_gpu)
   self%time%dt = huge(1._R8P)
   ! PERCHÈ NON ABBIAMO PARALLELIZATO DENTRO IL CUF ANCHE L'INDICE DI BLOCCO?
   do b=1, self%field%blocks_number
      call compute_umax_cuf(b, ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, ns=self%ns,           &
                            dx=self%field%dxyz(1,b), dy=self%field%dxyz(2,b), dz=self%field%dxyz(3,b), &
                            q_aux_gpu=self%q_aux_gpu, umax=umax, mu=self%physics%eos(1)%mu)
      self%time%dt = min(self%time%dt, self%time%CFL / umax)
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, self%time%dt, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, self%mpih%error)
   endsubroutine compute_dt

   subroutine compute_residuals(self)
   !< Compute residuals of equation.
   class(nasto_nvf_object), intent(inout) :: self         !< The equation.
   real(R8P)                              :: ib_eps       !< Tolerance ib delta ratio.
   integer(I4P)                           :: iercuda      !< Error trapping flag for CUDAFortran.
   type(dim3)                             :: grid, tBlock !< CUDA grid and block.

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk,                                                             &
             ngc=>self%ngc, ns=>self%ns, blocks_number=>self%blocks_number,                                     &
             dx_gpu=>self%base_gpu%dxyz_gpu(:,1),                                                               &
             dy_gpu=>self%base_gpu%dxyz_gpu(:,2),                                                               &
             dz_gpu=>self%base_gpu%dxyz_gpu(:,3),                                                               &
             q_aux_gpu=>self%q_aux_gpu, phi_gpu=>self%phi_gpu, fl_gpu=>self%fl_gpu,                             &
             flx_gpu=>self%flx_gpu, fly_gpu=>self%fly_gpu, flz_gpu=>self%flz_gpu,                               &
             cell_scheme_gpu=>self%cell_scheme_gpu, ror_stats_gpu=>self%ror_stats_gpu,                          &
             fc_coeff_gpu=>self%fc_coeff_gpu,                                                                   &
             gminus_x_gpu=>self%gminus_x_gpu, gminus_y_gpu=>self%gminus_y_gpu, gminus_z_gpu=>self%gminus_z_gpu, &
             gplus_x_gpu=>self%gplus_x_gpu, gplus_y_gpu=>self%gplus_y_gpu, gplus_z_gpu=>self%gplus_z_gpu,       &
             ror_schemes_gpu=>self%ror_schemes_gpu, ror_ivar_gpu=>self%ror_ivar_gpu,                            &
             ror_threshold=>self%schemes%ror_threshold, enable_ror_stats=>self%schemes%enable_ror_stats,        &
             lmax=>self%schemes%lmax, iweno=>self%schemes%iweno,                                                &
             cv=>self%physics%eos(1)%cv, g=>self%physics%eos(1)%g, R=>self%physics%eos(1)%R,                    &
             mu=>self%physics%eos(1)%mu, kd=>self%physics%eos(1)%kd, dha=>self%physics%eos(1)%dha)

   call self%check_cuda_error(error_code=-15, msg='CUDA error at start residuals computation')

   if (blocks_number > 0) then
      select case(self%schemes%fluxes_convective)
      case(SCHEME_FCONV_WENO_CENTRAL_2,SCHEME_FCONV_WENO_CENTRAL_4,SCHEME_FCONV_WENO_CENTRAL_6)
         call self%compute_cuda_dimensions(grid_x=blocks_number, grid_y=nj, grid=grid, tBlock=tBlock)
         call compute_flux_conv_x_central_kernel<<<grid, tBlock>>>(blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc,      &
                                                                   nv=ns+4,lmax=lmax,fc_coeff_gpu=fc_coeff_gpu,q_aux_gpu=q_aux_gpu,&
                                                                   dx_gpu=dx_gpu, flx_gpu=flx_gpu)

         call self%compute_cuda_dimensions(grid_x=blocks_number, grid_y=ni, grid=grid, tBlock=tBlock)
         call compute_flux_conv_y_central_kernel<<<grid, tBlock>>>(blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc,      &
                                                                   nv=ns+4,lmax=lmax,fc_coeff_gpu=fc_coeff_gpu,q_aux_gpu=q_aux_gpu,&
                                                                   dy_gpu=dy_gpu, fly_gpu=fly_gpu)

         call self%compute_cuda_dimensions(grid_x=blocks_number, grid_y=ni, grid=grid, tBlock=tBlock)
         call compute_flux_conv_z_central_kernel<<<grid, tBlock>>>(blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc,      &
                                                                   nv=ns+4,lmax=lmax,fc_coeff_gpu=fc_coeff_gpu,q_aux_gpu=q_aux_gpu,&
                                                                   dz_gpu=dz_gpu, flz_gpu=flz_gpu)
      case(SCHEME_FCONV_WENO_UPWIND)
         call self%compute_cuda_dimensions(grid_x=blocks_number, grid_y=nj, grid=grid, tBlock=tBlock)
         call compute_flux_conv_x_kernel<<<grid, tBlock>>>(blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=ns+4, &
                                                           iweno=iweno, dha=dha, g=g, R=R, cv=cv,                              &
                                                           ror_threshold=ror_threshold, enable_ror_stats=enable_ror_stats,     &
                                                           cell_scheme_gpu=cell_scheme_gpu, ror_ivar_gpu=ror_ivar_gpu,         &
                                                           ror_schemes_gpu=ror_schemes_gpu, q_aux_gpu=q_aux_gpu,               &
                                                           ror_stats_gpu=ror_stats_gpu,                                        &
                                                           gplus=gplus_x_gpu, gminus=gminus_x_gpu, flx_gpu=flx_gpu)

         call self%compute_cuda_dimensions(grid_x=blocks_number, grid_y=ni, grid=grid, tBlock=tBlock)
         call compute_flux_conv_y_kernel<<<grid, tBlock>>>(blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=ns+4, &
                                                           iweno=iweno, dha=dha, g=g, R=R, cv=cv,                              &
                                                           ror_threshold=ror_threshold, enable_ror_stats=enable_ror_stats,     &
                                                           cell_scheme_gpu=cell_scheme_gpu, ror_ivar_gpu=ror_ivar_gpu,         &
                                                           ror_schemes_gpu=ror_schemes_gpu, q_aux_gpu=q_aux_gpu,               &
                                                           ror_stats_gpu=ror_stats_gpu,                                        &
                                                           gplus=gplus_y_gpu, gminus=gminus_y_gpu, fly_gpu=fly_gpu)

         call self%compute_cuda_dimensions(grid_x=blocks_number, grid_y=ni, grid=grid, tBlock=tBlock)
         call compute_flux_conv_z_kernel<<<grid, tBlock>>>(blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=ns+4, &
                                                           iweno=iweno, dha=dha, g=g, R=R, cv=cv,                              &
                                                           ror_threshold=ror_threshold, enable_ror_stats=enable_ror_stats,     &
                                                           cell_scheme_gpu=cell_scheme_gpu, ror_ivar_gpu=ror_ivar_gpu,         &
                                                           ror_schemes_gpu=ror_schemes_gpu, q_aux_gpu=q_aux_gpu,               &
                                                           ror_stats_gpu=ror_stats_gpu,                                        &
                                                           gplus=gplus_z_gpu, gminus=gminus_z_gpu, flz_gpu=flz_gpu)
      endselect
   endif

   call self%check_cuda_error(error_code=-15, msg='CUDA error after convective fluxes computation')

   if (mu > 0.) call compute_fluxes_diffusive_cuf(blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=ns+4, &
                                                  mu=mu, kd=kd, q_aux_gpu=q_aux_gpu,                                  &
                                                  dx_gpu=dx_gpu, dy_gpu=dy_gpu, dz_gpu=dz_gpu,                        &
                                                  flx_gpu=flx_gpu, fly_gpu=fly_gpu, flz_gpu=flz_gpu)
   !@cuf iercuda=cudaDeviceSynchronize()

   call self%check_cuda_error(error_code=-15, msg='CUDA error after diffusive fluxes computation')

   ib_eps = 1.e-12_R8P
   call compute_fluxes_difference_cuf(blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=ns+4,                &
                                      fl_gpu=fl_gpu, flx_gpu=flx_gpu, fly_gpu=fly_gpu, flz_gpu=flz_gpu, phi_gpu=phi_gpu, &
                                      dx_gpu=dx_gpu, dy_gpu=dy_gpu, dz_gpu=dz_gpu, ib_eps=ib_eps)

   call self%check_cuda_error(error_code=-15, msg='CUDA error after fluxes difference computation')
   endassociate
   endsubroutine compute_residuals

   subroutine compute_rk_q_gpu(self, s)
   !< Compute RK approximation over q.
   class(nasto_nvf_object), intent(inout) :: self !< The equation.
   integer(I4P),            intent(in)    :: s    !< Current RK stage.

   call compute_rk_q_gpu_cuf(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, nv=self%nv, blocks_number=self%blocks_number, &
                             dt=self%time%dt, s=s, q_gpu=self%q_gpu, q_old_gpu=self%q_old_gpu,                               &
                             fl_gpu=self%fl_gpu, phi_gpu=self%phi_gpu,                                                       &
                             ark=self%schemes%ark(s), brk=self%schemes%brk(s), crk=self%schemes%crk(s))
   endsubroutine compute_rk_q_gpu

   subroutine integrate(self, t, do_ghost_syncro, residual)
   !< Perform one step integration.
   class(nasto_nvf_object), intent(inout)         :: self             !< The equation.
   real(R8P),               intent(in)            :: t                !< Time.
   logical,                 intent(in),  optional :: do_ghost_syncro  !< Flag to do syncrous ghost update.
   real(R8P),               intent(out), optional :: residual         !< Global residual.
   logical                                        :: do_ghost_syncro_ !< Flag to do syncrous ghost update, local var.
   integer(I4P)                                   :: s                !< Counter.
   integer(I4P)                                   :: i_eikonal        !< Counter.
   integer(I4P), parameter                        :: n_eikonal=2      !< Counter.

   do_ghost_syncro_ = .true. ; if (present(do_ghost_syncro)) do_ghost_syncro_ = do_ghost_syncro
   self%q_old_gpu = self%q_gpu
   do s=1, self%schemes%nrk
      if (self%ib%solids_number > 0) then
         call self%update_ghost_gpu(q_gpu=self%q_gpu)
         do i_eikonal=1,n_eikonal
            call MPI_Barrier(MPI_COMM_WORLD, self%mpih%error)
            call self%evolve_eikonal_q_gpu
            call self%update_ghost_gpu(q_gpu=self%q_gpu)
         enddo
         call self%invert_eikonal_q_gpu
      else
         ! added for restart debug
         self%q_invert_gpu = self%q_gpu
      endif
      call MPI_Barrier(MPI_COMM_WORLD, self%mpih%error)
      call self%compute_q_aux_gpu(q_gpu=self%q_invert_gpu, q_aux_gpu=self%q_aux_gpu)
      call self%compute_residuals
      call self%compute_rk_q_gpu(s=s)
   enddo
   endsubroutine integrate

   subroutine simulate(self, filename)
   !< Perform the simulation.
   class(nasto_nvf_object), intent(inout) :: self             !< The equation.
   character(*),            intent(in)    :: filename         !< Input file name.
   real(R8P)                              :: timing(1:2)      !< Tic toc timing.
   real(R8P)                              :: timing_step(1:2) !< Tic toc timing.
   integer(I4P)                           :: i                !< Counter.

   call self%initialize(filename=filename)
   if (self%io%restart) then
      print '(A)', self%mpih%myrankstr//'restart simulation from "'//trim(self%io%restart_basename)//'" files'
      call self%load_restart_files(t=self%time%it, time=self%time%time)
      print '(A)', self%mpih%myrankstr//'restart [t, time]: '//trim(str(self%time%it))//', '//trim(str(self%time%time))
   else
      do i=1, 10
         call self%set_initial_conditions
         if (self%ib%solids_number > 0) call self%update_phi()
         call self%amr_update()
      enddo
      call self%set_initial_conditions
      self%time%time = 0._R8P
      self%time%it = 0
   endif
   if (self%ib%solids_number > 0) call self%update_phi()

   call self%amr_update()

   call self%save_simulation_data

   call self%mpih%barrier(tictoc=.true., timing=timing(1), single=.true.)
   integration: do
      call self%mpih%barrier(tictoc=.true., timing=timing_step(1), single=.true.)
      self%time%it = self%time%it + 1

      if (self%io%save_memory_status) then
         call save_memory_cpu_status(file_name='memory_cpu-'//self%mpih%myrankstr//'.dat', tag=str(self%time%it,.true.))
         call save_memory_gpu_status(file_name='memory_gpu-'//self%mpih%myrankstr//'.dat', tag=str(self%time%it,.true.))
      endif

      if (mod(self%time%it,self%amr%frequency)==0) then
         call self%mpih%barrier(tictoc=.true.)
         call self%amr_update()
         call self%mpih%barrier(tictoc=.true.)
         ! print '(A, F18.10)', self%mpih%myrankstr//'step timing (AMR): ', self%mpih%tictoc_timing()
      endif

      call self%compute_dt()
      if ((self%time%it_max <= 0).and.(self%time%time+self%time%dt > self%time%time_max)) &
         self%time%dt=self%time%time_max-self%time%time

      call self%integrate(t=self%time%time)

      self%time%time = self%time%time + self%time%dt
      call self%time%print_progress(nodes_number=self%adam%tree%nodes_number)

      call self%save_simulation_data

      if (((self%time%it_max <= 0).and.(self%time%time >= self%time%time_max)).or.&
         ((self%time%it>=self%time%it_max).and.(self%time%it_max > 0))) exit integration

      call self%mpih%barrier(tictoc=.true., timing=timing_step(2), single=.true.)
      ! print '(A, F18.10)', self%mpih%myrankstr//'step timing: ', timing_step(2) - timing_step(1)
   enddo integration
   call self%mpih%barrier(tictoc=.true., timing=timing(2), single=.true.)
   ! print '(A, F18.10)', self%mpih%myrankstr//'averaged timing: ', (timing(2) - timing(1))/self%time%it

   call self%save_simulation_data
   endsubroutine simulate
endmodule adam_nasto_nvf_object
