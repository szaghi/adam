!< ADAM, PRISM (Plasma Research usIng Simulation Methods) equations system class definition, GPU (FNL) backend.

#include "fundal.H"

module adam_prism_fnl_object
!< ADAM, PRISM (Plasma Research usIng Simulation Methods) equations system class definition, GPU (FNL) backend.

! ADAM classes, libraries, parameters
use :: adam_common_library
! PRISM modules
use :: adam_prism_common_library
! PRISM FNL modules
use :: adam_prism_fnl_library
! third party modules
use :: fundal, save_memory_status_gpu=>save_memory_status
use :: penf,   save_memory_status_cpu=>save_memory_status
! sdk modules
use :: mpi

implicit none
private
public :: prism_fnl_object

integer(I4P), parameter :: FV_FLUX_VARIANT_MAXWELL          = 0_I4P !< Plain dimensional Maxwell fluxes.
integer(I4P), parameter :: FV_FLUX_VARIANT_MAXWELL_ADIM     = 1_I4P !< Plain adimensional Maxwell fluxes.
integer(I4P), parameter :: FV_FLUX_VARIANT_MAXWELL_DIV_D    = 2_I4P !< Dimensional hyperbolic fluxes with D cleaning.
integer(I4P), parameter :: FV_FLUX_VARIANT_MAXWELL_DIV_B    = 3_I4P !< Dimensional hyperbolic fluxes with B cleaning.
integer(I4P), parameter :: FV_FLUX_VARIANT_MAXWELL_DIV_D_B  = 4_I4P !< Dimensional hyperbolic fluxes with D/B cleaning.
integer(I4P), parameter :: FV_FLUX_VARIANT_MAXWELL_ADIM_DIV_D   = 5_I4P !< Adimensional hyperbolic fluxes with D cleaning.
integer(I4P), parameter :: FV_FLUX_VARIANT_MAXWELL_ADIM_DIV_B   = 6_I4P !< Adimensional hyperbolic fluxes with B cleaning.
integer(I4P), parameter :: FV_FLUX_VARIANT_MAXWELL_ADIM_DIV_D_B = 7_I4P !< Adimensional hyperbolic fluxes with D/B cleaning.
integer(I4P), parameter :: FD_RESIDUAL_VARIANT_PLAIN   = 0_I4P !< Centered-FD Maxwell residual without hyperbolic cleaning.
integer(I4P), parameter :: FD_RESIDUAL_VARIANT_PHI     = 1_I4P !< Centered-FD Maxwell residual with phi cleaning only.
integer(I4P), parameter :: FD_RESIDUAL_VARIANT_PSI     = 2_I4P !< Centered-FD Maxwell residual with psi cleaning only.
integer(I4P), parameter :: FD_RESIDUAL_VARIANT_PHI_PSI = 3_I4P !< Centered-FD Maxwell residual with phi/psi cleaning.
integer(I4P), parameter :: PML_FACE_X_M = 1_I4P
integer(I4P), parameter :: PML_FACE_X_P = 2_I4P
integer(I4P), parameter :: PML_FACE_Y_M = 3_I4P
integer(I4P), parameter :: PML_FACE_Y_P = 4_I4P
integer(I4P), parameter :: PML_FACE_Z_M = 5_I4P
integer(I4P), parameter :: PML_FACE_Z_P = 6_I4P
real(R8P),    parameter :: GRMS_3DB_RATIO = 10.0_R8P**(-3.0_R8P/20.0_R8P)

type, extends(prism_common_object) :: prism_fnl_object
   !< PRISM equations system class definition, GPU (FNL) backend.
   type(field_fnl_object)         :: field_fnl   !< GPU field handler.
   type(ib_fnl_object)            :: ib_fnl      !< GPU immersed boundary.
   type(rk_fnl_object)            :: rk_fnl      !< GPU Runge-Kutta integrator.
   type(weno_fnl_object)          :: weno_fnl    !< GPU WENO reconstructor.
   type(prism_fnl_coil_object)    :: coil_fnl    !< GPU coil source.
   type(prism_fnl_leapfrog_pic_object) :: leapfrog_pic_fnl !< GPU PIC leapfrog integrator.
   type(prism_fnl_pic_object)     :: pic_fnl     !< GPU PIC support state.
   type(prism_fnl_rk_pic_object)  :: rk_pic_fnl  !< GPU PIC RK integrator.
   type(prism_fnl_pml_object)     :: pml_fnl     !< GPU PML compact metadata/state.
   type(prism_fnl_rk_pml_object)  :: rk_pml_fnl  !< GPU PML SSP RK support.
   ! device data
   real(R8P), pointer :: q_gpu(:,:,:,:,:)=>null()           !< Field cell centered variables.
   real(R8P), pointer :: dq_gpu(:,:,:,:,:)=>null()          !< Residuals right hand side.
   real(R8P), pointer :: flxyz_c_gpu(:,:,:,:,:,:,:)=>null() !< Fluxes at cell center with +/- decomposition for all directions.
   real(R8P), pointer :: flx_f_gpu(:,:,:,:,:)=>null()       !< Fluxes along x at cell face.
   real(R8P), pointer :: fly_f_gpu(:,:,:,:,:)=>null()       !< Fluxes along y at cell face.
   real(R8P), pointer :: flz_f_gpu(:,:,:,:,:)=>null()       !< Fluxes along z at cell face.
   real(R8P), pointer :: curl_gpu(:,:,:,:,:)=>null()        !< Curl fields.
   real(R8P), pointer :: divergence_gpu(:,:,:,:,:)=>null()  !< Divergence fields.
   ! host-device data, transposing buffers
   integer(I4P)           :: db5(2,5)              !< Device data bounds (rank 5): bb(1,:)=lower, bb(2,:)=upper.
   integer(I4P)           :: hb5(2,5)              !< Host buffer data bounds (rank 5): bb(1,:)=lower, bb(2,:)=upper.
   real(R8P), allocatable :: buf_5D_R8P(:,:,:,:,:) !< Buffer (host memory, device shape), rank 5, R8P.
   integer(I4P)           :: db6(2,6)              !< Device data bounds (rank 6): bb(1,:)=lower, bb(2,:)=upper.
   integer(I4P)           :: hb6(2,6)              !< Host buffer data bounds (rank 6): bb(1,:)=lower, bb(2,:)=upper.
   real(R8P), allocatable :: buf_6D_R8P(:,:,:,:,:,:) !< Buffer (host memory, device shape), rank 6, R8P.
   integer(I4P)          :: fv_flux_variant    = FV_FLUX_VARIANT_MAXWELL !< FV Maxwell-flux variant selector.
   logical               :: fv_add_phi_damping = .false.                  !< Apply phi damping in FV residuals.
   logical               :: fv_add_psi_damping = .false.                  !< Apply psi damping in FV residuals.
   integer(I4P)          :: fv_ivar_phi        = 0_I4P                    !< Phi slot index in FV residuals.
   integer(I4P)          :: fv_ivar_psi        = 0_I4P                    !< Psi slot index in FV residuals.
   integer(I4P)          :: fd_residual_variant = FD_RESIDUAL_VARIANT_PLAIN !< Centered-FD residual selector.
   integer(I4P)          :: fd_ivar_phi         = 0_I4P                     !< Phi slot index in FD residuals.
   integer(I4P)          :: fd_ivar_psi         = 0_I4P                     !< Psi slot index in FD residuals.
   real(R8P)             :: fd_inv_mu_scale     = 1._R8P / MU0              !< Dimensional/adimensional scaling on curl(B).
   real(R8P)             :: fd_inv_eps_scale    = 1._R8P / EPS0             !< Dimensional/adimensional scaling on curl(D).
   real(R8P)             :: fd_chi_wave         = 0._R8P                    !< Hyperbolic transport speed for phi/psi equations.
   real(R8P)             :: fd_chi_damp         = 0._R8P                    !< Dedner damping speed for phi/psi equations.
   !< Pointer (abstract) TBP.
   procedure(compute_curl_interface_dev),       pass(self),pointer :: compute_curl_dev       =>null()!< Compute curl.
   procedure(compute_derivative1_interface_dev),pass(self),pointer :: compute_derivative1_dev=>null()!< Compute derivative1.
   procedure(compute_derivative2_interface_dev),pass(self),pointer :: compute_derivative2_dev=>null()!< Compute derivative2.
   procedure(compute_derivative4_interface_dev),pass(self),pointer :: compute_derivative4_dev=>null()!< Compute derivative4.
   procedure(compute_divergence_interface_dev), pass(self),pointer :: compute_divergence_dev =>null()!< Compute divergence.
   procedure(compute_gradient_interface_dev),   pass(self),pointer :: compute_gradient_dev   =>null()!< Compute gradient.
   procedure(compute_laplacian_interface_dev),  pass(self),pointer :: compute_laplacian_dev  =>null()!< Compute laplacian.
   procedure(compute_residuals_interface_dev),  pass(self),pointer :: compute_residuals_dev  =>null()!< Compute residuals.
   procedure(integrate_interface_dev),          pass(self),pointer :: integrate_dev          =>null()!< Integrate, time operator.
   contains
      ! auxiliary methods
      procedure, pass(self) :: allocate_gpu     !< Allocate GPU data.
      procedure, pass(self) :: copy_cpu_gpu     !< Copy data from CPU to GPU.
      procedure, pass(self) :: copy_gpu_cpu     !< Copy data from GPU to CPU.
      procedure, pass(self) :: destroy          !< Free device data owned by this realm.
      procedure, pass(self) :: initialize_prism !< Initialize PRISM equation.
      ! IO methods
      procedure, pass(self) :: load_restart_files   !< Load restart files.
      procedure, pass(self) :: save_residuals       !< Save residuals history.
      procedure, pass(self) :: save_simulation_data !< Save all simulation data.
      ! IC/BC/sources
      procedure, pass(self) :: apply_fwl_correction    !< Apply fWLayer correction (if present)
      procedure, pass(self) :: compute_coils_current   !< Compute current coils sources.
      procedure, pass(self) :: verify_no_pic_deposition_on_coils_dev !< Guard against PIC deposition on coil cells.
      procedure, pass(self) :: set_boundary_conditions !< Set boundary conditions of equation.
      procedure, pass(self) :: set_initial_conditions  !< Set initial conditions of equation.
      procedure, pass(self) :: update_ghost            !< Update ghost cells and set boundary conditions.
      procedure, pass(self) :: update_rk_ghost         !< Update RK stage ghost cells.
      ! numerical methods, FDV operators
      procedure, pass(self) :: compute_curl_fd_dev        !< Compute curl of vector field by finite difference.
      procedure, pass(self) :: compute_curl_fv_dev        !< Compute curl of vector field by finite volume.
      procedure, pass(self) :: compute_derivative1_fd_dev !< Compute derivative1 of scalar fields, finite difference schemes.
      procedure, pass(self) :: compute_derivative1_fv_dev !< Compute derivative1 of scalar fields, finite volume schemes.
      procedure, pass(self) :: compute_derivative2_fd_dev !< Compute derivative2 of scalar fields, finite difference schemes.
      procedure, pass(self) :: compute_derivative2_fv_dev !< Compute derivative2 of scalar fields, finite volume schemes.
      procedure, pass(self) :: compute_derivative4_fd_dev !< Compute derivative4 of scalar fields, finite difference schemes.
      procedure, pass(self) :: compute_divergence_fd_dev  !< Compute divergence of vector field by finite difference.
      procedure, pass(self) :: compute_divergence_fv_dev  !< Compute divergence of vector field by finite volume.
      procedure, pass(self) :: compute_gradient_fd_dev    !< Compute gradient of scalar field, finite difference schemes.
      procedure, pass(self) :: compute_gradient_fv_dev    !< Compute gradient of scalar field, finite volume schemes.
      procedure, pass(self) :: compute_laplacian_fd_dev   !< Compute laplacian of scalar field, finite difference schemes.
      procedure, pass(self) :: compute_laplacian_fv_dev   !< Compute laplacian of scalar field, finite volume schemes.
      ! numerical methods, space operators
      procedure, pass(self) :: apply_pml_fd_centered_ade_dev !< Add ADE-PML terms to FD-centered residuals.
      procedure, pass(self) :: compute_residuals_fd_centered_dev !< Compute residuals, centered finite difference schemes.
      ! procedure, pass(self) :: compute_residuals_fv_centered !< Compute residuals, centered finite volume schemes.
      ! procedure, pass(self) :: compute_residuals_weno        !< Compute residuals, WENO schemes.
      ! numerical methods, time operators
      procedure, pass(self) :: integrate_blanesmoan_dev !< Blanes and Moan scheme.
      procedure, pass(self) :: integrate_cfm_dev        !< Commutator-Free Magnus scheme.
      procedure, pass(self) :: integrate_leapfrog_dev   !< Leapfrog scheme.
      procedure, pass(self) :: integrate_rk_ls_dev      !< RK classical low storage schemes.
      procedure, pass(self) :: integrate_rk_ssp_dev     !< SSP RK schemes.
      procedure, pass(self) :: integrate_rk_yoshida_dev !< Yoshida schemes.
      ! forest orchestrator contract methods overridings
      procedure, pass(self) :: initialize_forest                 !< Invoked by forest%initialize per realm at startup.
      procedure, pass(self) :: compute_local_dt_forest           !< Invoked by forest%compute_global_dt during the min reduction.
      procedure, pass(self) :: advance_one_step_forest           !< Invoked by forest%evolve_one_step per realm per timestep.
      procedure, pass(self) :: stages_per_step_forest            !< Number of integrator stages this realm exposes per step.
      procedure, pass(self) :: open_step_forest                  !< Per-step prologue (multi-realm path).
      procedure, pass(self) :: begin_stage_forest                !< Begin one integrator stage on this realm (multi-realm path).
      procedure, pass(self) :: end_stage_forest                  !< End the stage: residuals + assignment (multi-realm path).
      procedure, pass(self) :: close_step_forest                 !< Per-step epilogue (multi-realm path).
      procedure, pass(self) :: post_step_forest                  !< Invoked by forest%post_step per realm per timestep.
      procedure, pass(self) :: is_done_forest                    !< Invoked by forest%is_done during the termination reduction.
      procedure, pass(self) :: finalize_forest                   !< Invoked by forest%finalize per realm at shutdown.
      procedure, pass(self) :: finalize_mpi_forest               !< Process-global MPI finalize (mpih_fnl).
      ! inter-realm seam ghost-fill contract (agnostic-dummy redesign)
      procedure, pass(self) :: fill_seam_from_peer_forest    !< OpenACC device-direct copy of peer's GPU interior into self's GPU ghosts.
      procedure, pass(self) :: after_topology_build_forest   !< Refresh device-resident seam maps from freshly-built host maps.
      procedure, pass(self) :: apply_reflux_to_stage_forest  !< Apply Berger-Colella reflux to self's RK stage buffer (FNL no-op).
      ! numerical methods, miscellanea
      procedure, pass(self) :: compute_dt             !< Compute time step.
      procedure, pass(self) :: compute_energy       	!< Compute energy.
      procedure, pass(self) :: compute_energy_error 	!< Compute energy error.
      procedure, pass(self) :: compute_grms           !< Compute Grms of the rotating magnetic-field amplitude.
      procedure, pass(self) :: compute_magnetic_field_at_center_domain !< Compute B at the domain center.
		procedure, pass(self) :: compute_max_divergence !< Compute divergence of D, B and J fields for diagnostics.
      procedure, pass(self) :: impose_ct_correction_dev  !< Device-side constrained-transport correction on q_gpu.
      procedure, pass(self) :: impose_div_free      	!< Impose divergence-free property.
      procedure, pass(self) :: simulate             	!< Perform the simulation.
endtype prism_fnl_object

interface
   subroutine compute_curl_interface_dev(self, ivar, q_gpu, curl_gpu)
   !< Compute curl of vector fields, div(q(ivar:ivar+2).
   import :: prism_fnl_object, I4P, R8P
   class(prism_fnl_object), intent(in)    :: self                                                !< The equation.
   integer(I4P),            intent(in)    :: ivar                                                !< Start index of variable of q.
   real(R8P),               intent(in)    :: q_gpu(   1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: curl_gpu(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Curl.
   endsubroutine compute_curl_interface_dev

   subroutine compute_derivative1_interface_dev(self, dir, ivar, q_gpu, dq_ds_gpu)
   !< Compute derivative1 of scalar fields, dq(ivar)/ds.
   import :: prism_fnl_object, I4P, R8P
   class(prism_fnl_object), intent(in)    :: self                                              !< The equation.
   integer(I4P),            intent(in)    :: dir                                               !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),            intent(in)    :: ivar                                              !< Start index of variable of q.
   real(R8P),               intent(in)    :: q_gpu(1:, 1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: dq_ds_gpu(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Derivative1, dq/ds.
   endsubroutine compute_derivative1_interface_dev

   subroutine compute_derivative2_interface_dev(self, dir, ivar, q_gpu, d2q_ds2_gpu)
   !< Compute derivative2 of scalar fields, d2q(ivar)/ds2.
   import :: prism_fnl_object, I4P, R8P
   class(prism_fnl_object), intent(in)    :: self                                                !< The equation.
   integer(I4P),            intent(in)    :: dir                                                 !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),            intent(in)    :: ivar                                                !< Start index of variable of q.
   real(R8P),               intent(in)    :: q_gpu(1:,   1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: d2q_ds2_gpu(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Derivative2, d2q/ds2.
   endsubroutine compute_derivative2_interface_dev

   subroutine compute_derivative4_interface_dev(self, dir, ivar, q_gpu, d4q_ds4_gpu)
   !< Compute derivative4 of scalar fields, d4q(ivar)/ds4.
   import :: prism_fnl_object, I4P, R8P
   class(prism_fnl_object), intent(in)    :: self                                                !< The equation.
   integer(I4P),            intent(in)    :: dir                                                 !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),            intent(in)    :: ivar                                                !< Start index of variable of q.
   real(R8P),               intent(in)    :: q_gpu(1:,   1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: d4q_ds4_gpu(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Derivative4, d4q/ds4.
   endsubroutine compute_derivative4_interface_dev

   subroutine compute_divergence_interface_dev(self, ivar, ovar, q_gpu, divergence_gpu, maxdiv)
   !< Compute divergence of vector fields, div(q(ivar:ivar+2).
   import :: prism_fnl_object, I4P, R8P
   class(prism_fnl_object), intent(in)             :: self               !< The equation.
   integer(I4P),            intent(in)             :: ivar               !< Start index of field of q.
   integer(I4P),            intent(in)             :: ovar               !< Output index in divergence.
   real(R8P),               intent(in)             :: q_gpu(1:,         &
                                                            1-self%ngc:,&
                                                            1-self%ngc:,&
                                                            1-self%ngc:,&
                                                            1:)          !< Field variables.
   real(R8P),               intent(inout)          :: divergence_gpu(1:,         &
                                                                     1-self%ngc:,&
                                                                     1-self%ngc:,&
                                                                     1-self%ngc:,&
                                                                     1:) !< Divergence.
   real(R8P),               intent(out), optional  :: maxdiv             !< Maximum divergence for diagnostics.
   endsubroutine compute_divergence_interface_dev

   subroutine compute_gradient_interface_dev(self, ivar, q_gpu, gradient_gpu)
   !< Compute gradient of scalar variable q(ivar).
   import :: prism_fnl_object, I4P, R8P
   class(prism_fnl_object), intent(in)    :: self                                                    !< The equation.
   integer(I4P),            intent(in)    :: ivar                                                    !< Index of scalar var of q.
   real(R8P),               intent(in)    :: q_gpu(       1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: gradient_gpu(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Gradient.
   endsubroutine compute_gradient_interface_dev

   subroutine compute_laplacian_interface_dev(self, ivar, q_gpu, laplacian_gpu)
   !< Compute laplacian of scalar variable q(ivar).
   import :: prism_fnl_object, I4P, R8P
   class(prism_fnl_object), intent(in)    :: self                                                  !< The equation.
   integer(I4P),            intent(in)    :: ivar                                                  !< Index of scalar variable of q.
   real(R8P),               intent(in)    :: q_gpu(     1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: laplacian_gpu(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Laplacian.
   endsubroutine compute_laplacian_interface_dev

   subroutine compute_residuals_interface_dev(self, q_gpu, dq_gpu, s, flux_register)
   !< Compute residuals of equation, space operator.
   !<
   !< Inter-realm seam ghost cells are filled by the forest's Phase 2 seam
   !< exchange BEFORE this routine fires, so the stencil reads valid halo
   !< data without any peer-realm access. The `realm(:)` optional that
   !< used to thread through this signature has been retired by the
   !< agnostic-dummy seam redesign.
   !<
   !< `flux_register` (issue #23 R3): the forest's flux register for FV seam
   !< reflux, threaded from end_stage_forest; consumed by the fv_centered
   !< implementation at the realm's final RK substage, ignored elsewhere.
   import :: prism_fnl_object, R8P, I4P, flux_register_object
   class(prism_fnl_object), intent(inout)             :: self       !< The equation.
   real(R8P),               intent(inout)             :: q_gpu(1:,         &
                                                               1-self%ngc:,&
                                                               1-self%ngc:,&
                                                               1-self%ngc:,&
                                                               1:)  !< Conservative variables.
   real(R8P),               intent(inout)             :: dq_gpu(1:,         &
                                                                1-self%ngc:,&
                                                                1-self%ngc:,&
                                                                1-self%ngc:,&
                                                                1:) !< Residuals.
   integer(I4P),            intent(in),    optional   :: s          !< Stage counter.
   class(flux_register_object), intent(inout), optional :: flux_register !< Forest's flux register for FV seam reflux.
   endsubroutine compute_residuals_interface_dev

   subroutine integrate_interface_dev(self)
   !< Integrate equation, time operator.
   import :: prism_fnl_object, R8P
   class(prism_fnl_object), intent(inout) :: self !< The equation.
   endsubroutine integrate_interface_dev
endinterface

contains
   subroutine destroy(self)
   !< Free device data owned by this PRISM FNL realm.
   class(prism_fnl_object), intent(inout) :: self !< The equation.

   call self%rk_pml_fnl%destroy()
   call self%pml_fnl%destroy()
   call self%rk_pic_fnl%destroy()
   call self%leapfrog_pic_fnl%destroy()
   call self%pic_fnl%destroy()
   call self%coil_fnl%destroy()
   call destroy_rk_fnl_device(self%rk_fnl)
   call destroy_weno_fnl_device(self%weno_fnl)
   call destroy_ib_fnl_device(self%ib_fnl)
   call destroy_field_fnl_device(self%field_fnl)
   call free_prism_core_gpu(self)
   nullify(self%compute_curl_dev)
   nullify(self%compute_derivative1_dev)
   nullify(self%compute_derivative2_dev)
   nullify(self%compute_derivative4_dev)
   nullify(self%compute_divergence_dev)
   nullify(self%compute_gradient_dev)
   nullify(self%compute_laplacian_dev)
   nullify(self%compute_residuals_dev)
   nullify(self%integrate_dev)
   endsubroutine destroy

   subroutine free_prism_core_gpu(self)
   !< Free raw device buffers allocated directly by prism_fnl_object%allocate_gpu.
   class(prism_fnl_object), intent(inout) :: self !< The equation.

   if (associated(self%q_gpu)) then
      call dev_free(self%q_gpu, mydev)
      nullify(self%q_gpu)
   endif
   if (associated(self%dq_gpu)) then
      call dev_free(self%dq_gpu, mydev)
      nullify(self%dq_gpu)
   endif
   if (associated(self%flxyz_c_gpu)) then
      call dev_free(self%flxyz_c_gpu, mydev)
      nullify(self%flxyz_c_gpu)
   endif
   if (associated(self%flx_f_gpu)) then
      call dev_free(self%flx_f_gpu, mydev)
      nullify(self%flx_f_gpu)
   endif
   if (associated(self%fly_f_gpu)) then
      call dev_free(self%fly_f_gpu, mydev)
      nullify(self%fly_f_gpu)
   endif
   if (associated(self%flz_f_gpu)) then
      call dev_free(self%flz_f_gpu, mydev)
      nullify(self%flz_f_gpu)
   endif
   if (associated(self%curl_gpu)) then
      call dev_free(self%curl_gpu, mydev)
      nullify(self%curl_gpu)
   endif
   if (associated(self%divergence_gpu)) then
      call dev_free(self%divergence_gpu, mydev)
      nullify(self%divergence_gpu)
   endif
   if (allocated(self%buf_5D_R8P)) deallocate(self%buf_5D_R8P)
   if (allocated(self%buf_6D_R8P)) deallocate(self%buf_6D_R8P)
   self%db5 = 0_I4P
   self%hb5 = 0_I4P
   self%db6 = 0_I4P
   self%hb6 = 0_I4P
   endsubroutine free_prism_core_gpu

   subroutine destroy_field_fnl_device(field_fnl)
   !< Free device data owned by the field FNL component used by this PRISM realm.
   type(field_fnl_object), intent(inout) :: field_fnl !< FNL field helper.

   call destroy_maps_fnl_device(field_fnl%maps)
   if (associated(field_fnl%fec_1_6_array_gpu)) then
      call dev_free(field_fnl%fec_1_6_array_gpu, mydev)
      nullify(field_fnl%fec_1_6_array_gpu)
   endif
   if (associated(field_fnl%x_cell_gpu)) then
      call dev_free(field_fnl%x_cell_gpu, mydev)
      nullify(field_fnl%x_cell_gpu)
   endif
   if (associated(field_fnl%y_cell_gpu)) then
      call dev_free(field_fnl%y_cell_gpu, mydev)
      nullify(field_fnl%y_cell_gpu)
   endif
   if (associated(field_fnl%z_cell_gpu)) then
      call dev_free(field_fnl%z_cell_gpu, mydev)
      nullify(field_fnl%z_cell_gpu)
   endif
   if (associated(field_fnl%dxyz_gpu)) then
      call dev_free(field_fnl%dxyz_gpu, mydev)
      nullify(field_fnl%dxyz_gpu)
   endif
   nullify(field_fnl%ngc)
   nullify(field_fnl%ni)
   nullify(field_fnl%nj)
   nullify(field_fnl%nk)
   nullify(field_fnl%nb)
   nullify(field_fnl%blocks_number)
   nullify(field_fnl%nv)
   endsubroutine destroy_field_fnl_device

   subroutine destroy_maps_fnl_device(maps_fnl)
   !< Free device data owned by the maps FNL component used by this PRISM realm.
   type(maps_fnl_object), intent(inout) :: maps_fnl !< FNL maps helper.

   if (associated(maps_fnl%local_map_ghost_cell_gpu)) then
      call dev_free(maps_fnl%local_map_ghost_cell_gpu, mydev)
      nullify(maps_fnl%local_map_ghost_cell_gpu)
   endif
   if (associated(maps_fnl%comm_map_recv_ghost_cell_gpu)) then
      call dev_free(maps_fnl%comm_map_recv_ghost_cell_gpu, mydev)
      nullify(maps_fnl%comm_map_recv_ghost_cell_gpu)
   endif
   if (associated(maps_fnl%comm_map_send_ghost_cell_gpu)) then
      call dev_free(maps_fnl%comm_map_send_ghost_cell_gpu, mydev)
      nullify(maps_fnl%comm_map_send_ghost_cell_gpu)
   endif
   if (associated(maps_fnl%send_buffer_ghost_gpu)) then
      call dev_free(maps_fnl%send_buffer_ghost_gpu, mydev)
      nullify(maps_fnl%send_buffer_ghost_gpu)
   endif
   if (associated(maps_fnl%recv_buffer_ghost_gpu)) then
      call dev_free(maps_fnl%recv_buffer_ghost_gpu, mydev)
      nullify(maps_fnl%recv_buffer_ghost_gpu)
   endif
   if (associated(maps_fnl%local_map_bc_crown_gpu)) then
      call dev_free(maps_fnl%local_map_bc_crown_gpu, mydev)
      nullify(maps_fnl%local_map_bc_crown_gpu)
   endif
   if (associated(maps_fnl%seam_local_map_ghost_cell_gpu)) then
      call dev_free(maps_fnl%seam_local_map_ghost_cell_gpu, mydev)
      nullify(maps_fnl%seam_local_map_ghost_cell_gpu)
   endif
   if (associated(maps_fnl%seam_local_send_buf_gpu)) then
      call dev_free(maps_fnl%seam_local_send_buf_gpu, mydev)
      nullify(maps_fnl%seam_local_send_buf_gpu)
   endif
   if (associated(maps_fnl%seam_local_recv_buf_gpu)) then
      call dev_free(maps_fnl%seam_local_recv_buf_gpu, mydev)
      nullify(maps_fnl%seam_local_recv_buf_gpu)
   endif
   if (associated(maps_fnl%seam_comm_map_send_ghost_cell_gpu)) then
      call dev_free(maps_fnl%seam_comm_map_send_ghost_cell_gpu, mydev)
      nullify(maps_fnl%seam_comm_map_send_ghost_cell_gpu)
   endif
   if (associated(maps_fnl%seam_comm_map_recv_ghost_cell_gpu)) then
      call dev_free(maps_fnl%seam_comm_map_recv_ghost_cell_gpu, mydev)
      nullify(maps_fnl%seam_comm_map_recv_ghost_cell_gpu)
   endif
   if (associated(maps_fnl%seam_mpi_send_buf_gpu)) then
      call dev_free(maps_fnl%seam_mpi_send_buf_gpu, mydev)
      nullify(maps_fnl%seam_mpi_send_buf_gpu)
   endif
   if (associated(maps_fnl%seam_mpi_recv_buf_gpu)) then
      call dev_free(maps_fnl%seam_mpi_recv_buf_gpu, mydev)
      nullify(maps_fnl%seam_mpi_recv_buf_gpu)
   endif
   endsubroutine destroy_maps_fnl_device

   subroutine destroy_ib_fnl_device(ib_fnl)
   !< Free device data owned by the IB FNL component used by this PRISM realm.
   type(ib_fnl_object), intent(inout) :: ib_fnl !< FNL IB helper.

   if (associated(ib_fnl%q_bcs_vars_gpu)) then
      call dev_free(ib_fnl%q_bcs_vars_gpu, mydev)
      nullify(ib_fnl%q_bcs_vars_gpu)
   endif
   if (associated(ib_fnl%phi_gpu)) then
      call dev_free(ib_fnl%phi_gpu, mydev)
      nullify(ib_fnl%phi_gpu)
   endif
   endsubroutine destroy_ib_fnl_device

   subroutine destroy_rk_fnl_device(rk_fnl)
   !< Free device data owned by the RK FNL component used by this PRISM realm.
   type(rk_fnl_object), intent(inout) :: rk_fnl !< FNL RK helper.

   if (associated(rk_fnl%alph_gpu)) then
      call dev_free(rk_fnl%alph_gpu, mydev)
      nullify(rk_fnl%alph_gpu)
   endif
   if (associated(rk_fnl%beta_gpu)) then
      call dev_free(rk_fnl%beta_gpu, mydev)
      nullify(rk_fnl%beta_gpu)
   endif
   if (associated(rk_fnl%gamm_gpu)) then
      call dev_free(rk_fnl%gamm_gpu, mydev)
      nullify(rk_fnl%gamm_gpu)
   endif
   if (associated(rk_fnl%q_rk_gpu)) then
      call dev_free(rk_fnl%q_rk_gpu, mydev)
      nullify(rk_fnl%q_rk_gpu)
   endif
   endsubroutine destroy_rk_fnl_device

   subroutine destroy_weno_fnl_device(weno_fnl)
   !< Free device data owned by the WENO FNL component used by this PRISM realm.
   type(weno_fnl_object), intent(inout) :: weno_fnl !< FNL WENO helper.

   if (associated(weno_fnl%a_gpu)) then
      call dev_free(weno_fnl%a_gpu, mydev)
      nullify(weno_fnl%a_gpu)
   endif
   if (associated(weno_fnl%p_gpu)) then
      call dev_free(weno_fnl%p_gpu, mydev)
      nullify(weno_fnl%p_gpu)
   endif
   if (associated(weno_fnl%d_gpu)) then
      call dev_free(weno_fnl%d_gpu, mydev)
      nullify(weno_fnl%d_gpu)
   endif
   if (associated(weno_fnl%ror_schemes_gpu)) then
      call dev_free(weno_fnl%ror_schemes_gpu, mydev)
      nullify(weno_fnl%ror_schemes_gpu)
   endif
   if (associated(weno_fnl%ror_ivar_gpu)) then
      call dev_free(weno_fnl%ror_ivar_gpu, mydev)
      nullify(weno_fnl%ror_ivar_gpu)
   endif
   if (associated(weno_fnl%ror_stats_gpu)) then
      call dev_free(weno_fnl%ror_stats_gpu, mydev)
      nullify(weno_fnl%ror_stats_gpu)
   endif
   if (associated(weno_fnl%cell_scheme_gpu)) then
      call dev_free(weno_fnl%cell_scheme_gpu, mydev)
      nullify(weno_fnl%cell_scheme_gpu)
   endif
   endsubroutine destroy_weno_fnl_device

   ! auxiliary methods
   subroutine allocate_gpu(self)
   !< Allocate GPU data.
   class(prism_fnl_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: ierr !< Error status.
   integer(I4P)                           :: nc   !< Number of coils.
   integer(I4P)                           :: nv_coil !< Number of coil vector components.

   call mpih_fnl%print_message('prism_fnl_object%allocate_gpu start')
   call free_prism_core_gpu(self)
   associate(nv=>self%nv, nb=>self%nb, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk)
   call dev_alloc(fptr_dev=self%q_gpu, &
                  ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,nv], lbounds=[1,1-ngc,1-ngc,1-ngc,1], init_value=0._R8P, ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate q_gpu in prism_fnl_object%allocate_gpu')
   call dev_alloc(fptr_dev=self%dq_gpu, &
                  ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,nv], lbounds=[1,1-ngc,1-ngc,1-ngc,1], init_value=0._R8P, ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate dq_gpu in prism_fnl_object%allocate_gpu')
   call dev_alloc(fptr_dev=self%flxyz_c_gpu, &
                  ubounds=[nb,3,3,ni+ngc,nj+ngc,nk+ngc,nv], lbounds=[1,1,1,1-ngc,1-ngc,1-ngc,1], init_value=0._R8P, ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate flxyz_c_gpu in prism_fnl_object%allocate_gpu')
   call dev_alloc(fptr_dev=self%flx_f_gpu, &
                  ubounds=[nb,ni,nj,nk,nv], lbounds=[1,0,1,1,1], init_value=0._R8P, ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate flx_f_gpu in prism_fnl_object%allocate_gpu')
   call dev_alloc(fptr_dev=self%fly_f_gpu, &
                  ubounds=[nb,ni,nj,nk,nv], lbounds=[1,1,0,1,1], init_value=0._R8P, ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate fly_f_gpu in prism_fnl_object%allocate_gpu')
   call dev_alloc(fptr_dev=self%flz_f_gpu, &
                  ubounds=[nb,ni,nj,nk,nv], lbounds=[1,1,1,0,1], init_value=0._R8P, ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate flz_f_gpu in prism_fnl_object%allocate_gpu')
   call dev_alloc(fptr_dev=self%curl_gpu, &
                  ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,nv], lbounds=[1,1-ngc,1-ngc,1-ngc,1], init_value=0._R8P, ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate curl_gpu in prism_fnl_object%allocate_gpu')
   call dev_alloc(fptr_dev=self%divergence_gpu, &
                  ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,nv], lbounds=[1,1-ngc,1-ngc,1-ngc,1], init_value=0._R8P, ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate divergence_gpu in prism_fnl_object%allocate_gpu')
   ! buffer for transposed copy
   call allocate_variable(var=self%buf_5D_R8P,       &
                          ulb=reshape([1,nb,         &
                                       1-ngc,ni+ngc, &
                                       1-ngc,nj+ngc, &
                                       1-ngc,nk+ngc, &
                                       1,nv],[2,5]), &
                          msg=mpih%myrankstr//'prism_fnl_object%allocate_common(buf_5D_R8P_hh) ', verbose=.true.)
   self%db5(1,:) = [1 ,1-ngc ,1-ngc ,1-ngc ,1 ] ! lower device bounds, rank 5D
   self%db5(2,:) = [nb,ni+ngc,nj+ngc,nk+ngc,nv] ! upper device bounds, rank 5D
   self%hb5(1,:) = [1 ,1-ngc ,1-ngc ,1-ngc ,1 ] ! lower host   bounds, rank 5D
   self%hb5(2,:) = [nv,ni+ngc,nj+ngc,nk+ngc,nb] ! upper host   bounds, rank 5D
   nc = self%coil%total_coils_number
   if (nc > 0_I4P) then
      nv_coil = int(size(self%coil%j_vec, dim=1), I4P)
      call allocate_variable(var=self%buf_6D_R8P,       &
                             ulb=reshape([1,nb,         &
                                          1-ngc,ni+ngc, &
                                          1-ngc,nj+ngc, &
                                          1-ngc,nk+ngc, &
                                          1,nv_coil,    &
                                          1,nc],[2,6]), &
                             msg=mpih%myrankstr//'prism_fnl_object%allocate_common(buf_6D_R8P_hh) ', verbose=.true.)
      self%db6(1,:) = [1 ,1-ngc ,1-ngc ,1-ngc ,1      ,1 ] ! lower device bounds, rank 6D
      self%db6(2,:) = [nb,ni+ngc,nj+ngc,nk+ngc,nv_coil,nc] ! upper device bounds, rank 6D
      self%hb6(1,:) = [1      ,1-ngc ,1-ngc ,1-ngc ,1 ,1 ] ! lower host   bounds, rank 6D
      self%hb6(2,:) = [nv_coil,ni+ngc,nj+ngc,nk+ngc,nb,nc] ! upper host   bounds, rank 6D
   endif
   endassociate
   call mpih_fnl%print_message('prism_fnl_object%allocate_gpu finish')
   endsubroutine allocate_gpu

   subroutine copy_cpu_gpu(self, verbose)
   !< Copy data from CPU to GPU.
   class(prism_fnl_object), intent(inout)        :: self    !< The equation.
   logical,                 intent(in), optional :: verbose !< Trigger verbose output.

   call dev_memcpy_to_device(bb=self%db5,ij=[1,5],tb=self%hb5,dst=self%q_gpu         ,src=self%q         ,buf=self%buf_5D_R8P)
   ! call dev_memcpy_to_device(bb=self%db5,ij=[1,5],tb=self%hb5,dst=self%curl_gpu      ,src=self%curl      ,buf=self%buf_5D_R8P)
   ! call dev_memcpy_to_device(bb=self%db5,ij=[1,5],tb=self%hb5,dst=self%divergence_gpu,src=self%divergence,buf=self%buf_5D_R8P)
   ! call dev_assign_to_device(src=self%q         ,dst=self%q_gpu            ,ij=[1,5])
   ! call dev_assign_to_device(src=self%curl      ,dst=self%curl_gpu         ,ij=[1,5])
   ! call dev_assign_to_device(src=self%divergence,dst=self%divergence_gpu   ,ij=[1,5])
   if (self%physics%physical_model == PIC_PHYSICAL_MODEL) &
      call self%pic_fnl%copy_cpu_gpu(pic=self%pic, q_pic=self%q_pic, pic_fields=self%pic_fields, &
                                     verbose=verbose)
   if (self%physics%physical_model == PIC_PHYSICAL_MODEL .and. self%pic%scheme_time == NUM_SCHEME_TIME_PIC_LEAPFROG) &
      call self%leapfrog_pic_fnl%copy_cpu_gpu(q_pic_old=self%leapfrog_pic%q_pic_old, verbose=verbose)
   if (allocated(self%buf_6D_R8P)) then
      call self%coil_fnl%copy_cpu_gpu(coil=self%coil, grid=self%adam%grid, buf6D=self%buf_6D_R8P, db6=self%db6, hb6=self%hb6)
   else
      call self%coil_fnl%copy_cpu_gpu(coil=self%coil, grid=self%adam%grid)
   endif
   call self%field_fnl%copy_cpu_gpu(field=self%adam%field, maps=self%adam%maps, verbose=verbose)
   endsubroutine copy_cpu_gpu

   subroutine copy_gpu_cpu(self, compute_copy_q_aux, copy_phi, verbose)
   !< Copy data from GPU to CPU.
   class(prism_fnl_object), intent(inout)        :: self               !< The equation.
   logical,                 intent(in), optional :: compute_copy_q_aux !< Flag to compute auxiliary variables.
   logical,                 intent(in), optional :: copy_phi           !< Copy also phi.
   logical,                 intent(in), optional :: verbose            !< Trigger verbose output.

   call dev_memcpy_from_device(bb=self%db5,ij=[1,5],tb=self%hb5,dst=self%q         ,src=self%q_gpu         ,buf=self%buf_5D_R8P)
   call dev_memcpy_from_device(bb=self%db5,ij=[1,5],tb=self%hb5,dst=self%dq        ,src=self%dq_gpu        ,buf=self%buf_5D_R8P)
   call dev_memcpy_from_device(bb=self%db5,ij=[1,5],tb=self%hb5,dst=self%curl      ,src=self%curl_gpu      ,buf=self%buf_5D_R8P)
   call dev_memcpy_from_device(bb=self%db5,ij=[1,5],tb=self%hb5,dst=self%divergence,src=self%divergence_gpu,buf=self%buf_5D_R8P)
   ! call dev_assign_from_device(src=self%q_gpu         ,dst=self%q         ,ij=[1,5])
   ! call dev_assign_from_device(src=self%curl_gpu      ,dst=self%curl      ,ij=[1,5])
   ! call dev_assign_from_device(src=self%divergence_gpu,dst=self%divergence,ij=[1,5])
   if (self%physics%physical_model == PIC_PHYSICAL_MODEL) &
      call self%pic_fnl%copy_gpu_cpu(pic=self%pic, q_pic=self%q_pic, pic_fields=self%pic_fields, &
                                     verbose=verbose)
   if (self%physics%physical_model == PIC_PHYSICAL_MODEL .and. self%pic%scheme_time == NUM_SCHEME_TIME_PIC_LEAPFROG) &
      call self%leapfrog_pic_fnl%copy_gpu_cpu(q_pic_old=self%leapfrog_pic%q_pic_old, verbose=verbose)
   if (allocated(self%buf_6D_R8P)) then
      call self%coil_fnl%copy_gpu_cpu(coil=self%coil, grid=self%adam%grid, buf6D=self%buf_6D_R8P, db6=self%db6, hb6=self%hb6)
   else
      call self%coil_fnl%copy_gpu_cpu(coil=self%coil, grid=self%adam%grid)
   endif
   endsubroutine copy_gpu_cpu

   subroutine initialize_prism(self, filename, realms_number)
   !< Initialize PRISM equation.
   class(prism_fnl_object), intent(inout), target :: self                !< The equation.
   character(*),            intent(in)            :: filename            !< Input file name.
   integer(I4P),            intent(in), optional  :: realms_number
                                                                         !< Forest realm count; divides the per-device budget
                                                                         !< (default 1).
   logical                                        :: is_mpih_initialized !< Flag to check if MPI has been inizialied.
   real(R8P)                                      :: memory_avail_       !< Per-realm device budget (GB) after the forest split.
   integer(I4P)                                   :: realms_number_      !< Local realm count (>=1).

   realms_number_ = 1_I4P ; if (present(realms_number)) realms_number_ = max(realms_number, 1_I4P)
   if (.not. mpih_fnl_is_initialized) then
      call MPI_INITIALIZED(is_mpih_initialized, mpih_fnl%error)
      call mpih_fnl%initialize(do_mpi_init=.not.is_mpih_initialized, do_device_init=.true., verbose=.true.)
      mpih_fnl_is_initialized = .true.
   endif
   call mpih_fnl%print_message('prism_fnl_object%initialize start')
   memory_avail_ = real(mpih_fnl%dev_memory_avail/1e9, R8P) / real(realms_number_, R8P)
   call self%prism_common_object%initialize(filename=filename, memory_avail=memory_avail_, verbose=.true.)
   call check_pml_configuration()
   call self%field_fnl%initialize(grid=self%adam%grid, field=self%adam%field, maps=self%adam%maps, verbose=.true.)
   call self%ib_fnl%initialize(grid=self%adam%grid, field=self%adam%field, ib=self%ib)
   call self%rk_fnl%initialize(grid=self%adam%grid, field=self%adam%field, rk=self%rk)
   call self%weno_fnl%initialize(weno=self%weno)
   call self%pml_fnl%initialize(pml=self%pml, grid=self%adam%grid, field=self%adam%field)
   if (self%pml_fnl%enabled) call self%rk_pml_fnl%initialize(rk=self%rk, pml_fnl=self%pml_fnl)
   call self%allocate_gpu
   if (allocated(self%buf_6D_R8P)) then
      call self%coil_fnl%initialize(coil=self%coil, field=self%adam%field, grid=self%adam%grid, &
                                    buf6D=self%buf_6D_R8P, db6=self%db6, hb6=self%hb6)
   else
      call self%coil_fnl%initialize(coil=self%coil, field=self%adam%field, grid=self%adam%grid)
   endif
   if (self%physics%physical_model == PIC_PHYSICAL_MODEL) then
      call self%pic_fnl%initialize(pic=self%pic, q_pic=self%q_pic, pic_fields=self%pic_fields)
      if (self%pic%scheme_time == NUM_SCHEME_TIME_PIC_LEAPFROG) &
         call self%leapfrog_pic_fnl%initialize(pic=self%pic, leapfrog_pic=self%leapfrog_pic)
      if (self%pic%scheme_time == NUM_SCHEME_TIME_PIC_RUNGE_KUTTA) &
         call self%rk_pic_fnl%initialize(pic=self%pic, rk_pic=self%rk_pic)
   endif

   ! set pointer (abstract) TBP
   if (self%physics%physical_model == EM_PHYSICAL_MODEL .or. self%physics%physical_model == ADIM_EM_PHYSICAL_MODEL) then
      select case(self%numerics%scheme_time)
      case(NUM_SCHEME_TIME_BLANES_MOAN)        ; self%integrate_dev => integrate_blanesmoan_dev
      case(NUM_SCHEME_TIME_CFM)                ; self%integrate_dev => integrate_cfm_dev
      case(NUM_SCHEME_TIME_LEAPFROG)           ; self%integrate_dev => integrate_leapfrog_dev
      case(NUM_SCHEME_TIME_RUNGE_KUTTA)
         select case(self%rk%scheme)
         case(RK_1, RK_2, RK_3)                ; self%integrate_dev => integrate_rk_ls_dev
         case(RK_SSP_11, RK_SSP_22, RK_SSP_33, RK_SSP_54) ; self%integrate_dev => integrate_rk_ssp_dev
         case(RK_YOSHIDA)                      ; self%integrate_dev => integrate_rk_yoshida_dev
         endselect
      endselect
   elseif (self%physics%physical_model == PIC_PHYSICAL_MODEL) then
      select case(self%numerics%scheme_time)
      case(NUM_SCHEME_TIME_LEAPFROG)
         select case(self%pic%scheme_time)
         case(NUM_SCHEME_TIME_PIC_LEAPFROG)
            self%integrate_dev => integrate_leapfrog_pic
         case default
            call mpih_fnl%error_stop(msg=': PIC time integration combination not ported to FNL backend')
         endselect
      case(NUM_SCHEME_TIME_RUNGE_KUTTA)
         select case(self%pic%scheme_time)
         case(NUM_SCHEME_TIME_PIC_LEAPFROG)
            self%integrate_dev => integrate_leapfrog_pic
         case(NUM_SCHEME_TIME_PIC_RUNGE_KUTTA)
            select case(self%rk_pic%scheme)
            case(RK_SSP_22, RK_SSP_33, RK_SSP_54)
               self%integrate_dev => integrate_rk_ssp_pic
            case default
               call mpih_fnl%error_stop(msg=': PIC RK scheme not ported to FNL backend')
            endselect
         case default
            call mpih_fnl%error_stop(msg=': PIC time integration combination not ported to FNL backend')
         endselect
      case default
         call mpih_fnl%error_stop(msg=': PIC time integration combination not ported to FNL backend')
      endselect
   endif

   self%fv_flux_variant    = FV_FLUX_VARIANT_MAXWELL
   self%fv_add_phi_damping = .false.
   self%fv_add_psi_damping = .false.
   self%fv_ivar_phi        = 0_I4P
   self%fv_ivar_psi        = 0_I4P
   select case(self%physics%physical_model)
   case(ADIM_EM_PHYSICAL_MODEL)
      select case(self%numerics%div_corr_var)
      case(DIV_CORR_VAR_HYPER)
         if (self%numerics%constrained_transport_D .and. .not.self%numerics%constrained_transport_B) then
            self%fv_flux_variant    = FV_FLUX_VARIANT_MAXWELL_ADIM_DIV_D
            self%fv_add_phi_damping = .true.
            self%fv_ivar_phi        = self%nv_c
         elseif (.not.self%numerics%constrained_transport_D .and. self%numerics%constrained_transport_B) then
            self%fv_flux_variant    = FV_FLUX_VARIANT_MAXWELL_ADIM_DIV_B
            self%fv_add_psi_damping = .true.
            self%fv_ivar_psi        = self%nv_c
         elseif (self%numerics%constrained_transport_D .and. self%numerics%constrained_transport_B) then
            self%fv_flux_variant    = FV_FLUX_VARIANT_MAXWELL_ADIM_DIV_D_B
            self%fv_add_phi_damping = .true.
            self%fv_add_psi_damping = .true.
            self%fv_ivar_phi        = self%nv_c - 1_I4P
            self%fv_ivar_psi        = self%nv_c
         else
            self%fv_flux_variant = FV_FLUX_VARIANT_MAXWELL_ADIM
         endif
      case default
         self%fv_flux_variant = FV_FLUX_VARIANT_MAXWELL_ADIM
      endselect
   case default
      select case(self%numerics%div_corr_var)
      case(DIV_CORR_VAR_HYPER)
         if (self%numerics%constrained_transport_D .and. .not.self%numerics%constrained_transport_B) then
            self%fv_flux_variant    = FV_FLUX_VARIANT_MAXWELL_DIV_D
            self%fv_add_phi_damping = .true.
            self%fv_ivar_phi        = self%nv_c
         elseif (.not.self%numerics%constrained_transport_D .and. self%numerics%constrained_transport_B) then
            self%fv_flux_variant    = FV_FLUX_VARIANT_MAXWELL_DIV_B
            self%fv_add_psi_damping = .true.
            self%fv_ivar_psi        = self%nv_c
         elseif (self%numerics%constrained_transport_D .and. self%numerics%constrained_transport_B) then
            self%fv_flux_variant    = FV_FLUX_VARIANT_MAXWELL_DIV_D_B
            self%fv_add_phi_damping = .true.
            self%fv_add_psi_damping = .true.
            self%fv_ivar_phi        = self%nv_c - 1_I4P
            self%fv_ivar_psi        = self%nv_c
         else
            self%fv_flux_variant = FV_FLUX_VARIANT_MAXWELL
         endif
      case default
         self%fv_flux_variant = FV_FLUX_VARIANT_MAXWELL
      endselect
   endselect

   self%fd_residual_variant = FD_RESIDUAL_VARIANT_PLAIN
   self%fd_ivar_phi         = 0_I4P
   self%fd_ivar_psi         = 0_I4P
   select case(self%physics%physical_model)
   case(ADIM_EM_PHYSICAL_MODEL)
      self%fd_inv_mu_scale  = 1._R8P
      self%fd_inv_eps_scale = 1._R8P
      self%fd_chi_wave      = self%physics%chi
      self%fd_chi_damp      = self%physics%chi
   case default
      self%fd_inv_mu_scale  = 1._R8P / MU0
      self%fd_inv_eps_scale = 1._R8P / EPS0
      self%fd_chi_wave      = self%physics%chi * C0
      self%fd_chi_damp      = self%physics%chi * C0
   endselect
   if (self%fv_add_phi_damping .and. self%fv_add_psi_damping) then
      self%fd_residual_variant = FD_RESIDUAL_VARIANT_PHI_PSI
      self%fd_ivar_phi         = self%fv_ivar_phi
      self%fd_ivar_psi         = self%fv_ivar_psi
   elseif (self%fv_add_phi_damping) then
      self%fd_residual_variant = FD_RESIDUAL_VARIANT_PHI
      self%fd_ivar_phi         = self%fv_ivar_phi
   elseif (self%fv_add_psi_damping) then
      self%fd_residual_variant = FD_RESIDUAL_VARIANT_PSI
      self%fd_ivar_psi         = self%fv_ivar_psi
   endif

   select case(self%numerics%scheme_space)
   case(NUM_SCHEME_SPACE_WENO)
      self%compute_curl_dev        => compute_curl_fv_dev
      self%compute_derivative1_dev => compute_derivative1_fv_dev
      self%compute_derivative2_dev => compute_derivative2_fv_dev
      self%compute_divergence_dev  => compute_divergence_fv_dev
      self%compute_gradient_dev    => compute_gradient_fv_dev
      self%compute_laplacian_dev   => compute_laplacian_fv_dev
      ! self%compute_residuals_dev   => compute_residuals_weno_dev
      ! issue #23 R1: the WENO residual path is not ported to FNL — refuse cleanly at
      ! initialization instead of leaving a null procedure pointer (segfault at step 1).
      call mpih_fnl%error_stop(msg=': scheme_space "weno" residual path is not ported to FNL — '// &
                                   'run this case on the CPU backend')
   case(NUM_SCHEME_SPACE_FD_CENTERED)
      self%compute_curl_dev        => compute_curl_fd_dev
      self%compute_derivative1_dev => compute_derivative1_fd_dev
      self%compute_derivative2_dev => compute_derivative2_fd_dev
      self%compute_divergence_dev  => compute_divergence_fd_dev
      self%compute_gradient_dev    => compute_gradient_fd_dev
      self%compute_laplacian_dev   => compute_laplacian_fd_dev
      self%compute_residuals_dev   => compute_residuals_fd_centered_dev
   case(NUM_SCHEME_SPACE_FV_CENTERED)
      self%compute_curl_dev        => compute_curl_fv_dev
      self%compute_derivative1_dev => compute_derivative1_fv_dev
      self%compute_derivative2_dev => compute_derivative2_fv_dev
      self%compute_divergence_dev  => compute_divergence_fv_dev
      self%compute_gradient_dev    => compute_gradient_fv_dev
      self%compute_laplacian_dev   => compute_laplacian_fv_dev
      self%compute_residuals_dev   => compute_residuals_fv_centered_dev
   endselect

   call external_fields_initialize_dev(external_fields=self%external_fields)

   call mpih_fnl%print_message('prism_fnl_object%initialize finish')
   contains
      subroutine check_pml_configuration()
      if (.not. self%pml%enabled) return
      if (self%numerics%scheme_space /= NUM_SCHEME_SPACE_FD_CENTERED) &
         call mpih_fnl%error_stop(msg=': FNL PML is currently available only with scheme_space = FD_centered')
      if (self%numerics%scheme_time /= NUM_SCHEME_TIME_RUNGE_KUTTA) &
         call mpih_fnl%error_stop(msg=': FNL PML is currently available only with scheme_time = Runge_Kutta')
      select case (trim(self%rk%scheme))
      case (RK_SSP_11, RK_SSP_22, RK_SSP_33, RK_SSP_54)
         continue
      case default
         call mpih_fnl%error_stop(msg=': FNL PML is currently available only with SSP RK schemes')
      endselect
      if (trim(self%numerics%div_corr_var) == DIV_CORR_VAR_HYPER .or. trim(self%numerics%div_corr_var) == DIV_CORR_VAR_POISS) &
         call mpih_fnl%error_stop(msg=': FNL PML is currently available only without divergence correction')
      if (self%numerics%constrained_transport_D .or. self%numerics%constrained_transport_B) &
         call mpih_fnl%error_stop(msg=': FNL PML is currently available only without constrained transport')
      endsubroutine check_pml_configuration
   endsubroutine initialize_prism

   ! IO methods
   subroutine load_restart_files(self, t, time)
   !< Save restart files.
   class(prism_fnl_object), intent(inout) :: self !< The equation.
   integer(I4P),            intent(out)   :: t    !< Time iteration.
   real(R8P),               intent(out)   :: time !< Time.

   call self%prism_common_object%load_restart_files(t=t, time=time)
   call self%copy_cpu_gpu
   endsubroutine load_restart_files

   subroutine save_residuals(self)
   !< Save residuals history.
   class(prism_fnl_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: v    !< Counter.

   if (self%time%is_to_save(it_save=self%io%residuals_save)) then
      call compute_normL2_residuals_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, nv=self%nv, &
                                        blocks_number=self%blocks_number, dq_gpu=self%dq_gpu, norm=self%adam%field%residuals)
      do v=1, self%nv
         call MPI_ALLREDUCE(MPI_IN_PLACE, self%adam%field%residuals(v), 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih_fnl%error)
         self%adam%field%residuals(v) = sqrt(self%adam%field%residuals(v))/sqrt(real(self%ni*self%nj*self%nk, R8P))
      enddo
      if (mpih_fnl%myrank==0) call self%io%save_residuals(it=self%time%it, time=self%time%time, &
                                                          blocks_number=self%blocks_number,     &
                                                          residuals=self%adam%field%residuals)
   endif
   endsubroutine save_residuals

   subroutine save_simulation_data(self)
   !< Save all simulation data.
   class(prism_fnl_object), intent(inout) :: self !< The equation.

   if ((self%time%is_to_save(it_save=self%io%it_save)).or.      &
       (self%time%is_to_save(it_save=self%io%restart_save)).or. &
       (self%slices%is_to_save(it=self%time%it,it_max=self%time%it_max,time=self%time%time,time_max=self%time%time_max))) then
      call self%update_ghost(q_gpu=self%q_gpu)
      call self%copy_gpu_cpu
      call self%compute_auxiliary_fields

      if (self%time%is_to_save(it_save=self%io%it_save)) call self%save_xh5f(with_ghost=.true.)
      if (mod(self%time%it,self%io%restart_save)==0) call self%save_restart_files
      ! if (self%slices%is_to_save(it=self%time%it,it_max=self%time%it_max,time=self%time%time,time_max=self%time%time_max)) then
      !    call self%slices%save_mat(basename=self%io%output_basename, &
      !                              it=self%time%it,                  &
      !                              it_max=self%time%it_max,          &
      !                              time=self%time%time,              &
      !                              time_max=self%time%time_max,      &
      !                              adam=self%adam,                   &
      !                              q=self%q,                         &
      !                              q_name=self%q_name)
      ! endif
   endif
   if (self%pic%problem_type == SINGLE_PARTICLE_TYPE_PROBLEM) then
      call self%pic_fnl%copy_q_pic_gpu_cpu(q_pic=self%q_pic)
      call write_single_particle_output(filename='single_particle_output.dat', time=self%time%time, q_pic=self%q_pic)
   endif
   endsubroutine save_simulation_data

   ! IC/BC/sources
   subroutine apply_fwl_correction(self, q_gpu)
   !< Apply correction if a fWL is present.
   class(prism_fnl_object), intent(inout) :: self                                !< The equation.
   real(R8P),               intent(inout) :: q_gpu(1:,         &
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1:)                           !< Conservative variables.
   integer(I4P)                           :: face, b                           !< Counter.

   associate(layer=>self%fWLayer%layer, C=>self%fWLayer%C, ni_fWL=>self%fWLayer%ni_fWL, &
            nj_fWL=>self%fWLayer%nj_fWL, nk_fWL=>self%fWLayer%nk_fWL, n=>self%fWLayer%n, &
            s2=>self%fWLayer%s2, alfa_D=>self%fWLayer%alfa_D, beta_D=>self%fWLayer%beta_D, &
            alfa_B=>self%fWLayer%alfa_B, beta_B=>self%fWLayer%beta_B, profile_extent=>self%fWLayer%profile_extent, &
            profile_cells=>self%fWLayer%profile_cells, dxyz_gpu=>self%field_fnl%dxyz_gpu, x_cell_gpu=>self%field_fnl%x_cell_gpu, &
            y_cell_gpu=>self%field_fnl%y_cell_gpu, z_cell_gpu=>self%field_fnl%z_cell_gpu)
   if (allocated(self%fWLayer%C)) then
      do face=1, 6
         if (.not. layer(face)) cycle
         do b=1, self%blocks_number
            if (C(b,face) <= 0_I4P) cycle
            call apply_fwl_correction_dev_kernel(block_idx=b, ngc=self%ngc, ni=self%ni, nj=self%nj, nk=self%nk, &
                                                 ni1=ni_fWL(1,b,face), ni2=ni_fWL(2,b,face), &
                                                 nj1=nj_fWL(1,b,face), nj2=nj_fWL(2,b,face), &
                                                 nk1=nk_fWL(1,b,face), nk2=nk_fWL(2,b,face), &
                                                 face=face,                                    &
                                                 n=n(face), s2=s2(face), alfa_D=alfa_D(face), beta_D=beta_D(face), &
                                                 alfa_B=alfa_B(face), beta_B=beta_B(face),    &
                                                 domain_emin_n=self%adam%grid%domain_emin(n(face)),                  &
                                                 domain_emax_n=self%adam%grid%domain_emax(n(face)),                  &
                                                 profile_extent=profile_extent(face), profile_cells=profile_cells(face), &
                                                 x_cell_gpu=x_cell_gpu, y_cell_gpu=y_cell_gpu, z_cell_gpu=z_cell_gpu, &
                                                 dxyz_gpu=dxyz_gpu, q_gpu=q_gpu)
         enddo
      enddo
   endif
   endassociate
   endsubroutine apply_fwl_correction

   subroutine compute_coils_current(self, q_gpu, gamm)
   !< Compute current coils sources.
   class(prism_fnl_object), intent(in)           :: self            !< The equation.
   real(R8P),               intent(inout)        :: q_gpu(1:,         &
                                                          1-self%ngc:,&
                                                          1-self%ngc:,&
                                                          1-self%ngc:,&
                                                          1:)       !< Conservative variables.
   real(R8P),               intent(in), optional :: gamm            !< RK coefficient.
   real(R8P)                                     :: time_s          !< Local time.
   real(R8P)                                     :: current_density !< Current density.
   real(R8P)                                     :: g               !< Starting polynomial transitory of coils.
   real(R8P)                                     :: s               !< Non dimensional time, clamp.
   integer(I4P)                                  :: coil_id         !< Uniq coild ID.
   integer(I4P)                                  :: i,j,k,b         !< Counter.
   real(R8P)                                     :: phi_rad         !< Phase in rads.
   real(R8P)                                     :: omega           !< Frequency.
   real(R8P)                                     :: theta           !< Phase.
   real(R8P)                                     :: f_abs           !< Absolute frequency.
   integer(I4P)                                  :: w_ac            !< Switch AC/DC, =1 AC, =0 DC (branchless).
   integer(I4P)                                  :: n               !< Coils number counter.
   real(R8P), parameter                          :: f_tol=1.e-30_R8P!< Tolerance on frequency for branchless switch AC/DC.

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, blocks_number=>self%blocks_number,    &
             time=>self%time%time, dt=>self%time%dt, td=>self%coil%td,                                                  &
             A=>self%coil%coil_amplitude, f=>self%coil%f, phase=>self%coil%phase,                                       &
             var_Jx=>self%physics%var_Jx, var_Jy=>self%physics%var_Jy, var_Jz=>self%physics%var_Jz)
   time_s = time ; if (present(gamm)) time_s = time + dt*gamm
   if (self%coil%total_coils_number >= 1_I4P) then
      ! Azzero termini sorgenti (NB: col PIC potresti voler accumulare in un buffer)
      call nullify_j_vec_vars_kernel(ni            = ni           ,&
                                     nj            = nj           ,&
                                     nk            = nk           ,&
                                     ngc           = ngc          ,&
                                     blocks_number = blocks_number,&
                                     coils_number  = self%coil%total_coils_number,&
                                     var_jx        = var_jx       ,&
                                     var_jy        = var_jy       ,&
                                     var_jz        = var_jz       ,&
                                     j_vec_gpu     = self%coil_fnl%j_vec_gpu,&
                                     q_gpu         = q_gpu)

      ! Envelope C^2: clamp(s) in [0,1], g(0)=0, g(1)=1, g'(0)=g'(1)=0, g''(0)=g''(1)=0
      s = 1._R8P ; if (td > 0._R8P) s = time_s / td ; s = max(0._R8P, min(1._R8P, s))
      g = 10._R8P*s**3 - 15._R8P*s**4 + 6._R8P*s**5

      do n=1, self%coil%total_coils_number
         coil_id = n

         phi_rad = phase(coil_id) * PI / 180._R8P
         omega   = 2._R8P * PI * f(coil_id)

         ! Se f ~ 0 -> DC (w_ac=0). Se f != 0 -> AC (w_ac=1). Branchless.
         f_abs = abs(f(coil_id))
         w_ac  = nint( (sign(1._R8P, f_abs - f_tol) + 1._R8P) * 0.5_R8P )

         ! Theta: DC -> theta = phi ; AC -> theta = omega*(t-td) + phi
         theta = w_ac * omega * (time_s - td) + phi_rad

         ! Unica formula: DC e AC
         current_density = A(coil_id) * g * cos(theta)

         call apply_j_vec_kernel(ni              = ni                     ,&
                                 nj              = nj                     ,&
                                 nk              = nk                     ,&
                                 ngc             = ngc                    ,&
                                 blocks_number   = blocks_number          ,&
                                 current_density = current_density        ,&
                                 n               = n                      ,&
                                 var_jx          = var_jx                 ,&
                                 var_jy          = var_jy                 ,&
                                 var_jz          = var_jz                 ,&
                                 j_vec_gpu       = self%coil_fnl%j_vec_gpu,&
                                 q_gpu           = q_gpu)
      enddo
   endif
   endassociate
   contains
      subroutine nullify_j_vec_vars_kernel(ni,nj,nk,ngc,blocks_number,coils_number,var_jx,var_jy,var_jz,j_vec_gpu,q_gpu)
      !< Nullify J_Vec vars in q only on the support of the analytic coil current.
      integer(I4P), intent(in)    :: ni,nj,nk,ngc,blocks_number,coils_number        !< Grids dimensions / coil count.
      integer(I4P), intent(in)    :: var_jx,var_jy,var_jz                            !< Indexes of J_vec variables.
      real(R8P),    intent(in)    :: j_vec_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:,1:)       !< Analytic coil support.
      real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)               !< Field cell centered variables.
      integer(I4P)                :: i,j,k,b,n                                       !< Counter.
      logical                     :: has_coil                                         !< Coil support flag.

      ! FULL ghost-inclusive extent (issue #26 G2): the CPU stamp covers 1-ngc..n+ngc;
      ! the former interior-only loops left the q J-row ghosts UNSTAMPED (stale exchange
      ! content) -- the divergence stencils near every block border read wrong J, the
      ! dominant term of the FNL div(J) untruthfulness (2.35E+01 vs 3.2E-05).
      !$acc parallel loop independent gang vector collapse(4) &
      !$acc DEVICEVAR(q_gpu,j_vec_gpu)                        &
      !$acc firstprivate(ni,nj,nk,ngc,blocks_number,coils_number,var_jx,var_jy,var_jz)
      !$omp OMPLOOP collapse(4) &
      !$omp DEVICEPTR(q_gpu,j_vec_gpu) &
      !$omp firstprivate(ni,nj,nk,ngc,blocks_number,coils_number,var_jx,var_jy,var_jz)
      do b=1, blocks_number
      do k=1-ngc, nk+ngc
      do j=1-ngc, nj+ngc
      do i=1-ngc, ni+ngc
         has_coil = .false.
         !$acc loop seq
         do n=1, coils_number
            has_coil = has_coil .or. any(j_vec_gpu(b,i,j,k,1:3,n) /= 0._R8P)
         enddo
         if (has_coil) then
            q_gpu(b,i,j,k,var_Jx) = 0._R8P
            q_gpu(b,i,j,k,var_Jy) = 0._R8P
            q_gpu(b,i,j,k,var_Jz) = 0._R8P
         endif
      enddo
      enddo
      enddo
      enddo
      endsubroutine nullify_j_vec_vars_kernel

      subroutine apply_j_vec_kernel(ni,nj,nk,ngc,blocks_number,current_density,n,var_jx,var_jy,var_jz,j_vec_gpu,q_gpu)
      !< Apply J_Vec to q, device kernel.
      integer(I4P), intent(in)    :: ni,nj,nk,ngc,blocks_number               !< Grids dimensions.
      real(R8P),    intent(in)    :: current_density                          !< Current density.
      integer(I4P), intent(in)    :: n                                        !< Current coil index.
      integer(I4P), intent(in)    :: var_jx,var_jy,var_jz                     !< Indexes of J_vec variables.
      real(R8P),    intent(in)    :: j_vec_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:,1:) !< Field cell centered variables.
      real(R8P),    intent(inout) :: q_gpu(    1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Field cell centered variables.
      integer(I4P)                :: i,j,k,b                                  !< Counter.

      ! FULL ghost-inclusive extent (issue #26 G2): twin of the nullify kernel above;
      ! j_vec_gpu carries ghost cells by allocation and matches the CPU J_vec content
      ! to round-off, so the analytic stamp is valid over the whole slab.
      !$acc parallel loop independent gang vector collapse(4) &
      !$acc DEVICEVAR(q_gpu,j_vec_gpu)                        &
      !$acc firstprivate(ni,nj,nk,ngc,blocks_number,current_density,n,var_jx,var_jy,var_jz)
      !$omp OMPLOOP collapse(4) &
      !$omp DEVICEPTR(q_gpu,j_vec_gpu) &
      !$omp firstprivate(ni,nj,nk,ngc,blocks_number,current_density,n,var_jx,var_jy,var_jz)
      do b=1, blocks_number
      do k=1-ngc, nk+ngc
      do j=1-ngc, nj+ngc
      do i=1-ngc, ni+ngc
         q_gpu(b,i,j,k,var_Jx) = q_gpu(b,i,j,k,var_Jx) + current_density * j_vec_gpu(b,i,j,k,1,n)
         q_gpu(b,i,j,k,var_Jy) = q_gpu(b,i,j,k,var_Jy) + current_density * j_vec_gpu(b,i,j,k,2,n)
         q_gpu(b,i,j,k,var_Jz) = q_gpu(b,i,j,k,var_Jz) + current_density * j_vec_gpu(b,i,j,k,3,n)
      enddo
      enddo
      enddo
      enddo
      endsubroutine apply_j_vec_kernel
   endsubroutine compute_coils_current

   subroutine verify_no_pic_deposition_on_coils_dev(self, q_gpu, check_current, check_charge, context)
   !< Ensure device-side PIC deposition does not populate cells carrying analytic coil current.
   class(prism_fnl_object), intent(in)              :: self                                               !< The equation.
   real(R8P),               intent(in)              :: q_gpu(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Deposited field.
   logical,                 intent(in),    optional :: check_current                                      !< Check Jx/Jy/Jz overlap.
   logical,                 intent(in),    optional :: check_charge                                       !< Check rho overlap.
   character(*),            intent(in),    optional :: context                                            !< Call-site label.
   logical                                         :: do_current                                          !< Check current overlap.
   logical                                         :: do_charge                                           !< Check charge overlap.
   character(len=:), allocatable                   :: context_                                            !< Context label, local copy.
   integer(I4P)                                    :: overlap_current                                     !< Reduction flag for current overlap.
   integer(I4P)                                    :: overlap_charge                                      !< Reduction flag for charge overlap.
   integer(I4P)                                    :: nv_q                                                !< Runtime q width.

   do_current = .false. ; if (present(check_current)) do_current = check_current
   do_charge  = .false. ; if (present(check_charge )) do_charge  = check_charge
   if ((.not. do_current) .and. (.not. do_charge)) return
   if (self%coil%total_coils_number <= 0_I4P) return

   context_ = 'PIC deposition'
   if (present(context)) context_ = trim(context)

   overlap_current = 0_I4P
   overlap_charge  = 0_I4P
   nv_q = int(size(q_gpu, dim=5), I4P)

   if (do_current .or. do_charge) then
      call verify_no_pic_deposition_on_coils_dev_kernel(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc,          &
                                                        blocks_number=self%blocks_number, coils_number=self%coil%total_coils_number, &
                                                        var_jx=self%physics%var_Jx, var_jy=self%physics%var_Jy, var_jz=self%physics%var_Jz, &
                                                        nv_q=nv_q, do_current=merge(1_I4P,0_I4P,do_current),         &
                                                        do_charge=merge(1_I4P,0_I4P,do_charge), j_vec_gpu=self%coil_fnl%j_vec_gpu, &
                                                        q_gpu=q_gpu, overlap_current=overlap_current, overlap_charge=overlap_charge)
   endif

   if (overlap_current /= 0_I4P) then
      call mpih_fnl%error_stop(msg=': '//trim(context_)//' deposited plasma current on a coil cell')
   endif
   if (overlap_charge /= 0_I4P) then
      call mpih_fnl%error_stop(msg=': '//trim(context_)//' deposited plasma charge on a coil cell')
   endif

   contains
      subroutine verify_no_pic_deposition_on_coils_dev_kernel(ni, nj, nk, ngc, blocks_number, coils_number,       &
                                                              var_jx, var_jy, var_jz, nv_q, do_current, do_charge, &
                                                              j_vec_gpu, q_gpu, overlap_current, overlap_charge)
      !< Check overlap between PIC deposition and coil support on device.
      integer(I4P), intent(in)    :: ni, nj, nk, ngc, blocks_number, coils_number       !< Grid / coil sizes.
      integer(I4P), intent(in)    :: var_jx, var_jy, var_jz, nv_q                        !< Variable indexes.
      integer(I4P), intent(in)    :: do_current, do_charge                               !< Integerized logicals.
      real(R8P),    intent(in)    :: j_vec_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:,1:)           !< Analytic coil support.
      real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)                  !< Deposited field.
      integer(I4P), intent(inout) :: overlap_current, overlap_charge                     !< Host reduction flags.
      integer(I4P)                :: b, i, j, k, n                                       !< Counters.
      logical                     :: has_coil                                             !< Coil support at cell.

      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(j_vec_gpu, q_gpu) &
      !$acc& firstprivate(ni, nj, nk, ngc, blocks_number, coils_number, var_jx, var_jy, var_jz, nv_q, do_current, do_charge) &
      !$acc& reduction(max: overlap_current, overlap_charge)
      !$omp OMPLOOP collapse(4) DEVICEPTR(j_vec_gpu, q_gpu) &
      !$omp& firstprivate(ni, nj, nk, ngc, blocks_number, coils_number, var_jx, var_jy, var_jz, nv_q, do_current, do_charge) &
      !$omp& reduction(max: overlap_current, overlap_charge)
      do b=1, blocks_number
      do k=1-ngc, nk+ngc
      do j=1-ngc, nj+ngc
      do i=1-ngc, ni+ngc
         has_coil = .false.
         !$acc loop seq
         do n=1, coils_number
            has_coil = has_coil .or. any(j_vec_gpu(b,i,j,k,1:3,n) /= 0._R8P)
         enddo
         if (.not. has_coil) cycle
         if (do_current /= 0_I4P) then
            if ((q_gpu(b,i,j,k,var_jx) /= 0._R8P) .or. (q_gpu(b,i,j,k,var_jy) /= 0._R8P) .or. &
                (q_gpu(b,i,j,k,var_jz) /= 0._R8P)) overlap_current = 1_I4P
         endif
         if (do_charge /= 0_I4P) then
            if (q_gpu(b,i,j,k,nv_q) /= 0._R8P) overlap_charge = 1_I4P
         endif
      enddo
      enddo
      enddo
      enddo
      endsubroutine verify_no_pic_deposition_on_coils_dev_kernel
   endsubroutine verify_no_pic_deposition_on_coils_dev

   subroutine set_boundary_conditions(self, q_gpu)
   !< Set boundary conditions of equation.
   class(prism_fnl_object), intent(in)    :: self                  !< The equation.
   real(R8P),               intent(inout) :: q_gpu(1:,         &
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1-self%ngc:,1:) !< Conservative variables.
   integer(I4P)                           :: crown                 !< Crown counter.
   if (associated(self%field_fnl%maps%local_map_bc_crown_gpu)) then
      do crown=1, self%ngc
         call set_boundary_conditions_kernel(ni                     = self%ni                                   ,&
                                             nj                     = self%nj                                   ,&
                                             nk                     = self%nk                                   ,&
                                             ngc                    = self%ngc                                  ,&
                                             nv                     = self%nv                                   ,&
                                             nv_c                   = self%physics%nv_c                         ,&
                                             nv_cl                  = self%physics%nv_cl                        ,&
                                             crown                  = crown                                     ,&
                                             local_map_bc_crown_gpu = self%field_fnl%maps%local_map_bc_crown_gpu,&
                                             q_gpu                  = q_gpu)
      enddo
      call enforce_silver_muller_normal_bc_fnl(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, nv=self%nv, &
                                               hs=self%fdv_half_stencils(1), has_rho=merge(1_I4P, 0_I4P,      &
                                               self%physics%physical_model == PIC_PHYSICAL_MODEL),             &
                                               var_rho=merge(self%nv, 0_I4P,                                  &
                                               self%physics%physical_model == PIC_PHYSICAL_MODEL),             &
                                               local_map_bc_crown_gpu=self%field_fnl%maps%local_map_bc_crown_gpu, &
                                               dxyz_gpu=self%field_fnl%dxyz_gpu, q_gpu=q_gpu)
   endif
   contains
      subroutine set_boundary_conditions_kernel(ni, nj, nk, ngc, nv, nv_c, nv_cl, crown, local_map_bc_crown_gpu, q_gpu)
      !< Set boundary conditions of equation, kernel device.
      integer(I4P), intent(in)    :: ni,nj,nk,ngc                      !< Grid dimensions.
      integer(I4P), intent(in)    :: nv                                !< Number of conservative variables.
      integer(I4P), intent(in)    :: nv_c                              !< Number of physical conservative variables.
      integer(I4P), intent(in)    :: nv_cl                             !< Number of cleaning variables.
      integer(I4P), intent(in)    :: crown                             !< Crown counter.
      integer(I8P), intent(in)    :: local_map_bc_crown_gpu(:,:,:)     !< Local map for face BC ghost cells, "crown" order.
      real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
      integer(I4P)                :: b, c, i, j, k, v                  !< Counter.
      integer(I4P)                :: idelta                            !< IJK i delta step for extrapolation.
      integer(I4P)                :: jdelta                            !< IJK j delta step for extrapolation.
      integer(I4P)                :: kdelta                            !< IJK k delta step for extrapolation.
      integer(I4P)                :: bc_type                           !< Boundary condition type.
      integer(I4P)                :: fec                               !< Boundary fec (1 to 26).
      integer(I4P)                :: fec_1_6                           !< Boundary fec (1 to 6).
      integer(I4P)                :: iref, jref, kref                  !< Interior reference indexes for face BCs.
      integer(I4P)                :: alfa_D, beta_D, gamma_D           !< Tangential/normal D components.
      integer(I4P)                :: alfa_B, beta_B, gamma_B           !< Tangential/normal B components.
      real(R8P)                   :: s1                                !< Face orientation sign.
      !$acc parallel loop independent gang vector &
      !$acc& DEVICEVAR(local_map_bc_crown_gpu, q_gpu) firstprivate(ni,nj,nk,nv,nv_c,nv_cl,crown)
      !$omp OMPLOOP &
      !$omp& DEVICEPTR(local_map_bc_crown_gpu, q_gpu) firstprivate(ni,nj,nk,nv,nv_c,nv_cl,crown)
      do c=1, size(local_map_bc_crown_gpu, dim=1)
         b = local_map_bc_crown_gpu(c, 1 ,crown)
         if (b>0) then
            i       = local_map_bc_crown_gpu(c, 2 ,crown)
            j       = local_map_bc_crown_gpu(c, 3 ,crown)
            k       = local_map_bc_crown_gpu(c, 4 ,crown)
            idelta  = local_map_bc_crown_gpu(c, 5 ,crown)
            jdelta  = local_map_bc_crown_gpu(c, 6 ,crown)
            kdelta  = local_map_bc_crown_gpu(c, 7 ,crown)
            bc_type = local_map_bc_crown_gpu(c, 8 ,crown)
            fec     = local_map_bc_crown_gpu(c, 9 ,crown)
            if (fec > 6_I4P) cycle
            fec_1_6 = fec_1_6_array(fec)
            if (bc_type == BC_EXTRAPOLATION) then
               do v=1, nv
                  q_gpu(b,i,j,k,v) = q_gpu(b,i-idelta,j-jdelta,k-kdelta,v)
               enddo
            elseif (bc_type == BC_NEUMANN) then
               call compute_face_mirror_indexes(face=fec_1_6, ni=ni, nj=nj, nk=nk, i_gc=i, j_gc=j, k_gc=k, &
                                                idelta=idelta, jdelta=jdelta, kdelta=kdelta,               &
                                                i_d=iref, j_d=jref, k_d=kref)
               do v=1, nv
                  q_gpu(b,i,j,k,v) = q_gpu(b,iref,jref,kref,v)
               enddo
            elseif (bc_type == BC_SILVER_MULLER) then
               ! With outward normal n, the Silver-Muller conditions are
               ! B_t = (n x E_d) / c, B_n = B_n,d, E_t = c (B_d x n), E_n = E_n,d.
               call compute_face_mirror_indexes(face=fec_1_6, ni=ni, nj=nj, nk=nk, i_gc=i, j_gc=j, k_gc=k, &
                                                idelta=idelta, jdelta=jdelta, kdelta=kdelta,               &
                                                i_d=iref, j_d=jref, k_d=kref)
               select case(fec_1_6)
               case(1)
                  s1 = -1.0_R8P
                  alfa_D = 2_I4P
                  beta_D = 3_I4P
                  gamma_D = 1_I4P
                  alfa_B = 5_I4P
                  beta_B = 6_I4P
                  gamma_B = 4_I4P
               case(2)
                  s1 = 1.0_R8P
                  alfa_D = 2_I4P
                  beta_D = 3_I4P
                  gamma_D = 1_I4P
                  alfa_B = 5_I4P
                  beta_B = 6_I4P
                  gamma_B = 4_I4P
               case(3)
                  s1 = -1.0_R8P
                  alfa_D = 3_I4P
                  beta_D = 1_I4P
                  gamma_D = 2_I4P
                  alfa_B = 6_I4P
                  beta_B = 4_I4P
                  gamma_B = 5_I4P
               case(4)
                  s1 = 1.0_R8P
                  alfa_D = 3_I4P
                  beta_D = 1_I4P
                  gamma_D = 2_I4P
                  alfa_B = 6_I4P
                  beta_B = 4_I4P
                  gamma_B = 5_I4P
               case(5)
                  s1 = -1.0_R8P
                  alfa_D = 1_I4P
                  beta_D = 2_I4P
                  gamma_D = 3_I4P
                  alfa_B = 4_I4P
                  beta_B = 5_I4P
                  gamma_B = 6_I4P
               case(6)
                  s1 = 1.0_R8P
                  alfa_D = 1_I4P
                  beta_D = 2_I4P
                  gamma_D = 3_I4P
                  alfa_B = 4_I4P
                  beta_B = 5_I4P
                  gamma_B = 6_I4P
               endselect
               q_gpu(b,i,j,k,alfa_D ) =  s1*C0*q_gpu(b,iref,jref,kref,beta_B )*EPS0
               q_gpu(b,i,j,k,beta_D ) = -s1*C0*q_gpu(b,iref,jref,kref,alfa_B)*EPS0
               q_gpu(b,i,j,k,gamma_D) =       q_gpu(b,iref,jref,kref,gamma_D)
               q_gpu(b,i,j,k,alfa_B ) = -s1/C0*q_gpu(b,iref,jref,kref,beta_D)/EPS0
               q_gpu(b,i,j,k,beta_B ) =  s1/C0*q_gpu(b,iref,jref,kref,alfa_D)/EPS0
               q_gpu(b,i,j,k,gamma_B) =       q_gpu(b,iref,jref,kref,gamma_B)
               do v=nv_c-nv_cl+1, nv
                  q_gpu(b,i,j,k,v) = q_gpu(b,iref,jref,kref,v)
               enddo
            elseif (bc_type == BC_PEC) then
               call compute_face_mirror_indexes(face=fec_1_6, ni=ni, nj=nj, nk=nk, i_gc=i, j_gc=j, k_gc=k, &
                                                idelta=idelta, jdelta=jdelta, kdelta=kdelta,               &
                                                i_d=iref, j_d=jref, k_d=kref)
               do v=1, nv
                  q_gpu(b,i,j,k,v) = q_gpu(b,iref,jref,kref,v)
               enddo
               select case(fec_1_6)
               case(1, 2)
                  q_gpu(b,i,j,k,VAR_DX) =  q_gpu(b,iref,jref,kref,VAR_DX)
                  q_gpu(b,i,j,k,VAR_DY) = -q_gpu(b,iref,jref,kref,VAR_DY)
                  q_gpu(b,i,j,k,VAR_DZ) = -q_gpu(b,iref,jref,kref,VAR_DZ)
                  q_gpu(b,i,j,k,VAR_BX) = -q_gpu(b,iref,jref,kref,VAR_BX)
                  q_gpu(b,i,j,k,VAR_BY) =  q_gpu(b,iref,jref,kref,VAR_BY)
                  q_gpu(b,i,j,k,VAR_BZ) =  q_gpu(b,iref,jref,kref,VAR_BZ)
               case(3, 4)
                  q_gpu(b,i,j,k,VAR_DX) = -q_gpu(b,iref,jref,kref,VAR_DX)
                  q_gpu(b,i,j,k,VAR_DY) =  q_gpu(b,iref,jref,kref,VAR_DY)
                  q_gpu(b,i,j,k,VAR_DZ) = -q_gpu(b,iref,jref,kref,VAR_DZ)
                  q_gpu(b,i,j,k,VAR_BX) =  q_gpu(b,iref,jref,kref,VAR_BX)
                  q_gpu(b,i,j,k,VAR_BY) = -q_gpu(b,iref,jref,kref,VAR_BY)
                  q_gpu(b,i,j,k,VAR_BZ) =  q_gpu(b,iref,jref,kref,VAR_BZ)
               case(5, 6)
                  q_gpu(b,i,j,k,VAR_DX) = -q_gpu(b,iref,jref,kref,VAR_DX)
                  q_gpu(b,i,j,k,VAR_DY) = -q_gpu(b,iref,jref,kref,VAR_DY)
                  q_gpu(b,i,j,k,VAR_DZ) =  q_gpu(b,iref,jref,kref,VAR_DZ)
                  q_gpu(b,i,j,k,VAR_BX) =  q_gpu(b,iref,jref,kref,VAR_BX)
                  q_gpu(b,i,j,k,VAR_BY) =  q_gpu(b,iref,jref,kref,VAR_BY)
                  q_gpu(b,i,j,k,VAR_BZ) = -q_gpu(b,iref,jref,kref,VAR_BZ)
               endselect
            elseif (bc_type == BC_DIRICHLET) then
               do v=1, nv
                  q_gpu(b,i,j,k,v) = 0._R8P
               enddo
            elseif (bc_type == BC_PERIOD) then
               do v=1, nv
                  q_gpu(b,i,j,k,v) = 0._R8P
               enddo
               select case(fec_1_6)
               case(1)
                  do v=1, nv_c
                     q_gpu(b,i,j,k,v) = q_gpu(b,ni+i,j   ,k   ,v)
                  enddo
               case(2)
                  do v=1, nv_c
                     q_gpu(b,i,j,k,v) = q_gpu(b,i-ni,j   ,k   ,v)
                  enddo
               case(3)
                  do v=1, nv_c
                     q_gpu(b,i,j,k,v) = q_gpu(b,i   ,nj+j,k   ,v)
                  enddo
               case(4)
                  do v=1, nv_c
                     q_gpu(b,i,j,k,v) = q_gpu(b,i   ,j-nj,k   ,v)
                  enddo
               case(5)
                  do v=1, nv_c
                     q_gpu(b,i,j,k,v) = q_gpu(b,i   ,j   ,nk+k,v)
                  enddo
               case(6)
                  do v=1, nv_c
                     q_gpu(b,i,j,k,v) = q_gpu(b,i   ,j   ,k-nk,v)
                  enddo
               endselect
            endif
         endif
      enddo
      endsubroutine set_boundary_conditions_kernel

      subroutine compute_face_mirror_indexes(face, ni, nj, nk, i_gc, j_gc, k_gc, idelta, jdelta, kdelta, i_d, j_d, k_d)
      !$acc routine seq
      !< Return the donor indexes mirrored across the selected boundary face.
      integer(I4P), intent(in)  :: face                       !< Boundary face index in [1, 6].
      integer(I4P), intent(in)  :: ni, nj, nk                !< Interior grid extents.
      integer(I4P), intent(in)  :: i_gc, j_gc, k_gc         !< Ghost-cell indexes.
      integer(I4P), intent(in)  :: idelta, jdelta, kdelta   !< One-step inward deltas.
      integer(I4P), intent(out) :: i_d, j_d, k_d            !< Mirrored donor indexes.

      i_d = i_gc - idelta
      j_d = j_gc - jdelta
      k_d = k_gc - kdelta

      select case(face)
      case(1)
         i_d = 1_I4P - i_gc
      case(2)
         i_d = 2_I4P * ni + 1_I4P - i_gc
      case(3)
         j_d = 1_I4P - j_gc
      case(4)
         j_d = 2_I4P * nj + 1_I4P - j_gc
      case(5)
         k_d = 1_I4P - k_gc
      case(6)
         k_d = 2_I4P * nk + 1_I4P - k_gc
      endselect
      endsubroutine compute_face_mirror_indexes

      subroutine enforce_silver_muller_normal_bc_fnl(ni, nj, nk, ngc, nv, hs, has_rho, var_rho, local_map_bc_crown_gpu, dxyz_gpu, q_gpu)
      integer(I4P), intent(in)    :: ni, nj, nk, ngc, nv, hs, has_rho, var_rho
      integer(I8P), intent(in)    :: local_map_bc_crown_gpu(:,:,:)
      real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)
      real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
      integer(I4P), parameter     :: sm_face_sweeps = 2_I4P
      integer(I4P)                :: c, face_stage, sweep

      if (hs <= 0_I4P) return
      if (ngc < hs) call mpih_fnl%error_stop(msg='Silver_Muller requires ngc >= fdv_half_stencils(1)')

      do sweep=1, sm_face_sweeps
         do face_stage=1, 6
            !$acc parallel loop independent gang vector &
            !$acc& DEVICEVAR(local_map_bc_crown_gpu, dxyz_gpu, q_gpu) firstprivate(ni, nj, nk, ngc, nv, hs, has_rho, var_rho, face_stage)
            !$omp OMPLOOP DEVICEPTR(local_map_bc_crown_gpu, dxyz_gpu, q_gpu) &
            !$omp& firstprivate(ni, nj, nk, ngc, nv, hs, has_rho, var_rho, face_stage)
            do c=1, size(local_map_bc_crown_gpu, dim=1)
               call enforce_silver_muller_normal_line_kernel(c=c, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, hs=hs, &
                                                             has_rho=has_rho, var_rho=var_rho,                    &
                                                             local_map_bc_crown_gpu=local_map_bc_crown_gpu,       &
                                                             dxyz_gpu=dxyz_gpu, face_filter=face_stage, q_gpu=q_gpu)
            enddo
         enddo
         do face_stage=6, 1, -1
            !$acc parallel loop independent gang vector &
            !$acc& DEVICEVAR(local_map_bc_crown_gpu, dxyz_gpu, q_gpu) firstprivate(ni, nj, nk, ngc, nv, hs, has_rho, var_rho, face_stage)
            !$omp OMPLOOP DEVICEPTR(local_map_bc_crown_gpu, dxyz_gpu, q_gpu) &
            !$omp& firstprivate(ni, nj, nk, ngc, nv, hs, has_rho, var_rho, face_stage)
            do c=1, size(local_map_bc_crown_gpu, dim=1)
               call enforce_silver_muller_normal_line_kernel(c=c, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, hs=hs, &
                                                             has_rho=has_rho, var_rho=var_rho,                    &
                                                             local_map_bc_crown_gpu=local_map_bc_crown_gpu,       &
                                                             dxyz_gpu=dxyz_gpu, face_filter=face_stage, q_gpu=q_gpu)
            enddo
         enddo
      enddo
      endsubroutine enforce_silver_muller_normal_bc_fnl

      subroutine enforce_silver_muller_normal_line_kernel(c, ni, nj, nk, ngc, nv, hs, has_rho, var_rho, local_map_bc_crown_gpu, dxyz_gpu, face_filter, q_gpu)
      !$acc routine seq
      integer(I4P), intent(in)    :: c, ni, nj, nk, ngc, nv, hs, has_rho, var_rho, face_filter
      integer(I8P), intent(in)    :: local_map_bc_crown_gpu(:,:,:)
      real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)
      real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
      integer(I4P)                :: b, i, j, k, face
      real(R8P)                   :: dxyz(3)

      b = local_map_bc_crown_gpu(c, 1, 1)
      if (b <= 0_I4P) return
      if (local_map_bc_crown_gpu(c, 8, 1) /= BC_SILVER_MULLER) return
      if (local_map_bc_crown_gpu(c, 9, 1) > 6_I4P) return

      i = local_map_bc_crown_gpu(c, 2, 1)
      j = local_map_bc_crown_gpu(c, 3, 1)
      k = local_map_bc_crown_gpu(c, 4, 1)
      face = local_map_bc_crown_gpu(c, 9, 1)
      if (face /= face_filter) return
      if (.not. is_face_line_seed_fnl(face=face, i=i, j=j, k=k, ni=ni, nj=nj, nk=nk)) return

      dxyz = dxyz_gpu(b, 1:3)
      call solve_silver_muller_normal_line_fnl(q_gpu=q_gpu, ngc=ngc, ni=ni, nj=nj, nk=nk, b=b, face=face, &
                                               i_seed=i, j_seed=j, k_seed=k, hs=hs, dxyz=dxyz,             &
                                               has_rho=has_rho, var_rho=var_rho)
      endsubroutine enforce_silver_muller_normal_line_kernel

      logical function is_face_line_seed_fnl(face, i, j, k, ni, nj, nk)
      !$acc routine seq
      integer(I4P), intent(in) :: face, i, j, k, ni, nj, nk

      is_face_line_seed_fnl = .false.
      select case(face)
      case(1)
         is_face_line_seed_fnl = i == 0_I4P      .and. j >= 1_I4P .and. j <= nj .and. k >= 1_I4P .and. k <= nk
      case(2)
         is_face_line_seed_fnl = i == ni + 1_I4P .and. j >= 1_I4P .and. j <= nj .and. k >= 1_I4P .and. k <= nk
      case(3)
         is_face_line_seed_fnl = j == 0_I4P      .and. i >= 1_I4P .and. i <= ni .and. k >= 1_I4P .and. k <= nk
      case(4)
         is_face_line_seed_fnl = j == nj + 1_I4P .and. i >= 1_I4P .and. i <= ni .and. k >= 1_I4P .and. k <= nk
      case(5)
         is_face_line_seed_fnl = k == 0_I4P      .and. i >= 1_I4P .and. i <= ni .and. j >= 1_I4P .and. j <= nj
      case(6)
         is_face_line_seed_fnl = k == nk + 1_I4P .and. i >= 1_I4P .and. i <= ni .and. j >= 1_I4P .and. j <= nj
      endselect
      endfunction is_face_line_seed_fnl

      subroutine solve_silver_muller_normal_line_fnl(q_gpu, ngc, ni, nj, nk, b, face, i_seed, j_seed, k_seed, hs, dxyz, has_rho, var_rho)
      !$acc routine seq
      integer(I4P), intent(in)    :: ngc, ni, nj, nk, b, face, i_seed, j_seed, k_seed, hs, has_rho, var_rho
      real(R8P),    intent(in)    :: dxyz(3)
      real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
      integer(I4P)                :: dir_t1, dir_t2
      integer(I4P)                :: var_d_t1, var_d_t2, var_d_n
      integer(I4P)                :: var_b_t1, var_b_t2, var_b_n

      call silver_muller_face_metadata_fnl(face=face, dir_t1=dir_t1, dir_t2=dir_t2,                              &
                                           var_d_t1=var_d_t1, var_d_t2=var_d_t2, var_d_n=var_d_n,               &
                                           var_b_t1=var_b_t1, var_b_t2=var_b_t2, var_b_n=var_b_n)

      call solve_silver_muller_normal_field_fnl(q_gpu=q_gpu, ngc=ngc, ni=ni, nj=nj, nk=nk, b=b, face=face,         &
                                                i_seed=i_seed, j_seed=j_seed, k_seed=k_seed, hs=hs, dxyz=dxyz,         &
                                                dir_t1=dir_t1, dir_t2=dir_t2, var_t1=var_d_t1, var_t2=var_d_t2,      &
                                                var_n=var_d_n, has_target=has_rho, target_var=var_rho)

      call solve_silver_muller_normal_field_fnl(q_gpu=q_gpu, ngc=ngc, ni=ni, nj=nj, nk=nk, b=b, face=face,         &
                                                i_seed=i_seed, j_seed=j_seed, k_seed=k_seed, hs=hs, dxyz=dxyz,         &
                                                dir_t1=dir_t1, dir_t2=dir_t2, var_t1=var_b_t1, var_t2=var_b_t2,      &
                                                var_n=var_b_n, has_target=0_I4P, target_var=0_I4P)
      endsubroutine solve_silver_muller_normal_line_fnl

      subroutine silver_muller_face_metadata_fnl(face, dir_t1, dir_t2, var_d_t1, var_d_t2, var_d_n, var_b_t1, var_b_t2, var_b_n)
      !$acc routine seq
      integer(I4P), intent(in)  :: face
      integer(I4P), intent(out) :: dir_t1, dir_t2
      integer(I4P), intent(out) :: var_d_t1, var_d_t2, var_d_n
      integer(I4P), intent(out) :: var_b_t1, var_b_t2, var_b_n

      select case(face)
      case(1, 2)
         dir_t1 = 2_I4P ; dir_t2 = 3_I4P
         var_d_t1 = VAR_DY ; var_d_t2 = VAR_DZ ; var_d_n = VAR_DX
         var_b_t1 = VAR_BY ; var_b_t2 = VAR_BZ ; var_b_n = VAR_BX
      case(3, 4)
         dir_t1 = 3_I4P ; dir_t2 = 1_I4P
         var_d_t1 = VAR_DZ ; var_d_t2 = VAR_DX ; var_d_n = VAR_DY
         var_b_t1 = VAR_BZ ; var_b_t2 = VAR_BX ; var_b_n = VAR_BY
      case(5, 6)
         dir_t1 = 1_I4P ; dir_t2 = 2_I4P
         var_d_t1 = VAR_DX ; var_d_t2 = VAR_DY ; var_d_n = VAR_DZ
         var_b_t1 = VAR_BX ; var_b_t2 = VAR_BY ; var_b_n = VAR_BZ
      case default
         dir_t1 = 0_I4P ; dir_t2 = 0_I4P
         var_d_t1 = 0_I4P ; var_d_t2 = 0_I4P ; var_d_n = 0_I4P
         var_b_t1 = 0_I4P ; var_b_t2 = 0_I4P ; var_b_n = 0_I4P
      endselect
      endsubroutine silver_muller_face_metadata_fnl

      subroutine solve_silver_muller_normal_field_fnl(q_gpu, ngc, ni, nj, nk, b, face, i_seed, j_seed, k_seed, hs, dxyz, &
                                                      dir_t1, dir_t2, var_t1, var_t2, var_n, has_target, target_var)
      !$acc routine seq
      integer(I4P), intent(in)    :: ngc, ni, nj, nk, b, face, i_seed, j_seed, k_seed, hs
      integer(I4P), intent(in)    :: dir_t1, dir_t2, var_t1, var_t2, var_n, has_target, target_var
      real(R8P),    intent(in)    :: dxyz(3)
      real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
      real(R8P)                   :: unknown(FDV_S_MAX)
      real(R8P)                   :: rhs, coeff, target
      integer(I4P)                :: eq, u, u_new, m, dir_n
      integer(I4P)                :: ic, jc, kc

      unknown = 0._R8P
      dir_n = face_normal_direction_fnl(face)

      do eq=hs, 1, -1
         call interior_cell_from_face_fnl(face=face, eq=eq, ni=ni, nj=nj, nk=nk, i_seed=i_seed, j_seed=j_seed, k_seed=k_seed, &
                                          ic=ic, jc=jc, kc=kc)

         if (has_target /= 0_I4P) then
            target = q_gpu(b, ic, jc, kc, target_var)
         else
            target = 0._R8P
         endif

         rhs = target
         rhs = rhs - tangential_divergence_at_cell_fnl(q_gpu=q_gpu, ngc=ngc, b=b, i=ic, j=jc, k=kc, hs=hs, dxyz=dxyz,           &
                                                       dir_t1=dir_t1, dir_t2=dir_t2, var_t1=var_t1, var_t2=var_t2)
         rhs = rhs - normal_known_contribution_fnl(q_gpu=q_gpu, ngc=ngc, ni=ni, nj=nj, nk=nk, b=b, face=face, i=ic, j=jc, k=kc, &
                                                   hs=hs, dxyz=dxyz, var_n=var_n)

         u_new = hs - eq + 1_I4P
         do u=1, u_new-1
            m = eq + u - 1_I4P
            coeff = silver_muller_unknown_coefficient_fnl(face=face, m=m, hs=hs, ds=dxyz(dir_n))
            rhs = rhs - coeff * unknown(u)
         enddo

         coeff = silver_muller_unknown_coefficient_fnl(face=face, m=hs, hs=hs, ds=dxyz(dir_n))
         unknown(u_new) = rhs / coeff
      enddo

      do u=1, hs
         call ghost_cell_from_face_fnl(face=face, ghost_layer=u, ni=ni, nj=nj, nk=nk, i_seed=i_seed, j_seed=j_seed, k_seed=k_seed, &
                                       ic=ic, jc=jc, kc=kc)
         q_gpu(b, ic, jc, kc, var_n) = unknown(u)
      enddo
      endsubroutine solve_silver_muller_normal_field_fnl

      integer(I4P) function face_normal_direction_fnl(face)
      !$acc routine seq
      integer(I4P), intent(in) :: face
      select case(face)
      case(1, 2)
         face_normal_direction_fnl = 1_I4P
      case(3, 4)
         face_normal_direction_fnl = 2_I4P
      case(5, 6)
         face_normal_direction_fnl = 3_I4P
      case default
         face_normal_direction_fnl = 0_I4P
      endselect
      endfunction face_normal_direction_fnl

      real(R8P) function silver_muller_unknown_coefficient_fnl(face, m, hs, ds)
      !$acc routine seq
      integer(I4P), intent(in) :: face, m, hs
      real(R8P),    intent(in) :: ds
      silver_muller_unknown_coefficient_fnl = FD1_CC(m, hs) / ds
      if (face == 1_I4P .or. face == 3_I4P .or. face == 5_I4P) then
         silver_muller_unknown_coefficient_fnl = -silver_muller_unknown_coefficient_fnl
      endif
      endfunction silver_muller_unknown_coefficient_fnl

      subroutine interior_cell_from_face_fnl(face, eq, ni, nj, nk, i_seed, j_seed, k_seed, ic, jc, kc)
      !$acc routine seq
      integer(I4P), intent(in)  :: face, eq, ni, nj, nk, i_seed, j_seed, k_seed
      integer(I4P), intent(out) :: ic, jc, kc
      select case(face)
      case(1)
         ic = eq ; jc = j_seed ; kc = k_seed
      case(2)
         ic = ni - eq + 1_I4P ; jc = j_seed ; kc = k_seed
      case(3)
         ic = i_seed ; jc = eq ; kc = k_seed
      case(4)
         ic = i_seed ; jc = nj - eq + 1_I4P ; kc = k_seed
      case(5)
         ic = i_seed ; jc = j_seed ; kc = eq
      case(6)
         ic = i_seed ; jc = j_seed ; kc = nk - eq + 1_I4P
      case default
         ic = 0_I4P ; jc = 0_I4P ; kc = 0_I4P
      endselect
      endsubroutine interior_cell_from_face_fnl

      subroutine ghost_cell_from_face_fnl(face, ghost_layer, ni, nj, nk, i_seed, j_seed, k_seed, ic, jc, kc)
      !$acc routine seq
      integer(I4P), intent(in)  :: face, ghost_layer, ni, nj, nk, i_seed, j_seed, k_seed
      integer(I4P), intent(out) :: ic, jc, kc
      select case(face)
      case(1)
         ic = 1_I4P - ghost_layer ; jc = j_seed ; kc = k_seed
      case(2)
         ic = ni + ghost_layer    ; jc = j_seed ; kc = k_seed
      case(3)
         ic = i_seed ; jc = 1_I4P - ghost_layer ; kc = k_seed
      case(4)
         ic = i_seed ; jc = nj + ghost_layer    ; kc = k_seed
      case(5)
         ic = i_seed ; jc = j_seed ; kc = 1_I4P - ghost_layer
      case(6)
         ic = i_seed ; jc = j_seed ; kc = nk + ghost_layer
      case default
         ic = 0_I4P ; jc = 0_I4P ; kc = 0_I4P
      endselect
      endsubroutine ghost_cell_from_face_fnl

      real(R8P) function tangential_divergence_at_cell_fnl(q_gpu, ngc, b, i, j, k, hs, dxyz, dir_t1, dir_t2, var_t1, var_t2)
      !$acc routine seq
      integer(I4P), intent(in) :: ngc, b, i, j, k, hs, dir_t1, dir_t2, var_t1, var_t2
      real(R8P),    intent(in) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
      real(R8P),    intent(in) :: dxyz(3)
      integer(I4P)             :: m
      tangential_divergence_at_cell_fnl = 0._R8P
      do m=1, hs
         tangential_divergence_at_cell_fnl = tangential_divergence_at_cell_fnl + FD1_CC(m, hs) * &
            (component_along_direction_fnl(q_gpu, ngc, b, var_t1, dir_t1, i, j, k,  m) - &
             component_along_direction_fnl(q_gpu, ngc, b, var_t1, dir_t1, i, j, k, -m)) / dxyz(dir_t1)
         tangential_divergence_at_cell_fnl = tangential_divergence_at_cell_fnl + FD1_CC(m, hs) * &
            (component_along_direction_fnl(q_gpu, ngc, b, var_t2, dir_t2, i, j, k,  m) - &
             component_along_direction_fnl(q_gpu, ngc, b, var_t2, dir_t2, i, j, k, -m)) / dxyz(dir_t2)
      enddo
      endfunction tangential_divergence_at_cell_fnl

      real(R8P) function normal_known_contribution_fnl(q_gpu, ngc, ni, nj, nk, b, face, i, j, k, hs, dxyz, var_n)
      !$acc routine seq
      integer(I4P), intent(in) :: ngc, ni, nj, nk, b, face, i, j, k, hs, var_n
      real(R8P),    intent(in) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
      real(R8P),    intent(in) :: dxyz(3)
      integer(I4P)             :: dir_n, m, depth

      dir_n = face_normal_direction_fnl(face)
      depth = normal_depth_from_cell_fnl(face=face, i=i, j=j, k=k, ni=ni, nj=nj, nk=nk)
      normal_known_contribution_fnl = 0._R8P

      select case(face)
      case(1, 3, 5)
         do m=1, hs
            normal_known_contribution_fnl = normal_known_contribution_fnl + FD1_CC(m, hs) * &
               component_along_direction_fnl(q_gpu, ngc, b, var_n, dir_n, i, j, k,  m) / dxyz(dir_n)
         enddo
         do m=1, depth-1
            normal_known_contribution_fnl = normal_known_contribution_fnl - FD1_CC(m, hs) * &
               component_along_direction_fnl(q_gpu, ngc, b, var_n, dir_n, i, j, k, -m) / dxyz(dir_n)
         enddo
      case(2, 4, 6)
         do m=1, depth-1
            normal_known_contribution_fnl = normal_known_contribution_fnl + FD1_CC(m, hs) * &
               component_along_direction_fnl(q_gpu, ngc, b, var_n, dir_n, i, j, k,  m) / dxyz(dir_n)
         enddo
         do m=1, hs
            normal_known_contribution_fnl = normal_known_contribution_fnl - FD1_CC(m, hs) * &
               component_along_direction_fnl(q_gpu, ngc, b, var_n, dir_n, i, j, k, -m) / dxyz(dir_n)
         enddo
      endselect
      endfunction normal_known_contribution_fnl

      integer(I4P) function normal_depth_from_cell_fnl(face, i, j, k, ni, nj, nk)
      !$acc routine seq
      integer(I4P), intent(in) :: face, i, j, k, ni, nj, nk
      select case(face)
      case(1)
         normal_depth_from_cell_fnl = i
      case(2)
         normal_depth_from_cell_fnl = ni - i + 1_I4P
      case(3)
         normal_depth_from_cell_fnl = j
      case(4)
         normal_depth_from_cell_fnl = nj - j + 1_I4P
      case(5)
         normal_depth_from_cell_fnl = k
      case(6)
         normal_depth_from_cell_fnl = nk - k + 1_I4P
      case default
         normal_depth_from_cell_fnl = 0_I4P
      endselect
      endfunction normal_depth_from_cell_fnl

      real(R8P) function component_along_direction_fnl(q_gpu, ngc, b, var, dir, i, j, k, offset)
      !$acc routine seq
      integer(I4P), intent(in) :: ngc, b, var, dir, i, j, k, offset
      real(R8P),    intent(in) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
      select case(dir)
      case(1)
         component_along_direction_fnl = q_gpu(b, i + offset, j, k, var)
      case(2)
         component_along_direction_fnl = q_gpu(b, i, j + offset, k, var)
      case(3)
         component_along_direction_fnl = q_gpu(b, i, j, k + offset, var)
      case default
         component_along_direction_fnl = 0._R8P
      endselect
      endfunction component_along_direction_fnl
   endsubroutine set_boundary_conditions

   subroutine set_initial_conditions(self, is_restart)
   !< Set initial conditions of field.
   class(prism_fnl_object), intent(inout) :: self       !< The equation.
   logical,                 intent(in)    :: is_restart !< Branching sentinel for restart/non restart path.

   if (.not.is_restart) call self%ic%set_initial_conditions(physics=self%physics, field=self%adam%field, grid=self%adam%grid, &
                                                            q=self%q)
   if (.not.is_restart) call self%initialize_pic_time_zero()

   call self%initialize_coils
   if (.not.is_restart) then
      call self%compute_coils_current_time_zero()
      if (self%physics%physical_model == PIC_PHYSICAL_MODEL) &
         call self%impose_pic_fields_time_zero(ivar=VAR_DX)
      if (maxval(abs(self%q(self%physics%var_Jx:self%physics%var_Jz,:,:,:,:))) > 0._R8P) &
         call self%impose_pic_fields_time_zero(ivar=VAR_BX)
   endif

   call self%copy_cpu_gpu

   if (.not.is_restart) then
      call self%apply_fWL_correction(q_gpu=self%q_gpu)
      if (self%external_fields%ef_type/=EF_TYPE_NONE) &
         call add_external_fields_dev(external_fields=self%external_fields, field_gpu=self%field_fnl, &
                                      dt=0._R8P, time=0._R8P, q_gpu=self%q_gpu)
      if (self%physics%physical_model == PIC_PHYSICAL_MODEL) then
         call dev_memcpy_from_device(bb=self%db5,ij=[1,5],tb=self%hb5,dst=self%q,src=self%q_gpu,buf=self%buf_5D_R8P)
         call self%weight_pic_fields_time_zero()
         call self%pic_fnl%copy_cpu_gpu(pic=self%pic, q_pic=self%q_pic, pic_fields=self%pic_fields)
      endif
   endif
   endsubroutine set_initial_conditions

   subroutine update_ghost(self, q_gpu, step, s)
   !< Update ghost cells (intra-realm only).
   !<
   !< Inter-realm seam ghost cells are filled by the forest's Phase 2 seam
   !< exchange (between assemble and evaluate) via
   !< `fill_seam_from_peer_forest`, NOT here. The legacy
   !< `realm(:)` optional dummy and the
   !< `exchange_inter_realm_halos_forest` call have been retired by the
   !< agnostic-dummy seam redesign.
   class(prism_fnl_object), intent(inout)                   :: self            !< The equation.
   real(R8P),               intent(inout)                   :: q_gpu(1:,         &
                                                                     1-self%ngc:,&
                                                                     1-self%ngc:,&
                                                                     1-self%ngc:,&
                                                                     1:)       !< Conservative variables.
   integer(I4P),            intent(in),    optional         :: step            !< Step to be performed in async comp.
   integer(I4P),            intent(in),    optional         :: s               !< Stage counter.
   logical                                                  :: do_local_update !< Flag for triggering local update.
   logical                                                  :: do_set_bc       !< Flag for triggering setting bc.

   ! perform local update if step is not specified or if first step is selected
   do_local_update = .false.
   do_set_bc       = .false.
   if (.not.present(step)) then
      do_local_update = .true.
      do_set_bc       = .true.
   else
      if (step==1) do_local_update = .true.
      if (step==3) do_set_bc       = .true.
   endif

   if (do_local_update) call self%field_fnl%update_ghost_local_gpu(q_gpu=q_gpu)
                        call self%field_fnl%update_ghost_mpi_gpu(comm_map_send_ptr_ghost=self%adam%maps%comm_map_send_ptr_ghost, &
                                                                 comm_map_recv_ptr_ghost=self%adam%maps%comm_map_recv_ptr_ghost, &
                                                                 q_gpu=q_gpu, step=step)
   if (do_set_bc) call self%set_boundary_conditions(q_gpu=q_gpu)
   ! Trailing coil re-stamp, CPU-parity (issue #26 G2): the CPU update_ghost ends with an
   ! unconditional compute_coils_current, so q's J rows are analytic-fresh over the FULL
   ! extent (interiors + ghosts) after every ghost refresh -- the ghost machinery never
   ! has the last word on J. The call was dropped on FNL in the fWL resequencing (the `s`
   ! dummy survived as "API symmetry"); its absence left the committed q with a
   ! one-dt-lagged J interior and stale J ghosts. The allocated() guard reproduces the
   ! CPU low-storage behavior exactly: CPU LS calls compute_residuals without `s`, so its
   ! stamp is the bare-time one; FNL LS threads `s` through, and gamm does not exist for
   ! LS schemes.
   if (present(s) .and. allocated(self%rk%gamm)) then
      call self%compute_coils_current(q_gpu=q_gpu, gamm=self%rk%gamm(s))
   else
      call self%compute_coils_current(q_gpu=q_gpu)
   endif
   endsubroutine update_ghost

   subroutine update_rk_ghost(self, dt, phi_gpu)
   !< Update RK ghost cells.
   class(prism_fnl_object), intent(inout)        :: self        !< RK object.
   real(R8P),               intent(in)           :: dt          !< Current time step.
   real(R8P),               intent(in), optional :: phi_gpu(1:,          &
                                                            1-self%ngc:, &
                                                            1-self%ngc:, &
                                                            1-self%ngc:, &
                                                            1:) !< IB distance.
   ! to be implemented
   endsubroutine update_rk_ghost

   ! numerical methods, FDV operators
   subroutine compute_curl_fd_dev(self, ivar, q_gpu, curl_gpu)
   !< Compute curl of vector fields, div(q(ivar:ivar+2), using finite difference schemes.
   class(prism_fnl_object), intent(in)    :: self                                                !< The equation.
   integer(I4P),            intent(in)    :: ivar                                                !< Start index of variable of q.
   real(R8P),               intent(in)    :: q_gpu(   1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: curl_gpu(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Curl.

   call  compute_curl_fd_dev_kernel(ni            = self%ni                  ,&
                                    nj            = self%nj                  ,&
                                    nk            = self%nk                  ,&
                                    ngc           = self%ngc                 ,&
                                    blocks_number = self%blocks_number       ,&
                                    ivar          = ivar                     ,&
                                    s1            = self%fdv_half_stencils(1),&
                                    dxyz_gpu      = self%field_fnl%dxyz_gpu  ,&
                                    q_gpu         = q_gpu                    ,&
                                    curl_gpu      = curl_gpu)
   contains
      subroutine compute_curl_fd_dev_kernel(ni,nj,nk,ngc,blocks_number, &
                                            ivar,s1,dxyz_gpu,q_gpu,curl_gpu)
      !< Compute curl, space operator, centered finite difference schemes, kernel device.
      integer(I4P), intent(in)    :: ni,nj,nk,ngc,blocks_number                          !< Grids dimensions.
      integer(I4P), intent(in)    :: ivar                                                !< Start index of variable of q.
      integer(I4P), intent(in)    :: s1                                                  !< Half FDV stencil length.
      real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)                                     !< Delta cells GPU [nb,3].
      real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)                   !< Field cell centered variables.
      real(R8P),    intent(inout) :: curl_gpu(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Curl.
      integer(I4P)                :: i,j,k,b,s                                           !< Counter
      ! rank 1D stencil for computations on device that contiguos memory is mandatory
      real(R8P) :: qsx_y(1-s1:1+s1) !< Y component of vector field over the x stencil.
      real(R8P) :: qsx_z(1-s1:1+s1) !< Z component of vector field over the x stencil.
      real(R8P) :: qsy_x(1-s1:1+s1) !< X component of vector field over the y stencil.
      real(R8P) :: qsy_z(1-s1:1+s1) !< Z component of vector field over the y stencil.
      real(R8P) :: qsz_x(1-s1:1+s1) !< X component of vector field over the z stencil.
      real(R8P) :: qsz_y(1-s1:1+s1) !< Y component of vector field over the z stencil.
      real(R8P) :: dxyz_b(3)        !< Per-block deltas, PRIVATE copy (no strided-section temp: issue #26 G1.b).
      real(R8P) :: curl_(3)         !< Curl, PRIVATE result buffer (no strided-section OUT temp: issue #26 G1.b).

      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu,curl_gpu) &
      !$acc& firstprivate(ni,nj,nk,blocks_number,ivar,s1)                                        &
      !$acc& private(qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y,dxyz_b,curl_)
      !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,q_gpu,curl_gpu) &
      !$omp& firstprivate(ni,nj,nk,blocks_number,ivar,s1) &
      !$omp& private(qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y,dxyz_b,curl_)
      do b=1,blocks_number
      do k=1,nk
      do j=1,nj
      do i=1,ni
         ! hoist per-block deltas into a PRIVATE vector and return the curl through a PRIVATE
         ! buffer (issue #26 G1.b, rule from #22 F1-bis): both the strided IN section
         ! dxyz_gpu(b,1:3) and the strided OUT section curl_gpu(b,i,j,k,ivar:) materialize
         ! compiler temporaries that are NOT privatized -- threads race on them; benign on
         ! uniform grids, live at 2:1 level mixes. This kernel runs on mixed-level AMR
         ! topologies (curl field saves, coil diagnostics).
         dxyz_b(1) = dxyz_gpu(b,1) ; dxyz_b(2) = dxyz_gpu(b,2) ; dxyz_b(3) = dxyz_gpu(b,3)
         !$acc loop seq
         do s=1-s1, 1+s1
            qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,ivar+1)
            qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,ivar+2)
            qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,ivar+0)
            qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,ivar+2)
            qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,ivar+0)
            qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,ivar+1)
         enddo
         call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_b,                   &
                                           qsx_y=qsx_y,qsx_z=qsx_z,qsy_x=qsy_x,&
                                           qsy_z=qsy_z,qsz_x=qsz_x,qsz_y=qsz_y,&
                                           curl=curl_)
         curl_gpu(b,i,j,k,ivar+0) = curl_(1)
         curl_gpu(b,i,j,k,ivar+1) = curl_(2)
         curl_gpu(b,i,j,k,ivar+2) = curl_(3)
      enddo
      enddo
      enddo
      enddo
      endsubroutine compute_curl_fd_dev_kernel
   endsubroutine compute_curl_fd_dev

   subroutine compute_curl_fv_dev(self, ivar, q_gpu, curl_gpu)
   !< Compute curl of vector fields, div(q(ivar:ivar+2), using finite volume schemes.
   class(prism_fnl_object), intent(in)    :: self                                                !< The equation.
   integer(I4P),            intent(in)    :: ivar                                                !< Start index of variable of q.
   real(R8P),               intent(in)    :: q_gpu(   1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: curl_gpu(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Curl.
   integer(I4P)                           :: i,j,k,b                                             !< Counter.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz_gpu=>self%field_fnl%dxyz_gpu,&
             hs=>self%fdv_half_stencils(1))

   endassociate
   endsubroutine compute_curl_fv_dev

   subroutine compute_derivative1_fd_dev(self, dir, ivar, q_gpu, dq_ds_gpu)
   !< Compute derivative1 of scalar fields, dq(ivar)/ds, using finite difference schemes.
   class(prism_fnl_object), intent(in)    :: self                                              !< The equation.
   integer(I4P),            intent(in)    :: dir                                               !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),            intent(in)    :: ivar                                              !< Start index of variable of q.
   real(R8P),               intent(in)    :: q_gpu(1:, 1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: dq_ds_gpu(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Derivative1, dq/ds.
   integer(I4P)                           :: i,j,k,b                                           !< Counter.
   integer(I4P)                           :: is,js,ks                                          !< Stencils.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz_gpu=>self%field_fnl%dxyz_gpu,&
             hs=>self%fdv_half_stencils(1))

   endassociate
   endsubroutine compute_derivative1_fd_dev

   subroutine compute_derivative1_fv_dev(self, dir, ivar, q_gpu, dq_ds_gpu)
   !< Compute derivative1 of scalar fields, dq(ivar)/ds, using finite volume schemes.
   class(prism_fnl_object), intent(in)    :: self                                              !< The equation.
   integer(I4P),            intent(in)    :: dir                                               !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),            intent(in)    :: ivar                                              !< Start index of variable of q.
   real(R8P),               intent(in)    :: q_gpu(1:, 1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: dq_ds_gpu(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Derivative1, dq/ds.
   integer(I4P)                           :: i,j,k,b                                           !< Counter.
   integer(I4P)                           :: is,js,ks                                          !< Stencils.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz_gpu=>self%field_fnl%dxyz_gpu,&
             hs=>self%fdv_half_stencils(1))

   endassociate
   endsubroutine compute_derivative1_fv_dev

   subroutine compute_derivative2_fd_dev(self, dir, ivar, q_gpu, d2q_ds2_gpu)
   !< Compute derivative2 of scalar fields, d2q(ivar)/ds2, using finite difference schemes.
   class(prism_fnl_object), intent(in)    :: self                                                !< The equation.
   integer(I4P),            intent(in)    :: dir                                                 !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),            intent(in)    :: ivar                                                !< Start index of variable of q.
   real(R8P),               intent(in)    :: q_gpu(1:, 1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)   !< Field variables.
   real(R8P),               intent(inout) :: d2q_ds2_gpu(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Derivative2, d2q/ds2.
   integer(I4P)                           :: i,j,k,b                                             !< Counter.
   integer(I4P)                           :: is,js,ks                                            !< Stencils.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz_gpu=>self%field_fnl%dxyz_gpu,&
             hs=>self%fdv_half_stencils(1))

   endassociate
   endsubroutine compute_derivative2_fd_dev

   subroutine compute_derivative2_fv_dev(self, dir, ivar, q_gpu, d2q_ds2_gpu)
   !< Compute derivative2 of scalar fields, d2q(ivar)/ds2, using finite volume schemes.
   class(prism_fnl_object), intent(in)    :: self                                                !< The equation.
   integer(I4P),            intent(in)    :: dir                                                 !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),            intent(in)    :: ivar                                                !< Start index of variable of q.
   real(R8P),               intent(in)    :: q_gpu(1:, 1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)   !< Field variables.
   real(R8P),               intent(inout) :: d2q_ds2_gpu(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Derivative2, d2q/ds2.
   integer(I4P)                           :: i,j,k,b                                             !< Counter.
   integer(I4P)                           :: is,js,ks                                            !< Stencils.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz_gpu=>self%field_fnl%dxyz_gpu,&
             hs=>self%fdv_half_stencils(1))

   endassociate
   endsubroutine compute_derivative2_fv_dev

   subroutine compute_derivative4_fd_dev(self, dir, ivar, q_gpu, d4q_ds4_gpu)
   !< Compute derivative4 of scalar fields, d4q(ivar)/ds4, using finite difference schemes.
   class(prism_fnl_object), intent(in)    :: self                                                !< The equation.
   integer(I4P),            intent(in)    :: dir                                                 !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),            intent(in)    :: ivar                                                !< Start index of variable of q.
   real(R8P),               intent(in)    :: q_gpu(1:, 1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)   !< Field variables.
   real(R8P),               intent(inout) :: d4q_ds4_gpu(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Derivative4, d4q/ds4.
   integer(I4P)                           :: i,j,k,b                                             !< Counter.
   integer(I4P)                           :: is,js,ks                                            !< Stencils.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz_gpu=>self%field_fnl%dxyz_gpu,&
             hs=>self%fdv_half_stencils(1))

   endassociate
   endsubroutine compute_derivative4_fd_dev

   subroutine compute_divergence_fd_dev(self, ivar, ovar, q_gpu, divergence_gpu, maxdiv)
   !< Compute divergence of vector fields, div(q(ivar:ivar+2), using finite difference schemes.
   !< Directly computes divergence from transposed GPU layout (b,i,j,k,v).
   class(prism_fnl_object), intent(in)            :: self                                                      !< The equation.
   integer(I4P),            intent(in)            :: ivar
                                                                                                               !< Start index of
                                                                                                               !< field of q.
   integer(I4P),            intent(in)            :: ovar
                                                                                                               !< Output index in
                                                                                                               !< div.
   real(R8P),               intent(in)            :: q_gpu(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)          !< Field variables.
   real(R8P),               intent(inout)         :: divergence_gpu(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Divergence.
   real(R8P),               intent(out), optional :: maxdiv
                                                                                                               !< Max divergence,
                                                                                                               !< for checking.

   call  compute_divergence_fd_dev_kernel(ni             = self%ni                  ,&
                                          nj             = self%nj                  ,&
                                          nk             = self%nk                  ,&
                                          ngc            = self%ngc                 ,&
                                          blocks_number  = self%blocks_number       ,&
                                          ivar           = ivar                     ,&
                                          ovar           = ovar                     ,&
                                          s1             = self%fdv_half_stencils(1),&
                                          dxyz_gpu       = self%field_fnl%dxyz_gpu       ,&
                                          q_gpu          = q_gpu                    ,&
                                          divergence_gpu = divergence_gpu           ,&
                                          maxdiv         = maxdiv)
   contains
      subroutine compute_divergence_fd_dev_kernel(ni,nj,nk,ngc,blocks_number, &
                                                  ivar,ovar,s1,dxyz_gpu,q_gpu,divergence_gpu,maxdiv)
      !< Compute divergence, centered finite difference schemes, kernel device.
      integer(I4P), intent(in)              :: ni,nj,nk,ngc,blocks_number                                !< Grids dimensions.
      integer(I4P), intent(in)              :: ivar
                                                                                                         !< Start index of variable
                                                                                                         !< of q.
      integer(I4P), intent(in)              :: ovar                                                      !< Output index in div.
      integer(I4P), intent(in)              :: s1                                                        !< Half FDV stencil length.
      real(R8P),    intent(in)              :: dxyz_gpu(1:,1:)                                           !< Delta cells GPU [nb,3].
      real(R8P),    intent(in)              :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
                                                                                                         !< Field cell centered
                                                                                                         !< variables.
      real(R8P),    intent(inout)           :: divergence_gpu(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Divergence.
      real(R8P),    intent(out), optional   :: maxdiv
                                                                                                         !< Max divergence, for
                                                                                                         !< checking.
      integer(I4P)                          :: i,j,k,b,s                                                 !< Counter
      ! rank 1D stencil for computations on device that contiguos memory is mandatory
      real(R8P) :: qsx(1-s1:1+s1) !< X component of vector field over the x stencil.
      real(R8P) :: qsy(1-s1:1+s1) !< Y component of vector field over the y stencil.
      real(R8P) :: qsz(1-s1:1+s1) !< Z component of vector field over the z stencil.
      real(R8P) :: dxyz_b(3)      !< Per-block deltas, PRIVATE copy (no strided-section temp: issue #26 G1.b).
      real(R8P) :: maxdiv_

      maxdiv_ = -huge(1._R8P)
      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu,divergence_gpu) &
      !$acc& firstprivate(ivar,ovar,s1)                                                                &
      !$acc& private(qsx,qsy,qsz,dxyz_b) reduction(max: maxdiv_)
      !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,q_gpu,divergence_gpu) &
      !$omp& firstprivate(ivar,ovar,s1) &
      !$omp& private(qsx,qsy,qsz,dxyz_b) reduction(max: maxdiv_)
      do b=1,blocks_number
      do k=1,nk
      do j=1,nj
      do i=1,ni
         ! hoist per-block deltas into a PRIVATE vector (issue #26 G1.b, rule from #22 F1-bis):
         ! passing the strided section dxyz_gpu(b,1:3) materializes a compiler temporary that
         ! is NOT privatized -- threads race on it; benign on uniform grids (equal values),
         ! live at 2:1 level mixes. This kernel runs on mixed-level AMR topologies (div(J)
         ! diagnostics, divergence field saves).
         dxyz_b(1) = dxyz_gpu(b,1) ; dxyz_b(2) = dxyz_gpu(b,2) ; dxyz_b(3) = dxyz_gpu(b,3)
         !$acc loop seq
         do s=1-s1, 1+s1
            qsx(s) = q_gpu(b,i+s-1,j    ,k    ,ivar+0)
            qsy(s) = q_gpu(b,i    ,j+s-1,k    ,ivar+1)
            qsz(s) = q_gpu(b,i    ,j    ,k+s-1,ivar+2)
         enddo
         call compute_divergence_fd_centered_dev(s=s1,dxyz=dxyz_b,        &
                                                 qsx=qsx,qsy=qsy,qsz=qsz, &
                                                 divergence=divergence_gpu(b,i,j,k,ovar))
         maxdiv_ = max(maxdiv_, abs(divergence_gpu(b,i,j,k,ovar)))
      enddo
      enddo
      enddo
      enddo
      if (present(maxdiv)) maxdiv = maxdiv_
      endsubroutine compute_divergence_fd_dev_kernel
   endsubroutine compute_divergence_fd_dev

   subroutine compute_divergence_fv_dev(self, ivar, ovar, q_gpu, divergence_gpu, maxdiv)
   !< Compute divergence of vector fields, div(q(ivar:ivar+2), using finite volume schemes.
   !< Directly computes divergence from transposed GPU layout (b,i,j,k,v).
   class(prism_fnl_object), intent(in)            :: self                                                      !< The equation.
   integer(I4P),            intent(in)            :: ivar
                                                                                                               !< Start index of
                                                                                                               !< field of q.
   integer(I4P),            intent(in)            :: ovar
                                                                                                               !< Output index in
                                                                                                               !< divergence.
   real(R8P),               intent(in)            :: q_gpu(1:,      1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)    !< Field variables.
   real(R8P),               intent(inout)         :: divergence_gpu(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Divergence.
   real(R8P),               intent(out), optional :: maxdiv
                                                                                                               !< Max divergence,
                                                                                                               !< for checking.
   integer(I4P)                                   :: i,j,k,b,m                                                 !< Counter.
   real(R8P)                                      :: div_x, div_y, div_z
                                                                                                               !< Partial
                                                                                                               !< derivatives.
   real(R8P)                                      :: q_line(1-4:1+4)
                                                                                                               !< 1D stencil for
                                                                                                               !< reconstruction.
   real(R8P)                                      :: ql, qr
                                                                                                               !< Left/right
                                                                                                               !< reconstructions.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz_gpu=>self%field_fnl%dxyz_gpu,&
             hs=>self%fdv_half_stencils(1))

   endassociate
   endsubroutine compute_divergence_fv_dev

   subroutine compute_gradient_fd_dev(self, ivar, q_gpu, gradient_gpu)
   !< Compute gradient of scalar variable q(ivar), finite difference schemes.
   class(prism_fnl_object), intent(in)    :: self                                                    !< The equation.
   integer(I4P),            intent(in)    :: ivar                                                    !< Index of scalar var of q.
   real(R8P),               intent(in)    :: q_gpu(       1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: gradient_gpu(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Gradient.
   integer(I4P)                           :: i, j, k, b                                              !< Counter.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz_gpu=>self%field_fnl%dxyz_gpu,&
             hs=>self%fdv_half_stencils(1))

   endassociate
   endsubroutine compute_gradient_fd_dev

   subroutine compute_gradient_fv_dev(self, ivar, q_gpu, gradient_gpu)
   !< Compute gradient of scalar variable q(ivar), finite volume schemes.
   class(prism_fnl_object), intent(in)    :: self                                                    !< The equation.
   integer(I4P),            intent(in)    :: ivar                                                    !< Index of scalar var of q.
   real(R8P),               intent(in)    :: q_gpu(       1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: gradient_gpu(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Gradient.
   integer(I4P)                           :: i, j, k, b                                              !< Counter.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz_gpu=>self%field_fnl%dxyz_gpu,&
             hs=>self%fdv_half_stencils(1))

   endassociate
   endsubroutine compute_gradient_fv_dev

   subroutine compute_laplacian_fd_dev(self, ivar, q_gpu, laplacian_gpu)
   !< Compute laplacian of scalar variable q(ivar), finite difference schemes.
   class(prism_fnl_object), intent(in)    :: self                                                  !< The equation.
   integer(I4P),            intent(in)    :: ivar                                                  !< Index of scalar variable of q.
   real(R8P),               intent(in)    :: q_gpu(     1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: laplacian_gpu(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Gradient.
   integer(I4P)                           :: i, j, k, b                                            !< Counter.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz_gpu=>self%field_fnl%dxyz_gpu,&
             hs=>self%fdv_half_stencils(2))

   endassociate
   endsubroutine compute_laplacian_fd_dev

   subroutine compute_laplacian_fv_dev(self, ivar, q_gpu, laplacian_gpu)
   !< Compute laplacian of scalar variable q(ivar), finite volume schemes.
   class(prism_fnl_object), intent(in)    :: self                                                  !< The equation.
   integer(I4P),            intent(in)    :: ivar                                                  !< Index of scalar variable of q.
   real(R8P),               intent(in)    :: q_gpu(     1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: laplacian_gpu(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Gradient.
   integer(I4P)                           :: i, j, k, b                                            !< Counter.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz_gpu=>self%field_fnl%dxyz_gpu,&
             hs=>self%fdv_half_stencils(2))

   endassociate
   endsubroutine compute_laplacian_fv_dev

   ! numerical methods, space operators
   subroutine compute_residuals_fd_centered_dev(self, q_gpu, dq_gpu, s, flux_register)
   !< Compute residuals of equation, space operator, centered finite difference schemes.
   class(prism_fnl_object), intent(inout)                :: self       !< The equation.
   real(R8P),               intent(inout)                :: q_gpu(1:,         &
                                                                  1-self%ngc:,&
                                                                  1-self%ngc:,&
                                                                  1-self%ngc:,&
                                                                  1:)  !< Conservative variables.
   real(R8P),               intent(inout)                :: dq_gpu(1:,         &
                                                                   1-self%ngc:,&
                                                                   1-self%ngc:,&
                                                                   1-self%ngc:,&
                                                                   1:) !< Residuals.
   integer(I4P),            intent(in),    optional      :: s          !< Stage counter.
   class(flux_register_object), intent(inout), optional  :: flux_register !< Forest's flux register (interface parity; the
                                                                          !< FD path does not accumulate seam fluxes).

   if (present(flux_register)) continue ! FV-only machinery; accepted for interface conformance
   if (self%blocks_number > 0) then
      !call self%apply_fwl_correction(q_gpu=q_gpu)
      if (present(s)) then
         call self%update_ghost(q_gpu=q_gpu, s=s)
      else
         call self%update_ghost(q_gpu=q_gpu)
      endif
      select case(self%fd_residual_variant)
      case(FD_RESIDUAL_VARIANT_PLAIN)
         call fd_centered_plain_dev_kernel(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, &
                                           blocks_number=self%blocks_number,                   &
                                           var_jx=self%physics%var_jx, var_jy=self%physics%var_jy, &
                                           var_jz=self%physics%var_jz, s1=self%fdv_half_stencils(1), &
                                           inv_mu_scale=self%fd_inv_mu_scale, inv_eps_scale=self%fd_inv_eps_scale, &
                                           dxyz_gpu=self%field_fnl%dxyz_gpu, q_gpu=q_gpu, dq_gpu=dq_gpu)
      case(FD_RESIDUAL_VARIANT_PHI)
         call fd_centered_phi_dev_kernel(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, &
                                         blocks_number=self%blocks_number,                 &
                                         var_jx=self%physics%var_jx, var_jy=self%physics%var_jy, &
                                         var_jz=self%physics%var_jz, s1=self%fdv_half_stencils(1), &
                                         inv_mu_scale=self%fd_inv_mu_scale, inv_eps_scale=self%fd_inv_eps_scale, &
                                         chi_wave=self%fd_chi_wave, chi_damp=self%fd_chi_damp, c_r=self%physics%c_r, &
                                         ivar_phi=self%fd_ivar_phi, dxyz_gpu=self%field_fnl%dxyz_gpu, q_gpu=q_gpu, dq_gpu=dq_gpu)
      case(FD_RESIDUAL_VARIANT_PSI)
         call fd_centered_psi_dev_kernel(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, &
                                         blocks_number=self%blocks_number,                 &
                                         var_jx=self%physics%var_jx, var_jy=self%physics%var_jy, &
                                         var_jz=self%physics%var_jz, s1=self%fdv_half_stencils(1), &
                                         inv_mu_scale=self%fd_inv_mu_scale, inv_eps_scale=self%fd_inv_eps_scale, &
                                         chi_wave=self%fd_chi_wave, chi_damp=self%fd_chi_damp, c_r=self%physics%c_r, &
                                         ivar_psi=self%fd_ivar_psi, dxyz_gpu=self%field_fnl%dxyz_gpu, q_gpu=q_gpu, dq_gpu=dq_gpu)
      case(FD_RESIDUAL_VARIANT_PHI_PSI)
         call fd_centered_phi_psi_dev_kernel(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, &
                                             blocks_number=self%blocks_number,                 &
                                             var_jx=self%physics%var_jx, var_jy=self%physics%var_jy, &
                                             var_jz=self%physics%var_jz, s1=self%fdv_half_stencils(1), &
                                             inv_mu_scale=self%fd_inv_mu_scale, inv_eps_scale=self%fd_inv_eps_scale, &
                                             chi_wave=self%fd_chi_wave, chi_damp=self%fd_chi_damp, c_r=self%physics%c_r, &
                                             ivar_phi=self%fd_ivar_phi, ivar_psi=self%fd_ivar_psi, &
                                             dxyz_gpu=self%field_fnl%dxyz_gpu, q_gpu=q_gpu, dq_gpu=dq_gpu)
      case default
         call mpih_fnl%error_stop(msg=': unknown FD residual variant in FNL backend')
      end select
      if (self%pml_fnl%enabled) call self%apply_pml_fd_centered_ade_dev(q_gpu=q_gpu, dq_gpu=dq_gpu, s=s)
   endif
   contains
      subroutine compute_residuals_fd_centered_dev_kernel(ni, nj, nk, ngc, blocks_number,        &
                                                          var_Jx, var_Jy, var_Jz, nv_c, chi, c_r, s1, &
                                                          dxyz_gpu, q_gpu, dq_gpu)
      !< Compute residuals of equation, space operator, centered finite difference schemes, kernel device.
      integer(I4P), intent(in)    :: ni,nj,nk,ngc,blocks_number         !< Grids dimensions.
      integer(I4P), intent(in)    :: var_jx,var_jy,var_jz               !< Indexes of J_vec variables.
      integer(I4P), intent(in)    :: nv_c                               !< Number of conservative variables.
      real(R8P),    intent(in)    :: chi                                !< Hyperbolic correction speed.
      real(R8P),    intent(in)    :: c_r                                !< Dedner GLM parabolic-damping ratio (issue #29).
      integer(I4P), intent(in)    :: s1                                 !< Half FDV stencil length.
      real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)                    !< Delta cells GPU [nb,3].
      real(R8P),    intent(in)    :: q_gpu( 1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field cell centered variables.
      real(R8P),    intent(inout) :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Residuals.
      integer(I4P)                :: i,j,k,b,s                          !< Counter
      real(R8P)                   :: curlD(3), curlB(3)                 !< Residuals components.
      real(R8P)                   :: divergenceD, divergenceB           !< Divergence for hyperbolic correction.
      real(R8P)   			   	 :: gradphi(3), gradpsi(3) 	         !< Gradient for hyperbolic correction.
      real(R8P)   			   	 :: dxyz_b(3) 	                     !< Per-block deltas, PRIVATE copy (no strided-section temp).
      real(R8P)                   :: damping_coeff                      !< Optional GLM parabolic damping coefficient.
      ! rank 1D stencil for curl computations on device that contiguos memory is mandatory
      real(R8P) :: qsx_y(1-FDV_S_MAX:1+FDV_S_MAX) !< Y component of vector field over the x stencil.
      real(R8P) :: qsx_z(1-FDV_S_MAX:1+FDV_S_MAX) !< Z component of vector field over the x stencil.
      real(R8P) :: qsy_x(1-FDV_S_MAX:1+FDV_S_MAX) !< X component of vector field over the y stencil.
      real(R8P) :: qsy_z(1-FDV_S_MAX:1+FDV_S_MAX) !< Z component of vector field over the y stencil.
      real(R8P) :: qsz_x(1-FDV_S_MAX:1+FDV_S_MAX) !< X component of vector field over the z stencil.
      real(R8P) :: qsz_y(1-FDV_S_MAX:1+FDV_S_MAX) !< Y component of vector field over the z stencil.
      ! rank 1D stencil for divergence computations on device that contiguos memory is mandatory
      real(R8P) :: qsx_x(1-FDV_S_MAX:1+FDV_S_MAX) !< X component of vector field over the x stencil.
      real(R8P) :: qsy_y(1-FDV_S_MAX:1+FDV_S_MAX) !< Y component of vector field over the y stencil.
      real(R8P) :: qsz_z(1-FDV_S_MAX:1+FDV_S_MAX) !< Z component of vector field over the z stencil.

      if (self%physics%physical_model == EM_PHYSICAL_MODEL .or.  self%physics%physical_model == PIC_PHYSICAL_MODEL) then
		   if (self%numerics%div_corr_var == DIV_CORR_VAR_HYPER .and. &
            self%numerics%constrained_transport_D .and. &
		      .not.self%numerics%constrained_transport_B) then
            ! RHS:
            ! dD/dt = curl(B/MU0) - grad(phi) - J
            ! dB/dt = -curl(D/EPS0)
            ! dphi/dt = -ch^2*div(D)
            !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu,dq_gpu) &
            !$acc& firstprivate(var_jx,var_jy,var_jz,nv_c,chi,s1)                                    &
            !$acc& private(curlD,curlB,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y,                          &
            !$acc&         qsx_x,qsy_y,qsz_z,divergenceD,gradphi,dxyz_b)
            !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,q_gpu,dq_gpu) &
            !$omp& firstprivate(var_jx,var_jy,var_jz,nv_c,chi,s1) &
            !$omp& private(curlD,curlB,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y, &
            !$omp&         qsx_x,qsy_y,qsz_z,divergenceD,gradphi,dxyz_b)
            do b=1,blocks_number
            do k=1,nk
            do j=1,nj
            do i=1,ni
               ! per-block deltas -> PRIVATE copy (issue #26 G1.c; rationale at the plain-EM branch of this kernel)
               dxyz_b(1) = dxyz_gpu(b,1) ; dxyz_b(2) = dxyz_gpu(b,2) ; dxyz_b(3) = dxyz_gpu(b,3)
               !$acc loop seq
               do s=1-s1, 1+s1
                  qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DY)
                  qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DZ)
                  qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DX)
                  qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DZ)
                  qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DX)
                  qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DY)
               enddo
               call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_b,                                           &
                                                 qsx_y=qsx_y(1-s1:1+s1),qsx_z=qsx_z(1-s1:1+s1),qsy_x=qsy_x(1-s1:1+s1),&
                                                 qsy_z=qsy_z(1-s1:1+s1),qsz_x=qsz_x(1-s1:1+s1),qsz_y=qsz_y(1-s1:1+s1),&
                                                 curl=curlD)
               !$acc loop seq
               do s=1-s1, 1+s1
                  qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BY)
                  qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BZ)
                  qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BX)
                  qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BZ)
                  qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BX)
                  qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BY)
               enddo
               call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_b,                                           &
                                                 qsx_y=qsx_y(1-s1:1+s1),qsx_z=qsx_z(1-s1:1+s1),qsy_x=qsy_x(1-s1:1+s1),&
                                                 qsy_z=qsy_z(1-s1:1+s1),qsz_x=qsz_x(1-s1:1+s1),qsz_y=qsz_y(1-s1:1+s1),&
                                                 curl=curlB)
               !$acc loop seq
               do s=1-s1, 1+s1
                  qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DX)
                  qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DY)
                  qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DZ)
               enddo
               call compute_divergence_fd_centered_dev(s=s1,dxyz=dxyz_b,                                      &
                                                       qsx=qsx_x(1-s1:1+s1),qsy=qsy_y(1-s1:1+s1),qsz=qsz_z(1-s1:1+s1), &
                                                       divergence=divergenceD)
               !$acc loop seq
               do s=1-s1,1+s1
                  qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,nv_c)
                  qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,nv_c)
                  qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,nv_c)
               enddo
               call compute_gradient_fd_centered_dev(s=s1,dxyz=dxyz_b,                                      &
                                                     qsx=qsx_x(1-s1:1+s1),qsy=qsy_y(1-s1:1+s1),qsz=qsz_z(1-s1:1+s1), &
                                                     gradient=gradphi)

               dq_gpu(b,i,j,k,VAR_DX) =  curlB(1)/MU0 - gradphi(1) - q_gpu(b,i,j,k,var_Jx)
               dq_gpu(b,i,j,k,VAR_DY) =  curlB(2)/MU0 - gradphi(2) - q_gpu(b,i,j,k,var_Jy)
               dq_gpu(b,i,j,k,VAR_DZ) =  curlB(3)/MU0 - gradphi(3) - q_gpu(b,i,j,k,var_Jz)
               dq_gpu(b,i,j,k,VAR_BX) = -curlD(1)/EPS0
               dq_gpu(b,i,j,k,VAR_BY) = -curlD(2)/EPS0
               dq_gpu(b,i,j,k,VAR_BZ) = -curlD(3)/EPS0
               dq_gpu(b,i,j,k,nv_c) = -(chi*C0)**2*divergenceD
               if (c_r > 0._R8P) then
                  damping_coeff = chi * C0 / (c_r * minval(dxyz_b))
                  dq_gpu(b,i,j,k,nv_c) = dq_gpu(b,i,j,k,nv_c) - damping_coeff * q_gpu(b,i,j,k,nv_c)
               endif
            enddo
            enddo
            enddo
            enddo
		   elseif (self%numerics%div_corr_var == DIV_CORR_VAR_HYPER .and. &
                 .not.self%numerics%constrained_transport_D .and.       &
		           self%numerics%constrained_transport_B) then
            ! RHS:
            ! dD/dt = curl(B/MU0) - J
            ! dB/dt = -curl(D/EPS0) -grad(psi)
            ! dpsi/dt = -ch^2*div(B)
            !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu,dq_gpu) &
            !$acc& firstprivate(var_jx,var_jy,var_jz,nv_c,chi,s1)                                    &
            !$acc& private(curlD,curlB,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y,                          &
            !$acc&         qsx_x,qsy_y,qsz_z,divergenceB,gradpsi,dxyz_b)
            !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,q_gpu,dq_gpu) &
            !$omp& firstprivate(var_jx,var_jy,var_jz,nv_c,chi,s1) &
            !$omp& private(curlD,curlB,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y, &
            !$omp&         qsx_x,qsy_y,qsz_z,divergenceB,gradpsi,dxyz_b)
            do b=1,blocks_number
            do k=1,nk
            do j=1,nj
            do i=1,ni
               ! per-block deltas -> PRIVATE copy (issue #26 G1.c; rationale at the plain-EM branch of this kernel)
               dxyz_b(1) = dxyz_gpu(b,1) ; dxyz_b(2) = dxyz_gpu(b,2) ; dxyz_b(3) = dxyz_gpu(b,3)
               !$acc loop seq
               do s=1-s1, 1+s1
                  qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DY)
                  qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DZ)
                  qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DX)
                  qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DZ)
                  qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DX)
                  qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DY)
               enddo
               call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_b,                                            &
                                                 qsx_y=qsx_y(1-s1:1+s1),qsx_z=qsx_z(1-s1:1+s1),qsy_x=qsy_x(1-s1:1+s1), &
                                                 qsy_z=qsy_z(1-s1:1+s1),qsz_x=qsz_x(1-s1:1+s1),qsz_y=qsz_y(1-s1:1+s1), &
                                                 curl=curlD)
               !$acc loop seq
               do s=1-s1, 1+s1
                  qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BY)
                  qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BZ)
                  qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BX)
                  qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BZ)
                  qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BX)
                  qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BY)
               enddo
               call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_b,                                            &
                                                 qsx_y=qsx_y(1-s1:1+s1),qsx_z=qsx_z(1-s1:1+s1),qsy_x=qsy_x(1-s1:1+s1), &
                                                 qsy_z=qsy_z(1-s1:1+s1),qsz_x=qsz_x(1-s1:1+s1),qsz_y=qsz_y(1-s1:1+s1), &
                                                 curl=curlB)
               !$acc loop seq
               do s=1-s1, 1+s1
                  qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BX)
                  qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BY)
                  qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BZ)
               enddo
               call compute_divergence_fd_centered_dev(s=s1,dxyz=dxyz_b,                                      &
                                                       qsx=qsx_x(1-s1:1+s1),qsy=qsy_y(1-s1:1+s1),qsz=qsz_z(1-s1:1+s1), &
                                                       divergence=divergenceB)
               !$acc loop seq
               do s=1-s1,1+s1
                  qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,nv_c)
                  qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,nv_c)
                  qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,nv_c)
               enddo
               call compute_gradient_fd_centered_dev(s=s1,dxyz=dxyz_b,                                      &
                                                     qsx=qsx_x(1-s1:1+s1),qsy=qsy_y(1-s1:1+s1),qsz=qsz_z(1-s1:1+s1), &
                                                     gradient=gradpsi)
               dq_gpu(b,i,j,k,VAR_DX) =  curlB(1)/MU0 - q_gpu(b,i,j,k,var_Jx)
               dq_gpu(b,i,j,k,VAR_DY) =  curlB(2)/MU0 - q_gpu(b,i,j,k,var_Jy)
               dq_gpu(b,i,j,k,VAR_DZ) =  curlB(3)/MU0 - q_gpu(b,i,j,k,var_Jz)
               dq_gpu(b,i,j,k,VAR_BX) = -curlD(1)/EPS0 - gradpsi(1)
               dq_gpu(b,i,j,k,VAR_BY) = -curlD(2)/EPS0 - gradpsi(2)
               dq_gpu(b,i,j,k,VAR_BZ) = -curlD(3)/EPS0 - gradpsi(3)
               dq_gpu(b,i,j,k,nv_c) = -(chi*C0)**2*divergenceB
               if (c_r > 0._R8P) then
                  damping_coeff = chi * C0 / (c_r * minval(dxyz_b))
                  dq_gpu(b,i,j,k,nv_c) = dq_gpu(b,i,j,k,nv_c) - damping_coeff * q_gpu(b,i,j,k,nv_c)
               endif
            enddo
            enddo
            enddo
            enddo
		   elseif (self%numerics%div_corr_var == DIV_CORR_VAR_HYPER .and. &
                 self%numerics%constrained_transport_D .and. &
		           self%numerics%constrained_transport_B) then
            ! RHS:
            ! dD/dt   = curl(B/MU0) - grad(phi) - J
            ! dB/dt   = -curl(D/EPS0) -grad(psi)
            ! dphi/dt = -ch^2*div(D)
            ! dpsi/dt = -ch^2*div(B)
            !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu,dq_gpu) &
            !$acc& firstprivate(var_jx,var_jy,var_jz,nv_c,chi,s1)                                    &
            !$acc& private(curlD,curlB,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y,                          &
            !$acc&         qsx_x,qsy_y,qsz_z,divergenceD,divergenceB,gradphi,gradpsi,dxyz_b)
            !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,q_gpu,dq_gpu) &
            !$omp& firstprivate(var_jx,var_jy,var_jz,nv_c,chi,s1) &
            !$omp& private(curlD,curlB,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y, &
            !$omp&         qsx_x,qsy_y,qsz_z,divergenceD,divergenceB,gradphi,gradpsi,dxyz_b)
            do b=1,blocks_number
            do k=1,nk
            do j=1,nj
            do i=1,ni
               ! per-block deltas -> PRIVATE copy (issue #26 G1.c; rationale at the plain-EM branch of this kernel)
               dxyz_b(1) = dxyz_gpu(b,1) ; dxyz_b(2) = dxyz_gpu(b,2) ; dxyz_b(3) = dxyz_gpu(b,3)
               !$acc loop seq
               do s=1-s1, 1+s1
                  qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DY)
                  qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DZ)
                  qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DX)
                  qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DZ)
                  qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DX)
                  qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DY)
               enddo
               call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_b,                                            &
                                                 qsx_y=qsx_y(1-s1:1+s1),qsx_z=qsx_z(1-s1:1+s1),qsy_x=qsy_x(1-s1:1+s1), &
                                                 qsy_z=qsy_z(1-s1:1+s1),qsz_x=qsz_x(1-s1:1+s1),qsz_y=qsz_y(1-s1:1+s1), &
                                                 curl=curlD)
               !$acc loop seq
               do s=1-s1, 1+s1
                  qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BY)
                  qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BZ)
                  qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BX)
                  qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BZ)
                  qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BX)
                  qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BY)
               enddo
               call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_b,                                            &
                                                 qsx_y=qsx_y(1-s1:1+s1),qsx_z=qsx_z(1-s1:1+s1),qsy_x=qsy_x(1-s1:1+s1), &
                                                 qsy_z=qsy_z(1-s1:1+s1),qsz_x=qsz_x(1-s1:1+s1),qsz_y=qsz_y(1-s1:1+s1), &
                                                 curl=curlB)
               !$acc loop seq
               do s=1-s1, 1+s1
                  qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DX)
                  qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DY)
                  qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DZ)
               enddo
               call compute_divergence_fd_centered_dev(s=s1,dxyz=dxyz_b,                                      &
                                                       qsx=qsx_x(1-s1:1+s1),qsy=qsy_y(1-s1:1+s1),qsz=qsz_z(1-s1:1+s1), &
                                                       divergence=divergenceD)
               !$acc loop seq
               do s=1-s1,1+s1
                  qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,nv_c-1_I4P)
                  qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,nv_c-1_I4P)
                  qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,nv_c-1_I4P)
               enddo
               call compute_gradient_fd_centered_dev(s=s1,dxyz=dxyz_b,                                      &
                                                     qsx=qsx_x(1-s1:1+s1),qsy=qsy_y(1-s1:1+s1),qsz=qsz_z(1-s1:1+s1), &
                                                     gradient=gradphi)
               !$acc loop seq
               do s=1-s1, 1+s1
                  qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BX)
                  qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BY)
                  qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BZ)
               enddo
               call compute_divergence_fd_centered_dev(s=s1,dxyz=dxyz_b,                                      &
                                                       qsx=qsx_x(1-s1:1+s1),qsy=qsy_y(1-s1:1+s1),qsz=qsz_z(1-s1:1+s1), &
                                                       divergence=divergenceB)
               !$acc loop seq
               do s=1-s1,1+s1
                  qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,nv_c)
                  qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,nv_c)
                  qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,nv_c)
               enddo
               call compute_gradient_fd_centered_dev(s=s1,dxyz=dxyz_b,                                      &
                                                     qsx=qsx_x(1-s1:1+s1),qsy=qsy_y(1-s1:1+s1),qsz=qsz_z(1-s1:1+s1), &
                                                     gradient=gradpsi)
               dq_gpu(b,i,j,k,VAR_DX    ) =  curlB(1)/MU0  - gradphi(1) - q_gpu(b,i,j,k,var_Jx)
               dq_gpu(b,i,j,k,VAR_DY    ) =  curlB(2)/MU0  - gradphi(2) - q_gpu(b,i,j,k,var_Jy)
               dq_gpu(b,i,j,k,VAR_DZ    ) =  curlB(3)/MU0  - gradphi(3) - q_gpu(b,i,j,k,var_Jz)
               dq_gpu(b,i,j,k,VAR_BX    ) = -curlD(1)/EPS0 - gradpsi(1)
               dq_gpu(b,i,j,k,VAR_BY    ) = -curlD(2)/EPS0 - gradpsi(2)
               dq_gpu(b,i,j,k,VAR_BZ    ) = -curlD(3)/EPS0 - gradpsi(3)
               dq_gpu(b,i,j,k,nv_c-1_I4P) = -(chi*C0)**2*divergenceD
               dq_gpu(b,i,j,k,nv_c)       = -(chi*C0)**2*divergenceB
               if (c_r > 0._R8P) then
                  damping_coeff = chi * C0 / (c_r * minval(dxyz_b))
                  dq_gpu(b,i,j,k,nv_c-1_I4P) = dq_gpu(b,i,j,k,nv_c-1_I4P) - damping_coeff * q_gpu(b,i,j,k,nv_c-1_I4P)
                  dq_gpu(b,i,j,k,nv_c)       = dq_gpu(b,i,j,k,nv_c)       - damping_coeff * q_gpu(b,i,j,k,nv_c)
               endif
            enddo
            enddo
            enddo
            enddo
         else
            ! RHS:
            ! dD/dt = curl(B/MU0) - J
            ! dB/dt = -curl(D/EPS0)
            !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu,dq_gpu) &
            !$acc& firstprivate(var_jx,var_jy,var_jz,s1)                                             &
            !$acc& private(curlD,curlB,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y,dxyz_b)
            !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,q_gpu,dq_gpu) &
            !$omp& firstprivate(var_jx,var_jy,var_jz,s1) &
            !$omp& private(curlD,curlB,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y,dxyz_b)
            do b=1,blocks_number
            do k=1,nk
            do j=1,nj
            do i=1,ni
               ! hoist per-block deltas into a PRIVATE vector (issue #22 F1-bis): passing the
               ! strided section dxyz_gpu(b,1:3) materializes a compiler temporary that is NOT
               ! privatized -- threads race on it; benign on uniform grids (equal values),
               ! live at 2:1 level mixes.
               dxyz_b(1) = dxyz_gpu(b,1) ; dxyz_b(2) = dxyz_gpu(b,2) ; dxyz_b(3) = dxyz_gpu(b,3)
               !$acc loop seq
               do s=1-s1, 1+s1
                  qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DY)
                  qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DZ)
                  qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DX)
                  qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DZ)
                  qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DX)
                  qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DY)
               enddo
               call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_b,                                                      &
                                                 qsx_y=qsx_y(1-s1:1+s1),qsx_z=qsx_z(1-s1:1+s1),qsy_x=qsy_x(1-s1:1+s1), &
                                                 qsy_z=qsy_z(1-s1:1+s1),qsz_x=qsz_x(1-s1:1+s1),qsz_y=qsz_y(1-s1:1+s1), &
                                                 curl=curlD)
               !$acc loop seq
               do s=1-s1, 1+s1
                  qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BY)
                  qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BZ)
                  qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BX)
                  qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BZ)
                  qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BX)
                  qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BY)
               enddo
               call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_b,                                                      &
                                                 qsx_y=qsx_y(1-s1:1+s1),qsx_z=qsx_z(1-s1:1+s1),qsy_x=qsy_x(1-s1:1+s1), &
                                                 qsy_z=qsy_z(1-s1:1+s1),qsz_x=qsz_x(1-s1:1+s1),qsz_y=qsz_y(1-s1:1+s1), &
                                                 curl=curlB)
               dq_gpu(b,i,j,k,VAR_DX) =  curlB(1)/MU0 - q_gpu(b,i,j,k,var_Jx)
               dq_gpu(b,i,j,k,VAR_DY) =  curlB(2)/MU0 - q_gpu(b,i,j,k,var_Jy)
               dq_gpu(b,i,j,k,VAR_DZ) =  curlB(3)/MU0 - q_gpu(b,i,j,k,var_Jz)
               dq_gpu(b,i,j,k,VAR_BX) = -curlD(1)/EPS0
               dq_gpu(b,i,j,k,VAR_BY) = -curlD(2)/EPS0
               dq_gpu(b,i,j,k,VAR_BZ) = -curlD(3)/EPS0
            enddo
            enddo
            enddo
            enddo
         endif
      elseif (self%physics%physical_model == ADIM_EM_PHYSICAL_MODEL) then
         	if (self%numerics%div_corr_var == DIV_CORR_VAR_HYPER .and. &
            self%numerics%constrained_transport_D .and. &
		      .not.self%numerics%constrained_transport_B) then
            ! RHS:
            ! dD/dt = curl(B) - grad(phi) - J
            ! dB/dt = -curl(D)
            ! dphi/dt = -ch^2*div(D)
            !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu,dq_gpu) &
            !$acc& firstprivate(var_jx,var_jy,var_jz,nv_c,chi,s1)                                    &
            !$acc& private(curlD,curlB,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y,                          &
            !$acc&         qsx_x,qsy_y,qsz_z,divergenceD,gradphi,dxyz_b)
            !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,q_gpu,dq_gpu) &
            !$omp& firstprivate(var_jx,var_jy,var_jz,nv_c,chi,s1) &
            !$omp& private(curlD,curlB,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y, &
            !$omp&         qsx_x,qsy_y,qsz_z,divergenceD,gradphi,dxyz_b)
            do b=1,blocks_number
            do k=1,nk
            do j=1,nj
            do i=1,ni
               ! per-block deltas -> PRIVATE copy (issue #26 G1.c; rationale at the plain-EM branch of this kernel)
               dxyz_b(1) = dxyz_gpu(b,1) ; dxyz_b(2) = dxyz_gpu(b,2) ; dxyz_b(3) = dxyz_gpu(b,3)
               !$acc loop seq
               do s=1-s1, 1+s1
                  qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DY)
                  qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DZ)
                  qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DX)
                  qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DZ)
                  qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DX)
                  qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DY)
               enddo
               call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_b,                                           &
                                                 qsx_y=qsx_y(1-s1:1+s1),qsx_z=qsx_z(1-s1:1+s1),qsy_x=qsy_x(1-s1:1+s1),&
                                                 qsy_z=qsy_z(1-s1:1+s1),qsz_x=qsz_x(1-s1:1+s1),qsz_y=qsz_y(1-s1:1+s1),&
                                                 curl=curlD)
               !$acc loop seq
               do s=1-s1, 1+s1
                  qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BY)
                  qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BZ)
                  qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BX)
                  qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BZ)
                  qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BX)
                  qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BY)
               enddo
               call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_b,                                           &
                                                 qsx_y=qsx_y(1-s1:1+s1),qsx_z=qsx_z(1-s1:1+s1),qsy_x=qsy_x(1-s1:1+s1),&
                                                 qsy_z=qsy_z(1-s1:1+s1),qsz_x=qsz_x(1-s1:1+s1),qsz_y=qsz_y(1-s1:1+s1),&
                                                 curl=curlB)
               !$acc loop seq
               do s=1-s1, 1+s1
                  qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DX)
                  qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DY)
                  qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DZ)
               enddo
               call compute_divergence_fd_centered_dev(s=s1,dxyz=dxyz_b,                                      &
                                                       qsx=qsx_x(1-s1:1+s1),qsy=qsy_y(1-s1:1+s1),qsz=qsz_z(1-s1:1+s1), &
                                                       divergence=divergenceD)
               !$acc loop seq
               do s=1-s1,1+s1
                  qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,nv_c)
                  qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,nv_c)
                  qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,nv_c)
               enddo
               call compute_gradient_fd_centered_dev(s=s1,dxyz=dxyz_b,                                      &
                                                     qsx=qsx_x(1-s1:1+s1),qsy=qsy_y(1-s1:1+s1),qsz=qsz_z(1-s1:1+s1), &
                                                     gradient=gradphi)

               dq_gpu(b,i,j,k,VAR_DX) =  curlB(1) - gradphi(1) - q_gpu(b,i,j,k,var_Jx)
               dq_gpu(b,i,j,k,VAR_DY) =  curlB(2) - gradphi(2) - q_gpu(b,i,j,k,var_Jy)
               dq_gpu(b,i,j,k,VAR_DZ) =  curlB(3) - gradphi(3) - q_gpu(b,i,j,k,var_Jz)
               dq_gpu(b,i,j,k,VAR_BX) = -curlD(1)
               dq_gpu(b,i,j,k,VAR_BY) = -curlD(2)
               dq_gpu(b,i,j,k,VAR_BZ) = -curlD(3)
               dq_gpu(b,i,j,k,nv_c) = -(chi)**2*divergenceD
               if (c_r > 0._R8P) then
                  damping_coeff = chi / (c_r * minval(dxyz_b))
                  dq_gpu(b,i,j,k,nv_c) = dq_gpu(b,i,j,k,nv_c) - damping_coeff * q_gpu(b,i,j,k,nv_c)
               endif
            enddo
            enddo
            enddo
            enddo
		   elseif (self%numerics%div_corr_var == DIV_CORR_VAR_HYPER .and. &
                 .not.self%numerics%constrained_transport_D .and.       &
		           self%numerics%constrained_transport_B) then
            ! RHS:
            ! dD/dt = curl(B) - J
            ! dB/dt = -curl(D) -grad(psi)
            ! dpsi/dt = -ch^2*div(B)
            !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu,dq_gpu) &
            !$acc& firstprivate(var_jx,var_jy,var_jz,nv_c,chi,s1)                                    &
            !$acc& private(curlD,curlB,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y,                          &
            !$acc&         qsx_x,qsy_y,qsz_z,divergenceB,gradpsi,dxyz_b)
            !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,q_gpu,dq_gpu) &
            !$omp& firstprivate(var_jx,var_jy,var_jz,nv_c,chi,s1) &
            !$omp& private(curlD,curlB,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y, &
            !$omp&         qsx_x,qsy_y,qsz_z,divergenceB,gradpsi,dxyz_b)
            do b=1,blocks_number
            do k=1,nk
            do j=1,nj
            do i=1,ni
               ! per-block deltas -> PRIVATE copy (issue #26 G1.c; rationale at the plain-EM branch of this kernel)
               dxyz_b(1) = dxyz_gpu(b,1) ; dxyz_b(2) = dxyz_gpu(b,2) ; dxyz_b(3) = dxyz_gpu(b,3)
               !$acc loop seq
               do s=1-s1, 1+s1
                  qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DY)
                  qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DZ)
                  qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DX)
                  qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DZ)
                  qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DX)
                  qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DY)
               enddo
               call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_b,                                            &
                                                 qsx_y=qsx_y(1-s1:1+s1),qsx_z=qsx_z(1-s1:1+s1),qsy_x=qsy_x(1-s1:1+s1), &
                                                 qsy_z=qsy_z(1-s1:1+s1),qsz_x=qsz_x(1-s1:1+s1),qsz_y=qsz_y(1-s1:1+s1), &
                                                 curl=curlD)
               !$acc loop seq
               do s=1-s1, 1+s1
                  qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BY)
                  qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BZ)
                  qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BX)
                  qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BZ)
                  qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BX)
                  qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BY)
               enddo
               call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_b,                                            &
                                                 qsx_y=qsx_y(1-s1:1+s1),qsx_z=qsx_z(1-s1:1+s1),qsy_x=qsy_x(1-s1:1+s1), &
                                                 qsy_z=qsy_z(1-s1:1+s1),qsz_x=qsz_x(1-s1:1+s1),qsz_y=qsz_y(1-s1:1+s1), &
                                                 curl=curlB)
               !$acc loop seq
               do s=1-s1, 1+s1
                  qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BX)
                  qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BY)
                  qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BZ)
               enddo
               call compute_divergence_fd_centered_dev(s=s1,dxyz=dxyz_b,                                      &
                                                       qsx=qsx_x(1-s1:1+s1),qsy=qsy_y(1-s1:1+s1),qsz=qsz_z(1-s1:1+s1), &
                                                       divergence=divergenceB)
               !$acc loop seq
               do s=1-s1,1+s1
                  qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,nv_c)
                  qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,nv_c)
                  qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,nv_c)
               enddo
               call compute_gradient_fd_centered_dev(s=s1,dxyz=dxyz_b,                                      &
                                                     qsx=qsx_x(1-s1:1+s1),qsy=qsy_y(1-s1:1+s1),qsz=qsz_z(1-s1:1+s1), &
                                                     gradient=gradpsi)
               dq_gpu(b,i,j,k,VAR_DX) =  curlB(1) - q_gpu(b,i,j,k,var_Jx)
               dq_gpu(b,i,j,k,VAR_DY) =  curlB(2) - q_gpu(b,i,j,k,var_Jy)
               dq_gpu(b,i,j,k,VAR_DZ) =  curlB(3) - q_gpu(b,i,j,k,var_Jz)
               dq_gpu(b,i,j,k,VAR_BX) = -curlD(1) - gradpsi(1)
               dq_gpu(b,i,j,k,VAR_BY) = -curlD(2) - gradpsi(2)
               dq_gpu(b,i,j,k,VAR_BZ) = -curlD(3) - gradpsi(3)
               dq_gpu(b,i,j,k,nv_c) = -(chi)**2*divergenceB
               if (c_r > 0._R8P) then
                  damping_coeff = chi / (c_r * minval(dxyz_b))
                  dq_gpu(b,i,j,k,nv_c) = dq_gpu(b,i,j,k,nv_c) - damping_coeff * q_gpu(b,i,j,k,nv_c)
               endif
            enddo
            enddo
            enddo
            enddo
		   elseif (self%numerics%div_corr_var == DIV_CORR_VAR_HYPER .and. &
                 self%numerics%constrained_transport_D .and. &
		           self%numerics%constrained_transport_B) then
            ! RHS:
            ! dD/dt   = curl(B) - grad(phi) - J
            ! dB/dt   = -curl(D) -grad(psi)
            ! dphi/dt = -ch^2*div(D)
            ! dpsi/dt = -ch^2*div(B)
            !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu,dq_gpu) &
            !$acc& firstprivate(var_jx,var_jy,var_jz,nv_c,chi,s1)                                    &
            !$acc& private(curlD,curlB,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y,                          &
            !$acc&         qsx_x,qsy_y,qsz_z,divergenceD,divergenceB,gradphi,gradpsi,dxyz_b)
            !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,q_gpu,dq_gpu) &
            !$omp& firstprivate(var_jx,var_jy,var_jz,nv_c,chi,s1) &
            !$omp& private(curlD,curlB,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y, &
            !$omp&         qsx_x,qsy_y,qsz_z,divergenceD,divergenceB,gradphi,gradpsi,dxyz_b)
            do b=1,blocks_number
            do k=1,nk
            do j=1,nj
            do i=1,ni
               ! per-block deltas -> PRIVATE copy (issue #26 G1.c; rationale at the plain-EM branch of this kernel)
               dxyz_b(1) = dxyz_gpu(b,1) ; dxyz_b(2) = dxyz_gpu(b,2) ; dxyz_b(3) = dxyz_gpu(b,3)
               !$acc loop seq
               do s=1-s1, 1+s1
                  qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DY)
                  qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DZ)
                  qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DX)
                  qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DZ)
                  qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DX)
                  qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DY)
               enddo
               call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_b,                                            &
                                                 qsx_y=qsx_y(1-s1:1+s1),qsx_z=qsx_z(1-s1:1+s1),qsy_x=qsy_x(1-s1:1+s1), &
                                                 qsy_z=qsy_z(1-s1:1+s1),qsz_x=qsz_x(1-s1:1+s1),qsz_y=qsz_y(1-s1:1+s1), &
                                                 curl=curlD)
               !$acc loop seq
               do s=1-s1, 1+s1
                  qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BY)
                  qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BZ)
                  qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BX)
                  qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BZ)
                  qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BX)
                  qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BY)
               enddo
               call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_b,                                            &
                                                 qsx_y=qsx_y(1-s1:1+s1),qsx_z=qsx_z(1-s1:1+s1),qsy_x=qsy_x(1-s1:1+s1), &
                                                 qsy_z=qsy_z(1-s1:1+s1),qsz_x=qsz_x(1-s1:1+s1),qsz_y=qsz_y(1-s1:1+s1), &
                                                 curl=curlB)
               !$acc loop seq
               do s=1-s1, 1+s1
                  qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DX)
                  qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DY)
                  qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DZ)
               enddo
               call compute_divergence_fd_centered_dev(s=s1,dxyz=dxyz_b,                                      &
                                                       qsx=qsx_x(1-s1:1+s1),qsy=qsy_y(1-s1:1+s1),qsz=qsz_z(1-s1:1+s1), &
                                                       divergence=divergenceD)
               !$acc loop seq
               do s=1-s1,1+s1
                  qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,nv_c-1_I4P)
                  qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,nv_c-1_I4P)
                  qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,nv_c-1_I4P)
               enddo
               call compute_gradient_fd_centered_dev(s=s1,dxyz=dxyz_b,                                      &
                                                     qsx=qsx_x(1-s1:1+s1),qsy=qsy_y(1-s1:1+s1),qsz=qsz_z(1-s1:1+s1), &
                                                     gradient=gradphi)
               !$acc loop seq
               do s=1-s1, 1+s1
                  qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BX)
                  qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BY)
                  qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BZ)
               enddo
               call compute_divergence_fd_centered_dev(s=s1,dxyz=dxyz_b,                                      &
                                                       qsx=qsx_x(1-s1:1+s1),qsy=qsy_y(1-s1:1+s1),qsz=qsz_z(1-s1:1+s1), &
                                                       divergence=divergenceB)
               !$acc loop seq
               do s=1-s1,1+s1
                  qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,nv_c)
                  qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,nv_c)
                  qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,nv_c)
               enddo
               call compute_gradient_fd_centered_dev(s=s1,dxyz=dxyz_b,                                      &
                                                     qsx=qsx_x(1-s1:1+s1),qsy=qsy_y(1-s1:1+s1),qsz=qsz_z(1-s1:1+s1), &
                                                     gradient=gradpsi)
               dq_gpu(b,i,j,k,VAR_DX    ) =  curlB(1) - gradphi(1) - q_gpu(b,i,j,k,var_Jx)
               dq_gpu(b,i,j,k,VAR_DY    ) =  curlB(2) - gradphi(2) - q_gpu(b,i,j,k,var_Jy)
               dq_gpu(b,i,j,k,VAR_DZ    ) =  curlB(3) - gradphi(3) - q_gpu(b,i,j,k,var_Jz)
               dq_gpu(b,i,j,k,VAR_BX    ) = -curlD(1) - gradpsi(1)
               dq_gpu(b,i,j,k,VAR_BY    ) = -curlD(2) - gradpsi(2)
               dq_gpu(b,i,j,k,VAR_BZ    ) = -curlD(3) - gradpsi(3)
               dq_gpu(b,i,j,k,nv_c-1_I4P) = -(chi)**2*divergenceD
               dq_gpu(b,i,j,k,nv_c)       = -(chi)**2*divergenceB
               if (c_r > 0._R8P) then
                  damping_coeff = chi / (c_r * minval(dxyz_b))
                  dq_gpu(b,i,j,k,nv_c-1_I4P) = dq_gpu(b,i,j,k,nv_c-1_I4P) - damping_coeff * q_gpu(b,i,j,k,nv_c-1_I4P)
                  dq_gpu(b,i,j,k,nv_c)       = dq_gpu(b,i,j,k,nv_c)       - damping_coeff * q_gpu(b,i,j,k,nv_c)
               endif
            enddo
            enddo
            enddo
            enddo
         else
            ! RHS:
            ! dD/dt = curl(B) - J
            ! dB/dt = -curl(D)
            !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu,dq_gpu) &
            !$acc& firstprivate(var_jx,var_jy,var_jz,s1)                                             &
            !$acc& private(curlD,curlB,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y,dxyz_b)
            !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,q_gpu,dq_gpu) &
            !$omp& firstprivate(var_jx,var_jy,var_jz,s1) &
            !$omp& private(curlD,curlB,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y,dxyz_b)
            do b=1,blocks_number
            do k=1,nk
            do j=1,nj
            do i=1,ni
               ! per-block deltas -> PRIVATE copy (issue #26 G1.c; rationale at the plain-EM branch of this kernel)
               dxyz_b(1) = dxyz_gpu(b,1) ; dxyz_b(2) = dxyz_gpu(b,2) ; dxyz_b(3) = dxyz_gpu(b,3)
               !$acc loop seq
               do s=1-s1, 1+s1
                  qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DY)
                  qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DZ)
                  qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DX)
                  qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DZ)
                  qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DX)
                  qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DY)
               enddo
               call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_b,                                            &
                                                 qsx_y=qsx_y(1-s1:1+s1),qsx_z=qsx_z(1-s1:1+s1),qsy_x=qsy_x(1-s1:1+s1), &
                                                 qsy_z=qsy_z(1-s1:1+s1),qsz_x=qsz_x(1-s1:1+s1),qsz_y=qsz_y(1-s1:1+s1), &
                                                 curl=curlD)
               !$acc loop seq
               do s=1-s1, 1+s1
                  qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BY)
                  qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BZ)
                  qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BX)
                  qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BZ)
                  qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BX)
                  qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BY)
               enddo
               call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_b,                                            &
                                                 qsx_y=qsx_y(1-s1:1+s1),qsx_z=qsx_z(1-s1:1+s1),qsy_x=qsy_x(1-s1:1+s1), &
                                                 qsy_z=qsy_z(1-s1:1+s1),qsz_x=qsz_x(1-s1:1+s1),qsz_y=qsz_y(1-s1:1+s1), &
                                                 curl=curlB)
               dq_gpu(b,i,j,k,VAR_DX) =  curlB(1) - q_gpu(b,i,j,k,var_Jx)
               dq_gpu(b,i,j,k,VAR_DY) =  curlB(2) - q_gpu(b,i,j,k,var_Jy)
               dq_gpu(b,i,j,k,VAR_DZ) =  curlB(3) - q_gpu(b,i,j,k,var_Jz)
               dq_gpu(b,i,j,k,VAR_BX) = -curlD(1)
               dq_gpu(b,i,j,k,VAR_BY) = -curlD(2)
               dq_gpu(b,i,j,k,VAR_BZ) = -curlD(3)
            enddo
            enddo
            enddo
            enddo
         endif
      endif
      endsubroutine compute_residuals_fd_centered_dev_kernel
   endsubroutine compute_residuals_fd_centered_dev

   subroutine apply_pml_fd_centered_ade_dev(self, q_gpu, dq_gpu, s)
   !< Add ADE-PML source terms to the FD-centered Maxwell residuals and build the device-side auxiliary RHSs.
   class(prism_fnl_object), intent(inout) :: self
   real(R8P),               intent(in)    :: q_gpu(1:,         &
                                                   1-self%ngc:, &
                                                   1-self%ngc:, &
                                                   1-self%ngc:, &
                                                   1:)
   real(R8P),               intent(inout) :: dq_gpu(1:,         &
                                                    1-self%ngc:, &
                                                    1-self%ngc:, &
                                                    1-self%ngc:, &
                                                    1:)
   integer(I4P),            intent(in), optional :: s
   integer(I4P)                           :: ni
   integer(I4P)                           :: nj
   integer(I4P)                           :: nk
   integer(I4P)                           :: s1
   real(R8P)                              :: inv_eps_scale
   real(R8P)                              :: inv_mu_scale
   real(R8P), pointer                     :: dxyz_gpu(:,:)

   if (.not. self%pml_fnl%enabled) return
   ni = self%ni
   nj = self%nj
   nk = self%nk
   s1 = self%fdv_half_stencils(1)
   dxyz_gpu => self%field_fnl%dxyz_gpu

   select case (self%physics%physical_model)
   case (ADIM_EM_PHYSICAL_MODEL)
      inv_eps_scale = 1._R8P
      inv_mu_scale  = 1._R8P
   case default
      inv_eps_scale = 1._R8P / EPS0
      inv_mu_scale  = 1._R8P / MU0
   endselect

   call self%rk_pml_fnl%reset_rhs()

   if (present(s)) then
      if (associated(self%rk_pml_fnl%q_pml_x_m_rk_gpu)) then
         call apply_x_face_field_correction_dev(q_face_gpu=self%rk_pml_fnl%q_pml_x_m_rk_gpu(:,:,:,:,:,s),           &
                                                blocks_gpu=self%pml_fnl%blocks_x_m_gpu, start_gpu=self%pml_fnl%start_x_m_gpu, &
                                                cells_gpu=self%pml_fnl%cells_x_m_gpu, kappa_gpu=self%pml_fnl%kappa_x_m_gpu)
         call compute_x_face_pml_rhs_dev(q_face_gpu=self%rk_pml_fnl%q_pml_x_m_rk_gpu(:,:,:,:,:,s),                   &
                                         dq_face_gpu=self%rk_pml_fnl%dq_pml_x_m_gpu, blocks_gpu=self%pml_fnl%blocks_x_m_gpu, &
                                         start_gpu=self%pml_fnl%start_x_m_gpu, cells_gpu=self%pml_fnl%cells_x_m_gpu,       &
                                         gamma_gpu=self%pml_fnl%gamma_x_m_gpu, alpha_gpu=self%pml_fnl%alpha_x_m_gpu,       &
                                         kappa_gpu=self%pml_fnl%kappa_x_m_gpu)
      endif
      if (associated(self%rk_pml_fnl%q_pml_x_p_rk_gpu)) then
         call apply_x_face_field_correction_dev(q_face_gpu=self%rk_pml_fnl%q_pml_x_p_rk_gpu(:,:,:,:,:,s),           &
                                                blocks_gpu=self%pml_fnl%blocks_x_p_gpu, start_gpu=self%pml_fnl%start_x_p_gpu, &
                                                cells_gpu=self%pml_fnl%cells_x_p_gpu, kappa_gpu=self%pml_fnl%kappa_x_p_gpu)
         call compute_x_face_pml_rhs_dev(q_face_gpu=self%rk_pml_fnl%q_pml_x_p_rk_gpu(:,:,:,:,:,s),                   &
                                         dq_face_gpu=self%rk_pml_fnl%dq_pml_x_p_gpu, blocks_gpu=self%pml_fnl%blocks_x_p_gpu, &
                                         start_gpu=self%pml_fnl%start_x_p_gpu, cells_gpu=self%pml_fnl%cells_x_p_gpu,       &
                                         gamma_gpu=self%pml_fnl%gamma_x_p_gpu, alpha_gpu=self%pml_fnl%alpha_x_p_gpu,       &
                                         kappa_gpu=self%pml_fnl%kappa_x_p_gpu)
      endif
      if (associated(self%rk_pml_fnl%q_pml_y_m_rk_gpu)) then
         call apply_y_face_field_correction_dev(q_face_gpu=self%rk_pml_fnl%q_pml_y_m_rk_gpu(:,:,:,:,:,s),           &
                                                blocks_gpu=self%pml_fnl%blocks_y_m_gpu, start_gpu=self%pml_fnl%start_y_m_gpu, &
                                                cells_gpu=self%pml_fnl%cells_y_m_gpu, kappa_gpu=self%pml_fnl%kappa_y_m_gpu)
         call compute_y_face_pml_rhs_dev(q_face_gpu=self%rk_pml_fnl%q_pml_y_m_rk_gpu(:,:,:,:,:,s),                   &
                                         dq_face_gpu=self%rk_pml_fnl%dq_pml_y_m_gpu, blocks_gpu=self%pml_fnl%blocks_y_m_gpu, &
                                         start_gpu=self%pml_fnl%start_y_m_gpu, cells_gpu=self%pml_fnl%cells_y_m_gpu,       &
                                         gamma_gpu=self%pml_fnl%gamma_y_m_gpu, alpha_gpu=self%pml_fnl%alpha_y_m_gpu,       &
                                         kappa_gpu=self%pml_fnl%kappa_y_m_gpu)
      endif
      if (associated(self%rk_pml_fnl%q_pml_y_p_rk_gpu)) then
         call apply_y_face_field_correction_dev(q_face_gpu=self%rk_pml_fnl%q_pml_y_p_rk_gpu(:,:,:,:,:,s),           &
                                                blocks_gpu=self%pml_fnl%blocks_y_p_gpu, start_gpu=self%pml_fnl%start_y_p_gpu, &
                                                cells_gpu=self%pml_fnl%cells_y_p_gpu, kappa_gpu=self%pml_fnl%kappa_y_p_gpu)
         call compute_y_face_pml_rhs_dev(q_face_gpu=self%rk_pml_fnl%q_pml_y_p_rk_gpu(:,:,:,:,:,s),                   &
                                         dq_face_gpu=self%rk_pml_fnl%dq_pml_y_p_gpu, blocks_gpu=self%pml_fnl%blocks_y_p_gpu, &
                                         start_gpu=self%pml_fnl%start_y_p_gpu, cells_gpu=self%pml_fnl%cells_y_p_gpu,       &
                                         gamma_gpu=self%pml_fnl%gamma_y_p_gpu, alpha_gpu=self%pml_fnl%alpha_y_p_gpu,       &
                                         kappa_gpu=self%pml_fnl%kappa_y_p_gpu)
      endif
      if (associated(self%rk_pml_fnl%q_pml_z_m_rk_gpu)) then
         call apply_z_face_field_correction_dev(q_face_gpu=self%rk_pml_fnl%q_pml_z_m_rk_gpu(:,:,:,:,:,s),           &
                                                blocks_gpu=self%pml_fnl%blocks_z_m_gpu, start_gpu=self%pml_fnl%start_z_m_gpu, &
                                                cells_gpu=self%pml_fnl%cells_z_m_gpu, kappa_gpu=self%pml_fnl%kappa_z_m_gpu)
         call compute_z_face_pml_rhs_dev(q_face_gpu=self%rk_pml_fnl%q_pml_z_m_rk_gpu(:,:,:,:,:,s),                   &
                                         dq_face_gpu=self%rk_pml_fnl%dq_pml_z_m_gpu, blocks_gpu=self%pml_fnl%blocks_z_m_gpu, &
                                         start_gpu=self%pml_fnl%start_z_m_gpu, cells_gpu=self%pml_fnl%cells_z_m_gpu,       &
                                         gamma_gpu=self%pml_fnl%gamma_z_m_gpu, alpha_gpu=self%pml_fnl%alpha_z_m_gpu,       &
                                         kappa_gpu=self%pml_fnl%kappa_z_m_gpu)
      endif
      if (associated(self%rk_pml_fnl%q_pml_z_p_rk_gpu)) then
         call apply_z_face_field_correction_dev(q_face_gpu=self%rk_pml_fnl%q_pml_z_p_rk_gpu(:,:,:,:,:,s),           &
                                                blocks_gpu=self%pml_fnl%blocks_z_p_gpu, start_gpu=self%pml_fnl%start_z_p_gpu, &
                                                cells_gpu=self%pml_fnl%cells_z_p_gpu, kappa_gpu=self%pml_fnl%kappa_z_p_gpu)
         call compute_z_face_pml_rhs_dev(q_face_gpu=self%rk_pml_fnl%q_pml_z_p_rk_gpu(:,:,:,:,:,s),                   &
                                         dq_face_gpu=self%rk_pml_fnl%dq_pml_z_p_gpu, blocks_gpu=self%pml_fnl%blocks_z_p_gpu, &
                                         start_gpu=self%pml_fnl%start_z_p_gpu, cells_gpu=self%pml_fnl%cells_z_p_gpu,       &
                                         gamma_gpu=self%pml_fnl%gamma_z_p_gpu, alpha_gpu=self%pml_fnl%alpha_z_p_gpu,       &
                                         kappa_gpu=self%pml_fnl%kappa_z_p_gpu)
      endif
   else
      if (associated(self%pml_fnl%q_pml_x_m_gpu)) then
         call apply_x_face_field_correction_dev(q_face_gpu=self%pml_fnl%q_pml_x_m_gpu, blocks_gpu=self%pml_fnl%blocks_x_m_gpu, &
                                                start_gpu=self%pml_fnl%start_x_m_gpu, cells_gpu=self%pml_fnl%cells_x_m_gpu,     &
                                                kappa_gpu=self%pml_fnl%kappa_x_m_gpu)
         call compute_x_face_pml_rhs_dev(q_face_gpu=self%pml_fnl%q_pml_x_m_gpu, dq_face_gpu=self%rk_pml_fnl%dq_pml_x_m_gpu, &
                                         blocks_gpu=self%pml_fnl%blocks_x_m_gpu, start_gpu=self%pml_fnl%start_x_m_gpu,      &
                                         cells_gpu=self%pml_fnl%cells_x_m_gpu, gamma_gpu=self%pml_fnl%gamma_x_m_gpu,        &
                                         alpha_gpu=self%pml_fnl%alpha_x_m_gpu, kappa_gpu=self%pml_fnl%kappa_x_m_gpu)
      endif
      if (associated(self%pml_fnl%q_pml_x_p_gpu)) then
         call apply_x_face_field_correction_dev(q_face_gpu=self%pml_fnl%q_pml_x_p_gpu, blocks_gpu=self%pml_fnl%blocks_x_p_gpu, &
                                                start_gpu=self%pml_fnl%start_x_p_gpu, cells_gpu=self%pml_fnl%cells_x_p_gpu,     &
                                                kappa_gpu=self%pml_fnl%kappa_x_p_gpu)
         call compute_x_face_pml_rhs_dev(q_face_gpu=self%pml_fnl%q_pml_x_p_gpu, dq_face_gpu=self%rk_pml_fnl%dq_pml_x_p_gpu, &
                                         blocks_gpu=self%pml_fnl%blocks_x_p_gpu, start_gpu=self%pml_fnl%start_x_p_gpu,      &
                                         cells_gpu=self%pml_fnl%cells_x_p_gpu, gamma_gpu=self%pml_fnl%gamma_x_p_gpu,        &
                                         alpha_gpu=self%pml_fnl%alpha_x_p_gpu, kappa_gpu=self%pml_fnl%kappa_x_p_gpu)
      endif
      if (associated(self%pml_fnl%q_pml_y_m_gpu)) then
         call apply_y_face_field_correction_dev(q_face_gpu=self%pml_fnl%q_pml_y_m_gpu, blocks_gpu=self%pml_fnl%blocks_y_m_gpu, &
                                                start_gpu=self%pml_fnl%start_y_m_gpu, cells_gpu=self%pml_fnl%cells_y_m_gpu,     &
                                                kappa_gpu=self%pml_fnl%kappa_y_m_gpu)
         call compute_y_face_pml_rhs_dev(q_face_gpu=self%pml_fnl%q_pml_y_m_gpu, dq_face_gpu=self%rk_pml_fnl%dq_pml_y_m_gpu, &
                                         blocks_gpu=self%pml_fnl%blocks_y_m_gpu, start_gpu=self%pml_fnl%start_y_m_gpu,      &
                                         cells_gpu=self%pml_fnl%cells_y_m_gpu, gamma_gpu=self%pml_fnl%gamma_y_m_gpu,        &
                                         alpha_gpu=self%pml_fnl%alpha_y_m_gpu, kappa_gpu=self%pml_fnl%kappa_y_m_gpu)
      endif
      if (associated(self%pml_fnl%q_pml_y_p_gpu)) then
         call apply_y_face_field_correction_dev(q_face_gpu=self%pml_fnl%q_pml_y_p_gpu, blocks_gpu=self%pml_fnl%blocks_y_p_gpu, &
                                                start_gpu=self%pml_fnl%start_y_p_gpu, cells_gpu=self%pml_fnl%cells_y_p_gpu,     &
                                                kappa_gpu=self%pml_fnl%kappa_y_p_gpu)
         call compute_y_face_pml_rhs_dev(q_face_gpu=self%pml_fnl%q_pml_y_p_gpu, dq_face_gpu=self%rk_pml_fnl%dq_pml_y_p_gpu, &
                                         blocks_gpu=self%pml_fnl%blocks_y_p_gpu, start_gpu=self%pml_fnl%start_y_p_gpu,      &
                                         cells_gpu=self%pml_fnl%cells_y_p_gpu, gamma_gpu=self%pml_fnl%gamma_y_p_gpu,        &
                                         alpha_gpu=self%pml_fnl%alpha_y_p_gpu, kappa_gpu=self%pml_fnl%kappa_y_p_gpu)
      endif
      if (associated(self%pml_fnl%q_pml_z_m_gpu)) then
         call apply_z_face_field_correction_dev(q_face_gpu=self%pml_fnl%q_pml_z_m_gpu, blocks_gpu=self%pml_fnl%blocks_z_m_gpu, &
                                                start_gpu=self%pml_fnl%start_z_m_gpu, cells_gpu=self%pml_fnl%cells_z_m_gpu,     &
                                                kappa_gpu=self%pml_fnl%kappa_z_m_gpu)
         call compute_z_face_pml_rhs_dev(q_face_gpu=self%pml_fnl%q_pml_z_m_gpu, dq_face_gpu=self%rk_pml_fnl%dq_pml_z_m_gpu, &
                                         blocks_gpu=self%pml_fnl%blocks_z_m_gpu, start_gpu=self%pml_fnl%start_z_m_gpu,      &
                                         cells_gpu=self%pml_fnl%cells_z_m_gpu, gamma_gpu=self%pml_fnl%gamma_z_m_gpu,        &
                                         alpha_gpu=self%pml_fnl%alpha_z_m_gpu, kappa_gpu=self%pml_fnl%kappa_z_m_gpu)
      endif
      if (associated(self%pml_fnl%q_pml_z_p_gpu)) then
         call apply_z_face_field_correction_dev(q_face_gpu=self%pml_fnl%q_pml_z_p_gpu, blocks_gpu=self%pml_fnl%blocks_z_p_gpu, &
                                                start_gpu=self%pml_fnl%start_z_p_gpu, cells_gpu=self%pml_fnl%cells_z_p_gpu,     &
                                                kappa_gpu=self%pml_fnl%kappa_z_p_gpu)
         call compute_z_face_pml_rhs_dev(q_face_gpu=self%pml_fnl%q_pml_z_p_gpu, dq_face_gpu=self%rk_pml_fnl%dq_pml_z_p_gpu, &
                                         blocks_gpu=self%pml_fnl%blocks_z_p_gpu, start_gpu=self%pml_fnl%start_z_p_gpu,      &
                                         cells_gpu=self%pml_fnl%cells_z_p_gpu, gamma_gpu=self%pml_fnl%gamma_z_p_gpu,        &
                                         alpha_gpu=self%pml_fnl%alpha_z_p_gpu, kappa_gpu=self%pml_fnl%kappa_z_p_gpu)
      endif
   endif

   contains
      subroutine apply_x_face_field_correction_dev(q_face_gpu, blocks_gpu, start_gpu, cells_gpu, kappa_gpu)
      real(R8P),    intent(in)    :: q_face_gpu(:,:,:,:,:)
      integer(I4P), intent(in)    :: blocks_gpu(:), start_gpu(:), cells_gpu(:)
      real(R8P),    intent(in)    :: kappa_gpu(:,:)
      integer(I4P)                :: b, cells, i, i0, j, k, lid, li, ss
      real(R8P)                   :: d_field, dxyz
      real(R8P)                   :: kappa
      real(R8P)                   :: q_line(1-FDV_S_MAX:1+FDV_S_MAX)
      !$acc routine(compute_derivative1_fd_centered)
      !$omp declare target(compute_derivative1_fd_centered)
      !$acc parallel loop gang collapse(3) independent DEVICEVAR(q_gpu, dq_gpu, q_face_gpu, blocks_gpu, start_gpu, cells_gpu, kappa_gpu, dxyz_gpu) &
      !$acc& firstprivate(s1, inv_eps_scale, inv_mu_scale) private(b,cells,i,i0,li,d_field,dxyz,kappa,ss,q_line)
      !$omp OMPLOOP collapse(3) DEVICEPTR(q_gpu, dq_gpu, q_face_gpu, blocks_gpu, start_gpu, cells_gpu, kappa_gpu, dxyz_gpu) &
      !$omp& firstprivate(s1, inv_eps_scale, inv_mu_scale) private(b,cells,i,i0,li,d_field,dxyz,kappa,ss,q_line)
      do lid = 1, size(blocks_gpu)
         do k = 1, nk
            do j = 1, nj
               b = blocks_gpu(lid)
               i0 = start_gpu(lid)
               cells = cells_gpu(lid)
               dxyz = dxyz_gpu(b,1)
               !$acc loop vector
               do li = 1, cells
                  i = i0 + li - 1_I4P
                  kappa = kappa_gpu(lid,li)
                  !$acc loop seq
                  do ss = 1 - s1, 1 + s1
                     q_line(ss) = q_gpu(b,i+ss-1_I4P,j,k,VAR_BZ)
                  enddo
                  call compute_derivative1_fd_centered(s=s1, ds=dxyz, q=q_line(1-s1:1+s1), dq_ds=d_field)
                  dq_gpu(b,i,j,k,VAR_DY) = dq_gpu(b,i,j,k,VAR_DY) + (1._R8P / kappa - 1._R8P) * (-d_field * inv_mu_scale) + q_face_gpu(lid,li,j,k,4)
                  !$acc loop seq
                  do ss = 1 - s1, 1 + s1
                     q_line(ss) = q_gpu(b,i+ss-1_I4P,j,k,VAR_BY)
                  enddo
                  call compute_derivative1_fd_centered(s=s1, ds=dxyz, q=q_line(1-s1:1+s1), dq_ds=d_field)
                  dq_gpu(b,i,j,k,VAR_DZ) = dq_gpu(b,i,j,k,VAR_DZ) + (1._R8P / kappa - 1._R8P) * (d_field * inv_mu_scale) - q_face_gpu(lid,li,j,k,3)
                  !$acc loop seq
                  do ss = 1 - s1, 1 + s1
                     q_line(ss) = q_gpu(b,i+ss-1_I4P,j,k,VAR_DZ)
                  enddo
                  call compute_derivative1_fd_centered(s=s1, ds=dxyz, q=q_line(1-s1:1+s1), dq_ds=d_field)
                  dq_gpu(b,i,j,k,VAR_BY) = dq_gpu(b,i,j,k,VAR_BY) + (1._R8P / kappa - 1._R8P) * (d_field * inv_eps_scale) - q_face_gpu(lid,li,j,k,2)
                  !$acc loop seq
                  do ss = 1 - s1, 1 + s1
                     q_line(ss) = q_gpu(b,i+ss-1_I4P,j,k,VAR_DY)
                  enddo
                  call compute_derivative1_fd_centered(s=s1, ds=dxyz, q=q_line(1-s1:1+s1), dq_ds=d_field)
                  dq_gpu(b,i,j,k,VAR_BZ) = dq_gpu(b,i,j,k,VAR_BZ) + (1._R8P / kappa - 1._R8P) * (-d_field * inv_eps_scale) + q_face_gpu(lid,li,j,k,1)
               enddo
            enddo
         enddo
      enddo
      endsubroutine apply_x_face_field_correction_dev

      subroutine compute_x_face_pml_rhs_dev(q_face_gpu, dq_face_gpu, blocks_gpu, start_gpu, cells_gpu, gamma_gpu, alpha_gpu, kappa_gpu)
      real(R8P),    intent(in)    :: q_face_gpu(:,:,:,:,:)
      real(R8P),    intent(inout) :: dq_face_gpu(:,:,:,:,:)
      integer(I4P), intent(in)    :: blocks_gpu(:), start_gpu(:), cells_gpu(:)
      real(R8P),    intent(in)    :: gamma_gpu(:,:), alpha_gpu(:,:), kappa_gpu(:,:)
      integer(I4P)                :: b, cells, i, i0, j, k, lid, li, ss
      real(R8P)                   :: alpha, d_field, dxyz, gamma, kappa
      real(R8P)                   :: q_line(1-FDV_S_MAX:1+FDV_S_MAX)
      !$acc routine(compute_derivative1_fd_centered)
      !$omp declare target(compute_derivative1_fd_centered)
      !$acc parallel loop gang collapse(3) independent DEVICEVAR(q_gpu, q_face_gpu, dq_face_gpu, blocks_gpu, start_gpu, cells_gpu, gamma_gpu, alpha_gpu, kappa_gpu, dxyz_gpu) &
      !$acc& firstprivate(s1, inv_eps_scale, inv_mu_scale) private(b,cells,i,i0,li,d_field,dxyz,gamma,alpha,kappa,ss,q_line)
      !$omp OMPLOOP collapse(3) DEVICEPTR(q_gpu, q_face_gpu, dq_face_gpu, blocks_gpu, start_gpu, cells_gpu, gamma_gpu, alpha_gpu, kappa_gpu, dxyz_gpu) &
      !$omp& firstprivate(s1, inv_eps_scale, inv_mu_scale) private(b,cells,i,i0,li,d_field,dxyz,gamma,alpha,kappa,ss,q_line)
      do lid = 1, size(blocks_gpu)
         do k = 1, nk
            do j = 1, nj
               b = blocks_gpu(lid)
               i0 = start_gpu(lid)
               cells = cells_gpu(lid)
               dxyz = dxyz_gpu(b,1)
               !$acc loop vector
               do li = 1, cells
                  i = i0 + li - 1_I4P
                  gamma = gamma_gpu(lid,li)
                  alpha = alpha_gpu(lid,li)
                  kappa = kappa_gpu(lid,li)
                  !$acc loop seq
                  do ss = 1 - s1, 1 + s1
                     q_line(ss) = q_gpu(b,i+ss-1_I4P,j,k,VAR_DY)
                  enddo
                  call compute_derivative1_fd_centered(s=s1, ds=dxyz, q=q_line(1-s1:1+s1), dq_ds=d_field)
                  dq_face_gpu(lid,li,j,k,1) = gamma / kappa**2 * d_field * inv_eps_scale - (alpha + gamma / kappa) * q_face_gpu(lid,li,j,k,1)
                  !$acc loop seq
                  do ss = 1 - s1, 1 + s1
                     q_line(ss) = q_gpu(b,i+ss-1_I4P,j,k,VAR_DZ)
                  enddo
                  call compute_derivative1_fd_centered(s=s1, ds=dxyz, q=q_line(1-s1:1+s1), dq_ds=d_field)
                  dq_face_gpu(lid,li,j,k,2) = gamma / kappa**2 * d_field * inv_eps_scale - (alpha + gamma / kappa) * q_face_gpu(lid,li,j,k,2)
                  !$acc loop seq
                  do ss = 1 - s1, 1 + s1
                     q_line(ss) = q_gpu(b,i+ss-1_I4P,j,k,VAR_BY)
                  enddo
                  call compute_derivative1_fd_centered(s=s1, ds=dxyz, q=q_line(1-s1:1+s1), dq_ds=d_field)
                  dq_face_gpu(lid,li,j,k,3) = gamma / kappa**2 * d_field * inv_mu_scale - (alpha + gamma / kappa) * q_face_gpu(lid,li,j,k,3)
                  !$acc loop seq
                  do ss = 1 - s1, 1 + s1
                     q_line(ss) = q_gpu(b,i+ss-1_I4P,j,k,VAR_BZ)
                  enddo
                  call compute_derivative1_fd_centered(s=s1, ds=dxyz, q=q_line(1-s1:1+s1), dq_ds=d_field)
                  dq_face_gpu(lid,li,j,k,4) = gamma / kappa**2 * d_field * inv_mu_scale - (alpha + gamma / kappa) * q_face_gpu(lid,li,j,k,4)
               enddo
            enddo
         enddo
      enddo
      endsubroutine compute_x_face_pml_rhs_dev

      subroutine apply_y_face_field_correction_dev(q_face_gpu, blocks_gpu, start_gpu, cells_gpu, kappa_gpu)
      real(R8P),    intent(in)    :: q_face_gpu(:,:,:,:,:)
      integer(I4P), intent(in)    :: blocks_gpu(:), start_gpu(:), cells_gpu(:)
      real(R8P),    intent(in)    :: kappa_gpu(:,:)
      integer(I4P)                :: b, cells, i, j, j0, k, lid, lj, ss
      real(R8P)                   :: d_field, dxyz, kappa
      real(R8P)                   :: q_line(1-FDV_S_MAX:1+FDV_S_MAX)
      !$acc routine(compute_derivative1_fd_centered)
      !$omp declare target(compute_derivative1_fd_centered)
      !$acc parallel loop gang collapse(3) independent DEVICEVAR(q_gpu, dq_gpu, q_face_gpu, blocks_gpu, start_gpu, cells_gpu, kappa_gpu, dxyz_gpu) &
      !$acc& firstprivate(s1, inv_eps_scale, inv_mu_scale) private(b,cells,j,j0,lj,d_field,dxyz,kappa,ss,q_line)
      !$omp OMPLOOP collapse(3) DEVICEPTR(q_gpu, dq_gpu, q_face_gpu, blocks_gpu, start_gpu, cells_gpu, kappa_gpu, dxyz_gpu) &
      !$omp& firstprivate(s1, inv_eps_scale, inv_mu_scale) private(b,cells,j,j0,lj,d_field,dxyz,kappa,ss,q_line)
      do lid = 1, size(blocks_gpu)
         do k = 1, nk
            do i = 1, ni
               b = blocks_gpu(lid)
               j0 = start_gpu(lid)
               cells = cells_gpu(lid)
               dxyz = dxyz_gpu(b,2)
               !$acc loop vector
               do lj = 1, cells
                  j = j0 + lj - 1_I4P
                  kappa = kappa_gpu(lid,lj)
                  !$acc loop seq
                  do ss = 1 - s1, 1 + s1
                     q_line(ss) = q_gpu(b,i,j+ss-1_I4P,k,VAR_BZ)
                  enddo
                  call compute_derivative1_fd_centered(s=s1, ds=dxyz, q=q_line(1-s1:1+s1), dq_ds=d_field)
                  dq_gpu(b,i,j,k,VAR_DX) = dq_gpu(b,i,j,k,VAR_DX) + (1._R8P / kappa - 1._R8P) * (d_field * inv_mu_scale) - q_face_gpu(lid,i,lj,k,4)
                  !$acc loop seq
                  do ss = 1 - s1, 1 + s1
                     q_line(ss) = q_gpu(b,i,j+ss-1_I4P,k,VAR_BX)
                  enddo
                  call compute_derivative1_fd_centered(s=s1, ds=dxyz, q=q_line(1-s1:1+s1), dq_ds=d_field)
                  dq_gpu(b,i,j,k,VAR_DZ) = dq_gpu(b,i,j,k,VAR_DZ) + (1._R8P / kappa - 1._R8P) * (-d_field * inv_mu_scale) + q_face_gpu(lid,i,lj,k,3)
                  !$acc loop seq
                  do ss = 1 - s1, 1 + s1
                     q_line(ss) = q_gpu(b,i,j+ss-1_I4P,k,VAR_DZ)
                  enddo
                  call compute_derivative1_fd_centered(s=s1, ds=dxyz, q=q_line(1-s1:1+s1), dq_ds=d_field)
                  dq_gpu(b,i,j,k,VAR_BX) = dq_gpu(b,i,j,k,VAR_BX) + (1._R8P / kappa - 1._R8P) * (-d_field * inv_eps_scale) + q_face_gpu(lid,i,lj,k,2)
                  !$acc loop seq
                  do ss = 1 - s1, 1 + s1
                     q_line(ss) = q_gpu(b,i,j+ss-1_I4P,k,VAR_DX)
                  enddo
                  call compute_derivative1_fd_centered(s=s1, ds=dxyz, q=q_line(1-s1:1+s1), dq_ds=d_field)
                  dq_gpu(b,i,j,k,VAR_BZ) = dq_gpu(b,i,j,k,VAR_BZ) + (1._R8P / kappa - 1._R8P) * (d_field * inv_eps_scale) - q_face_gpu(lid,i,lj,k,1)
               enddo
            enddo
         enddo
      enddo
      endsubroutine apply_y_face_field_correction_dev

      subroutine compute_y_face_pml_rhs_dev(q_face_gpu, dq_face_gpu, blocks_gpu, start_gpu, cells_gpu, gamma_gpu, alpha_gpu, kappa_gpu)
      real(R8P),    intent(in)    :: q_face_gpu(:,:,:,:,:)
      real(R8P),    intent(inout) :: dq_face_gpu(:,:,:,:,:)
      integer(I4P), intent(in)    :: blocks_gpu(:), start_gpu(:), cells_gpu(:)
      real(R8P),    intent(in)    :: gamma_gpu(:,:), alpha_gpu(:,:), kappa_gpu(:,:)
      integer(I4P)                :: b, cells, i, j, j0, k, lid, lj, ss
      real(R8P)                   :: alpha, d_field, dxyz, gamma, kappa
      real(R8P)                   :: q_line(1-FDV_S_MAX:1+FDV_S_MAX)
      !$acc routine(compute_derivative1_fd_centered)
      !$omp declare target(compute_derivative1_fd_centered)
      !$acc parallel loop gang collapse(3) independent DEVICEVAR(q_gpu, q_face_gpu, dq_face_gpu, blocks_gpu, start_gpu, cells_gpu, gamma_gpu, alpha_gpu, kappa_gpu, dxyz_gpu) &
      !$acc& firstprivate(s1, inv_eps_scale, inv_mu_scale) private(b,cells,j,j0,lj,d_field,dxyz,gamma,alpha,kappa,ss,q_line)
      !$omp OMPLOOP collapse(3) DEVICEPTR(q_gpu, q_face_gpu, dq_face_gpu, blocks_gpu, start_gpu, cells_gpu, gamma_gpu, alpha_gpu, kappa_gpu, dxyz_gpu) &
      !$omp& firstprivate(s1, inv_eps_scale, inv_mu_scale) private(b,cells,j,j0,lj,d_field,dxyz,gamma,alpha,kappa,ss,q_line)
      do lid = 1, size(blocks_gpu)
         do k = 1, nk
            do i = 1, ni
               b = blocks_gpu(lid)
               j0 = start_gpu(lid)
               cells = cells_gpu(lid)
               dxyz = dxyz_gpu(b,2)
               !$acc loop vector
               do lj = 1, cells
                  j = j0 + lj - 1_I4P
                  gamma = gamma_gpu(lid,lj)
                  alpha = alpha_gpu(lid,lj)
                  kappa = kappa_gpu(lid,lj)
                  !$acc loop seq
                  do ss = 1 - s1, 1 + s1
                     q_line(ss) = q_gpu(b,i,j+ss-1_I4P,k,VAR_DX)
                  enddo
                  call compute_derivative1_fd_centered(s=s1, ds=dxyz, q=q_line(1-s1:1+s1), dq_ds=d_field)
                  dq_face_gpu(lid,i,lj,k,1) = gamma / kappa**2 * d_field * inv_eps_scale - (alpha + gamma / kappa) * q_face_gpu(lid,i,lj,k,1)
                  !$acc loop seq
                  do ss = 1 - s1, 1 + s1
                     q_line(ss) = q_gpu(b,i,j+ss-1_I4P,k,VAR_DZ)
                  enddo
                  call compute_derivative1_fd_centered(s=s1, ds=dxyz, q=q_line(1-s1:1+s1), dq_ds=d_field)
                  dq_face_gpu(lid,i,lj,k,2) = gamma / kappa**2 * d_field * inv_eps_scale - (alpha + gamma / kappa) * q_face_gpu(lid,i,lj,k,2)
                  !$acc loop seq
                  do ss = 1 - s1, 1 + s1
                     q_line(ss) = q_gpu(b,i,j+ss-1_I4P,k,VAR_BX)
                  enddo
                  call compute_derivative1_fd_centered(s=s1, ds=dxyz, q=q_line(1-s1:1+s1), dq_ds=d_field)
                  dq_face_gpu(lid,i,lj,k,3) = gamma / kappa**2 * d_field * inv_mu_scale - (alpha + gamma / kappa) * q_face_gpu(lid,i,lj,k,3)
                  !$acc loop seq
                  do ss = 1 - s1, 1 + s1
                     q_line(ss) = q_gpu(b,i,j+ss-1_I4P,k,VAR_BZ)
                  enddo
                  call compute_derivative1_fd_centered(s=s1, ds=dxyz, q=q_line(1-s1:1+s1), dq_ds=d_field)
                  dq_face_gpu(lid,i,lj,k,4) = gamma / kappa**2 * d_field * inv_mu_scale - (alpha + gamma / kappa) * q_face_gpu(lid,i,lj,k,4)
               enddo
            enddo
         enddo
      enddo
      endsubroutine compute_y_face_pml_rhs_dev

      subroutine apply_z_face_field_correction_dev(q_face_gpu, blocks_gpu, start_gpu, cells_gpu, kappa_gpu)
      real(R8P),    intent(in)    :: q_face_gpu(:,:,:,:,:)
      integer(I4P), intent(in)    :: blocks_gpu(:), start_gpu(:), cells_gpu(:)
      real(R8P),    intent(in)    :: kappa_gpu(:,:)
      integer(I4P)                :: b, cells, i, j, k, k0, lid, lk, ss
      real(R8P)                   :: d_field, dxyz, kappa
      real(R8P)                   :: q_line(1-FDV_S_MAX:1+FDV_S_MAX)
      !$acc routine(compute_derivative1_fd_centered)
      !$omp declare target(compute_derivative1_fd_centered)
      !$acc parallel loop gang collapse(3) independent DEVICEVAR(q_gpu, dq_gpu, q_face_gpu, blocks_gpu, start_gpu, cells_gpu, kappa_gpu, dxyz_gpu) &
      !$acc& firstprivate(s1, inv_eps_scale, inv_mu_scale) private(b,cells,k,k0,lk,d_field,dxyz,kappa,ss,q_line)
      !$omp OMPLOOP collapse(3) DEVICEPTR(q_gpu, dq_gpu, q_face_gpu, blocks_gpu, start_gpu, cells_gpu, kappa_gpu, dxyz_gpu) &
      !$omp& firstprivate(s1, inv_eps_scale, inv_mu_scale) private(b,cells,k,k0,lk,d_field,dxyz,kappa,ss,q_line)
      do lid = 1, size(blocks_gpu)
         do j = 1, nj
            do i = 1, ni
               b = blocks_gpu(lid)
               k0 = start_gpu(lid)
               cells = cells_gpu(lid)
               dxyz = dxyz_gpu(b,3)
               !$acc loop vector
               do lk = 1, cells
                  k = k0 + lk - 1_I4P
                  kappa = kappa_gpu(lid,lk)
                  !$acc loop seq
                  do ss = 1 - s1, 1 + s1
                     q_line(ss) = q_gpu(b,i,j,k+ss-1_I4P,VAR_BY)
                  enddo
                  call compute_derivative1_fd_centered(s=s1, ds=dxyz, q=q_line(1-s1:1+s1), dq_ds=d_field)
                  dq_gpu(b,i,j,k,VAR_DX) = dq_gpu(b,i,j,k,VAR_DX) + (1._R8P / kappa - 1._R8P) * (-d_field * inv_mu_scale) + q_face_gpu(lid,i,j,lk,4)
                  !$acc loop seq
                  do ss = 1 - s1, 1 + s1
                     q_line(ss) = q_gpu(b,i,j,k+ss-1_I4P,VAR_BX)
                  enddo
                  call compute_derivative1_fd_centered(s=s1, ds=dxyz, q=q_line(1-s1:1+s1), dq_ds=d_field)
                  dq_gpu(b,i,j,k,VAR_DY) = dq_gpu(b,i,j,k,VAR_DY) + (1._R8P / kappa - 1._R8P) * (d_field * inv_mu_scale) - q_face_gpu(lid,i,j,lk,3)
                  !$acc loop seq
                  do ss = 1 - s1, 1 + s1
                     q_line(ss) = q_gpu(b,i,j,k+ss-1_I4P,VAR_DY)
                  enddo
                  call compute_derivative1_fd_centered(s=s1, ds=dxyz, q=q_line(1-s1:1+s1), dq_ds=d_field)
                  dq_gpu(b,i,j,k,VAR_BX) = dq_gpu(b,i,j,k,VAR_BX) + (1._R8P / kappa - 1._R8P) * (d_field * inv_eps_scale) - q_face_gpu(lid,i,j,lk,2)
                  !$acc loop seq
                  do ss = 1 - s1, 1 + s1
                     q_line(ss) = q_gpu(b,i,j,k+ss-1_I4P,VAR_DX)
                  enddo
                  call compute_derivative1_fd_centered(s=s1, ds=dxyz, q=q_line(1-s1:1+s1), dq_ds=d_field)
                  dq_gpu(b,i,j,k,VAR_BY) = dq_gpu(b,i,j,k,VAR_BY) + (1._R8P / kappa - 1._R8P) * (-d_field * inv_eps_scale) + q_face_gpu(lid,i,j,lk,1)
               enddo
            enddo
         enddo
      enddo
      endsubroutine apply_z_face_field_correction_dev

      subroutine compute_z_face_pml_rhs_dev(q_face_gpu, dq_face_gpu, blocks_gpu, start_gpu, cells_gpu, gamma_gpu, alpha_gpu, kappa_gpu)
      real(R8P),    intent(in)    :: q_face_gpu(:,:,:,:,:)
      real(R8P),    intent(inout) :: dq_face_gpu(:,:,:,:,:)
      integer(I4P), intent(in)    :: blocks_gpu(:), start_gpu(:), cells_gpu(:)
      real(R8P),    intent(in)    :: gamma_gpu(:,:), alpha_gpu(:,:), kappa_gpu(:,:)
      integer(I4P)                :: b, cells, i, j, k, k0, lid, lk, ss
      real(R8P)                   :: alpha, d_field, dxyz, gamma, kappa
      real(R8P)                   :: q_line(1-FDV_S_MAX:1+FDV_S_MAX)
      !$acc routine(compute_derivative1_fd_centered)
      !$omp declare target(compute_derivative1_fd_centered)
      !$acc parallel loop gang collapse(3) independent DEVICEVAR(q_gpu, q_face_gpu, dq_face_gpu, blocks_gpu, start_gpu, cells_gpu, gamma_gpu, alpha_gpu, kappa_gpu, dxyz_gpu) &
      !$acc& firstprivate(s1, inv_eps_scale, inv_mu_scale) private(b,cells,k,k0,lk,d_field,dxyz,gamma,alpha,kappa,ss,q_line)
      !$omp OMPLOOP collapse(3) DEVICEPTR(q_gpu, q_face_gpu, dq_face_gpu, blocks_gpu, start_gpu, cells_gpu, gamma_gpu, alpha_gpu, kappa_gpu, dxyz_gpu) &
      !$omp& firstprivate(s1, inv_eps_scale, inv_mu_scale) private(b,cells,k,k0,lk,d_field,dxyz,gamma,alpha,kappa,ss,q_line)
      do lid = 1, size(blocks_gpu)
         do j = 1, nj
            do i = 1, ni
               b = blocks_gpu(lid)
               k0 = start_gpu(lid)
               cells = cells_gpu(lid)
               dxyz = dxyz_gpu(b,3)
               !$acc loop vector
               do lk = 1, cells
                  k = k0 + lk - 1_I4P
                  gamma = gamma_gpu(lid,lk)
                  alpha = alpha_gpu(lid,lk)
                  kappa = kappa_gpu(lid,lk)
                  !$acc loop seq
                  do ss = 1 - s1, 1 + s1
                     q_line(ss) = q_gpu(b,i,j,k+ss-1_I4P,VAR_DX)
                  enddo
                  call compute_derivative1_fd_centered(s=s1, ds=dxyz, q=q_line(1-s1:1+s1), dq_ds=d_field)
                  dq_face_gpu(lid,i,j,lk,1) = gamma / kappa**2 * d_field * inv_eps_scale - (alpha + gamma / kappa) * q_face_gpu(lid,i,j,lk,1)
                  !$acc loop seq
                  do ss = 1 - s1, 1 + s1
                     q_line(ss) = q_gpu(b,i,j,k+ss-1_I4P,VAR_DY)
                  enddo
                  call compute_derivative1_fd_centered(s=s1, ds=dxyz, q=q_line(1-s1:1+s1), dq_ds=d_field)
                  dq_face_gpu(lid,i,j,lk,2) = gamma / kappa**2 * d_field * inv_eps_scale - (alpha + gamma / kappa) * q_face_gpu(lid,i,j,lk,2)
                  !$acc loop seq
                  do ss = 1 - s1, 1 + s1
                     q_line(ss) = q_gpu(b,i,j,k+ss-1_I4P,VAR_BX)
                  enddo
                  call compute_derivative1_fd_centered(s=s1, ds=dxyz, q=q_line(1-s1:1+s1), dq_ds=d_field)
                  dq_face_gpu(lid,i,j,lk,3) = gamma / kappa**2 * d_field * inv_mu_scale - (alpha + gamma / kappa) * q_face_gpu(lid,i,j,lk,3)
                  !$acc loop seq
                  do ss = 1 - s1, 1 + s1
                     q_line(ss) = q_gpu(b,i,j,k+ss-1_I4P,VAR_BY)
                  enddo
                  call compute_derivative1_fd_centered(s=s1, ds=dxyz, q=q_line(1-s1:1+s1), dq_ds=d_field)
                  dq_face_gpu(lid,i,j,lk,4) = gamma / kappa**2 * d_field * inv_mu_scale - (alpha + gamma / kappa) * q_face_gpu(lid,i,j,lk,4)
               enddo
            enddo
         enddo
      enddo
      endsubroutine compute_z_face_pml_rhs_dev
   endsubroutine apply_pml_fd_centered_ade_dev

   subroutine fd_centered_plain_dev_kernel(ni, nj, nk, ngc, blocks_number, var_jx, var_jy, var_jz, s1, &
                                           inv_mu_scale, inv_eps_scale, dxyz_gpu, q_gpu, dq_gpu)
   !< Centered-FD Maxwell residual without hyperbolic cleaning.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc, blocks_number, var_jx, var_jy, var_jz, s1
   real(R8P),    intent(in)    :: inv_mu_scale, inv_eps_scale
   real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(inout) :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   integer(I4P)                :: i, j, k, b, s
   real(R8P)                   :: curlD(3), curlB(3), dxyz_b(3)
   real(R8P)                   :: qsx_y(1-FDV_S_MAX:1+FDV_S_MAX), qsx_z(1-FDV_S_MAX:1+FDV_S_MAX)
   real(R8P)                   :: qsy_x(1-FDV_S_MAX:1+FDV_S_MAX), qsy_z(1-FDV_S_MAX:1+FDV_S_MAX)
   real(R8P)                   :: qsz_x(1-FDV_S_MAX:1+FDV_S_MAX), qsz_y(1-FDV_S_MAX:1+FDV_S_MAX)

   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu,dq_gpu) &
   !$acc& firstprivate(var_jx,var_jy,var_jz,s1,inv_mu_scale,inv_eps_scale)                 &
   !$acc& private(curlD,curlB,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y,dxyz_b)
   !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,q_gpu,dq_gpu) &
   !$omp& firstprivate(var_jx,var_jy,var_jz,s1,inv_mu_scale,inv_eps_scale) &
   !$omp& private(curlD,curlB,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y,dxyz_b)
   do b=1,blocks_number
   do k=1,nk
   do j=1,nj
   do i=1,ni
      dxyz_b(1) = dxyz_gpu(b,1) ; dxyz_b(2) = dxyz_gpu(b,2) ; dxyz_b(3) = dxyz_gpu(b,3)
      !$acc loop seq
      do s=1-s1, 1+s1
         qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DY)
         qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DZ)
         qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DX)
         qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DZ)
         qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DX)
         qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DY)
      enddo
      call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_b,                                            &
                                        qsx_y=qsx_y(1-s1:1+s1),qsx_z=qsx_z(1-s1:1+s1),qsy_x=qsy_x(1-s1:1+s1), &
                                        qsy_z=qsy_z(1-s1:1+s1),qsz_x=qsz_x(1-s1:1+s1),qsz_y=qsz_y(1-s1:1+s1), &
                                        curl=curlD)
      !$acc loop seq
      do s=1-s1, 1+s1
         qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BY)
         qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BZ)
         qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BX)
         qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BZ)
         qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BX)
         qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BY)
      enddo
      call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_b,                                            &
                                        qsx_y=qsx_y(1-s1:1+s1),qsx_z=qsx_z(1-s1:1+s1),qsy_x=qsy_x(1-s1:1+s1), &
                                        qsy_z=qsy_z(1-s1:1+s1),qsz_x=qsz_x(1-s1:1+s1),qsz_y=qsz_y(1-s1:1+s1), &
                                        curl=curlB)
      dq_gpu(b,i,j,k,VAR_DX) =  curlB(1) * inv_mu_scale - q_gpu(b,i,j,k,var_jx)
      dq_gpu(b,i,j,k,VAR_DY) =  curlB(2) * inv_mu_scale - q_gpu(b,i,j,k,var_jy)
      dq_gpu(b,i,j,k,VAR_DZ) =  curlB(3) * inv_mu_scale - q_gpu(b,i,j,k,var_jz)
      dq_gpu(b,i,j,k,VAR_BX) = -curlD(1) * inv_eps_scale
      dq_gpu(b,i,j,k,VAR_BY) = -curlD(2) * inv_eps_scale
      dq_gpu(b,i,j,k,VAR_BZ) = -curlD(3) * inv_eps_scale
   enddo
   enddo
   enddo
   enddo
   endsubroutine fd_centered_plain_dev_kernel

   subroutine fd_centered_phi_dev_kernel(ni, nj, nk, ngc, blocks_number, var_jx, var_jy, var_jz, s1, &
                                         inv_mu_scale, inv_eps_scale, chi_wave, chi_damp, c_r, ivar_phi, &
                                         dxyz_gpu, q_gpu, dq_gpu)
   !< Centered-FD Maxwell residual with phi cleaning only.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc, blocks_number, var_jx, var_jy, var_jz, s1, ivar_phi
   real(R8P),    intent(in)    :: inv_mu_scale, inv_eps_scale, chi_wave, chi_damp, c_r
   real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(inout) :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   integer(I4P)                :: i, j, k, b, s
   real(R8P)                   :: curlD(3), curlB(3), gradphi(3), divergenceD, dxyz_b(3)
   real(R8P)                   :: damping_coeff, min_h
   real(R8P)                   :: qsx_y(1-FDV_S_MAX:1+FDV_S_MAX), qsx_z(1-FDV_S_MAX:1+FDV_S_MAX)
   real(R8P)                   :: qsy_x(1-FDV_S_MAX:1+FDV_S_MAX), qsy_z(1-FDV_S_MAX:1+FDV_S_MAX)
   real(R8P)                   :: qsz_x(1-FDV_S_MAX:1+FDV_S_MAX), qsz_y(1-FDV_S_MAX:1+FDV_S_MAX)
   real(R8P)                   :: qsx_x(1-FDV_S_MAX:1+FDV_S_MAX), qsy_y(1-FDV_S_MAX:1+FDV_S_MAX)
   real(R8P)                   :: qsz_z(1-FDV_S_MAX:1+FDV_S_MAX)

   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu,dq_gpu) &
   !$acc& firstprivate(var_jx,var_jy,var_jz,s1,inv_mu_scale,inv_eps_scale,chi_wave,chi_damp,c_r,ivar_phi) &
   !$acc& private(curlD,curlB,gradphi,divergenceD,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y,qsx_x,qsy_y,qsz_z,dxyz_b,min_h,damping_coeff)
   !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,q_gpu,dq_gpu) &
   !$omp& firstprivate(var_jx,var_jy,var_jz,s1,inv_mu_scale,inv_eps_scale,chi_wave,chi_damp,c_r,ivar_phi) &
   !$omp& private(curlD,curlB,gradphi,divergenceD,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y,qsx_x,qsy_y,qsz_z,dxyz_b,min_h,damping_coeff)
   do b=1,blocks_number
   do k=1,nk
   do j=1,nj
   do i=1,ni
      dxyz_b(1) = dxyz_gpu(b,1) ; dxyz_b(2) = dxyz_gpu(b,2) ; dxyz_b(3) = dxyz_gpu(b,3)
      !$acc loop seq
      do s=1-s1, 1+s1
         qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DY)
         qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DZ)
         qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DX)
         qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DZ)
         qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DX)
         qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DY)
      enddo
      call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_b,                                            &
                                        qsx_y=qsx_y(1-s1:1+s1),qsx_z=qsx_z(1-s1:1+s1),qsy_x=qsy_x(1-s1:1+s1), &
                                        qsy_z=qsy_z(1-s1:1+s1),qsz_x=qsz_x(1-s1:1+s1),qsz_y=qsz_y(1-s1:1+s1), &
                                        curl=curlD)
      !$acc loop seq
      do s=1-s1, 1+s1
         qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BY)
         qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BZ)
         qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BX)
         qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BZ)
         qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BX)
         qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BY)
      enddo
      call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_b,                                            &
                                        qsx_y=qsx_y(1-s1:1+s1),qsx_z=qsx_z(1-s1:1+s1),qsy_x=qsy_x(1-s1:1+s1), &
                                        qsy_z=qsy_z(1-s1:1+s1),qsz_x=qsz_x(1-s1:1+s1),qsz_y=qsz_y(1-s1:1+s1), &
                                        curl=curlB)
      !$acc loop seq
      do s=1-s1, 1+s1
         qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DX)
         qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DY)
         qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DZ)
      enddo
      call compute_divergence_fd_centered_dev(s=s1,dxyz=dxyz_b,                                      &
                                              qsx=qsx_x(1-s1:1+s1),qsy=qsy_y(1-s1:1+s1),qsz=qsz_z(1-s1:1+s1), &
                                              divergence=divergenceD)
      !$acc loop seq
      do s=1-s1, 1+s1
         qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,ivar_phi)
         qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,ivar_phi)
         qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,ivar_phi)
      enddo
      call compute_gradient_fd_centered_dev(s=s1,dxyz=dxyz_b,                                        &
                                            qsx=qsx_x(1-s1:1+s1),qsy=qsy_y(1-s1:1+s1),qsz=qsz_z(1-s1:1+s1), &
                                            gradient=gradphi)
      dq_gpu(b,i,j,k,VAR_DX) =  curlB(1) * inv_mu_scale - gradphi(1) - q_gpu(b,i,j,k,var_jx)
      dq_gpu(b,i,j,k,VAR_DY) =  curlB(2) * inv_mu_scale - gradphi(2) - q_gpu(b,i,j,k,var_jy)
      dq_gpu(b,i,j,k,VAR_DZ) =  curlB(3) * inv_mu_scale - gradphi(3) - q_gpu(b,i,j,k,var_jz)
      dq_gpu(b,i,j,k,VAR_BX) = -curlD(1) * inv_eps_scale
      dq_gpu(b,i,j,k,VAR_BY) = -curlD(2) * inv_eps_scale
      dq_gpu(b,i,j,k,VAR_BZ) = -curlD(3) * inv_eps_scale
      dq_gpu(b,i,j,k,ivar_phi) = -(chi_wave * chi_wave) * divergenceD
      if (c_r > 0._R8P) then
         min_h = min(dxyz_b(1), min(dxyz_b(2), dxyz_b(3)))
         damping_coeff = chi_damp / (c_r * min_h)
         dq_gpu(b,i,j,k,ivar_phi) = dq_gpu(b,i,j,k,ivar_phi) - damping_coeff * q_gpu(b,i,j,k,ivar_phi)
      endif
   enddo
   enddo
   enddo
   enddo
   endsubroutine fd_centered_phi_dev_kernel

   subroutine fd_centered_psi_dev_kernel(ni, nj, nk, ngc, blocks_number, var_jx, var_jy, var_jz, s1, &
                                         inv_mu_scale, inv_eps_scale, chi_wave, chi_damp, c_r, ivar_psi, &
                                         dxyz_gpu, q_gpu, dq_gpu)
   !< Centered-FD Maxwell residual with psi cleaning only.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc, blocks_number, var_jx, var_jy, var_jz, s1, ivar_psi
   real(R8P),    intent(in)    :: inv_mu_scale, inv_eps_scale, chi_wave, chi_damp, c_r
   real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(inout) :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   integer(I4P)                :: i, j, k, b, s
   real(R8P)                   :: curlD(3), curlB(3), gradpsi(3), divergenceB, dxyz_b(3)
   real(R8P)                   :: damping_coeff, min_h
   real(R8P)                   :: qsx_y(1-FDV_S_MAX:1+FDV_S_MAX), qsx_z(1-FDV_S_MAX:1+FDV_S_MAX)
   real(R8P)                   :: qsy_x(1-FDV_S_MAX:1+FDV_S_MAX), qsy_z(1-FDV_S_MAX:1+FDV_S_MAX)
   real(R8P)                   :: qsz_x(1-FDV_S_MAX:1+FDV_S_MAX), qsz_y(1-FDV_S_MAX:1+FDV_S_MAX)
   real(R8P)                   :: qsx_x(1-FDV_S_MAX:1+FDV_S_MAX), qsy_y(1-FDV_S_MAX:1+FDV_S_MAX)
   real(R8P)                   :: qsz_z(1-FDV_S_MAX:1+FDV_S_MAX)

   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu,dq_gpu) &
   !$acc& firstprivate(var_jx,var_jy,var_jz,s1,inv_mu_scale,inv_eps_scale,chi_wave,chi_damp,c_r,ivar_psi) &
   !$acc& private(curlD,curlB,gradpsi,divergenceB,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y,qsx_x,qsy_y,qsz_z,dxyz_b,min_h,damping_coeff)
   !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,q_gpu,dq_gpu) &
   !$omp& firstprivate(var_jx,var_jy,var_jz,s1,inv_mu_scale,inv_eps_scale,chi_wave,chi_damp,c_r,ivar_psi) &
   !$omp& private(curlD,curlB,gradpsi,divergenceB,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y,qsx_x,qsy_y,qsz_z,dxyz_b,min_h,damping_coeff)
   do b=1,blocks_number
   do k=1,nk
   do j=1,nj
   do i=1,ni
      dxyz_b(1) = dxyz_gpu(b,1) ; dxyz_b(2) = dxyz_gpu(b,2) ; dxyz_b(3) = dxyz_gpu(b,3)
      !$acc loop seq
      do s=1-s1, 1+s1
         qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DY)
         qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DZ)
         qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DX)
         qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DZ)
         qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DX)
         qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DY)
      enddo
      call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_b,                                            &
                                        qsx_y=qsx_y(1-s1:1+s1),qsx_z=qsx_z(1-s1:1+s1),qsy_x=qsy_x(1-s1:1+s1), &
                                        qsy_z=qsy_z(1-s1:1+s1),qsz_x=qsz_x(1-s1:1+s1),qsz_y=qsz_y(1-s1:1+s1), &
                                        curl=curlD)
      !$acc loop seq
      do s=1-s1, 1+s1
         qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BY)
         qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BZ)
         qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BX)
         qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BZ)
         qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BX)
         qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BY)
      enddo
      call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_b,                                            &
                                        qsx_y=qsx_y(1-s1:1+s1),qsx_z=qsx_z(1-s1:1+s1),qsy_x=qsy_x(1-s1:1+s1), &
                                        qsy_z=qsy_z(1-s1:1+s1),qsz_x=qsz_x(1-s1:1+s1),qsz_y=qsz_y(1-s1:1+s1), &
                                        curl=curlB)
      !$acc loop seq
      do s=1-s1, 1+s1
         qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BX)
         qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BY)
         qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BZ)
      enddo
      call compute_divergence_fd_centered_dev(s=s1,dxyz=dxyz_b,                                      &
                                              qsx=qsx_x(1-s1:1+s1),qsy=qsy_y(1-s1:1+s1),qsz=qsz_z(1-s1:1+s1), &
                                              divergence=divergenceB)
      !$acc loop seq
      do s=1-s1, 1+s1
         qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,ivar_psi)
         qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,ivar_psi)
         qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,ivar_psi)
      enddo
      call compute_gradient_fd_centered_dev(s=s1,dxyz=dxyz_b,                                        &
                                            qsx=qsx_x(1-s1:1+s1),qsy=qsy_y(1-s1:1+s1),qsz=qsz_z(1-s1:1+s1), &
                                            gradient=gradpsi)
      dq_gpu(b,i,j,k,VAR_DX) =  curlB(1) * inv_mu_scale - q_gpu(b,i,j,k,var_jx)
      dq_gpu(b,i,j,k,VAR_DY) =  curlB(2) * inv_mu_scale - q_gpu(b,i,j,k,var_jy)
      dq_gpu(b,i,j,k,VAR_DZ) =  curlB(3) * inv_mu_scale - q_gpu(b,i,j,k,var_jz)
      dq_gpu(b,i,j,k,VAR_BX) = -curlD(1) * inv_eps_scale - gradpsi(1)
      dq_gpu(b,i,j,k,VAR_BY) = -curlD(2) * inv_eps_scale - gradpsi(2)
      dq_gpu(b,i,j,k,VAR_BZ) = -curlD(3) * inv_eps_scale - gradpsi(3)
      dq_gpu(b,i,j,k,ivar_psi) = -(chi_wave * chi_wave) * divergenceB
      if (c_r > 0._R8P) then
         min_h = min(dxyz_b(1), min(dxyz_b(2), dxyz_b(3)))
         damping_coeff = chi_damp / (c_r * min_h)
         dq_gpu(b,i,j,k,ivar_psi) = dq_gpu(b,i,j,k,ivar_psi) - damping_coeff * q_gpu(b,i,j,k,ivar_psi)
      endif
   enddo
   enddo
   enddo
   enddo
   endsubroutine fd_centered_psi_dev_kernel

   subroutine fd_centered_phi_psi_dev_kernel(ni, nj, nk, ngc, blocks_number, var_jx, var_jy, var_jz, s1, &
                                             inv_mu_scale, inv_eps_scale, chi_wave, chi_damp, c_r, ivar_phi, ivar_psi, &
                                             dxyz_gpu, q_gpu, dq_gpu)
   !< Centered-FD Maxwell residual with phi/psi cleaning.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc, blocks_number, var_jx, var_jy, var_jz, s1, ivar_phi, ivar_psi
   real(R8P),    intent(in)    :: inv_mu_scale, inv_eps_scale, chi_wave, chi_damp, c_r
   real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(inout) :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   integer(I4P)                :: i, j, k, b, s
   real(R8P)                   :: curlD(3), curlB(3), gradphi(3), gradpsi(3), divergenceD, divergenceB, dxyz_b(3)
   real(R8P)                   :: damping_coeff, min_h
   real(R8P)                   :: qsx_y(1-FDV_S_MAX:1+FDV_S_MAX), qsx_z(1-FDV_S_MAX:1+FDV_S_MAX)
   real(R8P)                   :: qsy_x(1-FDV_S_MAX:1+FDV_S_MAX), qsy_z(1-FDV_S_MAX:1+FDV_S_MAX)
   real(R8P)                   :: qsz_x(1-FDV_S_MAX:1+FDV_S_MAX), qsz_y(1-FDV_S_MAX:1+FDV_S_MAX)
   real(R8P)                   :: qsx_x(1-FDV_S_MAX:1+FDV_S_MAX), qsy_y(1-FDV_S_MAX:1+FDV_S_MAX)
   real(R8P)                   :: qsz_z(1-FDV_S_MAX:1+FDV_S_MAX)

   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu,dq_gpu) &
   !$acc& firstprivate(var_jx,var_jy,var_jz,s1,inv_mu_scale,inv_eps_scale,chi_wave,chi_damp,c_r,ivar_phi,ivar_psi) &
   !$acc& private(curlD,curlB,gradphi,gradpsi,divergenceD,divergenceB,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y, &
   !$acc&         qsx_x,qsy_y,qsz_z,dxyz_b,min_h,damping_coeff)
   !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,q_gpu,dq_gpu) &
   !$omp& firstprivate(var_jx,var_jy,var_jz,s1,inv_mu_scale,inv_eps_scale,chi_wave,chi_damp,c_r,ivar_phi,ivar_psi) &
   !$omp& private(curlD,curlB,gradphi,gradpsi,divergenceD,divergenceB,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y, &
   !$omp&         qsx_x,qsy_y,qsz_z,dxyz_b,min_h,damping_coeff)
   do b=1,blocks_number
   do k=1,nk
   do j=1,nj
   do i=1,ni
      dxyz_b(1) = dxyz_gpu(b,1) ; dxyz_b(2) = dxyz_gpu(b,2) ; dxyz_b(3) = dxyz_gpu(b,3)
      !$acc loop seq
      do s=1-s1, 1+s1
         qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DY)
         qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DZ)
         qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DX)
         qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DZ)
         qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DX)
         qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DY)
      enddo
      call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_b,                                            &
                                        qsx_y=qsx_y(1-s1:1+s1),qsx_z=qsx_z(1-s1:1+s1),qsy_x=qsy_x(1-s1:1+s1), &
                                        qsy_z=qsy_z(1-s1:1+s1),qsz_x=qsz_x(1-s1:1+s1),qsz_y=qsz_y(1-s1:1+s1), &
                                        curl=curlD)
      !$acc loop seq
      do s=1-s1, 1+s1
         qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BY)
         qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BZ)
         qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BX)
         qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BZ)
         qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BX)
         qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BY)
      enddo
      call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_b,                                            &
                                        qsx_y=qsx_y(1-s1:1+s1),qsx_z=qsx_z(1-s1:1+s1),qsy_x=qsy_x(1-s1:1+s1), &
                                        qsy_z=qsy_z(1-s1:1+s1),qsz_x=qsz_x(1-s1:1+s1),qsz_y=qsz_y(1-s1:1+s1), &
                                        curl=curlB)
      !$acc loop seq
      do s=1-s1, 1+s1
         qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DX)
         qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DY)
         qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DZ)
      enddo
      call compute_divergence_fd_centered_dev(s=s1,dxyz=dxyz_b,                                      &
                                              qsx=qsx_x(1-s1:1+s1),qsy=qsy_y(1-s1:1+s1),qsz=qsz_z(1-s1:1+s1), &
                                              divergence=divergenceD)
      !$acc loop seq
      do s=1-s1, 1+s1
         qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,ivar_phi)
         qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,ivar_phi)
         qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,ivar_phi)
      enddo
      call compute_gradient_fd_centered_dev(s=s1,dxyz=dxyz_b,                                        &
                                            qsx=qsx_x(1-s1:1+s1),qsy=qsy_y(1-s1:1+s1),qsz=qsz_z(1-s1:1+s1), &
                                            gradient=gradphi)
      !$acc loop seq
      do s=1-s1, 1+s1
         qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BX)
         qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BY)
         qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BZ)
      enddo
      call compute_divergence_fd_centered_dev(s=s1,dxyz=dxyz_b,                                      &
                                              qsx=qsx_x(1-s1:1+s1),qsy=qsy_y(1-s1:1+s1),qsz=qsz_z(1-s1:1+s1), &
                                              divergence=divergenceB)
      !$acc loop seq
      do s=1-s1, 1+s1
         qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,ivar_psi)
         qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,ivar_psi)
         qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,ivar_psi)
      enddo
      call compute_gradient_fd_centered_dev(s=s1,dxyz=dxyz_b,                                        &
                                            qsx=qsx_x(1-s1:1+s1),qsy=qsy_y(1-s1:1+s1),qsz=qsz_z(1-s1:1+s1), &
                                            gradient=gradpsi)
      dq_gpu(b,i,j,k,VAR_DX) =  curlB(1) * inv_mu_scale - gradphi(1) - q_gpu(b,i,j,k,var_jx)
      dq_gpu(b,i,j,k,VAR_DY) =  curlB(2) * inv_mu_scale - gradphi(2) - q_gpu(b,i,j,k,var_jy)
      dq_gpu(b,i,j,k,VAR_DZ) =  curlB(3) * inv_mu_scale - gradphi(3) - q_gpu(b,i,j,k,var_jz)
      dq_gpu(b,i,j,k,VAR_BX) = -curlD(1) * inv_eps_scale - gradpsi(1)
      dq_gpu(b,i,j,k,VAR_BY) = -curlD(2) * inv_eps_scale - gradpsi(2)
      dq_gpu(b,i,j,k,VAR_BZ) = -curlD(3) * inv_eps_scale - gradpsi(3)
      dq_gpu(b,i,j,k,ivar_phi) = -(chi_wave * chi_wave) * divergenceD
      dq_gpu(b,i,j,k,ivar_psi) = -(chi_wave * chi_wave) * divergenceB
      if (c_r > 0._R8P) then
         min_h = min(dxyz_b(1), min(dxyz_b(2), dxyz_b(3)))
         damping_coeff = chi_damp / (c_r * min_h)
         dq_gpu(b,i,j,k,ivar_phi) = dq_gpu(b,i,j,k,ivar_phi) - damping_coeff * q_gpu(b,i,j,k,ivar_phi)
         dq_gpu(b,i,j,k,ivar_psi) = dq_gpu(b,i,j,k,ivar_psi) - damping_coeff * q_gpu(b,i,j,k,ivar_psi)
      endif
   enddo
   enddo
   enddo
   enddo
   endsubroutine fd_centered_phi_psi_dev_kernel

   subroutine compute_residuals_fv_centered_dev(self, q_gpu, dq_gpu, s, flux_register)
   !< Compute residuals, space operator, centered finite volume scheme — FNL device
   !< backend (issue #23 R2). 1:1 structural mirror of
   !< `prism_cpu_object%compute_residuals_fv_centered`: ghost refresh on q,
   !< pointwise Maxwell fluxes at ALL cells (incl. ghosts) into `flxyz_c_gpu`,
   !< three staggered face-reconstruction sweeps into `fl{x,y,z}_f_gpu` (the m=0 SOTA
   !< primitive, already `acc routine seq`), flux difference + J source into `dq_gpu`.
   !< The Maxwell-flux variant is selected ONCE at host side, then a specialized
   !< device kernel runs with no per-cell variant branch. Inter-realm seam accumulation
   !< (issue #23 R3): after the
   !< face sweeps, at the realm's FINAL RK substage only (α.r1, CPU parity), the seam
   !< face skins are device-packed, D2H-copied (tiny slabs) and accumulated into the
   !< HOST-side flux register — see accumulate_seam_fluxes_fv_dev.
   !<
   !< Race discipline (CLAUDE-gpu): the kernels below are MODULE-LEVEL procedures
   !< (not contained here) with constant-bound private gathers — the certified
   !< #22 F2/F3 pattern.
   class(prism_fnl_object), intent(inout)           :: self       !< The equation.
   real(R8P),               intent(inout)           :: q_gpu(1:,         &
                                                            1-self%ngc:,&
                                                            1-self%ngc:,&
                                                            1-self%ngc:,&
                                                            1:)  !< Conservative variables.
   real(R8P),               intent(inout)           :: dq_gpu(1:,         &
                                                             1-self%ngc:,&
                                                             1-self%ngc:,&
                                                             1-self%ngc:,&
                                                             1:) !< Residuals.
   integer(I4P),            intent(in),    optional :: s          !< Stage counter (gates the seam-flux accumulation).
   class(flux_register_object), intent(inout), optional :: flux_register !< Forest's flux register for FV seam reflux.
   integer(I4P)                                     :: stage_idx  !< Captured stage index (CPU-parity gate).
   real(R8P)                                        :: chi_damp   !< Host-side damping-speed factor.

   stage_idx = 0_I4P ; if (present(s)) stage_idx = s
   if (self%blocks_number > 0) then
      if (present(s)) then
         call self%update_ghost(q_gpu=q_gpu, s=s)
      else
         call self%update_ghost(q_gpu=q_gpu)
      endif
      select case(self%fv_flux_variant)
      case(FV_FLUX_VARIANT_MAXWELL)
         call fv_cell_fluxes_maxwell_dev_kernel(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, &
                                                nv_c=self%nv_c, blocks_number=self%blocks_number, &
                                                chi=self%physics%chi, q_gpu=q_gpu, flxyz_c_gpu=self%flxyz_c_gpu)
      case(FV_FLUX_VARIANT_MAXWELL_ADIM)
         call fv_cell_fluxes_maxwell_adim_dev_kernel(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, &
                                                     nv_c=self%nv_c, blocks_number=self%blocks_number, &
                                                     chi=self%physics%chi, q_gpu=q_gpu, flxyz_c_gpu=self%flxyz_c_gpu)
      case(FV_FLUX_VARIANT_MAXWELL_DIV_D)
         call fv_cell_fluxes_maxwell_div_d_dev_kernel(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, &
                                                      nv_c=self%nv_c, blocks_number=self%blocks_number, &
                                                      chi=self%physics%chi, q_gpu=q_gpu, flxyz_c_gpu=self%flxyz_c_gpu)
      case(FV_FLUX_VARIANT_MAXWELL_DIV_B)
         call fv_cell_fluxes_maxwell_div_b_dev_kernel(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, &
                                                      nv_c=self%nv_c, blocks_number=self%blocks_number, &
                                                      chi=self%physics%chi, q_gpu=q_gpu, flxyz_c_gpu=self%flxyz_c_gpu)
      case(FV_FLUX_VARIANT_MAXWELL_DIV_D_B)
         call fv_cell_fluxes_maxwell_div_d_b_dev_kernel(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, &
                                                        nv_c=self%nv_c, blocks_number=self%blocks_number, &
                                                        chi=self%physics%chi, q_gpu=q_gpu, flxyz_c_gpu=self%flxyz_c_gpu)
      case(FV_FLUX_VARIANT_MAXWELL_ADIM_DIV_D)
         call fv_cell_fluxes_maxwell_adim_div_d_dev_kernel(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, &
                                                           nv_c=self%nv_c, blocks_number=self%blocks_number, &
                                                           chi=self%physics%chi, q_gpu=q_gpu, flxyz_c_gpu=self%flxyz_c_gpu)
      case(FV_FLUX_VARIANT_MAXWELL_ADIM_DIV_B)
         call fv_cell_fluxes_maxwell_adim_div_b_dev_kernel(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, &
                                                           nv_c=self%nv_c, blocks_number=self%blocks_number, &
                                                           chi=self%physics%chi, q_gpu=q_gpu, flxyz_c_gpu=self%flxyz_c_gpu)
      case(FV_FLUX_VARIANT_MAXWELL_ADIM_DIV_D_B)
         call fv_cell_fluxes_maxwell_adim_div_d_b_dev_kernel(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, &
                                                             nv_c=self%nv_c, blocks_number=self%blocks_number, &
                                                             chi=self%physics%chi, q_gpu=q_gpu, flxyz_c_gpu=self%flxyz_c_gpu)
      case default
         call mpih_fnl%error_stop(msg=': unknown FV flux variant in FNL backend')
      end select
      call fv_recon_x_dev_kernel(s1=self%fdv_half_stencils(1), ni=self%ni, nj=self%nj,        &
                                 nk=self%nk, ngc=self%ngc, nv_c=self%nv_c,                    &
                                 blocks_number=self%blocks_number,                            &
                                 flxyz_c_gpu=self%flxyz_c_gpu, flx_f_gpu=self%flx_f_gpu)
      call fv_recon_y_dev_kernel(s1=self%fdv_half_stencils(1), ni=self%ni, nj=self%nj,        &
                                 nk=self%nk, ngc=self%ngc, nv_c=self%nv_c,                    &
                                 blocks_number=self%blocks_number,                            &
                                 flxyz_c_gpu=self%flxyz_c_gpu, fly_f_gpu=self%fly_f_gpu)
      call fv_recon_z_dev_kernel(s1=self%fdv_half_stencils(1), ni=self%ni, nj=self%nj,        &
                                 nk=self%nk, ngc=self%ngc, nv_c=self%nv_c,                    &
                                 blocks_number=self%blocks_number,                            &
                                 flxyz_c_gpu=self%flxyz_c_gpu, flz_f_gpu=self%flz_f_gpu)
      ! Inter-realm seam flux accumulation (issue #23 R3): α.r1 end-of-step gate, CPU
      ! parity with compute_residuals_fv_centered's hook — accumulate ONLY at the
      ! realm's final RK substage, from the just-reconstructed face fluxes, BEFORE
      ! the conservative update consumes them.
      if (present(flux_register) .and. stage_idx == self%rk%nrk &
                                 .and. allocated(self%adam%maps%inter_realm_face_register_index)) then
         if (flux_register%nfaces > 0_I4P) then
            call accumulate_seam_fluxes_fv_dev(self, flux_register)
         endif
      endif
      chi_damp = self%physics%chi * C0
      if (self%physics%physical_model == ADIM_EM_PHYSICAL_MODEL) chi_damp = self%physics%chi
      if (self%fv_add_phi_damping .and. self%fv_add_psi_damping) then
         call fv_flux_diff_phi_psi_dev_kernel(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, &
                                              nv_c=self%nv_c, blocks_number=self%blocks_number, &
                                              var_jx=self%physics%var_jx, var_jy=self%physics%var_jy, &
                                              var_jz=self%physics%var_jz, chi_damp=chi_damp, c_r=self%physics%c_r, &
                                              fv_ivar_phi=self%fv_ivar_phi, fv_ivar_psi=self%fv_ivar_psi, &
                                              dxyz_gpu=self%field_fnl%dxyz_gpu, flx_f_gpu=self%flx_f_gpu, &
                                              fly_f_gpu=self%fly_f_gpu, flz_f_gpu=self%flz_f_gpu, q_gpu=q_gpu, dq_gpu=dq_gpu)
      elseif (self%fv_add_phi_damping) then
         call fv_flux_diff_phi_dev_kernel(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, &
                                          nv_c=self%nv_c, blocks_number=self%blocks_number, &
                                          var_jx=self%physics%var_jx, var_jy=self%physics%var_jy, &
                                          var_jz=self%physics%var_jz, chi_damp=chi_damp, c_r=self%physics%c_r, &
                                          fv_ivar_phi=self%fv_ivar_phi, dxyz_gpu=self%field_fnl%dxyz_gpu, &
                                          flx_f_gpu=self%flx_f_gpu, fly_f_gpu=self%fly_f_gpu, flz_f_gpu=self%flz_f_gpu, &
                                          q_gpu=q_gpu, dq_gpu=dq_gpu)
      elseif (self%fv_add_psi_damping) then
         call fv_flux_diff_psi_dev_kernel(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, &
                                          nv_c=self%nv_c, blocks_number=self%blocks_number, &
                                          var_jx=self%physics%var_jx, var_jy=self%physics%var_jy, &
                                          var_jz=self%physics%var_jz, chi_damp=chi_damp, c_r=self%physics%c_r, &
                                          fv_ivar_psi=self%fv_ivar_psi, dxyz_gpu=self%field_fnl%dxyz_gpu, &
                                          flx_f_gpu=self%flx_f_gpu, fly_f_gpu=self%fly_f_gpu, flz_f_gpu=self%flz_f_gpu, &
                                          q_gpu=q_gpu, dq_gpu=dq_gpu)
      else
         call fv_flux_diff_plain_dev_kernel(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, &
                                            nv_c=self%nv_c, blocks_number=self%blocks_number, &
                                            var_jx=self%physics%var_jx, var_jy=self%physics%var_jy, &
                                            var_jz=self%physics%var_jz, dxyz_gpu=self%field_fnl%dxyz_gpu, &
                                            flx_f_gpu=self%flx_f_gpu, fly_f_gpu=self%fly_f_gpu, flz_f_gpu=self%flz_f_gpu, &
                                            q_gpu=q_gpu, dq_gpu=dq_gpu)
      endif
   endif
   endsubroutine compute_residuals_fv_centered_dev

   subroutine fv_cell_fluxes_maxwell_dev_kernel(ni, nj, nk, ngc, nv_c, blocks_number, chi, q_gpu, flxyz_c_gpu)
   !< Pointwise plain dimensional Maxwell fluxes at every cell (incl. ghosts).
   integer(I4P), intent(in)    :: ni, nj, nk, ngc, blocks_number
   integer(I4P), intent(in)    :: nv_c
   real(R8P),    intent(in)    :: chi
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(inout) :: flxyz_c_gpu(1:,1:,1:,1-ngc:,1-ngc:,1-ngc:,1:)
   integer(I4P), parameter     :: NVC_CAP=8_I4P
   real(R8P)                   :: qv(1:NVC_CAP), fv(1:NVC_CAP), sir(1:3)
   integer(I4P)                :: i, j, k, b, d, v

   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu, flxyz_c_gpu) &
   !$acc& firstprivate(ni, nj, nk, ngc, nv_c, blocks_number, chi)                        &
   !$acc& private(qv, fv, sir)
   !$omp OMPLOOP collapse(4) DEVICEPTR(q_gpu, flxyz_c_gpu) &
   !$omp& firstprivate(ni, nj, nk, ngc, nv_c, blocks_number, chi) &
   !$omp& private(qv, fv, sir)
   do b=1, blocks_number
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
      !$acc loop seq
      do v=1, nv_c
         qv(v) = q_gpu(b,i,j,k,v)
      enddo
      !$acc loop seq
      do d=1, 3
         sir(1) = 0._R8P ; sir(2) = 0._R8P ; sir(3) = 0._R8P ; sir(d) = 1._R8P
         call compute_convective_fluxes_maxwell(sir=sir, q=qv(1:nv_c), f=fv(1:nv_c), chi=chi)
         do v=1, nv_c
            flxyz_c_gpu(b,1,d,i,j,k,v) = fv(v)
         enddo
      enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine fv_cell_fluxes_maxwell_dev_kernel

   subroutine fv_cell_fluxes_maxwell_adim_dev_kernel(ni, nj, nk, ngc, nv_c, blocks_number, chi, q_gpu, flxyz_c_gpu)
   !< Pointwise plain adimensional Maxwell fluxes at every cell (incl. ghosts).
   integer(I4P), intent(in)    :: ni, nj, nk, ngc, blocks_number
   integer(I4P), intent(in)    :: nv_c
   real(R8P),    intent(in)    :: chi
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(inout) :: flxyz_c_gpu(1:,1:,1:,1-ngc:,1-ngc:,1-ngc:,1:)
   integer(I4P), parameter     :: NVC_CAP=8_I4P
   real(R8P)                   :: qv(1:NVC_CAP), fv(1:NVC_CAP), sir(1:3)
   integer(I4P)                :: i, j, k, b, d, v

   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu, flxyz_c_gpu) &
   !$acc& firstprivate(ni, nj, nk, ngc, nv_c, blocks_number, chi)                        &
   !$acc& private(qv, fv, sir)
   !$omp OMPLOOP collapse(4) DEVICEPTR(q_gpu, flxyz_c_gpu) &
   !$omp& firstprivate(ni, nj, nk, ngc, nv_c, blocks_number, chi) &
   !$omp& private(qv, fv, sir)
   do b=1, blocks_number
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
      !$acc loop seq
      do v=1, nv_c
         qv(v) = q_gpu(b,i,j,k,v)
      enddo
      !$acc loop seq
      do d=1, 3
         sir(1) = 0._R8P ; sir(2) = 0._R8P ; sir(3) = 0._R8P ; sir(d) = 1._R8P
         call compute_convective_fluxes_maxwell_adim(sir=sir, q=qv(1:nv_c), f=fv(1:nv_c), chi=chi)
         do v=1, nv_c
            flxyz_c_gpu(b,1,d,i,j,k,v) = fv(v)
         enddo
      enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine fv_cell_fluxes_maxwell_adim_dev_kernel

   subroutine fv_cell_fluxes_maxwell_div_d_dev_kernel(ni, nj, nk, ngc, nv_c, blocks_number, chi, q_gpu, flxyz_c_gpu)
   !< Pointwise dimensional Maxwell fluxes with D hyperbolic cleaning.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc, blocks_number
   integer(I4P), intent(in)    :: nv_c
   real(R8P),    intent(in)    :: chi
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(inout) :: flxyz_c_gpu(1:,1:,1:,1-ngc:,1-ngc:,1-ngc:,1:)
   integer(I4P), parameter     :: NVC_CAP=8_I4P
   real(R8P)                   :: qv(1:NVC_CAP), fv(1:NVC_CAP), sir(1:3)
   integer(I4P)                :: i, j, k, b, d, v

   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu, flxyz_c_gpu) &
   !$acc& firstprivate(ni, nj, nk, ngc, nv_c, blocks_number, chi)                        &
   !$acc& private(qv, fv, sir)
   !$omp OMPLOOP collapse(4) DEVICEPTR(q_gpu, flxyz_c_gpu) &
   !$omp& firstprivate(ni, nj, nk, ngc, nv_c, blocks_number, chi) &
   !$omp& private(qv, fv, sir)
   do b=1, blocks_number
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
      !$acc loop seq
      do v=1, nv_c
         qv(v) = q_gpu(b,i,j,k,v)
      enddo
      !$acc loop seq
      do d=1, 3
         sir(1) = 0._R8P ; sir(2) = 0._R8P ; sir(3) = 0._R8P ; sir(d) = 1._R8P
         call compute_convective_fluxes_maxwell_div_d(sir=sir, q=qv(1:nv_c), f=fv(1:nv_c), chi=chi)
         do v=1, nv_c
            flxyz_c_gpu(b,1,d,i,j,k,v) = fv(v)
         enddo
      enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine fv_cell_fluxes_maxwell_div_d_dev_kernel

   subroutine fv_cell_fluxes_maxwell_div_b_dev_kernel(ni, nj, nk, ngc, nv_c, blocks_number, chi, q_gpu, flxyz_c_gpu)
   !< Pointwise dimensional Maxwell fluxes with B hyperbolic cleaning.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc, blocks_number
   integer(I4P), intent(in)    :: nv_c
   real(R8P),    intent(in)    :: chi
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(inout) :: flxyz_c_gpu(1:,1:,1:,1-ngc:,1-ngc:,1-ngc:,1:)
   integer(I4P), parameter     :: NVC_CAP=8_I4P
   real(R8P)                   :: qv(1:NVC_CAP), fv(1:NVC_CAP), sir(1:3)
   integer(I4P)                :: i, j, k, b, d, v

   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu, flxyz_c_gpu) &
   !$acc& firstprivate(ni, nj, nk, ngc, nv_c, blocks_number, chi)                        &
   !$acc& private(qv, fv, sir)
   !$omp OMPLOOP collapse(4) DEVICEPTR(q_gpu, flxyz_c_gpu) &
   !$omp& firstprivate(ni, nj, nk, ngc, nv_c, blocks_number, chi) &
   !$omp& private(qv, fv, sir)
   do b=1, blocks_number
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
      !$acc loop seq
      do v=1, nv_c
         qv(v) = q_gpu(b,i,j,k,v)
      enddo
      !$acc loop seq
      do d=1, 3
         sir(1) = 0._R8P ; sir(2) = 0._R8P ; sir(3) = 0._R8P ; sir(d) = 1._R8P
         call compute_convective_fluxes_maxwell_div_b(sir=sir, q=qv(1:nv_c), f=fv(1:nv_c), chi=chi)
         do v=1, nv_c
            flxyz_c_gpu(b,1,d,i,j,k,v) = fv(v)
         enddo
      enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine fv_cell_fluxes_maxwell_div_b_dev_kernel

   subroutine fv_cell_fluxes_maxwell_div_d_b_dev_kernel(ni, nj, nk, ngc, nv_c, blocks_number, chi, q_gpu, flxyz_c_gpu)
   !< Pointwise dimensional Maxwell fluxes with D/B hyperbolic cleaning.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc, blocks_number
   integer(I4P), intent(in)    :: nv_c
   real(R8P),    intent(in)    :: chi
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(inout) :: flxyz_c_gpu(1:,1:,1:,1-ngc:,1-ngc:,1-ngc:,1:)
   integer(I4P), parameter     :: NVC_CAP=8_I4P
   real(R8P)                   :: qv(1:NVC_CAP), fv(1:NVC_CAP), sir(1:3)
   integer(I4P)                :: i, j, k, b, d, v

   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu, flxyz_c_gpu) &
   !$acc& firstprivate(ni, nj, nk, ngc, nv_c, blocks_number, chi)                        &
   !$acc& private(qv, fv, sir)
   !$omp OMPLOOP collapse(4) DEVICEPTR(q_gpu, flxyz_c_gpu) &
   !$omp& firstprivate(ni, nj, nk, ngc, nv_c, blocks_number, chi) &
   !$omp& private(qv, fv, sir)
   do b=1, blocks_number
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
      !$acc loop seq
      do v=1, nv_c
         qv(v) = q_gpu(b,i,j,k,v)
      enddo
      !$acc loop seq
      do d=1, 3
         sir(1) = 0._R8P ; sir(2) = 0._R8P ; sir(3) = 0._R8P ; sir(d) = 1._R8P
         call compute_convective_fluxes_maxwell_div_d_b(sir=sir, q=qv(1:nv_c), f=fv(1:nv_c), chi=chi)
         do v=1, nv_c
            flxyz_c_gpu(b,1,d,i,j,k,v) = fv(v)
         enddo
      enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine fv_cell_fluxes_maxwell_div_d_b_dev_kernel

   subroutine fv_cell_fluxes_maxwell_adim_div_d_dev_kernel(ni, nj, nk, ngc, nv_c, blocks_number, chi, q_gpu, flxyz_c_gpu)
   !< Pointwise adimensional Maxwell fluxes with D hyperbolic cleaning.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc, blocks_number
   integer(I4P), intent(in)    :: nv_c
   real(R8P),    intent(in)    :: chi
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(inout) :: flxyz_c_gpu(1:,1:,1:,1-ngc:,1-ngc:,1-ngc:,1:)
   integer(I4P), parameter     :: NVC_CAP=8_I4P
   real(R8P)                   :: qv(1:NVC_CAP), fv(1:NVC_CAP), sir(1:3)
   integer(I4P)                :: i, j, k, b, d, v

   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu, flxyz_c_gpu) &
   !$acc& firstprivate(ni, nj, nk, ngc, nv_c, blocks_number, chi)                        &
   !$acc& private(qv, fv, sir)
   !$omp OMPLOOP collapse(4) DEVICEPTR(q_gpu, flxyz_c_gpu) &
   !$omp& firstprivate(ni, nj, nk, ngc, nv_c, blocks_number, chi) &
   !$omp& private(qv, fv, sir)
   do b=1, blocks_number
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
      !$acc loop seq
      do v=1, nv_c
         qv(v) = q_gpu(b,i,j,k,v)
      enddo
      !$acc loop seq
      do d=1, 3
         sir(1) = 0._R8P ; sir(2) = 0._R8P ; sir(3) = 0._R8P ; sir(d) = 1._R8P
         call compute_convective_fluxes_maxwell_adim_div_d(sir=sir, q=qv(1:nv_c), f=fv(1:nv_c), chi=chi)
         do v=1, nv_c
            flxyz_c_gpu(b,1,d,i,j,k,v) = fv(v)
         enddo
      enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine fv_cell_fluxes_maxwell_adim_div_d_dev_kernel

   subroutine fv_cell_fluxes_maxwell_adim_div_b_dev_kernel(ni, nj, nk, ngc, nv_c, blocks_number, chi, q_gpu, flxyz_c_gpu)
   !< Pointwise adimensional Maxwell fluxes with B hyperbolic cleaning.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc, blocks_number
   integer(I4P), intent(in)    :: nv_c
   real(R8P),    intent(in)    :: chi
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(inout) :: flxyz_c_gpu(1:,1:,1:,1-ngc:,1-ngc:,1-ngc:,1:)
   integer(I4P), parameter     :: NVC_CAP=8_I4P
   real(R8P)                   :: qv(1:NVC_CAP), fv(1:NVC_CAP), sir(1:3)
   integer(I4P)                :: i, j, k, b, d, v

   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu, flxyz_c_gpu) &
   !$acc& firstprivate(ni, nj, nk, ngc, nv_c, blocks_number, chi)                        &
   !$acc& private(qv, fv, sir)
   !$omp OMPLOOP collapse(4) DEVICEPTR(q_gpu, flxyz_c_gpu) &
   !$omp& firstprivate(ni, nj, nk, ngc, nv_c, blocks_number, chi) &
   !$omp& private(qv, fv, sir)
   do b=1, blocks_number
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
      !$acc loop seq
      do v=1, nv_c
         qv(v) = q_gpu(b,i,j,k,v)
      enddo
      !$acc loop seq
      do d=1, 3
         sir(1) = 0._R8P ; sir(2) = 0._R8P ; sir(3) = 0._R8P ; sir(d) = 1._R8P
         call compute_convective_fluxes_maxwell_adim_div_b(sir=sir, q=qv(1:nv_c), f=fv(1:nv_c), chi=chi)
         do v=1, nv_c
            flxyz_c_gpu(b,1,d,i,j,k,v) = fv(v)
         enddo
      enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine fv_cell_fluxes_maxwell_adim_div_b_dev_kernel

   subroutine fv_cell_fluxes_maxwell_adim_div_d_b_dev_kernel(ni, nj, nk, ngc, nv_c, blocks_number, chi, q_gpu, flxyz_c_gpu)
   !< Pointwise adimensional Maxwell fluxes with D/B hyperbolic cleaning.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc, blocks_number
   integer(I4P), intent(in)    :: nv_c
   real(R8P),    intent(in)    :: chi
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(inout) :: flxyz_c_gpu(1:,1:,1:,1-ngc:,1-ngc:,1-ngc:,1:)
   integer(I4P), parameter     :: NVC_CAP=8_I4P
   real(R8P)                   :: qv(1:NVC_CAP), fv(1:NVC_CAP), sir(1:3)
   integer(I4P)                :: i, j, k, b, d, v

   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu, flxyz_c_gpu) &
   !$acc& firstprivate(ni, nj, nk, ngc, nv_c, blocks_number, chi)                        &
   !$acc& private(qv, fv, sir)
   !$omp OMPLOOP collapse(4) DEVICEPTR(q_gpu, flxyz_c_gpu) &
   !$omp& firstprivate(ni, nj, nk, ngc, nv_c, blocks_number, chi) &
   !$omp& private(qv, fv, sir)
   do b=1, blocks_number
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
      !$acc loop seq
      do v=1, nv_c
         qv(v) = q_gpu(b,i,j,k,v)
      enddo
      !$acc loop seq
      do d=1, 3
         sir(1) = 0._R8P ; sir(2) = 0._R8P ; sir(3) = 0._R8P ; sir(d) = 1._R8P
         call compute_convective_fluxes_maxwell_adim_div_d_b(sir=sir, q=qv(1:nv_c), f=fv(1:nv_c), chi=chi)
         do v=1, nv_c
            flxyz_c_gpu(b,1,d,i,j,k,v) = fv(v)
         enddo
      enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine fv_cell_fluxes_maxwell_adim_div_d_b_dev_kernel

   subroutine fv_recon_x_dev_kernel(s1, ni, nj, nk, ngc, nv_c, blocks_number, flxyz_c_gpu, flx_f_gpu)
   !< Reconstruct x-fluxes at x-faces from cell-center fluxes (issue #23 R2): per face
   !< i+1/2 (i = 0..ni), gather the 2*s1 cell-flux stencil into a constant-bound private
   !< buffer and call the m=0 SOTA face-reconstruction primitive (already acc routine seq).
   integer(I4P), intent(in)    :: s1                                                 !< Half FDV stencil length.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc, blocks_number                     !< Grids dimensions.
   integer(I4P), intent(in)    :: nv_c                                               !< Conservative variables number.
   real(R8P),    intent(in)    :: flxyz_c_gpu(1:,1:,1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Cell-center fluxes (b,1,d,i,j,k,v).
   real(R8P),    intent(inout) :: flx_f_gpu(1:,0:,1:,1:,1:)                          !< X-face fluxes (b,0:ni,j,k,v).
   real(R8P)                   :: qs(1-FDV_S_MAX:FDV_S_MAX)                          !< Private stencil gather.
   integer(I4P)                :: i, j, k, b, v, m                                   !< Counters.

   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(flxyz_c_gpu, flx_f_gpu) &
   !$acc& firstprivate(s1, ni, nj, nk, ngc, nv_c, blocks_number)                             &
   !$acc& private(qs)
   !$omp OMPLOOP collapse(4) DEVICEPTR(flxyz_c_gpu, flx_f_gpu) &
   !$omp& firstprivate(s1, ni, nj, nk, ngc, nv_c, blocks_number) &
   !$omp& private(qs)
   do b=1, blocks_number
   do k=1, nk
   do j=1, nj
   do i=0, ni
      !$acc loop seq
      do v=1, nv_c
         do m=1-s1, s1
            qs(m) = flxyz_c_gpu(b,1,1,i+m,j,k,v)
         enddo
         call compute_reconstruction_r_fv_centered(s=s1, q=qs(1-s1:s1), qr=flx_f_gpu(b,i,j,k,v))
      enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine fv_recon_x_dev_kernel

   subroutine fv_recon_y_dev_kernel(s1, ni, nj, nk, ngc, nv_c, blocks_number, flxyz_c_gpu, fly_f_gpu)
   !< Reconstruct y-fluxes at y-faces (issue #23 R2); twin of fv_recon_x_dev_kernel.
   integer(I4P), intent(in)    :: s1                                                 !< Half FDV stencil length.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc, blocks_number                     !< Grids dimensions.
   integer(I4P), intent(in)    :: nv_c                                               !< Conservative variables number.
   real(R8P),    intent(in)    :: flxyz_c_gpu(1:,1:,1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Cell-center fluxes (b,1,d,i,j,k,v).
   real(R8P),    intent(inout) :: fly_f_gpu(1:,1:,0:,1:,1:)                          !< Y-face fluxes (b,i,0:nj,k,v).
   real(R8P)                   :: qs(1-FDV_S_MAX:FDV_S_MAX)                          !< Private stencil gather.
   integer(I4P)                :: i, j, k, b, v, m                                   !< Counters.

   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(flxyz_c_gpu, fly_f_gpu) &
   !$acc& firstprivate(s1, ni, nj, nk, ngc, nv_c, blocks_number)                             &
   !$acc& private(qs)
   !$omp OMPLOOP collapse(4) DEVICEPTR(flxyz_c_gpu, fly_f_gpu) &
   !$omp& firstprivate(s1, ni, nj, nk, ngc, nv_c, blocks_number) &
   !$omp& private(qs)
   do b=1, blocks_number
   do k=1, nk
   do j=0, nj
   do i=1, ni
      !$acc loop seq
      do v=1, nv_c
         do m=1-s1, s1
            qs(m) = flxyz_c_gpu(b,1,2,i,j+m,k,v)
         enddo
         call compute_reconstruction_r_fv_centered(s=s1, q=qs(1-s1:s1), qr=fly_f_gpu(b,i,j,k,v))
      enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine fv_recon_y_dev_kernel

   subroutine fv_recon_z_dev_kernel(s1, ni, nj, nk, ngc, nv_c, blocks_number, flxyz_c_gpu, flz_f_gpu)
   !< Reconstruct z-fluxes at z-faces (issue #23 R2); twin of fv_recon_x_dev_kernel.
   integer(I4P), intent(in)    :: s1                                                 !< Half FDV stencil length.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc, blocks_number                     !< Grids dimensions.
   integer(I4P), intent(in)    :: nv_c                                               !< Conservative variables number.
   real(R8P),    intent(in)    :: flxyz_c_gpu(1:,1:,1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Cell-center fluxes (b,1,d,i,j,k,v).
   real(R8P),    intent(inout) :: flz_f_gpu(1:,1:,1:,0:,1:)                          !< Z-face fluxes (b,i,j,0:nk,v).
   real(R8P)                   :: qs(1-FDV_S_MAX:FDV_S_MAX)                          !< Private stencil gather.
   integer(I4P)                :: i, j, k, b, v, m                                   !< Counters.

   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(flxyz_c_gpu, flz_f_gpu) &
   !$acc& firstprivate(s1, ni, nj, nk, ngc, nv_c, blocks_number)                             &
   !$acc& private(qs)
   !$omp OMPLOOP collapse(4) DEVICEPTR(flxyz_c_gpu, flz_f_gpu) &
   !$omp& firstprivate(s1, ni, nj, nk, ngc, nv_c, blocks_number) &
   !$omp& private(qs)
   do b=1, blocks_number
   do k=0, nk
   do j=1, nj
   do i=1, ni
      !$acc loop seq
      do v=1, nv_c
         do m=1-s1, s1
            qs(m) = flxyz_c_gpu(b,1,3,i,j,k+m,v)
         enddo
         call compute_reconstruction_r_fv_centered(s=s1, q=qs(1-s1:s1), qr=flz_f_gpu(b,i,j,k,v))
      enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine fv_recon_z_dev_kernel

   subroutine fv_flux_diff_plain_dev_kernel(ni, nj, nk, ngc, nv_c, blocks_number, var_jx, var_jy, var_jz, &
                                            dxyz_gpu, flx_f_gpu, fly_f_gpu, flz_f_gpu, q_gpu, dq_gpu)
   !< Conservative flux difference + J source, no hyperbolic damping fields.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc, blocks_number
   integer(I4P), intent(in)    :: nv_c
   integer(I4P), intent(in)    :: var_jx, var_jy, var_jz
   real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)
   real(R8P),    intent(in)    :: flx_f_gpu(1:,0:,1:,1:,1:)
   real(R8P),    intent(in)    :: fly_f_gpu(1:,1:,0:,1:,1:)
   real(R8P),    intent(in)    :: flz_f_gpu(1:,1:,1:,0:,1:)
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(inout) :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P)                   :: dxb, dyb, dzb
   integer(I4P)                :: i, j, k, b, v

   !$acc parallel loop independent gang vector collapse(4)                              &
   !$acc& DEVICEVAR(dxyz_gpu, flx_f_gpu, fly_f_gpu, flz_f_gpu, q_gpu, dq_gpu)            &
   !$acc& firstprivate(ni, nj, nk, nv_c, blocks_number, var_jx, var_jy, var_jz)
   !$omp OMPLOOP collapse(4) &
   !$omp& DEVICEPTR(dxyz_gpu, flx_f_gpu, fly_f_gpu, flz_f_gpu, q_gpu, dq_gpu) &
   !$omp& firstprivate(ni, nj, nk, nv_c, blocks_number, var_jx, var_jy, var_jz)
   do b=1, blocks_number
   do k=1, nk
   do j=1, nj
   do i=1, ni
      dxb = dxyz_gpu(b,1) ; dyb = dxyz_gpu(b,2) ; dzb = dxyz_gpu(b,3)
      !$acc loop seq
      do v=1, nv_c
         dq_gpu(b,i,j,k,v) = - (flx_f_gpu(b,i,j,k,v) - flx_f_gpu(b,i-1,j,k,v)) / dxb &
                             - (fly_f_gpu(b,i,j,k,v) - fly_f_gpu(b,i,j-1,k,v)) / dyb &
                             - (flz_f_gpu(b,i,j,k,v) - flz_f_gpu(b,i,j,k-1,v)) / dzb
      enddo
      dq_gpu(b,i,j,k,VAR_DX) = dq_gpu(b,i,j,k,VAR_DX) - q_gpu(b,i,j,k,var_jx)
      dq_gpu(b,i,j,k,VAR_DY) = dq_gpu(b,i,j,k,VAR_DY) - q_gpu(b,i,j,k,var_jy)
      dq_gpu(b,i,j,k,VAR_DZ) = dq_gpu(b,i,j,k,VAR_DZ) - q_gpu(b,i,j,k,var_jz)
   enddo
   enddo
   enddo
   enddo
   endsubroutine fv_flux_diff_plain_dev_kernel

   subroutine fv_flux_diff_phi_dev_kernel(ni, nj, nk, ngc, nv_c, blocks_number, var_jx, var_jy, var_jz, &
                                          chi_damp, c_r, fv_ivar_phi, dxyz_gpu, flx_f_gpu, fly_f_gpu, flz_f_gpu, q_gpu, dq_gpu)
   !< Conservative flux difference + J source + phi damping.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc, blocks_number
   integer(I4P), intent(in)    :: nv_c
   integer(I4P), intent(in)    :: var_jx, var_jy, var_jz, fv_ivar_phi
   real(R8P),    intent(in)    :: chi_damp, c_r
   real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)
   real(R8P),    intent(in)    :: flx_f_gpu(1:,0:,1:,1:,1:)
   real(R8P),    intent(in)    :: fly_f_gpu(1:,1:,0:,1:,1:)
   real(R8P),    intent(in)    :: flz_f_gpu(1:,1:,1:,0:,1:)
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(inout) :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P)                   :: dxb, dyb, dzb, min_h, damping_coeff
   integer(I4P)                :: i, j, k, b, v

   !$acc parallel loop independent gang vector collapse(4)                              &
   !$acc& DEVICEVAR(dxyz_gpu, flx_f_gpu, fly_f_gpu, flz_f_gpu, q_gpu, dq_gpu)            &
   !$acc& firstprivate(ni, nj, nk, nv_c, blocks_number, var_jx, var_jy, var_jz, chi_damp, c_r, fv_ivar_phi)
   !$omp OMPLOOP collapse(4) &
   !$omp& DEVICEPTR(dxyz_gpu, flx_f_gpu, fly_f_gpu, flz_f_gpu, q_gpu, dq_gpu) &
   !$omp& firstprivate(ni, nj, nk, nv_c, blocks_number, var_jx, var_jy, var_jz, chi_damp, c_r, fv_ivar_phi)
   do b=1, blocks_number
   do k=1, nk
   do j=1, nj
   do i=1, ni
      dxb = dxyz_gpu(b,1) ; dyb = dxyz_gpu(b,2) ; dzb = dxyz_gpu(b,3)
      !$acc loop seq
      do v=1, nv_c
         dq_gpu(b,i,j,k,v) = - (flx_f_gpu(b,i,j,k,v) - flx_f_gpu(b,i-1,j,k,v)) / dxb &
                             - (fly_f_gpu(b,i,j,k,v) - fly_f_gpu(b,i,j-1,k,v)) / dyb &
                             - (flz_f_gpu(b,i,j,k,v) - flz_f_gpu(b,i,j,k-1,v)) / dzb
      enddo
      dq_gpu(b,i,j,k,VAR_DX) = dq_gpu(b,i,j,k,VAR_DX) - q_gpu(b,i,j,k,var_jx)
      dq_gpu(b,i,j,k,VAR_DY) = dq_gpu(b,i,j,k,VAR_DY) - q_gpu(b,i,j,k,var_jy)
      dq_gpu(b,i,j,k,VAR_DZ) = dq_gpu(b,i,j,k,VAR_DZ) - q_gpu(b,i,j,k,var_jz)
      if (c_r > 0._R8P) then
         min_h = min(dxb, min(dyb, dzb))
         damping_coeff = chi_damp / (c_r * min_h)
         dq_gpu(b,i,j,k,fv_ivar_phi) = dq_gpu(b,i,j,k,fv_ivar_phi) - damping_coeff * q_gpu(b,i,j,k,fv_ivar_phi)
      endif
   enddo
   enddo
   enddo
   enddo
   endsubroutine fv_flux_diff_phi_dev_kernel

   subroutine fv_flux_diff_psi_dev_kernel(ni, nj, nk, ngc, nv_c, blocks_number, var_jx, var_jy, var_jz, &
                                          chi_damp, c_r, fv_ivar_psi, dxyz_gpu, flx_f_gpu, fly_f_gpu, flz_f_gpu, q_gpu, dq_gpu)
   !< Conservative flux difference + J source + psi damping.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc, blocks_number
   integer(I4P), intent(in)    :: nv_c
   integer(I4P), intent(in)    :: var_jx, var_jy, var_jz, fv_ivar_psi
   real(R8P),    intent(in)    :: chi_damp, c_r
   real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)
   real(R8P),    intent(in)    :: flx_f_gpu(1:,0:,1:,1:,1:)
   real(R8P),    intent(in)    :: fly_f_gpu(1:,1:,0:,1:,1:)
   real(R8P),    intent(in)    :: flz_f_gpu(1:,1:,1:,0:,1:)
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(inout) :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P)                   :: dxb, dyb, dzb, min_h, damping_coeff
   integer(I4P)                :: i, j, k, b, v

   !$acc parallel loop independent gang vector collapse(4)                              &
   !$acc& DEVICEVAR(dxyz_gpu, flx_f_gpu, fly_f_gpu, flz_f_gpu, q_gpu, dq_gpu)            &
   !$acc& firstprivate(ni, nj, nk, nv_c, blocks_number, var_jx, var_jy, var_jz, chi_damp, c_r, fv_ivar_psi)
   !$omp OMPLOOP collapse(4) &
   !$omp& DEVICEPTR(dxyz_gpu, flx_f_gpu, fly_f_gpu, flz_f_gpu, q_gpu, dq_gpu) &
   !$omp& firstprivate(ni, nj, nk, nv_c, blocks_number, var_jx, var_jy, var_jz, chi_damp, c_r, fv_ivar_psi)
   do b=1, blocks_number
   do k=1, nk
   do j=1, nj
   do i=1, ni
      dxb = dxyz_gpu(b,1) ; dyb = dxyz_gpu(b,2) ; dzb = dxyz_gpu(b,3)
      !$acc loop seq
      do v=1, nv_c
         dq_gpu(b,i,j,k,v) = - (flx_f_gpu(b,i,j,k,v) - flx_f_gpu(b,i-1,j,k,v)) / dxb &
                             - (fly_f_gpu(b,i,j,k,v) - fly_f_gpu(b,i,j-1,k,v)) / dyb &
                             - (flz_f_gpu(b,i,j,k,v) - flz_f_gpu(b,i,j,k-1,v)) / dzb
      enddo
      dq_gpu(b,i,j,k,VAR_DX) = dq_gpu(b,i,j,k,VAR_DX) - q_gpu(b,i,j,k,var_jx)
      dq_gpu(b,i,j,k,VAR_DY) = dq_gpu(b,i,j,k,VAR_DY) - q_gpu(b,i,j,k,var_jy)
      dq_gpu(b,i,j,k,VAR_DZ) = dq_gpu(b,i,j,k,VAR_DZ) - q_gpu(b,i,j,k,var_jz)
      if (c_r > 0._R8P) then
         min_h = min(dxb, min(dyb, dzb))
         damping_coeff = chi_damp / (c_r * min_h)
         dq_gpu(b,i,j,k,fv_ivar_psi) = dq_gpu(b,i,j,k,fv_ivar_psi) - damping_coeff * q_gpu(b,i,j,k,fv_ivar_psi)
      endif
   enddo
   enddo
   enddo
   enddo
   endsubroutine fv_flux_diff_psi_dev_kernel

   subroutine fv_flux_diff_phi_psi_dev_kernel(ni, nj, nk, ngc, nv_c, blocks_number, var_jx, var_jy, var_jz, &
                                              chi_damp, c_r, fv_ivar_phi, fv_ivar_psi, dxyz_gpu, flx_f_gpu, fly_f_gpu, &
                                              flz_f_gpu, q_gpu, dq_gpu)
   !< Conservative flux difference + J source + phi/psi damping.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc, blocks_number
   integer(I4P), intent(in)    :: nv_c
   integer(I4P), intent(in)    :: var_jx, var_jy, var_jz, fv_ivar_phi, fv_ivar_psi
   real(R8P),    intent(in)    :: chi_damp, c_r
   real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)
   real(R8P),    intent(in)    :: flx_f_gpu(1:,0:,1:,1:,1:)
   real(R8P),    intent(in)    :: fly_f_gpu(1:,1:,0:,1:,1:)
   real(R8P),    intent(in)    :: flz_f_gpu(1:,1:,1:,0:,1:)
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(inout) :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P)                   :: dxb, dyb, dzb, min_h, damping_coeff
   integer(I4P)                :: i, j, k, b, v

   !$acc parallel loop independent gang vector collapse(4)                              &
   !$acc& DEVICEVAR(dxyz_gpu, flx_f_gpu, fly_f_gpu, flz_f_gpu, q_gpu, dq_gpu)            &
   !$acc& firstprivate(ni, nj, nk, nv_c, blocks_number, var_jx, var_jy, var_jz, chi_damp, c_r, &
   !$acc&              fv_ivar_phi, fv_ivar_psi)
   !$omp OMPLOOP collapse(4) &
   !$omp& DEVICEPTR(dxyz_gpu, flx_f_gpu, fly_f_gpu, flz_f_gpu, q_gpu, dq_gpu) &
   !$omp& firstprivate(ni, nj, nk, nv_c, blocks_number, var_jx, var_jy, var_jz, chi_damp, c_r, &
   !$omp&              fv_ivar_phi, fv_ivar_psi)
   do b=1, blocks_number
   do k=1, nk
   do j=1, nj
   do i=1, ni
      dxb = dxyz_gpu(b,1) ; dyb = dxyz_gpu(b,2) ; dzb = dxyz_gpu(b,3)
      !$acc loop seq
      do v=1, nv_c
         dq_gpu(b,i,j,k,v) = - (flx_f_gpu(b,i,j,k,v) - flx_f_gpu(b,i-1,j,k,v)) / dxb &
                             - (fly_f_gpu(b,i,j,k,v) - fly_f_gpu(b,i,j-1,k,v)) / dyb &
                             - (flz_f_gpu(b,i,j,k,v) - flz_f_gpu(b,i,j,k-1,v)) / dzb
      enddo
      dq_gpu(b,i,j,k,VAR_DX) = dq_gpu(b,i,j,k,VAR_DX) - q_gpu(b,i,j,k,var_jx)
      dq_gpu(b,i,j,k,VAR_DY) = dq_gpu(b,i,j,k,VAR_DY) - q_gpu(b,i,j,k,var_jy)
      dq_gpu(b,i,j,k,VAR_DZ) = dq_gpu(b,i,j,k,VAR_DZ) - q_gpu(b,i,j,k,var_jz)
      if (c_r > 0._R8P) then
         min_h = min(dxb, min(dyb, dzb))
         damping_coeff = chi_damp / (c_r * min_h)
         dq_gpu(b,i,j,k,fv_ivar_phi) = dq_gpu(b,i,j,k,fv_ivar_phi) - damping_coeff * q_gpu(b,i,j,k,fv_ivar_phi)
         dq_gpu(b,i,j,k,fv_ivar_psi) = dq_gpu(b,i,j,k,fv_ivar_psi) - damping_coeff * q_gpu(b,i,j,k,fv_ivar_psi)
      endif
   enddo
   enddo
   enddo
   enddo
   endsubroutine fv_flux_diff_phi_psi_dev_kernel

   subroutine accumulate_seam_fluxes_fv_dev(self, flux_register)
   !< Accumulate end-of-step FV seam face fluxes into the forest's flux register —
   !< FNL twin of `prism_cpu_object%accumulate_seam_fluxes_fv` (issue #23 R3).
   !<
   !< Register crossing (grilled R0 default): the flux register stays HOST-side.
   !< Per seam (block, fec) hit, a tiny device kernel packs the face skin into a
   !< contiguous device slab (inner axis fastest — the register's cell order), ONE
   !< D2H copy (~nv*nface_cells doubles, tens of KB per step) brings it to the
   !< host, and the shared register machinery does the rest: coarse side →
   !< `accumulate_coarse_flux` directly; fine side → quadrant offset from
   !< `maps%amr_seam_quadrant` (Morton-code precompute, issue #28 D2) + the pure 2:1 restriction
   !< `restrict_fine_face_to_quadrant` (moved to adam_flux_register_object, shared
   !< with the CPU backend) → `accumulate_fine_flux`. Fires once per realm per
   !< step (final-substage gate at the call site). Third-axis register index
   !< hardcoded to 1 (α.r1 collapsed register), exactly as on CPU.
   class(prism_fnl_object),     intent(inout) :: self          !< The realm.
   class(flux_register_object), intent(inout) :: flux_register !< Forest's flux register.
   real(R8P), pointer                         :: skin_gpu(:,:)   !< Device face-skin slab (nface_cells, nv_c).
   real(R8P), allocatable                     :: skin(:,:)       !< Host copy of the packed skin.
   real(R8P), allocatable                     :: flux_slab(:,:)  !< Register-shaped contribution (nv_reg, nface_cells).
   real(R8P), allocatable                     :: fine_face(:,:,:)!< Fine face flux (nv_c, inner_n, outer_n).
   integer(I4P)                               :: sgn_idx         !< Signed register index for (b, fec).
   integer(I4P)                               :: face_idx        !< |sgn_idx| → register face.
   integer(I4P)                               :: nv_reg          !< Register state-vector width.
   integer(I4P)                               :: nface_cells     !< Coarse-face skin cell count.
   integer(I4P)                               :: fec, b          !< Face, block counters.
   integer(I4P)                               :: inner_n, outer_n   !< Cell counts along tangential axes.
   integer(I4P)                               :: ioff, joff      !< Fine-block quadrant offset ∈ {0,1}.
   integer(I4P)                               :: c, v, fi, fo    !< Packing counters.
   integer(I4P)                               :: ierr            !< Device allocation error flag.

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, nv_c=>self%nv_c)
   do b=1, self%blocks_number
      do fec=1, 6
         sgn_idx = self%adam%maps%inter_realm_face_register_index(b, fec)
         if (sgn_idx == 0_I4P) cycle
         face_idx = abs(sgn_idx)
         if (face_idx > flux_register%nfaces) cycle
         if (.not. allocated(flux_register%face(face_idx)%F_coarse)) cycle
         nv_reg      = int(size(flux_register%face(face_idx)%F_coarse, dim=1), I4P)
         nface_cells = flux_register%face(face_idx)%nface_cells
         ! Tangential extents per face (inner fastest in the linear index c) — CPU parity.
         select case (fec)
         case (1_I4P, 2_I4P) ; inner_n = nj ; outer_n = nk
         case (3_I4P, 4_I4P) ; inner_n = ni ; outer_n = nk
         case (5_I4P, 6_I4P) ; inner_n = ni ; outer_n = nj
         case default ; cycle
         end select

         call dev_alloc(fptr_dev=skin_gpu, lbounds=[1,1], ubounds=[inner_n*outer_n, nv_c], ierr=ierr)
         if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate skin_gpu in accumulate_seam_fluxes_fv')
         call fv_pack_face_skin_dev_kernel(fec=fec, ni=ni, nj=nj, nk=nk, nv_c=nv_c, b=b,       &
                                           inner_n=inner_n, outer_n=outer_n,                   &
                                           flx_f_gpu=self%flx_f_gpu, fly_f_gpu=self%fly_f_gpu, &
                                           flz_f_gpu=self%flz_f_gpu, skin_gpu=skin_gpu)
         allocate(skin(1:inner_n*outer_n, 1:nv_c))
         call dev_memcpy_from_device(dst=skin, src=skin_gpu)
         call dev_free(skin_gpu, mydev)

         allocate(flux_slab(1:nv_reg, 1:nface_cells))
         flux_slab = 0._R8P
         if (sgn_idx > 0_I4P) then
            ! Coarse side: the packed skin IS the register-order slab (transpose only).
            do c=1_I4P, nface_cells
               do v=1_I4P, nv_c
                  flux_slab(v, c) = skin(c, v)
               enddo
            enddo
            call flux_register%accumulate_coarse_flux(face_index=face_idx, stage=1_I4P, flux_face=flux_slab)
         else
            ! Fine side: reshape to (nv_c, inner, outer), 2:1-restrict into this
            ! block's quadrant of the coarse skin. Quadrant offsets are read from
            ! `maps%amr_seam_quadrant`, precomputed at registration from Morton
            ! codes (issue #28 D2, CPU parity): the register's `coarse_block` is
            ! an owner-rank-LOCAL index, so deriving the quadrant from its
            ! emin/emax here reads an unrelated local block's geometry whenever
            ! the coarse partner lives on another rank. Table allocated only by
            ! the intra-realm AMR pass; inter-realm mirror seams have no quadrant.
            if (allocated(self%adam%maps%amr_seam_quadrant)) then
               ioff = self%adam%maps%amr_seam_quadrant(1, b, fec)
               joff = self%adam%maps%amr_seam_quadrant(2, b, fec)
            else
               ioff = 0_I4P ; joff = 0_I4P
            endif
            allocate(fine_face(1:nv_c, 1:inner_n, 1:outer_n))
            do fo=1_I4P, outer_n
               do fi=1_I4P, inner_n
                  c = (fo - 1_I4P)*inner_n + fi
                  do v=1_I4P, nv_c
                     fine_face(v, fi, fo) = skin(c, v)
                  enddo
               enddo
            enddo
            call restrict_fine_face_to_quadrant(fine_face=fine_face, inner_n=inner_n, outer_n=outer_n, &
                                                ioff=ioff, joff=joff, slab=flux_slab)
            call flux_register%accumulate_fine_flux(face_index=face_idx, stage=1_I4P, flux_face=flux_slab)
            deallocate(fine_face)
         endif
         deallocate(skin, flux_slab)
      enddo
   enddo
   endassociate
   endsubroutine accumulate_seam_fluxes_fv_dev

   subroutine fv_pack_face_skin_dev_kernel(fec, ni, nj, nk, nv_c, b, inner_n, outer_n, &
                                           flx_f_gpu, fly_f_gpu, flz_f_gpu, skin_gpu)
   !< Pack one block's face flux skin into a contiguous device slab (issue #23 R3):
   !< `skin_gpu(c, v)` with `c = (outer-1)*inner_n + inner` — the register cell order,
   !< identical to the CPU pack (pack_coarse_face / fine_face_cell mappings). Scalar
   !< element reads/writes only: no private arrays, no section actuals. The `fec`
   !< branch is uniform per launch (firstprivate scalar). Tiny kernel: one launch per
   !< seam (block, fec) hit, once per step.
   integer(I4P), intent(in)    :: fec                        !< Face code 1..6 (-x,+x,-y,+y,-z,+z).
   integer(I4P), intent(in)    :: ni, nj, nk                 !< Interior cell counts.
   integer(I4P), intent(in)    :: nv_c                       !< Conservative variables number.
   integer(I4P), intent(in)    :: b                          !< Block index.
   integer(I4P), intent(in)    :: inner_n, outer_n           !< Tangential cell counts (inner fastest).
   real(R8P),    intent(in)    :: flx_f_gpu(1:,0:,1:,1:,1:)  !< X-face fluxes.
   real(R8P),    intent(in)    :: fly_f_gpu(1:,1:,0:,1:,1:)  !< Y-face fluxes.
   real(R8P),    intent(in)    :: flz_f_gpu(1:,1:,1:,0:,1:)  !< Z-face fluxes.
   real(R8P),    intent(inout) :: skin_gpu(1:,1:)            !< Packed skin (inner_n*outer_n, nv_c).
   integer(I4P)                :: ii, o, v                   !< Counters.

   !$acc parallel loop independent gang vector collapse(2)                    &
   !$acc& DEVICEVAR(flx_f_gpu, fly_f_gpu, flz_f_gpu, skin_gpu)                &
   !$acc& firstprivate(fec, ni, nj, nk, nv_c, b, inner_n, outer_n)
   !$omp OMPLOOP collapse(2) &
   !$omp& DEVICEPTR(flx_f_gpu, fly_f_gpu, flz_f_gpu, skin_gpu) &
   !$omp& firstprivate(fec, ni, nj, nk, nv_c, b, inner_n, outer_n)
   do o=1, outer_n
   do ii=1, inner_n
      select case (fec)
      case (1_I4P)
         do v=1, nv_c ; skin_gpu((o-1)*inner_n+ii, v) = flx_f_gpu(b, 0,  ii, o, v) ; enddo ! -x: inner=j, outer=k
      case (2_I4P)
         do v=1, nv_c ; skin_gpu((o-1)*inner_n+ii, v) = flx_f_gpu(b, ni, ii, o, v) ; enddo ! +x
      case (3_I4P)
         do v=1, nv_c ; skin_gpu((o-1)*inner_n+ii, v) = fly_f_gpu(b, ii, 0,  o, v) ; enddo ! -y: inner=i, outer=k
      case (4_I4P)
         do v=1, nv_c ; skin_gpu((o-1)*inner_n+ii, v) = fly_f_gpu(b, ii, nj, o, v) ; enddo ! +y
      case (5_I4P)
         do v=1, nv_c ; skin_gpu((o-1)*inner_n+ii, v) = flz_f_gpu(b, ii, o, 0,  v) ; enddo ! -z: inner=i, outer=j
      case (6_I4P)
         do v=1, nv_c ; skin_gpu((o-1)*inner_n+ii, v) = flz_f_gpu(b, ii, o, nk, v) ; enddo ! +z
      endselect
   enddo
   enddo
   endsubroutine fv_pack_face_skin_dev_kernel

   subroutine fv_apply_reflux_face_dev_kernel(axis, sgn, b, ni, nj, nk, ngc, nv_reg, nface_cells, &
                                              scale, delta_gpu, q_gpu)
   !< Apply the Berger-Colella end-of-step correction for ONE register face to the
   !< committed q_gpu (issue #23 R4): per skin cell c (register cell order — the
   !< SAME (axis, sgn) → (i,j,k) mapping as the CPU apply), add
   !< scale * delta(v, c) to every state row. Scalar element ops only: no private
   !< arrays, no section actuals; each c writes a distinct cell (disjoint), the
   !< inner v loop is per-thread seq. Tiny kernel — one launch per face per step.
   integer(I4P), intent(in)    :: axis, sgn                  !< Face normal axis (1..3) and sign (+-1).
   integer(I4P), intent(in)    :: b                          !< Coarse block index.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc            !< Grids dimensions.
   integer(I4P), intent(in)    :: nv_reg                     !< Register state-vector width.
   integer(I4P), intent(in)    :: nface_cells                !< Skin cell count.
   real(R8P),    intent(in)    :: scale                      !< sign * dt / dx_coarse.
   real(R8P),    intent(in)    :: delta_gpu(1:,1:)           !< Flux mismatch slab (nv_reg, nface_cells).
   real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Committed conservative variables.
   integer(I4P)                :: c, c0, i, j, k, v          !< Counters / decoded cell indexes.

   !$acc parallel loop independent gang vector          &
   !$acc& DEVICEVAR(delta_gpu, q_gpu)                   &
   !$acc& firstprivate(axis, sgn, b, ni, nj, nk, nv_reg, nface_cells, scale)
   !$omp OMPLOOP DEVICEPTR(delta_gpu, q_gpu) &
   !$omp& firstprivate(axis, sgn, b, ni, nj, nk, nv_reg, nface_cells, scale)
   do c=1, nface_cells
      c0 = c - 1
      select case (axis)
      case (1_I4P)  ! x-normal face: i fixed; tangentials (j, k) walk (mod nj, div nj).
         i = merge(ni, 1_I4P, sgn > 0_I4P) ; j = 1_I4P + mod(c0, nj) ; k = 1_I4P + c0/nj
      case (2_I4P)  ! y-normal face: j fixed; tangentials (i, k) walk (mod ni, div ni).
         i = 1_I4P + mod(c0, ni) ; j = merge(nj, 1_I4P, sgn > 0_I4P) ; k = 1_I4P + c0/ni
      case (3_I4P)  ! z-normal face: k fixed; tangentials (i, j) walk (mod ni, div ni).
         i = 1_I4P + mod(c0, ni) ; j = 1_I4P + c0/ni ; k = merge(nk, 1_I4P, sgn > 0_I4P)
      case default
         i = 0_I4P ; j = 0_I4P ; k = 0_I4P
      endselect
      if (i /= 0_I4P) then
         !$acc loop seq
         do v=1, nv_reg
            q_gpu(b,i,j,k,v) = q_gpu(b,i,j,k,v) + scale * delta_gpu(v,c)
         enddo
      endif
   enddo
   endsubroutine fv_apply_reflux_face_dev_kernel

   ! numerical methods, time operators
   subroutine integrate_blanesmoan_dev(self)
   !< Integrate equation, time operator, Blanes and Moan scheme.
   class(prism_fnl_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: s    !< Counter.

   ! associate(nc=>self%blanesmoan%nc,a=>self%blanesmoan%a,b=>self%blanesmoan%b)
   ! call self%compute_coils_current(q_gpu=self%q_gpu)
   ! do s=1, nc
   !    call self%compute_residuals_dev(q=self%q, dq=self%dq)
   !    if (s==1) call self%save_residuals
   !    self%q(VAR_BX,:,:,:,:) = self%q(VAR_BX,:,:,:,:) + b(s) * self%time%dt * self%dq(VAR_BX,:,:,:,:)
   !    self%q(VAR_BY,:,:,:,:) = self%q(VAR_BY,:,:,:,:) + b(s) * self%time%dt * self%dq(VAR_BY,:,:,:,:)
   !    self%q(VAR_BZ,:,:,:,:) = self%q(VAR_BZ,:,:,:,:) + b(s) * self%time%dt * self%dq(VAR_BZ,:,:,:,:)
   !    call self%compute_residuals_dev(q=self%q, dq=self%dq)
   !    self%q(VAR_DX,:,:,:,:) = self%q(VAR_DX,:,:,:,:) + a(s) * self%time%dt * self%dq(VAR_DX,:,:,:,:)
   !    self%q(VAR_DY,:,:,:,:) = self%q(VAR_DY,:,:,:,:) + a(s) * self%time%dt * self%dq(VAR_DY,:,:,:,:)
   !    self%q(VAR_DZ,:,:,:,:) = self%q(VAR_DZ,:,:,:,:) + a(s) * self%time%dt * self%dq(VAR_DZ,:,:,:,:)
   ! enddo
   ! call self%impose_div_free
   ! endassociate
   endsubroutine integrate_blanesmoan_dev

   subroutine integrate_cfm_dev(self)
   !< Integrate equation, time operator, Commutator-Free Magnus integrator.
   class(prism_fnl_object), intent(inout) :: self             !< The equation.
   real(R8P), parameter                   :: toll=1.0e-14_R8P !< CFM coefficients tollerance.
   integer(I4P)                           :: s,ss             !< Counter.

   ! call self%compute_coils_current(q_gpu=self%q_gpu)
   ! associate(dt=>self%time%dt,s_coeffs=>self%cfm%s_coeffs,e_coeffs=>self%cfm%e_coeffs)
   ! self%cfm%q = self%q
   ! call self%compute_residuals_dev(q=self%cfm%q, dq=self%cfm%dq(:,:,:,:,:,1))
   ! do s=2, self%cfm%n_stages
   !    do ss=1, s-1
   !       if (abs(s_coeffs(s,ss))>toll) &
   !          call self%cfm%compute_exponential_update(alpha=dt*s_coeffs(s,ss),dq=self%cfm%dq(:,:,:,:,:,ss),q=self%cfm%q)
   !    enddo
   !    call self%compute_residuals_dev(q=self%cfm%q, dq=self%cfm%dq(:,:,:,:,:,s))
   ! enddo
   ! self%cfm%q = self%q
   ! do s=1, self%cfm%n_stages
   !    if (abs(self%cfm%e_coeffs(s))>toll) &
   !       call self%cfm%compute_exponential_update(alpha=dt*e_coeffs(s),dq=self%cfm%dq(:,:,:,:,:,s),q=self%cfm%q)
   ! enddo
   ! self%q = self%cfm%q
   ! endassociate
   ! call self%impose_div_free
   endsubroutine integrate_cfm_dev

   subroutine integrate_leapfrog_dev(self)
   !< Integrate equation, time operator, leapfrog scheme.
   class(prism_fnl_object), intent(inout) :: self !< The equation.

   ! call self%compute_coils_current(q_gpu=self%q_gpu)
   ! call self%compute_residuals_dev(q=self%q, dq=self%dq)
   ! call self%save_residuals
   ! call self%leapfrog%integrate(dt=self%time%dt, q=self%q, dq=self%dq)
   ! call self%impose_div_free
   endsubroutine integrate_leapfrog_dev

   subroutine integrate_leapfrog_pic(self)
   !< Integrate equation, time operator, leapfrog scheme.
   class(prism_fnl_object), intent(inout) :: self !< The equation.

   !!Fai check su come si parlano questa subroutine e quella che calcola la corrente associata alle particelle
   !call self%compute_coils_current(q_gpu=self%q_gpu)

   !! qua ci va la chiamata alla subroutine che aggiorna la neighbour list delle particelle
   !! qua ci va la chiamata alla subroutine che calcola la corrente associata alle particelle

   !call self%compute_residuals_dev(q=self%q, dq=self%dq) !< Calcolo i residui relativi ai campi E e B
   !call self%save_residuals

   !call self%leapfrog%integrate(dt=self%time%dt, q=self%q, dq=self%dq)

   !!Qua ci va la chiamata alla subroutine che integra aggiornando le velocita e le posizioni delle particelle
   !            !decidi se fare qui all'inizio del tempo successivo l'aggiornamento della neighbour list delle particelle

   !call self%impose_div_free
   !!call self%apply_fWL_correction
   endsubroutine integrate_leapfrog_pic

   subroutine integrate_rk_ssp_pic(self)
   !< Integrate PIC equations with SSP RK on device.
   class(prism_fnl_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: s    !< Counter.

   call self%rk_fnl%initialize_stages(grid=self%adam%grid, field=self%adam%field, q_gpu=self%q_gpu)
   call self%rk_pic_fnl%initialize_stages(q_pic_gpu=self%pic_fnl%q_pic_gpu)
   if (self%pml_fnl%enabled) call self%rk_pml_fnl%initialize_stages(pml_fnl=self%pml_fnl)

   do s=1, self%rk%nrk
      if (self%ib%solids_number>0) then
         call self%rk_fnl%compute_stage(grid=self%adam%grid, field=self%adam%field, s=s, dt=self%time%dt, &
                                        phi_gpu=self%ib_fnl%phi_gpu)
      else
         call self%rk_fnl%compute_stage(grid=self%adam%grid, field=self%adam%field, s=s, dt=self%time%dt)
      endif
      call self%rk_pic_fnl%compute_stage(s=s, dt=self%time%dt)
      if (self%pml_fnl%enabled) call self%rk_pml_fnl%compute_stage(s=s, dt=self%time%dt)

      call self%pic_fnl%particle_cartesian_grid_index_dev(field_fnl=self%field_fnl, field=self%adam%field, &
                                                          grid=self%adam%grid, q_pic_gpu=self%rk_pic_fnl%q_pic_rk_gpu(:,:,s))
      call self%pic_fnl%current_weighting_dev(field_fnl=self%field_fnl, field=self%adam%field, grid=self%adam%grid, &
                                              q_gpu=self%rk_fnl%q_rk_gpu(:,:,:,:,:,s), &
                                              q_pic_gpu=self%rk_pic_fnl%q_pic_rk_gpu(:,:,s), nv=self%nv)
      call self%verify_no_pic_deposition_on_coils_dev(q_gpu=self%rk_fnl%q_rk_gpu(:,:,:,:,:,s), check_current=.true., &
                                                      context='integrate_rk_ssp_pic(stage current)')
      call self%compute_coils_current(q_gpu=self%rk_fnl%q_rk_gpu(:,:,:,:,:,s), gamm=self%rk%gamm(s))
      call self%compute_residuals_dev(q_gpu=self%rk_fnl%q_rk_gpu(:,:,:,:,:,s), dq_gpu=self%dq_gpu, s=s)
      if (s==1) call self%save_residuals

      call self%pic_fnl%field_weighting_dev(field_fnl=self%field_fnl, field=self%adam%field, grid=self%adam%grid, &
                                            pic_fields_gpu=self%pic_fnl%pic_fields_gpu, &
                                            q_gpu=self%rk_fnl%q_rk_gpu(:,:,:,:,:,s), &
                                            q_pic_gpu=self%rk_pic_fnl%q_pic_rk_gpu(:,:,s), nv=self%nv)

      if (self%ib%solids_number>0) then
         call self%rk_fnl%assign_stage(grid=self%adam%grid, field=self%adam%field, s=s, q_gpu=self%dq_gpu, &
                                       phi_gpu=self%ib_fnl%phi_gpu)
      else
         call self%rk_fnl%assign_stage(grid=self%adam%grid, field=self%adam%field, s=s, q_gpu=self%dq_gpu)
      endif
      call self%rk_pic_fnl%assign_stage(s=s, pic_fields_gpu=self%pic_fnl%pic_fields_gpu)
      if (self%pml_fnl%enabled) call self%rk_pml_fnl%assign_stage(s=s)
   enddo

   if (self%ib%solids_number>0) then
      call self%rk_fnl%update_q(grid=self%adam%grid, field=self%adam%field, rk=self%rk, dt=self%time%dt, &
                                phi_gpu=self%ib_fnl%phi_gpu, q_gpu=self%q_gpu)
   else
      call self%rk_fnl%update_q(grid=self%adam%grid, field=self%adam%field, rk=self%rk, dt=self%time%dt, q_gpu=self%q_gpu)
   endif
   call self%rk_pic_fnl%update_q_pic(dt=self%time%dt, q_pic_gpu=self%pic_fnl%q_pic_gpu)
   if (self%pml_fnl%enabled) call self%rk_pml_fnl%update_q_pml(dt=self%time%dt, pml_fnl=self%pml_fnl)

   call self%apply_fwl_correction(q_gpu=self%q_gpu)
   call self%impose_div_free
   call self%pic_fnl%particle_cartesian_grid_index_dev(field_fnl=self%field_fnl, field=self%adam%field, &
                                                       grid=self%adam%grid, q_pic_gpu=self%pic_fnl%q_pic_gpu)
   call self%pic_fnl%current_weighting_dev(field_fnl=self%field_fnl, field=self%adam%field, grid=self%adam%grid, &
                                           q_gpu=self%q_gpu, q_pic_gpu=self%pic_fnl%q_pic_gpu, nv=self%nv)
   call self%pic_fnl%particle_weighting_dev(field_fnl=self%field_fnl, field=self%adam%field, grid=self%adam%grid, &
                                            q_gpu=self%q_gpu, q_pic_gpu=self%pic_fnl%q_pic_gpu, nv=self%nv)
   call self%verify_no_pic_deposition_on_coils_dev(q_gpu=self%q_gpu, check_current=.true., check_charge=.true., &
                                                   context='integrate_rk_ssp_pic(final deposition)')
   call self%compute_coils_current(q_gpu=self%q_gpu)
   endsubroutine integrate_rk_ssp_pic

   subroutine integrate_rk_ls_dev(self)
   !< Integrate equation, time operator, RK classical low storage schemes.
   !< Low storage RK working on q_rk(:,:,:,:,:,1)/q as stages, update q in place.
   class(prism_fnl_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: s    !< Counter.

   call self%compute_coils_current(q_gpu=self%q_gpu)
   call self%rk_fnl%initialize_stages(grid=self%adam%grid, field=self%adam%field, q_gpu=self%q_gpu)
   do s=1, self%rk%nrk
      call self%compute_residuals_dev(q_gpu=self%q_gpu, dq_gpu=self%dq_gpu, s=s)
      if (s==1) call self%save_residuals
      if (self%ib%solids_number>0) then
         call self%rk_fnl%compute_stage_ls(grid=self%adam%grid, field=self%adam%field, rk=self%rk, s=s, dt=self%time%dt, &
                                           phi_gpu=self%ib_fnl%phi_gpu, &
                                           dq_gpu=self%dq_gpu, q_gpu=self%q_gpu)
      else
         call self%rk_fnl%compute_stage_ls(grid=self%adam%grid, field=self%adam%field, rk=self%rk, s=s, dt=self%time%dt, &
                                           dq_gpu=self%dq_gpu, q_gpu=self%q_gpu)
      endif
   enddo
   call self%apply_fwl_correction(q_gpu=self%q_gpu)
   call self%compute_coils_current(q_gpu=self%q_gpu)
   call self%impose_div_free
   endsubroutine integrate_rk_ls_dev

   subroutine integrate_rk_ssp_dev(self)
   !< Integrate equation, time operator, SSP RK schemes.
   !< SSP RK working on q_rk as stages.
   class(prism_fnl_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: s    !< Counter.

   if (self%external_fields%ef_type/=EF_TYPE_NONE) &
      call sub_external_fields_dev(external_fields=self%external_fields, field_gpu=self%field_fnl, &
                                   dt=self%time%dt, time=self%time%time, q_gpu=self%q_gpu)
   call self%rk_fnl%initialize_stages(grid=self%adam%grid, field=self%adam%field, q_gpu=self%q_gpu)
   if (self%pml_fnl%enabled) call self%rk_pml_fnl%initialize_stages(pml_fnl=self%pml_fnl)
   do s=1, self%rk%nrk
      if (self%ib%solids_number>0) then
         call self%rk_fnl%compute_stage(grid=self%adam%grid, field=self%adam%field, s=s, dt=self%time%dt, &
                                        phi_gpu=self%ib_fnl%phi_gpu)
      else
         call self%rk_fnl%compute_stage(grid=self%adam%grid, field=self%adam%field, s=s, dt=self%time%dt)
      endif
      if (self%pml_fnl%enabled) call self%rk_pml_fnl%compute_stage(s=s, dt=self%time%dt)
      call self%compute_coils_current(q_gpu=self%rk_fnl%q_rk_gpu(:,:,:,:,:,s), gamm=self%rk%gamm(s))
      call self%compute_residuals_dev(q_gpu=self%rk_fnl%q_rk_gpu(:,:,:,:,:,s), dq_gpu=self%dq_gpu, s=s)
      ! if (s==1) call self%save_residuals
      if (self%ib%solids_number>0) then
         call self%rk_fnl%assign_stage(grid=self%adam%grid, field=self%adam%field, s=s, q_gpu=self%dq_gpu, &
                                       phi_gpu=self%ib_fnl%phi_gpu)
      else
         call self%rk_fnl%assign_stage(grid=self%adam%grid, field=self%adam%field, s=s, q_gpu=self%dq_gpu)
      endif
      if (self%pml_fnl%enabled) call self%rk_pml_fnl%assign_stage(s=s)
   enddo
   if (self%ib%solids_number>0) then
      call self%rk_fnl%update_q(grid=self%adam%grid, field=self%adam%field, rk=self%rk, dt=self%time%dt, &
                                phi_gpu=self%ib_fnl%phi_gpu, q_gpu=self%q_gpu)
      ! call self%update_rk_ghost(dt=self%time%dt, phi_gpu=ib_fnl%phi_gpu)
   else
      call self%rk_fnl%update_q(grid=self%adam%grid, field=self%adam%field, rk=self%rk, dt=self%time%dt, q_gpu=self%q_gpu)
      ! call self%update_rk_ghost(dt=self%time%dt)
      call self%save_residuals
   endif
   if (self%pml_fnl%enabled) call self%rk_pml_fnl%update_q_pml(dt=self%time%dt, pml_fnl=self%pml_fnl)
   call self%apply_fwl_correction(q_gpu=self%q_gpu)
   call self%compute_coils_current(q_gpu=self%q_gpu)
   call self%impose_div_free
   if (self%external_fields%ef_type/=EF_TYPE_NONE) &
      call add_external_fields_dev(external_fields=self%external_fields, field_gpu=self%field_fnl, &
                                   dt=self%time%dt, time=self%time%time, q_gpu=self%q_gpu)
   endsubroutine integrate_rk_ssp_dev

   subroutine integrate_rk_yoshida_dev(self)
   !< Integrate equation, time operator, Yoshida RK scheme.
   class(prism_fnl_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: s    !< Counter.

   ! call self%compute_coils_current(q_gpu=self%q_gpu)
   ! do s=1, rk%nrk - 1
   !    call self%compute_residuals_dev(q=self%q, dq=self%dq)
   !    if (s==1) call self%save_residuals
   !    self%q(VAR_BX,:,:,:,:) = self%q(VAR_BX,:,:,:,:) + rk%ssa(s) * self%time%dt * self%dq(VAR_BX,:,:,:,:)
   !    self%q(VAR_BY,:,:,:,:) = self%q(VAR_BY,:,:,:,:) + rk%ssa(s) * self%time%dt * self%dq(VAR_BY,:,:,:,:)
   !    self%q(VAR_BZ,:,:,:,:) = self%q(VAR_BZ,:,:,:,:) + rk%ssa(s) * self%time%dt * self%dq(VAR_BZ,:,:,:,:)
   !    call self%compute_residuals_dev(q=self%q, dq=self%dq)
   !    self%q(VAR_DX,:,:,:,:) = self%q(VAR_DX,:,:,:,:) + rk%ssb(s) * self%time%dt * self%dq(VAR_DX,:,:,:,:)
   !    self%q(VAR_DY,:,:,:,:) = self%q(VAR_DY,:,:,:,:) + rk%ssb(s) * self%time%dt * self%dq(VAR_DY,:,:,:,:)
   !    self%q(VAR_DZ,:,:,:,:) = self%q(VAR_DZ,:,:,:,:) + rk%ssb(s) * self%time%dt * self%dq(VAR_DZ,:,:,:,:)
   ! enddo
   ! call self%compute_residuals_dev(q=self%q, dq=self%dq)
   ! self%q(VAR_BX,:,:,:,:) = self%q(VAR_BX,:,:,:,:) + rk%ssa(rk%nrk) * self%time%dt * self%dq(VAR_BX,:,:,:,:)
   ! self%q(VAR_BY,:,:,:,:) = self%q(VAR_BY,:,:,:,:) + rk%ssa(rk%nrk) * self%time%dt * self%dq(VAR_BY,:,:,:,:)
   ! self%q(VAR_BZ,:,:,:,:) = self%q(VAR_BZ,:,:,:,:) + rk%ssa(rk%nrk) * self%time%dt * self%dq(VAR_BZ,:,:,:,:)
   ! call self%impose_div_free
   endsubroutine integrate_rk_yoshida_dev

   subroutine is_done_forest(self, done)
   !< Decide whether this realm has reached its local termination criterion.
   !<
   !< Invoked by forest%is_done. PRISM-FNL override: matches the legacy
   !< condition inline in simulate — terminate when either the simulated
   !< time has reached self%time%time_max (time-driven mode, it_max <= 0) or
   !< the iteration count has reached self%time%it_max (iteration-driven mode).
   !< Today the test reads time-state from the `time` module singleton,
   !< so `self` is unused; once the forest takes over time bookkeeping
   !< the body will consume self%time%* instead.
   class(prism_fnl_object), intent(in)  :: self !< The realm.
   logical,                 intent(out) :: done !< True if this realm is done evolving.

   done = ((self%time%it_max <= 0).and.(self%time%time >= self%time%time_max)).or.&
          ((self%time%it >= self%time%it_max).and.(self%time%it_max > 0))
   endsubroutine is_done_forest

   subroutine finalize_forest(self)
   !< Shut this realm down: final state dump, close residuals file, finalize
   !< MPI handler.
   class(prism_fnl_object), intent(inout) :: self !< The realm.
   logical                                :: is_open

   !call self%compute_energy_error !Cazzo
   call self%save_simulation_data
   call self%io%close_file_residuals
   if (mpih_fnl%myrank == 0) then
      inquire(unit=self%io%energy_history_unit, opened=is_open)
      if (is_open) close(self%io%energy_history_unit)
      if (self%grms%history_unit > 0_I4P) then
         inquire(unit=self%grms%history_unit, opened=is_open)
         if (is_open) close(self%grms%history_unit)
      endif
      if (self%magnetic_field_at_center_domain%history_unit > 0_I4P) then
         inquire(unit=self%magnetic_field_at_center_domain%history_unit, opened=is_open)
         if (is_open) close(self%magnetic_field_at_center_domain%history_unit)
      endif
      inquire(unit=self%io%divergence_history_unit, opened=is_open)
      if (is_open) close(self%io%divergence_history_unit)
   endif
   call self%destroy()
   endsubroutine finalize_forest

   subroutine finalize_mpi_forest(self)
   !< Finalize the process-global FNL MPI handler (mpih_fnl).
   !<
   !< Overrides realm_object%finalize_mpi_forest to target the FNL `mpih_fnl`
   !< singleton instead of the CPU `mpih`. Called ONCE by forest%finalize
   !< after every realm's finalize_forest (see the base for the rationale).
   class(prism_fnl_object), intent(inout) :: self !< The realm (carries no MPI state).

   call mpih_fnl%finalize
   endsubroutine finalize_mpi_forest

   subroutine simulate(self, filename)
   !< Perform the simulation: legacy single-realm entry point.
   class(prism_fnl_object), intent(inout) :: self      !< The equation.
   character(*),            intent(in)    :: filename  !< Input file name.
   logical                                :: loop_done !< Termination predicate.

   call self%initialize_forest(filename=filename)
   integration: do
      call self%compute_dt
      call self%advance_one_step_forest(dt=self%time%dt)
      call self%post_step_forest(dt=self%time%dt, t=self%time%time, it=self%time%it)
      call self%is_done_forest(done=loop_done)
      if (loop_done) exit integration
   enddo integration
   call self%finalize_forest
   call self%finalize_mpi_forest ! finalize_forest no longer finalizes MPI; single-realm legacy path
   endsubroutine simulate

   ! forest orchestrator contract methods overridings
   subroutine initialize_forest(self, filename, realms_number, memory_avail, nv, verbose)
   !< Initialize this realm from scratch: PRISM init, IC injection (or restart load), initial ghost update,
   !< initial diagnostics dump, IO files open, plus PIC/leapfrog priming if those schemes are active.
   !<
   !< Invoked by forest%initialize. The forest writes `self%realm_index = is`
   !< BEFORE calling this routine, so the body can already read the 1-based
   !< forest position from `self%realm_index` if needed.
   class(prism_fnl_object), intent(inout)           :: self          !< The realm.
   character(*),            intent(in)              :: filename      !< Input parameters file name.
   integer(I4P),            intent(in),    optional :: realms_number !< Realm count; divides the per-device budget.
   real(R8P),               intent(in),    optional :: memory_avail  !< Per-process memory budget override.
   integer(I4P),            intent(in),    optional :: nv            !< Number of field variables override.
   logical,                 intent(in),    optional :: verbose       !< Trigger verbose output.
   integer(I4P)                                     :: i             !< Counter.
   integer(I4P)                                     :: n             !< Coil counter.
   integer(I4P)                                     :: b             !< Block counter.
   integer(I4P)                                     :: ind           !< Charge-density variable index for PIC diagnostics.

   call self%initialize_prism(filename=filename, realms_number=realms_number)
   if (self%io%restart) then
      call mpih_fnl%print_message('restart simulation from "'//trim(self%io%restart_basename)//'" files')
      call self%load_restart_files(t=self%time%it, time=self%time%time)
      call mpih_fnl%print_message('restart [t, time]: '//trim(str(self%time%it))//', '//trim(str(self%time%time)))
      call self%set_initial_conditions(is_restart=self%io%restart)
   else
      call mpih_fnl%print_message('impose initial conditions start')
      do i=1, self%ic%amr_iterations
         call mpih_fnl%print_message('  AMR/set IC iteration:'//trim(str(i,.true.)))
         call self%set_initial_conditions(is_restart=self%io%restart)
         call self%amr_update ! host-side (promoted to prism_common, issue #22 F0); device resync below
      enddo
      call self%set_initial_conditions(is_restart=self%io%restart)
      call self%adam%make_comm_local_maps_ghost_bc
      ! Device topology resync (issue #22 F0/GA2): amr_update rebuilt the host tree/field/maps and
      ! make_comm_local_maps_ghost_bc above rebuilt the ghost maps AFTER the last set_initial_conditions'
      ! copy_cpu_gpu push — re-push everything (q, coils, fWL, field coords/dxyz AND the maps) so the
      ! device sees the final refined topology. Idempotent (dev_assign_to_device reallocates dst).
      ! Verbose: this is the copy that establishes the FINAL device topology — the printed map row
      ! counts (incl. seam flag-4 rows, issue #22 F3) are the record of what the kernels will consume.
      call self%copy_cpu_gpu(verbose=.true.)
      self%time%time = 0._R8P
      self%time%it = 0
      call mpih_fnl%print_message('impose initial conditions finish')
   endif

   associate(hs => self%fdv_half_stencils(1))
   call self%compute_divergence(hs=hs, ivar=1_I4P, q=self%q(VAR_DX:VAR_DZ,:,:,:,:), divergence=self%divergence(1,:,:,:,:))
   call self%compute_divergence(hs=hs, ivar=1_I4P, q=self%q(VAR_BX:VAR_BZ,:,:,:,:), divergence=self%divergence(2,:,:,:,:))
   endassociate
   call mpih_fnl%print_message('Initial conditions setting completed')
   if (self%physics%physical_model == EM_PHYSICAL_MODEL .or. self%physics%physical_model == ADIM_EM_PHYSICAL_MODEL) then
      call mpih_fnl%print_message('   max div(D) at t0='//trim(str(maxval(abs(self%divergence(1,:,:,:,:))))))
   elseif (self%physics%physical_model == PIC_PHYSICAL_MODEL) then
      ind = size(self%q(:,1,1,1,1))
      call mpih_fnl%print_message('   max div(D)-rho at t0='//trim(str(maxval(abs(self%divergence(1,:,:,:,:)-self%q(ind,:,:,:,:))))))
   endif
   call mpih_fnl%print_message('   max div(B) at t0='//trim(str(maxval(abs(self%divergence(2,:,:,:,:))))))

   call self%update_ghost(q_gpu=self%q_gpu)

   call mpih_fnl%print_message('Coils initialization values')
   do n=1, self%coil%total_coils_number
      call self%compute_divergence(hs=self%fdv_half_stencils(1), ivar=1_I4P, q=self%coil%J_vec(1:3,:,:,:,:,n), &
                                    divergence=self%divergence(3,:,:,:,:))
      call mpih_fnl%print_message('Coil n='//trim(str(n,.true.)))
      call mpih_fnl%print_message('   max div(j_vec)='//trim(str(maxval(abs(self%divergence(3,:,:,:,:))))))
   enddo

   call mpih_fnl%print_message('assigned block number: '//trim(str(self%adam%field%blocks_number,.true.)))
   do b = 1, self%adam%field%blocks_number
      call mpih_fnl%print_message('  b='//trim(str(b,.true.))//' code='//trim(str(self%adam%field%code(b))))
   enddo

   if (.not.self%io%restart .and. self%pic%problem_type == SINGLE_PARTICLE_TYPE_PROBLEM) then
      call initialize_single_particle_output(filename='single_particle_output.dat')
   endif

   call self%copy_gpu_cpu
   associate(hs => self%fdv_half_stencils(1))
   call self%compute_divergence(hs=hs, ivar=1_I4P, q=self%q, divergence=self%divergence(1,:,:,:,:))
   call self%compute_divergence(hs=hs, ivar=4_I4P, q=self%q, divergence=self%divergence(2,:,:,:,:))
   call self%compute_divergence(hs=hs, ivar=self%physics%var_Jx, q=self%q, divergence=self%divergence(3,:,:,:,:))
   endassociate

   if (self%physics%physical_model == EM_PHYSICAL_MODEL .or. self%physics%physical_model == ADIM_EM_PHYSICAL_MODEL) then
      call mpih_fnl%print_message('   max div(D) at t0 after update_ghost='//trim(str(maxval(abs(self%divergence(1,:,:,:,:))))))
   elseif (self%physics%physical_model == PIC_PHYSICAL_MODEL) then
      ind = size(self%q(:,1,1,1,1))
      call mpih_fnl%print_message('   max div(D)-rho at t0 after update ghost='// &
                                  trim(str(maxval(abs(self%divergence(1,:,:,:,:)-self%q(ind,:,:,:,:))))))
   endif
   call mpih_fnl%print_message('   max div(B) at t0 after update_ghost='//trim(str(maxval(abs(self%divergence(2,:,:,:,:))))))
   call self%save_simulation_data
   call self%compute_energy
   if (self%grms%do_save_history) call self%compute_grms
   if (self%magnetic_field_at_center_domain%do_save_history) call self%compute_magnetic_field_at_center_domain
   !call self%save_energy_error(is_to_open=.true.)
   call self%save_energy_history(is_to_open=.true.) !Cazzo
   call self%save_grms_history(is_to_open=.true.)
   call self%save_magnetic_field_at_center_domain_history(is_to_open=.true.)
   call self%compute_max_divergence
   ! issue #22 F1: pass the maxima compute_max_divergence just stored — the former locals were never assigned
   call self%save_divergence_history(is_to_open=.true., div_D=self%max_divergence_D, div_B=self%max_divergence_B, &
                                     div_J=self%max_divergence_J)
   call self%io%open_file_residuals(nv=self%nv)

   if (self%numerics%scheme_time==NUM_SCHEME_TIME_LEAPFROG) then
      ! to be implemented leapfrog on device
      ! call self%leapfrog%assign_step(s=1, q=self%q)
      ! call self%compute_dt
      ! call self%compute_residuals_dev(q=self%q, dq=self%dq)
      ! self%q = self%q + self%time%dt * self%dq
   endif

   if (self%physics%physical_model == PIC_PHYSICAL_MODEL) then
      if (self%pic%scheme_time==NUM_SCHEME_TIME_PIC_LEAPFROG) then
         ! to be implemented
      endif
   endif

   ! Lock AMR (issue #22 GA6): the FNL device state cannot follow a regrid past this point —
   ! any later amr_update call error_stops instead of computing on stale device maps.
   self%amr_locked_ = .true.
   endsubroutine initialize_forest

   subroutine compute_local_dt_forest(self, dt_local)
   !< Compute this realm's local stability-limited dt (no MPI reduction).
   !<
   !< Invoked by forest%compute_global_dt during the min reduction across all realms in the forest. The reduction
   !< itself is the orchestrator's job; this method computes only the value local to `self`.
   !< Compute this realm's local stability-limited dt (no MPI reduction).
   class(prism_fnl_object), intent(in)  :: self     !< The realm.
   real(R8P),               intent(out) :: dt_local !< Local stability-limited dt.
   real(R8P)                            :: dxyz_min !< Minimum space step.

   call compute_dxyz_min_kernel(blocks_number=self%blocks_number, dxyz_gpu=self%field_fnl%dxyz_gpu, dxyz_min=dxyz_min)
   dxyz_min = dxyz_min * 0.5_R8P
   dt_local = self%time%CFL*dxyz_min / self%physics%evmax
   contains
      subroutine compute_dxyz_min_kernel(blocks_number, dxyz_gpu, dxyz_min)
      !< Compute minimum space step accordingly, kernel device.
      integer(I4P), intent(in)  :: blocks_number   !< Blocks number.
      real(R8P),    intent(in)  :: dxyz_gpu(1:,1:) !< Delta cells GPU [nb,3].
      real(R8P),    intent(out) :: dxyz_min        !< Minimum space step.
      integer(I4P)              :: b               !< Counter.

      dxyz_min = huge(0._R8P)
      !$acc parallel loop independent gang vector DEVICEVAR(dxyz_gpu) firstprivate (blocks_number) reduction(min: dxyz_min)
      !$omp OMPLOOP DEVICEPTR(dxyz_gpu) firstprivate(blocks_number) reduction(min: dxyz_min)
      do b=1, blocks_number
         dxyz_min = min(dxyz_min, dxyz_gpu(b,1), dxyz_gpu(b,2), dxyz_gpu(b,3))
      enddo
      endsubroutine compute_dxyz_min_kernel
   endsubroutine compute_local_dt_forest

   subroutine advance_one_step_forest(self, dt)
   !< Advance this realm by one full timestep of size `dt`.
   !<
   !< Invoked by forest%evolve_one_step once per realm per timestep. Owns the integration itself, i.e. everything that turns
   !< `q` at time `t` into `q` at time `t + dt`.
   class(prism_fnl_object), intent(inout) :: self    !< The realm.
   real(R8P),               intent(in)    :: dt      !< Timestep size from the forest's global reduction.
   real(R8P)                              :: dt_step !< Local copy, possibly capped for time_max.

   self%time%it = self%time%it + 1
   dt_step = dt
   if ((self%time%it_max <= 0).and.(self%time%time+dt_step > self%time%time_max)) dt_step = self%time%time_max - self%time%time
   self%time%dt = dt_step
   call self%integrate_dev
   self%time%time = self%time%time + dt_step
   call self%time%print_progress(nodes_number=self%adam%tree%nodes_number)
   endsubroutine advance_one_step_forest

   function stages_per_step_forest(self) result(K)
   !< Number of integrator stages this realm exposes per step (FNL).
   !<
   !< SSP-only contract, twin of `prism_cpu_object%stages_per_step_forest`
   !< (issue #25): the staged protocol reads `gamm(k)` per stage and
   !< `beta_gpu` in `close_step_forest` — for low-storage schemes those are
   !< never allocated (the FNL symptom was CUDA_ERROR_ILLEGAL_ADDRESS in
   !< `rk_update_q_dev` through a bogus `beta_gpu`). Refusing here leaves
   !< the fused fast path — where LS schemes legitimately run — untouched.
   class(prism_fnl_object), intent(in) :: self !< The realm.
   integer(I4P)                        :: K    !< Number of integrator stages per step.

   select case(self%rk%scheme)
   case(RK_SSP_11, RK_SSP_22, RK_SSP_33, RK_SSP_54)
      K = self%rk%nrk
   case default
      K = 0
      ! Routed through the ADAM (CPU) mpih handler, NOT mpih_fnl: FUNDAL's error_stop
      ! ends with a plain `stop` (exit code 0), so a refusal through it looks SUCCESSFUL
      ! to the calling shell/harness (upstream FUNDAL defect, flagged in issue #25).
      ! This is host code; both handlers wrap the same communicator.
      call mpih%error_stop(msg=': RK scheme "'//trim(adjustl(self%rk%scheme))//'" is not stage-splittable: '// &
                               'the staged forest path (multi-realm, or intra-realm AMR seam faces) requires '// &
                               'an SSP scheme (runge-kutta-ssp-*); low-storage schemes run only on the fused '// &
                               'single-realm/no-seam fast path')
   endselect
   endfunction stages_per_step_forest

   subroutine open_step_forest(self, dt)
   !< Per-step prologue on the multi-realm path: set dt, init RK stages.
   !<
   !< Mirrors the head of `integrate_rk_ssp`: external-field prelude (if active), `rk%initialize_stages(q=self%q)`, and
   !< the per-step time bookkeeping `advance_one_step_forest` does inline (it increment, dt cap for time_max, self%time%dt update).
   !< Time advance and progress print run in `close_step_forest`, mirroring the legacy ordering.
   class(prism_fnl_object), intent(inout) :: self    !< The realm.
   real(R8P),               intent(in)    :: dt      !< Timestep size from the forest.
   real(R8P)                              :: dt_step !< Local copy, possibly capped for time_max.

   self%time%it = self%time%it + 1
   dt_step = dt
   if ((self%time%it_max <= 0).and.(self%time%time+dt_step > self%time%time_max)) dt_step = self%time%time_max - self%time%time
   self%time%dt = dt_step
   if (self%external_fields%ef_type/=EF_TYPE_NONE) &
      call sub_external_fields_dev(external_fields=self%external_fields, field_gpu=self%field_fnl, &
                                   dt=self%time%dt, time=self%time%time, q_gpu=self%q_gpu)
   call self%rk_fnl%initialize_stages(grid=self%adam%grid, field=self%adam%field, q_gpu=self%q_gpu)
   if (self%pml_fnl%enabled) call self%rk_pml_fnl%initialize_stages(pml_fnl=self%pml_fnl)
   endsubroutine open_step_forest

   subroutine begin_stage_forest(self, k, K_total, dt, realm)
   !< Begin integrator stage `k` on the multi-realm path (FNL).
   !<
   !< Mirrors the FIRST line of `integrate_rk_ssp`'s substage body:
   !< `rk_fnl%compute_stage(s=k, dt=self%time%dt [, phi_gpu=ib_fnl%phi_gpu])`.
   !< This populates `rk_fnl%q_rk_gpu(:,:,:,:,:,k)` from previously computed
   !< substages and from `self%q_gpu`, and publishes the stage buffer to
   !< peers by setting `self%stage_active = k`. No ghost reads, no
   !< peer-realm access — peer realms may not yet have opened their
   !< stage-`k` buffer when this fires.
   class(prism_fnl_object), intent(inout)                   :: self     !< The realm.
   integer(I4P),            intent(in)                      :: k        !< Stage index (1..K_total).
   integer(I4P),            intent(in)                      :: K_total  !< Forest-wide stage count for this step.
   real(R8P),               intent(in)                      :: dt       !< Timestep size from the forest.
   class(realm_object),     intent(inout), optional, target :: realm(:) !< Sibling realms (contract parity).

   self%stage_active = k
   if (self%ib%solids_number>0) then
      call self%rk_fnl%compute_stage(grid=self%adam%grid, field=self%adam%field, s=k, dt=self%time%dt, phi_gpu=self%ib_fnl%phi_gpu)
   else
      call self%rk_fnl%compute_stage(grid=self%adam%grid, field=self%adam%field, s=k, dt=self%time%dt)
   endif
   if (self%pml_fnl%enabled) call self%rk_pml_fnl%compute_stage(s=k, dt=self%time%dt)
   call self%compute_coils_current(q_gpu=self%rk_fnl%q_rk_gpu(:,:,:,:,:,k), gamm=self%rk%gamm(k))
   endsubroutine begin_stage_forest

   subroutine end_stage_forest(self, k, K_total, dt, realm, flux_register)
   !< Close integrator stage `k`: residuals + stage assignment in one sweep (FNL).
   !<
   !< This is the FNL forest's Phase 3 TBP. Inter-realm seam ghost cells
   !< on `rk_fnl%q_rk_gpu(:,...,k)` are already filled by Phase 2 (which
   !< called `fill_seam_from_peer_forest` on this realm and its
   !< peers), so `compute_residuals_dev` reads valid halo data and
   !< `rk_fnl%assign_stage`'s in-place overwrite of the stage interior
   !< no longer races with peer reads.
   !<
   !< `realm` accepted on contract for parity but unused: the seam has
   !< already been refreshed by the forest. FV seam-flux ACCUMULATION is
   !< threaded through `flux_register` (issue #23 R3, CPU parity); the
   !< reflux APPLICATION (`apply_reflux_to_stage_forest`) remains a no-op
   !< until #23 R4.
   class(prism_fnl_object),     intent(inout)                   :: self          !< The realm.
   integer(I4P),                intent(in)                      :: k             !< Stage index (1..K_total).
   integer(I4P),                intent(in)                      :: K_total       !< Forest-wide stage count for this step.
   real(R8P),                   intent(in)                      :: dt            !< Timestep size from the forest.
   class(realm_object),         intent(inout), optional, target :: realm(:)      !< Sibling realms (FNL parity only; unused).
   class(flux_register_object), intent(inout), optional         :: flux_register !< Forest's flux register for FV reflux.

   if (present(realm)) continue
   if (present(flux_register)) then
      call self%compute_residuals_dev(q_gpu=self%rk_fnl%q_rk_gpu(:,:,:,:,:,k), dq_gpu=self%dq_gpu, s=k, &
                                      flux_register=flux_register)
   else
      call self%compute_residuals_dev(q_gpu=self%rk_fnl%q_rk_gpu(:,:,:,:,:,k), dq_gpu=self%dq_gpu, s=k)
   endif
   if (self%ib%solids_number>0) then
      call self%rk_fnl%assign_stage(grid=self%adam%grid, field=self%adam%field, s=k, q_gpu=self%dq_gpu, phi_gpu=self%ib_fnl%phi_gpu)
   else
      call self%rk_fnl%assign_stage(grid=self%adam%grid, field=self%adam%field, s=k, q_gpu=self%dq_gpu)
   endif
   if (self%pml_fnl%enabled) call self%rk_pml_fnl%assign_stage(s=k)
   endsubroutine end_stage_forest

   subroutine close_step_forest(self, dt)
   !< Per-step epilogue on the multi-realm path: q assembly, BC, div-clean,
   !< residual save, coil source refresh, time advance, progress print (FNL).
   class(prism_fnl_object), intent(inout) :: self !< The realm.
   real(R8P),               intent(in)    :: dt   !< Timestep size from the forest.

   if (self%ib%solids_number>0) then
      call self%rk_fnl%update_q(grid=self%adam%grid, field=self%adam%field, rk=self%rk, dt=self%time%dt, &
                                phi_gpu=self%ib_fnl%phi_gpu, q_gpu=self%q_gpu)
   else
      call self%rk_fnl%update_q(grid=self%adam%grid, field=self%adam%field, rk=self%rk, dt=self%time%dt, q_gpu=self%q_gpu)
      call self%save_residuals
   endif
   if (self%pml_fnl%enabled) call self%rk_pml_fnl%update_q_pml(dt=self%time%dt, pml_fnl=self%pml_fnl)
   call self%apply_fwl_correction(q_gpu=self%q_gpu)
   call self%compute_coils_current(q_gpu=self%q_gpu)
   call self%impose_div_free
   if (self%external_fields%ef_type/=EF_TYPE_NONE) &
      call add_external_fields_dev(external_fields=self%external_fields, field_gpu=self%field_fnl, &
                                   dt=self%time%dt, time=self%time%time, q_gpu=self%q_gpu)
   self%time%time = self%time%time + self%time%dt
   call self%time%print_progress(nodes_number=self%adam%tree%nodes_number)
   ! Clear active-stage marker (mirror of CPU). Inter-step get_cell /
   ! `fill_seam_from_peer_forest` queries read self%q_gpu (committed) instead of
   ! a stale q_rk_gpu slot from the previous step.
   self%stage_active = 0_I4P
   endsubroutine close_step_forest

   ! Inter-realm seam ghost-fill contract — PRISM FNL override.
   !
   ! Same signature as CPU (`fill_seam_from_peer_forest(self, peer, p_idx)`);
   ! only the body differs. Walks the DEVICE-resident seam map
   ! (`maps_fnl%seam_local_map_ghost_cell_gpu`) under an OpenACC
   ! `parallel loop` and copies peer's active GPU buffer cells into self's
   ! GPU ghost cells. No host staging; same-rank peers share device
   ! memory so the cross-realm read is a device-to-device copy on the
   ! same GPU. Cross-rank peers (not yet implemented) will instead route
   ! the same row range through the `seam_comm_map_*` MPI path
   ! (GPU-direct MPI_Isend/MPI_Irecv on the device-resident
   ! `seam_mpi_send/recv_buf_gpu`).

   subroutine fill_seam_from_peer_forest(self, peer, p_idx)
   !< FNL twin of `prism_cpu_object%fill_seam_from_peer_forest`. Walks the
   !< host-side row-range bookkeeping (`seam_local_peer_row_start/count`,
   !< populated by `build_seam_local_map` on the host), then dispatches
   !< an OpenACC kernel that reads peer's GPU buffer at peer-INTERIOR
   !< coords (FNL layout `(b,i,j,k,v)`) and writes self's GPU buffer at
   !< recv-side coords.
   class(prism_fnl_object), intent(inout)         :: self
   class(realm_object),     intent(in),    target :: peer
   integer(I4P),            intent(in)            :: p_idx
   integer(I4P)                                   :: c, v, row
   integer(I4P)                                   :: b_send, i_send, j_send, k_send
   integer(I4P)                                   :: b_recv, i_recv, j_recv, k_recv
   integer(I4P)                                   :: row_start, row_count
   integer(I4P)                                   :: nv_loc, ngc_loc
   integer(I4P)                                   :: self_s, peer_s
   integer(I4P), pointer                          :: rows_gpu(:,:)

   row_start = self%adam%maps%seam_local_peer_row_start(p_idx)
   row_count = self%adam%maps%seam_local_peer_row_count(p_idx)
   rows_gpu => self%field_fnl%maps%seam_local_map_ghost_cell_gpu
   nv_loc   = self%nv
   ngc_loc  = self%ngc
   self_s   = self%stage_active
   select type (peer)
   class is (prism_fnl_object)
      peer_s = peer%stage_active
      if (self_s > 0_I4P .and. peer_s > 0_I4P) then
         associate(self_q_rk_gpu => self%rk_fnl%q_rk_gpu, peer_q_rk_gpu => peer%rk_fnl%q_rk_gpu)
         !$acc parallel loop collapse(2) independent DEVICEVAR(rows_gpu, self_q_rk_gpu, peer_q_rk_gpu) &
         !$acc&                                     firstprivate(ngc_loc, nv_loc, row_start, row_count, self_s, peer_s)
         !$omp OMPLOOP collapse(2) DEVICEPTR(rows_gpu, self_q_rk_gpu, peer_q_rk_gpu) &
         !$omp& firstprivate(ngc_loc, nv_loc, row_start, row_count, self_s, peer_s)
         do c = 1, row_count
            do v = 1, nv_loc
               row    = row_start + c - 1
               b_send = rows_gpu(row, 2) ; i_send = rows_gpu(row, 4)
               j_send = rows_gpu(row, 5) ; k_send = rows_gpu(row, 6)
               b_recv = rows_gpu(row, 3) ; i_recv = rows_gpu(row, 7)
               j_recv = rows_gpu(row, 8) ; k_recv = rows_gpu(row, 9)
               self_q_rk_gpu(b_recv, i_recv, j_recv, k_recv, v, self_s) = &
                  peer_q_rk_gpu(b_send, i_send, j_send, k_send, v, peer_s)
            enddo
         enddo
         endassociate
      else if (self_s > 0_I4P) then
         associate(self_q_rk_gpu => self%rk_fnl%q_rk_gpu, peer_q_gpu => peer%q_gpu)
         !$acc parallel loop collapse(2) independent DEVICEVAR(rows_gpu, self_q_rk_gpu, peer_q_gpu) &
         !$acc&                                     firstprivate(ngc_loc, nv_loc, row_start, row_count, self_s)
         !$omp OMPLOOP collapse(2) DEVICEPTR(rows_gpu, self_q_rk_gpu, peer_q_gpu) &
         !$omp& firstprivate(ngc_loc, nv_loc, row_start, row_count, self_s)
         do c = 1, row_count
            do v = 1, nv_loc
               row    = row_start + c - 1
               b_send = rows_gpu(row, 2) ; i_send = rows_gpu(row, 4)
               j_send = rows_gpu(row, 5) ; k_send = rows_gpu(row, 6)
               b_recv = rows_gpu(row, 3) ; i_recv = rows_gpu(row, 7)
               j_recv = rows_gpu(row, 8) ; k_recv = rows_gpu(row, 9)
               self_q_rk_gpu(b_recv, i_recv, j_recv, k_recv, v, self_s) = &
                  peer_q_gpu(b_send, i_send, j_send, k_send, v)
            enddo
         enddo
         endassociate
      else if (peer_s > 0_I4P) then
         associate(self_q_gpu => self%q_gpu, peer_q_rk_gpu => peer%rk_fnl%q_rk_gpu)
         !$acc parallel loop collapse(2) independent DEVICEVAR(rows_gpu, self_q_gpu, peer_q_rk_gpu) &
         !$acc&                                     firstprivate(ngc_loc, nv_loc, row_start, row_count, peer_s)
         !$omp OMPLOOP collapse(2) DEVICEPTR(rows_gpu, self_q_gpu, peer_q_rk_gpu) &
         !$omp& firstprivate(ngc_loc, nv_loc, row_start, row_count, peer_s)
         do c = 1, row_count
            do v = 1, nv_loc
               row    = row_start + c - 1
               b_send = rows_gpu(row, 2) ; i_send = rows_gpu(row, 4)
               j_send = rows_gpu(row, 5) ; k_send = rows_gpu(row, 6)
               b_recv = rows_gpu(row, 3) ; i_recv = rows_gpu(row, 7)
               j_recv = rows_gpu(row, 8) ; k_recv = rows_gpu(row, 9)
               self_q_gpu(b_recv, i_recv, j_recv, k_recv, v) = &
                  peer_q_rk_gpu(b_send, i_send, j_send, k_send, v, peer_s)
            enddo
         enddo
         endassociate
      else
         associate(self_q_gpu => self%q_gpu, peer_q_gpu => peer%q_gpu)
         !$acc parallel loop collapse(2) independent DEVICEVAR(rows_gpu, self_q_gpu, peer_q_gpu) &
         !$acc&                                     firstprivate(ngc_loc, nv_loc, row_start, row_count)
         !$omp OMPLOOP collapse(2) DEVICEPTR(rows_gpu, self_q_gpu, peer_q_gpu) &
         !$omp& firstprivate(ngc_loc, nv_loc, row_start, row_count)
         do c = 1, row_count
            do v = 1, nv_loc
               row    = row_start + c - 1
               b_send = rows_gpu(row, 2) ; i_send = rows_gpu(row, 4)
               j_send = rows_gpu(row, 5) ; k_send = rows_gpu(row, 6)
               b_recv = rows_gpu(row, 3) ; i_recv = rows_gpu(row, 7)
               j_recv = rows_gpu(row, 8) ; k_recv = rows_gpu(row, 9)
               self_q_gpu(b_recv, i_recv, j_recv, k_recv, v) = &
                  peer_q_gpu(b_send, i_send, j_send, k_send, v)
            enddo
         enddo
         endassociate
      endif
   class default
      call mpih_fnl%error_stop(msg='prism_fnl_object%fill_seam_from_peer_forest: peer realm is not prism_fnl_object')
   end select
   endsubroutine fill_seam_from_peer_forest

   subroutine after_topology_build_forest(self)
   !< Propagate the freshly-built host seam maps and per-peer buffers to
   !< device-resident counterparts on `field_fnl%maps`. Invoked by the
   !< forest at the end of `populate_inter_realm_topology` AFTER
   !< `build_seam_local_map` has populated `self%adam%maps%seam_local_*`.
   class(prism_fnl_object), intent(inout) :: self

   call self%field_fnl%maps%copy_cpu_gpu(maps=self%adam%maps)
   endsubroutine after_topology_build_forest

   subroutine apply_reflux_to_stage_forest(self, stage, dt, flux_register)
   !< PRISM-FNL override of the Berger-Colella reflux correction TBP —
   !< the apply twin of `prism_cpu_object%apply_reflux_to_stage_forest`
   !< (issue #23 R4, M4 semantics).
   !<
   !< **α.r1 cadence: end-of-step barrier.** Fires exactly once per realm
   !< per step at `stage == self%rk%nrk`; earlier substages return
   !< immediately. For each register face whose coarse side this realm —
   !< and this RANK (#28 D4) — owns, the correction
   !< `sign*(dt/dx_coarse)*(F_coarse − F_fine_sum)(:,:,1)` is applied to
   !< the COMMITTED `q_gpu` on the one-cell-thick coarse seam skin (this
   !< TBP runs AFTER `close_step_forest`'s `update_q`). Full step weight
   !< `dt/dx` — NOT a stage RK coefficient, and NOT a `q_rk_gpu(...,stage)`
   !< write (both were the pre-M4 sketch's errors: `ark` is never
   !< allocated for the SSP family, and a stage-buffer write would
   !< entangle the stage beta weight — see the CPU apply's M4 note).
   !< The tiny host mismatch slab is H2D-copied per face and
   !< `fv_apply_reflux_face_dev_kernel` adds it to `q_gpu` (scalar ops
   !< only, disjoint cells).
   class(prism_fnl_object),     intent(inout) :: self          !< The realm.
   integer(I4P),                intent(in)    :: stage         !< Integrator stage 1..K_total.
   real(R8P),                   intent(in)    :: dt            !< Time step.
   class(flux_register_object), intent(in)    :: flux_register !< Forest's flux register.
   integer(I4P)                               :: f             !< Register face counter.
   integer(I4P)                               :: axis, sgn     !< Face normal axis and sign.
   integer(I4P)                               :: nv_reg        !< Register state-vector width.
   integer(I4P)                               :: nface_cells   !< Coarse-face skin cell count.
   integer(I4P)                               :: ierr          !< Device allocation error flag.
   real(R8P)                                  :: dx_coarse     !< Coarse spacing along the face normal.
   real(R8P)                                  :: scale_        !< Correction scale: sign * dt / dx_coarse.
   real(R8P), allocatable                     :: delta(:,:)    !< Host flux mismatch F_coarse - F_fine_sum (:,:,1).
   real(R8P), pointer                         :: delta_gpu(:,:)!< Device copy of the mismatch slab.

   ! α.r1 end-of-step gate: parity with PRISM-CPU M3. Returns immediately
   ! for any non-final substage.
   if (stage /= self%rk%nrk) return
   if (.not. flux_register%is_initialized_) return
   if (flux_register%nfaces == 0_I4P)       return
   if (.not. allocated(flux_register%face)) return

   ! Issue #23 R4 — the correction application, M4 semantics (CPU parity): for each
   ! face this realm owns the coarse side of, apply sign*(dt/dx_coarse)*(F_coarse -
   ! F_fine_sum)(:,:,1) to the COMMITTED q_gpu on the one-cell-thick coarse seam
   ! skin, once per realm per step (this TBP runs AFTER close_step_forest's
   ! update_q). Full step weight dt/dx — NOT a stage RK coefficient (see the CPU
   ! apply's M4 note). The tiny host mismatch slab is H2D-copied and a per-face
   ! device kernel adds it to q_gpu (scalar ops only, disjoint cells).
   do f=1_I4P, flux_register%nfaces
      associate(face_f => flux_register%face(f))
      if (face_f%coarse_realm == self%realm_index) then
      ! Issue #28 D4 (CPU parity): only the coarse block's OWNER rank applies
      ! (and prints) this face — `coarse_block` is an owner-rank-LOCAL field
      ! slot, an aliased unrelated block on any other rank. Accumulators are
      ! complete on every rank after reduce_fine_sums (#28 D3).
      if (face_f%coarse_rank == mpih%myrank) then
      if (allocated(face_f%F_coarse) .and. allocated(face_f%F_fine_sum)) then
         call face_axis_sign(face_f%coarse_face, axis, sgn)
         if (axis /= 0_I4P) then
            dx_coarse = self%adam%field%dxyz(axis, face_f%coarse_block)
            if (dx_coarse > 0._R8P) then
               ! Register-level diagnostic (issue #23 R3): format matched with the CPU
               ! apply — the two backends' register contents stay log-comparable.
               call mpih_fnl%print_message('reflux face '//trim(str(f, .true.))//' coarse_block '//                 &
                                           trim(str(face_f%coarse_block, .true.))//' max|F_coarse-F_fine_sum| = '// &
                                           trim(str(maxval(abs(face_f%F_coarse(:,:,1) - face_f%F_fine_sum(:,:,1))))))
               scale_      = real(sgn, R8P) * dt / dx_coarse
               nv_reg      = int(size(face_f%F_coarse, dim=1), I4P)
               nface_cells = face_f%nface_cells
               allocate(delta(1:nv_reg, 1:nface_cells))
               delta = face_f%F_coarse(:,:,1) - face_f%F_fine_sum(:,:,1)
               call dev_alloc(fptr_dev=delta_gpu, lbounds=[1,1], ubounds=[nv_reg,nface_cells], ierr=ierr)
               if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate delta_gpu in apply_reflux_to_stage_forest')
               call dev_memcpy_to_device(dst=delta_gpu, src=delta)
               call fv_apply_reflux_face_dev_kernel(axis=axis, sgn=sgn, b=face_f%coarse_block,          &
                                                    ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc,   &
                                                    nv_reg=nv_reg, nface_cells=nface_cells,             &
                                                    scale=scale_, delta_gpu=delta_gpu, q_gpu=self%q_gpu)
               call dev_free(delta_gpu, mydev)
               deallocate(delta)
            endif
         endif
      endif
      endif
      endif
      endassociate
   enddo
   endsubroutine apply_reflux_to_stage_forest

   subroutine post_step_forest(self, dt, t, it, do_save_state, do_save_residuals, do_save_restart, do_amr, realm)
   !< Run PRISM-FNL's per-timestep post-step work: state IO, energy
   !< diagnostics, max-divergence diagnostics.
   !<
   !< Invoked by forest%post_step. v1 implementation is the verbatim post-
   !< step block formerly inline in `simulate` — every action runs every
   !< step, since today's cadence is enforced inside the save_* routines
   !< themselves (e.g. save_simulation_data honours `io%it_save`). The
   !< `do_*` flags are signature-only for now: when the forest takes
   !< over cadence the flags will gate the individual calls. For now
   !< they are accepted but unused, preserving present-day behavior
   !< bit-for-bit.
   !<
   !< `dt`, `t`, `it` are not consumed by the current body; they are on
   !< the contract so the forest can supply them once it owns time-state
   !< (today they are still read from the `time` module singleton).
   !<
   !< Optional `realm(:)`: forwarded to save_simulation_data and
   !< update_ghost for the dummy-argument inter-realm halo refresh path.
   class(prism_fnl_object), intent(inout)                   :: self              !< The realm.
   real(R8P),               intent(in)                      :: dt                !< Timestep size just advanced.
   real(R8P),               intent(in)                      :: t                 !< Simulation time after the advance.
   integer(I4P),            intent(in)                      :: it                !< Iteration index after the advance.
   logical,                 intent(in),    optional         :: do_save_state     !< Save state output this step.
   logical,                 intent(in),    optional         :: do_save_residuals !< Save residuals output this step.
   logical,                 intent(in),    optional         :: do_save_restart   !< Save restart dump this step.
   logical,                 intent(in),    optional         :: do_amr            !< Run AMR update this step.
   class(realm_object),     intent(inout), optional, target :: realm(:)          !< Sibling realms for inter-realm halo refresh.

   if (self%io%save_memory_status) then
      call save_memory_status_cpu(file_name='memory_cpu-'//mpih_fnl%myrankstr//'.dat', tag=str(self%time%it,.true.))
      call save_memory_status_gpu(file_name='memory_gpu-'//mpih_fnl%myrankstr//'.dat', tag=str(self%time%it,.true.))
   endif
   if (mod(self%time%it,self%amr%frequency)==0) then
      call mpih_fnl%barrier(tictoc=.true.)
      !call self%amr_update
      call mpih_fnl%barrier(tictoc=.true.)
   endif
   call self%save_simulation_data
   call self%update_ghost(q_gpu=self%q_gpu) !Cazzo
   ! Issue #31 (CPU-parity): self%update_ghost fills the intra-realm ghosts + physical BCs
   ! but NOT the inter-realm seam (the forest owns that). The max-divergence diagnostic
   ! below reads an s1-deep stencil that, at seam-adjacent cells, needs the seam ghosts —
   ! which update_ghost leaves stale/BC-filled, producing a spurious div(B) at the seam skin
   ! (the evolved field is div-free). Re-establish the inter-realm seam from peers on the
   ! committed q_gpu (stage_active==0 → fill writes self%q_gpu) before the diagnostic.
   ! Mirrors the CPU fix in prism_cpu_object%post_step_forest.
   if (present(realm) .and. allocated(self%adam%maps%seam_local_map_ghost_cell) .and. &
       allocated(self%adam%maps%seam_local_peer_realm)) then
      block
         integer(I4P) :: p_s
         do p_s = 1_I4P, int(size(self%adam%maps%seam_local_peer_realm), I4P)
            call self%fill_seam_from_peer_forest(peer=realm(self%adam%maps%seam_local_peer_realm(p_s)), p_idx=p_s)
         enddo
      endblock
   endif
   call self%compute_energy
   if (self%grms%do_save_history) call self%compute_grms
   if (self%magnetic_field_at_center_domain%do_save_history) call self%compute_magnetic_field_at_center_domain
   !call self%save_energy_error !Cazzo
   call self%save_energy_history !Cazzo
   call self%save_grms_history
   call self%save_magnetic_field_at_center_domain_history
   call self%compute_max_divergence
   ! issue #22 F1: pass the maxima compute_max_divergence just stored — the former locals were never assigned
   call self%save_divergence_history(div_D=self%max_divergence_D, div_B=self%max_divergence_B, &
                                     div_J=self%max_divergence_J)
   endsubroutine post_step_forest

   ! numerical methods
   subroutine compute_dt(self)
   !< Compute the global stability-limited dt and store it on `self%time%dt`.
   !<
   !< Body delegates the local computation to compute_local_dt_forest
   !< (orchestrator contract method), then performs the legacy
   !< MPI_ALLREDUCE on MPI_COMM_WORLD for backward compatibility with
   !< simulate. The forest's `compute_global_dt` performs its own
   !< reduction, possibly on a per-realm sub-comm; the redundancy
   !< disappears once the legacy compute_dt is retired.
   class(prism_fnl_object), intent(inout) :: self !< The equation.

   call self%compute_local_dt_forest(dt_local=self%time%dt)
   call MPI_ALLREDUCE(MPI_IN_PLACE, self%time%dt, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, mpih_fnl%error)
   endsubroutine compute_dt

   subroutine compute_energy(self)
   !< Compute energy.
   class(prism_fnl_object), intent(inout) :: self          !< The equation.
   real(R8P)                              :: energy_D      !< Energy of D field.
   real(R8P)                              :: energy_B      !< Energy of B field.
   real(R8P)                              :: coil_power    !< Coil power.
   real(R8P)                              :: poynting_flux !< Total Poynting flux from boundary.

   call compute_e_dev_kernel(ni=self%ni,nj=self%nj,nk=self%nk,ngc=self%ngc,blocks_number=self%blocks_number,&
                             ivar=VAR_DX,dxyz_gpu=self%field_fnl%dxyz_gpu,q_gpu=self%q_gpu,energy=energy_D)
   call compute_e_dev_kernel(ni=self%ni,nj=self%nj,nk=self%nk,ngc=self%ngc,blocks_number=self%blocks_number,&
                             ivar=VAR_BX,dxyz_gpu=self%field_fnl%dxyz_gpu,q_gpu=self%q_gpu,energy=energy_B)
   if (self%coil%total_coils_number > 0_I4P) then
      call compute_coil_power_dev_kernel(ni=self%ni,nj=self%nj,nk=self%nk,ngc=self%ngc,blocks_number=self%blocks_number,&
                                         ivar=self%physics%var_Jx,            &
                                         dxyz_gpu=self%field_fnl%dxyz_gpu,q_gpu=self%q_gpu,coil_power=coil_power)
   else
      coil_power = 0._R8P
   endif
   call compute_poynting_flux_dev_kernel(ni=self%ni,nj=self%nj,nk=self%nk,ngc=self%ngc,blocks_number=self%blocks_number,&
                                         s=self%fdv_half_stencils(1),                                                   &
                                         dxyz_gpu=self%field_fnl%dxyz_gpu,q_gpu=self%q_gpu,poynting_flux=poynting_flux)
   call MPI_ALLREDUCE(MPI_IN_PLACE, energy_D,      1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih_fnl%error)
   call MPI_ALLREDUCE(MPI_IN_PLACE, energy_B,      1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih_fnl%error)
   call MPI_ALLREDUCE(MPI_IN_PLACE, coil_power,    1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih_fnl%error)
   call MPI_ALLREDUCE(MPI_IN_PLACE, poynting_flux, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih_fnl%error)
   if (allocated(self%energy_D  ).and.allocated(self%energy_B     ).and. &
       allocated(self%coil_power).and.allocated(self%Poynting_flux)) then
      self%energy_D      = [self%energy_D,      energy_D     ]
      self%energy_B      = [self%energy_B,      energy_B     ]
      self%coil_power    = [self%coil_power,    coil_power   ]
      self%poynting_flux = [self%poynting_flux, poynting_flux]
   else
      allocate(self%energy_D(     1:self%time%it))
      allocate(self%energy_B(     1:self%time%it))
      allocate(self%coil_power(   1:self%time%it))
      allocate(self%poynting_flux(1:self%time%it))
      self%energy_D      = energy_D
      self%energy_B      = energy_B
      self%coil_power    = coil_power
      self%poynting_flux = poynting_flux
   endif
   contains
      subroutine compute_e_dev_kernel(ni,nj,nk,ngc,blocks_number,ivar,dxyz_gpu,q_gpu,energy)
      !< Compute energy of vector field starting from ivar, device kernel.
      integer(I4P), intent(in)  :: ni,nj,nk,ngc,blocks_number        !< Grids dimensions.
      integer(I4P), intent(in)  :: ivar                              !< Starting position of vector field.
      real(R8P),    intent(in)  :: dxyz_gpu(1:,1:)                   !< Delta cells GPU [nb,3].
      real(R8P),    intent(in)  :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
      real(R8P),    intent(out) :: energy                            !< Energy of the vector field starting from ivar.
      real(R8P)                 :: const                             !< Costant for the energy computation.
      integer(I4P)              :: i,j,k,b                           !< Counter.

      if (ivar==VAR_DX) then
         const = EPS0
      elseif (ivar==VAR_BX) then
         const = MU0
      endif
      energy = 0.0_R8P
      !$acc parallel loop independent gang vector collapse(4) &
      !$acc& DEVICEVAR(q_gpu,dxyz_gpu) firstprivate(ni,nj,nk,blocks_number,ivar,const) reduction(+: energy)
      !$omp OMPLOOP collapse(4) &
      !$omp& DEVICEPTR(q_gpu,dxyz_gpu) firstprivate(ni,nj,nk,blocks_number,ivar,const) reduction(+: energy)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         energy = energy + 0.5_R8P * (q_gpu(b,i,j,k,ivar  )*q_gpu(b,i,j,k,ivar  ) + &
                                      q_gpu(b,i,j,k,ivar+1)*q_gpu(b,i,j,k,ivar+1) + &
                                      q_gpu(b,i,j,k,ivar+2)*q_gpu(b,i,j,k,ivar+2))/const*(dxyz_gpu(b,1)*dxyz_gpu(b,2)*dxyz_gpu(b,3))
      enddo
      enddo
      enddo
      enddo
      endsubroutine compute_e_dev_kernel

      subroutine compute_coil_power_dev_kernel(ni,nj,nk,ngc,blocks_number,ivar,dxyz_gpu,q_gpu,coil_power)
      !< Compute coil power of vector field starting from ivar, device kernel.
      integer(I4P), intent(in)  :: ni,nj,nk,ngc,blocks_number        !< Grids dimensions.
      integer(I4P), intent(in)  :: ivar                              !< Starting position of vector field.
      real(R8P),    intent(in)  :: dxyz_gpu(1:,1:)                   !< Delta cells GPU [nb,3].
      real(R8P),    intent(in)  :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
      real(R8P),    intent(out) :: coil_power                        !< Coil power of the vector field.
      integer(I4P)              :: i,j,k,b                           !< Counter.

      coil_power = 0.0_R8P
      !$acc parallel loop independent gang vector collapse(4) &
      !$acc& DEVICEVAR(q_gpu,dxyz_gpu) firstprivate(ni,nj,nk,blocks_number,ivar) reduction(+: coil_power)
      !$omp OMPLOOP collapse(4) &
      !$omp& DEVICEPTR(q_gpu,dxyz_gpu) firstprivate(ni,nj,nk,blocks_number,ivar) reduction(+: coil_power)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         coil_power = coil_power-(q_gpu(b,i,j,k,VAR_DX  )*q_gpu(b,i,j,k,ivar  ) + &
                                  q_gpu(b,i,j,k,VAR_DX+1)*q_gpu(b,i,j,k,ivar+1) + &
                                  q_gpu(b,i,j,k,VAR_DX+2)*q_gpu(b,i,j,k,ivar+2))/EPS0*(dxyz_gpu(b,1)*dxyz_gpu(b,2)*dxyz_gpu(b,3))
      enddo
      enddo
      enddo
      enddo
      endsubroutine compute_coil_power_dev_kernel

      subroutine compute_poynting_flux_dev_kernel(ni,nj,nk,ngc,blocks_number,s,dxyz_gpu,q_gpu,poynting_flux)
      !< Compute Poynting flux, device kernel.
      integer(I4P), intent(in)  :: ni,nj,nk,ngc,blocks_number        !< Grids dimensions.
      integer(I4P), intent(in)  :: s                                 !< FDV half stencil.
      real(R8P),    intent(in)  :: dxyz_gpu(1:,1:)                   !< Delta cells GPU [nb,3].
      real(R8P),    intent(in)  :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
      real(R8P),    intent(out) :: poynting_flux                     !< Power irradiated outside computational domain.
      integer(I4P)              :: i,j,k,b,v,m                       !< Counter.
      real(R8P)                 :: q_boundary(6)                     !< Variables at boundary for the Poynting flux computation.
      real(R8P)                 :: q_buff(1-FDV_S_MAX:FDV_S_MAX)     !< 1D contiguos buffer.
      real(R8P)                 :: cpb                               !< Cross product buffer, 1 component.

      poynting_flux = 0.0_R8P
      ! face -x
      !$acc parallel loop independent gang vector collapse(3)                 &
      !$acc& DEVICEVAR(q_gpu,dxyz_gpu) firstprivate(ni,nj,nk,blocks_number,s) &
      !$acc& private(v,m,cpb,q_boundary,q_buff) reduction(+: poynting_flux)
      !$omp OMPLOOP collapse(3) &
      !$omp& DEVICEPTR(q_gpu,dxyz_gpu) firstprivate(ni,nj,nk,blocks_number,s) &
      !$omp& private(v,m,cpb,q_boundary,q_buff) reduction(+: poynting_flux)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
         do v=1, 6
            do m = 1-s, s
               q_buff(m) = q_gpu(b, m, j, k, v)
            enddo
            call compute_reconstruction_r_fd_centered_dev(s=s,q=q_buff,qr=q_boundary(v))
         enddo
         cpb = (q_boundary(2) * q_boundary(6)) - (q_boundary(3) * q_boundary(5))
         poynting_flux = poynting_flux + (-cpb/MU0)*(dxyz_gpu(b,2)*dxyz_gpu(b,3))
      enddo
      enddo
      enddo
      ! face +x
      !$acc parallel loop independent gang vector collapse(3) &
      !$acc& DEVICEVAR(q_gpu,dxyz_gpu) firstprivate(s) private(cpb,q_boundary,q_buff) reduction(+: poynting_flux)
      !$omp OMPLOOP collapse(3) &
      !$omp& DEVICEPTR(q_gpu,dxyz_gpu) firstprivate(s) private(cpb,q_boundary,q_buff) reduction(+: poynting_flux)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
         do v=1, 6
            do m = 1-s, s
               q_buff(m) = q_gpu(b, ni+m, j, k, v)
            enddo
            call compute_reconstruction_r_fd_centered_dev(s=s,q=q_buff,qr=q_boundary(v))
         enddo
         cpb = (q_boundary(2) * q_boundary(6)) - (q_boundary(3) * q_boundary(5))
         poynting_flux = poynting_flux + (cpb/MU0)*(dxyz_gpu(b,2)*dxyz_gpu(b,3))
      enddo
      enddo
      enddo
      ! face -y
      !$acc parallel loop independent gang vector collapse(3) &
      !$acc& DEVICEVAR(q_gpu,dxyz_gpu) firstprivate(s) private(cpb,q_boundary,q_buff) reduction(+: poynting_flux)
      !$omp OMPLOOP collapse(3) &
      !$omp& DEVICEPTR(q_gpu,dxyz_gpu) firstprivate(s) private(cpb,q_boundary,q_buff) reduction(+: poynting_flux)
      do b=1, blocks_number
      do k=1, nk
      do i=1, ni
         do v=1, 6
            do m = 1-s, s
               q_buff(m) = q_gpu(b, i, m, k, v)
            enddo
            call compute_reconstruction_r_fd_centered_dev(s=s,q=q_buff,qr=q_boundary(v))
         enddo
         cpb = (q_boundary(3) * q_boundary(4)) - (q_boundary(1) * q_boundary(6))
         poynting_flux = poynting_flux + (-cpb/MU0)*(dxyz_gpu(b,1)*dxyz_gpu(b,3))
      enddo
      enddo
      enddo
      ! face +y
      !$acc parallel loop independent gang vector collapse(3) &
      !$acc& DEVICEVAR(q_gpu,dxyz_gpu) firstprivate(s) private(cpb,q_boundary,q_buff) reduction(+: poynting_flux)
      !$omp OMPLOOP collapse(3) &
      !$omp& DEVICEPTR(q_gpu,dxyz_gpu) firstprivate(s) private(cpb,q_boundary,q_buff) reduction(+: poynting_flux)
      do b=1, blocks_number
      do k=1, nk
      do i=1, ni
         do v=1, 6
            do m = 1-s, s
               q_buff(m) = q_gpu(b, i, nj+m, k, v)
            enddo
            call compute_reconstruction_r_fd_centered_dev(s=s,q=q_buff,qr=q_boundary(v))
         enddo
         cpb = (q_boundary(3) * q_boundary(4)) - (q_boundary(1) * q_boundary(6))
         poynting_flux = poynting_flux + (cpb/MU0)*(dxyz_gpu(b,1)*dxyz_gpu(b,3))
      enddo
      enddo
      enddo
      ! face -z
      !$acc parallel loop independent gang vector collapse(3) &
      !$acc& DEVICEVAR(q_gpu,dxyz_gpu) firstprivate(s) private(cpb,q_boundary,q_buff) reduction(+: poynting_flux)
      !$omp OMPLOOP collapse(3) &
      !$omp& DEVICEPTR(q_gpu,dxyz_gpu) firstprivate(s) private(cpb,q_boundary,q_buff) reduction(+: poynting_flux)
      do b=1, blocks_number
      do j=1, nj
      do i=1, ni
         do v=1, 6
            do m = 1-s, s
               q_buff(m) = q_gpu(b, i, j, m, v)
            enddo
            call compute_reconstruction_r_fd_centered_dev(s=s,q=q_buff,qr=q_boundary(v))
         enddo
         cpb = (q_boundary(1) * q_boundary(5)) - (q_boundary(2) * q_boundary(4))
         poynting_flux = poynting_flux + (-cpb/MU0)*(dxyz_gpu(b,1)*dxyz_gpu(b,2))
      enddo
      enddo
      enddo
      ! face +z
      !$acc parallel loop independent gang vector collapse(3) &
      !$acc& DEVICEVAR(q_gpu,dxyz_gpu) firstprivate(s) private(cpb,q_boundary,q_buff) reduction(+: poynting_flux)
      !$omp OMPLOOP collapse(3) &
      !$omp& DEVICEPTR(q_gpu,dxyz_gpu) firstprivate(s) private(cpb,q_boundary,q_buff) reduction(+: poynting_flux)
      do b=1, blocks_number
      do j=1, nj
      do i=1, ni
         do v=1, 6
            do m = 1-s, s
               q_buff(m) = q_gpu(b, i, j, nk+m, v)
            enddo
            call compute_reconstruction_r_fd_centered_dev(s=s,q=q_buff,qr=q_boundary(v))
         enddo
         cpb = (q_boundary(1) * q_boundary(5)) - (q_boundary(2) * q_boundary(4))
         poynting_flux = poynting_flux + (cpb/MU0)*(dxyz_gpu(b,1)*dxyz_gpu(b,2))
      enddo
      enddo
      enddo
      endsubroutine compute_poynting_flux_dev_kernel
   endsubroutine compute_energy

   subroutine compute_energy_error(self)
   !< Compute energy error.
   class(prism_fnl_object), intent(inout) :: self       !< The equation.
   real(R8P)                              :: energy_D0  !< Initial energy of D field.
   real(R8P)                              :: energy_B0  !< Initial energy of B field.
   real(R8P), allocatable                 :: error_D(:) !< Error (normalized) of energy of D field.
   real(R8P), allocatable                 :: error_B(:) !< Error (normalized) of energy of B field.

   if (allocated(self%energy_D).and.allocated(self%energy_B)) then
      energy_D0 = self%energy_D(1)
      energy_B0 = self%energy_B(1)

      error_D = abs(self%energy_D - energy_D0) / abs(energy_D0)
      error_B = abs(self%energy_B - energy_B0) / abs(energy_B0)

      self%rms_energy_error_D = sqrt(sum(error_D)/size(error_D))
      self%rms_energy_error_B = sqrt(sum(error_B)/size(error_B))
   endif
   endsubroutine compute_energy_error

   subroutine compute_grms(self)
   !< Compute the RMS gradient of the rotating magnetic-field amplitude over the selected domain
   !< and over its -3 dB subset.
   class(prism_fnl_object), intent(inout) :: self
   integer(I4P), allocatable              :: hi_i(:), hi_j(:), hi_k(:), lo_i(:), lo_j(:), lo_k(:)
   integer(I4P)                           :: b, i, j, k
   integer(I4P)                           :: b_ref, i_ref, j_ref, k_ref
   integer(I8P)                           :: cells_domain
   integer(I8P)                           :: cells_3db
   real(R8P)                              :: best_r2_global, best_r2_local
   real(R8P)                              :: bref_candidate, bref_global, bref_local
   real(R8P)                              :: cells_count_3db
   real(R8P)                              :: cells_count_domain
   real(R8P)                              :: domain_center(3)
   real(R8P)                              :: half_length, radius2
   real(R8P)                              :: measure_3db
   real(R8P)                              :: measure_domain
   real(R8P)                              :: threshold
   real(R8P)                              :: weighted_sum_3db
   real(R8P)                              :: weighted_sum_domain
   real(R8P)                              :: x, y, z
   logical                                :: use_cylinder

   use_cylinder = self%grms%use_cylindrical_region
   domain_center = [0._R8P, 0._R8P, 0._R8P]
   if (use_cylinder) domain_center = self%grms%center
   half_length = 0.5_R8P * self%grms%length
   radius2 = self%grms%radius * self%grms%radius
   allocate(lo_i(1:self%blocks_number), hi_i(1:self%blocks_number), lo_j(1:self%blocks_number), hi_j(1:self%blocks_number), &
            lo_k(1:self%blocks_number), hi_k(1:self%blocks_number))
   do b = 1, self%blocks_number
      call get_valid_window(self=self, hs=self%fdv_half_stencils(1), b=b, lo_i=lo_i(b), hi_i=hi_i(b), lo_j=lo_j(b), &
                            hi_j=hi_j(b), lo_k=lo_k(b), hi_k=hi_k(b))
   enddo

   best_r2_local = huge(1.0_R8P)
   bref_local = 0.0_R8P
   b_ref = 0_I4P ; i_ref = 1_I4P ; j_ref = 1_I4P ; k_ref = 1_I4P
   do b = 1, self%blocks_number
      if (lo_i(b) > hi_i(b) .or. lo_j(b) > hi_j(b) .or. lo_k(b) > hi_k(b)) cycle
      do k = lo_k(b), hi_k(b)
         z = self%adam%field%z_cell(k,b)
         do j = lo_j(b), hi_j(b)
            y = self%adam%field%y_cell(j,b)
            do i = lo_i(b), hi_i(b)
               x = self%adam%field%x_cell(i,b)
               if (.not. is_inside_selected_domain(x=x, y=y, z=z)) cycle
               if (reference_distance2(x=x, y=y, z=z) < best_r2_local) then
                  best_r2_local = reference_distance2(x=x, y=y, z=z)
                  b_ref = b ; i_ref = i ; j_ref = j ; k_ref = k
               endif
            enddo
         enddo
      enddo
   enddo

   if (best_r2_local < huge(1.0_R8P)) then
      call compute_bref_dev_kernel(b_ref=b_ref, i_ref=i_ref, j_ref=j_ref, k_ref=k_ref, ngc=self%ngc, q_gpu=self%q_gpu, bref=bref_local)
   endif
   call MPI_ALLREDUCE(best_r2_local, best_r2_global, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, mpih_fnl%error)
   bref_candidate = 0.0_R8P
   if (abs(best_r2_local - best_r2_global) <= 1.0E-12_R8P * max(1.0_R8P, abs(best_r2_global))) bref_candidate = bref_local
   bref_global = bref_candidate
   call MPI_ALLREDUCE(MPI_IN_PLACE, bref_global, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, mpih_fnl%error)

   threshold = GRMS_3DB_RATIO * bref_global
   call compute_grms_dev_kernel(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                                s1=self%fdv_half_stencils(1), threshold=threshold, use_cylinder=use_cylinder, &
                                center_x=self%grms%center(1), center_y=self%grms%center(2), center_z=self%grms%center(3), &
                                axis_x=self%grms%axis(1), axis_y=self%grms%axis(2), axis_z=self%grms%axis(3), &
                                half_length=half_length, radius2=radius2, lo_i=lo_i, hi_i=hi_i, lo_j=lo_j, hi_j=hi_j, &
                                lo_k=lo_k, hi_k=hi_k, x_cell_gpu=self%field_fnl%x_cell_gpu, y_cell_gpu=self%field_fnl%y_cell_gpu, &
                                z_cell_gpu=self%field_fnl%z_cell_gpu, dxyz_gpu=self%field_fnl%dxyz_gpu, q_gpu=self%q_gpu, &
                                weighted_sum_domain=weighted_sum_domain, measure_domain=measure_domain, &
                                cells_count_domain=cells_count_domain, weighted_sum_3db=weighted_sum_3db, &
                                measure_3db=measure_3db, cells_count_3db=cells_count_3db)
   call MPI_ALLREDUCE(MPI_IN_PLACE, weighted_sum_domain, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih_fnl%error)
   call MPI_ALLREDUCE(MPI_IN_PLACE, weighted_sum_3db, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih_fnl%error)
   call MPI_ALLREDUCE(MPI_IN_PLACE, measure_domain, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih_fnl%error)
   call MPI_ALLREDUCE(MPI_IN_PLACE, measure_3db, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih_fnl%error)
   call MPI_ALLREDUCE(MPI_IN_PLACE, cells_count_domain, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih_fnl%error)
   call MPI_ALLREDUCE(MPI_IN_PLACE, cells_count_3db, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih_fnl%error)

   cells_domain = int(nint(cells_count_domain), I8P)
   cells_3db = int(nint(cells_count_3db), I8P)
   self%grms%reference_B = bref_global
   self%grms%threshold_B = threshold
   self%grms%domain_measure = measure_domain
   self%grms%measure_3db = measure_3db
   self%grms%domain_cells_number = cells_domain
   self%grms%cells_number_3db = cells_3db
   self%grms%grms_domain_B = 0.0_R8P
   self%grms%grms_3db_B = 0.0_R8P
   if (measure_domain > 0.0_R8P) self%grms%grms_domain_B = sqrt(weighted_sum_domain / measure_domain)
   if (measure_3db > 0.0_R8P) self%grms%grms_3db_B = sqrt(weighted_sum_3db / measure_3db)
   deallocate(lo_i, hi_i, lo_j, hi_j, lo_k, hi_k)
   contains
      subroutine get_valid_window(self, hs, b, lo_i, hi_i, lo_j, hi_j, lo_k, hi_k)
      class(prism_fnl_object), intent(in)  :: self
      integer(I4P),            intent(in)  :: hs
      integer(I4P),            intent(in)  :: b
      integer(I4P),            intent(out) :: lo_i, hi_i, lo_j, hi_j, lo_k, hi_k

      lo_i = 1_I4P ; hi_i = self%ni
      lo_j = 1_I4P ; hi_j = self%nj
      lo_k = 1_I4P ; hi_k = self%nk
      if (allocated(self%fWLayer%ni_fWL)) then
         call exclude_face(lo=lo_i, hi=hi_i, face_first=self%fWLayer%ni_fWL(1,b,PML_FACE_X_M), &
                           face_last=self%fWLayer%ni_fWL(2,b,PML_FACE_X_M), is_minus=.true., hs=hs)
         call exclude_face(lo=lo_i, hi=hi_i, face_first=self%fWLayer%ni_fWL(1,b,PML_FACE_X_P), &
                           face_last=self%fWLayer%ni_fWL(2,b,PML_FACE_X_P), is_minus=.false., hs=hs)
         call exclude_face(lo=lo_j, hi=hi_j, face_first=self%fWLayer%nj_fWL(1,b,PML_FACE_Y_M), &
                           face_last=self%fWLayer%nj_fWL(2,b,PML_FACE_Y_M), is_minus=.true., hs=hs)
         call exclude_face(lo=lo_j, hi=hi_j, face_first=self%fWLayer%nj_fWL(1,b,PML_FACE_Y_P), &
                           face_last=self%fWLayer%nj_fWL(2,b,PML_FACE_Y_P), is_minus=.false., hs=hs)
         call exclude_face(lo=lo_k, hi=hi_k, face_first=self%fWLayer%nk_fWL(1,b,PML_FACE_Z_M), &
                           face_last=self%fWLayer%nk_fWL(2,b,PML_FACE_Z_M), is_minus=.true., hs=hs)
         call exclude_face(lo=lo_k, hi=hi_k, face_first=self%fWLayer%nk_fWL(1,b,PML_FACE_Z_P), &
                           face_last=self%fWLayer%nk_fWL(2,b,PML_FACE_Z_P), is_minus=.false., hs=hs)
      endif
      if (allocated(self%pml%ni_pml)) then
         call exclude_face(lo=lo_i, hi=hi_i, face_first=self%pml%ni_pml(1,b,PML_FACE_X_M), &
                           face_last=self%pml%ni_pml(2,b,PML_FACE_X_M), is_minus=.true., hs=hs)
         call exclude_face(lo=lo_i, hi=hi_i, face_first=self%pml%ni_pml(1,b,PML_FACE_X_P), &
                           face_last=self%pml%ni_pml(2,b,PML_FACE_X_P), is_minus=.false., hs=hs)
         call exclude_face(lo=lo_j, hi=hi_j, face_first=self%pml%nj_pml(1,b,PML_FACE_Y_M), &
                           face_last=self%pml%nj_pml(2,b,PML_FACE_Y_M), is_minus=.true., hs=hs)
         call exclude_face(lo=lo_j, hi=hi_j, face_first=self%pml%nj_pml(1,b,PML_FACE_Y_P), &
                           face_last=self%pml%nj_pml(2,b,PML_FACE_Y_P), is_minus=.false., hs=hs)
         call exclude_face(lo=lo_k, hi=hi_k, face_first=self%pml%nk_pml(1,b,PML_FACE_Z_M), &
                           face_last=self%pml%nk_pml(2,b,PML_FACE_Z_M), is_minus=.true., hs=hs)
         call exclude_face(lo=lo_k, hi=hi_k, face_first=self%pml%nk_pml(1,b,PML_FACE_Z_P), &
                           face_last=self%pml%nk_pml(2,b,PML_FACE_Z_P), is_minus=.false., hs=hs)
      endif
      endsubroutine get_valid_window

      pure subroutine exclude_face(lo, hi, face_first, face_last, is_minus, hs)
      integer(I4P), intent(inout) :: lo
      integer(I4P), intent(inout) :: hi
      integer(I4P), intent(in)    :: face_first
      integer(I4P), intent(in)    :: face_last
      logical,      intent(in)    :: is_minus
      integer(I4P), intent(in)    :: hs

      if (face_first <= 0_I4P .or. face_last <= 0_I4P) return
      if (is_minus) then
         lo = max(lo, face_last + hs + 1_I4P)
      else
         hi = min(hi, face_first - hs - 1_I4P)
      endif
      endsubroutine exclude_face

      subroutine compute_bref_dev_kernel(b_ref, i_ref, j_ref, k_ref, ngc, q_gpu, bref)
      integer(I4P), intent(in)  :: b_ref, i_ref, j_ref, k_ref, ngc
      real(R8P),    intent(in)  :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
      real(R8P),    intent(out) :: bref
      integer(I4P)              :: one

      bref = 0.0_R8P
      !$acc parallel loop independent DEVICEVAR(q_gpu) private(one) reduction(max:bref)
      !$omp OMPLOOP DEVICEPTR(q_gpu) private(one) reduction(max:bref)
      do one = 1, 1
         bref = sqrt(q_gpu(b_ref,i_ref,j_ref,k_ref,VAR_BX)**2 + q_gpu(b_ref,i_ref,j_ref,k_ref,VAR_BY)**2)
      enddo
      endsubroutine compute_bref_dev_kernel

      logical function is_inside_selected_domain(x, y, z)
      real(R8P), intent(in) :: x, y, z
      real(R8P)             :: axial_distance, radial2_local

      if (.not. use_cylinder) then
         is_inside_selected_domain = .true.
         return
      endif
      axial_distance = (x - self%grms%center(1)) * self%grms%axis(1) + &
                       (y - self%grms%center(2)) * self%grms%axis(2) + &
                       (z - self%grms%center(3)) * self%grms%axis(3)
      if (abs(axial_distance) > half_length) then
         is_inside_selected_domain = .false.
         return
      endif
      radial2_local = (x - self%grms%center(1))**2 + (y - self%grms%center(2))**2 + (z - self%grms%center(3))**2 - &
                      axial_distance * axial_distance
      is_inside_selected_domain = radial2_local <= radius2
      endfunction is_inside_selected_domain

      real(R8P) function reference_distance2(x, y, z)
      real(R8P), intent(in) :: x, y, z

      reference_distance2 = (x - domain_center(1))**2 + (y - domain_center(2))**2 + (z - domain_center(3))**2
      endfunction reference_distance2

      subroutine compute_grms_dev_kernel(ni, nj, nk, ngc, blocks_number, s1, threshold, use_cylinder, center_x, center_y, center_z, &
                                         axis_x, axis_y, axis_z, half_length, radius2, lo_i, hi_i, lo_j, hi_j, lo_k, hi_k, &
                                         x_cell_gpu, y_cell_gpu, z_cell_gpu, dxyz_gpu, q_gpu, weighted_sum_domain, measure_domain, &
                                         cells_count_domain, weighted_sum_3db, measure_3db, cells_count_3db)
      integer(I4P), intent(in)  :: ni, nj, nk, ngc, blocks_number, s1
      integer(I4P), intent(in)  :: lo_i(1:), hi_i(1:), lo_j(1:), hi_j(1:), lo_k(1:), hi_k(1:)
      real(R8P),    intent(in)  :: threshold
      logical,      intent(in)  :: use_cylinder
      real(R8P),    intent(in)  :: center_x, center_y, center_z, axis_x, axis_y, axis_z, half_length, radius2
      real(R8P),    intent(in)  :: x_cell_gpu(1:,1-ngc:)
      real(R8P),    intent(in)  :: y_cell_gpu(1:,1-ngc:)
      real(R8P),    intent(in)  :: z_cell_gpu(1:,1-ngc:)
      real(R8P),    intent(in)  :: dxyz_gpu(1:,1:)
      real(R8P),    intent(in)  :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
      real(R8P),    intent(out) :: weighted_sum_domain, measure_domain, cells_count_domain
      real(R8P),    intent(out) :: weighted_sum_3db, measure_3db, cells_count_3db
      real(R8P)                 :: axial_distance, b_amp, b_minus, b_plus, dBdx, dBdy, dBdz, grad2, radial2_local
      integer(I4P)              :: b, i, j, k, s

      weighted_sum_domain = 0.0_R8P
      weighted_sum_3db = 0.0_R8P
      measure_domain = 0.0_R8P
      measure_3db = 0.0_R8P
      cells_count_domain = 0.0_R8P
      cells_count_3db = 0.0_R8P
      !$acc parallel loop collapse(4) DEVICEVAR(x_cell_gpu,y_cell_gpu,z_cell_gpu,dxyz_gpu,q_gpu) copyin(lo_i,hi_i,lo_j,hi_j,lo_k,hi_k) &
      !$acc& firstprivate(ni,nj,nk,blocks_number,s1,threshold,use_cylinder,center_x,center_y,center_z,axis_x,axis_y,axis_z,half_length,radius2) &
      !$acc& private(axial_distance,b_amp,b_minus,b_plus,dBdx,dBdy,dBdz,grad2,radial2_local,s) &
      !$acc& reduction(+: weighted_sum_domain, measure_domain, cells_count_domain, weighted_sum_3db, measure_3db, cells_count_3db)
      !$omp OMPLOOP collapse(4) DEVICEPTR(x_cell_gpu,y_cell_gpu,z_cell_gpu,dxyz_gpu,q_gpu) map(to:lo_i,hi_i,lo_j,hi_j,lo_k,hi_k) &
      !$omp& firstprivate(ni,nj,nk,blocks_number,s1,threshold,use_cylinder,center_x,center_y,center_z,axis_x,axis_y,axis_z,half_length,radius2) &
      !$omp& private(axial_distance,b_amp,b_minus,b_plus,dBdx,dBdy,dBdz,grad2,radial2_local,s) &
      !$omp& reduction(+: weighted_sum_domain, measure_domain, cells_count_domain, weighted_sum_3db, measure_3db, cells_count_3db)
      do b = 1, blocks_number
      do k = 1, nk
      do j = 1, nj
      do i = 1, ni
         if (i < lo_i(b) .or. i > hi_i(b) .or. j < lo_j(b) .or. j > hi_j(b) .or. k < lo_k(b) .or. k > hi_k(b)) cycle
         if (use_cylinder) then
            axial_distance = (x_cell_gpu(b,i) - center_x) * axis_x + (y_cell_gpu(b,j) - center_y) * axis_y + &
                             (z_cell_gpu(b,k) - center_z) * axis_z
            if (abs(axial_distance) > half_length) cycle
            radial2_local = (x_cell_gpu(b,i) - center_x)**2 + (y_cell_gpu(b,j) - center_y)**2 + (z_cell_gpu(b,k) - center_z)**2 - &
                            axial_distance * axial_distance
            if (radial2_local > radius2) cycle
         endif
         b_amp = sqrt(q_gpu(b,i,j,k,VAR_BX)**2 + q_gpu(b,i,j,k,VAR_BY)**2)
         dBdx = 0.0_R8P
         dBdy = 0.0_R8P
         dBdz = 0.0_R8P
         !$acc loop seq
         do s = 1, s1
            b_plus = sqrt(q_gpu(b,i+s,j,k,VAR_BX)**2 + q_gpu(b,i+s,j,k,VAR_BY)**2)
            b_minus = sqrt(q_gpu(b,i-s,j,k,VAR_BX)**2 + q_gpu(b,i-s,j,k,VAR_BY)**2)
            dBdx = dBdx + FD1_CC(s,s1) * (b_plus - b_minus) / dxyz_gpu(b,1)
            b_plus = sqrt(q_gpu(b,i,j+s,k,VAR_BX)**2 + q_gpu(b,i,j+s,k,VAR_BY)**2)
            b_minus = sqrt(q_gpu(b,i,j-s,k,VAR_BX)**2 + q_gpu(b,i,j-s,k,VAR_BY)**2)
            dBdy = dBdy + FD1_CC(s,s1) * (b_plus - b_minus) / dxyz_gpu(b,2)
            b_plus = sqrt(q_gpu(b,i,j,k+s,VAR_BX)**2 + q_gpu(b,i,j,k+s,VAR_BY)**2)
            b_minus = sqrt(q_gpu(b,i,j,k-s,VAR_BX)**2 + q_gpu(b,i,j,k-s,VAR_BY)**2)
            dBdz = dBdz + FD1_CC(s,s1) * (b_plus - b_minus) / dxyz_gpu(b,3)
         enddo
         grad2 = dBdx*dBdx + dBdy*dBdy + dBdz*dBdz
         weighted_sum_domain = weighted_sum_domain + grad2 * (dxyz_gpu(b,1) * dxyz_gpu(b,2) * dxyz_gpu(b,3))
         measure_domain = measure_domain + (dxyz_gpu(b,1) * dxyz_gpu(b,2) * dxyz_gpu(b,3))
         cells_count_domain = cells_count_domain + 1.0_R8P
         if (b_amp >= threshold) then
            weighted_sum_3db = weighted_sum_3db + grad2 * (dxyz_gpu(b,1) * dxyz_gpu(b,2) * dxyz_gpu(b,3))
            measure_3db = measure_3db + (dxyz_gpu(b,1) * dxyz_gpu(b,2) * dxyz_gpu(b,3))
            cells_count_3db = cells_count_3db + 1.0_R8P
         endif
      enddo
      enddo
      enddo
      enddo
      endsubroutine compute_grms_dev_kernel
   endsubroutine compute_grms

   subroutine compute_magnetic_field_at_center_domain(self)
   !< Compute the magnetic field in the cell center closest to the geometrical domain center.
   class(prism_fnl_object), intent(inout) :: self
   integer(I4P)                          :: b, i, j, k
   integer(I4P)                          :: b_ref, i_ref, j_ref, k_ref
   integer(I4P)                          :: owner_rank
   real(R8P)                             :: best_r2_global(2)
   real(R8P)                             :: best_r2_local(2)
   real(R8P)                             :: center(3)
   real(R8P)                             :: distance2
   real(R8P)                             :: magnetic_field(3)
   real(R8P)                             :: sample_point(3)
   real(R8P)                             :: x, y, z

   center = 0.5_R8P * (self%adam%grid%domain_emin + self%adam%grid%domain_emax)
   best_r2_local = [huge(1.0_R8P), real(mpih_fnl%myrank, R8P)]
   magnetic_field = 0.0_R8P
   sample_point = center
   b_ref = 0_I4P ; i_ref = 1_I4P ; j_ref = 1_I4P ; k_ref = 1_I4P

   do b = 1, self%blocks_number
      do k = 1, self%nk
         z = self%adam%field%z_cell(k,b)
         do j = 1, self%nj
            y = self%adam%field%y_cell(j,b)
            do i = 1, self%ni
               x = self%adam%field%x_cell(i,b)
               distance2 = (x - center(1))**2 + (y - center(2))**2 + (z - center(3))**2
               if (distance2 < best_r2_local(1)) then
                  best_r2_local(1) = distance2
                  b_ref = b ; i_ref = i ; j_ref = j ; k_ref = k
                  sample_point = [x, y, z]
               endif
            enddo
         enddo
      enddo
   enddo

   if (b_ref > 0_I4P) then
      call compute_magnetic_field_at_center_domain_dev_kernel(b_ref=b_ref, i_ref=i_ref, j_ref=j_ref, k_ref=k_ref, &
                                                              ngc=self%ngc, q_gpu=self%q_gpu, bx=magnetic_field(1), &
                                                              by=magnetic_field(2), bz=magnetic_field(3))
   endif
   call MPI_ALLREDUCE(best_r2_local, best_r2_global, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_WORLD, mpih_fnl%error)
   owner_rank = nint(best_r2_global(2), I4P)
   call MPI_BCAST(magnetic_field, 3, MPI_REAL8, owner_rank, MPI_COMM_WORLD, mpih_fnl%error)
   call MPI_BCAST(sample_point, 3, MPI_REAL8, owner_rank, MPI_COMM_WORLD, mpih_fnl%error)

   self%magnetic_field_at_center_domain%center = center
   self%magnetic_field_at_center_domain%sample_point = sample_point
   self%magnetic_field_at_center_domain%magnetic_field = magnetic_field
   self%magnetic_field_at_center_domain%distance = sqrt(best_r2_global(1))
   contains
      subroutine compute_magnetic_field_at_center_domain_dev_kernel(b_ref, i_ref, j_ref, k_ref, ngc, q_gpu, bx, by, bz)
      integer(I4P), intent(in)  :: b_ref, i_ref, j_ref, k_ref, ngc
      real(R8P),    intent(in)  :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
      real(R8P),    intent(out) :: bx, by, bz
      integer(I4P)              :: one

      bx = 0.0_R8P ; by = 0.0_R8P ; bz = 0.0_R8P
      !$acc parallel loop independent DEVICEVAR(q_gpu) private(one) reduction(+:bx,by,bz)
      !$omp OMPLOOP DEVICEPTR(q_gpu) private(one) reduction(+:bx,by,bz)
      do one = 1, 1
         bx = q_gpu(b_ref,i_ref,j_ref,k_ref,VAR_BX)
         by = q_gpu(b_ref,i_ref,j_ref,k_ref,VAR_BY)
         bz = q_gpu(b_ref,i_ref,j_ref,k_ref,VAR_BZ)
      enddo
      endsubroutine compute_magnetic_field_at_center_domain_dev_kernel
   endsubroutine compute_magnetic_field_at_center_domain

   subroutine compute_max_divergence(self)
   !< Compute maximum divergence.
   class(prism_fnl_object), intent(inout) :: self       !< The equation.
   real(R8P)                              :: max_div(3) !< Maximum divergence.
   integer(I4P), allocatable              :: fwl_cells(:,:) !< Block-local fWLayer cell counts; zero when layer is disabled.
   integer(I4P)                           :: rho_ivar   !< Charge-density slot for PIC, appended at the end of q.
   integer(I4P)                           :: use_rho    !< Integerized PIC flag for OpenACC firstprivate handling.

   allocate(fwl_cells(1:self%blocks_number,1:6))
   fwl_cells = 0_I4P
   if (allocated(self%fWLayer%C)) fwl_cells = self%fWLayer%C(1:self%blocks_number,1:6)
   rho_ivar = self%nv
   use_rho = 0_I4P
   if (self%physics%physical_model == PIC_PHYSICAL_MODEL) use_rho = 1_I4P

	call compute_max_divergence_dev_kernel(ni            = self%ni                  ,&
                                          nj            = self%nj                  ,&
                                          nk            = self%nk                  ,&
                                          blocks_number = self%blocks_number       ,&
															ngc           = self%ngc                 ,&
                                          var_jx        = self%physics%var_jx      ,&
                                          var_jy        = self%physics%var_jy      ,&
                                          var_jz        = self%physics%var_jz      ,&
                                          rho_ivar      = rho_ivar                 ,&
                                          use_rho       = use_rho                  ,&
                                          s1            = self%fdv_half_stencils(1),&
                                          fwl_c         = fwl_cells                ,&
                                          dxyz_gpu      = self%field_fnl%dxyz_gpu  ,&
                                          q_gpu         = self%q_gpu               ,&
                                          max_div       = max_div)
   deallocate(fwl_cells)

	self%max_divergence_D = max_div(1)
	self%max_divergence_B = max_div(2)
	self%max_divergence_J = max_div(3)
   contains
      subroutine compute_max_divergence_dev_kernel(ni, nj, nk, blocks_number, ngc, &
                                                   var_Jx, var_Jy, var_Jz, rho_ivar, use_rho, s1, fwl_c, &
                                                   dxyz_gpu, q_gpu, max_div)

			!< Compute maximum divergence of D, B and J fields, device kernel.
			integer(I4P), intent(in)  :: ni, nj, nk, blocks_number, ngc        !< Grids dimensions.
			integer(I4P), intent(in)  :: var_Jx, var_Jy, var_Jz                !< Current variables indices.
			integer(I4P), intent(in)  :: rho_ivar                              !< Charge-density slot for PIC.
			integer(I4P), intent(in)  :: use_rho                               !< Integerized PIC flag.
			integer(I4P), intent(in)  :: s1                                    !< FDV half stencil.
			integer(I4P), intent(in)  :: fwl_c(1:,1:)                          !< fWLayer cell counts [nb,6].
			real(R8P),    intent(in)  :: dxyz_gpu(1:,1:)                       !< Delta cells GPU [nb,3].
			real(R8P),    intent(in)  :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Conservative variables.
			real(R8P),    intent(out) :: max_div(3)                            !< Maximum divergence of D, B and J fields.
			! Stencil buffers MUST have compile-time-constant bounds (FDV_S_MAX), like the
			! residual kernels: runtime-sized (automatic) private arrays inside the acc
			! collapse(4) region are mis-privatized by nvfortran -- scheduling-dependent
			! garbage that grows with gang count (issue #22 F1: the divergence history read
			! O(1e3) on a bit-perfect solution, nondeterministically, above ~32 blocks).
			real(R8P)                 :: divergenceD, divergenceB, divergenceJ !< Divergence of D, B and J fields.
			real(R8P)                 :: max_divD, max_divB, max_divJ			 !< Maximum divergence of D, B and J fields.
			real(R8P)                 :: dxyz_b(3)                             !< Per-block deltas, PRIVATE copy (no strided-section temp: issue #22 F1-bis).
			integer(I4P)              :: i,j,k,b,s                             !< Counter.
			integer(I4P)              :: lo_i, hi_i, lo_j, hi_j, lo_k, hi_k    !< fWLayer skin-exclusion bounds (CPU-parity).

			max_divD = 0.0_R8P
			max_divB = 0.0_R8P
			max_divJ = 0.0_R8P
		      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu) copyin(fwl_c) &
	         !$acc& firstprivate(var_jx,var_jy,var_jz,rho_ivar,use_rho,s1)                             &
	         !$acc& private(divergenceD,divergenceB,divergenceJ,dxyz_b,lo_i,hi_i,lo_j,hi_j,lo_k,hi_k) &
				!$acc& reduction(max: max_divD, max_divB, max_divJ)
		      !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,q_gpu) map(to:fwl_c) &
	         !$omp& firstprivate(var_jx,var_jy,var_jz,rho_ivar,use_rho,s1)                             &
	         !$omp& private(divergenceD,divergenceB,divergenceJ,dxyz_b,lo_i,hi_i,lo_j,hi_j,lo_k,hi_k) &
				!$omp& reduction(max: max_divD, max_divB, max_divJ)
	         do b=1,blocks_number
	         do k=1,nk
	         do j=1,nj
	         do i=1,ni
	            dxyz_b(1) = dxyz_gpu(b,1) ; dxyz_b(2) = dxyz_gpu(b,2) ; dxyz_b(3) = dxyz_gpu(b,3)
               lo_i = 1_I4P ; hi_i = ni
               lo_j = 1_I4P ; hi_j = nj
               lo_k = 1_I4P ; hi_k = nk
               if (fwl_c(b,1) > 0_I4P) lo_i = 1_I4P + fwl_c(b,1) + s1
               if (fwl_c(b,2) > 0_I4P) hi_i = ni - (fwl_c(b,2) + s1 - 1_I4P)
               if (fwl_c(b,3) > 0_I4P) lo_j = 1_I4P + fwl_c(b,3) + s1
               if (fwl_c(b,4) > 0_I4P) hi_j = nj - (fwl_c(b,4) + s1 - 1_I4P)
               if (fwl_c(b,5) > 0_I4P) lo_k = 1_I4P + fwl_c(b,5) + s1
               if (fwl_c(b,6) > 0_I4P) hi_k = nk - (fwl_c(b,6) + s1 - 1_I4P)
	            ! Buffer-free divergences (issue #22 F1-bis): the former private stencil
            ! buffers of this CONTAINED kernel were mis-privatized by nvfortran even
            ! with constant bounds (threads bled each other's fills: div(J) tracked
            ! div(D) with J identically zero, nondeterministically, seam-only because
            ! zero planes are race-invisible on the uniform cases). Scalars only:
            ! pair-form FD1_CC accumulation, no arrays, no callees.
            divergenceD = 0._R8P
            divergenceB = 0._R8P
            divergenceJ = 0._R8P
            !$acc loop seq
            do s=1, s1
               divergenceD = divergenceD + FD1_CC(s,s1)*((q_gpu(b,i+s,j,k,VAR_DX) - q_gpu(b,i-s,j,k,VAR_DX))/dxyz_b(1)  &
                                                       + (q_gpu(b,i,j+s,k,VAR_DY) - q_gpu(b,i,j-s,k,VAR_DY))/dxyz_b(2)  &
                                                       + (q_gpu(b,i,j,k+s,VAR_DZ) - q_gpu(b,i,j,k-s,VAR_DZ))/dxyz_b(3))
               divergenceB = divergenceB + FD1_CC(s,s1)*((q_gpu(b,i+s,j,k,VAR_BX) - q_gpu(b,i-s,j,k,VAR_BX))/dxyz_b(1)  &
                                                       + (q_gpu(b,i,j+s,k,VAR_BY) - q_gpu(b,i,j-s,k,VAR_BY))/dxyz_b(2)  &
                                                       + (q_gpu(b,i,j,k+s,VAR_BZ) - q_gpu(b,i,j,k-s,VAR_BZ))/dxyz_b(3))
               divergenceJ = divergenceJ + FD1_CC(s,s1)*((q_gpu(b,i+s,j,k,var_Jx) - q_gpu(b,i-s,j,k,var_Jx))/dxyz_b(1)  &
                                                       + (q_gpu(b,i,j+s,k,var_Jy) - q_gpu(b,i,j-s,k,var_Jy))/dxyz_b(2)  &
                                                       + (q_gpu(b,i,j,k+s,var_Jz) - q_gpu(b,i,j,k-s,var_Jz))/dxyz_b(3))
            enddo
            if (use_rho /= 0_I4P) divergenceD = divergenceD - q_gpu(b,i,j,k,rho_ivar)
	            ! fWLayer-skin exclusion (CPU-parity): only cells outside the local block's
	            ! layer skin contribute to the reported maximum.
	            if (i >= lo_i .and. i <= hi_i .and. j >= lo_j .and. j <= hi_j .and. k >= lo_k .and. k <= hi_k) then
	               max_divD = max(max_divD, abs(divergenceD))
	               max_divB = max(max_divB, abs(divergenceB))
	               max_divJ = max(max_divJ, abs(divergenceJ))
            endif
         enddo
         enddo
         enddo
         enddo
			max_div(1) = max_divD
			max_div(2) = max_divB
			max_div(3) = max_divJ
		endsubroutine compute_max_divergence_dev_kernel
	endsubroutine compute_max_divergence

   subroutine impose_ct_correction_dev(self, ivar)
   !< Impose Constrained Transport Correction on vectorial variable q(ivar:ivar+2).
   !< Note that self%divergence memory is used as buffer, be carefull.
   class(prism_fnl_object), intent(inout) :: self   !< The equation.
   integer(I4P),            intent(in)    :: ivar   !< Variable (start) index in q.
   real(R8P)                              :: dq_max !< Maximum residual.
   integer(I4P)                           :: iter   !< Counter.

   associate(blocks_number=>self%blocks_number)
   if (blocks_number>0) then
      ! call self%compute_divergence(ivar=ivar,ovar=4,q_gpu=q_gpu,divergence_gpu=buffer)
      do iter=1, self%flail%iterations
         ! call compute_smoothing_multigrid(ni=ni,nj=nj,nk=nk,ngc=ngc,nv=1_I4P,blocks_number=blocks_number, &
         !                                  dxyz=field%dxyz,                                                &
         !                                  f=-buffer(4:4,:,:,:,:),                                         &
         !                                  q=buffer(7:7,:,:,:,:),                                          &
         !                                  dq_max=dq_max,                                                  &
         !                                  dq=buffer(5:5,:,:,:,:),                                         &
         !                                  iterations_init=self%flail%iterations_init,                     &
         !                                  iterations_fine=self%flail%iterations_fine,                     &
         !                                  iterations_coarse=self%flail%iterations_coarse)
         if (dq_max < self%flail%tolerance) exit
      enddo
      call mpih_fnl%print_message('FLAIL convergence reached at iteration '//trim(str(iter,.true.)))
      ! call self%compute_gradient(ivar=1,q_gpu=buffer(:,:,:,:,7:7),gradient_gpu=buffer(:,:,:,:,4:6))
      call impose_ct_correction_kernel(ni            = self%ni,            &
                                       nj            = self%nj,            &
                                       nk            = self%nk,            &
                                       ngc           = self%ngc,           &
                                       blocks_number = blocks_number,      &
                                       ivar          = ivar,               &
                                       q_gpu         = self%q_gpu,         &
                                       buffer_gpu    = self%divergence_gpu)
   endif
   endassociate
   contains
      subroutine impose_ct_correction_kernel(ni, nj, nk, ngc, blocks_number, ivar, q_gpu, buffer_gpu)
      !< Apply the CT correction increment, device kernel.
      integer(I4P), intent(in)    :: ni, nj, nk, ngc, blocks_number          !< Grid dimensions.
      integer(I4P), intent(in)    :: ivar                                    !< Variable (start) index in q.
      real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)       !< Conservative variables.
      real(R8P),    intent(in)    :: buffer_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Correction buffer (self%divergence_gpu shape).
      integer(I4P)                :: i, j, k, b, v                           !< Counter.

      !$acc parallel loop independent gang vector collapse(5) DEVICEVAR(q_gpu,buffer_gpu) &
      !$acc& firstprivate(ni,nj,nk,blocks_number,ivar)
      !$omp OMPLOOP collapse(5) DEVICEPTR(q_gpu,buffer_gpu) &
      !$omp& firstprivate(ni,nj,nk,blocks_number,ivar)
      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do v=1, 3
                     q_gpu(b,i,j,k,ivar+v-1) = q_gpu(b,i,j,k,ivar+v-1) + buffer_gpu(b,i,j,k,3+v)
                  enddo
               enddo
            enddo
         enddo
      enddo
      endsubroutine impose_ct_correction_kernel
   endsubroutine impose_ct_correction_dev

   subroutine impose_div_free(self)
   !< Impose divergence-free property.
   class(prism_fnl_object), intent(inout) :: self !< The equation.

   associate(constrained_transport_D=>self%numerics%constrained_transport_D,&
             constrained_transport_B=>self%numerics%constrained_transport_B,div_corr_var=>self%numerics%div_corr_var)
   if (constrained_transport_D.and.div_corr_var==DIV_CORR_VAR_POISS) call self%impose_ct_correction_dev(ivar=1_I4P)
   if (constrained_transport_B.and.div_corr_var==DIV_CORR_VAR_POISS) call self%impose_ct_correction_dev(ivar=4_I4P)
   ! here should go also other corrections...
   endassociate
   endsubroutine impose_div_free
endmodule adam_prism_fnl_object
