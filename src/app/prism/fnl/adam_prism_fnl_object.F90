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

type, extends(prism_common_object) :: prism_fnl_object
   !< PRISM equations system class definition, GPU (FNL) backend.
   ! PRISM FNL C.3 closure (issue #13 D.4a follow-up): the 6 genuinely per-realm
   ! FNL objects (field/ib/rk/weno/coil/fwlayer) are value components reached
   ! through self%, replacing the former module-scope singletons + pointer shims.
   ! The MPI handler (mpih_fnl) is NOT per-realm — it is one-per-process (rank,
   ! device, communicator) and stays a true program-scope singleton in
   ! adam_fnl_mpih_global, mirroring the CPU adam_mpih_global.
   type(field_fnl_object)         :: field_fnl   !< GPU field handler.
   type(ib_fnl_object)            :: ib_fnl      !< GPU immersed boundary.
   type(rk_fnl_object)            :: rk_fnl      !< GPU Runge-Kutta integrator.
   type(weno_fnl_object)          :: weno_fnl    !< GPU WENO reconstructor.
   type(prism_fnl_coil_object)    :: coil_fnl    !< GPU coil source.
   type(prism_fnl_fwlayer_object) :: fwlayer_fnl !< GPU fWLayer.
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
      procedure, pass(self) :: initialize_prism !< Initialize PRISM equation.
      ! IO methods
      procedure, pass(self) :: load_restart_files   !< Load restart files.
      procedure, pass(self) :: save_residuals       !< Save residuals history.
      procedure, pass(self) :: save_simulation_data !< Save all simulation data.
      ! IC/BC/sources
      procedure, pass(self) :: apply_fwl_correction    !< Apply fWLayer correction (if present)
      procedure, pass(self) :: compute_coils_current   !< Compute current coils sources.
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
      ! numerical methods, miscellanea
      procedure, pass(self) :: compute_dt                        !< Compute time step.
      procedure, pass(self) :: initialize_forest                 !< Orchestrator contract; overrides realm_object default.
      procedure, pass(self) :: compute_local_dt_forest           !< Orchestrator contract; overrides realm_object default.
      procedure, pass(self) :: advance_one_step_forest           !< Orchestrator contract; overrides realm_object default (N=1 fast path).
      procedure, pass(self) :: nrk_forest                        !< Orchestrator contract; overrides realm_object default (multi-realm path).
      procedure, pass(self) :: prepare_step_forest               !< Orchestrator contract; overrides realm_object default (multi-realm path).
      procedure, pass(self) :: assemble_substage_forest          !< Orchestrator contract; overrides realm_object default (multi-realm path).
      procedure, pass(self) :: residuals_substage_forest         !< Orchestrator contract; overrides realm_object default (multi-realm 3-phase loop, issue #13).
      procedure, pass(self) :: assign_substage_forest            !< Orchestrator contract; overrides realm_object default (multi-realm 3-phase loop, issue #13).
      procedure, pass(self) :: evaluate_substage_forest          !< Orchestrator contract; overrides realm_object default (multi-realm path, legacy 2-phase).
      procedure, pass(self) :: finalize_step_forest              !< Orchestrator contract; overrides realm_object default (multi-realm path).
      procedure, pass(self) :: post_step_forest                  !< Orchestrator contract; overrides realm_object default.
      procedure, pass(self) :: is_done_forest                    !< Orchestrator contract; overrides realm_object default.
      procedure, pass(self) :: finalize_forest                   !< Orchestrator contract; overrides realm_object default.
      procedure, pass(self) :: finalize_mpi_forest               !< Process-global MPI finalize (mpih_fnl); overrides realm_object default.
      procedure, pass(self) :: exchange_inter_realm_halos_forest !< Orchestrator contract; overrides realm_object default.
      procedure, pass(self) :: compute_energy       	!< Compute energy.
      procedure, pass(self) :: compute_energy_error 	!< Compute energy error.
		procedure, pass(self) :: compute_max_divergence !< Compute divergence of D, B and J fields for diagnostics.
      procedure, pass(self) :: impose_ct_correction 	!< Impose Constrained Transport correction on q(ivar:ivar+2).
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
   class(prism_fnl_object), intent(in)             :: self                                                   !< The equation.
   integer(I4P),            intent(in)             :: ivar                                                   !< Start index of field of q.
   integer(I4P),            intent(in)             :: ovar                                                   !< Output index in divergence.
   real(R8P),               intent(in)             :: q_gpu(1:,      1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout)          :: divergence_gpu(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Divergence.
   real(R8P),               intent(out), optional  :: maxdiv                                                 !< Maximum divergence for diagnostics.
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

   subroutine compute_residuals_interface_dev(self, q_gpu, dq_gpu, s)
   !< Compute residuals of equation, space operator.
   import :: prism_fnl_object, R8P, I4P
   class(prism_fnl_object), intent(inout) :: self                                              !< The equation.
   real(R8P),               intent(inout) :: q_gpu( 1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Conservative variables.
   real(R8P),               intent(inout) :: dq_gpu(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Residuals.
   integer(I4P),  optional, intent(in)    :: s                                                 !< Stage counter.
   endsubroutine compute_residuals_interface_dev

   subroutine integrate_interface_dev(self)
   !< Integrate equation, time operator.
   import :: prism_fnl_object, R8P
   class(prism_fnl_object), intent(inout) :: self !< The equation.
   endsubroutine integrate_interface_dev
endinterface

contains
   ! auxiliary methods
   subroutine allocate_gpu(self)
   !< Allocate GPU data.
   class(prism_fnl_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: ierr !< Error status.

   call mpih_fnl%print_message('prism_fnl_object%allocate_gpu start')
   associate(nv=>self%nv, nb=>self%nb, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk)
   call dev_alloc(fptr_dev=self%q_gpu, &
                  ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,nv], lbounds=[1,1-ngc,1-ngc,1-ngc,1], init_value=0._R8P, ierr=ierr)
   call dev_alloc(fptr_dev=self%dq_gpu, &
                  ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,nv], lbounds=[1,1-ngc,1-ngc,1-ngc,1], init_value=0._R8P, ierr=ierr)
   call dev_alloc(fptr_dev=self%flxyz_c_gpu, &
                  ubounds=[nb,3,3,ni+ngc,nj+ngc,nk+ngc,nv], lbounds=[1,1,1,1-ngc,1-ngc,1-ngc,1], init_value=0._R8P, ierr=ierr)
   call dev_alloc(fptr_dev=self%flx_f_gpu, &
                  ubounds=[nb,ni,nj,nk,nv], lbounds=[1,0,1,1,1], init_value=0._R8P, ierr=ierr)
   call dev_alloc(fptr_dev=self%fly_f_gpu, &
                  ubounds=[nb,ni,nj,nk,nv], lbounds=[1,1,0,1,1], init_value=0._R8P, ierr=ierr)
   call dev_alloc(fptr_dev=self%flz_f_gpu, &
                  ubounds=[nb,ni,nj,nk,nv], lbounds=[1,1,1,0,1], init_value=0._R8P, ierr=ierr)
   call dev_alloc(fptr_dev=self%curl_gpu, &
                  ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,nv], lbounds=[1,1-ngc,1-ngc,1-ngc,1], init_value=0._R8P, ierr=ierr)
   call dev_alloc(fptr_dev=self%divergence_gpu, &
                  ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,nv], lbounds=[1,1-ngc,1-ngc,1-ngc,1], init_value=0._R8P, ierr=ierr)
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
   call self%coil_fnl%copy_cpu_gpu(coil=self%coil, grid=self%adam%grid)
   call self%fwlayer_fnl%copy_cpu_gpu(fwlayer=self%fWLayer, grid=self%adam%grid, buffer=self%buf_5D_R8P, verbose=verbose)
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
   call self%coil_fnl%copy_gpu_cpu(coil=self%coil, grid=self%adam%grid)
   call self%fwlayer_fnl%copy_gpu_cpu(fwlayer=self%fWLayer, grid=self%adam%grid, buffer=self%buf_5D_R8P, verbose=verbose)
   endsubroutine copy_gpu_cpu

   subroutine initialize_prism(self, filename, realms_number)
   !< Initialize PRISM equation.
   class(prism_fnl_object), intent(inout), target :: self                !< The equation.
   character(*),            intent(in)            :: filename            !< Input file name.
   integer(I4P),            intent(in), optional  :: realms_number       !< Forest realm count; divides the per-device budget (default 1).
   logical                                        :: is_mpih_initialized !< Flag to check if MPI has been inizialied.
   real(R8P)                                      :: memory_avail_       !< Per-realm device budget (GB) after the forest split.
   integer(I4P)                                   :: realms_number_      !< Local realm count (>=1).

   realms_number_ = 1_I4P ; if (present(realms_number)) realms_number_ = max(realms_number, 1_I4P)
   ! mpih_fnl is a program-scope singleton (one rank/device/communicator per process),
   ! not a per-realm component. Initialize it exactly once: FUNDAL's initialize is
   ! intent(out) and runs device init + communicator split, so a second call would
   ! wipe the handler and re-bind the GPU. The first realm initializes it; later
   ! realms skip via mpih_fnl_is_initialized.
   ! @NOTE: the MPI_INITIALIZED probe works around fundal not auto-checking MPI init;
   ! fundal should grow that auto-check as the adam CPU mpi handler already has.
   if (.not. mpih_fnl_is_initialized) then
      call MPI_INITIALIZED(is_mpih_initialized, mpih_fnl%error)
      call mpih_fnl%initialize(do_mpi_init=.not.is_mpih_initialized, do_device_init=.true., verbose=.true.)
      mpih_fnl_is_initialized = .true.
   endif
   call mpih_fnl%print_message('prism_fnl_object%initialize start')
   memory_avail_ = real(mpih_fnl%dev_memory_avail/1e9, R8P) / real(realms_number_, R8P)
   call self%prism_common_object%initialize(filename=filename, memory_avail=memory_avail_, verbose=.true.)
   call self%field_fnl%initialize(grid=self%adam%grid, field=self%adam%field, maps=self%adam%maps, verbose=.true.)
   call self%ib_fnl%initialize(grid=self%adam%grid, field=self%adam%field, ib=self%ib)
   call self%rk_fnl%initialize(grid=self%adam%grid, field=self%adam%field, rk=self%rk)
   call self%weno_fnl%initialize(weno=self%weno)
   call self%allocate_gpu
   call self%coil_fnl%initialize(coil=self%coil, field=self%adam%field, grid=self%adam%grid)
   call self%fwlayer_fnl%initialize(fwlayer=self%fWLayer, field=self%adam%field, grid=self%adam%grid)

   ! set pointer (abstract) TBP
   if (self%physics%physical_model == EM_PHYSICAL_MODEL) then
      select case(self%numerics%scheme_time)
      case(NUM_SCHEME_TIME_BLANES_MOAN)
         self%integrate_dev => integrate_blanesmoan_dev
      case(NUM_SCHEME_TIME_CFM)
         self%integrate_dev => integrate_cfm_dev
      case(NUM_SCHEME_TIME_LEAPFROG)
         self%integrate_dev => integrate_leapfrog_dev
      case(NUM_SCHEME_TIME_RUNGE_KUTTA)
         select case(self%rk%scheme)
         case(RK_1, RK_2, RK_3)
            self%integrate_dev => integrate_rk_ls_dev
         case(RK_SSP_22, RK_SSP_33, RK_SSP_54)
            self%integrate_dev => integrate_rk_ssp_dev
         case(RK_YOSHIDA)
            self%integrate_dev => integrate_rk_yoshida_dev
         endselect
      endselect
   !elseif (self%physics%physical_model == PIC_PHYSICAL_MODEL) then
   !   select case(self%numerics%scheme_time)
   !   case(NUM_SCHEME_TIME_LEAPFROG)
   !      select case(self%pic%scheme_time)
   !      case(NUM_SCHEME_TIME_PIC_LEAPFROG)
   !         self%integrate => integrate_leapfrog_pic
   !      case(NUM_SCHEME_TIME_PIC_RUNGE_KUTTA)
   !         !self%integrate =>
   !      endselect
   !   case(NUM_SCHEME_TIME_RUNGE_KUTTA)
   !      select case(self%pic%scheme_time)
   !      case(NUM_SCHEME_TIME_PIC_LEAPFROG)
   !         self%integrate => integrate_leapfrog_pic
   !      case(NUM_SCHEME_TIME_PIC_RUNGE_KUTTA)
   !         !self%integrate =>
   !      endselect
   !   endselect
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
      ! self%compute_residuals_dev   => compute_residuals_fv_centered_dev
   endselect

   call external_fields_initialize_dev(external_fields=self%external_fields)

   call mpih_fnl%print_message('prism_fnl_object%initialize finish')
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
                                                               blocks_number=self%blocks_number, residuals=self%adam%field%residuals)
   endif
   endsubroutine save_residuals

   subroutine save_simulation_data(self)
   !< Save all simulation data.
   class(prism_fnl_object), intent(inout) :: self      !< The equation.

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
   endsubroutine save_simulation_data

   ! IC/BC/sources
   subroutine apply_fwl_correction(self, q_gpu)
   !< Apply correction if a fWL is present.
   class(prism_fnl_object), intent(inout) :: self                    !< The equation.
   real(R8P),               intent(inout) :: q_gpu(1:,         &
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1:)               !< Conservative variables.
   integer(I4P)                           :: ni(2,6),nj(2,6),nk(2,6) !< Dimensions of FWL domain.
   real(R8P)                              :: s2(6)                   !< Side coefficient.
   integer(I4P)                           :: n(6)                    !< FWL f function index.
   integer(I4P)                           :: alfa_D(6), beta_D(6)    !< Corrected var index of D (Barbas' notation).
   integer(I4P)                           :: alfa_B(6), beta_B(6)    !< Corrected var index of D (Barbas' notation).
   integer(I4P)                           :: face                    !< Counter.

   associate(C=>self%fWLayer%C, layer=>self%fWLayer%layer)
   if (C>0) then
      ! below arrays should be initialized elsewhere...
        ni(1,1)=1_I4P   ;   ni(1,2)=self%ni-C ;   ni(1,3)=1_I4P   ;   ni(1,4)=1_I4P     ;   ni(1,5)=1_I4P   ;   ni(1,6)=1_I4P
        ni(2,1)=C       ;   ni(2,2)=self%ni   ;   ni(2,3)=self%ni ;   ni(2,4)=self%ni   ;   ni(2,5)=self%ni ;   ni(2,6)=self%ni
        nj(1,1)=1_I4P   ;   nj(1,2)=1_I4P     ;   nj(1,3)=1_I4P   ;   nj(1,4)=self%nj-C ;   nj(1,5)=1_I4P   ;   nj(1,6)=1_I4P
        nj(2,1)=self%nj ;   nj(2,2)=self%nj   ;   nj(2,3)=C       ;   nj(2,4)=self%nj   ;   nj(2,5)=self%nj ;   nj(2,6)=self%nj
        nk(1,1)=1_I4P   ;   nk(1,2)=1_I4P     ;   nk(1,3)=1_I4P   ;   nk(1,4)=1_I4P     ;   nk(1,5)=1_I4P   ;   nk(1,6)=self%nk-C
        nk(2,1)=self%nk ;   nk(2,2)=self%nk   ;   nk(2,3)=self%nk ;   nk(2,4)=self%nk   ;   nk(2,5)=C       ;   nk(2,6)=self%nk
           n(1)=1_I4P   ;      n(2)= 1_I4P    ;      n(3)=2_I4P   ;      n(4)= 2_I4P    ;      n(5)=3_I4P   ;      n(6)= 3_I4P
          s2(1)=1.0_R8P ;     s2(2)=-1.0_R8P  ;     s2(3)=1.0_R8P ;     s2(4)=-1.0_R8P  ;     s2(5)=1.0_R8P ;     s2(6)=-1.0_R8P
      alfa_D(1)=2_I4P   ; alfa_D(2)= 2_I4P    ; alfa_D(3)=3_I4P   ; alfa_D(4)= 3_I4P    ; alfa_D(5)=1_I4P   ; alfa_D(6)= 1_I4P
      beta_D(1)=3_I4P   ; beta_D(2)= 3_I4P    ; beta_D(3)=1_I4P   ; beta_D(4)= 1_I4P    ; beta_D(5)=2_I4P   ; beta_D(6)= 2_I4P
      alfa_B(1)=5_I4P   ; alfa_B(2)= 5_I4P    ; alfa_B(3)=6_I4P   ; alfa_B(4)= 6_I4P    ; alfa_B(5)=4_I4P   ; alfa_B(6)= 4_I4P
      beta_B(1)=6_I4P   ; beta_B(2)= 6_I4P    ; beta_B(3)=4_I4P   ; beta_B(4)= 4_I4P    ; beta_B(5)=5_I4P   ; beta_B(6)= 5_I4P
      do face=1, 6
         if (layer(face)) call apply_fwl_correction_dev_kernel(blocks_number=self%blocks_number,ngc=self%ngc,&
                                                               ni1   =  ni(1,face),                          &
                                                               ni2   =  ni(2,face),                          &
                                                               nj1   =  nj(1,face),                          &
                                                               nj2   =  nj(2,face),                          &
                                                               nk1   =  nk(1,face),                          &
                                                               nk2   =  nk(2,face),                          &
                                                               n     =     n(face),                          &
                                                               s2    =    s2(face),                          &
                                                               alfa_D=alfa_D(face),                          &
                                                               beta_D=beta_D(face),                          &
                                                               alfa_B=alfa_B(face),                          &
                                                               beta_B=beta_B(face),                          &
                                                               f_gpu=self%fwlayer_fnl%f_gpu,q_gpu=q_gpu)
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
                                     var_jx        = var_jx       ,&
                                     var_jy        = var_jy       ,&
                                     var_jz        = var_jz       ,&
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
      subroutine nullify_j_vec_vars_kernel(ni,nj,nk,ngc,blocks_number,var_jx,var_jy,var_jz,q_gpu)
      !< Nullify J_Vec vars in q, devide kernel.
      integer(I4P), intent(in)    :: ni,nj,nk,ngc,blocks_number        !< Grids dimensions.
      integer(I4P), intent(in)    :: var_jx,var_jy,var_jz              !< Indexes of J_vec variables.
      real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field cell centered variables.
      integer(I4P)                :: i,j,k,b                           !< Counter.

      !$acc parallel loop independent gang vector collapse(4) &
      !$acc DEVICEVAR(q_gpu)                                  &
      !$acc firstprivate(ni,nj,nk,blocks_number,var_jx,var_jy,var_jz)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         q_gpu(b,i,j,k,var_Jx) = 0._R8P
         q_gpu(b,i,j,k,var_Jy) = 0._R8P
         q_gpu(b,i,j,k,var_Jz) = 0._R8P
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

      !$acc parallel loop independent gang vector collapse(4) &
      !$acc DEVICEVAR(q_gpu,j_vec_gpu)                        &
      !$acc firstprivate(ni,nj,nk,blocks_number,current_density,n,var_jx,var_jy,var_jz)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         q_gpu(b,i,j,k,var_Jx) = q_gpu(b,i,j,k,var_Jx) + current_density * j_vec_gpu(b,i,j,k,1,n)
         q_gpu(b,i,j,k,var_Jy) = q_gpu(b,i,j,k,var_Jy) + current_density * j_vec_gpu(b,i,j,k,2,n)
         q_gpu(b,i,j,k,var_Jz) = q_gpu(b,i,j,k,var_Jz) + current_density * j_vec_gpu(b,i,j,k,3,n)
      enddo
      enddo
      enddo
      enddo
      endsubroutine apply_j_vec_kernel
   endsubroutine compute_coils_current

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
                                             crown                  = crown                                     ,&
                                             local_map_bc_crown_gpu = self%field_fnl%maps%local_map_bc_crown_gpu,&
                                             q_gpu                  = q_gpu)
      enddo
   endif
   contains
      subroutine set_boundary_conditions_kernel(ni, nj, nk, ngc, nv, crown, local_map_bc_crown_gpu, q_gpu)
      !< Set boundary conditions of equation, kernel device.
      integer(I4P), intent(in)    :: ni,nj,nk,ngc                      !< Grid dimensions.
      integer(I4P), intent(in)    :: nv                                !< Number of conservative variables.
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
      !$acc parallel loop independent gang vector &
      !$acc& DEVICEVAR(local_map_bc_crown_gpu, q_gpu) firstprivate(ni,nj,nk,nv,crown)
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
            fec_1_6 = fec_1_6_array(fec)
            if (bc_type == BC_EXTRAPOLATION) then
               do v=1, nv
                  q_gpu(b,i,j,k,v) = q_gpu(b,i-idelta,j-jdelta,k-kdelta,v)
               enddo
            elseif (bc_type == BC_NEUMANN) then
               do v=1, nv
                  q_gpu(b,i,j,k,v) = q_gpu(b,i+abs(idelta)*(-2*i+1+(idelta+1)*ni),&
                                             j+abs(jdelta)*(-2*j+1+(jdelta+1)*nj),&
                                             k+abs(kdelta)*(-2*k+1+(kdelta+1)*nk),v)
               enddo
            elseif (bc_type == BC_SILVER_MULLER) then
               ! to be impelmented
            elseif (bc_type == BC_DIRICHLET) then
               do v=1, nv
                  q_gpu(b,i,j,k,v) = 0._R8P
               enddo
            elseif (bc_type == BC_PERIOD) then
               ! to be implemented
            endif
         endif
      enddo
      endsubroutine set_boundary_conditions_kernel
   endsubroutine set_boundary_conditions

   subroutine set_initial_conditions(self, is_restart)
   !< Set initial conditions of field.
   class(prism_fnl_object), intent(inout) :: self       !< The equation.
   logical,                 intent(in)    :: is_restart !< Branching sentinel for restart/non restart path.

   if (.not.is_restart) call self%ic%set_initial_conditions(physics=self%physics, field=self%adam%field, grid=self%adam%grid, q=self%q)
   ! if (self%physics%physical_model == PIC_PHYSICAL_MODEL) then
   !    call self%particle_injection%set_particle_initial_injection(field=self%adam%field, pic=pic, q_pic=self%q_pic)
   !    call write_initial_injection_tab(filename='particle_injection.dat', q_pic=self%q_pic, np=self%pic%particle_number)
   !    call write_initial_injection_tab(filename='neighbour_list.dat', q_pic=real(self%pic%neighbour_list,R8P), &
   !                                     np=self%pic%particle_number)
   ! endif
   ! call self%coil%set_coils(physics=physics, field=self%adam%field)

   call self%initialize_coils

   ! if (self%physics%physical_model == PIC_PHYSICAL_MODEL) then
   !    call self%pic%current_weighting(field=self%adam%field, q=self%q, q_pic=self%q_pic, nv=self%nv)
   !    call self%pic%particle_weighting(field=self%adam%field, q=self%q, q_pic=self%q_pic, nv=self%nv)
   !    call self%pic%field_weighting(field=self%adam%field, q=self%q, q_pic=self%q_pic, pic_fields=pic_fields, nv=self%nv)
   ! endif

   call self%copy_cpu_gpu
   endsubroutine set_initial_conditions

   subroutine update_ghost(self, q_gpu, step, s)
   !< Update ghost cells.
   !< If not specified all steps are perfermod, syncronous computation
   class(prism_fnl_object), intent(inout)        :: self            !< The equation.
   real(R8P),               intent(inout)        :: q_gpu(1:,         &
                                                          1-self%ngc:,&
                                                          1-self%ngc:,&
                                                          1-self%ngc:,&
                                                          1:)       !< Conservative variables.
   integer(I4P),            intent(in), optional :: step            !< Step to be perfordmed in asyncronous comp.
   integer(I4P),            intent(in), optional :: s               !< Stage counter.
   logical                                       :: do_local_update !< Flag for triggering local update.
   logical                                       :: do_set_bc       !< Flag for triggering setting bc.

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

   if (do_local_update) call self%field_fnl%update_ghost_local_gpu(q_gpu=q_gpu)
                        call self%field_fnl%update_ghost_mpi_gpu(comm_map_send_ptr_ghost=self%adam%maps%comm_map_send_ptr_ghost, &
                                                                 comm_map_recv_ptr_ghost=self%adam%maps%comm_map_recv_ptr_ghost, &
                                                                 q_gpu=q_gpu, step=step)
   if (associated(forest_realm) .and. allocated(self%adam%maps%inter_realm_neighbors)) then
      call self%exchange_inter_realm_halos_forest(realm=forest_realm)
   endif
   if (do_set_bc)       call self%set_boundary_conditions(q_gpu=q_gpu)
   endsubroutine update_ghost

   subroutine update_rk_ghost(self, dt, phi_gpu)
   !< Update RK ghost cells.
   class(prism_fnl_object), intent(inout)        :: self                 !< RK object.
   real(R8P),               intent(in)           :: dt                   !< Current time step.
   real(R8P),               intent(in), optional :: phi_gpu(1:,          &
                                                            1-self%ngc:, &
                                                            1-self%ngc:, &
                                                            1-self%ngc:, &
                                                            1:)          !< IB distance.
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

      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu,curl_gpu) &
      !$acc& firstprivate(ni,nj,nk,blocks_number,ivar,s1)                                        &
      !$acc& private(qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y)
      do b=1,blocks_number
      do k=1,nk
      do j=1,nj
      do i=1,ni
         !$acc loop seq
         do s=1-s1, 1+s1
            qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,ivar+1)
            qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,ivar+2)
            qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,ivar+0)
            qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,ivar+2)
            qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,ivar+0)
            qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,ivar+1)
         enddo
         call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_gpu(b,1:3),          &
                                           qsx_y=qsx_y,qsx_z=qsx_z,qsy_x=qsy_x,&
                                           qsy_z=qsy_z,qsz_x=qsz_x,qsz_y=qsz_y,&
                                           curl=curl_gpu(b,i,j,k,ivar:))
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
   integer(I4P),            intent(in)            :: ivar                                                      !< Start index of field of q.
   integer(I4P),            intent(in)            :: ovar                                                      !< Output index in div.
   real(R8P),               intent(in)            :: q_gpu(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)          !< Field variables.
   real(R8P),               intent(inout)         :: divergence_gpu(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Divergence.
   real(R8P),               intent(out), optional :: maxdiv                                                    !< Max divergence, for checking.

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
      integer(I4P), intent(in)              :: ivar                                                      !< Start index of variable of q.
      integer(I4P), intent(in)              :: ovar                                                      !< Output index in div.
      integer(I4P), intent(in)              :: s1                                                        !< Half FDV stencil length.
      real(R8P),    intent(in)              :: dxyz_gpu(1:,1:)                                           !< Delta cells GPU [nb,3].
      real(R8P),    intent(in)              :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)                         !< Field cell centered variables.
      real(R8P),    intent(inout)           :: divergence_gpu(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Divergence.
      real(R8P),    intent(out), optional   :: maxdiv                                                    !< Max divergence, for checking.
      integer(I4P)                          :: i,j,k,b,s                                                 !< Counter
      ! rank 1D stencil for computations on device that contiguos memory is mandatory
      real(R8P) :: qsx(1-s1:1+s1) !< X component of vector field over the x stencil.
      real(R8P) :: qsy(1-s1:1+s1) !< Y component of vector field over the y stencil.
      real(R8P) :: qsz(1-s1:1+s1) !< Z component of vector field over the z stencil.
      real(R8P) :: maxdiv_

      maxdiv_ = -huge(1._R8P)
      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu,divergence_gpu) &
      !$acc& firstprivate(ivar,ovar,s1)                                                                &
      !$acc& private(qsx,qsy,qsz) reduction(max: maxdiv_)
      do b=1,blocks_number
      do k=1,nk
      do j=1,nj
      do i=1,ni
         !$acc loop seq
         do s=1-s1, 1+s1
            qsx(s) = q_gpu(b,i+s-1,j    ,k    ,ivar+0)
            qsy(s) = q_gpu(b,i    ,j+s-1,k    ,ivar+1)
            qsz(s) = q_gpu(b,i    ,j    ,k+s-1,ivar+2)
         enddo
         call compute_divergence_fd_centered_dev(s=s1,dxyz=dxyz_gpu(b,1:3), &
                                                 qsx=qsx,qsy=qsy,qsz=qsz,   &
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
   integer(I4P),            intent(in)            :: ivar                                                      !< Start index of field of q.
   integer(I4P),            intent(in)            :: ovar                                                      !< Output index in divergence.
   real(R8P),               intent(in)            :: q_gpu(1:,      1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)    !< Field variables.
   real(R8P),               intent(inout)         :: divergence_gpu(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Divergence.
   real(R8P),               intent(out), optional :: maxdiv                                                    !< Max divergence, for checking.
   integer(I4P)                                   :: i,j,k,b,m                                                 !< Counter.
   real(R8P)                                      :: div_x, div_y, div_z                                       !< Partial derivatives.
   real(R8P)                                      :: q_line(1-4:1+4)                                           !< 1D stencil for reconstruction.
   real(R8P)                                      :: ql, qr                                                    !< Left/right reconstructions.

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
   subroutine compute_residuals_fd_centered_dev(self, q_gpu, dq_gpu, s)
   !< Compute residuals of equation, space operator, centered finite difference schemes.
   class(prism_fnl_object), intent(inout) :: self       !< The equation.
   real(R8P),               intent(inout) :: q_gpu(1:,         &
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1:)  !< Conservative variables.
   real(R8P),               intent(inout) :: dq_gpu(1:,         &
                                                    1-self%ngc:,&
                                                    1-self%ngc:,&
                                                    1-self%ngc:,&
                                                    1:) !< Residuals.
   integer(I4P),  optional, intent(in)    :: s          !< Stage counter.

   if (self%blocks_number > 0) then
      call self%apply_fwl_correction(q_gpu=q_gpu)
      call self%update_ghost(q_gpu=q_gpu, s=s)
      call compute_residuals_fd_centered_dev_kernel(ni            = self%ni                  ,&
                                                    nj            = self%nj                  ,&
                                                    nk            = self%nk                  ,&
                                                    ngc           = self%ngc                 ,&
                                                    blocks_number = self%blocks_number       ,&
                                                    var_jx        = self%physics%var_jx           ,&
                                                    var_jy        = self%physics%var_jy           ,&
                                                    var_jz        = self%physics%var_jz           ,&
                                                    nv_c          = self%physics%nv_c             ,&
                                                    chi           = self%physics%chi              ,&
                                                    s1            = self%fdv_half_stencils(1),&
                                                    dxyz_gpu      = self%field_fnl%dxyz_gpu       ,&
                                                    q_gpu         = q_gpu                    ,&
                                                    dq_gpu        = dq_gpu)
   endif
   contains
      subroutine compute_residuals_fd_centered_dev_kernel(ni, nj, nk, ngc, blocks_number,        &
                                                          var_Jx, var_Jy, var_Jz, nv_c, chi, s1, &
                                                          dxyz_gpu, q_gpu, dq_gpu)
      !< Compute residuals of equation, space operator, centered finite difference schemes, kernel device.
      integer(I4P), intent(in)    :: ni,nj,nk,ngc,blocks_number         !< Grids dimensions.
      integer(I4P), intent(in)    :: var_jx,var_jy,var_jz               !< Indexes of J_vec variables.
      integer(I4P), intent(in)    :: nv_c                               !< Number of conservative variables.
      real(R8P),    intent(in)    :: chi                                !< Hyperbolic correction speed.
      integer(I4P), intent(in)    :: s1                                 !< Half FDV stencil length.
      real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)                    !< Delta cells GPU [nb,3].
      real(R8P),    intent(in)    :: q_gpu( 1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field cell centered variables.
      real(R8P),    intent(inout) :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Residuals.
      integer(I4P)                :: i,j,k,b,s                          !< Counter
      real(R8P)                   :: curlD(3), curlB(3)                 !< Residuals components.
      real(R8P)                   :: divergenceD, divergenceB           !< Divergence for hyperbolic correction.
      real(R8P)   			   	 :: gradphi(3), gradpsi(3) 	         !< Gradient for hyperbolic correction.

      ! rank 1D stencil for curl computations on device that contiguos memory is mandatory
      real(R8P) :: qsx_y(1-s1:1+s1) !< Y component of vector field over the x stencil.
      real(R8P) :: qsx_z(1-s1:1+s1) !< Z component of vector field over the x stencil.
      real(R8P) :: qsy_x(1-s1:1+s1) !< X component of vector field over the y stencil.
      real(R8P) :: qsy_z(1-s1:1+s1) !< Z component of vector field over the y stencil.
      real(R8P) :: qsz_x(1-s1:1+s1) !< X component of vector field over the z stencil.
      real(R8P) :: qsz_y(1-s1:1+s1) !< Y component of vector field over the z stencil.

      ! rank 1D stencil for divergence computations on device that contiguos memory is mandatory
      real(R8P) :: qsx_x(1-s1:1+s1) !< X component of vector field over the x stencil.
      real(R8P) :: qsy_y(1-s1:1+s1) !< Y component of vector field over the y stencil.
      real(R8P) :: qsz_z(1-s1:1+s1) !< Z component of vector field over the z stencil.

      ! rank 1D stencil for gradient computations on device that contiguos memory is mandatory
      real(R8P) :: qs(1-s1:1+s1,1-s1:1+s1,1-s1:1+s1) !< Scalar field over the stencil [1-s:1+s,1-s:1+s,1-s:1+s].

		if (self%numerics%div_corr_var == DIV_CORR_VAR_HYPER .and. &
         self%numerics%constrained_transport_D .and. &
		   .not.self%numerics%constrained_transport_B) then
         ! compute RHS dD/dt = curl(B/MU0) - grad(phi) - J, dB/dt = -curl(D/EPS0), dphi/dt = -ch^2*div(D)

         !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu,dq_gpu) &
         !$acc& firstprivate(var_jx,var_jy,var_jz,nv_c,chi,s1)                                    &
         !$acc& private(curlD,curlB,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y,                          &
         !$acc&         qsx_x,qsy_y,qsz_z,qs,divergenceD,gradphi)
         do b=1,blocks_number
         do k=1,nk
         do j=1,nj
         do i=1,ni

            !$acc loop seq
            do s=1-s1, 1+s1
               qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DY)
               qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DZ)
               qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DX)
               qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DZ)
               qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DX)
               qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DY)
            enddo
            call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_gpu(b,1:3),          &
                                              qsx_y=qsx_y,qsx_z=qsx_z,qsy_x=qsy_x,&
                                              qsy_z=qsy_z,qsz_x=qsz_x,qsz_y=qsz_y,&
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
            call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_gpu(b,1:3),          &
                                              qsx_y=qsx_y,qsx_z=qsx_z,qsy_x=qsy_x,&
                                              qsy_z=qsy_z,qsz_x=qsz_x,qsz_y=qsz_y,&
                                              curl=curlB)

            !$acc loop seq
            do s=1-s1, 1+s1
               qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DX)
               qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DY)
               qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DZ)
            enddo
            call compute_divergence_fd_centered_dev(s=s1,dxyz=dxyz_gpu(b,1:3),       &
                                                    qsx=qsx_x,qsy=qsy_y,qsz=qsz_z,   &
                                                    divergence=divergenceD)

            !$acc loop seq
            do s=1-s1,1+s1
               qs(s,1,1) = q_gpu(b,i+s-1,j,k,nv_c)
               qs(1,s,1) = q_gpu(b,i,j+s-1,k,nv_c)
               qs(1,1,s) = q_gpu(b,i,j,k+s-1,nv_c)
            enddo
            call compute_gradient_fd_centered_dev(s=s1,dxyz=dxyz_gpu(b,1:3),q=qs,gradient=gradphi)

            dq_gpu(b,i,j,k,VAR_DX) =  curlB(1)/MU0 - gradphi(1) - q_gpu(b,i,j,k,var_Jx)
            dq_gpu(b,i,j,k,VAR_DY) =  curlB(2)/MU0 - gradphi(2) - q_gpu(b,i,j,k,var_Jy)
            dq_gpu(b,i,j,k,VAR_DZ) =  curlB(3)/MU0 - gradphi(3) - q_gpu(b,i,j,k,var_Jz)
            dq_gpu(b,i,j,k,VAR_BX) = -curlD(1)/EPS0
            dq_gpu(b,i,j,k,VAR_BY) = -curlD(2)/EPS0
            dq_gpu(b,i,j,k,VAR_BZ) = -curlD(3)/EPS0
            dq_gpu(b,i,j,k,nv_c)	  = -(chi*C0)**2*divergenceD
         enddo
         enddo
         enddo
         enddo

		elseif (self%numerics%div_corr_var == DIV_CORR_VAR_HYPER .and. &
              .not.self%numerics%constrained_transport_D .and.       &
		        self%numerics%constrained_transport_B) then
         ! compute RHS dD/dt = curl(B/MU0) - J, dB/dt = -curl(D/EPS0) -grad(psi), dpsi/dt = -ch^2*div(B)

         !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu,dq_gpu) &
         !$acc& firstprivate(var_jx,var_jy,var_jz,nv_c,chi,s1)                                    &
         !$acc& private(curlD,curlB,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y,                          &
         !$acc&         qsx_x,qsy_y,qsz_z,qs,divergenceB,gradpsi)
         do b=1,blocks_number
         do k=1,nk
         do j=1,nj
         do i=1,ni
            !$acc loop seq
            do s=1-s1, 1+s1
               qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DY)
               qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DZ)
               qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DX)
               qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DZ)
               qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DX)
               qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DY)
            enddo
            call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_gpu(b,1:3),          &
                                              qsx_y=qsx_y,qsx_z=qsx_z,qsy_x=qsy_x,&
                                              qsy_z=qsy_z,qsz_x=qsz_x,qsz_y=qsz_y,&
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
            call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_gpu(b,1:3),          &
                                              qsx_y=qsx_y,qsx_z=qsx_z,qsy_x=qsy_x,&
                                              qsy_z=qsy_z,qsz_x=qsz_x,qsz_y=qsz_y,&
                                              curl=curlB)
            !$acc loop seq
            do s=1-s1, 1+s1
               qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BX)
               qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BY)
               qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BZ)
            enddo
            call compute_divergence_fd_centered_dev(s=s1,dxyz=dxyz_gpu(b,1:3),       &
                                                    qsx=qsx_x,qsy=qsy_y,qsz=qsz_z,   &
                                                    divergence=divergenceB)
            !$acc loop seq
            do s=1-s1,1+s1
               qs(s,1,1) = q_gpu(b,i+s-1,j,k,nv_c)
               qs(1,s,1) = q_gpu(b,i,j+s-1,k,nv_c)
               qs(1,1,s) = q_gpu(b,i,j,k+s-1,nv_c)
            enddo
            call compute_gradient_fd_centered_dev(s=s1,dxyz=dxyz_gpu(b,1:3),q=qs,gradient=gradpsi)
            dq_gpu(b,i,j,k,VAR_DX) =  curlB(1)/MU0 - q_gpu(b,i,j,k,var_Jx)
            dq_gpu(b,i,j,k,VAR_DY) =  curlB(2)/MU0 - q_gpu(b,i,j,k,var_Jy)
            dq_gpu(b,i,j,k,VAR_DZ) =  curlB(3)/MU0 - q_gpu(b,i,j,k,var_Jz)
            dq_gpu(b,i,j,k,VAR_BX) = -curlD(1)/EPS0 - gradpsi(1)
            dq_gpu(b,i,j,k,VAR_BY) = -curlD(2)/EPS0 - gradpsi(2)
            dq_gpu(b,i,j,k,VAR_BZ) = -curlD(3)/EPS0 - gradpsi(3)
            dq_gpu(b,i,j,k,nv_c)	  = -(chi*C0)**2*divergenceB
         enddo
         enddo
         enddo
         enddo

		elseif (self%numerics%div_corr_var == DIV_CORR_VAR_HYPER .and. &
              self%numerics%constrained_transport_D .and. &
		        self%numerics%constrained_transport_B) then
		! compute RHS dD/dt = curl(B/MU0) - grad(phi) - J, dB/dt = -curl(D/EPS0) -grad(psi),
		!             dphi/dt = -ch^2*div(D), dpsi/dt = -ch^2*div(B)

         !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu,dq_gpu) &
         !$acc& firstprivate(var_jx,var_jy,var_jz,nv_c,chi,s1)                                    &
         !$acc& private(curlD,curlB,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y,                          &
         !$acc&         qsx_x,qsy_y,qsz_z,qs,divergenceD,divergenceB,gradphi,gradpsi)
         do b=1,blocks_number
         do k=1,nk
         do j=1,nj
         do i=1,ni
            !$acc loop seq
            do s=1-s1, 1+s1
               qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DY)
               qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DZ)
               qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DX)
               qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DZ)
               qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DX)
               qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DY)
            enddo
            call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_gpu(b,1:3),          &
                                              qsx_y=qsx_y,qsx_z=qsx_z,qsy_x=qsy_x,&
                                              qsy_z=qsy_z,qsz_x=qsz_x,qsz_y=qsz_y,&
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
            call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_gpu(b,1:3),          &
                                              qsx_y=qsx_y,qsx_z=qsx_z,qsy_x=qsy_x,&
                                              qsy_z=qsy_z,qsz_x=qsz_x,qsz_y=qsz_y,&
                                              curl=curlB)
            !$acc loop seq
            do s=1-s1, 1+s1
               qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DX)
               qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DY)
               qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DZ)
            enddo
            call compute_divergence_fd_centered_dev(s=s1,dxyz=dxyz_gpu(b,1:3),       &
                                                    qsx=qsx_x,qsy=qsy_y,qsz=qsz_z,   &
                                                    divergence=divergenceD)
            !$acc loop seq
            do s=1-s1,1+s1
               qs(s,1,1) = q_gpu(b,i+s-1,j,k,nv_c-1_I4P)
               qs(1,s,1) = q_gpu(b,i,j+s-1,k,nv_c-1_I4P)
               qs(1,1,s) = q_gpu(b,i,j,k+s-1,nv_c-1_I4P)
            enddo
            call compute_gradient_fd_centered_dev(s=s1,dxyz=dxyz_gpu(b,1:3),q=qs,gradient=gradphi)
            !$acc loop seq
            do s=1-s1, 1+s1
               qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BX)
               qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BY)
               qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BZ)
            enddo
            call compute_divergence_fd_centered_dev(s=s1,dxyz=dxyz_gpu(b,1:3),       &
                                                    qsx=qsx_x,qsy=qsy_y,qsz=qsz_z,   &
                                                    divergence=divergenceB)
            !$acc loop seq
            do s=1-s1,1+s1
               qs(s,1,1) = q_gpu(b,i+s-1,j,k,nv_c)
               qs(1,s,1) = q_gpu(b,i,j+s-1,k,nv_c)
               qs(1,1,s) = q_gpu(b,i,j,k+s-1,nv_c)
            enddo
            call compute_gradient_fd_centered_dev(s=s1,dxyz=dxyz_gpu(b,1:3),q=qs,gradient=gradpsi)
            dq_gpu(b,i,j,k,VAR_DX)       =  curlB(1)/MU0  - gradphi(1) - q_gpu(b,i,j,k,var_Jx)
            dq_gpu(b,i,j,k,VAR_DY)       =  curlB(2)/MU0  - gradphi(2) - q_gpu(b,i,j,k,var_Jy)
            dq_gpu(b,i,j,k,VAR_DZ)       =  curlB(3)/MU0  - gradphi(3) - q_gpu(b,i,j,k,var_Jz)
            dq_gpu(b,i,j,k,VAR_BX)       = -curlD(1)/EPS0 - gradpsi(1)
            dq_gpu(b,i,j,k,VAR_BY)       = -curlD(2)/EPS0 - gradpsi(2)
            dq_gpu(b,i,j,k,VAR_BZ)       = -curlD(3)/EPS0 - gradpsi(3)
            dq_gpu(b,i,j,k,nv_c-1_I4P)	  = -(chi*C0)**2*divergenceD
            dq_gpu(b,i,j,k,nv_c)	        = -(chi*C0)**2*divergenceB
         enddo
         enddo
         enddo
         enddo
      else
         ! compute RHS dD/dt = curl(B/MU0) - J, dB/dt = -curl(D/EPS0)

         !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu,dq_gpu) &
         !$acc& firstprivate(var_jx,var_jy,var_jz,s1)                                             &
         !$acc& private(curlD,curlB,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y)
         do b=1,blocks_number
         do k=1,nk
         do j=1,nj
         do i=1,ni
            !$acc loop seq
            do s=1-s1, 1+s1
               qsx_y(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DY)
               qsx_z(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DZ)
               qsy_x(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DX)
               qsy_z(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DZ)
               qsz_x(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DX)
               qsz_y(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DY)
            enddo
            call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_gpu(b,1:3),          &
                                              qsx_y=qsx_y,qsx_z=qsx_z,qsy_x=qsy_x,&
                                              qsy_z=qsy_z,qsz_x=qsz_x,qsz_y=qsz_y,&
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
            call compute_curl_fd_centered_dev(s=s1,dxyz=dxyz_gpu(b,1:3),          &
                                              qsx_y=qsx_y,qsx_z=qsx_z,qsy_x=qsy_x,&
                                              qsy_z=qsy_z,qsz_x=qsz_x,qsz_y=qsz_y,&
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
      endsubroutine compute_residuals_fd_centered_dev_kernel
   endsubroutine compute_residuals_fd_centered_dev

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

   !!Qua ci va la chiamata alla subroutine che calcola i resiudi delle particellle dq_pic

   !call self%leapfrog%integrate(dt=self%time%dt, q=self%q, dq=self%dq)

   !!Qua ci va la chiamata alla subroutine che integra aggiornando le velocita e le posizioni delle particelle
   !            !decidi se fare qui all'inizio del tempo successivo l'aggiornamento della neighbour list delle particelle

   !call self%impose_div_free
   !!call self%apply_fWL_correction
   endsubroutine integrate_leapfrog_pic

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
         call self%rk_fnl%compute_stage_ls(grid=self%adam%grid, field=self%adam%field, rk=self%rk, s=s, dt=self%time%dt, phi_gpu=self%ib_fnl%phi_gpu, &
                                           dq_gpu=self%dq_gpu, q_gpu=self%q_gpu)
      else
         call self%rk_fnl%compute_stage_ls(grid=self%adam%grid, field=self%adam%field, rk=self%rk, s=s, dt=self%time%dt, dq_gpu=self%dq_gpu, q_gpu=self%q_gpu)
      endif
   enddo
   call self%impose_div_free
   call self%apply_fwl_correction(q_gpu=self%q_gpu)
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
   do s=1, self%rk%nrk
      if (self%ib%solids_number>0) then
         call self%rk_fnl%compute_stage(grid=self%adam%grid, field=self%adam%field, s=s, dt=self%time%dt, phi_gpu=self%ib_fnl%phi_gpu)
      else
         call self%rk_fnl%compute_stage(grid=self%adam%grid, field=self%adam%field, s=s, dt=self%time%dt)
      endif
      call self%compute_coils_current(q_gpu=self%rk_fnl%q_rk_gpu(:,:,:,:,:,s), gamm=self%rk%gamm(s))
      call self%compute_residuals_dev(q_gpu=self%rk_fnl%q_rk_gpu(:,:,:,:,:,s), dq_gpu=self%dq_gpu, s=s)
      ! if (s==1) call self%save_residuals
      if (self%ib%solids_number>0) then
         call self%rk_fnl%assign_stage(grid=self%adam%grid, field=self%adam%field, s=s, q_gpu=self%dq_gpu, phi_gpu=self%ib_fnl%phi_gpu)
      else
         call self%rk_fnl%assign_stage(grid=self%adam%grid, field=self%adam%field, s=s, q_gpu=self%dq_gpu)
      endif
   enddo
   if (self%ib%solids_number>0) then
      call self%rk_fnl%update_q(grid=self%adam%grid, field=self%adam%field, rk=self%rk, dt=self%time%dt, phi_gpu=self%ib_fnl%phi_gpu, q_gpu=self%q_gpu)
      ! call self%update_rk_ghost(dt=self%time%dt, phi_gpu=ib_fnl%phi_gpu)
   else
      call self%rk_fnl%update_q(grid=self%adam%grid, field=self%adam%field, rk=self%rk, dt=self%time%dt, q_gpu=self%q_gpu)
      ! call self%update_rk_ghost(dt=self%time%dt)
      call self%save_residuals
   endif
   call self%compute_coils_current(q_gpu=self%q_gpu)
   call self%impose_div_free
   ! call self%apply_fwl_correction  ! to be removed, probably
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

   subroutine initialize_forest(self, filename, realm_index, realms_number, memory_avail, nv, verbose)
   !< Initialize this realm from scratch: PRISM init, IC injection (or
   !< restart load), initial ghost update on device, initial diagnostics
   !< dump, IO files open, plus PIC/leapfrog priming if those schemes are
   !< active.
   !<
   !< Invoked by forest%initialize. v1 implementation wraps the legacy
   !< `initialize_prism(filename)` plus the verbatim post-init / pre-loop
   !< block formerly inline in `simulate`. Behavior unchanged.
   !<
   !< `realms_number` (passed by the forest as the realm count) divides the
   !< per-device memory budget so each realm sizes its block storage to its
   !< share, not the full device budget — without it, N realms each grab the
   !< whole VRAM and over-subscribe the GPU by Nx (issue #13 rmf-2realm OOM).
   !< The optional `realm_index`, `memory_avail`, `nv`, `verbose` are accepted;
   !< `memory_avail` stays a door-open placeholder (the budget source is
   !< mpih_fnl, only valid after device init inside initialize_prism).
   class(prism_fnl_object), intent(inout)        :: self          !< The realm.
   character(*),            intent(in)           :: filename      !< Input parameters file name.
   integer(I4P),            intent(in), optional :: realm_index   !< Index of this realm in the forest (Phase D).
   integer(I4P),            intent(in), optional :: realms_number !< Realm count; divides the per-device budget (Phase D).
   real(R8P),               intent(in), optional :: memory_avail  !< Per-process memory budget override (door-open placeholder).
   integer(I4P),            intent(in), optional :: nv            !< Number of field variables override.
   logical,                 intent(in), optional :: verbose       !< Trigger verbose output.
   integer(I4P)                                  :: i             !< Counter.

   if (present(realm_index)) continue
   if (present(memory_avail)) continue
   if (present(nv)) continue
   if (present(verbose)) continue

   call self%initialize_prism(filename=filename, realms_number=realms_number)
   if (self%io%restart) then
      call mpih_fnl%print_message('restart simulation from "'//trim(self%io%restart_basename)//'" files')
      call self%load_restart_files(t=self%time%it, time=self%time%time)
      call self%set_initial_conditions(is_restart=self%io%restart)
      call mpih_fnl%print_message('restart [t, time]: '//trim(str(self%time%it))//', '//trim(str(self%time%time)))
   else
      call mpih_fnl%print_message('impose initial conditions start')
      do i=1, self%ic%amr_iterations
         call mpih_fnl%print_message('  AMR/set IC iteration:'//trim(str(i,.true.)))
         call self%set_initial_conditions(is_restart=self%io%restart)
         !if (ib%solids_number > 0) call self%compute_phi()
         !call self%amr_update
      enddo
      call self%adam%make_comm_local_maps_ghost_bc
      call self%set_initial_conditions(is_restart=self%io%restart)
      self%time%time = 0._R8P
      self%time%it = 0
      call mpih_fnl%print_message('impose initial conditions finish')
   endif
   !if (ib%solids_number > 0) call self%compute_phi()
   ! call self%amr_update
   ! call self%update_ghost(q=self%q) ! Aggiunto da FN
   call self%update_ghost(q_gpu=self%q_gpu)

   call self%save_simulation_data
   call self%compute_energy
   !call self%save_energy_error(is_to_open=.true.)
   call self%save_energy_history(is_to_open=.true.) !Cazzo
   call self%compute_max_divergence
   call self%save_divergence_history(is_to_open=.true.) !Cazzo
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
   endsubroutine initialize_forest

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

   associate(self_unused => self)
   end associate
   done = ((self%time%it_max <= 0).and.(self%time%time >= self%time%time_max)).or.&
          ((self%time%it >= self%time%it_max).and.(self%time%it_max > 0))
   endsubroutine is_done_forest

   subroutine finalize_forest(self)
   !< Shut this realm down: final state dump, close residuals file, finalize
   !< MPI handler.
   !<
   !< Invoked by forest%finalize. v1 implementation is the verbatim post-
   !< loop block formerly inline in `simulate`. Behavior unchanged.
   class(prism_fnl_object), intent(inout) :: self !< The realm.

   !call self%compute_energy_error !Cazzo
   call self%save_simulation_data
   call self%io%close_file_residuals
   ! NB: MPI_FINALIZE is NOT called here — it is process-global and runs once via
   ! forest%finalize -> finalize_mpi_forest after ALL realms finish (issue #13).
   endsubroutine finalize_forest

   subroutine finalize_mpi_forest(self)
   !< Finalize the process-global FNL MPI handler (mpih_fnl).
   !<
   !< Overrides realm_object%finalize_mpi_forest to target the FNL `mpih_fnl`
   !< singleton instead of the CPU `mpih`. Called ONCE by forest%finalize
   !< after every realm's finalize_forest (see the base for the rationale).
   class(prism_fnl_object), intent(inout) :: self !< The realm (carries no MPI state).

   associate(self_unused => self); end associate
   call mpih_fnl%finalize
   endsubroutine finalize_mpi_forest

   subroutine simulate(self, filename)
   !< Perform the simulation: legacy single-realm entry point.
   !<
   !< Preserved (R5 of issue #10) so app developers can still drive a single
   !< realm without setting up a forest. Mirrors the forest's `simulate`
   !< orchestration: initialize → loop {compute_dt → advance → post → done}
   !< → finalize. Per-step bookkeeping (self%time%it increment, time_max dt cap,
   !< time advance, progress print, save_memory_status, AMR-update hook)
   !< now lives inside `advance_one_step_forest` and `post_step_forest`, so
   !< this body becomes a thin orchestration that matches the forest path
   !< bit-for-bit.
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
   call self%finalize_mpi_forest ! finalize_forest no longer finalizes MPI (issue #13); single-realm legacy path
   endsubroutine simulate

   ! numerical methods, miscellanea
   subroutine compute_dt(self)
   !< Compute the global stability-limited dt and store it on `self%time%dt`.
   !<
   !< Body delegates the local computation to compute_local_dt_forest
   !< (orchestrator contract method), then performs the legacy
   !< MPI_ALLREDUCE on MPI_COMM_WORLD for backward compatibility with
   !< simulate. The forest's `compute_global_dt` (Phase A.6 commit 8)
   !< will perform its own reduction, possibly on a per-realm sub-comm;
   !< the redundancy disappears once the legacy compute_dt is retired.
   class(prism_fnl_object), intent(inout) :: self !< The equation.

   call self%compute_local_dt_forest(dt_local=self%time%dt)
   call MPI_ALLREDUCE(MPI_IN_PLACE, self%time%dt, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, mpih_fnl%error)
   endsubroutine compute_dt

   subroutine compute_local_dt_forest(self, dt_local)
   !< Compute this realm's local stability-limited dt (no MPI reduction).
   !<
   !< Invoked by forest%compute_global_dt during the min reduction across
   !< all realms in the forest. The reduction itself is the orchestrator's
   !< job; this method computes only the value local to `self`.
   !<
   !< PRISM-FNL override: CFL on the minimum cell size divided by the
   !< maximum wave speed (light speed for Maxwell), local to this MPI
   !< rank's blocks. Mirrors the previous body of compute_dt minus the
   !< MPI_ALLREDUCE.
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
      do b=1, blocks_number
         dxyz_min = min(dxyz_min, dxyz_gpu(b,1), dxyz_gpu(b,2), dxyz_gpu(b,3))
      enddo
      endsubroutine compute_dxyz_min_kernel
   endsubroutine compute_local_dt_forest

   subroutine advance_one_step_forest(self, dt)
   !< Advance this realm by one full timestep of size `dt`.
   !<
   !< Invoked by forest%evolve_one_step once per realm per timestep. Owns
   !< the integration itself — RK substages, BC application, intra-realm
   !< ghost exchange, divergence cleaning — i.e. everything that turns
   !< `q` at time `t` into `q` at time `t + dt`. Also owns the per-step
   !< bookkeeping the legacy `simulate` loop did inline: self%time%it
   !< increment, time_max cap on `dt`, self%time%time advance, progress
   !< printing. These move here so the forest's evolve_one_step path is
   !< bit-identical to the legacy per-realm simulate loop.
   !<
   !< PRISM-FNL override: thin wrapper around the legacy `integrate_dev`
   !< procedure pointer dispatched by `self%numerics%scheme_time` at init.
   !< `integrate_dev`'s body still reads `self%time%dt`, so the wrapper
   !< propagates the (possibly capped) `dt` into the module-scope
   !< singleton before dispatching, then advances `self%time%time`. Once the
   !< forest takes over time bookkeeping (Phase D), the self%time%* updates
   !< move to the forest and this wrapper reduces to `call self%integrate_dev`.
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

   function nrk_forest(self) result(nrk)
   !< Number of RK substages used by this realm's integrator (FNL).
   !<
   !< See PRISM-CPU's `nrk_forest` for the multi-realm-path rationale. The
   !< PRISM FNL C.3 closure (issue #13 D.4a follow-up) promoted the 7 FNL
   !< singletons (mpih_fnl/field_fnl/ib_fnl/rk_fnl/weno_fnl/coil_fnl/fwlayer_fnl)
   !< to per-realm value components, so the substage state in `rk_fnl%q_rk_gpu`
   !< is no longer shared across realms and the multi-realm path is structurally
   !< viable. Bit-comparability with single-realm `rmf` is the D.4b acceptance
   !< criterion, gated on a workstation with MPI.
   class(prism_fnl_object), intent(in) :: self !< The realm.
   integer(I4P)                        :: nrk  !< Number of substages.

   associate(self_unused => self)
   end associate
   nrk = self%rk%nrk
   endfunction nrk_forest

   subroutine prepare_step_forest(self, dt)
   !< Per-step prologue on the multi-realm path: set dt, external-field prelude,
   !< RK GPU stage-array init. Mirror of `integrate_rk_ssp_dev`'s head.
   !<
   !< The per-step time bookkeeping (`self%time%it`, `self%time%dt`, time_max cap) that
   !< `advance_one_step_forest` does inline on the N=1 fast path moves here on
   !< the multi-realm path. `self%time%time` advance and progress print run in
   !< `finalize_step_forest`, mirroring the legacy SSP ordering.
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
   endsubroutine prepare_step_forest

   subroutine assemble_substage_forest(self, s, nrk, dt)
   !< Assemble RK substage `s` on the multi-realm path (FNL/GPU).
   !<
   !< Mirrors the FIRST two lines of `integrate_rk_ssp_dev`'s substage body:
   !< `rk_fnl%compute_stage(s, dt [, phi_gpu])` then `compute_coils_current(q_gpu=rk_fnl%q_rk_gpu(:,...,s), gamm=rk%gamm(s))`.
   !< The coil-current refresh happens here (substage prep) because it's a
   !< source-term update on `q_rk_gpu(:,...,s)` and is needed BEFORE the
   !< inter-realm exchange in the evaluate pass.
   class(prism_fnl_object), intent(inout) :: self !< The realm.
   integer(I4P),            intent(in)    :: s    !< Substage index (1..nrk).
   integer(I4P),            intent(in)    :: nrk  !< Total number of substages.
   real(R8P),               intent(in)    :: dt   !< Timestep size from the forest (unused; self%time%dt is canonical).

   associate(dt_unused => dt, nrk_unused => nrk)
   end associate
   if (self%ib%solids_number>0) then
      call self%rk_fnl%compute_stage(grid=self%adam%grid, field=self%adam%field, s=s, dt=self%time%dt, phi_gpu=self%ib_fnl%phi_gpu)
   else
      call self%rk_fnl%compute_stage(grid=self%adam%grid, field=self%adam%field, s=s, dt=self%time%dt)
   endif
   call self%compute_coils_current(q_gpu=self%rk_fnl%q_rk_gpu(:,:,:,:,:,s), gamm=self%rk%gamm(s))
   endsubroutine assemble_substage_forest

   subroutine residuals_substage_forest(self, s, nrk, dt)
   !< Compute residuals for RK substage `s` (FNL multi-realm path,
   !< 3-phase loop). Reads q_rk_gpu(s), writes dq_gpu. NO assign_stage.
   !< See [[prism_cpu_object]]%`residuals_substage_forest` for the
   !< BC_SEAM-coherence rationale of the residuals/assign split.
   class(prism_fnl_object), intent(inout) :: self !< The realm.
   integer(I4P),            intent(in)    :: s    !< Substage index (1..nrk).
   integer(I4P),            intent(in)    :: nrk  !< Total number of substages.
   real(R8P),               intent(in)    :: dt   !< Timestep size from the forest (unused; self%time%dt is canonical).

   associate(dt_unused => dt, nrk_unused => nrk)
   end associate
   call self%compute_residuals_dev(q_gpu=self%rk_fnl%q_rk_gpu(:,:,:,:,:,s), dq_gpu=self%dq_gpu, s=s)
   endsubroutine residuals_substage_forest

   subroutine assign_substage_forest(self, s, nrk, dt)
   !< Assign RK substage `s` state from dq_gpu (FNL multi-realm path,
   !< 3-phase loop). Overwrites q_rk_gpu(:, interior, :, b, s) with
   !< dq_gpu. Runs AFTER all peers have computed residuals.
   class(prism_fnl_object), intent(inout) :: self !< The realm.
   integer(I4P),            intent(in)    :: s    !< Substage index (1..nrk).
   integer(I4P),            intent(in)    :: nrk  !< Total number of substages.
   real(R8P),               intent(in)    :: dt   !< Timestep size from the forest (unused; self%time%dt is canonical).

   associate(dt_unused => dt, nrk_unused => nrk)
   end associate
   if (self%ib%solids_number>0) then
      call self%rk_fnl%assign_stage(grid=self%adam%grid, field=self%adam%field, s=s, q_gpu=self%dq_gpu, phi_gpu=self%ib_fnl%phi_gpu)
   else
      call self%rk_fnl%assign_stage(grid=self%adam%grid, field=self%adam%field, s=s, q_gpu=self%dq_gpu)
   endif
   endsubroutine assign_substage_forest

   subroutine evaluate_substage_forest(self, s, nrk, dt)
   !< [Deprecated] Combined residuals + stage assignment for RK substage `s`
   !< (FNL multi-realm path, legacy 2-phase). The forest now uses the
   !< residuals/assign split; this wrapper is kept for backward-compat.
   class(prism_fnl_object), intent(inout) :: self !< The realm.
   integer(I4P),            intent(in)    :: s    !< Substage index (1..nrk).
   integer(I4P),            intent(in)    :: nrk  !< Total number of substages.
   real(R8P),               intent(in)    :: dt   !< Timestep size from the forest (unused; self%time%dt is canonical).

   call self%residuals_substage_forest(s=s, nrk=nrk, dt=dt)
   call self%assign_substage_forest(s=s, nrk=nrk, dt=dt)
   endsubroutine evaluate_substage_forest

   subroutine finalize_step_forest(self, dt)
   !< Per-step epilogue on the multi-realm path: q assembly, BC, div-clean,
   !< residual save, coil source refresh, time advance, progress print (FNL).
   class(prism_fnl_object), intent(inout) :: self !< The realm.
   real(R8P),               intent(in)    :: dt   !< Timestep size from the forest.

   associate(dt_unused => dt)
   end associate
   if (self%ib%solids_number>0) then
      call self%rk_fnl%update_q(grid=self%adam%grid, field=self%adam%field, rk=self%rk, dt=self%time%dt, phi_gpu=self%ib_fnl%phi_gpu, q_gpu=self%q_gpu)
   else
      call self%rk_fnl%update_q(grid=self%adam%grid, field=self%adam%field, rk=self%rk, dt=self%time%dt, q_gpu=self%q_gpu)
      call self%save_residuals
   endif
   call self%compute_coils_current(q_gpu=self%q_gpu)
   call self%impose_div_free
   if (self%external_fields%ef_type/=EF_TYPE_NONE) &
      call add_external_fields_dev(external_fields=self%external_fields, field_gpu=self%field_fnl, &
                                   dt=self%time%dt, time=self%time%time, q_gpu=self%q_gpu)
   self%time%time = self%time%time + self%time%dt
   call self%time%print_progress(nodes_number=self%adam%tree%nodes_number)
   endsubroutine finalize_step_forest

   subroutine exchange_inter_realm_halos_forest(self, realm)
   !< Refresh `self`'s ghost cells from peer realms on GPU buffers (FNL).
   !<
   !< FNL twin of [[prism_cpu_object]]%`exchange_inter_realm_halos_forest`.
   !< Algorithm structure is identical (geometric block matching, mirror
   !< copy at the face slab); the only difference is the buffer location —
   !< host arrays on CPU vs device pointers on FNL.
   !<
   !< Substage selection: when `forest_active_substage > 0` the buffer copied
   !< is `rk_fnl%q_rk_gpu(:,:,:,:,:,s)`; outside any substage (once-per-step
   !< floor) it is `self%q_gpu`. This branch is taken inside `copy_peer_face_gpu`
   !< per face. The PRISM FNL C.3 closure (issue #13 D.4a follow-up) made
   !< `rk_fnl` a per-realm value component so the per-substage path is now
   !< well-defined; bit-comparability with single-realm `rmf` is the D.4b
   !< acceptance criterion (blocked on workstation MPI today).
   class(prism_fnl_object), intent(inout) :: self     !< This realm.
   class(realm_object),     intent(in)    :: realm(:) !< All realms in the forest.
   integer(I4P)                           :: n        !< Neighbour entry counter.
   integer(I4P)                           :: peer     !< Peer realm index.
   integer(I4P)                           :: my_face  !< Self face code.
   integer(I4P)                           :: peer_face !< Peer face code.
   integer(I4P)                           :: my_axis    !< Coupling axis (1=x, 2=y, 3=z).
   integer(I4P)                           :: my_sign    !< +1 (MAX face) or -1 (MIN face).
   integer(I4P)                           :: b          !< Self block index.
   integer(I4P)                           :: b_peer     !< Peer block index.
   logical                                :: found      !< Peer block resolution flag.

   if (.not. allocated(self%adam%maps%inter_realm_neighbors)) return
   do n = 1_I4P, int(size(self%adam%maps%inter_realm_neighbors), I4P)
      associate(entry => self%adam%maps%inter_realm_neighbors(n))
         if (entry%coupling /= COUPLING_MIRROR) &
            call mpih%error_stop(msg='prism_fnl_object%exchange_inter_realm_halos_forest: only COUPLING_MIRROR is implemented')
         peer      = entry%peer_realm
         my_face   = entry%my_face
         peer_face = entry%peer_face
         call face_axis_sign(my_face, my_axis, my_sign)
         do b = 1_I4P, int(self%adam%field%blocks_number, I4P)
            if (.not. block_face_on_realm_boundary(b=b, axis=my_axis, sgn=my_sign)) cycle
            call find_peer_block(peer=peer, peer_face=peer_face, my_block=b, &
                                  my_axis=my_axis, my_sign=my_sign,           &
                                  b_peer=b_peer, found=found, realm=realm)
            if (.not. found) cycle
            call copy_peer_face_gpu(b=b, my_axis=my_axis, my_sign=my_sign, &
                                    b_peer=b_peer, peer=peer, realm=realm)
         enddo
      end associate
   enddo
   contains
      pure subroutine face_axis_sign(face_code, axis, sgn)
      integer(I4P), intent(in)  :: face_code
      integer(I4P), intent(out) :: axis
      integer(I4P), intent(out) :: sgn

      select case (face_code)
      case (FACE_X_MAX); axis = 1_I4P; sgn = +1_I4P
      case (FACE_X_MIN); axis = 1_I4P; sgn = -1_I4P
      case (FACE_Y_MAX); axis = 2_I4P; sgn = +1_I4P
      case (FACE_Y_MIN); axis = 2_I4P; sgn = -1_I4P
      case (FACE_Z_MAX); axis = 3_I4P; sgn = +1_I4P
      case (FACE_Z_MIN); axis = 3_I4P; sgn = -1_I4P
      case default;      axis = 0_I4P; sgn = 0_I4P
      end select
      endsubroutine face_axis_sign

      function block_face_on_realm_boundary(b, axis, sgn) result(yes)
      integer(I4P), intent(in) :: b
      integer(I4P), intent(in) :: axis
      integer(I4P), intent(in) :: sgn
      logical                  :: yes
      real(R8P)                :: face_coord, target_coord, tol

      if (sgn > 0_I4P) then
         face_coord   = self%adam%field%emax(axis, b)
         target_coord = self%adam%grid%domain_emax(axis)
      else
         face_coord   = self%adam%field%emin(axis, b)
         target_coord = self%adam%grid%domain_emin(axis)
      endif
      tol = max(abs(target_coord), 1._R8P) * 1.0e-10_R8P
      yes = abs(face_coord - target_coord) <= tol
      endfunction block_face_on_realm_boundary

      subroutine find_peer_block(peer, peer_face, my_block, my_axis, my_sign, b_peer, found, realm)
      integer(I4P),        intent(in)  :: peer
      integer(I4P),        intent(in)  :: peer_face
      integer(I4P),        intent(in)  :: my_block
      integer(I4P),        intent(in)  :: my_axis
      integer(I4P),        intent(in)  :: my_sign
      integer(I4P),        intent(out) :: b_peer
      logical,             intent(out) :: found
      class(realm_object), intent(in)  :: realm(:)
      integer(I4P)                     :: bp, peer_axis, peer_sign
      integer(I4P)                     :: tax1, tax2
      real(R8P)                        :: my_face_coord
      real(R8P)                        :: my_tmin(2), my_tmax(2)
      real(R8P)                        :: peer_face_coord
      real(R8P)                        :: tol

      found  = .false.
      b_peer = 0_I4P
      call face_axis_sign(peer_face, peer_axis, peer_sign)
      if (peer_axis /= my_axis) return
      if (peer_sign == my_sign) return
      tax1 = mod(my_axis,        3_I4P) + 1_I4P
      tax2 = mod(my_axis + 1_I4P, 3_I4P) + 1_I4P
      if (my_sign > 0_I4P) then
         my_face_coord = self%adam%field%emax(my_axis, my_block)
      else
         my_face_coord = self%adam%field%emin(my_axis, my_block)
      endif
      my_tmin(1) = self%adam%field%emin(tax1, my_block)
      my_tmax(1) = self%adam%field%emax(tax1, my_block)
      my_tmin(2) = self%adam%field%emin(tax2, my_block)
      my_tmax(2) = self%adam%field%emax(tax2, my_block)
      tol = max(abs(my_face_coord), 1._R8P) * 1.0e-10_R8P
      select type (peer_realm => realm(peer))
      class is (prism_fnl_object)
         do bp = 1_I4P, int(peer_realm%adam%field%blocks_number, I4P)
            if (peer_sign > 0_I4P) then
               peer_face_coord = peer_realm%adam%field%emax(peer_axis, bp)
            else
               peer_face_coord = peer_realm%adam%field%emin(peer_axis, bp)
            endif
            if (abs(peer_face_coord - my_face_coord) > tol) cycle
            if (abs(peer_realm%adam%field%emin(tax1, bp) - my_tmin(1)) > tol) cycle
            if (abs(peer_realm%adam%field%emax(tax1, bp) - my_tmax(1)) > tol) cycle
            if (abs(peer_realm%adam%field%emin(tax2, bp) - my_tmin(2)) > tol) cycle
            if (abs(peer_realm%adam%field%emax(tax2, bp) - my_tmax(2)) > tol) cycle
            b_peer = bp
            found  = .true.
            return
         enddo
      class default
         call mpih%error_stop(msg='prism_fnl_object%exchange_inter_realm_halos_forest: peer realm is not prism_fnl_object')
      end select
      endsubroutine find_peer_block

      subroutine copy_peer_face_gpu(b, my_axis, my_sign, b_peer, peer, realm)
      !< Host-staged peer→self ghost-slab copy for FNL.
      !<
      !< Pulls the peer's interior slab from device → host, writes into
      !< self's ghost slab on host, pushes back to device. The buffer copied
      !< depends on the active substage: when `forest_active_substage > 0`,
      !< both sides are at substage `s_active` and the buffer is
      !< `rk_fnl%q_rk_gpu(:,:,:,:,:,s_active)`; outside the substage loop
      !< the canonical state is `self%q_gpu` / `peer_realm%q_gpu`. This is
      !< the GPU twin of CPU's per-substage selection in `copy_peer_face`.
      integer(I4P),        intent(in) :: b
      integer(I4P),        intent(in) :: my_axis
      integer(I4P),        intent(in) :: my_sign
      integer(I4P),        intent(in) :: b_peer
      integer(I4P),        intent(in) :: peer
      class(realm_object), intent(in) :: realm(:)
      integer(I4P)                    :: d, i, j, k, ng, ni_, nj_, nk_
      integer(I4P)                    :: peer_ni, peer_nj, peer_nk
      integer(I4P)                    :: s_active

      ng = self%ngc
      ni_ = self%ni
      nj_ = self%nj
      nk_ = self%nk
      s_active = forest_active_substage
      select type (peer_realm => realm(peer))
      class is (prism_fnl_object)
         peer_ni = peer_realm%ni
         peer_nj = peer_realm%nj
         peer_nk = peer_realm%nk
         select case (my_axis)
         case (1_I4P)
            if (my_sign > 0_I4P) then
               do d = 1_I4P, ng ; do k = 1_I4P, nk_ ; do j = 1_I4P, nj_
                  if (s_active > 0_I4P) then
                     !$acc update host(peer_realm%rk_fnl%q_rk_gpu(:, d, j, k, b_peer, s_active))
                     self%rk_fnl%q_rk_gpu(:, ni_ + d, j, k, b, s_active) = peer_realm%rk_fnl%q_rk_gpu(:, d, j, k, b_peer, s_active)
                     !$acc update device(self%rk_fnl%q_rk_gpu(:, ni_ + d, j, k, b, s_active))
                  else
                     !$acc update host(peer_realm%q_gpu(:, d, j, k, b_peer))
                     self%q_gpu(:, ni_ + d, j, k, b) = peer_realm%q_gpu(:, d, j, k, b_peer)
                     !$acc update device(self%q_gpu(:, ni_ + d, j, k, b))
                  endif
               enddo ; enddo ; enddo
            else
               do d = 1_I4P, ng ; do k = 1_I4P, nk_ ; do j = 1_I4P, nj_
                  if (s_active > 0_I4P) then
                     !$acc update host(peer_realm%rk_fnl%q_rk_gpu(:, peer_ni - d + 1_I4P, j, k, b_peer, s_active))
                     self%rk_fnl%q_rk_gpu(:, 1_I4P - d, j, k, b, s_active) = &
                        peer_realm%rk_fnl%q_rk_gpu(:, peer_ni - d + 1_I4P, j, k, b_peer, s_active)
                     !$acc update device(self%rk_fnl%q_rk_gpu(:, 1_I4P - d, j, k, b, s_active))
                  else
                     !$acc update host(peer_realm%q_gpu(:, peer_ni - d + 1_I4P, j, k, b_peer))
                     self%q_gpu(:, 1_I4P - d, j, k, b) = peer_realm%q_gpu(:, peer_ni - d + 1_I4P, j, k, b_peer)
                     !$acc update device(self%q_gpu(:, 1_I4P - d, j, k, b))
                  endif
               enddo ; enddo ; enddo
            endif
         case (2_I4P)
            if (my_sign > 0_I4P) then
               do d = 1_I4P, ng ; do k = 1_I4P, nk_ ; do i = 1_I4P, ni_
                  if (s_active > 0_I4P) then
                     !$acc update host(peer_realm%rk_fnl%q_rk_gpu(:, i, d, k, b_peer, s_active))
                     self%rk_fnl%q_rk_gpu(:, i, nj_ + d, k, b, s_active) = peer_realm%rk_fnl%q_rk_gpu(:, i, d, k, b_peer, s_active)
                     !$acc update device(self%rk_fnl%q_rk_gpu(:, i, nj_ + d, k, b, s_active))
                  else
                     !$acc update host(peer_realm%q_gpu(:, i, d, k, b_peer))
                     self%q_gpu(:, i, nj_ + d, k, b) = peer_realm%q_gpu(:, i, d, k, b_peer)
                     !$acc update device(self%q_gpu(:, i, nj_ + d, k, b))
                  endif
               enddo ; enddo ; enddo
            else
               do d = 1_I4P, ng ; do k = 1_I4P, nk_ ; do i = 1_I4P, ni_
                  if (s_active > 0_I4P) then
                     !$acc update host(peer_realm%rk_fnl%q_rk_gpu(:, i, peer_nj - d + 1_I4P, k, b_peer, s_active))
                     self%rk_fnl%q_rk_gpu(:, i, 1_I4P - d, k, b, s_active) = &
                        peer_realm%rk_fnl%q_rk_gpu(:, i, peer_nj - d + 1_I4P, k, b_peer, s_active)
                     !$acc update device(self%rk_fnl%q_rk_gpu(:, i, 1_I4P - d, k, b, s_active))
                  else
                     !$acc update host(peer_realm%q_gpu(:, i, peer_nj - d + 1_I4P, k, b_peer))
                     self%q_gpu(:, i, 1_I4P - d, k, b) = peer_realm%q_gpu(:, i, peer_nj - d + 1_I4P, k, b_peer)
                     !$acc update device(self%q_gpu(:, i, 1_I4P - d, k, b))
                  endif
               enddo ; enddo ; enddo
            endif
         case (3_I4P)
            if (my_sign > 0_I4P) then
               do d = 1_I4P, ng ; do j = 1_I4P, nj_ ; do i = 1_I4P, ni_
                  if (s_active > 0_I4P) then
                     !$acc update host(peer_realm%rk_fnl%q_rk_gpu(:, i, j, d, b_peer, s_active))
                     self%rk_fnl%q_rk_gpu(:, i, j, nk_ + d, b, s_active) = peer_realm%rk_fnl%q_rk_gpu(:, i, j, d, b_peer, s_active)
                     !$acc update device(self%rk_fnl%q_rk_gpu(:, i, j, nk_ + d, b, s_active))
                  else
                     !$acc update host(peer_realm%q_gpu(:, i, j, d, b_peer))
                     self%q_gpu(:, i, j, nk_ + d, b) = peer_realm%q_gpu(:, i, j, d, b_peer)
                     !$acc update device(self%q_gpu(:, i, j, nk_ + d, b))
                  endif
               enddo ; enddo ; enddo
            else
               do d = 1_I4P, ng ; do j = 1_I4P, nj_ ; do i = 1_I4P, ni_
                  if (s_active > 0_I4P) then
                     !$acc update host(peer_realm%rk_fnl%q_rk_gpu(:, i, j, peer_nk - d + 1_I4P, b_peer, s_active))
                     self%rk_fnl%q_rk_gpu(:, i, j, 1_I4P - d, b, s_active) = &
                        peer_realm%rk_fnl%q_rk_gpu(:, i, j, peer_nk - d + 1_I4P, b_peer, s_active)
                     !$acc update device(self%rk_fnl%q_rk_gpu(:, i, j, 1_I4P - d, b, s_active))
                  else
                     !$acc update host(peer_realm%q_gpu(:, i, j, peer_nk - d + 1_I4P, b_peer))
                     self%q_gpu(:, i, j, 1_I4P - d, b) = peer_realm%q_gpu(:, i, j, peer_nk - d + 1_I4P, b_peer)
                     !$acc update device(self%q_gpu(:, i, j, 1_I4P - d, b))
                  endif
               enddo ; enddo ; enddo
            endif
         end select
      class default
         call mpih%error_stop(msg='prism_fnl_object%exchange_inter_realm_halos_forest: peer realm is not prism_fnl_object')
      end select
      endsubroutine copy_peer_face_gpu
   endsubroutine exchange_inter_realm_halos_forest

   subroutine post_step_forest(self, dt, t, it, do_save_state, do_save_residuals, do_save_restart, do_amr)
   !< Run PRISM-FNL's per-timestep post-step work: state IO, energy
   !< diagnostics, max-divergence diagnostics.
   !<
   !< Invoked by forest%post_step. v1 implementation is the verbatim post-
   !< step block formerly inline in `simulate` — every action runs every
   !< step, since today's cadence is enforced inside the save_* routines
   !< themselves (e.g. save_simulation_data honours `io%it_save`). The
   !< `do_*` flags are signature-only for now: when the forest takes over
   !< cadence (Phase A.6 commit 8 / Phase C), the flags will gate the
   !< individual calls. For now they are accepted but unused, preserving
   !< present-day behavior bit-for-bit.
   !<
   !< `dt`, `t`, `it` are not consumed by the current body; they are on
   !< the contract so the forest can supply them once it owns time-state
   !< (today they are still read from the `time` module singleton).
   class(prism_fnl_object), intent(inout)        :: self              !< The realm.
   real(R8P),               intent(in)           :: dt                !< Timestep size just advanced.
   real(R8P),               intent(in)           :: t                 !< Simulation time after the advance.
   integer(I4P),            intent(in)           :: it                !< Iteration index after the advance.
   logical,                 intent(in), optional :: do_save_state     !< Save state output this step.
   logical,                 intent(in), optional :: do_save_residuals !< Save residuals output this step.
   logical,                 intent(in), optional :: do_save_restart   !< Save restart dump this step.
   logical,                 intent(in), optional :: do_amr            !< Run AMR update this step.

   associate(dt_unused => dt, t_unused => t, it_unused => it) ! v1: dt/t/it not consumed yet
   end associate
   if (present(do_save_state)) continue
   if (present(do_save_residuals)) continue
   if (present(do_save_restart)) continue
   if (present(do_amr)) continue

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
   call self%compute_energy
   !call self%save_energy_error !Cazzo
   call self%save_energy_history !Cazzo
   call self%compute_max_divergence
   call self%save_divergence_history !Cazzo
   endsubroutine post_step_forest

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
      integer(I4P), intent(in)  :: ni,nj,nk,ngc,blocks_number             !< Grids dimensions.
      integer(I4P), intent(in)  :: ivar                                   !< Starting position of vector field.
      real(R8P),    intent(in)  :: dxyz_gpu(1:,1:)                        !< Delta cells GPU [nb,3].
      real(R8P),    intent(in)  :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Conservative variables.
      real(R8P),    intent(out) :: coil_power                             !< Coil power of the vector field.
      integer(I4P)              :: i,j,k,b                                !< Counter.

      coil_power = 0.0_R8P
      !$acc parallel loop independent gang vector collapse(4) &
      !$acc& DEVICEVAR(q_gpu,dxyz_gpu) firstprivate(ni,nj,nk,blocks_number,ivar) reduction(+: coil_power)
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

   subroutine compute_max_divergence(self)
   !< Compute maximum divergence.
   class(prism_fnl_object), intent(inout) :: self       !< The equation.
   real(R8P)                              :: max_div(3) !< Maximum divergence.

	call compute_max_divergence_dev_kernel( ni            = self%ni                  ,&
                                           nj            = self%nj                  ,&
                                           nk            = self%nk                  ,&
                                           blocks_number = self%blocks_number       ,&
														 ngc           = self%ngc                 ,&
                                           var_jx        = self%physics%var_jx           ,&
                                           var_jy        = self%physics%var_jy           ,&
                                           var_jz        = self%physics%var_jz           ,&
                                           s1            = self%fdv_half_stencils(1),&
                                           dxyz_gpu      = self%field_fnl%dxyz_gpu       ,&
                                           q_gpu         = self%q_gpu               ,&
                                           max_div       = max_div)

	self%max_divergence_D = max_div(1)
	self%max_divergence_B = max_div(2)
	self%max_divergence_J = max_div(3)
   contains
      subroutine compute_max_divergence_dev_kernel(ni, nj, nk, blocks_number, ngc, &
                                                   var_Jx, var_Jy, var_Jz, s1, 	  &
                                                   dxyz_gpu, q_gpu, max_div)

			!< Compute maximum divergence of D, B and J fields, device kernel.
			integer(I4P), intent(in)  :: ni, nj, nk, blocks_number, ngc        !< Grids dimensions.
			integer(I4P), intent(in)  :: var_Jx, var_Jy, var_Jz                !< Current variables indices.
			integer(I4P), intent(in)  :: s1                                    !< FDV half stencil.
			real(R8P),    intent(in)  :: dxyz_gpu(1:,1:)                       !< Delta cells GPU [nb,3].
			real(R8P),    intent(in)  :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Conservative variables.
			real(R8P),    intent(out) :: max_div(3)                            !< Maximum divergence of D, B and J fields.
			real(R8P)                 :: qsx_x(1-s1:1+s1)                      !< Buffer for x-derivative in x direction (for divergence).
			real(R8P)                 :: qsy_y(1-s1:1+s1)                      !< Buffer for y-derivative in y direction (for divergence).
			real(R8P)                 :: qsz_z(1-s1:1+s1)                      !< Buffer for z-derivative in z direction (for divergence).
			real(R8P)                 :: divergenceD, divergenceB, divergenceJ !< Divergence of D, B and J fields.
			real(R8P)                 :: max_divD, max_divB, max_divJ			 !< Maximum divergence of D, B and J fields.
			integer(I4P)              :: i,j,k,b,s                             !< Counter.

			max_divD = 0.0_R8P
			max_divB = 0.0_R8P
			max_divJ = 0.0_R8P
	      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu) &
         !$acc& firstprivate(var_jx,var_jy,var_jz,s1)                                      &
         !$acc& private(qsx_x,qsy_y,qsz_z,divergenceD,divergenceB,divergenceJ)			 	 &
			!$acc& reduction(max: max_divD, max_divB, max_divJ)
         do b=1,blocks_number
         do k=1,nk
         do j=1,nj
         do i=1,ni
            !$acc loop seq
            do s=1-s1, 1+s1
               qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_DX)
               qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_DY)
               qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_DZ)
            enddo
            call compute_divergence_fd_centered_dev(s=s1,dxyz=dxyz_gpu(b,1:3),       &
                                                    qsx=qsx_x,qsy=qsy_y,qsz=qsz_z,   &
                                                    divergence=divergenceD)
				max_divD = max(max_divD, abs(divergenceD))
            !$acc loop seq
            do s=1-s1, 1+s1
               qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,VAR_BX)
               qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,VAR_BY)
               qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,VAR_BZ)
            enddo
            call compute_divergence_fd_centered_dev(s=s1,dxyz=dxyz_gpu(b,1:3),       &
                                                    qsx=qsx_x,qsy=qsy_y,qsz=qsz_z,   &
                                                    divergence=divergenceB)
				max_divB = max(max_divB, abs(divergenceB))
            !$acc loop seq
            do s=1-s1, 1+s1
               qsx_x(s) = q_gpu(b,i+s-1,j    ,k    ,var_Jx)
               qsy_y(s) = q_gpu(b,i    ,j+s-1,k    ,var_Jy)
               qsz_z(s) = q_gpu(b,i    ,j    ,k+s-1,var_Jz)
            enddo
            call compute_divergence_fd_centered_dev(s=s1,dxyz=dxyz_gpu(b,1:3),       &
                                                    qsx=qsx_x,qsy=qsy_y,qsz=qsz_z,   &
                                                    divergence=divergenceJ)
				max_divJ = max(max_divJ, abs(divergenceJ))
         enddo
         enddo
         enddo
         enddo
			max_div(1) = max_divD
			max_div(2) = max_divB
			max_div(3) = max_divJ
		endsubroutine compute_max_divergence_dev_kernel
	endsubroutine compute_max_divergence
   
   subroutine impose_ct_correction(self, ivar)
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
   endsubroutine impose_ct_correction

   subroutine impose_div_free(self)
   !< Impose divergence-free property.
   class(prism_fnl_object), intent(inout) :: self !< The equation.

   associate(constrained_transport_D=>self%numerics%constrained_transport_D,&
             constrained_transport_B=>self%numerics%constrained_transport_B,div_corr_var=>self%numerics%div_corr_var)
   if (constrained_transport_D.and.div_corr_var==DIV_CORR_VAR_POISS) call self%impose_ct_correction(ivar=1_I4P)
   if (constrained_transport_B.and.div_corr_var==DIV_CORR_VAR_POISS) call self%impose_ct_correction(ivar=4_I4P)
   ! here should go also other corrections...
   endassociate
   endsubroutine impose_div_free
endmodule adam_prism_fnl_object
