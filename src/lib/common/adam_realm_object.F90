!< ADAM, realm class definition — common per-realm state and orchestrator contract.
module adam_realm_object
!< ADAM, realm class definition — common per-realm state and orchestrator contract.
!<
!< A **realm** is the unit the [[forest_object]] orchestrator tends. The name is
!< the crasi of two aspects this type entails simultaneously:
!<
!<   1. A **physical-numerical model**: the PDE system, numerics configuration,
!<      time integrator, boundary conditions, etc. (math layer).
!<   2. A **spatial region modelization**: the grid, the AMR Morton tree, the
!<      field data, the communication maps (where layer).
!<
!< A realm has both **laws** (1) and **borders** (2). Neither name `equation_object`
!< (math-only) nor `subdomain_object` (region-only) captured both aspects — hence
!< `realm_object`, evoking a bounded territory governed by its own rules.
!<
!< Consumer apps (PRISM, NASTO, CHASE, ...) **extend** `realm_object` to define
!< their specific PDE family (`prism_common_object extends realm_object`, etc.).
!< Each extended type adds app-private data and overrides the orchestrator
!< contract TBPs.
!<
!< The orchestrator contract is the set of TBPs carrying the **`_forest`** suffix.
!< These are the methods [[forest_object]] may call on a realm; `grep "_forest"`
!< across the codebase reveals the entire orchestrator contract surface. Three
!< architectural rules govern these TBPs (see issue #10, Step 6 Phase A.2):
!<
!<   * O1 — signature uses only ADAM-lib-visible types (intrinsic kinds plus
!<     `class(realm_object)`, `grid_object`, etc.). Apps override the
!<     implementation; signatures are sacred.
!<   * O2 — app-specific dispatch (which integrator, which physics) lives in
!<     `self%`'s components, not in orchestrator-supplied arguments.
!<   * O3 — `q` is never public on `realm_object`. Halo exchange between realms
!<     (Phase D work) goes through realm-side TBPs that operate on q internally;
!<     the orchestrator schedules but never touches q.

! ADAM classes, libraries, parameters
use :: adam_adam_object
use :: adam_amr_object
use :: adam_blanes_moan_object
use :: adam_cfm_object
use :: adam_eos_ic_object
use :: adam_fdv_operators_library
use :: adam_field_object
use :: adam_flail_object
use :: adam_ib_object
use :: adam_io_object
use :: adam_leapfrog_object
use :: adam_maps_object
use :: adam_parameters
use :: adam_rk_object
use :: adam_riemann_euler_library
use :: adam_slices_object
use :: adam_tree_node_object
use :: adam_tree_bucket_object
use :: adam_tree_object
use :: adam_weno_object
! ADAM singleton objects
use :: adam_globals,    only : adam, field, grid, ib, mpih, rk, weno
use :: adam_adam_bind,  only : bind_globals_to_adam
! third party modules
use :: finer
use :: motion
use :: penf
use :: stringifor

implicit none
private
public :: realm_object

type :: realm_object
   !< Equation system class definition, common data to all backends and applications.
   ! ADAM library objects
   type(io_object)         :: io         !< IO handler.
   type(amr_object)        :: amr        !< AMR marker handler.
   type(slices_object)     :: slices     !< Slices handler.
   type(blanesmoan_object) :: blanesmoan !< Blanes-Moan integrator.
   type(cfm_object)        :: cfm        !< Commutator-Free Magnus integrator.
   type(leapfrog_object)   :: leapfrog   !< Leapfrog integrator.
   type(flail_object)      :: flail      !< Linear algebra methods handler.
   ! Owned sub-objects — Step 2 of forest-of-trees migration (issue #10).
   ! Value components: weno/ib/rk are now first-class state of the equation
   ! and live in the running solver instance (prism, etc.) rather than in
   ! parallel module-scope singletons. The adam_*_global modules become
   ! pointer shims aliased into these components by inline pointer
   ! assignments in initialize (no separate binder module — see note there).
   type(weno_object) :: weno !< WENO reconstructor.
   type(ib_object)   :: ib   !< Immersed boundary.
   type(rk_object)   :: rk   !< Runge-Kutta integrator.
   ! C.3 closure of issue #10 / D.3 follow-up of issue #13: adam was the last
   ! sub-object whose state lived in a module-scope singleton owner (formerly
   ! `adam_singleton`) rather than per realm. Promoted here to a value
   ! component so that for N>1 each realm has its own grid/field/tree/maps
   ! (carried inside its own adam) and inter-realm topology can be stored
   ! per-realm in `realm(is)%adam%maps%inter_realm_neighbors`. The legacy
   ! adam_adam_global module becomes a pointer shim aliasing realm(1)%adam,
   ! matching the seven other shims; the 28 consumer files that read `adam%...`
   ! through the shim continue to work unchanged for the N=1 case.
   type(adam_object) :: adam !< ADAM (grid + tree + field + maps) container.
   ! FDV data
   character(:), allocatable :: fdv_scheme                   !< FDV scheme, fd/fv.
   integer(I4P)              :: fdv_order=2_I4P              !< Order of finite difference/volume schemes, general order.
   integer(I4P)              :: fdv_half_stencil=1_I4P       !< Half stencil length of finite difference/volume schemes.
   integer(I4P)              :: fdv_half_stencils(6)=[1_I4P,&
                                                      1_I4P,&
                                                      1_I4P,&
                                                      1_I4P,&
                                                      1_I4P,&
                                                      1_I4P] !< Half stencil length of fdv schemes for each derivative up to 6.
   ! grid/field data replica for easy handling
   integer(I4P), pointer :: ngc=>null()           !< Number of ghost cells.
   integer(I4P), pointer :: ni=>null()            !< Number of cells in i direction.
   integer(I4P), pointer :: nj=>null()            !< Number of cells in j direction.
   integer(I4P), pointer :: nk=>null()            !< Number of cells in k direction.
   integer(I4P), pointer :: nb=>null()            !< Total blocks number for MPI.
   integer(I4P), pointer :: blocks_number=>null() !< Actual blocks number.
   integer(I4P), pointer :: nv=>null()            !< Number of variables in q vector.
   !< Procedure pointer TBPs for FDV operators (set at initialization by backend).
   procedure(compute_block_total_variation_interface), pass(self),pointer :: compute_block_total_variation=>null()!< Compute TV.
   procedure(compute_curl_interface),                  pass(self),pointer :: compute_curl                 =>null()!< Compute curl.
   procedure(compute_derivative1_interface),           pass(self),pointer :: compute_derivative1          =>null()!< Compute deriv1.
   procedure(compute_derivative2_interface),           pass(self),pointer :: compute_derivative2          =>null()!< Compute deriv2.
   procedure(compute_derivative4_interface),           pass(self),pointer :: compute_derivative4          =>null()!< Compute deriv4.
   procedure(compute_divergence_interface),            pass(self),pointer :: compute_divergence           =>null()!< Compute dive.
   procedure(compute_gradient_interface),              pass(self),pointer :: compute_gradient             =>null()!< Compute grad.
   procedure(compute_laplacian_interface),             pass(self),pointer :: compute_laplacian            =>null()!< Compute laplac.
   contains
      ! public methods
      procedure, pass(self) :: initialize         !< Initialize common data.
      procedure, pass(self) :: load_fdv_from_file !< Load FDV config from file.
      ! Orchestrator contract — methods carrying the `_forest` suffix are
      ! the surface forest_object may call on a realm. See issue #10 Step 6
      ! Phase A.4. Default implementations here error-stop with a
      ! "not overridden" message; each consumer app must override.
      procedure, pass(self) :: initialize_forest                !< Invoked by forest%initialize per realm at startup.
      procedure, pass(self) :: compute_local_dt_forest          !< Invoked by forest%compute_global_dt during the min reduction.
      procedure, pass(self) :: advance_one_step_forest          !< Invoked by forest%evolve_one_step per realm per timestep.
      procedure, pass(self) :: post_step_forest                 !< Invoked by forest%post_step per realm per timestep.
      procedure, pass(self) :: is_done_forest                   !< Invoked by forest%is_done during the termination reduction.
      procedure, pass(self) :: finalize_forest                  !< Invoked by forest%finalize per realm at shutdown.
      procedure, pass(self) :: exchange_inter_realm_halos_forest !< Invoked by forest%exchange_halos to refresh inter-realm ghosts.
      ! IO methods
      procedure, nopass     :: close_block_xh5f !< Close XH5F file block.
      procedure, nopass     :: close_file_xh5f  !< Close XH5F file.
      procedure, pass(self) :: open_block_xh5f  !< Open block file XH5F.
      procedure, pass(self) :: open_file_xh5f   !< Open file XH5F.
      procedure, pass(self) :: save_q_xh5f      !< Save in XH5F (XDMF/HDF5) format.
      ! private FDV operators
      procedure, pass(self), private :: compute_block_total_variation_fd !< Return the max of block total variation for a given var.
      procedure, pass(self), private :: compute_curl_fd                  !< Compute curl of vector field, finite difference.
      procedure, pass(self), private :: compute_curl_fv                  !< Compute curl of vector field, finite volume.
      procedure, pass(self), private :: compute_derivative1_fd           !< Compute derivative1 of scalar field, finite difference.
      procedure, pass(self), private :: compute_derivative1_fv           !< Compute derivative1 of scalar field, finite volume.
      procedure, pass(self), private :: compute_derivative2_fd           !< Compute derivative2 of scalar field, finite difference.
      procedure, pass(self), private :: compute_derivative2_fv           !< Compute derivative2 of scalar field, finite volume.
      procedure, pass(self), private :: compute_derivative4_fd           !< Compute derivative4 of scalar field, finite difference.
      procedure, pass(self), private :: compute_divergence_fd            !< Compute divergence of vector field, finite difference.
      procedure, pass(self), private :: compute_divergence_fv            !< Compute divergence of vector field, finite volume.
      procedure, pass(self), private :: compute_gradient_fd              !< Compute gradient of scalar field, finite difference.
      procedure, pass(self), private :: compute_gradient_fv              !< Compute gradient of scalar field, finite volume.
      procedure, pass(self), private :: compute_laplacian_fd             !< Compute laplacian of scalar field, finite difference.
      procedure, pass(self), private :: compute_laplacian_fv             !< Compute laplacian of scalar field, finite volume.
endtype realm_object

interface
   subroutine compute_block_total_variation_interface(self, hs, dxyz, ivar, q, tot_var_field, total_variation)
   !< Return the max of block total variation for a given var.
   import :: realm_object, I4P, R8P
   class(realm_object), intent(in)     :: self                                                 !< Coils.
   integer(I4P),           intent(in)     :: hs                                                   !< FDV half stencil length.
   real(R8P),              intent(in)     :: dxyz(3)                                              !< Space steps.
   integer(I4P),           intent(in)     :: ivar                                                 !< Index of first component of vec field.
   real(R8P),              intent(in)     :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:)            !< Field variables.
   real(R8P),              intent(inout)  :: tot_var_field(1-self%ngc:,1-self%ngc:,1-self%ngc:)   !< Total variation field on blocks.
   real(R8P),              intent(out)    :: total_variation                                      !< Max total variation on given block.
   endsubroutine compute_block_total_variation_interface

   subroutine compute_curl_interface(self, hs, ivar, q, curl)
   !< Compute curl of vector field, curl(q(ivar:ivar+2)).
   import :: realm_object, I4P, R8P
   class(realm_object), intent(in)    :: self                                            !< The equation.
   integer(I4P),           intent(in)    :: hs                                              !< FDV half stencil length.
   integer(I4P),           intent(in)    :: ivar                                            !< Start index of variable of q.
   real(R8P),              intent(in)    :: q(   1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),              intent(inout) :: curl(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Curl.
   endsubroutine compute_curl_interface

   subroutine compute_derivative1_interface(self, hs, dir, ivar, q, dq_ds)
   !< Compute first derivative of scalar field, dq(ivar)/ds.
   import :: realm_object, I4P, R8P
   class(realm_object), intent(in)    :: self                                             !< The equation.
   integer(I4P),           intent(in)    :: hs                                               !< FDV half stencil length.
   integer(I4P),           intent(in)    :: dir                                              !< Direction: 1=X, 2=Y, 3=Z.
   integer(I4P),           intent(in)    :: ivar                                             !< Index of variable of q.
   real(R8P),              intent(in)    :: q(    1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),              intent(inout) :: dq_ds(   1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< First derivative dq/ds.
   endsubroutine compute_derivative1_interface

   subroutine compute_derivative2_interface(self, hs, dir, ivar, q, d2q_ds2)
   !< Compute second derivative of scalar field, d2q(ivar)/ds2.
   import :: realm_object, I4P, R8P
   class(realm_object), intent(in)    :: self                                               !< The equation.
   integer(I4P),           intent(in)    :: hs                                                 !< FDV half stencil length.
   integer(I4P),           intent(in)    :: dir                                                !< Direction: 1=X, 2=Y, 3=Z.
   integer(I4P),           intent(in)    :: ivar                                               !< Index of variable of q.
   real(R8P),              intent(in)    :: q(      1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),              intent(inout) :: d2q_ds2(   1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Second derivative d2q/ds2.
   endsubroutine compute_derivative2_interface

   subroutine compute_derivative4_interface(self, hs, dir, ivar, q, d4q_ds4)
   !< Compute fourth derivative of scalar field, d4q(ivar)/ds4.
   import :: realm_object, I4P, R8P
   class(realm_object), intent(in)    :: self                                               !< The equation.
   integer(I4P),           intent(in)    :: hs                                                 !< FDV half stencil length.
   integer(I4P),           intent(in)    :: dir                                                !< Direction: 1=X, 2=Y, 3=Z.
   integer(I4P),           intent(in)    :: ivar                                               !< Index of variable of q.
   real(R8P),              intent(in)    :: q(      1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),              intent(inout) :: d4q_ds4(   1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Fourth derivative d4q/ds4.
   endsubroutine compute_derivative4_interface

   subroutine compute_divergence_interface(self, hs, ivar, q, divergence)
   !< Compute divergence of vector field, div(q(ivar:ivar+2)).
   import :: realm_object, I4P, R8P
   class(realm_object), intent(in)    :: self                                                  !< The equation.
   integer(I4P),           intent(in)    :: hs                                                    !< FDV half stencil length.
   integer(I4P),           intent(in)    :: ivar                                                  !< Start index of field of q.
   real(R8P),              intent(in)    :: q(         1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),              intent(inout) :: divergence(   1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Divergence.
   endsubroutine compute_divergence_interface

   subroutine compute_gradient_interface(self, hs, ivar, q, gradient)
   !< Compute gradient of scalar variable q(ivar).
   import :: realm_object, I4P, R8P
   class(realm_object), intent(in)    :: self                                                  !< The equation.
   integer(I4P),           intent(in)    :: hs                                                    !< FDV half stencil length.
   integer(I4P),           intent(in)    :: ivar                                                  !< Index of scalar variable of q.
   real(R8P),              intent(in)    :: q(        1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)  !< Field variables.
   real(R8P),              intent(inout) :: gradient( 1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)  !< Gradient.
   endsubroutine compute_gradient_interface

   subroutine compute_laplacian_interface(self, hs, ivar, q, laplacian)
   !< Compute laplacian of scalar variable q(ivar).
   import :: realm_object, I4P, R8P
   class(realm_object), intent(in)    :: self                                                  !< The equation.
   integer(I4P),           intent(in)    :: hs                                                    !< FDV half stencil length.
   integer(I4P),           intent(in)    :: ivar                                                  !< Index of scalar variable of q.
   real(R8P),              intent(in)    :: q(        1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)  !< Field variables.
   real(R8P),              intent(inout) :: laplacian(   1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)  !< Laplacian.
   endsubroutine compute_laplacian_interface
endinterface

contains
   ! orchestrator contract — see issue #10 Step 6 Phase A.4
   subroutine initialize_forest(self, filename, realm_index, memory_avail, nv, verbose)
   !< Initialize this realm from scratch: app-level initialize, IC injection
   !< (or restart load), initial AMR, IO file open, initial diagnostics dump.
   !<
   !< Invoked by forest%initialize once per realm at startup, before the
   !< time loop begins. This is the single per-realm setup entry point the
   !< orchestrator uses: it both runs the basic per-realm initialize
   !< (which the legacy simulate did via `realm%initialize_prism(filename)`
   !< or equivalent) and performs the post-init / pre-loop setup that has
   !< to happen before `forest%evolve_one_step` can be called.
   !<
   !< Default implementation error-stops: every consumer app MUST override
   !< this method (the initial-condition catalog, AMR strategy, and which
   !< IO files to open are app-specific). PRISM's override lives on
   !< prism_cpu_object/prism_fnl_object.
   !<
   !< Optional `realm_index`, `memory_avail`, `nv`, `verbose` are door-
   !< kept-open for Phase D: for N=1 they are unused, but the signature is
   !< locked in so adding them later does not break callers.
   class(realm_object), intent(inout)        :: self         !< The realm.
   character(*),        intent(in)           :: filename     !< Input parameters file name.
   integer(I4P),        intent(in), optional :: realm_index  !< Index of this realm in the forest (Phase D).
   real(R8P),           intent(in), optional :: memory_avail !< Per-process memory budget override.
   integer(I4P),        intent(in), optional :: nv           !< Number of field variables override.
   logical,             intent(in), optional :: verbose      !< Trigger verbose output.

   associate(filename_unused => filename) ! quiet "unused dummy" warnings before the stop
   end associate
   if (present(realm_index)) continue
   if (present(memory_avail)) continue
   if (present(nv)) continue
   if (present(verbose)) continue
   error stop 'realm_object%initialize_forest: not overridden by app extension'
   endsubroutine initialize_forest

   subroutine compute_local_dt_forest(self, dt_local)
   !< Compute this realm's local stability-limited dt (no MPI reduction).
   !<
   !< Invoked by forest%compute_global_dt during the min reduction across
   !< all realms in the forest. The reduction itself is the orchestrator's
   !< job; this method computes only the value local to `self`.
   !<
   !< Default implementation error-stops: every consumer app MUST override
   !< this method (the CFL criterion and umax are app-specific). PRISM's
   !< override lives on prism_cpu_object/prism_fnl_object.
   class(realm_object), intent(in)  :: self     !< The realm.
   real(R8P),           intent(out) :: dt_local !< Local stability-limited dt.

   dt_local = 0._R8P ! quiet "may be uninitialised" warnings before the stop
   error stop 'realm_object%compute_local_dt_forest: not overridden by app extension'
   endsubroutine compute_local_dt_forest

   subroutine advance_one_step_forest(self, dt)
   !< Advance this realm by one full timestep of size `dt`.
   !<
   !< Invoked by forest%evolve_one_step once per realm per timestep. The
   !< orchestrator owns global dt selection (compute_global_dt) and the
   !< termination check; this method owns the integration itself — RK
   !< substages, BC application, intra-realm ghost exchange, divergence
   !< cleaning — i.e. everything that turns `q` at time `t` into `q` at
   !< time `t + dt`.
   !<
   !< Default implementation error-stops: every consumer app MUST override
   !< this method (the integrator catalog and the substage layout are
   !< app-specific). PRISM's override lives on prism_cpu_object and
   !< prism_fnl_object as a thin wrapper around the legacy `integrate`
   !< / `integrate_dev` dispatch.
   class(realm_object), intent(inout) :: self !< The realm.
   real(R8P),           intent(in)    :: dt   !< Timestep size.

   associate(dt_unused => dt) ! quiet "unused dummy" warnings before the stop
   end associate
   error stop 'realm_object%advance_one_step_forest: not overridden by app extension'
   endsubroutine advance_one_step_forest

   subroutine post_step_forest(self, dt, t, it, do_save_state, do_save_residuals, do_save_restart, do_amr)
   !< Run this realm's per-timestep post-step work: diagnostics, IO, AMR.
   !<
   !< Invoked by forest%post_step once per realm per timestep, after the
   !< advance. The orchestrator owns cadence decisions (when to save state,
   !< when to refine the mesh, etc.) and conveys them via the optional
   !< `do_*` flags; this method executes the work that those flags enable.
   !< Cadence flags absent ⇒ caller wants the default cadence (today: do
   !< everything every step, since save_*_data routines internally respect
   !< the realm's own `it_save` settings).
   !<
   !< Default implementation error-stops: every consumer app MUST override
   !< this method (the diagnostic catalog and IO layout are app-specific).
   !< PRISM's override lives on prism_cpu_object and prism_fnl_object.
   class(realm_object), intent(inout)        :: self              !< The realm.
   real(R8P),           intent(in)           :: dt                !< Timestep size just advanced.
   real(R8P),           intent(in)           :: t                 !< Simulation time after the advance.
   integer(I4P),        intent(in)           :: it                !< Iteration index after the advance.
   logical,             intent(in), optional :: do_save_state     !< Save state output this step.
   logical,             intent(in), optional :: do_save_residuals !< Save residuals output this step.
   logical,             intent(in), optional :: do_save_restart   !< Save restart dump this step.
   logical,             intent(in), optional :: do_amr            !< Run AMR update this step.

   associate(dt_unused => dt, t_unused => t, it_unused => it) ! quiet "unused dummy" warnings
   end associate
   if (present(do_save_state)) continue
   if (present(do_save_residuals)) continue
   if (present(do_save_restart)) continue
   if (present(do_amr)) continue
   error stop 'realm_object%post_step_forest: not overridden by app extension'
   endsubroutine post_step_forest

   subroutine is_done_forest(self, done)
   !< Decide whether this realm's local termination criterion is met.
   !<
   !< Invoked by forest%is_done during the termination reduction across
   !< all realms in the forest. The reduction itself (typically a logical
   !< AND or OR depending on coupling semantics) is the orchestrator's
   !< job; this method computes only the predicate local to `self`.
   !<
   !< Default implementation error-stops: every consumer app MUST override
   !< this method (the termination criterion mixes iteration count,
   !< simulated time, and app-specific physics signals). PRISM's override
   !< lives on prism_cpu_object/prism_fnl_object.
   class(realm_object), intent(in)  :: self !< The realm.
   logical,             intent(out) :: done !< True if this realm is done evolving.

   done = .false. ! quiet "may be uninitialised" warnings before the stop
   associate(self_unused => self)
   end associate
   error stop 'realm_object%is_done_forest: not overridden by app extension'
   endsubroutine is_done_forest

   subroutine finalize_forest(self)
   !< Shut this realm down: final state dump, final diagnostics, close
   !< output files, finalize the realm's MPI handler.
   !<
   !< Invoked by forest%finalize once per realm at shutdown, after the
   !< time loop exits. Mirror of initialize_forest: every IO file opened
   !< in initialize_forest is closed here; every per-realm resource the
   !< orchestrator did not acquire is released here.
   !<
   !< Default implementation error-stops: every consumer app MUST override
   !< this method (the close-order and final-dump catalog are app-specific).
   !< PRISM's override lives on prism_cpu_object/prism_fnl_object.
   class(realm_object), intent(inout) :: self !< The realm.

   associate(self_unused => self) ! quiet "unused dummy" warnings before the stop
   end associate
   error stop 'realm_object%finalize_forest: not overridden by app extension'
   endsubroutine finalize_forest

   subroutine exchange_inter_realm_halos_forest(self, realm)
   !< Refresh THIS realm's ghost cells that depend on neighbour realms.
   !<
   !< Invoked by forest%exchange_halos once per realm per call. The peer
   !< realms are passed through as the `realm(:)` assumed-shape dummy so the
   !< override can index into `realm(n%peer_realm)%field%q` (host-side chain
   !< walk, in line with R2 of issue #10 — kernels never walk this chain).
   !<
   !< Default implementation is a **no-op**: the base class cannot read
   !< `q` (which is app-private per O3 of issue #10) and therefore cannot
   !< copy peer cells into self's ghosts. Apps that participate in a
   !< multi-realm forest MUST override this TBP; apps that only ever run as
   !< N=1 may leave it. The no-op default keeps single-realm rmf
   !< (`self%maps%inter_realm_neighbors` unallocated) bit-identical to its
   !< pre-Phase-D behaviour.
   !<
   !< Granularity note: this TBP is invoked by the forest *between* whole
   !< timesteps (after `advance_one_step_forest`). Bit-comparability with a
   !< single-realm reference also requires inter-realm refresh between RK
   !< substages — that finer integration point is a separate decision
   !< documented in issue #13 and is NOT addressed by this TBP. The host
   !< app's `advance_one_step_forest` override may also call into the same
   !< exchange machinery between substages if needed.
   class(realm_object), intent(inout) :: self     !< This realm.
   class(realm_object), intent(in)    :: realm(:) !< All realms in the forest (so peers are reachable as realm(n%peer_realm)).

   associate(self_unused => self, realm_unused => realm) ! no-op default; signature locked in for app overrides
   end associate
   endsubroutine exchange_inter_realm_halos_forest

   ! public methods
   subroutine initialize(self, filename, memory_avail, nv, verbose)
   !< Initialize common data: MPI, ADAM, grid, field and data replica pointers.
   class(realm_object), intent(inout), target :: self         !< The equation.
   character(*),           intent(in)            :: filename     !< Input parameters file name.
   real(R8P),              intent(in), value     :: memory_avail !< Memory available for single MPI process.
   integer(I4P),           intent(in), optional  :: nv           !< Number of field variables.
   logical,                intent(in), optional  :: verbose      !< Trigger verbose output.
   logical                                       :: verbose_     !< Trigger verbose output, local variable.
   integer(I8P)                                  :: nodes_number !< Allocated nodes on tree.
   integer(I4P)                                  :: nb           !< Number of allocated blocks.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   call mpih%initialize(verbose=verbose_)
   if (verbose_) call mpih%print_message('realm_object%initialize start')
   call self%io%initialize(filename=trim(filename),verbose=verbose_)
   associate(file_parameters=>self%io%file_parameters)
      ! Bind the five legacy shim singletons (adam, grid, field, tree, maps) into
      ! `self`'s value components BEFORE adam%initialize. The outer `adam`
      ! pointer aliases `self%adam`; the four inner shims (grid/field/tree/maps)
      ! then alias `self%adam%grid` etc. via bind_globals_to_adam (which reads
      ! `adam%...` and is therefore correct once the outer alias is set).
      ! Binding must happen here, not after: tree/field/maps sub-initializations
      ! called inside self%adam%initialize read `grid%...` through the shim, so
      ! the shim must be associated before they run. See issue #10 step 1 and
      ! the C.3 closure note on the `adam` value-component declaration above.
      adam => self%adam
      call bind_globals_to_adam
      call self%adam%initialize(file_parameters=file_parameters, memory_avail=memory_avail, nv=nv, verbose=verbose_)
      ! Step 2 of forest-of-trees migration (issue #10): bind the three legacy
      ! shim singletons (weno, ib, rk) into `self`'s value components BEFORE
      ! their sub-initializations. Same shell-then-populate trick as the adam
      ! shim binding above: the shim aliases the empty `self%weno` shell;
      ! self%weno%initialize then populates that shell; the shim transparently
      ! sees the populated value. Required because consumers (fdv operators,
      ! prism kernels, etc.) read `weno%...` / `ib%...` / `rk%...` through the
      ! shim across the codebase.
      !
      ! The binding is inlined here rather than delegated to a separate
      ! `adam_equation_bind` module: that module would need to `use ::
      ! adam_realm_object`, and this module would need to `use ::
      ! adam_equation_bind` to call it — a direct circular dependency. The
      ! adam binder works because adam_adam_object does not need the binder
      ! itself; the same is not true for realm_object.
      weno => self%weno
      ib   => self%ib
      rk   => self%rk
      call self%amr%initialize(file_parameters=file_parameters)
      call self%ib%initialize(file_parameters=file_parameters)
      call self%slices%initialize(file_parameters=file_parameters)
      ! call self%blanesmoan%initialize(file_parameters=file_parameters)
      ! call self%cfm%initialize(file_parameters=file_parameters)
      ! call self%leapfrog%initialize(file_parameters=file_parameters)
      call self%rk%initialize(file_parameters=file_parameters)
      self%ngc           => grid%ngc
      self%ni            => grid%ni
      self%nj            => grid%nj
      self%nk            => grid%nk
      self%nb            => field%nb
      self%blocks_number => field%blocks_number
      self%nv            => field%nv
      call self%weno%initialize(file_parameters=file_parameters, nb=self%nb, ngc=self%ngc, ni=self%ni, nj=self%nj, nk=self%nk)
      call self%flail%initialize(file_parameters=file_parameters)
      call self%load_fdv_from_file(file_parameters=file_parameters)
   endassociate
   endsubroutine initialize

   subroutine load_fdv_from_file(self, file_parameters, go_on_fail)
   !< Load FDV config from file.
   class(realm_object), intent(inout)        :: self                   !< The equation.
   type(file_ini),         intent(in)           :: file_parameters        !< Simulation parameters ini file handler.
   logical,                intent(in), optional :: go_on_fail             !< Go on if load fails.
   logical                                      :: go_on_fail_            !< Go on if load fails.
   integer(I4P)                                 :: error                  !< Error status.
   character(99)                                :: buff                   !< Character buffer.
   character(len=3), parameter                  :: INI_SECTION_NAME='fdv' !< INI file section name containing FDV config.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='fdv_scheme', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(fdv_scheme)')
   select case(trim(adjustl(buff)))
   case('FD', 'fd', 'Fd', 'fD')
      self%fdv_scheme = 'FD'
      self%compute_block_total_variation => compute_block_total_variation_fd
      self%compute_curl                  => compute_curl_fd
      self%compute_derivative1           => compute_derivative1_fd
      self%compute_derivative2           => compute_derivative2_fd
      self%compute_divergence            => compute_divergence_fd
      self%compute_gradient              => compute_gradient_fd
      self%compute_laplacian             => compute_laplacian_fd
   case('FV', 'fv', 'Fv', 'fV')
      self%fdv_scheme = 'FV'
      ! self%compute_block_total_variation => compute_block_total_variation_fd
      self%compute_curl                  => compute_curl_fv
      self%compute_derivative1           => compute_derivative1_fv
      self%compute_derivative2           => compute_derivative2_fv
      self%compute_divergence            => compute_divergence_fv
      self%compute_gradient              => compute_gradient_fv
      self%compute_laplacian             => compute_laplacian_fv
   case default
      ! implement error message
   endselect

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='fdv_order', val=self%fdv_order,error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(fdv_order)')
   ! valid only for centered schems
   self%fdv_half_stencil = self%fdv_order / 2
   self%fdv_half_stencils(1) = self%fdv_half_stencil + 0
   self%fdv_half_stencils(2) = self%fdv_half_stencil + 0
   self%fdv_half_stencils(3) = self%fdv_half_stencil + 1
   self%fdv_half_stencils(4) = self%fdv_half_stencil + 1
   self%fdv_half_stencils(5) = self%fdv_half_stencil + 2
   self%fdv_half_stencils(6) = self%fdv_half_stencil + 2
   endsubroutine load_fdv_from_file

   ! IO methods
   subroutine close_block_xh5f(xh5f)
   !< Close XH5F file block.
   type(xh5f_file_object), intent(inout) :: xh5f !< XH5F file handler.

   call xh5f%close_block
   endsubroutine close_block_xh5f

   subroutine close_file_xh5f(xh5f)
   !< Close XH5F file.
   type(xh5f_file_object), intent(inout) :: xh5f !< XH5F file handler.

   call xh5f%close_grid
   call xh5f%close_grid(grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION_ASYNC)
   call xh5f%close_file
   endsubroutine close_file_xh5f

   subroutine open_block_xh5f(self, xh5f, b, nijk, t, time)
   !< Open XH5F file block.
   class(realm_object), intent(inout)        :: self    !< The equation.
   type(xh5f_file_object), intent(inout)        :: xh5f    !< XH5F file handler.
   integer(I4P),           intent(in)           :: b       !< Block index.
   integer(I8P),           intent(in)           :: nijk(3) !< Blocks dimensions.
   integer(I4P),           intent(in), optional :: t       !< Time iteration.
   real(R8P),              intent(in), optional :: time    !< Time.
   integer(I4P)                                 :: t_      !< Time iteration, local var.
   real(R8P)                                    :: time_   !< Time, local var.
   real(R8P)                                    :: emin(3) !< Minimum abscissa of current block.
   character(:), allocatable                    :: bn      !< Block name.

   t_    = 0_I4P   ; if (present(t   )) t_    = t
   time_ = 0._R8P  ; if (present(time)) time_ = time
   emin = [field%emin(1,b)-self%ngc*field%dxyz(1,b), &
           field%emin(2,b)-self%ngc*field%dxyz(2,b), &
           field%emin(3,b)-self%ngc*field%dxyz(3,b)]
   bn = 'block_'//trim(strz(b,9))//'-proc'//trim(strz(mpih%myrank,6))
   call xh5f%open_block(block_type = XH5F_PARAMETERS%XH5F_BLOCK_CARTESIAN_UNIFORM, &
                        block_name = bn,                                           &
                        nijk       = nijk,                                         &
                        emin       = emin,                                         &
                        dxyz       = field%dxyz(:,b),                              &
                        time       = time_)
   call xh5f%save_block_field(xdmf_field_name = 'time_iteration',                                &
                              field           = t_,                                              &
                              field_center    = XDMF_PARAMETERS%XDMF_ATTR_CENTER_GRID,           &
                              field_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF, &
                              hdf5_field_name = bn//'-time_iteration')
   endsubroutine open_block_xh5f

   subroutine open_file_xh5f(self, basename, xh5f, directory)
   !< Open XH5F file.
   class(realm_object), intent(inout)        :: self          !< The equation.
   character(*),           intent(in)           :: basename      !< Base name of output files.
   type(xh5f_file_object), intent(out)          :: xh5f          !< XH5F file handler.
   character(*),           intent(in), optional :: directory     !< Directory name of output files.
   character(:), allocatable                    :: directory_    !< Directory name of output files, local var.
   character(:), allocatable                    :: filename_hdf5 !< HDF5 file name.
   character(:), allocatable                    :: filename_xdmf !< XDMF file name.

   directory_        = ''      ; if (present(directory       )) directory_        = trim(adjustl(directory))
   filename_hdf5 = directory_//trim(adjustl(basename))//'-proc'//trim(strz(mpih%myrank,6))//'.h5'
   filename_xdmf = directory_//trim(adjustl(basename))//'.xdmf'
   call xh5f%open_file(filename_hdf5=filename_hdf5, filename_xdmf=filename_xdmf, act=FILE_PARAMETERS%FILE_ACTION_OVERWRITE)
   call xh5f%open_grid(grid_name='adam',                                 grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION_ASYNC)
   call xh5f%open_grid(grid_name='mpi_'//trim(strz(mpih%myrank,6)), grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION)
   endsubroutine open_file_xh5f

   subroutine save_q_xh5f(self, basename, q, q_name, directory, with_ghost, with_cell_morton, t, time)
   !< Save ADAM in XH5F format.
   class(realm_object), intent(inout)        :: self                   !< The equation.
   character(*),           intent(in)           :: basename                !< Base name of output files.
   real(R8P),              intent(in)           :: q(1:,              &
                                                     1-grid%ngc:,&
                                                     1-grid%ngc:,&
                                                     1-grid%ngc:,&
                                                     1:)                   !< Q-vector variables [nv,ni,nj,nk,nb].
   character(*),           intent(in), optional :: q_name(:)               !< Q-vector variables names [nv].
   character(*),           intent(in), optional :: directory               !< Directory name of output files.
   logical,                intent(in), optional :: with_ghost              !< Flag to save ghost cells.
   logical,                intent(in), optional :: with_cell_morton        !< Flag to save Morton code also in cells.
   integer(I4P),           intent(in), optional :: t                       !< Time iteration.
   real(R8P),              intent(in), optional :: time                    !< Time.
   type(string), allocatable                    :: q_name_(:)              !< Q-vector variables names [nv].
   character(:), allocatable                    :: directory_              !< Directory name of output files, local var.
   logical                                      :: with_ghost_             !< Flag to save ghost cells, local var.
   logical                                      :: with_cell_morton_       !< Flag to save Morton code also in cells, local var.
   integer(I4P)                                 :: t_                      !< Time iteration, local var.
   real(R8P)                                    :: time_                   !< Time, local var.
   integer(I4P)                                 :: ngc                     !< Ghost cells saved.
   integer(I4P)                                 :: ijk(2,3)                !< Blocks extents.
   integer(I8P)                                 :: nijk(3)                 !< Blocks dimensions.
   real(R8P)                                    :: emin(3)                 !< Minimum abscissa of current block.
   character(:), allocatable                    :: filename_hdf5           !< HDF5 file name.
   character(:), allocatable                    :: filename_xdmf           !< XDMF file name.
   character(:), allocatable                    :: bn                      !< Block name.
   type(xh5f_file_object)                       :: xh5f                    !< XH5F file handler.
   integer(I4P)                                 :: i, b, v                 !< Counter.

   if (present(q_name)) then
      allocate(q_name_(size(q, dim=1)))
      do v=1, size(q, dim=1)
         q_name_(v) = trim(adjustl(q_name(v)))
      enddo
   else
      do v=1, size(q, dim=1)
         q_name_(v) = 'q-'//trim(strz(v,2))
      enddo
   endif
   directory_        = ''      ; if (present(directory       )) directory_        = trim(adjustl(directory))
   with_ghost_       = .false. ; if (present(with_ghost      )) with_ghost_       = with_ghost
   with_cell_morton_ = .false. ; if (present(with_cell_morton)) with_cell_morton_ = with_cell_morton
   t_                = 0_I4P   ; if (present(t               )) t_                = t
   time_             = 0._R8P  ; if (present(time            )) time_             = time
   if (with_ghost_) then
      ngc = grid%ngc
   else
      ngc = 0_I4P
   endif
   associate(ni=>grid%ni, nj=>grid%nj, nk=>grid%nk)
   ijk(:,1) = [1-ngc,ni+ngc]
   ijk(:,2) = [1-ngc,nj+ngc]
   ijk(:,3) = [1-ngc,nk+ngc]
   nijk = [ijk(2,1)-ijk(1,1)+1, &
           ijk(2,2)-ijk(1,2)+1, &
           ijk(2,3)-ijk(1,3)+1]
   endassociate
   filename_hdf5 = directory_//trim(adjustl(basename))//'-proc'//trim(strz(mpih%myrank,6))//'.h5'
   filename_xdmf = directory_//trim(adjustl(basename))//'.xdmf'
   call xh5f%open_file(filename_hdf5=filename_hdf5, filename_xdmf=filename_xdmf, act=FILE_PARAMETERS%FILE_ACTION_OVERWRITE)
   call xh5f%open_grid(grid_name='adam',                                 grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION_ASYNC)
   call xh5f%open_grid(grid_name='mpi_'//trim(strz(mpih%myrank,6)), grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION)
   do b=1, field%blocks_number
      emin = [field%emin(1,b)-ngc*field%dxyz(1,b), &
              field%emin(2,b)-ngc*field%dxyz(2,b), &
              field%emin(3,b)-ngc*field%dxyz(3,b)]
      bn = 'block_'//trim(strz(b,9))//'-proc'//trim(strz(mpih%myrank,6))
      call xh5f%open_block(block_type = XH5F_PARAMETERS%XH5F_BLOCK_CARTESIAN_UNIFORM, &
                           block_name = bn,                                           &
                           nijk       = nijk,                                         &
                           emin       = emin,                                         &
                           dxyz       = field%dxyz(:,b),                              &
                           time       = time_)
      call xh5f%save_block_field(xdmf_field_name = 'time_iteration',                                &
                                 field           = t_,                                              &
                                 field_center    = XDMF_PARAMETERS%XDMF_ATTR_CENTER_GRID,           &
                                 field_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF, &
                                 hdf5_field_name = bn//'-time_iteration')
      call self%io%save_field(xh5f=xh5f, block_name=bn, ijk=ijk, nijk=nijk, q=q(:,:,:,:,b), q_name=q_name_)
      call xh5f%close_block
   enddo
   call xh5f%close_grid
   call xh5f%close_grid(grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION_ASYNC)
   call xh5f%close_file
   endsubroutine save_q_xh5f

   ! FDV operators numerical methods
   subroutine compute_block_total_variation_fd(self, hs, dxyz, ivar, q, tot_var_field, total_variation)
   !< Return the max of block total variation for a given var.
   class(realm_object), intent(in)     :: self                                                 !< Coils.
   integer(I4P),           intent(in)     :: hs                                                   !< FDV half stencil length.
   real(R8P),              intent(in)     :: dxyz(3)                                              !< Space steps.
   integer(I4P),           intent(in)     :: ivar                                                 !< Index of first component of vec field.
   real(R8P),              intent(in)     :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:)            !< Field variables.
   real(R8P),              intent(inout)  :: tot_var_field(1-self%ngc:,1-self%ngc:,1-self%ngc:)   !< Total variation field on blocks.
   real(R8P),              intent(out)    :: total_variation                                      !< Max total variation on given block.
   real(R8P)                              :: gradient(3,3)                                        !< Gradient.
   real(R8P)                              :: tv                                                   !< Total variation buffer.
   integer(I4P)                           :: i,j,k                                                !< Counter.

   total_variation = -huge(1._R8P)
   do k=1, self%nk
   do j=1, self%nj
   do i=1, self%ni
      call compute_gradient_fd_centered(s=hs,dxyz=dxyz,q=q(ivar,  i-hs:i+hs,j-hs:j+hs,k-hs:k+hs),&
                                        gradient=gradient(:,1))
      call compute_gradient_fd_centered(s=hs,dxyz=dxyz,q=q(ivar+1,i-hs:i+hs,j-hs:j+hs,k-hs:k+hs),&
                                        gradient=gradient(:,2))
      call compute_gradient_fd_centered(s=hs,dxyz=dxyz,q=q(ivar+2,i-hs:i+hs,j-hs:j+hs,k-hs:k+hs),&
                                        gradient=gradient(:,3))
      tv = sqrt(gradient(1,1)*gradient(1,1) + gradient(2,1)*gradient(2,1) + gradient(3,1)*gradient(3,1) +&
                gradient(1,2)*gradient(1,2) + gradient(2,2)*gradient(2,2) + gradient(3,2)*gradient(3,2) +&
                gradient(1,3)*gradient(1,3) + gradient(2,3)*gradient(2,3) + gradient(3,3)*gradient(3,3))
      tot_var_field(i,j,k) = tv
      total_variation = max(total_variation,tv)
   enddo
   enddo
   enddo
   endsubroutine compute_block_total_variation_fd

   subroutine compute_curl_fd(self, hs, ivar, q, curl)
   !< Compute curl of vector fields, div(q(ivar:ivar+2), using finite difference schemes.
   class(realm_object), intent(in)    :: self                                            !< The equation.
   integer(I4P),           intent(in)    :: hs                                              !< FDV half stencil length.
   integer(I4P),           intent(in)    :: ivar                                            !< Start index of variable of q.
   real(R8P),              intent(in)    :: q(   1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),              intent(inout) :: curl(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Curl.
   integer(I4P)                          :: i,j,k,b                                         !< Counter.

   associate(dxyz=>field%dxyz)
   do b=1, self%blocks_number
   do k=1, self%nk
   do j=1, self%nj
   do i=1, self%ni
      call compute_curl_fd_centered(s=hs,dxyz=dxyz(:,b),q=q(ivar:ivar+2,i-hs:i+hs,j-hs:j+hs,k-hs:k+hs,b),curl=curl(ivar:,i,j,k,b))
   enddo
   enddo
   enddo
   enddo
   endassociate
   endsubroutine compute_curl_fd

   subroutine compute_curl_fv(self, hs, ivar, q, curl)
   !< Compute curl of vector fields, div(q(ivar:ivar+2), using finite volume schemes.
   class(realm_object), intent(in)    :: self                                            !< The equation.
   integer(I4P),           intent(in)    :: hs                                              !< FDV half stencil length.
   integer(I4P),           intent(in)    :: ivar                                            !< Start index of variable of q.
   real(R8P),              intent(in)    :: q(   1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),              intent(inout) :: curl(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Curl.
   integer(I4P)                          :: i,j,k,b                                         !< Counter.

   associate(dxyz=>field%dxyz)
   do b=1, self%blocks_number
   do k=1, self%nk
   do j=1, self%nj
   do i=1, self%ni
      call compute_curl_fv_centered(s=hs,dxyz=dxyz(:,b),q=q(ivar:ivar+2,i-hs:i+hs,j-hs:j+hs,k-hs:k+hs,b),curl=curl(ivar:,i,j,k,b))
   enddo
   enddo
   enddo
   enddo
   endassociate
   endsubroutine compute_curl_fv

   subroutine compute_derivative1_fd(self, hs, dir, ivar, q, dq_ds)
   !< Compute derivative1 of scalar fields, dq(ivar)/ds, using finite difference schemes.
   class(realm_object), intent(in)    :: self                                         !< The equation.
   integer(I4P),           intent(in)    :: hs                                           !< FDV half stencil length.
   integer(I4P),           intent(in)    :: dir                                          !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),           intent(in)    :: ivar                                         !< Start index of (vec.) variable of q.
   real(R8P),              intent(in)    :: q(1:, 1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Field variables.
   real(R8P),              intent(inout) :: dq_ds(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Derivative1, dq/ds.
   integer(I4P)                          :: i,j,k,b                                      !< Counter.
   integer(I4P)                          :: is,js,ks                                     !< Stencils.

   associate(dxyz=>field%dxyz)
   select case(dir)
   case(1)
      do b=1, self%blocks_number
      do k=1, self%nk
      do j=1, self%nj
      do i=1, self%ni
         call compute_derivative1_fd_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i-hs:i+hs,j,k,b),dq_ds=dq_ds(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   case(2)
      do b=1, self%blocks_number
      do k=1, self%nk
      do j=1, self%nj
      do i=1, self%ni
         call compute_derivative1_fd_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i,j-hs:j+hs,k,b),dq_ds=dq_ds(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   case(3)
      do b=1, self%blocks_number
      do k=1, self%nk
      do j=1, self%nj
      do i=1, self%ni
         call compute_derivative1_fd_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i,j,k-hs:k+hs,b),dq_ds=dq_ds(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   endselect
   endassociate
   endsubroutine compute_derivative1_fd

   subroutine compute_derivative1_fv(self, hs, dir, ivar, q, dq_ds)
   !< Compute derivative1 of scalar fields, dq(ivar)/ds, using finite volume schemes.
   class(realm_object), intent(in)    :: self                                         !< The equation.
   integer(I4P),           intent(in)    :: hs                                           !< FDV half stencil length.
   integer(I4P),           intent(in)    :: dir                                          !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),           intent(in)    :: ivar                                         !< Start index of (vec.) variable of q.
   real(R8P),              intent(in)    :: q(1:, 1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Field variables.
   real(R8P),              intent(inout) :: dq_ds(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Derivative1, dq/ds.
   integer(I4P)                          :: i,j,k,b                                      !< Counter.

   associate(dxyz=>field%dxyz)
   select case(dir)
   case(1)
      do b=1, self%blocks_number
      do k=1, self%nk
      do j=1, self%nj
      do i=1, self%ni
         call compute_derivative1_fv_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i-hs:i+hs,j,k,b),dq_ds=dq_ds(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   case(2)
      do b=1, self%blocks_number
      do k=1, self%nk
      do j=1, self%nj
      do i=1, self%ni
         call compute_derivative1_fv_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i,j-hs:j+hs,k,b),dq_ds=dq_ds(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   case(3)
      do b=1, self%blocks_number
      do k=1, self%nk
      do j=1, self%nj
      do i=1, self%ni
         call compute_derivative1_fv_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i,j,k-hs:k+hs,b),dq_ds=dq_ds(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   endselect
   endassociate
   endsubroutine compute_derivative1_fv

   subroutine compute_derivative2_fd(self, hs, dir, ivar, q, d2q_ds2)
   !< Compute derivative2 of scalar fields, d2q(ivar)/ds2, using finite difference schemes.
   class(realm_object), intent(in)    :: self                                           !< The equation.
   integer(I4P),           intent(in)    :: hs                                             !< FDV half stencil length.
   integer(I4P),           intent(in)    :: dir                                            !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),           intent(in)    :: ivar                                           !< Start index of vec variable of q.
   real(R8P),              intent(in)    :: q(1:,   1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Field variables.
   real(R8P),              intent(inout) :: d2q_ds2(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Derivative2, d2q/ds2.
   integer(I4P)                          :: i,j,k,b                                        !< Counter.
   integer(I4P)                          :: is,js,ks                                       !< Stencils.

   associate(dxyz=>field%dxyz)
   select case(dir)
   case(1)
      do b=1, self%blocks_number
      do k=1, self%nk
      do j=1, self%nj
      do i=1, self%ni
         call compute_derivative2_fd_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i-hs:i+hs,j,k,b),d2q_ds2=d2q_ds2(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   case(2)
      do b=1, self%blocks_number
      do k=1, self%nk
      do j=1, self%nj
      do i=1, self%ni
         call compute_derivative2_fd_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i,j-hs:j+hs,k,b),d2q_ds2=d2q_ds2(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   case(3)
      do b=1, self%blocks_number
      do k=1, self%nk
      do j=1, self%nj
      do i=1, self%ni
         call compute_derivative2_fd_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i,j,k-hs:k+hs,b),d2q_ds2=d2q_ds2(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   endselect
   endassociate
   endsubroutine compute_derivative2_fd

   subroutine compute_derivative2_fv(self, hs, dir, ivar, q, d2q_ds2)
   !< Compute derivative2 of scalar fields, d2q(ivar)/ds2, using finite volume schemes.
   class(realm_object), intent(in)    :: self                                            !< The equation.
   integer(I4P),           intent(in)    :: hs                                              !< FDV half stencil length.
   integer(I4P),           intent(in)    :: dir                                             !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),           intent(in)    :: ivar                                            !< Start index of variable of q.
   real(R8P),              intent(in)    :: q(1:,   1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),              intent(inout) :: d2q_ds2(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Derivative2, d2q/ds2.
   integer(I4P)                          :: i,j,k,b                                         !< Counter.
   integer(I4P)                          :: is,js,ks                                        !< Stencils.

   associate(dxyz=>field%dxyz)
   select case(dir)
   case(1)
      do b=1, self%blocks_number
      do k=1, self%nk
      do j=1, self%nj
      do i=1, self%ni
         ! call compute_derivative2_fv_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i-hs:i+hs,j,k,b),d2q_ds2=d2q_ds2(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   case(2)
      do b=1, self%blocks_number
      do k=1, self%nk
      do j=1, self%nj
      do i=1, self%ni
         ! call compute_derivative2_fv_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i,j-hs:j+hs,k,b),d2q_ds2=d2q_ds2(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   case(3)
      do b=1, self%blocks_number
      do k=1, self%nk
      do j=1, self%nj
      do i=1, self%ni
         ! call compute_derivative2_fv_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i,j,k-hs:k+hs,b),d2q_ds2=d2q_ds2(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   endselect
   endassociate
   endsubroutine compute_derivative2_fv

   subroutine compute_derivative4_fd(self, hs, dir, ivar, q, d4q_ds4)
   !< Compute derivative2 of scalar fields, d4q(ivar)/ds4, using finite difference schemes.
   class(realm_object), intent(in)    :: self                                            !< The equation.
   integer(I4P),           intent(in)    :: hs                                              !< FDV half stencil length.
   integer(I4P),           intent(in)    :: dir                                             !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),           intent(in)    :: ivar                                            !< Start index of variable of q.
   real(R8P),              intent(in)    :: q(1:,   1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),              intent(inout) :: d4q_ds4(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Derivative4, d4q/ds4.
   integer(I4P)                          :: i,j,k,b                                         !< Counter.
   integer(I4P)                          :: is,js,ks                                        !< Stencils.

   associate(dxyz=>field%dxyz)
   select case(dir)
   case(1)
      do b=1, self%blocks_number
      do k=1, self%nk
      do j=1, self%nj
      do i=1, self%ni
         call compute_derivative4_fd_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i-hs:i+hs,j,k,b),d4q_ds4=d4q_ds4(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   case(2)
      do b=1, self%blocks_number
      do k=1, self%nk
      do j=1, self%nj
      do i=1, self%ni
         call compute_derivative4_fd_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i,j-hs:j+hs,k,b),d4q_ds4=d4q_ds4(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   case(3)
      do b=1, self%blocks_number
      do k=1, self%nk
      do j=1, self%nj
      do i=1, self%ni
         call compute_derivative4_fd_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i,j,k-hs:k+hs,b),d4q_ds4=d4q_ds4(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   endselect
   endassociate
   endsubroutine compute_derivative4_fd

   subroutine compute_divergence_fd(self, hs, ivar, q, divergence)
   !< Compute divergence of vector fields, div(q(ivar:ivar+2), using finite difference schemes.
   class(realm_object), intent(in)    :: self                                               !< The equation.
   integer(I4P),           intent(in)    :: hs                                                 !< FDV half stencil length.
   integer(I4P),           intent(in)    :: ivar                                               !< Start index of field of q.
   real(R8P),              intent(in)    :: q(1:,      1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),              intent(inout) :: divergence(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Divergence.
   integer(I4P)                          :: i,j,k,b                                            !< Counter.

   associate(dxyz=>field%dxyz)
   do b=1, self%blocks_number
   do k=1, self%nk
   do j=1, self%nj
   do i=1, self%ni
      call compute_divergence_fd_centered(s=hs,dxyz=dxyz(:,b),q=q(ivar:ivar+2,i-hs:i+hs,j-hs:j+hs,k-hs:k+hs,b),&
                                          divergence=divergence(i,j,k,b))
   enddo
   enddo
   enddo
   enddo
   endassociate
   endsubroutine compute_divergence_fd

   subroutine compute_divergence_fv(self, hs, ivar, q, divergence)
   !< Compute divergence of vector fields, div(q(ivar:ivar+2), using finite volume schemes.
   class(realm_object), intent(in)    :: self                                               !< The equation.
   integer(I4P),           intent(in)    :: hs                                                 !< FDV half stencil length.
   integer(I4P),           intent(in)    :: ivar                                               !< Start index of field of q.
   real(R8P),              intent(in)    :: q(1:,      1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),              intent(inout) :: divergence(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Divergence.
   integer(I4P)                          :: i,j,k,b                                            !< Counter.

   associate(dxyz=>field%dxyz)
   do b=1, self%blocks_number
   do k=1, self%nk
   do j=1, self%nj
   do i=1, self%ni
      call compute_divergence_fv_centered(s=hs,dxyz=dxyz(:,b),q=q(ivar:ivar+2,i-hs:i+hs,j-hs:j+hs,k-hs:k+hs,b), &
                                          divergence=divergence(i,j,k,b))
   enddo
   enddo
   enddo
   enddo
   endassociate
   endsubroutine compute_divergence_fv

   subroutine compute_gradient_fd(self, hs, ivar, q, gradient)
   !< Compute gradient of scalar variable q(ivar), finite difference schemes.
   class(realm_object), intent(in)    :: self                                               !< The equation.
   integer(I4P),           intent(in)    :: hs                                                 !< FDV half stencil length.
   integer(I4P),           intent(in)    :: ivar                                               !< Index of scalar variable of q.
   real(R8P),              intent(in)    :: q(       1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Field variables.
   real(R8P),              intent(inout) :: gradient(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Gradient.
   integer(I4P)                          :: i, j, k, b                                         !< Counter.

   associate(dxyz=>field%dxyz)
   do b=1, self%blocks_number
   do k=1, self%nk
   do j=1, self%nj
   do i=1, self%ni
      call compute_gradient_fd_centered(s=hs,dxyz=dxyz(:,b),q=q(ivar,i-hs:i+hs,j-hs:j+hs,k-hs:k+hs,b),&
                                        gradient=gradient(1:3,i,j,k,b))
   enddo
   enddo
   enddo
   enddo
   endassociate
   endsubroutine compute_gradient_fd

   subroutine compute_gradient_fv(self, hs, ivar, q, gradient)
   !< Compute gradient of scalar variable q(ivar), finite volume schemes.
   class(realm_object), intent(in)    :: self                                               !< The equation.
   integer(I4P),           intent(in)    :: hs                                                 !< FDV half stencil length.
   integer(I4P),           intent(in)    :: ivar                                               !< Index of scalar variable of q.
   real(R8P),              intent(in)    :: q(       1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Field variables.
   real(R8P),              intent(inout) :: gradient(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Gradient.
   integer(I4P)                          :: i, j, k, b                                         !< Counter.

   associate(dxyz=>field%dxyz)
   do b=1, self%blocks_number
   do k=1, self%nk
   do j=1, self%nj
   do i=1, self%ni
      call compute_gradient_fv_centered(s=hs,dxyz=dxyz(:,b),q=q(ivar,i-hs:i+hs,j-hs:j+hs,k-hs:k+hs,b),&
                                        gradient=gradient(1:3,i,j,k,b))
   enddo
   enddo
   enddo
   enddo
   endassociate
   endsubroutine compute_gradient_fv

   subroutine compute_laplacian_fd(self, hs, ivar, q, laplacian)
   !< Compute laplacian of scalar variable q(ivar), finite difference schemes.
   class(realm_object), intent(in)    :: self                                              !< The equation.
   integer(I4P),           intent(in)    :: hs                                                !< FDV half stencil length.
   integer(I4P),           intent(in)    :: ivar                                              !< Index of scalar variable of q.
   real(R8P),              intent(in)    :: q(     1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),              intent(inout) :: laplacian(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Gradient.
   integer(I4P)                          :: i, j, k, b                                        !< Counter.

   associate(dxyz=>field%dxyz)
   do b=1, self%blocks_number
   do k=1, self%nk
   do j=1, self%nj
   do i=1, self%ni
      call compute_laplacian_fd_centered(s=hs,dxyz=dxyz(:,b),q=q(ivar,i-hs:i+hs,j-hs:j+hs,k-hs:k+hs,b),laplacian=laplacian(i,j,k,b))
   enddo
   enddo
   enddo
   enddo
   endassociate
   endsubroutine compute_laplacian_fd

   subroutine compute_laplacian_fv(self, hs, ivar, q, laplacian)
   !< Compute laplacian of scalar variable q(ivar), finite volume schemes.
   class(realm_object), intent(in)    :: self                                              !< The equation.
   integer(I4P),           intent(in)    :: hs                                                !< FDV half stencil length.
   integer(I4P),           intent(in)    :: ivar                                              !< Index of scalar variable of q.
   real(R8P),              intent(in)    :: q(     1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),              intent(inout) :: laplacian(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Gradient.
   integer(I4P)                          :: i, j, k, b                                        !< Counter.

   associate(dxyz=>field%dxyz)
   do b=1, self%blocks_number
   do k=1, self%nk
   do j=1, self%nj
   do i=1, self%ni
     !call compute_laplacian_fv_centered(s=hs,dxyz=dxyz(:,b),q=q(ivar,i-hs:i+hs,j-hs:j+hs,k-hs:k+hs,b),laplacian=laplacian(i,j,k,b))
   enddo
   enddo
   enddo
   enddo
   endassociate
   endsubroutine compute_laplacian_fv
endmodule adam_realm_object
