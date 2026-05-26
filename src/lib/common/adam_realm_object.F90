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
!< architectural rules govern these TBPs:
!<
!<   * O1 — signature uses only ADAM-lib-visible types (intrinsic kinds plus
!<     `class(realm_object)`, `grid_object`, etc.). Apps override the
!<     implementation; signatures are sacred.
!<   * O2 — app-specific dispatch (which integrator, which physics) lives in
!<     `self%`'s components, not in orchestrator-supplied arguments.
!<   * O3 — `q` is never public on `realm_object`. Halo exchange between realms
!<     goes through realm-side TBPs that operate on q internally;
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
use :: adam_flux_register_object
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
use :: adam_globals, only : mpih
! third party modules
use :: finer
use :: motion
use :: penf
use :: stringifor

implicit none
private
public :: realm_object

type :: realm_object
   !< Realm system class definition, common data to all backends and applications.
   !<
   !< @NOTE: the orchestrator contract methods carrying the `_forest` suffix are the surface
   !< forest_object may call on a realm. Default implementations here error-stop with
   !< a "not overridden" message; each consumer app must override.
   ! ADAM library objects
   type(io_object)         :: io         !< IO handler.
   type(amr_object)        :: amr        !< AMR marker handler.
   type(slices_object)     :: slices     !< Slices handler.
   type(blanesmoan_object) :: blanesmoan !< Blanes-Moan integrator.
   type(cfm_object)        :: cfm        !< Commutator-Free Magnus integrator.
   type(leapfrog_object)   :: leapfrog   !< Leapfrog integrator.
   type(flail_object)      :: flail      !< Linear algebra methods handler.
   type(weno_object)       :: weno       !< WENO reconstructor.
   type(ib_object)         :: ib         !< Immersed boundary.
   type(rk_object)         :: rk         !< Runge-Kutta integrator.
   type(adam_object)       :: adam       !< ADAM, grid + tree + field + maps orchestrator.
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
      ! forest orchestrator contract methods
      procedure, pass(self) :: initialize_forest                 !< Invoked by forest%initialize per realm at startup.
      procedure, pass(self) :: compute_local_dt_forest           !< Invoked by forest%compute_global_dt during the min reduction.
      procedure, pass(self) :: advance_one_step_forest           !< Invoked by forest%evolve_one_step per realm per timestep.
      procedure, pass(self) :: nrk_forest                        !< Number of substages this realm's integrator uses.
      procedure, pass(self) :: prepare_step_forest               !< Per-step prologue (multi-realm path).
      procedure, pass(self) :: assemble_substage_forest          !< One RK substage assembly (multi-realm path).
      procedure, pass(self) :: residuals_substage_forest         !< One RK substage residual computation.
      procedure, pass(self) :: assign_substage_forest            !< One RK substage stage assignment.
      procedure, pass(self) :: evaluate_substage_forest          !< [Deprecated] residuals + assignment in one call.
      procedure, pass(self) :: finalize_step_forest              !< Per-step epilogue (multi-realm path).
      procedure, pass(self) :: post_step_forest                  !< Invoked by forest%post_step per realm per timestep.
      procedure, pass(self) :: is_done_forest                    !< Invoked by forest%is_done during the termination reduction.
      procedure, pass(self) :: finalize_forest                   !< Invoked by forest%finalize per realm at shutdown.
      procedure, pass(self) :: finalize_mpi_forest               !< Process-global MPI finalize; forest calls it ONCE after all.
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
   class(realm_object), intent(in)     :: self                       !< Coils.
   integer(I4P),        intent(in)     :: hs                         !< FDV half stencil length.
   real(R8P),           intent(in)     :: dxyz(3)                    !< Space steps.
   integer(I4P),        intent(in)     :: ivar                       !< Index of first component of vec field.
   real(R8P),           intent(in)     :: q(1:,&
                                            1-self%ngc:,&
                                            1-self%ngc:,&
                                            1-self%ngc:)             !< Field variables.
   real(R8P),           intent(inout)  :: tot_var_field(1-self%ngc:,&
                                                        1-self%ngc:,&
                                                        1-self%ngc:) !< Total variation field on blocks.
   real(R8P),           intent(out)    :: total_variation            !< Max total variation on given block.
   endsubroutine compute_block_total_variation_interface

   subroutine compute_curl_interface(self, hs, ivar, q, curl)
   !< Compute curl of vector field, curl(q(ivar:ivar+2)).
   import :: realm_object, I4P, R8P
   class(realm_object), intent(in)    :: self                                            !< The equation.
   integer(I4P),        intent(in)    :: hs                                              !< FDV half stencil length.
   integer(I4P),        intent(in)    :: ivar                                            !< Start index of variable of q.
   real(R8P),           intent(in)    :: q(   1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),           intent(inout) :: curl(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Curl.
   endsubroutine compute_curl_interface

   subroutine compute_derivative1_interface(self, hs, dir, ivar, q, dq_ds)
   !< Compute first derivative of scalar field, dq(ivar)/ds.
   import :: realm_object, I4P, R8P
   class(realm_object), intent(in)    :: self                                             !< The equation.
   integer(I4P),        intent(in)    :: hs                                               !< FDV half stencil length.
   integer(I4P),        intent(in)    :: dir                                              !< Direction: 1=X, 2=Y, 3=Z.
   integer(I4P),        intent(in)    :: ivar                                             !< Index of variable of q.
   real(R8P),           intent(in)    :: q(    1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),           intent(inout) :: dq_ds(   1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< First derivative dq/ds.
   endsubroutine compute_derivative1_interface

   subroutine compute_derivative2_interface(self, hs, dir, ivar, q, d2q_ds2)
   !< Compute second derivative of scalar field, d2q(ivar)/ds2.
   import :: realm_object, I4P, R8P
   class(realm_object), intent(in)    :: self                                               !< The equation.
   integer(I4P),        intent(in)    :: hs                                                 !< FDV half stencil length.
   integer(I4P),        intent(in)    :: dir                                                !< Direction: 1=X, 2=Y, 3=Z.
   integer(I4P),        intent(in)    :: ivar                                               !< Index of variable of q.
   real(R8P),           intent(in)    :: q(      1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),           intent(inout) :: d2q_ds2(   1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Second derivative d2q/ds2.
   endsubroutine compute_derivative2_interface

   subroutine compute_derivative4_interface(self, hs, dir, ivar, q, d4q_ds4)
   !< Compute fourth derivative of scalar field, d4q(ivar)/ds4.
   import :: realm_object, I4P, R8P
   class(realm_object), intent(in)    :: self                                               !< The equation.
   integer(I4P),        intent(in)    :: hs                                                 !< FDV half stencil length.
   integer(I4P),        intent(in)    :: dir                                                !< Direction: 1=X, 2=Y, 3=Z.
   integer(I4P),        intent(in)    :: ivar                                               !< Index of variable of q.
   real(R8P),           intent(in)    :: q(      1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),           intent(inout) :: d4q_ds4(   1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Fourth derivative d4q/ds4.
   endsubroutine compute_derivative4_interface

   subroutine compute_divergence_interface(self, hs, ivar, q, divergence)
   !< Compute divergence of vector field, div(q(ivar:ivar+2)).
   import :: realm_object, I4P, R8P
   class(realm_object), intent(in)    :: self                                                  !< The equation.
   integer(I4P),        intent(in)    :: hs                                                    !< FDV half stencil length.
   integer(I4P),        intent(in)    :: ivar                                                  !< Start index of field of q.
   real(R8P),           intent(in)    :: q(         1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),           intent(inout) :: divergence(   1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Divergence.
   endsubroutine compute_divergence_interface

   subroutine compute_gradient_interface(self, hs, ivar, q, gradient)
   !< Compute gradient of scalar variable q(ivar).
   import :: realm_object, I4P, R8P
   class(realm_object), intent(in)    :: self                                                  !< The equation.
   integer(I4P),        intent(in)    :: hs                                                    !< FDV half stencil length.
   integer(I4P),        intent(in)    :: ivar                                                  !< Index of scalar variable of q.
   real(R8P),           intent(in)    :: q(        1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)  !< Field variables.
   real(R8P),           intent(inout) :: gradient( 1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)  !< Gradient.
   endsubroutine compute_gradient_interface

   subroutine compute_laplacian_interface(self, hs, ivar, q, laplacian)
   !< Compute laplacian of scalar variable q(ivar).
   import :: realm_object, I4P, R8P
   class(realm_object), intent(in)    :: self                                                  !< The equation.
   integer(I4P),        intent(in)    :: hs                                                    !< FDV half stencil length.
   integer(I4P),        intent(in)    :: ivar                                                  !< Index of scalar variable of q.
   real(R8P),           intent(in)    :: q(        1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)  !< Field variables.
   real(R8P),           intent(inout) :: laplacian(   1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)  !< Laplacian.
   endsubroutine compute_laplacian_interface
endinterface

contains
   ! public methods
   subroutine initialize(self, filename, memory_avail, nv, verbose, L0)
   !< Initialize common data: MPI, ADAM, grid, field and data replica pointers.
   class(realm_object), intent(inout), target :: self         !< The equation.
   character(*),        intent(in)            :: filename     !< Input parameters file name.
   real(R8P),           intent(in), value     :: memory_avail !< Memory available for single MPI process.
   integer(I4P),        intent(in), optional  :: nv           !< Number of field variables.
   logical,             intent(in), optional  :: verbose      !< Trigger verbose output.
   real(R8P),           intent(in), optional  :: L0           !< Adimensionalization parameter.
   logical                                    :: verbose_     !< Trigger verbose output, local variable.
   integer(I8P)                               :: nodes_number !< Allocated nodes on tree.
   integer(I4P)                               :: nb           !< Number of allocated blocks.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   call mpih%initialize(verbose=verbose_)
   if (verbose_) call mpih%print_message('realm_object%initialize start')
   call self%io%initialize(filename=trim(filename),verbose=verbose_)
   associate(file_parameters=>self%io%file_parameters)
      ! The legacy adam/grid/field/tree/maps shim singletons have been retired;
      ! adam%initialize and its sub-initializations now read their realm-local
      ! objects via self%... or threaded dummy arguments, so no pre-binding is
      ! needed here (see issue #15).
      if (present(L0)) then
         call self%adam%initialize(file_parameters=file_parameters, memory_avail=memory_avail, nv=nv, verbose=verbose_, L0=L0)
      else
         call self%adam%initialize(file_parameters=file_parameters, memory_avail=memory_avail, nv=nv, verbose=verbose_)
      endif
      call self%amr%initialize(file_parameters=file_parameters)
      call self%ib%initialize(field=self%adam%field, grid=self%adam%grid, file_parameters=file_parameters)
      call self%slices%initialize(file_parameters=file_parameters)
      ! call self%blanesmoan%initialize(file_parameters=file_parameters)
      ! call self%cfm%initialize(file_parameters=file_parameters)
      ! call self%leapfrog%initialize(file_parameters=file_parameters)
      call self%rk%initialize(field=self%adam%field, grid=self%adam%grid, file_parameters=file_parameters)
      self%ngc           => self%adam%grid%ngc
      self%ni            => self%adam%grid%ni
      self%nj            => self%adam%grid%nj
      self%nk            => self%adam%grid%nk
      self%nb            => self%adam%field%nb
      self%blocks_number => self%adam%field%blocks_number
      self%nv            => self%adam%field%nv
      call self%weno%initialize(file_parameters=file_parameters, nb=self%nb, ngc=self%ngc, ni=self%ni, nj=self%nj, nk=self%nk)
      call self%flail%initialize(file_parameters=file_parameters)
      call self%load_fdv_from_file(file_parameters=file_parameters)
   endassociate
   endsubroutine initialize

   subroutine load_fdv_from_file(self, file_parameters, go_on_fail)
   !< Load FDV config from file.
   class(realm_object), intent(inout)        :: self                   !< The equation.
   type(file_ini),      intent(in)           :: file_parameters        !< Simulation parameters ini file handler.
   logical,             intent(in), optional :: go_on_fail             !< Go on if load fails.
   logical                                   :: go_on_fail_            !< Go on if load fails.
   integer(I4P)                              :: error                  !< Error status.
   character(99)                             :: buff                   !< Character buffer.
   character(len=3), parameter               :: INI_SECTION_NAME='fdv' !< INI file section name containing FDV config.

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
      call mpih%error_stop(msg=': unknown FDV scheme "'//trim(adjustl(buff))//'"')
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

   ! forest orchestrator contract methods
   subroutine initialize_forest(self, filename, realm_index, realms_number, memory_avail, nv, verbose, realm)
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
   class(realm_object), intent(inout)                   :: self          !< The realm.
   character(*),        intent(in)                      :: filename      !< Input parameters file name.
   integer(I4P),        intent(in),    optional         :: realm_index   !< Index of this realm in the forest (Phase D).
   integer(I4P),        intent(in),    optional         :: realms_number !< Realm count in the forest.
   real(R8P),           intent(in),    optional         :: memory_avail  !< Per-process memory budget override.
   integer(I4P),        intent(in),    optional         :: nv            !< Number of field variables override.
   logical,             intent(in),    optional         :: verbose       !< Trigger verbose output.
   class(realm_object), intent(inout), optional, target :: realm(:)      !< Sibling realms (forest array) for inter-realm.

   call mpih%error_stop(msg='realm_object%initialize_forest: not overridden by app extension')
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

   call mpih%error_stop(msg='realm_object%compute_local_dt_forest: not overridden by app extension')
   endsubroutine compute_local_dt_forest

   subroutine advance_one_step_forest(self, dt)
   !< Advance this realm by one full timestep of size `dt`.
   !<
   !< Invoked by forest%evolve_one_step once per realm per timestep when the
   !< forest contains a single realm (N=1 fast path). For N>1 the forest
   !< drives the substage loop itself via `prepare_step_forest`,
   !< `compute_substage_forest`, `finalize_step_forest`, and this TBP is
   !< NOT invoked — see [[forest_object]]%`evolve_one_step` for the
   !< multi-realm path.
   !<
   !< The orchestrator owns global dt selection (compute_global_dt) and the
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

   call mpih%error_stop(msg='realm_object%advance_one_step_forest: not overridden by app extension')
   endsubroutine advance_one_step_forest

   function nrk_forest(self) result(nrk)
   !< Return the number of RK substages this realm's integrator uses.
   !<
   !< Invoked by forest%evolve_one_step on the multi-realm path to size the
   !< substage loop. The forest assumes every realm reports the same `nrk`
   !< (Phase D v1 — all realms run the same integrator scheme); a mismatch
   !< is flagged by the forest with `mpih%error_stop`.
   !<
   !< Default implementation error-stops: every realm that participates in a
   !< multi-realm forest MUST override this. Apps that only ever run as N=1
   !< may leave it (the N=1 fast path doesn't call this).
   class(realm_object), intent(in) :: self !< The realm.
   integer(I4P)                    :: nrk  !< Number of substages.

   call mpih%error_stop(msg='realm_object%nrk_forest: not overridden by app extension')
   endfunction nrk_forest

   subroutine prepare_step_forest(self, dt)
   !< Per-step prologue on the multi-realm path: pre-substage setup.
   !<
   !< Invoked once per realm per timestep, BEFORE the forest's substage loop
   !< starts. Body owns whatever the integrator does between substages once
   !< per step: external-field prelude, RK stage-array initialization, time
   !< bookkeeping (it increment, dt capping for time_max), progress prints.
   !<
   !< Default implementation error-stops: every multi-realm participant MUST
   !< override this. N=1 apps can leave it.
   class(realm_object), intent(inout) :: self !< The realm.
   real(R8P),           intent(in)    :: dt   !< Timestep size from the forest.

   call mpih%error_stop(msg='realm_object%prepare_step_forest: not overridden by app extension')
   endsubroutine prepare_step_forest

   subroutine assemble_substage_forest(self, s, nrk, dt, realm)
   !< Assemble substage q-buffer on the multi-realm path.
   !<
   !< Invoked once per (realm, substage) by the forest, with `s` iterating
   !< from 1 to `nrk`. Body assembles `rk%q_rk(:,...,s)` from the previously
   !< computed substages and from `self%q` — i.e. it runs whatever
   !< `rk%compute_stage(s)` does in the legacy integrator body — and stops
   !< there. Crucially it does NOT invoke `compute_residuals` or anything
   !< that touches ghost cells from peer realms; the peer realms may not
   !< yet have assembled their own substage-s buffer when this fires.
   !<
   !< The forest invokes `assemble_substage_forest` on EVERY realm before
   !< proceeding to `evaluate_substage_forest`, so by the time residuals
   !< evaluate, every peer realm has its substage-s buffer populated and
   !< the inter-realm ghost exchange returns coherent values.
   !<
   !< Default implementation error-stops: every multi-realm participant MUST
   !< override this.
   class(realm_object), intent(inout)                   :: self     !< The realm.
   integer(I4P),        intent(in)                      :: s        !< Substage index (1..nrk).
   integer(I4P),        intent(in)                      :: nrk      !< Total number of substages.
   real(R8P),           intent(in)                      :: dt       !< Timestep size from the forest.
   class(realm_object), intent(inout), optional, target :: realm(:) !< Sibling realms (contract parity).

   call mpih%error_stop(msg='realm_object%assemble_substage_forest: not overridden by app extension')
   endsubroutine assemble_substage_forest

   subroutine residuals_substage_forest(self, s, nrk, dt, realm, flux_register)
   !< Compute residuals for RK substage `s` on the multi-realm path.
   !<
   !< Invoked once per (realm, substage) AFTER every realm has run
   !< `assemble_substage_forest(s)`. Body calls `compute_residuals(q=
   !< rk%q_rk(:,...,s), dq=self%dq, s=s)` ONLY — it must NOT call
   !< `rk%assign_stage`, because doing so would overwrite
   !< `rk%q_rk(:, interior, :, b, s)` with `dq` IN-PLACE, corrupting the
   !< substage-state read by peers' subsequent residual evaluations at
   !< the seam.
   !<
   !< The forest invokes `residuals_substage_forest` on EVERY realm
   !< before any realm runs `assign_substage_forest(s)` — that is the
   !< three-phase substage loop that keeps `q_rk(:, interior, :, b, s)`
   !< stable as the assembled substage state for the duration of all
   !< realms' residual evaluations.
   !<
   !< Default implementation error-stops: every multi-realm participant
   !< MUST override this.
   !<
   !< Optional `realm(:)`: forest realm array
   !< for inter-realm halo refresh during update_ghost calls inside the
   !< residual computation. When present, the FNL backend forwards it
   !< through to `update_ghost(realm=realm)` so the seam exchange uses
   !< the dummy-argument path instead of the polymorphic `forest_realm`
   !< module pointer (works around an nvfortran-OpenACC runtime bug
   !< triggered by polymorphic-class-pointer dereferences followed by
   !< FUNDAL cross-realm device memcpy).
   class(realm_object),         intent(inout)                   :: self          !< The realm.
   integer(I4P),                intent(in)                      :: s             !< Substage index (1..nrk).
   integer(I4P),                intent(in)                      :: nrk           !< Total number of substages.
   real(R8P),                   intent(in)                      :: dt            !< Timestep size from the forest.
   class(realm_object),         intent(inout), optional, target :: realm(:)      !< Sibling realms for inter-realm halo refresh.
   class(flux_register_object), intent(inout), optional         :: flux_register !< Forest's flux register for FV reflux.

   call mpih%error_stop(msg='realm_object%residuals_substage_forest: not overridden by app extension')
   endsubroutine residuals_substage_forest

   subroutine assign_substage_forest(self, s, nrk, dt, realm)
   !< Assign RK substage `s` state from `self%dq` on the multi-realm path.
   !<
   !< Invoked once per (realm, substage) AFTER every realm has run
   !< `residuals_substage_forest(s)`. Body calls `rk%assign_stage(s=s,
   !< q=self%dq)`, which overwrites `rk%q_rk(:, interior, :, b, s)` with
   !< the residual buffer (the SSP-RK trick — q_rk(s) holds the dq of
   !< stage s for the next stage's linear combination).
   !<
   !< By the time this fires, every realm's residual evaluation has
   !< completed, so the in-place overwrite of `q_rk(s)` does not race
   !< with peer reads.
   !<
   !< Default implementation error-stops: every multi-realm participant
   !< MUST override this.
   class(realm_object), intent(inout)                   :: self     !< The realm.
   integer(I4P),        intent(in)                      :: s        !< Substage index (1..nrk).
   integer(I4P),        intent(in)                      :: nrk      !< Total number of substages.
   real(R8P),           intent(in)                      :: dt       !< Timestep size from the forest.
   class(realm_object), intent(inout), optional, target :: realm(:) !< Sibling realms (contract parity).

   call mpih%error_stop(msg='realm_object%assign_substage_forest: not overridden by app extension')
   endsubroutine assign_substage_forest

   subroutine evaluate_substage_forest(self, s, nrk, dt, realm)
   !< [Deprecated] Combined residuals + stage assignment.
   !<
   !< Was the single-phase substage entry point before the
   !< residuals/assign split required by the BC_SEAM-coherence fix
   !< (issue #13, 2026-05-24). Retained as a backward-compat wrapper:
   !< default implementation calls residuals then assign in sequence.
   !< The forest no longer invokes this — it drives the three-phase
   !< loop directly via the split TBPs. App overrides may still
   !< implement this for use outside the forest orchestrator.
   class(realm_object), intent(inout)                   :: self     !< The realm.
   integer(I4P),        intent(in)                      :: s        !< Substage index (1..nrk).
   integer(I4P),        intent(in)                      :: nrk      !< Total number of substages.
   real(R8P),           intent(in)                      :: dt       !< Timestep size from the forest.
   class(realm_object), intent(inout), optional, target :: realm(:) !< Sibling realms (forwarded to residuals/assign).

   if (present(realm)) then
      call self%residuals_substage_forest(s=s, nrk=nrk, dt=dt, realm=realm)
      call self%assign_substage_forest(s=s, nrk=nrk, dt=dt, realm=realm)
   else
      call self%residuals_substage_forest(s=s, nrk=nrk, dt=dt)
      call self%assign_substage_forest(s=s, nrk=nrk, dt=dt)
   endif
   endsubroutine evaluate_substage_forest

   subroutine finalize_step_forest(self, dt)
   !< Per-step epilogue on the multi-realm path: post-substage finalization.
   !<
   !< Invoked once per realm per timestep, AFTER the forest's substage loop
   !< completes. Body owns whatever the integrator does once per step at the
   !< end: q assembly from substage state (rk%update_q), BC update,
   !< divergence cleaning, coil current refresh, residual save, time advance.
   !<
   !< Default implementation error-stops: every multi-realm participant MUST
   !< override this.
   class(realm_object), intent(inout) :: self !< The realm.
   real(R8P),           intent(in)    :: dt   !< Timestep size from the forest.

   call mpih%error_stop(msg='realm_object%finalize_step_forest: not overridden by app extension')
   endsubroutine finalize_step_forest

   subroutine post_step_forest(self, dt, t, it, do_save_state, do_save_residuals, do_save_restart, do_amr, realm)
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
   !<
   !< Optional `realm(:)` (issue #13, 2026-05-25): forest realm array
   !< threaded for inter-realm halo refresh inside `save_simulation_data`
   !< (which calls update_ghost). FNL backend forwards it through;
   !< CPU backend ignores it.
   class(realm_object), intent(inout)                   :: self              !< The realm.
   real(R8P),           intent(in)                      :: dt                !< Timestep size just advanced.
   real(R8P),           intent(in)                      :: t                 !< Simulation time after the advance.
   integer(I4P),        intent(in)                      :: it                !< Iteration index after the advance.
   logical,             intent(in),    optional         :: do_save_state     !< Save state output this step.
   logical,             intent(in),    optional         :: do_save_residuals !< Save residuals output this step.
   logical,             intent(in),    optional         :: do_save_restart   !< Save restart dump this step.
   logical,             intent(in),    optional         :: do_amr            !< Run AMR update this step.
   class(realm_object), intent(inout), optional, target :: realm(:)          !< Sibling realms for inter-realm halo refresh.

   call mpih%error_stop(msg='realm_object%post_step_forest: not overridden by app extension')
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

   call mpih%error_stop(msg='realm_object%is_done_forest: not overridden by app extension')
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

   call mpih%error_stop(msg='realm_object%finalize_forest: not overridden by app extension')
   endsubroutine finalize_forest

   subroutine finalize_mpi_forest(self)
   !< Finalize the process-global MPI handler.
   !<
   !< MPI_FINALIZE is process-global and must run EXACTLY ONCE, after every
   !< realm has completed all MPI-using teardown. The forest therefore calls
   !< this on a SINGLE realm after the per-realm `finalize_forest` loop — not
   !< inside `finalize_forest` (which runs per realm; doing the finalize there
   !< tore MPI down while later realms still had IO/ghost work pending, e.g.
   !< issue #13 rmf-2realm: realm 2's update_ghost hit MPI_Type_f2c after
   !< realm 1 had already finalized). The default targets the CPU `mpih`
   !< singleton; the FNL backend overrides to finalize `mpih_fnl`.
   class(realm_object), intent(inout) :: self !< The realm (carries no MPI state; dispatch picks the backend handler).

   call mpih%finalize
   endsubroutine finalize_mpi_forest

   subroutine exchange_inter_realm_halos_forest(self, realm, s_active)
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
   !< timesteps (after `advance_one_step_forest`) AND between RK substages
   !< from inside the realm's `update_ghost` path. The optional `s_active`
   !< dummy carries the current substage index when called from substage
   !< depth (1..nrk), or 0 when called outside a substage (e.g. initial
   !< seam refresh during `initialize_forest`). The CPU/FNL overrides use
   !< `s_active` to select the peer's substage-`s` buffer (`q_rk(s)` or
   !< `q_rk_gpu(s)`) instead of the canonical `q`/`q_gpu`. This dummy
   !< replaces the retired `forest_active_substage` module-scope shared
   !< integer (issue #13, 2026-05-25).
   class(realm_object), intent(inout)           :: self     !< This realm.
   class(realm_object), intent(in)              :: realm(:) !< All realms in the forest.
   integer(I4P),        intent(in),    optional :: s_active !< Active RK substage (1..nrk) or 0/absent.

   if (int(size(realm), I4P) > 1_I4P) &
   call mpih%error_stop(msg='realm_object%exchange_inter_realm_halos_forest: not overridden by app extension')
   endsubroutine exchange_inter_realm_halos_forest

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
   class(realm_object),    intent(inout)        :: self    !< The equation.
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
   emin = [self%adam%field%emin(1,b)-self%ngc*self%adam%field%dxyz(1,b), &
           self%adam%field%emin(2,b)-self%ngc*self%adam%field%dxyz(2,b), &
           self%adam%field%emin(3,b)-self%ngc*self%adam%field%dxyz(3,b)]
   bn = 'block_'//trim(strz(b,9))//'-proc'//trim(strz(mpih%myrank,6))
   call xh5f%open_block(block_type = XH5F_PARAMETERS%XH5F_BLOCK_CARTESIAN_UNIFORM, &
                        block_name = bn,                                           &
                        nijk       = nijk,                                         &
                        emin       = emin,                                         &
                        dxyz       = self%adam%field%dxyz(:,b),                              &
                        time       = time_)
   call xh5f%save_block_field(xdmf_field_name = 'time_iteration',                                &
                              field           = t_,                                              &
                              field_center    = XDMF_PARAMETERS%XDMF_ATTR_CENTER_GRID,           &
                              field_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF, &
                              hdf5_field_name = bn//'-time_iteration')
   endsubroutine open_block_xh5f

   subroutine open_file_xh5f(self, basename, xh5f, directory)
   !< Open XH5F file.
   class(realm_object),    intent(inout)        :: self          !< The equation.
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
   class(realm_object), intent(inout)        :: self              !< The equation.
   character(*),        intent(in)           :: basename          !< Base name of output files.
   real(R8P),           intent(in)           :: q(1:,                   &
                                                  1-self%adam%grid%ngc:,&
                                                  1-self%adam%grid%ngc:,&
                                                  1-self%adam%grid%ngc:,&
                                                  1:)             !< Q-vector variables [nv,ni,nj,nk,nb].
   character(*),        intent(in), optional :: q_name(:)         !< Q-vector variables names [nv].
   character(*),        intent(in), optional :: directory         !< Directory name of output files.
   logical,             intent(in), optional :: with_ghost        !< Flag to save ghost cells.
   logical,             intent(in), optional :: with_cell_morton  !< Flag to save Morton code also in cells.
   integer(I4P),        intent(in), optional :: t                 !< Time iteration.
   real(R8P),           intent(in), optional :: time              !< Time.
   type(string), allocatable                 :: q_name_(:)        !< Q-vector variables names [nv].
   character(:), allocatable                 :: directory_        !< Directory name of output files, local var.
   logical                                   :: with_ghost_       !< Flag to save ghost cells, local var.
   logical                                   :: with_cell_morton_ !< Flag to save Morton code also in cells, local var.
   integer(I4P)                              :: t_                !< Time iteration, local var.
   real(R8P)                                 :: time_             !< Time, local var.
   integer(I4P)                              :: ngc               !< Ghost cells saved.
   integer(I4P)                              :: ijk(2,3)          !< Blocks extents.
   integer(I8P)                              :: nijk(3)           !< Blocks dimensions.
   real(R8P)                                 :: emin(3)           !< Minimum abscissa of current block.
   character(:), allocatable                 :: filename_hdf5     !< HDF5 file name.
   character(:), allocatable                 :: filename_xdmf     !< XDMF file name.
   character(:), allocatable                 :: bn                !< Block name.
   type(xh5f_file_object)                    :: xh5f              !< XH5F file handler.
   integer(I4P)                              :: i, b, v           !< Counter.

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
      ngc = self%adam%grid%ngc
   else
      ngc = 0_I4P
   endif
   associate(ni=>self%adam%grid%ni, nj=>self%adam%grid%nj, nk=>self%adam%grid%nk)
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
   do b=1, self%adam%field%blocks_number
      emin = [self%adam%field%emin(1,b)-ngc*self%adam%field%dxyz(1,b), &
              self%adam%field%emin(2,b)-ngc*self%adam%field%dxyz(2,b), &
              self%adam%field%emin(3,b)-ngc*self%adam%field%dxyz(3,b)]
      bn = 'block_'//trim(strz(b,9))//'-proc'//trim(strz(mpih%myrank,6))
      call xh5f%open_block(block_type = XH5F_PARAMETERS%XH5F_BLOCK_CARTESIAN_UNIFORM, &
                           block_name = bn,                                           &
                           nijk       = nijk,                                         &
                           emin       = emin,                                         &
                           dxyz       = self%adam%field%dxyz(:,b),                              &
                           time       = time_)
      call xh5f%save_block_field(xdmf_field_name = 'time_iteration',                                &
                                 field           = t_,                                              &
                                 field_center    = XDMF_PARAMETERS%XDMF_ATTR_CENTER_GRID,           &
                                 field_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF, &
                                 hdf5_field_name = bn//'-time_iteration')
      call self%io%save_field(xh5f=xh5f, grid=self%adam%grid, block_name=bn, ijk=ijk, nijk=nijk, q=q(:,:,:,:,b), q_name=q_name_)
      call xh5f%close_block
   enddo
   call xh5f%close_grid
   call xh5f%close_grid(grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION_ASYNC)
   call xh5f%close_file
   endsubroutine save_q_xh5f

   ! FDV operators numerical methods
   subroutine compute_block_total_variation_fd(self, hs, dxyz, ivar, q, tot_var_field, total_variation)
   !< Return the max of block total variation for a given var.
   class(realm_object), intent(in)    :: self                       !< The equation.
   integer(I4P),        intent(in)    :: hs                         !< FDV half stencil length.
   real(R8P),           intent(in)    :: dxyz(3)                    !< Space steps.
   integer(I4P),        intent(in)    :: ivar                       !< Index of first component of vec field.
   real(R8P),           intent(in)    :: q(1:,         &
                                           1-self%ngc:,&
                                           1-self%ngc:,&
                                           1-self%ngc:)             !< Field variables.
   real(R8P),           intent(inout) :: tot_var_field(1-self%ngc:,&
                                                       1-self%ngc:,&
                                                       1-self%ngc:) !< Total variation field on blocks.
   real(R8P),           intent(out)   :: total_variation            !< Max total variation on given block.
   real(R8P)                          :: gradient(3,3)              !< Gradient.
   real(R8P)                          :: tv                         !< Total variation buffer.
   integer(I4P)                       :: i,j,k                      !< Counter.

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
   integer(I4P),        intent(in)    :: hs                                              !< FDV half stencil length.
   integer(I4P),        intent(in)    :: ivar                                            !< Start index of variable of q.
   real(R8P),           intent(in)    :: q(   1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),           intent(inout) :: curl(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Curl.
   integer(I4P)                       :: i,j,k,b                                         !< Counter.

   associate(dxyz=>self%adam%field%dxyz)
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
   integer(I4P),        intent(in)    :: hs                                              !< FDV half stencil length.
   integer(I4P),        intent(in)    :: ivar                                            !< Start index of variable of q.
   real(R8P),           intent(in)    :: q(   1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),           intent(inout) :: curl(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Curl.
   integer(I4P)                       :: i,j,k,b                                         !< Counter.

   associate(dxyz=>self%adam%field%dxyz)
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
   integer(I4P),        intent(in)    :: hs                                           !< FDV half stencil length.
   integer(I4P),        intent(in)    :: dir                                          !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),        intent(in)    :: ivar                                         !< Start index of (vec.) variable of q.
   real(R8P),           intent(in)    :: q(1:, 1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Field variables.
   real(R8P),           intent(inout) :: dq_ds(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Derivative1, dq/ds.
   integer(I4P)                       :: i,j,k,b                                      !< Counter.

   associate(dxyz=>self%adam%field%dxyz)
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
   integer(I4P),        intent(in)    :: hs                                           !< FDV half stencil length.
   integer(I4P),        intent(in)    :: dir                                          !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),        intent(in)    :: ivar                                         !< Start index of (vec.) variable of q.
   real(R8P),           intent(in)    :: q(1:, 1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Field variables.
   real(R8P),           intent(inout) :: dq_ds(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Derivative1, dq/ds.
   integer(I4P)                       :: i,j,k,b                                      !< Counter.

   associate(dxyz=>self%adam%field%dxyz)
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
   integer(I4P),        intent(in)    :: hs                                             !< FDV half stencil length.
   integer(I4P),        intent(in)    :: dir                                            !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),        intent(in)    :: ivar                                           !< Start index of vec variable of q.
   real(R8P),           intent(in)    :: q(1:,   1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Field variables.
   real(R8P),           intent(inout) :: d2q_ds2(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Derivative2, d2q/ds2.
   integer(I4P)                       :: i,j,k,b                                        !< Counter.

   associate(dxyz=>self%adam%field%dxyz)
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
   integer(I4P),        intent(in)    :: hs                                              !< FDV half stencil length.
   integer(I4P),        intent(in)    :: dir                                             !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),        intent(in)    :: ivar                                            !< Start index of variable of q.
   real(R8P),           intent(in)    :: q(1:,   1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),           intent(inout) :: d2q_ds2(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Derivative2, d2q/ds2.
   integer(I4P)                       :: i,j,k,b                                         !< Counter.

   associate(dxyz=>self%adam%field%dxyz)
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
   integer(I4P),        intent(in)    :: hs                                              !< FDV half stencil length.
   integer(I4P),        intent(in)    :: dir                                             !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),        intent(in)    :: ivar                                            !< Start index of variable of q.
   real(R8P),           intent(in)    :: q(1:,   1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),           intent(inout) :: d4q_ds4(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Derivative4, d4q/ds4.
   integer(I4P)                       :: i,j,k,b                                         !< Counter.

   associate(dxyz=>self%adam%field%dxyz)
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
   integer(I4P),        intent(in)    :: hs                                                 !< FDV half stencil length.
   integer(I4P),        intent(in)    :: ivar                                               !< Start index of field of q.
   real(R8P),           intent(in)    :: q(1:,      1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),           intent(inout) :: divergence(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Divergence.
   integer(I4P)                       :: i,j,k,b                                            !< Counter.

   associate(dxyz=>self%adam%field%dxyz)
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
   integer(I4P),        intent(in)    :: hs                                                 !< FDV half stencil length.
   integer(I4P),        intent(in)    :: ivar                                               !< Start index of field of q.
   real(R8P),           intent(in)    :: q(1:,      1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),           intent(inout) :: divergence(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Divergence.
   integer(I4P)                       :: i,j,k,b                                            !< Counter.

   associate(dxyz=>self%adam%field%dxyz)
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
   integer(I4P),        intent(in)    :: hs                                                 !< FDV half stencil length.
   integer(I4P),        intent(in)    :: ivar                                               !< Index of scalar variable of q.
   real(R8P),           intent(in)    :: q(       1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Field variables.
   real(R8P),           intent(inout) :: gradient(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Gradient.
   integer(I4P)                       :: i, j, k, b                                         !< Counter.

   associate(dxyz=>self%adam%field%dxyz)
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
   integer(I4P),        intent(in)    :: hs                                                 !< FDV half stencil length.
   integer(I4P),        intent(in)    :: ivar                                               !< Index of scalar variable of q.
   real(R8P),           intent(in)    :: q(       1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Field variables.
   real(R8P),           intent(inout) :: gradient(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Gradient.
   integer(I4P)                       :: i, j, k, b                                         !< Counter.

   associate(dxyz=>self%adam%field%dxyz)
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
   integer(I4P),        intent(in)    :: hs                                                !< FDV half stencil length.
   integer(I4P),        intent(in)    :: ivar                                              !< Index of scalar variable of q.
   real(R8P),           intent(in)    :: q(     1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),           intent(inout) :: laplacian(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Gradient.
   integer(I4P)                       :: i, j, k, b                                        !< Counter.

   associate(dxyz=>self%adam%field%dxyz)
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
   integer(I4P),        intent(in)    :: hs                                                !< FDV half stencil length.
   integer(I4P),        intent(in)    :: ivar                                              !< Index of scalar variable of q.
   real(R8P),           intent(in)    :: q(     1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),           intent(inout) :: laplacian(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Gradient.
   integer(I4P)                       :: i, j, k, b                                        !< Counter.

   associate(dxyz=>self%adam%field%dxyz)
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
