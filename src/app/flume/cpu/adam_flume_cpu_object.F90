!< ADAM, FLUME CPU backend object.
module adam_flume_cpu_object
!< ADAM, FLUME CPU backend object.
!<
!< Implements the forest contract on the host (MPI + OpenMP). The dispatch procedure pointers are type components,
!< bound in `initialize_flume` by exhaustive `select case` with a fatal `case default` (issue #35, D-10).
!< Milestone 1 status: skeleton (P1), the space operator computes the auxiliary variables and a null residual.

! ADAM classes, libraries, parameters
use :: adam_flux_register_object, only : flux_register_object
use :: adam_parameters,           only : FEC_1_6_ARRAY
use :: adam_realm_object,         only : realm_object
use :: adam_rk_object,            only : RK_1, RK_2, RK_3, RK_SSP_11, RK_SSP_22, RK_SSP_33, RK_SSP_54
! ADAM singleton objects
use :: adam_mpih_global,          only : mpih
! FLUME modules
use :: adam_flume_common_library, only : flume_common_object, conservative_to_auxiliary, BC_EXTRAPOLATION, BC_INFLOW, &
                                         BC_WALL_INVISCID, IA_A, IA_U, IQ_RU, NV_AUX, SCHEME_SPACE_WENO
! third party modules
use :: mpi
use :: penf,                      only : I4P, R8P, str

implicit none
private
public :: flume_cpu_object

type, extends(flume_common_object) :: flume_cpu_object
   !< FLUME CPU backend object.
   ! fluxes data
   real(R8P), allocatable :: flx_f(:,:,:,:,:) !< X-face fluxes [nv, 0:ni, 1:nj, 1:nk, nb], face i = i+1/2.
   real(R8P), allocatable :: fly_f(:,:,:,:,:) !< Y-face fluxes [nv, 1:ni, 0:nj, 1:nk, nb], face j = j+1/2.
   real(R8P), allocatable :: flz_f(:,:,:,:,:) !< Z-face fluxes [nv, 1:ni, 1:nj, 0:nk, nb], face k = k+1/2.
   ! dispatch
   procedure(compute_residuals_interface), pass(self), pointer :: compute_residuals=>null() !< Space operator.
   procedure(integrate_interface),         pass(self), pointer :: integrate=>null()         !< Time operator.
   contains
      ! public methods
      procedure, pass(self) :: allocate_cpu            !< Allocate CPU data.
      procedure, pass(self) :: compute_conservation    !< Compute and save the conservation integrals.
      procedure, pass(self) :: compute_q_aux           !< Compute the auxiliary variables.
      procedure, pass(self) :: initialize_flume        !< Initialize the CPU backend.
      procedure, pass(self) :: save_residuals          !< Save residuals history.
      procedure, pass(self) :: save_simulation_data    !< Save fields, restart and diagnostics on their cadence.
      procedure, pass(self) :: set_boundary_conditions !< Set boundary conditions on the crown maps.
      procedure, pass(self) :: set_initial_conditions  !< Set initial conditions.
      procedure, pass(self) :: update_ghost            !< Update ghost cells: local, MPI, boundary conditions.
      ! forest methods
      procedure, pass(self) :: advance_one_step_forest      !< Advance one full step (fast path).
      procedure, pass(self) :: apply_reflux_to_stage_forest !< Apply the reflux correction.
      procedure, pass(self) :: begin_stage_forest           !< Begin an integrator stage (staged path).
      procedure, pass(self) :: close_step_forest            !< Close a step (staged path).
      procedure, pass(self) :: compute_local_dt_forest      !< Compute the local stability-limited time step.
      procedure, pass(self) :: end_stage_forest             !< End an integrator stage (staged path).
      procedure, pass(self) :: fill_seam_from_peer_forest   !< Fill inter-realm seam ghosts from a peer.
      procedure, pass(self) :: finalize_forest              !< Finalize the realm.
      procedure, pass(self) :: initialize_forest            !< Initialize the realm.
      procedure, pass(self) :: is_done_forest               !< Return true if the realm is done.
      procedure, pass(self) :: open_step_forest             !< Open a step (staged path).
      procedure, pass(self) :: post_step_forest             !< Post-step IO and diagnostics.
      procedure, pass(self) :: stages_per_step_forest       !< Return the integrator stages per step.
endtype flume_cpu_object

abstract interface
   subroutine compute_residuals_interface(self, q, dq, s, flux_register)
   !< Compute the residuals, space operator; `q` is `intent(inout)` because its ghost cells are filled inside.
   import :: flume_cpu_object, flux_register_object, I4P, R8P
   class(flume_cpu_object),     intent(inout)           :: self          !< The equation.
   real(R8P),                   intent(inout)           :: q(1:,         &
                                                             1-self%ngc:,&
                                                             1-self%ngc:,&
                                                             1-self%ngc:,&
                                                             1:)         !< Conservative variables.
   real(R8P),                   intent(inout)           :: dq(1:,         &
                                                              1-self%ngc:,&
                                                              1-self%ngc:,&
                                                              1-self%ngc:,&
                                                              1:)         !< Residuals.
   integer(I4P),                intent(in),    optional :: s             !< Runge-Kutta stage.
   class(flux_register_object), intent(inout), optional :: flux_register !< Forest's flux register for reflux.
   endsubroutine compute_residuals_interface

   subroutine integrate_interface(self)
   !< Integrate one time step, time operator.
   import :: flume_cpu_object
   class(flume_cpu_object), intent(inout) :: self !< The equation.
   endsubroutine integrate_interface
endinterface

contains
   ! public methods
   subroutine allocate_cpu(self)
   !< Allocate CPU data.
   class(flume_cpu_object), intent(inout) :: self       !< The equation.
   integer(I4P)                           :: alloc_stat !< Allocation status.
   character(999)                         :: alloc_msg  !< Allocation error message.

   associate(nv=>self%physics%nv, ni=>self%ni, nj=>self%nj, nk=>self%nk, nb=>self%nb)
   allocate(self%flx_f(1:nv,0:ni,1:nj,1:nk,1:nb), stat=alloc_stat, errmsg=alloc_msg)
   if (alloc_stat /= 0_I4P) call mpih%error_stop(msg=': failed to allocate flx_f: '//trim(alloc_msg))
   allocate(self%fly_f(1:nv,1:ni,0:nj,1:nk,1:nb), stat=alloc_stat, errmsg=alloc_msg)
   if (alloc_stat /= 0_I4P) call mpih%error_stop(msg=': failed to allocate fly_f: '//trim(alloc_msg))
   allocate(self%flz_f(1:nv,1:ni,1:nj,0:nk,1:nb), stat=alloc_stat, errmsg=alloc_msg)
   if (alloc_stat /= 0_I4P) call mpih%error_stop(msg=': failed to allocate flz_f: '//trim(alloc_msg))
   endassociate
   self%flx_f = 0._R8P
   self%fly_f = 0._R8P
   self%flz_f = 0._R8P
   endsubroutine allocate_cpu

   subroutine compute_conservation(self)
   !< Compute the volume integrals of the conservative variables and save them on the diagnostics cadence.
   class(flume_cpu_object), intent(inout) :: self         !< The equation.
   real(R8P)                              :: integrals(5) !< Volume integrals.
   real(R8P)                              :: volume       !< Cell volume.
   integer(I4P)                           :: b, i, j, k   !< Counters.

   if (.not.self%time%is_to_save(cadence=self%diagnostics%conservation_history_save)) return
   integrals = 0._R8P
   do b=1, self%blocks_number
      volume = product(self%adam%field%dxyz(:,b), mask=.not.self%adam%grid%null_xyz)
      do k=1, self%nk
         do j=1, self%nj
            do i=1, self%ni
               integrals = integrals + self%q(1:5,i,j,k,b) * volume
            enddo
         enddo
      enddo
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, integrals, 5, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih%error)
   call self%diagnostics%save_conservation_row(it=self%time%it, time=self%time%time, integrals=integrals)
   endsubroutine compute_conservation

   subroutine compute_q_aux(self, q)
   !< Compute the auxiliary variables on every cell, ghost cells included.
   class(flume_cpu_object), intent(inout) :: self              !< The equation.
   real(R8P),               intent(in)    :: q(1:,         &
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1:)             !< Conservative variables.
   real(R8P)                              :: gamma             !< Specific heats ratio.
   real(R8P)                              :: R                 !< Gas constant.
   integer(I4P)                           :: ngc, ni, nj, nk   !< Grid dimensions.
   integer(I4P)                           :: b, i, j, k        !< Counters.

   gamma = self%physics%gamma ; R = self%physics%R
   ngc = self%ngc ; ni = self%ni ; nj = self%nj ; nk = self%nk
   !$omp parallel do collapse(4) default(firstprivate) shared(self, q)
   do b=1, self%blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               call conservative_to_auxiliary(gamma=gamma, R=R, q=q(:,i,j,k,b), qa=self%q_aux(:,i,j,k,b))
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine compute_q_aux

   subroutine initialize_flume(self, filename, realms_number)
   !< Initialize the CPU backend: MPI, memory budget, common data, CPU data, dispatch.
   class(flume_cpu_object), intent(inout)        :: self           !< The equation.
   character(*),            intent(in)           :: filename       !< Input file name.
   integer(I4P),            intent(in), optional :: realms_number  !< Realm count; divides the per-process budget.
   integer(I4P)                                  :: realms_number_ !< Realm count, local variable.

   realms_number_ = 1_I4P ; if (present(realms_number)) realms_number_ = max(1_I4P, realms_number)
   call mpih%initialize(do_mpi_init=.true., verbose=.true.)
   call self%flume_common_object%initialize(filename=filename, memory_avail=mpih%memory_avail/realms_number_, &
                                            verbose=.true.)
   call self%allocate_cpu
   select case(self%numerics%scheme_space)
   case(SCHEME_SPACE_WENO)
      self%compute_residuals => compute_residuals_weno
   case default
      call mpih%error_stop(msg=': no CPU space operator for scheme_space "'//self%numerics%scheme_space//'"')
   endselect
   select case(self%rk%scheme)
   case(RK_1, RK_2, RK_3)
      self%integrate => integrate_rk_ls
   case(RK_SSP_11, RK_SSP_22, RK_SSP_33, RK_SSP_54)
      self%integrate => integrate_rk_ssp
   case default
      call mpih%error_stop(msg=': no CPU time integrator for [runge_kutta].(scheme) "'//trim(self%rk%scheme)//'"')
   endselect
   endsubroutine initialize_flume

   subroutine save_residuals(self)
   !< Save residuals history (L2 norm of dq, MPI-reduced, rank 0 writes).
   class(flume_cpu_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: v    !< Counter.

   if (.not.self%time%is_to_save(cadence=self%io%residuals_save)) return
   call self%adam%field%compute_normL2_residuals(grid=self%adam%grid, dq=self%dq, norm=self%adam%field%residuals)
   do v=1, self%nv
      call MPI_ALLREDUCE(MPI_IN_PLACE, self%adam%field%residuals(v), 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih%error)
      self%adam%field%residuals(v) = sqrt(self%adam%field%residuals(v)) / sqrt(real(self%ni*self%nj*self%nk, R8P))
   enddo
   if (mpih%myrank == 0) call self%io%save_residuals(it=self%time%it, time=self%time%time, &
                                                     blocks_number=self%blocks_number, residuals=self%adam%field%residuals)
   endsubroutine save_residuals

   subroutine save_simulation_data(self)
   !< Save fields, restart and conservation history, each on its own cadence (one predicate per output).
   class(flume_cpu_object), intent(inout) :: self !< The equation.

   if (self%time%is_to_save(cadence=self%io%it_save) .or. self%time%is_to_save(cadence=self%io%restart_save)) then
      call self%update_ghost(q=self%q)
      if (self%time%is_to_save(cadence=self%io%it_save)) call self%save_xh5f(with_ghost=.true.)
      if (self%time%is_to_save(cadence=self%io%restart_save)) call self%save_restart_files
   endif
   call self%compute_conservation
   endsubroutine save_simulation_data

   subroutine set_boundary_conditions(self, q)
   !< Set boundary conditions on the crown maps (face ghosts by kind, edge and corner ghosts by extrapolation).
   !<
   !< Periodic faces have no crown rows: their ghosts are filled by the ghost exchange. Edge and corner ghosts
   !< (fec > 6) are never read by the directional stencils; they are extrapolated so that every ghost holds a finite,
   !< deterministic state for the auxiliary variables computation.
   class(flume_cpu_object), intent(inout) :: self                   !< The equation.
   real(R8P),               intent(inout) :: q(1:,         &
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1:)                  !< Conservative variables.
   integer(I4P)                           :: b, c, i, j, k          !< Counters.
   integer(I4P)                           :: idelta, jdelta, kdelta !< IJK inward step.
   integer(I4P)                           :: bc_type                !< Boundary condition type.
   integer(I4P)                           :: crown                  !< Crown counter.
   integer(I4P)                           :: fec                    !< Boundary fec (1 to 26).
   integer(I4P)                           :: face                   !< Boundary face (1 to 6).
   integer(I4P)                           :: iref, jref, kref       !< Mirrored donor indexes.

   if (.not.allocated(self%adam%maps%local_map_bc_crown)) return
   associate(crown_map=>self%adam%maps%local_map_bc_crown, ni=>self%ni, nj=>self%nj, nk=>self%nk)
   do crown=1, self%ngc
      do c=1, size(crown_map, dim=1)
         b = int(crown_map(c,1,crown), I4P)
         if (b <= 0_I4P) cycle
         i       = int(crown_map(c,2,crown), I4P)
         j       = int(crown_map(c,3,crown), I4P)
         k       = int(crown_map(c,4,crown), I4P)
         idelta  = int(crown_map(c,5,crown), I4P)
         jdelta  = int(crown_map(c,6,crown), I4P)
         kdelta  = int(crown_map(c,7,crown), I4P)
         bc_type = int(crown_map(c,8,crown), I4P)
         fec     = int(crown_map(c,9,crown), I4P)
         if (fec > 6_I4P) then
            q(:,i,j,k,b) = q(:,i-idelta,j-jdelta,k-kdelta,b)
            cycle
         endif
         face = FEC_1_6_ARRAY(fec)
         select case(bc_type)
         case(BC_EXTRAPOLATION)
            q(:,i,j,k,b) = q(:,i-idelta,j-jdelta,k-kdelta,b)
         case(BC_INFLOW)
            q(:,i,j,k,b) = self%bc%q_inflow(:,face)
         case(BC_WALL_INVISCID)
            call compute_face_mirror_indexes(face=face, ni=ni, nj=nj, nk=nk, i_gc=i, j_gc=j, k_gc=k, &
                                             idelta=idelta, jdelta=jdelta, kdelta=kdelta, i_d=iref, j_d=jref, k_d=kref)
            q(:,i,j,k,b) = q(:,iref,jref,kref,b)
            q(IQ_RU+(face-1)/2,i,j,k,b) = -q(IQ_RU+(face-1)/2,i,j,k,b)
         case default
            call mpih%error_stop(msg=': unexpected boundary condition type '//trim(str(bc_type))//' on the crown map')
         endselect
      enddo
   enddo
   endassociate
   endsubroutine set_boundary_conditions

   subroutine set_initial_conditions(self)
   !< Set initial conditions.
   class(flume_cpu_object), intent(inout) :: self !< The equation.

   call self%ic%set_initial_conditions(field=self%adam%field, q=self%q)
   endsubroutine set_initial_conditions

   subroutine update_ghost(self, q)
   !< Update ghost cells: intra-realm local copies, MPI exchange, boundary conditions.
   !<
   !< Inter-realm seam ghosts are filled by the forest (`fill_seam_from_peer_forest`), not here.
   class(flume_cpu_object), intent(inout) :: self          !< The equation.
   real(R8P),               intent(inout) :: q(1:,         &
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1:)         !< Conservative variables.

   call self%adam%field%update_ghost_local(grid=self%adam%grid, maps=self%adam%maps, q=q)
   call self%adam%field%update_ghost_mpi(grid=self%adam%grid, maps=self%adam%maps, q=q)
   call self%set_boundary_conditions(q=q)
   endsubroutine update_ghost

   ! forest methods
   subroutine advance_one_step_forest(self, dt)
   !< Advance one full step of size `dt` (fast path: single realm, no AMR seam faces).
   class(flume_cpu_object), intent(inout) :: self !< The equation.
   real(R8P),               intent(in)    :: dt   !< Time step from the forest's global reduction.

   self%time%it = self%time%it + 1_I4P
   self%time%dt = dt
   if ((self%time%it_max <= 0_I4P) .and. (self%time%time + dt > self%time%time_max)) &
      self%time%dt = self%time%time_max - self%time%time
   call self%integrate
   self%time%time = self%time%time + self%time%dt
   call self%time%print_progress(nodes_number=self%adam%tree%nodes_number)
   endsubroutine advance_one_step_forest

   subroutine apply_reflux_to_stage_forest(self, stage, dt, flux_register)
   !< Apply the Berger-Colella reflux correction: not available before milestone 1 phase P5, so a run with AMR seam
   !< faces is refused rather than silently left non-conservative.
   class(flume_cpu_object),     intent(inout) :: self          !< The equation.
   integer(I4P),                intent(in)    :: stage         !< Integrator stage.
   real(R8P),                   intent(in)    :: dt            !< Time step.
   class(flux_register_object), intent(in)    :: flux_register !< Forest's flux register.

   if (.not.flux_register%is_initialized_) return
   if (flux_register%nfaces == 0_I4P) return
   call mpih%error_stop(msg=': AMR coarse-fine seams need reflux, not implemented yet (stage '//trim(str(stage))// &
                            ', dt '//trim(str(dt))//')')
   endsubroutine apply_reflux_to_stage_forest

   subroutine begin_stage_forest(self, k, K_total, dt, realm)
   !< Begin integrator stage `k` (staged path): publish the stage and compute its state.
   class(flume_cpu_object), intent(inout)                   :: self     !< The equation.
   integer(I4P),            intent(in)                      :: k        !< Stage index (1..K_total).
   integer(I4P),            intent(in)                      :: K_total  !< Forest-wide stage count for this step.
   real(R8P),               intent(in)                      :: dt       !< Time step from the forest.
   class(realm_object),     intent(inout), optional, target :: realm(:) !< Sibling realms (contract parity).

   self%stage_active = k
   call self%rk%compute_stage(field=self%adam%field, s=k, dt=self%time%dt)
   endsubroutine begin_stage_forest

   subroutine close_step_forest(self, dt)
   !< Close a step (staged path): assemble q, save residuals, advance time, clear the active stage.
   class(flume_cpu_object), intent(inout) :: self !< The equation.
   real(R8P),               intent(in)    :: dt   !< Time step from the forest (the local capped value is time%dt).

   call self%rk%update_q(field=self%adam%field, dt=self%time%dt, q=self%q, dq=self%dq)
   call self%save_residuals
   self%time%time = self%time%time + self%time%dt
   call self%time%print_progress(nodes_number=self%adam%tree%nodes_number)
   self%stage_active = 0_I4P
   endsubroutine close_step_forest

   subroutine compute_local_dt_forest(self, dt_local)
   !< Compute the local stability-limited time step, `dt = CFL / max(sum_d (|u_d| + a) / dx_d)` (no MPI reduction).
   !<
   !< The auxiliary variables are recomputed from the committed `q` (not read from `q_aux`, which holds the last
   !< stage state); null directions do not contribute.
   class(flume_cpu_object), intent(in)  :: self        !< The equation.
   real(R8P),               intent(out) :: dt_local    !< Local stability-limited time step.
   real(R8P)                            :: qa(NV_AUX)  !< Auxiliary variables of one cell.
   real(R8P)                            :: lambda_max  !< Maximum of sum_d (|u_d| + a) / dx_d.
   real(R8P)                            :: gamma       !< Specific heats ratio.
   real(R8P)                            :: R           !< Gas constant.
   logical                              :: is_null(3)  !< Null directions.
   integer(I4P)                         :: b, i, j, k  !< Counters.
   integer(I4P)                         :: d           !< Direction counter.

   gamma = self%physics%gamma ; R = self%physics%R ; is_null = self%adam%grid%null_xyz
   lambda_max = 0._R8P
   !$omp parallel do collapse(4) default(firstprivate) shared(self) reduction(max:lambda_max)
   do b=1, self%blocks_number
      do k=1, self%nk
         do j=1, self%nj
            do i=1, self%ni
               call conservative_to_auxiliary(gamma=gamma, R=R, q=self%q(:,i,j,k,b), qa=qa)
               lambda_max = max(lambda_max, sum([((abs(qa(IA_U+d-1)) + qa(IA_A)) / self%adam%field%dxyz(d,b), d=1, 3)], &
                                                mask=.not.is_null))
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   dt_local = huge(1._R8P)
   if (lambda_max > 0._R8P) dt_local = self%time%CFL / lambda_max
   endsubroutine compute_local_dt_forest

   subroutine end_stage_forest(self, k, K_total, dt, realm, flux_register)
   !< End integrator stage `k` (staged path): residuals on the stage state, then stage assignment.
   class(flume_cpu_object),     intent(inout)                   :: self          !< The equation.
   integer(I4P),                intent(in)                      :: k             !< Stage index (1..K_total).
   integer(I4P),                intent(in)                      :: K_total       !< Forest-wide stage count.
   real(R8P),                   intent(in)                      :: dt            !< Time step from the forest.
   class(realm_object),         intent(inout), optional, target :: realm(:)      !< Sibling realms (parity only).
   class(flux_register_object), intent(inout), optional         :: flux_register !< Forest's flux register.

   if (present(flux_register)) then
      call self%compute_residuals(q=self%rk%q_rk(:,:,:,:,:,k), dq=self%dq, s=k, flux_register=flux_register)
   else
      call self%compute_residuals(q=self%rk%q_rk(:,:,:,:,:,k), dq=self%dq, s=k)
   endif
   call self%rk%assign_stage(field=self%adam%field, s=k, q=self%dq)
   endsubroutine end_stage_forest

   subroutine fill_seam_from_peer_forest(self, peer, p_idx)
   !< Fill this realm's seam ghosts for peer slot `p_idx` from the peer's interior, on the active buffers (`q` when
   !< `stage_active == 0`, else the active stage of `q_rk`) of both realms.
   class(flume_cpu_object), intent(inout)         :: self                           !< The equation.
   class(realm_object),     intent(in),    target :: peer                           !< Peer realm.
   integer(I4P),            intent(in)            :: p_idx                          !< Peer slot.
   integer(I4P)                                   :: c, row_start, row_count, row   !< Counters and seam map rows.
   integer(I4P)                                   :: b_send, i_send, j_send, k_send !< Peer interior cell.
   integer(I4P)                                   :: b_recv, i_recv, j_recv, k_recv !< Own seam ghost cell.

   row_start = self%adam%maps%seam_local_peer_row_start(p_idx)
   row_count = self%adam%maps%seam_local_peer_row_count(p_idx)
   select type(peer)
   class is(flume_cpu_object)
      associate(rows=>self%adam%maps%seam_local_map_ghost_cell)
      do c=1_I4P, row_count
         row = row_start + c - 1_I4P
         b_send = rows(row,2) ; i_send = rows(row,4) ; j_send = rows(row,5) ; k_send = rows(row,6)
         b_recv = rows(row,3) ; i_recv = rows(row,7) ; j_recv = rows(row,8) ; k_recv = rows(row,9)
         if (self%stage_active > 0_I4P) then
            if (peer%stage_active > 0_I4P) then
               self%rk%q_rk(:,i_recv,j_recv,k_recv,b_recv,self%stage_active) = &
                  peer%rk%q_rk(:,i_send,j_send,k_send,b_send,peer%stage_active)
            else
               self%rk%q_rk(:,i_recv,j_recv,k_recv,b_recv,self%stage_active) = peer%q(:,i_send,j_send,k_send,b_send)
            endif
         else
            if (peer%stage_active > 0_I4P) then
               self%q(:,i_recv,j_recv,k_recv,b_recv) = peer%rk%q_rk(:,i_send,j_send,k_send,b_send,peer%stage_active)
            else
               self%q(:,i_recv,j_recv,k_recv,b_recv) = peer%q(:,i_send,j_send,k_send,b_send)
            endif
         endif
      enddo
      endassociate
   class default
      call mpih%error_stop(msg=': flume_cpu_object%fill_seam_from_peer_forest: peer realm is not a flume_cpu_object')
   endselect
   endsubroutine fill_seam_from_peer_forest

   subroutine finalize_forest(self)
   !< Finalize the realm: close the output files and free the data (MPI is finalized once by the forest).
   class(flume_cpu_object), intent(inout) :: self !< The equation.

   call self%io%close_file_residuals
   call self%diagnostics%close_file
   call self%destroy_common
   if (allocated(self%flx_f)) deallocate(self%flx_f)
   if (allocated(self%fly_f)) deallocate(self%fly_f)
   if (allocated(self%flz_f)) deallocate(self%flz_f)
   nullify(self%compute_residuals)
   nullify(self%integrate)
   endsubroutine finalize_forest

   subroutine initialize_forest(self, filename, realms_number, memory_avail, nv, verbose)
   !< Initialize the realm (issue #35, section 6.1): backend init, IC (or restart), initial AMR, ghost update, initial
   !< output, output files open, AMR lock.
   class(flume_cpu_object), intent(inout)           :: self          !< The equation.
   character(*),            intent(in)              :: filename      !< Input file name.
   integer(I4P),            intent(in),    optional :: realms_number !< Realm count; divides the per-process budget.
   real(R8P),               intent(in),    optional :: memory_avail  !< Unused: the budget comes from the MPI handler.
   integer(I4P),            intent(in),    optional :: nv            !< Unused: nv is decided by the physics.
   logical,                 intent(in),    optional :: verbose       !< Unused: initialization is always verbose.
   integer(I4P)                                     :: i             !< Counter.

   call self%initialize_flume(filename=filename, realms_number=realms_number)
   if (self%io%restart) then
      call mpih%print_message('restart simulation from "'//trim(self%io%restart_basename)//'" files')
      call self%load_restart_files(t=self%time%it, time=self%time%time)
   else
      do i=1, self%ic%amr_iterations
         call self%set_initial_conditions
         call self%amr_update
      enddo
      call self%set_initial_conditions
      call self%adam%make_comm_local_maps_ghost_bc
      self%time%time = 0._R8P
      self%time%it   = 0_I4P
   endif
   call self%update_ghost(q=self%q)
   call self%compute_q_aux(q=self%q)
   call self%diagnostics%open_file(output_basename=self%io%output_basename, q_name=self%q_name, &
                                   is_restart=self%io%restart)
   call self%save_simulation_data
   call self%io%open_file_residuals(nv=self%nv)
   self%amr_locked_ = .true.
   endsubroutine initialize_forest

   subroutine is_done_forest(self, done)
   !< Return true if the realm has reached its end.
   class(flume_cpu_object), intent(in)  :: self !< The equation.
   logical,                 intent(out) :: done !< True if the realm is done.

   done = self%time%is_done()
   endsubroutine is_done_forest

   subroutine open_step_forest(self, dt)
   !< Open a step (staged path): time bookkeeping and Runge-Kutta stages initialization.
   class(flume_cpu_object), intent(inout) :: self !< The equation.
   real(R8P),               intent(in)    :: dt   !< Time step from the forest.

   self%time%it = self%time%it + 1_I4P
   self%time%dt = dt
   if ((self%time%it_max <= 0_I4P) .and. (self%time%time + dt > self%time%time_max)) &
      self%time%dt = self%time%time_max - self%time%time
   call self%rk%initialize_stages(field=self%adam%field, q=self%q)
   endsubroutine open_step_forest

   subroutine post_step_forest(self, dt, t, it, do_save_state, do_save_residuals, do_save_restart, do_amr, realm)
   !< Post-step work: fields, restart and conservation history on their cadence.
   class(flume_cpu_object), intent(inout)                   :: self              !< The equation.
   real(R8P),               intent(in)                      :: dt                !< Time step just advanced.
   real(R8P),               intent(in)                      :: t                 !< Time after the advance.
   integer(I4P),            intent(in)                      :: it                !< Iteration after the advance.
   logical,                 intent(in),    optional         :: do_save_state     !< Unused: cadence is internal.
   logical,                 intent(in),    optional         :: do_save_residuals !< Unused: cadence is internal.
   logical,                 intent(in),    optional         :: do_save_restart   !< Unused: cadence is internal.
   logical,                 intent(in),    optional         :: do_amr            !< Unused: AMR is init-time only.
   class(realm_object),     intent(inout), optional, target :: realm(:)          !< Sibling realms.

   call self%save_simulation_data
   endsubroutine post_step_forest

   function stages_per_step_forest(self) result(K)
   !< Return the integrator stages per step: only SSP schemes are stage-splittable (staged path).
   class(flume_cpu_object), intent(in) :: self !< The equation.
   integer(I4P)                        :: K    !< Integrator stages per step.

   select case(self%rk%scheme)
   case(RK_SSP_11, RK_SSP_22, RK_SSP_33, RK_SSP_54)
      K = self%rk%nrk
   case default
      K = 0_I4P
      call mpih%error_stop(msg=': RK scheme "'//trim(self%rk%scheme)//'" is not stage-splittable: the staged forest '// &
                               'path (multi-realm, or AMR seam faces) requires an SSP scheme (runge-kutta-ssp-*)')
   endselect
   endfunction stages_per_step_forest

   ! private procedures
   subroutine compute_face_mirror_indexes(face, ni, nj, nk, i_gc, j_gc, k_gc, idelta, jdelta, kdelta, i_d, j_d, k_d)
   !< Return the donor indexes mirrored across a boundary face.
   integer(I4P), intent(in)  :: face                   !< Face index, 1 to 6.
   integer(I4P), intent(in)  :: ni, nj, nk             !< Grid dimensions.
   integer(I4P), intent(in)  :: i_gc, j_gc, k_gc       !< Ghost cell indexes.
   integer(I4P), intent(in)  :: idelta, jdelta, kdelta !< Inward steps.
   integer(I4P), intent(out) :: i_d, j_d, k_d          !< Mirrored donor indexes.

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
   case default
      call mpih%error_stop(msg=': compute_face_mirror_indexes: invalid face index '//trim(str(face)))
   endselect
   endsubroutine compute_face_mirror_indexes

   subroutine compute_residuals_weno(self, q, dq, s, flux_register)
   !< Compute the residuals with the WENO space operator.
   !<
   !< P1 skeleton: ghost update and auxiliary variables are computed, the residual is null. The WENO flux splitting
   !< lands in P3 (issue #35, section 3.4).
   class(flume_cpu_object),     intent(inout)           :: self          !< The equation.
   real(R8P),                   intent(inout)           :: q(1:,         &
                                                             1-self%ngc:,&
                                                             1-self%ngc:,&
                                                             1-self%ngc:,&
                                                             1:)         !< Conservative variables.
   real(R8P),                   intent(inout)           :: dq(1:,         &
                                                              1-self%ngc:,&
                                                              1-self%ngc:,&
                                                              1-self%ngc:,&
                                                              1:)         !< Residuals.
   integer(I4P),                intent(in),    optional :: s             !< Runge-Kutta stage.
   class(flux_register_object), intent(inout), optional :: flux_register !< Forest's flux register for reflux.

   call self%update_ghost(q=q)
   call self%compute_q_aux(q=q)
   dq = 0._R8P
   endsubroutine compute_residuals_weno

   subroutine integrate_rk_ls(self)
   !< Integrate one time step with a low-storage Runge-Kutta scheme.
   class(flume_cpu_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: s    !< Counter.

   call self%rk%initialize_stages(field=self%adam%field, q=self%q)
   do s=1, self%rk%nrk
      call self%compute_residuals(q=self%q, dq=self%dq)
      if (s == 1) call self%save_residuals
      call self%rk%compute_stage_ls(field=self%adam%field, s=s, dt=self%time%dt, dq=self%dq, q=self%q)
   enddo
   endsubroutine integrate_rk_ls

   subroutine integrate_rk_ssp(self)
   !< Integrate one time step with a strong stability preserving Runge-Kutta scheme.
   class(flume_cpu_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: s    !< Counter.

   call self%rk%initialize_stages(field=self%adam%field, q=self%q)
   do s=1, self%rk%nrk
      call self%rk%compute_stage(field=self%adam%field, s=s, dt=self%time%dt)
      call self%compute_residuals(q=self%rk%q_rk(:,:,:,:,:,s), dq=self%dq, s=s)
      call self%rk%assign_stage(field=self%adam%field, s=s, q=self%dq)
   enddo
   call self%rk%update_q(field=self%adam%field, dt=self%time%dt, q=self%q, dq=self%dq)
   call self%save_residuals
   endsubroutine integrate_rk_ssp
endmodule adam_flume_cpu_object
