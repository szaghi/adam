!< ADAM, FLUME CPU backend object.
module adam_flume_cpu_object
!< ADAM, FLUME CPU backend object.
!<
!< Implements the forest contract on the host (MPI + OpenMP). The dispatch procedure pointers are type components,
!< bound in `initialize_flume` by exhaustive `select case` with a fatal `case default` (issue #35, D-10).
!< Space operator: characteristic (or conservative) WENO flux splitting, face fluxes then flux difference (issue #35,
!< section 3.4); the per-model loops live in `adam_flume_cpu_<model>_kernels`, selected here by `select case` on the
!< physical model, never inside a loop (issue #41, section 4).

! ADAM classes, libraries, parameters
use :: adam_flux_register_object, only : face_tangential_ratios, flux_register_object
use :: adam_maps_object,          only : face_axis_sign
use :: adam_parameters,           only : BC_SEAM, FEC_1_6_ARRAY, FEC_TO_DELTA
use :: adam_realm_object,         only : realm_object
use :: adam_seam_exchange,        only : seam_fill_all, pack_seam_rows, unpack_seam_rows
use :: adam_rk_object,            only : RK_1, RK_2, RK_3, RK_SSP_11, RK_SSP_22, RK_SSP_33, RK_SSP_54
! ADAM singleton objects
use :: adam_mpih_global,          only : mpih
! FLUME modules
use :: adam_flume_common_library,      only : flume_common_object, flume_seam_sync_object, ib_cut_spacing,             &
                                              seam_face_cells, seam_fine_to_coarse, seam_skin_cell, seam_skin_index,   &
                                              BC_EXTRAPOLATION,                                                        &
                                              BC_INFLOW, BC_WALL_INVISCID, BC_WALL_ISOTHERMAL, BC_WALL_NOSLIP,          &
                                              realm_edge_donor, realm_edge_face, wall_noslip_ghost,                     &
                                              MODEL_EULER,                                                              &
                                              MODEL_MHD, MODEL_MHD_EGLM, MODEL_MHD_GLM, POSITIVITY_LIMITER_CELL,        &
                                              RECON_CHARACTERISTIC, RIEMANN_SOLVER_HLL,                                 &
                                              RIEMANN_SOLVER_HLLC, RIEMANN_SOLVER_HLLD, RIEMANN_SOLVER_LLF,            &
                                              SCHEME_SPACE_WENO,                                                        &
                                              SCHEME_SPACE_WENO_RIEMANN
use :: adam_flume_cpu_euler_hll_kernels,  only : compute_riemann_face_fluxes_euler_hll=>compute_riemann_face_fluxes
use :: adam_flume_cpu_euler_hllc_kernels, only : compute_riemann_face_fluxes_euler_hllc=>compute_riemann_face_fluxes
use :: adam_flume_cpu_euler_llf_kernels,  only : compute_riemann_face_fluxes_euler_llf=>compute_riemann_face_fluxes
use :: adam_flume_cpu_mhd_hll_kernels,          only : compute_riemann_face_fluxes_mhd_hll=>compute_riemann_face_fluxes
use :: adam_flume_cpu_mhd_hlld_kernels,         only : compute_riemann_face_fluxes_mhd_hlld=>compute_riemann_face_fluxes
use :: adam_flume_cpu_mhd_llf_kernels,          only : compute_riemann_face_fluxes_mhd_llf=>compute_riemann_face_fluxes
use :: adam_flume_cpu_mhd_glm_hll_kernels,      only : compute_riemann_face_fluxes_mhd_glm_hll=>compute_riemann_face_fluxes
use :: adam_flume_cpu_mhd_glm_hlld_kernels,     only : compute_riemann_face_fluxes_mhd_glm_hlld=>compute_riemann_face_fluxes
use :: adam_flume_cpu_mhd_glm_llf_kernels,      only : compute_riemann_face_fluxes_mhd_glm_llf=>compute_riemann_face_fluxes
use :: adam_flume_cpu_mhd_eglm_hll_kernels,     only : compute_riemann_face_fluxes_mhd_eglm_hll=>compute_riemann_face_fluxes
use :: adam_flume_cpu_mhd_eglm_hlld_kernels,    only : compute_riemann_face_fluxes_mhd_eglm_hlld=>compute_riemann_face_fluxes
use :: adam_flume_cpu_mhd_eglm_llf_kernels,     only : compute_riemann_face_fluxes_mhd_eglm_llf=>compute_riemann_face_fluxes
use :: adam_flume_cpu_dissipation_kernels, only : add_dissipative_fluxes_o2, add_dissipative_fluxes_o4
use :: adam_flume_cpu_euler_kernels,   only : blend_inadmissible_ghosts_euler=>blend_inadmissible_ghosts,            &
                                              blend_positivity_fluxes_euler=>blend_positivity_fluxes,                &
                                              compute_backbone_fluxes_euler=>compute_backbone_fluxes,                &
                                              compute_seam_positivity_factors_euler=>compute_seam_positivity_factors,&
                                              compute_face_fluxes_euler=>compute_face_fluxes,                        &
                                              compute_positivity_factors_euler=>compute_positivity_factors,          &
                                              compute_lambda_max_euler=>compute_lambda_max,                          &
                                              compute_lambda_max_dissipative_euler=>compute_lambda_max_dissipative,  &
                                              compute_q_aux_euler=>compute_q_aux,                                    &
                                              count_nonfinite_euler=>count_nonfinite
use :: adam_flume_cpu_mhd_kernels,     only : apply_floors_mhd=>apply_floors,                                   &
                                              blend_inadmissible_ghosts_mhd=>blend_inadmissible_ghosts,            &
                                              blend_positivity_fluxes_mhd=>blend_positivity_fluxes,                &
                                              compute_backbone_fluxes_mhd=>compute_backbone_fluxes,                &
                                              compute_seam_positivity_factors_mhd=>compute_seam_positivity_factors,&
                                              compute_positivity_factors_mhd=>compute_positivity_factors,          &
                                              compute_divb_norms_mhd=>compute_divb_norms,                          &
                                              compute_face_fluxes_mhd=>compute_face_fluxes,                        &
                                              compute_lambda_max_mhd=>compute_lambda_max,                            &
                                              compute_q_aux_mhd=>compute_q_aux,                                      &
                                              count_nonfinite_mhd=>count_nonfinite
use :: adam_flume_cpu_mhd_eglm_kernels, only : add_eglm_sources, add_eglm_sources_limited,                      &
                                               add_glm_damping_eglm=>add_glm_damping,                           &
                                               blend_inadmissible_ghosts_mhd_eglm=>blend_inadmissible_ghosts,   &
                                               blend_positivity_fluxes_mhd_eglm=>blend_positivity_fluxes,       &
                                               compute_backbone_fluxes_mhd_eglm=>compute_backbone_fluxes,       &
                                               compute_seam_positivity_factors_mhd_eglm=>                       &
                                                  compute_seam_positivity_factors,                             &
                                               compute_positivity_factors_mhd_eglm=>compute_positivity_factors, &
                                               apply_floors_mhd_eglm=>apply_floors,                             &
                                               compute_divb_norms_mhd_eglm=>compute_divb_norms,                 &
                                               compute_face_fluxes_mhd_eglm=>compute_face_fluxes,               &
                                               compute_lambda_max_mhd_eglm=>compute_lambda_max,                 &
                                               compute_q_aux_mhd_eglm=>compute_q_aux,                           &
                                               compute_speed_max_mhd_eglm=>compute_speed_max,                   &
                                               count_nonfinite_mhd_eglm=>count_nonfinite
use :: adam_flume_cpu_mhd_glm_kernels, only : add_glm_damping, apply_floors_mhd_glm=>apply_floors,              &
                                              blend_inadmissible_ghosts_mhd_glm=>blend_inadmissible_ghosts,        &
                                              compute_divb_norms_mhd_glm=>compute_divb_norms,                      &
                                              compute_face_fluxes_mhd_glm=>compute_face_fluxes,                      &
                                              compute_lambda_max_mhd_glm=>compute_lambda_max,                        &
                                              compute_q_aux_mhd_glm=>compute_q_aux,                                  &
                                              compute_speed_max_mhd_glm=>compute_speed_max,                          &
                                              count_nonfinite_mhd_glm=>count_nonfinite
! third party modules
use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
use :: mpi
use :: penf,                      only : I4P, I8P, R8P, str

implicit none
private
public :: flume_cpu_object

type, extends(flume_common_object) :: flume_cpu_object
   !< FLUME CPU backend object.
   ! fluxes data
   real(R8P), allocatable :: flx_f(:,:,:,:,:) !< X-face fluxes [nv, 0:ni, 1:nj, 1:nk, nb], face i = i+1/2.
   real(R8P), allocatable :: fly_f(:,:,:,:,:) !< Y-face fluxes [nv, 1:ni, 0:nj, 1:nk, nb], face j = j+1/2.
   real(R8P), allocatable :: flz_f(:,:,:,:,:) !< Z-face fluxes [nv, 1:ni, 1:nj, 0:nk, nb], face k = k+1/2.
   ! positivity limiter data
   type(flume_seam_sync_object) :: seam        !< Per-stage seam flux synchronisation of the limiter (issue #50).
   real(R8P), allocatable :: lam(:,:,:,:,:)   !< Cell factors of the positivity limiter (component 1, q-shaped for the
                                              !< ghost exchange) [nv, 1-ngc:ni+ngc, ..., nb]; allocated with the limiter.
   ! dispatch
   procedure(compute_residuals_interface), pass(self), pointer :: compute_residuals=>null() !< Space operator.
   procedure(integrate_interface),         pass(self), pointer :: integrate=>null()         !< Time operator.
   contains
      ! public methods
      procedure, pass(self) :: accumulate_seam_fluxes  !< Accumulate the weighted seam face fluxes of one stage.
      procedure, pass(self) :: add_dissipative_fluxes  !< Add the dissipative face fluxes (issue #65).
      procedure, pass(self) :: apply_floors            !< Apply the MHD positivity floors of a stage.
      procedure, pass(self) :: allocate_cpu            !< Allocate CPU data.
      procedure, pass(self) :: check_glm_ch            !< Check the GLM c_h against the fastest wave.
      procedure, pass(self) :: compute_conservation    !< Compute and save the conservation integrals.
      procedure, pass(self) :: compute_divb_history    !< Compute and save the div(B) norms (MHD).
      procedure, pass(self) :: compute_q_aux           !< Compute the auxiliary variables.
      procedure, pass(self) :: check_nonfinite         !< Stop on a non-finite committed state.
      procedure, pass(self) :: initialize_flume        !< Initialize the CPU backend.
      procedure, pass(self) :: limit_positivity        !< Apply the positivity limiter to the stage face fluxes.
      procedure, pass(self) :: seam_sync_blend         !< Set the seam face fluxes with the seam factors (issue #50).
      procedure, pass(self) :: seam_sync_factors       !< Recompute the seam cells' factors (issue #50).
      procedure, pass(self) :: seam_sync_fluxes        !< Publish the donor states and the fine means (issue #50).
      procedure, pass(self) :: seam_sync_theta         !< Compute the seam factors (issue #50).
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
      procedure, pass(self) :: pack_seam_cells_forest       !< Pack own cells of the cross-rank seam send rows.
      procedure, pass(self) :: unpack_seam_cells_forest     !< Unpack the cross-rank seam receive rows into own ghosts.
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
   subroutine accumulate_seam_fluxes(self, s, flux_register)
   !< Accumulate the seam face fluxes of stage `s`, weighted by its SSP coefficient, into the forest's flux register.
   !<
   !< Every stage contributes (issue #35, P5): the register holds `sum_s beta_s F_s`, the flux of the committed step.
   class(flume_cpu_object),     intent(inout) :: self          !< The equation.
   integer(I4P),                intent(in)    :: s             !< Runge-Kutta stage.
   class(flux_register_object), intent(inout) :: flux_register !< Forest's flux register.
   real(R8P), allocatable                     :: skin(:,:)     !< Face skin (nv, inner_n*outer_n).
   integer(I4P)                               :: b, fec        !< Block, face counters.
   integer(I4P)                               :: i, j, k, c    !< Cell counters.

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, nv=>self%nv)
   do b=1, self%blocks_number
      do fec=1, 6
         if (self%adam%maps%inter_realm_face_register_index(b, fec) == 0_I4P) cycle
         if (allocated(skin)) deallocate(skin)
         c = 0
         select case(fec)
         case(1_I4P, 2_I4P)
            allocate(skin(nv, nj*nk))
            i = merge(0_I4P, ni, fec == 1_I4P)
            do k=1, nk
               do j=1, nj
                  c = c + 1 ; skin(:,c) = self%flx_f(1:nv,i,j,k,b)
               enddo
            enddo
         case(3_I4P, 4_I4P)
            allocate(skin(nv, ni*nk))
            j = merge(0_I4P, nj, fec == 3_I4P)
            do k=1, nk
               do i=1, ni
                  c = c + 1 ; skin(:,c) = self%fly_f(1:nv,i,j,k,b)
               enddo
            enddo
         case default
            allocate(skin(nv, ni*nj))
            k = merge(0_I4P, nk, fec == 5_I4P)
            do j=1, nj
               do i=1, ni
                  c = c + 1 ; skin(:,c) = self%flz_f(1:nv,i,j,k,b)
               enddo
            enddo
         endselect
         call self%accumulate_seam_skin(flux_register=flux_register, b=b, fec=fec, weight=self%rk%beta(s), skin=skin)
      enddo
   enddo
   endassociate
   endsubroutine accumulate_seam_fluxes

   subroutine add_dissipative_fluxes(self)
   !< Add the dissipative face fluxes of the active directions to the inviscid ones (issue #65), from the auxiliary
   !< variables of the stage: the kernel of the order is selected here, never inside the loops.
   class(flume_cpu_object), intent(inout) :: self              !< The equation.
   real(R8P)                              :: tref              !< Reference temperature of the laws.
   real(R8P)                              :: omega_mu, omega_k !< Exponents of the laws.
   procedure(add_dissipative_fluxes_o4), pointer :: add_fluxes  !< Kernel of the order.

   call self%physics%dissipation%laws(tref=tref, omega_mu=omega_mu, omega_k=omega_k)
   select case(self%numerics%dissipative_order)
   case(2_I4P)
      add_fluxes => add_dissipative_fluxes_o2
   case(4_I4P)
      add_fluxes => add_dissipative_fluxes_o4
   case default
      call mpih%error_stop(msg=': no dissipative kernels of order '//trim(str(self%numerics%dissipative_order)))
   endselect
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nb=>self%blocks_number, &
             mu0=>self%physics%dissipation%mu, k0=>self%physics%dissipation%k, is_null=>self%adam%grid%null_xyz)
   if (.not.is_null(1)) call add_fluxes(d=1_I4P, di=1_I4P, dj=0_I4P, dk=0_I4P, ni=ni, nj=nj, nk=nk, ngc=ngc,            &
                                        blocks_number=nb, mu0=mu0, k0=k0, tref=tref, omega_mu=omega_mu, omega_k=omega_k, &
                                        dxyz=self%adam%field%dxyz, is_null=is_null, q_aux=self%q_aux, fl=self%flx_f)
   if (.not.is_null(2)) call add_fluxes(d=2_I4P, di=0_I4P, dj=1_I4P, dk=0_I4P, ni=ni, nj=nj, nk=nk, ngc=ngc,            &
                                        blocks_number=nb, mu0=mu0, k0=k0, tref=tref, omega_mu=omega_mu, omega_k=omega_k, &
                                        dxyz=self%adam%field%dxyz, is_null=is_null, q_aux=self%q_aux, fl=self%fly_f)
   if (.not.is_null(3)) call add_fluxes(d=3_I4P, di=0_I4P, dj=0_I4P, dk=1_I4P, ni=ni, nj=nj, nk=nk, ngc=ngc,            &
                                        blocks_number=nb, mu0=mu0, k0=k0, tref=tref, omega_mu=omega_mu, omega_k=omega_k, &
                                        dxyz=self%adam%field%dxyz, is_null=is_null, q_aux=self%q_aux, fl=self%flz_f)
   endassociate
   endsubroutine add_dissipative_fluxes

   subroutine apply_floors(self, q)
   !< Apply the MHD positivity floors to the interior of a stage state, before its ghost exchange (issue #41, 3.8).
   !<
   !< Euler has no floors (return before any work). MHD: the floored cells of the stage are logged by rank 0 when any;
   !< a non-positive density or pressure with the floors disabled (both zero) is fatal, reported with the global
   !< minimum density and pressure.
   class(flume_cpu_object), intent(inout) :: self      !< The equation.
   real(R8P),               intent(inout) :: q(1:,         &
                                                1-self%ngc:,&
                                                1-self%ngc:,&
                                                1-self%ngc:,&
                                                1:)         !< Conservative variables.
   integer(I4P)                           :: counts(2) !< Floored cells, non-positive cells.
   real(R8P)                              :: mins(2)   !< Minimum density and pressure.

   select case(self%physics%model)
   case(MODEL_EULER)
      return
   case(MODEL_MHD)
      call apply_floors_mhd(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number,    &
                                gamma=self%physics%gamma, R=self%physics%R, rho_floor=self%physics%mhd%rho_floor,   &
                                p_floor=self%physics%mhd%p_floor, q=q,               &
                                floored=counts(1), nonpositive=counts(2), rho_min=mins(1), p_min=mins(2))
   case(MODEL_MHD_GLM)
      call apply_floors_mhd_glm(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                                gamma=self%physics%gamma, R=self%physics%R, rho_floor=self%physics%mhd%rho_floor,   &
                                p_floor=self%physics%mhd%p_floor, q=q,                            &
                                floored=counts(1), nonpositive=counts(2), rho_min=mins(1), p_min=mins(2))
   case(MODEL_MHD_EGLM)
      call apply_floors_mhd_eglm(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                                gamma=self%physics%gamma, R=self%physics%R, rho_floor=self%physics%mhd%rho_floor,   &
                                p_floor=self%physics%mhd%p_floor, q=q,                            &
                                floored=counts(1), nonpositive=counts(2), rho_min=mins(1), p_min=mins(2))
   case default
      call mpih%error_stop(msg=': no floors for physical model "'//self%physics%physical_model//'"')
   endselect
   call MPI_ALLREDUCE(MPI_IN_PLACE, counts, 2, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpih%error)
   call MPI_ALLREDUCE(MPI_IN_PLACE, mins, 2, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, mpih%error)
   if (counts(2) > 0_I4P .and. .not.(self%physics%mhd%rho_floor > 0._R8P .or. self%physics%mhd%p_floor > 0._R8P)) &
      call mpih%error_stop(msg=': '//trim(str(counts(2)))//' cells with a non-positive density or pressure at step '// &
                               trim(str(self%time%it))//' (min rho '//trim(str(mins(1)))//', min p '//               &
                               trim(str(mins(2)))//'); the [mhd] floors rho_floor, p_floor are disabled')
   if (counts(1) > 0_I4P .and. mpih%myrank == 0) &
      print '(A)', mpih%myrankstr//'MHD floors: '//trim(str(counts(1)))//' cells floored at step '// &
                   trim(str(self%time%it))//' (min rho '//trim(str(mins(1)))//', min p '//trim(str(mins(2)))//')'
   endsubroutine apply_floors

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
   if (self%numerics%positivity_limiter == POSITIVITY_LIMITER_CELL) then
      allocate(self%lam(1:nv,1-self%ngc:ni+self%ngc,1-self%ngc:nj+self%ngc,1-self%ngc:nk+self%ngc,1:nb), &
               stat=alloc_stat, errmsg=alloc_msg)
      if (alloc_stat /= 0_I4P) call mpih%error_stop(msg=': failed to allocate lam: '//trim(alloc_msg))
      self%lam = 1._R8P
   endif
   endassociate
   self%flx_f = 0._R8P
   self%fly_f = 0._R8P
   self%flz_f = 0._R8P
   endsubroutine allocate_cpu

   subroutine check_glm_ch(self)
   !< Check the GLM cleaning speed against the fastest wave of the committed state (issue #41, section 3.5); a no-op
   !< without GLM.
   class(flume_cpu_object), intent(inout) :: self      !< The equation.
   real(R8P)                              :: speed_max !< Fastest wave speed of this rank.

   select case(self%physics%model)
   case(MODEL_EULER, MODEL_MHD)
      return
   case(MODEL_MHD_GLM)
      call compute_speed_max_mhd_glm(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                                     gamma=self%physics%gamma, R=self%physics%R, is_null=self%adam%grid%null_xyz,       &
                                     q=self%q, speed_max=speed_max)
   case(MODEL_MHD_EGLM)
      call compute_speed_max_mhd_eglm(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                                     gamma=self%physics%gamma, R=self%physics%R, is_null=self%adam%grid%null_xyz,       &
                                     q=self%q, speed_max=speed_max)
   case default
      call mpih%error_stop(msg=': no CPU kernels for physical model "'//self%physics%physical_model//'"')
   endselect
   call self%report_glm_speed(speed_max=speed_max)
   endsubroutine check_glm_ch

   subroutine check_nonfinite(self)
   !< Stop the run when the committed state holds a non-finite (NaN or infinite) value (issue #45): a run could otherwise
   !< finish with NaN fields and exit 0. One pass over the interior cells, the model selected outside the loops.
   class(flume_cpu_object), intent(inout) :: self !< The equation.
   integer(I8P)                           :: n    !< Non-finite values number of this rank.

   select case(self%physics%model)
   case(MODEL_EULER)
      call count_nonfinite_euler(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                                 q=self%q, n=n)
   case(MODEL_MHD)
      call count_nonfinite_mhd(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                               q=self%q, n=n)
   case(MODEL_MHD_GLM)
      call count_nonfinite_mhd_glm(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                                   q=self%q, n=n)
   case(MODEL_MHD_EGLM)
      call count_nonfinite_mhd_eglm(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                                   q=self%q, n=n)
   case default
      call mpih%error_stop(msg=': no CPU kernels for physical model "'//self%physics%physical_model//'"')
   endselect
   n = self%nonfinite_total(n_local=n)
   if (n > 0_I8P) call self%stop_nonfinite(n_total=n)
   endsubroutine check_nonfinite

   subroutine compute_conservation(self)
   !< Compute the volume integrals of the conservative variables and save them on the diagnostics cadence.
   !<
   !< The cell volume includes the null directions: the tree splits them too, so a refined block's cells are smaller
   !< along them, and a volume without them overweights the fine cells (issue #37).
   !<
   !< The sum is compensated (Neumaier): a plain running sum of a uniform field adds equal terms to a growing total, so
   !< its rounding errors do not cancel and grow with the cells number (M2-P7b: 1.3e-12 relative on the ~2.6e5 cells of
   !< a uniform-density AMR run, against a 1e-13 conservation criterion); the FNL twin reduces on the device (a tree).
   class(flume_cpu_object), intent(inout) :: self             !< The equation.
   real(R8P), allocatable                 :: integrals(:)     !< Volume integrals [nv].
   real(R8P), allocatable                 :: compensation(:)  !< Neumaier compensation of each integral [nv].
   real(R8P)                              :: volume           !< Cell volume.
   real(R8P)                              :: term, total      !< Cell contribution, updated sum.
   integer(I4P)                           :: b, i, j, k, v    !< Counters.

   if (.not.self%time%is_to_save(cadence=self%diagnostics%conservation_history_save)) return
   allocate(integrals(self%physics%nv), compensation(self%physics%nv))
   integrals    = 0._R8P
   compensation = 0._R8P
   do b=1, self%blocks_number
      volume = product(self%adam%field%dxyz(:,b))
      do k=1, self%nk
         do j=1, self%nj
            do i=1, self%ni
               do v=1, self%physics%nv
                  term  = self%q(v,i,j,k,b) * volume
                  total = integrals(v) + term
                  if (abs(integrals(v)) >= abs(term)) then
                     compensation(v) = compensation(v) + ((integrals(v) - total) + term)
                  else
                     compensation(v) = compensation(v) + ((term - total) + integrals(v))
                  endif
                  integrals(v) = total
               enddo
            enddo
         enddo
      enddo
   enddo
   integrals = integrals + compensation
   call MPI_ALLREDUCE(MPI_IN_PLACE, integrals, size(integrals), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih%error)
   call self%diagnostics%save_conservation_row(it=self%time%it, time=self%time%time*self%units%time_output(),        &
                                               integrals=integrals*self%output_factors(self%q_name)*            &
                                                         self%units%length_output()**3)
   endsubroutine compute_conservation

   subroutine compute_divb_history(self, realm)
   !< Compute and save the div(B) norms of the committed state on the diagnostics cadence (issue #41, D-10); a no-op
   !< for Euler. The ghost cells are refreshed first and, with sibling realms, the inter-realm seam ghosts are refilled
   !< from the peers (`update_ghost` does not fill them, issue #31): the stencils of the seam-local cells read them.
   class(flume_cpu_object), intent(inout)                   :: self     !< The equation.
   class(realm_object),     intent(inout), optional, target :: realm(:) !< Sibling realms.
   real(R8P)                                                :: norms(3) !< Norms of this rank.

   if (self%physics%model == MODEL_EULER) return
   if (.not.self%time%is_to_save(cadence=self%diagnostics%conservation_history_save)) return
   call self%update_ghost(q=self%q)
   if (present(realm)) call seam_fill_all(self=self, realm=realm) ! local and cross-rank seam rows (issue #40)
   call self%set_divb_seam
   associate(hs=>self%fdv_half_stencils(1))
   if (hs > self%ngc) call mpih%error_stop(msg=': the div(B) stencil ([fdv].(fdv_order)) exceeds the ghost cells')
   select case(self%physics%model)
   case(MODEL_MHD)
      call compute_divb_norms_mhd(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                                  hs=hs, band=self%ngc, dxyz=self%adam%field%dxyz, is_null=self%adam%grid%null_xyz,  &
                                  seam=self%divb_seam, q=self%q, divb_max=norms(1), divb_l1=norms(2),                &
                                  divb_seam_max=norms(3))
   case(MODEL_MHD_GLM)
      call compute_divb_norms_mhd_glm(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc,                             &
                                      blocks_number=self%blocks_number, hs=hs, band=self%ngc,                       &
                                      dxyz=self%adam%field%dxyz, is_null=self%adam%grid%null_xyz, seam=self%divb_seam, &
                                      q=self%q, divb_max=norms(1), divb_l1=norms(2), divb_seam_max=norms(3))
   case(MODEL_MHD_EGLM)
      call compute_divb_norms_mhd_eglm(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc,                             &
                                      blocks_number=self%blocks_number, hs=hs, band=self%ngc,                       &
                                      dxyz=self%adam%field%dxyz, is_null=self%adam%grid%null_xyz, seam=self%divb_seam, &
                                      q=self%q, divb_max=norms(1), divb_l1=norms(2), divb_seam_max=norms(3))
   case default
      call mpih%error_stop(msg=': no CPU kernels for physical model "'//self%physics%physical_model//'"')
   endselect
   endassociate
   call self%report_divb(norms=norms)
   endsubroutine compute_divb_history

   subroutine compute_q_aux(self, q)
   !< Compute the auxiliary variables on every cell, ghost cells included.
   class(flume_cpu_object), intent(inout) :: self              !< The equation.
   real(R8P),               intent(in)    :: q(1:,         &
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1:)             !< Conservative variables.

   select case(self%physics%model)
   case(MODEL_EULER)
      call compute_q_aux_euler(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                               gamma=self%physics%gamma, R=self%physics%R, q=q, q_aux=self%q_aux)
   case(MODEL_MHD)
      call compute_q_aux_mhd(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                             gamma=self%physics%gamma, R=self%physics%R, q=q, q_aux=self%q_aux)
   case(MODEL_MHD_GLM)
      call compute_q_aux_mhd_glm(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                                 gamma=self%physics%gamma, R=self%physics%R, q=q, q_aux=self%q_aux)
   case(MODEL_MHD_EGLM)
      call compute_q_aux_mhd_eglm(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                                 gamma=self%physics%gamma, R=self%physics%R, q=q, q_aux=self%q_aux)
   case default
      call mpih%error_stop(msg=': no CPU kernels for physical model "'//self%physics%physical_model//'"')
   endselect
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
   if (realms_number_ > 1_I4P .and. self%units%is_active) &
      call mpih%error_stop(msg=': [reference] is not supported on multi-realm runs (each realm converts its own '// &
                               'input and nothing yet checks that the realms share the references)')
   if (realms_number_ > 1_I4P .and. self%numerics%positivity_limiter == POSITIVITY_LIMITER_CELL) &
      call mpih%error_stop(msg=': [numerics].(positivity_limiter)=cell is not supported on multi-realm runs (the '// &
                               'inter-realm seam faces carry no limiting factor)')
   call self%allocate_cpu
   select case(self%numerics%scheme_space)
   case(SCHEME_SPACE_WENO)
      self%compute_residuals => compute_residuals_weno
   case(SCHEME_SPACE_WENO_RIEMANN)
      self%compute_residuals => compute_residuals_riemann
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

   subroutine limit_positivity(self, q, flux_register)
   !< Apply the positivity limiter to the face fluxes of the stage (issue #47, D-9): the cell factors of the model's
   !< kernel, their ghost exchange (intra-realm copies and MPI only: a ghost the exchange does not fill, at a physical
   !< boundary, keeps the factor 1 it was allocated with, so the face takes the interior cell's factor), the blending of
   !< the active directions. The limited faces and the cells with an inadmissible backbone are reduced over the ranks and
   !< logged. The forward-Euler step of the limiter is the step `time%dt` (every stage of an SSP scheme is a convex
   !< combination of forward-Euler steps of size at most `dt`).
   !<
   !< With 2:1 AMR seams (an intra-realm face in the forest's register) the seam flux is synchronised at this stage
   !< (issue #50, D1; `adam_flume_seam_sync_object`): the fine seam faces take the coarse donor state as backbone outer
   !< state, the coarse seam cells' factors use the means of the fine fluxes, and both sides of a seam face blend with
   !< `theta_s = min(Lambda_C, min Lambda_f)`, so the coarse flux is the mean of the fine ones and every cell stays in its
   !< corner box; the end-of-step reflux then corrects round-off only.
   class(flume_cpu_object),     intent(inout)           :: self          !< The equation.
   real(R8P),                   intent(in)              :: q(1:,         &
                                                             1-self%ngc:,&
                                                             1-self%ngc:,&
                                                             1-self%ngc:,&
                                                             1:)         !< Conservative variables of the stage.
   class(flux_register_object), intent(in),    optional :: flux_register !< Forest's flux register.
   integer(I4P)                                         :: counts(5)     !< Inadmissible backbones, limited faces per
                                                                         !< direction, cells with a non-finite flux.
   integer(I4P), allocatable                            :: ridx(:,:)     !< Register index when the map is unallocated.
   logical                                              :: seam          !< The stage synchronises seam faces.
   integer(I4P)                                         :: d             !< Counter.

   counts = 0_I4P
   seam = .false.
   if (present(flux_register)) then
      select type(flux_register)
      type is(flux_register_object)
         if (allocated(self%adam%maps%inter_realm_face_register_index)) then
            seam = self%seam%build(flux_register=flux_register,                                         &
                                   register_index=self%adam%maps%inter_realm_face_register_index,       &
                                   blocks_number=self%blocks_number, ni=self%ni, nj=self%nj, nk=self%nk, &
                                   nv=self%physics%nv)
         else
            allocate(ridx(self%blocks_number, 6)) ; ridx = 0_I4P
            seam = self%seam%build(flux_register=flux_register, register_index=ridx,                    &
                                   blocks_number=self%blocks_number, ni=self%ni, nj=self%nj, nk=self%nk, &
                                   nv=self%physics%nv)
         endif
      endselect
   endif
   if (seam) call self%seam_sync_fluxes(q=q)
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nb=>self%blocks_number, gamma=>self%physics%gamma, &
             ch=>self%physics%mhd%glm_ch, is_null=>self%adam%grid%null_xyz, dt=>self%time%dt, hs=>self%weno%S,      &
             dxyz=>self%adam%field%dxyz)
   select case(self%physics%model)
   case(MODEL_EULER)
      call compute_positivity_factors_euler(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, gamma=gamma, ch=ch,         &
                                            damping=0._R8P, hs=hs, dt=dt, dxyz=dxyz, is_null=is_null, q=q,             &
                                            q_aux=self%q_aux, flx=self%flx_f, fly=self%fly_f, flz=self%flz_f,          &
                                            lam=self%lam, bad=counts(1), nonfinite=counts(5))
   case(MODEL_MHD)
      call compute_positivity_factors_mhd(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, gamma=gamma, ch=ch,           &
                                          damping=0._R8P, hs=hs, dt=dt, dxyz=dxyz, is_null=is_null, q=q,               &
                                          q_aux=self%q_aux, flx=self%flx_f, fly=self%fly_f, flz=self%flz_f,            &
                                          lam=self%lam, bad=counts(1), nonfinite=counts(5))
   case(MODEL_MHD_EGLM)
      call compute_positivity_factors_mhd_eglm(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, gamma=gamma, ch=ch,      &
                                               damping=self%physics%mhd%glm_damping, hs=hs, dt=dt, dxyz=dxyz,          &
                                               is_null=is_null, q=q, q_aux=self%q_aux, flx=self%flx_f,                 &
                                               fly=self%fly_f, flz=self%flz_f, lam=self%lam, bad=counts(1),            &
                                               nonfinite=counts(5))
   case default
      call mpih%error_stop(msg=': no CPU positivity limiter for physical model "'//self%physics%physical_model//'"')
   endselect
   if (seam) call self%seam_sync_factors(q=q, dbad=counts(1), dnonfinite=counts(5))
   call self%adam%field%update_ghost_local(grid=self%adam%grid, maps=self%adam%maps, q=self%lam)
   call self%adam%field%update_ghost_mpi(grid=self%adam%grid, maps=self%adam%maps, q=self%lam)
   if (seam) call self%seam_sync_theta
   do d=1, 3
      if (is_null(d)) cycle
      select case(self%physics%model)
      case(MODEL_EULER)
         select case(d)
         case(1)
            call blend_positivity_fluxes_euler(d=d, di=1_I4P, dj=0_I4P, dk=0_I4P, ni=ni, nj=nj, nk=nk, ngc=ngc,      &
                                               blocks_number=nb, gamma=gamma, ch=ch, q=q, lam=self%lam,            &
                                               fl=self%flx_f, limited=counts(1+d))
         case(2)
            call blend_positivity_fluxes_euler(d=d, di=0_I4P, dj=1_I4P, dk=0_I4P, ni=ni, nj=nj, nk=nk, ngc=ngc,      &
                                               blocks_number=nb, gamma=gamma, ch=ch, q=q, lam=self%lam,            &
                                               fl=self%fly_f, limited=counts(1+d))
         case(3)
            call blend_positivity_fluxes_euler(d=d, di=0_I4P, dj=0_I4P, dk=1_I4P, ni=ni, nj=nj, nk=nk, ngc=ngc,      &
                                               blocks_number=nb, gamma=gamma, ch=ch, q=q, lam=self%lam,            &
                                               fl=self%flz_f, limited=counts(1+d))
         endselect
      case(MODEL_MHD)
         select case(d)
         case(1)
            call blend_positivity_fluxes_mhd(d=d, di=1_I4P, dj=0_I4P, dk=0_I4P, ni=ni, nj=nj, nk=nk, ngc=ngc,        &
                                             blocks_number=nb, gamma=gamma, ch=ch, q=q, lam=self%lam,              &
                                             fl=self%flx_f, limited=counts(1+d))
         case(2)
            call blend_positivity_fluxes_mhd(d=d, di=0_I4P, dj=1_I4P, dk=0_I4P, ni=ni, nj=nj, nk=nk, ngc=ngc,        &
                                             blocks_number=nb, gamma=gamma, ch=ch, q=q, lam=self%lam,              &
                                             fl=self%fly_f, limited=counts(1+d))
         case(3)
            call blend_positivity_fluxes_mhd(d=d, di=0_I4P, dj=0_I4P, dk=1_I4P, ni=ni, nj=nj, nk=nk, ngc=ngc,        &
                                             blocks_number=nb, gamma=gamma, ch=ch, q=q, lam=self%lam,              &
                                             fl=self%flz_f, limited=counts(1+d))
         endselect
      case(MODEL_MHD_EGLM)
         select case(d)
         case(1)
            call blend_positivity_fluxes_mhd_eglm(d=d, di=1_I4P, dj=0_I4P, dk=0_I4P, ni=ni, nj=nj, nk=nk, ngc=ngc,   &
                                                  blocks_number=nb, gamma=gamma, ch=ch, q=q, lam=self%lam,         &
                                                  fl=self%flx_f, limited=counts(1+d))
         case(2)
            call blend_positivity_fluxes_mhd_eglm(d=d, di=0_I4P, dj=1_I4P, dk=0_I4P, ni=ni, nj=nj, nk=nk, ngc=ngc,   &
                                                  blocks_number=nb, gamma=gamma, ch=ch, q=q, lam=self%lam,         &
                                                  fl=self%fly_f, limited=counts(1+d))
         case(3)
            call blend_positivity_fluxes_mhd_eglm(d=d, di=0_I4P, dj=0_I4P, dk=1_I4P, ni=ni, nj=nj, nk=nk, ngc=ngc,   &
                                                  blocks_number=nb, gamma=gamma, ch=ch, q=q, lam=self%lam,         &
                                                  fl=self%flz_f, limited=counts(1+d))
         endselect
      endselect
   enddo
   endassociate
   if (seam) call self%seam_sync_blend
   call MPI_ALLREDUCE(MPI_IN_PLACE, counts, 5, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpih%error)
   if (sum(counts) > 0_I4P .and. mpih%myrank == 0) &
      print '(A)', mpih%myrankstr//'positivity limiter: '//trim(str(sum(counts(2:4))))//' faces limited, '// &
                   trim(str(counts(1)))//' inadmissible backbones at step '//trim(str(self%time%it))
   if (counts(5) > 0_I4P .and. mpih%myrank == 0) &
      print '(A)', mpih%myrankstr//'positivity limiter (non-finite): '//trim(str(counts(5)))// &
                   ' cells with a non-finite high-order flux took the backbone at step '//trim(str(self%time%it))
   endsubroutine limit_positivity

   subroutine seam_sync_fluxes(self, q)
   !< Seam synchronisation, first phase (issue #50): publish the coarse donor states of the seam skins, then on every
   !< fine seam face compute the backbone with the donor state as outer state, keep it with the face's high-order flux,
   !< and accumulate the means of both into the coarse skins (2x2 fine faces per coarse face cell, 2x1 along the
   !< unrefined z of a quadtree, issue #46).
   class(flume_cpu_object), intent(inout) :: self                     !< The equation.
   real(R8P),               intent(in)    :: q(1:,         &
                                                 1-self%ngc:,&
                                                 1-self%ngc:,&
                                                 1-self%ngc:,&
                                                 1:)                    !< Conservative variables of the stage.
   real(R8P), allocatable                 :: qL(:,:), qR(:,:), fl(:,:) !< Face states, backbone fluxes.
   integer(I4P)                           :: b, fec, s, f, c, cc       !< Counters.
   integer(I4P)                           :: axis, sg, nc, off, foff   !< Face axis, side, cells, offsets.
   integer(I4P)                           :: ioff, joff                !< Fine quadrant offsets.
   integer(I4P)                           :: ratios(2)                 !< Tangential refinement ratios (inner, outer).
   real(R8P)                              :: w                         !< Mean weight, 1 / (ratios(1) ratios(2)).
   integer(I4P)                           :: i, j, k                   !< Cell indexes.

   if (.not.allocated(self%adam%maps%inter_realm_face_register_index)) then ! no block on this rank: collectives only
      call self%seam%reduce_states
      call self%seam%reduce_fluxes
      return
   endif
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, nv=>self%physics%nv, idx=>self%adam%maps%inter_realm_face_register_index)
   do b=1, self%blocks_number
      do fec=1, 6
         s = idx(b, fec)
         if (s <= 0_I4P) cycle
         off = self%seam%off(s)
         if (off < 0_I4P) cycle
         axis = (fec + 1_I4P) / 2_I4P ; sg = merge(-1_I4P, 1_I4P, mod(fec, 2_I4P) == 1_I4P)
         do c=1, seam_face_cells(fec=fec, ni=ni, nj=nj, nk=nk)
            call seam_skin_cell(axis=axis, sgn=sg, ni=ni, nj=nj, nk=nk, c=c, i=i, j=j, k=k)
            self%seam%qc(1:nv,off+c) = q(1:nv,i,j,k,b)
         enddo
      enddo
   enddo
   call self%seam%reduce_states
   do b=1, self%blocks_number
      do fec=1, 6
         foff = self%seam%fine_off(b, fec)
         if (foff < 0_I4P) cycle
         f = -idx(b, fec) ; off = self%seam%off(f)
         axis = (fec + 1_I4P) / 2_I4P ; sg = merge(-1_I4P, 1_I4P, mod(fec, 2_I4P) == 1_I4P)
         ioff = 0_I4P ; joff = 0_I4P
         if (allocated(self%adam%maps%amr_seam_quadrant)) then
            ioff = self%adam%maps%amr_seam_quadrant(1, b, fec) ; joff = self%adam%maps%amr_seam_quadrant(2, b, fec)
         endif
         ratios = face_tangential_ratios(fec=fec, refine_ratio=self%adam%maps%refine_ratio)
         w = 1._R8P / real(ratios(1) * ratios(2), R8P)
         nc = seam_face_cells(fec=fec, ni=ni, nj=nj, nk=nk)
         allocate(qL(nv,nc), qR(nv,nc), fl(nv,nc))
         do c=1, nc
            call seam_skin_cell(axis=axis, sgn=sg, ni=ni, nj=nj, nk=nk, c=c, i=i, j=j, k=k)
            cc = seam_fine_to_coarse(fec=fec, ni=ni, nj=nj, nk=nk, ioff=ioff, joff=joff, ri=ratios(1), ro=ratios(2), c=c)
            if (sg < 0_I4P) then
               qL(:,c) = self%seam%qc(1:nv,off+cc) ; qR(:,c) = q(1:nv,i,j,k,b)
            else
               qL(:,c) = q(1:nv,i,j,k,b) ; qR(:,c) = self%seam%qc(1:nv,off+cc)
            endif
            select case(axis)
            case(1)
               self%seam%fine_hi(1:nv,foff+c) = self%flx_f(1:nv,merge(0_I4P, ni, sg < 0_I4P),j,k,b)
            case(2)
               self%seam%fine_hi(1:nv,foff+c) = self%fly_f(1:nv,i,merge(0_I4P, nj, sg < 0_I4P),k,b)
            case default
               self%seam%fine_hi(1:nv,foff+c) = self%flz_f(1:nv,i,j,merge(0_I4P, nk, sg < 0_I4P),b)
            endselect
         enddo
         select case(self%physics%model)
         case(MODEL_EULER)
            call compute_backbone_fluxes_euler(n=nc, d=axis, gamma=self%physics%gamma, ch=self%physics%mhd%glm_ch, &
                                               qL=qL, qR=qR, flo=fl)
         case(MODEL_MHD)
            call compute_backbone_fluxes_mhd(n=nc, d=axis, gamma=self%physics%gamma, ch=self%physics%mhd%glm_ch, &
                                             qL=qL, qR=qR, flo=fl)
         case(MODEL_MHD_EGLM)
            call compute_backbone_fluxes_mhd_eglm(n=nc, d=axis, gamma=self%physics%gamma, ch=self%physics%mhd%glm_ch, &
                                                  qL=qL, qR=qR, flo=fl)
         endselect
         do c=1, nc
            cc = seam_fine_to_coarse(fec=fec, ni=ni, nj=nj, nk=nk, ioff=ioff, joff=joff, ri=ratios(1), ro=ratios(2), c=c)
            self%seam%fine_lo(1:nv,foff+c) = fl(:,c)
            self%seam%flo(1:nv,off+cc) = self%seam%flo(1:nv,off+cc) + w * fl(:,c)
            self%seam%fhi(1:nv,off+cc) = self%seam%fhi(1:nv,off+cc) + w * self%seam%fine_hi(1:nv,foff+c)
         enddo
         deallocate(qL, qR, fl)
      enddo
   enddo
   endassociate
   call self%seam%reduce_fluxes
   endsubroutine seam_sync_fluxes

   subroutine seam_sync_factors(self, q, dbad, dnonfinite)
   !< Seam synchronisation, second phase (issue #50): recompute the factor of every cell beside a seam face with that
   !< face's pair overridden, the coarse side by the fine means, the fine side by its donor-state backbone (and its own
   !< high-order flux), and add the change of the counts.
   class(flume_cpu_object), intent(inout) :: self                     !< The equation.
   real(R8P),               intent(in)    :: q(1:,         &
                                                 1-self%ngc:,&
                                                 1-self%ngc:,&
                                                 1-self%ngc:,&
                                                 1:)                    !< Conservative variables of the stage.
   integer(I4P),            intent(inout) :: dbad                      !< Inadmissible backbones count.
   integer(I4P),            intent(inout) :: dnonfinite                !< Non-finite cells count.
   integer(I4P), allocatable              :: cell(:,:)                 !< Seam cells (i, j, k, b).
   logical,      allocatable              :: omask(:,:)                !< Overridden faces.
   real(R8P),    allocatable              :: olo(:,:,:), ohi(:,:,:)    !< Override fluxes.
   logical                                :: mask(6)                   !< Seam faces of a cell.
   integer(I4P)                           :: n, m, b, i, j, k, fec, s  !< Counters.
   integer(I4P)                           :: c, db, dn                 !< Skin index, count changes.

   if (.not.allocated(self%adam%maps%inter_realm_face_register_index)) return
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, nv=>self%physics%nv, idx=>self%adam%maps%inter_realm_face_register_index)
   do m=1, 2 ! count, then fill
      n = 0
      do b=1, self%blocks_number
         if (all(idx(b,:) == 0_I4P)) cycle
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do fec=1, 6
                     s = idx(b, fec)
                     mask(fec) = .false.
                     if (s == 0_I4P .or. self%adam%grid%null_xyz((fec + 1) / 2)) cycle
                     if (self%seam%off(abs(s)) < 0_I4P) cycle
                     select case(fec)
                     case(1) ; mask(fec) = i == 1
                     case(2) ; mask(fec) = i == ni
                     case(3) ; mask(fec) = j == 1
                     case(4) ; mask(fec) = j == nj
                     case(5) ; mask(fec) = k == 1
                     case(6) ; mask(fec) = k == nk
                     endselect
                  enddo
                  if (.not.any(mask)) cycle
                  n = n + 1
                  if (m == 1) cycle
                  cell(:,n) = [i, j, k, b]
                  omask(:,n) = mask
                  do fec=1, 6
                     if (.not.mask(fec)) cycle
                     s = idx(b, fec)
                     c = seam_skin_index(fec=fec, ni=ni, nj=nj, i=i, j=j, k=k)
                     if (s > 0_I4P) then
                        olo(:,fec,n) = self%seam%flo(1:nv,self%seam%off(s)+c)
                        ohi(:,fec,n) = self%seam%fhi(1:nv,self%seam%off(s)+c)
                     else
                        olo(:,fec,n) = self%seam%fine_lo(1:nv,self%seam%fine_off(b,fec)+c)
                        ohi(:,fec,n) = self%seam%fine_hi(1:nv,self%seam%fine_off(b,fec)+c)
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
      if (m == 1) then
         allocate(cell(4,max(1,n)), omask(6,max(1,n)), olo(nv,6,max(1,n)), ohi(nv,6,max(1,n)))
         olo = 0._R8P ; ohi = 0._R8P
      endif
   enddo
   endassociate
   if (n == 0) return
   associate(ngc=>self%ngc, gamma=>self%physics%gamma, ch=>self%physics%mhd%glm_ch, is_null=>self%adam%grid%null_xyz, &
             dt=>self%time%dt, hs=>self%weno%S, dxyz=>self%adam%field%dxyz)
   select case(self%physics%model)
   case(MODEL_EULER)
      call compute_seam_positivity_factors_euler(ngc=ngc, gamma=gamma, ch=ch, damping=0._R8P, hs=hs, dt=dt, dxyz=dxyz, &
                                                 is_null=is_null, q=q, q_aux=self%q_aux, flx=self%flx_f, fly=self%fly_f, &
                                                 flz=self%flz_f, ncells=n, cell=cell, omask=omask, olo=olo, ohi=ohi,    &
                                                 lam=self%lam, dbad=db, dnonfinite=dn)
   case(MODEL_MHD)
      call compute_seam_positivity_factors_mhd(ngc=ngc, gamma=gamma, ch=ch, damping=0._R8P, hs=hs, dt=dt, dxyz=dxyz, &
                                               is_null=is_null, q=q, q_aux=self%q_aux, flx=self%flx_f, fly=self%fly_f, &
                                               flz=self%flz_f, ncells=n, cell=cell, omask=omask, olo=olo, ohi=ohi,    &
                                               lam=self%lam, dbad=db, dnonfinite=dn)
   case(MODEL_MHD_EGLM)
      call compute_seam_positivity_factors_mhd_eglm(ngc=ngc, gamma=gamma, ch=ch, damping=self%physics%mhd%glm_damping, &
                                                    hs=hs, dt=dt, dxyz=dxyz, is_null=is_null, q=q, q_aux=self%q_aux,   &
                                                    flx=self%flx_f, fly=self%fly_f, flz=self%flz_f, ncells=n,          &
                                                    cell=cell, omask=omask, olo=olo, ohi=ohi, lam=self%lam, dbad=db,   &
                                                    dnonfinite=dn)
   endselect
   endassociate
   dbad = dbad + db
   dnonfinite = dnonfinite + dn
   endsubroutine seam_sync_factors

   subroutine seam_sync_theta(self)
   !< Seam synchronisation, third phase (issue #50): the smallest factor of the fine cells under every coarse skin cell,
   !< then the seam factor `theta_s = min(Lambda_C, min Lambda_f)` in [0, 1] written by the coarse owner.
   class(flume_cpu_object), intent(inout) :: self               !< The equation.
   integer(I4P)                           :: b, fec, s, c, cc   !< Counters.
   integer(I4P)                           :: axis, sg, off, foff !< Face axis, side, offsets.
   integer(I4P)                           :: ioff, joff, i, j, k !< Fine quadrant offsets, cell indexes.
   integer(I4P)                           :: ratios(2)          !< Tangential refinement ratios (inner, outer).

   if (.not.allocated(self%adam%maps%inter_realm_face_register_index)) then ! no block on this rank: collectives only
      call self%seam%reduce_factors
      call self%seam%reduce_theta
      return
   endif
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, idx=>self%adam%maps%inter_realm_face_register_index)
   do b=1, self%blocks_number
      do fec=1, 6
         foff = self%seam%fine_off(b, fec)
         if (foff < 0_I4P) cycle
         off = self%seam%off(-idx(b, fec))
         axis = (fec + 1_I4P) / 2_I4P ; sg = merge(-1_I4P, 1_I4P, mod(fec, 2_I4P) == 1_I4P)
         ioff = 0_I4P ; joff = 0_I4P
         if (allocated(self%adam%maps%amr_seam_quadrant)) then
            ioff = self%adam%maps%amr_seam_quadrant(1, b, fec) ; joff = self%adam%maps%amr_seam_quadrant(2, b, fec)
         endif
         ratios = face_tangential_ratios(fec=fec, refine_ratio=self%adam%maps%refine_ratio)
         do c=1, seam_face_cells(fec=fec, ni=ni, nj=nj, nk=nk)
            call seam_skin_cell(axis=axis, sgn=sg, ni=ni, nj=nj, nk=nk, c=c, i=i, j=j, k=k)
            cc = seam_fine_to_coarse(fec=fec, ni=ni, nj=nj, nk=nk, ioff=ioff, joff=joff, ri=ratios(1), ro=ratios(2), &
                                     c=c)
            self%seam%lmin(off+cc) = min(self%seam%lmin(off+cc), self%lam(1,i,j,k,b))
         enddo
      enddo
   enddo
   call self%seam%reduce_factors
   do b=1, self%blocks_number
      do fec=1, 6
         s = idx(b, fec)
         if (s <= 0_I4P) cycle
         off = self%seam%off(s)
         if (off < 0_I4P) cycle
         axis = (fec + 1_I4P) / 2_I4P ; sg = merge(-1_I4P, 1_I4P, mod(fec, 2_I4P) == 1_I4P)
         do c=1, seam_face_cells(fec=fec, ni=ni, nj=nj, nk=nk)
            call seam_skin_cell(axis=axis, sgn=sg, ni=ni, nj=nj, nk=nk, c=c, i=i, j=j, k=k)
            self%seam%th(off+c) = max(0._R8P, min(1._R8P, self%lam(1,i,j,k,b), self%seam%lmin(off+c)))
         enddo
      enddo
   enddo
   endassociate
   call self%seam%reduce_theta
   endsubroutine seam_sync_theta

   subroutine seam_sync_blend(self)
   !< Seam synchronisation, last phase (issue #50): set the seam face fluxes with the seam factor, after the blend of
   !< the other faces. Coarse face `F_LF + theta_s (F_H - F_LF)` of the fine means, fine face the same of its own pair:
   !< the coarse flux is the mean of the fine ones (four, two on a quadtree face normal to x or y).
   class(flume_cpu_object), intent(inout) :: self                   !< The equation.
   real(R8P)                              :: fl(self%physics%nv)    !< Face flux.
   real(R8P)                              :: th                     !< Seam factor.
   integer(I4P)                           :: b, fec, s, c, cc       !< Counters.
   integer(I4P)                           :: axis, sg, off, foff    !< Face axis, side, offsets.
   integer(I4P)                           :: ioff, joff, i, j, k    !< Fine quadrant offsets, cell indexes.
   integer(I4P)                           :: ratios(2)              !< Tangential refinement ratios (inner, outer).

   if (.not.allocated(self%adam%maps%inter_realm_face_register_index)) return
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, nv=>self%physics%nv, idx=>self%adam%maps%inter_realm_face_register_index)
   do b=1, self%blocks_number
      do fec=1, 6
         s = idx(b, fec)
         if (s == 0_I4P) cycle
         off = self%seam%off(abs(s))
         if (off < 0_I4P) cycle
         axis = (fec + 1_I4P) / 2_I4P ; sg = merge(-1_I4P, 1_I4P, mod(fec, 2_I4P) == 1_I4P)
         foff = self%seam%fine_off(b, fec)
         ioff = 0_I4P ; joff = 0_I4P
         if (s < 0_I4P .and. allocated(self%adam%maps%amr_seam_quadrant)) then
            ioff = self%adam%maps%amr_seam_quadrant(1, b, fec) ; joff = self%adam%maps%amr_seam_quadrant(2, b, fec)
         endif
         ratios = face_tangential_ratios(fec=fec, refine_ratio=self%adam%maps%refine_ratio)
         do c=1, seam_face_cells(fec=fec, ni=ni, nj=nj, nk=nk)
            call seam_skin_cell(axis=axis, sgn=sg, ni=ni, nj=nj, nk=nk, c=c, i=i, j=j, k=k)
            if (s > 0_I4P) then
               th = self%seam%th(off+c)
               fl = self%seam%flo(1:nv,off+c)
               if (th > 0._R8P) fl = fl + th * (self%seam%fhi(1:nv,off+c) - fl)
            else
               cc = seam_fine_to_coarse(fec=fec, ni=ni, nj=nj, nk=nk, ioff=ioff, joff=joff, ri=ratios(1), ro=ratios(2), &
                                        c=c)
               th = self%seam%th(off+cc)
               fl = self%seam%fine_lo(1:nv,foff+c)
               if (th > 0._R8P) fl = fl + th * (self%seam%fine_hi(1:nv,foff+c) - fl)
            endif
            select case(axis)
            case(1)
               self%flx_f(1:nv,merge(0_I4P, ni, sg < 0_I4P),j,k,b) = fl
            case(2)
               self%fly_f(1:nv,i,merge(0_I4P, nj, sg < 0_I4P),k,b) = fl
            case default
               self%flz_f(1:nv,i,j,merge(0_I4P, nk, sg < 0_I4P),b) = fl
            endselect
         enddo
      enddo
   enddo
   endassociate
   endsubroutine seam_sync_blend

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
   if (mpih%myrank == 0) call self%io%save_residuals(it=self%time%it, time=self%time%time*self%units%time_output(), &
                                                     blocks_number=self%blocks_number,                             &
                                                     residuals=self%adam%field%residuals*self%output_factors(self%dq_name))
   endsubroutine save_residuals

   subroutine save_simulation_data(self, realm)
   !< Save fields, restart, slices and conservation history, each on its own cadence (one predicate per output).
   class(flume_cpu_object), intent(inout)                   :: self      !< The equation.
   class(realm_object),     intent(inout), optional, target :: realm(:)  !< Sibling realms.
   logical                                                  :: is_slices !< Slices save step.

   is_slices = self%slices%is_to_save(it=self%time%it, it_max=self%time%it_max, time=self%time%time, &
                                      time_max=self%time%time_max)
   if (self%time%is_to_save(cadence=self%io%it_save) .or. self%time%is_to_save(cadence=self%io%restart_save) .or. &
       is_slices) then
      ! the stage order: inter-realm seams first, then the intra-realm ghosts and the boundary conditions, so that the
      ! saved ghosts (edges and corners included) are those the residual stencils read; with [diagnostics]
      ! ghost_poison every ghost is NaN before the fill, so one the fill misses or reads too early shows (issue #65)
      if (self%diagnostics%ghost_poison) call poison_ghosts(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc,         &
                                                           blocks_number=self%blocks_number, q=self%q)
      if (present(realm)) call seam_fill_all(self=self, realm=realm)
      call self%update_ghost(q=self%q)
      if (self%time%is_to_save(cadence=self%io%it_save)) call self%save_xh5f(with_ghost=.true.)
      if (self%time%is_to_save(cadence=self%io%restart_save)) call self%save_restart_files
      if (is_slices) call self%save_slices
   endif
   call self%compute_conservation
   endsubroutine save_simulation_data

   subroutine set_boundary_conditions(self, q)
   !< Set boundary conditions on the crown maps, in three passes: the rows beyond one realm face, then beyond two (realm
   !< edges), then beyond three (realm corners).
   !<
   !< A row beyond one face applies that face's kind from the cell inward along its normal; the donor is interior or in
   !< a lower crown. This covers the face ghosts and the edge and corner ghosts whose other directions lead into
   !< another block (their tree boundary is a face). A row beyond several realm faces applies the kind of one physical
   !< face among them, inflow first, to the donor mirrored (wall) or clamped (extrapolation) along that face's axis.
   !< The donor lies beyond the other faces only, so it is a face ghost (for an edge row) or an edge ghost (for a
   !< corner row), already filled by the seams, the exchange or the earlier pass. The realm-edge ghosts then continue
   !< the face ghosts, as a cross derivative at a boundary face needs (issue #65 P0: the inward diagonal copy used
   !< before was off by O(1) at walls). A row beyond seams only lies outside the forest (a re-entrant corner of an
   !< L-shaped forest): it keeps a copy of its inward diagonal, finite, and read by no stencil of a fluid cell.
   !<
   !< Periodic faces have no crown rows: their ghosts are filled by the ghost exchange.
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
   integer(I4P)                           :: iref, jref, kref       !< Donor indexes.
   integer(I4P)                           :: pass                   !< Pass: realm faces the row lies beyond.
   integer(I4P)                           :: face_kind(6)           !< Kind of each realm face, BC_SEAM for a seam.
   real(R8P)                              :: psi_energy             !< 1 if the total energy holds psi^2 / 2 (EGLM).

   if (.not.allocated(self%adam%maps%local_map_bc_crown)) return
   psi_energy = merge(1._R8P, 0._R8P, self%physics%model == MODEL_MHD_EGLM)
   face_kind = self%bc%bc_type
   where (self%adam%maps%seam_face) face_kind = BC_SEAM
   associate(crown_map=>self%adam%maps%local_map_bc_crown, ni=>self%ni, nj=>self%nj, nk=>self%nk)
   do pass=1, 3
      do crown=1, self%ngc
         do c=1, size(crown_map, dim=1)
            b = int(crown_map(c,1,crown), I4P)
            if (b <= 0_I4P) cycle
            fec = int(crown_map(c,9,crown), I4P)
            if (count(FEC_TO_DELTA(:,fec) /= 0_I4P) /= pass) cycle
            i       = int(crown_map(c,2,crown), I4P)
            j       = int(crown_map(c,3,crown), I4P)
            k       = int(crown_map(c,4,crown), I4P)
            idelta  = int(crown_map(c,5,crown), I4P)
            jdelta  = int(crown_map(c,6,crown), I4P)
            kdelta  = int(crown_map(c,7,crown), I4P)
            bc_type = int(crown_map(c,8,crown), I4P)
            if (pass > 1_I4P) then
               face = realm_edge_face(fec=fec, face_kind=face_kind)
               if (face == 0_I4P) then ! beyond seams only: outside the forest
                  q(:,i,j,k,b) = q(:,i-idelta,j-jdelta,k-kdelta,b)
                  cycle
               endif
               call realm_edge_donor(face=face, face_kind=face_kind(face), ni=ni, nj=nj, nk=nk, i=i, j=j, k=k, &
                                     iref=iref, jref=jref, kref=kref)
               select case(face_kind(face))
               case(BC_INFLOW)
                  q(:,i,j,k,b) = self%bc%q_inflow(:,face)
               case(BC_WALL_INVISCID)
                  q(:,i,j,k,b) = self%bc%wall_sign(:,(face+1)/2) * q(:,iref,jref,kref,b)
               case(BC_WALL_NOSLIP, BC_WALL_ISOTHERMAL)
                  call wall_noslip_ghost(nv=self%nv, d=(face+1)/2, gamma=self%physics%gamma, R=self%physics%R,           &
                                         psi_energy=psi_energy, isothermal=face_kind(face) == BC_WALL_ISOTHERMAL,        &
                                         wall_velocity=self%bc%wall_velocity(:,face),                                   &
                                         wall_temperature=self%bc%wall_temperature(face), q=q(:,iref,jref,kref,b),      &
                                         qg=q(:,i,j,k,b))
               case default ! BC_EXTRAPOLATION
                  q(:,i,j,k,b) = q(:,iref,jref,kref,b)
               endselect
               cycle
            endif
            face = FEC_1_6_ARRAY(fec)
            select case(bc_type)
            case(BC_EXTRAPOLATION)
               ! the first interior cell along the normal, not the previous ghost: the donor then lies beyond block
               ! interfaces only (filled by the exchange), never in a row of this pass (issue #65 P1: the chain
               ! raced on the device, and on the host relied on the row order)
               call realm_edge_donor(face=face, face_kind=BC_EXTRAPOLATION, ni=ni, nj=nj, nk=nk, i=i, j=j, k=k, &
                                     iref=iref, jref=jref, kref=kref)
               q(:,i,j,k,b) = q(:,iref,jref,kref,b)
            case(BC_INFLOW)
               q(:,i,j,k,b) = self%bc%q_inflow(:,face)
            case(BC_WALL_INVISCID)
               call compute_face_mirror_indexes(face=face, ni=ni, nj=nj, nk=nk, i_gc=i, j_gc=j, k_gc=k,                     &
                                                idelta=idelta, jdelta=jdelta, kdelta=kdelta, i_d=iref, j_d=jref, k_d=kref)
               q(:,i,j,k,b) = self%bc%wall_sign(:,(face+1)/2) * q(:,iref,jref,kref,b)
            case(BC_WALL_NOSLIP, BC_WALL_ISOTHERMAL)
               call compute_face_mirror_indexes(face=face, ni=ni, nj=nj, nk=nk, i_gc=i, j_gc=j, k_gc=k,                     &
                                                idelta=idelta, jdelta=jdelta, kdelta=kdelta, i_d=iref, j_d=jref, k_d=kref)
               call wall_noslip_ghost(nv=self%nv, d=(face+1)/2, gamma=self%physics%gamma, R=self%physics%R,                &
                                      psi_energy=psi_energy, isothermal=bc_type == BC_WALL_ISOTHERMAL,                     &
                                      wall_velocity=self%bc%wall_velocity(:,face),                                        &
                                      wall_temperature=self%bc%wall_temperature(face), q=q(:,iref,jref,kref,b),           &
                                      qg=q(:,i,j,k,b))
            case(BC_SEAM)
               ! inter-realm seam face: filled by the forest (fill_seam_from_peer_forest), nothing to do here
            case default
               call mpih%error_stop(msg=': unexpected boundary condition type '//trim(str(bc_type))//' on the crown map')
            endselect
         enddo
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
   !< Update ghost cells: intra-realm local copies, MPI exchange, boundary conditions, then the positivity blend of the
   !< inadmissible face ghosts (issue #50, D2: a coarse-to-fine seam interpolant of admissible states may be inadmissible;
   !< every other ghost is a copy or a boundary state, admissible, and untouched).
   !<
   !< Inter-realm seam ghosts are filled by the forest (`fill_seam_from_peer_forest`), not here.
   class(flume_cpu_object), intent(inout) :: self          !< The equation.
   real(R8P),               intent(inout) :: q(1:,         &
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1:)         !< Conservative variables.
   integer(I4P)                           :: blended       !< Face ghosts blended toward the interior.

   call self%adam%field%update_ghost_local(grid=self%adam%grid, maps=self%adam%maps, q=q)
   call self%adam%field%update_ghost_mpi(grid=self%adam%grid, maps=self%adam%maps, q=q)
   call self%set_boundary_conditions(q=q)
   blended = 0_I4P
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nb=>self%blocks_number, is_null=>self%adam%grid%null_xyz)
   select case(self%physics%model)
   case(MODEL_EULER)
      call blend_inadmissible_ghosts_euler(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, is_null=is_null, q=q, &
                                           blended=blended)
   case(MODEL_MHD)
      call blend_inadmissible_ghosts_mhd(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, is_null=is_null, q=q, &
                                         blended=blended)
   case(MODEL_MHD_GLM)
      call blend_inadmissible_ghosts_mhd_glm(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, is_null=is_null, q=q, &
                                             blended=blended)
   case(MODEL_MHD_EGLM)
      call blend_inadmissible_ghosts_mhd_eglm(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, is_null=is_null, q=q, &
                                              blended=blended)
   endselect
   endassociate
   if (blended > 0_I4P) print '(A)', mpih%myrankstr//'ghost positivity: '//trim(str(blended))// &
                                     ' inadmissible face ghosts blended toward the interior at step '//trim(str(self%time%it))
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
   !< Apply the Berger-Colella reflux correction to the committed `q` (the forest calls it once per step, after
   !< `close_step_forest`, with the final stage).
   !<
   !< For every register face whose coarse side this realm and this rank own, the coarse skin cells receive
   !< `sgn dt / dx_coarse (F_coarse - F_fine_sum)`: the register holds the step fluxes `sum_s beta_s F_s` of both sides
   !< (`accumulate_seam_fluxes`), so the coarse flux the step used is replaced by the restricted fine one, and the
   !< coarse-fine interface conserves to round-off.
   !<
   !< The scale uses the step the update used, `time%dt`: on the last step of a time-driven run the realm caps it to
   !< land on `time_max`, and the forest's `dt` argument is the uncapped value (with it the correction was off by the
   !< ratio of the two, 1.3e-9 of the mass on sod-amr; issue #37).
   class(flume_cpu_object),     intent(inout) :: self          !< The equation.
   integer(I4P),                intent(in)    :: stage         !< Integrator stage.
   real(R8P),                   intent(in)    :: dt            !< Forest time step (unused: see the scale below).
   class(flux_register_object), intent(in)    :: flux_register !< Forest's flux register.
   real(R8P)                                  :: scale         !< Correction scale, sgn dt / dx_coarse.
   integer(I4P)                               :: f, c          !< Face, skin cell counters.
   integer(I4P)                               :: axis, sgn     !< Face normal axis and side.
   integer(I4P)                               :: i, j, k       !< Coarse cell indexes.

   if (.not.self%numerics%reflux) return
   if (.not.flux_register%is_initialized_) return
   if (flux_register%nfaces == 0_I4P) return
   if (.not.allocated(flux_register%face)) return
   if (stage /= self%rk%nrk) return
   do f=1, flux_register%nfaces
      associate(face=>flux_register%face(f))
      if (face%coarse_realm /= self%realm_index .or. face%coarse_rank /= mpih%myrank) cycle
      if (.not.(allocated(face%F_coarse) .and. allocated(face%F_fine_sum))) cycle
      call face_axis_sign(face_code=face%coarse_face, axis=axis, sgn=sgn)
      if (axis == 0_I4P) call mpih%error_stop(msg=': malformed coarse face code of register face '//trim(str(f)))
      scale = real(sgn, R8P) * self%time%dt / self%adam%field%dxyz(axis,face%coarse_block)
      do c=1, face%nface_cells
         call seam_skin_cell(axis=axis, sgn=sgn, ni=self%ni, nj=self%nj, nk=self%nk, c=c, i=i, j=j, k=k)
         self%q(:,i,j,k,face%coarse_block) = self%q(:,i,j,k,face%coarse_block) + &
                                             scale * (face%F_coarse(1:self%nv,c,1) - face%F_fine_sum(1:self%nv,c,1))
      enddo
      endassociate
   enddo
   endsubroutine apply_reflux_to_stage_forest

   subroutine begin_stage_forest(self, k, K_total, dt, realm)
   !< Begin integrator stage `k` (staged path): publish the stage and compute its state.
   class(flume_cpu_object), intent(inout)                   :: self     !< The equation.
   integer(I4P),            intent(in)                      :: k        !< Stage index (1..K_total).
   integer(I4P),            intent(in)                      :: K_total  !< Forest-wide stage count for this step.
   real(R8P),               intent(in)                      :: dt       !< Time step from the forest.
   class(realm_object),     intent(inout), optional, target :: realm(:) !< Sibling realms (contract parity).

   self%stage_active = k
   call self%rk%compute_stage(field=self%adam%field, s=k, dt=self%time%dt, phi=self%ib%phi)
   endsubroutine begin_stage_forest

   subroutine close_step_forest(self, dt)
   !< Close a step (staged path): assemble q, save residuals, advance time, clear the active stage.
   class(flume_cpu_object), intent(inout) :: self !< The equation.
   real(R8P),               intent(in)    :: dt   !< Time step from the forest (the local capped value is time%dt).

   call self%rk%update_q(field=self%adam%field, dt=self%time%dt, phi=self%ib%phi, q=self%q, dq=self%dq)
   if (allocated(self%ib%phi)) call compute_rk_ssp_residual(self)
   call self%save_residuals
   self%time%time = self%time%time + self%time%dt
   call self%time%print_progress(nodes_number=self%adam%tree%nodes_number)
   self%stage_active = 0_I4P
   endsubroutine close_step_forest

   subroutine compute_local_dt_forest(self, dt_local)
   !< Compute the local stability-limited time step, `dt = CFL / max(sum_d (|u_d| + a) / dx_d)` (no MPI reduction);
   !< with GLM, also `dt <= CFL / (c_h max sum_d 1 / dx_d)` (`glm_lambda`, issue #41, section 3.5). With dissipative
   !< terms the cell sum gains `2 nu sum_d 1 / dx_d^2`, `nu = max(4/3 mu / rho, gamma k / (rho cp), eta)` (issue #65,
   !< D-M4-3), and the first step logs which limit dominates.
   !<
   !< The auxiliary variables are recomputed from the committed `q` (not read from `q_aux`, which holds the last
   !< stage state); null directions do not contribute.
   class(flume_cpu_object), intent(in)  :: self        !< The equation.
   real(R8P),               intent(out) :: dt_local    !< Local stability-limited time step.
   real(R8P)                            :: lambda_max  !< Maximum of sum_d (|u_d| + a) / dx_d (+ the diffusive part).
   real(R8P)                            :: lambda_hyp  !< Maximum of the hyperbolic part (dissipative runs).
   real(R8P)                            :: lambda_dif  !< Maximum of the diffusive part (dissipative runs).
   real(R8P)                            :: re_cell_min !< Minimum cell Reynolds number (dissipative runs).
   real(R8P)                            :: tref        !< Reference temperature of the laws.
   real(R8P)                            :: omega_mu    !< Viscosity law exponent.
   real(R8P)                            :: omega_k     !< Conductivity law exponent.

   select case(self%physics%model)
   case(MODEL_EULER)
      if (self%physics%dissipation%is_active) then
         call self%physics%dissipation%laws(tref=tref, omega_mu=omega_mu, omega_k=omega_k)
         call compute_lambda_max_dissipative_euler(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc,                 &
                                                   blocks_number=self%blocks_number, gamma=self%physics%gamma,       &
                                                   R=self%physics%R, cp=self%physics%cp,                             &
                                                   mu0=self%physics%dissipation%mu, k0=self%physics%dissipation%k,   &
                                                   eta=self%physics%dissipation%eta, t_wall=self%wall_temperature_bound(), &
                                                   tref=tref, omega_mu=omega_mu,                                     &
                                                   omega_k=omega_k, dxyz=self%adam%field%dxyz,                       &
                                                   is_null=self%adam%grid%null_xyz, q=self%q, lambda_max=lambda_max, &
                                                   lambda_hyp=lambda_hyp, lambda_dif=lambda_dif,                     &
                                                   re_cell_min=re_cell_min)
         if (self%time%it == 0_I4P) call self%log_dt_limit(lambda_hyp=lambda_hyp, lambda_dif=lambda_dif, &
                                                            re_cell_min=re_cell_min)
      else
         call compute_lambda_max_euler(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc,                              &
                                       blocks_number=self%blocks_number, gamma=self%physics%gamma, R=self%physics%R, &
                                       dxyz=self%adam%field%dxyz, is_null=self%adam%grid%null_xyz, q=self%q,         &
                                       lambda_max=lambda_max)
      endif
   case(MODEL_MHD)
      call compute_lambda_max_mhd(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                                  gamma=self%physics%gamma, R=self%physics%R, dxyz=self%adam%field%dxyz,              &
                                  is_null=self%adam%grid%null_xyz, q=self%q, lambda_max=lambda_max)
   case(MODEL_MHD_GLM)
      call compute_lambda_max_mhd_glm(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc,                                &
                                      blocks_number=self%blocks_number, gamma=self%physics%gamma, R=self%physics%R, &
                                      dxyz=self%adam%field%dxyz, is_null=self%adam%grid%null_xyz, q=self%q,         &
                                      lambda_max=lambda_max)
      lambda_max = max(lambda_max, self%glm_lambda())
   case(MODEL_MHD_EGLM)
      call compute_lambda_max_mhd_eglm(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc,                                &
                                      blocks_number=self%blocks_number, gamma=self%physics%gamma, R=self%physics%R, &
                                      dxyz=self%adam%field%dxyz, is_null=self%adam%grid%null_xyz, q=self%q,         &
                                      lambda_max=lambda_max)
      lambda_max = max(lambda_max, self%glm_lambda())
   case default
      call mpih%error_stop(msg=': no CPU kernels for physical model "'//self%physics%physical_model//'"')
   endselect
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
   call self%rk%assign_stage(field=self%adam%field, s=k, q=self%dq, phi=self%ib%phi)
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

   subroutine pack_seam_cells_forest(self, p_idx, buf)
   !< Pack the ghost values of the cross-rank seam send rows of peer slot `p_idx` from the active buffer (`q` when
   !< `stage_active == 0`, else the active stage of `q_rk`): cell copies, or the 2:1 interpolations and restrictions of a
   !< refined seam (issues #40, #52; `adam_seam_exchange%pack_seam_rows`).
   class(flume_cpu_object), intent(in)  :: self   !< The realm.
   integer(I4P),            intent(in)  :: p_idx  !< Peer slot (the realm owning the ghosts).
   real(R8P),               intent(out) :: buf(:) !< Packed values.

   associate(maps=>self%adam%maps)
   if (self%stage_active > 0_I4P) then
      call pack_seam_rows(rows=maps%seam_mpi_send_cell, row_start=maps%seam_mpi_send_row_start(p_idx),            &
                          row_count=maps%seam_mpi_send_row_count(p_idx), regime=maps%seam_ghost_fill, ngc=self%ngc, &
                          q=self%rk%q_rk(:,:,:,:,:,self%stage_active), buf=buf)
   else
      call pack_seam_rows(rows=maps%seam_mpi_send_cell, row_start=maps%seam_mpi_send_row_start(p_idx),            &
                          row_count=maps%seam_mpi_send_row_count(p_idx), regime=maps%seam_ghost_fill, ngc=self%ngc, &
                          q=self%q, buf=buf)
   endif
   endassociate
   endsubroutine pack_seam_cells_forest

   subroutine unpack_seam_cells_forest(self, p_idx, buf)
   !< Unpack `buf` into this realm's seam ghosts of the cross-rank receive rows of peer slot `p_idx`, on the active
   !< buffer (issue #40; `adam_seam_exchange%unpack_seam_rows`).
   class(flume_cpu_object), intent(inout) :: self   !< The realm.
   integer(I4P),            intent(in)    :: p_idx  !< Peer slot (the realm owning the cells).
   real(R8P),               intent(in)    :: buf(:) !< Packed values.

   associate(maps=>self%adam%maps)
   if (self%stage_active > 0_I4P) then
      call unpack_seam_rows(rows=maps%seam_mpi_recv_cell, row_start=maps%seam_mpi_recv_row_start(p_idx),          &
                            row_count=maps%seam_mpi_recv_row_count(p_idx), ngc=self%ngc, buf=buf,                   &
                            q=self%rk%q_rk(:,:,:,:,:,self%stage_active))
   else
      call unpack_seam_rows(rows=maps%seam_mpi_recv_cell, row_start=maps%seam_mpi_recv_row_start(p_idx),          &
                            row_count=maps%seam_mpi_recv_row_count(p_idx), ngc=self%ngc, buf=buf, q=self%q)
   endif
   endassociate
   endsubroutine unpack_seam_cells_forest

   subroutine finalize_forest(self)
   !< Finalize the realm: close the output files and free the data (MPI is finalized once by the forest).
   class(flume_cpu_object), intent(inout) :: self !< The equation.

   call self%io%close_file_residuals
   call self%diagnostics%close_file
   call self%destroy_common
   if (allocated(self%flx_f)) deallocate(self%flx_f)
   if (allocated(self%fly_f)) deallocate(self%fly_f)
   if (allocated(self%flz_f)) deallocate(self%flz_f)
   if (allocated(self%lam)) deallocate(self%lam)
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
      call self%compute_phi
   else
      do i=1, self%ic%amr_iterations
         call self%set_initial_conditions
         call self%compute_phi
         call self%amr_update
      enddo
      call self%set_initial_conditions
      call self%compute_phi
      call self%adam%make_comm_local_maps_ghost_bc
      self%time%time = 0._R8P
      self%time%it   = 0_I4P
   endif
   call self%set_glm_damping
   call self%update_ghost(q=self%q)
   call self%compute_q_aux(q=self%q)
   call self%check_glm_ch
   call self%diagnostics%open_file(output_basename=self%io%output_basename, q_name=self%q_name, &
                                   is_restart=self%io%restart, with_divb=self%physics%model /= MODEL_EULER)
   ! a restarted run starts from a step its predecessor already saved: saving it again would duplicate the history rows
   if (.not.self%io%restart) then
      call self%compute_divb_history
      call self%save_simulation_data
   endif
   call self%io%open_file_residuals(nv=self%nv, is_restart=self%io%restart)
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
   !< Post-step work: the non-finite state check, the GLM c_h check, the div(B) history (MHD), fields, restart and
   !< conservation history on their cadence.
   class(flume_cpu_object), intent(inout)                   :: self              !< The equation.
   real(R8P),               intent(in)                      :: dt                !< Time step just advanced.
   real(R8P),               intent(in)                      :: t                 !< Time after the advance.
   integer(I4P),            intent(in)                      :: it                !< Iteration after the advance.
   logical,                 intent(in),    optional         :: do_save_state     !< Unused: cadence is internal.
   logical,                 intent(in),    optional         :: do_save_residuals !< Unused: cadence is internal.
   logical,                 intent(in),    optional         :: do_save_restart   !< Unused: cadence is internal.
   logical,                 intent(in),    optional         :: do_amr            !< Unused: AMR is init-time only.
   class(realm_object),     intent(inout), optional, target :: realm(:)          !< Sibling realms.

   call self%check_nonfinite
   call self%check_glm_ch
   call self%compute_divb_history(realm=realm)
   call self%save_simulation_data(realm=realm)
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

   subroutine compute_flux_difference(nv, ni, nj, nk, ngc, blocks_number, is_null, freeze, dxyz, flx, fly, flz, dq, phi)
   !< Compute the residuals from the face fluxes, `dq = -sum_d (F_{d,i+1/2} - F_{d,i-1/2}) / dx_d`.
   !<
   !< A null direction weighs zero, and the residual of the variable it freezes, `freeze(d)` (0: none), is zero: the
   !< normal momentum for Euler (CHASE semantics, issue #35, section 3.4), none for MHD (issue #41, M2-P3).
   !< With immersed solids (`phi` present, its last slot the all-solids summary), the spacing of a fluid cell is cut by
   !< the solid surface (`ib_cut_spacing`, CHASE semantics, D-9).
   integer(I4P), intent(in)           :: nv                              !< Conservative variables number.
   integer(I4P), intent(in)           :: ni, nj, nk, ngc                 !< Grid dimensions.
   integer(I4P), intent(in)           :: blocks_number                   !< Actual blocks number.
   logical,      intent(in)           :: is_null(3)                      !< Null directions.
   integer(I4P), intent(in)           :: freeze(3)                       !< Variable frozen by each null direction.
   real(R8P),    intent(in)           :: dxyz(1:,1:)                     !< Blocks space steps [3, nb].
   real(R8P),    intent(in)           :: flx(1:,0:,1:,1:,1:)             !< X-face fluxes.
   real(R8P),    intent(in)           :: fly(1:,1:,0:,1:,1:)             !< Y-face fluxes.
   real(R8P),    intent(in)           :: flz(1:,1:,1:,0:,1:)             !< Z-face fluxes.
   real(R8P),    intent(inout)        :: dq(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Residuals.
   real(R8P),    intent(in), optional :: phi(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Immersed solids distance function.
   real(R8P)                          :: wx, wy, wz                      !< Direction weights: 1 active, 0 null.
   real(R8P)                          :: dx, dy, dz                      !< Cell spacings.
   integer(I4P)                       :: ns                              !< All-solids summary slot of phi.
   integer(I4P)                       :: b, i, j, k, v                   !< Counters.

   wx = merge(0._R8P, 1._R8P, is_null(1))
   wy = merge(0._R8P, 1._R8P, is_null(2))
   wz = merge(0._R8P, 1._R8P, is_null(3))
   ns = 0_I4P ; if (present(phi)) ns = ubound(phi, dim=1)
   !$omp parallel do collapse(4) default(firstprivate) shared(dxyz, flx, fly, flz, dq, phi)
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               dx = dxyz(1,b) ; dy = dxyz(2,b) ; dz = dxyz(3,b)
               if (ns > 0_I4P) then
                  dx = ib_cut_spacing(phi_c=phi(ns,i,j,k,b), phi_m=phi(ns,i-1,j,k,b), phi_p=phi(ns,i+1,j,k,b), ds=dx)
                  dy = ib_cut_spacing(phi_c=phi(ns,i,j,k,b), phi_m=phi(ns,i,j-1,k,b), phi_p=phi(ns,i,j+1,k,b), ds=dy)
                  dz = ib_cut_spacing(phi_c=phi(ns,i,j,k,b), phi_m=phi(ns,i,j,k-1,b), phi_p=phi(ns,i,j,k+1,b), ds=dz)
               endif
               do v=1, nv
                  dq(v,i,j,k,b) = -(wx * (flx(v,i,j,k,b) - flx(v,i-1,j,k,b)) / dx + &
                                    wy * (fly(v,i,j,k,b) - fly(v,i,j-1,k,b)) / dy + &
                                    wz * (flz(v,i,j,k,b) - flz(v,i,j,k-1,b)) / dz)
               enddo
               if (freeze(1) > 0_I4P) dq(freeze(1),i,j,k,b) = 0._R8P
               if (freeze(2) > 0_I4P) dq(freeze(2),i,j,k,b) = 0._R8P
               if (freeze(3) > 0_I4P) dq(freeze(3),i,j,k,b) = 0._R8P
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine compute_flux_difference

   subroutine poison_ghosts(ni, nj, nk, ngc, blocks_number, q)
   !< Set every ghost cell to a quiet NaN (verification instrument, `[diagnostics] ghost_poison`, issue #65).
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                 !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                   !< Actual blocks number.
   real(R8P),    intent(inout) :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Conservative variables.
   real(R8P)                   :: nan                             !< Quiet NaN.
   integer(I4P)                :: b, i, j, k                      !< Counters.

   nan = ieee_value(nan, ieee_quiet_nan)
   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               if (i < 1 .or. i > ni .or. j < 1 .or. j > nj .or. k < 1 .or. k > nk) q(:,i,j,k,b) = nan
            enddo
         enddo
      enddo
   enddo
   endsubroutine poison_ghosts

   subroutine compute_rk_ssp_residual(self)
   !< Compute the residual of a strong stability preserving step, `dq = sum_s beta_s dq_s`, from the stored stages.
   !<
   !< Used with immersed solids only: the library `update_q` computes the step residual on its unmasked path but not on
   !< the masked one. Same summation order as the FNL kernel `compute_rk_ssp_residual_dev`.
   class(flume_cpu_object), intent(inout) :: self          !< The equation.
   integer(I4P)                           :: b, i, j, k, s !< Counters.

   self%dq = 0._R8P
   !$omp parallel do collapse(4) default(firstprivate) shared(self)
   do b=1, self%blocks_number
      do k=1, self%nk
         do j=1, self%nj
            do i=1, self%ni
               do s=1, self%rk%nrk
                  self%dq(:,i,j,k,b) = self%dq(:,i,j,k,b) + self%rk%beta(s) * self%rk%q_rk(:,i,j,k,b,s)
               enddo
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine compute_rk_ssp_residual

   subroutine compute_residuals_weno(self, q, dq, s, flux_register)
   !< Compute the residuals with the WENO space operator: ghost update, auxiliary variables, face fluxes of the active
   !< directions, flux difference, and with GLM the psi damping source.
   !<
   !< The face fluxes of a null direction are never computed: they keep their zero initialization. On the staged path
   !< (AMR seam faces), the seam face fluxes of every stage are accumulated into the forest's flux register. With
   !< immersed solids, the eikonal extrapolation fills the solid cells first (`n_eikonal` Jacobi sweeps, each followed by
   !< a ghost exchange, then the wall inversion), and the flux difference uses the spacing cut by the surface.
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
   logical                                              :: is_char       !< Characteristic reconstruction flag.
   integer(I4P)                                         :: e             !< Eikonal iterations counter.

   if (self%ib%solids_number > 0_I4P) then
      call self%update_ghost(q=q)
      do e=1, self%ib%n_eikonal
         call self%ib%evolve_eikonal(field=self%adam%field, grid=self%adam%grid, q=q, dq=dq)
         call self%update_ghost(q=q)
      enddo
      call self%ib%invert_eikonal(field=self%adam%field, grid=self%adam%grid, q=q)
   endif
   call self%apply_floors(q=q)
   call self%update_ghost(q=q)
   call self%compute_q_aux(q=q)
   is_char = self%numerics%reconstruction_variables == RECON_CHARACTERISTIC
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nb=>self%blocks_number, gamma=>self%physics%gamma, &
             is_null=>self%adam%grid%null_xyz)
   select case(self%physics%model)
   case(MODEL_EULER)
      if (.not.is_null(1)) call compute_face_fluxes_euler(d=1_I4P, di=1_I4P, dj=0_I4P, dk=0_I4P, ni=ni,          &
                                                          nj=nj, nk=nk, ngc=ngc, blocks_number=nb, gamma=gamma,  &
                                                          ch=self%physics%mhd%glm_ch, is_characteristic=is_char, &
                                                          weno=self%weno, q=q, q_aux=self%q_aux, fl=self%flx_f)
      if (.not.is_null(2)) call compute_face_fluxes_euler(d=2_I4P, di=0_I4P, dj=1_I4P, dk=0_I4P, ni=ni,          &
                                                          nj=nj, nk=nk, ngc=ngc, blocks_number=nb, gamma=gamma,  &
                                                          ch=self%physics%mhd%glm_ch, is_characteristic=is_char, &
                                                          weno=self%weno, q=q, q_aux=self%q_aux, fl=self%fly_f)
      if (.not.is_null(3)) call compute_face_fluxes_euler(d=3_I4P, di=0_I4P, dj=0_I4P, dk=1_I4P, ni=ni,          &
                                                          nj=nj, nk=nk, ngc=ngc, blocks_number=nb, gamma=gamma,  &
                                                          ch=self%physics%mhd%glm_ch, is_characteristic=is_char, &
                                                          weno=self%weno, q=q, q_aux=self%q_aux, fl=self%flz_f)
   case(MODEL_MHD)
      if (.not.is_null(1)) call compute_face_fluxes_mhd(d=1_I4P, di=1_I4P, dj=0_I4P, dk=0_I4P, ni=ni,          &
                                                        nj=nj, nk=nk, ngc=ngc, blocks_number=nb, gamma=gamma,  &
                                                        ch=self%physics%mhd%glm_ch, is_characteristic=is_char, &
                                                        weno=self%weno, q=q, q_aux=self%q_aux, fl=self%flx_f)
      if (.not.is_null(2)) call compute_face_fluxes_mhd(d=2_I4P, di=0_I4P, dj=1_I4P, dk=0_I4P, ni=ni,          &
                                                        nj=nj, nk=nk, ngc=ngc, blocks_number=nb, gamma=gamma,  &
                                                        ch=self%physics%mhd%glm_ch, is_characteristic=is_char, &
                                                        weno=self%weno, q=q, q_aux=self%q_aux, fl=self%fly_f)
      if (.not.is_null(3)) call compute_face_fluxes_mhd(d=3_I4P, di=0_I4P, dj=0_I4P, dk=1_I4P, ni=ni,          &
                                                        nj=nj, nk=nk, ngc=ngc, blocks_number=nb, gamma=gamma,  &
                                                        ch=self%physics%mhd%glm_ch, is_characteristic=is_char, &
                                                        weno=self%weno, q=q, q_aux=self%q_aux, fl=self%flz_f)
   case(MODEL_MHD_GLM)
      if (.not.is_null(1)) call compute_face_fluxes_mhd_glm(d=1_I4P, di=1_I4P, dj=0_I4P, dk=0_I4P, ni=ni,          &
                                                            nj=nj, nk=nk, ngc=ngc, blocks_number=nb, gamma=gamma,  &
                                                            ch=self%physics%mhd%glm_ch, is_characteristic=is_char, &
                                                            weno=self%weno, q=q, q_aux=self%q_aux, fl=self%flx_f)
      if (.not.is_null(2)) call compute_face_fluxes_mhd_glm(d=2_I4P, di=0_I4P, dj=1_I4P, dk=0_I4P, ni=ni,          &
                                                            nj=nj, nk=nk, ngc=ngc, blocks_number=nb, gamma=gamma,  &
                                                            ch=self%physics%mhd%glm_ch, is_characteristic=is_char, &
                                                            weno=self%weno, q=q, q_aux=self%q_aux, fl=self%fly_f)
      if (.not.is_null(3)) call compute_face_fluxes_mhd_glm(d=3_I4P, di=0_I4P, dj=0_I4P, dk=1_I4P, ni=ni,          &
                                                            nj=nj, nk=nk, ngc=ngc, blocks_number=nb, gamma=gamma,  &
                                                            ch=self%physics%mhd%glm_ch, is_characteristic=is_char, &
                                                            weno=self%weno, q=q, q_aux=self%q_aux, fl=self%flz_f)
   case(MODEL_MHD_EGLM)
      if (.not.is_null(1)) call compute_face_fluxes_mhd_eglm(d=1_I4P, di=1_I4P, dj=0_I4P, dk=0_I4P, ni=ni,          &
                                                            nj=nj, nk=nk, ngc=ngc, blocks_number=nb, gamma=gamma,  &
                                                            ch=self%physics%mhd%glm_ch, is_characteristic=is_char, &
                                                            weno=self%weno, q=q, q_aux=self%q_aux, fl=self%flx_f)
      if (.not.is_null(2)) call compute_face_fluxes_mhd_eglm(d=2_I4P, di=0_I4P, dj=1_I4P, dk=0_I4P, ni=ni,          &
                                                            nj=nj, nk=nk, ngc=ngc, blocks_number=nb, gamma=gamma,  &
                                                            ch=self%physics%mhd%glm_ch, is_characteristic=is_char, &
                                                            weno=self%weno, q=q, q_aux=self%q_aux, fl=self%fly_f)
      if (.not.is_null(3)) call compute_face_fluxes_mhd_eglm(d=3_I4P, di=0_I4P, dj=0_I4P, dk=1_I4P, ni=ni,          &
                                                            nj=nj, nk=nk, ngc=ngc, blocks_number=nb, gamma=gamma,  &
                                                            ch=self%physics%mhd%glm_ch, is_characteristic=is_char, &
                                                            weno=self%weno, q=q, q_aux=self%q_aux, fl=self%flz_f)
   case default
      call mpih%error_stop(msg=': no CPU kernels for physical model "'//self%physics%physical_model//'"')
   endselect
   if (self%physics%dissipation%is_active) call self%add_dissipative_fluxes
   if (self%numerics%positivity_limiter == POSITIVITY_LIMITER_CELL) call self%limit_positivity(q=q, &
                                                                                         flux_register=flux_register)
   if (present(flux_register) .and. present(s) .and. self%numerics%reflux) then
      if (flux_register%nfaces > 0_I4P) call self%accumulate_seam_fluxes(s=s, flux_register=flux_register)
   endif
   call compute_flux_difference(nv=self%physics%nv, ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, is_null=is_null,  &
                                freeze=self%null_freeze(),                                                           &
                                dxyz=self%adam%field%dxyz, flx=self%flx_f, fly=self%fly_f, flz=self%flz_f, dq=dq, &
                                phi=self%ib%phi)
   if (self%physics%model == MODEL_MHD_GLM) call add_glm_damping(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb,          &
                                                                  damping=self%physics%mhd%glm_damping, q=q, dq=dq)
   if (self%physics%model == MODEL_MHD_EGLM) then
      call add_glm_damping_eglm(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, damping=self%physics%mhd%glm_damping, &
                                q=q, dq=dq)
      if (self%numerics%positivity_limiter == POSITIVITY_LIMITER_CELL) then
         call add_eglm_sources_limited(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, hs=self%weno%S,              &
                                       dxyz=self%adam%field%dxyz, is_null=is_null, q=q, q_aux=self%q_aux,         &
                                       lam=self%lam, dq=dq)
      else
         call add_eglm_sources(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, hs=self%weno%S,                      &
                               dxyz=self%adam%field%dxyz, is_null=is_null, q=q, q_aux=self%q_aux, dq=dq)
      endif
   endif
   endassociate
   endsubroutine compute_residuals_weno

   subroutine compute_residuals_riemann(self, q, dq, s, flux_register)
   !< Compute the residuals with the `weno-riemann` space operator (issue #47): as `compute_residuals_weno`, with the
   !< face fluxes of the WENO interpolation, the Riemann solver and the high-order correction.
   !<
   !< The host selects the kernel of the (model, solver) pair; a pair without kernels is fatal. The pre-steps (immersed
   !< solids, floors, ghosts, auxiliary variables) and the post-steps (seam accumulation, flux difference, GLM damping,
   !< EGLM damping and sources) are those of `compute_residuals_weno`.
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
   logical                                              :: is_char       !< Characteristic interpolation flag.
   real(R8P)                                            :: cc(3)         !< Correction coefficients.
   real(R8P)                                            :: tau           !< Correction sensor threshold.
   integer(I4P)                                         :: e             !< Eikonal iterations counter.
   procedure(compute_riemann_face_fluxes_euler_llf), pointer :: face_fluxes !< Kernel of the (model, solver) pair.
   integer(I4P)                                         :: fallbacks(3)  !< Riemann solver fallbacks per direction.

   if (self%ib%solids_number > 0_I4P) then
      call self%update_ghost(q=q)
      do e=1, self%ib%n_eikonal
         call self%ib%evolve_eikonal(field=self%adam%field, grid=self%adam%grid, q=q, dq=dq)
         call self%update_ghost(q=q)
      enddo
      call self%ib%invert_eikonal(field=self%adam%field, grid=self%adam%grid, q=q)
   endif
   call self%apply_floors(q=q)
   call self%update_ghost(q=q)
   call self%compute_q_aux(q=q)
   is_char = self%numerics%reconstruction_variables == RECON_CHARACTERISTIC
   cc = self%numerics%correction_coefficients()
   tau = self%numerics%correction_threshold()
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nb=>self%blocks_number, gamma=>self%physics%gamma, &
             ch=>self%physics%mhd%glm_ch, is_null=>self%adam%grid%null_xyz)
   face_fluxes => null()
   if (self%physics%model == MODEL_EULER) then
      select case(self%numerics%riemann_solver)
      case(RIEMANN_SOLVER_LLF)
         face_fluxes => compute_riemann_face_fluxes_euler_llf
      case(RIEMANN_SOLVER_HLL)
         face_fluxes => compute_riemann_face_fluxes_euler_hll
      case(RIEMANN_SOLVER_HLLC)
         face_fluxes => compute_riemann_face_fluxes_euler_hllc
      endselect
   elseif (self%physics%model == MODEL_MHD) then
      select case(self%numerics%riemann_solver)
      case(RIEMANN_SOLVER_LLF)
         face_fluxes => compute_riemann_face_fluxes_mhd_llf
      case(RIEMANN_SOLVER_HLL)
         face_fluxes => compute_riemann_face_fluxes_mhd_hll
      case(RIEMANN_SOLVER_HLLD)
         face_fluxes => compute_riemann_face_fluxes_mhd_hlld
      endselect
   elseif (self%physics%model == MODEL_MHD_GLM) then
      select case(self%numerics%riemann_solver)
      case(RIEMANN_SOLVER_LLF)
         face_fluxes => compute_riemann_face_fluxes_mhd_glm_llf
      case(RIEMANN_SOLVER_HLL)
         face_fluxes => compute_riemann_face_fluxes_mhd_glm_hll
      case(RIEMANN_SOLVER_HLLD)
         face_fluxes => compute_riemann_face_fluxes_mhd_glm_hlld
      endselect
   elseif (self%physics%model == MODEL_MHD_EGLM) then
      select case(self%numerics%riemann_solver)
      case(RIEMANN_SOLVER_LLF)
         face_fluxes => compute_riemann_face_fluxes_mhd_eglm_llf
      case(RIEMANN_SOLVER_HLL)
         face_fluxes => compute_riemann_face_fluxes_mhd_eglm_hll
      case(RIEMANN_SOLVER_HLLD)
         face_fluxes => compute_riemann_face_fluxes_mhd_eglm_hlld
      endselect
   endif
   if (.not.associated(face_fluxes)) call mpih%error_stop(msg=': no CPU weno-riemann kernels yet for '//          &
                                                              '[physics].(physical_model)='//                     &
                                                              self%physics%physical_model//                       &
                                                              ' with [numerics].(riemann_solver)='//              &
                                                              self%numerics%riemann_solver)
   fallbacks = 0_I4P
   if (.not.is_null(1)) call face_fluxes(d=1_I4P, di=1_I4P, dj=0_I4P, dk=0_I4P, ni=ni, nj=nj, nk=nk, ngc=ngc,         &
                                         blocks_number=nb, gamma=gamma, ch=ch, is_characteristic=is_char,             &
                                         weno=self%weno, cc=cc, tau=tau, q=q, q_aux=self%q_aux, fl=self%flx_f,   &
                                         fallbacks=fallbacks(1))
   if (.not.is_null(2)) call face_fluxes(d=2_I4P, di=0_I4P, dj=1_I4P, dk=0_I4P, ni=ni, nj=nj, nk=nk, ngc=ngc,         &
                                         blocks_number=nb, gamma=gamma, ch=ch, is_characteristic=is_char,             &
                                         weno=self%weno, cc=cc, tau=tau, q=q, q_aux=self%q_aux, fl=self%fly_f,   &
                                         fallbacks=fallbacks(2))
   if (.not.is_null(3)) call face_fluxes(d=3_I4P, di=0_I4P, dj=0_I4P, dk=1_I4P, ni=ni, nj=nj, nk=nk, ngc=ngc,         &
                                         blocks_number=nb, gamma=gamma, ch=ch, is_characteristic=is_char,             &
                                         weno=self%weno, cc=cc, tau=tau, q=q, q_aux=self%q_aux, fl=self%flz_f,   &
                                         fallbacks=fallbacks(3))
   if (self%numerics%riemann_solver == RIEMANN_SOLVER_HLLD) then
      call MPI_ALLREDUCE(MPI_IN_PLACE, fallbacks, 3, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpih%error)
      if (sum(fallbacks) > 0_I4P .and. mpih%myrank == 0) &
         print '(A)', mpih%myrankstr//'HLLD fallbacks to HLL: '//trim(str(sum(fallbacks)))//' faces at step '// &
                      trim(str(self%time%it))
   endif
   if (self%physics%dissipation%is_active) call self%add_dissipative_fluxes
   if (self%numerics%positivity_limiter == POSITIVITY_LIMITER_CELL) call self%limit_positivity(q=q, &
                                                                                         flux_register=flux_register)
   if (present(flux_register) .and. present(s) .and. self%numerics%reflux) then
      if (flux_register%nfaces > 0_I4P) call self%accumulate_seam_fluxes(s=s, flux_register=flux_register)
   endif
   call compute_flux_difference(nv=self%physics%nv, ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, is_null=is_null,  &
                                freeze=self%null_freeze(),                                                           &
                                dxyz=self%adam%field%dxyz, flx=self%flx_f, fly=self%fly_f, flz=self%flz_f, dq=dq, &
                                phi=self%ib%phi)
   if (self%physics%model == MODEL_MHD_GLM) call add_glm_damping(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb,          &
                                                                  damping=self%physics%mhd%glm_damping, q=q, dq=dq)
   if (self%physics%model == MODEL_MHD_EGLM) then
      call add_glm_damping_eglm(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, damping=self%physics%mhd%glm_damping, &
                                q=q, dq=dq)
      if (self%numerics%positivity_limiter == POSITIVITY_LIMITER_CELL) then
         call add_eglm_sources_limited(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, hs=self%weno%S,              &
                                       dxyz=self%adam%field%dxyz, is_null=is_null, q=q, q_aux=self%q_aux,         &
                                       lam=self%lam, dq=dq)
      else
         call add_eglm_sources(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, hs=self%weno%S,                      &
                               dxyz=self%adam%field%dxyz, is_null=is_null, q=q, q_aux=self%q_aux, dq=dq)
      endif
   endif
   endassociate
   endsubroutine compute_residuals_riemann

   subroutine integrate_rk_ls(self)
   !< Integrate one time step with a low-storage Runge-Kutta scheme.
   class(flume_cpu_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: s    !< Counter.

   call self%rk%initialize_stages(field=self%adam%field, q=self%q)
   do s=1, self%rk%nrk
      call self%compute_residuals(q=self%q, dq=self%dq)
      if (s == 1) call self%save_residuals
      call self%rk%compute_stage_ls(field=self%adam%field, s=s, dt=self%time%dt, phi=self%ib%phi, dq=self%dq, q=self%q)
   enddo
   endsubroutine integrate_rk_ls

   subroutine integrate_rk_ssp(self)
   !< Integrate one time step with a strong stability preserving Runge-Kutta scheme.
   class(flume_cpu_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: s    !< Counter.

   call self%rk%initialize_stages(field=self%adam%field, q=self%q)
   do s=1, self%rk%nrk
      call self%rk%compute_stage(field=self%adam%field, s=s, dt=self%time%dt, phi=self%ib%phi)
      call self%compute_residuals(q=self%rk%q_rk(:,:,:,:,:,s), dq=self%dq, s=s)
      call self%rk%assign_stage(field=self%adam%field, s=s, q=self%dq, phi=self%ib%phi)
   enddo
   call self%rk%update_q(field=self%adam%field, dt=self%time%dt, phi=self%ib%phi, q=self%q, dq=self%dq)
   if (allocated(self%ib%phi)) call compute_rk_ssp_residual(self)
   call self%save_residuals
   endsubroutine integrate_rk_ssp
endmodule adam_flume_cpu_object
