!< ADAM, FLUME FNL (OpenACC / OpenMP offload) backend object.

#include "fundal.H"

module adam_flume_fnl_object
!< ADAM, FLUME FNL (OpenACC / OpenMP offload) backend object.
!<
!< Implements the forest contract on the device. The state lives on the device for the whole run in the transposed
!< layout `(b, i, j, k, v)`; full-field device-to-host copies happen only on save steps. The FNL helpers (field, IB,
!< RK, WENO) are per-realm components initialized after the common initialization from the realm's own objects;
!< `mpih_fnl` is the only FNL singleton and is initialized once per process.
!< Space operator: characteristic (or conservative) WENO flux splitting, face fluxes then flux difference (issue #35,
!< section 3.4); the per-face physics is the shared `adam_flume_euler_library`, so only the kernels are FNL-specific.

! ADAM classes, libraries, parameters
use :: adam_flux_register_object, only : face_tangential_ratios, flux_register_object
use :: adam_maps_object,          only : face_axis_sign
use :: adam_parameters,           only : BC_SEAM
use :: adam_realm_object,         only : realm_object
use :: adam_seam_exchange,        only : seam_fill_all
use :: adam_rk_object,            only : RK_1, RK_2, RK_3, RK_SSP_11, RK_SSP_22, RK_SSP_33, RK_SSP_54
! ADAM FNL classes, libraries
use :: adam_fnl_field_kernels,    only : compute_normL2_residuals_dev, pack_seam_rows_dev, unpack_seam_rows_dev
use :: adam_fnl_field_object,     only : field_fnl_object
use :: adam_fnl_ib_object,        only : ib_fnl_object
use :: adam_fnl_rk_object,        only : rk_fnl_object
use :: adam_fnl_weno_object,      only : weno_fnl_object
! ADAM singleton objects
use :: adam_fnl_mpih_global,      only : mpih_fnl, mpih_fnl_is_initialized
! FLUME modules
use :: adam_flume_common_library,      only : flume_common_object, flume_seam_sync_object, seam_face_cells,             &
                                              seam_fine_to_coarse, seam_skin_cell, seam_skin_index,                    &
                                              MODEL_EULER, MODEL_MHD, MODEL_MHD_EGLM,                                  &
                                              MODEL_MHD_GLM, POSITIVITY_LIMITER_CELL,                                &
                                              RECON_CHARACTERISTIC, RIEMANN_SOLVER_HLL, RIEMANN_SOLVER_HLLC,            &
                                              RIEMANN_SOLVER_HLLD, RIEMANN_SOLVER_LLF, SCHEME_SPACE_WENO,              &
                                              SCHEME_SPACE_WENO_RIEMANN
use :: adam_flume_fnl_euler_hll_kernels,  only : compute_riemann_face_fluxes_euler_hll_dev=>compute_riemann_face_fluxes_dev
use :: adam_flume_fnl_euler_hllc_kernels, only : compute_riemann_face_fluxes_euler_hllc_dev=>compute_riemann_face_fluxes_dev
use :: adam_flume_fnl_euler_llf_kernels,  only : compute_riemann_face_fluxes_euler_llf_dev=>compute_riemann_face_fluxes_dev
use :: adam_flume_fnl_mhd_hll_kernels,          only : compute_riemann_face_fluxes_mhd_hll_dev=>compute_riemann_face_fluxes_dev
use :: adam_flume_fnl_mhd_hlld_kernels,         only : compute_riemann_face_fluxes_mhd_hlld_dev=>compute_riemann_face_fluxes_dev
use :: adam_flume_fnl_mhd_llf_kernels,          only : compute_riemann_face_fluxes_mhd_llf_dev=>compute_riemann_face_fluxes_dev
use :: adam_flume_fnl_mhd_glm_hll_kernels,      only : compute_riemann_face_fluxes_mhd_glm_hll_dev=>compute_riemann_face_fluxes_dev
use :: adam_flume_fnl_mhd_glm_hlld_kernels,     only : compute_riemann_face_fluxes_mhd_glm_hlld_dev=>compute_riemann_face_fluxes_dev
use :: adam_flume_fnl_mhd_glm_llf_kernels,      only : compute_riemann_face_fluxes_mhd_glm_llf_dev=>compute_riemann_face_fluxes_dev
use :: adam_flume_fnl_mhd_eglm_hll_kernels,     only : compute_riemann_face_fluxes_mhd_eglm_hll_dev=>compute_riemann_face_fluxes_dev
use :: adam_flume_fnl_mhd_eglm_hlld_kernels,    only : &
                                                   compute_riemann_face_fluxes_mhd_eglm_hlld_dev=>compute_riemann_face_fluxes_dev
use :: adam_flume_fnl_mhd_eglm_llf_kernels,     only : compute_riemann_face_fluxes_mhd_eglm_llf_dev=>compute_riemann_face_fluxes_dev
use :: adam_flume_fnl_euler_kernels,   only : blend_inadmissible_ghosts_euler_dev=>blend_inadmissible_ghosts_dev,      &
                                              blend_positivity_fluxes_euler_dev=>blend_positivity_fluxes_dev,          &
                                              compute_backbone_fluxes_euler_host=>compute_backbone_fluxes_host,        &
                                              compute_seam_positivity_factors_euler_host=>                             &
                                                 compute_seam_positivity_factors_host,                                &
                                              compute_conservation_euler_dev=>compute_conservation_dev,                &
                                              compute_positivity_factors_euler_dev=>compute_positivity_factors_dev,    &
                                              compute_face_fluxes_euler_dev=>compute_face_fluxes_dev,                  &
                                              compute_lambda_max_euler_dev=>compute_lambda_max_dev,                    &
                                              compute_q_aux_euler_dev=>compute_q_aux_dev,                              &
                                              count_nonfinite_euler_dev=>count_nonfinite_dev
use :: adam_flume_fnl_mhd_kernels,     only : apply_floors_mhd_dev=>apply_floors_dev,                           &
                                              blend_inadmissible_ghosts_mhd_dev=>blend_inadmissible_ghosts_dev,       &
                                              blend_positivity_fluxes_mhd_dev=>blend_positivity_fluxes_dev,           &
                                              compute_backbone_fluxes_mhd_host=>compute_backbone_fluxes_host,         &
                                              compute_seam_positivity_factors_mhd_host=>                              &
                                                 compute_seam_positivity_factors_host,                               &
                                              compute_positivity_factors_mhd_dev=>compute_positivity_factors_dev,     &
                                              compute_divb_norms_mhd_dev=>compute_divb_norms_dev,                     &
                                              compute_conservation_mhd_dev=>compute_conservation_dev,                 &
                                              compute_face_fluxes_mhd_dev=>compute_face_fluxes_dev,                   &
                                              compute_lambda_max_mhd_dev=>compute_lambda_max_dev,                     &
                                              compute_q_aux_mhd_dev=>compute_q_aux_dev,                               &
                                              count_nonfinite_mhd_dev=>count_nonfinite_dev
use :: adam_flume_fnl_mhd_eglm_kernels, only : add_eglm_sources_dev, add_eglm_sources_limited_dev,                  &
                                               add_glm_damping_eglm_dev=>add_glm_damping_dev,                         &
                                               blend_inadmissible_ghosts_mhd_eglm_dev=>blend_inadmissible_ghosts_dev, &
                                               blend_positivity_fluxes_mhd_eglm_dev=>blend_positivity_fluxes_dev,     &
                                               compute_backbone_fluxes_mhd_eglm_host=>compute_backbone_fluxes_host,   &
                                               compute_seam_positivity_factors_mhd_eglm_host=>                        &
                                                  compute_seam_positivity_factors_host,                              &
                                               compute_positivity_factors_mhd_eglm_dev=>                              &
                                               compute_positivity_factors_dev,                                        &
                                               apply_floors_mhd_eglm_dev=>apply_floors_dev,                           &
                                               compute_divb_norms_mhd_eglm_dev=>compute_divb_norms_dev,               &
                                               compute_conservation_mhd_eglm_dev=>compute_conservation_dev,           &
                                               compute_face_fluxes_mhd_eglm_dev=>compute_face_fluxes_dev,             &
                                               compute_lambda_max_mhd_eglm_dev=>compute_lambda_max_dev,               &
                                               compute_q_aux_mhd_eglm_dev=>compute_q_aux_dev,                         &
                                               compute_speed_max_mhd_eglm_dev=>compute_speed_max_dev,                 &
                                               count_nonfinite_mhd_eglm_dev=>count_nonfinite_dev
use :: adam_flume_fnl_mhd_glm_kernels, only : add_glm_damping_dev, apply_floors_mhd_glm_dev=>apply_floors_dev,  &
                                              blend_inadmissible_ghosts_mhd_glm_dev=>blend_inadmissible_ghosts_dev,   &
                                              compute_divb_norms_mhd_glm_dev=>compute_divb_norms_dev,                 &
                                              compute_conservation_mhd_glm_dev=>compute_conservation_dev,             &
                                              compute_face_fluxes_mhd_glm_dev=>compute_face_fluxes_dev,               &
                                              compute_lambda_max_mhd_glm_dev=>compute_lambda_max_dev,                 &
                                              compute_q_aux_mhd_glm_dev=>compute_q_aux_dev,                           &
                                              compute_speed_max_mhd_glm_dev=>compute_speed_max_dev,                   &
                                              count_nonfinite_mhd_glm_dev=>count_nonfinite_dev
use :: adam_flume_fnl_kernels,         only : apply_reflux_face_dev,                                                   &
                                              compute_flux_difference_dev, compute_flux_difference_ib_dev,             &
                                              compute_rk_ssp_residual_dev, fill_seam_copy_dev, gather_seam_cells_dev,  &
                                              gather_seam_faces_dev, gather_seam_stencils_dev, pack_seam_skin_dev,     &
                                              scatter_seam_cells_dev, scatter_seam_faces_dev,                          &
                                              set_boundary_conditions_dev
! third party modules
use :: fundal,                    only : dev_alloc, dev_assign_to_device, dev_free, dev_memcpy_from_device,     &
                                         dev_memcpy_to_device, mydev
use :: mpi
use :: penf,                      only : I4P, I8P, R8P, str

implicit none
private
public :: flume_fnl_object

type, extends(flume_common_object) :: flume_fnl_object
   !< FLUME FNL backend object.
   ! FNL helpers
   type(field_fnl_object) :: field_fnl                    !< Field helper (coordinates, maps, ghost exchange).
   type(ib_fnl_object)    :: ib_fnl                       !< Immersed boundary helper.
   type(rk_fnl_object)    :: rk_fnl                       !< Runge-Kutta helper.
   type(weno_fnl_object)  :: weno_fnl                     !< WENO helper.
   ! device data
   real(R8P), pointer     :: q_gpu(:,:,:,:,:)=>null()     !< Conservative variables [nb, i, j, k, nv].
   real(R8P), pointer     :: dq_gpu(:,:,:,:,:)=>null()    !< Residuals [nb, i, j, k, nv].
   real(R8P), pointer     :: q_aux_gpu(:,:,:,:,:)=>null() !< Auxiliary variables [nb, i, j, k, nv_aux].
   real(R8P), pointer     :: flx_f_gpu(:,:,:,:,:)=>null() !< X-face fluxes [nb, 0:ni, 1:nj, 1:nk, nv].
   real(R8P), pointer     :: fly_f_gpu(:,:,:,:,:)=>null() !< Y-face fluxes [nb, 1:ni, 0:nj, 1:nk, nv].
   real(R8P), pointer     :: flz_f_gpu(:,:,:,:,:)=>null() !< Z-face fluxes [nb, 1:ni, 1:nj, 0:nk, nv].
   type(flume_seam_sync_object) :: seam                   !< Per-stage seam flux synchronisation of the limiter (issue #50).
   real(R8P), pointer     :: lam_gpu(:,:,:,:,:)=>null()   !< Positivity limiter cell factors (component 1, q-shaped for
                                                          !< the ghost exchange) [nb, i, j, k, nv]; with the limiter.
   real(R8P), pointer     :: q_inflow_gpu(:,:)=>null()    !< Conservative inflow state of each face [nv, 6].
   real(R8P), pointer     :: wall_sign_gpu(:,:)=>null()   !< Wall mirror sign per variable and direction [nv, 3].
   integer(I4P), pointer  :: divb_seam_gpu(:,:)=>null()   !< Seam faces flags of the div(B) history [nb, 6].
   ! host staging
   real(R8P), allocatable :: buf_5D_R8P(:,:,:,:,:)        !< Transposed copy buffer, extent identical to q_gpu.
   integer(I4P)           :: db5(2,5)=0_I4P               !< Device bounds of the transposed copies.
   integer(I4P)           :: hb5(2,5)=0_I4P               !< Host bounds of the transposed copies.
   ! dispatch
   procedure(compute_residuals_dev_interface), pass(self), pointer :: compute_residuals_dev=>null() !< Space operator.
   procedure(integrate_dev_interface),         pass(self), pointer :: integrate_dev=>null()         !< Time operator.
   contains
      ! public methods
      procedure, pass(self) :: accumulate_seam_fluxes  !< Accumulate the weighted seam face fluxes of one stage.
      procedure, pass(self) :: apply_floors            !< Apply the MHD positivity floors of a stage.
      procedure, pass(self) :: allocate_gpu            !< Allocate device data.
      procedure, pass(self) :: check_glm_ch            !< Check the GLM c_h against the fastest wave.
      procedure, pass(self) :: compute_conservation    !< Compute and save the conservation integrals.
      procedure, pass(self) :: compute_divb_history    !< Compute and save the div(B) norms (MHD).
      procedure, pass(self) :: compute_q_aux           !< Compute the auxiliary variables.
      procedure, pass(self) :: check_nonfinite         !< Stop on a non-finite committed state.
      procedure, pass(self) :: copy_cpu_gpu            !< Copy state and topology from host to device.
      procedure, pass(self) :: copy_gpu_cpu            !< Copy state from device to host.
      procedure, pass(self) :: copy_phi_gpu            !< Copy the immersed solids distance function to the device.
      procedure, pass(self) :: destroy                 !< Free device and host data.
      procedure, pass(self) :: initialize_flume        !< Initialize the FNL backend.
      procedure, pass(self) :: limit_positivity_dev    !< Apply the positivity limiter to the stage face fluxes.
      procedure, pass(self) :: seam_lists              !< Host lists of the seam skin cells and faces (issue #50).
      procedure, pass(self) :: seam_sync_blend_dev     !< Set the seam face fluxes with the seam factors (issue #50).
      procedure, pass(self) :: seam_sync_factors_dev   !< Recompute the seam cells' factors (issue #50).
      procedure, pass(self) :: seam_sync_fluxes_dev    !< Publish the donor states and the fine means (issue #50).
      procedure, pass(self) :: seam_sync_theta_dev     !< Compute the seam factors (issue #50).
      procedure, pass(self) :: save_residuals          !< Save residuals history.
      procedure, pass(self) :: save_simulation_data    !< Save fields, restart and diagnostics on their cadence.
      procedure, pass(self) :: set_boundary_conditions !< Set boundary conditions on the device crown maps.
      procedure, pass(self) :: update_ghost            !< Update ghost cells: local, MPI, boundary conditions.
      ! forest methods
      procedure, pass(self) :: advance_one_step_forest      !< Advance one full step (fast path).
      procedure, pass(self) :: after_topology_build_forest  !< Copy the forest-built seam and BC maps to the device.
      procedure, pass(self) :: apply_reflux_to_stage_forest !< Apply the reflux correction.
      procedure, pass(self) :: begin_stage_forest           !< Begin an integrator stage (staged path).
      procedure, pass(self) :: close_step_forest            !< Close a step (staged path).
      procedure, pass(self) :: compute_local_dt_forest      !< Compute the local stability-limited time step.
      procedure, pass(self) :: end_stage_forest             !< End an integrator stage (staged path).
      procedure, pass(self) :: fill_seam_from_peer_forest   !< Fill inter-realm seam ghosts from a peer.
      procedure, pass(self) :: pack_seam_cells_forest       !< Pack own cells of the cross-rank seam send rows.
      procedure, pass(self) :: unpack_seam_cells_forest     !< Unpack the cross-rank seam receive rows into own ghosts.
      procedure, pass(self) :: finalize_forest              !< Finalize the realm.
      procedure, pass(self) :: finalize_mpi_forest          !< Finalize the FNL MPI handler.
      procedure, pass(self) :: initialize_forest            !< Initialize the realm.
      procedure, pass(self) :: is_done_forest               !< Return true if the realm is done.
      procedure, pass(self) :: open_step_forest             !< Open a step (staged path).
      procedure, pass(self) :: post_step_forest             !< Post-step IO and diagnostics.
      procedure, pass(self) :: stages_per_step_forest       !< Return the integrator stages per step.
endtype flume_fnl_object

abstract interface
   subroutine compute_residuals_dev_interface(self, q_gpu, dq_gpu, s, flux_register)
   !< Compute the residuals on the device, space operator; ghost cells of `q_gpu` are filled inside.
   import :: flume_fnl_object, flux_register_object, I4P, R8P
   class(flume_fnl_object),     intent(inout)           :: self              !< The equation.
   real(R8P),                   intent(inout)           :: q_gpu(1:,         &
                                                                 1-self%ngc:,&
                                                                 1-self%ngc:,&
                                                                 1-self%ngc:,&
                                                                 1:)         !< Conservative variables.
   real(R8P),                   intent(inout)           :: dq_gpu(1:,         &
                                                                  1-self%ngc:,&
                                                                  1-self%ngc:,&
                                                                  1-self%ngc:,&
                                                                  1:)         !< Residuals.
   integer(I4P),                intent(in),    optional :: s                 !< Runge-Kutta stage.
   class(flux_register_object), intent(inout), optional :: flux_register     !< Forest's flux register for reflux.
   endsubroutine compute_residuals_dev_interface

   subroutine integrate_dev_interface(self)
   !< Integrate one time step on the device, time operator.
   import :: flume_fnl_object
   class(flume_fnl_object), intent(inout) :: self !< The equation.
   endsubroutine integrate_dev_interface
endinterface

contains
   ! public methods
   subroutine accumulate_seam_fluxes(self, s, flux_register)
   !< Accumulate the seam face fluxes of stage `s`, weighted by its SSP coefficient, into the forest's flux register.
   !<
   !< The register is host-side: each seam face skin is packed on the device and copied to the host (a few KB per face
   !< and stage), then routed by the shared `accumulate_seam_skin`, as on the CPU.
   class(flume_fnl_object),     intent(inout) :: self          !< The equation.
   integer(I4P),                intent(in)    :: s             !< Runge-Kutta stage.
   class(flux_register_object), intent(inout) :: flux_register !< Forest's flux register.
   real(R8P), pointer                         :: skin_gpu(:,:) !< Device face skin (nv, cells).
   real(R8P), allocatable                     :: skin(:,:)     !< Host face skin (nv, cells).
   integer(I4P)                               :: b, fec        !< Block, face counters.
   integer(I4P)                               :: cells         !< Skin cells number.
   integer(I4P)                               :: ierr          !< Error status.

   do b=1, self%blocks_number
      do fec=1, 6
         if (self%adam%maps%inter_realm_face_register_index(b, fec) == 0_I4P) cycle
         select case(fec)
         case(1_I4P, 2_I4P)
            cells = self%nj * self%nk
         case(3_I4P, 4_I4P)
            cells = self%ni * self%nk
         case default
            cells = self%ni * self%nj
         endselect
         call dev_alloc(fptr_dev=skin_gpu, lbounds=[1,1], ubounds=[self%nv,cells], ierr=ierr)
         if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate skin_gpu in accumulate_seam_fluxes')
         call pack_seam_skin_dev(fec=fec, b=b, ni=self%ni, nj=self%nj, nk=self%nk, nv=self%nv, flx_f_gpu=self%flx_f_gpu, &
                                 fly_f_gpu=self%fly_f_gpu, flz_f_gpu=self%flz_f_gpu, skin_gpu=skin_gpu)
         allocate(skin(self%nv,cells))
         call dev_memcpy_from_device(dst=skin, src=skin_gpu)
         call dev_free(skin_gpu, mydev)
         call self%accumulate_seam_skin(flux_register=flux_register, b=b, fec=fec, weight=self%rk%beta(s), skin=skin)
         deallocate(skin)
      enddo
   enddo
   endsubroutine accumulate_seam_fluxes

   subroutine apply_floors(self, q_gpu)
   !< Apply the MHD positivity floors to the interior of a stage state, before its ghost exchange (issue #41, 3.8).
   !<
   !< Euler has no floors (return before any work). MHD: the floored cells of the stage are logged by rank 0 when any;
   !< a non-positive density or pressure with the floors disabled (both zero) is fatal, reported with the global
   !< minimum density and pressure.
   class(flume_fnl_object), intent(inout) :: self      !< The equation.
   real(R8P),               intent(inout) :: q_gpu(1:,         &
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
      call apply_floors_mhd_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number,    &
                                gamma=self%physics%gamma, R=self%physics%R, rho_floor=self%physics%mhd%rho_floor,   &
                                p_floor=self%physics%mhd%p_floor, q_gpu=q_gpu,        &
                                floored=counts(1), nonpositive=counts(2), rho_min=mins(1), p_min=mins(2))
   case(MODEL_MHD_GLM)
      call apply_floors_mhd_glm_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                                    gamma=self%physics%gamma, R=self%physics%R, rho_floor=self%physics%mhd%rho_floor,   &
                                    p_floor=self%physics%mhd%p_floor, q_gpu=q_gpu,            &
                                    floored=counts(1), nonpositive=counts(2), rho_min=mins(1), p_min=mins(2))
   case(MODEL_MHD_EGLM)
      call apply_floors_mhd_eglm_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                                    gamma=self%physics%gamma, R=self%physics%R, rho_floor=self%physics%mhd%rho_floor,   &
                                    p_floor=self%physics%mhd%p_floor, q_gpu=q_gpu,            &
                                    floored=counts(1), nonpositive=counts(2), rho_min=mins(1), p_min=mins(2))
   case default
      call mpih_fnl%error_stop(msg=': no floors for physical model "'//self%physics%physical_model//'"')
   endselect
   call MPI_ALLREDUCE(MPI_IN_PLACE, counts, 2, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpih_fnl%error)
   call MPI_ALLREDUCE(MPI_IN_PLACE, mins, 2, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, mpih_fnl%error)
   if (counts(2) > 0_I4P .and. .not.(self%physics%mhd%rho_floor > 0._R8P .or. self%physics%mhd%p_floor > 0._R8P)) &
      call mpih_fnl%error_stop(msg=': '//trim(str(counts(2)))//' cells with a non-positive density or pressure at step '// &
                               trim(str(self%time%it))//' (min rho '//trim(str(mins(1)))//', min p '//               &
                               trim(str(mins(2)))//'); the [mhd] floors rho_floor, p_floor are disabled')
   if (counts(1) > 0_I4P .and. mpih_fnl%myrank == 0) &
      print '(A)', mpih_fnl%myrankstr//'MHD floors: '//trim(str(counts(1)))//' cells floored at step '// &
                   trim(str(self%time%it))//' (min rho '//trim(str(mins(1)))//', min p '//trim(str(mins(2)))//')'
   endsubroutine apply_floors

   subroutine allocate_gpu(self)
   !< Allocate device data (every `dev_alloc` checked) and the host staging buffer.
   class(flume_fnl_object), intent(inout) :: self       !< The equation.
   integer(I4P)                           :: ierr       !< Error status.
   integer(I4P)                           :: alloc_stat !< Allocation status.
   character(999)                         :: alloc_msg  !< Allocation error message.

   associate(nv=>self%physics%nv, nv_aux=>self%physics%nv_aux, nb=>self%nb, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, &
             nk=>self%nk)
   call dev_alloc(fptr_dev=self%q_gpu, ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,nv], lbounds=[1,1-ngc,1-ngc,1-ngc,1], &
                  init_value=0._R8P, ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate q_gpu in flume_fnl_object%allocate_gpu')
   call dev_alloc(fptr_dev=self%dq_gpu, ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,nv], lbounds=[1,1-ngc,1-ngc,1-ngc,1], &
                  init_value=0._R8P, ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate dq_gpu in flume_fnl_object%allocate_gpu')
   call dev_alloc(fptr_dev=self%q_aux_gpu, ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,nv_aux], lbounds=[1,1-ngc,1-ngc,1-ngc,1], &
                  init_value=0._R8P, ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate q_aux_gpu in flume_fnl_object%allocate_gpu')
   call dev_alloc(fptr_dev=self%flx_f_gpu, ubounds=[nb,ni,nj,nk,nv], lbounds=[1,0,1,1,1], init_value=0._R8P, ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate flx_f_gpu in flume_fnl_object%allocate_gpu')
   call dev_alloc(fptr_dev=self%fly_f_gpu, ubounds=[nb,ni,nj,nk,nv], lbounds=[1,1,0,1,1], init_value=0._R8P, ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate fly_f_gpu in flume_fnl_object%allocate_gpu')
   call dev_alloc(fptr_dev=self%flz_f_gpu, ubounds=[nb,ni,nj,nk,nv], lbounds=[1,1,1,0,1], init_value=0._R8P, ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate flz_f_gpu in flume_fnl_object%allocate_gpu')
   if (self%numerics%positivity_limiter == POSITIVITY_LIMITER_CELL) then
      call dev_alloc(fptr_dev=self%lam_gpu, ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,nv], lbounds=[1,1-ngc,1-ngc,1-ngc,1], &
                     init_value=1._R8P, ierr=ierr)
      if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate lam_gpu in flume_fnl_object%allocate_gpu')
   endif
   call dev_assign_to_device(src=self%bc%q_inflow,  dst=self%q_inflow_gpu)
   call dev_assign_to_device(src=self%bc%wall_sign, dst=self%wall_sign_gpu)
   allocate(self%buf_5D_R8P(1:nb,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nv), stat=alloc_stat, errmsg=alloc_msg)
   if (alloc_stat /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate buf_5D_R8P: '//trim(alloc_msg))
   self%db5(1,:) = [1 , 1-ngc , 1-ngc , 1-ngc , 1 ]
   self%db5(2,:) = [nb, ni+ngc, nj+ngc, nk+ngc, nv]
   self%hb5(1,:) = [1 , 1-ngc , 1-ngc , 1-ngc , 1 ]
   self%hb5(2,:) = [nv, ni+ngc, nj+ngc, nk+ngc, nb]
   endassociate
   endsubroutine allocate_gpu

   subroutine check_nonfinite(self)
   !< Stop the run when the committed state holds a non-finite (NaN or infinite) value (issue #45): a run could otherwise
   !< finish with NaN fields and exit 0. One device pass over the interior cells, the model selected outside the kernels;
   !< the state is copied to the host only to locate a failure.
   class(flume_fnl_object), intent(inout) :: self !< The equation.
   integer(I8P)                           :: n    !< Non-finite values number of this rank.

   select case(self%physics%model)
   case(MODEL_EULER)
      call count_nonfinite_euler_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                                     q_gpu=self%q_gpu, n=n)
   case(MODEL_MHD)
      call count_nonfinite_mhd_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                                   q_gpu=self%q_gpu, n=n)
   case(MODEL_MHD_GLM)
      call count_nonfinite_mhd_glm_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc,                  &
                                       blocks_number=self%blocks_number, q_gpu=self%q_gpu, n=n)
   case(MODEL_MHD_EGLM)
      call count_nonfinite_mhd_eglm_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc,                  &
                                       blocks_number=self%blocks_number, q_gpu=self%q_gpu, n=n)
   case default
      call mpih_fnl%error_stop(msg=': no FNL kernels for physical model "'//self%physics%physical_model//'"')
   endselect
   n = self%nonfinite_total(n_local=n)
   if (n > 0_I8P) then
      call self%copy_gpu_cpu
      call self%stop_nonfinite(n_total=n)
   endif
   endsubroutine check_nonfinite

   subroutine check_glm_ch(self)
   !< Check the GLM cleaning speed against the fastest wave of the committed state (issue #41, section 3.5); a no-op
   !< without GLM.
   class(flume_fnl_object), intent(inout) :: self      !< The equation.
   real(R8P)                              :: speed_max !< Fastest wave speed of this rank.

   select case(self%physics%model)
   case(MODEL_EULER, MODEL_MHD)
      return
   case(MODEL_MHD_GLM)
      call compute_speed_max_mhd_glm_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc,                              &
                                         blocks_number=self%blocks_number, gamma=self%physics%gamma, R=self%physics%R, &
                                         is_null=self%adam%grid%null_xyz, q_gpu=self%q_gpu, speed_max=speed_max)
   case(MODEL_MHD_EGLM)
      call compute_speed_max_mhd_eglm_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc,                              &
                                         blocks_number=self%blocks_number, gamma=self%physics%gamma, R=self%physics%R, &
                                         is_null=self%adam%grid%null_xyz, q_gpu=self%q_gpu, speed_max=speed_max)
   case default
      call mpih_fnl%error_stop(msg=': no FNL kernels for physical model "'//self%physics%physical_model//'"')
   endselect
   call self%report_glm_speed(speed_max=speed_max)
   endsubroutine check_glm_ch

   subroutine compute_conservation(self)
   !< Compute the volume integrals of the conservative variables on the device and save them on their cadence.
   class(flume_fnl_object), intent(inout) :: self         !< The equation.
   real(R8P), allocatable                 :: integrals(:) !< Volume integrals [nv].

   if (.not.self%time%is_to_save(cadence=self%diagnostics%conservation_history_save)) return
   allocate(integrals(self%physics%nv))
   select case(self%physics%model)
   case(MODEL_EULER)
      call compute_conservation_euler_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc,                   &
                                          blocks_number=self%blocks_number, dxyz_gpu=self%field_fnl%dxyz_gpu, &
                                          q_gpu=self%q_gpu, integrals=integrals)
   case(MODEL_MHD)
      call compute_conservation_mhd_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc,                   &
                                        blocks_number=self%blocks_number, dxyz_gpu=self%field_fnl%dxyz_gpu, &
                                        q_gpu=self%q_gpu, integrals=integrals)
   case(MODEL_MHD_GLM)
      call compute_conservation_mhd_glm_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc,                   &
                                            blocks_number=self%blocks_number, dxyz_gpu=self%field_fnl%dxyz_gpu, &
                                            q_gpu=self%q_gpu, integrals=integrals)
   case(MODEL_MHD_EGLM)
      call compute_conservation_mhd_eglm_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc,                   &
                                            blocks_number=self%blocks_number, dxyz_gpu=self%field_fnl%dxyz_gpu, &
                                            q_gpu=self%q_gpu, integrals=integrals)
   case default
      call mpih_fnl%error_stop(msg=': no FNL kernels for physical model "'//self%physics%physical_model//'"')
   endselect
   call MPI_ALLREDUCE(MPI_IN_PLACE, integrals, size(integrals), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih_fnl%error)
   call self%diagnostics%save_conservation_row(it=self%time%it, time=self%time%time*self%units%time_output(),        &
                                               integrals=integrals*self%output_factors(self%q_name)*            &
                                                         self%units%length_output()**3)
   endsubroutine compute_conservation

   subroutine compute_divb_history(self, realm)
   !< Compute and save the div(B) norms of the committed state on the diagnostics cadence (issue #41, D-10); a no-op
   !< for Euler. The ghost cells are refreshed first and, with sibling realms, the inter-realm seam ghosts are refilled
   !< from the peers (`update_ghost` does not fill them, issue #31): the stencils of the seam-local cells read them. The
   !< seam faces flags are rebuilt and copied to the device at every call (6 nb integers): at step 0 the forest may not
   !< have registered the seams yet.
   class(flume_fnl_object), intent(inout)                   :: self     !< The equation.
   class(realm_object),     intent(inout), optional, target :: realm(:) !< Sibling realms.
   real(R8P)                                                :: norms(3) !< Norms of this rank.

   if (self%physics%model == MODEL_EULER) return
   if (.not.self%time%is_to_save(cadence=self%diagnostics%conservation_history_save)) return
   call self%update_ghost(q_gpu=self%q_gpu)
   if (present(realm)) call seam_fill_all(self=self, realm=realm) ! local and cross-rank seam rows (issue #40)
   call self%set_divb_seam
   if (associated(self%divb_seam_gpu)) then
      call dev_free(self%divb_seam_gpu, mydev)
      nullify(self%divb_seam_gpu)
   endif
   call dev_assign_to_device(src=self%divb_seam, dst=self%divb_seam_gpu)
   associate(hs=>self%fdv_half_stencils(1))
   if (hs > min(self%ngc, 3_I4P)) &
      call mpih_fnl%error_stop(msg=': the div(B) stencil ([fdv].(fdv_order)) exceeds the ghost cells or sixth order')
   select case(self%physics%model)
   case(MODEL_MHD)
      call compute_divb_norms_mhd_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc,                              &
                                      blocks_number=self%blocks_number, hs=hs, band=self%ngc,                        &
                                      dxyz_gpu=self%field_fnl%dxyz_gpu, is_null=self%adam%grid%null_xyz,             &
                                      seam_gpu=self%divb_seam_gpu, q_gpu=self%q_gpu, divb_max=norms(1),               &
                                      divb_l1=norms(2), divb_seam_max=norms(3))
   case(MODEL_MHD_GLM)
      call compute_divb_norms_mhd_glm_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc,                          &
                                          blocks_number=self%blocks_number, hs=hs, band=self%ngc,                    &
                                          dxyz_gpu=self%field_fnl%dxyz_gpu, is_null=self%adam%grid%null_xyz,         &
                                          seam_gpu=self%divb_seam_gpu, q_gpu=self%q_gpu, divb_max=norms(1),           &
                                          divb_l1=norms(2), divb_seam_max=norms(3))
   case(MODEL_MHD_EGLM)
      call compute_divb_norms_mhd_eglm_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc,                          &
                                          blocks_number=self%blocks_number, hs=hs, band=self%ngc,                    &
                                          dxyz_gpu=self%field_fnl%dxyz_gpu, is_null=self%adam%grid%null_xyz,         &
                                          seam_gpu=self%divb_seam_gpu, q_gpu=self%q_gpu, divb_max=norms(1),           &
                                          divb_l1=norms(2), divb_seam_max=norms(3))
   case default
      call mpih_fnl%error_stop(msg=': no FNL kernels for physical model "'//self%physics%physical_model//'"')
   endselect
   endassociate
   call self%report_divb(norms=norms)
   endsubroutine compute_divb_history

   subroutine compute_q_aux(self, q_gpu)
   !< Compute the auxiliary variables on the device, ghost cells included.
   class(flume_fnl_object), intent(inout) :: self              !< The equation.
   real(R8P),               intent(in)    :: q_gpu(1:,         &
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1:)         !< Conservative variables.

   select case(self%physics%model)
   case(MODEL_EULER)
      call compute_q_aux_euler_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                                   gamma=self%physics%gamma, R=self%physics%R, q_gpu=q_gpu, q_aux_gpu=self%q_aux_gpu)
   case(MODEL_MHD)
      call compute_q_aux_mhd_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                                 gamma=self%physics%gamma, R=self%physics%R, q_gpu=q_gpu, q_aux_gpu=self%q_aux_gpu)
   case(MODEL_MHD_GLM)
      call compute_q_aux_mhd_glm_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc,                           &
                                     blocks_number=self%blocks_number, gamma=self%physics%gamma, R=self%physics%R, &
                                     q_gpu=q_gpu, q_aux_gpu=self%q_aux_gpu)
   case(MODEL_MHD_EGLM)
      call compute_q_aux_mhd_eglm_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc,                           &
                                     blocks_number=self%blocks_number, gamma=self%physics%gamma, R=self%physics%R, &
                                     q_gpu=q_gpu, q_aux_gpu=self%q_aux_gpu)
   case default
      call mpih_fnl%error_stop(msg=': no FNL kernels for physical model "'//self%physics%physical_model//'"')
   endselect
   endsubroutine compute_q_aux

   subroutine copy_cpu_gpu(self, verbose)
   !< Copy state and topology (coordinates, maps) from host to device.
   class(flume_fnl_object), intent(inout)        :: self    !< The equation.
   logical,                 intent(in), optional :: verbose !< Trigger verbose output.

   call dev_memcpy_to_device(bb=self%db5, ij=[1,5], tb=self%hb5, dst=self%q_gpu, src=self%q, buf=self%buf_5D_R8P)
   call self%field_fnl%copy_cpu_gpu(field=self%adam%field, maps=self%adam%maps, verbose=verbose)
   endsubroutine copy_cpu_gpu

   subroutine copy_gpu_cpu(self)
   !< Copy state and residuals from device to host.
   class(flume_fnl_object), intent(inout) :: self !< The equation.

   call dev_memcpy_from_device(bb=self%db5, ij=[1,5], tb=self%hb5, dst=self%q, src=self%q_gpu, buf=self%buf_5D_R8P)
   call dev_memcpy_from_device(bb=self%db5, ij=[1,5], tb=self%hb5, dst=self%dq, src=self%dq_gpu, buf=self%buf_5D_R8P)
   endsubroutine copy_gpu_cpu

   subroutine copy_phi_gpu(self)
   !< Copy the immersed solids distance function (computed on the host, the solids are static) to the device, transposed
   !< from `(s, i, j, k, b)` to `(b, i, j, k, s)`; the staging buffer has the extent of the device array (issue #31).
   class(flume_fnl_object), intent(inout) :: self           !< The equation.
   real(R8P), allocatable                 :: buf(:,:,:,:,:) !< Transposed staging buffer.
   integer(I4P)                           :: db(2,5)        !< Device bounds.
   integer(I4P)                           :: hb(2,5)        !< Host bounds.

   if (self%ib%solids_number == 0_I4P) return
   associate(ns=>self%ib%solids_number+1, nb=>self%nb, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk)
   db(1,:) = [1 , 1-ngc , 1-ngc , 1-ngc , 1 ]
   db(2,:) = [nb, ni+ngc, nj+ngc, nk+ngc, ns]
   hb(1,:) = [1 , 1-ngc , 1-ngc , 1-ngc , 1 ]
   hb(2,:) = [ns, ni+ngc, nj+ngc, nk+ngc, nb]
   allocate(buf(1:nb,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:ns))
   call dev_memcpy_to_device(bb=db, ij=[1,5], tb=hb, dst=self%ib_fnl%phi_gpu, src=self%ib%phi, buf=buf)
   endassociate
   endsubroutine copy_phi_gpu

   subroutine destroy(self)
   !< Free device and host data: own device buffers, the FNL helpers (library teardown) and the common data.
   class(flume_fnl_object), intent(inout) :: self !< The equation.

   call free_gpu(self%q_gpu)
   call free_gpu(self%dq_gpu)
   call free_gpu(self%q_aux_gpu)
   call free_gpu(self%flx_f_gpu)
   call free_gpu(self%fly_f_gpu)
   call free_gpu(self%flz_f_gpu)
   call free_gpu(self%lam_gpu)
   if (associated(self%q_inflow_gpu)) then
      call dev_free(self%q_inflow_gpu, mydev)
      nullify(self%q_inflow_gpu)
   endif
   if (associated(self%wall_sign_gpu)) then
      call dev_free(self%wall_sign_gpu, mydev)
      nullify(self%wall_sign_gpu)
   endif
   if (associated(self%divb_seam_gpu)) then
      call dev_free(self%divb_seam_gpu, mydev)
      nullify(self%divb_seam_gpu)
   endif
   if (allocated(self%buf_5D_R8P)) deallocate(self%buf_5D_R8P)
   call self%rk_fnl%destroy()
   call self%weno_fnl%destroy()
   call self%ib_fnl%destroy()
   call self%field_fnl%destroy()
   call self%destroy_common
   nullify(self%compute_residuals_dev)
   nullify(self%integrate_dev)
   contains
      subroutine free_gpu(fptr)
      !< Free one device buffer, if allocated.
      real(R8P), pointer, intent(inout) :: fptr(:,:,:,:,:) !< Device buffer.

      if (associated(fptr)) then
         call dev_free(fptr, mydev)
         nullify(fptr)
      endif
      endsubroutine free_gpu
   endsubroutine destroy

   subroutine initialize_flume(self, filename, realms_number)
   !< Initialize the FNL backend: MPI and device (once per process), device budget, common data, FNL helpers,
   !< device data, dispatch.
   class(flume_fnl_object), intent(inout), target :: self               !< The equation.
   character(*),            intent(in)            :: filename           !< Input file name.
   integer(I4P),            intent(in), optional  :: realms_number      !< Realm count; divides the device budget.
   logical                                        :: is_mpi_initialized !< MPI initialization status.
   integer(I4P)                                   :: realms_number_     !< Realm count, local variable.
   real(R8P)                                      :: memory_avail_      !< Per-realm device budget (GB).

   realms_number_ = 1_I4P ; if (present(realms_number)) realms_number_ = max(1_I4P, realms_number)
   if (.not.mpih_fnl_is_initialized) then
      call MPI_INITIALIZED(is_mpi_initialized, mpih_fnl%error)
      call mpih_fnl%initialize(do_mpi_init=.not.is_mpi_initialized, do_device_init=.true., verbose=.true.)
      mpih_fnl_is_initialized = .true.
   endif
   memory_avail_ = real(mpih_fnl%dev_memory_total, R8P) / 1e9_R8P / real(realms_number_, R8P)
   call self%flume_common_object%initialize(filename=filename, memory_avail=memory_avail_, verbose=.true.)
   if (realms_number_ > 1_I4P .and. self%units%is_active) &
      call mpih_fnl%error_stop(msg=': [reference] is not supported on multi-realm runs (each realm converts its own '// &
                                   'input and nothing yet checks that the realms share the references)')
   if (realms_number_ > 1_I4P .and. self%numerics%positivity_limiter == POSITIVITY_LIMITER_CELL) &
      call mpih_fnl%error_stop(msg=': [numerics].(positivity_limiter)=cell is not supported on multi-realm runs (the '// &
                                   'inter-realm seam faces carry no limiting factor)')
   call self%field_fnl%initialize(grid=self%adam%grid, field=self%adam%field, maps=self%adam%maps, verbose=.true.)
   call self%ib_fnl%initialize(grid=self%adam%grid, field=self%adam%field, ib=self%ib)
   call self%rk_fnl%initialize(grid=self%adam%grid, field=self%adam%field, rk=self%rk)
   call self%weno_fnl%initialize(weno=self%weno)
   if (self%numerics%scheme_space == SCHEME_SPACE_WENO_RIEMANN) call self%weno_fnl%initialize_interpolation(weno=self%weno)
   call self%allocate_gpu
   select case(self%numerics%scheme_space)
   case(SCHEME_SPACE_WENO)
      self%compute_residuals_dev => compute_residuals_weno_dev
   case(SCHEME_SPACE_WENO_RIEMANN)
      self%compute_residuals_dev => compute_residuals_riemann_dev
   case default
      call mpih_fnl%error_stop(msg=': no FNL space operator for scheme_space "'//self%numerics%scheme_space//'"')
   endselect
   select case(self%rk%scheme)
   case(RK_1, RK_2, RK_3)
      self%integrate_dev => integrate_rk_ls_dev
   case(RK_SSP_11, RK_SSP_22, RK_SSP_33, RK_SSP_54)
      self%integrate_dev => integrate_rk_ssp_dev
   case default
      call mpih_fnl%error_stop(msg=': no FNL time integrator for [runge_kutta].(scheme) "'//trim(self%rk%scheme)//'"')
   endselect
   endsubroutine initialize_flume

   subroutine limit_positivity_dev(self, q_gpu, flux_register)
   !< Apply the positivity limiter to the device face fluxes of the stage (issue #47, D-9; device twin of the CPU
   !< `limit_positivity`): the model's factors kernel, the ghost exchange of the factors on the device (intra-realm copies
   !< and GPU-direct MPI only: a physical-boundary ghost keeps the factor 1 it was allocated with), the blending of the
   !< active directions; the counters are reduced over the ranks and logged. With 2:1 AMR seams the seam flux is
   !< synchronised as on the CPU (issue #50, `adam_flume_seam_sync_object`): the seam passes run on the host over the
   !< seam cells, fed by device gathers, and scatter their factors and face fluxes back.
   class(flume_fnl_object),     intent(inout)           :: self          !< The equation.
   real(R8P),                   intent(in)              :: q_gpu(1:,         &
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
   integer(I4P)                                         :: di(3,3)       !< Unit steps of the directions.

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
   if (seam) call self%seam_sync_fluxes_dev(q_gpu=q_gpu)
   di = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nb=>self%blocks_number, gamma=>self%physics%gamma, &
             ch=>self%physics%mhd%glm_ch, is_null=>self%adam%grid%null_xyz, dt=>self%time%dt, hs=>self%weno%S,      &
             dxyz_gpu=>self%field_fnl%dxyz_gpu)
   select case(self%physics%model)
   case(MODEL_EULER)
      call compute_positivity_factors_euler_dev(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, gamma=gamma, ch=ch,    &
                                                damping=0._R8P, hs=hs, dt=dt, dxyz_gpu=dxyz_gpu, is_null=is_null,      &
                                                q_gpu=q_gpu, q_aux_gpu=self%q_aux_gpu, flx_gpu=self%flx_f_gpu,         &
                                                fly_gpu=self%fly_f_gpu, flz_gpu=self%flz_f_gpu, lam_gpu=self%lam_gpu,  &
                                                bad=counts(1), nonfinite=counts(5))
   case(MODEL_MHD)
      call compute_positivity_factors_mhd_dev(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, gamma=gamma, ch=ch,      &
                                              damping=0._R8P, hs=hs, dt=dt, dxyz_gpu=dxyz_gpu, is_null=is_null,        &
                                              q_gpu=q_gpu, q_aux_gpu=self%q_aux_gpu, flx_gpu=self%flx_f_gpu,           &
                                              fly_gpu=self%fly_f_gpu, flz_gpu=self%flz_f_gpu, lam_gpu=self%lam_gpu,    &
                                              bad=counts(1), nonfinite=counts(5))
   case(MODEL_MHD_EGLM)
      call compute_positivity_factors_mhd_eglm_dev(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, gamma=gamma, ch=ch, &
                                                   damping=self%physics%mhd%glm_damping, hs=hs, dt=dt,                 &
                                                   dxyz_gpu=dxyz_gpu, is_null=is_null, q_gpu=q_gpu,                    &
                                                   q_aux_gpu=self%q_aux_gpu, flx_gpu=self%flx_f_gpu,                   &
                                                   fly_gpu=self%fly_f_gpu, flz_gpu=self%flz_f_gpu,                     &
                                                   lam_gpu=self%lam_gpu, bad=counts(1), nonfinite=counts(5))
   case default
      call mpih_fnl%error_stop(msg=': no FNL positivity limiter for physical model "'//self%physics%physical_model//'"')
   endselect
   if (seam) call self%seam_sync_factors_dev(q_gpu=q_gpu, dbad=counts(1), dnonfinite=counts(5))
   call self%field_fnl%update_ghost_local_gpu(q_gpu=self%lam_gpu)
   call self%field_fnl%update_ghost_mpi_gpu(comm_map_send_ptr_ghost=self%adam%maps%comm_map_send_ptr_ghost, &
                                            comm_map_recv_ptr_ghost=self%adam%maps%comm_map_recv_ptr_ghost, &
                                            q_gpu=self%lam_gpu)
   if (seam) call self%seam_sync_theta_dev
   do d=1, 3
      if (is_null(d)) cycle
      select case(self%physics%model)
      case(MODEL_EULER)
         select case(d)
         case(1)
            call blend_positivity_fluxes_euler_dev(d=d, di=di(1,d), dj=di(2,d), dk=di(3,d), ni=ni, nj=nj, nk=nk,      &
                                                   ngc=ngc, blocks_number=nb, gamma=gamma, ch=ch, q_gpu=q_gpu,      &
                                                   lam_gpu=self%lam_gpu, fl_gpu=self%flx_f_gpu, limited=counts(1+d))
         case(2)
            call blend_positivity_fluxes_euler_dev(d=d, di=di(1,d), dj=di(2,d), dk=di(3,d), ni=ni, nj=nj, nk=nk,      &
                                                   ngc=ngc, blocks_number=nb, gamma=gamma, ch=ch, q_gpu=q_gpu,      &
                                                   lam_gpu=self%lam_gpu, fl_gpu=self%fly_f_gpu, limited=counts(1+d))
         case(3)
            call blend_positivity_fluxes_euler_dev(d=d, di=di(1,d), dj=di(2,d), dk=di(3,d), ni=ni, nj=nj, nk=nk,      &
                                                   ngc=ngc, blocks_number=nb, gamma=gamma, ch=ch, q_gpu=q_gpu,      &
                                                   lam_gpu=self%lam_gpu, fl_gpu=self%flz_f_gpu, limited=counts(1+d))
         endselect
      case(MODEL_MHD)
         select case(d)
         case(1)
            call blend_positivity_fluxes_mhd_dev(d=d, di=di(1,d), dj=di(2,d), dk=di(3,d), ni=ni, nj=nj, nk=nk,        &
                                                 ngc=ngc, blocks_number=nb, gamma=gamma, ch=ch, q_gpu=q_gpu,        &
                                                 lam_gpu=self%lam_gpu, fl_gpu=self%flx_f_gpu, limited=counts(1+d))
         case(2)
            call blend_positivity_fluxes_mhd_dev(d=d, di=di(1,d), dj=di(2,d), dk=di(3,d), ni=ni, nj=nj, nk=nk,        &
                                                 ngc=ngc, blocks_number=nb, gamma=gamma, ch=ch, q_gpu=q_gpu,        &
                                                 lam_gpu=self%lam_gpu, fl_gpu=self%fly_f_gpu, limited=counts(1+d))
         case(3)
            call blend_positivity_fluxes_mhd_dev(d=d, di=di(1,d), dj=di(2,d), dk=di(3,d), ni=ni, nj=nj, nk=nk,        &
                                                 ngc=ngc, blocks_number=nb, gamma=gamma, ch=ch, q_gpu=q_gpu,        &
                                                 lam_gpu=self%lam_gpu, fl_gpu=self%flz_f_gpu, limited=counts(1+d))
         endselect
      case(MODEL_MHD_EGLM)
         select case(d)
         case(1)
            call blend_positivity_fluxes_mhd_eglm_dev(d=d, di=di(1,d), dj=di(2,d), dk=di(3,d), ni=ni, nj=nj, nk=nk,   &
                                                      ngc=ngc, blocks_number=nb, gamma=gamma, ch=ch, q_gpu=q_gpu,   &
                                                      lam_gpu=self%lam_gpu, fl_gpu=self%flx_f_gpu, limited=counts(1+d))
         case(2)
            call blend_positivity_fluxes_mhd_eglm_dev(d=d, di=di(1,d), dj=di(2,d), dk=di(3,d), ni=ni, nj=nj, nk=nk,   &
                                                      ngc=ngc, blocks_number=nb, gamma=gamma, ch=ch, q_gpu=q_gpu,   &
                                                      lam_gpu=self%lam_gpu, fl_gpu=self%fly_f_gpu, limited=counts(1+d))
         case(3)
            call blend_positivity_fluxes_mhd_eglm_dev(d=d, di=di(1,d), dj=di(2,d), dk=di(3,d), ni=ni, nj=nj, nk=nk,   &
                                                      ngc=ngc, blocks_number=nb, gamma=gamma, ch=ch, q_gpu=q_gpu,   &
                                                      lam_gpu=self%lam_gpu, fl_gpu=self%flz_f_gpu, limited=counts(1+d))
         endselect
      endselect
   enddo
   endassociate
   if (seam) call self%seam_sync_blend_dev
   call MPI_ALLREDUCE(MPI_IN_PLACE, counts, 5, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpih_fnl%error)
   if (sum(counts) > 0_I4P .and. mpih_fnl%myrank == 0) &
      print '(A)', mpih_fnl%myrankstr//'positivity limiter: '//trim(str(sum(counts(2:4))))//' faces limited, '// &
                   trim(str(counts(1)))//' inadmissible backbones at step '//trim(str(self%time%it))
   if (counts(5) > 0_I4P .and. mpih_fnl%myrank == 0) &
      print '(A)', mpih_fnl%myrankstr//'positivity limiter (non-finite): '//trim(str(counts(5)))// &
                   ' cells with a non-finite high-order flux took the backbone at step '//trim(str(self%time%it))
   endsubroutine limit_positivity_dev

   subroutine seam_lists(self, coarse, ccell, cface, cdst, fine, fcell, fface, fdst, fcc, fw)
   !< Host lists of this rank's seam skin cells (issue #50): `coarse` (or `fine`) selects the coarse (fine) seam faces;
   !< each entry has its interior cell `(i, j, k, b)`, its face `(axis, i, j, k, b)` in the flux arrays' indexing and its
   !< skin destination (`cdst`: skin cell of the coarse skins; `fdst`: compact fine store, `fcc`: the coarse skin cell
   !< covering it, `fw`: its weight in the coarse mean, 1/4, or 1/2 along the unrefined z of a quadtree, issue #46). The
   !< order is block, face, skin cell, as on the CPU.
   class(flume_fnl_object),   intent(in)  :: self                       !< The equation.
   logical,                   intent(in)  :: coarse, fine               !< Lists to build.
   integer(I4P), allocatable, intent(out) :: ccell(:,:), cface(:,:)     !< Coarse cells [4, n], faces [5, n].
   integer(I4P), allocatable, intent(out) :: cdst(:)                    !< Coarse skin destinations.
   integer(I4P), allocatable, intent(out) :: fcell(:,:), fface(:,:)     !< Fine cells [4, n], faces [5, n].
   integer(I4P), allocatable, intent(out) :: fdst(:), fcc(:)            !< Fine store slots, covering coarse skin cells.
   real(R8P),    allocatable, intent(out), optional :: fw(:)            !< Fine weights in the coarse means.
   integer(I4P)                           :: b, fec, s, c, m, nc, mc, mf !< Counters.
   integer(I4P)                           :: axis, sg, ioff, joff        !< Face axis, side, fine quadrant offsets.
   integer(I4P)                           :: ratios(2)                   !< Tangential refinement ratios (inner, outer).
   integer(I4P)                           :: i, j, k                     !< Cell indexes.

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, idx=>self%adam%maps%inter_realm_face_register_index)
   do m=1, 2 ! count, then fill
      mc = 0 ; mf = 0
      do b=1, self%blocks_number
         do fec=1, 6
            s = idx(b, fec)
            if (s == 0_I4P) cycle
            if (self%seam%off(abs(s)) < 0_I4P) cycle
            if ((s > 0_I4P .and. .not.coarse) .or. (s < 0_I4P .and. .not.fine)) cycle
            axis = (fec + 1_I4P) / 2_I4P ; sg = merge(-1_I4P, 1_I4P, mod(fec, 2_I4P) == 1_I4P)
            ioff = 0_I4P ; joff = 0_I4P
            if (s < 0_I4P .and. allocated(self%adam%maps%amr_seam_quadrant)) then
               ioff = self%adam%maps%amr_seam_quadrant(1, b, fec) ; joff = self%adam%maps%amr_seam_quadrant(2, b, fec)
            endif
            nc = seam_face_cells(fec=fec, ni=ni, nj=nj, nk=nk)
            do c=1, nc
               call seam_skin_cell(axis=axis, sgn=sg, ni=ni, nj=nj, nk=nk, c=c, i=i, j=j, k=k)
               if (s > 0_I4P) then
                  mc = mc + 1
                  if (m == 1) cycle
                  ccell(:,mc) = [i, j, k, b]
                  cface(:,mc) = [axis, merge(merge(0_I4P, ni, sg < 0_I4P), i, axis == 1),   &
                                       merge(merge(0_I4P, nj, sg < 0_I4P), j, axis == 2),   &
                                       merge(merge(0_I4P, nk, sg < 0_I4P), k, axis == 3), b]
                  cdst(mc) = self%seam%off(s) + c
               else
                  mf = mf + 1
                  if (m == 1) cycle
                  fcell(:,mf) = [i, j, k, b]
                  fface(:,mf) = [axis, merge(merge(0_I4P, ni, sg < 0_I4P), i, axis == 1),   &
                                       merge(merge(0_I4P, nj, sg < 0_I4P), j, axis == 2),   &
                                       merge(merge(0_I4P, nk, sg < 0_I4P), k, axis == 3), b]
                  fdst(mf) = self%seam%fine_off(b, fec) + c
                  ratios = face_tangential_ratios(fec=fec, refine_ratio=self%adam%maps%refine_ratio)
                  fcc(mf) = self%seam%off(-s) + seam_fine_to_coarse(fec=fec, ni=ni, nj=nj, nk=nk, ioff=ioff, joff=joff, &
                                                                    ri=ratios(1), ro=ratios(2), c=c)
                  if (present(fw)) fw(mf) = 1._R8P / real(ratios(1) * ratios(2), R8P)
               endif
            enddo
         enddo
      enddo
      if (m == 1) allocate(ccell(4,mc), cface(5,mc), cdst(mc), fcell(4,mf), fface(5,mf), fdst(mf), fcc(mf))
      if (m == 1 .and. present(fw)) allocate(fw(mf))
   enddo
   endassociate
   endsubroutine seam_lists

   subroutine seam_sync_fluxes_dev(self, q_gpu)
   !< Seam synchronisation, first phase (issue #50; host twin of the CPU `seam_sync_fluxes`, fed by device gathers).
   class(flume_fnl_object), intent(inout) :: self                         !< The equation.
   real(R8P),               intent(in)    :: q_gpu(1:,         &
                                                     1-self%ngc:,&
                                                     1-self%ngc:,&
                                                     1-self%ngc:,&
                                                     1:)                    !< Conservative variables of the stage.
   integer(I4P), allocatable              :: ccell(:,:), cface(:,:), cdst(:) !< Coarse lists.
   integer(I4P), allocatable              :: fcell(:,:), fface(:,:), fdst(:), fcc(:) !< Fine lists.
   real(R8P),    allocatable              :: fw(:)                           !< Fine weights in the coarse means.
   real(R8P),    allocatable              :: buf(:,:), fhi(:,:)              !< Gathered states, fine face fluxes.
   real(R8P)                              :: qL(self%physics%nv,1), qR(self%physics%nv,1), fl(self%physics%nv,1) !< Face.
   integer(I4P)                           :: m                               !< Counter.

   associate(nv=>self%physics%nv)
   if (allocated(self%adam%maps%inter_realm_face_register_index)) then
      call self%seam_lists(coarse=.true., ccell=ccell, cface=cface, cdst=cdst, fine=.true., fcell=fcell, fface=fface, &
                           fdst=fdst, fcc=fcc, fw=fw)
      allocate(buf(nv,size(cdst)))
      call seam_gather_cells(cells=ccell, nv=nv, ngc=self%ngc, a_gpu=q_gpu, out=buf)
      do m=1, size(cdst)
         self%seam%qc(1:nv,cdst(m)) = buf(:,m)
      enddo
      deallocate(buf)
   endif
   call self%seam%reduce_states
   if (allocated(self%adam%maps%inter_realm_face_register_index)) then
      allocate(buf(nv,size(fdst)), fhi(nv,size(fdst)))
      call seam_gather_cells(cells=fcell, nv=nv, ngc=self%ngc, a_gpu=q_gpu, out=buf)
      call seam_gather_faces(faces=fface, nv=nv, flx_f_gpu=self%flx_f_gpu, fly_f_gpu=self%fly_f_gpu, &
                             flz_f_gpu=self%flz_f_gpu, out=fhi)
      do m=1, size(fdst)
         if (fface(2+fface(1,m)-1,m) == 0_I4P) then ! the fine block's minimum face: the donor is the left state
            qL(:,1) = self%seam%qc(1:nv,fcc(m)) ; qR(:,1) = buf(:,m)
         else
            qL(:,1) = buf(:,m) ; qR(:,1) = self%seam%qc(1:nv,fcc(m))
         endif
         select case(self%physics%model)
         case(MODEL_EULER)
            call compute_backbone_fluxes_euler_host(n=1, d=fface(1,m), gamma=self%physics%gamma,              &
                                                    ch=self%physics%mhd%glm_ch, qL=qL, qR=qR, flo=fl)
         case(MODEL_MHD)
            call compute_backbone_fluxes_mhd_host(n=1, d=fface(1,m), gamma=self%physics%gamma,                &
                                                  ch=self%physics%mhd%glm_ch, qL=qL, qR=qR, flo=fl)
         case(MODEL_MHD_EGLM)
            call compute_backbone_fluxes_mhd_eglm_host(n=1, d=fface(1,m), gamma=self%physics%gamma,           &
                                                       ch=self%physics%mhd%glm_ch, qL=qL, qR=qR, flo=fl)
         endselect
         self%seam%fine_lo(1:nv,fdst(m)) = fl(:,1)
         self%seam%fine_hi(1:nv,fdst(m)) = fhi(:,m)
         self%seam%flo(1:nv,fcc(m)) = self%seam%flo(1:nv,fcc(m)) + fw(m) * fl(:,1)
         self%seam%fhi(1:nv,fcc(m)) = self%seam%fhi(1:nv,fcc(m)) + fw(m) * fhi(:,m)
      enddo
   endif
   endassociate
   call self%seam%reduce_fluxes
   endsubroutine seam_sync_fluxes_dev

   subroutine seam_sync_factors_dev(self, q_gpu, dbad, dnonfinite)
   !< Seam synchronisation, second phase (issue #50; host twin of the CPU `seam_sync_factors`): gather the seam cells'
   !< stencils, auxiliary variables and face fluxes, recompute their factors with the seam faces overridden on the
   !< host, scatter them to the device factors.
   class(flume_fnl_object), intent(inout) :: self                      !< The equation.
   real(R8P),               intent(in)    :: q_gpu(1:,         &
                                                     1-self%ngc:,&
                                                     1-self%ngc:,&
                                                     1-self%ngc:,&
                                                     1:)                 !< Conservative variables of the stage.
   integer(I4P),            intent(inout) :: dbad                      !< Inadmissible backbones count.
   integer(I4P),            intent(inout) :: dnonfinite                !< Non-finite cells count.
   integer(I4P), allocatable              :: cell(:,:), face(:,:)      !< Seam cells (i, j, k, b), their faces.
   logical,      allocatable              :: omask(:,:)                !< Overridden faces.
   real(R8P),    allocatable              :: olo(:,:,:), ohi(:,:,:)    !< Override fluxes.
   real(R8P),    allocatable              :: qs(:,:,:,:), qa(:,:)      !< Gathered stencils, auxiliary variables.
   real(R8P),    allocatable              :: fh(:,:), ds(:,:), lam(:,:) !< Gathered face fluxes, steps; factors.
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
                  face(:,6*(n-1)+1) = [1, i-1, j, k, b] ; face(:,6*(n-1)+2) = [1, i, j, k, b]
                  face(:,6*(n-1)+3) = [2, i, j-1, k, b] ; face(:,6*(n-1)+4) = [2, i, j, k, b]
                  face(:,6*(n-1)+5) = [3, i, j, k-1, b] ; face(:,6*(n-1)+6) = [3, i, j, k, b]
                  ds(:,n) = self%adam%field%dxyz(:,b)
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
         if (n == 0) exit
         allocate(cell(4,n), face(5,6*n), omask(6,n), olo(nv,6,n), ohi(nv,6,n), ds(3,n))
         olo = 0._R8P ; ohi = 0._R8P
      endif
   enddo
   if (n > 0) then
      allocate(qs(nv,2*self%weno%S+1,3,n), qa(self%physics%nv_aux,n), fh(nv,6*n), lam(1,n))
      call seam_gather_stencils(cells=cell, nv=nv, ngc=self%ngc, s=self%weno%S, q_gpu=q_gpu, out=qs)
      call seam_gather_cells(cells=cell, nv=self%physics%nv_aux, ngc=self%ngc, a_gpu=self%q_aux_gpu, out=qa)
      call seam_gather_faces(faces=face, nv=nv, flx_f_gpu=self%flx_f_gpu, fly_f_gpu=self%fly_f_gpu, &
                             flz_f_gpu=self%flz_f_gpu, out=fh)
      associate(gamma=>self%physics%gamma, ch=>self%physics%mhd%glm_ch, is_null=>self%adam%grid%null_xyz, &
                dt=>self%time%dt, hs=>self%weno%S)
      select case(self%physics%model)
      case(MODEL_EULER)
         call compute_seam_positivity_factors_euler_host(ncells=n, hs=hs, gamma=gamma, ch=ch, damping=0._R8P, dt=dt, &
                                                         ds=ds, is_null=is_null, qs=qs, qa=qa,                     &
                                                         fh=reshape(fh, [nv, 6, n]), omask=omask, olo=olo,         &
                                                         ohi=ohi, lam=lam(1,:), dbad=db, dnonfinite=dn)
      case(MODEL_MHD)
         call compute_seam_positivity_factors_mhd_host(ncells=n, hs=hs, gamma=gamma, ch=ch, damping=0._R8P, dt=dt, &
                                                       ds=ds, is_null=is_null, qs=qs, qa=qa,                     &
                                                       fh=reshape(fh, [nv, 6, n]), omask=omask, olo=olo,         &
                                                       ohi=ohi, lam=lam(1,:), dbad=db, dnonfinite=dn)
      case(MODEL_MHD_EGLM)
         call compute_seam_positivity_factors_mhd_eglm_host(ncells=n, hs=hs, gamma=gamma, ch=ch,                    &
                                                            damping=self%physics%mhd%glm_damping, dt=dt, ds=ds,    &
                                                            is_null=is_null, qs=qs, qa=qa,                         &
                                                            fh=reshape(fh, [nv, 6, n]), omask=omask, olo=olo,      &
                                                            ohi=ohi, lam=lam(1,:), dbad=db, dnonfinite=dn)
      endselect
      endassociate
      call seam_scatter_cells(cells=cell, nv=1_I4P, ngc=self%ngc, in=lam, a_gpu=self%lam_gpu)
      dbad = dbad + db
      dnonfinite = dnonfinite + dn
   endif
   endassociate
   endsubroutine seam_sync_factors_dev

   subroutine seam_sync_theta_dev(self)
   !< Seam synchronisation, third phase (issue #50; host twin of the CPU `seam_sync_theta`).
   class(flume_fnl_object), intent(inout) :: self                         !< The equation.
   integer(I4P), allocatable              :: ccell(:,:), cface(:,:), cdst(:) !< Coarse lists.
   integer(I4P), allocatable              :: fcell(:,:), fface(:,:), fdst(:), fcc(:) !< Fine lists.
   real(R8P),    allocatable              :: lam(:,:)                        !< Gathered factors.
   integer(I4P)                           :: m                               !< Counter.

   if (allocated(self%adam%maps%inter_realm_face_register_index)) then
      call self%seam_lists(coarse=.true., ccell=ccell, cface=cface, cdst=cdst, fine=.true., fcell=fcell, fface=fface, &
                           fdst=fdst, fcc=fcc)
      allocate(lam(1,size(fdst)))
      call seam_gather_cells(cells=fcell, nv=1_I4P, ngc=self%ngc, a_gpu=self%lam_gpu, out=lam)
      do m=1, size(fdst)
         self%seam%lmin(fcc(m)) = min(self%seam%lmin(fcc(m)), lam(1,m))
      enddo
      deallocate(lam)
   endif
   call self%seam%reduce_factors
   if (allocated(self%adam%maps%inter_realm_face_register_index)) then
      allocate(lam(1,size(cdst)))
      call seam_gather_cells(cells=ccell, nv=1_I4P, ngc=self%ngc, a_gpu=self%lam_gpu, out=lam)
      do m=1, size(cdst)
         self%seam%th(cdst(m)) = max(0._R8P, min(1._R8P, lam(1,m), self%seam%lmin(cdst(m))))
      enddo
   endif
   call self%seam%reduce_theta
   endsubroutine seam_sync_theta_dev

   subroutine seam_sync_blend_dev(self)
   !< Seam synchronisation, last phase (issue #50; host twin of the CPU `seam_sync_blend`): the seam face fluxes with
   !< the seam factor, scattered to the device face fluxes.
   class(flume_fnl_object), intent(inout) :: self                         !< The equation.
   integer(I4P), allocatable              :: ccell(:,:), cface(:,:), cdst(:) !< Coarse lists.
   integer(I4P), allocatable              :: fcell(:,:), fface(:,:), fdst(:), fcc(:) !< Fine lists.
   real(R8P),    allocatable              :: fl(:,:)                         !< Seam face fluxes.
   real(R8P)                              :: th                              !< Seam factor.
   integer(I4P)                           :: m                               !< Counter.

   if (.not.allocated(self%adam%maps%inter_realm_face_register_index)) return
   associate(nv=>self%physics%nv)
   call self%seam_lists(coarse=.true., ccell=ccell, cface=cface, cdst=cdst, fine=.true., fcell=fcell, fface=fface, &
                        fdst=fdst, fcc=fcc)
   allocate(fl(nv,size(cdst)))
   do m=1, size(cdst)
      th = self%seam%th(cdst(m))
      fl(:,m) = self%seam%flo(1:nv,cdst(m))
      if (th > 0._R8P) fl(:,m) = fl(:,m) + th * (self%seam%fhi(1:nv,cdst(m)) - fl(:,m))
   enddo
   call seam_scatter_faces(faces=cface, nv=nv, in=fl, flx_f_gpu=self%flx_f_gpu, fly_f_gpu=self%fly_f_gpu, &
                           flz_f_gpu=self%flz_f_gpu)
   deallocate(fl)
   allocate(fl(nv,size(fdst)))
   do m=1, size(fdst)
      th = self%seam%th(fcc(m))
      fl(:,m) = self%seam%fine_lo(1:nv,fdst(m))
      if (th > 0._R8P) fl(:,m) = fl(:,m) + th * (self%seam%fine_hi(1:nv,fdst(m)) - fl(:,m))
   enddo
   call seam_scatter_faces(faces=fface, nv=nv, in=fl, flx_f_gpu=self%flx_f_gpu, fly_f_gpu=self%fly_f_gpu, &
                           flz_f_gpu=self%flz_f_gpu)
   endassociate
   endsubroutine seam_sync_blend_dev

   subroutine save_residuals(self)
   !< Save residuals history (L2 norm of dq on the device, MPI-reduced, rank 0 writes).
   class(flume_fnl_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: v    !< Counter.

   if (.not.self%time%is_to_save(cadence=self%io%residuals_save)) return
   call compute_normL2_residuals_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, nv=self%nv, &
                                     blocks_number=self%blocks_number, dq_gpu=self%dq_gpu, norm=self%adam%field%residuals)
   do v=1, self%nv
      call MPI_ALLREDUCE(MPI_IN_PLACE, self%adam%field%residuals(v), 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih_fnl%error)
      self%adam%field%residuals(v) = sqrt(self%adam%field%residuals(v)) / sqrt(real(self%ni*self%nj*self%nk, R8P))
   enddo
   if (mpih_fnl%myrank == 0) call self%io%save_residuals(it=self%time%it, time=self%time%time*self%units%time_output(), &
                                                         blocks_number=self%blocks_number,                             &
                                                         residuals=self%adam%field%residuals*                          &
                                                                   self%output_factors(self%dq_name))
   endsubroutine save_residuals

   subroutine save_simulation_data(self, realm)
   !< Save fields, restart, slices and conservation history, each on its own cadence; state copied to host only when
   !< saved.
   class(flume_fnl_object), intent(inout)                   :: self      !< The equation.
   class(realm_object),     intent(inout), optional, target :: realm(:)  !< Sibling realms.
   logical                                                  :: is_slices !< Slices save step.

   is_slices = self%slices%is_to_save(it=self%time%it, it_max=self%time%it_max, time=self%time%time, &
                                      time_max=self%time%time_max)
   if (self%time%is_to_save(cadence=self%io%it_save) .or. self%time%is_to_save(cadence=self%io%restart_save) .or. &
       is_slices) then
      ! the stage order: inter-realm seams first, then the intra-realm ghosts and the boundary conditions, so that the
      ! saved ghosts (edges and corners included) are those the residual stencils read
      if (present(realm)) call seam_fill_all(self=self, realm=realm)
      call self%update_ghost(q_gpu=self%q_gpu)
      call self%copy_gpu_cpu
      if (self%time%is_to_save(cadence=self%io%it_save)) call self%save_xh5f(with_ghost=.true.)
      if (self%time%is_to_save(cadence=self%io%restart_save)) call self%save_restart_files
      if (is_slices) call self%save_slices
   endif
   call self%compute_conservation
   endsubroutine save_simulation_data

   subroutine set_boundary_conditions(self, q_gpu)
   !< Set boundary conditions on the device crown maps: three passes (rows beyond one, two, three realm faces), crown by
   !< crown within each (see the CPU `set_boundary_conditions`, issue #65 P0).
   class(flume_fnl_object), intent(inout) :: self              !< The equation.
   real(R8P),               intent(inout) :: q_gpu(1:,         &
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1:)         !< Conservative variables.
   integer(I4P)                           :: crown             !< Crown counter.
   integer(I4P)                           :: pass              !< Pass: realm faces the rows lie beyond.
   integer(I4P)                           :: face_kind(6)      !< Kind of each realm face, BC_SEAM for a seam.

   if (.not.associated(self%field_fnl%maps%local_map_bc_crown_gpu)) return
   face_kind = self%bc%bc_type
   where (self%adam%maps%seam_face) face_kind = BC_SEAM
   do pass=1, 3
      do crown=1, self%ngc
         call set_boundary_conditions_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, nv=self%nv, crown=crown, &
                                          pass=pass, face_kind=face_kind,                                          &
                                          local_map_bc_crown_gpu=self%field_fnl%maps%local_map_bc_crown_gpu,       &
                                          q_inflow_gpu=self%q_inflow_gpu, wall_sign_gpu=self%wall_sign_gpu, q_gpu=q_gpu)
      enddo
   enddo
   endsubroutine set_boundary_conditions

   subroutine update_ghost(self, q_gpu)
   !< Update ghost cells on the device: intra-realm local copies, GPU-direct MPI exchange, boundary conditions, then the
   !< positivity blend of the inadmissible face ghosts (issue #50, D2; see the CPU `update_ghost`).
   class(flume_fnl_object), intent(inout) :: self              !< The equation.
   real(R8P),               intent(inout) :: q_gpu(1:,         &
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1:)         !< Conservative variables.
   integer(I4P)                           :: blended           !< Face ghosts blended toward the interior.

   call self%field_fnl%update_ghost_local_gpu(q_gpu=q_gpu)
   call self%field_fnl%update_ghost_mpi_gpu(comm_map_send_ptr_ghost=self%adam%maps%comm_map_send_ptr_ghost, &
                                            comm_map_recv_ptr_ghost=self%adam%maps%comm_map_recv_ptr_ghost, &
                                            q_gpu=q_gpu)
   call self%set_boundary_conditions(q_gpu=q_gpu)
   blended = 0_I4P
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nb=>self%blocks_number, is_null=>self%adam%grid%null_xyz)
   select case(self%physics%model)
   case(MODEL_EULER)
      call blend_inadmissible_ghosts_euler_dev(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, is_null=is_null, &
                                               q_gpu=q_gpu, blended=blended)
   case(MODEL_MHD)
      call blend_inadmissible_ghosts_mhd_dev(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, is_null=is_null, &
                                             q_gpu=q_gpu, blended=blended)
   case(MODEL_MHD_GLM)
      call blend_inadmissible_ghosts_mhd_glm_dev(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, is_null=is_null, &
                                                 q_gpu=q_gpu, blended=blended)
   case(MODEL_MHD_EGLM)
      call blend_inadmissible_ghosts_mhd_eglm_dev(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, is_null=is_null, &
                                                  q_gpu=q_gpu, blended=blended)
   endselect
   endassociate
   if (blended > 0_I4P) print '(A)', mpih_fnl%myrankstr//'ghost positivity: '//trim(str(blended))// &
                                     ' inadmissible face ghosts blended toward the interior at step '//trim(str(self%time%it))
   endsubroutine update_ghost

   ! forest methods
   subroutine advance_one_step_forest(self, dt)
   !< Advance one full step of size `dt` (fast path: single realm, no AMR seam faces).
   class(flume_fnl_object), intent(inout) :: self !< The equation.
   real(R8P),               intent(in)    :: dt   !< Time step from the forest's global reduction.

   self%time%it = self%time%it + 1_I4P
   self%time%dt = dt
   if ((self%time%it_max <= 0_I4P) .and. (self%time%time + dt > self%time%time_max)) &
      self%time%dt = self%time%time_max - self%time%time
   call self%integrate_dev
   self%time%time = self%time%time + self%time%dt
   call self%time%print_progress(nodes_number=self%adam%tree%nodes_number)
   endsubroutine advance_one_step_forest

   subroutine after_topology_build_forest(self)
   !< Copy the host maps the forest built after the realm initialization to the device: the inter-realm seam map and
   !< buffers (`maps%seam_local_*`) and the BC crown map, whose seam rows the forest rewrote to `BC_SEAM` (issue #37;
   !< without it the seam fill kernel reads an unset device map and the BC kernel extrapolates over the seam ghosts).
   class(flume_fnl_object), intent(inout) :: self !< The equation.

   call self%field_fnl%maps%copy_cpu_gpu(maps=self%adam%maps)
   endsubroutine after_topology_build_forest

   subroutine apply_reflux_to_stage_forest(self, stage, dt, flux_register)
   !< Apply the Berger-Colella reflux correction to the committed `q_gpu` (the forest calls it once per step, after
   !< `close_step_forest`, with the final stage); device twin of the CPU apply.
   !<
   !< For every register face whose coarse side this realm and this rank own, the host mismatch slab
   !< `F_coarse - F_fine_sum` (step fluxes, `sum_s beta_s F_s`) is copied to the device and added, scaled by
   !< `sgn dt / dx_coarse`, to the coarse skin cells.
   !<
   !< The scale uses the step the update used, `time%dt`: on the last step of a time-driven run the realm caps it to
   !< land on `time_max`, and the forest's `dt` argument is the uncapped value (with it the correction was off by the
   !< ratio of the two, 1.3e-9 of the mass on sod-amr; issue #37).
   class(flume_fnl_object),     intent(inout) :: self           !< The equation.
   integer(I4P),                intent(in)    :: stage          !< Integrator stage.
   real(R8P),                   intent(in)    :: dt             !< Forest time step (unused: see the scale below).
   class(flux_register_object), intent(in)    :: flux_register  !< Forest's flux register.
   real(R8P), allocatable                     :: delta(:,:)     !< Host flux mismatch (nv, cells).
   real(R8P), pointer                         :: delta_gpu(:,:) !< Device flux mismatch (nv, cells).
   integer(I4P)                               :: f              !< Face counter.
   integer(I4P)                               :: axis, sgn      !< Face normal axis and side.
   integer(I4P)                               :: ierr           !< Error status.

   if (.not.self%numerics%reflux) return
   if (.not.flux_register%is_initialized_) return
   if (flux_register%nfaces == 0_I4P) return
   if (.not.allocated(flux_register%face)) return
   if (stage /= self%rk%nrk) return
   do f=1, flux_register%nfaces
      associate(face=>flux_register%face(f))
      if (face%coarse_realm /= self%realm_index .or. face%coarse_rank /= mpih_fnl%myrank) cycle
      if (.not.(allocated(face%F_coarse) .and. allocated(face%F_fine_sum))) cycle
      call face_axis_sign(face_code=face%coarse_face, axis=axis, sgn=sgn)
      if (axis == 0_I4P) call mpih_fnl%error_stop(msg=': malformed coarse face code of register face '//trim(str(f)))
      allocate(delta(self%nv,face%nface_cells))
      delta = face%F_coarse(1:self%nv,:,1) - face%F_fine_sum(1:self%nv,:,1)
      call dev_alloc(fptr_dev=delta_gpu, lbounds=[1,1], ubounds=[self%nv,face%nface_cells], ierr=ierr)
      if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate delta_gpu in apply_reflux_to_stage_forest')
      call dev_memcpy_to_device(dst=delta_gpu, src=delta)
      call apply_reflux_face_dev(axis=axis, sgn=sgn, b=face%coarse_block, ni=self%ni, nj=self%nj, nk=self%nk,  &
                                 ngc=self%ngc, nv=self%nv, nface_cells=face%nface_cells,                       &
                                 scale=real(sgn, R8P) * self%time%dt / self%adam%field%dxyz(axis,face%coarse_block), &
                                 delta_gpu=delta_gpu, q_gpu=self%q_gpu)
      call dev_free(delta_gpu, mydev)
      deallocate(delta)
      endassociate
   enddo
   endsubroutine apply_reflux_to_stage_forest

   subroutine begin_stage_forest(self, k, K_total, dt, realm)
   !< Begin integrator stage `k` (staged path): publish the stage and compute its state.
   class(flume_fnl_object), intent(inout)                   :: self     !< The equation.
   integer(I4P),            intent(in)                      :: k        !< Stage index (1..K_total).
   integer(I4P),            intent(in)                      :: K_total  !< Forest-wide stage count for this step.
   real(R8P),               intent(in)                      :: dt       !< Time step from the forest.
   class(realm_object),     intent(inout), optional, target :: realm(:) !< Sibling realms (contract parity).

   self%stage_active = k
   call rk_compute_stage(self, s=k)
   endsubroutine begin_stage_forest

   subroutine close_step_forest(self, dt)
   !< Close a step (staged path): assemble q, save residuals, advance time, clear the active stage.
   class(flume_fnl_object), intent(inout) :: self !< The equation.
   real(R8P),               intent(in)    :: dt   !< Time step from the forest (the local capped value is time%dt).

   call compute_rk_ssp_residual_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, nv=self%nv,                     &
                                    blocks_number=self%blocks_number, nrk=self%rk%nrk, beta_gpu=self%rk_fnl%beta_gpu, &
                                    q_rk_gpu=self%rk_fnl%q_rk_gpu, dq_gpu=self%dq_gpu)
   call rk_update_q(self)
   call self%save_residuals
   self%time%time = self%time%time + self%time%dt
   call self%time%print_progress(nodes_number=self%adam%tree%nodes_number)
   self%stage_active = 0_I4P
   endsubroutine close_step_forest

   subroutine compute_local_dt_forest(self, dt_local)
   !< Compute the local stability-limited time step on the device, `dt = CFL / max(sum_d (|u_d| + a) / dx_d)`; with GLM,
   !< also `dt <= CFL / (c_h max sum_d 1 / dx_d)` (`glm_lambda`, host data, issue #41, section 3.5).
   class(flume_fnl_object), intent(in)  :: self       !< The equation.
   real(R8P),               intent(out) :: dt_local   !< Local stability-limited time step.
   real(R8P)                            :: lambda_max !< Maximum of sum_d (|u_d| + a) / dx_d.

   select case(self%physics%model)
   case(MODEL_EULER)
      call compute_lambda_max_euler_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                                        gamma=self%physics%gamma, R=self%physics%R, dxyz_gpu=self%field_fnl%dxyz_gpu,       &
                                        is_null=self%adam%grid%null_xyz, q_gpu=self%q_gpu, lambda_max=lambda_max)
   case(MODEL_MHD)
      call compute_lambda_max_mhd_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                                      gamma=self%physics%gamma, R=self%physics%R, dxyz_gpu=self%field_fnl%dxyz_gpu,       &
                                      is_null=self%adam%grid%null_xyz, q_gpu=self%q_gpu, lambda_max=lambda_max)
   case(MODEL_MHD_GLM)
      call compute_lambda_max_mhd_glm_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc,                              &
                                          blocks_number=self%blocks_number, gamma=self%physics%gamma, R=self%physics%R, &
                                          dxyz_gpu=self%field_fnl%dxyz_gpu, is_null=self%adam%grid%null_xyz,           &
                                          q_gpu=self%q_gpu, lambda_max=lambda_max)
      lambda_max = max(lambda_max, self%glm_lambda())
   case(MODEL_MHD_EGLM)
      call compute_lambda_max_mhd_eglm_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc,                              &
                                          blocks_number=self%blocks_number, gamma=self%physics%gamma, R=self%physics%R, &
                                          dxyz_gpu=self%field_fnl%dxyz_gpu, is_null=self%adam%grid%null_xyz,           &
                                          q_gpu=self%q_gpu, lambda_max=lambda_max)
      lambda_max = max(lambda_max, self%glm_lambda())
   case default
      call mpih_fnl%error_stop(msg=': no FNL kernels for physical model "'//self%physics%physical_model//'"')
   endselect
   dt_local = huge(1._R8P)
   if (lambda_max > 0._R8P) dt_local = self%time%CFL / lambda_max
   endsubroutine compute_local_dt_forest

   subroutine end_stage_forest(self, k, K_total, dt, realm, flux_register)
   !< End integrator stage `k` (staged path): residuals on the stage state, then stage assignment.
   class(flume_fnl_object),     intent(inout)                   :: self          !< The equation.
   integer(I4P),                intent(in)                      :: k             !< Stage index (1..K_total).
   integer(I4P),                intent(in)                      :: K_total       !< Forest-wide stage count.
   real(R8P),                   intent(in)                      :: dt            !< Time step from the forest.
   class(realm_object),         intent(inout), optional, target :: realm(:)      !< Sibling realms (parity only).
   class(flux_register_object), intent(inout), optional         :: flux_register !< Forest's flux register.

   if (present(flux_register)) then
      call self%compute_residuals_dev(q_gpu=self%rk_fnl%q_rk_gpu(:,:,:,:,:,k), dq_gpu=self%dq_gpu, s=k, &
                                      flux_register=flux_register)
   else
      call self%compute_residuals_dev(q_gpu=self%rk_fnl%q_rk_gpu(:,:,:,:,:,k), dq_gpu=self%dq_gpu, s=k)
   endif
   call rk_assign_stage(self, s=k)
   endsubroutine end_stage_forest

   subroutine fill_seam_from_peer_forest(self, peer, p_idx)
   !< Fill this realm's seam ghosts for peer slot `p_idx` from the peer's interior, on the active device buffers
   !< (`q_gpu` when `stage_active == 0`, else the active stage of `q_rk_gpu`) of both realms.
   class(flume_fnl_object), intent(inout)         :: self      !< The equation.
   class(realm_object),     intent(in),    target :: peer      !< Peer realm.
   integer(I4P),            intent(in)            :: p_idx     !< Peer slot.
   integer(I4P)                                   :: row_start !< First seam map row of the peer.
   integer(I4P)                                   :: row_count !< Seam map rows of the peer.
   integer(I4P)                                   :: s_self    !< Own active stage (0: committed state).
   integer(I4P)                                   :: s_peer    !< Peer active stage (0: committed state).

   row_start = self%adam%maps%seam_local_peer_row_start(p_idx)
   row_count = self%adam%maps%seam_local_peer_row_count(p_idx)
   s_self    = self%stage_active
   select type(peer)
   class is(flume_fnl_object)
      s_peer = peer%stage_active
      associate(rows_gpu=>self%field_fnl%maps%seam_local_map_ghost_cell_gpu)
      if (s_self > 0_I4P .and. s_peer > 0_I4P) then
         call fill_seam_copy_dev(row_start=row_start, row_count=row_count, nv=self%nv, ngc=self%ngc, &
                                 rows_gpu=rows_gpu, src_gpu=peer%rk_fnl%q_rk_gpu(:,:,:,:,:,s_peer),  &
                                 dst_gpu=self%rk_fnl%q_rk_gpu(:,:,:,:,:,s_self))
      elseif (s_self > 0_I4P) then
         call fill_seam_copy_dev(row_start=row_start, row_count=row_count, nv=self%nv, ngc=self%ngc, &
                                 rows_gpu=rows_gpu, src_gpu=peer%q_gpu,                              &
                                 dst_gpu=self%rk_fnl%q_rk_gpu(:,:,:,:,:,s_self))
      elseif (s_peer > 0_I4P) then
         call fill_seam_copy_dev(row_start=row_start, row_count=row_count, nv=self%nv, ngc=self%ngc, &
                                 rows_gpu=rows_gpu, src_gpu=peer%rk_fnl%q_rk_gpu(:,:,:,:,:,s_peer),  &
                                 dst_gpu=self%q_gpu)
      else
         call fill_seam_copy_dev(row_start=row_start, row_count=row_count, nv=self%nv, ngc=self%ngc, &
                                 rows_gpu=rows_gpu, src_gpu=peer%q_gpu, dst_gpu=self%q_gpu)
      endif
      endassociate
   class default
      call mpih_fnl%error_stop(msg=': flume_fnl_object%fill_seam_from_peer_forest: peer realm is not a flume_fnl_object')
   endselect
   endsubroutine fill_seam_from_peer_forest

   subroutine pack_seam_cells_forest(self, p_idx, buf)
   !< Pack the ghost values of the seam send rows of peer slot `p_idx` from the active device buffer (`q_gpu` when
   !< `stage_active == 0`, else the active stage of `q_rk_gpu`) into the host `buf`, `nv` values per row: cell copies, or
   !< the 2:1 interpolations and restrictions of a refined seam (issues #40, #52).
   class(flume_fnl_object), intent(in)  :: self       !< The equation.
   integer(I4P),            intent(in)  :: p_idx      !< Peer slot (the realm owning the ghosts).
   real(R8P),               intent(out) :: buf(:)     !< Packed values.
   real(R8P), pointer                   :: buf_gpu(:) !< Device packed values.
   integer(I4P)                         :: row_start  !< First send row.
   integer(I4P)                         :: row_count  !< Send rows.
   integer(I4P)                         :: ierr       !< Error status.

   row_start = self%adam%maps%seam_mpi_send_row_start(p_idx)
   row_count = self%adam%maps%seam_mpi_send_row_count(p_idx)
   if (row_count == 0_I4P) return
   call dev_alloc(fptr_dev=buf_gpu, lbounds=[1], ubounds=[self%nv*row_count], ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate buf_gpu in pack_seam_cells_forest')
   if (self%stage_active > 0_I4P) then
      call pack_seam_rows_dev(row_start=row_start, row_count=row_count, nv=self%nv, ngc=self%ngc,              &
                              regime=self%adam%maps%seam_ghost_fill,                                          &
                              rows_gpu=self%field_fnl%maps%seam_mpi_send_cell_gpu,                            &
                              q_gpu=self%rk_fnl%q_rk_gpu(:,:,:,:,:,self%stage_active), buf_gpu=buf_gpu)
   else
      call pack_seam_rows_dev(row_start=row_start, row_count=row_count, nv=self%nv, ngc=self%ngc,              &
                              regime=self%adam%maps%seam_ghost_fill,                                          &
                              rows_gpu=self%field_fnl%maps%seam_mpi_send_cell_gpu, q_gpu=self%q_gpu, buf_gpu=buf_gpu)
   endif
   call dev_memcpy_from_device(dst=buf, src=buf_gpu)
   call dev_free(buf_gpu, mydev)
   endsubroutine pack_seam_cells_forest

   subroutine unpack_seam_cells_forest(self, p_idx, buf)
   !< Unpack the host `buf` into this realm's seam ghosts of the cross-rank receive rows of peer slot `p_idx`, on the
   !< active device buffer, `nv` values per row (issue #40).
   class(flume_fnl_object), intent(inout) :: self       !< The equation.
   integer(I4P),            intent(in)    :: p_idx      !< Peer slot (the realm owning the cells).
   real(R8P),               intent(in)    :: buf(:)     !< Packed values.
   real(R8P), pointer                     :: buf_gpu(:) !< Device packed values.
   integer(I4P)                           :: row_start  !< First receive row.
   integer(I4P)                           :: row_count  !< Receive rows.
   integer(I4P)                           :: ierr       !< Error status.

   row_start = self%adam%maps%seam_mpi_recv_row_start(p_idx)
   row_count = self%adam%maps%seam_mpi_recv_row_count(p_idx)
   if (row_count == 0_I4P) return
   call dev_alloc(fptr_dev=buf_gpu, lbounds=[1], ubounds=[self%nv*row_count], ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate buf_gpu in unpack_seam_cells_forest')
   call dev_memcpy_to_device(dst=buf_gpu, src=buf)
   if (self%stage_active > 0_I4P) then
      call unpack_seam_rows_dev(row_start=row_start, row_count=row_count, nv=self%nv, ngc=self%ngc,            &
                                rows_gpu=self%field_fnl%maps%seam_mpi_recv_cell_gpu, buf_gpu=buf_gpu,          &
                                q_gpu=self%rk_fnl%q_rk_gpu(:,:,:,:,:,self%stage_active))
   else
      call unpack_seam_rows_dev(row_start=row_start, row_count=row_count, nv=self%nv, ngc=self%ngc,            &
                                rows_gpu=self%field_fnl%maps%seam_mpi_recv_cell_gpu, buf_gpu=buf_gpu, q_gpu=self%q_gpu)
   endif
   call dev_free(buf_gpu, mydev)
   endsubroutine unpack_seam_cells_forest

   subroutine finalize_forest(self)
   !< Finalize the realm: close the output files and free device and host data (MPI is finalized once by the forest).
   class(flume_fnl_object), intent(inout) :: self !< The equation.

   call self%io%close_file_residuals
   call self%diagnostics%close_file
   call self%destroy
   endsubroutine finalize_forest

   subroutine finalize_mpi_forest(self)
   !< Finalize the FNL MPI handler (called once by the forest after every realm is finalized).
   class(flume_fnl_object), intent(inout) :: self !< The equation (carries no MPI state).

   call mpih_fnl%finalize
   endsubroutine finalize_mpi_forest

   subroutine initialize_forest(self, filename, realms_number, memory_avail, nv, verbose)
   !< Initialize the realm (issue #35, section 6.1): backend init, IC (or restart) and initial AMR on the host, device
   !< topology and state sync, ghost update, initial output, output files open, AMR lock.
   class(flume_fnl_object), intent(inout)           :: self          !< The equation.
   character(*),            intent(in)              :: filename      !< Input file name.
   integer(I4P),            intent(in),    optional :: realms_number !< Realm count; divides the device budget.
   real(R8P),               intent(in),    optional :: memory_avail  !< Unused: the budget comes from the device.
   integer(I4P),            intent(in),    optional :: nv            !< Unused: nv is decided by the physics.
   logical,                 intent(in),    optional :: verbose       !< Unused: initialization is always verbose.
   integer(I4P)                                     :: i             !< Counter.

   call self%initialize_flume(filename=filename, realms_number=realms_number)
   if (self%io%restart) then
      call mpih_fnl%print_message('restart simulation from "'//trim(self%io%restart_basename)//'" files')
      call self%load_restart_files(t=self%time%it, time=self%time%time)
      call self%compute_phi
   else
      do i=1, self%ic%amr_iterations
         call self%ic%set_initial_conditions(field=self%adam%field, q=self%q)
         call self%compute_phi
         call self%amr_update
      enddo
      call self%ic%set_initial_conditions(field=self%adam%field, q=self%q)
      call self%compute_phi
      call self%adam%make_comm_local_maps_ghost_bc
      self%time%time = 0._R8P
      self%time%it   = 0_I4P
   endif
   call self%set_glm_damping
   call self%copy_cpu_gpu(verbose=.true.)
   call self%copy_phi_gpu
   call self%update_ghost(q_gpu=self%q_gpu)
   call self%compute_q_aux(q_gpu=self%q_gpu)
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
   class(flume_fnl_object), intent(in)  :: self !< The equation.
   logical,                 intent(out) :: done !< True if the realm is done.

   done = self%time%is_done()
   endsubroutine is_done_forest

   subroutine open_step_forest(self, dt)
   !< Open a step (staged path): time bookkeeping and Runge-Kutta stages initialization.
   class(flume_fnl_object), intent(inout) :: self !< The equation.
   real(R8P),               intent(in)    :: dt   !< Time step from the forest.

   self%time%it = self%time%it + 1_I4P
   self%time%dt = dt
   if ((self%time%it_max <= 0_I4P) .and. (self%time%time + dt > self%time%time_max)) &
      self%time%dt = self%time%time_max - self%time%time
   call self%rk_fnl%initialize_stages(grid=self%adam%grid, field=self%adam%field, q_gpu=self%q_gpu)
   endsubroutine open_step_forest

   subroutine post_step_forest(self, dt, t, it, do_save_state, do_save_residuals, do_save_restart, do_amr, realm)
   !< Post-step work: the non-finite state check, the GLM c_h check, the div(B) history (MHD), fields, restart and
   !< conservation history on their cadence.
   class(flume_fnl_object), intent(inout)                   :: self              !< The equation.
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
   class(flume_fnl_object), intent(in) :: self !< The equation.
   integer(I4P)                        :: K    !< Integrator stages per step.

   select case(self%rk%scheme)
   case(RK_SSP_11, RK_SSP_22, RK_SSP_33, RK_SSP_54)
      K = self%rk%nrk
   case default
      K = 0_I4P
      call mpih_fnl%error_stop(msg=': RK scheme "'//trim(self%rk%scheme)//'" is not stage-splittable: the staged '// &
                                   'forest path (multi-realm, or AMR seam faces) requires an SSP scheme')
   endselect
   endfunction stages_per_step_forest

   ! private procedures
   subroutine compute_residuals_weno_dev(self, q_gpu, dq_gpu, s, flux_register)
   !< Compute the residuals with the WENO space operator on the device: ghost update, auxiliary variables, face fluxes
   !< of the active directions, flux difference.
   !<
   !< The face fluxes of a null direction are never computed: they keep their zero initialization (`dev_alloc`). On the
   !< staged path (AMR seam faces), the seam face fluxes of every stage are accumulated into the forest's flux register.
   class(flume_fnl_object),     intent(inout)           :: self              !< The equation.
   real(R8P),                   intent(inout)           :: q_gpu(1:,         &
                                                                 1-self%ngc:,&
                                                                 1-self%ngc:,&
                                                                 1-self%ngc:,&
                                                                 1:)         !< Conservative variables.
   real(R8P),                   intent(inout)           :: dq_gpu(1:,         &
                                                                  1-self%ngc:,&
                                                                  1-self%ngc:,&
                                                                  1-self%ngc:,&
                                                                  1:)         !< Residuals.
   integer(I4P),                intent(in),    optional :: s                 !< Runge-Kutta stage.
   class(flux_register_object), intent(inout), optional :: flux_register     !< Forest's flux register for reflux.
   logical                                              :: is_char           !< Characteristic reconstruction flag.
   integer(I4P)                                         :: e                 !< Eikonal iterations counter.

   if (self%ib%solids_number > 0_I4P) then
      call self%update_ghost(q_gpu=q_gpu)
      do e=1, self%ib%n_eikonal
         call self%ib_fnl%evolve_eikonal(grid=self%adam%grid, field=self%adam%field, ib=self%ib, dq_gpu=dq_gpu, &
                                         q_gpu=q_gpu, dxyz_gpu=self%field_fnl%dxyz_gpu)
         call self%update_ghost(q_gpu=q_gpu)
      enddo
      call self%ib_fnl%invert_eikonal(grid=self%adam%grid, field=self%adam%field, ib=self%ib, q_gpu=q_gpu)
   endif
   call self%apply_floors(q_gpu=q_gpu)
   call self%update_ghost(q_gpu=q_gpu)
   call self%compute_q_aux(q_gpu=q_gpu)
   is_char = self%numerics%reconstruction_variables == RECON_CHARACTERISTIC
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nb=>self%blocks_number, gamma=>self%physics%gamma, &
             is_null=>self%adam%grid%null_xyz, weno_s=>self%weno%S, zeps=>self%weno%zeps, a_gpu=>self%weno_fnl%a_gpu, &
             sigma=>self%weno%sigma,                                                                                  &
             p_gpu=>self%weno_fnl%p_gpu, d_gpu=>self%weno_fnl%d_gpu)
   select case(self%physics%model)
   case(MODEL_EULER)
      if (.not.is_null(1)) call compute_face_fluxes_euler_dev(d=1_I4P, di=1_I4P, dj=0_I4P, dk=0_I4P, ni=ni, nj=nj,     &
                                                              nk=nk, ngc=ngc, blocks_number=nb, S=weno_s, gamma=gamma, &
                                                              ch=self%physics%mhd%glm_ch, is_characteristic=is_char,   &
                                                              weno_a_gpu=a_gpu, weno_p_gpu=p_gpu, weno_d_gpu=d_gpu,    &
                                                              weno_zeps=zeps, weno_sigma=sigma, q_gpu=q_gpu,           &
                                                              q_aux_gpu=self%q_aux_gpu, fl_gpu=self%flx_f_gpu)
      if (.not.is_null(2)) call compute_face_fluxes_euler_dev(d=2_I4P, di=0_I4P, dj=1_I4P, dk=0_I4P, ni=ni, nj=nj,     &
                                                              nk=nk, ngc=ngc, blocks_number=nb, S=weno_s, gamma=gamma, &
                                                              ch=self%physics%mhd%glm_ch, is_characteristic=is_char,   &
                                                              weno_a_gpu=a_gpu, weno_p_gpu=p_gpu, weno_d_gpu=d_gpu,    &
                                                              weno_zeps=zeps, weno_sigma=sigma, q_gpu=q_gpu,           &
                                                              q_aux_gpu=self%q_aux_gpu, fl_gpu=self%fly_f_gpu)
      if (.not.is_null(3)) call compute_face_fluxes_euler_dev(d=3_I4P, di=0_I4P, dj=0_I4P, dk=1_I4P, ni=ni, nj=nj,     &
                                                              nk=nk, ngc=ngc, blocks_number=nb, S=weno_s, gamma=gamma, &
                                                              ch=self%physics%mhd%glm_ch, is_characteristic=is_char,   &
                                                              weno_a_gpu=a_gpu, weno_p_gpu=p_gpu, weno_d_gpu=d_gpu,    &
                                                              weno_zeps=zeps, weno_sigma=sigma, q_gpu=q_gpu,           &
                                                              q_aux_gpu=self%q_aux_gpu, fl_gpu=self%flz_f_gpu)
   case(MODEL_MHD)
      if (.not.is_null(1)) call compute_face_fluxes_mhd_dev(d=1_I4P, di=1_I4P, dj=0_I4P, dk=0_I4P, ni=ni, nj=nj,     &
                                                            nk=nk, ngc=ngc, blocks_number=nb, S=weno_s, gamma=gamma, &
                                                            ch=self%physics%mhd%glm_ch, is_characteristic=is_char,   &
                                                            weno_a_gpu=a_gpu, weno_p_gpu=p_gpu, weno_d_gpu=d_gpu,    &
                                                            weno_zeps=zeps, weno_sigma=sigma, q_gpu=q_gpu,           &
                                                            q_aux_gpu=self%q_aux_gpu, fl_gpu=self%flx_f_gpu)
      if (.not.is_null(2)) call compute_face_fluxes_mhd_dev(d=2_I4P, di=0_I4P, dj=1_I4P, dk=0_I4P, ni=ni, nj=nj,     &
                                                            nk=nk, ngc=ngc, blocks_number=nb, S=weno_s, gamma=gamma, &
                                                            ch=self%physics%mhd%glm_ch, is_characteristic=is_char,   &
                                                            weno_a_gpu=a_gpu, weno_p_gpu=p_gpu, weno_d_gpu=d_gpu,    &
                                                            weno_zeps=zeps, weno_sigma=sigma, q_gpu=q_gpu,           &
                                                            q_aux_gpu=self%q_aux_gpu, fl_gpu=self%fly_f_gpu)
      if (.not.is_null(3)) call compute_face_fluxes_mhd_dev(d=3_I4P, di=0_I4P, dj=0_I4P, dk=1_I4P, ni=ni, nj=nj,     &
                                                            nk=nk, ngc=ngc, blocks_number=nb, S=weno_s, gamma=gamma, &
                                                            ch=self%physics%mhd%glm_ch, is_characteristic=is_char,   &
                                                            weno_a_gpu=a_gpu, weno_p_gpu=p_gpu, weno_d_gpu=d_gpu,    &
                                                            weno_zeps=zeps, weno_sigma=sigma, q_gpu=q_gpu,           &
                                                            q_aux_gpu=self%q_aux_gpu, fl_gpu=self%flz_f_gpu)
   case(MODEL_MHD_GLM)
      if (.not.is_null(1)) call compute_face_fluxes_mhd_glm_dev(d=1_I4P, di=1_I4P, dj=0_I4P, dk=0_I4P, ni=ni, nj=nj,     &
                                                                nk=nk, ngc=ngc, blocks_number=nb, S=weno_s, gamma=gamma, &
                                                                ch=self%physics%mhd%glm_ch, is_characteristic=is_char,   &
                                                                weno_a_gpu=a_gpu, weno_p_gpu=p_gpu, weno_d_gpu=d_gpu,    &
                                                                weno_zeps=zeps, weno_sigma=sigma, q_gpu=q_gpu,           &
                                                                q_aux_gpu=self%q_aux_gpu, fl_gpu=self%flx_f_gpu)
      if (.not.is_null(2)) call compute_face_fluxes_mhd_glm_dev(d=2_I4P, di=0_I4P, dj=1_I4P, dk=0_I4P, ni=ni, nj=nj,     &
                                                                nk=nk, ngc=ngc, blocks_number=nb, S=weno_s, gamma=gamma, &
                                                                ch=self%physics%mhd%glm_ch, is_characteristic=is_char,   &
                                                                weno_a_gpu=a_gpu, weno_p_gpu=p_gpu, weno_d_gpu=d_gpu,    &
                                                                weno_zeps=zeps, weno_sigma=sigma, q_gpu=q_gpu,           &
                                                                q_aux_gpu=self%q_aux_gpu, fl_gpu=self%fly_f_gpu)
      if (.not.is_null(3)) call compute_face_fluxes_mhd_glm_dev(d=3_I4P, di=0_I4P, dj=0_I4P, dk=1_I4P, ni=ni, nj=nj,     &
                                                                nk=nk, ngc=ngc, blocks_number=nb, S=weno_s, gamma=gamma, &
                                                                ch=self%physics%mhd%glm_ch, is_characteristic=is_char,   &
                                                                weno_a_gpu=a_gpu, weno_p_gpu=p_gpu, weno_d_gpu=d_gpu,    &
                                                                weno_zeps=zeps, weno_sigma=sigma, q_gpu=q_gpu,           &
                                                                q_aux_gpu=self%q_aux_gpu, fl_gpu=self%flz_f_gpu)
   case(MODEL_MHD_EGLM)
      if (.not.is_null(1)) call compute_face_fluxes_mhd_eglm_dev(d=1_I4P, di=1_I4P, dj=0_I4P, dk=0_I4P, ni=ni, nj=nj,     &
                                                                nk=nk, ngc=ngc, blocks_number=nb, S=weno_s, gamma=gamma, &
                                                                ch=self%physics%mhd%glm_ch, is_characteristic=is_char,   &
                                                                weno_a_gpu=a_gpu, weno_p_gpu=p_gpu, weno_d_gpu=d_gpu,    &
                                                                weno_zeps=zeps, weno_sigma=sigma, q_gpu=q_gpu,           &
                                                                q_aux_gpu=self%q_aux_gpu, fl_gpu=self%flx_f_gpu)
      if (.not.is_null(2)) call compute_face_fluxes_mhd_eglm_dev(d=2_I4P, di=0_I4P, dj=1_I4P, dk=0_I4P, ni=ni, nj=nj,     &
                                                                nk=nk, ngc=ngc, blocks_number=nb, S=weno_s, gamma=gamma, &
                                                                ch=self%physics%mhd%glm_ch, is_characteristic=is_char,   &
                                                                weno_a_gpu=a_gpu, weno_p_gpu=p_gpu, weno_d_gpu=d_gpu,    &
                                                                weno_zeps=zeps, weno_sigma=sigma, q_gpu=q_gpu,           &
                                                                q_aux_gpu=self%q_aux_gpu, fl_gpu=self%fly_f_gpu)
      if (.not.is_null(3)) call compute_face_fluxes_mhd_eglm_dev(d=3_I4P, di=0_I4P, dj=0_I4P, dk=1_I4P, ni=ni, nj=nj,     &
                                                                nk=nk, ngc=ngc, blocks_number=nb, S=weno_s, gamma=gamma, &
                                                                ch=self%physics%mhd%glm_ch, is_characteristic=is_char,   &
                                                                weno_a_gpu=a_gpu, weno_p_gpu=p_gpu, weno_d_gpu=d_gpu,    &
                                                                weno_zeps=zeps, weno_sigma=sigma, q_gpu=q_gpu,           &
                                                                q_aux_gpu=self%q_aux_gpu, fl_gpu=self%flz_f_gpu)
   case default
      call mpih_fnl%error_stop(msg=': no FNL kernels for physical model "'//self%physics%physical_model//'"')
   endselect
   if (self%numerics%positivity_limiter == POSITIVITY_LIMITER_CELL) call self%limit_positivity_dev(q_gpu=q_gpu, &
                                                                                             flux_register=flux_register)
   if (present(flux_register) .and. present(s) .and. self%numerics%reflux) then
      if (flux_register%nfaces > 0_I4P) call self%accumulate_seam_fluxes(s=s, flux_register=flux_register)
   endif
   if (self%ib%solids_number > 0_I4P) then
      call compute_flux_difference_ib_dev(nv=self%physics%nv, ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb,           &
                                          is_null=is_null, freeze=self%null_freeze(), dxyz_gpu=self%field_fnl%dxyz_gpu, &
                                          flx_f_gpu=self%flx_f_gpu,                                                      &
                                          fly_f_gpu=self%fly_f_gpu, flz_f_gpu=self%flz_f_gpu,                            &
                                          phi_gpu=self%ib_fnl%phi_gpu, dq_gpu=dq_gpu)
   else
      call compute_flux_difference_dev(nv=self%physics%nv, ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb,              &
                                       is_null=is_null, freeze=self%null_freeze(), dxyz_gpu=self%field_fnl%dxyz_gpu,    &
                                       flx_f_gpu=self%flx_f_gpu,                                                         &
                                       fly_f_gpu=self%fly_f_gpu, flz_f_gpu=self%flz_f_gpu, dq_gpu=dq_gpu)
   endif
   if (self%physics%model == MODEL_MHD_GLM) call add_glm_damping_dev(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb,      &
                                                                      damping=self%physics%mhd%glm_damping, q_gpu=q_gpu, &
                                                                      dq_gpu=dq_gpu)
   if (self%physics%model == MODEL_MHD_EGLM) then
      call add_glm_damping_eglm_dev(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb,                               &
                                    damping=self%physics%mhd%glm_damping, q_gpu=q_gpu, dq_gpu=dq_gpu)
      if (self%numerics%positivity_limiter == POSITIVITY_LIMITER_CELL) then
         call add_eglm_sources_limited_dev(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, hs=weno_s,              &
                                           dxyz_gpu=self%field_fnl%dxyz_gpu, is_null=is_null, q_gpu=q_gpu,         &
                                           q_aux_gpu=self%q_aux_gpu, lam_gpu=self%lam_gpu, dq_gpu=dq_gpu)
      else
         call add_eglm_sources_dev(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, hs=weno_s,                      &
                                   dxyz_gpu=self%field_fnl%dxyz_gpu, is_null=is_null, q_gpu=q_gpu,                 &
                                   q_aux_gpu=self%q_aux_gpu, dq_gpu=dq_gpu)
      endif
   endif
   endassociate
   endsubroutine compute_residuals_weno_dev

   subroutine compute_residuals_riemann_dev(self, q_gpu, dq_gpu, s, flux_register)
   !< Compute the residuals with the `weno-riemann` space operator on the device (issue #47): as
   !< `compute_residuals_weno_dev`, with the face fluxes of the WENO interpolation, the Riemann solver and the high-order
   !< correction. The host selects the kernel of the (model, solver) pair; a pair without kernels is fatal.
   class(flume_fnl_object),     intent(inout)           :: self              !< The equation.
   real(R8P),                   intent(inout)           :: q_gpu(1:,         &
                                                                 1-self%ngc:,&
                                                                 1-self%ngc:,&
                                                                 1-self%ngc:,&
                                                                 1:)         !< Conservative variables.
   real(R8P),                   intent(inout)           :: dq_gpu(1:,         &
                                                                  1-self%ngc:,&
                                                                  1-self%ngc:,&
                                                                  1-self%ngc:,&
                                                                  1:)         !< Residuals.
   integer(I4P),                intent(in),    optional :: s                 !< Runge-Kutta stage.
   class(flux_register_object), intent(inout), optional :: flux_register     !< Forest's flux register for reflux.
   logical                                              :: is_char           !< Characteristic interpolation flag.
   real(R8P)                                            :: cc(3)             !< Correction coefficients.
   real(R8P)                                            :: tau               !< Correction sensor threshold.
   integer(I4P)                                         :: e                 !< Eikonal iterations counter.
   procedure(compute_riemann_face_fluxes_euler_llf_dev), pointer :: face_fluxes !< Kernel of the (model, solver) pair.
   integer(I4P)                                         :: fallbacks(3)      !< Riemann solver fallbacks per direction.

   if (self%ib%solids_number > 0_I4P) then
      call self%update_ghost(q_gpu=q_gpu)
      do e=1, self%ib%n_eikonal
         call self%ib_fnl%evolve_eikonal(grid=self%adam%grid, field=self%adam%field, ib=self%ib, dq_gpu=dq_gpu, &
                                         q_gpu=q_gpu, dxyz_gpu=self%field_fnl%dxyz_gpu)
         call self%update_ghost(q_gpu=q_gpu)
      enddo
      call self%ib_fnl%invert_eikonal(grid=self%adam%grid, field=self%adam%field, ib=self%ib, q_gpu=q_gpu)
   endif
   call self%apply_floors(q_gpu=q_gpu)
   call self%update_ghost(q_gpu=q_gpu)
   call self%compute_q_aux(q_gpu=q_gpu)
   is_char = self%numerics%reconstruction_variables == RECON_CHARACTERISTIC
   cc = self%numerics%correction_coefficients()
   tau = self%numerics%correction_threshold()
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nb=>self%blocks_number, gamma=>self%physics%gamma, &
             ch=>self%physics%mhd%glm_ch, is_null=>self%adam%grid%null_xyz, weno_s=>self%weno%S,                  &
             zeps=>self%weno%zeps, a_gpu=>self%weno_fnl%a_interp_gpu, p_gpu=>self%weno_fnl%p_interp_gpu,           &
             sigma=>self%weno%sigma,                                                                               &
             d_gpu=>self%weno_fnl%d_gpu)
   face_fluxes => null()
   if (self%physics%model == MODEL_EULER) then
      select case(self%numerics%riemann_solver)
      case(RIEMANN_SOLVER_LLF)
         face_fluxes => compute_riemann_face_fluxes_euler_llf_dev
      case(RIEMANN_SOLVER_HLL)
         face_fluxes => compute_riemann_face_fluxes_euler_hll_dev
      case(RIEMANN_SOLVER_HLLC)
         face_fluxes => compute_riemann_face_fluxes_euler_hllc_dev
      endselect
   elseif (self%physics%model == MODEL_MHD) then
      select case(self%numerics%riemann_solver)
      case(RIEMANN_SOLVER_LLF)
         face_fluxes => compute_riemann_face_fluxes_mhd_llf_dev
      case(RIEMANN_SOLVER_HLL)
         face_fluxes => compute_riemann_face_fluxes_mhd_hll_dev
      case(RIEMANN_SOLVER_HLLD)
         face_fluxes => compute_riemann_face_fluxes_mhd_hlld_dev
      endselect
   elseif (self%physics%model == MODEL_MHD_GLM) then
      select case(self%numerics%riemann_solver)
      case(RIEMANN_SOLVER_LLF)
         face_fluxes => compute_riemann_face_fluxes_mhd_glm_llf_dev
      case(RIEMANN_SOLVER_HLL)
         face_fluxes => compute_riemann_face_fluxes_mhd_glm_hll_dev
      case(RIEMANN_SOLVER_HLLD)
         face_fluxes => compute_riemann_face_fluxes_mhd_glm_hlld_dev
      endselect
   elseif (self%physics%model == MODEL_MHD_EGLM) then
      select case(self%numerics%riemann_solver)
      case(RIEMANN_SOLVER_LLF)
         face_fluxes => compute_riemann_face_fluxes_mhd_eglm_llf_dev
      case(RIEMANN_SOLVER_HLL)
         face_fluxes => compute_riemann_face_fluxes_mhd_eglm_hll_dev
      case(RIEMANN_SOLVER_HLLD)
         face_fluxes => compute_riemann_face_fluxes_mhd_eglm_hlld_dev
      endselect
   endif
   if (.not.associated(face_fluxes)) call mpih_fnl%error_stop(msg=': no FNL weno-riemann kernels yet for '//      &
                                                                  '[physics].(physical_model)='//                 &
                                                                  self%physics%physical_model//                   &
                                                                  ' with [numerics].(riemann_solver)='//          &
                                                                  self%numerics%riemann_solver)
   fallbacks = 0_I4P
   if (.not.is_null(1)) call face_fluxes(d=1_I4P, di=1_I4P, dj=0_I4P, dk=0_I4P, ni=ni, nj=nj, nk=nk, ngc=ngc,         &
                                         blocks_number=nb, S=weno_s, gamma=gamma, ch=ch, is_characteristic=is_char,   &
                                         weno_a_gpu=a_gpu, weno_p_gpu=p_gpu, weno_d_gpu=d_gpu, weno_zeps=zeps,        &
                                         weno_sigma=sigma, cc=cc,                                                 &
                                         tau=tau, q_gpu=q_gpu, q_aux_gpu=self%q_aux_gpu, fl_gpu=self%flx_f_gpu,   &
                                         fallbacks=fallbacks(1))
   if (.not.is_null(2)) call face_fluxes(d=2_I4P, di=0_I4P, dj=1_I4P, dk=0_I4P, ni=ni, nj=nj, nk=nk, ngc=ngc,         &
                                         blocks_number=nb, S=weno_s, gamma=gamma, ch=ch, is_characteristic=is_char,   &
                                         weno_a_gpu=a_gpu, weno_p_gpu=p_gpu, weno_d_gpu=d_gpu, weno_zeps=zeps,        &
                                         weno_sigma=sigma, cc=cc,                                                 &
                                         tau=tau, q_gpu=q_gpu, q_aux_gpu=self%q_aux_gpu, fl_gpu=self%fly_f_gpu,   &
                                         fallbacks=fallbacks(2))
   if (.not.is_null(3)) call face_fluxes(d=3_I4P, di=0_I4P, dj=0_I4P, dk=1_I4P, ni=ni, nj=nj, nk=nk, ngc=ngc,         &
                                         blocks_number=nb, S=weno_s, gamma=gamma, ch=ch, is_characteristic=is_char,   &
                                         weno_a_gpu=a_gpu, weno_p_gpu=p_gpu, weno_d_gpu=d_gpu, weno_zeps=zeps,        &
                                         weno_sigma=sigma, cc=cc,                                                 &
                                         tau=tau, q_gpu=q_gpu, q_aux_gpu=self%q_aux_gpu, fl_gpu=self%flz_f_gpu,   &
                                         fallbacks=fallbacks(3))
   if (self%numerics%riemann_solver == RIEMANN_SOLVER_HLLD) then
      call MPI_ALLREDUCE(MPI_IN_PLACE, fallbacks, 3, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpih_fnl%error)
      if (sum(fallbacks) > 0_I4P .and. mpih_fnl%myrank == 0) &
         print '(A)', mpih_fnl%myrankstr//'HLLD fallbacks to HLL: '//trim(str(sum(fallbacks)))//' faces at step '// &
                      trim(str(self%time%it))
   endif
   if (self%numerics%positivity_limiter == POSITIVITY_LIMITER_CELL) call self%limit_positivity_dev(q_gpu=q_gpu, &
                                                                                             flux_register=flux_register)
   if (present(flux_register) .and. present(s) .and. self%numerics%reflux) then
      if (flux_register%nfaces > 0_I4P) call self%accumulate_seam_fluxes(s=s, flux_register=flux_register)
   endif
   if (self%ib%solids_number > 0_I4P) then
      call compute_flux_difference_ib_dev(nv=self%physics%nv, ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb,           &
                                          is_null=is_null, freeze=self%null_freeze(), dxyz_gpu=self%field_fnl%dxyz_gpu, &
                                          flx_f_gpu=self%flx_f_gpu,                                                      &
                                          fly_f_gpu=self%fly_f_gpu, flz_f_gpu=self%flz_f_gpu,                            &
                                          phi_gpu=self%ib_fnl%phi_gpu, dq_gpu=dq_gpu)
   else
      call compute_flux_difference_dev(nv=self%physics%nv, ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb,              &
                                       is_null=is_null, freeze=self%null_freeze(), dxyz_gpu=self%field_fnl%dxyz_gpu,    &
                                       flx_f_gpu=self%flx_f_gpu,                                                         &
                                       fly_f_gpu=self%fly_f_gpu, flz_f_gpu=self%flz_f_gpu, dq_gpu=dq_gpu)
   endif
   if (self%physics%model == MODEL_MHD_GLM) call add_glm_damping_dev(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb,      &
                                                                      damping=self%physics%mhd%glm_damping, q_gpu=q_gpu, &
                                                                      dq_gpu=dq_gpu)
   if (self%physics%model == MODEL_MHD_EGLM) then
      call add_glm_damping_eglm_dev(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb,                               &
                                    damping=self%physics%mhd%glm_damping, q_gpu=q_gpu, dq_gpu=dq_gpu)
      if (self%numerics%positivity_limiter == POSITIVITY_LIMITER_CELL) then
         call add_eglm_sources_limited_dev(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, hs=weno_s,              &
                                           dxyz_gpu=self%field_fnl%dxyz_gpu, is_null=is_null, q_gpu=q_gpu,         &
                                           q_aux_gpu=self%q_aux_gpu, lam_gpu=self%lam_gpu, dq_gpu=dq_gpu)
      else
         call add_eglm_sources_dev(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, hs=weno_s,                      &
                                   dxyz_gpu=self%field_fnl%dxyz_gpu, is_null=is_null, q_gpu=q_gpu,                 &
                                   q_aux_gpu=self%q_aux_gpu, dq_gpu=dq_gpu)
      endif
   endif
   endassociate
   endsubroutine compute_residuals_riemann_dev

   subroutine integrate_rk_ls_dev(self)
   !< Integrate one time step with a low-storage Runge-Kutta scheme on the device.
   class(flume_fnl_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: s    !< Counter.

   call self%rk_fnl%initialize_stages(grid=self%adam%grid, field=self%adam%field, q_gpu=self%q_gpu)
   do s=1, self%rk%nrk
      call self%compute_residuals_dev(q_gpu=self%q_gpu, dq_gpu=self%dq_gpu, s=s)
      if (s == 1) call self%save_residuals
      call rk_compute_stage_ls(self, s=s)
   enddo
   endsubroutine integrate_rk_ls_dev

   subroutine integrate_rk_ssp_dev(self)
   !< Integrate one time step with a strong stability preserving Runge-Kutta scheme on the device.
   class(flume_fnl_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: s    !< Counter.

   call self%rk_fnl%initialize_stages(grid=self%adam%grid, field=self%adam%field, q_gpu=self%q_gpu)
   do s=1, self%rk%nrk
      call rk_compute_stage(self, s=s)
      call self%compute_residuals_dev(q_gpu=self%rk_fnl%q_rk_gpu(:,:,:,:,:,s), dq_gpu=self%dq_gpu, s=s)
      call rk_assign_stage(self, s=s)
   enddo
   call compute_rk_ssp_residual_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, nv=self%nv,                     &
                                    blocks_number=self%blocks_number, nrk=self%rk%nrk, beta_gpu=self%rk_fnl%beta_gpu, &
                                    q_rk_gpu=self%rk_fnl%q_rk_gpu, dq_gpu=self%dq_gpu)
   call rk_update_q(self)
   call self%save_residuals
   endsubroutine integrate_rk_ssp_dev

   subroutine rk_assign_stage(self, s)
   !< Assign the residual of stage `s` to the Runge-Kutta stage buffer; the solid cells are masked with immersed solids.
   class(flume_fnl_object), intent(inout) :: self !< The equation.
   integer(I4P),            intent(in)    :: s    !< Stage.

   if (associated(self%ib_fnl%phi_gpu)) then
      call self%rk_fnl%assign_stage(grid=self%adam%grid, field=self%adam%field, s=s, q_gpu=self%dq_gpu, &
                                    phi_gpu=self%ib_fnl%phi_gpu)
   else
      call self%rk_fnl%assign_stage(grid=self%adam%grid, field=self%adam%field, s=s, q_gpu=self%dq_gpu)
   endif
   endsubroutine rk_assign_stage

   subroutine rk_compute_stage(self, s)
   !< Compute the state of stage `s`; the solid cells are masked with immersed solids.
   class(flume_fnl_object), intent(inout) :: self !< The equation.
   integer(I4P),            intent(in)    :: s    !< Stage.

   if (associated(self%ib_fnl%phi_gpu)) then
      call self%rk_fnl%compute_stage(grid=self%adam%grid, field=self%adam%field, s=s, dt=self%time%dt, &
                                     phi_gpu=self%ib_fnl%phi_gpu)
   else
      call self%rk_fnl%compute_stage(grid=self%adam%grid, field=self%adam%field, s=s, dt=self%time%dt)
   endif
   endsubroutine rk_compute_stage

   subroutine rk_compute_stage_ls(self, s)
   !< Advance low-storage stage `s`; the solid cells are masked with immersed solids.
   class(flume_fnl_object), intent(inout) :: self !< The equation.
   integer(I4P),            intent(in)    :: s    !< Stage.

   if (associated(self%ib_fnl%phi_gpu)) then
      call self%rk_fnl%compute_stage_ls(grid=self%adam%grid, field=self%adam%field, rk=self%rk, s=s, dt=self%time%dt, &
                                        phi_gpu=self%ib_fnl%phi_gpu, dq_gpu=self%dq_gpu, q_gpu=self%q_gpu)
   else
      call self%rk_fnl%compute_stage_ls(grid=self%adam%grid, field=self%adam%field, rk=self%rk, s=s, dt=self%time%dt, &
                                        dq_gpu=self%dq_gpu, q_gpu=self%q_gpu)
   endif
   endsubroutine rk_compute_stage_ls

   subroutine rk_update_q(self)
   !< Assemble the committed state of a strong stability preserving step; the solid cells are masked with immersed solids.
   class(flume_fnl_object), intent(inout) :: self !< The equation.

   if (associated(self%ib_fnl%phi_gpu)) then
      call self%rk_fnl%update_q(grid=self%adam%grid, field=self%adam%field, rk=self%rk, dt=self%time%dt, &
                                phi_gpu=self%ib_fnl%phi_gpu, q_gpu=self%q_gpu)
   else
      call self%rk_fnl%update_q(grid=self%adam%grid, field=self%adam%field, rk=self%rk, dt=self%time%dt, q_gpu=self%q_gpu)
   endif
   endsubroutine rk_update_q

   subroutine seam_gather_cells(cells, nv, ngc, a_gpu, out)
   !< Gather on the host `out(v, m)` = the first `nv` components of the device field `a_gpu` at `cells(:, m)` (issue #50).
   integer(I4P), intent(in)    :: cells(1:,1:)                      !< Cells (i, j, k, b) [4, n].
   integer(I4P), intent(in)    :: nv                                !< Components.
   integer(I4P), intent(in)    :: ngc                               !< Ghost cells number.
   real(R8P),    intent(in)    :: a_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Device field.
   real(R8P),    intent(inout) :: out(1:,1:)                        !< Gathered values [nv, n].
   integer(I4P), pointer       :: cells_gpu(:,:)                    !< Device cells.
   real(R8P),    pointer       :: out_gpu(:,:)                      !< Device values.
   integer(I4P)                :: n, ierr                           !< Cells number, error status.

   n = size(cells, dim=2)
   if (n == 0_I4P) return
   call dev_alloc(fptr_dev=cells_gpu, lbounds=[1,1], ubounds=[4,n], ierr=ierr)
   call dev_alloc(fptr_dev=out_gpu, lbounds=[1,1], ubounds=[nv,n], ierr=ierr)
   call dev_memcpy_to_device(dst=cells_gpu, src=cells)
   call gather_seam_cells_dev(n=n, nv=nv, ngc=ngc, cells_gpu=cells_gpu, a_gpu=a_gpu, out_gpu=out_gpu)
   call dev_memcpy_from_device(dst=out, src=out_gpu)
   call dev_free(cells_gpu, mydev)
   call dev_free(out_gpu, mydev)
   endsubroutine seam_gather_cells

   subroutine seam_gather_stencils(cells, nv, ngc, s, q_gpu, out)
   !< Gather on the host the axis stencils of half width `s` around `cells(:, m)` (issue #50).
   integer(I4P), intent(in)    :: cells(1:,1:)                      !< Cells (i, j, k, b) [4, n].
   integer(I4P), intent(in)    :: nv                                !< Variables number.
   integer(I4P), intent(in)    :: ngc                               !< Ghost cells number.
   integer(I4P), intent(in)    :: s                                 !< Stencil half width.
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Device conservative variables.
   real(R8P),    intent(inout) :: out(1:,1:,1:,1:)                  !< Stencils [nv, 2 s + 1, 3, n].
   integer(I4P), pointer       :: cells_gpu(:,:)                    !< Device cells.
   real(R8P),    pointer       :: out_gpu(:,:,:,:)                  !< Device stencils.
   integer(I4P)                :: n, ierr                           !< Cells number, error status.

   n = size(cells, dim=2)
   if (n == 0_I4P) return
   call dev_alloc(fptr_dev=cells_gpu, lbounds=[1,1], ubounds=[4,n], ierr=ierr)
   call dev_alloc(fptr_dev=out_gpu, lbounds=[1,1,1,1], ubounds=[nv,2*s+1,3,n], ierr=ierr)
   call dev_memcpy_to_device(dst=cells_gpu, src=cells)
   call gather_seam_stencils_dev(n=n, nv=nv, ngc=ngc, s=s, cells_gpu=cells_gpu, q_gpu=q_gpu, out_gpu=out_gpu)
   call dev_memcpy_from_device(dst=out, src=out_gpu)
   call dev_free(cells_gpu, mydev)
   call dev_free(out_gpu, mydev)
   endsubroutine seam_gather_stencils

   subroutine seam_gather_faces(faces, nv, flx_f_gpu, fly_f_gpu, flz_f_gpu, out)
   !< Gather on the host the device face fluxes at `faces(:, m) = (axis, i, j, k, b)` (issue #50).
   integer(I4P), intent(in)    :: faces(1:,1:)              !< Faces [5, n].
   integer(I4P), intent(in)    :: nv                        !< Variables number.
   real(R8P),    intent(in)    :: flx_f_gpu(1:,0:,1:,1:,1:) !< X-face fluxes.
   real(R8P),    intent(in)    :: fly_f_gpu(1:,1:,0:,1:,1:) !< Y-face fluxes.
   real(R8P),    intent(in)    :: flz_f_gpu(1:,1:,1:,0:,1:) !< Z-face fluxes.
   real(R8P),    intent(inout) :: out(1:,1:)                !< Gathered fluxes [nv, n].
   integer(I4P), pointer       :: faces_gpu(:,:)            !< Device faces.
   real(R8P),    pointer       :: out_gpu(:,:)              !< Device fluxes.
   integer(I4P)                :: n, ierr                   !< Faces number, error status.

   n = size(faces, dim=2)
   if (n == 0_I4P) return
   call dev_alloc(fptr_dev=faces_gpu, lbounds=[1,1], ubounds=[5,n], ierr=ierr)
   call dev_alloc(fptr_dev=out_gpu, lbounds=[1,1], ubounds=[nv,n], ierr=ierr)
   call dev_memcpy_to_device(dst=faces_gpu, src=faces)
   call gather_seam_faces_dev(n=n, nv=nv, faces_gpu=faces_gpu, flx_f_gpu=flx_f_gpu, fly_f_gpu=fly_f_gpu, &
                              flz_f_gpu=flz_f_gpu, out_gpu=out_gpu)
   call dev_memcpy_from_device(dst=out, src=out_gpu)
   call dev_free(faces_gpu, mydev)
   call dev_free(out_gpu, mydev)
   endsubroutine seam_gather_faces

   subroutine seam_scatter_faces(faces, nv, in, flx_f_gpu, fly_f_gpu, flz_f_gpu)
   !< Set the device face fluxes at `faces(:, m) = (axis, i, j, k, b)` to the host `in(:, m)` (issue #50).
   integer(I4P), intent(in)    :: faces(1:,1:)              !< Faces [5, n].
   integer(I4P), intent(in)    :: nv                        !< Variables number.
   real(R8P),    intent(in)    :: in(1:,1:)                 !< Fluxes [nv, n].
   real(R8P),    intent(inout) :: flx_f_gpu(1:,0:,1:,1:,1:) !< X-face fluxes.
   real(R8P),    intent(inout) :: fly_f_gpu(1:,1:,0:,1:,1:) !< Y-face fluxes.
   real(R8P),    intent(inout) :: flz_f_gpu(1:,1:,1:,0:,1:) !< Z-face fluxes.
   integer(I4P), pointer       :: faces_gpu(:,:)            !< Device faces.
   real(R8P),    pointer       :: in_gpu(:,:)               !< Device fluxes.
   integer(I4P)                :: n, ierr                   !< Faces number, error status.

   n = size(faces, dim=2)
   if (n == 0_I4P) return
   call dev_alloc(fptr_dev=faces_gpu, lbounds=[1,1], ubounds=[5,n], ierr=ierr)
   call dev_alloc(fptr_dev=in_gpu, lbounds=[1,1], ubounds=[nv,n], ierr=ierr)
   call dev_memcpy_to_device(dst=faces_gpu, src=faces)
   call dev_memcpy_to_device(dst=in_gpu, src=in)
   call scatter_seam_faces_dev(n=n, nv=nv, faces_gpu=faces_gpu, in_gpu=in_gpu, flx_f_gpu=flx_f_gpu, &
                               fly_f_gpu=fly_f_gpu, flz_f_gpu=flz_f_gpu)
   call dev_free(faces_gpu, mydev)
   call dev_free(in_gpu, mydev)
   endsubroutine seam_scatter_faces

   subroutine seam_scatter_cells(cells, nv, ngc, in, a_gpu)
   !< Set the first `nv` components of the device field `a_gpu` at `cells(:, m)` to the host `in(:, m)` (issue #50).
   integer(I4P), intent(in)    :: cells(1:,1:)                      !< Cells (i, j, k, b) [4, n].
   integer(I4P), intent(in)    :: nv                                !< Components.
   integer(I4P), intent(in)    :: ngc                               !< Ghost cells number.
   real(R8P),    intent(in)    :: in(1:,1:)                         !< Values [nv, n].
   real(R8P),    intent(inout) :: a_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Device field.
   integer(I4P), pointer       :: cells_gpu(:,:)                    !< Device cells.
   real(R8P),    pointer       :: in_gpu(:,:)                       !< Device values.
   integer(I4P)                :: n, ierr                           !< Cells number, error status.

   n = size(cells, dim=2)
   if (n == 0_I4P) return
   call dev_alloc(fptr_dev=cells_gpu, lbounds=[1,1], ubounds=[4,n], ierr=ierr)
   call dev_alloc(fptr_dev=in_gpu, lbounds=[1,1], ubounds=[nv,n], ierr=ierr)
   call dev_memcpy_to_device(dst=cells_gpu, src=cells)
   call dev_memcpy_to_device(dst=in_gpu, src=in)
   call scatter_seam_cells_dev(n=n, nv=nv, ngc=ngc, cells_gpu=cells_gpu, in_gpu=in_gpu, a_gpu=a_gpu)
   call dev_free(cells_gpu, mydev)
   call dev_free(in_gpu, mydev)
   endsubroutine seam_scatter_cells
endmodule adam_flume_fnl_object
