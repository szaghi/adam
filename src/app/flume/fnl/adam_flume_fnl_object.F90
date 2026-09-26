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
use :: adam_flux_register_object, only : flux_register_object
use :: adam_maps_object,          only : face_axis_sign
use :: adam_realm_object,         only : realm_object
use :: adam_rk_object,            only : RK_1, RK_2, RK_3, RK_SSP_11, RK_SSP_22, RK_SSP_33, RK_SSP_54
! ADAM FNL classes, libraries
use :: adam_fnl_field_kernels,    only : compute_normL2_residuals_dev
use :: adam_fnl_field_object,     only : field_fnl_object
use :: adam_fnl_ib_object,        only : ib_fnl_object
use :: adam_fnl_rk_object,        only : rk_fnl_object
use :: adam_fnl_weno_object,      only : weno_fnl_object
! ADAM singleton objects
use :: adam_fnl_mpih_global,      only : mpih_fnl, mpih_fnl_is_initialized
! FLUME modules
use :: adam_flume_common_library,      only : flume_common_object, MODEL_EULER, MODEL_MHD, MODEL_MHD_GLM,                &
                                              RECON_CHARACTERISTIC, SCHEME_SPACE_WENO
use :: adam_flume_fnl_euler_kernels,   only : compute_conservation_euler_dev=>compute_conservation_dev,                &
                                              compute_face_fluxes_euler_dev=>compute_face_fluxes_dev,                  &
                                              compute_lambda_max_euler_dev=>compute_lambda_max_dev,                    &
                                              compute_q_aux_euler_dev=>compute_q_aux_dev
use :: adam_flume_fnl_mhd_kernels,     only : compute_conservation_mhd_dev=>compute_conservation_dev,                 &
                                              compute_face_fluxes_mhd_dev=>compute_face_fluxes_dev,                   &
                                              compute_lambda_max_mhd_dev=>compute_lambda_max_dev,                     &
                                              compute_q_aux_mhd_dev=>compute_q_aux_dev
use :: adam_flume_fnl_mhd_glm_kernels, only : compute_conservation_mhd_glm_dev=>compute_conservation_dev,             &
                                              compute_face_fluxes_mhd_glm_dev=>compute_face_fluxes_dev,               &
                                              compute_lambda_max_mhd_glm_dev=>compute_lambda_max_dev,                 &
                                              compute_q_aux_mhd_glm_dev=>compute_q_aux_dev
use :: adam_flume_fnl_kernels,         only : apply_reflux_face_dev,                                                   &
                                              compute_flux_difference_dev, compute_flux_difference_ib_dev,             &
                                              compute_rk_ssp_residual_dev, fill_seam_copy_dev, pack_seam_skin_dev,     &
                                              set_boundary_conditions_dev
! third party modules
use :: fundal,                    only : dev_alloc, dev_assign_to_device, dev_free, dev_memcpy_from_device,     &
                                         dev_memcpy_to_device, mydev
use :: mpi
use :: penf,                      only : I4P, R8P, str

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
   real(R8P), pointer     :: q_inflow_gpu(:,:)=>null()    !< Conservative inflow state of each face [nv, 6].
   real(R8P), pointer     :: wall_sign_gpu(:,:)=>null()   !< Wall mirror sign per variable and direction [nv, 3].
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
      procedure, pass(self) :: allocate_gpu            !< Allocate device data.
      procedure, pass(self) :: compute_conservation    !< Compute and save the conservation integrals.
      procedure, pass(self) :: compute_q_aux           !< Compute the auxiliary variables.
      procedure, pass(self) :: copy_cpu_gpu            !< Copy state and topology from host to device.
      procedure, pass(self) :: copy_gpu_cpu            !< Copy state from device to host.
      procedure, pass(self) :: copy_phi_gpu            !< Copy the immersed solids distance function to the device.
      procedure, pass(self) :: destroy                 !< Free device and host data.
      procedure, pass(self) :: initialize_flume        !< Initialize the FNL backend.
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
   case default
      call mpih_fnl%error_stop(msg=': no FNL kernels for physical model "'//self%physics%physical_model//'"')
   endselect
   call MPI_ALLREDUCE(MPI_IN_PLACE, integrals, size(integrals), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih_fnl%error)
   call self%diagnostics%save_conservation_row(it=self%time%it, time=self%time%time, integrals=integrals)
   endsubroutine compute_conservation

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
   if (associated(self%q_inflow_gpu)) then
      call dev_free(self%q_inflow_gpu, mydev)
      nullify(self%q_inflow_gpu)
   endif
   if (associated(self%wall_sign_gpu)) then
      call dev_free(self%wall_sign_gpu, mydev)
      nullify(self%wall_sign_gpu)
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
   call self%field_fnl%initialize(grid=self%adam%grid, field=self%adam%field, maps=self%adam%maps, verbose=.true.)
   call self%ib_fnl%initialize(grid=self%adam%grid, field=self%adam%field, ib=self%ib)
   call self%rk_fnl%initialize(grid=self%adam%grid, field=self%adam%field, rk=self%rk)
   call self%weno_fnl%initialize(weno=self%weno)
   call self%allocate_gpu
   select case(self%numerics%scheme_space)
   case(SCHEME_SPACE_WENO)
      self%compute_residuals_dev => compute_residuals_weno_dev
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
   if (mpih_fnl%myrank == 0) call self%io%save_residuals(it=self%time%it, time=self%time%time, &
                                                         blocks_number=self%blocks_number,     &
                                                         residuals=self%adam%field%residuals)
   endsubroutine save_residuals

   subroutine save_simulation_data(self)
   !< Save fields, restart, slices and conservation history, each on its own cadence; state copied to host only when
   !< saved.
   class(flume_fnl_object), intent(inout) :: self      !< The equation.
   logical                                :: is_slices !< Slices save step.

   is_slices = self%slices%is_to_save(it=self%time%it, it_max=self%time%it_max, time=self%time%time, &
                                      time_max=self%time%time_max)
   if (self%time%is_to_save(cadence=self%io%it_save) .or. self%time%is_to_save(cadence=self%io%restart_save) .or. &
       is_slices) then
      call self%update_ghost(q_gpu=self%q_gpu)
      call self%copy_gpu_cpu
      if (self%time%is_to_save(cadence=self%io%it_save)) call self%save_xh5f(with_ghost=.true.)
      if (self%time%is_to_save(cadence=self%io%restart_save)) call self%save_restart_files
      if (is_slices) call self%save_slices
   endif
   call self%compute_conservation
   endsubroutine save_simulation_data

   subroutine set_boundary_conditions(self, q_gpu)
   !< Set boundary conditions on the device crown maps, crown by crown.
   class(flume_fnl_object), intent(inout) :: self              !< The equation.
   real(R8P),               intent(inout) :: q_gpu(1:,         &
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1:)         !< Conservative variables.
   integer(I4P)                           :: crown             !< Crown counter.

   if (.not.associated(self%field_fnl%maps%local_map_bc_crown_gpu)) return
   do crown=1, self%ngc
      call set_boundary_conditions_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, nv=self%nv, crown=crown, &
                                       local_map_bc_crown_gpu=self%field_fnl%maps%local_map_bc_crown_gpu,        &
                                       q_inflow_gpu=self%q_inflow_gpu, wall_sign_gpu=self%wall_sign_gpu, q_gpu=q_gpu)
   enddo
   endsubroutine set_boundary_conditions

   subroutine update_ghost(self, q_gpu)
   !< Update ghost cells on the device: intra-realm local copies, GPU-direct MPI exchange, boundary conditions.
   class(flume_fnl_object), intent(inout) :: self              !< The equation.
   real(R8P),               intent(inout) :: q_gpu(1:,         &
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1:)         !< Conservative variables.

   call self%field_fnl%update_ghost_local_gpu(q_gpu=q_gpu)
   call self%field_fnl%update_ghost_mpi_gpu(comm_map_send_ptr_ghost=self%adam%maps%comm_map_send_ptr_ghost, &
                                            comm_map_recv_ptr_ghost=self%adam%maps%comm_map_recv_ptr_ghost, &
                                            q_gpu=q_gpu)
   call self%set_boundary_conditions(q_gpu=q_gpu)
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
   !< Compute the local stability-limited time step on the device, `dt = CFL / max(sum_d (|u_d| + a) / dx_d)`.
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
   call self%copy_cpu_gpu(verbose=.true.)
   call self%copy_phi_gpu
   call self%update_ghost(q_gpu=self%q_gpu)
   call self%compute_q_aux(q_gpu=self%q_gpu)
   call self%diagnostics%open_file(output_basename=self%io%output_basename, q_name=self%q_name, &
                                   is_restart=self%io%restart)
   ! a restarted run starts from a step its predecessor already saved: saving it again would duplicate the history rows
   if (.not.self%io%restart) call self%save_simulation_data
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
   !< Post-step work: fields, restart and conservation history on their cadence.
   class(flume_fnl_object), intent(inout)                   :: self              !< The equation.
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
   call self%update_ghost(q_gpu=q_gpu)
   call self%compute_q_aux(q_gpu=q_gpu)
   is_char = self%numerics%reconstruction_variables == RECON_CHARACTERISTIC
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nb=>self%blocks_number, gamma=>self%physics%gamma, &
             is_null=>self%adam%grid%null_xyz, weno_s=>self%weno%S, zeps=>self%weno%zeps, a_gpu=>self%weno_fnl%a_gpu, &
             p_gpu=>self%weno_fnl%p_gpu, d_gpu=>self%weno_fnl%d_gpu)
   select case(self%physics%model)
   case(MODEL_EULER)
      if (.not.is_null(1)) call compute_face_fluxes_euler_dev(d=1_I4P, di=1_I4P, dj=0_I4P, dk=0_I4P, ni=ni, nj=nj,     &
                                                              nk=nk, ngc=ngc, blocks_number=nb, S=weno_s, gamma=gamma, &
                                                              ch=self%physics%mhd%glm_ch, is_characteristic=is_char,   &
                                                              weno_a_gpu=a_gpu, weno_p_gpu=p_gpu, weno_d_gpu=d_gpu,    &
                                                              weno_zeps=zeps, q_gpu=q_gpu, q_aux_gpu=self%q_aux_gpu,   &
                                                              fl_gpu=self%flx_f_gpu)
      if (.not.is_null(2)) call compute_face_fluxes_euler_dev(d=2_I4P, di=0_I4P, dj=1_I4P, dk=0_I4P, ni=ni, nj=nj,     &
                                                              nk=nk, ngc=ngc, blocks_number=nb, S=weno_s, gamma=gamma, &
                                                              ch=self%physics%mhd%glm_ch, is_characteristic=is_char,   &
                                                              weno_a_gpu=a_gpu, weno_p_gpu=p_gpu, weno_d_gpu=d_gpu,    &
                                                              weno_zeps=zeps, q_gpu=q_gpu, q_aux_gpu=self%q_aux_gpu,   &
                                                              fl_gpu=self%fly_f_gpu)
      if (.not.is_null(3)) call compute_face_fluxes_euler_dev(d=3_I4P, di=0_I4P, dj=0_I4P, dk=1_I4P, ni=ni, nj=nj,     &
                                                              nk=nk, ngc=ngc, blocks_number=nb, S=weno_s, gamma=gamma, &
                                                              ch=self%physics%mhd%glm_ch, is_characteristic=is_char,   &
                                                              weno_a_gpu=a_gpu, weno_p_gpu=p_gpu, weno_d_gpu=d_gpu,    &
                                                              weno_zeps=zeps, q_gpu=q_gpu, q_aux_gpu=self%q_aux_gpu,   &
                                                              fl_gpu=self%flz_f_gpu)
   case(MODEL_MHD)
      if (.not.is_null(1)) call compute_face_fluxes_mhd_dev(d=1_I4P, di=1_I4P, dj=0_I4P, dk=0_I4P, ni=ni, nj=nj,     &
                                                            nk=nk, ngc=ngc, blocks_number=nb, S=weno_s, gamma=gamma, &
                                                            ch=self%physics%mhd%glm_ch, is_characteristic=is_char,   &
                                                            weno_a_gpu=a_gpu, weno_p_gpu=p_gpu, weno_d_gpu=d_gpu,    &
                                                            weno_zeps=zeps, q_gpu=q_gpu, q_aux_gpu=self%q_aux_gpu,   &
                                                            fl_gpu=self%flx_f_gpu)
      if (.not.is_null(2)) call compute_face_fluxes_mhd_dev(d=2_I4P, di=0_I4P, dj=1_I4P, dk=0_I4P, ni=ni, nj=nj,     &
                                                            nk=nk, ngc=ngc, blocks_number=nb, S=weno_s, gamma=gamma, &
                                                            ch=self%physics%mhd%glm_ch, is_characteristic=is_char,   &
                                                            weno_a_gpu=a_gpu, weno_p_gpu=p_gpu, weno_d_gpu=d_gpu,    &
                                                            weno_zeps=zeps, q_gpu=q_gpu, q_aux_gpu=self%q_aux_gpu,   &
                                                            fl_gpu=self%fly_f_gpu)
      if (.not.is_null(3)) call compute_face_fluxes_mhd_dev(d=3_I4P, di=0_I4P, dj=0_I4P, dk=1_I4P, ni=ni, nj=nj,     &
                                                            nk=nk, ngc=ngc, blocks_number=nb, S=weno_s, gamma=gamma, &
                                                            ch=self%physics%mhd%glm_ch, is_characteristic=is_char,   &
                                                            weno_a_gpu=a_gpu, weno_p_gpu=p_gpu, weno_d_gpu=d_gpu,    &
                                                            weno_zeps=zeps, q_gpu=q_gpu, q_aux_gpu=self%q_aux_gpu,   &
                                                            fl_gpu=self%flz_f_gpu)
   case(MODEL_MHD_GLM)
      if (.not.is_null(1)) call compute_face_fluxes_mhd_glm_dev(d=1_I4P, di=1_I4P, dj=0_I4P, dk=0_I4P, ni=ni, nj=nj,     &
                                                                nk=nk, ngc=ngc, blocks_number=nb, S=weno_s, gamma=gamma, &
                                                                ch=self%physics%mhd%glm_ch, is_characteristic=is_char,   &
                                                                weno_a_gpu=a_gpu, weno_p_gpu=p_gpu, weno_d_gpu=d_gpu,    &
                                                                weno_zeps=zeps, q_gpu=q_gpu, q_aux_gpu=self%q_aux_gpu,   &
                                                                fl_gpu=self%flx_f_gpu)
      if (.not.is_null(2)) call compute_face_fluxes_mhd_glm_dev(d=2_I4P, di=0_I4P, dj=1_I4P, dk=0_I4P, ni=ni, nj=nj,     &
                                                                nk=nk, ngc=ngc, blocks_number=nb, S=weno_s, gamma=gamma, &
                                                                ch=self%physics%mhd%glm_ch, is_characteristic=is_char,   &
                                                                weno_a_gpu=a_gpu, weno_p_gpu=p_gpu, weno_d_gpu=d_gpu,    &
                                                                weno_zeps=zeps, q_gpu=q_gpu, q_aux_gpu=self%q_aux_gpu,   &
                                                                fl_gpu=self%fly_f_gpu)
      if (.not.is_null(3)) call compute_face_fluxes_mhd_glm_dev(d=3_I4P, di=0_I4P, dj=0_I4P, dk=1_I4P, ni=ni, nj=nj,     &
                                                                nk=nk, ngc=ngc, blocks_number=nb, S=weno_s, gamma=gamma, &
                                                                ch=self%physics%mhd%glm_ch, is_characteristic=is_char,   &
                                                                weno_a_gpu=a_gpu, weno_p_gpu=p_gpu, weno_d_gpu=d_gpu,    &
                                                                weno_zeps=zeps, q_gpu=q_gpu, q_aux_gpu=self%q_aux_gpu,   &
                                                                fl_gpu=self%flz_f_gpu)
   case default
      call mpih_fnl%error_stop(msg=': no FNL kernels for physical model "'//self%physics%physical_model//'"')
   endselect
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
   endassociate
   endsubroutine compute_residuals_weno_dev

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
endmodule adam_flume_fnl_object
