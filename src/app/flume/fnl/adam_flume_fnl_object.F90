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
use :: adam_flume_common_library, only : flume_common_object, RECON_CHARACTERISTIC, SCHEME_SPACE_WENO
use :: adam_flume_fnl_kernels,    only : compute_conservation_dev, compute_face_fluxes_dev, compute_flux_difference_dev, &
                                         compute_lambda_max_dev, compute_q_aux_dev, compute_rk_ssp_residual_dev,       &
                                         fill_seam_copy_dev, set_boundary_conditions_dev
! third party modules
use :: fundal,                    only : dev_alloc, dev_free, dev_memcpy_from_device, dev_memcpy_to_device, mydev
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
   ! host staging
   real(R8P), allocatable :: buf_5D_R8P(:,:,:,:,:)        !< Transposed copy buffer, extent identical to q_gpu.
   integer(I4P)           :: db5(2,5)=0_I4P               !< Device bounds of the transposed copies.
   integer(I4P)           :: hb5(2,5)=0_I4P               !< Host bounds of the transposed copies.
   ! dispatch
   procedure(compute_residuals_dev_interface), pass(self), pointer :: compute_residuals_dev=>null() !< Space operator.
   procedure(integrate_dev_interface),         pass(self), pointer :: integrate_dev=>null()         !< Time operator.
   contains
      ! public methods
      procedure, pass(self) :: allocate_gpu            !< Allocate device data.
      procedure, pass(self) :: compute_conservation    !< Compute and save the conservation integrals.
      procedure, pass(self) :: compute_q_aux           !< Compute the auxiliary variables.
      procedure, pass(self) :: copy_cpu_gpu            !< Copy state and topology from host to device.
      procedure, pass(self) :: copy_gpu_cpu            !< Copy state from device to host.
      procedure, pass(self) :: destroy                 !< Free device and host data.
      procedure, pass(self) :: initialize_flume        !< Initialize the FNL backend.
      procedure, pass(self) :: save_residuals          !< Save residuals history.
      procedure, pass(self) :: save_simulation_data    !< Save fields, restart and diagnostics on their cadence.
      procedure, pass(self) :: set_boundary_conditions !< Set boundary conditions on the device crown maps.
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
   real(R8P)                              :: integrals(5) !< Volume integrals.

   if (.not.self%time%is_to_save(cadence=self%diagnostics%conservation_history_save)) return
   call compute_conservation_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                                 dxyz_gpu=self%field_fnl%dxyz_gpu, is_null=self%adam%grid%null_xyz, q_gpu=self%q_gpu, &
                                 integrals=integrals)
   call MPI_ALLREDUCE(MPI_IN_PLACE, integrals, 5, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih_fnl%error)
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

   call compute_q_aux_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                          gamma=self%physics%gamma, R=self%physics%R, q_gpu=q_gpu, q_aux_gpu=self%q_aux_gpu)
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

   subroutine destroy(self)
   !< Free device and host data: own device buffers, the FNL helpers (library teardown) and the common data.
   class(flume_fnl_object), intent(inout) :: self !< The equation.

   call free_gpu(self%q_gpu)
   call free_gpu(self%dq_gpu)
   call free_gpu(self%q_aux_gpu)
   call free_gpu(self%flx_f_gpu)
   call free_gpu(self%fly_f_gpu)
   call free_gpu(self%flz_f_gpu)
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
   !< Save fields, restart and conservation history, each on its own cadence; state copied to host only when saved.
   class(flume_fnl_object), intent(inout) :: self !< The equation.

   if (self%time%is_to_save(cadence=self%io%it_save) .or. self%time%is_to_save(cadence=self%io%restart_save)) then
      call self%update_ghost(q_gpu=self%q_gpu)
      call self%copy_gpu_cpu
      if (self%time%is_to_save(cadence=self%io%it_save)) call self%save_xh5f(with_ghost=.true.)
      if (self%time%is_to_save(cadence=self%io%restart_save)) call self%save_restart_files
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
                                       q_inflow=self%bc%q_inflow, q_gpu=q_gpu)
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

   subroutine apply_reflux_to_stage_forest(self, stage, dt, flux_register)
   !< Apply the Berger-Colella reflux correction: not available before milestone 1 phase P5, so a run with AMR seam
   !< faces is refused rather than silently left non-conservative.
   class(flume_fnl_object),     intent(inout) :: self          !< The equation.
   integer(I4P),                intent(in)    :: stage         !< Integrator stage.
   real(R8P),                   intent(in)    :: dt            !< Time step.
   class(flux_register_object), intent(in)    :: flux_register !< Forest's flux register.

   if (.not.flux_register%is_initialized_) return
   if (flux_register%nfaces == 0_I4P) return
   call mpih_fnl%error_stop(msg=': AMR coarse-fine seams need reflux, not implemented yet (stage '//trim(str(stage))// &
                                ', dt '//trim(str(dt))//')')
   endsubroutine apply_reflux_to_stage_forest

   subroutine begin_stage_forest(self, k, K_total, dt, realm)
   !< Begin integrator stage `k` (staged path): publish the stage and compute its state.
   class(flume_fnl_object), intent(inout)                   :: self     !< The equation.
   integer(I4P),            intent(in)                      :: k        !< Stage index (1..K_total).
   integer(I4P),            intent(in)                      :: K_total  !< Forest-wide stage count for this step.
   real(R8P),               intent(in)                      :: dt       !< Time step from the forest.
   class(realm_object),     intent(inout), optional, target :: realm(:) !< Sibling realms (contract parity).

   self%stage_active = k
   call self%rk_fnl%compute_stage(grid=self%adam%grid, field=self%adam%field, s=k, dt=self%time%dt)
   endsubroutine begin_stage_forest

   subroutine close_step_forest(self, dt)
   !< Close a step (staged path): assemble q, save residuals, advance time, clear the active stage.
   class(flume_fnl_object), intent(inout) :: self !< The equation.
   real(R8P),               intent(in)    :: dt   !< Time step from the forest (the local capped value is time%dt).

   call compute_rk_ssp_residual_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, nv=self%nv,                     &
                                    blocks_number=self%blocks_number, nrk=self%rk%nrk, beta_gpu=self%rk_fnl%beta_gpu, &
                                    q_rk_gpu=self%rk_fnl%q_rk_gpu, dq_gpu=self%dq_gpu)
   call self%rk_fnl%update_q(grid=self%adam%grid, field=self%adam%field, rk=self%rk, dt=self%time%dt, q_gpu=self%q_gpu)
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

   call compute_lambda_max_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                               gamma=self%physics%gamma, R=self%physics%R, dxyz_gpu=self%field_fnl%dxyz_gpu,       &
                               is_null=self%adam%grid%null_xyz, q_gpu=self%q_gpu, lambda_max=lambda_max)
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
   call self%rk_fnl%assign_stage(grid=self%adam%grid, field=self%adam%field, s=k, q_gpu=self%dq_gpu)
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
   else
      do i=1, self%ic%amr_iterations
         call self%ic%set_initial_conditions(field=self%adam%field, q=self%q)
         call self%amr_update
      enddo
      call self%ic%set_initial_conditions(field=self%adam%field, q=self%q)
      call self%adam%make_comm_local_maps_ghost_bc
      self%time%time = 0._R8P
      self%time%it   = 0_I4P
   endif
   call self%copy_cpu_gpu(verbose=.true.)
   call self%update_ghost(q_gpu=self%q_gpu)
   call self%compute_q_aux(q_gpu=self%q_gpu)
   call self%diagnostics%open_file(output_basename=self%io%output_basename, q_name=self%q_name, &
                                   is_restart=self%io%restart)
   call self%save_simulation_data
   call self%io%open_file_residuals(nv=self%nv)
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
   !< The face fluxes of a null direction are never computed: they keep their zero initialization (`dev_alloc`).
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

   call self%update_ghost(q_gpu=q_gpu)
   call self%compute_q_aux(q_gpu=q_gpu)
   is_char = self%numerics%reconstruction_variables == RECON_CHARACTERISTIC
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nb=>self%blocks_number, gamma=>self%physics%gamma, &
             is_null=>self%adam%grid%null_xyz, S=>self%weno%S, zeps=>self%weno%zeps, a_gpu=>self%weno_fnl%a_gpu,     &
             p_gpu=>self%weno_fnl%p_gpu, d_gpu=>self%weno_fnl%d_gpu)
   if (.not.is_null(1)) call compute_face_fluxes_dev(d=1_I4P, di=1_I4P, dj=0_I4P, dk=0_I4P, ni=ni, nj=nj, nk=nk, ngc=ngc, &
                                                     blocks_number=nb, S=S, gamma=gamma, is_characteristic=is_char,        &
                                                     weno_a_gpu=a_gpu, weno_p_gpu=p_gpu, weno_d_gpu=d_gpu, weno_zeps=zeps,  &
                                                     q_gpu=q_gpu, q_aux_gpu=self%q_aux_gpu, fl_gpu=self%flx_f_gpu)
   if (.not.is_null(2)) call compute_face_fluxes_dev(d=2_I4P, di=0_I4P, dj=1_I4P, dk=0_I4P, ni=ni, nj=nj, nk=nk, ngc=ngc, &
                                                     blocks_number=nb, S=S, gamma=gamma, is_characteristic=is_char,        &
                                                     weno_a_gpu=a_gpu, weno_p_gpu=p_gpu, weno_d_gpu=d_gpu, weno_zeps=zeps,  &
                                                     q_gpu=q_gpu, q_aux_gpu=self%q_aux_gpu, fl_gpu=self%fly_f_gpu)
   if (.not.is_null(3)) call compute_face_fluxes_dev(d=3_I4P, di=0_I4P, dj=0_I4P, dk=1_I4P, ni=ni, nj=nj, nk=nk, ngc=ngc, &
                                                     blocks_number=nb, S=S, gamma=gamma, is_characteristic=is_char,        &
                                                     weno_a_gpu=a_gpu, weno_p_gpu=p_gpu, weno_d_gpu=d_gpu, weno_zeps=zeps,  &
                                                     q_gpu=q_gpu, q_aux_gpu=self%q_aux_gpu, fl_gpu=self%flz_f_gpu)
   call compute_flux_difference_dev(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, is_null=is_null,                   &
                                    dxyz_gpu=self%field_fnl%dxyz_gpu, flx_f_gpu=self%flx_f_gpu, fly_f_gpu=self%fly_f_gpu, &
                                    flz_f_gpu=self%flz_f_gpu, dq_gpu=dq_gpu)
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
      call self%rk_fnl%compute_stage_ls(grid=self%adam%grid, field=self%adam%field, rk=self%rk, s=s, dt=self%time%dt, &
                                        dq_gpu=self%dq_gpu, q_gpu=self%q_gpu)
   enddo
   endsubroutine integrate_rk_ls_dev

   subroutine integrate_rk_ssp_dev(self)
   !< Integrate one time step with a strong stability preserving Runge-Kutta scheme on the device.
   class(flume_fnl_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: s    !< Counter.

   call self%rk_fnl%initialize_stages(grid=self%adam%grid, field=self%adam%field, q_gpu=self%q_gpu)
   do s=1, self%rk%nrk
      call self%rk_fnl%compute_stage(grid=self%adam%grid, field=self%adam%field, s=s, dt=self%time%dt)
      call self%compute_residuals_dev(q_gpu=self%rk_fnl%q_rk_gpu(:,:,:,:,:,s), dq_gpu=self%dq_gpu, s=s)
      call self%rk_fnl%assign_stage(grid=self%adam%grid, field=self%adam%field, s=s, q_gpu=self%dq_gpu)
   enddo
   call compute_rk_ssp_residual_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, nv=self%nv,                     &
                                    blocks_number=self%blocks_number, nrk=self%rk%nrk, beta_gpu=self%rk_fnl%beta_gpu, &
                                    q_rk_gpu=self%rk_fnl%q_rk_gpu, dq_gpu=self%dq_gpu)
   call self%rk_fnl%update_q(grid=self%adam%grid, field=self%adam%field, rk=self%rk, dt=self%time%dt, q_gpu=self%q_gpu)
   call self%save_residuals
   endsubroutine integrate_rk_ssp_dev
endmodule adam_flume_fnl_object
