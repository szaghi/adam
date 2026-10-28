!< RIPPLE, the FNL (OpenACC) backend realm.

#include "fundal.H"

module adam_ripple_fnl_object
!< RIPPLE, the FNL (OpenACC) backend realm.
!<
!< The CPU backend's twin with the state on the device: a sibling of `ripple_cpu_object`, extending the same common
!< realm. The initial condition and the initial AMR run on the host, then the state goes to the device and stays
!< there; it comes back only for the outputs. The library's FNL helpers do the device side of the library objects:
!< `field_fnl` the ghost exchange (local copies, 2:1 transfers, GPU-direct MPI), `rk_fnl` the Runge-Kutta stages,
!< `weno_fnl` the WENO coefficients. Device arrays are transposed, `(b, i, j, k, v)`.

! ADAM classes, libraries, parameters
use :: adam_flux_register_object, only : flux_register_object
use :: adam_maps_object,          only : face_axis_sign
use :: adam_realm_object,         only : realm_pointer
use :: adam_fnl_field_kernels,    only : pack_blocks_transposed_dev, unpack_blocks_transposed_dev
use :: adam_fnl_field_object,     only : field_fnl_object
use :: adam_fnl_rk_object,        only : rk_fnl_object
use :: adam_fnl_weno_object,      only : weno_fnl_object
! ADAM singleton objects
use :: adam_fnl_mpih_global,      only : count_device_sharers, mpih_fnl, mpih_fnl_device_sharers, mpih_fnl_is_initialized
! RIPPLE modules
use :: adam_ripple_common_object, only : ripple_common_object
use :: adam_ripple_fnl_kernels,   only : apply_reflux_face_dev, compute_face_fluxes_dev, compute_flux_difference_dev, &
                                         pack_face_skin_dev
! third party modules
use :: fundal,                    only : dev_alloc, dev_free, dev_memcpy_from_device, dev_memcpy_to_device, mydev
use :: mpi
use :: penf,                      only : I4P, R8P

implicit none
private
public :: ripple_fnl_object

type, extends(ripple_common_object) :: ripple_fnl_object
   !< RIPPLE FNL realm.
   type(field_fnl_object) :: field_fnl                !< Field helper: coordinates, maps, ghost exchange.
   type(rk_fnl_object)    :: rk_fnl                   !< Runge-Kutta helper: the stages on the device.
   type(weno_fnl_object)  :: weno_fnl                 !< WENO helper: the coefficients on the device.
   real(R8P), pointer     :: q_gpu(:,:,:,:,:)=>null() !< State [nb, 1-ngc:ni+ngc, ..., 1].
   real(R8P), pointer     :: dq_gpu(:,:,:,:,:)=>null()  !< Residuals, same shape.
   real(R8P), pointer     :: flx_gpu(:,:,:,:,:)=>null() !< X-face fluxes [nb, 0:ni, 1:nj, 1:nk, 1].
   real(R8P), pointer     :: fly_gpu(:,:,:,:,:)=>null() !< Y-face fluxes [nb, 1:ni, 0:nj, 1:nk, 1].
   real(R8P), pointer     :: flz_gpu(:,:,:,:,:)=>null() !< Z-face fluxes [nb, 1:ni, 1:nj, 0:nk, 1].
   real(R8P), pointer     :: face_gpu(:,:)=>null()      !< One face's skin or reflux mismatch [1, largest face].
   real(R8P), allocatable :: face(:,:)                  !< Its host copy, the same extent.
   contains
      ! private methods
      procedure, pass(self) :: allocate_gpu            !< Allocate the device data.
      procedure, pass(self) :: copy_state              !< Copy q between host and device.
      procedure, pass(self) :: update_ghost            !< Fill the ghost cells on the device.
      procedure, pass(self) :: compute_residuals       !< Residuals of one device state.
      procedure, pass(self) :: accumulate_seam_fluxes  !< Seam face fluxes of one stage to the flux register.
      procedure, pass(self) :: compute_total           !< Total of u dV on this rank, from the device state.
      ! forest methods: core layer
      procedure, pass(self) :: initialize_forest       !< Initialize the realm.
      procedure, pass(self) :: advance_one_step_forest !< Advance one whole step.
      procedure, pass(self) :: post_step_forest        !< Save the outputs on their cadence.
      procedure, pass(self) :: finalize_forest         !< Finalize the realm.
      procedure, pass(self) :: finalize_mpi_forest     !< Finalize the FNL MPI handler.
      ! forest methods: staged layer
      procedure, pass(self) :: open_step_forest             !< Open a step.
      procedure, pass(self) :: begin_stage_forest           !< Begin a stage.
      procedure, pass(self) :: end_stage_forest             !< End a stage.
      procedure, pass(self) :: close_step_forest            !< Close a step.
      procedure, pass(self) :: apply_reflux_to_stage_forest !< Reflux the committed step on the device.
endtype ripple_fnl_object

contains
   ! private methods
   subroutine allocate_gpu(self)
   !< Allocate the device data, every `dev_alloc` checked.
   class(ripple_fnl_object), intent(inout) :: self !< The realm.
   integer(I4P)                            :: ierr !< Error status.

   associate(nb=>self%nb, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk)
   call dev_alloc(fptr_dev=self%q_gpu, lbounds=[1,1-ngc,1-ngc,1-ngc,1], ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,1], &
                  init_value=0._R8P, ierr=ierr) ; call check('q_gpu')
   call dev_alloc(fptr_dev=self%dq_gpu, lbounds=[1,1-ngc,1-ngc,1-ngc,1], ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,1], &
                  init_value=0._R8P, ierr=ierr) ; call check('dq_gpu')
   call dev_alloc(fptr_dev=self%flx_gpu, lbounds=[1,0,1,1,1], ubounds=[nb,ni,nj,nk,1], init_value=0._R8P, ierr=ierr)
   call check('flx_gpu')
   call dev_alloc(fptr_dev=self%fly_gpu, lbounds=[1,1,0,1,1], ubounds=[nb,ni,nj,nk,1], init_value=0._R8P, ierr=ierr)
   call check('fly_gpu')
   call dev_alloc(fptr_dev=self%flz_gpu, lbounds=[1,1,1,0,1], ubounds=[nb,ni,nj,nk,1], init_value=0._R8P, ierr=ierr)
   call check('flz_gpu')
   call dev_alloc(fptr_dev=self%face_gpu, lbounds=[1,1], ubounds=[1,max(nj*nk, ni*nk, ni*nj)], init_value=0._R8P, &
                  ierr=ierr) ; call check('face_gpu')
   ! the host copy has the device buffer's extent: a copy between arrays of different extents goes through a host
   ! temporary and fails on the device (issue #31)
   allocate(self%face(1,max(nj*nk, ni*nk, ni*nj)), source=0._R8P)
   endassociate
   contains
      subroutine check(name)
      !< Stop on a failed device allocation.
      character(*), intent(in) :: name !< Array name.

      if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate '//name//' in ripple_fnl_object')
      endsubroutine check
   endsubroutine allocate_gpu

   subroutine copy_state(self, to_device)
   !< Copy the blocks in use of the state between the host `q(v, i, j, k, b)` and the device `q_gpu(b, i, j, k, v)`: a
   !< device kernel transposes them into a scratch buffer of the host slab's shape, moved by one contiguous copy.
   class(ripple_fnl_object), intent(inout) :: self      !< The realm.
   logical,                  intent(in)    :: to_device !< Host to device, else back.
   real(R8P), pointer                      :: buf_gpu(:,:,:,:,:) !< Device packed blocks.
   integer(I4P)                            :: ierr      !< Error status.

   if (self%blocks_number < 1_I4P) return
   associate(nb=>self%blocks_number, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk)
   call dev_alloc(fptr_dev=buf_gpu, lbounds=[1,1-ngc,1-ngc,1-ngc,1], ubounds=[1,ni+ngc,nj+ngc,nk+ngc,nb], ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate buf_gpu in ripple_fnl_object%copy_state')
   if (to_device) then
      call dev_memcpy_to_device(dst=buf_gpu, src=self%q(:,:,:,:,1:nb))
      call unpack_blocks_transposed_dev(blocks_number=nb, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=1, buf_gpu=buf_gpu, &
                                        a_gpu=self%q_gpu)
   else
      call pack_blocks_transposed_dev(blocks_number=nb, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=1, a_gpu=self%q_gpu, &
                                      buf_gpu=buf_gpu)
      call dev_memcpy_from_device(dst=self%q(:,:,:,:,1:nb), src=buf_gpu)
   endif
   endassociate
   call dev_free(buf_gpu, mydev)
   endsubroutine copy_state

   subroutine update_ghost(self, q_gpu)
   !< Fill the ghost cells on the device: local copies and 2:1 transfers, then the GPU-direct MPI exchange.
   class(ripple_fnl_object), intent(inout) :: self                                            !< The realm.
   real(R8P),                intent(inout) :: q_gpu(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< State.

   call self%field_fnl%update_ghost_local_gpu(q_gpu=q_gpu)
   call self%field_fnl%update_ghost_mpi_gpu(comm_map_send_ptr_ghost=self%adam%maps%comm_map_send_ptr_ghost, &
                                            comm_map_recv_ptr_ghost=self%adam%maps%comm_map_recv_ptr_ghost, q_gpu=q_gpu)
   endsubroutine update_ghost

   subroutine compute_residuals(self, q_gpu, s, flux_register)
   !< Compute the residuals `dq_gpu` of the device state `q_gpu` (the CPU `compute_residuals` on the device).
   class(ripple_fnl_object),    intent(inout)           :: self          !< The realm.
   real(R8P),                   intent(inout)           :: q_gpu(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< State.
   integer(I4P),                intent(in),    optional :: s             !< Runge-Kutta stage.
   class(flux_register_object), intent(inout), optional :: flux_register !< Forest's flux register.

   call self%update_ghost(q_gpu=q_gpu)
   associate(null=>self%adam%grid%null_xyz, ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, &
             nb=>self%blocks_number, ws=>self%weno%S, w=>self%weno_fnl)
   if (.not.null(1)) call compute_face_fluxes_dev(di=1, dj=0, dk=0, ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, S=ws, &
                                                  a=self%a(1), weno_a_gpu=w%a_gpu, weno_p_gpu=w%p_gpu,                  &
                                                  weno_d_gpu=w%d_gpu, weno_zeps=self%weno%zeps, q_gpu=q_gpu,            &
                                                  fl_gpu=self%flx_gpu)
   if (.not.null(2)) call compute_face_fluxes_dev(di=0, dj=1, dk=0, ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, S=ws, &
                                                  a=self%a(2), weno_a_gpu=w%a_gpu, weno_p_gpu=w%p_gpu,                  &
                                                  weno_d_gpu=w%d_gpu, weno_zeps=self%weno%zeps, q_gpu=q_gpu,            &
                                                  fl_gpu=self%fly_gpu)
   if (.not.null(3)) call compute_face_fluxes_dev(di=0, dj=0, dk=1, ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, S=ws, &
                                                  a=self%a(3), weno_a_gpu=w%a_gpu, weno_p_gpu=w%p_gpu,                  &
                                                  weno_d_gpu=w%d_gpu, weno_zeps=self%weno%zeps, q_gpu=q_gpu,            &
                                                  fl_gpu=self%flz_gpu)
   if (present(flux_register) .and. present(s)) then
      if (flux_register%nfaces > 0_I4P) call self%accumulate_seam_fluxes(s=s, flux_register=flux_register)
   endif
   call compute_flux_difference_dev(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=nb, w=merge(0._R8P, 1._R8P, null), &
                                    dxyz_gpu=self%field_fnl%dxyz_gpu, flx_gpu=self%flx_gpu, fly_gpu=self%fly_gpu,  &
                                    flz_gpu=self%flz_gpu, dq_gpu=self%dq_gpu)
   endassociate
   endsubroutine compute_residuals

   subroutine accumulate_seam_fluxes(self, s, flux_register)
   !< Pack the face flux skin of every registered seam face on the device, bring it to the host and route it to the
   !< (host) register, weighted by the SSP coefficient of stage `s`. One small copy per face: simple, not fast.
   class(ripple_fnl_object),    intent(inout) :: self          !< The realm.
   integer(I4P),                intent(in)    :: s             !< Runge-Kutta stage.
   class(flux_register_object), intent(inout) :: flux_register !< Forest's flux register.
   integer(I4P)                               :: b, fec, n     !< Block, face, skin cells.

   do b=1, self%blocks_number
      do fec=1, 6
         if (self%adam%maps%inter_realm_face_register_index(b, fec) == 0_I4P) cycle
         call pack_face_skin_dev(fec=fec, b=b, ni=self%ni, nj=self%nj, nk=self%nk, flx_gpu=self%flx_gpu, &
                                 fly_gpu=self%fly_gpu, flz_gpu=self%flz_gpu, skin_gpu=self%face_gpu)
         call dev_memcpy_from_device(dst=self%face, src=self%face_gpu)
         select case(fec)
         case(1, 2) ; n = self%nj * self%nk
         case(3, 4) ; n = self%ni * self%nk
         case default ; n = self%ni * self%nj
         endselect
         call self%accumulate_seam_skin(flux_register=flux_register, b=b, fec=fec, weight=self%rk%beta(s), &
                                        skin=self%face(:,1:n))
      enddo
   enddo
   endsubroutine accumulate_seam_fluxes

   function compute_total(self) result(total)
   !< Return the total of u dV over the interior cells of this rank's blocks, reduced on the device.
   class(ripple_fnl_object), intent(in) :: self                 !< The realm.
   real(R8P)                            :: total                !< Total of u dV.
   real(R8P), pointer                   :: q_gpu(:,:,:,:,:)     !< State.
   real(R8P), pointer                   :: dxyz_gpu(:,:)        !< Block spacings.
   real(R8P)                            :: t                    !< Reduction accumulator.
   real(R8P)                            :: wx, wy, wz           !< Axis weights of the cell volume.
   integer(I4P)                         :: ni, nj, nk, nb       !< Grid dimensions, blocks in use.
   integer(I4P)                         :: b, i, j, k           !< Counters.

   q_gpu => self%q_gpu ; dxyz_gpu => self%field_fnl%dxyz_gpu
   ni = self%ni ; nj = self%nj ; nk = self%nk ; nb = self%blocks_number
   wx = merge(0._R8P, 1._R8P, self%adam%grid%null_xyz(1))
   wy = merge(0._R8P, 1._R8P, self%adam%grid%null_xyz(2))
   wz = merge(0._R8P, 1._R8P, self%adam%grid%null_xyz(3))
   t = 0._R8P
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,dxyz_gpu) &
   !$acc& firstprivate(ni,nj,nk,nb,wx,wy,wz) reduction(+:t)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, nb
      t = t + q_gpu(b,i,j,k,1) * (wx * dxyz_gpu(b,1) + 1._R8P - wx) * (wy * dxyz_gpu(b,2) + 1._R8P - wy) * &
                                 (wz * dxyz_gpu(b,3) + 1._R8P - wz)
   enddo
   enddo
   enddo
   enddo
   total = t
   endfunction compute_total

   ! forest methods: core layer
   subroutine initialize_forest(self, filename, realms_number, memory_avail, nv, verbose)
   !< Initialize the realm: MPI and the device (once per process), the common data within the device budget, the FNL
   !< helpers, the device data; initial AMR and initial condition on the host, then the state to the device.
   class(ripple_fnl_object), intent(inout)        :: self          !< The realm.
   character(*),             intent(in)           :: filename      !< Input file name.
   integer(I4P),             intent(in), optional :: realms_number !< Unused: RIPPLE is single-realm.
   real(R8P),                intent(in), optional :: memory_avail  !< Unused: the budget comes from the device.
   integer(I4P),             intent(in), optional :: nv            !< Unused: nv = 1.
   logical,                  intent(in), optional :: verbose       !< Unused.
   logical                                        :: is_mpi        !< MPI already initialized.

   if (.not.mpih_fnl_is_initialized) then
      call MPI_INITIALIZED(is_mpi, mpih_fnl%error)
      call mpih_fnl%initialize(do_mpi_init=.not.is_mpi, do_device_init=.true., verbose=.true.)
      mpih_fnl_is_initialized = .true.
      mpih_fnl_device_sharers = count_device_sharers()
   endif
   ! the blocks budget is the device memory, shared by the ranks bound to the same device
   call self%initialize_ripple(filename=filename, memory_avail=real(mpih_fnl%dev_memory_total, R8P) / 1e9_R8P / &
                                                               real(mpih_fnl_device_sharers, R8P))
   call self%field_fnl%initialize(grid=self%adam%grid, field=self%adam%field, maps=self%adam%maps, verbose=.true.)
   call self%rk_fnl%initialize(grid=self%adam%grid, field=self%adam%field, rk=self%rk)
   call self%weno_fnl%initialize(weno=self%weno)
   call self%allocate_gpu
   if (self%amr%markers_number > 0) then
      call self%set_initial_conditions ! a registered criterion may mark on the initial state (issue #110)
      call self%amr_update
   endif
   call self%set_initial_conditions
   call self%field_fnl%copy_cpu_gpu(field=self%adam%field, maps=self%adam%maps, verbose=.true.)
   call self%copy_state(to_device=.true.)
   call self%save_output(fields=.true.)
   endsubroutine initialize_forest

   subroutine advance_one_step_forest(self, dt)
   !< Advance one whole SSP Runge-Kutta step on the device (fast path: no 2:1 face).
   class(ripple_fnl_object), intent(inout) :: self !< The realm.
   real(R8P),                intent(in)    :: dt   !< Time step from the forest.
   integer(I4P)                            :: s    !< Stage counter.

   call self%advance_time(dt=dt)
   call self%rk_fnl%initialize_stages(grid=self%adam%grid, field=self%adam%field, q_gpu=self%q_gpu)
   do s=1, self%rk%nrk
      call self%rk_fnl%compute_stage(grid=self%adam%grid, field=self%adam%field, s=s, dt=self%dt)
      call self%compute_residuals(q_gpu=self%rk_fnl%q_rk_gpu(:,:,:,:,:,s))
      call self%rk_fnl%assign_stage(grid=self%adam%grid, field=self%adam%field, s=s, q_gpu=self%dq_gpu)
   enddo
   call self%rk_fnl%update_q(grid=self%adam%grid, field=self%adam%field, rk=self%rk, dt=self%dt, q_gpu=self%q_gpu)
   self%time = self%time + self%dt
   endsubroutine advance_one_step_forest

   subroutine post_step_forest(self, dt, t, it, do_save_state, do_save_residuals, do_save_restart, do_amr, realm)
   !< Append the history row every step (the total reduced on the device); bring q to the host and save the fields
   !< every `[IO] it_save` steps and at the end.
   class(ripple_fnl_object), intent(inout)                   :: self              !< The realm.
   real(R8P),                intent(in)                      :: dt                !< Unused.
   real(R8P),                intent(in)                      :: t                 !< Unused.
   integer(I4P),             intent(in)                      :: it                !< Unused.
   logical,                  intent(in),    optional         :: do_save_state     !< Unused.
   logical,                  intent(in),    optional         :: do_save_residuals !< Unused.
   logical,                  intent(in),    optional         :: do_save_restart   !< Unused.
   logical,                  intent(in),    optional         :: do_amr            !< Unused.
   type(realm_pointer),      intent(inout), optional, target :: realm(:)          !< Unused.
   logical                                                   :: done, fields      !< The run is over, save the fields.

   call self%is_done_forest(done=done)
   fields = done .or. (self%io%it_save > 0 .and. mod(self%it, max(self%io%it_save, 1)) == 0)
   if (fields) call self%copy_state(to_device=.false.)
   call self%save_output(fields=fields)
   endsubroutine post_step_forest

   subroutine finalize_forest(self)
   !< Finalize the realm: close the history, free the device data (the forest finalizes MPI once).
   class(ripple_fnl_object), intent(inout) :: self !< The realm.

   if (self%history_unit /= 0) close(self%history_unit)
   call dev_free(self%q_gpu, mydev)
   call dev_free(self%dq_gpu, mydev)
   call dev_free(self%flx_gpu, mydev)
   call dev_free(self%fly_gpu, mydev)
   call dev_free(self%flz_gpu, mydev)
   call dev_free(self%face_gpu, mydev)
   call self%weno_fnl%destroy
   call self%rk_fnl%destroy
   call self%field_fnl%destroy
   endsubroutine finalize_forest

   subroutine finalize_mpi_forest(self)
   !< Finalize the FNL MPI handler (the forest calls it once, after every realm is finalized).
   class(ripple_fnl_object), intent(inout) :: self !< The realm (carries no MPI state).

   call mpih_fnl%finalize
   endsubroutine finalize_mpi_forest

   ! forest methods: staged layer
   subroutine open_step_forest(self, dt)
   !< Open a step: count it, cap dt, copy q into every stage.
   class(ripple_fnl_object), intent(inout) :: self !< The realm.
   real(R8P),                intent(in)    :: dt   !< Time step from the forest.

   call self%advance_time(dt=dt)
   call self%rk_fnl%initialize_stages(grid=self%adam%grid, field=self%adam%field, q_gpu=self%q_gpu)
   endsubroutine open_step_forest

   subroutine begin_stage_forest(self, k, K_total, dt, realm)
   !< Begin stage `k`: publish it (`stage_active`) and build its state.
   class(ripple_fnl_object), intent(inout)                   :: self     !< The realm.
   integer(I4P),             intent(in)                      :: k        !< Stage index.
   integer(I4P),             intent(in)                      :: K_total  !< Unused.
   real(R8P),                intent(in)                      :: dt       !< Unused: the capped dt is self%dt.
   type(realm_pointer),      intent(inout), optional, target :: realm(:) !< Unused.

   self%stage_active = k
   call self%rk_fnl%compute_stage(grid=self%adam%grid, field=self%adam%field, s=k, dt=self%dt)
   endsubroutine begin_stage_forest

   subroutine end_stage_forest(self, k, K_total, dt, realm, flux_register)
   !< End stage `k`: its residuals, the seam fluxes to the register, the stage assignment.
   class(ripple_fnl_object),    intent(inout)                   :: self          !< The realm.
   integer(I4P),                intent(in)                      :: k             !< Stage index.
   integer(I4P),                intent(in)                      :: K_total       !< Unused.
   real(R8P),                   intent(in)                      :: dt            !< Unused.
   type(realm_pointer),         intent(inout), optional, target :: realm(:)      !< Unused.
   class(flux_register_object), intent(inout), optional         :: flux_register !< Forest's flux register.

   call self%compute_residuals(q_gpu=self%rk_fnl%q_rk_gpu(:,:,:,:,:,k), s=k, flux_register=flux_register)
   call self%rk_fnl%assign_stage(grid=self%adam%grid, field=self%adam%field, s=k, q_gpu=self%dq_gpu)
   endsubroutine end_stage_forest

   subroutine close_step_forest(self, dt)
   !< Close a step: commit q, advance the time, clear the active stage.
   class(ripple_fnl_object), intent(inout) :: self !< The realm.
   real(R8P),                intent(in)    :: dt   !< Unused: the capped dt is self%dt.

   call self%rk_fnl%update_q(grid=self%adam%grid, field=self%adam%field, rk=self%rk, dt=self%dt, q_gpu=self%q_gpu)
   self%time = self%time + self%dt
   self%stage_active = 0_I4P
   endsubroutine close_step_forest

   subroutine apply_reflux_to_stage_forest(self, stage, dt, flux_register)
   !< Reflux the committed device state (the forest calls it after `close_step_forest`, with the last stage): the host
   !< mismatch `F_coarse - F_fine` of every coarse face this rank owns goes to the device and is added, scaled by
   !< `sgn dt / dx`, to its coarse skin cells.
   class(ripple_fnl_object),    intent(inout) :: self          !< The realm.
   integer(I4P),                intent(in)    :: stage         !< Stage: the correction acts at the last one.
   real(R8P),                   intent(in)    :: dt            !< Unused: the capped dt is self%dt.
   class(flux_register_object), intent(in)    :: flux_register !< Forest's flux register.
   integer(I4P)                               :: f, axis, sgn  !< Face counter, face normal and side.

   if (stage /= self%rk%nrk .or. flux_register%nfaces == 0_I4P .or. .not.allocated(flux_register%face)) return
   do f=1, flux_register%nfaces
      associate(face=>flux_register%face(f))
      if (face%coarse_realm /= self%realm_index .or. face%coarse_rank /= mpih_fnl%myrank) cycle
      call face_axis_sign(face_code=face%coarse_face, axis=axis, sgn=sgn)
      self%face = 0._R8P
      self%face(1,1:face%nface_cells) = face%F_coarse(1,:,1) - face%F_fine_sum(1,:,1)
      call dev_memcpy_to_device(dst=self%face_gpu, src=self%face)
      call apply_reflux_face_dev(axis=axis, sgn=sgn, b=face%coarse_block, ni=self%ni, nj=self%nj, nk=self%nk,   &
                                 ngc=self%ngc, nface_cells=face%nface_cells,                                    &
                                 scale=sgn * self%dt / self%adam%field%dxyz(axis,face%coarse_block),          &
                                 delta_gpu=self%face_gpu, q_gpu=self%q_gpu)
      endassociate
   enddo
   endsubroutine apply_reflux_to_stage_forest
endmodule adam_ripple_fnl_object
