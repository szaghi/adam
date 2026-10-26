!< RIPPLE, the CPU backend realm.
module adam_ripple_cpu_object
!< RIPPLE, the CPU backend realm.
!<
!< The residuals (ghost update, WENO face fluxes, flux difference) and the forest contract on the host. Two paths reach
!< the integration: the fast path (`advance_one_step_forest`, a whole SSP Runge-Kutta step at once) when the grid has no
!< 2:1 face, the staged path (`open_step`, `begin_stage`, `end_stage`, `close_step`, then the reflux) when it has one:
!< there the forest's flux register collects the seam face fluxes of every stage, and the end-of-step reflux makes the
!< coarse and fine sides of each 2:1 face exchange the same flux.

! ADAM classes, libraries, parameters
use :: adam_flux_register_object, only : flux_register_object
use :: adam_realm_object,         only : realm_pointer
use :: adam_weno_object,          only : weno_reconstruct_upwind, S_MAX
! ADAM singleton objects
use :: adam_mpih_global,          only : mpih
! RIPPLE modules
use :: adam_ripple_common_object, only : ripple_common_object
! third party modules
use :: penf

implicit none
private
public :: ripple_cpu_object

type, extends(ripple_common_object) :: ripple_cpu_object
   !< RIPPLE CPU realm.
   real(R8P), allocatable :: flx(:,:,:,:,:) !< X-face fluxes [1, 0:ni, 1:nj, 1:nk, nb].
   real(R8P), allocatable :: fly(:,:,:,:,:) !< Y-face fluxes [1, 1:ni, 0:nj, 1:nk, nb].
   real(R8P), allocatable :: flz(:,:,:,:,:) !< Z-face fluxes [1, 1:ni, 1:nj, 0:nk, nb].
   contains
      ! private methods
      procedure, pass(self) :: update_ghost            !< Fill the ghost cells.
      procedure, pass(self) :: compute_residuals       !< Residuals of one state.
      procedure, pass(self) :: accumulate_seam_fluxes  !< Seam face fluxes of one stage to the flux register.
      ! forest methods: core layer (single realm, no 2:1 face)
      procedure, pass(self) :: initialize_forest       !< Initialize the realm.
      procedure, pass(self) :: advance_one_step_forest !< Advance one whole step.
      procedure, pass(self) :: post_step_forest        !< Save the outputs on their cadence.
      procedure, pass(self) :: finalize_forest         !< Finalize the realm.
      ! forest methods: staged layer (2:1 AMR faces)
      procedure, pass(self) :: open_step_forest             !< Open a step.
      procedure, pass(self) :: begin_stage_forest           !< Begin a stage.
      procedure, pass(self) :: end_stage_forest             !< End a stage.
      procedure, pass(self) :: close_step_forest            !< Close a step.
      procedure, pass(self) :: apply_reflux_to_stage_forest !< Reflux the committed step.
endtype ripple_cpu_object

contains
   ! private methods
   subroutine update_ghost(self, q)
   !< Fill the ghost cells: copies (and 2:1 restriction/prolongation) between local blocks, then the MPI exchange. Every
   !< face is periodic, so there is no boundary condition to apply.
   class(ripple_cpu_object), intent(inout) :: self                                       !< The realm.
   real(R8P),                intent(inout) :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< State.

   call self%adam%field%update_ghost_local(grid=self%adam%grid, maps=self%adam%maps, q=q)
   call self%adam%field%update_ghost_mpi(grid=self%adam%grid, maps=self%adam%maps, q=q)
   endsubroutine update_ghost

   subroutine compute_residuals(self, q, s, flux_register)
   !< Compute the residuals `dq` of state `q`: ghost update, WENO face fluxes of each non-null axis, flux difference. On
   !< the staged path (`flux_register` present), the seam face fluxes of stage `s` go to the register.
   class(ripple_cpu_object),    intent(inout)           :: self          !< The realm.
   real(R8P),                   intent(inout)           :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< State.
   integer(I4P),                intent(in),    optional :: s             !< Runge-Kutta stage.
   class(flux_register_object), intent(inout), optional :: flux_register !< Forest's flux register.
   integer(I4P)                                         :: b, i, j, k    !< Counters.

   call self%update_ghost(q=q)
   associate(null=>self%adam%grid%null_xyz, dxyz=>self%adam%field%dxyz)
   if (.not.null(1)) call face_fluxes(di=1, dj=0, dk=0, a=self%a(1), fl=self%flx)
   if (.not.null(2)) call face_fluxes(di=0, dj=1, dk=0, a=self%a(2), fl=self%fly)
   if (.not.null(3)) call face_fluxes(di=0, dj=0, dk=1, a=self%a(3), fl=self%flz)
   if (present(flux_register) .and. present(s)) then
      if (flux_register%nfaces > 0_I4P) call self%accumulate_seam_fluxes(s=s, flux_register=flux_register)
   endif
   do b=1, self%blocks_number
      do k=1, self%nk
         do j=1, self%nj
            do i=1, self%ni
               self%dq(1,i,j,k,b) = -(merge(0._R8P, (self%flx(1,i,j,k,b) - self%flx(1,i-1,j,k,b)) / dxyz(1,b), null(1)) + &
                                      merge(0._R8P, (self%fly(1,i,j,k,b) - self%fly(1,i,j-1,k,b)) / dxyz(2,b), null(2)) + &
                                      merge(0._R8P, (self%flz(1,i,j,k,b) - self%flz(1,i,j,k-1,b)) / dxyz(3,b), null(3)))
            enddo
         enddo
      enddo
   enddo
   endassociate
   contains
      subroutine face_fluxes(di, dj, dk, a, fl)
      !< Face fluxes along the axis of unit step (di, dj, dk): the split fluxes `a+ u` (cells `m = 1-S ... S-1`) and
      !< `a- u` (cells `m+1`) of the stencil of face `i+1/2`, in the layout of the WENO upwind primitive, reconstructed
      !< and summed.
      integer(I4P), intent(in)    :: di, dj, dk                    !< Unit step.
      real(R8P),    intent(in)    :: a                             !< Velocity component.
      real(R8P),    intent(inout) :: fl(1:,1-di:,1-dj:,1-dk:,1:)   !< Face fluxes.
      real(R8P)                   :: v(2,2*S_MAX-1)                !< Packed stencil.
      real(R8P)                   :: vr(2)                         !< Reconstructed split fluxes.
      integer(I4P)                :: S, m                          !< WENO half stencil, counter.

      S = self%weno%S
      do b=1, self%blocks_number
         do k=1-dk, self%nk
            do j=1-dj, self%nj
               do i=1-di, self%ni
                  do m=1-S, S-1
                     v(2,m+S) = max(a, 0._R8P) * q(1,i+m*di,j+m*dj,k+m*dk,b)
                     v(1,m+S) = min(a, 0._R8P) * q(1,i+(m+1)*di,j+(m+1)*dj,k+(m+1)*dk,b)
                  enddo
                  call weno_reconstruct_upwind(S=S, weno_a=self%weno%a, weno_p=self%weno%p, weno_d=self%weno%d, &
                                               weno_zeps=self%weno%zeps, weno_rmu=1._R8P, v=v, vr=vr)
                  fl(1,i,j,k,b) = vr(1) + vr(2)
               enddo
            enddo
         enddo
      enddo
      endsubroutine face_fluxes
   endsubroutine compute_residuals

   subroutine accumulate_seam_fluxes(self, s, flux_register)
   !< Pack the face flux skin of every registered seam face of every block and route it to the register, weighted by
   !< the SSP coefficient of stage `s`.
   class(ripple_cpu_object),    intent(inout) :: self          !< The realm.
   integer(I4P),                intent(in)    :: s             !< Runge-Kutta stage.
   class(flux_register_object), intent(inout) :: flux_register !< Forest's flux register.
   integer(I4P)                               :: b, fec, f     !< Counters.

   do b=1, self%blocks_number
      do fec=1, 6
         if (self%adam%maps%inter_realm_face_register_index(b, fec) == 0_I4P) cycle
         f = merge(0_I4P, 1_I4P, mod(fec, 2_I4P) == 1_I4P) ! the min face is face 0, the max face is face n
         select case(fec)
         case(1, 2)
            call self%accumulate_seam_skin(flux_register=flux_register, b=b, fec=fec, weight=self%rk%beta(s), &
                                           skin=reshape(self%flx(1:1,f*self%ni,1:self%nj,1:self%nk,b), [1, self%nj*self%nk]))
         case(3, 4)
            call self%accumulate_seam_skin(flux_register=flux_register, b=b, fec=fec, weight=self%rk%beta(s), &
                                           skin=reshape(self%fly(1:1,1:self%ni,f*self%nj,1:self%nk,b), [1, self%ni*self%nk]))
         case default
            call self%accumulate_seam_skin(flux_register=flux_register, b=b, fec=fec, weight=self%rk%beta(s), &
                                           skin=reshape(self%flz(1:1,1:self%ni,1:self%nj,f*self%nk,b), [1, self%ni*self%nj]))
         endselect
      enddo
   enddo
   endsubroutine accumulate_seam_fluxes

   ! forest methods: core layer
   subroutine initialize_forest(self, filename, realms_number, memory_avail, nv, verbose)
   !< Initialize the realm: common data, face fluxes, initial AMR, initial condition, first output.
   class(ripple_cpu_object), intent(inout)        :: self          !< The realm.
   character(*),             intent(in)           :: filename      !< Input file name.
   integer(I4P),             intent(in), optional :: realms_number !< Unused: RIPPLE is single-realm.
   real(R8P),                intent(in), optional :: memory_avail  !< Unused: the budget comes from the MPI handler.
   integer(I4P),             intent(in), optional :: nv            !< Unused: nv = 1.
   logical,                  intent(in), optional :: verbose       !< Unused.

   call mpih%initialize(do_mpi_init=.true., verbose=.true.)
   call self%initialize_ripple(filename=filename, memory_avail=mpih%memory_avail)
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, nb=>self%nb)
      allocate(self%flx(1,0:ni,1:nj,1:nk,nb), self%fly(1,1:ni,0:nj,1:nk,nb), self%flz(1,1:ni,1:nj,0:nk,nb), source=0._R8P)
   endassociate
   if (self%amr%markers_number > 0) call self%amr_update
   call self%set_initial_conditions
   call self%save_output(fields=.true.)
   endsubroutine initialize_forest

   subroutine advance_one_step_forest(self, dt)
   !< Advance one whole SSP Runge-Kutta step (fast path: no 2:1 face).
   class(ripple_cpu_object), intent(inout) :: self !< The realm.
   real(R8P),                intent(in)    :: dt   !< Time step from the forest.
   integer(I4P)                            :: s    !< Stage counter.

   call self%advance_time(dt=dt)
   call self%rk%initialize_stages(field=self%adam%field, q=self%q)
   do s=1, self%rk%nrk
      call self%rk%compute_stage(field=self%adam%field, s=s, dt=self%dt)
      call self%compute_residuals(q=self%rk%q_rk(:,:,:,:,:,s))
      call self%rk%assign_stage(field=self%adam%field, s=s, q=self%dq)
   enddo
   call self%rk%update_q(field=self%adam%field, dt=self%dt, q=self%q, dq=self%dq)
   self%time = self%time + self%dt
   endsubroutine advance_one_step_forest

   subroutine post_step_forest(self, dt, t, it, do_save_state, do_save_residuals, do_save_restart, do_amr, realm)
   !< Append the history row every step; save the fields every `[IO] it_save` steps and at the end.
   class(ripple_cpu_object), intent(inout)                   :: self              !< The realm.
   real(R8P),                intent(in)                      :: dt                !< Unused: the realm keeps its own.
   real(R8P),                intent(in)                      :: t                 !< Unused.
   integer(I4P),             intent(in)                      :: it                !< Unused.
   logical,                  intent(in),    optional         :: do_save_state     !< Unused.
   logical,                  intent(in),    optional         :: do_save_residuals !< Unused.
   logical,                  intent(in),    optional         :: do_save_restart   !< Unused.
   logical,                  intent(in),    optional         :: do_amr            !< Unused.
   type(realm_pointer),      intent(inout), optional, target :: realm(:)          !< Unused.
   logical                                                   :: done              !< The run is over.

   call self%is_done_forest(done=done)
   call self%save_output(fields=done .or. (self%io%it_save > 0 .and. mod(self%it, max(self%io%it_save, 1)) == 0))
   endsubroutine post_step_forest

   subroutine finalize_forest(self)
   !< Finalize the realm (the forest finalizes MPI once).
   class(ripple_cpu_object), intent(inout) :: self !< The realm.

   if (self%history_unit /= 0) close(self%history_unit)
   endsubroutine finalize_forest

   ! forest methods: staged layer
   subroutine open_step_forest(self, dt)
   !< Open a step: count it, cap dt, copy q into every stage.
   class(ripple_cpu_object), intent(inout) :: self !< The realm.
   real(R8P),                intent(in)    :: dt   !< Time step from the forest.

   call self%advance_time(dt=dt)
   call self%rk%initialize_stages(field=self%adam%field, q=self%q)
   endsubroutine open_step_forest

   subroutine begin_stage_forest(self, k, K_total, dt, realm)
   !< Begin stage `k`: publish it (`stage_active`) and build its state.
   class(ripple_cpu_object), intent(inout)                   :: self     !< The realm.
   integer(I4P),             intent(in)                      :: k        !< Stage index.
   integer(I4P),             intent(in)                      :: K_total  !< Unused: the forest's stage count.
   real(R8P),                intent(in)                      :: dt       !< Unused: the capped dt is self%dt.
   type(realm_pointer),      intent(inout), optional, target :: realm(:) !< Unused.

   self%stage_active = k
   call self%rk%compute_stage(field=self%adam%field, s=k, dt=self%dt)
   endsubroutine begin_stage_forest

   subroutine end_stage_forest(self, k, K_total, dt, realm, flux_register)
   !< End stage `k`: its residuals, the seam fluxes to the register, the stage assignment.
   class(ripple_cpu_object),    intent(inout)                   :: self          !< The realm.
   integer(I4P),                intent(in)                      :: k             !< Stage index.
   integer(I4P),                intent(in)                      :: K_total       !< Unused.
   real(R8P),                   intent(in)                      :: dt            !< Unused.
   type(realm_pointer),         intent(inout), optional, target :: realm(:)      !< Unused.
   class(flux_register_object), intent(inout), optional         :: flux_register !< Forest's flux register.

   call self%compute_residuals(q=self%rk%q_rk(:,:,:,:,:,k), s=k, flux_register=flux_register)
   call self%rk%assign_stage(field=self%adam%field, s=k, q=self%dq)
   endsubroutine end_stage_forest

   subroutine close_step_forest(self, dt)
   !< Close a step: commit q, advance the time, clear the active stage.
   class(ripple_cpu_object), intent(inout) :: self !< The realm.
   real(R8P),                intent(in)    :: dt   !< Unused: the capped dt is self%dt.

   call self%rk%update_q(field=self%adam%field, dt=self%dt, q=self%q, dq=self%dq)
   self%time = self%time + self%dt
   self%stage_active = 0_I4P
   endsubroutine close_step_forest

   subroutine apply_reflux_to_stage_forest(self, stage, dt, flux_register)
   !< Reflux the committed step (the forest calls it after `close_step_forest`, with the last stage).
   class(ripple_cpu_object),    intent(inout) :: self          !< The realm.
   integer(I4P),                intent(in)    :: stage         !< Stage: the correction acts at the last one.
   real(R8P),                   intent(in)    :: dt            !< Unused: the capped dt is self%dt.
   class(flux_register_object), intent(in)    :: flux_register !< Forest's flux register.

   if (stage == self%rk%nrk) call self%apply_reflux(flux_register=flux_register)
   endsubroutine apply_reflux_to_stage_forest
endmodule adam_ripple_cpu_object
