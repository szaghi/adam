!< RIPPLE, the common realm: what the CPU and FNL backends share.
module adam_ripple_common_object
!< RIPPLE, the common realm: what the CPU and FNL backends share.
!<
!< RIPPLE solves the scalar linear advection `u_t + a . grad(u) = 0` on a periodic box: the worked example of the
!< "build your own solver" tutorial (issue #59). A RIPPLE realm extends the library's `realm_object`, which already owns
!< the grid, the tree, the field, the maps, the WENO and the Runge-Kutta objects; this module adds the state `q`, the
!< advection velocity, the time bookkeeping, the initial condition (a smooth profile with an exact translated solution),
!< the initial AMR and the parts of the forest contract that do not depend on where `q` lives.
!<
!< The scheme is the Shu-Osher conservative finite difference: point values at the cell centres, the upwind split flux
!< `a+ u` / `a- u` reconstructed at every face by the library's WENO, `dq = -sum_d (F_{i+1/2} - F_{i-1/2}) / dx_d`. The
!< backends own the residuals and the integration.

! ADAM classes, libraries, parameters
use :: adam_amr_criterion_object, only : amr_criterion_object, amr_criterion_pointer
use :: adam_amr_object,           only : amr_marker_object, AMR_GEO, AMR_GEO_PRIMITIVE_BOX, AMR_USER
use :: adam_flux_register_object, only : flux_register_object, face_tangential_ratios, restrict_fine_face_to_quadrant, &
                                         SEAM_KIND_INTRA_REALM_AMR
use :: adam_maps_object,          only : face_axis_sign
use :: adam_parameters,           only : BC_PERIODIC, TO_BE_DEREFINED, TO_BE_REFINED, TO_NOT_TOUCH
use :: adam_realm_object,         only : realm_object, realm_pointer
use :: adam_rk_object,            only : rk_stored_stages_number, RK_SSP_11, RK_SSP_22, RK_SSP_33, RK_SSP_54
use :: adam_weno_object,          only : WENO_WEIGHTS_JS
! ADAM singleton objects
use :: adam_mpih_global,          only : mpih
! third party modules
use :: mpi
use :: penf

implicit none
private
public :: ripple_common_object
public :: exact_solution

real(R8P), parameter :: PI2 = 2._R8P * acos(-1._R8P) !< 2 pi.

type, extends(realm_object) :: ripple_common_object
   !< RIPPLE common realm.
   real(R8P), allocatable :: q(:,:,:,:,:)            !< State, point values u [1, 1-ngc:ni+ngc, ..., nb].
   real(R8P), allocatable :: dq(:,:,:,:,:)           !< Residuals, same shape.
   real(R8P)              :: a(3)=0._R8P             !< Advection velocity.
   real(R8P)              :: cfl=0.5_R8P             !< CFL number.
   real(R8P)              :: time=0._R8P             !< Current time.
   real(R8P)              :: time_max=0._R8P         !< Final time.
   real(R8P)              :: dt=0._R8P               !< Current time step.
   integer(I4P)           :: it=0_I4P                !< Current iteration.
   integer(I4P)           :: it_max=0_I4P            !< Iterations number; <= 0: run to time_max.
   integer(I4P)           :: history_unit=0_I4P      !< Unit of the history file (rank 0).
   type(amr_criterion_pointer), allocatable :: amr_criteria(:) !< Registered refinement criteria (issue #110).
   contains
      ! common methods
      procedure, pass(self) :: initialize_ripple        !< Initialize the common data.
      procedure, pass(self) :: check_input              !< Refuse the inputs RIPPLE cannot honour.
      procedure, pass(self) :: set_initial_conditions   !< Set the initial profile on the interior cells.
      procedure, pass(self) :: amr_update               !< Initial AMR: box markers and registered criteria.
      procedure, pass(self) :: register_amr_criterion   !< Register a refinement criterion by name (issue #110).
      procedure, pass(self) :: compute_total            !< Total of u dV on this rank.
      procedure, pass(self) :: save_output              !< Save one history row and the state (XH5F).
      procedure, pass(self) :: advance_time             !< Count the step and cap dt on time_max.
      procedure, pass(self) :: accumulate_seam_skin     !< Route one seam face skin to the flux register.
      procedure, pass(self) :: apply_reflux             !< Berger-Colella correction of the coarse seam cells of q.
      ! forest methods (backend-independent)
      procedure, pass(self) :: compute_local_dt_forest  !< Local CFL time step.
      procedure, pass(self) :: is_done_forest           !< True when the run is over.
      procedure, pass(self) :: stages_per_step_forest   !< Runge-Kutta stages per step (staged path).
endtype ripple_common_object

contains
   ! common methods
   subroutine initialize_ripple(self, filename, memory_avail)
   !< Initialize the common data: the library objects (`realm_object%initialize`), periodic faces, `[ripple]` and
   !< `[time]`, the state, the uniform refinement. MPI is initialized by the backend, before.
   class(ripple_common_object), intent(inout) :: self          !< The realm.
   character(*),                intent(in)    :: filename      !< Input file name.
   real(R8P),                   intent(in)    :: memory_avail  !< Memory (GB) the blocks budget may use.
   integer(I4P)                               :: stages        !< Runge-Kutta stage fields.
   character(99)                              :: rk_scheme     !< Runge-Kutta scheme name.
   integer(I4P)                               :: error         !< Error status.

   call self%io%initialize(filename=trim(filename), verbose=.true.)
   ! the blocks budget counts the block-sized fields per block: q, dq, the three face fluxes, the stored stages
   call self%io%file_parameters%get(section_name='runge_kutta', option_name='scheme', val=rk_scheme, error=error)
   stages = rk_stored_stages_number(scheme=rk_scheme)
   if (error > 0 .or. stages < 0_I4P) call mpih%error_stop(msg=': failed to load [runge_kutta].(scheme)')
   call self%realm_object%initialize(filename=filename, memory_avail=memory_avail, nv=1_I4P, &
                                     fields_number=5_I4P + stages, verbose=.true.)
   ! every face periodic; the tree and the maps were built before the faces were known: rebuild them (issue #103)
   call self%adam%grid%set_bc_type(bc_type=[BC_PERIODIC, BC_PERIODIC, BC_PERIODIC, BC_PERIODIC, BC_PERIODIC, BC_PERIODIC])
   call self%adam%tree%rebuild_neighborhood(grid=self%adam%grid)
   call self%adam%make_comm_local_maps_ghost_bc
   associate(fp=>self%io%file_parameters)
      call fp%get(section_name='ripple', option_name='velocity_x', val=self%a(1), error=error) ; call check('velocity_x')
      call fp%get(section_name='ripple', option_name='velocity_y', val=self%a(2), error=error) ; call check('velocity_y')
      call fp%get(section_name='ripple', option_name='velocity_z', val=self%a(3), error=error) ; call check('velocity_z')
      call fp%get(section_name='time', option_name='time_max', val=self%time_max, error=error) ; call check('time_max')
      call fp%get(section_name='time', option_name='it_max', val=self%it_max, error=error)     ; call check('it_max')
      call fp%get(section_name='time', option_name='CFL', val=self%cfl, error=error)           ; call check('CFL')
   endassociate
   call self%check_input
   where (self%adam%grid%null_xyz) self%a = 0._R8P
   associate(ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk, nb=>self%nb)
      allocate(self%q(1,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,nb), source=0._R8P)
      allocate(self%dq(1,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,nb), source=0._R8P)
   endassociate
   if (self%adam%tree%iu_ref_levels > 0) &
      call self%adam%refine_uniform(refinement_levels=self%adam%tree%iu_ref_levels, do_mpi_redistribute=.true., &
                                    do_blocks_reorder=.false., q=self%q)
   contains
      subroutine check(key)
      !< Stop on a missing key.
      character(*), intent(in) :: key !< Key name.

      if (error > 0) call mpih%error_stop(msg=': failed to load the RIPPLE input key "'//key//'"')
      endsubroutine check
   endsubroutine initialize_ripple

   subroutine check_input(self)
   !< Refuse the inputs RIPPLE cannot honour, at initialization: a library key that the app silently ignored would run
   !< a different scheme than the input asks for.
   class(ripple_common_object), intent(in) :: self !< The realm.
   integer(I4P)                            :: m, c !< Counters.

   select case(self%rk%scheme)
   case(RK_SSP_11, RK_SSP_22, RK_SSP_33, RK_SSP_54)
   case default
      call mpih%error_stop(msg=': [runge_kutta].(scheme) must be an SSP scheme (runge-kutta-ssp-11/22/33/54): the '// &
                               'step combines the stages with the SSP alpha/beta tableau')
   endselect
   if (self%weno%is_centered) call mpih%error_stop(msg=': [weno].(scheme) must be upwind (weno-u-*): the split '// &
                                                       'fluxes a+ u and a- u are upwinded')
   if (self%weno%weights /= WENO_WEIGHTS_JS) call mpih%error_stop(msg=': [weno].(weights) must be js')
   if (self%weno%S > self%ngc) call mpih%error_stop(msg=': [grid].(ngc)='//trim(str(self%ngc, .true.))// &
                                                    ' is smaller than the WENO stencil half-width '//trim(str(self%weno%S, .true.)))
   if (self%amr%frequency > 0_I4P) call mpih%error_stop(msg=': [amr].(frequency) must be 0: RIPPLE refines at '// &
                                                            'initialization only')
   if (self%ib%solids_number > 0_I4P) call mpih%error_stop(msg=': [solids].(number) must be 0: RIPPLE has no '// &
                                                               'immersed boundary')
   do m=1, self%amr%markers_number
      associate(marker=>self%amr%markers(m))
      if (marker%mode == AMR_GEO .and. marker%geo_type == AMR_GEO_PRIMITIVE_BOX) cycle
      if (marker%mode /= AMR_USER) call mpih%error_stop(msg=': RIPPLE supports the primitive-box AMR marker and the '// &
                                                            'registered criteria (mode 6) only')
      c = criterion_index(marker%criterion)
      if (c == 0) call mpih%error_stop(msg=': [amr_marker_'//trim(str(m, .true.))//'] name = "'//marker%criterion//'"'// &
                                           ' is not a registered criterion: register it before the run')
      if (self%amr_criteria(c)%ptr%runtime_only) &
         call mpih%error_stop(msg=': the criterion "'//marker%criterion//'" acts at runtime regrids only, and RIPPLE '// &
                                  'refines at initialization only: set its runtime_only to .false.')
      endassociate
   enddo
   contains
      function criterion_index(name) result(c)
      !< The index of the criterion registered under `name`, 0 if none.
      character(*), intent(in) :: name !< Criterion name.
      integer(I4P)             :: c    !< Its index.

      if (allocated(self%amr_criteria)) then
         do c=1, size(self%amr_criteria, kind=I4P)
            if (self%amr_criteria(c)%ptr%name == name) return
         enddo
      endif
      c = 0_I4P
      endfunction criterion_index
   endsubroutine check_input

   subroutine set_initial_conditions(self)
   !< Set the initial profile `u0` (`exact_solution` at t = 0) on the interior cells of every block.
   class(ripple_common_object), intent(inout) :: self !< The realm.
   integer(I4P)                               :: b, i, j, k !< Counters.

   associate(emin=>self%adam%field%emin, dxyz=>self%adam%field%dxyz, grid=>self%adam%grid)
   do b=1, self%blocks_number
      do k=1, self%nk
         do j=1, self%nj
            do i=1, self%ni
               self%q(1,i,j,k,b) = exact_solution(x=emin(:,b) + ([i, j, k] - 0.5_R8P) * dxyz(:,b), t=0._R8P, a=self%a, &
                                                  emin=grid%domain_emin, emax=grid%domain_emax, is_null=grid%null_xyz)
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine set_initial_conditions

   subroutine amr_update(self)
   !< Initial AMR, `[amr] iters` sweeps over the `[amr_marker_N]` sections: a box (mode 1, `primitive-box`) refines the
   !< blocks whose centroid lies in it until they reach its `target_level`; a registered criterion (mode 6, issue #110)
   !< marks with its own votes (`mark` with `initial = .true.`). Other markers are refused by `check_input`.
   class(ripple_common_object), intent(inout) :: self     !< The realm.
   type(amr_marker_object)                    :: marker   !< Current marker.
   real(R8P)                                  :: c(3)     !< Block centroid.
   logical                                    :: changed  !< The grid changed.
   integer(I4P), allocatable                  :: votes(:) !< A criterion's votes, one per local block.
   integer(I4P)                               :: it, m, b, k !< Counters.

   do it=1, self%amr%iters
      do m=1, self%amr%markers_number
         marker = self%amr%markers(m)
         associate(field=>self%adam%field, tree=>self%adam%tree)
            field%refinements_needed = [(TO_NOT_TOUCH, b=1, self%blocks_number)]
            if (marker%mode == AMR_USER) then
               ! a registered criterion: its votes, checked
               do k=1, size(self%amr_criteria, kind=I4P)
                  if (self%amr_criteria(k)%ptr%name == marker%criterion) exit
               enddo
               allocate(votes(self%blocks_number))
               call self%amr_criteria(k)%ptr%mark(realm=self, votes=votes, initial=.true.)
               if (any(votes /= TO_BE_REFINED .and. votes /= TO_NOT_TOUCH .and. votes /= TO_BE_DEREFINED)) &
                  call mpih%error_stop(msg=': the criterion "'//marker%criterion//'" voted something that is not '// &
                                           '-1, 0 or 1')
               field%refinements_needed(1:self%blocks_number) = votes
               deallocate(votes)
            else
               do b=1, self%blocks_number
                  c = 0.5_R8P * (field%emin(:,b) + field%emax(:,b))
                  if (all(c >= marker%box_emin) .and. all(c <= marker%box_emax) .and. &
                      tree%level(field%code(b)) < marker%target_level) field%refinements_needed(b) = TO_BE_REFINED
               enddo
            endif
         endassociate
         call self%adam%amr_update(is_marked_by_field=.true., do_blocks_reorder=.false., is_grid_changed=changed, &
                                   q=self%q)
      enddo
   enddo
   call self%adam%make_comm_local_maps_ghost_bc
   endsubroutine amr_update

   subroutine register_amr_criterion(self, criterion)
   !< Register a refinement criterion (issue #110): `[amr_marker_N] mode = 6, name = <criterion%name>` selects it. RIPPLE
   !< refines at initialization only, so the criterion must act there (`runtime_only = .false.`). The realm keeps a
   !< pointer: the program owns the instance (a TARGET of the main program, outliving the run).
   class(ripple_common_object),          intent(inout) :: self      !< The realm.
   class(amr_criterion_object), pointer, intent(in)    :: criterion !< The criterion (a TARGET actual is pointed to).
   type(amr_criterion_pointer), allocatable           :: grown(:)  !< The list, one longer.
   integer(I4P)                                       :: c         !< Counter.

   if (.not.allocated(criterion%name)) call mpih%error_stop(msg=': register_amr_criterion: the criterion has no name')
   if (.not.allocated(self%amr_criteria)) allocate(self%amr_criteria(0))
   do c=1, size(self%amr_criteria, kind=I4P)
      if (self%amr_criteria(c)%ptr%name == criterion%name) &
         call mpih%error_stop(msg=': register_amr_criterion: "'//criterion%name//'" is registered twice')
   enddo
   allocate(grown(size(self%amr_criteria) + 1))
   grown(1:size(self%amr_criteria)) = self%amr_criteria
   grown(size(grown))%ptr => criterion
   call move_alloc(grown, self%amr_criteria)
   endsubroutine register_amr_criterion

   function compute_total(self) result(total)
   !< Return the total of u dV over the interior cells of this rank's blocks (host `q`).
   class(ripple_common_object), intent(in) :: self  !< The realm.
   real(R8P)                               :: total !< Total of u dV.
   integer(I4P)                            :: b     !< Counter.

   total = 0._R8P
   do b=1, self%blocks_number
      total = total + sum(self%q(1,1:self%ni,1:self%nj,1:self%nk,b)) * product(self%adam%field%dxyz(:,b), &
                                                                               mask=.not.self%adam%grid%null_xyz)
   enddo
   endfunction compute_total

   subroutine save_output(self, fields)
   !< Append one row (it, time, total of u dV over the domain) to the history file and, when `fields`, save the host `q`
   !< in XH5F (one HDF5 file per rank, interior cells, `<output_basename>-<it>`).
   class(ripple_common_object), intent(inout) :: self   !< The realm.
   logical,                     intent(in)    :: fields !< Save the fields too.
   real(R8P)                                  :: total  !< Total of u dV.
   integer                                    :: ierr   !< MPI error.

   total = self%compute_total()
   call MPI_ALLREDUCE(MPI_IN_PLACE, total, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
   if (mpih%myrank == 0) then
      if (self%history_unit == 0) then
         open(newunit=self%history_unit, file=trim(self%io%output_basename)//'-history.dat', action='write')
         write(self%history_unit, '(A)') '# it time total(u dV)'
      endif
      write(self%history_unit, '(I9,2(1X,ES23.15E3))') self%it, self%time, total
      flush(self%history_unit)
   endif
   if (fields) call self%save_q_xh5f(basename=trim(self%io%output_basename)//'-'//trim(strz(self%it, 9)), q=self%q, &
                                     q_name=['u'], t=self%it, time=self%time)
   endsubroutine save_output

   subroutine advance_time(self, dt)
   !< Count the step and take the forest's dt, capped so that a time-driven run lands on `time_max`.
   class(ripple_common_object), intent(inout) :: self !< The realm.
   real(R8P),                   intent(in)    :: dt   !< Time step from the forest.

   self%it = self%it + 1_I4P
   self%dt = dt
   if (self%it_max <= 0_I4P .and. self%time + dt > self%time_max) self%dt = self%time_max - self%time
   endsubroutine advance_time

   subroutine accumulate_seam_skin(self, flux_register, b, fec, weight, skin)
   !< Route one face skin of block `b`, face `fec`, weighted by `weight`, to the forest's flux register.
   !<
   !< `skin(1, c)` is the face flux with `c` running over the two tangential axes, inner fastest (x faces: j, k; y faces:
   !< i, k; z faces: i, j). A positive register index marks the coarse side of a 2:1 face: the skin is the coarse face. A
   !< negative one marks a fine block: its skin is restricted (2x2 average, 2x1 on a quadtree) into its quadrant of the
   !< coarse face. Weighting every stage by its SSP coefficient makes the register hold `sum_s beta_s F_s`, the flux of
   !< the committed step.
   class(ripple_common_object), intent(in)    :: self          !< The realm.
   class(flux_register_object), intent(inout) :: flux_register !< Forest's flux register.
   integer(I4P),                intent(in)    :: b, fec        !< Block, face (1..6: -x, +x, -y, +y, -z, +z).
   real(R8P),                   intent(in)    :: weight        !< Stage weight.
   real(R8P),                   intent(in)    :: skin(1:,1:)   !< Face skin (1, inner_n*outer_n).
   real(R8P), allocatable                     :: slab(:,:)     !< Register-shaped contribution.
   integer(I4P)                               :: idx, n(2)     !< Signed register index, tangential cell counts.
   integer(I4P)                               :: ratios(2)     !< Tangential refinement ratios.

   idx = self%adam%maps%inter_realm_face_register_index(b, fec)
   if (idx == 0_I4P) return
   if (flux_register%face(abs(idx))%seam_kind /= SEAM_KIND_INTRA_REALM_AMR) &
      call mpih%error_stop(msg=': RIPPLE is single-realm: an inter-realm register face is not expected')
   select case(fec)
   case(1_I4P, 2_I4P) ; n = [self%nj, self%nk]
   case(3_I4P, 4_I4P) ; n = [self%ni, self%nk]
   case default       ; n = [self%ni, self%nj]
   endselect
   allocate(slab(1, n(1)*n(2)), source=0._R8P)
   if (idx > 0_I4P) then
      slab = weight * skin
      call flux_register%accumulate_coarse_flux(face_index=idx, stage=1_I4P, flux_face=slab)
   else
      ratios = face_tangential_ratios(fec=fec, refine_ratio=self%adam%maps%refine_ratio)
      call restrict_fine_face_to_quadrant(fine_face=weight * reshape(skin, [1, n(1), n(2)]), inner_n=n(1), outer_n=n(2), &
                                          ioff=self%adam%maps%amr_seam_quadrant(1, b, fec),                            &
                                          joff=self%adam%maps%amr_seam_quadrant(2, b, fec),                            &
                                          slab=slab, inner_ratio=ratios(1), outer_ratio=ratios(2))
      call flux_register%accumulate_fine_flux(face_index=-idx, stage=1_I4P, flux_face=slab)
   endif
   endsubroutine accumulate_seam_skin

   subroutine apply_reflux(self, flux_register)
   !< Berger-Colella reflux on the host `q`: every coarse seam cell this rank owns receives
   !< `sgn dt / dx (F_coarse - F_fine)`, the coarse flux the step used replaced by the restricted fine one, so that the
   !< total of u dV is conserved across the 2:1 faces to round-off.
   class(ripple_common_object), intent(inout) :: self          !< The realm.
   class(flux_register_object), intent(in)    :: flux_register !< Forest's flux register.
   integer(I4P)                               :: f, c, axis, sgn, i, j, k !< Counters, face normal and side.

   if (flux_register%nfaces == 0_I4P .or. .not.allocated(flux_register%face)) return
   do f=1, flux_register%nfaces
      associate(face=>flux_register%face(f))
      if (face%coarse_realm /= self%realm_index .or. face%coarse_rank /= mpih%myrank) cycle
      call face_axis_sign(face_code=face%coarse_face, axis=axis, sgn=sgn)
      do c=1, face%nface_cells
         select case(axis)
         case(1_I4P) ; i = merge(self%ni, 1_I4P, sgn > 0) ; j = 1 + mod(c - 1, self%nj) ; k = 1 + (c - 1) / self%nj
         case(2_I4P) ; i = 1 + mod(c - 1, self%ni) ; j = merge(self%nj, 1_I4P, sgn > 0) ; k = 1 + (c - 1) / self%ni
         case default; i = 1 + mod(c - 1, self%ni) ; j = 1 + (c - 1) / self%ni ; k = merge(self%nk, 1_I4P, sgn > 0)
         endselect
         self%q(1,i,j,k,face%coarse_block) = self%q(1,i,j,k,face%coarse_block) +                              &
            sgn * self%dt / self%adam%field%dxyz(axis,face%coarse_block) * (face%F_coarse(1,c,1) - face%F_fine_sum(1,c,1))
      enddo
      endassociate
   enddo
   endsubroutine apply_reflux

   ! forest methods
   subroutine compute_local_dt_forest(self, dt_local)
   !< Return the local CFL time step, `CFL / max_b sum_d |a_d| / dx_d` (the forest takes the minimum over the ranks).
   class(ripple_common_object), intent(in)  :: self     !< The realm.
   real(R8P),                   intent(out) :: dt_local !< Local time step.
   integer(I4P)                             :: b        !< Counter.

   dt_local = huge(1._R8P)
   do b=1, self%blocks_number
      dt_local = min(dt_local, self%cfl / max(sum(abs(self%a) / self%adam%field%dxyz(:,b)), tiny(1._R8P)))
   enddo
   endsubroutine compute_local_dt_forest

   subroutine is_done_forest(self, done)
   !< Return true when `it_max` steps are done or, with `it_max <= 0`, when `time_max` is reached.
   class(ripple_common_object), intent(in)  :: self !< The realm.
   logical,                     intent(out) :: done !< True if the realm is done.

   if (self%it_max > 0_I4P) then
      done = self%it >= self%it_max
   else
      done = self%time >= self%time_max * (1._R8P - 10._R8P * epsilon(1._R8P))
   endif
   endsubroutine is_done_forest

   function stages_per_step_forest(self) result(K)
   !< Return the Runge-Kutta stages per step: the staged path (2:1 AMR faces) needs an SSP scheme.
   class(ripple_common_object), intent(in) :: self !< The realm.
   integer(I4P)                            :: K    !< Stages per step.

   select case(self%rk%scheme)
   case(RK_SSP_11, RK_SSP_22, RK_SSP_33, RK_SSP_54)
      K = self%rk%nrk
   case default
      K = 0_I4P
      call mpih%error_stop(msg=': the staged path (2:1 AMR faces) needs an SSP Runge-Kutta scheme')
   endselect
   endfunction stages_per_step_forest

   ! non TBP
   pure function exact_solution(x, t, a, emin, emax, is_null) result(u)
   !< Return the exact solution at point `x`, time `t`: the profile `u0 = prod_d sin(2 pi (x_d - emin_d) / L_d)` over
   !< the non-null axes, translated by `a t` (periodic box).
   real(R8P), intent(in) :: x(3)          !< Point.
   real(R8P), intent(in) :: t             !< Time.
   real(R8P), intent(in) :: a(3)          !< Advection velocity.
   real(R8P), intent(in) :: emin(3)       !< Domain minimum corner.
   real(R8P), intent(in) :: emax(3)       !< Domain maximum corner.
   logical,   intent(in) :: is_null(3)    !< Null axes.
   real(R8P)             :: u             !< Solution.
   integer(I4P)          :: d             !< Counter.

   u = 1._R8P
   do d=1, 3
      if (.not.is_null(d)) u = u * sin(PI2 * (x(d) - a(d) * t - emin(d)) / (emax(d) - emin(d)))
   enddo
   endfunction exact_solution
endmodule adam_ripple_common_object
