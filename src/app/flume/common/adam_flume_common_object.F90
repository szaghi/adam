!< ADAM, FLUME common object: data and methods shared by all backends.
module adam_flume_common_object
!< ADAM, FLUME common object: data and methods shared by all backends.
!<
!< `flume_common_object` extends the library `realm_object` (which owns grid, tree, field, maps, IO, AMR, IB, RK and
!< WENO objects) with the FLUME physics layer. The backends (`flume_cpu_object`, `flume_fnl_object`) extend it and
!< implement the forest contract. Initialization order (issue #35, section 6.1): IO, numerics, physics (decides nv),
!< blocks budget with the real per-block fields count, realm, BCs (before the first map build, so periodicity is seen
!< by the tree), time, IC, diagnostics, fields allocation, uniform refinement.

! ADAM classes, libraries, parameters
use :: adam_amr_object,               only : amr_marker_object, AMR_GEO, AMR_GEO_PRIMITIVE_BOX, AMR_GEO_STL, AMR_GRAD
use :: adam_parameters,               only : TO_BE_REFINED, TO_NOT_TOUCH
use :: adam_realm_object,             only : realm_object
use :: adam_rk_object,                only : rk_stored_stages_number
! ADAM singleton objects
use :: adam_mpih_global,              only : mpih
! FLUME modules
use :: adam_flume_bc_object,          only : flume_bc_object
use :: adam_flume_diagnostics_object, only : flume_diagnostics_object
use :: adam_flume_ic_object,          only : flume_ic_object
use :: adam_flume_numerics_object,    only : flume_numerics_object
use :: adam_flume_physics_object,     only : flume_physics_object
use :: adam_flume_time_object,        only : flume_time_object
! third party modules
use :: finer,                         only : file_ini
use :: motion,                        only : xh5f_file_object
use :: penf,                          only : I4P, I8P, R8P, str, strz
use :: stringifor,                    only : string

implicit none
private
public :: flume_common_object

character(len=11), parameter :: SCHEME_TIME_TAG="runge_kutta" !< Time-integration family tag (forest admissibility).

type, extends(realm_object) :: flume_common_object
   !< FLUME common object: data and methods shared by all backends.
   ! AMR
   logical                        :: amr_locked_=.false. !< Runtime AMR locked after initialization.
   ! fields data
   real(R8P),         allocatable :: q(:,:,:,:,:)        !< Conservative variables [nv, 1-ngc:ni+ngc, ..., nb].
   real(R8P),         allocatable :: dq(:,:,:,:,:)       !< Residuals [nv, 1-ngc:ni+ngc, ..., nb].
   real(R8P),         allocatable :: q_aux(:,:,:,:,:)    !< Auxiliary variables [nv_aux, 1-ngc:ni+ngc, ..., nb].
   type(string),      allocatable :: q_name(:)           !< Conservative variables names.
   type(string),      allocatable :: dq_name(:)          !< Residuals names.
   type(string),      allocatable :: q_aux_name(:)       !< Auxiliary variables names.
   ! FLUME classes
   type(flume_bc_object)          :: bc                  !< Boundary conditions.
   type(flume_diagnostics_object) :: diagnostics         !< Diagnostics.
   type(flume_ic_object)          :: ic                  !< Initial conditions.
   type(flume_numerics_object)    :: numerics            !< Numerics.
   type(flume_physics_object)     :: physics             !< Physics.
   type(flume_time_object)        :: time                !< Time handler.
   contains
      ! AMR methods
      procedure, pass(self) :: amr_update       !< Do AMR update (initialization-time only).
      procedure, pass(self) :: mark_by_geometry !< Mark blocks to be refined by a primitive geometric box.
      procedure, pass(self) :: mark_by_gradient !< Mark blocks by the gradient of an auxiliary variable (backend override).
      ! public methods
      procedure, pass(self) :: allocate_common       !< Allocate common data.
      procedure, pass(self) :: compute_fields_number !< Compute the block-sized fields allocated per block.
      procedure, pass(self) :: destroy_common        !< Free common data.
      procedure, pass(self) :: initialize            !< Initialize the common data.
      procedure, pass(self) :: load_restart_files    !< Load restart files.
      procedure, pass(self) :: save_restart_files    !< Save restart files.
      procedure, pass(self) :: save_xh5f             !< Save fields in XH5F format.
      ! forest methods
      procedure, pass(self) :: coupling_descriptor_forest !< Return the realm coupling descriptor.
      ! private methods
      procedure, pass(self), private :: check_ngc_number !< Check the ghost cells number against the stencils.
      procedure, pass(self), private :: io_initialize    !< Build the variables names.
endtype flume_common_object

contains
   ! AMR methods
   subroutine amr_update(self)
   !< Do AMR update: `amr%iters` sweeps over the markers until the grid stabilizes (initialization-time only).
   class(flume_common_object), intent(inout) :: self                !< The equation.
   logical                                   :: is_grid_changed     !< Flag to check grid changes for each marker.
   logical                                   :: is_grid_changed_all !< Flag to check grid changes for each iter.
   integer(I4P)                              :: i, i_marker         !< Counters.
   type(amr_marker_object)                   :: amr_marker          !< Current AMR marker.

   if (self%amr_locked_) &
      call mpih%error_stop(msg=': runtime AMR regrid is not supported, AMR is initialization-time only')
   amr: do i=1, self%amr%iters
      is_grid_changed_all = .false.
      do i_marker=1, self%amr%markers_number
         amr_marker = self%amr%markers(i_marker)
         select case(amr_marker%mode)
         case(AMR_GEO)
            select case(amr_marker%geo_type)
            case(AMR_GEO_PRIMITIVE_BOX)
               call self%mark_by_geometry(box_emin=amr_marker%box_emin, box_emax=amr_marker%box_emax, &
                                          target_level=amr_marker%target_level)
            case(AMR_GEO_STL)
               call mpih%error_stop(msg=': AMR marker geo_type STL is not supported by FLUME')
            case default
               call mpih%error_stop(msg=': AMR marker geo_type solid is not supported by FLUME yet')
            endselect
         case(AMR_GRAD)
            call self%mark_by_gradient(ivar=amr_marker%ivar, tol=amr_marker%tol)
         case default
            call mpih%error_stop(msg=': AMR marker mode '//trim(str(amr_marker%mode))//' is not supported by FLUME')
         endselect
         call self%adam%amr_update(is_marked_by_field=.true., do_blocks_reorder=.false., is_grid_changed=is_grid_changed, &
                                   q=self%q)
         is_grid_changed_all = is_grid_changed_all .or. is_grid_changed
      enddo
      if (.not.is_grid_changed_all) then
         call mpih%print_message('AMR grid stabilized after '//trim(str(i))//' AMR iterations')
         exit amr
      elseif (i == self%amr%iters) then
         call mpih%print_message('AMR grid is NOT stabilized after '//trim(str(i))//' AMR iterations')
      endif
   enddo amr
   endsubroutine amr_update

   subroutine mark_by_geometry(self, box_emin, box_emax, target_level, do_init)
   !< Mark blocks to be refined by a primitive axis-aligned box (deterministic, solution-independent).
   !<
   !< A block is flagged `TO_BE_REFINED` iff its centroid lies inside `[box_emin, box_emax]` and its refinement level
   !< is below `target_level`; every other block is left untouched (the marker is additive).
   class(flume_common_object), intent(inout)        :: self         !< The equation.
   real(R8P),                  intent(in)           :: box_emin(3)  !< Box minimum corner.
   real(R8P),                  intent(in)           :: box_emax(3)  !< Box maximum corner.
   integer(I4P),               intent(in)           :: target_level !< Refine blocks below this level.
   logical,                    intent(in), optional :: do_init      !< Re-initialize refinements queries.
   logical                                          :: do_init_     !< Re-initialize refinements queries, local var.
   real(R8P)                                        :: centroid(3)  !< Block centroid.
   integer(I4P)                                     :: b            !< Counter.

   do_init_ = .true. ; if (present(do_init)) do_init_ = do_init
   if (do_init_) self%adam%field%refinements_needed = [(TO_NOT_TOUCH, b=1, self%blocks_number)]
   associate(emin=>self%adam%field%emin, emax=>self%adam%field%emax, code=>self%adam%field%code, &
             tree=>self%adam%tree, refinements_needed=>self%adam%field%refinements_needed)
   do b=1, self%blocks_number
      centroid = 0.5_R8P * (emin(:,b) + emax(:,b))
      if (all(centroid >= box_emin) .and. all(centroid <= box_emax) .and. tree%level(code(b)) < target_level) &
         refinements_needed(b) = TO_BE_REFINED
   enddo
   endassociate
   endsubroutine mark_by_geometry

   subroutine mark_by_gradient(self, ivar, tol)
   !< Mark blocks by the gradient of an auxiliary variable: backends override, the common default is fatal.
   class(flume_common_object), intent(inout) :: self !< The equation.
   integer(I4P),               intent(in)    :: ivar !< Auxiliary variable index.
   real(R8P),                  intent(in)    :: tol  !< Refinement tolerance.

   call mpih%error_stop(msg=': mark_by_gradient is not implemented by this backend (ivar='//trim(str(ivar))// &
                            ', tol='//trim(str(tol))//')')
   endsubroutine mark_by_gradient

   ! public methods
   subroutine allocate_common(self)
   !< Allocate common data.
   class(flume_common_object), intent(inout) :: self       !< The equation.
   integer(I4P)                              :: alloc_stat !< Allocation status.
   character(999)                            :: alloc_msg  !< Allocation error message.

   associate(nv=>self%physics%nv, nv_aux=>self%physics%nv_aux, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk, &
             nb=>self%nb)
   allocate(self%q(1:nv,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb), stat=alloc_stat, errmsg=alloc_msg)
   if (alloc_stat /= 0_I4P) call mpih%error_stop(msg=': failed to allocate q: '//trim(alloc_msg))
   allocate(self%dq(1:nv,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb), stat=alloc_stat, errmsg=alloc_msg)
   if (alloc_stat /= 0_I4P) call mpih%error_stop(msg=': failed to allocate dq: '//trim(alloc_msg))
   allocate(self%q_aux(1:nv_aux,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb), stat=alloc_stat, errmsg=alloc_msg)
   if (alloc_stat /= 0_I4P) call mpih%error_stop(msg=': failed to allocate q_aux: '//trim(alloc_msg))
   endassociate
   self%q     = 0._R8P
   self%dq    = 0._R8P
   self%q_aux = 0._R8P
   endsubroutine allocate_common

   subroutine compute_fields_number(self, file_parameters, fields_number)
   !< Compute the block-sized fields FLUME allocates per block, the `fields_number` of the blocks budget.
   !<
   !< Both backends allocate the same nb-sized arrays, on the host (CPU) or on the device (FNL): `q`, `dq` (2 nv),
   !< `q_aux` (nv_aux), the three face fluxes (3 nv, counted as full block fields) and the Runge-Kutta stages
   !< (`rk_stored_stages_number` nv). Requires `physics` initialized.
   class(flume_common_object), intent(inout) :: self            !< The equation.
   type(file_ini),             intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   integer(I4P),               intent(out)   :: fields_number   !< Block-sized fields allocated per block.
   character(99)                             :: rk_scheme       !< Runge-Kutta scheme name.
   integer(I4P)                              :: stages_number   !< Runge-Kutta stage fields.
   integer(I4P)                              :: error           !< Error status.

   call file_parameters%get(section_name='runge_kutta', option_name='scheme', val=rk_scheme, error=error)
   if (error > 0) call mpih%error_stop(msg=': failed to load [runge_kutta].(scheme)')
   stages_number = rk_stored_stages_number(scheme=rk_scheme)
   if (stages_number < 0_I4P) &
      call mpih%error_stop(msg=': unknown Runge-Kutta scheme "'//trim(adjustl(rk_scheme))//'" in [runge_kutta].(scheme)')
   fields_number = self%physics%nv * (2_I4P + 3_I4P + stages_number) + self%physics%nv_aux
   endsubroutine compute_fields_number

   subroutine destroy_common(self)
   !< Free common data.
   class(flume_common_object), intent(inout) :: self !< The equation.

   if (allocated(self%q))          deallocate(self%q)
   if (allocated(self%dq))         deallocate(self%dq)
   if (allocated(self%q_aux))      deallocate(self%q_aux)
   if (allocated(self%q_name))     deallocate(self%q_name)
   if (allocated(self%dq_name))    deallocate(self%dq_name)
   if (allocated(self%q_aux_name)) deallocate(self%q_aux_name)
   endsubroutine destroy_common

   subroutine initialize(self, filename, memory_avail, nv, fields_number, verbose, L0)
   !< Initialize the common data (issue #35, section 6.1, step 3).
   class(flume_common_object), intent(inout), target :: self           !< The equation.
   character(*),               intent(in)            :: filename       !< Input file name.
   real(R8P),                  intent(in), value     :: memory_avail   !< Memory available for single MPI process.
   integer(I4P),               intent(in), optional  :: nv             !< Unused: nv is decided by the physics.
   integer(I4P),               intent(in), optional  :: fields_number  !< Block-sized fields per block (default: computed).
   logical,                    intent(in), optional  :: verbose        !< Trigger verbose output.
   real(R8P),                  intent(in), optional  :: L0             !< Unused: FLUME is dimensional.
   logical                                           :: verbose_       !< Trigger verbose output, local variable.
   integer(I4P)                                      :: fields_number_ !< Block-sized fields per block, local variable.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   call mpih%initialize(verbose=verbose_)
   if (verbose_) call mpih%print_message('flume_common_object%initialize start')
   call self%io%initialize(filename=trim(filename), verbose=verbose_)
   associate(file_parameters=>self%io%file_parameters)
   call self%numerics%initialize(file_parameters=file_parameters)
   call self%physics%initialize(file_parameters=file_parameters)
   if (present(fields_number)) then
      fields_number_ = fields_number
   else
      call self%compute_fields_number(file_parameters=file_parameters, fields_number=fields_number_)
   endif
   if (verbose_) call mpih%print_message('flume_common_object%initialize fields_number: '//trim(str(fields_number_)))
   call self%realm_object%initialize(filename=filename, memory_avail=memory_avail, nv=self%physics%nv, &
                                     fields_number=fields_number_, verbose=verbose_)
   call self%bc%initialize(file_parameters=file_parameters, physics=self%physics)
   call self%adam%grid%set_bc_type(bc_type=self%bc%bc_type)
   call self%time%initialize(file_parameters=file_parameters)
   call self%ic%initialize(file_parameters=file_parameters, physics=self%physics)
   call self%diagnostics%initialize(file_parameters=file_parameters)
   call self%check_ngc_number
   call self%allocate_common
   call self%io_initialize
   if (self%adam%tree%iu_ref_levels > 0) &
      call self%adam%refine_uniform(refinement_levels=self%adam%tree%iu_ref_levels, do_mpi_redistribute=.true., &
                                    do_blocks_reorder=.false., q=self%q)
   endassociate
   if (verbose_) call mpih%print_message('flume_common_object%initialize finish')
   endsubroutine initialize

   subroutine load_restart_files(self, t, time)
   !< Load restart files.
   class(flume_common_object), intent(inout) :: self !< The equation.
   integer(I4P),               intent(out)   :: t    !< Time iteration.
   real(R8P),                  intent(out)   :: time !< Time.

   call self%adam%load_restart_files(basename=self%io%restart_basename, t=t, time=time, q=self%q)
   call self%adam%make_comm_local_maps_ghost_bc
   endsubroutine load_restart_files

   subroutine save_restart_files(self)
   !< Save restart files.
   class(flume_common_object), intent(inout) :: self !< The equation.

   call mpih%barrier(tictoc=.true.)
   call mpih%print_message('save restart files t: '//trim(str(self%time%it, .true.))//', time: '// &
                           trim(str(self%time%time, .true.)))
   call self%adam%save_restart_files(basename=self%io%restart_basename, t=self%time%it, time=self%time%time, q=self%q)
   call self%save_xh5f(output_basename=self%io%restart_basename)
   call mpih%barrier(tictoc=.true.)
   endsubroutine save_restart_files

   subroutine save_xh5f(self, output_basename, with_ghost)
   !< Save fields in XH5F format: `q` always, `dq` when `[IO].(save_residual_fields)`.
   class(flume_common_object), intent(inout)        :: self             !< The equation.
   character(*),               intent(in), optional :: output_basename  !< Output basename.
   logical,                    intent(in), optional :: with_ghost       !< Flag to save ghost cells.
   character(:), allocatable                        :: output_basename_ !< Output basename, local var.
   logical                                          :: with_ghost_      !< Flag to save ghost cells, local var.
   type(xh5f_file_object)                           :: xh5f             !< XH5F file handler.
   integer(I4P)                                     :: ngc              !< Ghost cells saved.
   integer(I4P)                                     :: ijk(2,3)         !< Blocks extents.
   integer(I8P)                                     :: nijk(3)          !< Blocks dimensions.
   character(:), allocatable                        :: bn               !< Block name.
   integer(I4P)                                     :: b                !< Counter.

   call mpih%barrier(tictoc=.true.)
   call mpih%print_message('save HDF5 files t: '//trim(str(self%time%it, .true.))//', time: '// &
                           trim(str(self%time%time, .true.)))
   output_basename_ = trim(self%io%output_basename)//'-'//trim(strz(self%time%it, 9))
   if (present(output_basename)) output_basename_ = trim(output_basename)
   with_ghost_ = .false. ; if (present(with_ghost)) with_ghost_ = with_ghost
   ngc = 0_I4P ; if (with_ghost_) ngc = self%adam%grid%ngc
   associate(ni=>self%adam%grid%ni, nj=>self%adam%grid%nj, nk=>self%adam%grid%nk)
   ijk(:,1) = [1-ngc, ni+ngc]
   ijk(:,2) = [1-ngc, nj+ngc]
   ijk(:,3) = [1-ngc, nk+ngc]
   nijk = [ijk(2,1)-ijk(1,1)+1, ijk(2,2)-ijk(1,2)+1, ijk(2,3)-ijk(1,3)+1]
   endassociate
   call self%open_file_xh5f(basename=trim(output_basename_), xh5f=xh5f)
   do b=1, self%adam%field%blocks_number
      bn = 'block_'//trim(strz(b, 9))//'-proc'//trim(strz(mpih%myrank, 6))
      call self%open_block_xh5f(xh5f=xh5f, b=b, nijk=nijk, t=self%time%it, time=self%time%time)
      call self%io%save_field(xh5f=xh5f, grid=self%adam%grid, block_name=bn, ijk=ijk, nijk=nijk, &
                              q=self%q(:,:,:,:,b), q_name=self%q_name)
      if (self%io%save_residual_fields) &
         call self%io%save_field(xh5f=xh5f, grid=self%adam%grid, block_name=bn, ijk=ijk, nijk=nijk, &
                                 q=self%dq(:,:,:,:,b), q_name=self%dq_name)
      call self%close_block_xh5f(xh5f=xh5f)
   enddo
   call self%close_file_xh5f(xh5f=xh5f)
   call mpih%barrier(tictoc=.true.)
   endsubroutine save_xh5f

   ! forest methods
   subroutine coupling_descriptor_forest(self, scheme_time, rk_scheme, nv)
   !< Return the realm coupling descriptor checked by the forest for stage-coincident admissibility.
   class(flume_common_object), intent(in)  :: self        !< The equation.
   character(:), allocatable,  intent(out) :: scheme_time !< Time-integration family tag.
   character(:), allocatable,  intent(out) :: rk_scheme   !< Within-family scheme tag.
   integer(I4P),               intent(out) :: nv          !< Number of conserved variables on this realm.

   scheme_time = SCHEME_TIME_TAG
   rk_scheme   = trim(self%rk%scheme)
   nv          = self%physics%nv
   endsubroutine coupling_descriptor_forest

   ! private methods
   subroutine check_ngc_number(self)
   !< Check the ghost cells number against the WENO stencil half-width.
   class(flume_common_object), intent(in) :: self !< The equation.

   if (self%weno%S > self%ngc) &
      call mpih%error_stop(msg=': [grid].(ngc)='//trim(str(self%ngc))//' is smaller than the WENO stencil half-width '// &
                               trim(str(self%weno%S)))
   endsubroutine check_ngc_number

   subroutine io_initialize(self)
   !< Build the variables names from the physical model (the same predicate that decided nv).
   class(flume_common_object), intent(inout) :: self !< The equation.
   integer(I4P)                              :: v    !< Counter.

   allocate(self%q_name(1:self%physics%nv), self%dq_name(1:self%physics%nv), self%q_aux_name(1:self%physics%nv_aux))
   self%q_name(1) = 'r'
   self%q_name(2) = 'ru'
   self%q_name(3) = 'rv'
   self%q_name(4) = 'rw'
   self%q_name(5) = 'rE'
   do v=1, self%physics%nv
      self%dq_name(v) = 'dq_'//self%q_name(v)%chars()
   enddo
   self%q_aux_name(1) = 'rho'
   self%q_aux_name(2) = 'u'
   self%q_aux_name(3) = 'v'
   self%q_aux_name(4) = 'w'
   self%q_aux_name(5) = 'p'
   self%q_aux_name(6) = 'T'
   self%q_aux_name(7) = 'H'
   self%q_aux_name(8) = 'a'
   endsubroutine io_initialize
endmodule adam_flume_common_object
