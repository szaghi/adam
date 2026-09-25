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
use :: adam_amr_object,               only : amr_marker_object, AMR_DELTA_T_MAX, AMR_DELTA_T_X, AMR_DELTA_T_Y, AMR_DELTA_T_Z, &
                                             AMR_GEO, AMR_GEO_PRIMITIVE_BOX, AMR_GEO_SOLID, AMR_GEO_STL, AMR_GRAD
use :: adam_flux_register_object,     only : flux_register_object, restrict_fine_face_to_quadrant, SEAM_KIND_INTER_REALM
use :: adam_parameters,               only : TO_BE_DEREFINED, TO_BE_REFINED, TO_NOT_TOUCH
use :: adam_realm_object,             only : realm_object
use :: adam_rk_object,                only : rk_stored_stages_number
! ADAM singleton objects
use :: adam_mpih_global,              only : mpih
! FLUME modules
use :: adam_flume_bc_object,          only : flume_bc_object
use :: adam_flume_diagnostics_object, only : flume_diagnostics_object
use :: adam_flume_euler_library,      only : conservative_to_auxiliary
use :: adam_flume_ic_object,          only : flume_ic_object
use :: adam_flume_numerics_object,    only : flume_numerics_object
use :: adam_flume_parameters,         only : NV_AUX
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
public :: ib_cut_spacing
public :: seam_skin_cell

character(len=11), parameter :: SCHEME_TIME_TAG="runge_kutta" !< Time-integration family tag (forest admissibility).

type, extends(realm_object) :: flume_common_object
   !< FLUME common object: data and methods shared by all backends.
   ! AMR
   logical                        :: amr_locked_=.false. !< Runtime AMR locked after initialization.
   ! IO
   logical                        :: save_auxiliary_fields=.false. !< Save the auxiliary variables with the fields.
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
      procedure, pass(self) :: mark_by_gradient !< Mark blocks by the gradient of a conservative or auxiliary variable.
      procedure, pass(self) :: mark_by_solid    !< Mark blocks crossed by the surface of an immersed solid.
      ! public methods
      procedure, pass(self) :: accumulate_seam_skin  !< Route one weighted seam face skin to the forest's flux register.
      procedure, pass(self) :: allocate_common       !< Allocate common data.
      procedure, pass(self) :: compute_fields_number !< Compute the block-sized fields allocated per block.
      procedure, pass(self) :: compute_phi           !< Compute the immersed solids distance function (host).
      procedure, pass(self) :: destroy_common        !< Free common data.
      procedure, pass(self) :: initialize            !< Initialize the common data.
      procedure, pass(self) :: load_restart_files    !< Load restart files.
      procedure, pass(self) :: save_restart_files    !< Save restart files.
      procedure, pass(self) :: save_slices           !< Save the slices on their cadence.
      procedure, pass(self) :: save_xh5f             !< Save fields in XH5F format.
      ! forest methods
      procedure, pass(self) :: coupling_descriptor_forest !< Return the realm coupling descriptor.
      ! private methods
      procedure, pass(self), private :: block_spacing    !< Return the spacing of a block by a delta criterion.
      procedure, pass(self), private :: check_amr_block_cells !< Check the block cells numbers against the 2:1 refinement.
      procedure, pass(self), private :: check_ngc_number      !< Check the ghost cells number against the stencils.
      procedure, pass(self), private :: check_slices     !< Check the slices interpolation types.
      procedure, pass(self), private :: compute_q_aux_host !< Compute the auxiliary variables of the host q.
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
            case(AMR_GEO_SOLID)
               call self%compute_phi
               call self%mark_by_solid(solid=amr_marker%solid, delta_type=amr_marker%delta_type, &
                                       delta_fine=amr_marker%delta_fine, delta_coarse=amr_marker%delta_coarse)
            case(AMR_GEO_STL)
               call mpih%error_stop(msg=': AMR marker geo_type STL is not supported by FLUME')
            case default
               call mpih%error_stop(msg=': unknown AMR marker geo_type '//trim(str(amr_marker%geo_type)))
            endselect
         case(AMR_GRAD)
            call self%mark_by_gradient(field=amr_marker%field, ivar=amr_marker%ivar, tol=amr_marker%tol,           &
                                       delta_type=amr_marker%delta_type, delta_fine=amr_marker%delta_fine, &
                                       delta_coarse=amr_marker%delta_coarse)
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

   subroutine mark_by_gradient(self, field, ivar, tol, delta_type, delta_fine, delta_coarse)
   !< Mark blocks by the gradient magnitude of a conservative (`field = 1`) or auxiliary (`field = 2`) variable.
   !<
   !< CHASE semantics: the admissible cell spacing of a block is `delta_fine` where `max |grad var| > tol`, else
   !< `delta_coarse`; a block coarser than admissible is refined, a block whose parent (spacing doubled) would still be
   !< admissible is derefined, any other block is left untouched. The spacing of a block is chosen by `delta_type`
   !< (`x`, `y`, `z` or `max`). The gradient uses the interior cells only (centred differences, one-sided at the block
   !< edges), so the marker does not depend on the ghost cells: it runs on the host state before any ghost exchange,
   !< for both backends (AMR is initialization-time only).
   class(flume_common_object), intent(inout) :: self         !< The equation.
   integer(I4P),               intent(in)    :: field        !< Marker field: 1 conservative, 2 auxiliary variables.
   integer(I4P),               intent(in)    :: ivar         !< Variable index in the marker field.
   real(R8P),                  intent(in)    :: tol          !< Gradient magnitude tolerance.
   character(*),               intent(in)    :: delta_type   !< Block spacing criterion: x, y, z, max.
   real(R8P),                  intent(in)    :: delta_fine   !< Admissible spacing where the gradient exceeds tol.
   real(R8P),                  intent(in)    :: delta_coarse !< Admissible spacing elsewhere.
   real(R8P), allocatable                    :: var(:,:,:)   !< Marker variable of one block, interior cells.
   real(R8P)                                 :: qa(NV_AUX)   !< Auxiliary variables of one cell.
   real(R8P)                                 :: grad(3)      !< Gradient of one cell.
   real(R8P)                                 :: grad_max     !< Maximum gradient magnitude of one block.
   real(R8P)                                 :: dc           !< Block spacing.
   real(R8P)                                 :: delta        !< Admissible spacing.
   integer(I4P)                              :: b, i, j, k   !< Counters.

   if ((field == 1_I4P .and. (ivar < 1_I4P .or. ivar > self%physics%nv))     .or. &
       (field == 2_I4P .and. (ivar < 1_I4P .or. ivar > self%physics%nv_aux)) .or. (field < 1_I4P .or. field > 2_I4P)) &
      call mpih%error_stop(msg=': AMR gradient marker: invalid field '//trim(str(field))//' / ivar '//trim(str(ivar)))
   self%adam%field%refinements_needed = [(TO_NOT_TOUCH, b=1, self%blocks_number)]
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, dxyz=>self%adam%field%dxyz, is_null=>self%adam%grid%null_xyz, &
             refinements_needed=>self%adam%field%refinements_needed)
   allocate(var(ni,nj,nk))
   do b=1, self%blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               if (field == 1_I4P) then
                  var(i,j,k) = self%q(ivar,i,j,k,b)
               else
                  call conservative_to_auxiliary(gamma=self%physics%gamma, R=self%physics%R, q=self%q(:,i,j,k,b), qa=qa)
                  var(i,j,k) = qa(ivar)
               endif
            enddo
         enddo
      enddo
      grad_max = 0._R8P
      do k=1, nk
         do j=1, nj
            do i=1, ni
               grad = 0._R8P
               if (.not.is_null(1) .and. ni > 1) grad(1) = interior_derivative(var(:,j,k), i, dxyz(1,b))
               if (.not.is_null(2) .and. nj > 1) grad(2) = interior_derivative(var(i,:,k), j, dxyz(2,b))
               if (.not.is_null(3) .and. nk > 1) grad(3) = interior_derivative(var(i,j,:), k, dxyz(3,b))
               grad_max = max(grad_max, norm2(grad))
            enddo
         enddo
      enddo
      dc = self%block_spacing(b=b, delta_type=delta_type)
      delta = merge(delta_fine, delta_coarse, grad_max > tol)
      refinements_needed(b) = refinement_by_spacing(dc=dc, delta=delta)
   enddo
   endassociate
   contains
      pure function interior_derivative(v, n, ds) result(dv)
      !< Return the derivative of a line of interior cells at cell `n`: centred, one-sided at the ends.
      real(R8P),    intent(in) :: v(:) !< Line of interior values.
      integer(I4P), intent(in) :: n    !< Cell index.
      real(R8P),    intent(in) :: ds   !< Cell spacing.
      real(R8P)                :: dv   !< Derivative.

      if (n == 1) then
         dv = (v(2) - v(1)) / ds
      elseif (n == size(v)) then
         dv = (v(n) - v(n-1)) / ds
      else
         dv = 0.5_R8P * (v(n+1) - v(n-1)) / ds
      endif
      endfunction interior_derivative
   endsubroutine mark_by_gradient

   subroutine mark_by_solid(self, solid, delta_type, delta_fine, delta_coarse)
   !< Mark blocks by the surface of immersed solid `solid`: CHASE semantics, with the refine/derefine rule of the
   !< gradient marker.
   !<
   !< A block is crossed by the surface when its distance function (interior and ghost cells) changes sign; its
   !< admissible spacing is then `delta_fine`, `delta_coarse` otherwise. The distance function must be current (the
   !< caller computes it on the present grid). A run without solids, or a solid index out of range, is fatal (CHASE
   !< read `phi` unallocated in that case, issue #35 C-8).
   class(flume_common_object), intent(inout) :: self         !< The equation.
   integer(I4P),               intent(in)    :: solid        !< Solid index.
   character(*),               intent(in)    :: delta_type   !< Block spacing criterion: x, y, z, max.
   real(R8P),                  intent(in)    :: delta_fine   !< Admissible spacing across the surface.
   real(R8P),                  intent(in)    :: delta_coarse !< Admissible spacing elsewhere.
   real(R8P)                                 :: delta        !< Admissible spacing.
   integer(I4P)                              :: b            !< Counter.

   if (solid < 1_I4P .or. solid > self%ib%solids_number) &
      call mpih%error_stop(msg=': AMR solid marker: solid '//trim(str(solid))//' does not exist ([solids].(number) = '// &
                               trim(str(self%ib%solids_number))//')')
   self%adam%field%refinements_needed = [(TO_NOT_TOUCH, b=1, self%blocks_number)]
   do b=1, self%blocks_number
      delta = delta_coarse
      if (maxval(self%ib%phi(solid,:,:,:,b)) * minval(self%ib%phi(solid,:,:,:,b)) < 0._R8P) delta = delta_fine
      self%adam%field%refinements_needed(b) = refinement_by_spacing(dc=self%block_spacing(b=b, delta_type=delta_type), &
                                                                    delta=delta)
   enddo
   endsubroutine mark_by_solid

   ! public methods
   subroutine accumulate_seam_skin(self, flux_register, b, fec, weight, skin)
   !< Route one seam face skin of block `b`, face `fec`, weighted by `weight`, to the forest's flux register.
   !<
   !< Shared by the backends (the register is host-side): the CPU packs `skin` from its face fluxes, the FNL backend
   !< packs it on the device and copies it to the host. `skin(v, c)` is the face flux with the cell index `c` running
   !< over the two tangential axes, inner fastest (x faces: j, k; y faces: i, k; z faces: i, j), the register order.
   !< Coarse side (positive register index): the skin is the coarse face. Fine side (negative index): on an intra-realm
   !< AMR seam the skin is 2:1-restricted (2x2 average) into this block's quadrant of the coarse face, the quadrant offset
   !< precomputed by the forest from the Morton codes (`maps%amr_seam_quadrant`); on an inter-realm mirror seam (same
   !< resolution) it covers the coarse face 1:1 and is accumulated unrestricted.
   !<
   !< FLUME weighs every stage by its SSP coefficient, `weight = beta_s`: the register then holds the flux of the
   !< whole step, `sum_s beta_s F_s`, exactly the flux the committed update `q + dt sum_s beta_s dq_s` used, and the
   !< end-of-step correction restores conservation to round-off (issue #35, P5; AMReX non-subcycled flux register).
   class(flume_common_object),  intent(in)    :: self                !< The equation.
   class(flux_register_object), intent(inout) :: flux_register       !< Forest's flux register.
   integer(I4P),                intent(in)    :: b                   !< Block index.
   integer(I4P),                intent(in)    :: fec                 !< Face (1..6: -x, +x, -y, +y, -z, +z).
   real(R8P),                   intent(in)    :: weight              !< Stage weight.
   real(R8P),                   intent(in)    :: skin(1:,1:)         !< Face skin (nv, inner_n*outer_n).
   real(R8P), allocatable                     :: slab(:,:)           !< Register-shaped contribution (nv_reg, nface_cells).
   real(R8P), allocatable                     :: fine_face(:,:,:)    !< Weighted fine face (nv, inner_n, outer_n).
   integer(I4P)                               :: sgn_idx, face_idx   !< Signed and absolute register face index.
   integer(I4P)                               :: inner_n, outer_n    !< Tangential cell counts.
   integer(I4P)                               :: ioff, joff          !< Fine-block quadrant offset.
   integer(I4P)                               :: nv                  !< Skin variables number.

   sgn_idx = self%adam%maps%inter_realm_face_register_index(b, fec)
   if (sgn_idx == 0_I4P) return
   face_idx = abs(sgn_idx)
   if (face_idx > flux_register%nfaces) return
   if (.not.allocated(flux_register%face(face_idx)%F_coarse)) return
   select case(fec)
   case(1_I4P, 2_I4P)
      inner_n = self%nj ; outer_n = self%nk
   case(3_I4P, 4_I4P)
      inner_n = self%ni ; outer_n = self%nk
   case default
      inner_n = self%ni ; outer_n = self%nj
   endselect
   nv = size(skin, dim=1)
   if (size(skin, dim=2) /= inner_n * outer_n .or. flux_register%face(face_idx)%nface_cells /= inner_n * outer_n) &
      call mpih%error_stop(msg=': accumulate_seam_skin: skin size differs from the register face (block '// &
                               trim(str(b))//', face '//trim(str(fec))//')')
   allocate(slab(size(flux_register%face(face_idx)%F_coarse, dim=1), inner_n * outer_n))
   slab = 0._R8P
   if (sgn_idx > 0_I4P) then
      slab(1:nv,:) = weight * skin
      call flux_register%accumulate_coarse_flux(face_index=face_idx, stage=1_I4P, flux_face=slab)
   elseif (flux_register%face(face_idx)%seam_kind == SEAM_KIND_INTER_REALM) then
      ! same-resolution (mirror) inter-realm seam: the fine skin covers the coarse face 1:1, no restriction (issue #37)
      slab(1:nv,:) = weight * skin
      call flux_register%accumulate_fine_flux(face_index=face_idx, stage=1_I4P, flux_face=slab)
   else
      ioff = 0_I4P ; joff = 0_I4P
      if (allocated(self%adam%maps%amr_seam_quadrant)) then
         ioff = self%adam%maps%amr_seam_quadrant(1, b, fec)
         joff = self%adam%maps%amr_seam_quadrant(2, b, fec)
      endif
      allocate(fine_face(nv, inner_n, outer_n))
      fine_face = weight * reshape(skin, [nv, inner_n, outer_n])
      call restrict_fine_face_to_quadrant(fine_face=fine_face, inner_n=inner_n, outer_n=outer_n, ioff=ioff, joff=joff, &
                                          slab=slab)
      call flux_register%accumulate_fine_flux(face_index=face_idx, stage=1_I4P, flux_face=slab)
   endif
   endsubroutine accumulate_seam_skin

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

   subroutine compute_phi(self)
   !< Compute the distance function of the immersed solids on the host (a no-op without solids); the solids are
   !< static, so the backends compute it once per grid (initialization, restart, initial AMR).
   class(flume_common_object), intent(inout) :: self !< The equation.

   if (self%ib%solids_number > 0_I4P) call self%ib%compute_phi(field=self%adam%field, grid=self%adam%grid, verbose=.true.)
   endsubroutine compute_phi

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
   integer(I4P)                                      :: error          !< Error status.

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
   call file_parameters%get(section_name='IO', option_name='save_auxiliary_fields', val=self%save_auxiliary_fields, &
                            error=error)
   if (error > 0) call mpih%error_stop(msg=': failed to load [IO].(save_auxiliary_fields)')
   call self%check_slices
   call self%check_ngc_number
   call self%check_amr_block_cells
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
   !< Save fields in XH5F format: `q` always, `dq` when `[IO].(save_residual_fields)`, the auxiliary variables (computed
   !< from the saved `q`, ghost cells included) when `[IO].(save_auxiliary_fields)`.
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
   if (self%save_auxiliary_fields) call self%compute_q_aux_host
   call self%open_file_xh5f(basename=trim(output_basename_), xh5f=xh5f)
   do b=1, self%adam%field%blocks_number
      bn = 'block_'//trim(strz(b, 9))//'-proc'//trim(strz(mpih%myrank, 6))
      call self%open_block_xh5f(xh5f=xh5f, b=b, nijk=nijk, t=self%time%it, time=self%time%time)
      call self%io%save_field(xh5f=xh5f, grid=self%adam%grid, block_name=bn, ijk=ijk, nijk=nijk, &
                              q=self%q(:,:,:,:,b), q_name=self%q_name)
      if (self%io%save_residual_fields) &
         call self%io%save_field(xh5f=xh5f, grid=self%adam%grid, block_name=bn, ijk=ijk, nijk=nijk, &
                                 q=self%dq(:,:,:,:,b), q_name=self%dq_name)
      if (self%save_auxiliary_fields) &
         call self%io%save_field(xh5f=xh5f, grid=self%adam%grid, block_name=bn, ijk=ijk, nijk=nijk, &
                                 q=self%q_aux(:,:,:,:,b), q_name=self%q_aux_name)
      call self%close_block_xh5f(xh5f=xh5f)
   enddo
   call self%close_file_xh5f(xh5f=xh5f)
   call mpih%barrier(tictoc=.true.)
   endsubroutine save_xh5f

   subroutine save_slices(self)
   !< Save the slices (library `slices_object`, `[slices]` / `[slice_N]`) of the conservative variables on their
   !< cadence; the caller has refreshed the ghost cells of the host `q` (the interpolation stencils read them).
   class(flume_common_object), intent(inout) :: self    !< The equation.
   character(len=8), allocatable             :: name(:) !< Variables names.
   integer(I4P)                              :: v       !< Counter.

   if (.not.self%slices%is_to_save(it=self%time%it, it_max=self%time%it_max, time=self%time%time, &
                                   time_max=self%time%time_max)) return
   allocate(name(size(self%q_name)))
   do v=1, size(self%q_name)
      name(v) = self%q_name(v)%chars()
   enddo
   call self%slices%save_mat(basename=self%io%output_basename, it=self%time%it, it_max=self%time%it_max, &
                             time=self%time%time, time_max=self%time%time_max, adam=self%adam, q=self%q, q_name=name)
   endsubroutine save_slices

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

   ! public procedures
   pure function ib_cut_spacing(phi_c, phi_m, phi_p, ds, eps) result(ds_cut)
   !< Return the spacing of a fluid cell along one direction, shortened where the solid surface crosses its stencil
   !< (CHASE semantics, issue #35 D-9).
   !<
   !< When the fluid cell (`phi_c < 0`) has one neighbour inside the solid (`phi_m phi_p < 0`), the surface lies at the
   !< distance `delta = -phi_c / (phi_s - phi_c + eps) ds` from the cell centre towards the solid neighbour `s`, and
   !< the spacing becomes `ds / 2 + delta`; otherwise it is `ds`.
   real(R8P), intent(in) :: phi_c  !< Distance function of the cell (negative in the fluid).
   real(R8P), intent(in) :: phi_m  !< Distance function of the minus neighbour.
   real(R8P), intent(in) :: phi_p  !< Distance function of the plus neighbour.
   real(R8P), intent(in) :: ds     !< Spacing.
   real(R8P), intent(in) :: eps    !< Guard against a vanishing denominator.
   real(R8P)             :: ds_cut !< Spacing, cut by the surface.
   !$acc routine seq
   !$omp declare target

   ds_cut = ds
   if (phi_c < 0._R8P .and. phi_m * phi_p < 0._R8P) then
      if (phi_p > 0._R8P) then
         ds_cut = 0.5_R8P * ds - phi_c / (phi_p - phi_c + eps) * ds
      else
         ds_cut = 0.5_R8P * ds - phi_c / (phi_m - phi_c + eps) * ds
      endif
   endif
   endfunction ib_cut_spacing

   pure function refinement_by_spacing(dc, delta) result(refinement)
   !< Return the refinement query of a block of spacing `dc` against the admissible spacing `delta`: refine when too
   !< coarse, derefine when its parent (spacing doubled) would still be admissible, untouched otherwise.
   real(R8P), intent(in) :: dc         !< Block spacing.
   real(R8P), intent(in) :: delta      !< Admissible spacing.
   integer(I4P)          :: refinement !< Refinement query.

   if (dc > delta) then
      refinement = TO_BE_REFINED
   elseif (2._R8P * dc <= delta) then
      refinement = TO_BE_DEREFINED
   else
      refinement = TO_NOT_TOUCH
   endif
   endfunction refinement_by_spacing

   pure subroutine seam_skin_cell(axis, sgn, ni, nj, nk, c, i, j, k)
   !< Return the interior cell `(i, j, k)` of skin cell `c` on the face of normal `axis` and side `sgn` (register order,
   !< inner tangential axis fastest). Shared by the host and device reflux applications.
   integer(I4P), intent(in)  :: axis       !< Face normal axis, 1..3.
   integer(I4P), intent(in)  :: sgn        !< Face side, +1 maximum, -1 minimum.
   integer(I4P), intent(in)  :: ni, nj, nk !< Grid dimensions.
   integer(I4P), intent(in)  :: c          !< Skin cell index.
   integer(I4P), intent(out) :: i, j, k    !< Interior cell (0 on a malformed axis).
   !$acc routine seq
   !$omp declare target

   select case(axis)
   case(1_I4P)
      i = merge(ni, 1_I4P, sgn > 0_I4P) ; j = 1_I4P + mod(c - 1_I4P, nj) ; k = 1_I4P + (c - 1_I4P) / nj
   case(2_I4P)
      i = 1_I4P + mod(c - 1_I4P, ni) ; j = merge(nj, 1_I4P, sgn > 0_I4P) ; k = 1_I4P + (c - 1_I4P) / ni
   case(3_I4P)
      i = 1_I4P + mod(c - 1_I4P, ni) ; j = 1_I4P + (c - 1_I4P) / ni ; k = merge(nk, 1_I4P, sgn > 0_I4P)
   case default
      i = 0_I4P ; j = 0_I4P ; k = 0_I4P
   endselect
   endsubroutine seam_skin_cell

   ! private methods
   subroutine check_slices(self)
   !< Check the interpolation type of every slice: the library interpolation leaves the value undefined for an unknown
   !< type, so an unknown one is fatal here.
   class(flume_common_object), intent(in) :: self !< The equation.
   integer(I4P)                           :: s    !< Counter.

   do s=1, self%slices%slices_number
      select case(trim(self%slices%slice(s)%itype))
      case('trilinear', 'inverse_distance')
      case default
         call mpih%error_stop(msg=': unknown [slice_'//trim(str(s, .true.))//'].(itype) "'// &
                                  trim(self%slices%slice(s)%itype)//'"; expected one of trilinear, inverse_distance')
      endselect
   enddo
   endsubroutine check_slices

   subroutine compute_q_aux_host(self)
   !< Compute the auxiliary variables of the host `q` on every cell, ghost cells included (output only: the backends
   !< compute their own auxiliary variables in the space operator).
   class(flume_common_object), intent(inout) :: self       !< The equation.
   integer(I4P)                              :: b, i, j, k !< Counters.

   do b=1, self%blocks_number
      do k=1-self%ngc, self%nk+self%ngc
         do j=1-self%ngc, self%nj+self%ngc
            do i=1-self%ngc, self%ni+self%ngc
               call conservative_to_auxiliary(gamma=self%physics%gamma, R=self%physics%R, q=self%q(:,i,j,k,b), &
                                              qa=self%q_aux(:,i,j,k,b))
            enddo
         enddo
      enddo
   enddo
   endsubroutine compute_q_aux_host

   function block_spacing(self, b, delta_type) result(dc)
   !< Return the spacing of block `b` by the AMR delta criterion: `x`, `y`, `z`, or `max` over the active directions.
   class(flume_common_object), intent(in) :: self       !< The equation.
   integer(I4P),               intent(in) :: b          !< Block index.
   character(*),               intent(in) :: delta_type !< Delta criterion.
   real(R8P)                              :: dc         !< Block spacing.

   dc = 0._R8P
   select case(delta_type)
   case(AMR_DELTA_T_X)
      dc = self%adam%field%dxyz(1,b)
   case(AMR_DELTA_T_Y)
      dc = self%adam%field%dxyz(2,b)
   case(AMR_DELTA_T_Z)
      dc = self%adam%field%dxyz(3,b)
   case(AMR_DELTA_T_MAX)
      dc = maxval(self%adam%field%dxyz(:,b), mask=.not.self%adam%grid%null_xyz)
   case default
      call mpih%error_stop(msg=': unknown AMR marker delta_type "'//delta_type//'"; expected one of '// &
                               AMR_DELTA_T_X//', '//AMR_DELTA_T_Y//', '//AMR_DELTA_T_Z//', '//AMR_DELTA_T_MAX)
   endselect
   endfunction block_spacing

   subroutine check_amr_block_cells(self)
   !< Check that a run with init-time refinement has an even number of block cells along every non-null axis.
   !<
   !< A 2:1 child covers half of its parent: with an odd cell count its boundary falls in the middle of a parent cell,
   !< the coarse-fine ghost fill reads undefined values and the residual is NaN from the first stage (issue #37,
   !< measured with ni = 25). Uniform refinement (`iu_ref_levels`) has no coarse-fine faces and is not affected.
   class(flume_common_object), intent(in) :: self                     !< The equation.
   character(len=1), parameter            :: AXIS(3)=['i', 'j', 'k'] !< Axes names.
   integer(I4P)                           :: n(3)                     !< Block cells numbers.
   integer(I4P)                           :: d                        !< Axis counter.

   if (self%ic%amr_iterations <= 0_I4P) return
   n = [self%ni, self%nj, self%nk]
   do d=1, 3
      if (self%adam%grid%null_xyz(d)) cycle
      if (mod(n(d), 2_I4P) /= 0_I4P) &
         call mpih%error_stop(msg=': [grid].(n'//AXIS(d)//')='//trim(str(n(d)))//' is odd: the init-time 2:1 '// &
                                  'refinement ([initial_conditions].(amr_iterations) > 0) needs an even number of '// &
                                  'block cells along every non-null axis')
   enddo
   endsubroutine check_amr_block_cells

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
