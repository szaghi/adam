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
use :: adam_fdv_operators_library,    only : compute_derivative1_fd_centered
use :: adam_flux_register_object,     only : flux_register_object, restrict_fine_face_to_quadrant, SEAM_KIND_INTER_REALM
use :: adam_parameters,               only : TO_BE_DEREFINED, TO_BE_REFINED, TO_NOT_TOUCH
use :: adam_realm_object,             only : realm_object
use :: adam_rk_object,                only : rk_stored_stages_number, RK_SSP_11, RK_SSP_22, RK_SSP_33, RK_SSP_54
! ADAM singleton objects
use :: adam_mpih_global,              only : mpih
! FLUME modules
use :: adam_flume_bc_object,          only : flume_bc_object
use :: adam_flume_diagnostics_object, only : flume_diagnostics_object
use :: adam_flume_euler_library,      only : conservative_to_auxiliary
use :: adam_flume_ic_object,          only : flume_ic_object
use :: adam_flume_numerics_object,    only : flume_numerics_object
use :: adam_flume_mhd_library,        only : mhd_conservative_to_auxiliary, mhd_eglm_conservative_to_auxiliary
use :: adam_flume_parameters,         only : GLM_CH_CHECK_ERROR, IA_BX, IA_BY, IA_BZ, IA_P, IQ_BX, IQ_BY, IQ_BZ, IQ_RU, &
                                             MODEL_EULER, MODEL_MHD, MODEL_MHD_EGLM, MODEL_MHD_GLM,                   &
                                             POSITIVITY_LIMITER_CELL, RECON_PRIMITIVE, RIEMANN_SOLVER_HLL,            &
                                             RIEMANN_SOLVER_HLLC, RIEMANN_SOLVER_HLLD, RIEMANN_SOLVER_LLF,            &
                                             SCHEME_SPACE_WENO_RIEMANN
use :: adam_weno_object,              only : WENO_WEIGHTS_SI
use :: adam_flume_physics_object,     only : flume_physics_object
use :: adam_flume_reference_object,   only : flume_reference_object
use :: adam_flume_time_object,        only : flume_time_object
! third party modules
use :: finer,                         only : file_ini
use :: motion,                        only : xh5f_file_object
use :: mpi
use :: penf,                          only : I4P, I8P, R8P, str, strz
use :: stringifor,                    only : string

implicit none
private
public :: flume_common_object
public :: ib_cut_spacing
public :: seam_skin_cell

character(len=11), parameter :: SCHEME_TIME_TAG="runge_kutta" !< Time-integration family tag (forest admissibility).
character(len=4),  parameter :: MHD_DERIVED_NAME(4)=['pt  ', 'beta', 'bmag', 'divb'] !< MHD derived fields names.

type, extends(realm_object) :: flume_common_object
   !< FLUME common object: data and methods shared by all backends.
   ! AMR
   logical                        :: amr_locked_=.false. !< Runtime AMR locked after initialization.
   ! IO
   logical                        :: save_auxiliary_fields=.false. !< Save the auxiliary variables with the fields.
   ! GLM and div(B) monitors
   real(R8P)                      :: glm_speed_reported=0._R8P !< Largest wave speed above c_h reported so far.
   real(R8P)                      :: divb_reported=0._R8P      !< Largest max|div B| above divb_tol reported so far.
   integer(I4P),      allocatable :: divb_seam(:,:)            !< Seam faces flags of the div(B) history [nb, 6].
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
   type(flume_reference_object)   :: units               !< Reference layer (dimensional input, issue #49).
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
      procedure, pass(self) :: glm_lambda            !< Return the GLM bound of the local dt, c_h max sum_d 1/dx_d.
      procedure, pass(self) :: initialize            !< Initialize the common data.
      procedure, pass(self) :: null_freeze           !< Return the variable each null direction freezes.
      procedure, pass(self) :: output_factors        !< Return the output factors of written variables.
      procedure, pass(self) :: output_names          !< Return the names of every written variable.
      procedure, pass(self) :: report_divb           !< Reduce and save the div(B) norms, apply the divb_tol monitor.
      procedure, pass(self) :: report_glm_speed      !< Check c_h against the fastest wave (warning or stop).
      procedure, pass(self) :: nonfinite_total       !< Sum the non-finite values count over the ranks.
      procedure, pass(self) :: stop_nonfinite        !< Stop on a non-finite state, locating it on every rank.
      procedure, pass(self) :: load_restart_files    !< Load restart files.
      procedure, pass(self) :: save_restart_files    !< Save restart files.
      procedure, pass(self) :: save_slices           !< Save the slices on their cadence.
      procedure, pass(self) :: save_xh5f             !< Save fields in XH5F format.
      procedure, pass(self) :: set_divb_seam         !< Set the seam faces flags of the div(B) history.
      procedure, pass(self) :: set_glm_damping       !< Set the GLM damping once the grid exists.
      ! forest methods
      procedure, pass(self) :: coupling_descriptor_forest !< Return the realm coupling descriptor.
      ! private methods
      procedure, pass(self), private :: block_spacing    !< Return the spacing of a block by a delta criterion.
      procedure, pass(self), private :: check_ngc_number      !< Check the ghost cells number against the stencils.
      procedure, pass(self), private :: check_positivity_limiter !< Refuse the limiter where it cannot work.
      procedure, pass(self), private :: check_slices     !< Check the slices interpolation types.
      procedure, pass(self), private :: check_weno_scheme     !< Refuse the centred WENO schemes.
      procedure, pass(self), private :: initialize_riemann_scheme !< Check and set up the weno-riemann scheme.
      procedure, pass(self), private :: compute_mhd_derived !< Compute the MHD derived output fields of one block.
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
   real(R8P)                                 :: grad(3)      !< Gradient of one cell.
   real(R8P)                                 :: grad_max     !< Maximum gradient magnitude of one block.
   real(R8P)                                 :: dc           !< Block spacing.
   real(R8P)                                 :: delta        !< Admissible spacing.
   integer(I4P)                              :: b, i, j, k   !< Counters.

   if ((field == 1_I4P .and. (ivar < 1_I4P .or. ivar > self%physics%nv))     .or. &
       (field == 2_I4P .and. (ivar < 1_I4P .or. ivar > self%physics%nv_aux)) .or. (field < 1_I4P .or. field > 2_I4P)) &
      call mpih%error_stop(msg=': AMR gradient marker: invalid field '//trim(str(field))//' / ivar '//trim(str(ivar)))
   self%adam%field%refinements_needed = [(TO_NOT_TOUCH, b=1, self%blocks_number)]
   if (field == 2_I4P) call self%compute_q_aux_host
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
                  var(i,j,k) = self%q_aux(ivar,i,j,k,b)
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
   !< resolution) it is accumulated unrestricted into the register faces it overlaps (`maps%seam_overlap`, issue #51:
   !< the blocks of the two realms need not line up; one full-skin overlap when they do).
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
   if (sgn_idx < 0_I4P .and. flux_register%face(face_idx)%seam_kind == SEAM_KIND_INTER_REALM) then
      ! same-resolution (mirror) inter-realm seam: the skin goes unrestricted into every register face it overlaps
      ! (issues #37, #51)
      associate(mp => self%adam%maps)
         call flux_register%accumulate_fine_overlaps(                                                                &
            overlaps=mp%seam_overlap(:, mp%seam_overlap_start(b, fec):mp%seam_overlap_start(b, fec) +                 &
                                        mp%seam_overlap_count(b, fec) - 1_I4P),                                     &
            inner_n=inner_n, skin=skin, weight=weight)
      endassociate
      return
   endif
   if (size(skin, dim=2) /= inner_n * outer_n .or. flux_register%face(face_idx)%nface_cells /= inner_n * outer_n) &
      call mpih%error_stop(msg=': accumulate_seam_skin: skin size differs from the register face (block '// &
                               trim(str(b))//', face '//trim(str(fec))//')')
   allocate(slab(size(flux_register%face(face_idx)%F_coarse, dim=1), inner_n * outer_n))
   slab = 0._R8P
   if (sgn_idx > 0_I4P) then
      slab(1:nv,:) = weight * skin
      call flux_register%accumulate_coarse_flux(face_index=face_idx, stage=1_I4P, flux_face=slab)
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
   !< (`rk_stored_stages_number` nv), plus the positivity limiter's cell factors (nv, q-shaped for the ghost exchange).
   !< Requires `physics` and `numerics` initialized.
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
   if (self%numerics%positivity_limiter == POSITIVITY_LIMITER_CELL) fields_number = fields_number + self%physics%nv
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

   function glm_lambda(self) result(lambda)
   !< Return the GLM bound of the local time step, `c_h max_b sum_{d active} 1 / dx_d` (issue #41, section 3.5): the
   !< `(B_n, psi)` waves travel at `c_h` along every active direction, the multi-dimensional form of `dt <= CFL dx / c_h`
   !< consistent with the fluid bound `sum_d (|u_d| + c_{f,d}) / dx_d`. Zero without GLM; local (the forest reduces dt).
   class(flume_common_object), intent(in) :: self   !< The equation.
   real(R8P)                              :: lambda !< GLM bound, c_h max sum_d 1 / dx_d.
   real(R8P)                              :: w(3)   !< Direction weights: 1 active, 0 null.
   integer(I4P)                           :: b      !< Counter.

   lambda = 0._R8P
   if (.not.self%physics%mhd%has_glm) return
   w = merge(0._R8P, 1._R8P, self%adam%grid%null_xyz)
   do b=1, self%blocks_number
      lambda = max(lambda, sum(w / self%adam%field%dxyz(:,b)))
   enddo
   lambda = self%physics%mhd%glm_ch * lambda
   endfunction glm_lambda

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
   ! dimensional input: convert it to code units in the loaded file, before any parser reads it (the IO options, read
   ! by io%initialize, are read again)
   call self%units%initialize(file_parameters=self%io%file_parameters)
   if (self%units%is_active) call self%io%load_from_file(file_parameters=self%io%file_parameters)
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
   if (self%physics%model /= MODEL_EULER .and. self%ib%solids_number > 0_I4P) &
      call mpih%error_stop(msg=': immersed solids are not supported with [physics].(physical_model) = '// &
                               self%physics%physical_model)
   call self%bc%initialize(file_parameters=file_parameters, physics=self%physics)
   call self%adam%grid%set_bc_type(bc_type=self%bc%bc_type)
   call self%time%initialize(file_parameters=file_parameters)
   call self%ic%initialize(file_parameters=file_parameters, physics=self%physics)
   call self%diagnostics%initialize(file_parameters=file_parameters)
   call file_parameters%get(section_name='IO', option_name='save_auxiliary_fields', val=self%save_auxiliary_fields, &
                            error=error)
   if (error > 0) call mpih%error_stop(msg=': failed to load [IO].(save_auxiliary_fields)')
   call self%check_slices
   call self%check_weno_scheme
   call self%initialize_riemann_scheme
   call self%check_positivity_limiter
   call self%check_ngc_number
   call self%allocate_common
   call self%io_initialize
   if (self%units%is_active) call self%units%save_units(basename=trim(self%io%output_basename), names=self%output_names())
   if (self%adam%tree%iu_ref_levels > 0) &
      call self%adam%refine_uniform(refinement_levels=self%adam%tree%iu_ref_levels, do_mpi_redistribute=.true., &
                                    do_blocks_reorder=.false., q=self%q)
   endassociate
   if (verbose_) call mpih%print_message('flume_common_object%initialize finish')
   endsubroutine initialize

   function output_factors(self, names) result(factor)
   !< Return the output factors of written variables (1 with code output units, issue #49 N2c).
   class(flume_common_object), intent(in) :: self      !< The equation.
   type(string),               intent(in) :: names(1:) !< Variables names.
   real(R8P), allocatable                 :: factor(:) !< Output factors.
   integer(I4P)                           :: v         !< Counter.

   allocate(factor(size(names, dim=1)))
   do v=1, size(names, dim=1)
      factor(v) = self%units%variable_output(names(v)%chars())
   enddo
   endfunction output_factors

   function output_names(self) result(names)
   !< Return the names of every variable FLUME may write: conservative, residuals, auxiliary and, for MHD, derived.
   class(flume_common_object), intent(in) :: self     !< The equation.
   type(string), allocatable              :: names(:) !< Variables names.
   integer(I4P)                           :: n, v     !< Names number, counter.

   n = size(self%q_name) + size(self%dq_name) + size(self%q_aux_name)
   if (self%physics%model /= MODEL_EULER) n = n + size(MHD_DERIVED_NAME)
   allocate(names(n))
   n = 0
   do v=1, size(self%q_name)
      n = n + 1 ; names(n) = self%q_name(v)
   enddo
   do v=1, size(self%dq_name)
      n = n + 1 ; names(n) = self%dq_name(v)
   enddo
   do v=1, size(self%q_aux_name)
      n = n + 1 ; names(n) = self%q_aux_name(v)
   enddo
   if (self%physics%model /= MODEL_EULER) then
      do v=1, size(MHD_DERIVED_NAME)
         n = n + 1 ; names(n) = trim(MHD_DERIVED_NAME(v))
      enddo
   endif
   endfunction output_names

   function null_freeze(self) result(freeze)
   !< Return, per direction, the conservative variable whose residual a null direction freezes (0: none).
   !<
   !< Euler: the momentum along the null direction (CHASE semantics, issue #35, section 3.4). MHD: none, the 1-D MHD
   !< Riemann problems evolve the transverse momentum and field through the fluxes of the active directions (issue #41,
   !< M2-P3). Host data passed to the model-independent flux difference: no model branch in the kernels.
   class(flume_common_object), intent(in) :: self      !< The equation.
   integer(I4P)                           :: freeze(3) !< Frozen variable of each direction (0: none).
   integer(I4P)                           :: d         !< Direction counter.

   freeze = 0_I4P
   select case(self%physics%model)
   case(MODEL_EULER)
      do d=1, 3
         if (self%adam%grid%null_xyz(d)) freeze(d) = IQ_RU + d - 1
      enddo
   case(MODEL_MHD, MODEL_MHD_GLM, MODEL_MHD_EGLM)
   case default
      call mpih%error_stop(msg=': no null-direction rule for physical model "'//self%physics%physical_model//'"')
   endselect
   endfunction null_freeze

   subroutine report_divb(self, norms)
   !< Reduce the div(B) norms of this rank over all ranks, save the div(B) history row and apply the monitor (issue #41,
   !< D-10): `max|div B| > [mhd].(divb_tol)` (> 0) is fatal with `divb_error`, otherwise a warning logged by rank 0 each
   !< time the maximum exceeds the largest one reported so far.
   class(flume_common_object), intent(inout) :: self     !< The equation.
   real(R8P),                  intent(in)    :: norms(3) !< max|div B|, sum |div B| dV, seam-local max|div B| (rank).
   real(R8P)                                 :: g(3)     !< Global norms.

   g = norms
   call MPI_ALLREDUCE(MPI_IN_PLACE, g(1), 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, mpih%error)
   call MPI_ALLREDUCE(MPI_IN_PLACE, g(2), 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih%error)
   call MPI_ALLREDUCE(MPI_IN_PLACE, g(3), 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, mpih%error)
   call self%diagnostics%save_divb_row(it=self%time%it, time=self%time%time*self%units%time_output(),          &
                                       norms=[g(1)*self%units%variable_output('divb'),                         &
                                              g(2)*self%units%variable_output('bmag')*self%units%length_output()**2, &
                                              g(3)*self%units%variable_output('divb')])
   if (.not.(self%physics%mhd%divb_tol > 0._R8P .and. g(1) > self%physics%mhd%divb_tol)) return
   if (self%physics%mhd%divb_error) &
      call mpih%error_stop(msg=': max|div B| = '//trim(str(g(1)))//' > [mhd].(divb_tol) = '//             &
                               trim(str(self%physics%mhd%divb_tol))//' at step '//trim(str(self%time%it))// &
                               ' ([mhd].(divb_error) = .true.)')
   if (g(1) > self%divb_reported) then
      self%divb_reported = g(1)
      if (mpih%myrank == 0) print '(A)', mpih%myrankstr//'warning: max|div B| = '//trim(str(g(1)))//                &
                                         ' > [mhd].(divb_tol) = '//trim(str(self%physics%mhd%divb_tol))//' at step '// &
                                         trim(str(self%time%it))
   endif
   endsubroutine report_divb

   subroutine report_glm_speed(self, speed_max)
   !< Check the GLM cleaning speed against the fastest wave (issue #41, section 3.5, D-9): `max(|u_d| + c_{f,d}) > c_h`
   !< is fatal with `[mhd].(glm_ch_check) = error`, otherwise a warning logged by rank 0 each time the speed exceeds the
   !< largest one reported so far (GLM stays stable, dt includes c_h, but the cleaning is slower than the fastest wave).
   class(flume_common_object), intent(inout) :: self      !< The equation.
   real(R8P),                  intent(in)    :: speed_max !< Fastest wave speed of this rank.
   real(R8P)                                 :: speed     !< Fastest wave speed of all ranks.

   speed = speed_max
   call MPI_ALLREDUCE(MPI_IN_PLACE, speed, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, mpih%error)
   if (.not.(speed > self%physics%mhd%glm_ch)) return
   if (self%physics%mhd%glm_ch_check == GLM_CH_CHECK_ERROR) &
      call mpih%error_stop(msg=': max(|u| + c_f) = '//trim(str(speed))//' > [mhd].(glm_ch) = '//            &
                               trim(str(self%physics%mhd%glm_ch))//' at step '//trim(str(self%time%it))// &
                               ' ([mhd].(glm_ch_check) = error)')
   if (speed > self%glm_speed_reported) then
      self%glm_speed_reported = speed
      if (mpih%myrank == 0) print '(A)', mpih%myrankstr//'warning: max(|u| + c_f) = '//trim(str(speed))//            &
                                         ' > [mhd].(glm_ch) = '//trim(str(self%physics%mhd%glm_ch))//' at step '// &
                                         trim(str(self%time%it))//': the cleaning is slower than the fastest wave'
   endif
   endsubroutine report_glm_speed

   function nonfinite_total(self, n_local) result(n_total)
   !< Return the number of non-finite values of the committed state over every rank (issue #45).
   class(flume_common_object), intent(in) :: self    !< The equation.
   integer(I8P),               intent(in) :: n_local !< Non-finite values number of this rank.
   integer(I8P)                           :: n_total !< Non-finite values number of all ranks.

   n_total = n_local
   call MPI_ALLREDUCE(MPI_IN_PLACE, n_total, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, mpih%error)
   endfunction nonfinite_total

   subroutine stop_nonfinite(self, n_total)
   !< Stop the run on a non-finite (NaN or infinite) committed state (issue #45), reporting on every rank the first
   !< non-finite value of the host `q` (the caller has made it current): block Morton code, cell centre, variable. Every
   !< rank calls it after the reduction, so the stop is collective (the first `error_stop` aborts the others, #43).
   class(flume_common_object), intent(in) :: self                   !< The equation.
   integer(I8P),               intent(in) :: n_total                !< Non-finite values number of all ranks.
   integer(I8P), parameter                :: EXPONENT_MASK=2047_I8P !< The 11 exponent bits of a binary64.
   character(:), allocatable              :: site                   !< First non-finite value of this rank.
   integer(I4P)                           :: b, i, j, k, v          !< Counters.

   site = 'none on this rank'
   search: do b=1, self%blocks_number
      do k=1, self%nk
         do j=1, self%nj
            do i=1, self%ni
               do v=1, self%physics%nv
                  if (iand(ishft(transfer(self%q(v,i,j,k,b), 0_I8P), -52), EXPONENT_MASK) == EXPONENT_MASK) then
                     site = 'first at block code '//trim(str(self%adam%field%code(b)))//', cell ('//            &
                            trim(str(i))//', '//trim(str(j))//', '//trim(str(k))//') centre ('//                    &
                            trim(str(self%adam%field%x_cell(i,b)))//', '//trim(str(self%adam%field%y_cell(j,b)))// &
                            ', '//trim(str(self%adam%field%z_cell(k,b)))//'), variable '//self%q_name(v)%chars()
                     exit search
                  endif
               enddo
            enddo
         enddo
      enddo
   enddo search
   call mpih%error_stop(msg=': '//trim(str(n_total))//' non-finite (NaN or infinite) values in the state at step '// &
                            trim(str(self%time%it))//'; '//site)
   endsubroutine stop_nonfinite

   subroutine load_restart_files(self, t, time)
   !< Load restart files.
   class(flume_common_object), intent(inout) :: self !< The equation.
   integer(I4P),               intent(out)   :: t    !< Time iteration.
   real(R8P),                  intent(out)   :: time !< Time.

   call self%units%check_restart(basename=self%io%restart_basename)
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
   call self%units%save_restart(basename=self%io%restart_basename)
   call self%save_xh5f(output_basename=self%io%restart_basename)
   call mpih%barrier(tictoc=.true.)
   endsubroutine save_restart_files

   subroutine save_xh5f(self, output_basename, with_ghost)
   !< Save fields in XH5F format: `q` always, `dq` when `[IO].(save_residual_fields)`, the auxiliary variables (computed
   !< from the saved `q`, ghost cells included) when `[IO].(save_auxiliary_fields)`. With `[reference] output_units =
   !< dimensional` every field, the grid and the time are multiplied by their output factors (a scaled copy of each block,
   !< the state itself is untouched).
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
   real(R8P),    allocatable                        :: derived(:,:,:,:) !< MHD derived fields of one block.
   type(string)                                     :: derived_name(4)  !< MHD derived fields names.
   logical                                          :: with_derived     !< Save the MHD derived fields.
   real(R8P),    allocatable                        :: factor_q(:)      !< Output factors of q.
   real(R8P),    allocatable                        :: factor_dq(:)     !< Output factors of dq.
   real(R8P),    allocatable                        :: factor_aux(:)    !< Output factors of q_aux.
   real(R8P),    allocatable                        :: factor_derived(:) !< Output factors of the MHD derived fields.
   integer(I4P)                                     :: b, v             !< Counters.

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
   with_derived = self%save_auxiliary_fields .and. self%physics%model /= MODEL_EULER
   if (with_derived) then
      allocate(derived(4,1-self%ngc:self%ni+self%ngc,1-self%ngc:self%nj+self%ngc,1-self%ngc:self%nk+self%ngc))
      do v=1, 4
         derived_name(v) = trim(MHD_DERIVED_NAME(v))
      enddo
      factor_derived = self%output_factors(derived_name)
   endif
   factor_q   = self%output_factors(self%q_name)
   factor_dq  = self%output_factors(self%dq_name)
   factor_aux = self%output_factors(self%q_aux_name)
   call self%open_file_xh5f(basename=trim(output_basename_), xh5f=xh5f)
   do b=1, self%adam%field%blocks_number
      bn = 'block_'//trim(strz(b, 9))//'-proc'//trim(strz(mpih%myrank, 6))
      call self%open_block_xh5f(xh5f=xh5f, b=b, nijk=nijk, t=self%time%it, time=self%time%time*self%units%time_output(), &
                                length_scale=self%units%length_output())
      call save_block(q=self%q(:,:,:,:,b), q_name=self%q_name, factor=factor_q)
      if (self%io%save_residual_fields) call save_block(q=self%dq(:,:,:,:,b), q_name=self%dq_name, factor=factor_dq)
      if (self%save_auxiliary_fields) call save_block(q=self%q_aux(:,:,:,:,b), q_name=self%q_aux_name, factor=factor_aux)
      if (with_derived) then
         call self%compute_mhd_derived(b=b, derived=derived)
         call save_block(q=derived, q_name=derived_name, factor=factor_derived)
      endif
      call self%close_block_xh5f(xh5f=xh5f)
   enddo
   call self%close_file_xh5f(xh5f=xh5f)
   call mpih%barrier(tictoc=.true.)
   contains
      subroutine save_block(q, q_name, factor)
      !< Save the variables of one block, multiplied by their output factors (as they are with code output units).
      real(R8P),    intent(in) :: q(1:,1:,1:,1:) !< Variables of the block [nv,ni,nj,nk], ghost cells included.
      type(string), intent(in) :: q_name(1:)     !< Variables names [nv].
      real(R8P),    intent(in) :: factor(1:)     !< Output factors [nv].
      real(R8P), allocatable   :: q_out(:,:,:,:) !< Variables in output units.
      integer(I4P)             :: v_             !< Counter.

      if (.not.self%units%dimensional_output) then
         call self%io%save_field(xh5f=xh5f, grid=self%adam%grid, block_name=bn, ijk=ijk, nijk=nijk, q=q, q_name=q_name)
         return
      endif
      q_out = q
      do v_=1, size(q_out, dim=1)
         q_out(v_,:,:,:) = q_out(v_,:,:,:) * factor(v_)
      enddo
      call self%io%save_field(xh5f=xh5f, grid=self%adam%grid, block_name=bn, ijk=ijk, nijk=nijk, q=q_out, q_name=q_name)
      endsubroutine save_block
   endsubroutine save_xh5f

   subroutine set_divb_seam(self)
   !< Set the seam faces flags of the div(B) history (issue #41, D-10): a block face is a seam face when the forest
   !< registered it in the flux register (a 2:1 AMR or inter-realm seam, `inter_realm_face_register_index /= 0`).
   class(flume_common_object), intent(inout) :: self !< The equation.

   if (allocated(self%divb_seam)) deallocate(self%divb_seam)
   allocate(self%divb_seam(self%blocks_number,6))
   self%divb_seam = 0_I4P
   if (.not.allocated(self%adam%maps%inter_realm_face_register_index)) return
   if (size(self%adam%maps%inter_realm_face_register_index, dim=1) < self%blocks_number) return
   where (self%adam%maps%inter_realm_face_register_index(1:self%blocks_number,1:6) /= 0_I4P) self%divb_seam = 1_I4P
   endsubroutine set_divb_seam

   subroutine set_glm_damping(self)
   !< Set the GLM damping once the grid exists (issue #41, section 3.5): the minimum cell spacing of the realm over the
   !< active directions (MPI-reduced) is the `min-cell` damping length. A no-op without GLM.
   class(flume_common_object), intent(inout) :: self     !< The equation.
   real(R8P)                                 :: min_cell !< Minimum cell spacing.
   integer(I4P)                              :: b, d     !< Counters.

   if (.not.self%physics%mhd%has_glm) return
   min_cell = huge(1._R8P)
   do b=1, self%blocks_number
      do d=1, 3
         if (.not.self%adam%grid%null_xyz(d)) min_cell = min(min_cell, self%adam%field%dxyz(d,b))
      enddo
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, min_cell, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, mpih%error)
   call self%physics%mhd%set_glm_damping(min_cell=min_cell)
   endsubroutine set_glm_damping

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
                             time=self%time%time, time_max=self%time%time_max, adam=self%adam, q=self%q, q_name=name, &
                             length_scale=self%units%length_output(), q_scale=self%output_factors(self%q_name))
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
   pure function ib_cut_spacing(phi_c, phi_m, phi_p, ds) result(ds_cut)
   !< Return the spacing of a fluid cell along one direction, shortened where the solid surface crosses its stencil
   !< (CHASE semantics, issue #35 D-9).
   !<
   !< When the fluid cell (`phi_c < 0`) has one neighbour inside the solid (`phi_m phi_p < 0`), the surface lies at the
   !< distance `delta = -phi_c / (phi_s - phi_c) ds` from the cell centre towards the solid neighbour `s`, and the
   !< spacing becomes `ds / 2 + delta`; otherwise it is `ds`. `phi_s > 0 > phi_c` there, so `phi_s - phi_c > 0` needs no
   !< guard: CHASE added an absolute `1e-12`, a length that broke the scaling covariance of the scheme (issue #49).
   real(R8P), intent(in) :: phi_c  !< Distance function of the cell (negative in the fluid).
   real(R8P), intent(in) :: phi_m  !< Distance function of the minus neighbour.
   real(R8P), intent(in) :: phi_p  !< Distance function of the plus neighbour.
   real(R8P), intent(in) :: ds     !< Spacing.
   real(R8P)             :: ds_cut !< Spacing, cut by the surface.
   !$acc routine seq
   !$omp declare target

   ds_cut = ds
   if (phi_c < 0._R8P .and. phi_m * phi_p < 0._R8P) then
      if (phi_p > 0._R8P) then
         ds_cut = 0.5_R8P * ds - phi_c / (phi_p - phi_c) * ds
      else
         ds_cut = 0.5_R8P * ds - phi_c / (phi_m - phi_c) * ds
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

   subroutine compute_mhd_derived(self, b, derived)
   !< Compute the MHD derived output fields of block `b` (issue #41, section 9) from the host `q` and `q_aux`: total
   !< pressure `pt = p + |B|^2 / 2`, plasma beta `2 p / |B|^2` (`huge` where `B = 0`), `|B|` (every cell), and `div B`
   !< (interior cells, zero on the ghost cells) by the centred finite difference of the div(B) history (the library
   !< derivative, half stencil `[fdv]`, null directions weighted zero). The caller has refreshed the ghost cells.
   class(flume_common_object), intent(in)  :: self                            !< The equation.
   integer(I4P),               intent(in)  :: b                               !< Block index.
   real(R8P),                  intent(out) :: derived(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:) !< pt, beta, bmag, divb.
   real(R8P)                               :: b2                              !< |B|^2.
   real(R8P)                               :: w(3)                            !< Direction weights: 1 active, 0 null.
   real(R8P)                               :: db(3)                           !< dBx/dx, dBy/dy, dBz/dz.
   integer(I4P)                            :: hs                              !< Finite difference half stencil.
   integer(I4P)                            :: i, j, k                         !< Counters.

   hs = self%fdv_half_stencils(1)
   w = merge(0._R8P, 1._R8P, self%adam%grid%null_xyz)
   derived = 0._R8P
   do k=1-self%ngc, self%nk+self%ngc
      do j=1-self%ngc, self%nj+self%ngc
         do i=1-self%ngc, self%ni+self%ngc
            b2 = self%q_aux(IA_BX,i,j,k,b)**2 + self%q_aux(IA_BY,i,j,k,b)**2 + self%q_aux(IA_BZ,i,j,k,b)**2
            derived(1,i,j,k) = self%q_aux(IA_P,i,j,k,b) + 0.5_R8P * b2
            derived(2,i,j,k) = huge(1._R8P)
            if (b2 > 0._R8P) derived(2,i,j,k) = 2._R8P * self%q_aux(IA_P,i,j,k,b) / b2
            derived(3,i,j,k) = sqrt(b2)
         enddo
      enddo
   enddo
   do k=1, self%nk
      do j=1, self%nj
         do i=1, self%ni
            call compute_derivative1_fd_centered(s=hs, ds=self%adam%field%dxyz(1,b), q=self%q(IQ_BX,i-hs:i+hs,j,k,b), &
                                                 dq_ds=db(1))
            call compute_derivative1_fd_centered(s=hs, ds=self%adam%field%dxyz(2,b), q=self%q(IQ_BY,i,j-hs:j+hs,k,b), &
                                                 dq_ds=db(2))
            call compute_derivative1_fd_centered(s=hs, ds=self%adam%field%dxyz(3,b), q=self%q(IQ_BZ,i,j,k-hs:k+hs,b), &
                                                 dq_ds=db(3))
            derived(4,i,j,k) = w(1) * db(1) + w(2) * db(2) + w(3) * db(3)
         enddo
      enddo
   enddo
   endsubroutine compute_mhd_derived

   subroutine compute_q_aux_host(self)
   !< Compute the auxiliary variables of the host `q` on every cell, ghost cells included (output and AMR marking only:
   !< the backends compute their own auxiliary variables in the space operator). The model is selected outside the loops.
   class(flume_common_object), intent(inout) :: self       !< The equation.
   integer(I4P)                              :: b, i, j, k !< Counters.

   select case(self%physics%model)
   case(MODEL_EULER)
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
   case(MODEL_MHD, MODEL_MHD_GLM)
      do b=1, self%blocks_number
         do k=1-self%ngc, self%nk+self%ngc
            do j=1-self%ngc, self%nj+self%ngc
               do i=1-self%ngc, self%ni+self%ngc
                  call mhd_conservative_to_auxiliary(gamma=self%physics%gamma, R=self%physics%R, q=self%q(:,i,j,k,b), &
                                                     qa=self%q_aux(:,i,j,k,b))
               enddo
            enddo
         enddo
      enddo
   case(MODEL_MHD_EGLM)
      do b=1, self%blocks_number
         do k=1-self%ngc, self%nk+self%ngc
            do j=1-self%ngc, self%nj+self%ngc
               do i=1-self%ngc, self%ni+self%ngc
                  call mhd_eglm_conservative_to_auxiliary(gamma=self%physics%gamma, R=self%physics%R,             &
                                                          q=self%q(:,i,j,k,b), qa=self%q_aux(:,i,j,k,b))
               enddo
            enddo
         enddo
      enddo
   case default
      call mpih%error_stop(msg=': no host auxiliary variables for physical model "'//self%physics%physical_model//'"')
   endselect
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

   subroutine check_ngc_number(self)
   !< Check the ghost cells number against the WENO stencil half-width (`weno-riemann`: at least 2, the correction reads
   !< cells i-1 ... i+2).
   class(flume_common_object), intent(in) :: self !< The equation.

   if (self%weno%S > self%ngc) &
      call mpih%error_stop(msg=': [grid].(ngc)='//trim(str(self%ngc))//' is smaller than the WENO stencil half-width '// &
                               trim(str(self%weno%S)))
   if (self%numerics%scheme_space == SCHEME_SPACE_WENO_RIEMANN .and. self%ngc < 2_I4P) &
      call mpih%error_stop(msg=': [grid].(ngc)='//trim(str(self%ngc))//' is smaller than 2, the stencil of the '// &
                               'weno-riemann face flux correction')
   endsubroutine check_ngc_number

   subroutine check_positivity_limiter(self)
   !< Refuse the positivity limiter (issue #47, D-9, D-10) where it cannot work: with the mixed GLM (psi changes B_n
   !< outside the energy, so no backbone is admissible, M3-P0), with a Runge-Kutta scheme that is not SSP (the limiter
   !< makes each forward-Euler step of size dt admissible, and only an SSP scheme is a convex combination of such steps,
   !< all of size at most dt), and with immersed solids (the cut-cell flux difference is not the backbone's). The backends
   !< refuse it on multi-realm runs (the inter-realm seam carries no limiting factor).
   class(flume_common_object), intent(in) :: self !< The equation.

   if (self%numerics%positivity_limiter /= POSITIVITY_LIMITER_CELL) return
   if (self%physics%model == MODEL_MHD_GLM) &
      call mpih%error_stop(msg=': [numerics].(positivity_limiter)=cell is refused with [mhd].(divergence_control)=glm: '//&
                               'use eglm (issue #47, D-10)')
   if (.not.any(self%rk%scheme == [character(18) :: RK_SSP_11, RK_SSP_22, RK_SSP_33, RK_SSP_54])) &
      call mpih%error_stop(msg=': [numerics].(positivity_limiter)=cell needs an SSP Runge-Kutta scheme, not '// &
                               '[runge_kutta].(scheme)='//self%rk%scheme)
   if (self%ib%solids_number > 0_I4P) &
      call mpih%error_stop(msg=': [numerics].(positivity_limiter)=cell is not supported with immersed solids')
   endsubroutine check_positivity_limiter

   subroutine check_weno_scheme(self)
   !< Refuse the centred WENO schemes: the flux splitting calls the upwind primitive only, so a `weno-c-*` scheme would
   !< silently run it at order 2S-1 (issue #47). Refuse the scale-invariant weights on primitive variables: their
   !< descaler is the magnitude of each field, and a velocity or a magnetic field component crossing zero has none of
   !< its own, so the weights would lose accuracy on smooth data (issue #49; the characteristic and conservative variables
   !< take the magnitude of the projected state).
   class(flume_common_object), intent(in) :: self !< The equation.

   if (self%weno%is_centered) &
      call mpih%error_stop(msg=': [weno].(scheme)='//self%weno%scheme//' is a centred scheme: FLUME accepts only the '// &
                               'upwind schemes weno-u-1, weno-u-3, weno-u-5, weno-u-7, weno-u-9')
   if (self%weno%weights == WENO_WEIGHTS_SI .and. self%numerics%scheme_space == SCHEME_SPACE_WENO_RIEMANN .and. &
       self%numerics%reconstruction_variables == RECON_PRIMITIVE)                                                &
      call mpih%error_stop(msg=': [weno].(weights)=si is not available with [numerics].(reconstruction_variables)='// &
                               'primitive: use characteristic variables, or the default weights js')
   endsubroutine check_weno_scheme

   subroutine initialize_riemann_scheme(self)
   !< Check the Riemann solver against the physical model and build the WENO interpolation tables (`weno-riemann` only,
   !< issue #47): the Euler model takes `llf`, `hll`, `hllc`; the MHD models `llf`, `hll`, `hlld`.
   class(flume_common_object), intent(inout) :: self !< The equation.
   logical                                   :: ok   !< Model/solver pair accepted.

   if (self%numerics%scheme_space /= SCHEME_SPACE_WENO_RIEMANN) return
   select case(self%physics%model)
   case(MODEL_EULER)
      ok = any(self%numerics%riemann_solver == [character(4) :: RIEMANN_SOLVER_LLF, RIEMANN_SOLVER_HLL, RIEMANN_SOLVER_HLLC])
   case default
      ok = any(self%numerics%riemann_solver == [character(4) :: RIEMANN_SOLVER_LLF, RIEMANN_SOLVER_HLL, RIEMANN_SOLVER_HLLD])
   endselect
   if (.not.ok) call mpih%error_stop(msg=': [numerics].(riemann_solver)='//self%numerics%riemann_solver//' is not '// &
                                         'available for [physics].(physical_model)='//self%physics%physical_model)
   call self%weno%initialize_interpolation
   endsubroutine initialize_riemann_scheme

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
   self%q_aux_name(1) = 'rho'
   self%q_aux_name(2) = 'u'
   self%q_aux_name(3) = 'v'
   self%q_aux_name(4) = 'w'
   self%q_aux_name(5) = 'p'
   self%q_aux_name(6) = 'T'
   self%q_aux_name(7) = 'H'
   self%q_aux_name(8) = 'a'
   select case(self%physics%model)
   case(MODEL_MHD, MODEL_MHD_GLM, MODEL_MHD_EGLM)
      self%q_name(6) = 'bx'
      self%q_name(7) = 'by'
      self%q_name(8) = 'bz'
      if (self%physics%model /= MODEL_MHD) self%q_name(9) = 'psi'
      ! the auxiliary copies of B share the XH5F file with the conservative bx, by, bz: distinct names
      self%q_aux_name(9)  = 'Bx'
      self%q_aux_name(10) = 'By'
      self%q_aux_name(11) = 'Bz'
   endselect
   do v=1, self%physics%nv
      self%dq_name(v) = 'dq_'//self%q_name(v)%chars()
   enddo
   endsubroutine io_initialize
endmodule adam_flume_common_object
