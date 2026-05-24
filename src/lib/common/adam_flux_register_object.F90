!< ADAM, flux register class definition
module adam_flux_register_object
!< ADAM, flux register class definition
!<
!< Berger-Colella reflux machinery for coarse-fine interfaces — Phase A of
!< [issue #13]. Accumulates per-substage fluxes on registered seam faces from
!< both the coarse and the fine side; at end-of-substage, applies the
!< stage-weighted correction `(Δt/Δx)·(F_coarse_used − Σ F_fine)` to the
!< coarse-side conserved variables.
!<
!< Two seam kinds share the same accumulator topology:
!<
!<   - `SEAM_KIND_INTRA_REALM_AMR` — fine-coarse adjacency between two AMR
!<     levels inside one realm. Fine and coarse sides typically co-resident
!<     on the same rank (ADAM keeps levels co-resident by construction).
!<   - `SEAM_KIND_INTER_REALM` — face coupling between two realms in a forest.
!<     Fine and coarse sides may live on different ranks (the cross-rank case
!<     defect A of [issue #13]).
!<
!< Lifecycle (driven by `forest_object`):
!<
!<   - init / regrid:    `register_face` populates `face(:)`
!<   - per step (top):   `reset` zeroes accumulators
!<   - per substage:     `accumulate_*_flux` collect contributions from
!<                       `compute_residuals_*`
!<   - end of substage:  `apply_reflux` reduces (MPI) and corrects the coarse `q`
!<
!< Memory cost is O(seam-area × nv × nrk), bounded by the surface area of the
!< AMR boundaries — independent of the interior cell count.
!<
!< Architectural notes (constraints from [issue #10] R1/R2):
!<   - All accumulators are intrinsic-typed arrays; no derived-type pointer
!<     components reachable from kernels.
!<   - The register owns no field data — only the per-face flux skin. Field
!<     pointers (`q`, `dq`) are threaded in by callers, not stored.
!<
!< This is the **skeleton commit** of Phase A: every TBP body is a no-op
!< (or, for `reset`, a trivial allocation guard). The wiring into
!< `forest_object` (see #13 §3.2) and into PRISM's `compute_residuals_*`
!< (option α, #13 §3.2) lands in follow-up commits.

! ADAM classes, libraries, parameters
use :: adam_parameters
! ADAM singleton objects
use :: adam_mpih_global, only : mpih
! third party modules
use :: penf

implicit none
private
public :: flux_register_object
public :: flux_register_face_t
public :: SEAM_KIND_INTRA_REALM_AMR, SEAM_KIND_INTER_REALM

!< Seam-kind codes (see [[flux_register_face_t]]%seam_kind).
!<
!< Encoded as small integers (rather than a Fortran enum) for the same reason
!< the FACE_* codes in `adam_maps_object` are: they participate in checkpoint
!< dumps and INI-level diagnostics, so stable integer values are wanted.
integer(I4P), parameter :: SEAM_KIND_INTRA_REALM_AMR = 1_I4P !< Coarse-fine AMR jump inside one realm.
integer(I4P), parameter :: SEAM_KIND_INTER_REALM    = 2_I4P !< Inter-realm seam between two realms in a forest.

type :: flux_register_face_t
   !< Per-seam-face flux accumulator.
   !<
   !< One entry per coarse face touching a fine-side neighbourhood: either a
   !< fine-level patch (`SEAM_KIND_INTRA_REALM_AMR`) or another realm
   !< (`SEAM_KIND_INTER_REALM`). For a 2:1 refinement in 3D, the coarse face
   !< is covered by 4 fine faces; `fine_block(:)` lists them.
   !<
   !< Both accumulators (`F_coarse`, `F_fine_sum`) are sized
   !< `(nv, nface_cells, nrk)` — the per-substage flux on the **coarse-face
   !< skin** (the nface_cells coarse cells that lie under this coarse face).
   !< The fine side accumulates into the same coarse-face skin via 2:1
   !< averaging at the time of accumulation, so the two arrays are
   !< pointwise-comparable.
   integer(I4P) :: seam_kind    = 0_I4P  !< SEAM_KIND_* (set by `register_face`).
   integer(I4P) :: coarse_realm = 0_I4P  !< Realm owning the coarse side (1-based; 0 sentinel = unset).
   integer(I4P) :: coarse_block = 0_I4P  !< Block index on the coarse side.
   integer(I4P) :: coarse_face  = 0_I4P  !< Face code (FACE_X_MAX..FACE_Z_MIN from adam_maps_object).
   integer(I4P) :: fine_realm   = 0_I4P  !< Realm owning the fine side (1-based; 0 sentinel = unset).
   integer(I4P), allocatable :: fine_block(:)        !< Fine blocks covering this coarse face (size 4 for 2:1 refinement in 3D).
   integer(I4P) :: nface_cells  = 0_I4P  !< Coarse cells covered by this face (coarse-side count).
   real(R8P), allocatable :: F_coarse(:,:,:)         !< Coarse-side flux per substage: (nv, nface_cells, nrk).
   real(R8P), allocatable :: F_fine_sum(:,:,:)       !< Fine-side accumulator per substage: (nv, nface_cells, nrk).
endtype flux_register_face_t

type :: flux_register_object
   !< Flux register class definition.
   !<
   !< Owns the collection of registered seam faces. Lives on the forest (one
   !< instance per forest, shared across realms) so that the cross-realm
   !< reduce in `apply_reflux` has a single ownership point. Intra-realm AMR
   !< faces and inter-realm faces share the same `face(:)` array — they
   !< differ only in their `seam_kind` tag.
   logical :: is_initialized_ = .false. !< Flag: register has been initialized.
   integer(I4P) :: nfaces            = 0_I4P !< Number of registered faces (size of `face(:)` after `initialize`).
   type(flux_register_face_t), allocatable :: face(:) !< All registered seam faces; unallocated until `initialize`.
   contains
      ! public methods
      procedure, pass(self) :: initialize             !< Allocate `face(:)` with `nfaces` capacity (called once after topology is known).
      procedure, pass(self) :: register_face          !< Populate one entry in `face(:)`; called by topology pass at init/regrid.
      procedure, pass(self) :: accumulate_coarse_flux !< Add a coarse-side per-substage flux to F_coarse on the matching face (called by compute_residuals_* on coarse-side seam faces).
      procedure, pass(self) :: accumulate_fine_flux   !< Add a fine-side per-substage flux to F_fine_sum on the matching face (called by compute_residuals_* on fine-side seam faces).
      procedure, pass(self) :: apply_reflux           !< Reduce F_fine_sum across ranks if needed, then correct the coarse-side q with the (F_coarse − F_fine_sum) imbalance.
      procedure, pass(self) :: reset                  !< Zero F_coarse and F_fine_sum on every face; called at top-of-step.
      procedure, pass(self) :: destroy                !< Deallocate `face(:)` and reset state to uninitialized.
endtype flux_register_object

contains
   ! public methods
   subroutine initialize(self, nfaces)
   !< Allocate `face(:)` with `nfaces` capacity.
   !<
   !< Called once after the topology pass has counted how many seam faces
   !< exist. `register_face` is then called `nfaces` times to populate the
   !< entries.
   !<
   !< Idempotent: a re-initialization with the same forest topology yields
   !< the same allocation; under a regrid that changes block layouts, the
   !< caller is expected to call `destroy` first, then `initialize` with
   !< the new face count. (Phase A v1 has no regrid hook; this is the
   !< schema-reserved path for AMR-active follow-ups.)
   !<
   !< Calling with `nfaces == 0` is valid: it flips `is_initialized_` so
   !< the lifecycle hooks (`reset`, `apply_reflux`) become safe no-ops on
   !< the empty register. This is the N=1 single-realm path.
   class(flux_register_object), intent(inout) :: self    !< The register.
   integer(I4P),                intent(in)    :: nfaces  !< Number of seam faces to register.

   if (self%is_initialized_) call self%destroy
   self%nfaces = nfaces
   if (nfaces > 0_I4P) allocate(self%face(1:nfaces))
   self%is_initialized_ = .true.
   endsubroutine initialize

   subroutine register_face(self, face_index, seam_kind, coarse_realm, coarse_block, coarse_face, &
                            fine_realm, fine_block, nface_cells, nv, nrk)
   !< Populate one entry in `face(:)`; called by topology pass at init/regrid.
   !<
   !< For each coarse-fine adjacency discovered during topology population
   !< (in `forest_object%populate_inter_realm_topology` for inter-realm seams,
   !< and in the analogous AMR-jump walk for intra-realm seams), the topology
   !< pass calls this routine to fill in the corresponding `face(face_index)`
   !< entry. The accumulator arrays `F_coarse` and `F_fine_sum` are allocated
   !< here with shape (nv, nface_cells, nrk) and zeroed.
   !<
   !< Same-resolution case (current rmf-2realm with COUPLING_MIRROR and no
   !< intra-realm AMR jumps): `coarse_realm` and `fine_realm` are the two
   !< sides of the seam; the labels are conventional — both sides carry the
   !< same resolution, both fluxes are nominally equal, and the resulting
   !< correction is expected to be round-off zero in expectation. The
   !< accumulator structure is identical to the true coarse-fine AMR case;
   !< the same-resolution case just exercises the value-zero branch.
   class(flux_register_object), intent(inout) :: self         !< The register.
   integer(I4P),                intent(in)    :: face_index   !< Index into `face(:)` to populate (1-based, 1..nfaces).
   integer(I4P),                intent(in)    :: seam_kind    !< SEAM_KIND_*.
   integer(I4P),                intent(in)    :: coarse_realm !< Realm owning the coarse side (1-based).
   integer(I4P),                intent(in)    :: coarse_block !< Block index on the coarse side.
   integer(I4P),                intent(in)    :: coarse_face  !< Face code (FACE_X_MAX..FACE_Z_MIN).
   integer(I4P),                intent(in)    :: fine_realm   !< Realm owning the fine side (1-based).
   integer(I4P),                intent(in)    :: fine_block(:)!< Fine blocks covering the coarse face.
   integer(I4P),                intent(in)    :: nface_cells  !< Coarse cells covered by this face.
   integer(I4P),                intent(in)    :: nv           !< Number of variables (state-vector width).
   integer(I4P),                intent(in)    :: nrk          !< Number of RK substages.

   if (.not. self%is_initialized_) &
      call mpih%error_stop(msg='flux_register_object%register_face: called before initialize')
   if (face_index < 1_I4P .or. face_index > self%nfaces) &
      call mpih%error_stop(msg='flux_register_object%register_face: face_index out of range')
   if (nface_cells < 1_I4P) &
      call mpih%error_stop(msg='flux_register_object%register_face: nface_cells must be >= 1')
   if (nv < 1_I4P) &
      call mpih%error_stop(msg='flux_register_object%register_face: nv must be >= 1')
   if (nrk < 1_I4P) &
      call mpih%error_stop(msg='flux_register_object%register_face: nrk must be >= 1')

   associate(slot => self%face(face_index))
   slot%seam_kind    = seam_kind
   slot%coarse_realm = coarse_realm
   slot%coarse_block = coarse_block
   slot%coarse_face  = coarse_face
   slot%fine_realm   = fine_realm
   slot%nface_cells  = nface_cells
   if (allocated(slot%fine_block)) deallocate(slot%fine_block)
   allocate(slot%fine_block(1:size(fine_block, dim=1)))
   slot%fine_block = fine_block
   if (allocated(slot%F_coarse  )) deallocate(slot%F_coarse  )
   if (allocated(slot%F_fine_sum)) deallocate(slot%F_fine_sum)
   allocate(slot%F_coarse  (1:nv, 1:nface_cells, 1:nrk))
   allocate(slot%F_fine_sum(1:nv, 1:nface_cells, 1:nrk))
   slot%F_coarse   = 0._R8P
   slot%F_fine_sum = 0._R8P
   end associate
   endsubroutine register_face

   subroutine accumulate_coarse_flux(self, face_index, substage, flux_face)
   !< Add a coarse-side per-substage flux contribution to `F_coarse(:,:,substage)`.
   !<
   !< Called by `compute_residuals_*` (option α, #13 §3.2) when a face being
   !< processed has been identified as a coarse-side seam face. The `face_index`
   !< is looked up via a precomputed `is_seam_face(b, face)` mask owned by the
   !< residual routine; finding the index is O(1).
   !<
   !< **Skeleton commit:** body is a no-op pending the accumulation-wiring
   !< commit.
   class(flux_register_object), intent(inout) :: self          !< The register.
   integer(I4P),                intent(in)    :: face_index    !< Index into `face(:)` of the seam face being updated.
   integer(I4P),                intent(in)    :: substage      !< RK substage 1..nrk.
   real(R8P),                   intent(in)    :: flux_face(:,:)!< Per-face flux contribution shaped (nv, nface_cells).

   associate(unused_self => self, unused_face_index => face_index, &
             unused_substage => substage, unused_flux => flux_face )
   endassociate
   endsubroutine accumulate_coarse_flux

   subroutine accumulate_fine_flux(self, face_index, substage, flux_face)
   !< Add a fine-side per-substage flux contribution to `F_fine_sum(:,:,substage)`.
   !<
   !< Called by `compute_residuals_*` when a face being processed has been
   !< identified as a fine-side seam face. The fine-side contribution is
   !< already mapped to the coarse-face skin (i.e. the fine routine sums its 4
   !< fine-face fluxes onto a single coarse-face slot via 2:1 averaging before
   !< calling this routine, or the routine itself does the mapping — design
   !< decision deferred to the accumulation-wiring commit).
   !<
   !< **Skeleton commit:** body is a no-op pending the accumulation-wiring
   !< commit.
   class(flux_register_object), intent(inout) :: self          !< The register.
   integer(I4P),                intent(in)    :: face_index    !< Index into `face(:)` of the seam face being updated.
   integer(I4P),                intent(in)    :: substage      !< RK substage 1..nrk.
   real(R8P),                   intent(in)    :: flux_face(:,:)!< Per-face flux contribution shaped (nv, nface_cells); already mapped to the coarse-face skin.

   associate(unused_self => self, unused_face_index => face_index, &
             unused_substage => substage, unused_flux => flux_face )
   endassociate
   endsubroutine accumulate_fine_flux

   subroutine apply_reflux(self, substage, dt, dx_coarse, weight)
   !< Reduce `F_fine_sum` across ranks (if needed) and apply the Berger-Colella
   !< correction to the coarse-side conserved variables.
   !<
   !< Formula (Berger & Colella 1989 §4; Wang et al. 2018 Eq. 23 for the
   !< RK-stage-weighted form):
   !<
   !<     q_coarse(:, coarse cell on seam face)
   !<         +=  weight * (dt / dx_coarse) * ( F_coarse(:, cell, s) - F_fine_sum(:, cell, s) )
   !<
   !< where `weight` is the SSP-RK accumulation weight for substage `s`, read
   !< from `adam_rk_object%ark(s)`. For SSP-RK 3 (`adam_rk_object%nrk = 3`)
   !< the weights decompose to 1 per substage; for SSP-RK 5 they vary.
   !<
   !< Called by `forest_object%evolve_one_step` between substages, after
   !< `compute_residuals_*` has produced the fluxes but before the next ghost
   !< exchange. The cross-rank MPI reduce of `F_fine_sum` happens inside this
   !< routine (separate communicator pattern from `comm_map_send_ghost`).
   !<
   !< **Skeleton commit:** body is a no-op. The reflux correction is wired in
   !< the accumulation-wiring + correction-application follow-up commit.
   class(flux_register_object), intent(inout) :: self      !< The register.
   integer(I4P),                intent(in)    :: substage  !< RK substage 1..nrk.
   real(R8P),                   intent(in)    :: dt        !< Time step.
   real(R8P),                   intent(in)    :: dx_coarse !< Coarse cell width in the face-normal direction.
   real(R8P),                   intent(in)    :: weight    !< RK substage accumulation weight (from adam_rk_object%ark).

   associate(unused_self => self, unused_substage => substage, &
             unused_dt   => dt,   unused_dx       => dx_coarse,&
             unused_w    => weight                              )
   endassociate
   endsubroutine apply_reflux

   subroutine reset(self)
   !< Zero `F_coarse` and `F_fine_sum` on every registered face.
   !<
   !< Called by `forest_object%evolve_one_step` at top-of-step. After this
   !< call, every face's accumulators are zero; substage flux contributions
   !< land into clean slots.
   !<
   !< **Skeleton commit:** safe no-op when `face(:)` is unallocated (the
   !< case in this commit, since `initialize` does not allocate yet).
   class(flux_register_object), intent(inout) :: self !< The register.
   integer(I4P)                               :: i    !< Counter over faces.

   if (.not. allocated(self%face)) return
   do i=1, self%nfaces
      if (allocated(self%face(i)%F_coarse))   self%face(i)%F_coarse   = 0._R8P
      if (allocated(self%face(i)%F_fine_sum)) self%face(i)%F_fine_sum = 0._R8P
   enddo
   endsubroutine reset

   subroutine destroy(self)
   !< Deallocate `face(:)` and reset state to uninitialized.
   class(flux_register_object), intent(inout) :: self !< The register.
   integer(I4P)                               :: i    !< Counter over faces.

   if (allocated(self%face)) then
      do i=1, size(self%face, dim=1)
         if (allocated(self%face(i)%fine_block)) deallocate(self%face(i)%fine_block)
         if (allocated(self%face(i)%F_coarse  )) deallocate(self%face(i)%F_coarse  )
         if (allocated(self%face(i)%F_fine_sum)) deallocate(self%face(i)%F_fine_sum)
      enddo
      deallocate(self%face)
   endif
   self%nfaces = 0_I4P
   self%is_initialized_ = .false.
   endsubroutine destroy
endmodule adam_flux_register_object
