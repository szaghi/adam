!< ADAM, flux register class definition
module adam_flux_register_object
!< ADAM, flux register class definition
!<
!< Berger-Colella reflux machinery for coarse-fine interfaces — Phase A of
!< [issue #13]. Accumulates per-stage fluxes on registered seam faces from
!< both the coarse and the fine side; at end-of-stage, applies the
!< stage-weighted correction `(Δt/Δx)·(F_coarse_used − Σ F_fine)` to the
!< coarse-side conserved variables.
!<
!< Integrator-agnostic vocabulary: the per-face accumulator's third axis is
!< nominally sized by `n_stages` (the realm's `stages_per_step_forest()`
!< return value). For RK realms this happens to equal `rk%nrk`, but the
!< register API uses the neutral name so a non-RK integrator can register /
!< accumulate / correct against the same machinery without any vocabulary
!< mismatch.
!<
!< **α.r1 (current operating mode): third-axis collapsed to 1.** Under the
!< α end-of-step barrier seam coupling (PRD #16), per-stage RK-weighted
!< reflux refinement (Wang 2018) is dropped to admit asymmetric per-realm
!< stage counts (`K_realm(is) /= K_max`). The accumulator collapses to a
!< single end-of-step bucket: `F_coarse(:,:,1)` and `F_fine_sum(:,:,1)` hold
!< the realm-final-stage face flux, and the correction is applied once per
!< realm per step at `stage == stages_per_step_forest()`. The `n_stages`
!< argument on `register_face` is retained for API stability but **the
!< third-axis allocation is forced to 1 regardless of the value passed**;
!< the per-stage discipline at the call site lives in the caller (the forest
!< passes `n_stages=1`, PRISM's `accumulate_seam_fluxes_fv` passes `stage=1`
!< and gates the accumulation on `stage_idx == rk%nrk`). This is the AMReX
!< `Reflux` cadence; same-K and asymmetric-K both reduce to the same
!< Berger-Colella end-of-step correction.
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
!<   - per stage:        `accumulate_*_flux` collect contributions from
!<                       `compute_residuals_*`
!<   - end of stage:     `apply_reflux` reduces (MPI) and corrects the coarse `q`
!<
!< Memory cost is O(seam-area × nv × n_stages), bounded by the surface area
!< of the AMR boundaries — independent of the interior cell count.
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
   !< `(nv, nface_cells, 1)` under the α.r1 single-stage end-of-step contract
   !< (see module-level note). The third axis is retained dimensionally so
   !< the `accumulate_*_flux(face_index, stage, flux_face)` API survives
   !< unchanged; α-mode callers always pass `stage = 1`. The fine side
   !< accumulates into the same coarse-face skin via 2:1 averaging at the
   !< time of accumulation, so the two arrays are pointwise-comparable.
   integer(I4P) :: seam_kind    = 0_I4P  !< SEAM_KIND_* (set by `register_face`).
   integer(I4P) :: coarse_realm = 0_I4P  !< Realm owning the coarse side (1-based; 0 sentinel = unset).
   integer(I4P) :: coarse_block = 0_I4P  !< Block index on the coarse side.
   integer(I4P) :: coarse_face  = 0_I4P  !< Face code (FACE_X_MAX..FACE_Z_MIN from adam_maps_object).
   integer(I4P) :: fine_realm   = 0_I4P  !< Realm owning the fine side (1-based; 0 sentinel = unset).
   integer(I4P), allocatable :: fine_block(:)        !< Fine blocks covering this coarse face (size 4 for 2:1 refinement in 3D).
   integer(I4P) :: nface_cells  = 0_I4P  !< Coarse cells covered by this face (coarse-side count).
   real(R8P), allocatable :: F_coarse(:,:,:)         !< Coarse-side end-of-step flux: (nv, nface_cells, 1) under α.r1.
   real(R8P), allocatable :: F_fine_sum(:,:,:)       !< Fine-side end-of-step accumulator: (nv, nface_cells, 1) under α.r1.
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
      procedure, pass(self) :: initialize
                                                      !< Allocate `face(:)` with `nfaces` capacity (called once after topology is
                                                      !< known).
      procedure, pass(self) :: register_face          !< Populate one entry in `face(:)`; called by topology pass at init/regrid.
      procedure, pass(self) :: accumulate_coarse_flux
                                                      !< Add a coarse-side per-stage flux to F_coarse on the matching face
                                                      !< (called by compute_residuals_* on coarse-side seam faces).
      procedure, pass(self) :: accumulate_fine_flux
                                                      !< Add a fine-side per-stage flux to F_fine_sum on the matching face
                                                      !< (called by compute_residuals_* on fine-side seam faces).
      procedure, pass(self) :: reduce_fine_sums
                                                      !< MPI-reduce F_fine_sum across ranks so coarse owner has complete sum (Phase
                                                      !< A v1: no-op for replicated forest).
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
                            fine_realm, fine_block, nface_cells, nv, n_stages)
   !< Populate one entry in `face(:)`; called by topology pass at init/regrid.
   !<
   !< For each coarse-fine adjacency discovered during topology population
   !< (in `forest_object%populate_inter_realm_topology` for inter-realm seams,
   !< and in the analogous AMR-jump walk for intra-realm seams), the topology
   !< pass calls this routine to fill in the corresponding `face(face_index)`
   !< entry. The accumulator arrays `F_coarse` and `F_fine_sum` are allocated
   !< here with shape (nv, nface_cells, n_stages) and zeroed.
   !<
   !< Same-resolution case (current rmf-2realm with COUPLING_MIRROR and no
   !< intra-realm AMR jumps): `coarse_realm` and `fine_realm` are the two
   !< sides of the seam; the labels are conventional — both sides carry the
   !< same resolution, both fluxes are nominally equal, and the resulting
   !< correction is expected to be round-off zero in expectation. The
   !< accumulator structure is identical to the true coarse-fine AMR case;
   !< the same-resolution case just exercises the value-zero branch.
   !<
   !< **α.r1 storage convention.** The `n_stages` argument is validated
   !< (`>= 1`) for API hygiene but the third-axis allocation is **forced to
   !< `1`** regardless of the value passed: the single end-of-step Berger-
   !< Colella correction is the same expression whether the realm uses
   !< SSP-RK-3 or SSP-RK-5, and same-K Wang-2018 per-stage refinement is
   !< deferred (see module-level note). Callers under α should pass
   !< `n_stages = 1_I4P` to make intent explicit, but passing the realm's
   !< `stages_per_step_forest()` is also accepted (the surplus depth is
   !< silently dropped). When per-stage storage is reintroduced (γ research
   !< RFC #17 or a future Wang-2018 path), this routine grows a sizing
   !< branch and the contract goes back to per-stage.
   class(flux_register_object), intent(inout)           :: self          !< The register.
   integer(I4P),                intent(in)              :: face_index    !< Index into `face(:)` to populate (1-based, 1..nfaces).
   integer(I4P),                intent(in)              :: seam_kind     !< SEAM_KIND_*.
   integer(I4P),                intent(in)              :: coarse_realm  !< Realm owning the coarse side (1-based).
   integer(I4P),                intent(in)              :: coarse_block  !< Block index on the coarse side.
   integer(I4P),                intent(in)              :: coarse_face   !< Face code (FACE_X_MAX..FACE_Z_MIN).
   integer(I4P),                intent(in)              :: fine_realm    !< Realm owning the fine side (1-based).
   integer(I4P),                intent(in),    optional :: fine_block(:) !< Fine blocks covering the coarse face;
                                                                         !< absent ↔ no fine-side blocks recorded
                                                                         !< (same-resolution mirror seam, Phase A v1).
                                                                         !< Caller must NOT pass an empty array literal
                                                                         !< `[integer(I4P) ::]`: nvfortran 26.x propagates
                                                                         !< the null data pointer through the descriptor
                                                                         !< (debug: "Null pointer for zs$array"; release:
                                                                         !< silent corruption surfacing N kernels later as
                                                                         !< CUDA_ERROR_INVALID_VALUE).
   integer(I4P),                intent(in)              :: nface_cells   !< Coarse cells covered by this face.
   integer(I4P),                intent(in)              :: nv            !< Number of variables (state-vector width).
   integer(I4P),                intent(in)              :: n_stages      !< Realm's `stages_per_step_forest()` (advisory under α.r1:
                                                                          !< the third-axis allocation is forced to `1` regardless).

   if (.not. self%is_initialized_) &
      call mpih%error_stop(msg='flux_register_object%register_face: called before initialize')
   if (face_index < 1_I4P .or. face_index > self%nfaces) &
      call mpih%error_stop(msg='flux_register_object%register_face: face_index out of range')
   if (nface_cells < 1_I4P) &
      call mpih%error_stop(msg='flux_register_object%register_face: nface_cells must be >= 1')
   if (nv < 1_I4P) &
      call mpih%error_stop(msg='flux_register_object%register_face: nv must be >= 1')
   if (n_stages < 1_I4P) &
      call mpih%error_stop(msg='flux_register_object%register_face: n_stages must be >= 1')

   associate(slot => self%face(face_index))
   slot%seam_kind    = seam_kind
   slot%coarse_realm = coarse_realm
   slot%coarse_block = coarse_block
   slot%coarse_face  = coarse_face
   slot%fine_realm   = fine_realm
   slot%nface_cells  = nface_cells
   if (allocated(slot%fine_block)) deallocate(slot%fine_block)
   if (present(fine_block)) then
      if (size(fine_block) > 0_I4P) then
         allocate(slot%fine_block(1:size(fine_block)))
         slot%fine_block = fine_block
      endif
   endif
   ! No allocation ↔ no fine-side blocks recorded for this face.
   if (allocated(slot%F_coarse  )) deallocate(slot%F_coarse  )
   if (allocated(slot%F_fine_sum)) deallocate(slot%F_fine_sum)
   ! α.r1: third axis forced to 1 — see register_face docstring. The
   ! `n_stages` arg above is validated for API hygiene only.
   allocate(slot%F_coarse  (1:nv, 1:nface_cells, 1:1))
   allocate(slot%F_fine_sum(1:nv, 1:nface_cells, 1:1))
   slot%F_coarse   = 0._R8P
   slot%F_fine_sum = 0._R8P
   end associate
   endsubroutine register_face

   subroutine accumulate_coarse_flux(self, face_index, stage, flux_face)
   !< Add a coarse-side per-stage flux contribution to `F_coarse(:,:,stage)`.
   !<
   !< Called by PRISM's `compute_residuals_fv_centered` (option α, #13 §3.2)
   !< when processing a face that has been identified as a coarse-side seam
   !< face. The `face_index` is looked up via a precomputed
   !< `seam_face_index(b, face_code)` table owned by the residual routine;
   !< the lookup is O(1).
   !<
   !< Accumulation is additive: a face may receive contributions from
   !< multiple `accumulate_*_flux` calls within one stage (e.g. if the
   !< FV scheme processes the same face from both the cell-i and cell-(i+1)
   !< side); the register holds the running sum.
   class(flux_register_object), intent(inout) :: self          !< The register.
   integer(I4P),                intent(in)    :: face_index    !< Index into `face(:)` of the seam face being updated.
   integer(I4P),                intent(in)    :: stage         !< Integrator stage 1..n_stages.
   real(R8P),                   intent(in)    :: flux_face(:,:)!< Per-face flux contribution shaped (nv, nface_cells).

   if (face_index < 1_I4P .or. face_index > self%nfaces) &
      call mpih%error_stop(msg='flux_register_object%accumulate_coarse_flux: face_index out of range')
   if (.not. allocated(self%face(face_index)%F_coarse)) &
      call mpih%error_stop(msg='flux_register_object%accumulate_coarse_flux: F_coarse not allocated; was register_face called?')
   if (stage < 1_I4P .or. stage > size(self%face(face_index)%F_coarse, dim=3)) &
      call mpih%error_stop(msg='flux_register_object%accumulate_coarse_flux: stage out of range')
   self%face(face_index)%F_coarse(:,:,stage) = self%face(face_index)%F_coarse(:,:,stage) + flux_face(:,:)
   endsubroutine accumulate_coarse_flux

   subroutine accumulate_fine_flux(self, face_index, stage, flux_face)
   !< Add a fine-side per-stage flux contribution to `F_fine_sum(:,:,stage)`.
   !<
   !< Called by PRISM's `compute_residuals_fv_centered` when processing a
   !< face identified as a fine-side seam face. The contribution is assumed
   !< to be already mapped to the coarse-face skin (same-resolution case:
   !< identity mapping; AMR case: 2:1 averaging of 4 fine faces into one
   !< coarse slot, performed by the caller before this call). Storing the
   !< coarse-skin-shaped contribution directly keeps the register topology
   !< simple at the cost of caller-side averaging discipline for AMR.
   class(flux_register_object), intent(inout) :: self          !< The register.
   integer(I4P),                intent(in)    :: face_index    !< Index into `face(:)` of the seam face being updated.
   integer(I4P),                intent(in)    :: stage         !< Integrator stage 1..n_stages.
   real(R8P),                   intent(in)    :: flux_face(:,:)
                                                               !< Per-face flux contribution shaped (nv, nface_cells); already
                                                               !< mapped to the coarse-face skin.

   if (face_index < 1_I4P .or. face_index > self%nfaces) &
      call mpih%error_stop(msg='flux_register_object%accumulate_fine_flux: face_index out of range')
   if (.not. allocated(self%face(face_index)%F_fine_sum)) &
      call mpih%error_stop(msg='flux_register_object%accumulate_fine_flux: F_fine_sum not allocated; was register_face called?')
   if (stage < 1_I4P .or. stage > size(self%face(face_index)%F_fine_sum, dim=3)) &
      call mpih%error_stop(msg='flux_register_object%accumulate_fine_flux: stage out of range')
   self%face(face_index)%F_fine_sum(:,:,stage) = self%face(face_index)%F_fine_sum(:,:,stage) + flux_face(:,:)
   endsubroutine accumulate_fine_flux

   subroutine reduce_fine_sums(self)
   !< MPI-reduce `F_fine_sum` across ranks so the coarse-side rank has the
   !< full sum for the correction. Phase A v1 of [issue #13] uses the
   !< replicated-forest layout (both ranks own both realms) where every
   !< accumulator is already complete on every rank — this routine is
   !< a no-op in that case. Cross-rank reduce is left to a follow-up
   !< commit that adds disjoint-rank carve-out.
   !<
   !< The reduce happens here (not in `apply_reflux_corrections` on
   !< `forest_object`) because it is purely register-internal: it
   !< exchanges accumulator data between ranks owning copies of the
   !< same face entry, without any realm-side q access.
   class(flux_register_object), intent(inout) :: self !< The register.

   if (.not. self%is_initialized_) return
   if (self%nfaces == 0_I4P)        return
   ! Replicated-forest layout: no cross-rank reduce needed. The skeleton
   ! is in place for the disjoint-rank follow-up.
   endsubroutine reduce_fine_sums

   subroutine reset(self)
   !< Zero `F_coarse` and `F_fine_sum` on every registered face.
   !<
   !< Called by `forest_object%evolve_one_step` at top-of-step. After this
   !< call, every face's accumulators are zero; per-stage flux contributions
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
