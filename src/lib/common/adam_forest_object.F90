!< ADAM, forest class definition — orchestrator of a forest of realms.
module adam_forest_object
!< ADAM, forest class definition — orchestrator of a forest of realms.
!<
!< The **forest** tends an array of realms ([[realm_object]] extensions, e.g.
!< `prism_cpu_object`, `prism_fnl_object`). It is a **behavior-only** class:
!< it owns no derived-type state. Each realm lives in the program driver as a
!< concrete monomorphic array (`type(prism_cpu_object) :: realm(N)`); the
!< forest receives the array as `class(realm_object), intent(inout) :: realm(:)`
!< and orchestrates inter-realm operations:
!<
!<   * sequence per-realm initialize / finalize calls
!<   * reduce per-realm dt to a global dt (min reduction)
!<   * iterate per-realm timestep advance
!<   * iterate per-realm post-step diagnostics / IO
!<   * reduce per-realm termination predicate to a global done (AND reduction)
!<
!< The forest NEVER reaches inside a realm's private state; it only invokes
!< the realm-side TBPs that carry the **`_forest`** suffix. Together these
!< TBPs form the orchestrator contract (see [[adam_realm_object]]).
!<
!< Class-with-TBPs (not module-of-routines) so future forest-level
!< configuration (MPI sub-communicator topology, inter-realm
!< connectivity descriptor) can be added as intrinsic-typed state
!< without a module API break.
!<
!< See `docs/guide/forest.md` for the conceptual overview of the
!< multi-realm machinery (manifest schema, α/β cadence trade-offs, phase
!< cycle) and `src/lib/common/README.md` → "Forest orchestration" for the
!< library-developer contract surface.

use :: adam_realm_object,         only : realm_object
use :: adam_maps_object,          only : inter_realm_neighbor_t,                                                  &
                                         FACE_X_MAX, FACE_X_MIN, FACE_Y_MAX, FACE_Y_MIN, FACE_Z_MAX, FACE_Z_MIN, &
                                         CADENCE_END_OF_STEP, CADENCE_STAGE_COINCIDENT,                          &
                                         COUPLING_MIRROR, COUPLING_REFINED, face_axis_sign
use :: adam_tree_object,          only : tree_iterator_object, NODE_MORE_REFINED, NODE_LESS_REFINED
use :: adam_tree_node_object,     only : tree_node_object
use :: adam_forest_manifest,      only : forest_manifest_t, forest_face_pair_t
use :: adam_flux_register_object, only : flux_register_object, SEAM_KIND_INTER_REALM, SEAM_KIND_INTER_REALM_REFINED, &
                                         SEAM_KIND_INTRA_REALM_AMR
use :: adam_parameters,           only : BC_SEAM, FEC_1_6_ARRAY
use :: adam_seam_exchange,        only : seam_fill, SEAM_ROW_COPY, SEAM_ROW_INTERPOLATE, SEAM_ROW_RESTRICT
use :: adam_seam_interpolation_library, only : seam_meta_pack, seam_shift_anchor_pos, seam_tricubic_centered_pos, &
                                               seam_compatible_centered_pos
use :: adam_globals,              only : mpih
use :: mpi
use :: penf

implicit none
private
public :: forest_object

type :: forest_object
   !< Behavior-only orchestrator of an array of realms.
   integer(I4P)               :: n = 0_I4P     !< Number of realms in the forest (set by initialize from size(realm)).
   type(flux_register_object) :: flux_register !< Coarse-fine interface reflux machinery.
   contains
      ! public methods
      ! initialize/finalize
      procedure, pass(self) :: initialize               !< Sequence each realm's initialize_forest at startup (single shared INI).
      procedure, pass(self) :: initialize_from_manifest !< Like initialize, but each realm reads its own INI from a forest manifest.
      procedure, pass(self) :: finalize                 !< Sequence each realm's finalize_forest at shutdown.
      ! orchestrating methods
      procedure, pass(self) :: compute_global_dt      !< Min-reduce each realm's compute_local_dt_forest across the forest.
      procedure, pass(self) :: evolve_one_step        !< Iterate realm(:)%advance_one_step_forest(dt) for one global timestep.
      ! Inter-realm seam refresh runs inside evolve_one_step: β seams at every
      ! substage (Phase 2), α seams once at end-of-step (Phase 5).
      procedure, pass(self) :: is_done                !< AND-reduce each realm's is_done_forest across the forest.
      procedure, pass(self) :: post_step              !< Iterate realm(:)%post_step_forest for the per-step diagnostics/IO block.
      procedure, pass(self) :: simulate               !< Main entry point (single shared INI): drive the full simulation.
      procedure, pass(self) :: simulate_from_manifest !< Main entry point (per-realm INIs via forest manifest).
      ! private methods
      procedure, pass(self), private :: populate_inter_realm_topology   !< Translate manifest face-pairs into maps of neighbors.
      procedure, pass(self), private :: register_intra_realm_amr_seams !< Register intra-realm AMR coarse-fine faces in the flux
         !< register.
      procedure, pass(self), private :: apply_reflux_corrections        !< Apply Berger-Colella reflux to coarse-side.
endtype forest_object

type :: seam_rows_t
   !< Growable list of integer rows, used while the seam rows are enumerated.
   integer(I4P), allocatable :: row(:,:) !< Rows (n, columns).
   integer(I4P)              :: n = 0    !< Rows in use.
endtype seam_rows_t

contains
   ! public methods

   ! initialize/finalize
   subroutine initialize(self, realm, filename)
   !< Initialize the forest and every realm it tends.
   class(forest_object), intent(inout) :: self     !< The forest.
   class(realm_object),  intent(inout) :: realm(:) !< The realms to initialize.
   character(*),         intent(in)    :: filename !< Input parameters file name (shared across realms for v1).
   integer(I4P)                        :: is       !< Realm index.

   if (int(size(realm), I4P) > 1_I4P) &
      call mpih%error_stop(msg='forest_object%initialize: multi-realm forest requires initialize_from_manifest')
   self%n = int(size(realm), I4P)
   do is = 1, self%n
      ! Set the realm's self-aware forest position BEFORE invoking the per-
      ! realm initialize; downstream forest TBPs (apply_reflux_to_stage_forest,
      ! diagnostic prefixes) read self%realm_index instead of having `is`
      ! plumbed through every dispatch.
      realm(is)%realm_index = is
      call realm(is)%initialize_forest(filename=filename, realms_number=self%n)
   enddo
   ! Register intra-realm AMR coarse-fine faces. In the manifest-less (N=1) path
   ! there is no inter-realm seam pass, so this also owns the flux-register
   ! initialization (sizing it to the intra-realm AMR face count, or to zero when
   ! the realm has no coarse-fine jumps). Issue #13 §7.5 M2.
   call self%register_intra_realm_amr_seams(realm=realm)
   endsubroutine initialize

   subroutine initialize_from_manifest(self, realm, manifest)
   !< Initialize the forest and every realm using per-realm INIs from a manifest.
   !<
   !< Like `initialize` but each realm receives its OWN INI path (`manifest%realm_ini(is)`) instead of a shared filename. After all
   !< realms are initialized, translates the manifest's face-pair list into per-realm `maps%inter_realm_neighbors` entries.
   !<
   !< The driver MUST allocate `realm(size = manifest%realms_number)` with the concrete app type before calling this — the forest
   !< does not !< allocate the realm array (it cannot, since realm_object is abstract and each app has its own extension).
   class(forest_object),     intent(inout) :: self     !< The forest.
   class(realm_object),      intent(inout) :: realm(:) !< The realms to initialize.
   type(forest_manifest_t),  intent(in)    :: manifest !< Parsed manifest (per-realm INI paths + topology).
   integer(I4P)                            :: is       !< Realm index.

   if (size(realm) /= manifest%realms_number) &
      call mpih%error_stop(msg='forest_object%initialize_from_manifest: size(realm) /= manifest%realms_number')
   self%n = int(size(realm), I4P)
   do is = 1, self%n
      ! Set the realm's self-aware forest position BEFORE invoking the per-
      ! realm initialize; downstream forest TBPs (apply_reflux_to_stage_forest,
      ! diagnostic prefixes) read self%realm_index instead of having `is`
      ! plumbed through every dispatch.
      realm(is)%realm_index = is
      call realm(is)%initialize_forest(filename=trim(manifest%realm_ini(is)), realms_number=self%n)
   enddo
   call self%populate_inter_realm_topology(realm, manifest)
   call check_beta_admissibility(realm=realm, manifest=manifest)
   endsubroutine initialize_from_manifest

   subroutine finalize(self, realm)
   !< Shut down the forest and every realm it tends.
   !<
   !< Iterates `realm(is)%finalize_forest` in increasing index order; each realm closes its IO files, releases its
   !< resources, and finalizes its MPI handler.
   class(forest_object), intent(in)    :: self     !< The forest.
   class(realm_object),  intent(inout) :: realm(:) !< The realms to shut down.
   integer(I4P)                        :: is       !< Realm index.

   do is = 1, int(size(realm), I4P)
      call realm(is)%finalize_forest
   enddo
   ! End-of-run marker. Every realm has run finalize_forest above, so all
   ! checkpoints and residual/history files are flushed to disk. Emit a stable
   ! sentinel to stdout *before* the process-global MPI_Finalize below, which on
   ! some GPU-aware-MPI stacks (ROCm/UCX) can deadlock in teardown *after* the
   ! run is otherwise complete. Reaching here means the simulation finished
   ! cleanly; the regression harness keys on this line to distinguish "done,
   ! then hung in finalize" (terminate early) from "still running" (wait).
   call mpih%print_message('ADAM run complete: all realms finalized, entering MPI_Finalize')
   ! MPI_FINALIZE is process-global: run it ONCE here, after every realm has done its
   ! MPI-using teardown above — not per realm inside finalize_forest, which would tear
   ! MPI down while later realms still need it.
   if (size(realm) >= 1) call realm(1)%finalize_mpi_forest
   endsubroutine finalize

   ! orchestrating methods
   subroutine compute_global_dt(self, realm, dt)
   !< Compute the forest-global stability-limited dt.
   !<
   !< Each realm reports its local dt via `compute_local_dt_forest`; the forest takes the min across all realms
   !< (intra-process) and then across all MPI ranks (`MPI_ALLREDUCE` on `MPI_COMM_WORLD`). Returns
   !< the bit-identical global min every rank should advance by.
   class(forest_object), intent(in)    :: self     !< The forest.
   class(realm_object),  intent(inout) :: realm(:) !< The realms to query.
   real(R8P),            intent(out)   :: dt       !< Global stability-limited dt.
   real(R8P)                           :: dt_local !< Per-realm local dt.
   integer(I4P)                        :: is       !< Realm index.
   integer(I4P)                        :: ierr     !< MPI error code.

   dt = huge(0._R8P)
   do is = 1, int(size(realm), I4P)
      ! For N>1 each realm needs its singleton shims pointing at its own
      ! components before compute_local_dt_forest reads through them.
      call realm(is)%compute_local_dt_forest(dt_local=dt_local)
      dt = min(dt, dt_local)
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, dt, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
   endsubroutine compute_global_dt

   subroutine evolve_one_step(self, realm)
   !< Advance every realm by one global timestep — α end-of-step barrier semantics.
   !<
   !< N=1 fast path: realm owns the whole step (advance_one_step_forest).
   !< No seam machinery is exercised (there are no peer seams to refresh).
   !<
   !< N>1 multi-realm path: per-seam coupling cadence (α / β).
   !<
   !< The forest drives an integrator-agnostic K-stage clock with
   !< `K_realm(is) = realm(is)%stages_per_step_forest()` and
   !< `K_max = max(K_realm)`. Asymmetric per-realm K is a first-class
   !< operating mode under α; the per-stage TBPs (begin/end_stage_forest)
   !< are gated behind `k <= K_realm(is)` so a realm with K < K_max
   !< no-ops the trailing stages.
   !<
   !< Each inter-realm seam carries a `coupling_cadence` from the manifest
   !< (cached on the realm as `seam_local_cadence(p)`):
   !<
   !< α (CADENCE_END_OF_STEP, default; AMReX-aligned coarse-fine convention)
   !<   Mid-step peer ghosts are INTENTIONALLY STALE-BY-ONE-STEP. During
   !<   stages 1..K, each realm reads the peer ghosts established by the
   !<   previous end-of-step exchange (or by the initial-condition seam
   !<   fill at populate_inter_realm_topology time, for the first step).
   !<   This is structurally identical to AMReX's `FillCoarsePatch`
   !<   reading coarse `t^n` data during fine sub-steps (Berger-Oliger
   !<   1984; AMReX_Amr.cpp::timeStep). Phase 5 (below) is the α seam
   !<   coherence boundary: once per step, after `close_step_forest`,
   !<   every realm's `stage_active == 0` and the receiver reads peer's
   !<   committed `q`. α admits asymmetric per-realm K (the K-equality
   !<   guard is removed). Cost: first-order seam coupling in time.
   !<
   !< β (CADENCE_STAGE_COINCIDENT, opt-in)
   !<   Peer ghosts are refreshed once per RK substage, inside the K loop
   !<   (Phase 2 below), before `end_stage_forest` reads them. Admissible
   !<   only when both endpoint realms agree on (scheme_time, rk_scheme,
   !<   nv, K) — enforced at forest init by `check_beta_admissibility`.
   !<   Recovers bit-equivalence to a monolithic single-realm run on the
   !<   union grid when admissible. Seam coupling order matches the
   !<   per-realm interior order. Spatial operator MAY differ per realm.
   !<
   !< Per-seam selection: the same forest may carry α seams and β seams
   !< simultaneously. Phase 2 iterates seams and fires only for β; Phase 5
   !< iterates seams and fires only for α. A seam is filled exactly once
   !< per step under either cadence (Phase 2 may fire K times for β, but
   !< at successive substages, never duplicating the same substage).
   !<
   !< Reflux cadence (α.r1):
   !<
   !<   The flux register's third axis is collapsed to 1. PRISM realms
   !<   gate `apply_reflux_to_stage_forest` and `accumulate_seam_fluxes_fv`
   !<   on the realm's own final RK substage (`stage == rk%nrk`). Earlier
   !<   substages no-op. The mid-step `apply_reflux_corrections` call below
   !<   therefore fires for every k, but only the k == K_realm(is)
   !<   invocation does real work on realm `is`. Independent of α/β: β
   !<   does not restore Wang 2018 per-stage RK-weighted reflux (deferred).
   !<
   !< Phase outline:
   !<
   !<   Phase 0 — open_step_forest (per-realm prologue)
   !<   For k = 1..K_max:
   !<     Phase 1 — begin_stage_forest (per-realm, gated by k <= K_realm(is))
   !<     Phase 2 — fill_seam_from_peer_forest (per-seam, β-gated)
   !<     Phase 3 — end_stage_forest (per-realm, gated by k <= K_realm(is))
   !<               + reduce_fine_sums + apply_reflux_corrections
   !<               (reflux body is α.r1 end-of-step gated inside the realm)
   !<   Phase 4 — close_step_forest (per-realm epilogue)
   !<   Phase 5 — fill_seam_from_peer_forest (per-seam, α-gated)
   !<
   !< LOAD-BEARING INVARIANT (β): Phase 2 must complete on ALL realms
   !< before Phase 3 starts on ANY realm. The orchestrator's serial inner
   !< loops within a rank give this for free under the Phase-A
   !< replicated-forest layout. If a future refactor interleaves Phase 2
   !< and Phase 3, the read-after-overwrite race returns.
   !<
   !< Future per-realm-dt subcycling (Berger-Colella "case 4") and γ
   !< (dense-output peer reads) are deferred. β with asymmetric K would
   !< require γ-class interpolation; not in scope.
   class(forest_object), intent(inout)         :: self        !< The forest.
   class(realm_object),  intent(inout), target :: realm(:)    !< The realms to advance.
   real(R8P)                                   :: dt          !< Global timestep size.
   integer(I4P)                                :: is, p, k    !< Realm, peer, stage indices.
   integer(I4P), allocatable                   :: K_realm(:)  !< Per-realm stage counts (1..size(realm)).
   integer(I4P)                                :: K_max       !< Forest-wide max stage count (multi-realm path).

   call self%flux_register%reset
   call self%compute_global_dt(realm=realm, dt=dt)
   ! N=1 fast path: a single realm with NO coarse-fine seam faces owns its whole
   ! step via the fused `advance_one_step_forest` (= self%integrate). A single
   ! realm that DOES carry intra-realm AMR seam faces must go through the staged
   ! loop below instead, so the FV residual accumulates seam fluxes into the
   ! register (end_stage_forest threads it in) and the end-of-step reflux fires
   ! (apply_reflux_corrections). The `flux_register%nfaces > 0` test keeps every
   ! existing non-AMR single-realm case on the bit-identical fast path. #13 §7.5 M3.
   if (int(size(realm), I4P) == 1_I4P .and. self%flux_register%nfaces == 0_I4P) then
      call realm(1)%advance_one_step_forest(dt=dt)
   else
      do is = 1_I4P, int(size(realm), I4P)
         call realm(is)%open_step_forest(dt=dt)
      enddo
      ! Query each realm's stage count ONCE and cache it in K_realm(:); the
      ! per-stage gates below test against this cache, not against repeated
      ! TBP dispatches.
      allocate(K_realm(size(realm)))
      do is = 1_I4P, int(size(realm), I4P)
         K_realm(is) = realm(is)%stages_per_step_forest()
      enddo
      K_max = maxval(K_realm)
      ! α: no K-equality guard. Asymmetric K is a first-class mode; the
      ! per-stage gates below let a realm with K < K_max no-op trailing stages.
      do k = 1_I4P, K_max
         ! Phase 1 — open stage k on each participating realm; each realm sets stage_active=k.
         do is = 1_I4P, int(size(realm), I4P)
            if (k > K_realm(is)) cycle
            call realm(is)%begin_stage_forest(k=k, K_total=K_max, dt=dt)
         enddo

         ! Phase 2 — per-seam mid-step inter-realm seam fill (β).
         !
         ! Fires only for seams whose manifest `coupling_cadence` is
         ! `CADENCE_STAGE_COINCIDENT`. At this point all participating
         ! realms have completed `begin_stage_forest(k)` (their stage-k
         ! interior buffer is written), so reading peer's stage-k interior
         ! is well-defined. The TBP's buffer-selection logic in
         ! `fill_seam_from_peer_forest` reads peer's `rk%q_rk(:,...,
         ! peer%stage_active)` = peer's stage-k slice and writes self's
         ! stage-k ghosts. Race-free under the Phase-A replicated-forest
         ! layout (serial inner loops within a rank); the
         ! Phase 2 → Phase 3 ordering invariant the pre-α three-phase
         ! split mitigated is re-established.
         !
         ! α seams keep `CADENCE_END_OF_STEP` and skip this loop; their
         ! peer ghosts continue to hold the previous end-of-step exchange
         ! (Berger-Oliger 1984 / AMReX FillCoarsePatch).
         !
         ! `associate` wrapper: nvfortran 26.1 workaround for the polymorphic
         ! array element dispatch bug (same fix shape as Phase 3 below;
         ! commit 0062a237). Pre-emptive: this is a new dispatch site that
         ! the workaround should cover from the start.
         ! The slots come from the manifest (the same on every rank, whatever rows a rank holds): `seam_fill` moves the
         ! local and the cross-rank rows (issue #40), and every rank must reach the same (realm, slot) pairs in order.
         do is = 1_I4P, int(size(realm), I4P)
            if (.not. allocated(realm(is)%adam%maps%seam_local_peer_realm)) cycle
            do p = 1_I4P, int(size(realm(is)%adam%maps%seam_local_peer_realm), I4P)
               if (realm(is)%adam%maps%seam_local_cadence(p) /= CADENCE_STAGE_COINCIDENT) cycle
               associate(r=>realm(is)) ! nvfortran 26.1 polymorphic array element dispatch workaround (0062a237)
                  call seam_fill(self=r, realm=realm, p=p)
               endassociate
            enddo
         enddo

         ! Phase 3 — residuals (with flux_register) + assign on each participating realm.
         do is = 1_I4P, int(size(realm), I4P)
            if (k > K_realm(is)) cycle
            ! NVFORTRAN 26.1 WORKAROUND: dispatching the polymorphic array
            ! element `realm(is)%end_stage_forest(...)` directly segfaults
            ! inside libnvf's `pgf90_copy_f90_argl_i8` when marshalling the
            ! explicit-bound `q_gpu(1:, 1-self%ngc:, ...)` actual passed to
            ! `compute_residuals_dev` deep in the body (rmf-2realm/fnl
            ! reproducer; line 1895 of adam_prism_fnl_object.F90). Binding
            ! `realm(is)` to a polymorphic scalar via `associate` before
            ! dispatch lets the marshaller resolve the element's concrete-
            ! type stride correctly. Other dispatch sites in this routine
            ! happen not to trigger the bug because their bodies do not
            ! reach an equally complex arg-marshalling path; if a future
            ! refactor exposes them, apply the same wrapper. CPU build is
            ! unaffected (gfortran does not exhibit the issue).
            associate(r => realm(is))
               call r%end_stage_forest(k=k, K_total=K_max, dt=dt, flux_register=self%flux_register)
            end associate
         enddo

      enddo
      ! Cross-rank reduce of the fine-side accumulators (after the final stage's
      ! accumulation, which happened inside compute_residuals at k = nrk).
      call self%flux_register%reduce_fine_sums
      do is = 1_I4P, int(size(realm), I4P)
         call realm(is)%close_step_forest(dt=dt)
      enddo
      ! α.r1 reflux: applied AFTER close_step_forest's update_q has committed the
      ! step's q. The realm-side body writes the end-of-step Berger-Colella
      ! correction directly into self%q (full dt/dx weight, no stage RK
      ! coefficient), gated to each realm's final stage. Running it post-update_q
      ! is required: the correction is to the committed solution, not to a stage
      ! residual buffer (the pre-update_q q_rk path silently no-op'd for SSP and
      ! entangled the stage beta weight — see apply_reflux_to_stage_forest).
      call self%apply_reflux_corrections(realm=realm, dt=dt)

      ! Phase 5 — end-of-step inter-realm seam fill (α coherence barrier).
      !
      ! After close_step_forest, every realm has stage_active == 0 and its
      ! committed `q` is the peer-visible state. fill_seam_from_peer_forest's
      ! buffer-selection logic reads peer%q when peer%stage_active == 0, so
      ! no new TBP or call-site contract is needed. Same dispatch as the
      ! former per-stage Phase 2: receiver walks its own seam map, copies
      ! peer-INTERIOR cells into self's GHOST cells; backend dispatch via
      ! `select type(peer)` inside the receiver's override.
      !
      ! Per-seam gating: only seams with `coupling_cadence ==
      ! CADENCE_END_OF_STEP` (the α default) fire here. β seams
      ! (`CADENCE_STAGE_COINCIDENT`) were already filled at every substage
      ! in Phase 2 inside the K loop; firing again here would double-write
      ! and waste work.
      !
      ! `associate` wrapper: same nvfortran 26.1 workaround as Phase 2 and
      ! Phase 3 (commit 0062a237). Pre-emptive coverage.
      do is = 1_I4P, int(size(realm), I4P)
         if (.not. allocated(realm(is)%adam%maps%seam_local_peer_realm)) cycle
         do p = 1_I4P, int(size(realm(is)%adam%maps%seam_local_peer_realm), I4P)
            if (realm(is)%adam%maps%seam_local_cadence(p) /= CADENCE_END_OF_STEP) cycle
            associate(r=>realm(is)) ! nvfortran 26.1 polymorphic array element dispatch workaround (0062a237)
               call seam_fill(self=r, realm=realm, p=p)
            endassociate
         enddo
      enddo

      deallocate(K_realm)
   endif
   endsubroutine evolve_one_step

   subroutine is_done(self, realm, done)
   !< Decide whether the whole forest has finished evolving.
   !<
   !< Each realm reports its local predicate via `is_done_forest`; the forest AND-reduces across all realms (intra-process)
   !< and then across all MPI ranks (`MPI_ALLREDUCE` on `MPI_COMM_WORLD`). AND-reduction means the forest keeps evolving as
   !< long as ANY realm wants to — matching the legacy single-realm semantics for v1 (with one realm the global predicate
   !< equals that realm's local one).
   class(forest_object), intent(in)    :: self       !< The forest.
   class(realm_object),  intent(inout) :: realm(:)   !< The realms to query.
   logical,              intent(out)   :: done       !< Forest-global termination predicate.
   logical                             :: done_local !< Per-realm local predicate.
   integer(I4P)                        :: is         !< Realm index.
   integer(I4P)                        :: ierr       !< MPI error code.

   done = .true.
   do is = 1, int(size(realm), I4P)
      call realm(is)%is_done_forest(done=done_local)
      done = done .and. done_local
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, done, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
   endsubroutine is_done

   subroutine post_step(self, realm)
   !< Run every realm's post-step diagnostics / IO / AMR block.
   class(forest_object), intent(in)    :: self     !< The forest.
   class(realm_object),  intent(inout) :: realm(:) !< The realms to query.
   integer(I4P)                        :: is       !< Realm index.

   do is = 1, int(size(realm), I4P)
      call realm(is)%post_step_forest(dt=0._R8P, t=0._R8P, it=0_I4P, realm=realm)
   enddo
   endsubroutine post_step

   subroutine simulate(self, realm, filename)
   !< Drive the full simulation: initialize, time-loop, finalize.
   !<
   !< Top-level entry point the program driver calls. The time loop is: initialize → loop {evolve_one_step → post_step
   !< → is_done} → finalize. Each step invokes the orchestrator-contract TBPs on every realm; the per-realm body decides what
   !< app-specific work to do internally.
   class(forest_object), intent(inout) :: self     !< The forest.
   class(realm_object),  intent(inout) :: realm(:) !< The realms to evolve.
   character(*),         intent(in)    :: filename !< Input parameters file name.
   logical                             :: done     !< Forest-global termination predicate.

   call self%initialize(realm, filename=filename)
   done = .false.
   do
      call self%evolve_one_step(realm=realm)
      call self%post_step(realm=realm)
      call self%is_done(realm=realm, done=done)
      if (done) exit
   enddo
   call self%finalize(realm=realm)
   endsubroutine simulate

   subroutine simulate_from_manifest(self, realm, manifest)
   !< Drive the full simulation using per-realm INIs from a manifest.
   !<
   !< Like `simulate` but uses `initialize_from_manifest` to populate each realm from its own INI file, and (via that
   !< initialize) wires the inter-realm topology from the manifest. The time-loop body is identical to `simulate`.
   class(forest_object),     intent(inout) :: self     !< The forest.
   class(realm_object),      intent(inout) :: realm(:) !< The realms to evolve.
   type(forest_manifest_t),  intent(in)    :: manifest !< Parsed manifest.
   logical                                 :: done     !< Forest-global termination predicate.

   call self%initialize_from_manifest(realm=realm, manifest=manifest)
   done = .false.
   do
      call self%evolve_one_step(realm=realm)
      call self%post_step(realm=realm)
      call self%is_done(realm=realm, done=done)
      if (done) exit
   enddo
   call self%finalize(realm=realm)
   endsubroutine simulate_from_manifest

   ! private methods
   subroutine check_beta_admissibility(realm, manifest)
   !< Enforce β admissibility contract on every `stage_coincident` seam.
   !<
   !< For each manifest face-pair with `coupling_cadence == CADENCE_STAGE_COINCIDENT`,
   !< query both endpoint realms' `coupling_descriptor_forest` and verify
   !< they agree on (scheme_time, rk_scheme, nv). On any mismatch
   !< `error_stop` with a precise diagnostic naming the offending pair and
   !< the specific descriptor field that disagrees. Per-realm K participation
   !< (`stages_per_step_forest()`) is verified too: under β both endpoints
   !< must report the same K (asymmetric-K is α's domain, not β's).
   !<
   !< Runs once at `initialize_from_manifest` time, after every realm's
   !< `initialize_forest` has populated `numerics`, `rk`, `physics`.
   class(realm_object),     intent(inout) :: realm(:) !< The initialized realms.
   type(forest_manifest_t), intent(in)    :: manifest !< Parsed manifest.
   integer(I4P)                           :: f        !< Face-pair index.
   integer(I4P)                           :: ra, rb   !< Endpoint realm indices.
   character(:), allocatable              :: st_a, rk_a, st_b, rk_b !< Descriptor strings.
   integer(I4P)                           :: nv_a, nv_b, k_a, k_b   !< Descriptor integers.

   if (.not. allocated(manifest%face_pairs)) return
   do f = 1_I4P, int(size(manifest%face_pairs), I4P)
      if (manifest%face_pairs(f)%coupling_cadence /= CADENCE_STAGE_COINCIDENT) cycle
      ra = manifest%face_pairs(f)%realm_a
      rb = manifest%face_pairs(f)%realm_b
      call realm(ra)%coupling_descriptor_forest(scheme_time=st_a, rk_scheme=rk_a, nv=nv_a)
      call realm(rb)%coupling_descriptor_forest(scheme_time=st_b, rk_scheme=rk_b, nv=nv_b)
      k_a = realm(ra)%stages_per_step_forest()
      k_b = realm(rb)%stages_per_step_forest()
      if (nv_a < 0_I4P .or. nv_b < 0_I4P) &
         call mpih%error_stop(msg='forest_object%check_beta_admissibility: face_pair '//trim(str(f, .true.))// &
            ' has stage_coincident cadence but realm '//trim(str(ra, .true.))//' or '//                       &
            trim(str(rb, .true.))//' does not implement coupling_descriptor_forest')
      if (st_a /= st_b) &
         call mpih%error_stop(msg='forest_object%check_beta_admissibility: face_pair '//trim(str(f, .true.))// &
            ' stage_coincident requires equal scheme_time between realm '//trim(str(ra, .true.))//' ("'//     &
            st_a//'") and realm '//trim(str(rb, .true.))//' ("'//st_b//'")')
      if (rk_a /= rk_b) &
         call mpih%error_stop(msg='forest_object%check_beta_admissibility: face_pair '//trim(str(f, .true.))// &
            ' stage_coincident requires equal rk_scheme between realm '//trim(str(ra, .true.))//' ("'//       &
            rk_a//'") and realm '//trim(str(rb, .true.))//' ("'//rk_b//'")')
      if (nv_a /= nv_b) &
         call mpih%error_stop(msg='forest_object%check_beta_admissibility: face_pair '//trim(str(f, .true.))// &
            ' stage_coincident requires equal physics nv between realm '//trim(str(ra, .true.))//' ('//       &
            trim(str(nv_a, .true.))//') and realm '//trim(str(rb, .true.))//' ('//trim(str(nv_b, .true.))//')')
      if (k_a /= k_b) &
         call mpih%error_stop(msg='forest_object%check_beta_admissibility: face_pair '//trim(str(f, .true.))// &
            ' stage_coincident requires equal K (no stalling) between realm '//trim(str(ra, .true.))//' ('// &
            trim(str(k_a, .true.))//') and realm '//trim(str(rb, .true.))//' ('//trim(str(k_b, .true.))//')')
   enddo
   endsubroutine check_beta_admissibility

   subroutine populate_inter_realm_topology(self, realm, manifest)
   !< Translate manifest face-pairs into per-realm maps%inter_realm_neighbors.
   !<
   !< Each manifest face-pair becomes TWO entries — one in realm_a's
   !< `adam%maps%inter_realm_neighbors` (looking outward toward realm_b)
   !< and one in realm_b's (looking back). Both entries carry the SAME
   !< coupling kind. Block indices are NOT resolved at this stage; the
   !< manifest declares realm-level coupling (which realm-face touches
   !< which realm-face), and the realm-side override of
   !< `exchange_inter_realm_halos_forest` is responsible for enumerating
   !< per-block face cells at exchange time. This keeps the manifest small
   !< and avoids encoding block layouts that depend on AMR / decomposition
   !< state not known at INI parse time.
   class(forest_object),     intent(inout) :: self                !< The forest.
   class(realm_object),      intent(inout) :: realm(:)            !< Initialized realms whose adam%maps gets populated.
   type(forest_manifest_t),  intent(in)    :: manifest            !< Parsed manifest.
   integer(I4P), allocatable               :: per_realm_count(:)  !< How many neighbour entries each realm gets.
   integer(I4P), allocatable               :: per_realm_cursor(:) !< Write cursor per realm.
   integer(I4P)                            :: f, is               !< Face-pair and realm index counters.
   integer(I4P)                            :: n_intra             !< Intra-realm AMR faces in the flux register.
   type(forest_face_pair_t)                :: pair                !< Loop alias.

   if (.not. allocated(manifest%face_pairs)) then
      ! No inter-realm topology: the flux register holds the intra-realm AMR faces only.
      call self%register_intra_realm_amr_seams(realm=realm)
      return
   endif

   ! Pass 1: count entries per realm.
   allocate(per_realm_count(self%n))
   per_realm_count = 0_I4P
   do f = 1_I4P, int(size(manifest%face_pairs), I4P)
      pair = manifest%face_pairs(f)
      if (pair%realm_a < 1_I4P .or. pair%realm_a > self%n) &
         call mpih%error_stop(msg='forest_object%populate_inter_realm_topology: face pair realm_a out of range')
      if (pair%realm_b < 1_I4P .or. pair%realm_b > self%n) &
         call mpih%error_stop(msg='forest_object%populate_inter_realm_topology: face pair realm_b out of range')
      per_realm_count(pair%realm_a) = per_realm_count(pair%realm_a) + 1_I4P
      per_realm_count(pair%realm_b) = per_realm_count(pair%realm_b) + 1_I4P
   enddo
   ! Validate every coupling against the realm geometry before any seam data is built (issue #52).
   call check_couplings(realm=realm, manifest=manifest)
   do is = 1_I4P, self%n
      if (allocated(realm(is)%adam%maps%inter_realm_neighbors)) deallocate(realm(is)%adam%maps%inter_realm_neighbors)
      if (per_realm_count(is) > 0_I4P) allocate(realm(is)%adam%maps%inter_realm_neighbors(per_realm_count(is)))
   enddo
   ! Pass 2: write entries into each realm's array.
   allocate(per_realm_cursor(self%n))
   per_realm_cursor = 0_I4P
   do f = 1_I4P, int(size(manifest%face_pairs), I4P)
      pair = manifest%face_pairs(f)
      ! entry on realm_a's array: my=a, peer=b
      per_realm_cursor(pair%realm_a) = per_realm_cursor(pair%realm_a) + 1_I4P
      call set_neighbor(slot=realm(pair%realm_a)%adam%maps%inter_realm_neighbors(per_realm_cursor(pair%realm_a)), &
                        my_realm=pair%realm_a,my_face=pair%face_a,peer_realm=pair%realm_b,peer_face=pair%face_b,  &
                        coupling=pair%coupling)
      ! entry on realm_b's array: my=b, peer=a
      per_realm_cursor(pair%realm_b) = per_realm_cursor(pair%realm_b) + 1_I4P
      call set_neighbor(slot=realm(pair%realm_b)%adam%maps%inter_realm_neighbors(per_realm_cursor(pair%realm_b)), &
                        my_realm=pair%realm_b,my_face=pair%face_b,peer_realm=pair%realm_a,peer_face=pair%face_a,  &
                        coupling=pair%coupling)
   enddo
   ! Register inter-realm seams with the program-scope flux register
   ! For each (face-pair, block-on-coarse-side) tuple, one entry is added to
   ! the register. The "coarse" / "fine" labels follow the manifest's a/b
   ! ordering; for the current same-resolution (COUPLING_MIRROR) case the
   ! labels are conventional and the accumulator values are nominally equal
   ! on both sides — the reflux correction will be round-off zero in
   ! expectation. The structural cost (allocated registers, populated
   ! topology) is the same as for the true coarse-fine AMR case that will
   ! exercise these accumulators non-trivially in follow-up commits.
   ! Flux register: the intra-realm AMR faces first (cursors 1..n_intra), then the inter-realm seam faces after them
   ! (issue #37); both lists are walked on the replicated trees, the same on every rank (issue #40).
   call self%register_intra_realm_amr_seams(realm=realm, extra_faces=count_inter_realm_seam_faces(realm, manifest), &
                                            nfaces_intra=n_intra)
   call register_inter_realm_seams(realm=realm, manifest=manifest, flux_register=self%flux_register, first_cursor=n_intra)
   ! Build the seam ghost rows of every realm (issue #40): from the replicated trees, every rank enumerates every seam
   ! ghost of every realm in one canonical order and keeps the rows it takes part in: local rows (ghost and peer cell on
   ! this rank, `seam_local_map_ghost_cell`), receive rows (ghost here, cell on another rank, `seam_mpi_recv_cell`) and
   ! send rows (cell here, ghost on another rank, `seam_mpi_send_cell`). The peer slots and their cadence come from the
   ! manifest, so they are the same on every rank. `adam_seam_exchange` consumes them.
   call build_seam_rows(realm=realm, manifest=manifest)
   ! Override the BC crown's bc_type column to BC_SEAM for entries that
   ! lie on an inter-realm seam face.
   !
   ! Why this is needed: each realm parses its INI in isolation and
   ! declares physical BCs on all 6 faces (bc_x_max, bc_x_min, ...). For
   ! a realm whose +x face is glued to another realm's -x face by the
   ! manifest, the INI's bc_x_max = "Neumann" declaration is wrong —
   ! that face is a SEAM, not a physical boundary. PRISM's
   ! `make_local_maps_bc` (which runs during each realm's
   ! `initialize_forest`, well before this point) has already populated
   ! `local_map_bc_crown` with BC_NEUMANN entries for the seam face's
   ! cells. Without this override, `set_boundary_conditions` would then
   ! extrapolate Neumann values into those cells at every stage,
   ! overwriting the peer-interior values written by
   ! `exchange_inter_realm_halos_forest`.
   !
   ! Mechanism: walk each realm's BC crown post-hoc, find entries whose
   ! block-face matches a manifest-declared seam, and flip column 8
   ! (`bc_type`) from BC_NEUMANN/EXTRAPOLATION/etc. to BC_SEAM.
   ! `set_boundary_conditions` has no dispatch branch for BC_SEAM →
   ! those entries are silently no-oped → the seam-exchange-written
   ! ghosts survive.
   !
   ! The manifest is the authoritative source of truth about realm
   ! topology; this override applies that authority over the realm's
   ! own INI declarations at the right semantic layer.
   call override_seam_bc_in_crown(realm=realm, manifest=manifest)
   ! Backend hook: each realm propagates the freshly-built host seam maps
   ! to whatever device-side / backend-specific structures it owns. CPU
   ! realms no-op; FNL realms refresh maps_fnl%seam_local_* device pointers.
   block
      integer(I4P) :: is_tb
      do is_tb = 1_I4P, int(size(realm), I4P)
         call realm(is_tb)%after_topology_build_forest
      enddo
   endblock
   contains
      subroutine set_neighbor(slot, my_realm, my_face, peer_realm, peer_face, coupling)
      !< Set inter-realm neighbor.
      type(inter_realm_neighbor_t), intent(out) :: slot       !< Inter-realm neighbor slot.
      integer(I4P),                 intent(in)  :: my_realm   !< My realm.
      integer(I4P),                 intent(in)  :: my_face    !< My face.
      integer(I4P),                 intent(in)  :: peer_realm !< Peer realm.
      integer(I4P),                 intent(in)  :: peer_face  !< Peer face
      integer(I4P),                 intent(in)  :: coupling   !< Coupling type.

      slot%my_realm   = my_realm
      slot%my_block   = 0_I4P
      slot%my_face    = my_face
      slot%peer_realm = peer_realm
      slot%peer_block = 0_I4P
      slot%peer_face  = peer_face
      slot%coupling   = coupling
      endsubroutine set_neighbor

      function count_inter_realm_seam_faces(realm, manifest) result(nfaces)
      !< Count the inter-realm seam register faces: one per (face-pair, leaf of the coarse side on its seam face; realm_a
      !< for a mirror seam), every leaf of the replicated tree whatever rank owns it (the same count on every rank), the
      !< pass 1 of `register_inter_realm_seams`.
      class(realm_object),     intent(inout) :: realm(:)     !< Initialized realms.
      type(forest_manifest_t), intent(in)    :: manifest     !< Parsed manifest.
      integer(I4P)                           :: nfaces       !< Inter-realm seam faces.
      integer(I4P)                           :: f            !< Face-pair counter.
      integer(I4P)                           :: rc, fc       !< Coarse-side realm and face.
      integer(I4P)                           :: rf, ff       !< Fine-side realm and face.
      integer(I4P)                           :: axis, sgn    !< Coarse face axis and side.
      logical                                :: refined      !< 2:1 seam.
      type(tree_iterator_object)             :: iter         !< Tree traversal cursor.
      type(tree_node_object), pointer        :: node_ptr     !< Current leaf.

      nfaces = 0_I4P
      if (.not. allocated(manifest%face_pairs)) return
      do f = 1_I4P, int(size(manifest%face_pairs), I4P)
         call seam_sides(realm, manifest%face_pairs(f), rc, fc, rf, ff, refined)
         call face_axis_sign(fc, axis, sgn)
         iter%b = 1_I4P ; iter%p => null()
         do while (realm(rc)%adam%tree%loop(iter, node_ptr=node_ptr))
            if (leaf_on_realm_face(realm(rc), node_ptr%code, axis, sgn)) nfaces = nfaces + 1_I4P
         enddo
      enddo
      endfunction count_inter_realm_seam_faces

      subroutine register_inter_realm_seams(realm, manifest, flux_register, first_cursor)
      !< Populate `flux_register` from the manifest face-pairs.
      !<
      !< One register face per (face-pair, leaf of the coarse side on its seam face), enumerated on the replicated tree in
      !< its own order: every rank registers the same faces under the same cursors, whichever rank owns the leaves (issue
      !< #40), like the intra-realm AMR faces registered before them, so `reduce_fine_sums` completes them with the same
      !< collective sequence on every rank. `coarse_rank`/`coarse_block` are the owner and owner-local block of the coarse
      !< leaf; `nface_cells` comes from the coarse realm's grid, `nv` and the per-stage depth (`stages_per_step_forest()`)
      !< too.
      !<
      !< * `mirror` (`SEAM_KIND_INTER_REALM`, realm_a coarse by convention): the realm_b leaf containing the first cell
      !<   centre across the face must meet the realm_a leaf face to face (seams between blocks that do not line up: #51);
      !<   its skin covers the coarse face 1:1.
      !< * `refined` (`SEAM_KIND_INTER_REALM_REFINED`, issue #52): the coarse side is the realm with the larger cells; the
      !<   four fine leaves covering a coarse face are found at its quadrant centres, and each fine block records its 2:1
      !<   quadrant offsets in `maps%amr_seam_quadrant`, so the apps restrict its skin into that quadrant exactly as for an
      !<   intra-realm 2:1 face.
      !<
      !< The signed lookup `inter_realm_face_register_index(block, bc_fec)` (+cursor on the coarse block, -cursor on the
      !< fine blocks) and the quadrant table are written only for the blocks this rank owns: `block_index` is an owner-local
      !< slot (issue #28 D1).
      class(realm_object),        intent(inout) :: realm(:)       !< Initialized realms.
      type(forest_manifest_t),    intent(in)    :: manifest       !< Parsed manifest.
      type(flux_register_object), intent(inout) :: flux_register  !< Berger-Colella reflux accumulator owned by the forest.
      integer(I4P), optional,     intent(in)    :: first_cursor   !< Append after this cursor to a register (and index)
                                                                  !< already initialized by the intra-realm AMR pass.
      integer(I4P)                              :: f              !< Face-pair counter.
      integer(I4P)                              :: rc, fc         !< Coarse-side realm and face.
      integer(I4P)                              :: rf, ff         !< Fine-side realm and face.
      logical                                   :: refined        !< 2:1 seam.
      integer(I4P)                              :: c_axis, c_sign !< Coarse face axis and side.
      integer(I4P)                              :: t1, t2         !< Tangential axes (inner, outer: the skin order).
      integer(I4P)                              :: qi, qj         !< Quadrant offsets.
      integer(I4P)                              :: nfaces_total   !< Total register entries.
      integer(I4P)                              :: cursor         !< Write cursor into the register.
      integer(I4P)                              :: nface_cells    !< Cell count on the coarse-face skin.
      integer(I4P)                              :: bc_fec_c       !< BC fec of the coarse face.
      integer(I4P)                              :: bc_fec_f       !< BC fec of the fine face.
      integer(I8P)                              :: code_f         !< Morton code of a fine-side leaf.
      real(R8P)                                 :: emin_c(3)      !< Coarse leaf origin.
      real(R8P)                                 :: emax_c(3)      !< Coarse leaf upper corner.
      real(R8P)                                 :: d_c(3)         !< Coarse cell size.
      real(R8P)                                 :: emin_f(3)      !< Fine-side leaf origin.
      real(R8P)                                 :: emax_f(3)      !< Fine-side leaf upper corner.
      real(R8P)                                 :: d_f(3)         !< Fine-side cell size.
      real(R8P)                                 :: xq(3)          !< Lookup point across the face.
      type(tree_iterator_object)                :: iter           !< Tree traversal cursor.
      type(tree_node_object), pointer           :: node_c         !< Coarse-side leaf.
      type(tree_node_object), pointer           :: node_f         !< Fine-side leaf.

      if (.not. allocated(manifest%face_pairs)) then
         ! No inter-realm topology — initialize with zero faces so the
         ! register's `is_initialized_` flag flips and the per-step `reset`
         ! call becomes a safe no-op on the empty register.
         if (.not. present(first_cursor)) call flux_register%initialize(nfaces=0_I4P)
         return
      endif

      ! Pass 1: count register entries.
      nfaces_total = count_inter_realm_seam_faces(realm, manifest)
      if (.not. present(first_cursor)) call flux_register%initialize(nfaces=nfaces_total)
      if (nfaces_total == 0_I4P) return

      ! Allocate the per-realm (block, face_1_6) → register_index lookup.
      ! Sized (nb, 6); zero means "not a seam face", positive = 1-based
      ! index into flux_register%face(:). Consumed by the FV residual
      ! routines to know where to accumulate fluxes.
      ! When appending (`first_cursor` present) the lookup already holds the intra-realm AMR faces: keep it.
      do is = 1_I4P, int(size(realm), I4P)
         if (present(first_cursor)) exit
         block
            integer(I4P) :: nb_realm
            nb_realm = int(realm(is)%adam%field%blocks_number, I4P)
            if (allocated(realm(is)%adam%maps%inter_realm_face_register_index)) &
               deallocate(realm(is)%adam%maps%inter_realm_face_register_index)
            if (nb_realm > 0_I4P) then
               allocate(realm(is)%adam%maps%inter_realm_face_register_index(1:nb_realm, 1:6))
               realm(is)%adam%maps%inter_realm_face_register_index = 0_I4P
            endif
         endblock
      enddo

      ! Pass 2: register one entry per (face-pair, coarse leaf on its seam face) and fill the signed lookups.
      ! Sign convention (consumed by the FV reflux hooks):
      !   +cursor stored on the coarse-side realm
      !   -cursor stored on the fine-side   realm
      !       0   stored anywhere = "not a seam face".
      cursor = 0_I4P
      if (present(first_cursor)) cursor = first_cursor
      do f = 1_I4P, int(size(manifest%face_pairs), I4P)
         call seam_sides(realm, manifest%face_pairs(f), rc, fc, rf, ff, refined)
         call face_axis_sign(fc, c_axis, c_sign)
         t1 = merge(2_I4P, 1_I4P, c_axis == 1_I4P)
         t2 = merge(2_I4P, 3_I4P, c_axis == 3_I4P)
         nface_cells = tangential_cell_count(realm(rc), c_axis)
         bc_fec_c = face_code_to_bc_fec(fc)
         bc_fec_f = face_code_to_bc_fec(ff)
         iter%b = 1_I4P ; iter%p => null()
         do while (realm(rc)%adam%tree%loop(iter, node_ptr=node_c))
            if (.not. leaf_on_realm_face(realm(rc), node_c%code, c_axis, c_sign)) cycle
            cursor = cursor + 1_I4P
            ! `fine_block` deliberately omitted: nothing reads it, and passing an empty array literal
            ! `[integer(I4P) ::]` poisons the descriptor copy on nvfortran 26.x.
            call flux_register%register_face(face_index=cursor,                                                      &
                                             seam_kind=merge(SEAM_KIND_INTER_REALM_REFINED, SEAM_KIND_INTER_REALM,   &
                                                             refined),                                               &
                                             coarse_realm=rc,                                                        &
                                             coarse_rank=node_c%myrank,                                              &
                                             coarse_block=int(node_c%block_index, I4P),                              &
                                             coarse_face=fc,                                                         &
                                             fine_realm=rf,                                                          &
                                             nface_cells=nface_cells,                                                &
                                             nv=int(realm(rc)%adam%field%nv, I4P),                                   &
                                             n_stages=realm(rc)%stages_per_step_forest())
            if (node_c%myrank == mpih%myrank .and. bc_fec_c > 0_I4P .and. bc_fec_c <= 6_I4P) &
               realm(rc)%adam%maps%inter_realm_face_register_index(int(node_c%block_index, I4P), bc_fec_c) = +cursor
            call leaf_metrics(realm(rc), node_c%code, emin_c, d_c, emax_c)
            do qj = 0_I4P, merge(1_I4P, 0_I4P, refined)
               do qi = 0_I4P, merge(1_I4P, 0_I4P, refined)
                  ! lookup point: the fine-side cell centre across the face (in the quadrant (qi, qj) when refined)
                  xq = 0.5_R8P * (emin_c + emax_c)
                  if (refined) then
                     xq(t1) = emin_c(t1) + (real(qi, R8P) + 0.5_R8P) * 0.5_R8P * (emax_c(t1) - emin_c(t1))
                     xq(t2) = emin_c(t2) + (real(qj, R8P) + 0.5_R8P) * 0.5_R8P * (emax_c(t2) - emin_c(t2))
                     xq(c_axis) = merge(emax_c(c_axis) + 0.25_R8P * d_c(c_axis), emin_c(c_axis) - 0.25_R8P * d_c(c_axis), &
                                        c_sign > 0_I4P)
                  else
                     xq(c_axis) = merge(emax_c(c_axis) + 0.5_R8P * d_c(c_axis), emin_c(c_axis) - 0.5_R8P * d_c(c_axis), &
                                        c_sign > 0_I4P)
                  endif
                  code_f = realm(rf)%adam%tree%get_closest_block(grid=realm(rf)%adam%grid, point=xq)
                  node_f => realm(rf)%adam%tree%node(code=code_f)
                  call leaf_metrics(realm(rf), code_f, emin_f, d_f, emax_f)
                  if (.not. refined) then
                     if (.not. meet_face_to_face(emin_c, emax_c, emin_f, emax_f, c_axis, c_sign))                       &
                        call mpih%error_stop(msg='forest_object%populate_inter_realm_topology: the seam block '//        &
                           trim(str(int(node_c%block_index, I4P), .true.))//' of realm '//trim(str(rc, .true.))//        &
                           ' (rank '//trim(str(node_c%myrank, .true.))//') does not meet a block of realm '//             &
                           trim(str(rf, .true.))//' face to face (face_pair '//trim(str(f, .true.))//'): seams '//      &
                           'between blocks that do not line up are not implemented yet (issue #51)')
                  endif
                  if (node_f%myrank /= mpih%myrank .or. bc_fec_f < 1_I4P .or. bc_fec_f > 6_I4P) cycle
                  if (allocated(realm(rf)%adam%maps%inter_realm_face_register_index))                                    &
                     realm(rf)%adam%maps%inter_realm_face_register_index(int(node_f%block_index, I4P), bc_fec_f) = -cursor
                  ! the fine block's 2:1 quadrant within the coarse skin (inner, outer offsets), as for intra-realm faces
                  if (refined .and. allocated(realm(rf)%adam%maps%amr_seam_quadrant))                                    &
                     realm(rf)%adam%maps%amr_seam_quadrant(1:2, int(node_f%block_index, I4P), bc_fec_f) = [qi, qj]
               enddo
            enddo
         enddo
      enddo
      endsubroutine register_inter_realm_seams

      function block_face_on_realm_boundary(this_realm, b, axis, sgn) result(yes)
      !< Return .true. iff block `b`'s face on (axis, sgn) lies on the realm boundary.
      !<
      !< Geometric test against `grid%domain_emin/emax` and
      !< `field%emin/emax`, with a small absolute tolerance. This mirrors
      !< the identically-named helper inside the PRISM-CPU realm; a future
      !< refactor should lift the canonical version into `adam_maps_object`.
      class(realm_object), intent(in) :: this_realm   !< Realm to query.
      integer(I4P),        intent(in) :: b            !< Block index.
      integer(I4P),        intent(in) :: axis         !< 1=x, 2=y, 3=z.
      integer(I4P),        intent(in) :: sgn          !< +1 if checking MAX face, -1 if MIN.
      logical                         :: yes          !< Test result.
      real(R8P)                       :: face_coord   !< Face coordinate.
      real(R8P)                       :: target_coord !< Taget coordinate.
      real(R8P)                       :: tol          !< Tolerance.

      if (sgn > 0_I4P) then
         face_coord   = this_realm%adam%field%emax(axis, b)
         target_coord = this_realm%adam%grid%domain_emax(axis)
      else
         face_coord   = this_realm%adam%field%emin(axis, b)
         target_coord = this_realm%adam%grid%domain_emin(axis)
      endif
      tol = max(abs(target_coord), 1._R8P) * 1.0e-10_R8P
      yes = abs(face_coord - target_coord) <= tol
      endfunction block_face_on_realm_boundary

      pure function tangential_cell_count(this_realm, axis) result(n)
      !< Return the product of the two cell counts tangential to `axis`.
      !<
      !< For axis=1 (x-normal face): n = nj * nk.
      !< For axis=2 (y-normal face): n = ni * nk.
      !< For axis=3 (z-normal face): n = ni * nj.
      !<
      !< This is the cell count for a single block's face skin (NOT the
      !< whole-realm face skin); the realm-level face skin is the sum over
      !< the seam blocks, each contributing this count.
      class(realm_object), intent(in) :: this_realm !< Realm to query.
      integer(I4P),        intent(in) :: axis       !< 1=x, 2=y, 3=z.
      integer(I4P)                    :: n          !< Counter.

      associate(g => this_realm%adam%grid)
         select case (axis)
         case (1_I4P); n = g%nj * g%nk
         case (2_I4P); n = g%ni * g%nk
         case (3_I4P); n = g%ni * g%nj
         case default; n = 0_I4P
         endselect
      endassociate
      endfunction tangential_cell_count



      subroutine override_seam_bc_in_crown(realm, manifest)
      !< Overwrite the bc_type column of `local_map_bc_crown` to BC_SEAM
      !< for every entry whose (block, face) pair lies on a manifest-
      !< declared inter-realm seam.
      !<
      !< The walk uses `face_code_to_bc_fec` (FACE_X_MAX/MIN/... → BC
      !< face fec 1..6 in the FEC_TO_DELTA / FEC_1_6_ARRAY space) and
      !< `block_face_on_realm_boundary` (geometric test against the
      !< realm's `domain_emin/emax`) — the two helpers already defined
      !< below as siblings.
      !<
      !< The crown layout: `local_map_bc_crown(c, 1..9, crown)` =
      !<   [b, i, j, k, idelta, jdelta, kdelta, bc_type, fec]
      !<
      !< For each crown entry whose:
      !<   * block `b` lies on the realm's geometric seam face
      !<     (per `block_face_on_realm_boundary`), AND
      !<   * `fec_1_6_array(fec)` matches the seam face's BC fec code
      !<     (per `face_code_to_bc_fec`),
      !< the entry's `bc_type` is rewritten to `BC_SEAM`. The
      !< `set_boundary_conditions` dispatch ladder in each app's realm
      !< extension (currently `prism_cpu_object`) has no branch for
      !< `BC_SEAM`, so those entries become silent no-ops at consumption
      !< time. The cells they target retain the ghost values written by
      !< `exchange_inter_realm_halos_forest`.
      !<
      !< Edge and corner entries that share an axis with the seam face
      !< are ALSO overridden (because `fec_1_6_array(fec)` maps them to
      !< the seam-face code). Those entries don't fire any Neumann/etc.
      !< write in the current PRISM ladder (the dispatch is `case(fec)
      !< 1..6`, not 7..26), so the override is a no-op for them today —
      !< but it future-proofs the contract: any future BC kind that
      !< writes edge/corner cells will respect the seam override
      !< automatically.
      class(realm_object),     intent(inout) :: realm(:)
      type(forest_manifest_t), intent(in)    :: manifest
      integer(I4P)                           :: f
      type(forest_face_pair_t)               :: pair

      if (.not. allocated(manifest%face_pairs)) return
      do f = 1_I4P, int(size(manifest%face_pairs), I4P)
         pair = manifest%face_pairs(f)
         call mark_seam_in_crown(realm, my_realm_idx=pair%realm_a, my_face=pair%face_a)
         call mark_seam_in_crown(realm, my_realm_idx=pair%realm_b, my_face=pair%face_b)
      enddo
      endsubroutine override_seam_bc_in_crown

      subroutine mark_seam_in_crown(realm, my_realm_idx, my_face)
      !< Helper for `override_seam_bc_in_crown`: walk one realm's BC
      !< crown and rewrite bc_type → BC_SEAM for the entries lying on
      !< the (my_face) seam.
      class(realm_object), intent(inout) :: realm(:)
      integer(I4P),        intent(in)    :: my_realm_idx, my_face
      integer(I4P)                       :: axis, sgn, bc_fec_seam
      integer(I4P)                       :: ngc, c, crown
      integer(I8P)                       :: b_i8, fec_i8

      call face_axis_sign(my_face, axis, sgn)
      bc_fec_seam = face_code_to_bc_fec(my_face)
      if (bc_fec_seam == 0_I4P) return  ! malformed face code; defensive.
      if (.not. allocated(realm(my_realm_idx)%adam%maps%local_map_bc_crown)) return
      ngc = realm(my_realm_idx)%adam%grid%ngc
      associate(crown_map => realm(my_realm_idx)%adam%maps%local_map_bc_crown)
      do crown = 1_I4P, ngc
         do c = 1_I4P, int(size(crown_map, dim=1), I4P)
            b_i8 = crown_map(c, 1, crown)
            if (b_i8 <= 0_I8P) cycle  ! sentinel for unused slot
            fec_i8 = crown_map(c, 9, crown)
            if (FEC_1_6_ARRAY(int(fec_i8, I4P)) /= bc_fec_seam) cycle
            ! Geometric test: is this block on the realm's seam-face boundary?
            if (.not. block_face_on_realm_boundary(realm(my_realm_idx), int(b_i8, I4P), axis, sgn)) cycle
            ! Override the bc_type column (8).
            crown_map(c, 8, crown) = int(BC_SEAM, I8P)
         enddo
      enddo
      end associate
      endsubroutine mark_seam_in_crown

      pure function face_code_to_bc_fec(face_code) result(bc_fec)
      !< Translate adam_maps_object FACE_X_MAX..FACE_Z_MIN (1..6) into the
      !< BC routine's fec_1_6 numbering (1..6 via FEC_1_6_ARRAY).
      !<
      !< Table:
      !<   FACE_X_MAX (1, +x) → fec 2
      !<   FACE_X_MIN (2, -x) → fec 1
      !<   FACE_Y_MAX (3, +y) → fec 4
      !<   FACE_Y_MIN (4, -y) → fec 3
      !<   FACE_Z_MAX (5, +z) → fec 6
      !<   FACE_Z_MIN (6, -z) → fec 5
      !< (pairwise swap of MAX↔MIN within each axis).
      integer(I4P), intent(in) :: face_code
      integer(I4P)             :: bc_fec

      select case (face_code)
      case (1_I4P); bc_fec = 2_I4P
      case (2_I4P); bc_fec = 1_I4P
      case (3_I4P); bc_fec = 4_I4P
      case (4_I4P); bc_fec = 3_I4P
      case (5_I4P); bc_fec = 6_I4P
      case (6_I4P); bc_fec = 5_I4P
      case default; bc_fec = 0_I4P
      end select
      endfunction face_code_to_bc_fec

   endsubroutine populate_inter_realm_topology

   subroutine check_couplings(realm, manifest)
   !< Validate the coupling of every manifest face pair against the realm geometry (issue #52); error_stop otherwise.
   !<
   !< Walks the seam leaves of both sides on the replicated trees (the same result on every rank):
   !<
   !<   * `mirror`: one cell size on each side's seam face, equal across the seam;
   !<   * `refined`: one cell size per side with ratio exactly 2 (either side may be the coarse one), equal block cell
   !<     counts along the seam, the same `[amr] seam_ghost_fill`, and nested blocks: every fine seam leaf covers one 2:1
   !<     quadrant of the face of one coarse seam leaf, and even block cell counts in the fine realm;
   !<   * `periodic`, `interpolate`: schema-reserved, refused (before issue #52 they were accepted and run as `mirror`).
   class(realm_object),     intent(inout) :: realm(:)       !< Forest realms.
   type(forest_manifest_t), intent(in)    :: manifest       !< Parsed manifest.
   integer(I4P)                           :: f              !< Face-pair counter.
   real(R8P)                              :: d_a(3), d_b(3) !< Seam cell sizes of the two sides.
   character(:), allocatable              :: pfx            !< Face-pair label (error message prefix).
   integer(I4P)                           :: ic, jf         !< Coarse and fine realm of a refined pair.
   integer(I4P)                           :: face_f         !< Seam face of the fine realm.
   integer(I4P)                           :: axis, sgn      !< Seam face axis and side of the fine realm.
   integer(I4P)                           :: ntan_c(3)      !< Coarse block cell counts, normal one zeroed.
   integer(I4P)                           :: ntan_f(3)      !< Fine block cell counts, normal one zeroed.

   if (.not. allocated(manifest%face_pairs)) return
   do f = 1_I4P, int(size(manifest%face_pairs), I4P)
      associate(pair => manifest%face_pairs(f))
      pfx = 'forest_object%populate_inter_realm_topology: face_pair '//trim(str(f, .true.))//' (realms '// &
              trim(str(pair%realm_a, .true.))//' and '//trim(str(pair%realm_b, .true.))//'): '
      select case(pair%coupling)
      case(COUPLING_MIRROR, COUPLING_REFINED)
      case default
         call mpih%error_stop(msg=pfx//'only coupling = mirror or refined is implemented (periodic and interpolate '// &
                              'are reserved)')
      endselect
      d_a = seam_cell_size(realm(pair%realm_a), pair%face_a, pfx//'realm '//trim(str(pair%realm_a, .true.)))
      d_b = seam_cell_size(realm(pair%realm_b), pair%face_b, pfx//'realm '//trim(str(pair%realm_b, .true.)))
      if (pair%coupling == COUPLING_MIRROR) then
         if (any(abs(d_a - d_b) > 1.0e-6_R8P * d_a)) &
            call mpih%error_stop(msg=pfx//'coupling = mirror joins cells of the same size, the seam cells differ ('// &
                                 trim(str(d_a))//' and '//trim(str(d_b))//'): use coupling = refined for a 2:1 jump '// &
                                 '(issue #52)')
         cycle
      endif
      ! refined: identify the coarse side
      if (all(abs(d_a - 2._R8P * d_b) <= 1.0e-6_R8P * d_a)) then
         ic = pair%realm_a ; jf = pair%realm_b ; face_f = pair%face_b
      elseif (all(abs(d_b - 2._R8P * d_a) <= 1.0e-6_R8P * d_b)) then
         ic = pair%realm_b ; jf = pair%realm_a ; face_f = pair%face_a
      else
         call mpih%error_stop(msg=pfx//'coupling = refined needs a cell-size ratio of exactly 2 along every axis ('// &
                              trim(str(d_a))//' and '//trim(str(d_b))//')')
      endif
      ! the coarse face skin and the four fine skins under it must index the same cells: equal block cell counts along
      ! the two tangential axes (the normal count is free, the realms may have different extents across the seam)
      call face_axis_sign(face_f, axis, sgn)
      ! `associate` on the polymorphic array elements: nvfortran 26.1 reads their components with a wrong element stride
      ! (measured: realm(2)%adam%grid%nj = 0, nk = 3394; the same defect as the TBP dispatch workaround 0062a237)
      associate(r_c => realm(ic), r_f => realm(jf))
      ntan_c = [r_c%adam%grid%ni, r_c%adam%grid%nj, r_c%adam%grid%nk]
      ntan_f = [r_f%adam%grid%ni, r_f%adam%grid%nj, r_f%adam%grid%nk]
      ntan_c(axis) = 0_I4P ; ntan_f(axis) = 0_I4P
      if (any(ntan_c /= ntan_f)) &
         call mpih%error_stop(msg=pfx//'coupling = refined needs the same block cell counts along the seam (tangential '// &
                              'ni, nj, nk) in both realms')
      if (r_c%adam%maps%seam_ghost_fill /= r_f%adam%maps%seam_ghost_fill)                                             &
         call mpih%error_stop(msg=pfx//'coupling = refined needs the same [amr] seam_ghost_fill in both realms')
      ! a coarse ghost is the mean of 2x2x2 fine cells: even fine block cell counts keep them in one fine block
      if (any(mod([r_f%adam%grid%ni, r_f%adam%grid%nj, r_f%adam%grid%nk], 2_I4P) /= 0_I4P))                        &
         call mpih%error_stop(msg=pfx//'coupling = refined needs even block cell counts (ni, nj, nk) in the fine realm')
      endassociate
      call check_nested(coarse=realm(ic), fine=realm(jf), face_f=face_f, pfx=pfx)
      call mpih%print_message('forest: face_pair '//trim(str(f, .true.))//' is a 2:1 seam, realm '//               &
                              trim(str(ic, .true.))//' coarse, realm '//trim(str(jf, .true.))//' fine (issue #52)')
      endassociate
   enddo
   endsubroutine check_couplings

   subroutine seam_sides(realm, pair, rc, fc, rf, ff, refined)
   !< The coarse and fine sides of a face pair: for a `refined` seam the realm with the larger seam cells is the coarse
   !< side (issue #52), for a `mirror` seam realm_a by convention.
   class(realm_object),      intent(in)  :: realm(:)       !< Forest realms.
   type(forest_face_pair_t), intent(in)  :: pair           !< Face pair.
   integer(I4P),             intent(out) :: rc, fc         !< Coarse-side realm and face.
   integer(I4P),             intent(out) :: rf, ff         !< Fine-side realm and face.
   logical,                  intent(out) :: refined        !< 2:1 seam.
   real(R8P)                             :: d_a(3), d_b(3) !< Seam cell sizes of realm_a and realm_b.

   refined = pair%coupling == COUPLING_REFINED
   rc = pair%realm_a ; fc = pair%face_a ; rf = pair%realm_b ; ff = pair%face_b
   if (.not. refined) return
   d_a = seam_cell_size(realm(pair%realm_a), pair%face_a, 'forest_object%seam_sides: realm_a')
   d_b = seam_cell_size(realm(pair%realm_b), pair%face_b, 'forest_object%seam_sides: realm_b')
   if (d_b(1) > d_a(1)) then
      rc = pair%realm_b ; fc = pair%face_b ; rf = pair%realm_a ; ff = pair%face_a
   endif
   endsubroutine seam_sides

   function seam_cell_size(this_realm, face, pfx) result(dxyz)
   !< The cell size of the leaves of a realm on its seam face: one size for all of them, error_stop otherwise.
   class(realm_object), intent(in) :: this_realm !< Realm.
   integer(I4P),        intent(in) :: face       !< Seam face (FACE_* code).
   character(*),        intent(in) :: pfx      !< Error message prefix.
   real(R8P)                       :: dxyz(3)    !< Seam cell size.
   type(tree_iterator_object)      :: iter       !< Tree traversal cursor.
   type(tree_node_object), pointer :: node_ptr   !< Current leaf.
   integer(I4P)                    :: axis, sgn  !< Face axis and side.
   real(R8P)                       :: emin(3)    !< Leaf origin.
   real(R8P)                       :: d(3)       !< Leaf cell size.
   logical                         :: found      !< A seam leaf has been seen.

   call face_axis_sign(face, axis, sgn)
   found = .false.
   dxyz = 0._R8P
   iter%b = 1_I4P ; iter%p => null()
   do while (this_realm%adam%tree%loop(iter, node_ptr=node_ptr))
      if (.not. leaf_on_realm_face(this_realm, node_ptr%code, axis, sgn)) cycle
      call leaf_metrics(this_realm, node_ptr%code, emin, d)
      if (.not. found) then
         dxyz = d ; found = .true.
      elseif (any(abs(d - dxyz) > 1.0e-6_R8P * dxyz)) then
         call mpih%error_stop(msg=pfx//' has seam cells of different sizes (refined blocks along the seam); a seam '// &
                              'joins one cell size per side')
      endif
   enddo
   if (.not. found) call mpih%error_stop(msg=pfx//' has no block on its seam face')
   endfunction seam_cell_size

   subroutine check_nested(coarse, fine, face_f, pfx)
   !< Every fine seam leaf must cover one 2:1 quadrant of the seam face of one coarse leaf: half its tangential extent,
   !< starting at its origin or at its midpoint along each tangential axis.
   class(realm_object), intent(inout) :: coarse     !< Coarse realm.
   class(realm_object), intent(inout) :: fine       !< Fine realm.
   integer(I4P),        intent(in)    :: face_f     !< Fine seam face.
   character(*),        intent(in)    :: pfx      !< Error message prefix.
   type(tree_iterator_object)         :: iter       !< Tree traversal cursor.
   type(tree_node_object), pointer    :: node_f     !< Fine leaf.
   integer(I8P)                       :: code_c     !< Coarse leaf across the face.
   integer(I4P)                       :: axis, sgn  !< Fine face axis and side.
   integer(I4P)                       :: d          !< Axis counter.
   real(R8P)                          :: emin_f(3), emax_f(3), dx_f(3) !< Fine leaf.
   real(R8P)                          :: emin_c(3), emax_c(3), dx_c(3) !< Coarse leaf.
   real(R8P)                          :: xq(3)      !< First coarse cell centre across the fine face.
   real(R8P)                          :: tol        !< Geometric tolerance.
   logical                            :: ok         !< Nesting test.

   call face_axis_sign(face_f, axis, sgn)
   iter%b = 1_I4P ; iter%p => null()
   do while (fine%adam%tree%loop(iter, node_ptr=node_f))
      if (.not. leaf_on_realm_face(fine, node_f%code, axis, sgn)) cycle
      call leaf_metrics(fine, node_f%code, emin_f, dx_f, emax_f)
      xq = 0.5_R8P * (emin_f + emax_f)
      if (sgn > 0_I4P) then
         xq(axis) = emax_f(axis) + dx_f(axis)
      else
         xq(axis) = emin_f(axis) - dx_f(axis)
      endif
      ok = inside_domain(coarse, xq, 0.25_R8P * dx_f)
      if (ok) then
         code_c = coarse%adam%tree%get_closest_block(grid=coarse%adam%grid, point=xq)
         call leaf_metrics(coarse, code_c, emin_c, dx_c, emax_c)
         tol = max(maxval(abs(emax_c)), 1._R8P) * 1.0e-10_R8P
         ok = merge(abs(emin_c(axis) - emax_f(axis)), abs(emax_c(axis) - emin_f(axis)), sgn > 0_I4P) <= tol
         do d = 1_I4P, 3_I4P
            if (d == axis) cycle
            ! half the coarse extent, starting at the coarse origin or at its midpoint (one 2:1 quadrant)
            ok = ok .and. abs(2._R8P * (emax_f(d) - emin_f(d)) - (emax_c(d) - emin_c(d))) <= tol .and.        &
                 (abs(emin_f(d) - emin_c(d)) <= tol .or. abs(emin_f(d) - 0.5_R8P * (emin_c(d) + emax_c(d))) <= tol)
         enddo
      endif
      if (.not. ok) call mpih%error_stop(msg=pfx//'coupling = refined needs nested blocks: every fine seam block '// &
                                         'must cover a quarter of one coarse seam block face')
   enddo
   endsubroutine check_nested

   subroutine build_seam_rows(realm, manifest)
   !< Build the seam ghost rows of every realm from the replicated trees (issue #40).
   !<
   !< Every rank walks, for every face pair (side A then side B), every leaf of the realm owning the ghosts that lies on
   !< the seam face, owned or not, and every ghost of its slab (tangential edge and corner ghosts included), in the
   !< tree's own order: the same sequence on every rank. The peer cell is the interior cell of the leaf of the peer realm
   !< containing the ghost centre (`get_closest_block` on the peer's replicated tree); ghosts outside the peer domain are
   !< corners where the seam meets a physical boundary and are left to the physical BC. Each rank keeps the rows it
   !< takes part in:
   !<
   !<   * local (ghost and cell here): `[peer, b_cell, b_ghost, i_cell, j_cell, k_cell, i_ghost, j_ghost, k_ghost]`;
   !<   * receive (ghost here, cell on rank r): `[r, b_ghost, i_ghost, j_ghost, k_ghost]` on the ghost realm;
   !<   * send (cell here, ghost on rank r): `[r, b_cell, i_cell, j_cell, k_cell]` on the cell realm.
   !<
   !< Since the order is shared, the rows rank A sends to rank B for a seam are, in order, the rows B receives from A.
   !< The peer cell centre must coincide with the ghost centre: mirror seams join cells of the same size.
   class(realm_object),     intent(inout) :: realm(:) !< Forest realms.
   type(forest_manifest_t), intent(in)    :: manifest !< Parsed manifest.
   type(seam_rows_t), allocatable         :: local(:) !< Local rows, per ghost realm.
   type(seam_rows_t), allocatable         :: recv(:)  !< Receive rows, per ghost realm.
   type(seam_rows_t), allocatable         :: send(:)  !< Send rows, per cell realm.
   integer(I4P)                           :: f, is    !< Face-pair and realm counters.

   allocate(local(size(realm)), recv(size(realm)), send(size(realm)))
   if (allocated(manifest%face_pairs)) then
      do f=1_I4P, int(size(manifest%face_pairs), I4P)
         associate(pair=>manifest%face_pairs(f))
            call enumerate_seam_side(realm=realm, is=pair%realm_a, ip=pair%realm_b, face=pair%face_a, &
                                     local=local(pair%realm_a), recv=recv(pair%realm_a), send=send(pair%realm_b))
            call enumerate_seam_side(realm=realm, is=pair%realm_b, ip=pair%realm_a, face=pair%face_b, &
                                     local=local(pair%realm_b), recv=recv(pair%realm_b), send=send(pair%realm_a))
         endassociate
      enddo
   endif
   do is=1_I4P, int(size(realm), I4P)
      call store_seam_rows(realm=realm(is), is=is, manifest=manifest, local=local(is), recv=recv(is), send=send(is))
   enddo
   endsubroutine build_seam_rows

   subroutine enumerate_seam_side(realm, is, ip, face, local, recv, send)
   !< Enumerate the seam ghosts of realm `is` on `face`, whose peer is realm `ip`; append this rank's rows.
   class(realm_object), intent(inout) :: realm(:)                   !< Forest realms.
   integer(I4P),        intent(in)    :: is                         !< Realm owning the ghosts.
   integer(I4P),        intent(in)    :: ip                         !< Realm owning the cells.
   integer(I4P),        intent(in)    :: face                       !< Seam face of realm `is` (FACE_* code).
   type(seam_rows_t),   intent(inout) :: local                      !< Local rows of realm `is`.
   type(seam_rows_t),   intent(inout) :: recv                       !< Receive rows of realm `is`.
   type(seam_rows_t),   intent(inout) :: send                       !< Send rows of realm `ip`.
   type(tree_iterator_object)         :: iter                       !< Tree traversal cursor.
   type(tree_node_object), pointer    :: node_g                     !< Leaf owning the ghosts.
   type(tree_node_object), pointer    :: node_c                     !< Leaf owning the cell.
   integer(I8P)                       :: code_c                     !< Morton code of the cell leaf.
   integer(I4P)                       :: axis, sgn                  !< Face axis and side.
   integer(I4P)                       :: imin, imax, jmin, jmax     !< Ghost slab bounds.
   integer(I4P)                       :: kmin, kmax                 !< Ghost slab bounds.
   integer(I4P)                       :: i, j, k                    !< Ghost cell indices.
   integer(I4P)                       :: ijk_c(3)                   !< Peer cell indices.
   integer(I4P)                       :: b_g, b_c                   !< Owner-local block indices.
   integer(I4P)                       :: rank_g, rank_c             !< Owner ranks of ghost and cell.
   real(R8P)                          :: emin_g(3), emax_g(3), d_g(3) !< Ghost leaf extent and cell size.
   real(R8P)                          :: emin_c(3), d_c(3)          !< Cell leaf origin and cell size.
   real(R8P)                          :: xg(3), xc(3)               !< Ghost and cell centres.
   real(R8P)                          :: tol                        !< Face tolerance.
   integer(I4P)                       :: n_c(3)                     !< Block cell counts of the cell realm.
   integer(I4P)                       :: kind                       !< Row kind (SEAM_ROW_*).
   integer(I4P)                       :: meta                       !< Packed interpolation metadata.
   integer(I4P)                       :: sub(3), p4(3), p3(3)       !< Octant and anchor positions.
   integer(I4P)                       :: d                          !< Axis counter.
   logical                            :: ok                         !< Cell match.

   call face_axis_sign(face, axis, sgn)
   call ghost_slab_extents(realm(is), axis, sgn, imin, imax, jmin, jmax, kmin, kmax)
   iter%b = 1_I4P ; iter%p => null()
   do while (realm(is)%adam%tree%loop(iter, node_ptr=node_g))
      call leaf_metrics(realm(is), node_g%code, emin_g, d_g, emax_g)
      if (sgn > 0_I4P) then
         tol = max(abs(realm(is)%adam%grid%domain_emax(axis)), 1._R8P) * 1.0e-10_R8P
         if (abs(emax_g(axis) - realm(is)%adam%grid%domain_emax(axis)) > tol) cycle
      else
         tol = max(abs(realm(is)%adam%grid%domain_emin(axis)), 1._R8P) * 1.0e-10_R8P
         if (abs(emin_g(axis) - realm(is)%adam%grid%domain_emin(axis)) > tol) cycle
      endif
      rank_g = node_g%myrank
      b_g    = int(node_g%block_index, I4P)
      do k=kmin, kmax
         do j=jmin, jmax
            do i=imin, imax
               xg = emin_g + ([real(i, R8P), real(j, R8P), real(k, R8P)] - 0.5_R8P) * d_g
               if (.not.inside_domain(realm(ip), xg, 0.25_R8P * d_g)) cycle ! a corner at a physical boundary
               code_c = realm(ip)%adam%tree%get_closest_block(grid=realm(ip)%adam%grid, point=xg)
               call leaf_metrics(realm(ip), code_c, emin_c, d_c)
               ! `associate` on the polymorphic array element: nvfortran 26.1 builds this array constructor from a wrong
               ! element stride (garbage counts, every mirror ghost rejected on FNL; the 0062a237 defect)
               associate(r_p => realm(ip))
                  n_c = [r_p%adam%grid%ni, r_p%adam%grid%nj, r_p%adam%grid%nk]
               endassociate
               meta = 0_I4P
               if (all(abs(d_c - d_g) <= 1.0e-6_R8P * d_g)) then
                  ! same cell size: copy the cell whose centre coincides with the ghost centre
                  kind = SEAM_ROW_COPY
                  ijk_c = nint((xg - emin_c) / d_c + 0.5_R8P, I4P)
                  xc = emin_c + (real(ijk_c, R8P) - 0.5_R8P) * d_c
                  ok = all(abs(xc - xg) <= 1.0e-6_R8P * d_g) .and. all(ijk_c >= 1_I4P) .and. all(ijk_c <= n_c)
               elseif (all(abs(d_c - 2._R8P * d_g) <= 1.0e-6_R8P * d_c)) then
                  ! fine ghost, coarse cells (issue #52): the intra-realm coarse->fine interpolant around the coarse
                  ! cell containing the ghost centre, the octant from the half of it the centre falls in
                  kind = SEAM_ROW_INTERPOLATE
                  ijk_c = floor((xg - emin_c) / d_c, I4P) + 1_I4P
                  ok = all(ijk_c >= 1_I4P) .and. all(ijk_c <= n_c)
                  if (ok) then
                     do d = 1_I4P, 3_I4P
                        sub(d) = merge(1_I4P, 2_I4P, xg(d) < emin_c(d) + (real(ijk_c(d), R8P) - 0.5_R8P) * d_c(d))
                        p4(d) = seam_shift_anchor_pos(anchor=ijk_c(d), n_cells=n_c(d),                       &
                                                      p_centered=seam_tricubic_centered_pos(sub(d)), footprint_n=4_I4P)
                        p3(d) = seam_shift_anchor_pos(anchor=ijk_c(d), n_cells=n_c(d),                       &
                                                      p_centered=seam_compatible_centered_pos(sub(d)), footprint_n=3_I4P)
                     enddo
                     meta = seam_meta_pack(sub=sub, p4=p4, p3=p3)
                  endif
               elseif (all(abs(2._R8P * d_c - d_g) <= 1.0e-6_R8P * d_g)) then
                  ! coarse ghost, fine cells (issue #52): the mean of the 2x2x2 fine cells under the ghost, from the
                  ! base (lowest) one; nested blocks with even cell counts keep the eight in one fine block
                  kind = SEAM_ROW_RESTRICT
                  xc = xg - 0.5_R8P * d_c
                  code_c = realm(ip)%adam%tree%get_closest_block(grid=realm(ip)%adam%grid, point=xc)
                  call leaf_metrics(realm(ip), code_c, emin_c, d_c)
                  ijk_c = nint((xc - emin_c) / d_c + 0.5_R8P, I4P)
                  ok = all(ijk_c >= 1_I4P) .and. all(ijk_c + 1_I4P <= n_c) .and. &
                       all(abs(emin_c + (real(ijk_c, R8P) - 0.5_R8P) * d_c - xc) <= 1.0e-6_R8P * d_c)
               else
                  ok = .false.
               endif
               if (.not. ok) call mpih%error_stop(msg='forest_object%populate_inter_realm_topology: the seam '//  &
                  ghost_label(is, b_g, i, j, k)//' does not match the cells of realm '//trim(str(ip, .true.))//    &
                  ' (issues #40, #52): a seam joins cells of the same size (mirror) or of ratio 2 (refined)')
               node_c => realm(ip)%adam%tree%node(code=code_c)
               rank_c = node_c%myrank
               b_c    = int(node_c%block_index, I4P)
               if (kind == SEAM_ROW_COPY .and. rank_g == mpih%myrank .and. rank_c == mpih%myrank) then
                  call push_row(local, [ip, b_c, b_g, ijk_c, i, j, k])
               else
                  ! cross-rank rows, and same-rank interpolate/restrict rows as messages to self
                  if (rank_g == mpih%myrank) call push_row(recv, [ip, rank_c, b_g, i, j, k])
                  if (rank_c == mpih%myrank) call push_row(send, [is, rank_g, b_c, ijk_c, kind, meta])
               endif
            enddo
         enddo
      enddo
   enddo
   endsubroutine enumerate_seam_side

   subroutine store_seam_rows(realm, is, manifest, local, recv, send)
   !< Store the enumerated seam rows of realm `is` into its maps, grouped by peer slot (and by rank inside a slot).
   class(realm_object),     intent(inout) :: realm    !< Realm `is`.
   integer(I4P),            intent(in)    :: is       !< Realm index.
   type(forest_manifest_t), intent(in)    :: manifest !< Parsed manifest.
   type(seam_rows_t),       intent(in)    :: local    !< Local rows (col 1 = peer realm).
   type(seam_rows_t),       intent(in)    :: recv     !< Receive rows (col 1 = peer realm).
   type(seam_rows_t),       intent(in)    :: send     !< Send rows (col 1 = ghost realm).
   integer(I4P)                           :: peers(2*max(1, size_pairs(manifest))) !< Peer realms, manifest order.
   integer(I4P)                           :: n_peers  !< Number of peer slots.
   integer(I4P)                           :: f, p, c  !< Counters.
   integer(I4P)                           :: other    !< Realm on the other side of a face pair.
   integer(I4P)                           :: cadence  !< Cadence of a slot.

   associate(maps=>realm%adam%maps)
   if (allocated(maps%seam_local_map_ghost_cell)) deallocate(maps%seam_local_map_ghost_cell)
   if (allocated(maps%seam_local_peer_realm))     deallocate(maps%seam_local_peer_realm)
   if (allocated(maps%seam_local_peer_row_start)) deallocate(maps%seam_local_peer_row_start)
   if (allocated(maps%seam_local_peer_row_count)) deallocate(maps%seam_local_peer_row_count)
   if (allocated(maps%seam_local_cadence))        deallocate(maps%seam_local_cadence)
   if (allocated(maps%seam_local_send_buf))       deallocate(maps%seam_local_send_buf)
   if (allocated(maps%seam_local_recv_buf))       deallocate(maps%seam_local_recv_buf)
   if (allocated(maps%seam_mpi_recv_cell))        deallocate(maps%seam_mpi_recv_cell)
   if (allocated(maps%seam_mpi_recv_row_start))   deallocate(maps%seam_mpi_recv_row_start)
   if (allocated(maps%seam_mpi_recv_row_count))   deallocate(maps%seam_mpi_recv_row_count)
   if (allocated(maps%seam_mpi_send_cell))        deallocate(maps%seam_mpi_send_cell)
   if (allocated(maps%seam_mpi_send_row_start))   deallocate(maps%seam_mpi_send_row_start)
   if (allocated(maps%seam_mpi_send_row_count))   deallocate(maps%seam_mpi_send_row_count)
   ! peer slots: the realms glued to `is`, in manifest order (the same on every rank)
   n_peers = 0_I4P
   if (allocated(manifest%face_pairs)) then
      do f=1_I4P, int(size(manifest%face_pairs), I4P)
         other = 0_I4P
         if (manifest%face_pairs(f)%realm_a == is) other = manifest%face_pairs(f)%realm_b
         if (manifest%face_pairs(f)%realm_b == is) other = manifest%face_pairs(f)%realm_a
         if (other == 0_I4P) cycle
         if (any(peers(1:n_peers) == other)) cycle
         n_peers = n_peers + 1_I4P
         peers(n_peers) = other
      enddo
   endif
   if (n_peers == 0_I4P) return
   allocate(maps%seam_local_peer_realm(n_peers), source=peers(1:n_peers))
   ! cadence of each slot: a property of the realm pair, every face pair joining them must agree
   allocate(maps%seam_local_cadence(n_peers))
   do p=1_I4P, n_peers
      cadence = -1_I4P
      do f=1_I4P, int(size(manifest%face_pairs), I4P)
         associate(pair=>manifest%face_pairs(f))
         if (.not.((pair%realm_a == is .and. pair%realm_b == peers(p)) .or. &
                   (pair%realm_b == is .and. pair%realm_a == peers(p)))) cycle
         if (cadence == -1_I4P) then
            cadence = pair%coupling_cadence
         elseif (cadence /= pair%coupling_cadence) then
            call mpih%error_stop(msg='forest_object%populate_inter_realm_topology: conflicting coupling_cadence '// &
                                 'between realm '//trim(str(is, .true.))//' and realm '//trim(str(peers(p), .true.))//&
                                 ' across multiple face_pairs')
         endif
         endassociate
      enddo
      maps%seam_local_cadence(p) = cadence
   enddo
   ! local rows, grouped by slot
   allocate(maps%seam_local_peer_row_start(n_peers), maps%seam_local_peer_row_count(n_peers))
   maps%seam_local_peer_row_count = 0_I4P
   if (local%n > 0_I4P) then
      do p=1_I4P, n_peers
         maps%seam_local_peer_row_count(p) = count(local%row(1:local%n, 1) == peers(p))
      enddo
   endif
   maps%seam_local_peer_row_start(1) = 1_I4P
   do p=2_I4P, n_peers
      maps%seam_local_peer_row_start(p) = maps%seam_local_peer_row_start(p-1) + maps%seam_local_peer_row_count(p-1)
   enddo
   if (local%n > 0_I4P) then
      allocate(maps%seam_local_map_ghost_cell(local%n, 9))
      c = 0_I4P
      do p=1_I4P, n_peers
         do f=1_I4P, local%n
            if (local%row(f, 1) /= peers(p)) cycle
            c = c + 1_I4P
            maps%seam_local_map_ghost_cell(c, :) = local%row(f, 1:9)
         enddo
      enddo
      allocate(maps%seam_local_send_buf(realm%nv * maxval(maps%seam_local_peer_row_count), n_peers), source=0._R8P)
      allocate(maps%seam_local_recv_buf(realm%nv * maxval(maps%seam_local_peer_row_count), n_peers), source=0._R8P)
   endif
   ! cross-rank rows, grouped by slot and by rank
   call group_rows(rows=recv, peers=peers(1:n_peers), cell=maps%seam_mpi_recv_cell,      &
                   row_start=maps%seam_mpi_recv_row_start, row_count=maps%seam_mpi_recv_row_count)
   call group_rows(rows=send, peers=peers(1:n_peers), cell=maps%seam_mpi_send_cell,      &
                   row_start=maps%seam_mpi_send_row_start, row_count=maps%seam_mpi_send_row_count)
   call mpih%print_message('forest: realm '//trim(str(is, .true.))//' seam rows: '//trim(str(local%n, .true.))// &
                           ' local, '//trim(str(recv%n, .true.))//' received, '//trim(str(send%n, .true.))//' sent')
   endassociate
   endsubroutine store_seam_rows

   subroutine group_rows(rows, peers, cell, row_start, row_count)
   !< Group cross-rank rows `[peer, rank, b, i, j, k, ...]` by peer slot and, inside a slot, by ascending rank, keeping
   !< the enumeration order within a rank; store them without the peer column (`[rank, b, i, j, k]`, plus `kind, meta`
   !< for send rows).
   type(seam_rows_t),         intent(in)  :: rows         !< Enumerated rows.
   integer(I4P),              intent(in)  :: peers(:)     !< Peer realm of each slot.
   integer(I4P), allocatable, intent(out) :: cell(:,:)    !< Grouped rows.
   integer(I4P), allocatable, intent(out) :: row_start(:) !< First row of each slot.
   integer(I4P), allocatable, intent(out) :: row_count(:) !< Rows of each slot.
   integer(I4P)                           :: p, r, c      !< Counters.
   integer(I4P)                           :: rank, next   !< Current and next rank of a slot.

   allocate(row_start(size(peers)), row_count(size(peers)))
   row_count = 0_I4P ; row_start = 1_I4P
   if (rows%n > 0_I4P) allocate(cell(rows%n, size(rows%row, dim=2) - 1))
   c = 0_I4P
   do p=1_I4P, int(size(peers), I4P)
      row_start(p) = c + 1_I4P
      rank = -1_I4P
      do
         next = huge(1_I4P)
         do r=1_I4P, rows%n
            if (rows%row(r, 1) == peers(p) .and. rows%row(r, 2) > rank) next = min(next, rows%row(r, 2))
         enddo
         if (next == huge(1_I4P)) exit
         rank = next
         do r=1_I4P, rows%n
            if (rows%row(r, 1) /= peers(p) .or. rows%row(r, 2) /= rank) cycle
            c = c + 1_I4P
            cell(c, :) = rows%row(r, 2:)
         enddo
      enddo
      row_count(p) = c - row_start(p) + 1_I4P
   enddo
   endsubroutine group_rows

   subroutine push_row(list, row)
   !< Append a row to a growable row list (the capacity doubles when full).
   type(seam_rows_t), intent(inout) :: list    !< Row list.
   integer(I4P),      intent(in)    :: row(:)  !< Row to append.
   integer(I4P), allocatable        :: tmp(:,:) !< Grown storage.

   if (.not.allocated(list%row)) allocate(list%row(64, size(row)))
   if (list%n == size(list%row, dim=1)) then
      allocate(tmp(2 * list%n, size(row)))
      tmp(1:list%n, :) = list%row(1:list%n, :)
      call move_alloc(from=tmp, to=list%row)
   endif
   list%n = list%n + 1_I4P
   list%row(list%n, :) = row
   endsubroutine push_row

   pure function size_pairs(manifest) result(n)
   !< Number of face pairs of a manifest.
   type(forest_manifest_t), intent(in) :: manifest !< Parsed manifest.
   integer(I4P)                        :: n        !< Face pairs.

   n = 0_I4P
   if (allocated(manifest%face_pairs)) n = int(size(manifest%face_pairs), I4P)
   endfunction size_pairs

   subroutine leaf_metrics(this_realm, code, emin, dxyz, emax)
   !< Origin, cell size and (optionally) upper corner of the leaf `code` of a realm, from its tree and grid (valid on
   !< every rank, whichever rank owns the leaf).
   class(realm_object), intent(in)            :: this_realm !< Realm.
   integer(I8P),        intent(in)            :: code       !< Leaf Morton code.
   real(R8P),           intent(out)           :: emin(3)    !< Leaf origin.
   real(R8P),           intent(out)           :: dxyz(3)    !< Cell size.
   real(R8P),           intent(out), optional :: emax(3)    !< Leaf upper corner.
   integer(I4P)                               :: ijkl(4)    !< Leaf coordinates and level.

   call this_realm%adam%tree%morton_to_coordinates(code=code, i=ijkl(1), j=ijkl(2), k=ijkl(3), l=ijkl(4))
   call this_realm%adam%grid%compute_metrics(coordinates=ijkl, emin=emin, emax=emax, dx=dxyz(1), dy=dxyz(2), dz=dxyz(3))
   endsubroutine leaf_metrics

   function leaf_on_realm_face(this_realm, code, axis, sgn) result(yes)
   !< Return .true. iff the leaf `code` of a realm has its (axis, sgn) face on the realm boundary (any owner rank).
   class(realm_object), intent(in) :: this_realm !< Realm.
   integer(I8P),        intent(in) :: code       !< Leaf Morton code.
   integer(I4P),        intent(in) :: axis, sgn  !< Face axis (1..3) and side (+1 max, -1 min).
   logical                         :: yes        !< Test result.
   real(R8P)                       :: emin(3)    !< Leaf origin.
   real(R8P)                       :: emax(3)    !< Leaf upper corner.
   real(R8P)                       :: dxyz(3)    !< Cell size.
   real(R8P)                       :: bound      !< Realm boundary coordinate.

   call leaf_metrics(this_realm, code, emin, dxyz, emax)
   if (sgn > 0_I4P) then
      bound = this_realm%adam%grid%domain_emax(axis)
      yes = abs(emax(axis) - bound) <= max(abs(bound), 1._R8P) * 1.0e-10_R8P
   else
      bound = this_realm%adam%grid%domain_emin(axis)
      yes = abs(emin(axis) - bound) <= max(abs(bound), 1._R8P) * 1.0e-10_R8P
   endif
   endfunction leaf_on_realm_face

   pure function meet_face_to_face(emin_a, emax_a, emin_b, emax_b, axis, sgn) result(yes)
   !< Return .true. iff box b lies across the (axis, sgn) face of box a and covers exactly that face: same face plane,
   !< same extent on the two tangential axes.
   real(R8P),    intent(in) :: emin_a(3), emax_a(3) !< Box a.
   real(R8P),    intent(in) :: emin_b(3), emax_b(3) !< Box b.
   integer(I4P), intent(in) :: axis, sgn            !< Face of box a: axis (1..3) and side (+1 max, -1 min).
   logical                  :: yes                  !< Test result.
   real(R8P)                :: tol                  !< Geometric tolerance.
   integer(I4P)             :: d                    !< Axis counter.

   tol = max(maxval(abs(emin_a)), maxval(abs(emax_a)), 1._R8P) * 1.0e-10_R8P
   if (sgn > 0_I4P) then
      yes = abs(emin_b(axis) - emax_a(axis)) <= tol
   else
      yes = abs(emax_b(axis) - emin_a(axis)) <= tol
   endif
   do d=1_I4P, 3_I4P
      if (d == axis) cycle
      yes = yes .and. abs(emin_b(d) - emin_a(d)) <= tol .and. abs(emax_b(d) - emax_a(d)) <= tol
   enddo
   endfunction meet_face_to_face

   pure subroutine ghost_slab_extents(this_realm, axis, sgn, imin_g, imax_g, jmin_g, jmax_g, kmin_g, kmax_g)
   !< The (i, j, k) bounds of a block's ghost slab on one face: `ngc` cells deep on the (axis, sgn) face, the two
   !< tangential axes over the full ghost-extended range `[1-ngc, n+ngc]` (edge and corner ghosts included).
   class(realm_object), intent(in)  :: this_realm                                     !< Realm.
   integer(I4P),        intent(in)  :: axis, sgn                                      !< Face axis and side.
   integer(I4P),        intent(out) :: imin_g, imax_g, jmin_g, jmax_g, kmin_g, kmax_g !< Slab bounds.

   associate(g => this_realm%adam%grid)
   imin_g = 1_I4P - g%ngc ; imax_g = g%ni + g%ngc
   jmin_g = 1_I4P - g%ngc ; jmax_g = g%nj + g%ngc
   kmin_g = 1_I4P - g%ngc ; kmax_g = g%nk + g%ngc
   select case (axis)
   case (1_I4P)
      if (sgn > 0_I4P) then ; imin_g = g%ni + 1_I4P ; imax_g = g%ni + g%ngc
      else                  ; imin_g = 1_I4P - g%ngc ; imax_g = 0_I4P
      endif
   case (2_I4P)
      if (sgn > 0_I4P) then ; jmin_g = g%nj + 1_I4P ; jmax_g = g%nj + g%ngc
      else                  ; jmin_g = 1_I4P - g%ngc ; jmax_g = 0_I4P
      endif
   case (3_I4P)
      if (sgn > 0_I4P) then ; kmin_g = g%nk + 1_I4P ; kmax_g = g%nk + g%ngc
      else                  ; kmin_g = 1_I4P - g%ngc ; kmax_g = 0_I4P
      endif
   endselect
   endassociate
   endsubroutine ghost_slab_extents

   pure function inside_domain(this_realm, xc, tol) result(yes)
   !< Return .true. iff `xc` lies inside the domain of `this_realm`, at least `tol` away from its boundary.
   class(realm_object), intent(in) :: this_realm !< Realm to query.
   real(R8P),           intent(in) :: xc(3)      !< Point.
   real(R8P),           intent(in) :: tol(3)     !< Per-axis margin.
   logical                         :: yes        !< Test result.

   yes = all(xc > this_realm%adam%grid%domain_emin + tol) .and. all(xc < this_realm%adam%grid%domain_emax - tol)
   endfunction inside_domain

   function ghost_label(is, b, i, j, k) result(label)
   !< Name a seam ghost cell in an error message.
   integer(I4P), intent(in)  :: is         !< Realm index.
   integer(I4P), intent(in)  :: b, i, j, k !< Owner-local block and ghost cell indices.
   character(:), allocatable :: label      !< Label.

   label = 'ghost ('//trim(str(i))//','//trim(str(j))//','//trim(str(k))//') of block '//trim(str(b, .true.))// &
           ' of realm '//trim(str(is, .true.))
   endfunction ghost_label

   subroutine register_intra_realm_amr_seams(self, realm, extra_faces, nfaces_intra)
   !< Register every intra-realm AMR coarse-fine face in the forest flux register.
   !<
   !< Walks each realm's tree node neighborhood (already built by
   !< `make_neighborhood` at realm init) and, for every block whose face neighbor
   !< is MORE refined (`NODE_MORE_REFINED` — this block is the COARSE side of a
   !< 2:1 jump), adds one register entry. The fine-side blocks (the `ratio/2`
   !< finer neighbors covering that coarse face) are recorded in `fine_block(:)`.
   !<
   !< Only the 6 faces (`fec = 1..6`) are registered — reflux is a face-flux
   !< correction; edge/corner adjacencies (`fec = 7..26`) carry no conserved
   !< face flux and are skipped.
   !<
   !< Two-pass, mirroring `register_inter_realm_seams`:
   !<   * Pass 1 counts coarse-side faces across all realms → register size.
   !<   * Pass 2 calls `register_face(seam_kind = SEAM_KIND_INTRA_REALM_AMR, ...)`
   !<     and fills the signed `inter_realm_face_register_index` lookup: `+cursor`
   !<     on the coarse (block, bc_fec), `-cursor` on each fine neighbor block's
   !<     opposite face. For intra-realm jumps coarse and fine are the SAME realm,
   !<     so `coarse_realm = fine_realm = is`.
   !<
   !< **Register ownership.** This routine is the sole initializer of the flux
   !< register: it calls `flux_register%initialize` with the intra-realm face
   !< count (possibly 0, which still flips the register's `is_initialized_` so
   !< the per-step `reset`/reflux hooks are safe no-ops) plus `extra_faces`.
   !<
   !< **Composition with inter-realm seams (issue #37).** The manifest path
   !< calls it first with `extra_faces` = the inter-realm seam face count,
   !< then `register_inter_realm_seams` appends those faces after the
   !< `nfaces_intra` intra-realm ones. Both lists are walked on the
   !< replicated trees (issue #40), so every rank holds the same faces under
   !< the same cursors, which is what the register's cross-rank reduction of
   !< the fine sums relies on.
   class(forest_object), intent(inout)         :: self           !< The forest.
   class(realm_object),  intent(inout)         :: realm(:)       !< Initialized realms whose trees are walked.
   integer(I4P),         intent(in),  optional :: extra_faces    !< Register slots reserved for inter-realm faces.
   integer(I4P),         intent(out), optional :: nfaces_intra   !< Intra-realm AMR faces registered (cursors 1..n).
   integer(I4P)                                :: is             !< Realm index.
   integer(I4P)                                :: fec            !< Face/edge/corner direction (only 1..6 used).
   integer(I4P)                                :: nfaces_total   !< Total coarse-side AMR faces across the forest.
   integer(I4P)                                :: cursor         !< Write cursor into the register.
   integer(I4P)                                :: axis           !< Coarse-face axis (1=x,2=y,3=z) from the tree fec.
   integer(I4P)                                :: nface_cells    !< Coarse-face skin cell count for one block.
   integer(I4P)                                :: n_fine         !< Number of fine neighbor blocks on a coarse face.
   integer(I4P)                                :: kf             !< Fine-neighbor counter.
   type(tree_iterator_object)                  :: iter           !< Tree traversal cursor.
   type(tree_node_object), pointer             :: node_ptr       !< Current node.
   type(tree_node_object), pointer             :: fine_ptr       !< Fine neighbor node.
   integer(I4P), allocatable                   :: fine_blocks(:) !< Fine-side block indices on a coarse face.

   ! Pass 1: count coarse-side AMR faces.
   nfaces_total = 0_I4P
   do is = 1_I4P, int(size(realm), I4P)
      iter%b = 1_I4P ; iter%p => null()
      do while (realm(is)%adam%tree%loop(iter, node_ptr=node_ptr))
         do fec = 1_I4P, 6_I4P
            if (.not. allocated(node_ptr%neighbor(fec)%codes)) cycle
            if (node_ptr%neighbor(fec)%ntype == NODE_MORE_REFINED) nfaces_total = nfaces_total + 1_I4P
         enddo
      enddo
   enddo

   if (present(nfaces_intra)) nfaces_intra = nfaces_total
   if (present(extra_faces)) then
      call self%flux_register%initialize(nfaces=nfaces_total + extra_faces)
   else
      call self%flux_register%initialize(nfaces=nfaces_total)
   endif
   call mpih%print_message('forest: registered intra-realm AMR seam faces: '//trim(str(nfaces_total)))

   ! Allocate the per-realm (block, bc_fec) → signed register-index lookup. Same
   ! array the inter-realm pass uses and the FV reflux hook consumes; here it
   ! carries intra-realm AMR seam faces. Zero = not a seam face. The companion
   ! quadrant table (issue #28 D2) is allocated ONLY by this intra-realm pass —
   ! its allocation status is what the FV accumulation hooks use to discriminate
   ! intra-realm AMR seams (2:1 quadrant restriction) from inter-realm mirror
   ! seams (no quadrant).
   do is = 1_I4P, int(size(realm), I4P)
      block
         integer(I4P) :: nb_realm
         nb_realm = int(realm(is)%adam%field%blocks_number, I4P)
         if (allocated(realm(is)%adam%maps%inter_realm_face_register_index)) &
            deallocate(realm(is)%adam%maps%inter_realm_face_register_index)
         if (allocated(realm(is)%adam%maps%amr_seam_quadrant)) &
            deallocate(realm(is)%adam%maps%amr_seam_quadrant)
         if (nb_realm > 0_I4P) then
            allocate(realm(is)%adam%maps%inter_realm_face_register_index(1:nb_realm, 1:6))
            realm(is)%adam%maps%inter_realm_face_register_index = 0_I4P
            allocate(realm(is)%adam%maps%amr_seam_quadrant(1:2, 1:nb_realm, 1:6))
            realm(is)%adam%maps%amr_seam_quadrant = 0_I4P
         endif
      endblock
   enddo
   if (nfaces_total == 0_I4P) return

   ! Pass 2: register each coarse-side AMR face and fill the signed lookup.
   cursor = 0_I4P
   do is = 1_I4P, int(size(realm), I4P)
      iter%b = 1_I4P ; iter%p => null()
      do while (realm(is)%adam%tree%loop(iter, node_ptr=node_ptr))
         do fec = 1_I4P, 6_I4P
            if (.not. allocated(node_ptr%neighbor(fec)%codes)) cycle
            if (node_ptr%neighbor(fec)%ntype /= NODE_MORE_REFINED) cycle
            cursor = cursor + 1_I4P
            ! Tree fec 1..6 are faces ±x,±y,±z in order (FEC_TO_DELTA): axis = (fec+1)/2.
            axis = (fec + 1_I4P) / 2_I4P
            nface_cells = tangential_cells(realm(is), axis)
            ! Resolve the fine-side block indices (the ratio/2 finer neighbors).
            n_fine = int(size(node_ptr%neighbor(fec)%codes), I4P)
            if (allocated(fine_blocks)) deallocate(fine_blocks)
            allocate(fine_blocks(1:n_fine))
            do kf = 1_I4P, n_fine
               fine_ptr => realm(is)%adam%tree%node(code=node_ptr%neighbor(fec)%codes(kf))
               fine_blocks(kf) = int(fine_ptr%block_index, I4P)
            enddo
            call self%flux_register%register_face(face_index=cursor,                            &
                                                  seam_kind=SEAM_KIND_INTRA_REALM_AMR,          &
                                                  coarse_realm=is,                              &
                                                  coarse_rank=node_ptr%myrank,                  &
                                                  coarse_block=int(node_ptr%block_index, I4P),  &
                                                  coarse_face=fec_to_face_code(fec),            &
                                                  fine_realm=is,                                &
                                                  fine_block=fine_blocks,                       &
                                                  nface_cells=nface_cells,                      &
                                                  nv=int(realm(is)%adam%field%nv, I4P),         &
                                                  n_stages=realm(is)%stages_per_step_forest())
            ! Signed lookup, keyed by the TREE fec (1=-x,2=+x,3=-y,4=+y,5=-z,6=+z)
            ! — the SAME numbering `accumulate_seam_fluxes_fv`'s `select case (fec)`
            ! uses to pack the face skin. (The node's `%bc_fec` field is populated
            ! only for boundary-condition neighbors, not interior AMR neighbors, so
            ! it is unusable here.) +cursor on the coarse (block, fec); -cursor on
            ! each fine neighbor's opposite face (fec 1↔2, 3↔4, 5↔6).
            !
            ! Issue #28 D1: both writes are OWNERSHIP-GATED. `block_index` is an
            ! owner-rank-LOCAL field slot (assigned per-proc by the tree
            ! redistribution), while this walk covers the REPLICATED tree — every
            ! rank sees every node. Writing a remote node's entry would land the
            ! ±cursor on an unrelated LOCAL block that happens to share the local
            ! index (registration aliasing: at np2 on rmf-amr the later
            ! Morton-order write won and two physical faces swapped register
            ! entries) and is an out-of-bounds write whenever the remote rank
            ! stores more blocks than this one. The register face LIST stays
            ! replicated (same cursor order on every rank); only the rank-local
            ! LOOKUP is ownership-filtered.
            if (node_ptr%myrank == mpih%myrank) &
               realm(is)%adam%maps%inter_realm_face_register_index(int(node_ptr%block_index, I4P), fec) = +cursor
            do kf = 1_I4P, n_fine
               fine_ptr => realm(is)%adam%tree%node(code=node_ptr%neighbor(fec)%codes(kf))
               if (fine_ptr%myrank /= mpih%myrank) cycle
               realm(is)%adam%maps%inter_realm_face_register_index(int(fine_ptr%block_index, I4P), &
                                                                   opposite_fec(fec)) = -cursor
               ! Issue #28 D2: precompute this fine block's 2:1 quadrant offset
               ! within the coarse face skin from the two Morton codes — pure
               ! tree topology, valid regardless of which rank owns the coarse
               ! block (whose emin/emax the accumulation hooks previously read
               ! through an owner-local index, garbage when it is remote).
               block
                  integer(I4P) :: ioff, joff
                  call seam_quadrant_from_codes(realm(is), node_ptr%code, fine_ptr%code, fec, ioff, joff)
                  realm(is)%adam%maps%amr_seam_quadrant(1:2, int(fine_ptr%block_index, I4P), &
                                                        opposite_fec(fec)) = [ioff, joff]
               endblock
            enddo
         enddo
      enddo
   enddo
   if (allocated(fine_blocks)) deallocate(fine_blocks)


   contains
      pure function tangential_cells(this_realm, axis) result(n)
      !< Single-block face-skin cell count tangential to `axis` (nj*nk, ni*nk, ni*nj).
      class(realm_object), intent(in) :: this_realm !< Realm to query.
      integer(I4P),        intent(in) :: axis       !< 1=x, 2=y, 3=z.
      integer(I4P)                    :: n          !< Cell count.

      associate(g => this_realm%adam%grid)
         select case (axis)
         case (1_I4P); n = g%nj * g%nk
         case (2_I4P); n = g%ni * g%nk
         case (3_I4P); n = g%ni * g%nj
         case default; n = 0_I4P
         endselect
      endassociate
      endfunction tangential_cells

      pure function opposite_fec(fec) result(opp)
      !< Opposite tree fec within the same axis (1↔2, 3↔4, 5↔6).
      integer(I4P), intent(in) :: fec !< Tree fec (1..6).
      integer(I4P)             :: opp !< Opposite tree fec.

      select case (fec)
      case (1_I4P); opp = 2_I4P
      case (2_I4P); opp = 1_I4P
      case (3_I4P); opp = 4_I4P
      case (4_I4P); opp = 3_I4P
      case (5_I4P); opp = 6_I4P
      case (6_I4P); opp = 5_I4P
      case default; opp = 0_I4P
      endselect
      endfunction opposite_fec

      pure function fec_to_face_code(fec) result(face_code)
      !< Translate a tree fec (1=-x, 2=+x, 3=-y, 4=+y, 5=-z, 6=+z) into the FACE_* code the flux register stores and
      !< `face_axis_sign` decodes (FACE_X_MAX=1 is +x, FACE_X_MIN=2 is -x, ...).
      !<
      !< The two numberings are swapped within each axis: registering the raw tree fec as `coarse_face` (as this pass did
      !< before) makes every reflux application decode the opposite side, so the correction had the wrong sign
      !< (doubling the coarse-fine conservation defect instead of removing it) and landed on the far cell layer of the
      !< coarse block.
      integer(I4P), intent(in) :: fec       !< Tree fec (1..6).
      integer(I4P)             :: face_code !< FACE_* code.

      select case (fec)
      case (1_I4P); face_code = FACE_X_MIN
      case (2_I4P); face_code = FACE_X_MAX
      case (3_I4P); face_code = FACE_Y_MIN
      case (4_I4P); face_code = FACE_Y_MAX
      case (5_I4P); face_code = FACE_Z_MIN
      case (6_I4P); face_code = FACE_Z_MAX
      case default; face_code = 0_I4P
      endselect
      endfunction fec_to_face_code

      subroutine seam_quadrant_from_codes(this_realm, coarse_code, fine_code, fec, ioff, joff)
      !< Fine block's quadrant offset (inner, outer) ∈ {0,1}² within the coarse
      !< face skin, from the two nodes' Morton codes (issue #28 D2).
      !<
      !< At the fine level the coarse block spans integer coordinates
      !< `2c .. 2c+1` along every axis, so the offset along a tangential axis is
      !< `coord_fine − 2·coord_coarse`. Pure tree topology — no geometry reads,
      !< valid regardless of which ranks own the two blocks. The tangential
      !< (inner, outer) axis pair per face matches the accumulation hooks'
      !< convention: x-faces → (y, z); y-faces → (x, z); z-faces → (x, y).
      class(realm_object), intent(in)  :: this_realm  !< Realm owning the seam (both sides).
      integer(I8P),        intent(in)  :: coarse_code !< Coarse node Morton code.
      integer(I8P),        intent(in)  :: fine_code   !< Fine node Morton code.
      integer(I4P),        intent(in)  :: fec         !< Coarse-face tree fec (1..6).
      integer(I4P),        intent(out) :: ioff        !< Quadrant offset along the inner tangential axis.
      integer(I4P),        intent(out) :: joff        !< Quadrant offset along the outer tangential axis.
      integer(I4P)                     :: cc(3)       !< Coarse node (i,j,k) at its own level.
      integer(I4P)                     :: fc(3)       !< Fine node (i,j,k) at its own level.
      integer(I4P)                     :: lc, lf      !< Coarse / fine node levels.
      integer(I4P)                     :: inner_ax    !< Inner tangential axis (1=x,2=y,3=z).
      integer(I4P)                     :: outer_ax    !< Outer tangential axis (1=x,2=y,3=z).

      cc = 0_I4P ; fc = 0_I4P
      associate(tree => this_realm%adam%tree)
         select case (tree%ratio)
         case (2_I4P)
            call tree%morton_to_coordinates(code=coarse_code, i=cc(1), l=lc)
            call tree%morton_to_coordinates(code=fine_code,   i=fc(1), l=lf)
         case (4_I4P)
            call tree%morton_to_coordinates(code=coarse_code, i=cc(1), j=cc(2), l=lc)
            call tree%morton_to_coordinates(code=fine_code,   i=fc(1), j=fc(2), l=lf)
         case default
            call tree%morton_to_coordinates(code=coarse_code, i=cc(1), j=cc(2), k=cc(3), l=lc)
            call tree%morton_to_coordinates(code=fine_code,   i=fc(1), j=fc(2), k=fc(3), l=lf)
         endselect
      endassociate
      if (lf /= lc + 1_I4P) &
         call mpih%error_stop(msg='register_intra_realm_amr_seams: NODE_MORE_REFINED neighbor is not one level finer')
      select case ((fec + 1_I4P) / 2_I4P)
      case (1_I4P)  ; inner_ax = 2_I4P ; outer_ax = 3_I4P
      case (2_I4P)  ; inner_ax = 1_I4P ; outer_ax = 3_I4P
      case default  ; inner_ax = 1_I4P ; outer_ax = 2_I4P
      endselect
      ioff = fc(inner_ax) - 2_I4P * cc(inner_ax)
      joff = fc(outer_ax) - 2_I4P * cc(outer_ax)
      if (ioff < 0_I4P .or. ioff > 1_I4P .or. joff < 0_I4P .or. joff > 1_I4P) &
         call mpih%error_stop(msg='register_intra_realm_amr_seams: fine-block quadrant offset outside {0,1} — '// &
                                  'tree topology inconsistent with a 2:1 face jump')
      endsubroutine seam_quadrant_from_codes
   endsubroutine register_intra_realm_amr_seams

   subroutine apply_reflux_corrections(self, realm, dt)
   !< Dispatch the Berger-Colella reflux correction to every realm.
   !<
   !< Called once per step AFTER close_step_forest's update_q. The forest's
   !< role is purely orchestration: it iterates the realm array and invokes
   !< each realm's `apply_reflux_to_stage_forest` TBP, passing that realm's
   !< OWN final stage (`stages_per_step_forest()`) so the realm-side
   !< end-of-step gate (`stage == self%rk%nrk`) fires for every realm,
   !< asymmetric-K included. The realm-side body filters `flux_register%face(:)`
   !< by `face%coarse_realm == self%realm_index` and adds the per-cell
   !< end-of-step correction directly into its committed `self%q` (full
   !< `dt/dx` weight, no stage RK coefficient — AMReX Berger-Colella cadence).
   !<
   !< The forest never reaches `realm%rk`/`realm%q` directly: the per-cell
   !< write lives realm-side. This is the integrator-agnostic split — see
   !< [[realm_object]]%`apply_reflux_to_stage_forest` for the contract,
   !< and [[prism_cpu_object]] for the RK-specific override.
   !<
   !< Empty register fast path: when the register has no faces the per-realm
   !< calls each short-circuit on `flux_register%nfaces == 0`.
   class(forest_object), intent(in)    :: self     !< The forest (holds the flux register).
   class(realm_object),  intent(inout) :: realm(:) !< Realms; each realm's reflux TBP fires once.
   real(R8P),            intent(in)    :: dt       !< Time step.
   integer(I4P)                        :: is       !< Realm index.

   if (.not. self%flux_register%is_initialized_) return
   if (self%flux_register%nfaces == 0_I4P)        return
   if (.not. allocated(self%flux_register%face))  return

   ! Diagnostic: report the seam-flux mismatch the reflux is about to correct.
   ! For a same-resolution mirror seam this is round-off zero (no-op correction);
   ! for a true 2:1 AMR jump it is the non-trivial Berger-Colella delta. Emitted
   ! only when non-zero, so non-AMR runs stay silent. #13 §7.5 M3.
   block
      integer(I4P) :: ff
      real(R8P)    :: dmax
      dmax = 0._R8P
      do ff = 1_I4P, self%flux_register%nfaces
         if (allocated(self%flux_register%face(ff)%F_coarse) .and. &
             allocated(self%flux_register%face(ff)%F_fine_sum)) &
            dmax = max(dmax, maxval(abs(self%flux_register%face(ff)%F_coarse - &
                                        self%flux_register%face(ff)%F_fine_sum)))
      enddo
      if (dmax > 0._R8P) &
         call mpih%print_message('forest: reflux max|F_coarse-F_fine_sum| = '//trim(str(dmax)))
   endblock

   do is = 1_I4P, int(size(realm), I4P)
      call realm(is)%apply_reflux_to_stage_forest(stage=realm(is)%stages_per_step_forest(), dt=dt, &
                                                  flux_register=self%flux_register)
   enddo
   endsubroutine apply_reflux_corrections
endmodule adam_forest_object
