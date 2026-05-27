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
!< Class-with-TBPs (not module-of-routines) because Phase D (issue #10
!< Step 7) will need intrinsic-typed configuration on the forest itself —
!< MPI sub-communicator topology, inter-realm connectivity descriptor.
!< Adding that state later via a class extension is cheaper than via a
!< module API break.

use :: adam_realm_object,         only : realm_object
use :: adam_maps_object,          only : inter_realm_neighbor_t,                                                  &
                                         FACE_X_MAX, FACE_X_MIN, FACE_Y_MAX, FACE_Y_MIN, FACE_Z_MAX, FACE_Z_MIN, &
                                         face_axis_sign
use :: adam_forest_manifest,      only : forest_manifest_t, forest_face_pair_t
use :: adam_flux_register_object, only : flux_register_object, SEAM_KIND_INTER_REALM
use :: adam_parameters,           only : BC_SEAM, FEC_1_6_ARRAY
use :: adam_globals,              only : mpih
use :: mpi
use :: penf

implicit none
private
public :: forest_object

type :: forest_object
   !< Behavior-only orchestrator of an array of realms.
   integer(I4P)               :: n = 0_I4P     !< Number of realms in the forest (set by initialize from size(realm)).
   type(flux_register_object) :: flux_register !< Coarse-fine interface reflux machinery (Phase A of issue #13).
   contains
      ! public methods
      ! initialize/finalize
      procedure, pass(self) :: initialize               !< Sequence each realm's initialize_forest at startup (single shared INI).
      procedure, pass(self) :: initialize_from_manifest !< Like initialize, but each realm reads its own INI from a forest manifest.
      procedure, pass(self) :: finalize                 !< Sequence each realm's finalize_forest at shutdown.
      ! orchestrating methods
      procedure, pass(self) :: compute_global_dt      !< Min-reduce each realm's compute_local_dt_forest across the forest.
      procedure, pass(self) :: evolve_one_step        !< Iterate realm(:)%advance_one_step_forest(dt) for one global timestep.
      ! exchange_halos retired with `exchange_inter_realm_halos_forest` (#13);
      ! Phase 2 of evolve_one_step handles inter-realm seam refresh.
      procedure, pass(self) :: is_done                !< AND-reduce each realm's is_done_forest across the forest.
      procedure, pass(self) :: post_step              !< Iterate realm(:)%post_step_forest for the per-step diagnostics/IO block.
      procedure, pass(self) :: simulate               !< Main entry point (single shared INI): drive the full simulation.
      procedure, pass(self) :: simulate_from_manifest !< Main entry point (per-realm INIs via forest manifest).
      ! private methods
      procedure, pass(self), private :: populate_inter_realm_topology !< Translate manifest face-pairs into maps of neighbors.
      procedure, pass(self), private :: apply_reflux_corrections      !< Apply Berger-Colella reflux to coarse-side.
endtype forest_object

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
      call realm(is)%initialize_forest(filename=filename, realm_index=is, realms_number=self%n)
   enddo
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
      call realm(is)%initialize_forest(filename=trim(manifest%realm_ini(is)), realm_index=is, realms_number=self%n)
   enddo
   call self%populate_inter_realm_topology(realm, manifest)
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
   !< Advance every realm by one global timestep.
   !<
   !< N=1 fast path: realm owns the whole step (advance_one_step_forest).
   !<
   !< N>1 multi-realm path: forest drives an integrator-agnostic K-stage
   !< clock. `K_realm(is) = realm(is)%stages_per_step_forest()` is each
   !< realm's own stride — for RK realms `K_realm = rk%nrk`; other
   !< integrators report their own K. The forest takes
   !< `K_max = max(K_realm)` and gates each per-stage TBP behind
   !< `k <= K_realm(is)` so realms with a shorter stride participate
   !< only in the stages they have.
   !<
   !< Per-stage participation guard (Phase A v2): all realms MUST report
   !< the same K. The K_realm/K_max mechanics are in place so a future
   !< PR can lift the guard without touching the loop shape, but the
   !< temporal-coherence policy for sparse stages (what does a peer with
   !< `stage_active == K_realm < k` publish during stages k > K_realm?
   !< its committed `q`, or its last stage buffer?) is an open semantic
   !< decision and the safe Phase-A default is "no sparsity".
   !<
   !<   Phase 0 — `do is: open_step_forest(dt)`
   !<     Per-step prologue on each realm (external-field prelude,
   !<     integrator-private stage-buffer init, time bookkeeping).
   !<
   !<   For k = 1..K_max:
   !<
   !<     Phase 1 — `do is where k <= K_realm(is): begin_stage_forest(k)`
   !<       Each participating realm opens stage `k` on its integrator-
   !<       private buffer and sets `self%stage_active = k`. No ghost
   !<       reads, no peer-realm access.
   !<
   !<     Phase 2 — `do is, p: fill_seam_from_peer_forest`
   !<       Forest-driven inter-realm seam ghost fill. Each receiver realm
   !<       `is` iterates its peer realms `p`; for each peer the receiver
   !<       walks its own seam-map slice and copies peer-INTERIOR cells
   !<       (peer-side columns 2/4/5/6) into self's GHOST cells (recv-side
   !<       columns 3/7/8/9). Single TBP roundtrip per (receiver, peer);
   !<       same call site for CPU and FNL, backend dispatch via
   !<       `select type(peer)` inside the receiver's override. Buffer
   !<       selection on each side is driven by `stage_active`, so the
   !<       sparse-stage extension (peer with `K_realm < k`) is a no-op
   !<       once the guard above is lifted.
   !<
   !<     Phase 3 — `do is where k <= K_realm(is): end_stage_forest(k, flux_register)`
   !<       The race that originally required the begin/residuals/assign
   !<       three-phase split is GONE: by the time the stage-`k` buffer's
   !<       interior is overwritten, every receiver has already read its
   !<       peer's interior in Phase 2. Residuals + assignment can run in
   !<       a single sweep per realm, preserving today's reflux semantics
   !<       (flux_register accumulation happens inside
   !<       compute_residuals_fv_centered).
   !<
   !<   Phase 4 — `do is: close_step_forest`
   !<     Per-step epilogue (assembly of committed state from stage
   !<     buffers, BC, div-clean, time advance). Each realm clears its
   !<     `stage_active`.
   !<
   !< LOAD-BEARING INVARIANT: Phase 2 must complete on ALL realms before
   !< Phase 3 starts on ANY realm. If a future refactor interleaves Phase
   !< 2 and Phase 3 the read-after-overwrite race returns.
   class(forest_object), intent(inout)         :: self        !< The forest.
   class(realm_object),  intent(inout), target :: realm(:)    !< The realms to advance.
   real(R8P)                                   :: dt          !< Global timestep size.
   integer(I4P)                                :: is, p, k    !< Realm, peer, stage indices.
   integer(I4P), allocatable                   :: K_realm(:)  !< Per-realm stage counts (1..size(realm)).
   integer(I4P)                                :: K_max       !< Forest-wide max stage count (multi-realm path).
   integer(I4P)                                :: peer_idx    !< Peer realm index for a seam exchange.

   call self%flux_register%reset
   call self%compute_global_dt(realm=realm, dt=dt)
   if (int(size(realm), I4P) == 1_I4P) then
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
      ! Phase A v2 guard: sparse stages (K_realm /= K_max) are structurally
      ! supported by the loop below but semantically not yet defined; reject.
      do is = 1_I4P, int(size(realm), I4P)
         if (K_realm(is) /= K_max) &
            call mpih%error_stop(msg='forest_object%evolve_one_step: realms disagree on stages_per_step_forest')
      enddo
      do k = 1_I4P, K_max
         ! Phase 1 — open stage k on each participating realm; each realm sets stage_active=k.
         do is = 1_I4P, int(size(realm), I4P)
            if (k > K_realm(is)) cycle
            call realm(is)%begin_stage_forest(k=k, K_total=K_max, dt=dt)
         enddo

         ! Phase 2 — forest-driven inter-realm seam ghost fill.
         !
         ! Backend-agnostic dispatch via a single TBP on the RECEIVER:
         ! `realm(is)%fill_seam_from_peer_forest(peer=realm(peer_idx), p_idx=p)`.
         ! The receiver walks its own seam map (host on CPU, device on FNL),
         ! `select type`s the peer to access its active buffer, and copies
         ! peer-INTERIOR cells into self's GHOST cells. Single TBP, same
         ! call site for CPU and FNL; backend differences confined to one
         ! method per backend. Runs unconditionally regardless of K_realm —
         ! a stalled realm (K_realm < k) keeps `stage_active = K_realm`
         ! until close_step, so peers naturally pick up its last stage
         ! buffer; if the receiver itself is stalled it writes into its
         ! own last-stage buffer, which `close_step_forest` will overwrite.
         do is = 1_I4P, int(size(realm), I4P)
            if (.not. allocated(realm(is)%adam%maps%seam_local_map_ghost_cell)) cycle
            ! Phase-A guard: cross-rank seam path is not implemented.
            if (allocated(realm(is)%adam%maps%seam_comm_map_send_ghost_cell)) &
               call mpih%error_stop(msg='evolve_one_step: cross-rank seam not implemented (Phase B)')
            do p = 1_I4P, int(size(realm(is)%adam%maps%seam_local_peer_realm), I4P)
               peer_idx = realm(is)%adam%maps%seam_local_peer_realm(p)
               call realm(is)%fill_seam_from_peer_forest(peer=realm(peer_idx), p_idx=p)
            enddo
         enddo

         ! Phase 3 — residuals (with flux_register) + assign on each participating realm.
         ! Race-free: Phase 2 has already filled seam ghosts, and the
         ! intra-realm residuals→assign ordering is the integrator's normal
         ! flow within a single realm.
         do is = 1_I4P, int(size(realm), I4P)
            if (k > K_realm(is)) cycle
            call realm(is)%end_stage_forest(k=k, K_total=K_max, dt=dt, flux_register=self%flux_register)
         enddo

         ! Reflux corrections — preserved placement; integrator-weighted inside.
         call self%flux_register%reduce_fine_sums
         call self%apply_reflux_corrections(realm=realm, stage=k, dt=dt)
      enddo
      do is = 1_I4P, int(size(realm), I4P)
         call realm(is)%close_step_forest(dt=dt)
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
   type(forest_face_pair_t)                :: pair                !< Loop alias.

   associate(self_unused => self) ! method takes self for TBP-dispatch symmetry; uses no forest state
   end associate
   if (.not. allocated(manifest%face_pairs)) return  ! no inter-realm topology declared

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
   call register_inter_realm_seams(realm=realm, manifest=manifest, flux_register=self%flux_register)
   ! Build the per-cell inter-realm ghost map. Per-realm: enumerate every ghost cell in self's seam-block
   ! ghost region and resolve the (peer_realm, peer_block, peer_interior_cell)
   ! tuple. The runtime exchange then becomes a flat indexed loop, replacing
   ! the per-stage geometric find_peer_block + face-slab copy that misses
   ! corner / edge ghosts.
   call build_inter_realm_ghost_cell_map(realm=realm, manifest=manifest)
   ! Derive the sorted-by-peer seam_local map + per-peer index arrays + per-peer
   ! pack/unpack buffers from the just-built inter_realm_ghost_cell map. The new
   ! arrays are what the forest's Phase 2 seam fill (in evolve_one_step) consumes;
   ! inter_realm_ghost_cell is kept during migration so the legacy
   ! exchange_inter_realm_halos_forest TBP can be deleted in a single pass.
   ! Every entry is required to be same-rank under Phase A (replicated forest);
   ! cross-rank entries error_stop here to flag the unimplemented Phase B path.
   call build_seam_local_map(realm=realm)
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

      subroutine register_inter_realm_seams(realm, manifest, flux_register)
      !< Populate `flux_register` from the manifest face-pairs.
      !<
      !< Two-pass algorithm:
      !<
      !<   * Pass 1 — count the total number of register entries:
      !<     one entry per (face-pair, block-on-realm_a-side-of-the-seam).
      !<     The seam-block enumeration uses a geometric test against the
      !<     realm's domain extent (`grid%domain_emin/emax`) and the block's
      !<     face position (`field%emin/emax(axis, block)`); a block lies on
      !<     the realm's `face_a` iff its face coordinate matches the realm
      !<     boundary to within a small absolute tolerance.
      !<
      !<   * Pass 2 — call `register_face` once per (face-pair, seam-block),
      !<     supplying `coarse_realm = pair%realm_a`, `fine_realm =
      !<     pair%realm_b`, `nface_cells` from the realm-a side's grid (the
      !<     two tangential axes' cell counts), `nv` from the realm-a side's
      !<     `field%nv`, and the per-stage register depth from realm-a's
      !<     `stages_per_step_forest()` — the integrator-agnostic contract
      !<     TBP, NOT a direct `rk%nrk` reach.
      !<
      !< For the current same-resolution (COUPLING_MIRROR) case, realm_a and
      !< realm_b carry identical `nv` / stage-count / `nface_cells` by
      !< construction (validated upstream by the manifest's structural
      !< checks and by the `K_realm` equality guard in
      !< `forest_object%evolve_one_step`). A coarse-fine AMR case would
      !< resolve `fine_block(:)` by enumerating the realm_b-side blocks
      !< that geometrically cover the realm_a-side block face; for
      !< same-resolution that resolves to a single fine block, populated
      !< by `find_peer_block` at exchange time.
      !<
      !< The keyword arg on `flux_register%register_face` is `n_stages`
      !< (renamed from the original `nrk` together with `accumulate_*_flux`'s
      !< `stage`); both sides of the wiring now speak integrator-neutral
      !< vocabulary end-to-end.
      class(realm_object),        intent(inout) :: realm(:)       !< Initialized realms.
      type(forest_manifest_t),    intent(in)    :: manifest       !< Parsed manifest.
      type(flux_register_object), intent(inout) :: flux_register  !< Berger-Colella reflux accumulator owned by the forest.
      integer(I4P)                              :: f, b           !< Face-pair, block counters.
      integer(I4P)                              :: a_realm        !< Coarse-side realm index alias.
      integer(I4P)                              :: a_axis, a_sign !< Coarse-face axis and sign.
      integer(I4P)                              :: nfaces_total   !< Total register entries.
      integer(I4P)                              :: cursor         !< Write cursor into the register.
      integer(I4P)                              :: nface_cells    !< Cell count on the coarse-face skin.
      type(forest_face_pair_t)                  :: pair           !< Manifest face-pair alias.

      if (.not. allocated(manifest%face_pairs)) then
         ! No inter-realm topology — initialize with zero faces so the
         ! register's `is_initialized_` flag flips and the per-step `reset`
         ! call becomes a safe no-op on the empty register.
         call flux_register%initialize(nfaces=0_I4P)
         return
      endif

      ! Pass 1: count register entries.
      nfaces_total = 0_I4P
      do f = 1_I4P, int(size(manifest%face_pairs), I4P)
         pair    = manifest%face_pairs(f)
         a_realm = pair%realm_a
         call face_axis_sign(pair%face_a, a_axis, a_sign)
         do b = 1_I4P, int(realm(a_realm)%adam%field%blocks_number, I4P)
            if (block_face_on_realm_boundary(realm(a_realm), b, a_axis, a_sign)) &
               nfaces_total = nfaces_total + 1_I4P
         enddo
      enddo

      call flux_register%initialize(nfaces=nfaces_total)
      if (nfaces_total == 0_I4P) return

      ! Allocate the per-realm (block, face_1_6) → register_index lookup.
      ! Sized (nb, 6); zero means "not a seam face", positive = 1-based
      ! index into flux_register%face(:). Consumed by PRISM's
      ! compute_residuals_fv_centered to know where to accumulate fluxes.
      do is = 1_I4P, int(size(realm), I4P)
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

      ! Pass 2: register one entry per (face-pair, seam-block-on-coarse-side),
      ! AND populate both the coarse-side and fine-side
      ! inter_realm_face_register_index lookups so each realm's residual
      ! routine can find the right entry in O(1).
      cursor = 0_I4P
      do f = 1_I4P, int(size(manifest%face_pairs), I4P)
         pair    = manifest%face_pairs(f)
         a_realm = pair%realm_a
         call face_axis_sign(pair%face_a, a_axis, a_sign)
         ! Coarse-face skin cell count: product of the two tangential cell counts.
         ! For axis=x (a_axis=1) the tangential axes are y and z; etc.
         nface_cells = tangential_cell_count(realm(a_realm), a_axis)
         do b = 1_I4P, int(realm(a_realm)%adam%field%blocks_number, I4P)
            if (.not. block_face_on_realm_boundary(realm(a_realm), b, a_axis, a_sign)) cycle
            cursor = cursor + 1_I4P
            ! fine_block(:) is left empty (size 0) here — Phase A v1 uses
            ! same-resolution mirror seams where the fine-side block is
            ! resolved via the per-cell ghost map (inter_realm_ghost_cell).
            ! True AMR coarse-fine adjacency would populate fine_block(:)
            ! with the 4 (in 3D) finer blocks covering this coarse face.
            ! `fine_block` deliberately omitted: same-resolution mirror seam
            ! has no fine-side blocks to record. Passing an empty array literal
            ! `[integer(I4P) ::]` here was the original Phase A shortcut; it
            ! poisoned the descriptor copy on nvfortran 26.x. Absent ↔ unallocated
            ! `slot%fine_block`, which every consumer must guard with allocated().
            call flux_register%register_face(face_index=cursor,                                 &
                                             seam_kind=SEAM_KIND_INTER_REALM,                   &
                                             coarse_realm=pair%realm_a,                         &
                                             coarse_block=b,                                    &
                                             coarse_face=pair%face_a,                           &
                                             fine_realm=pair%realm_b,                           &
                                             nface_cells=nface_cells,                           &
                                             nv=int(realm(a_realm)%adam%field%nv, I4P),         &
                                             n_stages=realm(a_realm)%stages_per_step_forest())

            ! Index lookups: coarse side knows its own (b, face_a → fec)
            ! immediately. Fine side requires a geometric peer lookup to
            ! find the matching realm_b block; same-resolution mirror →
            ! identical (tangential extents, opposite-axis face) match.
            !
            ! Sign convention (consumed by PRISM's FV reflux hook):
            !   +cursor stored on the coarse-side realm (pair%realm_a)
            !   -cursor stored on the fine-side   realm (pair%realm_b)
            !       0   stored anywhere = "not a seam face".
            ! The consumer recovers the register index via `abs()` and
            ! the coarse/fine role via the sign — no realm-identity field
            ! is needed on the realm object.
            block
               integer(I4P) :: bc_fec_a, bc_fec_b
               integer(I4P) :: b_peer, b_peer_axis, b_peer_sign
               bc_fec_a = face_code_to_bc_fec(pair%face_a)
               bc_fec_b = face_code_to_bc_fec(pair%face_b)
               if (bc_fec_a > 0_I4P .and. bc_fec_a <= 6_I4P) &
                  realm(a_realm)%adam%maps%inter_realm_face_register_index(b, bc_fec_a) = +cursor
               call face_axis_sign(pair%face_b, b_peer_axis, b_peer_sign)
               call find_seam_peer_block(realm, my_realm_idx=pair%realm_a, my_block=b, &
                                         my_axis=a_axis, my_sign=a_sign,               &
                                         peer_realm_idx=pair%realm_b,                  &
                                         peer_axis=b_peer_axis, peer_sign=b_peer_sign, &
                                         b_peer=b_peer)
               if (b_peer > 0_I4P .and. bc_fec_b > 0_I4P .and. bc_fec_b <= 6_I4P) then
                  if (allocated(realm(pair%realm_b)%adam%maps%inter_realm_face_register_index)) &
                     realm(pair%realm_b)%adam%maps%inter_realm_face_register_index(b_peer, bc_fec_b) = -cursor
               endif
            endblock
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

      subroutine find_seam_peer_block(realm, my_realm_idx, my_block, my_axis, my_sign, &
                                      peer_realm_idx, peer_axis, peer_sign, b_peer)
      !< Find the peer-realm block whose face geometrically matches the
      !< given (my_realm, my_block, my_face) tuple — same-resolution mirror
      !< version. Returns `b_peer = 0` if no match (e.g. peer lives on
      !< another rank).
      !<
      !< Matching criteria (same-resolution mirror, COUPLING_MIRROR):
      !<   * peer block's face-coordinate along `peer_axis` equals my
      !<     face-coordinate along `my_axis` (the two coupled faces lie
      !<     in the same plane);
      !<   * peer block's extents in the two tangential axes equal mine.
      !<
      !< This mirrors the now-obsolete `find_peer_block` helper that lived
      !< inside `prism_cpu_object%exchange_inter_realm_halos_forest` before
      !< the seam comm-map commit; the canonical version lives here in
      !< `forest_object`'s contains scope and is used by the register
      !< topology pass to populate `inter_realm_face_register_index`.
      class(realm_object), intent(in)  :: realm(:)
      integer(I4P),        intent(in)  :: my_realm_idx, my_block
      integer(I4P),        intent(in)  :: my_axis, my_sign
      integer(I4P),        intent(in)  :: peer_realm_idx
      integer(I4P),        intent(in)  :: peer_axis, peer_sign
      integer(I4P),        intent(out) :: b_peer
      integer(I4P)                     :: bp
      integer(I4P)                     :: tax1, tax2
      real(R8P)                        :: my_face_coord, my_tmin(2), my_tmax(2)
      real(R8P)                        :: peer_face_coord
      real(R8P)                        :: tol

      b_peer = 0_I4P
      if (peer_axis /= my_axis) return
      if (peer_sign == my_sign) return  ! must be opposite (MIN ↔ MAX)
      tax1 = mod(my_axis,        3_I4P) + 1_I4P
      tax2 = mod(my_axis + 1_I4P, 3_I4P) + 1_I4P
      if (my_sign > 0_I4P) then
         my_face_coord = realm(my_realm_idx)%adam%field%emax(my_axis, my_block)
      else
         my_face_coord = realm(my_realm_idx)%adam%field%emin(my_axis, my_block)
      endif
      my_tmin(1) = realm(my_realm_idx)%adam%field%emin(tax1, my_block)
      my_tmax(1) = realm(my_realm_idx)%adam%field%emax(tax1, my_block)
      my_tmin(2) = realm(my_realm_idx)%adam%field%emin(tax2, my_block)
      my_tmax(2) = realm(my_realm_idx)%adam%field%emax(tax2, my_block)
      tol = max(abs(my_face_coord), 1._R8P) * 1.0e-10_R8P

      do bp = 1_I4P, int(realm(peer_realm_idx)%adam%field%blocks_number, I4P)
         if (peer_sign > 0_I4P) then
            peer_face_coord = realm(peer_realm_idx)%adam%field%emax(peer_axis, bp)
         else
            peer_face_coord = realm(peer_realm_idx)%adam%field%emin(peer_axis, bp)
         endif
         if (abs(peer_face_coord - my_face_coord) > tol) cycle
         if (abs(realm(peer_realm_idx)%adam%field%emin(tax1, bp) - my_tmin(1)) > tol) cycle
         if (abs(realm(peer_realm_idx)%adam%field%emax(tax1, bp) - my_tmax(1)) > tol) cycle
         if (abs(realm(peer_realm_idx)%adam%field%emin(tax2, bp) - my_tmin(2)) > tol) cycle
         if (abs(realm(peer_realm_idx)%adam%field%emax(tax2, bp) - my_tmax(2)) > tol) cycle
         b_peer = bp
         return
      enddo
      endsubroutine find_seam_peer_block

      subroutine build_inter_realm_ghost_cell_map(realm, manifest)
      !< Populate each realm's `adam%maps%inter_realm_ghost_cell` map.
      !<
      !< This is the per-cell seam ghost map driving
      !< `exchange_inter_realm_halos_forest` (Phase A of issue #13).
      !< Algorithm, per realm:
      !<
      !<   * Two-pass over (face_pair, seam_block, ghost_cell).
      !<   * Pass 1 counts entries: for every face pair whose `realm_a`
      !<     matches this realm (and symmetrically `realm_b`), every block
      !<     of this realm lying on `face_a`/`face_b`, every ghost cell in
      !<     the full ghost slab on that face (INCLUDING tangential corner
      !<     and edge ghosts — this is the defect-B fix).
      !<   * Pass 2 populates: for each ghost cell, find the peer-realm
      !<     block whose INTERIOR contains the global coordinates that
      !<     correspond to this ghost cell. Geometric match in 3D against
      !<     `field%emin/emax(:, bp)`. If no peer block matches, the
      !<     ghost is in the physical exterior of the peer realm (a
      !<     corner where the seam meets a physical boundary) and is
      !<     skipped — `set_boundary_conditions` on self fills it.
      !<
      !< Coordinate convention: a ghost cell at local indices
      !< (i_g, j_g, k_g) of block b in this realm has cell-center
      !< coordinates
      !<     x = field%emin(d, b) + (idx - 0.5) * field%dxyz(d, b)
      !< where (idx, d) ranges over (i_g, 1), (j_g, 2), (k_g, 3). The
      !< match in the peer realm finds a peer block bp such that for
      !< each axis d:
      !<     emin_peer(d, bp) < x_d < emax_peer(d, bp)
      !< (open interval — we want INTERIOR cells, not on-the-boundary
      !< ones, since cell centers are strictly inside their cells). The
      !< peer local cell index is:
      !<     i_peer = nint((x - emin_peer(d, bp)) / dxyz_peer(d, bp) + 0.5)
      !<
      !< For the rmf-2realm same-resolution case `field%dxyz` is uniform
      !< across both realms, so the match is integer-clean. For the
      !< coarse-fine AMR case (not exercised in Phase A v1) the per-cell
      !< resolution differs and the same geometric match degrades to a
      !< nearest-cell mapping — the `one_or_eight` column reserves the
      !< value 8 for that future case; v1 always writes 1.
      class(realm_object),     intent(inout) :: realm(:)  !< Forest realms (their maps get populated).
      type(forest_manifest_t), intent(in)    :: manifest  !< Parsed manifest.
      integer(I4P)                           :: f, is, ip !< Counters: face-pair, self-realm, peer-realm.
      integer(I4P)                           :: my_face
      integer(I4P)                           :: my_realm_idx
      integer(I4P)                           :: my_axis
      integer(I4P)                           :: my_sign
      integer(I4P)                           :: peer_face
      integer(I4P)                           :: peer_realm_idx
      integer(I4P)                           :: count_total
      integer(I4P), allocatable              :: per_realm_count(:)
      type(forest_face_pair_t)               :: pair

      if (.not. allocated(manifest%face_pairs)) then
         do is = 1_I4P, int(size(realm), I4P)
            if (allocated(realm(is)%adam%maps%inter_realm_ghost_cell)) &
               deallocate(realm(is)%adam%maps%inter_realm_ghost_cell)
         enddo
         return
      endif

      allocate(per_realm_count(size(realm)))
      per_realm_count = 0_I4P

      ! Pass 1: per-realm row counts
      do f = 1_I4P, int(size(manifest%face_pairs), I4P)
         pair = manifest%face_pairs(f)
         ! Side A: self = realm_a, peer = realm_b, looking through face_a.
         call count_seam_ghost_cells(realm=realm, per_realm_count=per_realm_count, &
                                     my_realm_idx=pair%realm_a, peer_realm_idx=pair%realm_b, &
                                     my_face=pair%face_a)
         ! Side B: symmetric pair entry: self = realm_b, peer = realm_a.
         call count_seam_ghost_cells(realm=realm, per_realm_count=per_realm_count, &
                                     my_realm_idx=pair%realm_b, peer_realm_idx=pair%realm_a, &
                                     my_face=pair%face_b)
      enddo

      count_total = 0_I4P
      do is = 1_I4P, int(size(realm), I4P)
         count_total = count_total + per_realm_count(is)
      enddo

      ! Allocate per-realm maps
      do is = 1_I4P, int(size(realm), I4P)
         if (allocated(realm(is)%adam%maps%inter_realm_ghost_cell)) &
            deallocate(realm(is)%adam%maps%inter_realm_ghost_cell)
         if (per_realm_count(is) > 0_I4P) &
            allocate(realm(is)%adam%maps%inter_realm_ghost_cell(1:per_realm_count(is), 1:10))
      enddo

      ! Pass 2: populate
      ! Per-realm cursor — re-use per_realm_count as the running write index,
      ! but reset to zero before reuse.
      per_realm_count = 0_I4P
      do f = 1_I4P, int(size(manifest%face_pairs), I4P)
         pair = manifest%face_pairs(f)
         call populate_seam_ghost_cells(realm=realm, per_realm_count=per_realm_count, &
                                        my_realm_idx=pair%realm_a, peer_realm_idx=pair%realm_b, &
                                        my_face=pair%face_a)
         call populate_seam_ghost_cells(realm=realm, per_realm_count=per_realm_count, &
                                        my_realm_idx=pair%realm_b, peer_realm_idx=pair%realm_a, &
                                        my_face=pair%face_b)
      enddo
      endsubroutine build_inter_realm_ghost_cell_map

      subroutine build_seam_local_map(realm)
      !< Derive `seam_local_map_ghost_cell` (sorted by `peer_realm`) +
      !< per-peer index arrays + per-peer pack/unpack buffers from the
      !< already-populated `inter_realm_ghost_cell` map.
      !<
      !< The output is what the forest's Phase 2 seam fill consumes via the
      !< realm `pack_seam_cells` / `unpack_seam_cells` TBPs. Rows are sorted
      !< by `peer_realm` so the forest can extract per-peer row ranges in
      !< O(1) via the index arrays.
      !<
      !< Invariant (Phase A): every entry must be same-rank — under the
      !< replicated-forest layout (`rmf-2realm`) both ranks own both realms,
      !< so the rank that owns `b_send` in the peer realm equals
      !< `mpih%myrank`. Cross-rank entries are detected via the peer's
      !< `comm_map_recv` (which lists who owns each block of the peer
      !< realm); a single cross-rank entry triggers an `error_stop` flagging
      !< the unimplemented Phase B `update_ghost_seam_mpi` path.
      class(realm_object), intent(inout) :: realm(:)
      integer(I4P)                       :: is, c, nrows, n_peers, p, peer
      integer(I4P)                       :: nv_is, max_rows_per_peer
      integer(I4P), allocatable          :: peer_list(:), peer_count(:), peer_cursor(:)
      integer(I4P)                       :: write_idx
      integer(I4P)                       :: src_col_peer  ! col 1 of inter_realm_ghost_cell

      src_col_peer = 1_I4P

      do is = 1_I4P, int(size(realm), I4P)
         associate(maps => realm(is)%adam%maps)
         ! Deallocate prior state.
         if (allocated(maps%seam_local_map_ghost_cell)) deallocate(maps%seam_local_map_ghost_cell)
         if (allocated(maps%seam_local_peer_realm))     deallocate(maps%seam_local_peer_realm)
         if (allocated(maps%seam_local_peer_row_start)) deallocate(maps%seam_local_peer_row_start)
         if (allocated(maps%seam_local_peer_row_count)) deallocate(maps%seam_local_peer_row_count)
         if (allocated(maps%seam_local_send_buf))       deallocate(maps%seam_local_send_buf)
         if (allocated(maps%seam_local_recv_buf))       deallocate(maps%seam_local_recv_buf)

         if (.not. allocated(maps%inter_realm_ghost_cell)) cycle

         nrows = int(size(maps%inter_realm_ghost_cell, dim=1), I4P)
         if (nrows == 0_I4P) cycle

         ! Pass 1: enumerate distinct peer realms appearing in col 1.
         ! Realm count is small (handful); a linear-scan dedup is cheap.
         allocate(peer_list(nrows))   ; peer_list = 0_I4P
         allocate(peer_count(nrows))  ; peer_count = 0_I4P
         n_peers = 0_I4P
         do c = 1_I4P, nrows
            peer = int(maps%inter_realm_ghost_cell(c, src_col_peer), I4P)
            p = 0_I4P
            do p = 1_I4P, n_peers
               if (peer_list(p) == peer) exit
            enddo
            if (p > n_peers) then
               n_peers = n_peers + 1_I4P
               peer_list(n_peers) = peer
               peer_count(n_peers) = 1_I4P
            else
               peer_count(p) = peer_count(p) + 1_I4P
            endif
         enddo

         ! Allocate index arrays sized to actual peer count.
         allocate(maps%seam_local_peer_realm(n_peers))
         allocate(maps%seam_local_peer_row_start(n_peers))
         allocate(maps%seam_local_peer_row_count(n_peers))
         maps%seam_local_peer_realm     = peer_list(1:n_peers)
         maps%seam_local_peer_row_count = peer_count(1:n_peers)

         maps%seam_local_peer_row_start(1) = 1_I4P
         do p = 2_I4P, n_peers
            maps%seam_local_peer_row_start(p) =                  &
               maps%seam_local_peer_row_start(p - 1_I4P)       + &
               maps%seam_local_peer_row_count(p - 1_I4P)
         enddo

         ! Allocate the sorted map (same row count, 9 cols — drop the
         ! one_or_eight column reserved for AMR coarse-fine).
         allocate(maps%seam_local_map_ghost_cell(nrows, 9))

         ! Pass 2: fill sorted map using per-peer write cursors.
         allocate(peer_cursor(n_peers))
         peer_cursor = maps%seam_local_peer_row_start
         do c = 1_I4P, nrows
            peer = int(maps%inter_realm_ghost_cell(c, src_col_peer), I4P)
            do p = 1_I4P, n_peers
               if (peer_list(p) == peer) exit
            enddo
            write_idx = peer_cursor(p)
            peer_cursor(p) = peer_cursor(p) + 1_I4P
            ! Columns 1..9 = [peer_realm, b_send, b_recv, i_send, j_send, k_send, i_recv, j_recv, k_recv].
            ! Source columns are the same 1..9; col 10 (one_or_eight) is intentionally dropped.
            maps%seam_local_map_ghost_cell(write_idx, 1:9) = &
               int(maps%inter_realm_ghost_cell(c, 1:9), I4P)
         enddo
         deallocate(peer_cursor)

         ! Phase A invariant: every entry is same-rank. Under the replicated-
         ! forest layout (rmf-2realm) all ranks own all realms' blocks, so the
         ! invariant holds by construction. The defect-A cross-rank case lands
         ! with update_ghost_seam_mpi (Phase B); until then the forest Phase 2
         ! loop error_stops if seam_comm_map_send_ghost_cell becomes allocated.
         !
         ! Sizing for per-peer buffers: nv × max(row_count_per_peer), one column per peer.
         ! Realm exposes nv as a pointer component (initialize binds self%nv => adam%field%nv).
         nv_is = realm(is)%nv
         max_rows_per_peer = maxval(maps%seam_local_peer_row_count)
         allocate(maps%seam_local_send_buf(nv_is * max_rows_per_peer, n_peers))
         allocate(maps%seam_local_recv_buf(nv_is * max_rows_per_peer, n_peers))
         maps%seam_local_send_buf = 0.0_R8P
         maps%seam_local_recv_buf = 0.0_R8P

         deallocate(peer_list)
         deallocate(peer_count)
         endassociate
      enddo
      endsubroutine build_seam_local_map

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

      subroutine count_seam_ghost_cells(realm, per_realm_count, my_realm_idx, peer_realm_idx, my_face)
      !< First-pass row counter for the inter-realm ghost-cell map.
      class(realm_object), intent(in)    :: realm(:)
      integer(I4P),        intent(inout) :: per_realm_count(:)
      integer(I4P),        intent(in)    :: my_realm_idx, peer_realm_idx, my_face
      integer(I4P) :: b, axis, sgn, i_g, j_g, k_g
      integer(I4P) :: imin_g, imax_g, jmin_g, jmax_g, kmin_g, kmax_g
      real(R8P)    :: xg(3)
      integer(I4P) :: bp, ip_dummy, jp_dummy, kp_dummy

      call face_axis_sign(my_face, axis, sgn)
      do b = 1_I4P, int(realm(my_realm_idx)%adam%field%blocks_number, I4P)
         if (.not. block_face_on_realm_boundary(realm(my_realm_idx), b, axis, sgn)) cycle
         call ghost_slab_extents(realm(my_realm_idx), axis, sgn, imin_g, imax_g, jmin_g, jmax_g, kmin_g, kmax_g)
         do k_g = kmin_g, kmax_g
            do j_g = jmin_g, jmax_g
               do i_g = imin_g, imax_g
                  call ghost_cell_center(realm(my_realm_idx), b, i_g, j_g, k_g, xg)
                  call find_peer_cell(realm(peer_realm_idx), xg, bp, ip_dummy, jp_dummy, kp_dummy)
                  if (bp > 0_I4P) per_realm_count(my_realm_idx) = per_realm_count(my_realm_idx) + 1_I4P
               enddo
            enddo
         enddo
      enddo
      endsubroutine count_seam_ghost_cells

      subroutine populate_seam_ghost_cells(realm, per_realm_count, my_realm_idx, peer_realm_idx, my_face)
      !< Second-pass row populator for the inter-realm ghost-cell map.
      class(realm_object), intent(inout) :: realm(:)
      integer(I4P),        intent(inout) :: per_realm_count(:)
      integer(I4P),        intent(in)    :: my_realm_idx, peer_realm_idx, my_face
      integer(I4P) :: b, axis, sgn, i_g, j_g, k_g
      integer(I4P) :: imin_g, imax_g, jmin_g, jmax_g, kmin_g, kmax_g
      real(R8P)    :: xg(3)
      integer(I4P) :: bp, ip, jp, kp, c

      call face_axis_sign(my_face, axis, sgn)
      do b = 1_I4P, int(realm(my_realm_idx)%adam%field%blocks_number, I4P)
         if (.not. block_face_on_realm_boundary(realm(my_realm_idx), b, axis, sgn)) cycle
         call ghost_slab_extents(realm(my_realm_idx), axis, sgn, imin_g, imax_g, jmin_g, jmax_g, kmin_g, kmax_g)
         do k_g = kmin_g, kmax_g
            do j_g = jmin_g, jmax_g
               do i_g = imin_g, imax_g
                  call ghost_cell_center(realm(my_realm_idx), b, i_g, j_g, k_g, xg)
                  call find_peer_cell(realm(peer_realm_idx), xg, bp, ip, jp, kp)
                  if (bp <= 0_I4P) cycle
                  per_realm_count(my_realm_idx) = per_realm_count(my_realm_idx) + 1_I4P
                  c = per_realm_count(my_realm_idx)
                  realm(my_realm_idx)%adam%maps%inter_realm_ghost_cell(c, 1)  = peer_realm_idx
                  realm(my_realm_idx)%adam%maps%inter_realm_ghost_cell(c, 2)  = bp
                  realm(my_realm_idx)%adam%maps%inter_realm_ghost_cell(c, 3)  = b
                  realm(my_realm_idx)%adam%maps%inter_realm_ghost_cell(c, 4)  = ip
                  realm(my_realm_idx)%adam%maps%inter_realm_ghost_cell(c, 5)  = jp
                  realm(my_realm_idx)%adam%maps%inter_realm_ghost_cell(c, 6)  = kp
                  realm(my_realm_idx)%adam%maps%inter_realm_ghost_cell(c, 7)  = i_g
                  realm(my_realm_idx)%adam%maps%inter_realm_ghost_cell(c, 8)  = j_g
                  realm(my_realm_idx)%adam%maps%inter_realm_ghost_cell(c, 9)  = k_g
                  realm(my_realm_idx)%adam%maps%inter_realm_ghost_cell(c, 10) = 1_I4P  ! same-resolution mirror
               enddo
            enddo
         enddo
      enddo
      endsubroutine populate_seam_ghost_cells

      pure subroutine ghost_slab_extents(this_realm, axis, sgn, imin_g, imax_g, jmin_g, jmax_g, kmin_g, kmax_g)
         !< Compute the (i,j,k) loop bounds for a block's ghost slab on one face.
         !<
         !< The slab is `ngc` cells deep on the (axis, sgn) face. The two
         !< tangential axes run the FULL ghost-extended range
         !< `[1-ngc..n+ngc]` — this is the defect-B fix.
         class(realm_object), intent(in)  :: this_realm
         integer(I4P),        intent(in)  :: axis, sgn
         integer(I4P),        intent(out) :: imin_g, imax_g, jmin_g, jmax_g, kmin_g, kmax_g

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
         end select
         end associate
         endsubroutine ghost_slab_extents

      pure subroutine ghost_cell_center(this_realm, b, i, j, k, xc)
      !< Cell-center coordinates of (i, j, k) in block b, regardless of
      !< whether (i, j, k) is an interior or ghost cell — same formula.
      class(realm_object), intent(in)  :: this_realm
      integer(I4P),        intent(in)  :: b, i, j, k
      real(R8P),           intent(out) :: xc(3)

      associate(field => this_realm%adam%field)
      xc(1) = field%emin(1, b) + (real(i, R8P) - 0.5_R8P) * field%dxyz(1, b)
      xc(2) = field%emin(2, b) + (real(j, R8P) - 0.5_R8P) * field%dxyz(2, b)
      xc(3) = field%emin(3, b) + (real(k, R8P) - 0.5_R8P) * field%dxyz(3, b)
      end associate
      endsubroutine ghost_cell_center

      subroutine find_peer_cell(peer_realm, xc, bp, ip, jp, kp)
      !< Find the peer-realm block + interior cell whose cell-center
      !< coincides (within tolerance) with the global point `xc`.
      !<
      !< Returns `bp = 0` if no peer block contains the point — that
      !< means the ghost cell maps outside the peer's physical extent
      !< (typical at corners where the seam meets a physical boundary)
      !< and the consumer should skip this entry (the physical BC on
      !< self will fill the ghost).
      class(realm_object), intent(in)  :: peer_realm
      real(R8P),           intent(in)  :: xc(3)
      integer(I4P),        intent(out) :: bp, ip, jp, kp
      integer(I4P) :: b
      real(R8P)    :: tol(3)
      real(R8P)    :: emin_b(3), emax_b(3), dxyz_b(3)
      integer(I4P) :: ii, jj, kk
      logical      :: inside_block

      bp = 0_I4P ; ip = 0_I4P ; jp = 0_I4P ; kp = 0_I4P
      do b = 1_I4P, int(peer_realm%adam%field%blocks_number, I4P)
         emin_b = peer_realm%adam%field%emin(:, b)
         emax_b = peer_realm%adam%field%emax(:, b)
         dxyz_b = peer_realm%adam%field%dxyz(:, b)
         tol = max(abs(emin_b), abs(emax_b), 1._R8P) * 1.0e-10_R8P
         inside_block = .true.
         if (xc(1) < emin_b(1) - tol(1) .or. xc(1) > emax_b(1) + tol(1)) inside_block = .false.
         if (xc(2) < emin_b(2) - tol(2) .or. xc(2) > emax_b(2) + tol(2)) inside_block = .false.
         if (xc(3) < emin_b(3) - tol(3) .or. xc(3) > emax_b(3) + tol(3)) inside_block = .false.
         if (.not. inside_block) cycle
         ! Block found — translate coordinates to local cell index.
         ii = nint((xc(1) - emin_b(1)) / dxyz_b(1) + 0.5_R8P, I4P)
         jj = nint((xc(2) - emin_b(2)) / dxyz_b(2) + 0.5_R8P, I4P)
         kk = nint((xc(3) - emin_b(3)) / dxyz_b(3) + 0.5_R8P, I4P)
         ! Reject if rounding hit an out-of-interior index.
         if (ii < 1_I4P .or. ii > peer_realm%adam%grid%ni) cycle
         if (jj < 1_I4P .or. jj > peer_realm%adam%grid%nj) cycle
         if (kk < 1_I4P .or. kk > peer_realm%adam%grid%nk) cycle
         bp = b ; ip = ii ; jp = jj ; kp = kk
         return
      enddo
      endsubroutine find_peer_cell
   endsubroutine populate_inter_realm_topology

   subroutine apply_reflux_corrections(self, realm, stage, dt)
   !< Dispatch the Berger-Colella reflux correction to every realm.
   !<
   !< The forest's role here is purely an orchestrator: it iterates the
   !< realm array and invokes each realm's `apply_reflux_to_stage_forest`
   !< TBP, passing the realm's 1-based forest index as `my_index`. The
   !< realm-side body filters `flux_register%face(:)` by
   !< `face%coarse_realm == my_index` and writes the per-cell correction
   !< into its OWN integrator-private stage buffer (for RK realms:
   !< `self%rk%q_rk(:, ..., stage)`, weighted by `self%rk%ark(stage)`).
   !<
   !< The forest never reaches `realm%rk` directly: the integrator-specific
   !< weight pickup, the buffer name, and the per-cell write all live
   !< realm-side. This is the integrator-agnostic split — see
   !< [[realm_object]]%`apply_reflux_to_stage_forest` for the contract,
   !< and [[prism_cpu_object]] for the RK-specific override.
   !<
   !< Empty register fast path: when the register has no faces (single-
   !< realm forest, no seams declared) the per-realm calls each short-
   !< circuit on `flux_register%nfaces == 0` and the dispatch is a true
   !< no-op for the N=1 path.
   class(forest_object), intent(in)    :: self     !< The forest (holds the flux register).
   class(realm_object),  intent(inout) :: realm(:) !< Realms; each realm's reflux TBP fires once.
   integer(I4P),         intent(in)    :: stage    !< Integrator stage 1..K_total.
   real(R8P),            intent(in)    :: dt       !< Time step.
   integer(I4P)                        :: is       !< Realm index.

   if (.not. self%flux_register%is_initialized_) return
   if (self%flux_register%nfaces == 0_I4P)        return
   if (.not. allocated(self%flux_register%face))  return

   do is = 1_I4P, int(size(realm), I4P)
      call realm(is)%apply_reflux_to_stage_forest(my_index=is, stage=stage, dt=dt, &
                                                  flux_register=self%flux_register)
   enddo
   endsubroutine apply_reflux_corrections
endmodule adam_forest_object
