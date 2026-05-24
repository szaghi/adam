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
!<
!< Architectural invariants:
!<
!<   * **No derived-type-pointer components** on `forest_object`, ever (R1
!<     of issue #10). The chain-resolution bug class is excluded by
!<     construction.
!<   * **No app-specific imports** in this module. The forest's body uses
!<     only ADAM-lib-visible types (intrinsics with explicit kinds,
!<     `realm_object`, MPI). App-specific dispatch lives inside the
!<     overridden realm TBPs (O2).
!<   * **`q` is never touched here.** All cell-centered variable access is
!<     realm-side, via the realm's own `*_forest` TBPs (O3).

use :: adam_realm_object,    only : realm_object
use :: adam_maps_object,     only : inter_realm_neighbor_t, &
                                    FACE_X_MAX, FACE_X_MIN, FACE_Y_MAX, FACE_Y_MIN, FACE_Z_MAX, FACE_Z_MIN
use :: adam_forest_manifest, only : forest_manifest_t, forest_face_pair_t
use :: adam_forest_global,   only : forest_realm, forest_active_substage, forest_flux_register
use :: adam_flux_register_object, only : SEAM_KIND_INTER_REALM
use :: adam_globals,         only : mpih
use :: mpi
use :: penf

implicit none
private
public :: forest_object

type :: forest_object
   !< Behavior-only orchestrator of an array of realms.
   integer(I4P) :: n = 0_I4P !< Number of realms in the forest (set by initialize from size(realm)).
   ! NO derived-type-pointer components, ever (R1).
   ! Phase D will likely add: MPI sub-communicator handle, inter-realm
   ! connectivity descriptor (both intrinsic-typed).
   contains
      procedure, pass(self) :: initialize               !< Sequence each realm's initialize_forest at startup (single shared INI).
      procedure, pass(self) :: initialize_from_manifest !< Like initialize, but each realm reads its own INI from a forest manifest.
      procedure, pass(self) :: simulate                 !< Main entry point (single shared INI): drive the full simulation.
      procedure, pass(self) :: simulate_from_manifest   !< Main entry point (per-realm INIs via forest manifest).
      procedure, pass(self) :: compute_global_dt        !< Min-reduce each realm's compute_local_dt_forest across the forest.
      procedure, pass(self) :: evolve_one_step          !< Iterate realm(:)%advance_one_step_forest(dt) for one global timestep.
      procedure, pass(self) :: exchange_halos           !< Iterate realm(:)%exchange_inter_realm_halos_forest to refresh inter-realm ghosts.
      procedure, pass(self) :: post_step                !< Iterate realm(:)%post_step_forest for the per-step diagnostics/IO block.
      procedure, pass(self) :: is_done                  !< AND-reduce each realm's is_done_forest across the forest.
      procedure, pass(self) :: finalize                 !< Sequence each realm's finalize_forest at shutdown.
      procedure, pass(self), private :: populate_inter_realm_topology !< Translate manifest face-pairs into per-realm maps%inter_realm_neighbors.
endtype forest_object

contains
   subroutine initialize(self, realm, filename)
   !< Initialize the forest and every realm it tends.
   !<
   !< Records `n = size(realm)`, then iterates `realm(is)%initialize_forest`
   !< in increasing index order. Each realm receives its index via the
   !< `realm_index` optional argument — useful for per-realm IO basenames
   !< or rank carve-outs (Phase D).
   class(forest_object), intent(inout) :: self     !< The forest.
   class(realm_object),  intent(inout) :: realm(:) !< The realms to initialize.
   character(*),         intent(in)    :: filename !< Input parameters file name (shared across realms for v1).
   integer(I4P)                        :: is       !< Realm index.

   self%n = int(size(realm), I4P)
   do is = 1, self%n
      call realm(is)%initialize_forest(filename=filename, realm_index=is, realms_number=self%n)
   enddo
   endsubroutine initialize

   subroutine initialize_from_manifest(self, realm, manifest)
   !< Initialize the forest and every realm using per-realm INIs from a manifest.
   !<
   !< Like `initialize` but each realm receives its OWN INI path
   !< (`manifest%realm_ini(is)`) instead of a shared filename. After all
   !< realms are initialized, translates the manifest's face-pair list into
   !< per-realm `maps%inter_realm_neighbors` entries.
   !<
   !< The driver MUST allocate `realm(size = manifest%realms_number)` with
   !< the concrete app type before calling this — the forest does not
   !< allocate the realm array (it cannot, since realm_object is abstract
   !< and each app has its own extension).
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

   subroutine simulate(self, realm, filename)
   !< Drive the full simulation: initialize, time-loop, finalize.
   !<
   !< Top-level entry point the program driver calls instead of the legacy
   !< per-realm `realm%simulate`. The time loop is:
   !<
   !<   initialize → loop { evolve_one_step → post_step → is_done } → finalize
   !<
   !< Each step invokes the orchestrator-contract TBPs on every realm; the
   !< per-realm body decides what app-specific work to do internally.
   class(forest_object), intent(inout) :: self     !< The forest.
   class(realm_object),  intent(inout) :: realm(:) !< The realms to evolve.
   character(*),         intent(in)    :: filename !< Input parameters file name.
   logical                             :: done     !< Forest-global termination predicate.

   call self%initialize(realm, filename=filename)
   done = .false.
   do
      call self%evolve_one_step(realm)
      call self%post_step(realm)
      call self%is_done(realm, done=done)
      if (done) exit
   enddo
   call self%finalize(realm)
   endsubroutine simulate

   subroutine simulate_from_manifest(self, realm, manifest)
   !< Drive the full simulation using per-realm INIs from a manifest.
   !<
   !< Like `simulate` but uses `initialize_from_manifest` to populate each
   !< realm from its own INI file, and (via that initialize) wires the
   !< inter-realm topology from the manifest. The time-loop body is
   !< identical to `simulate`.
   class(forest_object),     intent(inout) :: self     !< The forest.
   class(realm_object),      intent(inout) :: realm(:) !< The realms to evolve.
   type(forest_manifest_t),  intent(in)    :: manifest !< Parsed manifest.
   logical                                 :: done     !< Forest-global termination predicate.

   call self%initialize_from_manifest(realm, manifest=manifest)
   done = .false.
   do
      call self%evolve_one_step(realm)
      call self%post_step(realm)
      call self%is_done(realm, done=done)
      if (done) exit
   enddo
   call self%finalize(realm)
   endsubroutine simulate_from_manifest

   subroutine compute_global_dt(self, realm, dt)
   !< Compute the forest-global stability-limited dt.
   !<
   !< Each realm reports its local dt via `compute_local_dt_forest`; the
   !< forest takes the min across all realms (intra-process) and then
   !< across all MPI ranks (`MPI_ALLREDUCE` on `MPI_COMM_WORLD`). Returns
   !< the bit-identical global min every rank should advance by.
   !<
   !< Phase D may replace `MPI_COMM_WORLD` with a forest-specific
   !< sub-communicator once rank carve-outs land.
   class(forest_object), intent(in)    :: self     !< The forest.
   class(realm_object),  intent(inout) :: realm(:) !< The realms to query (inout for the multi-realm path's shim re-bind side effect).
   real(R8P),            intent(out)   :: dt       !< Global stability-limited dt.
   real(R8P)                           :: dt_local !< Per-realm local dt.
   integer(I4P)                        :: is       !< Realm index.
   integer(I4P)                        :: ierr     !< MPI error code.

   dt = huge(0._R8P)
   do is = 1, int(size(realm), I4P)
      ! For N>1 each realm needs its singleton shims pointing at its own
      ! components before compute_local_dt_forest reads through them.
      if (int(size(realm), I4P) > 1_I4P) call realm(is)%bind_my_globals_forest
      call realm(is)%compute_local_dt_forest(dt_local=dt_local)
      dt = min(dt, dt_local)
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, dt, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
   endsubroutine compute_global_dt

   subroutine evolve_one_step(self, realm)
   !< Advance every realm by one global timestep.
   !<
   !< Two paths, selected by `size(realm)`:
   !<
   !<   * **N=1 fast path**: invokes the realm's `advance_one_step_forest`
   !<     (the legacy entry point). The substage loop is internal to that
   !<     method; the once-per-step inter-realm exchange runs afterwards.
   !<     For single-realm forests `exchange_halos` iterates an empty
   !<     neighbour list and the path is bit-identical to pre-Phase-D.
   !<
   !<   * **N>1 multi-realm path**: drives the substage loop itself. Each
   !<     realm contributes a prologue (`prepare_step_forest`), `nrk`
   !<     substages (`compute_substage_forest`), and an epilogue
   !<     (`finalize_step_forest`); the forest invokes inter-realm halo
   !<     exchange via `exchange_halos` BETWEEN substages so realm A's
   !<     substage-s ghosts see realm B's substage-s interior values —
   !<     this is what bit-comparability with single-realm rmf requires.
   !<     All realms are assumed to use the same integrator (same `nrk`);
   !<     a mismatch is flagged with `mpih%error_stop`.
   !<
   !< Per-substage inter-realm exchange: the `forest_realm` module pointer
   !< is bound to `realm(:)` for the duration of this method so the realm
   !< code at substage depth (inside `compute_residuals → update_ghost`)
   !< can refresh inter-realm ghosts without threading the realm array
   !< through the legacy signature chain. See [[adam_forest_global]] for
   !< the rationale.
   class(forest_object), intent(in)            :: self     !< The forest.
   class(realm_object),  intent(inout), target :: realm(:) !< The realms to advance.
   real(R8P)                                   :: dt       !< Global timestep size.
   integer(I4P)                                :: is, s    !< Realm and substage indices.
   integer(I4P)                                :: nrk      !< Number of substages (multi-realm path).
   integer(I4P)                                :: nrk_chk  !< Per-realm consistency check.

   forest_realm => realm
   ! Zero flux register accumulators at top of step (Phase A of [issue #13]).
   ! Skeleton commit: safe no-op when face(:) is unallocated, which is the
   ! case until the topology-registration follow-up commit populates it.
   call forest_flux_register%reset
   call self%compute_global_dt(realm, dt=dt)
   if (int(size(realm), I4P) == 1_I4P) then
      ! N=1 fast path — bit-identical to pre-Phase-D.
      call realm(1)%advance_one_step_forest(dt=dt)
      call self%exchange_halos(realm)
   else
      ! N>1 multi-realm path — forest drives the substage loop. Before each
      ! per-realm TBP we re-bind the legacy singleton shims to that realm's
      ! value components; the shims would otherwise alias the last-initialized
      ! realm and feed wrong geometry / RK state to every other realm's
      ! per-step code. See [[realm_object]]%`bind_my_globals_forest`.
      do is = 1_I4P, int(size(realm), I4P)
         call realm(is)%bind_my_globals_forest
         call realm(is)%prepare_step_forest(dt=dt)
      enddo
      call realm(1)%bind_my_globals_forest
      nrk = realm(1)%nrk_forest()
      do is = 2_I4P, int(size(realm), I4P)
         call realm(is)%bind_my_globals_forest
         nrk_chk = realm(is)%nrk_forest()
         if (nrk_chk /= nrk) call mpih%error_stop(msg='forest_object%evolve_one_step: realms disagree on nrk_forest')
      enddo
      do s = 1_I4P, nrk
         forest_active_substage = s
         ! Phase 1: every realm assembles its q_rk(:,...,s) (no ghost reads).
         do is = 1_I4P, int(size(realm), I4P)
            call realm(is)%bind_my_globals_forest
            call realm(is)%assemble_substage_forest(s=s, nrk=nrk, dt=dt)
         enddo
         ! Phase 2: every realm evaluates residuals; update_ghost inside
         ! compute_residuals refreshes inter-realm ghosts via the
         ! forest_realm pointer, which now sees peer substage-s buffers.
         do is = 1_I4P, int(size(realm), I4P)
            call realm(is)%bind_my_globals_forest
            call realm(is)%evaluate_substage_forest(s=s, nrk=nrk, dt=dt)
         enddo
         ! Phase 3: Berger-Colella reflux at coarse-fine interfaces
         ! (Phase A of [issue #13]). Skeleton commit: apply_reflux is a
         ! no-op pending the accumulation + correction follow-up commit.
         ! The dt/dx/weight arguments are placeholders; the real signature
         ! will read dx per-face from the topology, weight from
         ! adam_rk_object%ark(s).
         call forest_flux_register%apply_reflux(substage=s, dt=dt, dx_coarse=0._R8P, weight=0._R8P)
      enddo
      forest_active_substage = 0_I4P
      do is = 1_I4P, int(size(realm), I4P)
         call realm(is)%bind_my_globals_forest
         call realm(is)%finalize_step_forest(dt=dt)
      enddo
   endif
   forest_realm => null()
   endsubroutine evolve_one_step

   subroutine exchange_halos(self, realm)
   !< Refresh inter-realm ghost cells across all realms.
   !<
   !< Iterates `realm(is)%exchange_inter_realm_halos_forest(realm)` in
   !< increasing index order. For single-realm forests (N=1, current rmf)
   !< the inter-realm neighbour list is empty and each iteration is a
   !< no-op; the call exists so the forest's evolve loop has a uniform
   !< shape regardless of N.
   !<
   !< Granularity: this method is invoked once per global timestep, AFTER
   !< all realms have completed `advance_one_step_forest`. For
   !< bit-comparability with a single-realm reference, additional refresh
   !< points are typically required between RK substages (inside each
   !< realm's `advance_one_step_forest` body) — see issue #13 for the
   !< per-substage design that the first concrete N>1 use case will need.
   !< This method's once-per-step invocation is the floor, not the ceiling.
   class(forest_object), intent(in)    :: self     !< The forest.
   class(realm_object),  intent(inout) :: realm(:) !< The realms whose ghosts to refresh.
   integer(I4P)                        :: is       !< Realm index.

   associate(self_unused => self) ! method takes self for TBP-dispatch symmetry
   end associate
   do is = 1, int(size(realm), I4P)
      if (int(size(realm), I4P) > 1_I4P) call realm(is)%bind_my_globals_forest
      call realm(is)%exchange_inter_realm_halos_forest(realm=realm)
   enddo
   endsubroutine exchange_halos

   subroutine post_step(self, realm)
   !< Run every realm's post-step diagnostics / IO / AMR block.
   !<
   !< Iterates `realm(is)%post_step_forest` without overriding cadence —
   !< for v1, each realm uses its own internal cadence (the `do_*` flags
   !< stay at their default). The `dt`, `t`, `it` arguments are NOT
   !< passed here because the forest does not yet own canonical time state
   !< (realm-side time module singletons still drive them). Once the forest
   !< takes over time bookkeeping (Phase D), this method will pass them
   !< through.
   class(forest_object), intent(in)    :: self     !< The forest.
   class(realm_object),  intent(inout) :: realm(:) !< The realms to query.
   integer(I4P)                        :: is       !< Realm index.

   do is = 1, int(size(realm), I4P)
      if (int(size(realm), I4P) > 1_I4P) call realm(is)%bind_my_globals_forest
      call realm(is)%post_step_forest(dt=0._R8P, t=0._R8P, it=0_I4P)
   enddo
   endsubroutine post_step

   subroutine is_done(self, realm, done)
   !< Decide whether the whole forest has finished evolving.
   !<
   !< Each realm reports its local predicate via `is_done_forest`; the
   !< forest AND-reduces across all realms (intra-process) and then across
   !< all MPI ranks (`MPI_ALLREDUCE` on `MPI_COMM_WORLD`). AND-reduction
   !< means the forest keeps evolving as long as ANY realm wants to —
   !< matching the legacy single-realm semantics for v1 (with one realm
   !< the global predicate equals that realm's local one).
   class(forest_object), intent(in)    :: self     !< The forest.
   class(realm_object),  intent(inout) :: realm(:) !< The realms to query (inout for the multi-realm path's shim re-bind side effect).
   logical,              intent(out)   :: done     !< Forest-global termination predicate.
   logical                             :: done_local !< Per-realm local predicate.
   integer(I4P)                        :: is         !< Realm index.
   integer(I4P)                        :: ierr       !< MPI error code.

   done = .true.
   do is = 1, int(size(realm), I4P)
      if (int(size(realm), I4P) > 1_I4P) call realm(is)%bind_my_globals_forest
      call realm(is)%is_done_forest(done=done_local)
      done = done .and. done_local
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, done, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
   endsubroutine is_done

   subroutine finalize(self, realm)
   !< Shut down the forest and every realm it tends.
   !<
   !< Iterates `realm(is)%finalize_forest` in increasing index order; each
   !< realm closes its IO files, releases its resources, and finalizes its
   !< MPI handler.
   class(forest_object), intent(in)    :: self     !< The forest.
   class(realm_object),  intent(inout) :: realm(:) !< The realms to shut down.
   integer(I4P)                        :: is       !< Realm index.

   do is = 1, int(size(realm), I4P)
      if (int(size(realm), I4P) > 1_I4P) call realm(is)%bind_my_globals_forest
      call realm(is)%finalize_forest
   enddo
   ! MPI_FINALIZE is process-global: run it ONCE here, after every realm has done its
   ! MPI-using teardown above — not per realm inside finalize_forest, which would tear
   ! MPI down while later realms still need it (issue #13 rmf-2realm MPI_Type_f2c abort).
   if (size(realm) >= 1) call realm(1)%finalize_mpi_forest
   endsubroutine finalize

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
   !<
   !< Writes per-realm into `realm(is)%adam%maps%inter_realm_neighbors`.
   !< This is the post-C.3-closure shape (each realm owns its own adam value
   !< component). For N=1 the array on `realm(1)%adam%maps` is the same data
   !< the singleton `maps%inter_realm_neighbors` aliases through the shim;
   !< for N>1 each realm has its own array with only the face-pair entries
   !< whose `my_realm` matches that realm's index.
   class(forest_object),     intent(in)    :: self     !< The forest.
   class(realm_object),      intent(inout) :: realm(:) !< Initialized realms whose adam%maps gets populated.
   type(forest_manifest_t),  intent(in)    :: manifest !< Parsed manifest.
   integer(I4P), allocatable               :: per_realm_count(:)  !< How many neighbour entries each realm gets.
   integer(I4P), allocatable               :: per_realm_cursor(:) !< Write cursor per realm.
   integer(I4P)                            :: f, is    !< Face-pair and realm index counters.
   type(forest_face_pair_t)                :: pair     !< Loop alias.

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
      call set_neighbor(realm(pair%realm_a)%adam%maps%inter_realm_neighbors(per_realm_cursor(pair%realm_a)), &
                        my_realm=pair%realm_a, my_face=pair%face_a, &
                        peer_realm=pair%realm_b, peer_face=pair%face_b, &
                        coupling=pair%coupling)
      ! entry on realm_b's array: my=b, peer=a
      per_realm_cursor(pair%realm_b) = per_realm_cursor(pair%realm_b) + 1_I4P
      call set_neighbor(realm(pair%realm_b)%adam%maps%inter_realm_neighbors(per_realm_cursor(pair%realm_b)), &
                        my_realm=pair%realm_b, my_face=pair%face_b, &
                        peer_realm=pair%realm_a, peer_face=pair%face_a, &
                        coupling=pair%coupling)
   enddo
   ! Register inter-realm seams with the program-scope flux register
   ! (Phase A of [issue #13], step 2: topology registration).
   !
   ! For each (face-pair, block-on-coarse-side) tuple, one entry is added to
   ! the register. The "coarse" / "fine" labels follow the manifest's a/b
   ! ordering; for the current same-resolution (COUPLING_MIRROR) case the
   ! labels are conventional and the accumulator values are nominally equal
   ! on both sides — the reflux correction will be round-off zero in
   ! expectation. The structural cost (allocated registers, populated
   ! topology) is the same as for the true coarse-fine AMR case that will
   ! exercise these accumulators non-trivially in follow-up commits.
   call register_inter_realm_seams(realm, manifest)
   ! Build the per-cell inter-realm ghost map (Phase A of issue #13 — seam
   ! comm-map). Per-realm: enumerate every ghost cell in self's seam-block
   ! ghost region and resolve the (peer_realm, peer_block, peer_interior_cell)
   ! tuple. The runtime exchange then becomes a flat indexed loop, replacing
   ! the per-substage geometric find_peer_block + face-slab copy that misses
   ! corner / edge ghosts (defect B of issue #13).
   call build_inter_realm_ghost_cell_map(realm, manifest)
   contains
      subroutine set_neighbor(slot, my_realm, my_face, peer_realm, peer_face, coupling)
      type(inter_realm_neighbor_t), intent(out) :: slot
      integer(I4P),                 intent(in)  :: my_realm, my_face, peer_realm, peer_face, coupling
      slot%my_realm   = my_realm
      slot%my_block   = 0_I4P              ! resolved at exchange time by the realm-side override
      slot%my_face    = my_face
      slot%peer_realm = peer_realm
      slot%peer_block = 0_I4P
      slot%peer_face  = peer_face
      slot%coupling   = coupling
      endsubroutine set_neighbor

      subroutine register_inter_realm_seams(realm, manifest)
      !< Populate `forest_flux_register` from the manifest face-pairs.
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
      !<     `field%nv`, and `nrk` from realm-a's `rk%nrk`.
      !<
      !< For the current same-resolution (COUPLING_MIRROR) case, realm_a and
      !< realm_b carry identical `nv`/`nrk`/`nface_cells` by construction
      !< (validated upstream by the manifest's structural checks). A
      !< coarse-fine AMR case would resolve `fine_block(:)` by enumerating
      !< the realm_b-side blocks that geometrically cover the realm_a-side
      !< block face; for same-resolution that resolves to a single fine
      !< block, populated by `find_peer_block` at exchange time.
      class(realm_object),     intent(in) :: realm(:) !< Initialized realms.
      type(forest_manifest_t), intent(in) :: manifest !< Parsed manifest.
      integer(I4P)                        :: f, b    !< Face-pair, block counters.
      integer(I4P)                        :: a_realm !< Coarse-side realm index alias.
      integer(I4P)                        :: a_axis, a_sign !< Coarse-face axis and sign.
      integer(I4P)                        :: nfaces_total   !< Total register entries.
      integer(I4P)                        :: cursor         !< Write cursor into the register.
      integer(I4P)                        :: nface_cells    !< Cell count on the coarse-face skin.
      type(forest_face_pair_t)            :: pair    !< Manifest face-pair alias.

      if (.not. allocated(manifest%face_pairs)) then
         ! No inter-realm topology — initialize with zero faces so the
         ! register's `is_initialized_` flag flips and the per-step `reset`
         ! call becomes a safe no-op on the empty register.
         call forest_flux_register%initialize(nfaces=0_I4P)
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

      call forest_flux_register%initialize(nfaces=nfaces_total)
      if (nfaces_total == 0_I4P) return

      ! Pass 2: register one entry per (face-pair, seam-block-on-coarse-side).
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
            ! fine_block(:) is left empty (size 0) in this commit — the
            ! per-block peer resolution lands in the accumulation-wiring
            ! follow-up commit (option α, #13 §3.2). The accumulators are
            ! allocated and ready to receive fluxes; only the
            ! coarse↔fine block mapping is deferred.
            call forest_flux_register%register_face(face_index=cursor,                &
                                                    seam_kind=SEAM_KIND_INTER_REALM,  &
                                                    coarse_realm=pair%realm_a,        &
                                                    coarse_block=b,                   &
                                                    coarse_face=pair%face_a,          &
                                                    fine_realm=pair%realm_b,          &
                                                    fine_block=[integer(I4P) ::],     &
                                                    nface_cells=nface_cells,          &
                                                    nv=int(realm(a_realm)%adam%field%nv, I4P), &
                                                    nrk=realm(a_realm)%rk%nrk)
         enddo
      enddo
      endsubroutine register_inter_realm_seams

      pure subroutine face_axis_sign(face_code, axis, sgn)
      !< Translate FACE_X_MAX / FACE_X_MIN / ... into (axis 1..3, sign ±1).
      !<
      !< Duplicate of the same helper inside
      !< `prism_cpu_object%exchange_inter_realm_halos_forest`; a future
      !< refactor should lift the canonical version into `adam_maps_object`
      !< and have both call sites use it. Kept local here to keep this
      !< Phase A topology-registration commit minimal in scope.
      integer(I4P), intent(in)  :: face_code !< Face code from inter_realm_neighbor_t.
      integer(I4P), intent(out) :: axis      !< 1=x, 2=y, 3=z.
      integer(I4P), intent(out) :: sgn       !< +1 if MAX, -1 if MIN.

      select case (face_code)
      case (FACE_X_MAX); axis = 1_I4P; sgn = +1_I4P
      case (FACE_X_MIN); axis = 1_I4P; sgn = -1_I4P
      case (FACE_Y_MAX); axis = 2_I4P; sgn = +1_I4P
      case (FACE_Y_MIN); axis = 2_I4P; sgn = -1_I4P
      case (FACE_Z_MAX); axis = 3_I4P; sgn = +1_I4P
      case (FACE_Z_MIN); axis = 3_I4P; sgn = -1_I4P
      case default;      axis = 0_I4P; sgn = 0_I4P
      end select
      endsubroutine face_axis_sign

      function block_face_on_realm_boundary(this_realm, b, axis, sgn) result(yes)
      !< Return .true. iff block `b`'s face on (axis, sgn) lies on the realm boundary.
      !<
      !< Geometric test against `grid%domain_emin/emax` and
      !< `field%emin/emax`, with a small absolute tolerance. This mirrors
      !< the identically-named helper inside the PRISM-CPU realm; a future
      !< refactor should lift the canonical version into `adam_maps_object`.
      class(realm_object), intent(in) :: this_realm  !< Realm to query.
      integer(I4P),        intent(in) :: b           !< Block index.
      integer(I4P),        intent(in) :: axis        !< 1=x, 2=y, 3=z.
      integer(I4P),        intent(in) :: sgn         !< +1 if checking MAX face, -1 if MIN.
      logical                         :: yes         !< Test result.
      real(R8P)                       :: face_coord, target_coord, tol

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
      class(realm_object), intent(in) :: this_realm
      integer(I4P),        intent(in) :: axis
      integer(I4P)                    :: n

      associate(g => this_realm%adam%grid)
      select case (axis)
      case (1_I4P); n = g%nj * g%nk
      case (2_I4P); n = g%ni * g%nk
      case (3_I4P); n = g%ni * g%nj
      case default; n = 0_I4P
      end select
      end associate
      endfunction tangential_cell_count

      ! ---------------------------------------------------------------------
      ! Phase A of issue #13 — inter-realm seam comm-map construction.
      ! All routines below are siblings inside this contains block (Fortran
      ! 2008 forbids contains nesting deeper than one level), and pass the
      ! realm array + per_realm_count explicitly rather than via host
      ! association to keep the dependency structure obvious.
      ! ---------------------------------------------------------------------

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
      class(realm_object),     intent(inout) :: realm(:) !< Forest realms (their maps get populated).
      type(forest_manifest_t), intent(in)    :: manifest !< Parsed manifest.
      integer(I4P)                           :: f, is, ip  !< Counters: face-pair, self-realm, peer-realm.
      integer(I4P)                           :: my_face, peer_face, my_realm_idx, peer_realm_idx
      integer(I4P)                           :: my_axis, my_sign
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

      ! ---- Pass 1: per-realm row counts --------------------------------
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

      ! ---- Allocate per-realm maps ------------------------------------
      do is = 1_I4P, int(size(realm), I4P)
         if (allocated(realm(is)%adam%maps%inter_realm_ghost_cell)) &
            deallocate(realm(is)%adam%maps%inter_realm_ghost_cell)
         if (per_realm_count(is) > 0_I4P) &
            allocate(realm(is)%adam%maps%inter_realm_ghost_cell(1:per_realm_count(is), 1:10))
      enddo

      ! ---- Pass 2: populate -------------------------------------------
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
endmodule adam_forest_object
