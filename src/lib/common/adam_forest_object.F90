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
use :: adam_maps_object,     only : inter_realm_neighbor_t
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
   ! (Phase A of [issue #13]). Skeleton commit: `initialize` is a no-op
   ! beyond setting the count; per-face `register_face` and per-block
   ! expansion land in the topology-registration follow-up commit.
   call forest_flux_register%initialize(nfaces=int(size(manifest%face_pairs), I4P))
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
   endsubroutine populate_inter_realm_topology
endmodule adam_forest_object
