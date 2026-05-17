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

use :: adam_realm_object, only : realm_object
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
      procedure, pass(self) :: initialize        !< Sequence each realm's initialize_forest at startup.
      procedure, pass(self) :: simulate          !< Main entry point: drive the full simulation.
      procedure, pass(self) :: compute_global_dt !< Min-reduce each realm's compute_local_dt_forest across the forest.
      procedure, pass(self) :: evolve_one_step   !< Iterate realm(:)%advance_one_step_forest(dt) for one global timestep.
      procedure, pass(self) :: post_step         !< Iterate realm(:)%post_step_forest for the per-step diagnostics/IO block.
      procedure, pass(self) :: is_done           !< AND-reduce each realm's is_done_forest across the forest.
      procedure, pass(self) :: finalize          !< Sequence each realm's finalize_forest at shutdown.
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
      call realm(is)%initialize_forest(filename=filename, realm_index=is)
   enddo
   endsubroutine initialize

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
   class(forest_object), intent(in)  :: self     !< The forest.
   class(realm_object),  intent(in)  :: realm(:) !< The realms to query.
   real(R8P),            intent(out) :: dt       !< Global stability-limited dt.
   real(R8P)                         :: dt_local !< Per-realm local dt.
   integer(I4P)                      :: is       !< Realm index.
   integer(I4P)                      :: ierr     !< MPI error code.

   dt = huge(0._R8P)
   do is = 1, int(size(realm), I4P)
      call realm(is)%compute_local_dt_forest(dt_local=dt_local)
      dt = min(dt, dt_local)
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, dt, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
   endsubroutine compute_global_dt

   subroutine evolve_one_step(self, realm)
   !< Advance every realm by one global timestep.
   !<
   !< Computes the global dt via `compute_global_dt`, then iterates
   !< `realm(is)%advance_one_step_forest(dt)` in increasing index order.
   !< For v1 (N=1) this is exactly one realm-side timestep.
   class(forest_object), intent(in)    :: self     !< The forest.
   class(realm_object),  intent(inout) :: realm(:) !< The realms to advance.
   real(R8P)                           :: dt       !< Global timestep size.
   integer(I4P)                        :: is       !< Realm index.

   call self%compute_global_dt(realm, dt=dt)
   do is = 1, int(size(realm), I4P)
      call realm(is)%advance_one_step_forest(dt=dt)
   enddo
   endsubroutine evolve_one_step

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
   class(forest_object), intent(in)  :: self     !< The forest.
   class(realm_object),  intent(in)  :: realm(:) !< The realms to query.
   logical,              intent(out) :: done     !< Forest-global termination predicate.
   logical                           :: done_local !< Per-realm local predicate.
   integer(I4P)                      :: is         !< Realm index.
   integer(I4P)                      :: ierr       !< MPI error code.

   done = .true.
   do is = 1, int(size(realm), I4P)
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
      call realm(is)%finalize_forest
   enddo
   endsubroutine finalize
endmodule adam_forest_object
