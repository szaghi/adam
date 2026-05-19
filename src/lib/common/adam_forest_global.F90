!< ADAM, global forest-realm pointer shim — used by realm-side per-substage hooks.
module adam_forest_global
!< ADAM, global forest-realm pointer shim — used by realm-side per-substage hooks.
!<
!< Phase D of issue #10 (design in issue #13). The forest's
!< `exchange_inter_realm_halos_forest` TBP is invoked once per global timestep by
!< the orchestrator (see [[forest_object]]%`exchange_halos`). For bit-comparability
!< with single-realm runs, the inter-realm halo MUST also be refreshed between
!< RK substages — that finer schedule is realm-side and runs inside the existing
!< `update_ghost` path. The realm code at that depth does NOT have the
!< `realm(:)` array in scope, because threading it through the legacy
!< `compute_residuals → update_ghost` chain would touch every app backend and
!< every per-stage helper.
!<
!< This module exposes a single program-scope pointer that the forest sets at
!< the start of `evolve_one_step` and clears at the end. The realm's
!< `update_ghost` checks the pointer and, if associated, invokes the
!< `exchange_inter_realm_halos_forest` override on every call — automatically
!< giving per-substage granularity for every substage shape the integrators
!< implement, without surgery in the substage loops themselves.
!<
!< Semantics:
!<
!<   * **Unassociated** outside `forest%evolve_one_step` — the realm-side hook
!<     becomes a no-op, preserving Phase A/B behaviour (single-step path).
!<   * **Associated** during `forest%evolve_one_step` — the realm-side hook
!<     iterates the realm's inter-realm neighbour list and refreshes ghosts.
!<     For single-realm forests the neighbour list is unallocated, so the
!<     iteration is empty and the call is bit-equivalent to a no-op.
!<
!< The pointer is a `class(realm_object)` pointer array so it can alias the
!< driver's monomorphic concrete array (e.g. `type(prism_cpu_object) :: realm(N)`)
!< by polymorphic dispatch — the forest passes the same array under
!< `class(realm_object), intent(inout) :: realm(:)` already, so the pointer
!< inherits the same dispatch rules.
!<
!< Architectural notes:
!<
!<   * **Not a long-lived alias.** Unlike the seven `adam_*_global` shims,
!<     which are bound once at startup, this pointer is set/cleared by the
!<     time loop. The "global" name keeps the family naming but the lifetime
!<     is per-step.
!<   * **No data ownership.** The pointer aliases the driver's realm array;
!<     it does not allocate or deallocate. Setting and clearing are O(1).
!<   * **Single-INI compatibility.** Even for the legacy single-INI path,
!<     `forest%simulate` calls `evolve_one_step`, so the pointer is set;
!<     the realm-side hook iterates an empty neighbour list and the
!<     behaviour is bit-identical to pre-Phase-D.

use :: adam_realm_object, only : realm_object
use :: penf,              only : I4P

implicit none
private
public :: forest_realm
public :: forest_active_substage

class(realm_object), pointer :: forest_realm(:) => null() !< Set by forest%evolve_one_step for the duration of the step.
integer(I4P)                 :: forest_active_substage = 0_I4P !< Set by the forest's substage loop on the multi-realm path; 0 elsewhere (N=1 fast path uses self%q directly).
endmodule adam_forest_global
