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

use :: adam_flux_register_object, only : flux_register_object
use :: adam_realm_object,         only : realm_object
use :: penf,                      only : I4P

implicit none
private
public :: forest_realm
public :: forest_active_substage
public :: forest_flux_register

class(realm_object), pointer :: forest_realm(:) => null() !< Set by forest%evolve_one_step for the duration of the step.

!< Berger-Colella flux register — coarse-fine interface machinery (Phase A of [issue #13]).
!<
!< Program-scope singleton (like `mpih`) rather than a `forest_object`
!< component, because:
!<
!<   * `forest_object` is behavior-only by invariant R1 of [issue #10] —
!<     "no derived-type-pointer components, ever" — and the register's
!<     internal `face(:)` is exactly such a component;
!<   * the register is indexed by (realm pair, block, face), not by
!<     a single realm, so it has no natural single-realm owner;
!<   * cross-rank reduce in `apply_reflux` wants a single rank-scope
!<     owner, mirroring the way `mpih` carries per-rank MPI state.
!<
!< Lifecycle is driven by `forest_object`:
!<
!<   * populated by `populate_inter_realm_topology` at init / regrid time;
!<   * `reset` at top of `evolve_one_step`;
!<   * `accumulate_*_flux` from PRISM's `compute_residuals_*` during substage;
!<   * `apply_reflux` between substages, before the next ghost exchange.
!<
!< In the skeleton commit of Phase A (this commit), every operation on the
!< register is a no-op or near-no-op: the type compiles, the lifecycle hooks
!< exist in `forest_object`, the wiring is in place. Real accumulation /
!< correction follow in subsequent commits of Phase A.
type(flux_register_object) :: forest_flux_register !< Coarse-fine interface reflux machinery.

!< INVARIANT — all realms advance their RK substages in LOCKSTEP.
!<
!< `forest_active_substage` is a SINGLE program-scope index shared by every
!< realm. It is correct only because `forest_object%evolve_one_step` drives one
!< substage loop over ALL realms together: it sets `forest_active_substage = s`
!< ONCE per outer substage and processes every realm at that `s` (assemble then
!< evaluate) before advancing — see adam_forest_object.F90 (the `do s = 1, nrk`
!< loop). The forest also asserts up front that every realm agrees on `nrk`
!< (`error_stop` otherwise, same routine), which forbids the most likely way to
!< break lockstep. Reader semantics: `> 0` selects the peer's substage-s buffer
!< in the inter-realm halo copy; `0` is the N=1 fast path (use `self%q`).
!<
!< This single global would become SILENTLY WRONG under realm-local subcycling
!< (different realms at different substages or different `nrk`) or any overlap
!< of two realms' substage phases in time. Neither is the current model. If
!< subcycling is ever scoped, this MUST become per-realm state (e.g. carried on
!< the realm via `bind_my_globals_forest`) — do not keep it as a shared scalar.
integer(I4P)                 :: forest_active_substage = 0_I4P !< Active RK substage shared by all realms; see the lockstep invariant above.
endmodule adam_forest_global
