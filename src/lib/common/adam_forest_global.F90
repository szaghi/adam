!< ADAM, forest module-scope state (substage index + flux register).
module adam_forest_global
!< ADAM, forest module-scope state used by realm-side per-substage hooks.
!<
!< Phase D of issue #10 (design in issue #13). This module retains
!< program-scope state that is NOT a derived-type pointer:
!<
!<   * `forest_active_substage` — integer index, lockstep across realms.
!<   * `forest_flux_register`   — `type(flux_register_object)` (not a
!<     pointer), Berger-Colella reflux accumulator owned by the forest.
!<
!< Earlier versions of this module also exposed a
!< `class(realm_object), pointer :: forest_realm(:)` shim that the forest
!< set at the start of `evolve_one_step` and cleared at the end, so
!< realm-side `update_ghost` could refresh inter-realm ghosts at substage
!< depth without threading the realm array through the call chain. That
!< pattern is now retired (2026-05-25, issue #13): pointer-to-class is
!< forbidden by the CLAUDE.md rule because nvfortran/OpenACC mishandles
!< it under offloading. The realm array now flows exclusively as the
!< optional `realm(:)` dummy argument through the `_forest` contract
!< TBPs; the legacy single-INI path leaves the dummy absent and the
!< realm-side inter-realm branch becomes a no-op (single realm has no
!< peers — `inter_realm_neighbors` is unallocated).

use :: adam_flux_register_object, only : flux_register_object
use :: penf,                      only : I4P

implicit none
private
public :: forest_active_substage
public :: forest_flux_register

! Historical note (issue #13, 2026-05-25): an earlier version of this
! module exposed a `class(realm_object), pointer :: forest_realm(:)`
! shim that the forest set/cleared each step and realm-side code
! dereferenced at substage depth. That pointer-to-class pattern is
! forbidden by the CLAUDE.md rule (nvfortran/OpenACC mishandles it).
! The realm array now flows exclusively as a dummy argument through
! the `_forest` contract TBPs; this module retains only the scalar
! `forest_active_substage` and the program-scope `forest_flux_register`,
! both of which are NOT derived-type pointers and are safe.

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
