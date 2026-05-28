# rmf-2realm regression case

Two-realm x-split version of the single-realm `rmf` regression case (Phase D
of issue #10 / D.4 of issue #13, ratified under α — PRD #16). Domain is split
at x=0 into two adjacent realms; the union of cell centers reconstructs the
single-realm rmf grid bit-for-bit.

The case anchors the multi-realm forest orchestrator's CPU path. The
companion case `rmf-2realm-asymK` covers the asymmetric-K mode introduced by
α (realm_1: SSP-RK-3, realm_2: SSP-RK-5).

## Files

```
rmf-2realm/
   input.ini       # forest manifest (auto-detected as a manifest by PRISM)
   realm_1.ini     # complete PRISM INI for realm 1 (x ∈ [-0.06, 0],   ni=8, nj=16, nk=16)
   realm_2.ini     # complete PRISM INI for realm 2 (x ∈ [ 0,  0.06], ni=8, nj=16, nk=16)
   golden/
      cpu/         # CPU-backend α golden (digest.txt + per-realm residuals.dat)
```

The companion case `../rmf-2realm-asymK/` extends this geometry to
asymmetric per-realm K (realm_1: SSP-RK-3, realm_2: SSP-RK-5); see PRD #16
M7 for the rationale.

The two `realm_*.ini` files differ only in:

- `[grid] emin_x`, `emax_x` (each realm spans half of the original x extent)
- `[IO] output_basename` (`rmf_2realm_r1` / `rmf_2realm_r2`)
- `[IO] restart_basename`

Every other section — `[numerics]`, `[physics]`, `[fdv]`, `[runge_kutta]`,
`[weno]`, `[time]`, `[coils_input]`, all four `[coils_input_coil_<n>]`
blocks, the BC sections — is byte-identical between realms and to the
single-realm `rmf/input.ini`.

There is no FNL golden yet — see "FNL status" below.

## Grid-alignment rationale

Single-realm `rmf` uses `[grid] ni=16, nj=16, nk=16` with `[amr] iu_ref_levels=2,
ratio=8`. At level 2 each direction has `2^2 = 4` blocks, giving block_dxyz =
0.03 and cell_dxyz = 0.001875.

For each two-realm realm to share the same cell-center positions:

- `ni = 8`, `nj = 16`, `nk = 16` and the same `iu_ref_levels = 2`, `ratio = 8`
- `emax_x - emin_x = 0.06` (half the single-realm extent in x)
- `emax_y - emin_y = emax_z - emin_z = 0.12` (full original extent in y, z)

At level 2 per realm: `nb_max(2) = 4` blocks per direction, giving block_dxyz =
(0.015, 0.03, 0.03) and cell_dxyz = (0.001875, 0.001875, 0.001875) in all
three directions — bit-equal to the single-realm cell_dxyz. The union of
realm_1 ∪ realm_2 cell centers reconstructs single-realm rmf exactly.

## Capturing the goldens

Both the single-realm `rmf` and the two-realm `rmf-2realm` CPU goldens are
already committed under the current digest format (step-index keying +
ghost-cell stripping; see `digest.py` header). To recapture either:

```bash
# Recapture single-realm rmf golden
./run.sh cpu                                              # produces rmf/work-cpu/{digest.txt, *-residuals.dat}
cp rmf/work-cpu/digest.txt                   rmf/golden/cpu/
cp rmf/work-cpu/rmf_regression-residuals.dat rmf/golden/cpu/

# Recapture two-realm golden
REGRESSION_RUN_GOLDENLESS=1 ./run.sh cpu                  # bypasses skip-when-no-golden
cp rmf-2realm/work-cpu/digest.txt                  rmf-2realm/golden/cpu/
cp rmf-2realm/work-cpu/rmf_2realm_r*-residuals.dat rmf-2realm/golden/cpu/
```

The FNL backend cannot be rebaselined locally (WSL2 runtime failure; see
"FNL status" below). Cluster validation under `run-fnl-local.sh` will
follow the same recipe once the runtime issue is resolved.

## α end-of-step barrier semantics

PRD #16 (M1–M7) restructured `forest_object%evolve_one_step` to follow the
AMReX coarse-fine interface convention (`AMReX_Amr.cpp::timeStep`,
`FillBoundary`, `Reflux`). Under α the seam coupling between realms is
once-per-step, not once-per-stage:

- **Mid-step peer ghosts are intentionally stale-by-one-step.** During RK
  stages 1..K each realm reads the seam ghosts established by the previous
  step's end-of-step exchange. This mirrors AMReX's `FillCoarsePatch`
  reading coarse `t^n` data during fine sub-steps (Berger-Oliger 1984) and
  is well-understood numerically: first-order seam coupling in time, full
  per-realm RK order in the interior.
- **End-of-step seam fill** fires once per global step, after every realm
  has completed `close_step_forest`. At that point every realm's
  `stage_active == 0` and `fill_seam_from_peer_forest` reads peer's
  committed `q` automatically (no new TBP).
- **Reflux at α.r1** is single-stage: the flux register's third axis is
  collapsed to 1; `apply_reflux_to_stage_forest` and the FV
  `accumulate_seam_fluxes_fv` call are gated on `stage == rk%nrk` so real
  work fires once per realm per step at its own final substage.

A consequence the regression validates: bit-identity with the single-realm
`rmf` digest is **no longer expected**. The α discretization deliberately
differs from the pre-α (and from single-realm) discretization at the seam.
This case's golden was rebaselined under α in commit `2df013a0` (PRD #16 M6).

### Seam localisation invariant

The α seam coupling error must remain spatially localised within one
block-width of x=0. The structural test (run before M6 ratified the
golden):

| distance from seam | max |Bz| | |Bz| / |By| |
|---|---|---|
| 0 – 0.015 (seam-adjacent block) | 1.0e-07 | 6.4% |
| 0.015 – 0.030 (next block out)  | 5.4e-10 | 0.04% |
| 0.030 – 0.045 (mid-interior)    | 4.3e-10 | 0.03% |
| 0.045 – 0.060 (far-interior)    | 4.3e-10 | 0.03% |

A ≥200× drop crossing one block (8 cells) into the interior confirms the
coupling is structurally seam-local, not interior-leaking. If a future
refactor of `forest_object` or PRISM's intra-stage routines spreads `|Bz|`
beyond the seam-adjacent slab, that is a hidden ordering assumption — see
PRD #16 risk clause and re-run the spatial check before rebaselining.

### Why the manifest still overrides the seam face BC

The manifest's `[forest.topology.face_*]` declaration overrides each realm's
INI-declared BC for the seam face (BC_SEAM sentinel, commit `6987c20`).
Without this override, PRISM's `set_boundary_conditions` would extrapolate
Neumann into the seam ghost cells and overwrite the peer-interior values
the inter-realm exchange wrote. This is unchanged by α — the BC override is
orthogonal to the per-stage vs end-of-step cadence change.

## FNL status

The architectural blocker that prevented FNL multi-realm has been resolved
on develop: `field_fnl`, `rk_fnl`, `ib_fnl`, `weno_fnl` (and PRISM-specific
`coil_fnl`, `fwlayer_fnl`) are now per-realm value components on
`prism_fnl_object` instead of program-scope singletons (commits
`8000ae7b`, `d2fb7bc8`, `90d0343e`, `76760724`). The α end-of-step gate is
in place on the FNL twin (commit `ebd5024d`, PRD #16 M5); the body remains
a no-op pending the FNL FV residual port (PRD #16 out-of-scope, Phase B).

A runtime crash remains on WSL2: `rmf-2realm/fnl` aborts at
`receive_recv_buffer_ghost_gpu_dev:139` with `CUDA_ERROR_INVALID_VALUE`
after the first cross-realm seam exchange. Standalone FUNDAL tests
(including `mpirun -np 2`) all pass; single-realm `rmf/fnl` passes. The
failure is not reproducible outside PRISM and is hypothesised to be a
WSL-UCX × OpenACC × FUNDAL device-memcpy interaction (the WSL `/dev/dxg`
shim does not retrieve the GPU primary context through UCX rendezvous;
see CLAUDE.md "Development Environment: WSL2 GPU+MPI Caveats").

Cluster validation is required before the FNL rmf-2realm golden can be
captured. Until then the case runs only under CPU.

## Per-realm INI duplication

`realm_1.ini` and `realm_2.ini` are essentially identical except for the
three keys listed above. This duplication is accepted as the cost of
per-realm autonomy (every PRISM INI key can be varied independently across
realms without growing a per-section override schema). If duplication
becomes painful in later use cases an `[include]` mechanism can be added
to FiNeR without touching the manifest design.
