# rmf-2realm regression case

Two-realm x-split version of the single-realm `rmf` regression case (Phase D
of issue #10 / D.4 of issue #13). Domain is split at x=0 into two adjacent
realms; the union of cell centers reconstructs the single-realm rmf grid
bit-for-bit, so the same physics integrated under the multi-realm forest
orchestrator MUST produce digests that match the single-realm rmf golden
within the project's standard regression tolerances.

## Files

```
rmf-2realm/
   input.ini       # forest manifest (auto-detected as a manifest by PRISM)
   realm_1.ini     # complete PRISM INI for realm 1 (x ∈ [-0.06, 0],   ni=8, nj=16, nk=16)
   realm_2.ini     # complete PRISM INI for realm 2 (x ∈ [ 0,  0.06], ni=8, nj=16, nk=16)
   golden/
      cpu/         # CPU-backend golden (digest.txt + per-realm residuals.dat)
      fnl/         # FNL-backend golden (see "FNL caveat" below)
```

The two `realm_*.ini` files differ only in:

- `[grid] emin_x`, `emax_x` (each realm spans half of the original x extent)
- `[IO] output_basename` (`rmf_2realm_r1` / `rmf_2realm_r2`)
- `[IO] restart_basename`

Every other section — `[numerics]`, `[physics]`, `[fdv]`, `[runge_kutta]`,
`[weno]`, `[time]`, `[coils_input]`, all four `[coils_input_coil_<n>]`
blocks, the BC sections — is byte-identical between realms and to the
single-realm `rmf/input.ini`.

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

The FNL backend follows the same recipe via `run-fnl-local.sh` instead of
`run.sh`. See the FNL caveat below for current status.

## Cross-config sanity check

In addition to passing against its own golden (the automatic harness
check), the `rmf-2realm` digest should also match the `rmf` single-realm
digest for every physics-meaningful FIELD variable — that is the deeper
correctness oracle for the seam machinery. The check is a manual one-liner
because the harness compares each case against its own golden only:

```bash
exe/.regression-venv/bin/python src/tests/prism/regression/digest.py compare \
    src/tests/prism/regression/rmf-2realm/work-cpu/digest.txt \
    src/tests/prism/regression/rmf/golden/cpu/digest.txt
```

Expected outcome:

```
>> 9 per-block-metadata discrepancies skipped (benign across block-partitioning changes; ...)
  SKIP_METADATA 000000000 / dxdydz — count 384 != 192
  SKIP_METADATA 000000000 / origin — count 384 != 192
  SKIP_METADATA 000000000 / time_iteration — count 128 != 64
  ... (3 more per output iteration)
>> digest match: 102 rows within rtol=1e-06 atol=0.001
```

The 9 SKIP_METADATA lines are benign per-block-metadata mismatches:
`dxdydz`, `origin`, `time_iteration` are per-block attributes (3 / 3 / 1
values per block respectively) whose digest reductions scale with the block
count, and the 2-realm config has 2× the blocks of the 1-realm config.
These rows are downgraded from MISMATCH to SKIP_METADATA by `digest.py
compare` automatically (see `_PER_BLOCK_METADATA_VARS` in `digest.py`); the
overall comparison passes when every FIELD variable matches.

This cross-config equivalence holds because:

1. Cell centers in realm_1 ∪ realm_2 = cell centers in single-realm rmf (by
   construction — see "Grid-alignment rationale" above).
2. The manifest's `[forest.topology.face_*]` declaration overrides each
   realm's INI-declared BC for the seam face (commit 6987c20, BC_SEAM
   sentinel) — without this override, PRISM's `set_boundary_conditions`
   would extrapolate Neumann into the seam ghost cells and overwrite the
   peer-interior values the inter-realm exchange wrote.
3. The forest's three-phase stage loop (same commit) keeps each realm's
   `q_rk(:, interior, :, b, k)` consistent across all peers' seam reads —
   the legacy two-phase loop (`compute_residuals` + `rk%assign_stage`
   back-to-back per realm) raced because `assign_stage` overwrites
   `q_rk(k)` with `dq` in-place.
4. Ghost cells are stripped from the digest, so the inter-realm-mirror
   ghost cells (which differ from single-realm's neighbouring-block reads
   only in how they are populated, not in what data they hold
   post-exchange) don't pollute the comparison.

## FNL caveat (Phase D.4)

The FNL backend has a known limitation in the multi-realm path: `rk_fnl`
(and other `*_fnl` singletons) is a program-scope singleton, not a per-realm
component. With two realms the stage state on the device is
overwritten on each realm's `begin_stage_forest`, so the per-stage
inter-realm exchange cannot reach a bit-comparable state. The FNL exchange
override therefore error-stops on the per-stage entry, surfacing the
limitation rather than silently producing wrong results.

Resolving this is a follow-up that promotes the FNL singletons (rk_fnl,
field_fnl, ib_fnl, weno_fnl, coil_fnl, fwlayer_fnl) to per-realm value
components on `prism_fnl_object`, paralleling the C.3-closure pattern
applied to the CPU side. Until that lands the FNL rmf-2realm golden cannot
be captured.

## Per-realm INI duplication

`realm_1.ini` and `realm_2.ini` are essentially identical except for the
three keys listed above. This duplication is accepted as the cost of
per-realm autonomy (every PRISM INI key can be varied independently across
realms without growing a per-section override schema). If duplication
becomes painful in later use cases an `[include]` mechanism can be added
to FiNeR without touching the manifest design.
