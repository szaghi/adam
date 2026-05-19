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

> **Important**: capturing the rmf-2realm golden also requires re-capturing
> the single-realm rmf golden because Phase D introduced a digest format
> change (step-index keying + ghost-cell stripping; see `digest.py` header).
> The old `rmf/golden/<backend>/digest.txt` files are NOT comparable with
> digests produced by the post-Phase-D `digest.py`.

```bash
# 1. Recapture single-realm rmf golden under the new digest format
./run.sh cpu                                            # produces rmf/work-cpu/{digest.txt, *-residuals.dat}
cp rmf/work-cpu/digest.txt                rmf/golden/cpu/
cp rmf/work-cpu/rmf_regression-residuals.dat rmf/golden/cpu/

# 2. Capture two-realm golden
./run.sh cpu                                            # produces rmf-2realm/work-cpu/{digest.txt, per-realm *-residuals.dat}
cp rmf-2realm/work-cpu/digest.txt           rmf-2realm/golden/cpu/
cp rmf-2realm/work-cpu/rmf_2realm_r*-residuals.dat rmf-2realm/golden/cpu/

# 3. Sanity-check the cross-config equivalence — the two digests SHOULD match
#    within the regression tolerance after ghost-cell stripping.
python3 digest.py compare rmf-2realm/golden/cpu/digest.txt rmf/golden/cpu/digest.txt
```

The FNL backend follows the same recipe via `run-fnl-local.sh` instead of
`run.sh`.

## Cross-config equivalence

The acceptance criterion is that **the rmf-2realm digest matches the rmf
single-realm digest** (within the same `rtol=1e-6, atol=1e-3` tolerance the
suite uses for cross-compiler matching). This holds because:

1. Cell centers in realm_1 ∪ realm_2 = cell centers in single-realm rmf (by
   construction — see "Grid-alignment rationale" above).
2. The inter-realm halo exchange at every RK substage propagates peer-realm
   interior values across the x=0 interface, so WENO/FDV stencils at the
   interface see the SAME data they would see in single-realm rmf (which
   has no interface — the cells straddling x=0 are interior cells reading
   from neighbouring blocks).
3. Ghost cells are stripped from the digest, so the inter-realm ghost cells
   (which differ from single-realm's neighbouring-block reads only in how
   they are populated, not in what data they hold post-exchange) don't
   pollute the comparison.

## FNL caveat (Phase D.4)

The FNL backend has a known limitation in the multi-realm path: `rk_fnl`
(and other `*_fnl` singletons) is a program-scope singleton, not a per-realm
component. With two realms the substage state on the device is
overwritten on each realm's `assemble_substage_forest`, so the per-substage
inter-realm exchange cannot reach a bit-comparable state. The FNL exchange
override therefore error-stops on the per-substage entry, surfacing the
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
