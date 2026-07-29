# PRISM regression suite

The `src/tests/prism/regression/` suite is the **structural-change regression baseline** for the PRISM Maxwell/EM solver — every step of the [forest-of-trees migration](https://github.com/szaghi/adam/issues/10) must leave it green. Each case runs in roughly a minute per backend and validates one backend against committed reference outputs (or, for the AMR-seam cases, against a physics oracle).

This page is the conceptual overview. The authoritative operational reference — golden-capture workflow, CI-as-authority rules, tolerance rationale, adding a case — is the suite's own [`README.md`](https://github.com/szaghi/adam/blob/develop/src/tests/prism/regression/README.md); it is not duplicated here. To stand the suite up on a new machine, compiler, or GPU architecture, follow the [bring-up tutorial](./prism-regression-tutorial).

::: warning Current status — two open items (as of `cf16e20d`)
The suite runs again after the `[fWLayer]` input migration, but is not fully green: **6 PASS, 1 FAIL, 3 SKIP** on a local CPU sweep.

1. **`rmf-amr` digest mismatch.** It is the only case on `fv_centered`; its golden dates from `d15fed4c` (2026-07-04) and `28625dbf` has since reworked FV coil initialisation. The input changed only in comments, so this is a source-behaviour change on the FV path — previously masked by the fWLayer breakage. It needs adjudication (correct behaviour → golden bump; incorrect → fix), **not a reflexive golden refresh**.
2. **`rmf-fwl` is still un-migrated** and still `error_stop`s: its layer is active, so translating `C = 6` to a physical `width` is not a guaranteed round-trip and its FNL golden must be re-verified.

Details and evidence in the [suite README](https://github.com/szaghi/adam/blob/develop/src/tests/prism/regression/README.md).
:::

## Design goals

- **Catch structural regressions, not physics changes.** The suite is derived from research cases but CI-tuned: small grids (`16³` or `8×16×16` per realm), `it_max ≈ 5–10`, per-step residual logging. It exists to catch a refactor that silently changes a reduction, a kernel chain, or a seam fill — not to validate the physics (that is the research cases' job).
- **Two backends, one baseline.** Every case runs on `prism-cpu-gnu` (the CI-gated GNU/CPU baseline) and `prism-fnl-nvf` (the workstation-only OpenACC/GPU backend). A case is fully covered only when goldens exist for both.
- **Fixed two-rank runs.** Every case runs at `mpirun -np 2`, so the golden is a two-rank golden and the suite exercises the MPI ghost-exchange and decomposition paths on every run. (A one-rank rerun will not match a two-rank golden.)
- **Deterministic across compilers.** AMR runtime refinement is disabled (`amr_iterations = 0`, `[amr] frequency = 999999`) so no float-comparison-sensitive refinement decision happens mid-run; initial uniform refinement (`iu_ref_levels = 2`) still gives 64 leaf blocks to distribute across the two ranks.

## The harness (`run.sh`)

```bash
# CPU backend (also the CI job)
src/tests/prism/regression/run.sh cpu

# FNL backend — workstation only, via the wrapper that loads nvhpc + WSL UCX knobs
src/tests/prism/regression/run-fnl-local.sh
```

`run.sh {cpu|fnl} [--no-build] [--varset <name>]` builds the backend, then **auto-discovers** cases: every subdirectory with an `input.ini` is a case. A case with no `golden/<backend>/` directory is **skipped** (it is not yet an anchor); set `REGRESSION_RUN_GOLDENLESS=1` to run it anyway for the initial golden capture. Before any case runs, a [GPU race-shape lint gate](https://github.com/szaghi/adam/issues/26) checks the FNL kernels.

### Two reference files per case

A full PRISM checkpoint is ~100 MB/rank, far too large to commit. Instead each case commits two compact references, both compared **tolerance-aware**:

| Reference | What it is | Comparison |
|---|---|---|
| `golden/<backend>/digest.txt` | per-variable field digest | `digest.py compare` — `rtol=1e-6, atol=1e-3` |
| `golden/<backend>/*-residuals.dat` | per-iteration residuals log | `digest.py compare-residuals` — same tolerance |

The **field digest** reduces every HDF5 checkpoint to `count, min, max, sum, sum_sq` per field variable, aggregated across all blocks, ranks, and (for forests) realms at each saved iteration. `div_*` diagnostics are excluded (they are ~0 constraint residues, pure round-off). The tolerance is calibrated for cross-compiler runs: `rtol=1e-6` absorbs FMA / libm / MPI-reduction-order noise (~9 significant digits still checked); `atol=1e-3` is the cancellation-residue floor for antisymmetric fields whose `sum` is pure noise. Per-block metadata (`dxdydz`, `origin`, `time_iteration`) that differs only by block partitioning is auto-downgraded to `SKIP_METADATA` — this is what lets the β cross-config oracle compare a two-realm run against a single-realm golden.

::: tip Golden provenance
The **CPU** golden is captured in CI (the single source of truth for the CPU backend); a dev-box golden on a different compiler would drift. The **FNL** golden is captured on the GPU workstation. See the suite README for the exact capture/refresh workflow.
:::

## Case inventory

All cases share `physical_model = EM`, `nv = 9`, `fdv_order = 6`, `ngc = 3`, `constrained_transport = no`, `divergence_correction = no` (so the [issue #11](https://github.com/szaghi/adam/issues/11) `nv`/CT hazard is not exercised).

### Digest/residuals-goldened cases

| Case | Realm | Cadence | Scheme | Exercises |
|---|---|---|---|---|
| `rmf`                  | single | — | fd | baseline EM case (4 AC coils); the digest anchor and the β cross-config reference target |
| `rmf-2realm`           | inter (x-split) | α | fd | single-realm `rmf` split at x=0; aggregated digest vs the `rmf` golden ([#13](https://github.com/szaghi/adam/issues/13) Phase D) |
| `rmf-2realm-asymK`     | inter | α | fd | the asymmetric-K validation: realm_1 SSP-33 (K=3) ∥ realm_2 SSP-54 (K=5), the α K-equality-guard removal ([#16](https://github.com/szaghi/adam/issues/16)) |
| `rmf-2realm-stagesync` | inter | β | fd | the load-bearing β oracle — digest must match the single-realm `rmf` golden bit-for-bit ([#18](https://github.com/szaghi/adam/issues/18)) |
| `rmf-fwl`              | single | — | fd | the only fWLayer case (six faces) + 4 coils — guards the FNL fWLayer host→device transfer fix ([#31](https://github.com/szaghi/adam/issues/31)); **FNL-goldened, CPU golden pending CI** |

### Inter-realm 1:1 seam anchor

`rmf-2realm-fd-pulse` — a source-free Gaussian pulse split at x=0 into two same-resolution realms glued by a **β** mirror seam. It isolates the seam *mechanism* (inter-realm peer copy) from any resolution *jump*: a 1:1 seam must be div-free exactly like a 1:1 intra-block interface, and it holds `div(B) = div(D) = 0` — **but only under β**. Under α the seam ghosts go unfilled during RK substages and div(B) leaks (this is the [#31](https://github.com/szaghi/adam/issues/31) diagnosis). For a 1:1 same-`dt` seam, β is correct and required.

::: warning Not automatically enforced
This case carries **no golden and no `check.sh`**, so `run.sh` skips it and nothing in CI or the local sweep checks the div-free property it documents. It is a **manual reproducer** for [#31](https://github.com/szaghi/adam/issues/31), not an enforced anchor. To make it a real guard it needs either a committed digest golden (both backends) or a `check.sh` asserting `max|div(B)|` and `max|div(D)|` at round-off, in the style of the AMR-seam oracles.
:::

### AMR-seam cases (digest-goldened **and** `check.sh`-driven)

Three single-realm cases carry a static intra-realm **2:1 AMR jump** (`markers_number = 1`, an AMR_GEO box marker covering the x<0 half). They are doubly covered: digest + residuals goldens on **both** backends (captured in [#24](https://github.com/szaghi/adam/issues/24), `d15fed4c`), which `run.sh` checks like any other case, **plus** a bespoke `check.sh` asserting a physics oracle with a marker-off control leg for contrast. The oracles are *not* invoked by `run.sh` — run them directly (they are the only thing that checks the seam `div(B)` behaviour).

| Case | Scheme | `check.sh` oracle |
|---|---|---|
| `rmf-amr`          | `fv_centered` | refinement happens → 2:1 jump exists; the AMR seam registers (`SEAM_KIND_INTRA_REALM_AMR`); Berger-Colella reflux fires non-zero; np1 ≡ np2 register parity ([#28](https://github.com/szaghi/adam/issues/28)) |
| `rmf-amr-fd`       | `fd_centered` | control `max\|div(B)\| ≤ 1e-13` (interior identity intact); seam `max\|div(B)\|` within ±5% of a pinned baseline; div(J) truthfulness within a factor-3 band ([#26](https://github.com/szaghi/adam/issues/26)) |
| `rmf-amr-fd-pulse` | `fd_centered` | source-free two-invariant oracle (`div(D)` structural-zero, `div(B)` seam baseline ±5%); a `--convergence` `p_obs ≥ 0.8` refinement ladder; and the RK-contract legs — leap-splitting refusal + SSP-11 instrument ([#25](https://github.com/szaghi/adam/issues/25)) |

The pulse case is the valid convergence instrument: with `J = 0` the Maxwell system conserves both discrete divergences, so any seam `div(B)` is attributable to the 2:1 jump alone (no coil `h⁻²` contamination). Its `--convergence` leg is what pins the seam `div(B)` source as O(h^p) truncation-order (see [`div(B)` at 2:1 seams](/guide/forest#the-2-1-seam-is-not-div-free-accept-truncation)).

## Reproducers

`src/tests/prism/regression/reproducers/` holds standalone minimal reproducers (not cases): `fnl_mpi_rndv.sh` for the FNL UCX-rendezvous SIGABRT on WSL ([#12](https://github.com/szaghi/adam/issues/12), the `UCX_RNDV_THRESH=inf` workaround) and a gfortran-16 debug-mode `-fcheck=bounds` false-positive reproducer.

## See also

- [Forest (multi-realm)](/guide/forest) — the seam machinery these cases exercise.
- Suite [`README.md`](https://github.com/szaghi/adam/blob/develop/src/tests/prism/regression/README.md) — operational reference (golden capture, tolerances, adding a case).
- `kernel_inventory.md` — the separate PRISM-FNL OpenACC kernel-class audit ([#10](https://github.com/szaghi/adam/issues/10) Step 3).
