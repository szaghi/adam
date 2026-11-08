# PRISM regression suite

The `src/tests/prism/regression/` suite is the **structural-change regression baseline** for the PRISM Maxwell/EM solver — every step of the [forest-of-trees migration](https://github.com/szaghi/adam/issues/10) must leave it green. Each case runs in roughly a minute per backend and validates one backend against committed reference outputs (or, for the AMR-seam cases, against a physics oracle).

This page is the conceptual overview. The authoritative operational reference — golden-capture workflow, CI-as-authority rules, tolerance rationale, adding a case — is the suite's own [`README.md`](https://github.com/szaghi/adam/blob/develop/src/tests/prism/regression/README.md); it is not duplicated here. To stand the suite up on a new machine, compiler, or GPU architecture, follow the [bring-up tutorial](./prism-regression-tutorial).

::: tip Coverage
The suite has 24 case directories (an `input.ini` each). 15 carry a digest golden for the CPU backend and 16 for FNL
(`rmf-fwl` is FNL-goldened only), so a green sweep reports 15 PASS on CPU and 16 on FNL. A SKIP is a case with no
golden for that backend, or one of the two non-case directories `reproducers/` and `slices/`. The other cases are
checked by their own `check.sh`, `scaling.sh` or `coils.sh`.
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

`run.sh {cpu|fnl} [--no-build] [--varset <name>]` builds the backend, then **auto-discovers** cases: every subdirectory with an `input.ini` is a case. A case with no `golden/<backend>/` directory is **skipped** (it is not yet an anchor); set `REGRESSION_RUN_GOLDENLESS=1` to run it anyway for the initial golden capture. `REGRESSION_CASES="rmf rmf-amr"` restricts the run to the named cases (as in the FLUME harness); a name that is not a case directory stops the harness with exit 2 before the build, so a misspelt filter cannot pass by running nothing. Before any case runs, a [GPU race-shape lint gate](https://github.com/szaghi/adam/issues/26) checks the FNL kernels.

A case fails when its outputs differ from the goldens, and also when its run exits non-zero (a crash or the hard
timeout), since a crash after the last checkpoint is still a crash. The harness's own termination of a run that hung in
`MPI_Finalize` after its end-of-run marker is exempt; `REGRESSION_ALLOW_MPIRUN_FAILURE=1` turns the failure into a
warning on an environment whose crash is known.

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

All cases use `fdv_order = 6`. The `rmf*`, `periodic-pulse-blocks` and `memory-scaling*` cases are `physical_model = EM`; the four `pic-*` cases are `physical_model = PIC`. The `pic-*`, `rmf-open-t0` and `rmf-pec-t0` cases use `ngc = 6`, every other case `ngc = 3`. All but `rmf-hyperbolic` run without divergence correction (`constrained_transport = no`, `divergence_correction = no`); `rmf-hyperbolic` is the only case with hyperbolic (GLM) cleaning, on D and B (`nv = 11`).

### Digest/residuals-goldened cases

| Case | Realm | Cadence | Scheme | Exercises |
|---|---|---|---|---|
| `rmf`                  | single | — | fd | baseline EM case (4 AC coils); the digest anchor and the β cross-config reference target |
| `rmf-2realm`           | inter (x-split) | α | fd | single-realm `rmf` split at x=0; aggregated digest vs the `rmf` golden ([#13](https://github.com/szaghi/adam/issues/13) Phase D) |
| `rmf-2realm-asymK`     | inter | α | fd | the asymmetric-K validation: realm_1 SSP-33 (K=3) ∥ realm_2 SSP-54 (K=5), the α K-equality-guard removal ([#16](https://github.com/szaghi/adam/issues/16)) |
| `rmf-2realm-stagesync` | inter | β | fd | the load-bearing β oracle — digest must match the single-realm `rmf` golden bit-for-bit ([#18](https://github.com/szaghi/adam/issues/18)) |
| `rmf-fwl`              | single | — | fd | the only fWLayer case (six faces) + 4 coils — guards the FNL fWLayer host→device transfer fix ([#31](https://github.com/szaghi/adam/issues/31)); **FNL-goldened, CPU golden pending CI** |
| `rmf-fwl-amr`          | single | — | fd | `rmf-fwl` with an init-time box marker refining the x-minus slab to level 3: the fWLayer ranges rebuilt on the refined blocks ([#95](https://github.com/szaghi/adam/issues/95) P5.1) |
| `rmf-pml-amr`          | single | — | fd | the same grid with a CFS PML instead of the layer: the PML support and its RK storage rebuilt on the refined blocks ([#95](https://github.com/szaghi/adam/issues/95) P5.1) |
| `rmf-hyperbolic`       | single | — | fd | `rmf` with hyperbolic cleaning on D and B (`chi = 1`, damping `alpha = 1e9` 1/s): the GLM state layout, φ and ψ in slots 7 and 8, J in 9..11 ([#61](https://github.com/szaghi/adam/issues/61)). Its two init refusals (`hyperbolic` without constrained transport, `poisson`) are legs of `refusals.sh` |
| `periodic-pulse-blocks` | single | — | fd | a Gaussian pulse along +x in a box periodic on every face, 48³ cells as 64 blocks split 32/32 over the two ranks, 192 steps = one period: library periodicity across blocks and ranks |
| `rmf-pec-t0`           | single | — | fd | the four `rmf` coils carrying their current at t0 inside PEC walls: FLAIL's t0 solves of A give a divergence-free B ([#84](https://github.com/szaghi/adam/issues/84) E3b); its `check.sh` asserts `max\|div B\| h / max\|B\| ≤ 1e-12` and np 1 = np 2 to 1e-12 |
| `rmf-open-t0`          | single | — | fd | the same coils inside exact/open faces (`Neumann`): the free-space potential in the ghosts, summed on rank 0 ([#84](https://github.com/szaghi/adam/issues/84) E3c); same `check.sh` acceptance, away from the faces |
| `pic-pec-t0`           | single | — | fd | the t0 electric field of a PIC plasma box between PEC walls: `div D = rho` up to the checkerboard mode ([#84](https://github.com/szaghi/adam/issues/84) E3b); `check.sh` bound `1e-12 max\|rho\|`, np 1 = np 2 |
| `pic-open-t0`          | single | — | fd | the same plasma inside exact/open faces ([#84](https://github.com/szaghi/adam/issues/84) E3c); `check.sh` bound `1e-11 max\|rho\|` beyond 3 cells of the faces, np 1 = np 2 |

### PIC cases (`check.sh`-driven, goldens in `golden-check/`)

Both run the production PIC path (TSC, conserving current `esirkepov_Q10` with the coherent gather, the t0 FLAIL solves,
FD6, SSP(5,4), fixed seed) through the shared driver `pic_check.sh`, as legs: CPU np 1 on one block for the first 4
steps, FNL np 1, and FNL np 2 on the same grid as 64 blocks. Every leg also asserts Gauss's law at t0 and charge
conservation to round-off. Their goldens sit in `golden-check/`, not `golden/`, because `run.sh` runs every goldened case
at `-np 2`, which the single-rank CPU leg refuses; `run.sh` reports them as SKIP.

| Case | Exercises |
|---|---|
| `pic-periodic`   | uniform plasma periodic on every face, plane wave of one box length, 32×24×24 cells, 40 steps: periodic current fold and particle wrap on three axes |
| `pic-slab-pulse` | a plasma slab 3 c/ω_pe thick, PEC along x and periodic in y, z, `gaussian_derivative_pulse` through the slab, 384×24×24 cells, 1152 steps |

### Capacity probes (no golden)

`memory-scaling` runs a minimal EM case (`it_max = 1`) at `-np 1` and `-np 2` and asserts, through `scaling.sh`, how the
per-rank block capacity `nb` scales: it halves on the CPU (host memory shared by the ranks of a node) and stays constant
on FNL (one device per rank). `memory-scaling-coils` is the same case with the four `rmf` coils; `coils.sh` asserts that
the capacity budget counts 3 block-sized fields per coil.

### Inter-realm 1:1 seam anchor

`rmf-2realm-fd-pulse` — a source-free Gaussian pulse split at x=0 into two same-resolution realms glued by a **β** mirror seam. It isolates the seam *mechanism* (inter-realm peer copy) from any resolution *jump*: a 1:1 seam must be div-free exactly like a 1:1 intra-block interface, and it holds `div(B) = div(D) = 0` — **but only under β**. Under α the seam ghosts go unfilled during RK substages and div(B) leaks (this is the [#31](https://github.com/szaghi/adam/issues/31) diagnosis). For a 1:1 same-`dt` seam, β is correct and required.

It now carries a **`check.sh` oracle** (three legs, `-np 1`, CPU by default, `PRISM_EXE=` to override):

1. **β cadence, as committed** — `max|div(D)|` and `max|div(B)|` must both sit at round-off (`1.0E-13`) on **both** realms. Measured: **exactly `0.0`** for all four quantities.
2. **Health** — both realms reach 100% with no error/abort/NaN and both divergence histories exist.
3. **α cadence, negative control** — rewriting `coupling_cadence` to `end_of_step` must push `div(B)` *above* tolerance. Measured: **`4.56E+02`**, identically on both realms. Without this leg the first two could pass vacuously (a seam that stopped being exercised would also read zero), so it is what gives the oracle teeth. It asserts only *that* α leaks, not a magnitude — α's drift is not a quantity this suite has reason to hold fixed.

::: warning Still not run by `run.sh`
Like the AMR-seam oracles, this `check.sh` is invoked **by hand**, not by the sweep — and the case still has no digest golden, so `run.sh` continues to report it as a SKIP. Run it directly after touching `fill_seam_from_peer_forest`, `post_step_forest`, or any forest seam path.
:::

### FV reflux on a capped last step

`rmf-amr-fv-pulse-timed` ([#38](https://github.com/szaghi/adam/issues/38)) — the source-free pulse of `rmf-amr-fd-pulse` on the FV path, a band x ∈ [-0.03, 0.03] refined 2:1 (two seam planes, both x boundaries coarse) and **time-driven**: `time_max` is 4.5 steps, so PRISM caps the fifth step at half a step. With J = 0 the FV update with reflux keeps the domain integrals of D and B constant, and the Berger-Colella correction must be scaled by the step the update used, not by the forest's uncapped `dt` (fixed in `4c1dd1a0`). Its `check.sh` asserts that the last step is capped, that the reflux diagnostic is printed, and that every integral drifts by at most `1e-9` (relative to the largest D, or B, integral). The y extent is twice that of `rmf-amr-fd-pulse` (12σ each side): with the narrower box the Gaussian tails leak through the y faces and even the grid without seams drifts by 1e-8 per step. The reflux-sensitive integral is ∫Dz (its x-flux is Hy); ∫Bx has no x-flux and cannot see the defect. Measured on CPU and FNL: ∫Dz drifts 1.0e-11 over the four full steps (the seams' disturbance reaching the x boundaries) and 2.2e-11 after the capped one. **Negative control** (the fix reverted): the same drift up to step 4, then **1.4e-7** on the capped step, 6500 times larger; the tolerance sits 45× above the fixed run and 140× below the reverted one.

### Inter-realm misaligned-block anchor

`rmf-2realm-fv-pulse-blocks` ([#51](https://github.com/szaghi/adam/issues/51)) — the 1:1 mirror pulse on the FV path (`fv_centered`, 20 steps) with the seam blocks of the two realms not lined up: one realm on 4 blocks per axis (the same cells), first realm 1 (the register side: each register face holds 2x2 realm-2 faces), then realm 2 (each realm-2 face scattered into 4 register faces). Its `check.sh` asserts that the register closes exactly, `max|F_coarse-F_fine_sum|` = 0 on every face and step (320 and 80 lines), and that the fields match the single-realm FV pulse within `1e-13`. They are not compared bitwise: PRISM FV is not invariant under a change of block layout at round-off (the single realm alone on 4 blocks per axis instead of 2 differs by 6.6e-16). Measured on CPU and FNL: mismatch 0, fields within 3.7e-15.

### Inter-realm 2:1 seam anchor

`rmf-2realm-fd-pulse-refined` ([#52](https://github.com/szaghi/adam/issues/52)) — `rmf-amr-fd-pulse` split at its 2:1 face x=0 into a fine realm (x<0, one level finer, no marker) and a coarse realm (x>0) glued by `coupling = refined`, β cadence. The union is the single-realm cell set, so the split must reproduce `rmf-amr-fd-pulse` **bit for bit**, including its [#29](https://github.com/szaghi/adam/issues/29) seam div(B) source: the inter-realm 2:1 seam reuses the intra-realm formulas, so the source must be neither smaller nor larger. Its `check.sh` (`-np 2` by default, `--np N`, `PRISM_EXE=` to override) runs both and asserts:

1. **Health** — both runs reach 100% with no error/abort/NaN.
2. **Fields** — every checkpoint field (B, D, J, `div_*`, `res_*`) bitwise on every cell (`multirealm_oracle.py --fields-only --tol 0`). Measured: 147456 cells, 27 fields, CPU and FNL.
3. **Divergence histories** — at every step the single-realm `max|div(D)|` and `max|div(B)|` equal the larger of the two realms' values. Measured: `max|div(B)|` 6.86 → 12.65 in both, the coarse realm holding the maximum at steps 1–2 and the fine one at steps 3–5.
4. **Seam div(B) monitor** — `[IO].seam_divB_tol = 1` (warn-only) fires on the split.

The FNL run found a pre-existing deadlock: the FNL residuals skipped the (collective) ghost update on a rank that owns no block of a realm, which the refined split produces at `-np 2`; the ghost update now runs on every rank, as on the CPU. Like the 1:1 anchor, the case has no digest golden and is run by hand.

### AMR-seam cases (digest-goldened **and** `check.sh`-driven)

Three single-realm cases carry a static intra-realm **2:1 AMR jump** (`markers_number = 1`, an AMR_GEO box marker covering the x<0 half). They are doubly covered: digest + residuals goldens on **both** backends (captured in [#24](https://github.com/szaghi/adam/issues/24), `d15fed4c`), which `run.sh` checks like any other case, **plus** a bespoke `check.sh` asserting a physics oracle with a marker-off control leg for contrast. The oracles are *not* invoked by `run.sh` — run them directly (they are the only thing that checks the seam `div(B)` behaviour).

| Case | Scheme | `check.sh` oracle |
|---|---|---|
| `rmf-amr`          | `fv_centered` | refinement happens → 2:1 jump exists; the AMR seam registers (`SEAM_KIND_INTRA_REALM_AMR`); Berger-Colella reflux fires non-zero; np1 ≡ np2 register parity ([#28](https://github.com/szaghi/adam/issues/28)) |
| `rmf-amr-fd`       | `fd_centered` | control `max\|div(B)\| ≤ 1e-13` (interior identity intact); seam `max\|div(B)\|` within ±5% of a pinned baseline; div(J) truthfulness within a factor-3 band ([#26](https://github.com/szaghi/adam/issues/26)) |
| `rmf-amr-fd-pulse` | `fd_centered` | source-free two-invariant oracle (`div(D)` structural-zero, `div(B)` seam baseline ±5%); a `--convergence` `p_obs ≥ 0.8` refinement ladder; and the RK-contract legs — leap-splitting refusal + SSP-11 instrument ([#25](https://github.com/szaghi/adam/issues/25)) |

The pulse case is the valid convergence instrument: with `J = 0` the Maxwell system conserves both discrete divergences, so any seam `div(B)` is attributable to the 2:1 jump alone (no coil `h⁻²` contamination). Its `--convergence` leg is what pins the seam `div(B)` source as O(h^p) truncation-order (see [`div(B)` at 2:1 seams](/concepts/forest#the-2-1-seam-is-not-div-free-accept-truncation)).

## Reproducers

`src/tests/prism/regression/reproducers/` holds standalone minimal reproducers (not cases): `fnl_mpi_rndv.sh` for the FNL UCX-rendezvous SIGABRT on WSL ([#12](https://github.com/szaghi/adam/issues/12), the `UCX_RNDV_THRESH=inf` workaround) and a gfortran-16 debug-mode `-fcheck=bounds` false-positive reproducer.

## See also

- [Forest (multi-realm)](/concepts/forest) — the seam machinery these cases exercise.
- Suite [`README.md`](https://github.com/szaghi/adam/blob/develop/src/tests/prism/regression/README.md) — operational reference (golden capture, tolerances, adding a case).
- [`kernel_inventory.md`](https://github.com/szaghi/adam/blob/develop/src/tests/prism/regression/kernel_inventory.md) — the separate PRISM-FNL OpenACC kernel-class audit ([#10](https://github.com/szaghi/adam/issues/10) Step 3).
