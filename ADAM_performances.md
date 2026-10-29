# ADAM / FLUME — the FNL face-flux kernel on MI300X

Investigation of `compute_face_fluxes_dev`, the dominant kernel of the FLUME OpenMP-offload (`fnl-amd`) backend,
on MI300X. Analysis October 2026; **all experiments measured 2026-10-27**, Slurm job 73669 on `TheraC18`.

- Case: `src/tests/flume/speedup/profile-couette3d` (3-D compressible Couette, WENO-5, SSP-RK5)
- Commit for every build and run: `8cb8d339`; no build flagged `src/ modified` except where stated
- AMD records: `results/runs-amd.csv`, run dirs `runs/`, counters `runs/counters-face/`
- NVIDIA comparison: medooza, one A30, `fnl-nvf` leg + `lto` / `r128` / `B0-lto` (jobs 1357614–1357635)
- Builds compared: **leg** (unmodified), **`noagpr`**, **`smax3`**, and the earlier **`B0`**

---

## 0. Executive summary

1. The face kernel is **5× slower per face on an MI300X than on an A30**, on a device with ~16× the FP64 peak.
   It is latency-bound at **1 wave/SIMD**, running at **0.216 % of vector peak**, with four of every five vector
   instructions being address/scratch shuffling rather than arithmetic.
2. The original diagnosis — *large private state → registers → occupancy 1* — is **wrong in its mechanism**.
   Across three builds the registers never moved (428 VGPR / 172 AGPR / 0 VGPR spills) and occupancy never left
   1 wave/SIMD, yet the kernel got **21.6 % faster**. The private arrays live in **scratch**, not registers.
3. The prescription — *cut per-thread private state* — **survives, and is confirmed by measurement**. The right
   ranking metric is **bytes of scratch removed per thread**, not registers. Measured rate: roughly
   **−20 % for −900 B**.
4. **There is no flag-level register or occupancy lever on this toolchain.** `-mllvm` on `amdflang` reaches the
   *host* pipeline only; every device-only channel was rejected. The remaining routes are all source changes.
5. **Recommended route: make the WENO order a build-time constant** (§6). ~5 files, keeps the runtime path as the
   default, captures both halves of the win, and is justified by the fact that **138 of 138 committed FLUME
   inputs use `weno-u-5`**.
6. **Open question:** what holds the 428 registers. Not the stencil arrays. Until answered, no variant of this
   work should be sold as an occupancy fix.

**→ What to do next is in [Call for action](#call-for-action--the-next-three-steps): one measurement run (A1),
one small 3-file code probe (A2), and only then the real implementation (A3).**

---

## Call for action — the next three steps

Ordered so that each step is cheap relative to the one it de-risks, and each has a stated decision criterion.
**Nothing below has been run or implemented.**

> **Update 2026-10-27 (job 73716):** A1 and A2 have since been run — results, code and environment notes in [§10](#10-results-of-a1-and-a2-2026-10-27-slurm-job-73716-on-therac18). A3 has not been started.

### A1 — Measurement, no code change: do the two levels compose? (~45 min, one build)

`smax3` (declarations, −21.6 %) and `B0` (loop bounds, −15.4 %) attack different things. Nobody has built them
together, so it is unknown whether the gains add, overlap, or saturate. This single number decides whether §6
needs **both** levels or only level 1 — i.e. whether the macro has to reach into the kernel loop bounds at all.

```bash
# from the repository root
git apply src/tests/flume/speedup/profile-couette3d/patches/face-kernel-s3.patch \
          src/tests/flume/speedup/profile-couette3d/patches/stencil-smax3.patch
cd src/tests/flume/speedup && ./build.sh fnl-amd --varset "thera_amd gfx942" --jobs 32 --tag s3b0
git -C ../../.. apply -R src/tests/flume/speedup/profile-couette3d/patches/face-kernel-s3.patch \
                         src/tests/flume/speedup/profile-couette3d/patches/stencil-smax3.patch
cd profile-couette3d
./profile.sh fnl-amd --size small --np 1 --tool rocprof --tag s3b0      # registers + face kernel
./profile.sh fnl-amd --size large --np 1 --tool none    --tag s3b0      # the number that matters
./profile.sh fnl-amd --size small --np 1 --steps 10 --io --tool none --tag s3b0   # then h5diff -c vs leg
```

Compare the `large` loop against **`smax3` 588.27 s** and **leg 750.36 s** (both `--tool none`).

| Outcome | Reading | Consequence |
|---|---|---|
| materially below 588 s | the levels compose | §6 must implement **both** levels — the `WENO_S` macro as sketched |
| ≈ 588 s | level 2 adds nothing once declarations shrink | §6 needs **only** `ADAM_S_MAX` + the guard: much smaller change |
| above 588 s | they interfere | investigate before building anything on either |

Optional, for a clean three-way: re-measure `B0` alone with `--tool none` (its 629.25 s was taken under
rocprof, so it is not directly comparable to the others). One more `large` run, ~13 min.

### A2 — Small code change to test the big one: a "lean" Euler face kernel (~3 files, one build)

Route #1 removes **1040 B** of scratch — more than the 896 B that bought −21.6 % — but costs ~110 files. **Do not
write that yet.** The same hypothesis can be tested on the Euler path alone, which is the only path this case
runs, by letting the Euler module include its *own* face-kernel body while the other twelve models keep the
shared one:

1. **`src/app/flume/common/adam_flume_euler_library_body.INC`** — add, alongside the existing routine (purely
   additive; every current caller including the whole CPU backend is untouched):
   ```fortran
   compute_face_split_fluxes_lean(gamma, ch, d, S, is_characteristic, &
                                  q_gpu, q_aux_gpu, ngc, b, i, j, k, di, dj, dk, fsplit, er, mu)
   ```
   Identical body, except it reads `q_m(NV_K)` and `qa_m(NV_AUX_K)` one stencil point at a time straight from
   `q_gpu` / `q_aux_gpu` instead of receiving gathered `qs` / `qas`. It already walks `m` twice (once for
   `alpha`, once for the split), so the stencil is re-read from L2 twice — the stated cost of #1. **No `R`
   argument and no recomputation of the auxiliary variables are needed**, so the arithmetic is untouched and the
   result is bitwise identical by construction.
2. **new `src/app/flume/fnl/adam_flume_fnl_face_kernels_euler_lean.INC`** — a copy of the shared face-kernel
   body with the `qs` / `qas` declarations and the gather loop deleted, calling the lean split routine.
   Keep `add_gravity_edge_correction_dev` unchanged in the copy.
3. **`src/app/flume/fnl/adam_flume_fnl_euler_kernels.F90`** — include the lean `.INC` instead of the shared one.

Expected: scratch **6048 → ~5008 B (−1040 B)**, so by the measured rate (≈ −20 % per −900 B) a loop around or
below `smax3`'s 588 s — and, unlike `smax3`, **legitimately**, with no cap on the WENO order.

| Outcome | Consequence |
|---|---|
| ≥ `smax3`'s gain | #1 is confirmed at the predicted rate → roll out to the other 12 models, and combine with §6 (they touch different arrays) |
| materially less | the scratch-bytes model is incomplete → re-open §9.1 before spending ~110 files |
| no gain | #1 is dead; §6 alone is the route |

The duplicated kernel body is deliberate and temporary — the shared `.INC` exists precisely to avoid it, and the
rollout in A3 deletes the copy.

### A3 — Only after A1 and A2: implement for real

Whichever of §6 and #1 the two probes justify, implemented properly rather than as a patch:

- **§6**: `ADAM_S_MAX` macro + the init-time guard (`error stop` when the input's `S` exceeds the compiled bound),
  plus the `WENO_S` macro if and only if A1 says level 2 pays. Default undefined → today's behaviour bit for bit.
  First verify the three items in §6 "Checks required before implementing" — in particular PRISM's `weno-u-1`
  and `weno-c-4`, and that NASTO's ROR path keeps `S_MAX` as a true runtime bound.
- **#1**: the lean split rolled out from Euler to the remaining 12 models (~110 files, mechanical signature
  changes), deleting the temporary duplicate body from A2.

Gate for both: the `--steps 10 --io` + `h5diff -c` check, then the **FLUME and PRISM FNL regression suites** —
neither probe above is a substitute for those, since the profiling case's step-10 regrid moves no data.

### Also outstanding

- **Multi-rank**: everything here is one rank on one GPU. `--np 8` on `large` would show whether the kernel
  saving survives alongside communication (§9.5).
- **The register mystery** (§9.1) and **why `noagpr` is faster at all** (§9.2) remain unexplained. `noagpr`
  should **not** be adopted on its current evidence.
- Delete the stray `adam_flume_fnl_euler_kernels.mod` from the repository root (§8, trap 3).

---

## 1. Measured results

All runs: one rank, one MI300X, commit `8cb8d339`. `small` numbers are under rocprofv3, **alternated** between
variants to remove ordering and thermal effects; `large` numbers are `--tool none`.

### Face kernel, `small` (300 dispatches, workgroup 256, grid 311296 — identical in every build)

| Build | Face kernel, mean per dispatch | vs leg |
|---|---|---|
| leg | 62.33 (job 73623), 62.26, 61.29, 61.39, 61.32, 61.37 ms | — |
| `noagpr` | 51.56, 51.87, 51.86 ms | **−16 %** |
| **`smax3`** | **48.79, 48.79 ms** | **−20.4 %** |

Per-dispatch spreads do not overlap (leg 60.64–62.67, `noagpr` 51.25–52.50 ms).

### Forest loop

| Build | `small` loop | `large` loop | `large` wall |
|---|---|---|---|
| leg | 21.16, 21.03, 21.06 s | 751.92 (job 73623), 744.66, 754.51 s — **mean 750.36 s, spread 1.3 %** | 767.05, 776.72 s |
| `noagpr` | 18.09 s (−14.5 %) | 633.31 s (**−15.6 %**) | 655.61 s |
| **`smax3`** | **17.21, 17.23 s (−18.2 %)** | **588.27 s (−21.6 %)** | **610.61 s** |
| `B0` (earlier, under rocprof) | — | 629.25 s vs leg-under-rocprof 743.83 s (−15.4 %) | 651.99 s |

The leg's 1.3 % spread across three `large` runs is well clear of both effects.

### Kernel resource records — the decisive negative

From `llvm-readelf --notes` on the device ELF (authoritative; see §8 for the extraction recipe). Euler
`compute_face_fluxes_dev`:

| Field | leg | `noagpr` | `smax3` |
|---|---|---|---|
| `vgpr_count` (unified total) | 428 | 428 | 428 |
| `agpr_count` | 172 | 172 | 172 |
| `sgpr_count` | 108 | 108 | 108 |
| `vgpr_spill_count` | 0 | 0 | 0 |
| `sgpr_spill_count` | 151 | 151 | 151 |
| **occupancy** `floor(512/428)` | **1 wave/SIMD** | **1** | **1** |
| `private_segment_fixed_size` (scratch) | 6048 B | 6048 B | **5152 B** |

### Correctness

`small`, refined walls, 10 steps, checkpoints at steps 0 and 10 (`--steps 10 --io`), `h5diff -c` against the leg:
**bitwise identical, both steps, for `noagpr` and `smax3`.**

Weak check, as in the A30 campaign: the step-10 regrid leaves the grid unchanged (geometric markers), so no
regrid data transfer is exercised. The FLUME and PRISM FNL regression suites remain the gate before adoption.

### Build times

leg 995 s, `B0` 996 s, `noagpr` 994 s, `smax3` 991 s (32 jobs, device LTO is the default on AMD and dominates
the tail of the build).

---

## 2. Why the kernel is slow — what holds up

### Hardware counters (`runs/counters-face/`, `small`, 3 dispatches)

| Counter | Value | Meaning |
|---|---|---|
| `MeanOccupancyPerActiveCU` | 1.0 | 1 wave per SIMD |
| `OccupancyPercent` | 11.8 % | out of 8 waves/SIMD max — the minimum |
| `SQ_WAVES` | 4864 | = 1216 workgroups × 4 waves, matches the grid |
| `SQ_INSTS_VALU` | 551.7 M | total vector instructions |
| FP64 arithmetic (FMA+ADD+MUL) | 114.1 M | only 20.7 % of VALU |
| Achieved FP64 | 176 GFLOP/s | 0.216 % of the 81.7 TFLOP/s vector peak |

1. **Four of five vector instructions are not arithmetic** — ~437 M of 551.7 M are address arithmetic and
   register/scratch shuffling, consistent with private arrays in scratch, dynamically indexed.
2. **The VALUs are idle ~98.6 %** — 551.7 M instructions in 61.3 ms is 9.0 G instr/s against ~638 G instr/s of
   capacity (1216 SIMDs, one wave64 VALU op per 4 cycles at ~2.1 GHz). Latency-bound, not compute-bound, with no
   second resident wave to cover a stall.

### Cross-vendor framing: 5× slower per face than an A30

Normalised per face; the two AMD sizes agree to 1 %, so it is not a size artefact:

| Machine / build | Case | Face kernel per launch | Faces | **ns per face** |
|---|---|---|---|---|
| MI300X, leg | `large` refined, 46 M | 2410.5 ms | 46 M | **52.4** |
| MI300X, leg | `small` refined, 1.2 M | 62.3 ms | 1.2 M | **51.9** |
| A30, `fnl-nvf` leg | `medium` uniform, 16.8 M | 172.4 ms | 16.8 M | **10.3** |
| A30, `fnl-nvf` `lto` | `medium` uniform | 103.1 ms | 16.8 M | **6.1** |

**5× slower than an A30, 8.5× slower than the A30's LTO build**, on ~16× the FP64 peak. Not an algorithmic cost:
the private state is handled far worse by this backend. On MI300X the face kernel is **97 % of GPU time**
(723 s of 737 s on `large`) against 67–71 % on the A30.

### Corroboration from the NVIDIA side

Nsight Compute on WSL (`small` uniform, one launch) of the A30 leg build: the kernel writes **5.75 GB to DRAM for
89 MB of fluxes** — the rest is local memory — and executes 139 M local loads/stores against 33 M global loads.
`lto` halves the local traffic; `B0-lto` cuts it 83 % with 65 % fewer instructions. Same private-state story,
different tool, different vendor.

---

## 3. The private state budget

`compute_face_fluxes_dev`, Euler model (`NV_K = 5`, `NV_AUX_K = 8`), doubles per thread. Sources:
`src/app/flume/fnl/adam_flume_fnl_face_kernels_agnostic.INC`,
`src/app/flume/common/adam_flume_euler_library_body.INC` (`compute_face_split_fluxes`),
`src/lib/fnl/adam_fnl_weno_kernels.F90`.

| Array | Declared in | `S_MAX = 5` | `S_MAX = 3` | Δ |
|---|---|---|---|---|
| `fsplit(2,1-S_MAX:S_MAX-1,NV_K)` | face kernel `:36` | 90 | 50 | −40 |
| `qas(NV_AUX_K,1-S_MAX:S_MAX)` | face kernel `:35` | 80 | 48 | −32 |
| `qs(NV_K,1-S_MAX:S_MAX)` | face kernel `:34` | 50 | 30 | −20 |
| `VP, w, IS, a` — each `(1:2,0:S_MAX_M1)` | `adam_fnl_weno_kernels.F90:30,31,120,121` | 40 | 24 | −16 |
| `v(2,2*S_MAX-1)` | face kernel `:39` | 18 | 10 | −8 |
| *S-dependent subtotal* | | *278* | *162* | **−116 (−928 B)** |
| `el(5,5)` + `er(5,5)` | split + kernel | 50 | 50 | — |
| `roe, lambda, alpha, f, qm, mu, vr, vr_, flux` | both | ~63 | ~63 | — |
| **peak live** | | **~391** | **~275** | |

**Predicted −928 B against a measured −896 B** (6048 → 5152). The 32 B residual is frame packing/alignment. The
array model is therefore sound — which is exactly what makes the register result below so informative.

The peak live set is **inside `compute_face_split_fluxes`** (`qs` 50 + `qas` 80 + `fsplit` 90 + `el` 25 + `er` 25
+ scalars ≈ 303), not in the WENO reconstruction loop (≈192).

`add_gravity_edge_correction_dev` shares the same body and declarations; every change below applies to it too.

---

## 4. The correction: scratch traffic, not occupancy

**Superseded reasoning** (recorded because it drove the experiments): gfx942 has a unified 512-register budget,
so 428 registers → `floor(512/428) = 1` wave/SIMD, matching the measured 11.8 %. It followed that halving the
live set to ≤256 would buy 2 waves/SIMD, and that the AGPRs were a spill area the compiler fell back on.

**What the measurements show instead:**

- `smax3` removed **896 B of scratch** — 38 % of the S-dependent private state — and **not one register moved**.
  `vgpr_count`, `agpr_count`, `sgpr_count` and both spill counts are identical to the leg. Occupancy stayed at
  1 wave/SIMD. The kernel still got **21.6 % faster**.
- `vgpr_spill_count` is **0** in every build, so the 172 AGPRs are *not* a VGPR spill area. The real spill story
  is **151 SGPR spills** plus the 6048 B of scratch, and the scratch is the private-array traffic of §3.
- `noagpr` did not touch the face kernel at all (identical record, identical disassembly: 2694 instructions,
  260 scratch accesses, 220 `s_waitcnt`), yet was 16 % faster — see §5.

**Conclusion.** The private arrays live in scratch, not registers. Shrinking them cuts **scratch traffic**, which
is what the speedup buys. Occupancy is not the thing to optimise, and the register count is not the metric to
rank changes by — **bytes of scratch removed per thread** is. Measured rate: **≈ −20 % per −900 B**.

### A caution against over-reading even the scratch story

`add_dissipative_fluxes_o4_dev` carries **21744 B** of scratch — 3.6× the face kernel's — with the same register
record, and runs **70× faster** per launch (34.1 ms against 2410.5 ms on `large`). Scratch *size* alone does not
predict cost. What makes the face kernel pathological is scratch traffic **inside a long dependent chain at one
resident wave**; the dissipative kernel has far less work and far less dependent reuse per thread.

---

## 5. Experiment `noagpr`: a real speedup of unknown mechanism

`-mllvm -amdgpu-spill-vgpr-to-agpr=0`, clean build, 994 s. Designed to discriminate "occupancy binds" from
"spill traffic binds" by dropping the face kernel's register count.

**It never dropped it.** The flag changed the device image (17,608,120 → 17,588,672 bytes) and **281 of 383
kernels** — all small ones, in the pattern `agpr_count 2 → 0`, `vgpr_count 118 → 114` — but the face kernel was
not among them. Same ISA, same resources, same grid, same work, and yet reproducibly **−16 % / −15.6 %**.

So: the experiment **did not test its hypothesis**, and the speedup is **real but unexplained**. Remaining
candidates are second-order — different host code and allocation order changing HBM page or scratch-pool
placement for a latency-bound kernel. Treat it as **a measured speedup of unknown mechanism, not a fix**: it
could evaporate on another size or another node. Do not adopt it without understanding it.

### Why it could not have worked: no device-flag channel exists

`-mllvm` on `amdflang` feeds the **host** pipeline. `-amdgpu-spill-vgpr-to-agpr=0` is merely *accepted* there
because the AMDGPU backend registers its `cl::opt` globally. Every device-only channel was rejected:

| Attempt | Result |
|---|---|
| `-mllvm -amdgpu-spill-vgpr-to-agpr=0` | accepted (host pipeline; no effect on the face kernel) |
| `-mllvm -amdgpu-waves-per-eu=2` | `Unknown command line argument` — identical to a bogus control flag |
| `-mllvm -amdgpu-num-vgpr=256` | `Unknown command line argument` |
| `-mllvm -amdgpu-function-attrs` | `Unknown command line argument` |
| `-Xopenmp-target=amdgcn-amd-amdhsa -mllvm …` | `flang-23: error: unknown argument` |
| `-Xoffload-linker --device-compiler=-mllvm=…` | `ld.lld: error: unknown argument` |

`amdgpu-waves-per-eu` is a **function attribute**, not a `cl::opt`, so there is no AMD analogue of nvfortran's
`maxregcount` as a flag — it would need `ompx_attribute(__launch_bounds__(…))` on the `target` directive, i.e. a
source change in FUNDAL's `OMPLOOP` (`src/third_party/FUNDAL/src/lib/fundal.H`, which emits a bare
`target teams distribute parallel do` with no clauses).

For the record, on the A30 the equivalent cap made things **worse**: `r128` (`-gpu=maxregcount:128`) cost the
face kernel 16 % (200.4 ms against 172.4 ms) because the spills it added outweighed the warps it gained.

---

## 6. Recommended route: make the WENO order a build-time constant

### The evidence for it

- **138 of 138 committed FLUME inputs use `weno-u-5` (S = 3).** None uses `weno-u-7` or `weno-u-9`. `S_MAX = 5`
  is a ~40 % private-state tax paid for an order the test suite never exercises.
- In FLUME, the only `S` reaching the reconstruction is the single `weno%S`, fixed at initialisation. The
  runtime-varying `S` in the library (**ROR**, the cell-by-cell reduced-order fallback, `ror_schemes_gpu(r)`)
  exists only in **NASTO** and **CHASE** — both unmaintained. No FLUME or PRISM path uses it.
- `smax3` already measured what level 1 buys: **−21.6 % on `large`**, bitwise identical.
- `B0` already measured what level 2 buys: **−15.4 %** on `large` (under rocprof), independently.

### Why the `smax3` patch itself is unsafe, and the fix

`patches/stencil-smax3.patch` shrinks the declarations but leaves `S` a runtime value, so a `weno-u-7` input
would loop `do m=1-4,4` over arrays declared `1-3:3` and corrupt memory **silently on the device** — no bounds
check fires there. It is valid only because this case runs WENO-5.

Make `S_MAX` a build macro and add an init-time guard, and the failure mode inverts from silently-wrong to
refused-at-init — the pattern the codebase already uses for `poisson`, unknown `PML_type` and `radiative`:

```fortran
! src/lib/common/adam_weno_object.F90
#ifndef ADAM_S_MAX
#define ADAM_S_MAX 5
#endif
integer(I4P), parameter :: S_MAX    = ADAM_S_MAX
integer(I4P), parameter :: S_MAX_M1 = ADAM_S_MAX - 1
```
```fortran
! in weno%initialize, right after the scheme is parsed
if (self%S > S_MAX) error stop 'weno: scheme '//trim(scheme)//' needs S=... but this build has S_MAX=... (-DADAM_S_MAX)'
```

Default stays 5, so a build that does not opt in is unchanged bit for bit.

### Two levels — the second is where the rest of the win is

| Level | What becomes compile-time | Measured independently |
|---|---|---|
| **1 — declarations** | private array *bounds* (`S_MAX`) | `smax3`: **−21.6 % on `large`** |
| **2 — loop bounds** | the stencil loops use the parameter, not the dummy `S` | `B0`: **−15.4 % on `large`**; −67 % on A30 under LTO |

Level 2 is only *legitimate* once the order is fixed at build time. Today `B0` is a cheat — it hard-codes
`SC = 3` while the input may still say otherwise. With the guard above it becomes correct. One macro gets both
from a single source, with no duplicated kernels:

```fortran
#ifdef ADAM_S
   integer(I4P), parameter :: SC = ADAM_S     ! fixed-order build
#  define WENO_S SC
#else
#  define WENO_S S                            ! runtime order: today's behaviour
#endif
```
Kernels then write `do m=1-WENO_S, WENO_S`. Macro undefined → current code, unchanged.

### Cost and alternatives

The cost is **one executable per WENO order**. Against the alternative — one kernel per `S` generated from the
same `.INC`, dispatched on the host with a `select case (S)` (what the case README calls "making it a feature") —
that is **4× the device code in 13 model modules** on top of a build that already takes ~994 s. Given 138/138
inputs at a single order, the build-time constant is far cheaper, and it is **~5 files against ~110** for the
`qs`/`qas` refactor below.

### Checks required before implementing

- **PRISM uses `weno-u-1` (S = 1) ×160 and `weno-c-4` ×6.** S = 1 would shrink the arrays dramatically —
  potentially a larger win than FLUME's — but it is a different build mode and the **centred** schemes may not
  size through this path. Verify before assuming.
- **NASTO genuinely needs runtime-varying `S`** (ROR). The shared library must keep `S_MAX` as a true bound for
  it; NASTO builds simply would not define the macro. (NASTO is unmaintained in any case.)
- Confirm the `weno-c-*` and interpolation paths share these declarations.

---

## 7. The other source routes, re-assessed

Originally ranked by registers saved; **re-ranked here by scratch bytes removed**, which is what the
measurements say matters. All are bitwise-identical by construction — they move or recompute values and reorder
no floating-point operation — so the `--steps 10 --io` + `h5diff` gate applies, with the FNL regression suites
before adoption.

| # | Change | Scratch removed | Cost | Status |
|---|---|---|---|---|
| **6** | **Build-time WENO order** (§6) | **928 B** (level 1) + level 2 | ~5 files; one exe per order | **recommended** |
| **1** | **Stop materialising `qs` and `qas`.** They exist only to feed `compute_face_split_fluxes`; give it a device variant taking `q_gpu, q_aux_gpu, b, i, j, k, di, dj, dk` and scalar `q_m(5)`, `qa_m(8)` buffers. It already walks `m` twice, so the stencil is re-read from L2 twice instead of 1040 B being held per thread. | **1040 B** | **~110 files**: 2 library bodies, 4 agnostic `.INC`, ~104 per-model `face_split_fluxes` adapters | viable, expensive |
| **3** | **Drop `qas` alone.** `compute_q_aux_dev` shows `q_aux` is exactly `conservative_to_auxiliary(q)` on every cell *including ghosts*, so recomputing per stencil point is bitwise identical: ~12 extra calls per face (one divide, one sqrt) against 640 B, at 0.216 % of FP64 peak. | 640 B | subset of #1 | viable |
| **4** | **Eliminate `fsplit`** by fusing the split with the reconstruction field by field (`k` outer, writing straight into `v(2,2S-1)`), caching `f_s` or recomputing `compute_flux` per field. | ~320 B net | moderate | viable |
| **2** | **Per-`S` kernel specialisation**, host `select case` | 928 B (same as #6 level 1) | 4× kernel code, longer builds | superseded by #6 |
| ~~**5**~~ | ~~Defer `er` past the WENO loop~~ | — | — | **withdrawn** |

**Why #5 is withdrawn:** it assumed the peak live set is during the WENO loop. It is not — it is inside
`compute_face_split_fluxes` (§3), and `compute_eigenvectors` produces `el` and `er` *together*, so `er` is live
at the peak regardless of when it is consumed. It removes nothing from the binding constraint.

#6 and #1 are complementary (different arrays), but #6 should be measured first: it is 20× less work and
already has its level-1 number.

---

## 8. Running on thera

### Environment

The toolchain is **not** the `rocm-afar-drop` of the system module tree; it is mbycklin's spack instance, exposed
by `spack-setup-env.sh`, which resets `MODULEPATH` to `/home/mbycklin/modules/Core` (the line exists but is
commented out at `~/.bashrc:77`):

```bash
source /home/mbycklin/bin/spack-setup-env.sh
module load rocm-afar-drop/23.2.1-gfx942 rocm-openmpi rocm-openblas rocm-hdf5-ompi-5.0.10/1.14.6
export PATH=/opt/rocm-7.2.1/bin:$PATH    # rocprofv3 1.1.0; the drop's 1.2.3 aborts this binary (case README)
```

This reproduces the leg's recorded toolchain exactly: `amdflang` from
`/home/mbycklin/amd/therock-23.2.1-gfx94X-7.13.0-7357b5084b/bin` ("AMD AFAR drop #23.2.0 04/18/26 flang version
23.0.0git … 35849413"), `mpif90`/`mpirun` from the spack `openmpi-5.0.10-3kcys5q…` prefix — the paths the
executable carries in its RPATH.

### Four traps, each of which cost a round trip

1. **The reservation carried no GPU gres.** Job 73669 had `AllocTRES=cpu=256,mem=…` and no `gpu=`. Steps saw no
   device (`rocm-smi`: "No AMD GPUs specified") and the application **segfaulted in `field_object%initialize`
   at start-up, exit 139**. The *leg* executable crashed identically — that is what identified it as
   environmental rather than a regression. Every step needs:
   ```bash
   srun --jobid=<id> --mpi=none --gres=gpu:8 -n1 bash -lc '…'
   ```
   `--mpi=none` is required too: the default `pmix_v6` plugin fails to load for a plain command.
2. **Never pipe `module load`.** It is a shell function; `module load X | tail` runs it in a subshell and the
   environment is silently lost, with no error.
3. **A stale `.mod` in the repository root.** An untracked `adam_flume_fnl_euler_kernels.mod` (left by the B0
   build) made the first `smax3` build fail with `Cannot use module file … is not the right module file` for
   every module the Euler kernels import. `fobis clean` wipes `exe/mod/` but not the repo root, and the compiler
   searches `.` first. It only bites when a change alters a module interface — which is why `noagpr` was
   unaffected and `smax3` was not. Moved to `~/scratch-flagtest/`; **it should be deleted from the tree.**
4. **The `free/total memory` log lines are not a device check.** They print free 1.58 GB against total 579 MB on
   a 192 GB MI300X, in both working and failing runs. Use `rocm-smi`.

### `build.sh --flags`

`build.sh` previously refused extra flags on any leg but `fnl-nvf` (`--gpu` prepends nvfortran's `-gpu=`
umbrella). It now also takes **`--flags "OPTS"`** — compiler-agnostic, appended verbatim to the compile and link
lines, same guards as `--gpu` (needs `--tag`, needs a clean build). This is the AMD flag-variant path the case
README recorded as missing.

### Reading kernel resources without running anything

`-Rpass-analysis=kernel-resource-usage` works and prints TotalSGPRs, VGPRs, AGPRs, ScratchSize,
`Occupancy [waves/SIMD]` and spill counts per kernel — but because device LTO is the default on AMD, the device
ISA is produced by the linker wrapper, so the remark appears at **link** time. A `-c` compile of one module
prints nothing.

For an already-linked executable, the authoritative record is the kernel metadata:

```bash
D=/home/mbycklin/amd/therock-23.2.1-gfx94X-7.13.0-7357b5084b
$D/llvm/bin/llvm-objcopy --dump-section=.llvm.offloading=img <exe> /dev/null   # ~17.6 MB device image
python3 -c "d=open('img','rb').read(); i=d.find(b'\x7fELF'); open('dev.elf','wb').write(d[i:])"
$D/llvm/bin/llvm-readelf --notes dev.elf                                       # .vgpr_count, .agpr_count, …
$D/llvm/bin/llvm-objdump -d dev.elf                                            # per-kernel ISA
```
Split the notes on `.agpr_count:` to get one record per kernel — a window-based search around a kernel name
picks up its neighbours' fields and will mislead you.

**rocprofv3's `VGPR_Count = 128` / `Accum_VGPR_Count = 304` do not map onto the metadata fields** (`vgpr_count`
428 / `agpr_count` 172). Prefer the metadata.

### Reproducing

```bash
cd src/tests/flume/speedup
./build.sh fnl-amd --varset "thera_amd gfx942" --jobs 32                     # the leg
./build.sh fnl-amd --varset "thera_amd gfx942" --jobs 32 --tag noagpr \
           --flags "-mllvm -amdgpu-spill-vgpr-to-agpr=0"
git apply profile-couette3d/patches/stencil-smax3.patch                      # from the repo root
./build.sh fnl-amd --varset "thera_amd gfx942" --jobs 32 --tag smax3
git apply -R profile-couette3d/patches/stencil-smax3.patch

cd profile-couette3d
./profile.sh fnl-amd --size small --np 1 --tool rocprof [--tag T]            # registers + kernel times
./profile.sh fnl-amd --size large --np 1 --tool none    [--tag T]            # loop timing
./profile.sh fnl-amd --size small --np 1 --steps 10 --io --tool none [--tag T]   # then h5diff -c
```

Alternate variants when timing — all `small` numbers above were taken leg-first and variant-first to prove
order-independence.

### Patches

Both are **measurements, not fixes**: apply to a clean tree, build a tagged variant, revert immediately, never
commit.

- `patches/face-kernel-s3.patch` (`B0`) — compile-time `SC = 3` in the **loop bounds** of the face kernel.
- `patches/stencil-smax3.patch` (`smax3`) — `S_MAX` 5 → 3 in `adam_flume_parameters.F90:144` and
  `adam_weno_object.F90:72-73`, shrinking every private stencil **declaration**.

Together they are the two halves of §6; separately, each is silently wrong for any scheme above WENO-5.

---

## 9. Open questions

1. **What holds the 428 registers?** Not the stencil arrays — `smax3` removed 38 % of them and moved nothing.
   Until this is answered, no occupancy gain should be promised from any of this work. Candidates: the scalar
   live values across the inlined call graph, and the 151 SGPR spills.
2. **Why is `noagpr` 16 % faster** with a byte-identical face kernel? Suspect allocation/placement effects;
   unconfirmed.
3. **Do levels 1 and 2 compose?** `smax3` (−21.6 %) and `B0` (−15.4 %) attack different things (declarations vs
   loop bounds). A combined build would measure it — and is exactly what §6 proposes to make legitimate.
4. **PRISM at S = 1** could gain more than FLUME. Unverified, and the centred schemes may differ.
5. **Does the gain hold at scale?** Everything here is one rank on one GPU. Multi-rank `large` runs
   (`--np 8`) would show whether the kernel saving survives alongside communication.

---

## 10. Results of A1 and A2 (2026-10-27, Slurm job 73716 on `TheraC18`)

Same case, same commit `8cb8d339`, one rank on one MI300X, same toolchain as §8 (AFAR drop's `amdflang`, spack Open MPI
5.0.10). `large` numbers are `--tool none`, **one run each**; the leg's spread across three runs was 1.3 %. Run directories:
`runs/fnl-amd-*-s3b0/`, `runs/fnl-amd-*-lean/`.

### A1 — do `smax3` and `B0` compose? **No: they interfere.**

Build `s3b0` = `face-kernel-s3.patch` + `stencil-smax3.patch` applied together, reverted right after the build.

| Build | `large` loop | vs leg | `small` face kernel (rocprof) |
|---|---|---|---|
| leg | 750.36 s (mean of 3) | — | 62.3 ms |
| `noagpr` | 633.31 s | −15.6 % | 51.6–51.9 ms |
| `smax3` | 588.27 s | −21.6 % | 48.8 ms |
| **`s3b0`** | **647.72 s** | **−13.7 %** | **53.8 ms** |

- Outcome: the "above 588 s — they interfere" row of the A1 table. The combined build is 10 % **slower** than `smax3`
  alone. `s3b0` is also slower than `noagpr`, whose cause is unexplained (§5).
- Correctness: `h5diff -c` against the leg is bitwise identical at step 0 and step 10 (the weak check of §1).
- Euler `compute_face_fluxes_dev` record: 428 VGPR / 172 AGPR / 0 VGPR spills, occupancy 1 (unchanged); scratch
  **5152 B** (the same as `smax3`, so `B0` removed no scratch); SGPR spills **127** (leg: 151).
- Consequence for §6: level 2 (the `WENO_S` macro in the loop bounds) does **not** add to level 1 once the
  declarations have shrunk, and looks harmful. **A3 should implement level 1 only** (`ADAM_S_MAX` + the init guard) and
  leave `WENO_S` out until the interference is understood. The optional `B0` alone with `--tool none` was not run.

### A2 — the lean Euler face kernel: partial gain, bitwise identical

Build `lean`: the changes below, left in the working tree (uncommitted; `run.info` records `src_modified=yes`).

| Build | `large` loop | vs leg | `small` face kernel (rocprof) | scratch |
|---|---|---|---|---|
| leg | 750.36 s | — | 62.3 ms | 6048 B |
| `smax3` | 588.27 s | −21.6 % | 48.8 ms | 5152 B |
| **`lean`** | **654.35 s** (wall 676.7 s) | **−12.8 %** | **54.4 ms** (−12.7 %) | **5280 B** |

- Correctness: `h5diff -c` against the leg is bitwise identical at step 0 and step 10 (the weak check of §1; the FLUME
  and PRISM FNL regression suites have **not** been run).
- Record: 428 VGPR / 172 AGPR (unchanged), occupancy 1, **VGPR spills 8** (leg: 0), **SGPR spills 223** (leg: 151).
- Scratch fell by **768 B**, not the predicted 1040 B (6048 → ~5008 B): about 270 B of the expected saving did not
  materialise. At the §4 rate (≈ −20 % per −900 B) 768 B should have bought ≈ −17 %; it bought −12.8 %.
- Outcome: the "materially less than `smax3`" row of the A2 table. The scratch-bytes model is incomplete: re-open §9.1
  before spending the ~110 files of #1. It is **not** the "no gain" row, so #1 is not dead.
- Caveat: `noagpr` was 15.6 % faster with a byte-identical face kernel (§5). That unexplained effect is as large as the
  gap between `lean` and `smax3`, so the per-byte rate above should not be trusted until it is understood.

#### What the lean change does to the code

Shared kernel (`adam_flume_fnl_face_kernels_agnostic.INC`), what every model uses. Each thread (one per face) gathers the
whole stencil into private arrays, then hands them to the split routine:
```fortran
real(R8P) :: qs(NV_K,1-S_MAX:S_MAX)        ! 5 x 10 doubles = 400 B  \  1040 B of private state per thread,
real(R8P) :: qas(NV_AUX_K,1-S_MAX:S_MAX)   ! 8 x 10 doubles = 640 B  /  held in scratch

do m=1-S, S                                  ! gather: copy the stencil from global memory
   do w=1, NV_K;     qs(w,m)  = q_gpu(b,i+m*di,j+m*dj,k+m*dk,w);      enddo
   do w=1, NV_AUX_K; qas(w,m) = q_aux_gpu(b,i+m*di,j+m*dj,k+m*dk,w);  enddo
enddo
call face_split_fluxes(gamma=gamma, ch=ch, d=d, S=S, is_characteristic=is_characteristic, qs=qs, qas=qas, &
                       fsplit=fsplit, er=er, mu=mu)
```

Lean kernel (`adam_flume_fnl_face_kernels_euler_lean.INC`). No `qs` / `qas`, no gather loop, `private(qs,qas,...)`
becomes `private(...)`; the thread passes the global arrays and its position instead:
```fortran
call compute_face_split_fluxes_lean(gamma=gamma, d=d, S=S, is_characteristic=is_characteristic, q_gpu=q_gpu, &
                                    q_aux_gpu=q_aux_gpu, ngc=ngc, b=b, i=i, j=j, k=k, di=di, dj=dj, dk=dk, &
                                    fsplit=fsplit, er=er, mu=mu)
```

New routine `compute_face_split_fluxes_lean` (`adam_flume_euler_library_body.INC`, purely additive). Each stencil cell is
read when needed into one-cell buffers `q_m(NV_K)` / `qa_m(NV_AUX_K)`, used and discarded:
```fortran
do m=1-S, S
   do n=1, NV_K;     q_m(n)  = q_gpu(b,i+m*di,j+m*dj,k+m*dk,n);      enddo   ! 5 doubles
   do n=1, NV_AUX_K; qa_m(n) = q_aux_gpu(b,i+m*di,j+m*dj,k+m*dk,n);  enddo   ! 8 doubles
   call compute_flux(d=d, q=q_m, qa=qa_m, f=f)
   call compute_state_magnitudes(q=q_m, qa=qa_m, qm=qm)
   ... split, accumulate mu (unchanged) ...
enddo
```
The routine walks the stencil twice (wave speeds, then split fields), so each cell is read from cache twice instead of being
held for the thread's lifetime. Every arithmetic operation and its order is that of `compute_face_split_fluxes`, hence
the bitwise-identical output. The Roe average of cells 0 and 1 reads their auxiliary values through `qa_m` / `qaR`.

Files touched (all under `src/app/flume/`):

| File | Change |
|---|---|
| `common/adam_flume_euler_library_body.INC` | + `compute_face_split_fluxes_lean` (96 lines, after `compute_face_split_fluxes`) |
| `common/adam_flume_euler_library.F90` | + `public :: compute_face_split_fluxes_lean` (plain Euler library only; the de, de_ps and ps modules that include the same body compile the routine but do not export it) |
| `fnl/adam_flume_fnl_face_kernels_euler_lean.INC` | **new**: the shared face-kernel body minus the gather; `add_gravity_edge_correction_dev` is a verbatim copy and still uses the shared path |
| `fnl/adam_flume_fnl_euler_kernels.F90` | includes the lean `.INC` instead of the shared one; imports the lean split |

Deviation from the A2 text: no `ch` argument on the lean routine (the existing split routine has none and the Euler
adapter drops it). The copy of the face-kernel body is deliberate and temporary (A2); the other twelve models still use
the shared one.

### Environment notes (additions to §8)

1. **`fobis` is not on `PATH`** in a non-interactive shell; it lives in `~/miniforge3/bin`. **Append** that directory
   (`export PATH=$PATH:~/miniforge3/bin`); prepending it would shadow the toolchain.
2. **Do not put `/opt/rocm-7.2.1/bin` first on `PATH` for the build.** It holds its own `amdflang`, which would replace
   the AFAR drop's and break comparability with the leg. Put it first only for the profiling runs (it is there for
   `rocprofv3` 1.1.0).
3. **`profile.sh` needs Python ≥ 3.7** (`make_profile.py` uses `from __future__ import annotations`); the system
   `python3` is 3.6.8. Put a directory holding only a `python3` symlink to `~/miniforge3/bin/python3` first on `PATH`.
4. **An allocation ends at its time limit.** Job 73716 ran 09:29–12:29 CDT; `srun --jobid` after that fails
   with "Invalid job id". Budget: one build ≈ 16.5 min, a `large` run ≈ 11 min.
5. The Euler kernel record can be read from the device ELF as in §8; split on `- .agpr_count:` and select the record
   whose text contains `euler_kernelsPcompute_face_fluxes_dev`.

### Next (not run)

- **`smax3` + lean together**: they touch different files, so they stack (`stencil-smax3.patch` on top of the lean working
  tree, tag `s3lean`). Not built: the reservation ended first. The driver is `~/a2b_s3lean.sh`; it needs ≈ 30 min
  (build, `--steps 10 --io` + `h5diff -c`, `small` rocprof, `large`).
- `B0` alone with `--tool none`, for a clean three-way (A1, optional).
- The register mystery (§9.1) and `noagpr` (§9.2) are now the blocking questions for the scratch-bytes model.

---

## 11. Re-run on the per-width tree and the new toolchain (2026-10-28, Slurm job 74031 on `TheraC18`)

Everything above (§1–§10) is commit `8cb8d339` with the spack AFAR drop (`amdflang` 23.2.0, HDF5 1.14.6). This section is
commit `faa434ab` (= upstream `5b39d420` + the AMD notes) with the system modules `rocm/10.1.0-afar-24.3.0`
(`amdflang` **24.0.0git**), `openmpi/5.0.10-ucx1.22.0-ucc1.9.0`, `openblas/0.3.34`, `hdf5/2.2.0`. **The two sets of numbers are
not comparable**: the commit and the compiler both changed, and the old `leg` was not rebuilt.

### 11.1 What changed in the code

Commit `421cd773` (the "yesterday" commit) is **not** "no S in the loops": the host still reads the WENO order at run time,
but `compute_face_fluxes_dev` and `add_gravity_edge_correction_dev` became dispatchers (`select case (S)`) onto kernels
compiled with `S` a **named constant** (`_s1_dev` … `_s5_dev`), the generic `_sn_dev` (S a dummy argument) being the
`case default`. All kernels include one body (`adam_flume_fnl_face_fluxes_kernel_body.INC`). The private arrays are still
declared with `S_MAX = 5` (`qs`, `qas`, `fsplit`, `v`), so only **level 2** (loop bounds) of §6 is in the code; level 1
(declarations) is still `stencil-smax3.patch`. Consequences for the tests:

- experiment `B0` is now the code: `patches/face-kernel-s3.patch` no longer applies (kept as the record of the experiment);
- the old `leg` ≈ the new `generic` variant (`patches/generic-kernel.patch`: `select case (-S)`, both dispatchers);
- `stencil-smax3.patch` still applies cleanly; the old "s3b0" build is now `head` + `smax3`;
- the untracked lean kernel (`adam_flume_fnl_face_kernels_euler_lean.INC`) copies the shared body that was since moved into the
  body `.INC` files: **stale, not rebuilt**.

The commit message says the AMD build was never compiled. It compiles and links on this toolchain (below), and the three
builds are bitwise identical at step 0 and 10.

### 11.2 Environment differences (what had to change)

| Item | Before | Now |
|---|---|---|
| toolchain | `source spack-setup-env.sh; module load rocm-afar-drop/…` | `module load rocm/10.1.0-afar-24.3.0 openmpi/5.0.10-ucx1.22.0-ucc1.9.0 openblas/0.3.34 hdf5/2.2.0` |
| HDF5 | 1.14.6, `$HDF5_PREFIX/lib` | 2.2.0 in **`lib64`**; `fobos` reads `$HDF5_PREFIX/lib` and `/include` → varset **`thera_amd2`** (local, uncommitted) with `$HDF5_PREFIX=~/hdf5-2.2.0-shim`, a directory of two symlinks (`include`, `lib`→`lib64`) |
| device link | worked | `clang-linker-wrapper: error: unable to find library -lamdhip64`, after a full compile → **`export LIBRARY_PATH=$ROCM_PATH/lib`**; `build.sh --no-clean` then relinks (454 s) from the surviving objects |
| `rocprofv3` | 1.2.3 aborts, needs `/opt/rocm-7.2.1` (1.1.0) | the module's **1.4.1 traces the executable to the end**, no workaround |
| parallel builds | one tree | one `git worktree` per variant; **`src/third_party` must be a copy** (FoBiS does not follow symlinks: `penf.mod` not found otherwise) |

`build.sh` records `src/ modified` only for `src/`, so the local `fobos` varset does not mark a build as modified. A clean build
is ~9 min to compile + ~7.5 min device-LTO link (≈ 16.5 min, as before), three builds at 32 jobs each ran at once.

### 11.3 Results

All runs: one rank, one MI300X, `profile.sh fnl-amd --np 1`, run dirs in `runs/`, rows in `results/runs-amd.csv` (the previous
file is `results/runs-amd-before-2026-10-28.csv`).

**Small, `rocprofv3`, face kernel per dispatch (300 dispatches), alternated `head`/`generic`/`smax3`, three rounds, one run at a time on GPU 0:**

| Build | Kernel | ms per dispatch (3 runs) | mean | loop (3 runs) | vs `head` | VGPR / AGPR | scratch |
|---|---|---|---|---|---|---|---|
| `head` | `compute_face_fluxes_s3_dev` | 53.82, 52.65, 52.72 | 53.07 | 22.13, 21.48, 21.67 s | — | 128 / 224 | 9088 B |
| `generic` | `compute_face_fluxes_sn_dev` | 56.53, 58.59, 57.58 | 57.57 | 22.91, 23.54, 23.22 s | **+8.5 %** | 128 / 224 | 9216 B |
| `smax3` | `compute_face_fluxes_s3_dev` | 50.06, 50.07, 50.08 | 50.07 | 20.80, 20.80, 20.83 s | **−5.6 %** | 128 / 224 | 8192 B |

The ranges do not overlap. The face kernel is 97 % of GPU time in `head` (16.3 s of 16.8 s).

**Large, refined walls, `--tool none`, forest loop (20 steps):**

| Build | alone on GPU 0 | at the same time as the other two (GPU 0 / 1 / 2) |
|---|---|---|
| `head` | **631.84 s** (wall 652.5) | 636.25 s (GPU 0) |
| `generic` | not run (no time) | 699.15 s (GPU 1) |
| `smax3` | **611.19 s** (wall 631.9), **−3.3 %** | 690.38 s (GPU 2) |

- Running `head` alone changes its loop by 0.7 % (GPU 0 both times), so the concurrency itself is not the problem.
- **The GPU (or the contention among the three concurrent runs) matters:** `smax3` took 690.4 s on GPU 2 and 611.2 s alone on
  GPU 0, +13 %; which of the two causes it was is not separated here. The concurrent large column therefore
  compares GPUs as much as builds; its `generic` and `smax3` rows (+9.9 %, +8.5 % over `head`) must not be read as build
  effects. Always pin the GPU (`ROCR_VISIBLE_DEVICES`) when comparing variants, and use the same one.
- Correctness: `h5diff -c` of `head` against `generic` and `smax3`, small refined walls, `--steps 10 --io`: bitwise identical at
  step 0 and 10 (positive control, step 0 against step 10 of `head`: differs). Same weak check as §1 (the step-10 regrid moves no
  data); the FLUME and PRISM FNL regression suites were **not** run.

### 11.4 What this says about §0–§10

1. **The per-width kernel is worth ≈ 8 % on MI300X** (face kernel, small, alternated), against −46 % at WENO-5 on the A30 (commit
   message). That is smaller than the −15.4 % the hand-edited `B0` gave on `large` in §1 (under rocprof); the compiler and the
   baseline changed, so it is a different measurement, and the large-size number for it (`generic` alone) is missing (§11.5).
2. **The interference of §10 A1 is not reproduced.** There, `smax3` + `B0` (`s3b0`) was 10 % *slower* than `smax3` alone. Here
   `head` is B0 in the code and `smax3` on top of it is **faster** (−5.6 % small, −3.3 % large, both alone on GPU 0). With the
   old toolchain and the hand-edited B0 they interfered; with the new one they compose, weakly.
3. **The scratch-bytes rate shrank.** `smax3` removes 896 B (9088 → 8192) and buys −5.6 % / −3.3 %, against ≈ −20 % per −900 B in §4
   (compiler 23.2.0). The registers moved too: 128 VGPR / 224 AGPR now (428 / 172 before), scratch is 9 KB (6 KB before). §3's
   private-state budget and §4's metric belong to the old compiler and have not been re-derived.
4. §6's recommendation (level 1 as a build-time `S_MAX`) is worth less than it looked (−3 to −6 %, not −21 %) once level 2 is in the code.
   It is still bitwise identical and still unsafe as a patch (silent beyond WENO-5), so it needs the init guard if adopted.

### 11.5 Not done

- `generic` large alone on GPU 0 (the only way to turn the 8 % small result into a loop number); `noagpr`, `lean` and
  `s3lean` were not rebuilt on the new tree.
- `--np 8`, uniform sizes, other WENO orders (`--scheme weno-u-<2S-1>`: only `weno-u-5` was run; the AMD builds of S = 1, 2, 4, 5 are compiled in `head` but untested).
- The FLUME/PRISM FNL regression suites on AMD. The FNL gate for the per-width commit was run on the A30 only.
- `make_figures.py figures` draws the AMD figures from large `rocprof` runs of the variants `leg`/`B0`; none of this section's runs is of
  that form (large runs were `--tool none`, small ones `rocprof`), so no new MI300X figure was drawn.

### 11.6 What is left on the machine

- Local, uncommitted: `fobos` (varset `thera_amd2`), `~/hdf5-2.2.0-shim/`, `profile-couette3d/README.md` (the patches, the
  system-module toolchain and the rocprofv3 note), `patches/generic-kernel.patch`, `results/runs-amd.csv` (+ the backup).
- Worktrees `~/Projects/adam-generic` (patch `generic-kernel.patch` applied) and `~/Projects/adam-smax3` (patch `stencil-smax3.patch`
  applied) hold the `generic` and `smax3` builds and their original run dirs; `git worktree remove --force` deletes them. The three
  executables are `exe/fnl-amd-head` (main tree), and `exe/fnl-amd-{generic,smax3}` in the worktrees.
- Drivers in `~`: `perf_build.sh`, `perf_build2.sh` (adds `LIBRARY_PATH`), `perf_run.sh`, `perf_suite.sh`, `perf_alt.sh`,
  `perf_large_seq2.sh`; logs in `~/perf_logs/`.
- The first sequential large attempt (07:48) was killed by a `pkill` that matched its own ssh shell and was restarted; its partial
  run dir is `runs/fnl-amd-large-np1-none-head.old-20261028-075152` (ignored by `collect`).

---

## 12. The cost was in the WENO routine: a per-width reconstruction takes the face kernel from 52.6 to 3.0 ms (2026-10-28, Slurm job 74047 on `TheraC18`)

Same tree (`faa434ab`), same toolchain, same one-rank-one-GPU setup as §11; every number below was measured **alone on GPU 0**, one run
at a time (the concurrent large runs of §11.3 are superseded by §12.3). Run dirs in `runs/`, rows in `results/runs-amd.csv` (70 runs).

### 12.1 Why the per-width kernels did so little on AMD (the device code)

The device ELF of each executable (`llvm-objcopy --dump-section=.llvm.offloading=…`, then `llvm-objdump -d`, as in §8) shows, for the
Euler `compute_face_fluxes_s3_dev` kernel of the tree `head`: 3715 instructions, 66 scratch loads and 112 stores, and **6 call sites**:
two are `__llvm_libc::gpu::allocate` (a **per-thread heap allocation of 40 B** at kernel entry, and one sized 8 B times the threads
of the workgroup behind a null test, so once per workgroup: read from the instructions around the calls, not from the source), two are `compute_face_split_fluxes` (not inlined, called with pointers to the thread's scratch),
and two are unresolved (register set elsewhere; by their arguments, the WENO reconstruction, which is also a separate function in
the image). The generic `sn` kernel is the same (3608 instructions, 6 calls). With S a named constant only the thin outer kernel knows
it: the callees still take S as a runtime argument and the WENO tables are indexed by it. The generic `weno_reconstruct_upwind_dev`
itself has 1226 instructions, **6 calls** (its three helpers are not inlined into it) and 94 scratch accesses.

### 12.2 The change (`patches/weno-width.patch`, variant `wenow`, 4 files, nothing else touched)

- new `src/lib/fnl/adam_fnl_weno_width_kernels.F90` and `adam_fnl_weno_width_body.INC`: `weno_reconstruct_upwind_s1_dev` … `s5_dev`,
  S a named constant, **no S argument**, the polynomials, the weights and the convolution of `adam_fnl_weno_kernels` merged into one
  routine **in the same order of operations**, the private `VP`, `w`, `IS`, `a` sized by S (not by `S_MAX`);
- `adam_flume_fnl_face_fluxes_kernel_body.INC`: `#ifdef FNL_WENO_WIDTH` call the specialised routine, else the generic one;
- `adam_flume_fnl_face_kernels_agnostic.INC`: each per-width face kernel `use`s its routine and defines `FNL_WENO_WIDTH` around the
  body (the generic `_sn_dev` and the gravity-edge kernels are unchanged).

Compiled device code of `compute_face_fluxes_s3_dev`: the specialised reconstruction has **0 call sites, 20 scratch accesses**, 1555
instructions of which 1122 FP64 (fully unrolled); the kernel now unrolls the 5-field loop, so it calls it 5 times per branch (14 call
sites in the kernel against 6, 3937 instructions against 3715); `compute_face_split_fluxes` and the two allocations are still calls.
Build: 946 s (16 jobs, a clean build, link included), no change to the environment of §11.2.

### 12.3 Results

**Small, `rocprofv3`, face kernel per dispatch (300 dispatches), alone on GPU 0, alternated `head`/`wenow`, 3 rounds** (`generic` of §11.3 for reference):

| WENO (S) | `generic` | `head` (per-width kernel) | **`wenow`** | `wenow` vs `head` | loop (generic / head / wenow) | `wenow` VGPR / AGPR / scratch |
|---|---|---|---|---|---|---|
| `weno-u-1` (1) | 52.71 ms | 49.25 | **2.53** | −94.9 % | 21.6 / 20.7 / **6.0** s | 84 / 132 / 2288 B |
| `weno-u-3` (2) | 54.98 | 50.85 | **3.61** | −92.9 % | 22.3 / 21.0 / **6.4** s | 84 / 132 / 2784 B |
| `weno-u-5` (3) | 57.57 | 52.62 (52.61, 52.64, 52.61) | **3.03** (2.99, 3.04, 3.05) | **−94.2 %** | 23.2 / 21.5 / **6.0** s | 84 / 132 / 2784 B |
| `weno-u-7` (4) | 58.95 | 53.84 | **3.50** | −93.5 % | 24.7 / 23.0 / **7.2** s | 100 / 132 / 2784 B |
| `weno-u-9` (5) | 59.15 | 54.27 | **4.05** | −92.5 % | 25.7 / 24.1 / **8.3** s | 120 / 128 / 2784 B |

(`head` and `generic`: 128 / 224 / 9088 B and 9216 B at every width. One run per width except `weno-u-5`, three; the `weno-u-5`
rounds agree to 2 %.)

**Large (refined walls, 20 steps), `--tool none`, forest loop, alone on GPU 0:**

| Build | runs (s) | mean | vs `head` |
|---|---|---|---|
| `head` | 631.8, 632.1, 632.7 | 632.2 | — |
| `generic` | 693.8, 705.2 | 699.5 | +10.6 % |
| `smax3` | 611.2, 609.0, 607.9 | 609.4 | −3.6 % |
| **`wenow`** | **50.56, 50.79** | **50.7** | **−92.0 % (12.5×)** |

The `head` runs are 75 minutes apart and spread by 0.1 %. These settle the §11.3 open items: the per-width kernel is **−9.6 %** at the
loop level (not +9.9 % as in the concurrent column), and `smax3` on top of it **−3.6 %** (not +8.5 %); the concurrent large column of
§11.3 is not to be used (its builds ran on different GPUs at the same time: the GPU or the contention, not the build, moved them).

**Where the large `wenow` loop goes** (one rocprofv3 run, loop 56.6 s with the profiler against 50.6 s without; kernels 44.1 s):
`compute_face_fluxes_s3_dev` 31.6 s (105.2 ms per dispatch, 71.5 % of kernel time; **2410.5 ms in §2**), `add_dissipative_fluxes_o4_dev`
8.9 s (29.6 ms, 20.1 %), `update_ghost_local_gpu_dev` 1.9 s (4.2 %), the rest under 1 % each; about 12 s of the loop have no kernel
running. By §2's normalisation (46 M faces) the face kernel is **2.3 ns per face** (52.4 before), against 10.3 / 6.1 ns for the A30 leg /
`lto` of §2 and its best build `B0-lto` (≈ 3.3); different cases and builds, indicative only.

### 12.4 Correctness

`h5diff -c` of `head` against `wenow`, `small`, 10 steps, checkpoints at steps 0 and 10 (`--steps 10 --io`): **bitwise identical for
all five widths (`weno-u-1`, 3, 5, 7, 9) and for the uniform variant: 12 comparisons, all exit 0** (the control of §11.3, step 0
against step 10, differs). It is the same weak check as §1 and §11.3: the step-10 regrid of this case moves no data, and only the Euler
model's WENO-5 is the physics of the case. The A30 gate of commit 421cd773 saw 3 of 172 cases not bitwise identical with S constant; that
is not seen here, but the FLUME and PRISM FNL regression suites have **not** been run on AMD with this patch.

### 12.5 What this changes in the picture

1. The face kernel was **not** slow because of the private stencil arrays of the kernel (§3, §4, §6, §7): `smax3` removed 896 B and gave
   −3.6 %, this patch removes 6304 B of scratch (9088 → 2784) and gives −94 %, so scratch bytes were not what ranked the changes. What
   the device code shows is the cost of the **calls**: a non-inlined WENO routine whose helpers are not inlined either, with tables
   indexed by a runtime S and loops of runtime bounds, fed through pointers into scratch. Constant S **inside the callee** removes the
   calls, unrolls the loops and drops the scratch traffic; constant S in the caller alone (421cd773) did almost nothing on AMD. This is
   consistent with the measurements, not proven by a counter run (no rocprof-compute pass was made on `wenow`).
2. The route of §6 (a build-time WENO order, `ADAM_S_MAX`) is not needed for this gain, and §7's #1, #3, #4 (`lean`, no `qas`, no
   `fsplit`) are now second-order for the face kernel; they are untested on the new tree.
3. Register occupancy is still 1 wave per SIMD at WENO-5 (84 + 132 = 216 unified registers is below 256, so 2 waves are possible in
   principle; the occupancy was **not** read from the kernel record here, the registers were). Not examined.
4. The same trap probably exists elsewhere: `compute_face_split_fluxes` is still an unresolved call per face (and the two heap
   allocations at kernel entry), the dissipative kernel (20 % of the kernel time now), and every other caller of the generic
   `weno_reconstruct_upwind_dev` (NASTO, RIPPLE, the gravity-edge kernels, the interpolation unit test) keeps the old path.

### 12.6 Not done

- The FLUME and PRISM FNL regression suites on AMD; a check with a case that moves data at a regrid and with models other than Euler (all
  models' per-width kernels compiled and linked: the executable holds them; only Euler was run).
- Multi-rank (`--np 8`), uniform sizes and `rocprof-compute` counters on `wenow`; the gravity-edge kernels (same pattern, not converted).
- `wenow` on the other compilers: the new module is compiled by every FNL mode (NVF OpenACC too) but was built with `amdflang` only.

### 12.7 On the machine

- `~/Projects/adam-wenow` (worktree, patch applied, `exe/fnl-amd-wenow`); the patch is `patches/weno-width.patch` in the case directory
  (applies to the main tree, checked with `git apply --check`); `git worktree remove --force` deletes the worktrees
  (`adam-generic`, `adam-smax3`, `adam-wenow`).
- Drivers in `~` (`perf_*.sh`), logs in `~/perf_logs/`, device images and disassembly in `~/isa/`.
- Nothing committed; `fobos` still has the local varset `thera_amd2`.

---

## 13. The general version, branch `WENO_refactor` (2026-10-28, Slurm job 74047 on `TheraC18`)

Branch `WENO_refactor` of `origin` (github.com:szaghi/adam), one commit, `b9775885`, cut from `origin/develop` (= `faa434ab`); the patch of §12
(`wenow`) was flux-splitting face kernels only. The general version, in 9 files (715 insertions, 172 deletions):

- `src/lib/fnl/adam_fnl_weno_width_kernels.F90` + `adam_fnl_weno_width_body.INC`: `weno_reconstruct_upwind_s<N>_dev` **and**
  `weno_reconstruct_upwind_wratio_s<N>_dev` (the flavour the Riemann kernels need), N = 1…5, S a named constant; one body, the
  macro `FNL_WENO_WRATIO` only selects the flavour; the three helpers of the generic routine merged in the same order of operations;
- flux splitting: `..._face_fluxes_kernel_body.INC` and `..._gravity_edge_kernel_body.INC` call the routine of their width through
  the macro `FNL_WENO_WIDTH` (one `#ifdef` each); `adam_flume_fnl_face_kernels_agnostic.INC` uses and defines it in every per-width
  kernel, face **and** gravity edge;
- Riemann (`weno-riemann`, every model and solver, ~40 modules share the file): `adam_flume_fnl_riemann_face_kernels_agnostic.INC`
  restructured as the flux-splitting ones were in 421cd773: a host dispatcher, kernels `_s1` … `_s5`, the generic `_sn`; the bodies
  moved unchanged but for S leaving the firstprivate lists into two new `..._riemann_face_fluxes_kernel_body.INC`,
  `..._riemann_gravity_edge_kernel_body.INC` (generated mechanically from the old file by a script, then read);
- `test_flume_weno_interpolation_fnl.F90`: per-width routines (plain and with the ratio) against the generic ones, bit for bit, on the
  device, every scheme, both weights, 4096 stencils.

Left alone: the generic `weno_reconstruct_upwind_dev` / `..._wratio_dev` (the `_sn` kernels, NASTO's reduced-order fallback, RIPPLE and
the other run-time-S callers), the CPU backend, PRISM. Macros: two, each names a procedure (no math, no loop bound), documented at
the use, per CLAUDE.md ("document WHY").

### 13.1 What was run (AMD, `thera_amd2`, §11.2 environment; nothing on NVHPC)

| Check | Result |
|---|---|
| clean build of `flume-fnl-omp-amd`, every model and solver | **1474 s, exit 0** (with the unit-test build running at the same time; §12: 946 s for the smaller change at 16 jobs) |
| `test_flume_weno_interpolation_fnl` (device vs host, per-width vs generic) | **PASSED** (build `test-flume-weno-interpolation-fnl-omp-amd`): per-width routines bit for bit the generic ones, 5 schemes × 2 weights × 4096 stencils, plain and ratio |
| `profile-couette3d` small, `--steps 10 --io`, `h5diff -c` against `head` (the tree without the change) | steps 0 and 10 **identical** |
| `profile-couette3d` small, rocprofv3, WENO-5, face kernel | **3.31 ms** per dispatch (the first version, alone on the node: 3.03; this one ran at the same time as the regression runs on other GPUs), 84 / 132 VGPR/AGPR, 2784 B scratch (as the first version), loop 6.2 s |
| `profile-couette3d` large, `--tool none` | loop **51.4 s** (first version alone: 50.7; `head`: 632.2), same caveat |
| 12 regression inputs (`src/tests/flume/regression`), `mpirun -np 2`, old (`head`) against new, every checkpoint with `h5diff -c`, the residual and conservation logs with `cmp` | 11 of 12 **bitwise identical**: `couette`, `brio-wu`, `orszag-tang` (MHD + GLM), `noh-de` (dual energy), `sod-hllc`, `rj2a-hlld`, `baro-collide-hllc`, `iso-collide-hllc`, `iso-tube`, `baro-tube`, `blast-limiter`; the twelfth, `baro-mhd-glm-hlld`, differs at its last checkpoint (step 641; a few cells of `psi`, `r`, `ru`, `bz` per block) **but the old executable differs from itself between two runs at that checkpoint** (step 0 identical both ways) and the residual and conservation-history logs are identical old vs new: non-determinism of the baseline in that case, not of the change |

The cases of the regression suite that use the Riemann solvers (`sod-hllc`, `rj2a-hlld`, `baro-collide-hllc`, `iso-collide-hllc`,
`baro-mhd-glm-hlld`) are the ones that exercise the restructured `weno-riemann` kernels (WENO-5, so the `_s3` kernels of the dispatcher).

### 13.2 Not checked

- **NVHPC / OpenACC**: the new module, the macros and the `use` inside the kernels are compiled by every FNL mode; only `amdflang` built them.
- The committed goldens of the suite (`golden/fnl`, from nvfortran) were **not** compared: the comparison is the previous tree against the new one, on AMD, bitwise;
  12 of the 57 cases only (self-gravity cases need the FFT executable and were skipped; the multi-realm ones were not run); 2 ranks, one block set each.
- WENO orders other than 5 in the Riemann kernels (the dispatcher's `_s1`, `_s2`, `_s4`, `_s5`, `_sn` kernels are compiled, not run); the gravity-edge kernels
  (flux splitting and Riemann) are compiled, not run.
- The speed of the Riemann kernels: no profile of a `weno-riemann` case was taken, so the gain there is unmeasured (it is the same mechanism: the wratio routine has the same
  helper chain), and the commit message says only what was measured (the flux-splitting face kernel).
- The commit message was written before the builds ran and says the generalised version was untested at that moment; this section is what they found.
  Its numbers are not in the commit; amending it would be a force-push of a branch that exists only on `origin`.

### 13.3 On the machine

`~/Projects/adam-weno-refactor` (the worktree of the branch; `fobos` there still has the local varset `thera_amd2`, **uncommitted** and not in the commit);
the executable `exe/fnl-amd-wrefactor/`; the regression outputs in `~/perf_regr2/{old,new,old2}/` (`old2` is the determinism rerun, the first `~/perf_regr/` is the failed
attempt: Open MPI had one slot per Slurm task, `--oversubscribe` fixed it); `results/runs-amd.csv` has 73 runs.

---

## 14. What is left after the per-width WENO reconstruction: the other four ideas and three observations (2026-10-28)

Where §12/§13 leave the face kernel (WENO-5, `profile-couette3d`, one MI300X, amdflang 24.0.0git): 3.0 ms per dispatch on `small` (52.6 before),
105 ms on `large` (2410 ms in §2), 84 VGPR / 132 AGPR / 2784 B of scratch, and in the `large` rocprofv3 run (loop 56.6 s, kernels 44.1 s)
`compute_face_fluxes_s3_dev` 71.5 % of the kernel time, `add_dissipative_fluxes_o4_dev` 20.1 %, `update_ghost_local_gpu_dev` 4.2 %, about
12 s (21 % of that loop) with no kernel running. Options 1 and 2 are the first two of the five ideas that came out of "how do we get the A30 gain on AMD";
2 is `WENO_refactor` (§13); this section writes down 1 (what is done and what is open), 3, 4, 5, and three observations from the device code
that are not among the five. Everything marked **inference** has not been measured.

### 14.1 Idea 1 — find out why the callees are not inlined (diagnosis done, the cause still open)

**What was done** (CPU only, on the executables; recipe of §8): the device ELF of `head`, `generic`, `wenow` was extracted
(`llvm-objcopy --dump-section=.llvm.offloading=… <exe> /dev/null`, then `llvm-objdump -d`; the images and disassembly are in `~/isa/`, the counting scripts on thera are `/tmp/isa2.py`, `/tmp/calls.py`, `/tmp/fn.py`, `~/perf_diss.py`, `~/perf_alloc.py`),
and for the Euler `compute_face_fluxes_s3_dev` kernel the instructions, the scratch accesses and the call sites (`s_swappc_b64`, the target resolved from the
`s_getpc_b64` + `s_add_u32` pair when it is a plain pc-relative call) were counted:

| Kernel | instructions | scratch loads / stores | call sites |
|---|---|---|---|
| `head` `_s3` (constant S in the kernel only) | 3715 | 66 / 112 | 6: 2 × `__llvm_libc::gpu::allocate`, 2 × `compute_face_split_fluxes`, 2 unresolved (by their arguments the WENO reconstruction) |
| `generic` `_sn` | 3608 | 60 / 103 | 6, same |
| `wenow` `_s3` | 3937 | 108 / 98 | 14: 2 × allocate, 2 × `compute_face_split_fluxes`, **10 unresolved** (5 per branch; **inference**: the 5-field loop is now unrolled, one WENO call per field) |
| generic `weno_reconstruct_upwind_dev` (a separate function in the image) | 1226 | 94 accesses | 6 (its three helpers are not inlined into it either) |
| `weno_reconstruct_upwind_s3_dev` (`wenow`) | 1555 (1122 FP64, fully unrolled) | 20 accesses | **0** |

The conclusion that drove idea 2 (and is **inference**, consistent with the timings): with S a constant only in the kernel, the callees see it as a run-time argument, take the stencil by
pointer into the thread's scratch, and index the tables by it; a constant **inside** the callee removes the calls and the scratch traffic.

**What is open, in the order to try:**
1. Resolve the 10 call sites of `wenow` (my resolver handles `s_getpc_b64` + `s_add_u32` into a plain register; those are set by a different sequence, e.g. through a scalar load of a function pointer). Expected: the WENO s3 routine. 5 minutes.
2. Why does the compiler not inline `compute_face_split_fluxes` and the WENO routine into the kernel? The remarks (`-Rpass-missed=inline -Rpass-analysis=inline`) appear at **link** time here (device LTO runs in the linker wrapper, §8): a relink from the objects, 7.5 min with `LIBRARY_PATH=$ROCM_PATH/lib`, with the flag in the link line (`build.sh --flags` needs a clean build, so call `fobis build --mode flume-fnl-omp-amd … --lflags "…"` by hand after a normal build). Expected: a "too costly to inline (cost=…, threshold=…)" line per call.
3. Try forcing it, in this order of intrusiveness: the Fortran directive `!dir$ forceinline` on the callee (a toy file compiled with this amdflang accepted it silently; the `-O2` IR showed no `alwaysinline` attribute, which proves nothing because the tiny routine may simply have been inlined — **not established that it is honoured**); an `ompx_attribute(__attribute__((always_inline)))` is not supported by flang as far as I know; textual inlining through an include (as the per-width kernels do for the body), which always works and costs code size.
**What it could buy:** the 2784 B of scratch are the arrays the callee receives by pointer; if the split routine and the WENO routine are inlined the compiler can keep them in registers, at the price of a larger kernel and more registers (the `wenow` kernel already has 3937 instructions, 14 call sites).
**Decision criterion:** a relink with the remarks says which threshold is missed and by how much; if the cost is within ×2 of the threshold a directive may be enough, if it is far above only textual inlining will do and the registers must be watched (84 + 132 now).

### 14.2 Idea 3 — stop holding `qs` and `qas` (the "lean" kernel of §10 A2)

**The idea.** The kernel gathers the whole stencil into two private arrays, `qs(NV_K, 1-S_MAX:S_MAX)` and `qas(NV_AUX_K, 1-S_MAX:S_MAX)` (Euler: 5 × 10 and 8 × 10 doubles, **1040 B**), and hands them to `compute_face_split_fluxes`. The lean kernel passes the global arrays and the thread's position instead, and the split routine reads each stencil cell
into one-cell buffers `q_m(NV_K)`, `qa_m(NV_AUX_K)` when it needs it: it already walks the stencil twice (wave speeds, then the split fields), so each cell is read twice from cache instead of held for the thread's lifetime. Every arithmetic operation and its order stay, so the result is bitwise identical by construction.

**What was measured** (§10, old tree and old compiler): scratch 6048 → 5280 B (−768 B, not the −1040 B predicted), face kernel 62.3 → 54.4 ms (−12.7 %), `large` loop −12.8 %, `h5diff -c` identical, but VGPR spills 0 → 8 and SGPR spills 151 → 223.

**State now.** The untracked `src/app/flume/fnl/adam_flume_fnl_face_kernels_euler_lean.INC` is still in the main tree; it is a copy of the face-kernel body **before** `421cd773` moved that body into `adam_flume_fnl_face_fluxes_kernel_body.INC`, so it is stale. The edits to tracked files that A2 listed (`compute_face_split_fluxes_lean` in `common/adam_flume_euler_library_body.INC`, its `public` line, the include in `adam_flume_fnl_euler_kernels.F90`) are not in the working tree as far as `git status` shows (only `fobos` and the README are modified): redo from §10's description.

**What it could buy now (inference from the declarations, not a measurement).** At `S_MAX = 5` the declared private arrays of the kernel are `qs` 400 B + `qas` 640 B + `fsplit(2,1-S_MAX:S_MAX-1,NV_K)` 720 B = 1760 B of the 2784 B of scratch of `wenow`, so the lean kernel touches 37 % of it (`qs` + `qas`) and removing `fsplit` as well (§7 #4, fuse the split with the reconstruction field by field) another 26 %. But the lesson of §12 is that scratch **bytes** were not what ranked the changes (−896 B gave −3.6 %, −6304 B gave −94 %), and the face kernel is now a different kernel (17× faster); there is no basis to carry over the old −12.8 %.
**How to test it cheaply:** Euler only, as A2 did: a lean copy of the new body INC (a macro-free copy is fine for a measurement), one build (~16 min), `small` rocprof alternated three times against `wenow`, then the 12-case regression (`~/perf_regr2.sh`). **Decision criterion:** worth rolling out to the ~13 model adapters (§7: ~110 files) only above ≈ 10 % on the `large` loop; below that, not.

### 14.3 Idea 4 — occupancy (waves per SIMD)

**The idea.** In §1/§4 the face kernel ran at 1 wave per SIMD (428 unified registers, `floor(512/428)`), latency-bound, nothing to hide a stall behind. If the kernel could fit twice, the SIMD would overlap two waves.

**What is known now.** rocprofv3 reports 84 VGPR and 132 AGPR for the `wenow` WENO-5 kernel (100 / 132 at WENO-7, 120 / 128 at WENO-9): the two files add up to about 216 (allocation granules round it up a little), below the 256 that gives 2 waves, and `head` was 128 + 224 = 352 (1 wave). The **occupancy itself was not read** from the kernel record in this work: `floor(512 / unified count)` predicts 2 for `wenow`, so the first thing to do is to read it (the notes of the device ELF, `llvm-readelf --notes dev-wenow.elf`, split on `.agpr_count`; or `-Rpass-analysis=kernel-resource-usage` at link). **It may already be 2 and this idea may already be spent.**

**If it is still 1** (or to get 3): the lever is the `amdgpu-waves-per-eu` function attribute, which is **not** reachable by a flag (§5: every `-mllvm` device channel was rejected). The place is FUNDAL's `OMPLOOP` macro (`src/third_party/FUNDAL/src/lib/fundal.H`, emits a bare `target teams distribute parallel do`): add an `ompx_attribute(...)` clause or a `thread_limit`, which means a change in the FUNDAL repository (a third-party dependency of this one) and a check that `flang` accepts the clause (unverified here). A runtime route exists without a rebuild, `OMP_NUM_TEAMS` / `OMP_TEAMS_THREAD_LIMIT` (pass them through the launcher, README of the case), but it changes the launch shape, not the register count.
**Risk:** a cap makes the compiler spill: on the A30 the equivalent (`-gpu=maxregcount:128`) made the face kernel 16 % **slower** (§5). **Decision criterion:** only if the occupancy read says 1; then try the cap at the register count that gives 2 and keep it if the kernel gets faster alternated against the uncapped one.

### 14.4 Idea 5 — compute the cell-centred quantities once per cell, not once per face

**The idea.** For each face the kernel evaluates, at each cell of the stencil (`2S` = 6 cells at WENO-5), the physical flux `f_d(q)` (`compute_flux`), the state magnitudes (`compute_state_magnitudes`) and, in the characteristic variant, the wave speeds `|lambda_k(m)|` for the maximum `alpha_k = max_m |lambda_k(m)|`; the Roe eigenvectors are for the two cells around the face only. Those quantities depend on the **cell** and the direction, not on the face, yet every cell is in the stencil of `2S` faces along each direction, so they are computed 6 times (at WENO-5; 10 at WENO-9). Precompute them in one cheap cell-centred kernel (the way `compute_q_aux_dev` already produces `q_aux` per cell, ghosts included, §7 #3) into global arrays and let the face kernel gather them: fewer divisions and square roots per face, fewer values to derive in the thread, a shorter dependent chain.

**Cost.** A new cell kernel per direction and per model, extra global arrays (`f`: `NV_K` doubles per cell; for `large`, 46 M cells, about 1.8 GB per direction if stored for all cells, reusable between directions), and changes in the per-model library bodies: the same ~13 `compute_face_split_fluxes` bodies (Euler, MHD, GLM, dual energy, barotropic, isothermal, passive scalars) and their ~100 adapters, i.e. the biggest change of the five. It can stay bitwise identical: the same operations produce the same bits whether they are stored and reloaded or recomputed. **Risk:** the face kernel would read more from global memory (HBM, not cache) than now; whether the kernel is bound by arithmetic or by the latency of its chain after §12 is **not known**: the old counters (§2: 0.216 % of FP64 peak, 1 wave) no longer describe this kernel, and **no counter run was made on `wenow`** (`rocprof-compute` would say: VALU utilisation, scratch, L2 hit rate).
**Decision criterion:** take the counters first. If the new kernel is arithmetic-bound or dependency-bound on the cell quantities, prototype for Euler only (flux and magnitudes), measure as above; if it is bound by memory, do not.

### 14.5 Three observations that are not among the five

**(a) Two device heap allocations at the entry of (almost) every kernel.** The `_s3` kernel starts with two calls to `__llvm_libc::gpu::allocate` (the device allocator of the LLVM libc: a global slab allocator): one of **40 B executed per thread** (under the thread's exec mask) and one sized `8 B × the threads of the workgroup` (`s_min(workgroup size, grid − block × size) << 3`) **behind a null test on an LDS slot**, so by the instructions once per workgroup (**inference**: read from the code around the calls, not from the source). They are **not** specific to the face kernel: the dissipative kernel has exactly the same two; of the **512** offload kernels in the `wenow` executable **none has 0 call sites, 255 have exactly 2, 257 more than 2**, and the 2-call ones include the small kernels (`compute_q_aux_dev`, `rk_compute_stage_dev`, `blend_inadmissible_ghosts_dev`: no private array to allocate). I did not resolve the targets of all 255, so "the same two allocate calls" is the natural reading, not a proof.
Why it matters: a device `malloc` per thread is an atomic on a global structure; for the kernels that run 1854 times at ~28 µs on `small` and ~99 µs on `large` (`blend_inadmissible_ghosts_dev`, 3.7 % of the `small` and 0.4 % of the `large` kernel time) it may be a large share of their own time, and it is inside the kernel time (the 12 s per `large` loop without a kernel are host time, a separate matter). **What to do:** find what the 40 B are (5 doubles: `mu(NV_K)` or `flux(NV_K)` are the only 40 B objects of the face kernel, but the tiny kernels have no such array, so more likely a piece of the OpenMP device-runtime state of the `collapse(n)` loop nest); look at the LLVM IR of one tiny kernel (`-save-temps` of one module with the kernel flags, or `-emit-llvm` for the device) for the call and its origin; a microbenchmark of an empty `target teams distribute parallel do` kernel under rocprofv3 gives the cost directly. This is a **compiler / runtime** question first (a report to the AMD toolchain team with the reproducer is the likely outcome), not a source change of this project.

**(b) `compute_face_split_fluxes` is still a real call, twice in the kernel** (the Euler one, `_QMadam_flume_euler_libraryPcompute_face_split_fluxes`, both in `head` and in `wenow`): it receives pointers to `qs`, `qas`, `fsplit`, `er`, `mu` in scratch (the constants loaded before the call are scratch offsets) and takes S as a run-time argument, with `S_MAX`-sized dummies. The same per-width treatment as the WENO routine applies in principle (a body included five times with S a parameter and the arrays sized by S), per model: 13 library bodies × 5 widths of device code (§6: "4× the device code in 13 model modules"), a heavier build (the build is already 16–25 min). What its share of the 3 ms is **is not known**; to measure it cheaply, a copy of the Euler routine only with S a parameter for WENO-5 (the profile case's only width), swapped in a worktree. This is the follow-up with the most direct line from §12's result, and it interacts with idea 3 (they touch the same arrays).

**(c) The dissipative kernel is the next one, but it is not the same disease.** `add_dissipative_fluxes_o4_dev` (WENO fourth-order viscous fluxes) is 8.9 s of the 44.1 s of kernel time of the `large` run (29.6 ms per dispatch, 300 dispatches), 20 %, the second after the face kernel. Its device code is **already inlined**: 8732 instructions, 2420 FP64, **no call except the two allocations**, scratch loads/stores 96 / 16, **272 B of scratch** (rocprofv3), 104 VGPR + 128 AGPR (§4 had 21744 B of scratch, 428 / 172 and 70× less time than the face kernel, with the old compiler and the older tree). So the "non-inlined chain with a run-time S" explanation of the face kernel does **not** carry over to it, and I have no hypothesis for its 29.6 ms; the counters (`rocprof-compute`) are the way: if it is memory-bound (it reads a wide stencil of `q` and `q_aux` and writes the flux), a fused or cache-blocked variant is the lever, not specialisation. It also runs 300 times with the face kernel (3 directions × 5 stages × 20 steps), so the pair together is 91.6 % of the kernel time (71.5 + 20.1).
Also seen: **about 12 s of the `large` loop (56.6 s under rocprofv3, kernels 44.1 s) have no kernel running**: host time (launches, the regrids at steps 10 and 20, the ghost updates `update_ghost_local_gpu_dev` is kernel time, 1.9 s). After the face-kernel fix this is 21 % of the loop and the next thing a timeline would show (rocprofv3 `--sys-trace` / the kernel trace gaps); it was never looked at, since before it was 2 %.

### 14.6 Order, and the protocol to measure any of them

| # | What | Cost | Needs a GPU | First action | Go / no-go |
|---|---|---|---|---|---|
| 1 | read the occupancy of `wenow` (idea 4, step 1) | minutes, CPU | no | `llvm-readelf --notes` on `~/isa/dev-wenow.elf` | decides whether idea 4 exists |
| 2 | inline remarks at link (idea 1, steps 1–3) | 7.5 min relink | no | `-Rpass-missed=inline` on the link line | decides forceinline vs textual inlining |
| 3 | counters on `wenow`: face and dissipative kernel (`rocprof-compute`) | 20 min | yes | one `small` run under `rocprof-compute` | decides ideas 5 and 14.5 (c) |
| 4 | per-width `compute_face_split_fluxes`, Euler only, WENO-5 (14.5 b) | one build (16 min) + alternated runs | yes | worktree from `WENO_refactor` | ≥ 10 % on the `large` loop → roll out |
| 5 | lean kernel on the new body (idea 3) | one build | yes | same worktree, same protocol | ≥ 10 % → roll out |
| 6 | allocations (14.5 a) | microbenchmark | yes | an empty offload kernel under rocprofv3 | report upstream |
| 7 | cell fluxes once per cell (idea 5) | the largest | yes | only after the counters of row 3 | memory-bound → drop |

**Protocol (what made the numbers of §11–§13 trustworthy, and what spoiled the others):** pin one GPU (`ROCR_VISIBLE_DEVICES`) and use the same one for every build; `small` rocprofv3 three rounds **alternated** build by build, one run at a time (the node otherwise idle: concurrent runs inflated the `small` loop ×1.6 and moved `large` by 13 %); `large` `--tool none`, alone, at least twice; bitwise check `h5diff -c` at steps 0 and 10 for every width; then the 12-case regression old against new (`~/perf_regr2.sh`: needs `mpirun --oversubscribe`, the Slurm step has one slot), remembering that `baro-mhd-glm-hlld` differs from itself between runs; the FLUME and PRISM FNL regression suites against their goldens are still to be run before any of this is adopted.

---

## 15. Follow-up (2026-10-28 evening, 2026-10-29 night): counters, the Riemann kernels, forceinline, and the MPI transport

Same tree as §13 (`WENO_refactor`, `f8bca7e8`), run through two unattended Slurm jobs (74182 on TheraC17, 4 h; 76074 on TheraC17, 90 min) whose
stage scripts are in `~/perf_queue/`, results in `~/perf_results/`, run directories under the worktrees (`runs/`). "Old" is `origin/develop`
(`faa434ab`), "new" is `WENO_refactor`; unless said otherwise, one MI300X and the **system modules** of §11.2 (amdflang 24.0.0git).

### 15.1 What happened to the machine

- 2026-10-28 14:21:30: TheraC16–19 and most other nodes of the cluster went `Not responding` at the same second (cluster-wide event, not
  mine). TheraC18 was then drained by root ("Kill task failed (JobId=74154)", the previous user's interactive job could not be killed) and stayed drained
  until after 15:50. My job 74182 started at 16:27:41, on **TheraC17** (back in service with 8 GPUs), not on C18, and ran 16:27–20:04 (the runner
  released the node when its queue was empty). Job 76074 ran from 02:50:46 (its last stage ended at 03:50).
- Two mistakes of mine, both without consequence for the numbers: the CPU-lane runner exited at 17:14 when its queue was empty (a stage added later did not run until I restarted it),
  and the first 8-rank retry (`G70`) never ran because I gave `profile.sh` a tag as if it were an executable (rc 2, no MPI at all). I also did not
  watch the stages finish, so results sat unread for hours.

### 15.2 Hardware counters (rocprofv3, small case, 3 steps, 45 dispatches each, GPU 0)

| | face kernel old | face kernel new | dissipative kernel (both) |
|---|---|---|---|
| `MeanOccupancyPerActiveCU` | **1.00** | **1.97** | 1.92 |
| `OccupancyPercent` | 11.7 % | 22.7 % | 21.4 % |
| `GRBM_GUI_ACTIVE` (cycles per dispatch) | 927 M | **53 M** | 9.6 M |
| `SQ_INSTS_VALU` | 521 M | **185 M** | 103 M |
| `SQ_INSTS_VALU_{ADD,FMA,MUL}_F64` | 43.4 / 54.8 / 15.8 M | **the same** | 17.4 / 15.4 / 5.5 M |
| `SQ_WAVE_CYCLES` | 33.1 G | 3.54 G | 0.63 G |
| `SQ_WAIT_ANY` / `SQ_WAVE_CYCLES` | 96 % | 73 % | 58 % |

The occupancy of §14.3 is now **measured**: 1 → 2 waves per SIMD (the device ELF record agrees: `vgpr_count` 346 → 210). The FP64 instructions are identical (the arithmetic
is untouched), the other vector instructions fell from 407 M to 71 M (−83 %: the scratch and address traffic; the whole `SQ_INSTS_VALU` −65 %), and the waves wait 23 points less of their time (96 % → 73 %); the kernel is still latency-bound (73 % of the wave cycles waiting).

### 15.3 The Riemann kernels (`weno-riemann`)

Scaled-up `sod-hllc` (32³ cells per block, 512 blocks, 20 steps, rocprofv3, GPU 0, alternated), `compute_riemann_face_fluxes*` per dispatch (60 dispatches):

| round | old | new |
|---|---|---|
| 1 | 743.5 ms | 43.2 ms |
| 2 | 839.6 ms | 42.9 ms |
| 3 | 841.3 ms | 41.7 ms |

**17× to 20×**; registers 128/224 → 108/132, scratch 9616 → 3072 B. This is the first measurement of the Riemann path (§13.2 listed it as not measured).

Other WENO orders in the Riemann kernels (`sod-hllc`, `rj2a-hlld`, `baro-collide-hllc` × `weno-u-1,3,7,9`, `ngc = 5`, old against new, bitwise): `rj2a-hlld` identical for the four orders; `sod-hllc` and
`baro-collide-hllc` differ at their last checkpoint for the four orders. The **control** (old twice, new twice, old against new, `sod-hllc` u-1, u-7 and `baro-collide-hllc` u-3, GPUs 6,7): **old differs from its own rerun,
and new from its own rerun, at the same checkpoints**: non-determinism of the baseline in those cases, not a difference between old and new.

### 15.4 Regression inputs, all of them

54 case directories; 7 self-gravity cases skipped (no FFT executable on AMD); of the 47 run, mpirun -np 2, old against new, every checkpoint (`h5diff -c`) and the residual and conservation logs (`cmp`):
**43 bitwise identical**; 4 differ (`baro-mhd-glm-hlld`, `iso-mhd-tube2`, `sod-amr`, `sod-amr-2realm`) and **all four differ from their own old-executable rerun at the same checkpoint** (non-determinism of the baseline).
The golden references of the suite (`golden/fnl`, nvfortran) were not compared.

### 15.5 `!dir$ forceinline` (variant `inl`, 4 directive lines, worktree `adam-inl`, **not committed**)

- The directive works **on the call statement** in this compiler: a toy with a large callee and 6 call sites keeps 6 calls by default and with `!dir$ inline`, and has 0 with `!dir$ forceinline` before each call; placed inside the callee it does nothing.
- Placed before the call of the model adapter, the two WENO calls (plain and macro form) of the shared face-flux body, and the adapter's call of `compute_face_split_fluxes` (Euler only).
- Checkpoints identical to old (step 0 and 10, small). Face kernel, small, 3 alternated rounds: 3.37 / 3.25 / 3.38 → **2.86 / 2.91 / 2.92 ms (−13 %)**. Large loop, alone, 3 alternated rounds, nothing else running:
  new 52.63 / 52.66 / 51.36 s, inl 50.38 / 50.17 / 50.10 s (**−3.8 %**). (The first pair of large rounds of the same stage, run while a relink was going, was inconclusive.)
- Registers 128 + 152; the ELF record says `vgpr_count` **276** (above 256): back to **1 wave per SIMD** (new: 210, 2 waves). It is faster in spite of that, so capping the registers of this variant at 256 is an idea not tried.

Link-time inliner remarks (relink of `WENO_refactor` with `-Rpass-missed=inline`, 6.9 MB log): the device LTO threshold is **3712**; the per-width WENO routines cost **1730 (S=1, inlined), 4080 (S=2), 6960 (S=3), 10640 (S=4), 16140 (S=5)**:
not inlined from S = 2, which is why a directive is needed. The **generic** `weno_reconstruct_upwind_dev` and `..._wratio_dev` carry, in the log, **6 `_FortranACopyInAssign`, 6 `_FortranACopyOutAssign` and 2 `__kmpc_alloc_shared`** each:
Fortran array-section temporaries (copy-in / copy-out through the Fortran runtime, with the temporary allocated by the OpenMP device runtime, `malloc` → `gpu::allocate`).
That is a mechanism for the heap allocations of §14.5 (a) **in the generic WENO routine** and a further reason why merging its helpers into one routine pays. **No such remark** exists for `compute_face_split_fluxes`, the dissipative kernel, the per-width kernels or the per-width routines, so the two
`gpu::allocate` calls that every kernel of the image still has (255 of 512 kernels with exactly two call sites) are **not** explained by this; toys (a trivial kernel, `has_device_addr`, assumed-shape dummies, a private array passed to a non-inlined callee) did not reproduce them.

### 15.6 Multi-rank, and the MPI transport (the larger finding of this section)

**New (system-module) environment** (Open MPI 5.0.10-ucx1.22.0-ucc1.9.0): small case, 2–8 ranks, `profile.sh`: **2, 3, 4 ranks run, 5, 6, 7, 8 all abort** with `MPI_ERR_INTERN` in `MPI_Waitall`
(after ~24 s, ranks all started and loaded); with the cluster's UCX settings (`OMPI_MCA_pml=ucx`, `UCX_TLS=self,sm,rocm`, `HSA_ENABLE_IPC_MODE_LEGACY=1`) and with `UCX_TLS=self,sm,rocm_copy,rocm_ipc` the same. The threshold is at the number of GPUs of one socket (4 of the 8; **hypothesis**, not tested).
One run of the old tree on 2 ranks in that environment hung at start-up (killed by the 20-minute timeout), the same environment ran the new tree on 2 ranks.
Where it runs, large case, GPUs 0–3: **4 ranks: new 15.08 s, old 173.47 s (11.5×)**; new, 2 ranks 27.31 s; new, 1 rank 52.6 s (speedups 1.9×, 3.5×).

**Original environment** (`spack-setup-env.sh`; `rocm-afar-drop/23.2.1-gfx942`, `rocm-openmpi`, `rocm-hdf5-ompi-5.0.10/1.14.6`; both trees rebuilt there with `thera_amd`, executables `exe/fnl-amd-oldenv`, worktrees `adam-oldhead`, `adam-oldwref`): 8 ranks run.
One rank, large: **old 483.2 s, new 53.8 s (9.0×)** (amdflang 23.2.0: the gain is not the property of one compiler).
Eight ranks, large, Open MPI **default** transport: old 280.9 / 292.3 s, new **208.7 / 209.0 s**; the forest timing shows imbalance 0.67 s, reflux 0.35 + 0.14 s, seam 0: about 200 s of the 204 s of steps are in no reported category.
Eight ranks with the **ROCm-aware UCX** transport (`OMPI_MCA_pml=ucx`, `OMPI_MCA_osc=ucx`, `UCX_TLS=self,sm,rocm`, `HSA_ENABLE_IPC_MODE_LEGACY=1`, `UCX_MAX_RNDV_RAILS=1`; or `UCX_TLS=self,sm,rocm_copy,rocm_ipc`):

| 8 ranks, large | loop |
|---|---|
| new, default transport | 208.7 s, 209.0 s |
| new, UCX, `rocm` | **9.20 s** |
| new, UCX, `rocm_copy,rocm_ipc` | **9.23 s** |
| new, default transport, `HSA_ENABLE_SDMA=0` | 212.7 s |
| **old, UCX** (2 rounds) | **68.97 s, 68.99 s** |
| **new, UCX** (2 more rounds) | **9.34 s, 9.25 s** |

So **the default transport costs ~200 s of a 209 s loop (23×)**, in line with the README of the speedup matrix ("the FNL halo exchange posts MPI directly on device buffers ... the MPI must be ROCm-aware, Open MPI + UCX with ROCm"):
without the UCX settings the run is correct and slow. With them, across the node the new tree is **7.4×** faster than the old one (69.0 → 9.3 s), and the parallel efficiency against one rank (original environment) is
old 483 → 69 s (7.0×, 88 %), new 53.8 → 9.3 s (5.8×, 72 %): with 8 ranks each rank holds an eighth of the cells, so the kernel time per rank shrinks about eightfold while the host gaps between launches and the regrids (§15.7: about 12 s of the loop at one rank) are not expected to shrink as fast; that is a **hypothesis** for the 72 %, no 8-rank trace was taken.

### 15.7 Where the large loop's time goes outside the kernels (rocprofv3 kernel trace of the large `wenow` run, one rank)

Trace 72.8 s, kernels 44.1 s, 6037 gaps totalling 28.6 s: one 15.9 s gap before the first update_ghost (initialisation, outside the loop), two ~1.95 s gaps before `count_nonfinite_dev` (the regrids at steps 10 and 20), the rest
**about 8.8 s in ~6000 small gaps** (4770 of 0.1–1 ms = 3.6 s, 1100 of 1–10 ms = 3.0 s, 161 of 10–100 ms = 2.0 s), before kernels such as `apply_reflux_face_dev`, `blend_inadmissible_ghosts_dev`, `set_boundary_conditions_dev`, `compute_flux_difference_dev`: host work between launches, 15 % of the loop now.

### 15.8 What this changes in §14

- §14.3 (occupancy): **settled**, the per-width WENO change already gave 2 waves per SIMD (1.97 measured); no launch-bounds attribute is needed for it. It reopens for the `inl` variant (276 registers, 1 wave).
- §14.1 (inlining): the cause is the **cost against a threshold of 3712** and a directive on the call statement removes it (toy, and `inl` in the real kernel: −13 % small, −3.8 % large); textual inlining is not needed.
- §14.5 (a) (allocations): partly explained (Fortran copy-in/out temporaries of the generic WENO routine); the two allocations every kernel still has are **not** explained.
- §14.6 table: rows 1, 2, 3 (counters), 6 (partly), and the multi-rank question of §9.5 are done; rows 4 (per-width `compute_face_split_fluxes`), 5 (lean) and 7 (cell fluxes) are **not** done and are now second-order next to the host gaps (§15.7) and, at 8 ranks, the transport settings.

### 15.9 Not done

- The per-width `compute_face_split_fluxes`, the lean kernel, the cell-flux precomputation (ideas 3 and 5 of the five; no experiment was made).
- Whether the socket boundary (4 GPUs) is the cause of the new environment's failure at 5+ ranks; whether a different Open MPI of the system tree runs 8 ranks.
- The NVHPC/OpenACC build; the golden references; the 7 self-gravity cases; an FFT executable on AMD.
- The origin of the two allocations in every kernel; the host gaps of §15.7; the gravity-edge kernels (still only compiled).
- The `inl` variant is **not committed** (worktree `adam-inl`, 4 directive lines); the notes of §11–§15 and the stage scripts are in `~` and the repository root only.

---

## 16. Back to the original environment: the dissipative kernel's memory traffic, the split-flux experiments (per-width and lean), and corrections to §15

From 2026-10-29 the work runs in the **original environment** (§8: `spack-setup-env.sh`, `rocm-afar-drop/23.2.1-gfx942`, `rocm-openmpi`, `rocm-hdf5-ompi-5.0.10/1.14.6`, amdflang **23.2.0**, `--varset "thera_amd gfx942"`),
with the ROCm-aware UCX settings of §15.6 for any multi-rank run; the baseline of every comparison is `WENO_refactor` (`a622e71b`) built in that environment (worktree `adam-oldwref`, `exe/fnl-amd-oldenv`).
Slurm job 76080 (04:44 to about 07:25, TheraC17; the runner released the node when its queue was empty); scripts `~/perf_queue`, results `~/perf_results2/`.

### 16.1 Corrections to §15 (found while preparing this section)

- **`shock-cylinder-ib` was never compared.** In every regression sweep of §15.4 (and in the 12-case lists before it) the case exits with code 137 and writes no checkpoint; my script printed "0 differences" for it
  (no file compared is not a difference). The cause is **not** the 400 s timeout I assumed and **not** any change of this branch: the case dies with `Memory allocation failed` (`adam_ib_object.F90:467`, `ib_object%compute_phi`) **on origin/develop too**,
  at 2 ranks (each rank holds 500–570 GB of host memory after a minute, 1443 GB of the node's 1511 GB in use within 90 s, then SIGKILL), alone, and at 8 ranks (peak 1488 GB). It cannot run on this node. The count of §15.4 is therefore
  **42 bitwise identical, 4 non-deterministic baseline cases, 1 not run (`shock-cylinder-ib`), 7 skipped (self-gravity)**, not "43 of 47". The commit message of `a622e71b` (and section 15.4 as pushed) carry the old count.
- The −2.9 % of the large loop that the flag-switched lean build (v5) first showed (§16.3, table) was inflated by a noisier baseline; the clean lean-only build measures **−0.3 %** (§16.3).
- §14.5 (c) and §14.6 row 7 said the dissipative kernel "does not have the face kernel's problem": true of amdflang 24.0.0git (272 B of scratch, inlined), **not** of 23.2.0 (§16.2).

### 16.2 The dissipative kernel and the face kernel: memory counters (amdflang 23.2.0, small case, 3 steps, 45 dispatches each, GPU 0, one counter per run)

| | face kernel | dissipative kernel |
|---|---|---|
| scratch per thread (grid 311296 threads) | 2972 B (925 MB) | **21744 B (6.8 GB)** |
| fetched / written per dispatch | 3.63 / 3.15 GB | 3.25 / 1.83 GB |
| duration under counters | 3.78 ms | 1.14 ms |
| achieved traffic | 1.8 TB/s | **4.4 TB/s** (MI300X peak ≈ 5.3 TB/s) |
| bytes read / written per thread | 11.7 / 10.1 KB | 10.4 / 5.9 KB |
| L2 hit rate (`TCC_HIT/(HIT+MISS)`) | 37.9 % | 40.3 % |
| `VALUUtilization` | 97 % | 91 % |
| `SQ_WAIT_INST_ANY` | 1.51 G | 0.12 G |

The useful output of one dispatch is about 12.5 MB (40 B per thread): both kernels move two orders of magnitude more data, mostly scratch. **With 23.2.0 the dissipative kernel is memory-bound on its own spills** (21.7 KB of scratch per thread, 4.4 TB/s);
with 24.0.0git it has 272 B and no calls (§14.5 c). `MemUnitBusy`, `L2CacheHit` and `TA_BUSY_sum` returned no values. The durations are those of the counter-collection runs (serialised dispatches), the bandwidths are estimates.
Which private arrays make up the 21.7 KB (the stencil arrays of `add_dissipative_fluxes_o4_dev`, sized by their maxima) is **not** examined.

### 16.3 Experiments 4 and 5 (Euler only; `compute_face_split_fluxes`)

**4, per-width split:** the body of the Euler split routine in an include, instantiated for the generic routine and for S = 1 … 5 (S a named constant), per-width adapters, the `_s<N>` face kernels call their width.
**5, lean:** `compute_face_split_fluxes_lean` reads each stencil cell from `q_gpu`/`q_aux_gpu` into one-cell buffers when needed (the stencil is walked twice) instead of receiving the gathered `qs`, `qas` (1040 B per thread);
same operations in the same order. Both switched by build flags in one tree (`FNL_EXP_SPLIT_WIDTH`, `FNL_EXP_LEAN`; worktrees `adam-v4`, `adam-v5`, `adam-v45`), a clean lean-only version in `adam-lean` (4 files, +107 lines, additive).

| | face kernel, small, 3 alternated rounds (ms) | mean | small loop (s) | large loop, 3 alternated rounds (s) | scratch |
|---|---|---|---|---|---|
| baseline | 3.683 / 3.681 / 3.855 | 3.74 | 3.466 / 3.416 / 3.398 | 56.49 / 56.03 / 54.44 (55.66) | 2972 B |
| v4 per-width | 3.963 / 3.938 / 3.945 | 3.95 (**+5.6 %**) | 3.569 / 3.413 / 3.531 | 58.21 / 58.30 / 57.26 (57.92, **+4.1 %**) | 3260 B |
| v5 lean | 3.287 / 3.119 / 3.256 | 3.22 (**−13.9 %**) | 3.359 / 3.260 / 3.349 | 54.54 / 54.75 / 52.84 (54.04, −2.9 %) | 2140 B |
| v45 both | 3.048 / 3.072 / 3.079 | 3.07 (**−18.0 %**) | 3.247 / 3.299 / 3.222 | 54.95 / 52.36 / 52.46 (53.26, −4.3 %) | 2140 B |

**Clean lean-only build** (alone on GPU 0, alternated with the baseline, bitwise identical at steps 0 and 10): face kernel, small, 3 rounds: baseline 3.621 / 3.637 / 3.637, lean 3.245 / 3.236 / 3.232 ms (**−10.8 %**); small loop 3.426 → 3.254 s (−5.0 %);
**large loop, 4 rounds: baseline 54.243 / 54.333 / 54.368 / 54.275 s, lean 54.130 / 54.182 / 54.111 / 54.162 s: 54.30 → 54.15 s, −0.29 %** (−0.11 to −0.26 s per round, same sign every round).
The baseline of the first large-loop comparison (55.66 s, rounds 56.5 / 56.0 / 54.4) was 2.5 % slower and noisier than that of the clean comparison (54.30 s, spread 0.2 %) twenty minutes of node time apart: the −2.9 % / −4.3 % of the table are partly that.

**Reading.** Lean removes 11–14 % of the face-kernel time in the small case and 832 B of scratch, but about **0.3 % of the large loop**: at that size the face kernel is no longer where the time goes (the dissipative kernel under 23.2.0, §16.2, the ghost updates and the host gaps, §15.7).
The per-width split is slower on its own (+5.6 % kernel, +4.1 % large loop) and on top of lean adds ≈ 1.4 % at the large size for five routines and five adapters. **Neither is committed.** The forceinline variant of §15.5 is in `stash@{0}` of the repository.
Patches of the experiments: `~/perf_patches/lean.patch` (clean, the candidate), `~/perf_patches/split-width-and-lean.patch` (flag-switched, v4/v5/v45).

**Bitwise checks.** Small case, 10 steps: v4, v5, v45 and the clean lean identical to the baseline at steps 0 and 10. Regression inputs, `mpirun -np 2` with the UCX settings, baseline against **v45**: of 54 case directories, 7 self-gravity skipped,
**42 identical (every checkpoint and the residual and conservation logs)**, 4 differing (`baro-mhd-glm-hlld`, `iso-mhd-tube2`, `sod-amr`, `sod-amr-2realm`) that **also differ from a rerun of the baseline** (non-determinism), 1 not run (`shock-cylinder-ib`, §16.1);
v4 and v5 alone identical on `couette`, `sod-x`, `vortex-periodic`. Clean lean, 12-case subset (`couette`, `sod-x`, `sod-y`, `sod-z`, `vortex-periodic`, `vortex-rk3`, `amr-periodic-reflux`, `becker-shock`, `blast-limiter`, `noh-de`, `brio-wu` and `shock-cylinder-ib`): the 11 compared cases identical, `shock-cylinder-ib` not run.

### 16.4 Where this leaves the list of §14 / §15

- Done: idea 3 (lean) and the per-width split of the split routine (14.5 b): **no worthwhile gain at the large size**, in either form. The kernel-level gain of lean (−11 %) does not reach the loop.
- Open and now better targeted: the dissipative kernel's 21.7 KB of scratch under 23.2.0 (§16.2), the host gaps and regrids (§15.7), the unexplained allocations per kernel (§15.5), the parallel efficiency at 8 ranks (§15.6).
- Not done: idea 5 (cell fluxes once per cell); the NVHPC build; the golden references; self-gravity cases; and `shock-cylinder-ib` (host-memory exhaustion of the case itself, origin/develop included).
