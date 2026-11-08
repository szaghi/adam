# Regression testing

Recipes for the two test layers: the **regression suites**, run on every change, which detect any change of a
committed result, and the **verification legs**, run for the part a change touches, which check the physics against
an oracle. The suites themselves (harness, goldens, tolerances, case inventories) are described in the
[regression reference](/reference/regression); this page gives the tasks.

## T1. Run a regression suite

**Problem.** Check that a change leaves every committed result unchanged.

**Commands:**

```bash
./src/tests/flume/regression/run.sh cpu                              # builds flume-cpu-gnu and flume-cpu-gnu-fft
./src/tests/flume/regression/run.sh fnl --varset local_nvf           # builds flume-fnl-nvf and flume-fnl-nvf-fft
./src/tests/prism/regression/run.sh cpu
./src/tests/prism/regression/run.sh fnl --varset local_nvf
./src/tests/flume/regression/run.sh cpu --no-build                   # reuse exe/
REGRESSION_CASES="sod-x brio-wu" ./src/tests/flume/regression/run.sh cpu --no-build   # only these cases
```

On a WSL2 workstation run the FNL suites through `run-fnl-local.sh [--no-build]` (it calls `run.sh fnl --varset
local_nvf` itself), which adds the WSL-only MPI settings ([GPU5](/cookbook/gpu#gpu5-run-on-a-development-box-with-wsl2-not-on-a-cluster)).

**What you get.** One line per case and a final `PASS: n FAIL: m SKIP: k`; the exit status is non-zero on any failure.
A case **fails** when its outputs differ from its goldens or when its run exits non-zero (a crash, or the hard timeout;
`REGRESSION_ALLOW_MPIRUN_FAILURE=1` downgrades that to a warning). A case **skips** when it has no golden for the
backend; in PRISM's count the two directories that are not cases (`reproducers/`, `slices/`) also show as SKIP. Each
case runs in its own `work-<backend>/` directory. Both suites start with the lints: the GPU race-shape check
(`check-gpu-race-shapes.sh`), then the removed-keys lint (T4).

**Reference.** [The harness](/reference/regression#the-harness); the usage comments at the top of each `run.sh`.

## T2. What a golden is, and what "pass" means

**Problem.** Know what the suite compares.

**Goldens** (`<case>/golden/<backend>/`, committed; the raw HDF5 checkpoints are not):
- `digest.txt`: for every checkpoint and every variable, `count, min, max, sum, sum_sq` over all blocks and ranks
  (`src/tests/prism/regression/digest.py`, shared by both suites);
- `*-residuals.dat`: the residual history;
- FLUME also `*-conservation_history.dat`: the volume integrals of the conserved variables.

**Pass** is a tolerance, not bitwise: `rtol = 1e-6` absorbs compiler, libm, FMA and reduction-order differences
between machines, `atol = 1e-3` is the floor for sums that are pure cancellation noise; integer columns and headers
must match exactly. Bitwise identity is checked separately where it is the contract: FLUME's OpenMP build against the
serial one (`run-omp-bitwise.sh`, in CI), RIPPLE's rank and CRLF legs (T6). A multi-realm case declared
`equivalent_to` a single-realm case is compared with it at the golden tolerance, not bitwise.

**Reference.** The tolerance comment in `digest.py`; [goldens and tolerance](/reference/regression#goldens-and-tolerance).

## T3. Add a case and capture its goldens

**Problem.** Turn an input into a regression anchor.

**Steps:**
1. Create `src/tests/<app>/regression/<case>/input.ini` (a frozen copy of a verification input, kept small: the FLUME
   suite runs every case on 2 ranks). The harness discovers every directory with an `input.ini`.
2. Run it without a golden and copy the outputs:

```bash
REGRESSION_RUN_GOLDENLESS=1 REGRESSION_CASES="<case>" ./src/tests/flume/regression/run.sh cpu --no-build
mkdir -p src/tests/flume/regression/<case>/golden/cpu
cp src/tests/flume/regression/<case>/work-cpu/digest.txt \
   src/tests/flume/regression/<case>/work-cpu/*-residuals.dat \
   src/tests/flume/regression/<case>/work-cpu/*-conservation_history.dat \
   src/tests/flume/regression/<case>/golden/cpu/
```

3. Repeat on the FNL backend for `golden/fnl/`, and add the case to the list in the `run.sh` header.
4. Check with `git status` that every golden file is staged (use `git add -f` if an ignore rule hides one): a golden
   that is not committed leaves the harness comparing the digest only.

**Refresh an existing golden** only after an intentional change, adjudicated case by case (which fields move, when, by
how much, and why) and recorded in the commit message. The run fails against the old golden but leaves its outputs in
`work-cpu/` (`work-fnl/` on the GPU):

```bash
REGRESSION_CASES=<case> ./src/tests/flume/regression/run.sh cpu --no-build
cp src/tests/flume/regression/<case>/work-cpu/{digest.txt,*-residuals.dat,*-conservation_history.dat} \
   src/tests/flume/regression/<case>/golden/cpu/
```

A failing case on a new machine is first a question about the machine (T7), not a reason to refresh.

**Reference.** [Goldens and tolerance](/reference/regression#goldens-and-tolerance), golden authority.

## T4. The removed-keys lint

**Problem.** Stop committed inputs from carrying keys the parser no longer reads.

**Command:** `./src/tests/lint/check-removed-input-keys.sh` (the regression suites run it first).

**What you get.** A failure for any committed PRISM, FLUME or AMR test input that holds a removed key:
`[fWLayer] C` (replaced by `width`), `[PIC] elliptic_correction`, the `[linear-algebra]` section. When you remove an
input key, add its rule here.

**Reference.** [Library input reference, parsing rules](/reference/library-input#parsing-rules).

## T5. Check that bad inputs stop at initialization

**Problem.** Keep a refused combination refused.

**Command:** `PRISM_EXE=$PWD/exe/adam_prism_cpu ./src/tests/prism/regression/refusals.sh` (or the FNL executable).

**What you get.** One leg per refusal: each rewrites the `rmf` input with `sed` and requires a non-zero exit **and**
the expected message in the log, so a crash is not mistaken for a refusal. Add a leg with every new refusal.

**Reference.** [PRISM hazards](/applications/prism/hazards#refused-combinations).

## T6. Verification legs and RIPPLE's oracles

**Problem.** Check the physics of the part a change touches.

**Commands:**

```bash
./src/tests/flume/verification/sod/check.sh --legs all          # every verification script takes --legs
./src/tests/ripple/check.sh                                      # RIPPLE: r1 r2 r3 r5 r6 (r6 also on test_ripple_crest_criterion)
RIPPLE_R5_EXE=$PWD/exe/adam_ripple_fnl ./src/tests/ripple/check.sh --legs r5   # the CRLF leg on FNL
```

**What you get.** `PASS`/`FAIL` per leg against an exact or reference solution. In the FLUME verification scripts a
selection that names no real leg exits 2 (`src/tests/flume/verification/legs.sh`); RIPPLE's `check.sh` stops on an
unknown leg too. Run the legs of the
part you changed, not all of them; the regression suite is the layer that runs on every change.

**Reference.** [FLUME verification](/applications/flume/verification/); RIPPLE's oracles in
`src/app/ripple/README.md`.

## T7. Bring up a new machine

**Problem.** Check that ADAM builds and runs correctly on a machine, compiler or GPU it has never run on.

**Steps:**
1. Toolchain, HDF5 and the `[varset:*]` block: [install, new machine](/start/install#new-machine).
2. The CPU gate, against the committed goldens:

```bash
./src/tests/prism/regression/run.sh cpu --varset mymachine_gnu     # 15 PASS, 11 SKIP
./src/tests/flume/regression/run.sh cpu --varset mymachine_gnu
```

3. The GPU gate. On a cluster, load the site's modules and call `run.sh` directly; `run-fnl-local.sh` is for the WSL2
   workstation only ([GPU5](/cookbook/gpu#gpu5-run-on-a-development-box-with-wsl2-not-on-a-cluster)). The FNL backend
   posts MPI directly on device buffers and has no host-staging fallback, so GPU-aware MPI is required.

```bash
module load nvhpc                                                  # or the site's modules
./src/tests/prism/regression/run.sh fnl --varset mymachine_nvf     # 16 PASS, 10 SKIP
```

4. The PRISM seam oracles, which `run.sh` does not invoke (`--build` builds `prism-cpu-gnu` first; `PRISM_EXE=` selects
   another executable, whose environment the caller sets):

```bash
cd src/tests/prism/regression
./rmf-amr/check.sh          --build     # refinement, seam registration, reflux
./rmf-amr-fd/check.sh       --build     # seam div(B) baseline, div(J) band
./rmf-amr-fd-pulse/check.sh --build     # source-free two-invariant oracle
./rmf-amr-fd-pulse/check.sh --convergence
PRISM_EXE=$PWD/../../../../exe/adam_prism_fnl ./rmf-amr-fd/check.sh
```

**Pass.** The committed goldens pass within tolerance: that is the port validation, and a new machine does not mint
goldens. The oracle bands are pinned (±5% on the seam div(B), a factor 3 on div(J)); the `p_obs` convergence assertion
is compiler-independent. A small miss is a tolerance question; a miss by orders of magnitude is a defect. Never widen a
band to make a machine pass.

**Troubleshooting** (the build-layer symptoms are in [install](/start/install#new-machine)):

| Symptom | Cause |
|---|---|
| `ERROR: mpirun not on PATH` | MPI module not loaded before `run.sh` |
| `ERROR: executable exe/adam_prism_* not found after build` | the build failed: read the build log; do not use `--no-build` |
| venv creation fails, or `h5py` import error | system Python lacks `venv`; delete `exe/.regression-venv` and retry |
| `no '*-*.h5' checkpoints produced` | the solver aborted: read the run output, try a debug build |
| SIGABRT in `ucp_proto_rndv_send_start` (FNL, np ≥ 2) | WSL2 rendezvous abort ([GPU5](/cookbook/gpu#gpu5-run-on-a-development-box-with-wsl2-not-on-a-cluster)) |
| digest mismatch around 1e-9 on stable fields | cross-compiler noise, inside `rtol = 1e-6`: not a failure |
| digest mismatch by orders of magnitude | a real regression: do not refresh the golden |

**Debug a crash.** `run.sh` builds the release modes. Rebuild the debug mode and run the case by hand on one rank (a
one-rank run isolates a crash; it does not match the two-rank golden):

```bash
fobis build --mode prism-fnl-nvf-debug --varset local_nvf
cd src/tests/prism/regression/<case>/work-fnl
mpirun -np 1 ../../../../../../exe/adam_prism_fnl input.ini
```

Debug adds `-Mbounds -Mchkptr -Ktrap=fp -traceback -O0` (nvfortran) or `-fcheck=all -ffpe-trap=invalid,zero,overflow
-fbacktrace` (gfortran). FoBiS does not always recompile a `.F90` when an `#include`d `.INC` changes (the FUNDAL
templates): if a rebuild ignores an edit, delete the stale objects under `exe/obj/prism/nvf/debug/` and
`exe/mod/prism/...`. Some failures appear only in release builds: a passing debug run does not mean fixed.

**Reference.** [Regression reference](/reference/regression); [install, new machine](/start/install#new-machine).
