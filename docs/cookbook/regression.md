# Regression testing

Recipes for the two test layers: the **regression suites**, run on every change, which detect any change of a
committed result, and the **verification legs**, run for the part a change touches, which check the physics against
an oracle. The suites are documented in full on the [PRISM regression page](/tests/prism-regression), its
[bring-up tutorial](/tests/prism-regression-tutorial) and the [FLUME regression page](/applications/flume/verification/regression);
this page gives the tasks.

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
case runs in its own `work-<backend>/` directory. Both suites start with the lints: the removed-keys lint (T4), and in FLUME's suite first the GPU race-shape check
(`check-gpu-race-shapes.sh`).

**Reference.** The usage comments at the top of each `run.sh`.

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

**Reference.** The tolerance comment in `digest.py`; the [PRISM regression page](/tests/prism-regression#two-reference-files-per-case).

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

**Rules.** Refreshing an existing golden is a deliberate act tied to an intended numerical change, adjudicated (what
moved, by how much, why) and recorded in the commit message; a failing case on a new machine is first a question about
the machine ([bring-up tutorial, step 7](/tests/prism-regression-tutorial#step-7-—-goldens-on-a-new-machine)).

**Reference.** [Refreshing a golden](/applications/flume/verification/regression#refreshing-a-golden).

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

**Reference.** [PRISM input hazards](/applications/prism/input#hazards).

## T6. Verification legs and RIPPLE's oracles

**Problem.** Check the physics of the part a change touches.

**Commands:**

```bash
./src/tests/flume/verification/sod/check.sh --legs all          # every verification script takes --legs
./src/tests/ripple/check.sh                                      # RIPPLE: r1 r2 r3 r5 on the CPU executable
RIPPLE_R5_EXE=$PWD/exe/adam_ripple_fnl ./src/tests/ripple/check.sh --legs r5   # the CRLF leg on FNL
```

**What you get.** `PASS`/`FAIL` per leg against an exact or reference solution. In the FLUME verification scripts a
selection that names no real leg exits 2 (`src/tests/flume/verification/legs.sh`); RIPPLE's `check.sh` stops on an
unknown leg too. Run the legs of the
part you changed, not all of them; the regression suite is the layer that runs on every change.

**Reference.** [FLUME verification](/applications/flume/verification/); RIPPLE's oracles in
`src/app/ripple/README.md`.
