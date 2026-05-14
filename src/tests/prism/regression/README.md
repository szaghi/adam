# PRISM regression suite

CI-tuned regression tests for the PRISM (Maxwell-equation) solver. Each case
runs in roughly one minute and validates one backend against a committed golden
output.

This suite is the **structural-change regression baseline** referenced by
[issue #10][issue10] (forest-of-trees migration plan). Every step of that plan
must leave this suite green.

[issue10]: https://github.com/szaghi/adam/issues/10

## What is and isn't a regression case

The existing test trees under `src/tests/prism/{cpu,fnl}/` are **research /
development** cases — long integration times, full AMR, sized for physics
validation. They are **not** regression anchors and are not consumed by
`run.sh`.

Cases under this directory (`src/tests/prism/regression/`) are derived from
research cases but **tuned for CI**:

- Wall-time budget: **soft target ~1 minute per case per backend**. Not a
  hard cap — CI runners are slower than dev workstations, and a per-step
  residuals log (`residuals_save = 1`) is a better regression signal than
  an aggressively short run. Tune `it_max` down only if a case starts
  dominating CI wall-time.
- Grid size kept small (e.g. `16³`).
- `it_max` reduced to roughly 10–20.
- `it_save` set so 2–3 checkpoints fall within `it_max` (e.g. 0, 5, 10).
- AMR runtime refinement **disabled** (`amr_iterations = 0`,
  `[amr] frequency = 999999`). Refinement decisions are
  float-comparison-sensitive across compilers and would break
  golden-output matching.
- **Initial uniform refinement is enabled** (`[amr] iu_ref_levels = 2`,
  `max_level = 2`). Setting `max_level = 0` would produce one root tree
  node total and leave each MPI rank with zero owned blocks, so the
  HDF5 writer emits an empty file. Uniform initial refinement gives
  64 leaf blocks (`ratio = 8`, level 2) — enough to distribute across
  the fixed two-rank run and exercise the kernels and exchange path.
  Because `frequency = 999999`, no refinement decisions occur during
  the run itself, so the case remains deterministic across compilers.
- Fixed `mpirun -np 2` — same across backends and across runs.

## Running the harness

The harness is `run.sh` in this directory — a backend-agnostic core. The CPU
backend runs in CI; the FNL backend is workstation-only (see "CPU vs FNL
backends" below) and has a convenience wrapper, `run-fnl-local.sh`, that
sets up the NVHPC environment for you.

```bash
# CPU backend — uses fobos's default varset (local_gnu)
./src/tests/prism/regression/run.sh cpu              # build, run, diff CPU golden
./src/tests/prism/regression/run.sh cpu --no-build   # skip the fobis build step

# FNL backend — GPU workstation only; the wrapper loads `nvhpc`, sets MPI
# tuning, and calls run.sh fnl --varset local_nvf.
./src/tests/prism/regression/run-fnl-local.sh
./src/tests/prism/regression/run-fnl-local.sh --no-build
```

Exit code 0 means every case passed; non-zero means at least one case failed.

### HDF5, NVF_CC and the fobos varset

HDF5 library paths and `NVF_CC` are **fobos variables**, defined per
`[varset:*]` in the repo `fobos` — not shell environment variables. The
harness does not set them and does not need them in the environment.

- **CPU** (`run.sh cpu`) builds with fobos's *default* varset, `local_gnu`,
  which resolves HDF5 to `lib/hdf5/develop/gnu/14.2.0`. CI relies on this
  default — no varset flag, no `HDF5_PREFIX`.
- **FNL** (`run.sh fnl`) needs `--varset local_nvf`, which resolves the
  nvf-built HDF5 and `NVF_CC`. `run-fnl-local.sh` supplies that flag.

`run.sh` also accepts a generic `--varset <name>` if a different fobos
varset is ever needed.

### Other prerequisites

| Tool / runtime | Required for | Notes                                                      |
|----------------|--------------|------------------------------------------------------------|
| `mpirun`       | both         | On `PATH`; for FNL the wrapper loads it via `nvhpc`        |
| `fobis` 3.8+   | both         | Long-form CLI; legacy `FoBiS.py -mode` not accepted        |
| `python3`     | both         | A venv with `h5py` is created on first run                 |
| `nvhpc` module | `fnl` only   | Loaded by `run-fnl-local.sh`; provides nvfortran + mpirun  |

On first invocation `run.sh` creates a private Python venv at
`exe/.regression-venv/` and installs `h5py` + `numpy` into it — `digest.py`
needs them, and WSL2 / CI system Python is externally managed. The venv is
gitignored (it sits under the already-ignored `exe/`) and reused across runs.

It deliberately lives under `exe/`, **not** under `src/`: `fobis` (the build)
and `formal` (the API-doc generator) both scan `src/` and would otherwise
compile or document the Fortran test fixtures shipped inside `numpy` —
breaking the build link step and producing junk doc pages.

## Reference data

A full PRISM checkpoint is ~100 MB per rank — six checkpoints across three
saved iterations is ~640 MB, far too large to commit. So the suite does **not**
commit raw HDF5. Each case commits two compact references instead:

| File                            | What it is                                   | Compared with        |
|---------------------------------|----------------------------------------------|----------------------|
| `golden/<backend>/digest.txt`   | Per-variable field digest (see `digest.py`)  | tolerance-aware      |
| `golden/<backend>/*-residuals.dat` | Per-iteration residuals log               | byte-exact (`diff`)  |

**The field digest** (`digest.py`) reduces every HDF5 checkpoint to a small set
of point-wise reductions — `count, min, max, sum, sum_sq` — aggregated per
field variable (all 32 blocks of `Bx` together, etc.). This is:

- **small** — tens of KB for the `rmf` case vs ~640 MB of raw HDF5;
- **point-wise sensitive** — any single changed cell shifts at least one
  reduction, which the aggregate residual norm could mask;
- **ownership-invariant** — aggregating across blocks/ranks means a pure
  block-redistribution (not a physics change) does not trip the digest. The
  residuals log is ownership-invariant for the same reason.

**Excluded variables.** The `div_*` divergence diagnostics (`div_D`, `div_B`,
`div_J`, and the numbered `div_05`..`div_12`) are *not* in the digest.
They are physical-constraint quantities that should be ~0 — their stored
values are pure round-off, so no reduction of them is reproducible across
compilers (the first CI run disagreed with the dev-box golden by ~24% on
`div_J` because both were noise). The exclusion is a name-prefix rule in
`digest.py`. If a divergence-cleaning regression is ever wanted it needs a
dedicated metric (e.g. ‖div‖ below a threshold), not a golden-value diff.

The raw HDF5 in `work-<backend>/` is the source the digest is computed from;
it is scratch and never committed.

## Adding a new case

```bash
mkdir -p src/tests/prism/regression/<case_name>
cp <a-research-case>/input.ini src/tests/prism/regression/<case_name>/input.ini
# tune per the rules above:
#   - reduce [grid] ni/nj/nk to 16
#   - [time] it_max = 10 (or 20 max)
#   - [IO] it_save = 5
#   - [IO] residuals_save = 1
#   - [IO] output_basename = <case_name>     # also the checkpoint file prefix
#   - [initial_conditions] amr_iterations = 0
#   - [amr] frequency = 999999, max_level = 2, iu_ref_levels = 2
```

`output_basename` doubles as the checkpoint file prefix the harness digests —
it must not collide with `restart_basename`, or the restart dump would be
swept into the digest.

The `work-<backend>/` directory (raw HDF5, the freshly computed digest, the
restart dump) is regenerated on every run and is gitignored — only
`golden/<backend>/digest.txt` and `golden/<backend>/*-residuals.dat` are
committed.

## Capturing and updating golden outputs

**CI is the golden authority for the CPU backend**, not a dev workstation.
The residuals log is compared byte-exact (see "Tolerance rationale"), so the
golden must be produced by the *same* toolchain that checks it — and the
checker is CI (Ubuntu, GCC 14). A golden captured on a dev box with a
different GCC will fail the residuals diff in CI even with no real change.
The FNL golden, by the same logic, is captured on the GPU workstation that
runs `run-fnl-local.sh` (see "Capturing the FNL golden").

To capture or refresh the CPU golden:

1. Push the branch (or open the PR). The `regression-prism-cpu` CI job runs
   `run.sh cpu` and — via an `always()` step — uploads the produced
   `work-cpu/digest.txt` and `work-cpu/*-residuals.dat` as the
   `prism-regression-cpu-work` artifact, even when the regression step fails.
2. Download that artifact from the run's summary page.
3. Copy its contents into the case's `golden/cpu/`:
   ```bash
   mkdir -p src/tests/prism/regression/<case_name>/golden/cpu
   cp <downloaded>/digest.txt <downloaded>/*-residuals.dat \
      src/tests/prism/regression/<case_name>/golden/cpu/
   ```
4. Commit. Re-run CI — the regression job now diffs against the committed
   golden and must pass.

**Goldens never change silently.** A refresh of an *existing* golden is by
definition a behaviour change and must be a deliberate, reviewer-approved act:

1. Identify the source change that *intentionally* alters numerical output —
   a bug fix, a numerical-scheme tightening, an FMA-flag change.
2. Open a PR titled `regression: bump goldens — <one-line reason>`, its
   description referencing the source-change PR.
3. Capture the new golden from that PR's CI run (steps 1–3 above).
4. Commit `golden/<backend>/digest.txt` and `golden/<backend>/*-residuals.dat`.
   The digest is plain text, so the PR diff makes every changed reduction
   explicit; a reviewer must agree the numerical change is correct.

Never bypass this with `--no-verify` or silent rewrites.

## Tolerance rationale

Two references, two comparison strategies. Both are calibrated for
**cross-compiler** runs — the golden is captured in CI, and the suite is
re-run in CI on the same image but also, manually, on dev workstations with
different GCC versions. Same-machine reruns are bit-identical; the tolerances
exist for the machine-to-machine case.

**Field digest** (`digest.py compare`): `count` must match exactly; `min`,
`max`, `sum`, `sum_sq` are compared with `numpy.isclose` at `rtol = 1e-6`,
`atol = 1e-3`.

- `rtol = 1e-6` absorbs floating-point association noise from differing
  compiler codegen, libm implementations, FMA contraction and MPI reduction
  order. The first CI run (gcc-14) vs a gcc-16 dev golden disagreed by
  `~1e-9` on stable fields — well inside `1e-6`, but far outside the
  original, naively tight `1e-11`. `1e-6` is still ~9 significant digits:
  any real algorithmic change moves a reduction far more than that.
- `atol = 1e-3` is the floor for **cancellation-residue** reductions. An
  antisymmetric field like `Jz` (±3.7e3 per cell, summing to ~4e-5 over
  340 k cells) has a `sum` that is pure cancellation noise — `1e-3` treats
  it as indistinguishable from zero. A genuine non-cancelling aggregate
  like `Jx` `sum` (~4.6e8) is utterly unaffected by a `1e-3` floor.

**Residuals log** (`diff -q`): compared **byte-exact**. Residuals are written
from a single rank, formatted decimally, and not subject to reduction
reordering. Byte-exact only holds **within one toolchain** — the residuals
golden must therefore be captured on the *same* toolchain that checks it
(CI captures and CI checks). Across GCC versions the low-order digits drift;
a dev-workstation run with a different compiler is expected to fail the
residuals diff and that is not a regression. If cross-compiler residual
checking is ever needed, the residuals comparison must move to a
tolerance-aware diff like the digest — it is byte-exact today only because
CI is the single source of truth.

If the digest tolerance proves too tight or too loose on a specific case,
`digest.py compare` accepts `--rtol` / `--atol` overrides; wire a per-case
override into `run.sh` when the second case demands it.

## CPU vs FNL backends

The two backends serve complementary purposes for the
[forest-of-trees migration plan][issue10]:

- **`prism-cpu-gnu`** — GNU/CPU baseline. Stable, fast, runs in CI on every
  PR. Catches host-side composition regressions.
- **`prism-fnl-nvf`** — OpenACC on nvfortran. Catches kernel
  chain-resolution regressions (the bug class the singleton refactor was
  built to avoid).

A case is fully covered only if its goldens exist for **both** backends.

### Why FNL is workstation-only

GitHub-hosted runners have **no GPU** — not just "no multi-GPU", zero GPU
on any plan — so `prism-fnl-nvf` cannot run as a CI job. The FNL backend is
therefore a **manually-invoked, workstation-local gate**: run it deliberately
before an "important" push, not on every PR. There is no automation deciding
what counts as "important" — that is the operator's judgement.

> **FNL status (as of Step 0): golden intentionally deferred.** The
> `prism-fnl-nvf` backend has known bugs under debugging, so `golden/fnl/`
> is deliberately *not* captured yet — pinning a golden to broken output
> would just enshrine the bug. The harness and `run-fnl-local.sh` are
> wired and ready; capture the FNL golden once the backend is fixed (it
> becomes load-bearing at Step 4 of the migration plan, the FNL kernel
> refactors). An empty `golden/fnl/` is expected, not an oversight.

`run-fnl-local.sh` is the trigger. It does only what the FNL run needs beyond
the CPU run: loads the `nvhpc` module (nvfortran + bundled MPI), sets the
OpenMPI/UCX tuning (the `nvhpc` module's `mpirun` alias does not expand in a
non-interactive script, so the flags are passed as `OMPI_MCA_*` / `UCX_TLS`
env vars instead), then calls `run.sh fnl --varset local_nvf`. The nvf HDF5
prefix and `NVF_CC` come from the `local_nvf` fobos varset — the wrapper does
not export them. `run.sh` itself stays backend-agnostic.

The case runs `mpirun -np 2`; on a multi-GPU box that is one rank per GPU.
The grid is small enough that two ranks could also share a single GPU if
needed — the FNL hardware floor is one GPU, not two.

### Capturing the FNL golden

The FNL golden does not exist until captured once on a GPU workstation:

```bash
./src/tests/prism/regression/run-fnl-local.sh   # runs, writes rmf/work-fnl/digest.txt
mkdir -p src/tests/prism/regression/rmf/golden/fnl
cp src/tests/prism/regression/rmf/work-fnl/digest.txt \
   src/tests/prism/regression/rmf/work-fnl/*-residuals.dat \
   src/tests/prism/regression/rmf/golden/fnl/
```

The digest is ~37 KB of plain text — committable like the CPU golden. The
"Updating golden outputs" rules above apply equally to the FNL golden: it
changes only by a deliberate, reviewer-approved act.
