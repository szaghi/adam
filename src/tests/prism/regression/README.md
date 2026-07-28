# PRISM regression suite

CI-tuned regression tests for the PRISM (Maxwell-equation) solver. Each case
runs in roughly one minute and validates one backend against a committed golden
output.

> ## ⚠️ THE SUITE IS CURRENTLY BROKEN — every case fails at init
>
> As of `96420ae4` (2026-07-27) **all nine cases abort during initialisation**
> with:
>
> ```
> [mpi-00000]prism_fWLayer_object%initialize start
> [mpi-00000]error stop : failed to load [fWLayer].(width)
> ```
>
> **Cause.** Commit `8e05d363` ("debug of fWL and divergence correction")
> reworked `prism_fWLayer_object` to take a **physical** layer width
> (`[fWLayer] width`, a real) instead of the old **cell-count** `C` (an
> integer); `C` became a *derived* per-block/face array
> (`C_face = min(ni, ceiling(width/ds))`,
> `src/app/prism/common/adam_prism_fWLayer_object.F90:144`). The parser now
> reads `width` and **only** `width` (`:189`), and `load_from_file` is called
> unconditionally from `initialize` (`:86`) with `go_on_fail` defaulting to
> `.false.` — so a missing `width` is a hard `error_stop` **even when the case
> uses no layer at all**.
>
> None of the nine regression `input.ini` files were migrated: all still carry
> `C = ...` and none defines `width`. The CI job `regression-prism-cpu` has been
> red on every commit since.
>
> **Fix** (not yet applied — it changes case inputs and therefore golden
> validity, so it is a deliberate, reviewer-approved act):
>
> - For the eight cases with `C = 0` (no layer), replace it with `width = 0.0`.
>   The layer is inactive either way, so the physics — and the committed
>   goldens — are unchanged.
> - For `rmf-fwl` the layer is **active** (`C = 6` cells at `ni = 16` on
>   `[-0.16, 0.16]`, i.e. `ds = 0.32/16 = 0.02`), so the equivalent physical
>   width is `6 * 0.02 = 0.12`. Because `C` is now derived via `ceiling`, an
>   exact round-trip is not guaranteed on refined blocks — **the `rmf-fwl` FNL
>   golden must be re-verified, and probably recaptured**, after the migration.
>
> Until this is fixed, treat every statement below about "green" as describing
> the suite's *design*, not its current state.

This suite is the **structural-change regression baseline** referenced by
[issue #10][issue10] (forest-of-trees migration plan). Every step of that plan
must leave this suite green.

[issue10]: https://github.com/szaghi/adam/issues/10

## What is and isn't a regression case

The existing test trees under `src/tests/prism/{cpu,fnl}/` are **research /
development** cases — long integration times, full AMR, sized for physics
validation. They are **not** regression anchors and are not consumed by
`run.sh`.

Since the 2026-07 reorganization those trees are split by spatial scheme,
`src/tests/prism/<backend>/{fd,fv}/<family>/`, where `fd` is the finite-difference
(`fd_centered`) path and `fv` the finite-volume (`fv_centered`) one. Current
families:

| Family | Where | What it is |
|---|---|---|
| `RMF` / `rmf` | `cpu/{fd,fv}`, `fnl/{fd,fv}` | the 4-AC-coil rotating-magnetic-field workhorse the regression cases derive from |
| `BC_tests` | `cpu/fd`, `fnl/fd` | one directory per boundary kind — `Neumann`, `Silver-Muller`, `fW_Layer`, and the combined `fWLayer_Silver_Muller` |
| `PIC` | `cpu/fd`, `fnl/fd` | `physical_model = PIC` particle cases: single particle (RK and leapfrog), plasma column/cylinder, and the `elliptic_solver_tests/` sub-family (incl. the `PEC` validation case) |
| `divergence_correction` | all four | the only cases exercising `constrained_transport = DB` + `divergence_correction = hyperbolic` (the [issue #11](https://github.com/szaghi/adam/issues/11) `nv`/CT hazard the regression suite deliberately avoids) |
| `fWL` | all four | the isolated absorbing-layer case (`[fWLayer] width`, all six faces) |
| `plane_wave` | `cpu/fd` | plane wave, Gaussian pulse, and `gaussian_pulse_PEC` (1-D PEC reflection) |
| `Gaussian_pulse` | `cpu/fv`, `fnl/fv` | the FV twin of the pulse case (`[weno] scheme = weno-c-4`) |
| `Helicon_antenna`, `magnetic_nozzle` | `cpu/fd` | long-integration application cases (helicon coil, solenoid nozzle) |
| `infinite_wire` | `fnl/fd` | single rectangular coil, `fdv_order = 4` |

None of these has a `check.sh`, none is goldened, and `run.sh` never descends
into them — it iterates only the immediate subdirectories of *this* directory.
The `fv` subtrees and the `divergence_correction`, `fWL`, PIC-on-FNL and
PEC/Silver-Müller families are all **new since 2026-07** and are not yet
represented by any regression anchor.

### Boundary-condition kinds

`[bc_*] type` accepts (`src/app/prism/common/adam_prism_bc_object.F90:84-101`):
`extrapolation`, `Neumann`, `Dirichlet`, `Silver_Muller`, `periodic`,
`radiative`, `PEC` (also lowercase `pec`). Note the **underscore** spelling
`Silver_Muller` — the hyphen appears only in a directory name. There is no
`default` branch in the `select case`, so an unrecognised string leaves the BC
type undefined rather than erroring: spelling mistakes fail silently.

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

### Case inventory notes

Most cases are auto-discovered (any subdirectory with an `input.ini` and a
`golden/<backend>/`). Two are worth calling out because they exercise a
device-path a plain field case does not:

- **`rmf-fwl`** — the only case that turns on the **fWLayer** (six faces)
  together with four AC coils. It is the regression anchor for the
  fWLayer host→device transfer: on the FNL backend the fWLayer field is stored
  transposed (`(nb,i,j,k,3)`) and copied through FUNDAL's transposed HtoD path,
  which — unlike the plain q-field copy — is sensitive to how the device
  destination pointer is passed. A `C=30`/`ni=32` research variant crashed at
  `np>1` with `cuMemcpyHtoDAsync → CUDA_ERROR_INVALID_VALUE` (an nvfortran
  copy-in temporary on the lbound-remapped device dummy handed a host address to
  the async HtoD); this regression-sized case guards that fix.

  **Layer width.** The layer used to be configured as a cell count `C`, which
  had to stay `< ni` — the layer stamps cells `ni-C+1 .. ni`, so `C ≥ ni` drove
  the index negative. Since `8e05d363` the input key is a **physical** width
  (`[fWLayer] width`) and the per-block/face cell count is derived and clamped:
  `C_face = min(ni, ceiling(width/ds))`. The clamp makes the old negative-index
  footgun unreachable, but a `width` wider than a block still silently saturates
  the layer to the whole block — keep `width` well under the block extent so the
  case exercises a genuine partial layer.

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
| `golden/<backend>/*-residuals.dat` | Per-iteration residuals log               | tolerance-aware (`digest.py compare-residuals`) |

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

**Residuals log** (`digest.py compare-residuals`): compared **tolerance-aware**,
with the *same* `(rtol = 1e-6, atol = 1e-3)` calibration as the field digest —
the header line must match verbatim, integer columns exact, float columns within
tolerance (`run.sh` calls `compare-residuals`, not `diff`). Residuals are written
from a single rank and are not subject to reduction reordering, so within one
toolchain they are in fact bit-stable; the tolerance exists so a dev-workstation
run on a **different** compiler (this box currently runs gfortran-16-experimental,
CI runs the Ubuntu default) is checked on its significant digits rather than
failing on last-digit drift. The golden is still captured in CI (the single
source of truth for the CPU backend); the tolerance simply means a same-case
cross-compiler rerun is a meaningful check rather than a guaranteed failure.

> **Historical note.** Earlier revisions of this suite compared residuals
> byte-exact with `diff -q` and required the golden to be captured on the exact
> checking toolchain. The comparison has since moved to the tolerance-aware
> `compare-residuals` path; this section and the reference table above reflect
> the current code.

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

`run-fnl-local.sh` is the trigger. It does only what the FNL run needs beyond
the CPU run: loads the `nvhpc` module (nvfortran + bundled MPI), sets the
OpenMPI/UCX tuning (the `nvhpc` module's `mpirun` alias does not expand in a
non-interactive script, so the flags are passed as `OMPI_MCA_*` / `UCX_TLS`
env vars instead), then calls `run.sh fnl --varset local_nvf`. The nvf HDF5
prefix and `NVF_CC` come from the `local_nvf` fobos varset — the wrapper does
not export them. `run.sh` itself stays backend-agnostic.

`run.sh` runs every case under `mpirun -np 2`, the same rank count for both
backends — the golden is therefore a two-rank golden. On a multi-GPU box that
is one rank per GPU; the grid is small enough that two ranks can also share a
single GPU if needed (the FNL hardware floor is one GPU, not two). When
diagnosing FNL by hand it is common to run `mpirun -np 1` directly against the
exe — useful for isolating a crash, but note a one-rank run will **not** match
the committed two-rank golden.

> **FNL dependency note.** The FNL backend reaches the GPU through the
> vendored **FUNDAL** library (`src/third_party/FUNDAL`). FUNDAL is a
> `fobis fetch` dependency and its tree is gitignored — after any
> `fobis fetch --update`, re-verify the FNL run, because a dependency bump
> can re-introduce device-layer regressions that only the FNL backend
> exercises.

### Capturing / refreshing the FNL golden

The FNL golden is captured on a GPU workstation, never in CI. To capture a
new case's golden, or to refresh an existing one after a reviewer-approved
behaviour change:

```bash
./src/tests/prism/regression/run-fnl-local.sh   # runs, writes rmf/work-fnl/digest.txt
mkdir -p src/tests/prism/regression/<case>/golden/fnl
cp src/tests/prism/regression/<case>/work-fnl/digest.txt \
   src/tests/prism/regression/<case>/work-fnl/*-residuals.dat \
   src/tests/prism/regression/<case>/golden/fnl/
```

The digest is plain text and committable like the CPU golden. The "Capturing
and updating golden outputs" rules above apply equally to the FNL golden — a
refresh of an existing golden is a deliberate, reviewer-approved act, not a
silent rewrite.

Unlike the CPU golden (CI is the authority), the FNL golden's authority **is**
the workstation that runs `run-fnl-local.sh` — there is no GPU CI to capture
from. Whoever refreshes it must run on a known-good tree and review the digest
diff before committing.

## AMR seam cases (goldenless, `check.sh`-driven)

`rmf-amr`, `rmf-amr-fd` and `rmf-amr-fd-pulse` carry a deterministic
intra-realm 2:1 AMR jump (the `AMR_GEO` primitive-box marker) and are driven
by their own `check.sh`, not by `run.sh` goldens. `rmf-amr` asserts the seam
*structure* (registration, restriction, reflux plumbing, fv path);
`rmf-amr-fd` and `rmf-amr-fd-pulse` assert the seam *divergence* behaviour of
the fd_centered path.

**The historical acceptance "seam max|div(B)| ≤ 1e-13" is RETIRED** (issue
#21 §1, premise correction): no surveyed cell-centered method achieves
machine-zero seam div(B) without staggered storage, and the E0/E0′ gate
established that the residual seam div(B) under the tricubic ghost fill
(`[amr] seam_ghost_fill`, default) is an ordinary truncation quantity —
convergent under refinement at p_obs ≈ +1.2…+1.6 and bounded in time
(#21 N3.5). The pointwise round-off criterion applies only to the
**controls** (no seam) and to the pulse case's **div(D)**, which is a
structural zero (partition-of-unity fill of a z-invariant field).

Acceptance now in force (implemented in the two `check.sh`):

| assertion | case | meaning |
|---|---|---|
| control max\|div(D,B)\| ≤ 1e-13 | both | interior matched-stencil identity intact |
| seam max\|div(D)\| ≤ 1e-13 | pulse | structural invariant of the fill |
| seam max\|div(B)\| = pinned baseline ± 5% | both | golden-style: catches regressions AND silent improvements |
| `--convergence`: p_obs ≥ 0.8 on the 16→32 matched-time ladder | pulse only | the seam leak stays truncation-order |

Convergence is asserted only on the **pulse** case: the rmf coil filament has
div(J) ~ h⁻², which contaminates matched-time seam metrics on `rmf-amr-fd`
(#21 N3 — the case remains a valuable regression anchor, but not an order
measurement). Baselines are pinned in each `check.sh` header with provenance;
moving them is a deliberate, reviewer-approved rebaseline, exactly like a
golden refresh.

## Debug builds for diagnosis

`run.sh` builds the **release** modes (`prism-cpu-gnu`, `prism-fnl-nvf`) — that
is what the suite validates, because release is what ships. The debug modes
(`prism-cpu-gnu-debug`, `prism-fnl-nvf-debug`) are **not** driven by the
harness; they are a manual diagnosis tool.

When a regression run crashes and the release backtrace is unhelpful, rebuild
the debug mode and run the case by hand:

```bash
# FNL example — debug mode adds -Mbounds, -Mchkptr, -Ktrap=fp, -traceback, -O0
fobis build --mode prism-fnl-nvf-debug --varset local_nvf
cd src/tests/prism/regression/<case>/work-fnl   # or run from a scratch dir
mpirun -np 1 ../../../../../../exe/adam_prism_fnl input.ini
```

The debug mode catches several bug classes the release `-fast` build silently
tolerates — uninitialised pointers (`-Mchkptr`), out-of-bounds access
(`-Mbounds`), FP exceptions (`-Ktrap=fp`). It is the right tool for "release
crashes somewhere in a device routine" — but two caveats:

- **Stale objects.** FoBiS does not always recompile a `.F90` when an
  `#include`d `.INC` it depends on changes (notably the vendored FUNDAL
  templates). If a debug rebuild seems not to pick up a source edit, delete
  the specific stale objects under `exe/obj/nvf/debug/` (and `exe/mod/...`)
  and rebuild.
- **`-fast`-only bugs.** Some failures appear *only* in the release build —
  the optimiser exposes undefined behaviour the `-O0` debug build steps past.
  If debug passes but release crashes, the bug is real and optimisation-
  sensitive; do not assume "debug works" means "fixed".
