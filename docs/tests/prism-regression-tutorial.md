# Running the PRISM regression suite on a new machine

A bring-up guide for standing the PRISM regression suite up on a machine it has
never run on — a new workstation, a new cluster, a new compiler, or a new GPU
architecture.

The [suite overview](./prism-regression) explains *what* the cases assert and
*why*; the suite's [`README.md`](https://github.com/szaghi/adam/blob/develop/src/tests/prism/regression/README.md)
is the operational reference for golden capture and tolerances. **This page is
the bring-up path**: from a bare checkout to a green `run.sh cpu`, then
optionally to a green FNL/GPU run.

::: danger Known-broken as of `96420ae4` (2026-07-27)
Every case currently aborts at initialisation with
`error stop : failed to load [fWLayer].(width)`, because commit `8e05d363`
changed `[fWLayer]` from a cell-count `C` to a physical `width` without
migrating the nine case `input.ini` files. `regression-prism-cpu` has been red
in CI since.

**This does not block a port.** Use the suite as your build/toolchain
validation — reaching the `[fWLayer].(width)` error already proves the
toolchain, HDF5, MPI and launch path are correct. Digest comparison only
becomes meaningful once the inputs are migrated.
:::

## The dependency stack

The suite sits on four layers. A failure at layer *n* usually shows up as a
confusing error at layer *n+1*.

```
4. regression harness  run.sh -> mpirun -np 2 -> digest.py (h5py venv)
                                      |
3. ADAM build          fobis build --mode prism-cpu-gnu [--varset <name>]
                                      |
2. varset (fobos)      $HDF5_PREFIX, $NVF_CC     <- the porting seam
                                      |
1. toolchain           gfortran|nvfortran + MPI + parallel HDF5 + python3
```

Layer 2 is the only machine-specific one. Porting is, almost entirely, adding
one `[varset:*]` block to `fobos`.

## Step 0 — Prerequisites

| Tool | Required for | Minimum | Check |
|---|---|---|---|
| `gfortran` | CPU backend | 9.x (CI uses **14**) | `gfortran --version` |
| `nvfortran` | FNL backend | NVHPC 22.x+ | `nvfortran --version` |
| MPI | both | OpenMPI ≥4.0 / MPICH ≥3.3 | `mpirun --version` |
| parallel HDF5 | both | ≥1.10, **same MPI as the app** | Step 1 |
| `fobis` | both | **3.8+** (long-form CLI) | `fobis --version` |
| `python3` | both | 3.8+, with `venv` | `python3 -m venv --help` |

```bash
pip install FoBiS.py          # the build tool; CLI binary is `fobis`
git clone https://github.com/szaghi/adam && cd adam
fobis fetch                   # vendored deps (PENF, FiNeR, FUNDAL, ...)
```

`fobis fetch` is mandatory on a fresh clone: `src/third_party/` is gitignored,
and the FNL backend reaches the GPU through the **FUNDAL** library fetched here.

::: warning FoBiS 3.8+ CLI form
Use `fobis build --mode X --varset Y`. The legacy `FoBiS.py build -mode X`
short-dash form is **not accepted** by 3.8+; every script in this repo emits the
long form.
:::

### Environment modules

On a module-managed box (including the reference WSL2 dev box), load the
compiler and MPI **before** anything else — `fobis` shells out to `mpif90`, and
a missing module surfaces as
`FileNotFoundError: [Errno 2] No such file or directory: 'mpif90'`:

```bash
module avail                          # see what the site provides
module load openmpi/5.0.7-gnu14.2.0   # CPU: MPI + matching GCC
module load nvhpc                     # GPU: nvfortran + bundled MPI
```

Match the MPI module to the compiler the HDF5 build used — that pairing is what
`$HDF5_PREFIX` encodes.

## Step 1 — Parallel HDF5

The most common bring-up failure. HDF5 **must be built against the same MPI**
the application links, or you get link-time symbol errors or run-time hangs.

Check for an existing parallel HDF5 first:

```bash
h5pfc -show 2>/dev/null || module avail hdf5 2>&1 | head
```

If one exists (`module load hdf5-parallel/...` on most clusters), note its
prefix and go to Step 2. Otherwise build one with the bundled helper, which also
fetches and builds szip and zlib:

```bash
./scripts/hdf5_build.sh -build \
   -hdf5 $PWD/lib/hdf5/develop/gnu/14.2.0 \
   -lsrc $PWD/untracked/hdf5-src \
   -abs-paths
```

`hdf5_build.sh` compiles with `CC=mpicc FC=mpif90`, so **whatever MPI is on
`PATH` when you run it is the MPI the app must use**. Load your modules first.
For an NVHPC build, `module load nvhpc` so `mpicc`/`mpif90` resolve to the
NVHPC-bundled MPI. `./scripts/hdf5_build.sh -h` lists all options
(`-use-autot` swaps CMake for autotools).

The prefix you pick here is `$HDF5_PREFIX` in Step 2. The in-tree convention is
`lib/hdf5/develop/<sdk>/<version>`, but any path works.

## Step 2 — Add a varset (the porting seam)

`$HDF5_PREFIX` and `$NVF_CC` are **fobos variables**, resolved from the active
`[varset:*]` in the repo `fobos`. They are *not* shell environment variables —
exporting them does nothing. This is the one file you edit to port.

Current varsets (`fobos`; `[varsets] default = local_gnu`):

| Varset | `$HDF5_PREFIX` | `$NVF_CC` | Machine |
|---|---|---|---|
| `local_gnu` | `lib/hdf5/develop/gnu/14.2.0` | — | dev workstation, GNU (**default**) |
| `local_nvf` | `lib/hdf5/develop/nvf/26.1` | `cc89` | dev workstation, Ada GPU |
| `leonardo` | Spack path, NVHPC 24.5 | `cc80` | CINECA Leonardo, A100 |
| `iac_gnu` | `lib/hdf5/develop/iac/gnu/12.2.1/openmpi-4.1.6` | — | IAC cluster, GNU |
| `iac_nvf` | `lib/hdf5/develop/iac/nvf/25.5` | `cc80` | IAC cluster, A100 |
| `spacehpc` | `/opt/cray/pe/hdf5/1.14.3.3/nvidia/23.3` | `cc90` | SpaceHPC (Cray), H100 |

Append a block for your machine:

```ini
[varset:mymachine_gnu]
$HDF5_PREFIX = /path/to/parallel/hdf5      # from Step 1

[varset:mymachine_nvf]
$HDF5_PREFIX = /path/to/nvf/parallel/hdf5
$NVF_CC      = cc90                        # your GPU's compute capability
```

### Picking `$NVF_CC`

`$NVF_CC` expands into `-gpu=$NVF_CC` in the nvfortran templates. Read it off
the device:

```bash
nvidia-smi --query-gpu=name,compute_cap --format=csv
# "NVIDIA GeForce RTX 4090, 8.9"  ->  cc89
```

Map `X.Y` → `ccXY`: A100 `8.0`→`cc80`, Ada/RTX-40 `8.9`→`cc89`, H100
`9.0`→`cc90`, Blackwell `10.0`→`cc100`.

::: danger NVF modes require an explicit varset
`$NVF_CC` is defined **only** in NVF varsets, never in the `local_gnu` default.
Building any `*-fnl-nvf` mode without `--varset <an nvf varset>` leaves
`-gpu=$NVF_CC` unexpanded and the build dies on the first nvfortran call.
:::

## Step 3 — Build

```bash
fobis build --lmodes                                    # list every mode
fobis build --mode prism-cpu-gnu                        # default varset
fobis build --mode prism-cpu-gnu --varset mymachine_gnu # explicit varset
```

PRISM modes (`fobos.d/prism.fobos`):

| Mode | Backend | Notes |
|---|---|---|
| `prism-cpu-gnu` | CPU / gfortran | what `run.sh cpu` builds; the CI-gated baseline |
| `prism-cpu-gnu-debug` | CPU / gfortran | `-O0 -fcheck=all -ffpe-trap` — manual diagnosis |
| `prism-fnl-nvf` | OpenACC / nvfortran | what `run.sh fnl` builds |
| `prism-fnl-nvf-debug` | OpenACC / nvfortran | `-O0 -Mbounds -Mchkptr -Ktrap=fp` — manual diagnosis |
| `prism-fnl-nvf-spacehpc` | OpenACC / nvfortran | Cray variant: links `hdf5_fortran hdf5 z` as `ext_libs` |

A successful build produces `exe/adam_prism_cpu` (or `exe/adam_prism_fnl`). A
harmless `requires executable stack` linker warning on
`adam_tree_bucket_object.o` is expected.

::: tip Cray / system-HDF5 machines
`prism-fnl-nvf-spacehpc` exists because Cray's HDF5 must be named explicitly
(`ext_libs = hdf5_fortran hdf5 z`) rather than discovered under a prefix. If
your cluster ships HDF5 as a system module and the normal mode fails to link,
copy that mode in `fobos.d/prism.fobos` and adapt `ext_libs`.
:::

## Step 4 — Run the CPU suite

```bash
./src/tests/prism/regression/run.sh cpu
```

It will:

1. run the GPU race-shape lint gate (`src/tests/lint/check-gpu-race-shapes.sh`),
   which fails the sweep if the forbidden strided device section
   `dxyz_gpu(b,1:3)` reappears in any `.F90`;
2. create a private Python venv at `exe/.regression-venv/` on first run and
   install `h5py` + `numpy` (needed by `digest.py`; WSL2 and CI system Python
   are externally managed, hence the venv);
3. `fobis build --mode prism-cpu-gnu`;
4. auto-discover every case (any immediate subdirectory with an `input.ini`),
   run it under `mpirun -np 2`, and compare the digest and residuals against
   `golden/cpu/`.

```bash
./src/tests/prism/regression/run.sh cpu --no-build            # reuse exe/
./src/tests/prism/regression/run.sh cpu --varset mymachine_gnu
```

Exit code 0 means every case passed:

```
== Summary [cpu]
PASS: 5
FAIL: 0
SKIP: 4
```

::: warning The venv lives under `exe/`, not `src/`
Deliberately. `fobis` (build) and the API-doc generator both scan `src/` and
would otherwise compile or document the Fortran test fixtures shipped inside
`numpy`/`h5py` — breaking the link step and generating junk doc pages.
:::

### Expected skips (not failures)

`SKIP` is normal. Cases without a `golden/<backend>/` directory are not anchors
and are skipped by design:

- `rmf-fwl` — FNL-goldened only; its CPU golden is pending a CI capture.
- `rmf-2realm-fd-pulse` — no golden and no `check.sh`; a manual reproducer only.
- `rmf-amr`, `rmf-amr-fd`, `rmf-amr-fd-pulse` — goldenless by design, driven by
  their own `check.sh` oracle (Step 6), not by `run.sh`.

## Step 5 — Run the FNL (GPU) suite

GPU runs are **local only** — GitHub-hosted runners have no GPU, so there is no
FNL CI job. On the WSL2 dev box:

```bash
./src/tests/prism/regression/run-fnl-local.sh
```

That wrapper loads `nvhpc`, exports the WSL2 MPI/UCX workarounds, then execs
`run.sh fnl --varset local_nvf`.

**On any other machine, do not use the wrapper** — its environment is
WSL2-specific. Set up your own and call `run.sh` directly:

```bash
module load nvhpc                                   # or your site's modules
./src/tests/prism/regression/run.sh fnl --varset mymachine_nvf
```

### WSL2 vs a real cluster

The wrapper's env vars are WSL2 workarounds and are a **performance regression
or a no-op on real InfiniBand + GPUDirect RDMA**. Never copy them into a Slurm
job script.

| Variable | WSL2 value | Why | On a real cluster |
|---|---|---|---|
| `UCX_RNDV_THRESH` | `inf` | forces every message EAGER; UCX rendezvous device transports cannot get the GPU primary context through WSL's `/dev/dxg` shim → SIGABRT ([#12](https://github.com/szaghi/adam/issues/12)) | **Do not set** — rendezvous *is* the fast path |
| `UCX_TLS` | `^cma` | WSL's cross-memory-attach shmem transport is broken (`ptrace_scope`) | prefer `rc_x,sm,cuda_copy,cuda_ipc` |
| `OMPI_MCA_coll_hcoll_enable` | `0` | hcoll needs Mellanox IB hardware this box lacks | leave enabled if you have IB |
| `OMPI_MCA_pml` | `ucx` | select the UCX point-to-point layer | usually correct as-is |

The FNL backend posts MPI directly on device-resident ghost buffers (GPU-direct
by design — there is **no host-staging fallback in the code**), so GPU-aware MPI
is a hard requirement, not an optimization.

`run.sh` uses `mpirun -np 2` for both backends. On a multi-GPU box that is one
rank per GPU; the grids are small enough that two ranks can share one GPU (the
hardware floor is one GPU, not two).

## Step 6 — Run the AMR-seam `check.sh` oracles

Three cases are goldenless by design and assert a physics oracle instead of a
committed digest. `run.sh` does **not** invoke them:

```bash
cd src/tests/prism/regression
./rmf-amr/check.sh            --build     # structural: refinement/registration/reflux
./rmf-amr-fd/check.sh         --build     # seam div(B) baseline + div(J) band
./rmf-amr-fd-pulse/check.sh   --build     # source-free two-invariant oracle
./rmf-amr-fd-pulse/check.sh --convergence # + the 16->32 refinement ladder
```

`--build` builds `prism-cpu-gnu` first. To check another backend's executable,
override `PRISM_EXE` — the caller owns the matching environment:

```bash
PRISM_EXE=$PWD/../../../../exe/adam_prism_fnl ./rmf-amr-fd/check.sh
```

::: warning Baselines are CPU-pinned
Acceptance thresholds are pinned numeric baselines with ±5% (div(B)) or factor-3
(div(J)) bands, calibrated for cross-compiler noise; the `p_obs` convergence
assertion is compiler-independent. A *small* miss on a new machine is a
tolerance question; an orders-of-magnitude miss is a real defect. Never widen a
band to make a new machine pass without understanding why.
:::

## Step 7 — Goldens on a new machine

A new machine does **not** mint goldens. The authority rules are fixed:

- **CPU golden authority is CI** (Ubuntu, GCC 14). A golden captured on a dev
  box with a different GCC can fail in CI with no real change.
- **FNL golden authority is the GPU workstation** — there is no GPU CI.

So the correct outcome on a newly-ported machine is: **the committed goldens
pass within tolerance**. That *is* the port validation. The tolerances
(`rtol=1e-6, atol=1e-3`) are calibrated precisely for the cross-machine,
cross-compiler case.

If a case fails, the question is "what did this machine do differently", not
"let me refresh the golden". Refreshing an existing golden is a deliberate,
reviewer-approved act tied to an intentional numerical change — see the suite
README.

To capture a golden for a case that has none on your backend:

```bash
REGRESSION_RUN_GOLDENLESS=1 ./src/tests/prism/regression/run.sh fnl --varset mymachine_nvf
mkdir -p src/tests/prism/regression/<case>/golden/fnl
cp src/tests/prism/regression/<case>/work-fnl/digest.txt \
   src/tests/prism/regression/<case>/work-fnl/*-residuals.dat \
   src/tests/prism/regression/<case>/golden/fnl/
```

## Troubleshooting by layer

| Symptom | Layer | Cause / fix |
|---|---|---|
| `FileNotFoundError: ... 'mpif90'` | 1 | MPI module not loaded before `fobis` |
| `-gpu=$NVF_CC` unexpanded; nvfortran dies immediately | 2 | built an NVF mode without `--varset <nvf varset>` |
| HDF5 link errors (`h5open_f` undefined) | 1–2 | `$HDF5_PREFIX` wrong, or HDF5 built with a different MPI |
| Run hangs at first I/O | 1 | serial HDF5, or HDF5's MPI ≠ the app's MPI |
| `ERROR: mpirun not on PATH` | 1 | load your MPI module before `run.sh` |
| `ERROR: executable exe/adam_prism_* not found after build` | 3 | the build failed — scroll up; don't trust `--no-build` |
| venv creation fails / `h5py` import error | 4 | system Python lacks `venv`; delete `exe/.regression-venv` and retry |
| `error stop : failed to load [fWLayer].(width)` | 4 | the known suite breakage — see the banner at the top |
| `no '*-*.h5' checkpoints produced` | 4 | the solver aborted — read the run output; try a debug build |
| SIGABRT in `ucp_proto_rndv_send_start` (FNL, np≥2) | 4 | WSL2 rendezvous bug — `UCX_RNDV_THRESH=inf` ([#12](https://github.com/szaghi/adam/issues/12)) |
| Digest mismatch ~1e-9 on stable fields | — | expected cross-compiler noise, inside `rtol=1e-6`; not a failure |
| Digest mismatch orders of magnitude | — | a real regression; do **not** refresh the golden |

### Debug builds

`run.sh` builds the **release** modes, because release is what ships. When a
release run crashes with an unhelpful backtrace, rebuild debug and run the case
by hand:

```bash
fobis build --mode prism-fnl-nvf-debug --varset local_nvf
cd src/tests/prism/regression/<case>/work-fnl
mpirun -np 1 ../../../../../../exe/adam_prism_fnl input.ini
```

Debug adds `-Mbounds -Mchkptr -Ktrap=fp -traceback -O0` (nvfortran) or
`-fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace` (gfortran). Two
caveats:

- **Stale objects.** FoBiS does not always recompile a `.F90` when an
  `#include`d `.INC` it depends on changes (notably the vendored FUNDAL
  templates). If a rebuild seems not to pick up an edit, delete the stale
  objects under `exe/obj/nvf/debug/` (and `exe/mod/...`) and rebuild.
- **`-fast`-only bugs.** Some failures appear only in release — the optimiser
  exposes undefined behaviour that `-O0` steps past. "Debug passes" does not
  mean "fixed".

A one-rank run will **not** match the committed two-rank golden; `-np 1` is for
isolating a crash, not for checking.

## Porting checklist

```
[ ] fobis 3.8+ installed; `fobis fetch` run on a fresh clone
[ ] compiler + MPI modules loaded (BEFORE building HDF5 and BEFORE fobis)
[ ] parallel HDF5 available, built against THAT MPI
[ ] [varset:<machine>] added to fobos with $HDF5_PREFIX (+ $NVF_CC for GPU)
[ ] $NVF_CC matches `nvidia-smi --query-gpu=compute_cap`
[ ] fobis build --mode prism-cpu-gnu --varset <machine> succeeds
[ ] run.sh cpu --varset <machine> reaches the run phase
[ ] (GPU) GPU-aware MPI confirmed; run.sh fnl --varset <machine>_nvf
[ ] (GPU) WSL2 UCX vars NOT copied into the cluster job script
[ ] check.sh oracles run and pass within their pinned bands
```

## See also

- [PRISM regression suite](./prism-regression) — what the cases assert and why.
- Suite [`README.md`](https://github.com/szaghi/adam/blob/develop/src/tests/prism/regression/README.md) — golden capture, tolerance rationale, adding a case.
- [Forest (multi-realm)](../guide/forest) — the seam machinery the cases exercise.
