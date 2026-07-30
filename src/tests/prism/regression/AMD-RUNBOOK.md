# PRISM regression on Thera (AMD) — interactive compute-node runbook

Step-by-step for building and running the PRISM regression suite with the AMD
`amdflang` toolchain on the **Thera** cluster: the CPU backend (`amdflang`) and
the ROCm **OpenMP target-offload GPU** backend.

Follow this after logging into `TheraS02` (the login node). The build works on
the login node, but **the run must happen on a GPU compute node**: the spack
OpenMPI has hard `librdmacm.so.1` / `libibverbs.so.1` (OFED) dependencies that
only exist on compute nodes — on the login node the executable dies with
`librdmacm.so.1: cannot open shared object file`.

---

## 0. The toolchain (must stay self-consistent)

Fortran `.mod` files are compiler-specific and HDF5 must link the same MPI the
app links, so these three MUST agree — they are the exact set the `thera_amd`
HDF5 in `fobos` was built with:

| Piece | Path | Notes |
|---|---|---|
| `amdflang`/`amdclang` | `/home/mbycklin/amd/therock-23.2.1-gfx94X-7.13.0-7357b5084b` | AFAR flang 23.0.0, gfx94X |
| OpenMPI 5.0.10 | `/home/mbycklin/code/spack/install/linux-zen4/openmpi-5.0.10-3rvtw7bdkzbm3tkapamurnizsr5v7hjn` | `mpif90` → the amdflang above |
| parallel HDF5 1.14.6 (Fortran ON) | `fobos` varset `thera_amd` | built with that amdflang + OpenMPI |
| ROCm | `/opt/rocm-6.4.1` | provides `libamdhip64` for the offload link |

The `fobos` build modes and varset already exist:

- `prism-cpu-amd`      — `amdflang` CPU build
- `prism-fnl-omp-amd`  — `amdflang` + `-fopenmp --offload-arch=$GPU_ARCH -lamdhip64`
  (preproc `-D_FNL -D_MPI_ -DDEV_OMP -DDEV_HIP`)
- varset `thera_amd`   — supplies `$HDF5_PREFIX`
- `$GPU_ARCH`          — `fobos [common-variables]`, default **`gfx942`** (MI300)

---

## 1. Reserve an interactive GPU node

Pick the partition whose GPU matches your `$GPU_ARCH`:

| Partition | GPU | `--offload-arch` / `$GPU_ARCH` | Nodes |
|---|---|---|---|
| **`MI300x`** | MI300X | **`gfx942`** (default; matches loaded drop) | TheraC[16-19] |
| `MI350x` | MI350X | `gfx950` | TheraC[70-77] |
| `MI355x` | MI355X | `gfx950` | TheraC[78-93] |
| `MI250` | MI250 | `gfx90a` | TheraC[53-67] |
| `MI210` | MI210 | `gfx90a` | TheraC[12-13,15] |

The loaded AFAR drop is **gfx94X**, so the matching partition is **`MI300x`
(gfx942)** — use that unless you switch drops.

```bash
export PATH=/share/opt/slurm/current/bin:$PATH      # if slurm not already on PATH

# Interactive allocation on one MI300X node (2 GPUs is enough; run.sh uses -np 2)
salloc -p MI300x -N 1 --gres=gpu:2 --cpus-per-task=16 -t 02:00:00
# (add '--account=<proj>' if the site requires it)
```

`salloc` drops you into a shell on the allocation. If it leaves you on the login
node instead, hop onto the compute node with:

```bash
srun --pty --gres=gpu:2 bash -i
```

Confirm you are on the GPU node and see the devices:

```bash
hostname                 # -> TheraC16..19
/opt/rocm-6.4.1/bin/rocminfo | grep -m2 'Name:.*gfx'   # -> gfx942
/opt/rocm-6.4.1/bin/rocm-smi
```

---

## 2. Set up the AMD toolchain environment

Paste this once on the compute node (it is exactly what `run-amd-local.sh`
does). These are shell env vars; `$HDF5_PREFIX`/`$GPU_ARCH` are **fobos**
variables and come from the varset, not from here.

```bash
AMD_OMPI=/home/mbycklin/code/spack/install/linux-zen4/openmpi-5.0.10-3rvtw7bdkzbm3tkapamurnizsr5v7hjn
AMD_FLANG=/home/mbycklin/amd/therock-23.2.1-gfx94X-7.13.0-7357b5084b
ROCM=/opt/rocm-6.4.1

export PATH="$AMD_OMPI/bin:$AMD_FLANG/bin:$ROCM/bin:$PATH"
export LD_LIBRARY_PATH="$AMD_OMPI/lib:$AMD_FLANG/lib:$ROCM/lib:${LD_LIBRARY_PATH:-}"

# Make the OpenMPI wrappers drive the AMD compilers (defensive; the wrapper
# defaults already point here).
export OMPI_FC="$AMD_FLANG/bin/amdflang"
export OMPI_CC="$AMD_FLANG/bin/amdclang"
export OMPI_CXX="$AMD_FLANG/bin/amdclang++"

# Sanity: mpif90 must resolve to amdflang, and mpirun must exist.
mpif90 -show | awk '{print "FC ->", $1}'
which mpirun
```

> On a compute node OFED is in the default loader path, so no `librdmacm`
> fallback is needed. (On the login node only, add
> `export LD_LIBRARY_PATH=/share/modules/mlnxofed/5.2-2.2.0.0/lib64:$LD_LIBRARY_PATH`.)

---

## 3. Compile

Run from the repo root (`/home/giarossi/Projects/adam`).

### 3a. CPU backend (amdflang)

```bash
cd /home/giarossi/Projects/adam
fobis build --mode prism-cpu-amd --varset thera_amd
# -> exe/adam_prism_cpu
```

### 3b. GPU backend — OpenMP target offload (the AMD-enabling target)

```bash
cd /home/giarossi/Projects/adam
fobis build --mode prism-fnl-omp-amd --varset thera_amd
# -> exe/adam_prism_fnl
```

What that mode expands to (from `fobos.d/prism.fobos` + `fobos.d/templates.fobos`):

```
compiler : amd            # FoBiS drives mpif90 (MPI=on), which wraps amdflang
preproc  : -D_FNL -D_MPI_ -DDEV_OMP -DDEV_HIP -w
cflags   : -cpp -c -O2 -fopenmp --offload-arch=gfx942
lflags   :         -O2 -fopenmp --offload-arch=gfx942 -lamdhip64 \
                   -Wl,-rpath,$HDF5_PREFIX/lib -Wl,-z,execstack
```

So the underlying per-file compile is effectively:

```bash
mpif90 -cpp -c -O2 -fopenmp --offload-arch=gfx942 \
       -D_FNL -D_MPI_ -DDEV_OMP -DDEV_HIP -w \
       -I$HDF5_PREFIX/include <file>.F90
```

**Different GPU?** Override `$GPU_ARCH` by adding its varset (and load the
matching AFAR drop): e.g. on MI350x/MI355x

```bash
fobis build --mode prism-fnl-omp-amd --varset "thera_amd gfx950"
```

**Debug build** (bounds/traps, `-O0`): append `-debug` to the mode, i.e.
`prism-cpu-amd-debug` or `prism-fnl-omp-amd-debug`.

> Stale-object caveat: FoBiS does not always recompile a `.F90` when an
> `#include`d FUNDAL `.INC` changes. If a rebuild seems to ignore an edit, wipe
> `exe/obj/amd/` (and `exe/mod/amd/`) and rebuild.

---

## 4. Run the regression suite

The harness `run.sh` learned two AMD backends; `run-amd-local.sh` wraps §2's env
and hands off with `--varset thera_amd`:

| Command | fobis mode | exe | diffed against |
|---|---|---|---|
| `./run-amd-local.sh amd` | `prism-cpu-amd` | `adam_prism_cpu` | `golden/cpu/` |
| `./run-amd-local.sh amd-omp` | `prism-fnl-omp-amd` | `adam_prism_fnl` | `golden/fnl/` |

```bash
cd /home/giarossi/Projects/adam/src/tests/prism/regression

# CPU (validates amdflang numerics vs the committed CPU golden)
./run-amd-local.sh amd

# GPU OpenMP offload (validates AMD offload vs the committed FNL golden)
./run-amd-local.sh amd-omp

# reuse an existing build:
./run-amd-local.sh amd-omp --no-build
```

If you prefer to skip the wrapper and drive `run.sh` directly (after §2's env):

```bash
./run.sh amd     --varset thera_amd
./run.sh amd-omp --varset thera_amd
```

Both AMD backends reuse an existing golden on purpose: they exercise the *same*
Fortran source as an already-anchored backend (CPU-gnu / FNL-nvf) and only swap
the compiler/offload runtime, so the committed golden is the correct
cross-compiler / cross-runtime reference. `digest.py` tolerances
(`rtol=1e-6, atol=1e-3`) are calibrated for exactly this.

Expected shape (per the bring-up tutorial): the seven `fd` cases are the ones to
judge; `rmf-amr` (fv path) and `rmf-fwl` (un-migrated input) are known issues.

### GPU binding note

`run.sh` launches `mpirun -np 2` with no explicit GPU binding, so both ranks land
on GPU 0 and share it — fine for these small correctness grids. To pin one rank
per GPU instead, run by hand:

```bash
cd <case>/work-amd-omp
mpirun -np 2 --map-by numa \
  bash -c 'export ROCR_VISIBLE_DEVICES=$OMPI_COMM_WORLD_LOCAL_RANK; exec ../../../../../../exe/adam_prism_fnl input.ini'
```

---

## 5. Quick reference (copy-paste, MI300x)

```bash
# on login node
export PATH=/share/opt/slurm/current/bin:$PATH
salloc -p MI300x -N 1 --gres=gpu:2 --cpus-per-task=16 -t 02:00:00
# (now on the compute node)
hostname; /opt/rocm-6.4.1/bin/rocminfo | grep -m1 'Name:.*gfx'

AMD_OMPI=/home/mbycklin/code/spack/install/linux-zen4/openmpi-5.0.10-3rvtw7bdkzbm3tkapamurnizsr5v7hjn
AMD_FLANG=/home/mbycklin/amd/therock-23.2.1-gfx94X-7.13.0-7357b5084b
ROCM=/opt/rocm-6.4.1
export PATH="$AMD_OMPI/bin:$AMD_FLANG/bin:$ROCM/bin:$PATH"
export LD_LIBRARY_PATH="$AMD_OMPI/lib:$AMD_FLANG/lib:$ROCM/lib:${LD_LIBRARY_PATH:-}"
export OMPI_FC=$AMD_FLANG/bin/amdflang OMPI_CC=$AMD_FLANG/bin/amdclang OMPI_CXX=$AMD_FLANG/bin/amdclang++

cd /home/giarossi/Projects/adam
fobis build --mode prism-fnl-omp-amd --varset thera_amd
cd src/tests/prism/regression
./run-amd-local.sh amd-omp --no-build
```

---

## Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `librdmacm.so.1: cannot open shared object file` | You are on the **login node** — run on a compute node (§1), or add the OFED fallback lib dir to `LD_LIBRARY_PATH`. |
| `mpif90` FC is not amdflang | `$AMD_OMPI/bin` not first on PATH, or a site `OMPI_FC` overrides — re-run §2. |
| HDF5 link errors (`h5open_f` undefined) / `.mod` version error | amdflang version ≠ the one that built the `thera_amd` HDF5 — use the drop in §0, don't mix AFAR drops. |
| `-gpu`/offload errors, `cannot find -lamdhip64` | ROCm not on `LD_LIBRARY_PATH`, or `--offload-arch` ≠ the node's GPU (`gfx942` on MI300x). |
| `error stop : failed to load [fWLayer].(width)` | expected for `rmf-fwl` (un-migrated input); not an AMD issue. |
| offload runs but is slow / 2 ranks on 1 GPU | expected default; pin per-rank GPUs (§4 GPU binding note). |
