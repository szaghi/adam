# GPU (FNL)

Recipes for building and running the GPU backend, FNL, and for writing device code that follows its rules. FNL is
OpenACC with nvfortran on NVIDIA GPUs, or OpenMP offloading with amdflang on AMD GPUs; the device memory is managed by
FUNDAL. The library side is described in [library on FNL](/concepts/library-fnl).

## GPU1. Build an FNL mode

**Problem.** Compile FLUME, PRISM or RIPPLE for the GPU.

**Commands** (modes from the [build-modes page](/reference/build-modes)):

```bash
fobis build --mode flume-fnl-nvf --varset local_nvf       # FLUME, OpenACC, nvfortran
fobis build --mode flume-fnl-nvf-fft --varset local_nvf   # FLUME with cuFFT (self-gravity)
fobis build --mode prism-fnl-nvf --varset local_nvf       # PRISM, OpenACC
fobis build --mode ripple-fnl-nvf --varset local_nvf      # RIPPLE, OpenACC
fobis build --mode flume-fnl-omp-amd --varset <an AMD varset>   # FLUME, OpenMP offload, amdflang
```

**Rules.**
- Every NVF mode needs an NVF varset: it supplies `$NVF_CC` (the `-gpu=ccXX` compute capability) and the HDF5 built by
  the same nvfortran. Without one, the build stops at the first nvfortran call on an unexpanded `-gpu=$NVF_CC`.
  `local_nvf` is the workstation's; add a `[varset:*]` block in `fobos` for a new machine
  ([install](/start/install#varsets-one-per-machine-and-toolchain)).
- Build clean after switching compilers: `.mod` files work only with the compiler that wrote them.
- The macros the modes set: `_FNL` (the FNL backend), `_MPI_`, and `DEV_OAC` + `COMPILER_NVF` (OpenACC on nvfortran) or
  `DEV_OMP` + `DEV_HIP` (OpenMP offload on amdflang). Code tests them only to select the backend.

**Reference.** [Install](/start/install), [build modes](/reference/build-modes).

## GPU2. Allocate device memory and check it

**Problem.** Allocate a device array so that a failure stops the run with a message instead of crashing later.

**Code** (`src/app/ripple/fnl/adam_ripple_fnl_object.F90`):

```fortran
call dev_alloc(fptr_dev=self%q_gpu, lbounds=[1,1-ngc,1-ngc,1-ngc,1], ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,1], &
               init_value=0._R8P, ierr=ierr) ; call check('q_gpu')
...
subroutine check(name)
character(*), intent(in) :: name
if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate '//name//' in ripple_fnl_object')
endsubroutine check
```

**What you get.** On a full device the run stops at initialization with `failed to allocate <name>`. FUNDAL's code for
a pointer it could not allocate is `101` (`FUNDAL_ERR_FPTR_DEV_NOT_ALLOCATED`); in practice it means the device memory
is exhausted, by this run or by another process on a shared GPU.

**The budget.** The number of blocks a rank can hold, `nb`, comes from the memory it is given:
$n_b = 0.85\,M / (F\,B)$, with $M$ the memory, $B = 8\,(n_i+2n_g)(n_j+2n_g)(n_k+2n_g)$ the bytes of **one scalar
variable** on one block (ghosts included, `grid%block_weight`) and $F$ = `fields_number` the number of such **scalar**
arrays the application allocates per block: FLUME counts `nv*(5 + stages) + nv_aux + ...` (state, residuals, stages
and fluxes, each `nv` wide). On FNL $M$ is the device
memory divided by the ranks that share the device (`count_device_sharers`) and, in a forest, by the realms. An
application that under-counts `fields_number` over-commits the device.

**Reference.** `adam_object%initialize` (`src/lib/common/adam_adam_object.F90`); FUNDAL's `dev_alloc`.

## GPU3. Copy between host and device without a hidden temporary

**Problem.** Copy a host array to its transposed device twin.

**Code** (`src/app/flume/fnl/adam_flume_fnl_object.F90`):

```fortran
allocate(self%buf_5D_R8P(1:nb,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nv))   ! the extent of q_gpu
...
if (any(shape(self%buf_5D_R8P) /= shape(self%q_gpu))) &
   call mpih_fnl%error_stop(msg=': flume_fnl_object%copy_cpu_gpu: buf_5D_R8P does not match q_gpu (issue #63)')
call dev_memcpy_to_device(bb=self%db5, ij=[1,5], tb=self%hb5, dst=self%q_gpu, src=self%q, buf=self%buf_5D_R8P)
```

**The rule.** Pass `buf=` (FUNDAL's buffered transposed copy) **only with a staging buffer of exactly the destination's
shape**. With a different extent the device destination reaches an lbound-remapped dummy, nvfortran materialises a host
copy-in temporary, and the copy fails with `CUDA_ERROR_INVALID_VALUE` (issue
[#31](https://github.com/szaghi/adam/issues/31)); otherwise use the whole-array `dev_assign_to_device(..., ij=[1,5])`.
The same holds for any host buffer that mirrors a device one: give it the device buffer's extent.

**Reference.** [Library on FNL](/concepts/library-fnl).

## GPU4. Write a device kernel for both GPU vendors

**Problem.** One kernel source for OpenACC (NVIDIA) and OpenMP offload (AMD).

**Code** (`src/app/flume/fnl/adam_flume_fnl_kernels.F90`): every offload `!$acc` directive is followed by its `!$omp`
twin, and the macros of `fundal.H` pick the clause spelling per backend:

```fortran
!$acc parallel loop independent gang vector DEVICEVAR(delta_gpu,q_gpu) &
!$acc& firstprivate(axis,sgn,b,ni,nj,nk,nv,nface_cells,scale,limit,mag_w,psi_w,cs2_inert) private(i,j,k,qc,dc,theta) &
!$acc& reduction(+:nclip)
!$omp OMPLOOP DEVICEPTR(delta_gpu,q_gpu) &
!$omp& firstprivate(axis,sgn,b,ni,nj,nk,nv,nface_cells,scale,limit,mag_w,psi_w,cs2_inert) private(i,j,k,qc,dc,theta) &
!$omp& reduction(+:nclip)
do c=1, nface_cells
```

`DEVICEVAR` is `deviceptr` on nvfortran; `OMPLOOP` is `target teams distribute parallel do` under `DEV_OMP`. The device
arrays are transposed, `(nb, ni, nj, nk, nv)`, with the block index stride-1. RIPPLE's kernels carry only the `!$acc`
directives: it is NVIDIA-only by decision, and FLUME is the example of the twin pattern.

**Rules.** No configuration branching inside a kernel (one kernel per model, the `select case` in the host
controller); a reduction goes into a procedure-local variable.

**Reference.** [Library on FNL](/concepts/library-fnl).

## GPU5. Run on a development box with WSL2, not on a cluster

**Problem.** Run an FNL executable with several ranks on a WSL2 workstation.

**Settings** (`src/tests/prism/regression/run-fnl-local.sh`, WSL2 only):

```bash
export OMPI_MCA_pml="ucx"
export UCX_TLS="^cma"
export OMPI_MCA_coll_hcoll_enable="0"
export UCX_RNDV_THRESH="inf"
```

**Why.** The FNL ghost exchange posts MPI directly on device buffers (GPU-direct; the code has no host-staging
fallback). With two or more ranks UCX moves them by rendezvous, whose device transports cannot reach the GPU context
through WSL's `/dev/dxg`, and the run aborts at the first exchange ([#12](https://github.com/szaghi/adam/issues/12));
`UCX_RNDV_THRESH=inf` forces every message eager. CMA is restricted under WSL and hcoll needs the InfiniBand hardware the box lacks; the `nvhpc` modulefile of the box
also sets `UCX_MEMTYPE_CACHE=n`, as WSL's memory-type detection is unreliable.

::: danger Never copy these settings into a cluster job
On a cluster with GPUDirect RDMA, rendezvous **is** the fast path: `UCX_RNDV_THRESH=inf` is a performance regression
there, and the others are no-ops or worse. WSL runs check correctness only; never quote their timings.
:::

**Reference.** [Clusters](/cookbook/clusters);
[library on FNL](/concepts/library-fnl#gpu-direct-ghost-exchange-and-the-wsl2-rendezvous-abort-issue-12).
