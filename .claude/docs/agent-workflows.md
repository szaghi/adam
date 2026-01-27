# Agent Workflows and Debugging

Structured protocols for common tasks and emergency debugging.

## Task: GPU Kernel Optimization

**Input**: User reports kernel is underperforming (low occupancy, low utilization)

**Response Protocol**:
1. Request profiling data: occupancy, achieved memory bandwidth, warp execution efficiency
2. Analyze memory access pattern: coalesced vs. strided, transactions per request
3. Check for: bank conflicts (shared memory), atomic contention, excessive register usage
4. Propose optimization with quantitative prediction: "Expected to improve bandwidth from X to Y GB/s"
5. Suggest verification: rerun profiler, compare before/after metrics

**Example**: If atomic operations detected -> suggest graph coloring (as used for vertex-based AMR operations in this codebase)

## Task: MPI Load Balancing

**Input**: User observes poor strong scaling or high MPI wait time

**Response Protocol**:
1. Request timing breakdown: computation vs. communication time per rank
2. Check domain decomposition: are blocks evenly distributed?
3. For AMR: analyze refinement pattern - clustered refinement causes imbalance
4. Propose solutions:
   - Dynamic load balancing (redistribute blocks between MPI ranks)
   - Space-filling curve (Morton/Hilbert) reordering for locality
   - Overdecomposition with task-based parallelism (MPI+X)

## Task: Fortran Modernization

**Input**: User has legacy Fortran 77/90 code to modernize

**Response Protocol**:
1. Identify obsolescent features: `common` blocks, `equivalence`, assumed-size arrays, `goto`
2. Propose replacements: `module` data, derived types, assumed-shape arrays, structured control flow
3. Introduce `submodule` for large modules (reduces compilation dependencies)
4. Convert to assumed-shape arrays: `real :: a(:,:,:)` instead of `real :: a(nx,ny,nz)`
5. Use `iso_fortran_env` for `real64`, `int32` instead of `selected_real_kind`

## Task: Debugging Numerical Instability

**Input**: Simulation crashes or produces non-physical results

**Response Protocol**:
1. Check CFL condition: `dt < CFL * min(dx, dy, dz) / max_wavespeed`
2. Verify boundary conditions: wall, symmetry, inflow/outflow implementations
3. Check for NaN/Inf: enable floating-point exception trapping (`-ffpe-trap=invalid,zero,overflow`)
4. Review discretization: upwind/central bias, dissipation terms, limiter functions
5. Suggest diagnostic output: write fields at crash point, check conservation errors

## Task: Adding a New Feature

1. Understand existing architecture (which backend? CPU/GPU/both?)
2. Propose design: module structure, type extensions, interfaces
3. Identify test cases for verification
4. Consider multi-compiler compatibility
5. Document with references (papers, standards)

---

## Emergency Debugging Checklist

### Compilation Fails
- [ ] Check module dependencies (`FoBiS.py build -lmodes` to verify mode exists)
- [ ] Verify HDF5/MPI library paths in fobos configuration
- [ ] Load correct module environment (compiler, MPI, HDF5 versions must match)
- [ ] Check preprocessor macros (only one of `_NVF`, `_FNL`, `_GMP` should be defined)

### Runtime Crash
- [ ] Enable bounds checking: `-fbounds-check` (GNU), `-check bounds` (Intel)
- [ ] Enable FP exceptions: `-ffpe-trap=invalid,zero,overflow` (GNU)
- [ ] Run with debugger: `gdb --args mpirun -np 4 ./exe/nasto`
- [ ] Check stack size: `ulimit -s unlimited` (Fortran recursive calls)
- [ ] Verify MPI library: `ldd exe/nasto | grep mpi` (should match loaded module)

### GPU Errors
- [ ] CUDA error codes: check return values, run with `cuda-memcheck`
- [ ] OpenACC: set `ACC_SYNCHRONOUS=1` for immediate error reporting
- [ ] Device memory: verify sufficient GPU memory (`nvidia-smi`)
- [ ] Data residency: check `!$acc update device/host` for stale data

### Performance Regression
- [ ] Profile before/after: compare wall time, memory bandwidth, kernel launch overhead
- [ ] Check compiler optimization flags: should be `-O3` or equivalent
- [ ] Verify no debug flags in production build (`-g`, `-fbounds-check`)
- [ ] Test on same hardware (CPU model, GPU model, memory size)
- [ ] Review recent changes: `git diff <last-good-commit>`

### MPI Hangs
- [ ] Deadlock detection: check for unmatched Send/Recv pairs
- [ ] Barrier analysis: ensure all ranks reach collective operations
- [ ] Set `MPICH_ASYNC_PROGRESS=1` or `OMPI_MCA_opal_progress_threads=1`
- [ ] Reduce process count to isolate issue (e.g., test with 2 ranks)
- [ ] Enable MPI debug: `export MPICH_DEBUG=1` or `--mca btl_base_verbose 30`
