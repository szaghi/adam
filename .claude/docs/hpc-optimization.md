# HPC Optimization Guide

Performance-driven development workflow for GPU, MPI, and memory optimization.

## Profile Before Action

- **Always ask first**: "Do you have profiling data?" (Nsight Compute, nvprof, gprof, Intel VTune)
- Identify bottlenecks: memory bandwidth, compute intensity, communication overhead, load imbalance
- Never optimize without quantitative evidence of the problem

## Correctness Verification Strategy

- Propose changes with explicit verification plan
- Specify test cases: unit tests, integration tests (Riemann problems), regression tests
- Multi-compiler validation: gfortran, Intel ifort/ifx, NVIDIA nvfortran
- Enable runtime checks during development: `-fbounds-check`, `-fcheck=all` (GNU), `-check all` (Intel)
- Memory debugging: Valgrind (CPU), cuda-memcheck (GPU)

## Performance Optimization Patterns

### GPU-Specific (OpenACC, CUDA Fortran)
- **Data movement**: Minimize host-device transfers, maximize data reuse on device
- **Loop directives**: Explicitly specify `gang`, `vector`, `seq` collapse depth
- **Atomic operations**: Red flag for serialization - consider alternatives (graph coloring, privatization, reduction)
- **Memory coalescing**: Ensure contiguous access patterns in innermost loops
- **Kernel fusion**: Combine multiple kernels to reduce launch overhead and improve data locality

### MPI Best Practices
- **Communication/computation overlap**: Use non-blocking MPI (Isend/Irecv, Iallreduce)
- **Domain decomposition**: Check load balancing across ranks - propose block-splitting if needed
- **Collective optimization**: Use MPI-3 neighborhood collectives for AMR stencil operations
- **One-sided communication**: Consider RMA (Put/Get) for ghost cell exchange

### Memory Layout (Fortran-specific)
- **Column-major awareness**: Inner loops on first index for cache efficiency
- **Array contiguity**: Use `contiguous` attribute, avoid unnecessary copies
- **Derived types**: Watch for alignment, padding, and GPU transfer efficiency
- **Avoid temporary arrays**: Use in-place operations or explicit buffers

## Benchmarking Requirements

For any performance-related change, specify:
- **Metrics**: Wall time, speedup, scaling efficiency, memory bandwidth utilization
- **Weak scaling**: Fixed problem size per process/GPU
- **Strong scaling**: Fixed total problem size, varying parallelism
- **Target**: Justify "acceptable performance" (e.g., >80% efficiency to 1024 GPUs)

## Documentation of Changes

- Comment optimizations with: rationale, expected impact, profiling evidence
- Reference specific sections of standards/papers when using advanced features
- Note compiler-specific behavior (especially GCC vs Intel vs NVIDIA differences)

## Common Pitfalls

1. **Atomic Operations on GPUs**:
   - Cause severe serialization (observed: 10-100x slowdown in vertex-based operations)
   - Solution: Graph coloring, data reordering, or reduction strategies

2. **Array Reallocation with Derived Types**:
   - GCC gfortran bounds checking differs from Intel/NVIDIA (observed: silent corruption)
   - Always verify bounds after `allocate`/`deallocate` with multiple compilers

3. **OpenACC Data Directives**:
   - Unstructured data regions (`enter data`/`exit data`) require careful lifecycle management
   - Missing `update device/host` causes stale data bugs (hard to debug)
   - Use `present` clause to verify data residency assumptions

4. **MPI Datatypes**:
   - Derived type MPI commits can fail silently - always check error codes
   - Padding/alignment issues when transferring structs across heterogeneous nodes

5. **Floating-Point Reproducibility**:
   - Non-associative reductions in MPI may cause answer changes with different process counts
   - GPU atomics introduce non-determinism - document or eliminate

## Performance Anti-Patterns

**Do NOT**:
- Optimize without profiling data
- Use `!$acc kernels` without loop directives (compiler may not parallelize optimally)
- Ignore compiler warnings (especially `-Wunused-variable`, `-Wimplicit-interface`)
- Mix MPI library versions or HDF5 versions across build
- Assume performance portability (GPU code often needs per-architecture tuning)

**DO**:
- Start with `!$acc parallel loop` with explicit clauses
- Use compiler feedback: `-Minfo=accel` (NVIDIA), `-qopt-report` (Intel)
- Validate correctness with sanitizers: AddressSanitizer (CPU), cuda-memcheck (GPU)
- Profile with multiple problem sizes (small, medium, large)
- Test on target architecture before claiming "optimization complete"
