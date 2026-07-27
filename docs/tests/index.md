# Tests

ADAM uses a multi-level testing strategy to ensure correctness across backends and configurations.

## Unit Tests

Finite difference operator tests with trigonometric analytical solutions:

```bash
# Build a specific test
fobis build --mode test-fdv-gradient-trigonometric-gnu
fobis build --mode test-fdv-divergence-trigonometric-gnu
fobis build --mode test-fdv-curl-trigonometric-gnu
fobis build --mode test-fdv-laplacian-trigonometric-gnu
fobis build --mode test-fdv-operators-trigonometric-gnu
fobis build --mode test-fdv-operators-step-gnu

# Run with MPI
mpirun -np <N> exe/<test_executable>
```

## Regression Tests

The [PRISM regression suite](./prism-regression) is the structural-change baseline for the Maxwell/EM solver: CI-tuned two-rank cases run on both the CPU and FNL (GPU) backends, compared tolerance-aware against committed field-digest and residuals goldens. It covers the forest multi-realm seam machinery (α/β cadence, inter-realm 1:1 mirror seams), the intra-realm 2:1 AMR seams (with `check.sh`-driven `div(B)` oracles), reflux, and the fWLayer. See the [Forest guide](../guide/forest) for the seam machinery itself.

## Integration Tests

NASTO integration tests with known CFD benchmarks:

| Test | Description |
|------|-------------|
| [Sod-X](./sod-x) | 1D Riemann Problem of Sod along X axis |
| [Sod-Y](./sod-y) | 1D Riemann Problem of Sod along Y axis |
| [Sod-Z](./sod-z) | 1D Riemann Problem of Sod along Z axis |
| [Shock-Sphere](./shock-sphere) | Shock-sphere interaction simulation |

## Compiler Testing Matrix

All changes affecting core numerics or parallel logic must pass:

- **GNU gfortran** — baseline CPU, strict bounds checking
- **NVIDIA nvfortran** — GPU targets (CUDA Fortran, OpenACC)
- **Intel ifort/ifx** — aggressive optimization, vectorization reports

Use `-fbounds-check -fcheck=all` (GNU) or `-check all -traceback` (Intel) during development.
