# Tests

ADAM uses a multi-level testing strategy to ensure correctness across backends and configurations.

## Unit Tests

Finite difference operator tests with trigonometric analytical solutions:

```bash
# Build a specific test
FoBiS.py build -mode test-fdv-gradient-trigonometric-gnu
FoBiS.py build -mode test-fdv-divergence-trigonometric-gnu
FoBiS.py build -mode test-fdv-curl-trigonometric-gnu
FoBiS.py build -mode test-fdv-laplacian-trigonometric-gnu
FoBiS.py build -mode test-fdv-operators-trigonometric-gnu
FoBiS.py build -mode test-fdv-operators-step-gnu

# Run with MPI
mpirun -np <N> exe/<test_executable>
```

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
