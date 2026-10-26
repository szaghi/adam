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

The [PRISM regression suite](./prism-regression) is the structural-change baseline for the Maxwell/EM solver: CI-tuned two-rank cases run on both the CPU and FNL (GPU) backends, compared tolerance-aware against committed field-digest and residuals goldens. It covers the forest multi-realm seam machinery (α/β cadence, inter-realm 1:1 mirror seams), the intra-realm 2:1 AMR seams (with `check.sh`-driven `div(B)` oracles), reflux, and the fWLayer. See the [Forest guide](/concepts/forest) for the seam machinery itself.

To stand the suite up on a new machine, compiler, or GPU architecture, follow the [regression bring-up tutorial](./prism-regression-tutorial).

::: tip Current status (`1c6ebee5`, medooza gates)
PRISM: **15 PASS, 9 SKIP** on the CPU backend and **16 PASS, 8 SKIP** on FNL; a SKIP is a case with no golden for that
backend, or one of the two non-case directories `reproducers/` and `slices/`. FLUME: **54/54 PASS** on both backends.
:::

The [FLUME regression suite](https://github.com/szaghi/adam/blob/develop/src/tests/flume/regression/run.sh)
(`src/tests/flume/regression/`, 54 cases) is the same harness for FLUME: each case is a frozen copy of a verification
input, compared against committed field digests, residual logs and conservation histories on each backend. The
physics is asserted by the verification oracles (`src/tests/flume/verification/`, the
[verification gallery](../applications/flume/verification/)); this suite detects changes.

```bash
src/tests/flume/regression/run.sh cpu                      # builds flume-cpu-gnu and flume-cpu-gnu-fft
src/tests/flume/regression/run.sh fnl --varset local_nvf   # builds flume-fnl-nvf and flume-fnl-nvf-fft
REGRESSION_CASES="sod-x brio-wu" src/tests/flume/regression/run.sh cpu --no-build
```

Self-gravity cases run on the `-fft` executable (FFTW on the CPU, cuFFT on FNL).

Research and development cases live separately under `src/tests/prism/<backend>/{fd,fv}/<family>/` — long integration times, full AMR, sized for physics validation. They are **not** regression anchors and `run.sh` never descends into them.

## Compiler Testing Matrix

Changes to core numerics or parallel logic are gated on both maintained toolchains (versions in
[Install](/start/install)):

- **GNU gfortran** — the CPU backend; the `*-debug` modes add bounds checking and FPE traps
- **NVIDIA nvfortran** — the FNL backend (OpenACC)
